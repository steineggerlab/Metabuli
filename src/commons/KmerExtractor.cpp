#include "KmerExtractor.h"
#include <unordered_map>

KmerExtractor::KmerExtractor(
    const LocalParameters &par,
    const GeneticCode & geneticCode,
    int kmerFormat) 
    : par(par), geneticCode(geneticCode), kmerScanners(new KmerScanner*[par.threads]) {
    // Initialize k-mer scanners for each thread
    kmerLen = 8;
    for (int i = 0; i < par.threads; ++i) {
        if (kmerFormat == 1) {
            if (i == 0) {
                std::cout << "Using OldMetamerScanner." << std::endl;
            }
            kmerScanners[i] = new OldMetamerScanner(geneticCode);
        } else if (kmerFormat == 2) {
            if (par.syncmer) {
                kmerScanners[i] = new SyncmerScanner(par.smerLen, geneticCode);
            } else {
                kmerScanners[i] = new MetamerScanner(geneticCode);
            }
        } else if (kmerFormat == 3) {
            kmerScanners[i] = new KmerScanner_dna2aa(geneticCode, 12);
            kmerLen = 12;
        } else if (kmerFormat == 4) {
            kmerScanners[i] = new KmerScanner_aa2aa(12);
        } else if (kmerFormat == 5) {
            kmerScanners[i] = new SyncmerScanner_dna2aa(geneticCode, 12, par.smerLen);
            kmerLen = 12;
        } else {
            std::cerr << "Error: Invalid k-mer format specified." << std::endl;
            exit(EXIT_FAILURE);
        }
    }
    spaceNum = 0;
    maskMode = par.maskMode;
    maskProb = par.maskProb;
    subMat = new NucleotideMatrix(par.scoringMatrixFile.values.nucleotide().c_str(), 1.0, 0.0);
    probMatrix = new ProbabilityMatrix(*(subMat));
}

KmerExtractor::~KmerExtractor() {
    delete probMatrix;
    delete subMat;
    for (int i = 0; i < par.threads; ++i) {
        delete kmerScanners[i];
    }
    delete[] kmerScanners;
}

void KmerExtractor::extractQueryKmers(Buffer<Kmer> &kmerBuffer,
                                      vector<Query> & queryList,
                                      const QuerySplit & currentSplit,
                                      const LocalParameters &par,
                                      KSeqWrapper* kseq1,
                                      KSeqWrapper* kseq2) {
    time_t beforeKmerExtraction = time(nullptr);
    std::cout << "Query k-mer extraction : " << flush;
    if (par.seqMode == 1 || par.seqMode == 3) { // Single-end short-read sequence or long-read sequence
        fillQueryKmerBufferParallel(kseq1,
                                    kmerBuffer,
                                    queryList,
                                    currentSplit,
                                    par);
    } else if (par.seqMode == 2) {
        fillQueryKmerBufferParallel_paired(kseq1,
                                           kseq2,
                                           kmerBuffer,
                                           queryList,
                                           currentSplit,
                                           par);
    }
    cout << double(time(nullptr) - beforeKmerExtraction) << " s" << endl;

    // Sort query k-mer
    time_t beforeQueryKmerSort = time(nullptr);
    std::cout << "Query k-mer sorting    : " << flush;
    SORT_PARALLEL(kmerBuffer.buffer, kmerBuffer.buffer + kmerBuffer.startIndexOfReserve, Kmer::compareQueryKmer);
    cout << double(time(nullptr) - beforeQueryKmerSort) << " s" << endl;
}

void KmerExtractor::fillQueryKmerBufferParallel(KSeqWrapper *kseq,
                                                Buffer<Kmer> &kmerBuffer,
                                                std::vector<Query> &queryList,
                                                const QuerySplit &currentSplit,
                                                const LocalParameters &par) {   
    size_t readLength = 1000;
    size_t processedQueryNum = 0;
    size_t chunkSize = 1000;
 
    // Reserve memory for each read of the chunk for each thread
    std::vector<std::vector<string>> chunkReads_thread(par.threads);
    for (int i = 0; i < par.threads; ++i) {
        chunkReads_thread[i].resize(chunkSize);
        for (size_t j = 0; j < chunkSize; ++j) {
            chunkReads_thread[i][j].reserve(readLength);
        }
    }

    // Vector to check empty reads
    std::vector<std::vector<bool>> emptyReads(par.threads);
    for (int i = 0; i < par.threads; ++i) {
        emptyReads[i].resize(chunkSize);
        for (size_t j = 0; j < chunkSize; ++j) {
            emptyReads[i][j] = false;
        }
    }

    // Initialize atomic variable for active tasks
    std::vector<atomic<bool>> busyThreads(par.threads);
    for (int i = 0; i < par.threads; ++i) {
        busyThreads[i].store(false);
    }

    // OpenMP parallel region with tasks
#pragma omp parallel default(none) shared(par, readLength, kmerBuffer, \
queryList, currentSplit, processedQueryNum, kseq, chunkSize, chunkReads_thread, busyThreads, emptyReads)
    {
        char *maskedSeq = new char[readLength];
        char *seq = nullptr;     
#pragma omp single nowait
        {
            int masterThread = 0;
            #ifdef OPENMP
               masterThread = omp_get_thread_num();
            #endif
            size_t count = 0;
            while (processedQueryNum < currentSplit.readCnt) {
                // Find an idle thread
                int threadId;
                if (par.threads == 1) {
                    threadId = 0;
                } else {
                    for (threadId = 0; threadId < par.threads; ++threadId) {
                        if (threadId == masterThread) { continue; }
                        if (!busyThreads[threadId].load()) {
                            busyThreads[threadId].store(true);
                            break;
                        }
                    }
                }
                
                if (threadId == par.threads) {
                    continue;
                }

                size_t chunkEnd = min(processedQueryNum + chunkSize, currentSplit.readCnt);
                count = 0;

                loadChunkOfReads(kseq,
                                 queryList,
                                 processedQueryNum,
                                 chunkSize,
                                 chunkEnd,
                                 chunkReads_thread[threadId],
                                 emptyReads[threadId],
                                 count,
                                 false);

                if (par.threads == 1) {
                    processSequence(count, processedQueryNum, chunkReads_thread[threadId], emptyReads[threadId], seq, maskedSeq, readLength, kmerBuffer, queryList, false);
                } else {
                    #pragma omp task firstprivate(count, processedQueryNum, threadId)
                    {
                        processSequence(count, processedQueryNum, chunkReads_thread[threadId], emptyReads[threadId], seq, maskedSeq, readLength, kmerBuffer, queryList, false);
                        busyThreads[threadId].store(false);
                    }
                }  
            }
        }
        delete[] maskedSeq;
    }
}

void KmerExtractor::fillQueryKmerBufferParallel_paired(KSeqWrapper *kseq1,
                                                       KSeqWrapper *kseq2,
                                                       Buffer<Kmer> &kmerBuffer,
                                                       vector<Query> &queryList,
                                                       const QuerySplit &currentSplit,
                                                       const LocalParameters &par) {
    size_t readLength = 1000;
    size_t processedQueryNum = 0;
    size_t chunkSize = 1000;
    // Reserve memory for each read of the chunk for each thread
    std::vector<std::vector<string>> chunkReads1_thread(par.threads);
    std::vector<std::vector<string>> chunkReads2_thread(par.threads);
    for (int i = 0; i < par.threads; ++i) {
        chunkReads1_thread[i].resize(chunkSize);
        chunkReads2_thread[i].resize(chunkSize);
        for (size_t j = 0; j < chunkSize; ++j) {
            chunkReads1_thread[i][j].reserve(readLength);
            chunkReads2_thread[i][j].reserve(readLength);
        }
    }
    // Vector to check empty reads
    std::vector<std::vector<bool>> emptyReads(par.threads);
    for (int i = 0; i < par.threads; ++i) {
        emptyReads[i].resize(chunkSize);
        for (size_t j = 0; j < chunkSize; ++j) {
            emptyReads[i][j] = false;
        }
    }

    // Initialize atomic variable for active tasks
    std::vector<std::atomic<bool>> busyThreads(par.threads);
    for (int i = 0; i < par.threads; ++i) {
        busyThreads[i].store(false);
    }

    // OpenMP parallel region with tasks
#pragma omp parallel default(none) shared(par, readLength, \
kmerBuffer, queryList, currentSplit, processedQueryNum, kseq1, kseq2, \
chunkSize, chunkReads1_thread, chunkReads2_thread, busyThreads, cout, emptyReads)
    {
        size_t maxReadLength1 = 1000;
        size_t maxReadLength2 = 1000;
        char *maskedSeq1 = new char[maxReadLength1];
        char *maskedSeq2 = new char[maxReadLength2];   
        char *seq1 = nullptr;
        char *seq2 = nullptr;  
#pragma omp single
        {
            int masterThread = 0;
            #ifdef OPENMP
               masterThread = omp_get_thread_num();
            #endif
           
            while (processedQueryNum < currentSplit.readCnt) {
                // Find an idle thread
                int threadId;
                if (par.threads == 1) {
                    threadId = 0;
                } else {
                    for (threadId = 0; threadId < par.threads; ++threadId) {
                        if (threadId == masterThread) { continue; }
                        if (!busyThreads[threadId].load()) {
                            busyThreads[threadId].store(true);
                            break;
                        }
                    }
                }
                if (threadId == par.threads) {
                    continue;
                }

                size_t chunkEnd = min(processedQueryNum + chunkSize, currentSplit.readCnt);
                
                size_t count = 0;
                size_t processedQueryNumCopy = processedQueryNum;
                loadChunkOfReads(kseq1,
                                 queryList,
                                 processedQueryNum,
                                 chunkSize,
                                 chunkEnd,
                                 chunkReads1_thread[threadId],
                                 emptyReads[threadId],
                                 count,
                                 false);

                count = 0;
                processedQueryNum = processedQueryNumCopy;
                loadChunkOfReads(kseq2,
                                 queryList,
                                 processedQueryNum,
                                 chunkSize,
                                 chunkEnd,
                                 chunkReads2_thread[threadId],
                                 emptyReads[threadId],
                                 count,
                                 true);

                if (par.threads == 1) {
                    processSequence(count, processedQueryNum, chunkReads1_thread[threadId], emptyReads[threadId], seq1, maskedSeq1, maxReadLength1, kmerBuffer, queryList, false);
                    processSequence(count, processedQueryNum, chunkReads2_thread[threadId], emptyReads[threadId], seq2, maskedSeq2, maxReadLength2, kmerBuffer, queryList, true);
                } 
                else {
                    #pragma omp task firstprivate(count, processedQueryNum, threadId)
                    {
                        processSequence(count, processedQueryNum, chunkReads1_thread[threadId], emptyReads[threadId], seq1, maskedSeq1, maxReadLength1, kmerBuffer, queryList, false);
                        processSequence(count, processedQueryNum, chunkReads2_thread[threadId], emptyReads[threadId], seq2, maskedSeq2, maxReadLength2, kmerBuffer, queryList, true);
                        busyThreads[threadId].store(false);
                    }
                }   
            }   
        }   
        delete[] maskedSeq1;
        delete[] maskedSeq2;
    }
}

void KmerExtractor::processSequence(
    size_t count,
    size_t processedQueryNum,
    const vector<string> & reads,
    const vector<bool> & emptyReads,
    char *seq,
    char *maskedSeq,
    size_t & maxReadLength,
    Buffer<Kmer> &kmerBuffer,
    const vector<Query> & queryList,
    bool isReverse) 
{
    for (size_t i = 0; i < count; ++i) {
        size_t queryIdx = processedQueryNum - count + i;
        if (emptyReads[i]) { continue; }
        // Get masked sequence
        if (maskMode) {
            if (maxReadLength < reads[i].length() + 1) {
                maxReadLength = reads[i].length() + 1;
                delete[] maskedSeq;
                maskedSeq = new char[maxReadLength];
            }
            SeqIterator::maskLowComplexityRegions((unsigned char *) reads[i].c_str(), (unsigned char *) maskedSeq, *probMatrix, maskProb, subMat);
            seq = maskedSeq;
        } else {
            seq = const_cast<char *>(reads[i].c_str());
        }

        size_t posToWrite = 0;
        if (isReverse) {
            posToWrite = kmerBuffer.reserveMemory(queryList[queryIdx].kmerCnt2);
            fillQueryKmerBuffer(
                seq,
                (int) reads[i].length(),
                kmerBuffer,
                posToWrite, 
                (uint32_t) queryIdx+1,
                queryList[queryIdx].queryLength+3);
        } else {
            posToWrite = kmerBuffer.reserveMemory(queryList[queryIdx].kmerCnt);
            fillQueryKmerBuffer(
                seq, 
                (int) reads[i].length(), 
                kmerBuffer, 
                posToWrite, 
                (uint32_t) queryIdx+1);                             
        }
    }
}

void KmerExtractor::fillQueryKmerBuffer(
    const char *seq,
    int seqLen, 
    Buffer<Kmer> &kmerBuffer, 
    size_t &posToWrite, 
    uint32_t seqID, 
    uint32_t offset) 
{
    int usedLen = LocalUtil::getMaxCoveredLength(seqLen);
#ifdef OPENMP
    size_t threadID = omp_get_thread_num();
#else
    size_t threadID = 0; // Single-threaded mode
#endif
    for (int frame = 0; frame < 6; frame++) {
        bool isForward = frame < 3;
        int begin = 0;
        if (isForward) {
            begin = frame % 3;
        } else {
            begin = (seqLen % 3) - (frame % 3);
            if (begin < 0) {
                begin += 3;
            }
        }
        kmerScanners[threadID]->initScanner(seq, begin, begin + usedLen - 1, isForward);
        Kmer kmer;
        while ((kmer = kmerScanners[threadID]->next()).value != UINT64_MAX) {
            kmerBuffer.buffer[posToWrite++] = {kmer.value, seqID, kmer.pos + offset, (uint8_t) frame};
        }
    }
}


int KmerExtractor::extractTargetKmers(
    const char *seq,
    Buffer<Kmer> &kmerBuffer,
    size_t &posToWrite,
    int taxId,
    int spTaxId,
    SequenceBlock block) 
{
#ifdef OPENMP
    size_t threadID = omp_get_thread_num();
#else
    size_t threadID = 0; // Single-threaded mode
#endif
    kmerScanners[threadID]->initScanner(seq, block.start, block.end, (block.strand > -1));
    Kmer kmer;
    while ((kmer = kmerScanners[threadID]->next()).value != UINT64_MAX) {
        kmerBuffer.buffer[posToWrite++] = { kmer.value, taxId, spTaxId };
    }
    return 0;
}


void KmerExtractor::loadChunkOfReads(KSeqWrapper *kseq,
                                     vector<Query> & queryList,
                                     size_t & processedQueryNum,
                                     size_t chunkSize,
                                     size_t chunkEnd,
                                     vector<string> & reads,
                                     vector<bool> & emptyReads,
                                     size_t & count,
                                     bool isReverse) {
    if (isReverse) {
        for (size_t i = 0; i < chunkSize && processedQueryNum < chunkEnd; ++i, ++processedQueryNum) {
            kseq->ReadEntry();
            queryList[processedQueryNum].queryLength2 = LocalUtil::getMaxCoveredLength((int) kseq->entry.sequence.l);   
            
            if (emptyReads[i]) { 
                count ++;
                continue; 
            }

            // Check if the read is too short
            int kmerCnt = LocalUtil::getQueryKmerNumber<int>((int) kseq->entry.sequence.l, spaceNum, this->kmerLen);
            if (kmerCnt < 1) {
                reads[i] = "";
                emptyReads[i] = true;
                queryList[processedQueryNum].kmerCnt2 = 0;
            } else {
                reads[i] = string(kseq->entry.sequence.s);
                // emptyReads[i] = false; // already set to false
                queryList[processedQueryNum].kmerCnt2 = kmerCnt;
            }
            count++;
        }
    } else {
        for (size_t i = 0; i < chunkSize && processedQueryNum < chunkEnd; ++i, ++processedQueryNum) {
            kseq->ReadEntry();
            queryList[processedQueryNum].name = string(kseq->entry.name.s);
            queryList[processedQueryNum].queryLength = LocalUtil::getMaxCoveredLength((int) kseq->entry.sequence.l);

            // Check if the read is too short
            int kmerCnt = LocalUtil::getQueryKmerNumber<int>((int) kseq->entry.sequence.l, spaceNum, this->kmerLen);
            if (kmerCnt < 1) {
                reads[i] = "";
                emptyReads[i] = true;
                queryList[processedQueryNum].kmerCnt = 0;
            } else {
                reads[i] = string(kseq->entry.sequence.s);
                emptyReads[i] = false;
                queryList[processedQueryNum].kmerCnt = kmerCnt;
            }
            count++;
        }
    }
}


bool KmerExtractor::extractKmers(
    KSeqWrapper *kseq,
    Buffer<Kmer> &kmerBuffer,
    std::unordered_map<string, uint32_t> & accession2index,
    uint32_t & idOffset,
    SeqEntry & savedSeq)
{
    if (!savedSeq.name.empty()) {
        if (savedSeq.s.length() >= 12) {
            size_t writePos = kmerBuffer.reserveMemory(savedSeq.s.length() - 11);
            kmerScanners[0]->initScanner(savedSeq.s.c_str(), 0, savedSeq.s.length() - 1, true);
            Kmer kmer;
            while((kmer = kmerScanners[0]->next()).value != UINT64_MAX) {
                kmerBuffer.buffer[writePos++] = {kmer.value, idOffset};
            }
            idOffset += 1;
        }
        savedSeq.s.clear(); 
        savedSeq.name.clear();
    }

    size_t seqLen = 1000;
    size_t chunkSize = 100;
    std::vector<std::vector<string>> readsPerThread(par.threads);
    std::vector<std::vector<string>> seqNamesPerThread(par.threads);
    for (int i = 0; i < par.threads; ++i) {
        readsPerThread[i].resize(chunkSize);
        seqNamesPerThread[i].resize(chunkSize);
        for (size_t j = 0; j < chunkSize; ++j) {
            readsPerThread[i][j].reserve(seqLen);
            seqNamesPerThread[i][j].reserve(128);
        }
    }

    BlockingQueue<int> freeBufferQueue;
    BlockingQueue<WorkItem> filledBufferQueue;

    for (int i = 1; i < par.threads; ++i) {
        freeBufferQueue.push(i); 
    }
    int producerBufferIndex = 0;
    bool moreData = true;
    #pragma omp parallel default(none) shared(par, seqLen, kmerBuffer, kseq, std::cout, \
    readsPerThread, seqNamesPerThread, chunkSize, filledBufferQueue, moreData, accession2index, savedSeq, producerBufferIndex, freeBufferQueue, idOffset)
    {
        int threadId = omp_get_thread_num();

        if (threadId == 0) {
            bool bufferNotFull = true;
            while (moreData && bufferNotFull) {
                auto & seqChunk = readsPerThread[producerBufferIndex];
                auto & nameChunk = seqNamesPerThread[producerBufferIndex];
                size_t seqCnt = 0;
                size_t kmerCnt = 0;
                for (; seqCnt < chunkSize; ++seqCnt) {
                    if (!kseq->ReadEntry()) {
                        // std::cout << "No more sequences to read" << std::endl;
                        moreData = false; // No more data to read
                        break;
                    }
                    size_t currentKmerCnt 
                        = (kseq->entry.sequence.l > 11) ? (kseq->entry.sequence.l - 11) : 0;
                    
                    if (kmerBuffer.startIndexOfReserve + kmerCnt + currentKmerCnt >= kmerBuffer.bufferSize) {
                        savedSeq.s = kseq->entry.sequence.s;
                        savedSeq.name = kseq->entry.name.s;
                        bufferNotFull = false;
                        // std::cout << "Buffer is full" << std::endl;
                        break;
                    }

                    kmerCnt += currentKmerCnt;
                    nameChunk[seqCnt] = kseq->entry.name.s;
                    seqChunk[seqCnt] = (currentKmerCnt > 0) ? kseq->entry.sequence.s : "";
                }
                if (seqCnt > 0) {
                    filledBufferQueue.push({producerBufferIndex, seqCnt, kmerBuffer.startIndexOfReserve, idOffset});
                    kmerBuffer.reserveMemory(kmerCnt);
                    idOffset += seqCnt;
                    if (moreData && bufferNotFull) {
                        freeBufferQueue.wait_and_pop(producerBufferIndex); // Wait for an idle thread
                    }
                }
            }

            for (int i = 0; i < par.threads - 1; ++i) {
                filledBufferQueue.push({-1, 0, 0, 0}); // Signal to other threads that no more data will be produced
            }
        } else {
            WorkItem workItem;
            std::unordered_map<string, uint32_t> localAccession2Index;
            while (true) {
                filledBufferQueue.wait_and_pop(workItem); 
                if (workItem.bufferIndex == -1) { break; }
                uint32_t currentId = workItem.sequenceIdOffset;
                size_t posToWrite = workItem.writePos;
                const auto & seqChunk = readsPerThread[workItem.bufferIndex];
                const auto & nameChunk = seqNamesPerThread[workItem.bufferIndex];
                for (size_t i = 0; i < workItem.sequenceCount; ++i) {
                    if (seqChunk[i].empty()) {
                        continue; // Skip empty reads
                    }
                    localAccession2Index[nameChunk[i]] = currentId;
                    kmerScanners[threadId]->initScanner(seqChunk[i].c_str(), 0, seqChunk[i].length() - 1, true);
                    Kmer kmer;
                    while((kmer = kmerScanners[threadId]->next()).value != UINT64_MAX) {
                        kmerBuffer.buffer[posToWrite++] = {kmer.value, currentId};
                    }
                    currentId++;
                }
                freeBufferQueue.push(workItem.bufferIndex);
            }
            // Merge local accession2index into global map
            #pragma omp critical
            {
                for (const auto & entry : localAccession2Index) {
                    accession2index[entry.first] = entry.second;
                }
            }
        }
    }
    return moreData;
}

bool KmerExtractor::extractUnirefKmers(
    KSeqWrapper *kseq,
    Buffer<Kmer> &kmerBuffer,
    std::unordered_map<string, uint32_t> & unirefName2Id,
    uint32_t & processedSeqCnt,
    SeqEntry & savedSeq)
{
    if (!savedSeq.name.empty()) {
        if (savedSeq.s.length() >= 12) {
            size_t writePos = kmerBuffer.reserveMemory(savedSeq.s.length() - 11);
            TaxID currentId = unirefName2Id[savedSeq.name];
            kmerScanners[0]->initScanner(savedSeq.s.c_str(), 0, savedSeq.s.length() - 1, true);
            Kmer kmer;
            while((kmer = kmerScanners[0]->next()).value != UINT64_MAX) {
                kmerBuffer.buffer[writePos++] = {kmer.value, currentId};
            }
        }
        savedSeq.s.clear(); 
        savedSeq.name.clear();
        ++processedSeqCnt;
    }
    size_t seqLen = 1000;
    size_t chunkSize = 100;
    std::vector<std::vector<string>> readsPerThread(par.threads);
    std::vector<std::vector<string>> seqNamesPerThread(par.threads);
    for (int i = 0; i < par.threads; ++i) {
        readsPerThread[i].resize(chunkSize);
        seqNamesPerThread[i].resize(chunkSize);
        for (size_t j = 0; j < chunkSize; ++j) {
            readsPerThread[i][j].reserve(seqLen);
            seqNamesPerThread[i][j].reserve(128);
        }
    }

    BlockingQueue<int> freeBufferQueue;
    BlockingQueue<WorkItem> filledBufferQueue;
    for (int i = 1; i < par.threads; ++i) {
        freeBufferQueue.push(i); 
    }
    int producerBufferIndex = 0;
    bool moreData = true;
    #pragma omp parallel default(none) shared(par, kmerBuffer, kseq, std::cout, processedSeqCnt, \
    readsPerThread, seqNamesPerThread, chunkSize, filledBufferQueue, unirefName2Id, moreData, savedSeq, producerBufferIndex, freeBufferQueue)
    {
        int threadId = omp_get_thread_num();

        if (threadId == 0) {
            bool bufferNotFull = true;
            while (moreData && bufferNotFull) {
                auto & seqChunk = readsPerThread[producerBufferIndex];
                auto & nameChunk = seqNamesPerThread[producerBufferIndex];
                size_t seqCnt = 0;
                size_t kmerCnt = 0;
                for (; seqCnt < chunkSize; ++seqCnt) {
                    if (!kseq->ReadEntry()) {
                        moreData = false; // No more data to read
                        break;
                    }
                    size_t currentKmerCnt 
                        = (kseq->entry.sequence.l > 11) ? (kseq->entry.sequence.l - 11) : 0;
                    
                    if (kmerBuffer.startIndexOfReserve + kmerCnt + currentKmerCnt >= kmerBuffer.bufferSize) {
                        savedSeq.s = kseq->entry.sequence.s;
                        savedSeq.name = kseq->entry.name.s;
                        bufferNotFull = false;
                        break;
                    }

                    kmerCnt += currentKmerCnt;
                    nameChunk[seqCnt] = kseq->entry.name.s;
                    seqChunk[seqCnt] = (currentKmerCnt > 0) ? kseq->entry.sequence.s : "";
                }
                if (seqCnt > 0) {
                    filledBufferQueue.push({producerBufferIndex, seqCnt, kmerBuffer.startIndexOfReserve, 0});
                    processedSeqCnt += seqCnt;
                    kmerBuffer.reserveMemory(kmerCnt);
                    if (moreData && bufferNotFull) {
                        freeBufferQueue.wait_and_pop(producerBufferIndex); // Wait for an idle thread
                    }
                }
            }

            for (int i = 0; i < par.threads - 1; ++i) {
                filledBufferQueue.push({-1, 0, 0, 0}); // Signal to other threads that no more data will be produced
            }
        } else {
            WorkItem workItem;
            while (true) {
                filledBufferQueue.wait_and_pop(workItem); 
                if (workItem.bufferIndex == -1) { break; }
                // uint32_t currentId = workItem.sequenceIdOffset;
                size_t posToWrite = workItem.writePos;
                const auto & seqChunk = readsPerThread[workItem.bufferIndex];
                const auto & nameChunk = seqNamesPerThread[workItem.bufferIndex];
                for (size_t i = 0; i < workItem.sequenceCount; ++i) {
                    if (seqChunk[i].empty()) {
                        continue; // Skip empty reads
                    }
                    TaxID currentId = unirefName2Id[nameChunk[i]];
                    // std::cout << nameChunk[i] << std::endl;
                    kmerScanners[threadId]->initScanner(seqChunk[i].c_str(), 0, seqChunk[i].length() - 1, true);
                    Kmer kmer;
                    while((kmer = kmerScanners[threadId]->next()).value != UINT64_MAX) {
                        kmerBuffer.buffer[posToWrite++] = {kmer.value, currentId};
                    }
                }
                freeBufferQueue.push(workItem.bufferIndex);
            }
        }
    }
    return moreData;
}

bool KmerExtractor::extractQueryKmers_aa2aa(
    Buffer<Kmer> &kmerBuffer,
    std::vector<ProteinQuery> & queryList,
    KSeqWrapper* kseq,
    uint64_t & processedSeqCnt,
    SeqEntry & savedSeq
) {
    uint32_t queryIdx = 1;
    if (!savedSeq.name.empty()) {
        queryList[queryIdx].queryId = queryIdx;
        queryList[queryIdx].name = savedSeq.name;
        queryList[queryIdx].length = savedSeq.s.length();
        queryList[queryIdx].kmerCnt = (savedSeq.s.length() >= 12) ? (savedSeq.s.length() - 11) : 0;

        if (savedSeq.s.length() >= 12) {
            size_t writePos = kmerBuffer.reserveMemory(savedSeq.s.length() - 11);
            kmerScanners[0]->initScanner(savedSeq.s.c_str(), 0, savedSeq.s.length() - 1, true);
            Kmer kmer;
            while((kmer = kmerScanners[0]->next()).value != UINT64_MAX) {
                kmerBuffer.buffer[writePos++] = {kmer.value, queryIdx}; // The first query of each iteration has ID 1
            }
        }
        savedSeq.s.clear(); 
        savedSeq.name.clear();
        ++processedSeqCnt;
        queryIdx ++;
    }
    

    size_t seqLen = 1000;
    size_t chunkSize = 100;
    std::vector<std::vector<string>> readsPerThread(par.threads);
    // std::vector<std::vector<string>> seqNamesPerThread(par.threads);
    for (int i = 0; i < par.threads; ++i) {
        readsPerThread[i].resize(chunkSize);
        // seqNamesPerThread[i].resize(chunkSize);
        for (size_t j = 0; j < chunkSize; ++j) {
            readsPerThread[i][j].reserve(seqLen);
            // seqNamesPerThread[i][j].reserve(128);
        }
    }

    BlockingQueue<int> freeBufferQueue;
    BlockingQueue<WorkItem> filledBufferQueue;
    for (int i = 1; i < par.threads; ++i) {
        freeBufferQueue.push(i); 
    }
    int producerBufferIndex = 0;
    bool moreData = true;
    #pragma omp parallel default(none) shared(par, kmerBuffer, kseq, std::cout, processedSeqCnt, queryIdx, \
    readsPerThread, chunkSize, filledBufferQueue, moreData, savedSeq, producerBufferIndex, queryList, freeBufferQueue)
    {
        int threadId = omp_get_thread_num();

        if (threadId == 0) {
            bool bufferNotFull = true;
            while (moreData && bufferNotFull) {
                auto & seqChunk = readsPerThread[producerBufferIndex];
                // auto & nameChunk = seqNamesPerThread[producerBufferIndex];
                size_t seqCnt = 0;
                size_t kmerCnt = 0;
                while (queryIdx + chunkSize - 1 >= queryList.size()) {
                    queryList.resize(queryList.size() * 2);
                }
                for (; seqCnt < chunkSize; ++seqCnt) {
                    if (!kseq->ReadEntry()) {
                        moreData = false; // No more data to read
                        break;
                    }
                    size_t currentKmerCnt 
                        = (kseq->entry.sequence.l > 11) ? (kseq->entry.sequence.l - 11) : 0;
                    
                    if (kmerBuffer.startIndexOfReserve + kmerCnt + currentKmerCnt >= kmerBuffer.bufferSize) {
                        savedSeq.s = kseq->entry.sequence.s;
                        savedSeq.name = kseq->entry.name.s;
                        bufferNotFull = false;
                        break;
                    }
                    kmerCnt += currentKmerCnt;
                    queryList[queryIdx + seqCnt] = {queryIdx + seqCnt,
                                                    (uint32_t) kseq->entry.sequence.l,
                                                    (uint32_t) currentKmerCnt,
                                                    0, 0.0,
                                                    string(kseq->entry.name.s)};
                    seqChunk[seqCnt] = (currentKmerCnt > 0) ? kseq->entry.sequence.s : "";
                }
                if (seqCnt > 0) {
                    filledBufferQueue.push({producerBufferIndex, seqCnt, kmerBuffer.startIndexOfReserve, queryIdx});
                    queryIdx += seqCnt;
                    processedSeqCnt += seqCnt;
                    kmerBuffer.reserveMemory(kmerCnt);
                    if (moreData && bufferNotFull) {
                        freeBufferQueue.wait_and_pop(producerBufferIndex); // Wait for an idle thread
                    }
                }
            }

            for (int i = 0; i < par.threads - 1; ++i) {
                filledBufferQueue.push({-1, 0, 0, 0}); // Signal to other threads that no more data will be produced
            }
        } else {
            WorkItem workItem;
            while (true) {
                filledBufferQueue.wait_and_pop(workItem); 

                if (workItem.bufferIndex == -1) { break; }
                // uint32_t idOffset = workItem.sequenceIdOffset;
                // size_t posToWrite = workItem.writePos;
                // const auto & seqChunk = readsPerThread[workItem.bufferIndex];

                processSequenceChunk_aa2aa(
                    kmerBuffer,
                    workItem.writePos,
                    workItem.sequenceIdOffset,
                    readsPerThread[workItem.bufferIndex],
                    workItem.sequenceCount,
                    threadId
                );
                freeBufferQueue.push(workItem.bufferIndex);
            }
        }
    }
    return moreData;
}

void KmerExtractor::processSequenceChunk_aa2aa(
    Buffer<Kmer> &kmerBuffer,
    size_t writePos,
    uint32_t queryOffset,
    const std::vector<std::string> & aaSeqs, 
    size_t seqNum,
    int threadID
) {
    for (size_t i = 0; i < seqNum; ++i) {
        if (aaSeqs[i].empty()) {continue;}
        kmerScanners[threadID]->initScanner(aaSeqs[i].c_str(), 0, aaSeqs[i].length() - 1, true);
        uint32_t id = queryOffset + i;
        Kmer kmer;
        while((kmer = kmerScanners[threadID]->next()).value != UINT64_MAX) {
            kmerBuffer.buffer[writePos++] = {kmer.value, id, kmer.pos, 0};
            // cout << kmer.pos << "\t"; kmer.printAA(geneticCode, 12); cout <<endl;
        }
    }
}
