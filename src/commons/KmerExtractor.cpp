#include "KmerExtractor.h"
#include <unordered_map>

KmerExtractor::KmerExtractor(
    const LocalParameters &par,
    const GeneticCode & geneticCode,
    int kmerFormat) 
    : par(par), kmerScanners(new KmerScanner*[par.threads]) {
    // Initialize k-mer scanners for each thread
    for (int i = 0; i < par.threads; ++i) {
        if (kmerFormat == 1) {
            kmerScanners[i] = new OldKmerScanner(geneticCode);
        } else {
            if (par.syncmer) {
                kmerScanners[i] = new SyncmerScanner(par.smerLen, geneticCode);
            } else {
                kmerScanners[i] = new KmerScanner(geneticCode);
            }
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

void KmerExtractor::extractQueryKmers(Buffer<QueryKmer> &kmerBuffer,
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
    SORT_PARALLEL(kmerBuffer.buffer, kmerBuffer.buffer + kmerBuffer.startIndexOfReserve, compareForLinearSearch);
    cout << double(time(nullptr) - beforeQueryKmerSort) << " s" << endl;
}

void KmerExtractor::fillQueryKmerBufferParallel(KSeqWrapper *kseq,
                                                Buffer<QueryKmer> &kmerBuffer,
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
                                                       Buffer<QueryKmer> &kmerBuffer,
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
#pragma omp single nowait
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
    Buffer<QueryKmer> &kmerBuffer,
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
    Buffer<QueryKmer> &kmerBuffer, 
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
    Buffer<TargetKmer> &kmerBuffer,
    size_t &posToWrite,
    int seqID,
    int taxIdAtRank,
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
        kmerBuffer.buffer[posToWrite++] = { kmer.value, taxIdAtRank, seqID };
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
            int kmerCnt = LocalUtil::getQueryKmerNumber<int>((int) kseq->entry.sequence.l, spaceNum);
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
            int kmerCnt = LocalUtil::getQueryKmerNumber<int>((int) kseq->entry.sequence.l, spaceNum);
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