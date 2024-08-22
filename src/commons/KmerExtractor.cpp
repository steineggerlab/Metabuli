#include "KmerExtractor.h"
#include <unordered_map>

KmerExtractor::KmerExtractor(const LocalParameters &par) {
    seqIterator = new SeqIterator(par);
    spaceNum = 0;
    maskMode = par.maskMode;
    maskProb = par.maskProb;
    subMat = new NucleotideMatrix(par.scoringMatrixFile.values.nucleotide().c_str(), 1.0, 0.0);
    probMatrix = new ProbabilityMatrix(*(subMat));
}

KmerExtractor::~KmerExtractor() {
    delete seqIterator;
    delete probMatrix;
    delete subMat;
}

void KmerExtractor::extractQueryKmers(QueryKmerBuffer &kmerBuffer,
                                      vector<Query> & queryList,
                                      const QuerySplit & currentSplit,
                                      const LocalParameters &par,
                                      KSeqWrapper* kseq1,
                                      KSeqWrapper* kseq2) {
    time_t beforeKmerExtraction = time(nullptr);
    std::cout << "Extracting query metamers ... " << endl;
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
    cout << "Time spent for metamer extraction: " << double(time(nullptr) - beforeKmerExtraction) << endl;

    // Sort query k-mer
    time_t beforeQueryKmerSort = time(nullptr);
    cout << "Sorting query metamer list ..." << endl;
    SORT_PARALLEL(kmerBuffer.buffer, kmerBuffer.buffer + kmerBuffer.startIndexOfReserve, compareForLinearSearch);
    cout << "Time spent for sorting query metamer list: " << double(time(nullptr) - beforeQueryKmerSort) << endl;
}

void KmerExtractor::fillQueryKmerBufferParallel(KSeqWrapper *kseq,
                                                QueryKmerBuffer &kmerBuffer,
                                                std::vector<Query> &queryList,
                                                const QuerySplit &currentSplit,
                                                const LocalParameters &par) {   
    size_t readLength = 1000;
    size_t processedQueryNum = 0;
    size_t chunkSize = 1000;
 
    // Reserve memory for each read of the chunk for each thread
    std::vector<std::vector<string>> chunkReads_thread(par.threads);
    for (size_t i = 0; i < par.threads; ++i) {
        chunkReads_thread[i].resize(chunkSize);
        for (size_t j = 0; j < chunkSize; ++j) {
            chunkReads_thread[i][j].reserve(readLength);
        }
    }

    // Initialize atomic variable for active tasks
    std::vector<std::atomic<bool>> busyThreads(par.threads);
    for (size_t i = 0; i < par.threads; ++i) {
        busyThreads[i].store(false);
    }

    // OpenMP parallel region with tasks
#pragma omp parallel default(none) shared(par, readLength, kmerBuffer, \
queryList, currentSplit, processedQueryNum, kseq, chunkSize, chunkReads_thread, busyThreads)
    {
        char *maskedSeq = new char[readLength];
        char *seq = nullptr;
        vector<int> aaFrames_read[6];
        // reserve memory for each frame
        for (int i = 0; i < 6; ++i) {
            aaFrames_read[i].reserve(readLength/3);
        }     
#pragma omp single nowait
        {
            int masterThread = omp_get_thread_num();
            size_t count = 0;
            while (processedQueryNum < currentSplit.readCnt) {
                // Find an idle thread
                int threadId;
                for (threadId = 0; threadId < par.threads; ++threadId) {
                    if (threadId == masterThread) { continue; }
                    if (!busyThreads[threadId].load()) {
                        busyThreads[threadId].store(true);
                        break;
                    }
                }
                
                if (threadId == par.threads) {
                    continue;
                }

                size_t chunkEnd = min(processedQueryNum + chunkSize, currentSplit.readCnt);
                count = 0;

                for (size_t i = 0; i < chunkSize && processedQueryNum < chunkEnd; ++i, ++processedQueryNum) {
                    kseq->ReadEntry();
                    queryList[processedQueryNum].name = string(kseq->entry.name.s);
                    queryList[processedQueryNum].queryLength = LocalUtil::getMaxCoveredLength((int) kseq->entry.sequence.l);
                    
                    // Check if the read is too short
                    int kmerCnt = LocalUtil::getQueryKmerNumber<int>((int) kseq->entry.sequence.l, spaceNum);
                    if (kmerCnt < 1) {
                        chunkReads_thread[threadId][i] = "";
                        queryList[processedQueryNum].kmerCnt = 0;
                    } else {
                        chunkReads_thread[threadId][i] = string(kseq->entry.sequence.s);
                        queryList[processedQueryNum].kmerCnt = kmerCnt;
                    }
                    count++;
                }
                
                // Process each chunk by idle threads
                #pragma omp task firstprivate(count, processedQueryNum, threadId)
                {
                    for (size_t i = 0; i < count; ++i) {
                        size_t queryIdx = processedQueryNum - count + i;
                        if (chunkReads_thread[threadId][i].empty()) { continue; }
                        
                        // Get masked sequence
                        if (maskMode) {
                            if (readLength < chunkReads_thread[threadId][i].length() + 1) {
                                readLength = chunkReads_thread[threadId][i].length() + 1;
                                delete[] maskedSeq;
                                maskedSeq = new char[readLength];
                            }
                            SeqIterator::maskLowComplexityRegions(chunkReads_thread[threadId][i].c_str(), maskedSeq, *probMatrix, maskProb, subMat);
                            seq = maskedSeq;
                        } else {
                            seq = const_cast<char *>(chunkReads_thread[threadId][i].c_str());
                        }

                        size_t posToWrite = kmerBuffer.reserveMemory(queryList[queryIdx].kmerCnt);
                        
                        // Process Read 1
                        seqIterator->sixFrameTranslation(seq, (int) chunkReads_thread[threadId][i].length(), aaFrames_read);
                        seqIterator->fillQueryKmerBuffer(seq, (int) chunkReads_thread[threadId][i].length(), kmerBuffer, posToWrite, 
                                                        (uint32_t) queryIdx, aaFrames_read);
                    }
                    busyThreads[threadId].store(false);
                }   
            }
        }
        delete[] maskedSeq;
    }
}

void KmerExtractor::fillQueryKmerBufferParallel_paired(KSeqWrapper *kseq1,
                                                       KSeqWrapper *kseq2,
                                                       QueryKmerBuffer &kmerBuffer,
                                                       vector<Query> &queryList,
                                                       const QuerySplit &currentSplit,
                                                       const LocalParameters &par) {
    size_t readLength = 1000;
    size_t processedQueryNum = 0;
    size_t chunkSize = 1000;
    // Reserve memory for each read of the chunk for each thread
    std::vector<std::vector<string>> chunkReads1_thread(par.threads);
    std::vector<std::vector<string>> chunkReads2_thread(par.threads);
    for (size_t i = 0; i < par.threads; ++i) {
        chunkReads1_thread[i].resize(chunkSize);
        chunkReads2_thread[i].resize(chunkSize);
        for (size_t j = 0; j < chunkSize; ++j) {
            chunkReads1_thread[i][j].reserve(readLength);
            chunkReads2_thread[i][j].reserve(readLength);
        }
    }

    // Initialize atomic variable for active tasks
    std::vector<std::atomic<bool>> busyThreads(par.threads);
    for (size_t i = 0; i < par.threads; ++i) {
        busyThreads[i].store(false);
    }

    // OpenMP parallel region with tasks
#pragma omp parallel default(none) shared(par, readLength, \
kmerBuffer, queryList, currentSplit, processedQueryNum, kseq1, kseq2, \
chunkSize, chunkReads1_thread, chunkReads2_thread, busyThreads, cout)
    {
        char *maskedSeq1 = new char[readLength];
        char *maskedSeq2 = new char[readLength];   
        char *seq1 = nullptr;
        char *seq2 = nullptr;
        vector<int> aaFrames_read1[6];
        vector<int> aaFrames_read2[6];
        // reserve memory for each frame
        for (int i = 0; i < 6; ++i) {
            aaFrames_read1[i].reserve(readLength/3);
            aaFrames_read2[i].reserve(readLength/3);
        }     
#pragma omp single nowait
        {
            int masterThread = omp_get_thread_num();
            size_t count = 0;
            while (processedQueryNum < currentSplit.readCnt) {
                // Find an idle thread
                int threadId;
                for (threadId = 0; threadId < par.threads; ++threadId) {
                    if (threadId == masterThread) { continue; }
                    if (!busyThreads[threadId].load()) {
                        busyThreads[threadId].store(true);
                        break;
                    }
                }

                if (threadId == par.threads) {
                    continue;
                }

                size_t chunkEnd = min(processedQueryNum + chunkSize, currentSplit.readCnt);
                count = 0;

                for (size_t i = 0; i < chunkSize && processedQueryNum < chunkEnd; ++i, ++processedQueryNum) {
                    kseq1->ReadEntry();
                    kseq2->ReadEntry();
                    queryList[processedQueryNum].name = string(kseq1->entry.name.s);
                    queryList[processedQueryNum].queryLength = LocalUtil::getMaxCoveredLength((int) kseq1->entry.sequence.l);
                    queryList[processedQueryNum].queryLength2 = LocalUtil::getMaxCoveredLength((int) kseq2->entry.sequence.l);
                    
                    // Check if the read is too short
                    int kmerCnt1 = LocalUtil::getQueryKmerNumber<int>((int) kseq1->entry.sequence.l, spaceNum);
                    int kmerCnt2 = LocalUtil::getQueryKmerNumber<int>((int) kseq2->entry.sequence.l, spaceNum);
                    if (kmerCnt1 < 1 || kmerCnt2 < 1) {
                        chunkReads1_thread[threadId][i] = "";
                        chunkReads2_thread[threadId][i] = "";
                        queryList[processedQueryNum].kmerCnt = 0;
                    } else {
                        // cout << threadId << endl
                        chunkReads1_thread[threadId][i] = string(kseq1->entry.sequence.s);
                        chunkReads2_thread[threadId][i] = string(kseq2->entry.sequence.s);
                        queryList[processedQueryNum].kmerCnt = kmerCnt1 + kmerCnt2;
                    }
                    count++;
                }
                

                // Process each chunk by idle threads
                #pragma omp task firstprivate(count, processedQueryNum, threadId)
                {
                    for (size_t i = 0; i < count; ++i) {
                        size_t queryIdx = processedQueryNum - count + i;
                        if (chunkReads1_thread[threadId][i].empty() || chunkReads2_thread[threadId][i].empty()) { continue; }
                        
                        // Get masked sequence
                        if (maskMode) {
                            if (readLength < chunkReads1_thread[threadId][i].length() + 1 || readLength < chunkReads2_thread[threadId][i].length() + 1) {
                                readLength = max(chunkReads1_thread[threadId][i].length() + 1, chunkReads2_thread[threadId][i].length() + 1);
                                delete[] maskedSeq1;
                                delete[] maskedSeq2;
                                maskedSeq1 = new char[readLength];
                                maskedSeq2 = new char[readLength];
                            }
                            SeqIterator::maskLowComplexityRegions(chunkReads1_thread[threadId][i].c_str(),maskedSeq1, *probMatrix, maskProb, subMat);
                            SeqIterator::maskLowComplexityRegions(chunkReads2_thread[threadId][i].c_str(),maskedSeq2, *probMatrix, maskProb, subMat);
                            seq1 = maskedSeq1;
                            seq2 = maskedSeq2;
                        } else {
                            seq1 = const_cast<char *>(chunkReads1_thread[threadId][i].c_str());
                            seq2 = const_cast<char *>(chunkReads2_thread[threadId][i].c_str());
                        }

                        size_t posToWrite = kmerBuffer.reserveMemory(queryList[queryIdx].kmerCnt);
                        
                        // Process Read 1
                        seqIterator->sixFrameTranslation(seq1, (int) chunkReads1_thread[threadId][i].length(), aaFrames_read1);
                        seqIterator->fillQueryKmerBuffer(seq1, (int) chunkReads1_thread[threadId][i].length(), kmerBuffer, posToWrite, 
                                                        (uint32_t) queryIdx, aaFrames_read1);

                        // Process Read 2
                        seqIterator->sixFrameTranslation(seq2, (int) chunkReads2_thread[threadId][i].length(), aaFrames_read2);
                        seqIterator->fillQueryKmerBuffer(seq2, (int) chunkReads2_thread[threadId][i].length(), kmerBuffer, posToWrite,
                                                        (uint32_t) queryIdx, aaFrames_read2, queryList[queryIdx].queryLength+3);
                    }
                    busyThreads[threadId].store(false);
                }   
            }
        }
        delete[] maskedSeq1;
        delete[] maskedSeq2;
    }
}
    


// void KmerExtractor::fillQueryKmerBufferParallel_paired2(KSeqWrapper *kseq1,
//                                                        KSeqWrapper *kseq2,
//                                                        QueryKmerBuffer &kmerBuffer,
//                                                        vector<Query> &queryList,
//                                                        const QuerySplit &currentSplit,
//                                                        const LocalParameters &par) {
//     size_t processedQueryNum = 0;

//     // Array to store reads of thread number
//     vector<string> reads1(par.threads);
//     vector<string> reads2(par.threads);

//     while (processedQueryNum < currentSplit.readCnt) {
//         size_t currentQueryNum = min(currentSplit.readCnt - processedQueryNum, (size_t) par.threads);
//         size_t count = 0;

//         // Fill reads in sequential
//         while (count < currentQueryNum) {
//             // Read query
//             kseq1->ReadEntry();
//             kseq2->ReadEntry();
//             const KSeqWrapper::KSeqEntry & e1 = kseq1->entry;
//             const KSeqWrapper::KSeqEntry & e2 = kseq2->entry;

//             // Get k-mer count
//             int kmerCnt = LocalUtil::getQueryKmerNumber<int>((int) e1.sequence.l, spaceNum);
//             int kmerCnt2 = LocalUtil::getQueryKmerNumber<int>((int) e2.sequence.l, spaceNum);

//             // Query Info
//             queryList[processedQueryNum].queryLength = LocalUtil::getMaxCoveredLength((int) e1.sequence.l);
//             queryList[processedQueryNum].queryLength2 = LocalUtil::getMaxCoveredLength((int) e2.sequence.l);
//             queryList[processedQueryNum].name = string(e1.name.s);
//             queryList[processedQueryNum].kmerCnt = (int) (kmerCnt + kmerCnt2);

//             // Store reads
//             reads1[count] = string(kseq1->entry.sequence.s);
//             reads2[count] = string(kseq2->entry.sequence.s);

//             processedQueryNum ++;
//             count ++;
//         }

//         // Process reads in parallel
// #pragma omp parallel default(none), shared(par, kmerBuffer, cout, processedQueryNum, queryList, currentQueryNum, currentSplit, count, reads1, reads2)
//         {
//             size_t posToWrite;
//             size_t readLength = 1000;
//             char *maskedSeq1 = new char[readLength];
//             char *maskedSeq2 = new char[readLength];
//             char *seq1 = nullptr;
//             char *seq2 = nullptr;            
// #pragma omp for schedule(dynamic, 1)
//             for (size_t i = 0; i < currentQueryNum; i ++) {
//                 size_t queryIdx = processedQueryNum - currentQueryNum + i;
//                 // Get k-mer count
//                 int kmerCnt = LocalUtil::getQueryKmerNumber<int>(reads1[i].length(), spaceNum);
//                 int kmerCnt2 = LocalUtil::getQueryKmerNumber<int>(reads2[i].length(), spaceNum);

//                 // Ignore short read
//                 if (kmerCnt2 < 1 || kmerCnt < 1) { continue; }

//                 // Get masked sequence
//                 if (maskMode) {
//                     if (readLength < reads1[i].length() + 1 || readLength < reads2[i].length() + 1) {
//                         readLength = max(reads1[i].length() + 1, reads2[i].length() + 1);
//                         delete[] maskedSeq1;
//                         delete[] maskedSeq2;
//                         maskedSeq1 = new char[readLength];
//                         maskedSeq2 = new char[readLength];
//                     }
//                     SeqIterator::maskLowComplexityRegions(reads1[i].c_str(),maskedSeq1, *probMatrix, maskProb, subMat);
//                     SeqIterator::maskLowComplexityRegions(reads2[i].c_str(),maskedSeq2, *probMatrix, maskProb, subMat);
//                     seq1 = maskedSeq1;
//                     seq2 = maskedSeq2;
//                 } else {
//                     seq1 = const_cast<char *>(reads1[i].c_str());
//                     seq2 = const_cast<char *>(reads2[i].c_str());
//                 }

//                 posToWrite = kmerBuffer.reserveMemory(kmerCnt + kmerCnt2);

//                 // Process Read 1
//                 seqIterator->sixFrameTranslation(seq1, (int) reads1[i].length());
//                 seqIterator->fillQueryKmerBuffer(seq1, (int) reads1[i].length(), kmerBuffer, posToWrite,
//                                                 (uint32_t) queryIdx);

//                 // Process Read 2
//                 seqIterator->sixFrameTranslation(seq2, (int) reads2[i].length());
//                 seqIterator->fillQueryKmerBuffer(seq2, (int) reads2[i].length(), kmerBuffer, posToWrite,
//                                                  (uint32_t) queryIdx, queryList[queryIdx].queryLength+3);
//             }
//             delete[] maskedSeq1;
//             delete[] maskedSeq2;
//         }
//     }
// }

