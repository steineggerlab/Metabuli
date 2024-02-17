#include "KmerExtractor.h"

KmerExtractor::KmerExtractor(const LocalParameters &par) {
    spaceNum = par.spaceMask.length() - 8;
    maskMode = par.maskMode;
    maskProb = par.maskProb;
    subMat = new NucleotideMatrix(par.scoringMatrixFile.values.nucleotide().c_str(), 1.0, 0.0);
    probMatrix = new ProbabilityMatrix(*(subMat));
}

KmerExtractor::~KmerExtractor() {
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

void KmerExtractor::fillQueryKmerBufferParallel(KSeqWrapper *kseq1,
                                                QueryKmerBuffer &kmerBuffer,
                                                vector<Query> &queryList,
                                                const QuerySplit &currentSplit,
                                                const LocalParameters &par) {                                                   
    size_t processedQueryNum = 0;
 
     // Array to store reads of thread number
     vector<string> reads1(par.threads);
 
    while (processedQueryNum < currentSplit.readCnt) {
        size_t currentQueryNum = min(currentSplit.readCnt - processedQueryNum, (size_t) par.threads);
        size_t count = 0;
        while (count < currentQueryNum) {
            // Read query
            kseq1->ReadEntry();
            const KSeqWrapper::KSeqEntry & e1 = kseq1->entry;

            // Get k-mer count
            int kmerCnt = LocalUtil::getQueryKmerNumber<int>((int) e1.sequence.l, spaceNum);

            // Query Info
            queryList[processedQueryNum].queryLength = LocalUtil::getMaxCoveredLength((int) e1.sequence.l);
            queryList[processedQueryNum].name = string(e1.name.s);
            queryList[processedQueryNum].kmerCnt = (int) (kmerCnt);

            // Store reads
            reads1[count] = string(kseq1->entry.sequence.s);

            processedQueryNum ++;
            count ++;
        }
#pragma omp parallel default(none), shared(par, kmerBuffer, cout, processedQueryNum, queryList, currentQueryNum, currentSplit, count, reads1)
        {
            SeqIterator seqIterator(par);
            size_t posToWrite;
#pragma omp for schedule(dynamic, 1)
            for (size_t i = 0; i < currentQueryNum; i ++) {
                size_t queryIdx = processedQueryNum - currentQueryNum + i;
                // Get k-mer count
                int kmerCnt = LocalUtil::getQueryKmerNumber<int>(reads1[i].length(), spaceNum);
                
                // Ignore short read
                if (kmerCnt < 1) { continue; }

                // Get masked sequence
                char *maskedSeq1 = nullptr;
                if (maskMode) {
                    maskedSeq1 = new char[reads1[i].length() + 1];
                    SeqIterator::maskLowComplexityRegions(reads1[i].c_str(),maskedSeq1, *probMatrix, maskProb, subMat);
                } else {
                    maskedSeq1 = const_cast<char *>(reads1[i].c_str());
                }

                posToWrite = kmerBuffer.reserveMemory(kmerCnt);

                // Process Read 1
                seqIterator.sixFrameTranslation(maskedSeq1, (int) reads1[i].length());
                seqIterator.fillQueryKmerBuffer(maskedSeq1, (int) reads1[i].length(), kmerBuffer, posToWrite,
                                                (uint32_t) queryIdx);

                if (maskMode) {
                    delete[] maskedSeq1;
                }
            }
        }
    }
}

void KmerExtractor::fillQueryKmerBufferParallel_paired(KSeqWrapper *kseq1,
                                                       KSeqWrapper *kseq2,
                                                       QueryKmerBuffer &kmerBuffer,
                                                       vector<Query> &queryList,
                                                       const QuerySplit &currentSplit,
                                                       const LocalParameters &par) {
    size_t processedQueryNum = 0;

    // Array to store reads of thread number
    vector<string> reads1(par.threads);
    vector<string> reads2(par.threads);

    while (processedQueryNum < currentSplit.readCnt) {
        size_t currentQueryNum = min(currentSplit.readCnt - processedQueryNum, (size_t) par.threads);
        size_t count = 0;

        // Fill reads in sequential
        while (count < currentQueryNum) {
            // Read query
            kseq1->ReadEntry();
            kseq2->ReadEntry();
            const KSeqWrapper::KSeqEntry & e1 = kseq1->entry;
            const KSeqWrapper::KSeqEntry & e2 = kseq2->entry;

            // Get k-mer count
            int kmerCnt = LocalUtil::getQueryKmerNumber<int>((int) e1.sequence.l, spaceNum);
            int kmerCnt2 = LocalUtil::getQueryKmerNumber<int>((int) e2.sequence.l, spaceNum);

            // Query Info
            queryList[processedQueryNum].queryLength = LocalUtil::getMaxCoveredLength((int) e1.sequence.l);
            queryList[processedQueryNum].queryLength2 = LocalUtil::getMaxCoveredLength((int) e2.sequence.l);
            queryList[processedQueryNum].name = string(e1.name.s);
            queryList[processedQueryNum].kmerCnt = (int) (kmerCnt + kmerCnt2);

            // Store reads
            reads1[count] = string(kseq1->entry.sequence.s);
            reads2[count] = string(kseq2->entry.sequence.s);

            processedQueryNum ++;
            count ++;
        }

        // Process reads in parallel
#pragma omp parallel default(none), shared(par, kmerBuffer, cout, processedQueryNum, queryList, currentQueryNum, currentSplit, count, reads1, reads2)
        {
            SeqIterator seqIterator(par);
            SeqIterator seqIterator2(par);
            size_t posToWrite;
#pragma omp for schedule(dynamic, 1)
            for (size_t i = 0; i < currentQueryNum; i ++) {
                size_t queryIdx = processedQueryNum - currentQueryNum + i;
                // Get k-mer count
                auto kmerCnt = LocalUtil::getQueryKmerNumber<size_t>(reads1[i].length(), spaceNum);
                auto kmerCnt2 = LocalUtil::getQueryKmerNumber<size_t>(reads2[i].length(), spaceNum);

                // Ignore short read
                if (kmerCnt2 < 1 || kmerCnt < 1) { continue; }

                // Get masked sequence
                char *maskedSeq1 = nullptr;
                char *maskedSeq2 = nullptr;
                if (maskMode) {
                    maskedSeq1 = new char[reads1[i].length() + 1];
                    maskedSeq2 = new char[reads2[i].length() + 1];
                    SeqIterator::maskLowComplexityRegions(reads1[i].c_str(),maskedSeq1, *probMatrix, maskProb, subMat);
                    SeqIterator::maskLowComplexityRegions(reads2[i].c_str(),maskedSeq2, *probMatrix, maskProb, subMat);
                } else {
                    maskedSeq1 = const_cast<char *>(reads1[i].c_str());
                    maskedSeq2 = const_cast<char *>(reads2[i].c_str());
                }

                posToWrite = kmerBuffer.reserveMemory(kmerCnt + kmerCnt2);

                // Process Read 1
                seqIterator.sixFrameTranslation(maskedSeq1, (int) reads1[i].length());
                seqIterator.fillQueryKmerBuffer(maskedSeq1, (int) reads1[i].length(), kmerBuffer, posToWrite,
                                                (uint32_t) queryIdx);

                // Process Read 2
                seqIterator2.sixFrameTranslation(maskedSeq2, (int) reads2[i].length());
                seqIterator2.fillQueryKmerBuffer(maskedSeq2, (int) reads2[i].length(), kmerBuffer, posToWrite,
                                                 (uint32_t) queryIdx, queryList[queryIdx].queryLength+3);

                if (maskMode) {
                    delete[] maskedSeq1;
                    delete[] maskedSeq2;
                }
            }
        }
    }
}

