#include "Classifier.h"
#include "LocalParameters.h"
#include <ctime>

Classifier::Classifier(LocalParameters & par) {
    MARKER = 16777215;
    MARKER = ~ MARKER;
    bitsForCodon = 3;
    numOfSplit = 0;
    minConsCnt = par.minConsCnt;
    minSpScore = par.minSpScore;

    // Mask for spaced k-mer
    size_t maskLen = par.spaceMask.length();
    mask = new uint32_t[maskLen];
    spaceNum = 0;
    spaceNum_int = 0;
    for(size_t i = 0, j = 0; i < maskLen; i++){
        mask[i] = par.spaceMask[i] - 48;
        spaceNum += (mask[i] == 0);
        spaceNum_int += (mask[i] == 0);
        if(mask[i]==1){
            unmaskedPos[j] = (int) i;
            j++;
        }
    }

    // Hamming Dist. margin
    hammingMargin = (uint8_t) par.hammingMargin;
}

Classifier::~Classifier() {
    delete[] mask;
}

void Classifier::startClassify(const char *queryFileName,
                               const char *targetDiffIdxFileName,
                               const char *targetInfoFileName,
                               const char *diffIdxSplitFileName,
                               vector<int> &taxIdList,
                               const LocalParameters &par,
                               NcbiTaxonomy &taxonomy) {

    vector<int> speciesTaxIdList;
    vector<TaxID> genusTaxIdList;
    taxonomy.createTaxIdListAtRank(taxIdList, speciesTaxIdList, "species");
    taxonomy.createTaxIdListAtRank(taxIdList, genusTaxIdList, "genus");

    //output file
//    char matchFileName[300];
//    sprintf(matchFileName, "%s_match2", queryFileName);
//    FILE *matchFile = fopen(matchFileName, "wb");

    // Allocate memory for buffers
    QueryKmerBuffer kmerBuffer(kmerBufSize);
    Buffer<Match> matchBuffer(size_t(kmerBufSize) * size_t(20));

    // Load query file
    cout << "Indexing query file ...";
    MmapedData<char> queryFile{};
    MmapedData<char> queryFile2{};
    vector<Sequence> sequences;
    vector<Sequence> sequences2;
    Query *queryList;
    size_t numOfSeq;
    size_t numOfSeq2;
    if (par.seqMode == 1 || par.seqMode == 3) {
        queryFile = mmapData<char>(par.filenames[0].c_str());
        IndexCreator::getSeqSegmentsWithHead(sequences, queryFile);
        numOfSeq = sequences.size();
        queryList = new Query[numOfSeq];
    } else if (par.seqMode == 2) {
        string queryFileName1 = par.filenames[0] + "_1";
        string queryFileName2 = par.filenames[0] + "_2";
        queryFile = mmapData<char>(queryFileName1.c_str());
        queryFile2 = mmapData<char>(queryFileName2.c_str());
        IndexCreator::getSeqSegmentsWithHead(sequences, queryFile);
        IndexCreator::getSeqSegmentsWithHead(sequences2, queryFile2);
        numOfSeq = sequences.size();
        numOfSeq2 = sequences2.size();
        if (numOfSeq > numOfSeq2) {
            queryList = new Query[numOfSeq];
        } else {
            numOfSeq = numOfSeq2;
            queryList = new Query[numOfSeq];
        }
    }
    cout << "Done" << endl;

    // Checker for multi-threading
    bool *processedSeqChecker = new bool[numOfSeq];
    fill_n(processedSeqChecker, numOfSeq, false);
    size_t processedSeqCnt = 0;

    //
    size_t numOfTatalQueryKmerCnt = 0;
    size_t totalMatchCnt = 0;
    // Extract k-mers from query sequences and compare them to target k-mer DB
    omp_set_num_threads(par.threads);
    while (processedSeqCnt < numOfSeq) {
        time_t beforeKmerExtraction = time(nullptr);

        // Initialize query k-mer buffer and match buffer
        kmerBuffer.startIndexOfReserve = 0;
        matchBuffer.startIndexOfReserve = 0;

        // Extract query k-mer
        cout << "K-mer extraction ... " << endl;
        if (par.seqMode == 1 || par.seqMode == 3) { // Single-end short-read sequence or long-read sequence
            fillQueryKmerBufferParallel(kmerBuffer,queryFile,sequences,processedSeqChecker,processedSeqCnt,
                                        queryList, par);
        } else if (par.seqMode == 2) {
            fillQueryKmerBufferParallel_paired(kmerBuffer,queryFile,queryFile2,sequences,sequences2,
                                               processedSeqChecker, processedSeqCnt, queryList, numOfSeq, par);
        }
        numOfTatalQueryKmerCnt += kmerBuffer.startIndexOfReserve;
        cout << "Time spent for k-mer extraction: " << double(time(nullptr) - beforeKmerExtraction) << endl;

        // Sort query k-mer
        time_t beforeQueryKmerSort = time(nullptr);
        SORT_PARALLEL(kmerBuffer.buffer, kmerBuffer.buffer + kmerBuffer.startIndexOfReserve,
                      Classifier::compareForLinearSearch);
        cout << "Time spent for sorting query k-mer list: " << double(time(nullptr) - beforeQueryKmerSort) << endl;

        // Search matches between query and target k-mers
        linearSearchParallel(kmerBuffer.buffer, kmerBuffer.startIndexOfReserve, targetDiffIdxFileName,
                             targetInfoFileName,
                             diffIdxSplitFileName, matchBuffer, taxIdList, speciesTaxIdList, genusTaxIdList,
                             par);

        // Sort matches
        time_t beforeSortMatches = time(nullptr);
        totalMatchCnt += matchBuffer.startIndexOfReserve;
        SORT_PARALLEL(matchBuffer.buffer, matchBuffer.buffer + matchBuffer.startIndexOfReserve,
                      Classifier::sortByGenusAndSpecies2);
        cout << "Time spent for sorting matches: " << double(time(nullptr) - beforeSortMatches) << endl;

        // Classify queries based on the matches
        time_t beforeAnalyze = time(nullptr);
        analyseResultParallel(taxonomy, matchBuffer.buffer, matchBuffer.startIndexOfReserve, (int) numOfSeq, queryList,
                              par);
        cout << "Time spent for analyzing: " << double(time(nullptr) - beforeAnalyze) << endl;
        cout << "The number of processed sequences: " << processedSeqCnt << " (" << processedSeqCnt / numOfSeq << ")" << endl;
    }
    cout << "Number of query k-mers: " << numOfTatalQueryKmerCnt << endl;
    cout << "The number of matches: " << totalMatchCnt << endl;

    // Write report files
    ofstream readClassificationFile;
    readClassificationFile.open(par.filenames[3] + "/" + par.filenames[4] + "_ReadClassification.tsv");
    writeReadClassification(queryList, (int) numOfSeq, readClassificationFile);
    readClassificationFile.close();
    writeReportFile(par.filenames[3] + "/" + par.filenames[4] + "_CompositionReport.tsv", taxonomy, numOfSeq);

    //Below is for developing
//    ofstream wr;
//    ofstream wr2;
//    vector<int> wrongClassifications;
//    sequences.clear();
//    IndexCreator::getSeqSegmentsWithHead(sequences, queryFile);
//    wr.open(par.filenames[0]+"_wrong_1");
//    if(par.seqMode == 2) {
//        sequences2.clear();
//        IndexCreator::getSeqSegmentsWithHead(sequences2, queryFile2);
//        wr2.open(par.filenames[0] + "_wrong_2");
//    }
//    performanceTest(taxonomy, queryList, numOfSeq, wrongClassifications);
//    for (size_t i = 0; i < wrongClassifications.size(); i++) {
//        kseq_buffer_t buffer(const_cast<char *>(&queryFile.data[sequences[wrongClassifications[i]].start]), sequences[wrongClassifications[i]].length);
//        kseq_t *seq = kseq_init(&buffer);
//        kseq_read(seq);
//        wr<<">"<<seq->name.s<<endl;
//        wr<<seq->seq.s<<endl;
//        kseq_destroy(seq);
//        if(par.seqMode == 2) {
//            kseq_buffer_t buffer2(const_cast<char *>(&queryFile2.data[sequences[wrongClassifications[i]].start]),
//                                  sequences2[wrongClassifications[i]].length);
//            kseq_t *seq2 = kseq_init(&buffer2);
//            kseq_read(seq2);
//            wr2 << ">" << seq2->name.s << endl;
//            wr2 << seq2->seq.s << endl;
//            kseq_destroy(seq2);
//        }
//    }
//    wr.close();
//    wr2.close();

    munmap(queryFile.data, queryFile.fileSize + 1);
    if (par.seqMode == 2) {
        munmap(queryFile2.data, queryFile2.fileSize + 1);
    }
}

void Classifier::fillQueryKmerBufferParallel(QueryKmerBuffer &kmerBuffer,
                                             MmapedData<char> &seqFile,
                                             vector<Sequence> &seqs,
                                             bool *checker,
                                             size_t &processedSeqCnt,
                                             Query *queryList,
                                             const LocalParameters &par) {
    bool hasOverflow = false;
    omp_set_num_threads(*(int *) par.PARAM_THREADS.value);
#pragma omp parallel default(none), shared(par, checker, hasOverflow, processedSeqCnt, kmerBuffer, seqFile, seqs, cout, queryList)
    {
        SeqIterator seqIterator(par);
        size_t posToWrite;
#pragma omp for schedule(dynamic, 1)
        for (size_t i = 0; i < seqs.size(); i++) {
            if (checker[i] == false && !hasOverflow) {
                kseq_buffer_t buffer(const_cast<char *>(&seqFile.data[seqs[i].start]), seqs[i].length);
                kseq_t *seq = kseq_init(&buffer);
                kseq_read(seq);
                seqIterator.sixFrameTranslation(seq->seq.s);
                int kmerCnt = getQueryKmerNumber((int) strlen(seq->seq.s));
                posToWrite = kmerBuffer.reserveMemory(kmerCnt);

                // Ignore short read
                if (kmerCnt < 1) continue;
                if (posToWrite + kmerCnt < kmerBuffer.bufferSize) {
                    seqIterator.fillQueryKmerBuffer(seq->seq.s, kmerBuffer, posToWrite, i);
                    checker[i] = true;
                    queryList[i].queryLength = getMaxCoveredLength((int) strlen(seq->seq.s));
                    queryList[i].queryId = i;
                    queryList[i].name = string(seq->name.s);
                    queryList[i].kmerCnt = kmerCnt;
#pragma omp atomic
                    processedSeqCnt++;
                } else {
#pragma omp atomic
                    kmerBuffer.startIndexOfReserve -= kmerCnt;
                    hasOverflow = true;
                }
                kseq_destroy(seq);
            }
        }
    }
}


int Classifier::getMaxCoveredLength(int queryLength) {
    if (queryLength % 3 == 2) {
        return queryLength - 2; // 2
    } else if (queryLength % 3 == 1) {
        return queryLength - 4; // 4
    } else {
        return queryLength - 3; // 3
    }
}

int Classifier::getQueryKmerNumber(int queryLength) {
    return (getMaxCoveredLength(queryLength) / 3 - kmerLength - spaceNum_int + 1) * 6;
}

void Classifier::fillQueryKmerBufferParallel_paired(QueryKmerBuffer &kmerBuffer,
                                                    MmapedData<char> &seqFile1,
                                                    MmapedData<char> &seqFile2,
                                                    vector<Sequence> &seqs,
                                                    vector<Sequence> &seqs2,
                                                    bool *checker,
                                                    size_t &processedSeqCnt,
                                                    Query *queryList,
                                                    size_t numOfSeq,
                                                    const LocalParameters &par) {
    bool hasOverflow = false;

#pragma omp parallel default(none), shared(par, checker, hasOverflow, processedSeqCnt, kmerBuffer, seqFile1, seqFile2, seqs, seqs2, cout, queryList, numOfSeq)
    {
        SeqIterator seqIterator(par);
        SeqIterator seqIterator2(par);
        size_t posToWrite;
#pragma omp for schedule(dynamic, 1)
        for (size_t i = 0; i < numOfSeq; i++) {
            if (!checker[i] && !hasOverflow) {
                // Read 1
                kseq_buffer_t buffer(const_cast<char *>(&seqFile1.data[seqs[i].start]), seqs[i].length);
                kseq_t *seq = kseq_init(&buffer);
                kseq_read(seq);
                int kmerCnt = getQueryKmerNumber((int) strlen(seq->seq.s));

                // Read 2
                kseq_buffer_t buffer2(const_cast<char *>(&seqFile2.data[seqs2[i].start]), seqs2[i].length);
                kseq_t *seq2 = kseq_init(&buffer2);
                kseq_read(seq2);
                int kmerCnt2 = getQueryKmerNumber((int) strlen(seq2->seq.s));

                // Ignore short read
                if (kmerCnt2 < 1 || kmerCnt < 1) {
                    processedSeqCnt++;
                    continue;
                }

                posToWrite = kmerBuffer.reserveMemory(kmerCnt + kmerCnt2);
                if (posToWrite + kmerCnt + kmerCnt2 < kmerBuffer.bufferSize) {
                    checker[i] = true;
                    // Read 1
                    seqIterator.sixFrameTranslation(seq->seq.s);
                    seqIterator.fillQueryKmerBuffer(seq->seq.s, kmerBuffer, posToWrite, (int) i);
                    queryList[i].queryLength = getMaxCoveredLength((int) strlen(seq->seq.s));

                    // Read 2
                    seqIterator2.sixFrameTranslation(seq2->seq.s);
                    seqIterator2.fillQueryKmerBuffer(seq2->seq.s, kmerBuffer, posToWrite, (int) i,
                                                     queryList[i].queryLength);

                    // Query Info
                    queryList[i].queryLength2 = getMaxCoveredLength((int) strlen(seq2->seq.s));
                    queryList[i].queryId = (int) i;
                    queryList[i].name = string(seq->name.s);
                    queryList[i].kmerCnt = kmerCnt + kmerCnt2;
#pragma omp atomic
                    processedSeqCnt++;
                } else {
#pragma omp atomic
                    kmerBuffer.startIndexOfReserve -= kmerCnt + kmerCnt2;
                    hasOverflow = true;
                }
                kseq_destroy(seq);
                kseq_destroy(seq2);
            }
        }
    }
}

void Classifier::linearSearchParallel(QueryKmer *queryKmerList, size_t &queryKmerCnt, const char *targetDiffIdxFileName,
                                      const char *targetInfoFileName, const char *diffIdxSplitsFileName,
                                      Buffer<Match> &matchBuffer, const vector<int> &taxIdList,
                                      const vector<int> &spTaxIdList, const vector<TaxID> &genusTaxIdList,
                                      const LocalParameters &par) {

    struct MmapedData<uint16_t> targetDiffIdxList = mmapData<uint16_t>(targetDiffIdxFileName, 2);
    struct MmapedData<TargetKmerInfo> targetInfoList = mmapData<TargetKmerInfo>(targetInfoFileName, 2);
    struct MmapedData<DiffIdxSplit> diffIdxSplits = mmapData<DiffIdxSplit>(diffIdxSplitsFileName);

    cout << "linearSearch start..." << endl;
    SeqIterator seqIterator1(par);
    // Find the first index of garbage query k-mer (UINT64_MAX) and discard from there
    for (size_t checkN = queryKmerCnt - 1; checkN > 0; checkN--) {
        if (queryKmerList[checkN].ADkmer != UINT64_MAX) {
            queryKmerCnt = checkN + 1;
            break;
        }
    }

    // Filter out meaningless target splits
    size_t numOfDiffIdxSplits = diffIdxSplits.fileSize / sizeof(DiffIdxSplit);
    size_t numOfDiffIdxSplits_use = numOfDiffIdxSplits;
    for (size_t i = 1; i < numOfDiffIdxSplits; i++) {
        if (diffIdxSplits.data[i].ADkmer == 0 || diffIdxSplits.data[i].ADkmer == UINT64_MAX) {
            diffIdxSplits.data[i] = {UINT64_MAX, UINT64_MAX, UINT64_MAX};
            numOfDiffIdxSplits_use--;
        }
    }

    cout << "Filtering out meaningless target splits ... done" << endl;

    // Divide query k-mer list into blocks for multi threading.
    // Each split has start and end points of query list + proper offset point of target k-mer list
    vector<QueryKmerSplit> querySplits;

    int threadNum = par.threads;
    uint64_t queryAA;
    if (threadNum == 1) { //Single thread
        querySplits.emplace_back(0, queryKmerCnt - 1, queryKmerCnt, diffIdxSplits.data[0]);
    } else if (threadNum == 2) { //Two threads
        size_t splitWidth = queryKmerCnt / 2;
        querySplits.emplace_back(0, splitWidth - 1, splitWidth, diffIdxSplits.data[0]);
        for (size_t tSplitCnt = 0; tSplitCnt < numOfDiffIdxSplits_use; tSplitCnt++) {
            queryAA = AminoAcidPart(queryKmerList[splitWidth].ADkmer);
            if (queryAA <= AminoAcidPart(diffIdxSplits.data[tSplitCnt].ADkmer)) {
                tSplitCnt = tSplitCnt - (tSplitCnt != 0);
                querySplits.emplace_back(splitWidth, queryKmerCnt - 1, queryKmerCnt - splitWidth,
                                         diffIdxSplits.data[tSplitCnt]);
                break;
            }
        }
    } else { //More than two threads
        size_t splitWidth = queryKmerCnt / (threadNum - 1);
        querySplits.emplace_back(0, splitWidth - 1, splitWidth, diffIdxSplits.data[0]);
        for (int i = 1; i < threadNum; i++) {
            queryAA = AminoAcidPart(queryKmerList[splitWidth * i].ADkmer);
            bool needLastTargetBlock = true;
            for (size_t j = 0; j < numOfDiffIdxSplits_use; j++) {
                if (queryAA <= AminoAcidPart(diffIdxSplits.data[j].ADkmer)) {
                    j = j - (j != 0);
                    if (i != threadNum - 1) {
                        querySplits.emplace_back(splitWidth * i, splitWidth * (i + 1) - 1, splitWidth,
                                                 diffIdxSplits.data[j]);
                    } else {
                        querySplits.emplace_back(splitWidth * i, queryKmerCnt - 1, queryKmerCnt - splitWidth * i,
                                                 diffIdxSplits.data[j]);
                    }
                    needLastTargetBlock = false;
                    break;
                }
            }
            if (needLastTargetBlock) {
                if (i != threadNum - 1) {
                    querySplits.emplace_back(splitWidth * i, splitWidth * (i + 1) - 1, splitWidth,
                                             diffIdxSplits.data[numOfDiffIdxSplits_use - 1]);
                } else {
                    querySplits.emplace_back(splitWidth * i, queryKmerCnt - 1, queryKmerCnt - splitWidth * i,
                                             diffIdxSplits.data[numOfDiffIdxSplits_use - 1]);
                }
            }
        }
    }

    bool *splitCheckList = (bool *) malloc(sizeof(bool) * threadNum);
    fill_n(splitCheckList, threadNum, false);
    int completedSplitCnt = 0;
    size_t numOfTargetKmer = targetInfoList.fileSize / sizeof(TargetKmerInfo);
    size_t numOfDiffIdx = targetDiffIdxList.fileSize / sizeof(uint16_t);

    cout << "The number of target k-mers: " << numOfTargetKmer << endl;

    time_t beforeSearch = time(nullptr);

    vector<vector<TaxID>> sspOrSp;
    sspOrSp.push_back(move(taxIdList));
    sspOrSp.push_back(move(spTaxIdList));
    while (completedSplitCnt < threadNum) {
        bool hasOverflow = false;
#pragma omp parallel default(none), shared(numOfDiffIdx, completedSplitCnt, splitCheckList, numOfTargetKmer, hasOverflow, \
querySplits, queryKmerList, targetDiffIdxList, targetInfoList, matchBuffer, cout, genusTaxIdList, taxIdList, spTaxIdList, sspOrSp, par)
        {
            //query variables
            uint64_t currentQuery = UINT64_MAX;
            uint64_t currentQueryAA = UINT64_MAX;

            //target variables
            size_t diffIdxPos = 0;
            size_t targetInfoIdx = 0;
            vector<uint64_t> candidateTargetKmers; //vector for candidate target k-mer, some of which are selected after based on hamming distance
            uint64_t currentTargetKmer;

            //Match buffer for each thread
            int localBufferSize = 2000000; // 64 Mb
            auto *matches = new Match[localBufferSize];
            int matchCnt = 0;

            // For debug
            SeqIterator seqIterator(par);

            //vectors for selected target k-mers
            vector<uint8_t> selectedHammingSum;
            vector<size_t> selectedMatches;
            vector<uint16_t> selectedHammings;
            size_t startIdxOfAAmatch = 0;
            size_t posToWrite;

            int currMatchNum;
            bool red;
            size_t idx;
            size_t range;
#pragma omp for schedule(dynamic, 1)
            for (size_t i = 0; i < querySplits.size(); i++) {
                if (hasOverflow || splitCheckList[i]) {
                    continue;
                }

                targetInfoIdx = querySplits[i].diffIdxSplit.infoIdxOffset - (i != 0);
                diffIdxPos = querySplits[i].diffIdxSplit.diffIdxOffset;
                currentTargetKmer = querySplits[i].diffIdxSplit.ADkmer;
                if (i == 0) {
                    currentTargetKmer = getNextTargetKmer(currentTargetKmer, targetDiffIdxList.data, diffIdxPos);
                }
                currentQuery = UINT64_MAX;
                currentQueryAA = UINT64_MAX;

                size_t lastMovedQueryIdx;
                for (size_t j = querySplits[i].start; j < querySplits[i].end + 1; j++) {
                    querySplits[i].start++;

                    // Reuse the comparison data if queries are exactly identical
                    if (currentQuery == queryKmerList[j].ADkmer) {
                        currMatchNum = selectedMatches.size();
                        // If local buffer is full, copy them to the shared buffer.
                        if (matchCnt + currMatchNum > localBufferSize) {
                            // Check if the shared buffer is full.
                            posToWrite = matchBuffer.reserveMemory(matchCnt);
                            if (posToWrite + matchCnt >=
                                matchBuffer.bufferSize) { // full -> write matches to file first
                                hasOverflow = true;

                                querySplits[i].start = lastMovedQueryIdx + 1;
                                __sync_fetch_and_sub(& matchBuffer.startIndexOfReserve, matchCnt);
                                break;
                            } else { // not full -> copy matches to the shared buffer
                                moveMatches(matchBuffer.buffer + posToWrite, matches, matchCnt);
                                lastMovedQueryIdx = j;
                            }
                        }
                        for (int k = 0; k < currMatchNum; k++) {
                            idx = selectedMatches[k];
                            red = targetInfoList.data[idx].redundancy;
                            matches[matchCnt] = {queryKmerList[j].info.sequenceID,
                                                 sspOrSp[red][targetInfoList.data[idx].sequenceID],
                                                 spTaxIdList[targetInfoList.data[idx].sequenceID],
                                                 genusTaxIdList[targetInfoList.data[idx].sequenceID],
                                                 queryKmerList[j].info.pos, selectedHammings[k],
                                                 selectedHammingSum[k], queryKmerList[j].info.frame};
                            matchCnt++;
                        }

                        if (currMatchNum != 0) {
                        }
                        continue;
                    }
                    selectedMatches.clear();
                    selectedHammingSum.clear();
                    selectedHammings.clear();

                    // Reuse the candidate target k-mers to compare in DNA level if queries are the same at amino acid level but not at DNA level
                    if (currentQueryAA == AminoAcidPart(queryKmerList[j].ADkmer)) {
                        compareDna(queryKmerList[j].ADkmer, candidateTargetKmers, startIdxOfAAmatch, selectedMatches,
                                   selectedHammingSum, selectedHammings, i);
                        currMatchNum = selectedMatches.size();

                        // If local buffer is full, copy them to the shared buffer.
                        if (matchCnt + currMatchNum > localBufferSize) {
                            // Check if the shared buffer is full.
                            posToWrite = matchBuffer.reserveMemory(matchCnt);
                            if (posToWrite + matchCnt >=
                                matchBuffer.bufferSize) { // full -> write matches to file first
                                hasOverflow = true;
                                querySplits[i].start = lastMovedQueryIdx + 1;
                                __sync_fetch_and_sub(& matchBuffer.startIndexOfReserve, matchCnt);
                                break;
                            } else { // not full -> copy matches to the shared buffer
                                moveMatches(matchBuffer.buffer + posToWrite, matches, matchCnt);
                                lastMovedQueryIdx = j;
                            }
                        }

                        for (int k = 0; k < currMatchNum; k++) {
                            idx = selectedMatches[k];
                            red = targetInfoList.data[idx].redundancy;
                            matches[matchCnt] = {queryKmerList[j].info.sequenceID,
                                                 sspOrSp[red][targetInfoList.data[idx].sequenceID],
                                                 spTaxIdList[targetInfoList.data[idx].sequenceID],
                                                 genusTaxIdList[targetInfoList.data[idx].sequenceID],
                                                 queryKmerList[j].info.pos, selectedHammings[k],
                                                 selectedHammingSum[k], queryKmerList[j].info.frame};
                            matchCnt++;
                        }
                        continue;
                    }
                    candidateTargetKmers.clear();

                    // Get next query, and start to find
                    currentQuery = queryKmerList[j].ADkmer;
                    currentQueryAA = AminoAcidPart(currentQuery);

                    // Skip target k-mers that are not matched in amino acid level
                    while (AminoAcidPart(currentQuery) > AminoAcidPart(currentTargetKmer) &&
                           (targetInfoIdx < numOfTargetKmer) && (diffIdxPos != numOfDiffIdx)) {
                        currentTargetKmer = getNextTargetKmer(currentTargetKmer, targetDiffIdxList.data, diffIdxPos);
                        targetInfoIdx++;
                    }

                    if (AminoAcidPart(currentQuery) !=
                            AminoAcidPart(currentTargetKmer)) // Move to next query k-mer if there isn't any match.
                        continue;
                    else
                        startIdxOfAAmatch = targetInfoIdx;

                    // Load target k-mers that are matched in amino acid level
                    while (AminoAcidPart(currentQuery) == AminoAcidPart(currentTargetKmer) &&
                           (targetInfoIdx < numOfTargetKmer) && (diffIdxPos != numOfDiffIdx)) {
                        candidateTargetKmers.push_back(currentTargetKmer);
                        currentTargetKmer = getNextTargetKmer(currentTargetKmer, targetDiffIdxList.data, diffIdxPos);
                        targetInfoIdx++;
                    }

                    // Compare the current query and the loaded target k-mers and select
                    compareDna(currentQuery, candidateTargetKmers, startIdxOfAAmatch, selectedMatches,
                               selectedHammingSum, selectedHammings, i);

                    // If local buffer is full, copy them to the shared buffer.
                    currMatchNum = selectedMatches.size();
                    if (matchCnt + currMatchNum > localBufferSize) {
                        // Check if the shared buffer is full.
                        posToWrite = matchBuffer.reserveMemory(matchCnt);
                        if (posToWrite + matchCnt >= matchBuffer.bufferSize) { // full -> write matches to file first
                            hasOverflow = true;
                            querySplits[i].start = lastMovedQueryIdx + 1;
                            __sync_fetch_and_sub(& matchBuffer.startIndexOfReserve, matchCnt);
                            break;
                        } else { // not full -> copy matches to the shared buffer
                            moveMatches(matchBuffer.buffer + posToWrite, matches, matchCnt);
                            lastMovedQueryIdx = j;
                        }
                    }

                    for (int k = 0; k < currMatchNum; k++) {
                        idx = selectedMatches[k];
                        red = targetInfoList.data[idx].redundancy;
                        matches[matchCnt] = {queryKmerList[j].info.sequenceID,
                                             sspOrSp[red][targetInfoList.data[idx].sequenceID],
                                             spTaxIdList[targetInfoList.data[idx].sequenceID],
                                             genusTaxIdList[targetInfoList.data[idx].sequenceID],
                                             queryKmerList[j].info.pos, selectedHammings[k],
                                             selectedHammingSum[k], queryKmerList[j].info.frame};
                        matchCnt++;
                    }
                } // End of one split

                // Move matches in the local buffer to the shared buffer
                posToWrite = matchBuffer.reserveMemory(matchCnt);
                if (posToWrite + matchCnt >= matchBuffer.bufferSize) {
                    hasOverflow = true;
                    querySplits[i].start = lastMovedQueryIdx + 1;
                    __sync_fetch_and_sub(& matchBuffer.startIndexOfReserve, matchCnt);
                } else {
                    moveMatches(matchBuffer.buffer + posToWrite, matches, matchCnt);
                }

                // Check whether current split is completed or not
                if (querySplits[i].start - 1 == querySplits[i].end) {
                    splitCheckList[i] = true;
                    __sync_fetch_and_add(&completedSplitCnt, 1);
                }
            } // End of omp for (Iterating for splits)
            delete[] matches;
        } // end of omp parallel
        if (hasOverflow) {
            cout << "overflow!!!" << endl;
            break;
        }
    } // end of while(completeSplitCnt < threadNum)
    cout << "Time spent for linearSearch: " << double(time(nullptr) - beforeSearch) << endl;

    munmap(targetDiffIdxList.data, targetDiffIdxList.fileSize + 1);
    munmap(targetInfoList.data, targetInfoList.fileSize + 1);
    munmap(diffIdxSplits.data, diffIdxSplits.fileSize + 1);
    free(splitCheckList);
    queryKmerCnt = 0;
}

void Classifier::moveMatches(Match *dest, Match *src, int &matchNum) {
    memcpy(dest, src, sizeof(Match) * matchNum);
    matchNum = 0;
}

// It compares query k-mers to target k-mers.
// If a query has matches, the matches with the smallest hamming distance will be selected
void Classifier::compareDna(uint64_t query, vector<uint64_t> &targetKmersToCompare, size_t startIdx,
                            vector<size_t> &selectedMatches, vector<uint8_t> &selectedHammingSum,
                            vector<uint16_t> &selectedHammings, int i2) {

    size_t size = targetKmersToCompare.size();
    uint8_t *hammingSums = new uint8_t[size + 1];
    uint8_t currentHammingSum;
    uint8_t minHammingSum = UINT8_MAX;

    // Calculate hamming distance
    for (size_t i = 0; i < size; i++) {
        currentHammingSum = getHammingDistanceSum(query, targetKmersToCompare[i]);
        if (currentHammingSum < minHammingSum) {
            minHammingSum = currentHammingSum;
        }
        hammingSums[i] = currentHammingSum;
    }

    // Select target k-mers that passed hamming criteria
    for (size_t h = 0; h < size; h++) {
        if (hammingSums[h] <= minHammingSum + hammingMargin) {
            selectedMatches.push_back(startIdx + h);
            selectedHammingSum.push_back(hammingSums[h]);
            selectedHammings.push_back(getHammings(query, targetKmersToCompare[h]));
        }
    }
    delete[] hammingSums;
}

// It analyses the result of linear search.
void Classifier::analyseResultParallel(NcbiTaxonomy &ncbiTaxonomy,
                                       Match *matchList,
                                       size_t numOfMatches,
                                       int seqNum,
                                       Query *queryList,
                                       const LocalParameters &par) {

    // Devide matches into blocks for multi threading
    cout << "Devide matches into blocks for multi threading" << endl;
    MatchBlock *matchBlocks = new MatchBlock[seqNum];
    size_t matchIdx = 0;
    size_t blockIdx = 0;
    uint32_t currentQuery;
    while (matchIdx < numOfMatches) {
        currentQuery = matchList[matchIdx].queryId;
        matchBlocks[blockIdx].id = currentQuery;
        matchBlocks[blockIdx].start = matchIdx;
        while ((currentQuery == matchList[matchIdx].queryId) && (matchIdx < numOfMatches)) ++matchIdx;
        matchBlocks[blockIdx].end = matchIdx - 1;
        blockIdx++;
    }

    if (PRINT) {
        omp_set_num_threads(1);
    } else {
        omp_set_num_threads(par.threads);
    }

    // Process each block
#pragma omp parallel default(none), shared(cout, matchBlocks, matchList, seqNum, ncbiTaxonomy, queryList, blockIdx, par)
    {
#pragma omp for schedule(dynamic, 1)
        for (size_t i = 0; i < blockIdx; ++i) {
            chooseBestTaxon(ncbiTaxonomy,
                            matchBlocks[i].id,
                            matchBlocks[i].start,
                            matchBlocks[i].end,
                            matchList,
                            queryList,
                            par);
        }
    }

    for (size_t i = 0; i < blockIdx; i++) {
        ++taxCounts[queryList[matchBlocks[i].id].classification];
    }
    delete[] matchBlocks;
    cout << "End of analyseResultParallel" << endl;
}


void Classifier::chooseBestTaxon(NcbiTaxonomy &ncbiTaxonomy, uint32_t currentQuery,
                                 size_t offset, size_t end, Match *matchList, Query *queryList,
                                 const LocalParameters &par) {
    int queryLength = queryList[currentQuery].queryLength; //queryList[currentQuery].queryLength; 13497
    TaxID selectedTaxon;
    if (PRINT) {
        cout << "# " << currentQuery << endl;
        for (size_t i = offset; i < end + 1; i++) {
            cout << matchList[i].genusTaxID << " " << matchList[i].speciesTaxID << " " << matchList[i].taxID << " " <<
                 matchList[i].position << " " << int(matchList[i].hamming) << endl;
        }
    }

    // Get the best genus for current query
    vector<Match> matchesForLCA;
    matchesForLCA.reserve(end - offset + 1);
    float highRankScore;
    int res;
    if (par.seqMode == 2) {
        res = getMatchesOfTheBestGenus_paired(matchesForLCA,
                                              matchList,
                                              end,
                                              offset,
                                              queryList[currentQuery].queryLength,
                                              queryList[currentQuery].queryLength2,
                                              highRankScore);
    } else {
        res = getMatchesOfTheBestGenus(matchesForLCA, matchList, end, offset, queryLength, highRankScore);
    }

    if (PRINT) {
        cout << "# " << currentQuery << " filtered" << endl;
        for (size_t i = 0; i < matchesForLCA.size(); i++) {
            cout << matchesForLCA[i].genusTaxID << " " << matchesForLCA[i].speciesTaxID << " " << matchesForLCA[i].taxID
                 << " " << matchesForLCA[i].position << " " << int(matchesForLCA[i].hamming)
                 << endl;
        }
    }

    // If there is no proper genus for current query, it is un-classified.
    if (res == 4) {
        queryList[currentQuery].isClassified = false;
        queryList[currentQuery].classification = 0;
        queryList[currentQuery].score = 0;
        queryList[currentQuery].newSpecies = false;
        return;
    }

    // If the score is too low, it is un-classified
    if (highRankScore < par.minScore) {
        queryList[currentQuery].isClassified = false;
        queryList[currentQuery].classification = 0;
        queryList[currentQuery].score = 0;
        queryList[currentQuery].newSpecies = false;
        return;
    }

    for (size_t i = 0; i < matchesForLCA.size(); i++) {
        queryList[currentQuery].taxCnt[matchesForLCA[i].taxID]++;
    }

    // If there are two or more good genus level candidates, find the LCA.
    if (res == 2) {
        vector<TaxID> taxIdList;
        taxIdList.reserve(matchesForLCA.size());
        for (size_t i = 0; i < matchesForLCA.size(); i++) {
            taxIdList.push_back(matchesForLCA[i].taxID);
        }
        selectedTaxon = ncbiTaxonomy.LCA(taxIdList)->taxId;
        queryList[currentQuery].isClassified = true;
        queryList[currentQuery].classification = selectedTaxon;
        queryList[currentQuery].score = highRankScore;
        if (PRINT) {
            cout << "# " << currentQuery << " " << res << endl;
            for (size_t i = 0; i < taxIdList.size(); i++) {
                cout << i << " " << matchesForLCA[i].position << " " <<
                     matchesForLCA[i].taxID << " " << int(matchesForLCA[i].hamming) << " "
                     << endl;
            }
            cout << "Score: " << highRankScore << " " << selectedTaxon << " "
                 << ncbiTaxonomy.taxonNode(selectedTaxon)->rank << endl;
        }
        return;
    }

    // Choose the species with the highest coverage.
    TaxID selectedSpecies;
    ScrCov speciesScrCov(0.f, 0.f);
    vector<TaxID> species;
    if (par.seqMode == 2) {
        classifyFurther_paired(matchesForLCA,
                               ncbiTaxonomy,
                               queryLength,
                               queryList[currentQuery].queryLength2,
                               speciesScrCov,
                               species);
    } else {
        chooseSpecies(matchesForLCA,
                      ncbiTaxonomy,
                      queryLength,
                      speciesScrCov,
                      species);
    }

    // Classify at the genus rank if more than one species are selected or the score at species level is not enough.
    if (species.size() > 1
        || (speciesScrCov.score < minSpScore && !ncbiTaxonomy.IsAncestor(par.virusTaxId, matchesForLCA[0].taxID))) {
        queryList[currentQuery].isClassified = true;
        queryList[currentQuery].classification = ncbiTaxonomy.getTaxIdAtRank(matchesForLCA[0].taxID, "genus");
        queryList[currentQuery].score = highRankScore;
        return;
    }

    selectedSpecies = species[0];

    // Check if it can be classified at the subspecies rank.
    int numOfstrains = 0;
    TaxID strainID = 0;
    int count = 1;
    int minStrainSpecificCnt = 1;
    if (par.seqMode == 1) {
        minStrainSpecificCnt = 1;
    } else if (par.seqMode == 2) {
        minStrainSpecificCnt = 2;
    } else if (par.seqMode == 3) {
        minStrainSpecificCnt = 3;
        if (queryLength > 3000) {
            minStrainSpecificCnt = queryLength / 1000;
        }
    }
    if (NcbiTaxonomy::findRankIndex(ncbiTaxonomy.taxonNode(selectedSpecies)->rank) == 4) {
        unordered_map<TaxID, int> strainMatchCnt;
        for (size_t i = 0; i < matchesForLCA.size(); i++) {
            if (selectedSpecies != matchesForLCA[i].taxID
                && ncbiTaxonomy.IsAncestor(selectedSpecies, matchesForLCA[i].taxID)) {
                strainMatchCnt[matchesForLCA[i].taxID]++;
            }
        }
        for (auto strainIt = strainMatchCnt.begin(); strainIt != strainMatchCnt.end(); strainIt++) {
            if (strainIt->second > minStrainSpecificCnt) {
                strainID = strainIt->first;
                numOfstrains++;
                count = strainIt->second;
            }
        }
        if (numOfstrains == 1 && count > minStrainSpecificCnt + 1) {
            selectedSpecies = strainID;
        }
    }

    // Store classification results
    queryList[currentQuery].isClassified = true;
    queryList[currentQuery].classification = selectedSpecies;
    queryList[currentQuery].score = speciesScrCov.score;
    queryList[currentQuery].newSpecies = false;

    if (PRINT) {
        cout << "# " << currentQuery << endl;
        for (size_t i = 0; i < matchesForLCA.size(); i++) {
            cout << i << " " << matchesForLCA[i].position << " " <<
                 matchesForLCA[i].taxID << " " << int(matchesForLCA[i].hamming) << endl;
        }
        cout << "Score: " << speciesScrCov.score << "  " << selectedSpecies << " "
             << ncbiTaxonomy.taxonNode(selectedSpecies)->rank
             << endl;
    }
}

int Classifier::getMatchesOfTheBestGenus_paired(vector<Match> &matchesForMajorityLCA, Match *matchList, size_t end,
                                                size_t offset, int readLength1, int readLength2, float &bestScore) {
    int conCnt;
    uint32_t hammingSum;
    float hammingMean;
    TaxID currentGenus;
    TaxID currentSpecies;

    vector<Match> tempMatchContainer;
    vector<Match> filteredMatches;
    vector<vector<Match>> matchesForEachGenus;
    vector<bool> conservedWithinGenus;
    vector<float> scoreOfEachGenus;
    size_t i = offset;
    size_t offsetIdx;
    bool newOffset;
    bool lastIn;
    size_t speciesMatchCnt;
    size_t speciesDiffPosCnt;
    size_t consecutiveCnt;
    int lastPos;
    uint8_t curFrame;
    while (i < end + 1) {
        currentGenus = matchList[i].genusTaxID;
        // For current genus
        while (currentGenus == matchList[i].genusTaxID && (i < end + 1)) {
            currentSpecies = matchList[i].speciesTaxID;
            // For current species
            // Filter un-consecutive matches (probably random matches)
            speciesMatchCnt = 0;
            speciesDiffPosCnt = 0;
            consecutiveCnt = 0;
            lastPos = -1;
            lastIn = false;
            while (currentSpecies == matchList[i + 1].speciesTaxID && (i < end + 1)) {
                if (matchList[i].position + 3 >= matchList[i + 1].position) {
                    //filteredMatches.push_back(matchList[i]);
                    tempMatchContainer.push_back(matchList[i]);
                    speciesMatchCnt++;
                    if (matchList[i].position / 3 != lastPos) {
                        lastPos = matchList[i].position / 3;
                        speciesDiffPosCnt++;
                        consecutiveCnt++;
                    }
                    lastIn = true;
                } else if (lastIn) {
                    lastIn = false;
                    //filteredMatches.push_back(matchList[i]);
                    tempMatchContainer.push_back(matchList[i]);
                    speciesMatchCnt++;
                    if (matchList[i].position / 3 != lastPos) {
                        lastPos = matchList[i].position / 3;
                        speciesDiffPosCnt++;
                        consecutiveCnt++;
                    }
//                    if (consecutiveCnt < minConsCnt) {
//                        for (size_t j = 0; j < speciesMatchCnt; j++) {
//                            filteredMatches.pop_back();
//                        }
//                    }
                    if (consecutiveCnt >= minConsCnt) {
                        filteredMatches.insert(filteredMatches.end(), tempMatchContainer.begin(), tempMatchContainer.end());
                    }
                    consecutiveCnt = 0;
                    speciesMatchCnt = 0;
                    tempMatchContainer.clear();
                }
                i++;
            }
            if (lastIn) {
                //filteredMatches.push_back(matchList[i]);
                tempMatchContainer.push_back(matchList[i]);
                speciesMatchCnt++;
                if (matchList[i].position / 3 != lastPos) {
                    lastPos = matchList[i].position / 3;
                    speciesDiffPosCnt++;
                    consecutiveCnt++;
                }
                //                    if (consecutiveCnt < minConsCnt) {
//                        for (size_t j = 0; j < speciesMatchCnt; j++) {
//                            filteredMatches.pop_back();
//                        }
//                    }
                if (consecutiveCnt >= minConsCnt) {
                    filteredMatches.insert(filteredMatches.end(), tempMatchContainer.begin(), tempMatchContainer.end());
                }
                tempMatchContainer.clear();
            }
            i++;
        }
        // Construct a match combination using filtered matches of current genus
        // so that it can best cover the query, and score the combination
        if (!filteredMatches.empty()) {
            constructMatchCombination_paired(filteredMatches, matchesForEachGenus, scoreOfEachGenus, readLength1,
                                             readLength2);
        }
        filteredMatches.clear();
    }

    // If there are no meaningful genus
    if (scoreOfEachGenus.empty()) {
        bestScore = 0;
        return 4;
    }

    float maxScore = *max_element(scoreOfEachGenus.begin(), scoreOfEachGenus.end());

    vector<size_t> maxIdx;
    for (size_t g = 0; g < scoreOfEachGenus.size(); g++) {
        if (scoreOfEachGenus[g] > maxScore * 0.95f) {
            maxIdx.push_back(g);
        }
    }
    bestScore = maxScore;

    for (size_t g = 0; g < maxIdx.size(); g++) {
        matchesForMajorityLCA.insert(matchesForMajorityLCA.end(), matchesForEachGenus[maxIdx[g]].begin(),
                                     matchesForEachGenus[maxIdx[g]].end());
    }

    if (maxIdx.size() > 1) {
        return 2;
    }
    return 1;

    //Three cases
    //1. one genus
    //2. more than one genus
    //3. no genus
}

int Classifier::getMatchesOfTheBestGenus(vector<Match> &matchesForMajorityLCA, Match *matchList, size_t end,
                                         size_t offset, int queryLength, float &bestScore) {
    int conCnt;
    uint32_t hammingSum;
    float hammingMean;
    TaxID currentGenus;
    TaxID currentSpecies;

    vector<Match> filteredMatches;
    vector<vector<Match>> matchesForEachGenus;
    vector<bool> conservedWithinGenus;
    vector<float> scoreOfEachGenus;
    size_t i = offset;
    size_t offsetIdx;
    bool newOffset;
    bool lastIn;
    size_t speciesMatchCnt;
    size_t speciesDiffPosCnt;
    size_t consecutiveCnt;
    int lastPos;
    uint8_t curFrame;
    while (i < end + 1) {
        currentGenus = matchList[i].genusTaxID;
        // For current genus
        while ((i < end + 1) && currentGenus == matchList[i].genusTaxID) {
            currentSpecies = matchList[i].speciesTaxID;
            // For current species
            // Filter un-consecutive matches (probably random matches)
            speciesMatchCnt = 0;
            speciesDiffPosCnt = 0;
            consecutiveCnt = 0;
            lastPos = -1;
            lastIn = false;
            while (currentSpecies == matchList[i + 1].speciesTaxID && (i < end + 1)) {
                if (matchList[i].position + 3 >= matchList[i + 1].position) {
                    filteredMatches.push_back(matchList[i]);
                    speciesMatchCnt++;
                    if (matchList[i].position / 3 != lastPos) {
                        lastPos = matchList[i].position / 3;
                        speciesDiffPosCnt++;
                        consecutiveCnt++;
                    }
                    lastIn = true;
                } else if (lastIn) {
                    lastIn = false;
                    filteredMatches.push_back(matchList[i]);
                    speciesMatchCnt++;
                    if (matchList[i].position / 3 != lastPos) {
                        lastPos = matchList[i].position / 3;
                        speciesDiffPosCnt++;
                        consecutiveCnt++;
                    }
                    if (consecutiveCnt < minConsCnt) {
                        for (size_t j = 0; j < speciesMatchCnt; j++) {
                            filteredMatches.pop_back();
                        }
                    }
                    consecutiveCnt = 0;
                    speciesMatchCnt = 0;
                }
                i++;
            }
            if (lastIn) {
                filteredMatches.push_back(matchList[i]);
                speciesMatchCnt++;
                if (matchList[i].position / 3 != lastPos) {
                    lastPos = matchList[i].position / 3;
                    speciesDiffPosCnt++;
                    consecutiveCnt++;
                }
                if (consecutiveCnt < minConsCnt) {
                    for (size_t j = 0; j < speciesMatchCnt; j++) {
                        filteredMatches.pop_back();
                    }
                }
            }
            i++;
        }

        // Construct a match combination using filtered matches of current genus
        // so that it can best cover the query, and score the combination
        if (!filteredMatches.empty()) {
            constructMatchCombination(filteredMatches, matchesForEachGenus, scoreOfEachGenus, queryLength);
        }
        filteredMatches.clear();
    }

    // If there are no meaningful genus
    if (scoreOfEachGenus.empty()) {
        bestScore = 0;
        return 3;
    }

    float maxScore = *max_element(scoreOfEachGenus.begin(), scoreOfEachGenus.end());

    vector<size_t> maxIdx;
    for (size_t g = 0; g < scoreOfEachGenus.size(); g++) {
        if (scoreOfEachGenus[g] > maxScore * 0.95f) {
            maxIdx.push_back(g);
        }
    }
    bestScore = maxScore;

    for (size_t g = 0; g < maxIdx.size(); g++) {
        matchesForMajorityLCA.insert(matchesForMajorityLCA.end(), matchesForEachGenus[maxIdx[g]].begin(),
                                     matchesForEachGenus[maxIdx[g]].end());
    }

    if (maxIdx.size() > 1) {
        return 2;
    }
    return 1;

    //Three cases
    //1. one genus
    //2. more than one genus
    //3. no genus
}

void Classifier::constructMatchCombination(vector<Match> &filteredMatches,
                                           vector<vector<Match>> &matchesForEachGenus,
                                           vector<float> &scoreOfEachGenus,
                                           int queryLength) {
    // Do not allow overlaps between the same species
    vector<Match> matches;
    size_t walker = 0;
    size_t numOfFitMat = filteredMatches.size();
    Match &currentMatch = filteredMatches[0];
    while (walker < numOfFitMat) {
        TaxID currentSpecies = filteredMatches[walker].speciesTaxID;
        int currentPosition = filteredMatches[walker].position / 3;
        currentMatch = filteredMatches[walker];
        // Look through overlaps within a species
        while (walker < numOfFitMat
               && filteredMatches[walker].speciesTaxID == currentSpecies
               && filteredMatches[walker].position / 3 == currentPosition) {
            // Take the match with lower hamming distance
            if (filteredMatches[walker].hamming < currentMatch.hamming) {
                currentMatch = filteredMatches[walker];
            }
                // Overlapping with same hamming distance but different subspecies taxonomy ID -> species level
            else if (currentMatch.taxID != filteredMatches[walker].taxID &&
                     currentMatch.hamming == filteredMatches[walker].hamming) {
                currentMatch.taxID = currentMatch.speciesTaxID;
            }
            walker++;
        }
        matches.push_back(currentMatch);
    }

    // Calculate Hamming distance & covered length
    int coveredPosCnt = 0;
    uint16_t currHammings;
    int aminoAcidNum = (int) queryLength / 3;
    int currPos;
    size_t matchNum = matches.size();
    size_t f = 0;

    // Get the largest hamming distance at each position of query
    auto *hammingsAtEachPos = new signed char[aminoAcidNum + 1];
    memset(hammingsAtEachPos, -1, (aminoAcidNum + 1));
    while (f < matchNum) {
        currPos = matches[f].position / 3;
        currHammings = matches[f].rightEndHamming;
        if (GET_2_BITS(currHammings) > hammingsAtEachPos[currPos])
            hammingsAtEachPos[currPos + unmaskedPos[0]] = GET_2_BITS(currHammings);
        if (GET_2_BITS(currHammings >> 2) > hammingsAtEachPos[currPos + 1])
            hammingsAtEachPos[currPos + unmaskedPos[1]] = GET_2_BITS(currHammings >> 2);
        if (GET_2_BITS(currHammings >> 4) > hammingsAtEachPos[currPos + 2])
            hammingsAtEachPos[currPos + unmaskedPos[2]] = GET_2_BITS(currHammings >> 4);
        if (GET_2_BITS(currHammings >> 6) > hammingsAtEachPos[currPos + 3])
            hammingsAtEachPos[currPos + unmaskedPos[3]] = GET_2_BITS(currHammings >> 6);
        if (GET_2_BITS(currHammings >> 8) > hammingsAtEachPos[currPos + 4])
            hammingsAtEachPos[currPos + unmaskedPos[4]] = GET_2_BITS(currHammings >> 8);
        if (GET_2_BITS(currHammings >> 10) > hammingsAtEachPos[currPos + 5])
            hammingsAtEachPos[currPos + unmaskedPos[5]] = GET_2_BITS(currHammings >> 10);
        if (GET_2_BITS(currHammings >> 12) > hammingsAtEachPos[currPos + 6])
            hammingsAtEachPos[currPos + unmaskedPos[6]] = GET_2_BITS(currHammings >> 12);
        if (GET_2_BITS(currHammings >> 14) > hammingsAtEachPos[currPos + 7])
            hammingsAtEachPos[currPos + unmaskedPos[7]] = GET_2_BITS(currHammings >> 14);
        f++;
    }

    // Sum up hamming distances and count the number of position covered by the matches.
    float hammingSum = 0;
    for (int h = 0; h < aminoAcidNum; h++) {
        if (hammingsAtEachPos[h] == 0) { // Add 0 for 0 hamming dist.
            coveredPosCnt++;
        } else if (hammingsAtEachPos[h] != -1) { // Add 1.5, 2, 2.5 for 1, 2, 3 hamming dist. respectively
            hammingSum += 1.0f + (0.5f * hammingsAtEachPos[h]);
            coveredPosCnt++;
        }
    }
    delete[] hammingsAtEachPos;

    // Ignore too few matches
    if (coveredPosCnt < 2) return;

    // Score current genus
    int coveredLength = coveredPosCnt * 3;
    if (coveredLength > queryLength) coveredLength = queryLength;
    float score = ((float) coveredLength - hammingSum) / (float) queryLength;
    scoreOfEachGenus.push_back(score);
    matchesForEachGenus.push_back(matches);

    if (PRINT) {
        cout << filteredMatches[0].genusTaxID << " " << coveredLength << " " << hammingSum << " "
             << ((float) coveredLength - hammingSum) / (float) queryLength <<
             " " << matches.size()
             << endl;
    }
}

void Classifier::constructMatchCombination_paired(vector<Match> &filteredMatches,
                                                  vector<vector<Match>> &matchesForEachGenus,
                                                  vector<float> &scoreOfEachGenus,
                                                  int readLength1,
                                                  int readLength2) {
    // Do not allow overlaps between the same species
    vector<Match> matches;
    size_t walker = 0;
    size_t numOfFitMat = filteredMatches.size();
    Match &currentMatch = filteredMatches[0];
    while (walker < numOfFitMat) {
        TaxID currentSpecies = filteredMatches[walker].speciesTaxID;
        int currentPosition = filteredMatches[walker].position / 3;
        currentMatch = filteredMatches[walker];
        // Look through overlaps within a species
        while (walker < numOfFitMat && filteredMatches[walker].speciesTaxID == currentSpecies
               && filteredMatches[walker].position / 3 == currentPosition
               ) {
            // Take the match with lower hamming distance
            if (filteredMatches[walker].hamming < currentMatch.hamming) {
                currentMatch = filteredMatches[walker];
            }
                // Overlapping with same hamming distance but different subspecies taxonomy ID -> species level
            else if (currentMatch.taxID != filteredMatches[walker].taxID &&
                     currentMatch.hamming == filteredMatches[walker].hamming) {
                currentMatch.taxID = currentMatch.speciesTaxID;
            }
            walker++;
        }
        matches.push_back(currentMatch);
    }

    // Calculate Hamming distance & covered length
    int coveredPosCnt_read1 = 0;
    int coveredPosCnt_read2 = 0;
    uint16_t currHammings;
    int aminoAcidNum_total = ((int) readLength1 / 3) + ((int) readLength2 / 3);
    int aminoAcidNum_read1 = ((int) readLength1 / 3);
    int currPos;
    size_t matchNum = matches.size();
    size_t f = 0;

    // Get the largest hamming distance at each position of query
    auto *hammingsAtEachPos = new signed char[aminoAcidNum_total + 3];
    memset(hammingsAtEachPos, -1, (aminoAcidNum_total + 3));
    while (f < matchNum) {
        currPos = matches[f].position / 3;
        currHammings = matches[f].rightEndHamming;
        if (GET_2_BITS(currHammings) > hammingsAtEachPos[currPos])
            hammingsAtEachPos[currPos + unmaskedPos[0]] = GET_2_BITS(currHammings);
        if (GET_2_BITS(currHammings >> 2) > hammingsAtEachPos[currPos + 1])
            hammingsAtEachPos[currPos + unmaskedPos[1]] = GET_2_BITS(currHammings >> 2);
        if (GET_2_BITS(currHammings >> 4) > hammingsAtEachPos[currPos + 2])
            hammingsAtEachPos[currPos + unmaskedPos[2]] = GET_2_BITS(currHammings >> 4);
        if (GET_2_BITS(currHammings >> 6) > hammingsAtEachPos[currPos + 3])
            hammingsAtEachPos[currPos + unmaskedPos[3]] = GET_2_BITS(currHammings >> 6);
        if (GET_2_BITS(currHammings >> 8) > hammingsAtEachPos[currPos + 4])
            hammingsAtEachPos[currPos + unmaskedPos[4]] = GET_2_BITS(currHammings >> 8);
        if (GET_2_BITS(currHammings >> 10) > hammingsAtEachPos[currPos + 5])
            hammingsAtEachPos[currPos + unmaskedPos[5]] = GET_2_BITS(currHammings >> 10);
        if (GET_2_BITS(currHammings >> 12) > hammingsAtEachPos[currPos + 6])
            hammingsAtEachPos[currPos + unmaskedPos[6]] = GET_2_BITS(currHammings >> 12);
        if (GET_2_BITS(currHammings >> 14) > hammingsAtEachPos[currPos + 7])
            hammingsAtEachPos[currPos + unmaskedPos[7]] = GET_2_BITS(currHammings >> 14);
        f++;
    }

    // Sum up hamming distances and count the number of position covered by the matches.
    float hammingSum = 0;
    for (int h = 0; h < aminoAcidNum_total; h++) {
        // Read 1
        if (h < aminoAcidNum_read1) {
            if (hammingsAtEachPos[h] == 0) { // Add 0 for 0 hamming dist.
                coveredPosCnt_read1++;
            } else if (hammingsAtEachPos[h] != -1) { // Add 1.5, 2, 2.5 for 1, 2, 3 hamming dist. respectively
                hammingSum += 1.0f + (0.5f * hammingsAtEachPos[h]);
                coveredPosCnt_read1++;
            }
        }
            // Read 2
        else {
            if (hammingsAtEachPos[h] == 0) { // Add 0 for 0 hamming dist.
                coveredPosCnt_read2++;
            } else if (hammingsAtEachPos[h] != -1) { // Add 1.5, 2, 2.5 for 1, 2, 3 hamming dist. respectively
                hammingSum += 1.0f + (0.5f * hammingsAtEachPos[h]);
                coveredPosCnt_read2++;
            }
        }
    }
    delete[] hammingsAtEachPos;

    // Ignore too few matches
    if (coveredPosCnt_read1 + coveredPosCnt_read2 < 2) return;

    // Score current genus
    int coveredLength_read1 = coveredPosCnt_read1 * 3;
    int coveredLength_read2 = coveredPosCnt_read2 * 3;
    if (coveredLength_read1 > readLength1) coveredLength_read1 = readLength1;
    if (coveredLength_read2 > readLength2) coveredLength_read2 = readLength2;
    float score =
            ((float) (coveredLength_read1 + coveredLength_read2) - hammingSum) / (float) (readLength1 + readLength2);
    scoreOfEachGenus.push_back(score);
    matchesForEachGenus.push_back(matches);

    if (PRINT) {
        cout << filteredMatches[0].genusTaxID << " " << coveredLength_read1 + coveredLength_read2 << " " << hammingSum
             << " " << score <<
             " " << matches.size()
             << endl;
    }
}

TaxID Classifier::classifyFurther2(const vector<Match> &matches,
                                   NcbiTaxonomy &taxonomy,
                                   float possibleKmerNum) {

    std::unordered_map<TaxID, int> taxIdCounts;
    float majorityCutoff = 0.8;
    float coverageThreshold = 0.8;

    // Subspecies && Species
    for (Match match: matches) {
        taxIdCounts[match.taxID] += 1;
    }

    // Subspecies -> Species
    for (Match match: matches) {
        if (taxIdCounts[match.taxID] > 1) {
            taxIdCounts[match.speciesTaxID] += (match.speciesTaxID != match.taxID);
        }
    }

    // Max count
    int maxCnt2 = (*max_element(taxIdCounts.begin(), taxIdCounts.end(),
                                [](const pair<TaxID, int> &p1, const pair<TaxID, int> &p2) {
                                    return p1.second < p2.second;
                                })
    ).second;
    if (maxCnt2 > (int) possibleKmerNum) { maxCnt2 = (int) possibleKmerNum; }

    float currentCoverage;
    float currnetPercentage;
    bool haveMetCovThr = false;
    bool haveMetMajorityThr = false;
    size_t matchNum = matches.size();
    TaxID currRank;
    vector<TaxID> ties;
    int minRank = INT_MAX;
    float selectedPercent = 0;
    TaxID selectedTaxon;
    for (auto it = taxIdCounts.begin(); it != taxIdCounts.end(); it++) {
        currnetPercentage = (float) it->second / matchNum;

        if (it->second >= maxCnt2) {
            it->second = maxCnt2 - 1;
        }
        currentCoverage = (float) it->second / (possibleKmerNum - 1);

        currRank = NcbiTaxonomy::findRankIndex(taxonomy.taxonNode(it->first)->rank);

        if (currentCoverage > coverageThreshold && (it->second >= maxCnt2 - 1)) {
            haveMetCovThr = true;
            ties.push_back(it->first);
        } else if (currnetPercentage >= majorityCutoff && (!haveMetCovThr)) {
            haveMetMajorityThr = true;
            if ((currRank < minRank) || ((currRank == minRank) && (currnetPercentage > selectedPercent))) {
                selectedTaxon = it->first;
                minRank = currRank;
                selectedPercent = currnetPercentage;
            }
        }
    }


    if (haveMetCovThr) {
        if (ties.size() > 1) {
            return taxonomy.LCA(ties)->taxId;
        } else {
            return ties[0];
        }
    } else if (haveMetMajorityThr) {
        return selectedTaxon;
    }
    return matches[0].genusTaxID;
}

void Classifier::chooseSpecies(const vector<Match> &matches,
                               NcbiTaxonomy &taxonomy,
                               int queryLength,
                               ScrCov &speciesScrCov,
                               vector<TaxID> &species) {
    // Score each species
    std::unordered_map<TaxID, ScrCov> speciesScrCovs;
    size_t i = 0;
    TaxID currentSpeices;
    size_t numOfMatch = matches.size();
    size_t speciesBegin, speciesEnd;
    float currScore, currCoverage;
    while (i < numOfMatch) {
        currentSpeices = matches[i].speciesTaxID;
        speciesBegin = i;
        while (currentSpeices == matches[i].speciesTaxID && (i < numOfMatch)) {
            i++;
        }
        speciesEnd = i;
        speciesScrCovs[currentSpeices] = scoreTaxon(matches, speciesBegin, speciesEnd, queryLength);
    }

    // Get the best species
    float bestCoverage = 0.f;
    for (auto sp = speciesScrCovs.begin(); sp != speciesScrCovs.end(); sp++) {
        if (sp->second.coverage > bestCoverage) {
            species.clear();
            species.push_back(sp->first);
            bestCoverage = sp->second.coverage;
            speciesScrCov.coverage = sp->second.coverage;
            speciesScrCov.score = sp->second.score;
        } else if (sp->second.coverage == bestCoverage) {
            species.push_back(sp->first);
        }
    }
}

void Classifier::classifyFurther_paired(const std::vector<Match> &matches,
                                        NcbiTaxonomy &taxonomy,
                                        int read1Length,
                                        int read2Length,
                                        ScrCov &speciesScrCov,
                                        vector<TaxID> &species) {
    // Score each species
    std::unordered_map<TaxID, ScrCov> speciesScrCovs;
    size_t i = 0;
    TaxID currentSpeices;
    size_t numOfMatch = matches.size();
    size_t speciesBegin, speciesEnd;
    while (i < numOfMatch) {
        currentSpeices = matches[i].speciesTaxID;
        speciesBegin = i;
        while (currentSpeices == matches[i].speciesTaxID && (i < numOfMatch)) {
            i++;
        }
        speciesEnd = i;
        speciesScrCovs[currentSpeices] = scoreTaxon_paired(matches, speciesBegin, speciesEnd, read1Length, read2Length);
    }

    // Get the best species
    float bestCovergae = 0.f;
    for (auto sp = speciesScrCovs.begin(); sp != speciesScrCovs.end(); sp++) {
        if (sp->second.coverage > bestCovergae) {
            species.clear();
            species.push_back(sp->first);
            bestCovergae = sp->second.coverage;
            speciesScrCov.coverage = sp->second.coverage;
            speciesScrCov.score = sp->second.score;
        } else if (sp->second.coverage == bestCovergae) {
            species.push_back(sp->first);
        }
    }
}

Classifier::ScrCov Classifier::scoreTaxon(const vector<Match> &matches,
                                          size_t begin,
                                          size_t end,
                                          int queryLength) {
    // Get the largest hamming distance at each position of query
    int aminoAcidNum = queryLength / 3;
    auto *hammingsAtEachPos = new signed char[aminoAcidNum + 1];
    memset(hammingsAtEachPos, -1, (aminoAcidNum + 1));
    int currPos;
    size_t walker = begin;
    uint16_t currHammings;
    while (walker < end) {
        currPos = matches[walker].position / 3;
        currHammings = matches[walker].rightEndHamming;
        if (GET_2_BITS(currHammings) > hammingsAtEachPos[currPos])
            hammingsAtEachPos[currPos + unmaskedPos[0]] = GET_2_BITS(currHammings);
        if (GET_2_BITS(currHammings >> 2) > hammingsAtEachPos[currPos + 1])
            hammingsAtEachPos[currPos + unmaskedPos[1]] = GET_2_BITS(currHammings >> 2);
        if (GET_2_BITS(currHammings >> 4) > hammingsAtEachPos[currPos + 2])
            hammingsAtEachPos[currPos + unmaskedPos[2]] = GET_2_BITS(currHammings >> 4);
        if (GET_2_BITS(currHammings >> 6) > hammingsAtEachPos[currPos + 3])
            hammingsAtEachPos[currPos + unmaskedPos[3]] = GET_2_BITS(currHammings >> 6);
        if (GET_2_BITS(currHammings >> 8) > hammingsAtEachPos[currPos + 4])
            hammingsAtEachPos[currPos + unmaskedPos[4]] = GET_2_BITS(currHammings >> 8);
        if (GET_2_BITS(currHammings >> 10) > hammingsAtEachPos[currPos + 5])
            hammingsAtEachPos[currPos + unmaskedPos[5]] = GET_2_BITS(currHammings >> 10);
        if (GET_2_BITS(currHammings >> 12) > hammingsAtEachPos[currPos + 6])
            hammingsAtEachPos[currPos + unmaskedPos[6]] = GET_2_BITS(currHammings >> 12);
        if (GET_2_BITS(currHammings >> 14) > hammingsAtEachPos[currPos + 7])
            hammingsAtEachPos[currPos + unmaskedPos[7]] = GET_2_BITS(currHammings >> 14);
        walker++;
    }

    // Sum up hamming distances and count the number of position covered by the matches.
    float hammingSum = 0;
    int coveredPosCnt = 0;
    for (int h = 0; h < aminoAcidNum; h++) {
        if (hammingsAtEachPos[h] == 0) { // Add 0 for 0 hamming dist.
            coveredPosCnt++;
        } else if (hammingsAtEachPos[h] != -1) { // Add 1.5, 2, 2.5 for 1, 2, 3 hamming dist. respectively
            hammingSum += 1.0f + (0.5f * hammingsAtEachPos[h]);
            coveredPosCnt++;
        }
    }
    delete[] hammingsAtEachPos;
    // Score
    int coveredLength = coveredPosCnt * 3;
    if (coveredLength >= queryLength) coveredLength = queryLength;
    return {((float) coveredLength - hammingSum) / (float) queryLength, (float) coveredLength / (float) queryLength};
}

Classifier::ScrCov
Classifier::scoreTaxon_paired(const vector<Match> &matches, size_t begin, size_t end, int queryLength,
                              int queryLength2) {

    // Get the largest hamming distance at each position of query
    int aminoAcidNum_total = queryLength / 3 + queryLength2 / 3;
    int aminoAcidNum_read1 = queryLength / 3;
    auto *hammingsAtEachPos = new signed char[aminoAcidNum_total + 3];
    memset(hammingsAtEachPos, -1, (aminoAcidNum_total + 3));

    int currPos;
    size_t walker = begin;
    uint16_t currHammings;
    while (walker < end) {
        currPos = matches[walker].position / 3;
        currHammings = matches[walker].rightEndHamming;
        if (GET_2_BITS(currHammings) > hammingsAtEachPos[currPos])
            hammingsAtEachPos[currPos + unmaskedPos[0]] = GET_2_BITS(currHammings);
        if (GET_2_BITS(currHammings >> 2) > hammingsAtEachPos[currPos + 1])
            hammingsAtEachPos[currPos + unmaskedPos[1]] = GET_2_BITS(currHammings >> 2);
        if (GET_2_BITS(currHammings >> 4) > hammingsAtEachPos[currPos + 2])
            hammingsAtEachPos[currPos + unmaskedPos[2]] = GET_2_BITS(currHammings >> 4);
        if (GET_2_BITS(currHammings >> 6) > hammingsAtEachPos[currPos + 3])
            hammingsAtEachPos[currPos + unmaskedPos[3]] = GET_2_BITS(currHammings >> 6);
        if (GET_2_BITS(currHammings >> 8) > hammingsAtEachPos[currPos + 4])
            hammingsAtEachPos[currPos + unmaskedPos[4]] = GET_2_BITS(currHammings >> 8);
        if (GET_2_BITS(currHammings >> 10) > hammingsAtEachPos[currPos + 5])
            hammingsAtEachPos[currPos + unmaskedPos[5]] = GET_2_BITS(currHammings >> 10);
        if (GET_2_BITS(currHammings >> 12) > hammingsAtEachPos[currPos + 6])
            hammingsAtEachPos[currPos + unmaskedPos[6]] = GET_2_BITS(currHammings >> 12);
        if (GET_2_BITS(currHammings >> 14) > hammingsAtEachPos[currPos + 7])
            hammingsAtEachPos[currPos + unmaskedPos[7]] = GET_2_BITS(currHammings >> 14);
        walker++;
    }

    // Sum up hamming distances and count the number of position covered by the matches.
    float hammingSum = 0;
    int coveredPosCnt_read1 = 0;
    int coveredPosCnt_read2 = 0;
    for (int h = 0; h < aminoAcidNum_total; h++) {
        // Read 1
        if (h < aminoAcidNum_read1) {
            if (hammingsAtEachPos[h] == 0) { // Add 0 for 0 hamming dist.
                coveredPosCnt_read1++;
            } else if (hammingsAtEachPos[h] != -1) { // Add 1.5, 2, 2.5 for 1, 2, 3 hamming dist. respectively
                hammingSum += 1.0f + (0.5f * hammingsAtEachPos[h]);
                coveredPosCnt_read1++;
            }
        }
        // Read 2
        else {
            if (hammingsAtEachPos[h] == 0) { // Add 0 for 0 hamming dist.
                coveredPosCnt_read2++;
            } else if (hammingsAtEachPos[h] != -1) { // Add 1.5, 2, 2.5 for 1, 2, 3 hamming dist. respectively
                hammingSum += 1.0f + (0.5f * hammingsAtEachPos[h]);
                coveredPosCnt_read2++;
            }
        }
    }
    delete[] hammingsAtEachPos;

    // Score
    hammingSum = 0;
    int coveredLength_read1 = coveredPosCnt_read1 * 3;
    int coveredLength_read2 = coveredPosCnt_read2 * 3;
    if (coveredLength_read1 >= queryLength) coveredLength_read1 = queryLength;
    if (coveredLength_read2 >= queryLength2) coveredLength_read2 = queryLength2;

    return {((float) coveredLength_read1 + coveredLength_read2 - hammingSum) / ((float) queryLength + queryLength2),
            ((float) coveredLength_read1 + coveredLength_read2) / ((float) queryLength + queryLength2)};
}

bool Classifier::compareForLinearSearch(const QueryKmer &a, const QueryKmer &b) {
    if (a.ADkmer < b.ADkmer) {
        return true;
    } else if (a.ADkmer == b.ADkmer) {
        return (a.info.sequenceID < b.info.sequenceID);
    }
    return false;
}

bool Classifier::sortByGenusAndSpecies2(const Match &a, const Match &b) {
    if (a.queryId < b.queryId) return true;
    else if (a.queryId == b.queryId) {
        if (a.genusTaxID < b.genusTaxID) return true;
        else if (a.genusTaxID == b.genusTaxID) {
            if (a.speciesTaxID < b.speciesTaxID) return true;
            else if (a.speciesTaxID == b.speciesTaxID) {
                if (a.position < b.position) return true;
                else if (a.position == b.position) {
                    return a.hamming < b.hamming;
                }
            }
        }
    }
    return false;
}

bool Classifier::sortMatchesByPos(const Match &a, const Match &b) {
    if (a.position / 3 < b.position / 3) return true;
    else if (a.position / 3 == b.position / 3) {
        if (a.speciesTaxID < b.speciesTaxID) return true;
        else if (a.speciesTaxID == b.speciesTaxID) {
            return a.hamming < b.hamming;
        }
    }
    return false;
}


void Classifier::writeReadClassification(Query *queryList, int queryNum, ofstream &readClassificationFile) {
    for (int i = 0; i < queryNum; i++) {
        readClassificationFile << queryList[i].isClassified << "\t" << queryList[i].name << "\t"
                               << queryList[i].classification << "\t"
                               << queryList[i].queryLength + queryList[i].queryLength2 << "\t"
                               << queryList[i].score << "\t";
        for (auto it = queryList[i].taxCnt.begin(); it != queryList[i].taxCnt.end(); ++it) {
            readClassificationFile << it->first << ":" << it->second << " ";
        }
        readClassificationFile << "\n";
    }
}

void Classifier::writeReportFile(const string &reportFileName, NcbiTaxonomy &ncbiTaxonomy, int numOfQuery) {
    unordered_map<TaxID, TaxonCounts> cladeCounts = ncbiTaxonomy.getCladeCounts(taxCounts);
    FILE *fp;
    fp = fopen(reportFileName.c_str(), "w");
    writeReport(fp, ncbiTaxonomy, cladeCounts, numOfQuery);
    fclose(fp);
}

void Classifier::writeReport(FILE *fp, const NcbiTaxonomy &ncbiTaxonomy,
                             const unordered_map<TaxID, TaxonCounts> &cladeCounts, unsigned long totalReads,
                             TaxID taxID, int depth) {
    auto it = cladeCounts.find(taxID);
    unsigned int cladeCount = (it == cladeCounts.end() ? 0 : it->second.cladeCount);
    unsigned int taxCount = (it == cladeCounts.end() ? 0 : it->second.taxCount);
    if (taxID == 0) {
        if (cladeCount > 0) {
            fprintf(fp, "%.2f\t%i\t%i\t0\tno rank\tunclassified\n", 100 * cladeCount / double(totalReads), cladeCount,
                    taxCount);
        }
        writeReport(fp, ncbiTaxonomy, cladeCounts, totalReads, 1);
    } else {
        if (cladeCount == 0) {
            return;
        }
        const TaxonNode *taxon = ncbiTaxonomy.taxonNode(taxID);
        fprintf(fp, "%.2f\t%i\t%i\t%i\t%s\t%s%s\n", 100 * cladeCount / double(totalReads), cladeCount, taxCount, taxID,
                taxon->rank.c_str(), string(2 * depth, ' ').c_str(), taxon->name.c_str());
        vector<TaxID> children = it->second.children;
        sort(children.begin(), children.end(),
             [&](int a, int b) { return cladeCountVal(cladeCounts, a) > cladeCountVal(cladeCounts, b); });
        for (TaxID childTaxId: children) {
            if (cladeCounts.count(childTaxId)) {
                writeReport(fp, ncbiTaxonomy, cladeCounts, totalReads, childTaxId, depth + 1);
            } else {
                break;
            }
        }
    }
}

unsigned int Classifier::cladeCountVal(const std::unordered_map<TaxID, TaxonCounts> &map, TaxID key) {
    typename std::unordered_map<TaxID, TaxonCounts>::const_iterator it = map.find(key);
    if (it == map.end()) {
        return 0;
    } else {
        return it->second.cladeCount;
    }
}

