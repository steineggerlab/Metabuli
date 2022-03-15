//
// Created by KJB on 01/09/2020.
//

#include "Classifier.h"
#include "LocalParameters.h"
#include <ctime>

Classifier::Classifier() {
    //seqIterator = new SeqIterator();
    numOfSplit = 0;
    selectedMatchCount = 0;
    queryCount = 0;
    perfectMatchCount = 0;

    correctCnt = 0;
    perfectCnt = 0;
    classifiedCnt = 0;
    speciesCnt = 0;
    subspCnt = 0;
    genusCnt = 0;
    phylumCnt = 0;
    classCnt = 0;
    orderCnt = 0;
    familyCnt = 0;
    superCnt = 0;
}

Classifier::~Classifier() {}//delete seqIterator; }

void Classifier::startClassify(const char *queryFileName,
                               const char *targetDiffIdxFileName,
                               const char *targetInfoFileName,
                               const char *diffIdxSplitFileName,
                               vector<int> &taxIdList,
                               const LocalParameters &par,
                               NcbiTaxonomy &taxonomy) {
//    ///-----------------------------------------------------
//    unordered_map<int,int> genus;
//    for(size_t i = 0 ; i < ncbiTaxonomy.taxonNodes.size(); i++){
//        if(ncbiTaxonomy.taxonNodes[i].rank == "species"){
//            genus[ncbiTaxonomy.taxonNodes[i].parentTaxId]++;
//        }
//    }
//    auto it = genus.begin();
//    size_t cnt = 0;
//    int archeaCnt = 0;
//    int bacteriaCnt = 0;
//    for(auto it = genus.begin(); it != genus.end(); it++){
//        if(it->second == 1){
//            if(it->first < 5839)
//                archeaCnt ++;
//            else
//                bacteriaCnt ++;
//        }
//        cnt ++;
//    }
//    cout<<archeaCnt<<" "<<bacteriaCnt<<" "<<cnt<<endl;
//    return;
    ///-------------------------------------------------------
    vector<int> speciesTaxIdList;
    vector<TaxID> genusTaxIdList;
    taxonomy.createTaxIdListAtRank(taxIdList, speciesTaxIdList, "species");
    taxonomy.createTaxIdListAtRank(taxIdList, genusTaxIdList, "genus");

    //output file
    char matchFileName[300];
    sprintf(matchFileName, "%s_match2", queryFileName);
    FILE *matchFile = fopen(matchFileName, "wb");

    // Allocate memory for buffers
    QueryKmerBuffer kmerBuffer(kmerBufSize);
    Buffer<Match> matchBuffer(kmerBufSize);

    // Load query file
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
//        queryList = new Query[numOfSeq + numOfSeq2];
        if (numOfSeq > numOfSeq2) {
            queryList = new Query[numOfSeq];
        } else {
            numOfSeq = numOfSeq2;
            queryList = new Query[numOfSeq];
        }
    }

    // Checker for multi-threading
    bool *processedSeqChecker = new bool[numOfSeq];
    fill_n(processedSeqChecker, numOfSeq, false);
    size_t processedSeqCnt = 0;

    // Timer
    size_t numOfTatalQueryKmerCnt = 0;

    // Extract k-mers from query sequences and compare them to target k-mer DB
    omp_set_num_threads(par.threads);
    while (processedSeqCnt < numOfSeq) {
        time_t beforeKmerExtraction = time(nullptr);
        if (par.seqMode == 1 || par.seqMode == 3) { // Single-end short-read sequence or long-read sequence
            fillQueryKmerBufferParallel(kmerBuffer,
                                        queryFile,
                                        sequences,
                                        processedSeqChecker,
                                        processedSeqCnt,
                                        queryList,
                                        par);
        } else if (par.seqMode == 2) {
            fillQueryKmerBufferParallel_paired(kmerBuffer,
                                               queryFile,
                                               queryFile2,
                                               sequences,
                                               sequences2,
                                               processedSeqChecker,
                                               processedSeqCnt,
                                               queryList,
                                               numOfSeq,
                                               par);
        }
        numOfTatalQueryKmerCnt += kmerBuffer.startIndexOfReserve;
        cout << "Time spent for k-mer extraction: " << double(time(nullptr) - beforeKmerExtraction) << endl;

        time_t beforeQueryKmerSort = time(nullptr);
        SORT_PARALLEL(kmerBuffer.buffer, kmerBuffer.buffer + kmerBuffer.startIndexOfReserve,
                      Classifier::compareForLinearSearch);
        cout << "Time spent for sorting query k-mer list: " << double(time(nullptr) - beforeQueryKmerSort) << endl;

        //time_t beforeSearch = time(nullptr);
        linearSearchParallel(kmerBuffer.buffer, kmerBuffer.startIndexOfReserve, targetDiffIdxFileName,
                             targetInfoFileName,
                             diffIdxSplitFileName, matchBuffer, taxIdList, speciesTaxIdList, genusTaxIdList, matchFile,
                             par);
        // cout<<"Time spent for linearSearch: " << double(time(nullptr) - beforeSearch) << endl;
        cout << "The number of matches: " << matchBuffer.startIndexOfReserve << endl;
    }
    cout << "Number of query k-mers: " << numOfTatalQueryKmerCnt << endl;

    if (par.memoryMode == 1) {
        writeMatches(matchBuffer, matchFile);
        fclose(matchFile);
        free(matchBuffer.buffer);
        struct MmapedData<Match> matchList = mmapData<Match>(matchFileName);
        size_t numOfMatches = matchList.fileSize / sizeof(Match);
        time_t beforeSortMatches = time(nullptr);
        SORT_PARALLEL(matchList.data, matchList.data + numOfMatches, Classifier::sortByGenusAndSpecies2);
        cout << "Time spent for sorting matches: " << double(time(nullptr) - beforeSortMatches) << endl;
        time_t beforeAnalyze = time(nullptr);
        analyseResultParallel(taxonomy, matchList.data, numOfMatches, (int) numOfSeq, queryList, par);
        cout << "Time spent for analyzing: " << double(time(nullptr) - beforeAnalyze) << endl;
    } else {
//        fclose(matchFile);
        time_t beforeSortMatches = time(nullptr);
        SORT_PARALLEL(matchBuffer.buffer, matchBuffer.buffer + matchBuffer.startIndexOfReserve,
                      Classifier::sortByGenusAndSpecies2);
        cout << "Time spent for sorting matches: " << double(time(nullptr) - beforeSortMatches) << endl;
        time_t beforeAnalyze = time(nullptr);
        analyseResultParallel(taxonomy, matchBuffer.buffer, matchBuffer.startIndexOfReserve, (int) numOfSeq, queryList,
                              par);
        cout << "Time spent for analyzing: " << double(time(nullptr) - beforeAnalyze) << endl;
        free(matchBuffer.buffer);
    }


    // Write report files
    ofstream readClassificationFile;
    readClassificationFile.open(par.filenames[3] + "/" + par.filenames[4] + "_ReadClassification.tsv");
    writeReadClassification(queryList, (int) numOfSeq, readClassificationFile);
    readClassificationFile.close();
    writeReportFile(par.filenames[3] + "/" + par.filenames[4] + "_CompositionReport.tsv", taxonomy, numOfSeq);

    // Below is for developing
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

    free(kmerBuffer.buffer);

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
#pragma omp parallel default(none), shared(checker, hasOverflow, processedSeqCnt, kmerBuffer, seqFile, seqs, cout, queryList)
    {
        SeqIterator seqIterator;
        size_t posToWrite;
#pragma omp for schedule(dynamic, 1)
        for (size_t i = 0; i < seqs.size(); i++) {
            if (checker[i] == false && !hasOverflow) {
                kseq_buffer_t buffer(const_cast<char *>(&seqFile.data[seqs[i].start]), seqs[i].length);
                kseq_t *seq = kseq_init(&buffer);
                kseq_read(seq);
                seqIterator.sixFrameTranslation(seq->seq.s);
                size_t kmerCnt = getQueryKmerNumber((int) strlen(seq->seq.s));
                posToWrite = kmerBuffer.reserveMemory(kmerCnt);
                if (posToWrite + kmerCnt < kmerBuffer.bufferSize) {
                    seqIterator.fillQueryKmerBuffer(seq->seq.s, kmerBuffer, posToWrite, i);
                    checker[i] = true;
                    queryList[i].queryLength = getMaxCoveredLength((int) strlen(seq->seq.s));
                    queryList[i].queryId = i;
                    queryList[i].name = string(seq->name.s);
                    queryList[i].kmerCnt = (int) kmerCnt;
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
    return (getMaxCoveredLength(queryLength) / 3 - kmerLength + 1) * 6;
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

#pragma omp parallel default(none), shared(checker, hasOverflow, processedSeqCnt, kmerBuffer, seqFile1, seqFile2, seqs, seqs2, cout, queryList, numOfSeq)
    {
        SeqIterator seqIterator;
        SeqIterator seqIterator2;
        size_t posToWrite;
#pragma omp for schedule(dynamic, 1)
        for (size_t i = 0; i < numOfSeq; i++) {
            if (checker[i] == false && !hasOverflow) {
                // Read 1
                kseq_buffer_t buffer(const_cast<char *>(&seqFile1.data[seqs[i].start]), seqs[i].length);
                kseq_t *seq = kseq_init(&buffer);
                kseq_read(seq);
                size_t kmerCnt = getQueryKmerNumber((int) strlen(seq->seq.s));

                // Read 2
                kseq_buffer_t buffer2(const_cast<char *>(&seqFile2.data[seqs2[i].start]), seqs2[i].length);
                kseq_t *seq2 = kseq_init(&buffer2);
                kseq_read(seq2);
                kmerCnt += getQueryKmerNumber((int) strlen(seq2->seq.s));

                posToWrite = kmerBuffer.reserveMemory(kmerCnt);
                if (posToWrite + kmerCnt < kmerBuffer.bufferSize) {
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
                    // queryList[i].queryLength = getMaxCoveredLength((int) strlen(seq->seq.s));
//                                                + getMaxCoveredLength((int) strlen(seq2->seq.s));
                    queryList[i].queryLength2 = getMaxCoveredLength((int) strlen(seq2->seq.s));
                    queryList[i].queryId = (int) i;
                    queryList[i].name = string(seq->name.s);
                    queryList[i].kmerCnt = (int) kmerCnt;
#pragma omp atomic
                    processedSeqCnt++;
                } else {
#pragma omp atomic
                    kmerBuffer.startIndexOfReserve -= kmerCnt;
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
                                      FILE *matchFile, const LocalParameters &par) {

    struct MmapedData<uint16_t> targetDiffIdxList = mmapData<uint16_t>(targetDiffIdxFileName, 2);
    struct MmapedData<TargetKmerInfo> targetInfoList = mmapData<TargetKmerInfo>(targetInfoFileName, 2);
    struct MmapedData<DiffIdxSplit> diffIdxSplits = mmapData<DiffIdxSplit>(diffIdxSplitsFileName);

    cout << "linearSearch start..." << endl;
    // Find the first index of garbage query k-mer (UINT64_MAX) and discard from there
    for (size_t checkN = queryKmerCnt - 1; checkN > 0; checkN--) {
        if (queryKmerList[checkN].ADkmer != UINT64_MAX) {
            queryKmerCnt = checkN + 1;
            break;
        }
    }

    // Filter out meaningless target querySplits
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
            queryAA = AminoAcid(queryKmerList[splitWidth].ADkmer);
            if (queryAA <= AminoAcid(diffIdxSplits.data[tSplitCnt].ADkmer)) {
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
            queryAA = AminoAcid(queryKmerList[splitWidth * i].ADkmer);
            bool needLastTargetBlock = true;
            for (size_t j = 0; j < numOfDiffIdxSplits_use; j++) {
                if (queryAA <= AminoAcid(diffIdxSplits.data[j].ADkmer)) {
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
                cout << "needLastTargetBlock" << endl;
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

    munmap(targetDiffIdxList.data, targetDiffIdxList.fileSize + 1);
    munmap(targetInfoList.data, targetInfoList.fileSize + 1);
    munmap(diffIdxSplits.data, diffIdxSplits.fileSize + 1);

    cout << "The number of target k-mers: " << numOfTargetKmer << endl;

    struct MmapedData<uint16_t> targetDiffIdxList2 = mmapData<uint16_t>(targetDiffIdxFileName, 2);
    struct MmapedData<TargetKmerInfo> targetInfoList2 = mmapData<TargetKmerInfo>(targetInfoFileName, 2);

    time_t beforeSearch = time(nullptr);

    while (completedSplitCnt < threadNum) {
        bool hasOverflow = false;
#pragma omp parallel default(none), shared(numOfDiffIdx, completedSplitCnt, splitCheckList, numOfTargetKmer, hasOverflow, \
querySplits, queryKmerList, targetDiffIdxList2, targetInfoList2, matchBuffer, cout, genusTaxIdList, taxIdList, spTaxIdList, matchFile)
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
            int matchBufferSize = 1000000; // 32 Mb
            Match *matches = new Match[1000000];
            int matchCnt = 0;

            //vectors for selected target k-mers
            vector<uint8_t> selectedHammingSum;
            vector<size_t> selectedMatches;
            vector<uint16_t> selectedHammings;
            size_t startIdxOfAAmatch = 0;
            size_t posToWrite;
            int currMatchNum;
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
                    currentTargetKmer = getNextTargetKmer(currentTargetKmer, targetDiffIdxList2.data, diffIdxPos);
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
                        if (matchCnt + currMatchNum > matchBufferSize) {
                            // Check if the shared buffer is full.
                            posToWrite = matchBuffer.reserveMemory(matchCnt);
                            if (posToWrite + matchCnt >=
                                matchBuffer.bufferSize) { // full -> write matches to file first
                                hasOverflow = true;
                                querySplits[i].start = lastMovedQueryIdx + 1; // TODO
#pragma omp atomic
                                matchBuffer.startIndexOfReserve -= matchCnt;
                                break;
                            } else { // not full -> copy matches to the shared buffer
                                moveMatches(matchBuffer.buffer + posToWrite, matches, matchCnt);
                                lastMovedQueryIdx = j;
                            }
                        }

                        for (int k = 0; k < currMatchNum; k++) {
                            if (targetInfoList2.data[selectedMatches[k]].redundancy) {
                                matches[matchCnt] = {queryKmerList[j].info.sequenceID,
                                                     spTaxIdList[targetInfoList2.data[selectedMatches[k]].sequenceID],
                                                     spTaxIdList[targetInfoList2.data[selectedMatches[k]].sequenceID],
                                                     genusTaxIdList[targetInfoList2.data[selectedMatches[k]].sequenceID],
                                                     queryKmerList[j].info.pos,
                                                     queryKmerList[j].info.frame,
                                                     selectedHammingSum[k], 1, selectedHammings[k]};
                            } else {
                                matches[matchCnt] = {queryKmerList[j].info.sequenceID,
                                                     taxIdList[targetInfoList2.data[selectedMatches[k]].sequenceID],
                                                     spTaxIdList[targetInfoList2.data[selectedMatches[k]].sequenceID],
                                                     genusTaxIdList[targetInfoList2.data[selectedMatches[k]].sequenceID],
                                                     queryKmerList[j].info.pos,
                                                     queryKmerList[j].info.frame,
                                                     selectedHammingSum[k], 0, selectedHammings[k]};
                            }
                            matchCnt++;
                        }
                        continue;
                    }
                    selectedMatches.clear();
                    selectedHammingSum.clear();
                    selectedHammings.clear();

                    // Reuse the candidate target k-mers to compare in DNA level if queries are the same at amino acid level but not at DNA level
                    if (currentQueryAA == AminoAcid(queryKmerList[j].ADkmer)) {
                        compareDna(queryKmerList[j].ADkmer, candidateTargetKmers, startIdxOfAAmatch, selectedMatches,
                                   selectedHammingSum, selectedHammings);
                        currMatchNum = selectedMatches.size();
                        // If local buffer is full, copy them to the shared buffer.
                        if (matchCnt + currMatchNum > matchBufferSize) {
                            // Check if the shared buffer is full.
                            posToWrite = matchBuffer.reserveMemory(matchCnt);
                            if (posToWrite + matchCnt >=
                                matchBuffer.bufferSize) { // full -> write matches to file first
                                hasOverflow = true;
                                querySplits[i].start = lastMovedQueryIdx + 1;
#pragma omp atomic
                                matchBuffer.startIndexOfReserve -= matchCnt;
                                break;
                            } else { // not full -> copy matches to the shared buffer
                                moveMatches(matchBuffer.buffer + posToWrite, matches, matchCnt);
                                lastMovedQueryIdx = j;
                            }
                        }

                        for (int k = 0; k < currMatchNum; k++) {
                            if (targetInfoList2.data[selectedMatches[k]].redundancy) {
                                matches[matchCnt] = {queryKmerList[j].info.sequenceID,
                                                     spTaxIdList[targetInfoList2.data[selectedMatches[k]].sequenceID],
                                                     spTaxIdList[targetInfoList2.data[selectedMatches[k]].sequenceID],
                                                     genusTaxIdList[targetInfoList2.data[selectedMatches[k]].sequenceID],
                                                     queryKmerList[j].info.pos,
                                                     queryKmerList[j].info.frame,
                                                     selectedHammingSum[k], 1, selectedHammings[k]};
                            } else {
                                matches[matchCnt] = {queryKmerList[j].info.sequenceID,
                                                     taxIdList[targetInfoList2.data[selectedMatches[k]].sequenceID],
                                                     spTaxIdList[targetInfoList2.data[selectedMatches[k]].sequenceID],
                                                     genusTaxIdList[targetInfoList2.data[selectedMatches[k]].sequenceID],
                                                     queryKmerList[j].info.pos,
                                                     queryKmerList[j].info.frame,
                                                     selectedHammingSum[k], 0, selectedHammings[k]};
                            }
                            matchCnt++;
                        }
                        continue;
                    }
                    candidateTargetKmers.clear();

                    // Get next query, and start to find
                    currentQuery = queryKmerList[j].ADkmer;
                    currentQueryAA = AminoAcid(currentQuery);

                    // Skip target k-mers that are not matched in amino acid level
                    while (AminoAcid(currentQuery) > AminoAcid(currentTargetKmer) &&
                           (targetInfoIdx < numOfTargetKmer) && (diffIdxPos != numOfDiffIdx)) {
                        currentTargetKmer = getNextTargetKmer(currentTargetKmer, targetDiffIdxList2.data, diffIdxPos);
                        targetInfoIdx++;
                    }

                    if (AminoAcid(currentQuery) !=
                        AminoAcid(currentTargetKmer)) // Move to next query k-mer if there isn't any match.
                        continue;
                    else
                        startIdxOfAAmatch = targetInfoIdx;

                    // Load target k-mers that are matched in amino acid level
                    while (AminoAcid(currentQuery) == AminoAcid(currentTargetKmer) &&
                           (targetInfoIdx < numOfTargetKmer) && (diffIdxPos != numOfDiffIdx)) {
                        candidateTargetKmers.push_back(currentTargetKmer);
                        currentTargetKmer = getNextTargetKmer(currentTargetKmer, targetDiffIdxList2.data, diffIdxPos);
                        targetInfoIdx++;
                    }

                    // Compare the current query and the loaded target k-mers and select
                    compareDna(currentQuery, candidateTargetKmers, startIdxOfAAmatch, selectedMatches,
                               selectedHammingSum, selectedHammings);

                    // If local buffer is full, copy them to the shared buffer.
                    currMatchNum = selectedMatches.size();
                    if (matchCnt + currMatchNum > matchBufferSize) {
                        // Check if the shared buffer is full.
                        posToWrite = matchBuffer.reserveMemory(matchCnt);
                        if (posToWrite + matchCnt >= matchBuffer.bufferSize) { // full -> write matches to file first
                            hasOverflow = true;
                            querySplits[i].start = lastMovedQueryIdx + 1;
#pragma omp atomic
                            matchBuffer.startIndexOfReserve -= matchCnt;
                            break;
                        } else { // not full -> copy matches to the shared buffer
                            moveMatches(matchBuffer.buffer + posToWrite, matches, matchCnt);
                            lastMovedQueryIdx = j;
                        }
                    }

                    for (int k = 0; k < currMatchNum; k++) {
                        if (targetInfoList2.data[selectedMatches[k]].redundancy) {
                            matches[matchCnt] = {queryKmerList[j].info.sequenceID,
                                                 spTaxIdList[targetInfoList2.data[selectedMatches[k]].sequenceID],
                                                 spTaxIdList[targetInfoList2.data[selectedMatches[k]].sequenceID],
                                                 genusTaxIdList[targetInfoList2.data[selectedMatches[k]].sequenceID],
                                                 queryKmerList[j].info.pos,
                                                 queryKmerList[j].info.frame,
                                                 selectedHammingSum[k], 1, selectedHammings[k]};
                        } else {
                            matches[matchCnt] = {queryKmerList[j].info.sequenceID,
                                                 taxIdList[targetInfoList2.data[selectedMatches[k]].sequenceID],
                                                 spTaxIdList[targetInfoList2.data[selectedMatches[k]].sequenceID],
                                                 genusTaxIdList[targetInfoList2.data[selectedMatches[k]].sequenceID],
                                                 queryKmerList[j].info.pos,
                                                 queryKmerList[j].info.frame,
                                                 selectedHammingSum[k], 0, selectedHammings[k]};
                        }
                        matchCnt++;
                    }
                }

                // Move matches in the local buffer to the shared buffer
                posToWrite = matchBuffer.reserveMemory(matchCnt);
                if (posToWrite + matchCnt >= matchBuffer.bufferSize) {
                    querySplits[i].start = lastMovedQueryIdx + 1;
#pragma omp atomic
                    matchBuffer.startIndexOfReserve -= matchCnt;
                } else {
                    moveMatches(matchBuffer.buffer, matches, matchCnt);
                }
                delete[] matches;

                // Check whether current split is completed or not
                if (querySplits[i].start - 1 == querySplits[i].end) {
                    splitCheckList[i] = true;
#pragma omp atomic
                    completedSplitCnt++; //sync~~
                }
            }
        } // end of omp parallel
        writeMatches(matchBuffer, matchFile);
    } // end of while(completeSplitCnt < threadNum)
    cout << "Time spent for linearSearch: " << double(time(nullptr) - beforeSearch) << endl;

    munmap(targetDiffIdxList2.data, targetDiffIdxList2.fileSize + 1);
    munmap(targetInfoList2.data, targetInfoList2.fileSize + 1);
    free(splitCheckList);
    queryKmerCnt = 0;
    cout << "end of linear seach parallel" << endl;
}

void Classifier::moveMatches(Match *dest, Match *src, int &matchNum) {
    memcpy(dest, src, sizeof(Match) * matchNum);
    matchNum = 0;
}

void
Classifier::linearSearchParallel2(QueryKmer *queryKmerList, size_t &queryKmerCnt, const char *targetDiffIdxFileName,
                                  const char *targetInfoFileName, const char *diffIdxSplitsFileName,
                                  Buffer<Match> &matchBuffer, const vector<int> &taxIdList,
                                  const vector<int> &spTaxIdList, const vector<TaxID> &genusTaxIdList,
                                  FILE *matchFile, const LocalParameters &par) {

    struct MmapedData<uint16_t> targetDiffIdxList = mmapData<uint16_t>(targetDiffIdxFileName, 2);
    struct MmapedData<TargetKmerInfo> targetInfoList = mmapData<TargetKmerInfo>(targetInfoFileName, 2);
    struct MmapedData<DiffIdxSplit> diffIdxSplits = mmapData<DiffIdxSplit>(diffIdxSplitsFileName);

    cout << "linearSearch start..." << endl;
    // Find the first index of garbage query k-mer (UINT64_MAX) and discard from there
    for (size_t checkN = queryKmerCnt - 1; checkN > 0; checkN--) {
        if (queryKmerList[checkN].ADkmer != UINT64_MAX) {
            queryKmerCnt = checkN + 1;
            break;
        }
    }

    // Filter out meaningless target querySplits
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
    vector<QueryKmer *> querySplits2;
    int threadNum = par.threads;
    uint64_t queryAA;
    if (threadNum == 1) { //Single thread
        querySplits.emplace_back(0, queryKmerCnt - 1, queryKmerCnt, diffIdxSplits.data[0]);
    } else if (threadNum == 2) { //Two threads
        size_t splitWidth = queryKmerCnt / 2;
        querySplits.emplace_back(0, splitWidth - 1, splitWidth, diffIdxSplits.data[0]);
        for (size_t tSplitCnt = 0; tSplitCnt < numOfDiffIdxSplits_use; tSplitCnt++) {
            queryAA = AminoAcid(queryKmerList[splitWidth].ADkmer);
            if (queryAA <= AminoAcid(diffIdxSplits.data[tSplitCnt].ADkmer)) {
                tSplitCnt = tSplitCnt - (tSplitCnt != 0);
                querySplits.emplace_back(splitWidth, queryKmerCnt - 1, queryKmerCnt - splitWidth,
                                         diffIdxSplits.data[tSplitCnt]);
                break;
            }
        }
    } else { //More than two threads
        size_t splitWidth = queryKmerCnt / (threadNum - 1);
        querySplits.emplace_back(0, splitWidth - 1, splitWidth, diffIdxSplits.data[0]);
        querySplits2.push_back(new QueryKmer[splitWidth]);
        memcpy(querySplits2[0], queryKmerList, sizeof(QueryKmer) * splitWidth);

        for (int i = 1; i < threadNum; i++) {
            queryAA = AminoAcid(queryKmerList[splitWidth * i].ADkmer);
            bool needLastTargetBlock = true;
            for (size_t j = 0; j < numOfDiffIdxSplits_use; j++) {
                if (queryAA <= AminoAcid(diffIdxSplits.data[j].ADkmer)) {
                    j = j - (j != 0);
                    if (i != threadNum - 1) {
                        querySplits.emplace_back(splitWidth * i, splitWidth * (i + 1) - 1, splitWidth,
                                                 diffIdxSplits.data[j]);
                        querySplits2.push_back(new QueryKmer[splitWidth]);
                        memcpy(querySplits2[i], queryKmerList + splitWidth * i, sizeof(QueryKmer) * splitWidth);
                    } else {
                        querySplits.emplace_back(splitWidth * i, queryKmerCnt - 1, queryKmerCnt - splitWidth * i,
                                                 diffIdxSplits.data[j]);
                        querySplits2.push_back(new QueryKmer[splitWidth]);
                        memcpy(querySplits2[i], queryKmerList + splitWidth * i,
                               sizeof(QueryKmer) * (queryKmerCnt - splitWidth * i));
                    }
                    needLastTargetBlock = false;
                    break;
                }
            }
            if (needLastTargetBlock) {
                cout << "needLastTargetBlock" << endl;
                if (i != threadNum - 1) {
                    querySplits.emplace_back(splitWidth * i, splitWidth * (i + 1) - 1, splitWidth,
                                             diffIdxSplits.data[numOfDiffIdxSplits_use - 1]);
                    querySplits2.push_back(new QueryKmer[splitWidth]);
                    memcpy(querySplits2[i], queryKmerList + splitWidth * i, sizeof(QueryKmer) * splitWidth);
                } else {
                    querySplits.emplace_back(splitWidth * i, queryKmerCnt - 1, queryKmerCnt - splitWidth * i,
                                             diffIdxSplits.data[numOfDiffIdxSplits_use - 1]);
                    querySplits2.push_back(new QueryKmer[splitWidth]);
                    memcpy(querySplits2[i], queryKmerList + splitWidth * i,
                           sizeof(QueryKmer) * (queryKmerCnt - splitWidth * i));
                }
            }
        }
    }

    bool *splitCheckList = (bool *) malloc(sizeof(bool) * threadNum);
    fill_n(splitCheckList, threadNum, false);
    int completedSplitCnt = 0;
    size_t numOfTargetKmer = targetInfoList.fileSize / sizeof(TargetKmerInfo);
    size_t numOfDiffIdx = targetDiffIdxList.fileSize / sizeof(uint16_t);

    munmap(targetDiffIdxList.data, targetDiffIdxList.fileSize + 1);
    munmap(targetInfoList.data, targetInfoList.fileSize + 1);
    munmap(diffIdxSplits.data, diffIdxSplits.fileSize + 1);

    cout << "The number of target k-mers: " << numOfTargetKmer << endl;

    struct MmapedData<uint16_t> targetDiffIdxList2 = mmapData<uint16_t>(targetDiffIdxFileName, 2);
    struct MmapedData<TargetKmerInfo> targetInfoList2 = mmapData<TargetKmerInfo>(targetInfoFileName, 2);

    time_t beforeSearch = time(nullptr);

    while (completedSplitCnt < threadNum) {
        bool hasOverflow = false;
#pragma omp parallel default(none), shared(numOfDiffIdx, completedSplitCnt, splitCheckList, numOfTargetKmer, hasOverflow, \
querySplits, queryKmerList, targetDiffIdxList2, targetInfoList2, matchBuffer, cout, genusTaxIdList, taxIdList, spTaxIdList)
        {
            //query variables
            uint64_t currentQuery = UINT64_MAX;
            uint64_t currentQueryAA = UINT64_MAX;

            //target variables
            size_t diffIdxPos = 0;
            size_t targetInfoIdx = 0;
            vector<uint64_t> candidateTargetKmers; //vector for candidate target k-mer, some of which are selected after based on hamming distance
            uint64_t currentTargetKmer;

            //vectors for selected target k-mers
            vector<uint8_t> selectedHammingSum;
            vector<size_t> selectedMatches;
            vector<uint16_t> selectedHammings;
            size_t startIdxOfAAmatch = 0;
            size_t posToWrite;
            size_t range;
            QueryKmer *currQuerySplit;
#pragma omp for schedule(dynamic, 1)
            for (size_t i = 0; i < querySplits.size(); i++) {
                currQuerySplit = new QueryKmer[querySplits[i].length];
                memcpy(currQuerySplit, queryKmerList + querySplits[i].start, sizeof(QueryKmer) * querySplits[i].length);
                if (hasOverflow || splitCheckList[i]) {
                    continue;
                }
                targetInfoIdx = querySplits[i].diffIdxSplit.infoIdxOffset - (i != 0);
                diffIdxPos = querySplits[i].diffIdxSplit.diffIdxOffset;
                currentTargetKmer = querySplits[i].diffIdxSplit.ADkmer;
                if (i == 0) {
                    currentTargetKmer = getNextTargetKmer(currentTargetKmer, targetDiffIdxList2.data, diffIdxPos);
                }
                currentQuery = UINT64_MAX;
                currentQueryAA = UINT64_MAX;

                for (size_t j = 0; j < querySplits[i].length; j++) {
                    querySplits[i].start++;
                    // Reuse the comparison data if queries are exactly identical
                    if (currentQuery == currQuerySplit[j].ADkmer) {
                        posToWrite = matchBuffer.reserveMemory(selectedMatches.size());
                        if (posToWrite + selectedMatches.size() >= matchBuffer.bufferSize) {
                            hasOverflow = true;
                            querySplits[i].start = j;
#pragma omp atomic
                            matchBuffer.startIndexOfReserve -= selectedMatches.size();
                            break;
                        } else {
                            range = selectedMatches.size();
                            for (size_t k = 0; k < range; k++) {
                                if (targetInfoList2.data[selectedMatches[k]].redundancy) {
                                    matchBuffer.buffer[posToWrite] = {currQuerySplit[j].info.sequenceID,
                                                                      spTaxIdList[targetInfoList2.data[selectedMatches[k]].sequenceID],
                                                                      spTaxIdList[targetInfoList2.data[selectedMatches[k]].sequenceID],
                                                                      genusTaxIdList[targetInfoList2.data[selectedMatches[k]].sequenceID],
                                                                      currQuerySplit[j].info.pos,
                                                                      currQuerySplit[j].info.frame,
                                                                      selectedHammingSum[k], 1, selectedHammings[k]};
                                } else {
                                    matchBuffer.buffer[posToWrite] = {currQuerySplit[j].info.sequenceID,
                                                                      taxIdList[targetInfoList2.data[selectedMatches[k]].sequenceID],
                                                                      spTaxIdList[targetInfoList2.data[selectedMatches[k]].sequenceID],
                                                                      genusTaxIdList[targetInfoList2.data[selectedMatches[k]].sequenceID],
                                                                      currQuerySplit[j].info.pos,
                                                                      currQuerySplit[j].info.frame,
                                                                      selectedHammingSum[k], 0, selectedHammings[k]};
                                }
                                posToWrite++;
                            }
                        }
                        continue;
                    }
                    selectedMatches.clear();
                    selectedHammingSum.clear();
                    selectedHammings.clear();

                    ///Reuse the candidate target k-mers to compare in DNA level if queries are the same at amino acid level but not at DNA level
                    if (currentQueryAA == AminoAcid(currQuerySplit[j].ADkmer)) {
                        compareDna(currQuerySplit[j].ADkmer, candidateTargetKmers, startIdxOfAAmatch, selectedMatches,
                                   selectedHammingSum, selectedHammings);
                        posToWrite = matchBuffer.reserveMemory(selectedMatches.size());
                        if (posToWrite + selectedMatches.size() >= matchBuffer.bufferSize) {
                            hasOverflow = true;
                            querySplits[i].start = j;
#pragma omp atomic
                            matchBuffer.startIndexOfReserve -= selectedMatches.size();
                            break;
                        } else {
                            range = selectedMatches.size();
                            for (size_t k = 0; k < range; k++) {
                                if (targetInfoList2.data[selectedMatches[k]].redundancy) {
                                    matchBuffer.buffer[posToWrite] = {currQuerySplit[j].info.sequenceID,
                                                                      spTaxIdList[targetInfoList2.data[selectedMatches[k]].sequenceID],
                                                                      spTaxIdList[targetInfoList2.data[selectedMatches[k]].sequenceID],
                                                                      genusTaxIdList[targetInfoList2.data[selectedMatches[k]].sequenceID],
                                                                      currQuerySplit[j].info.pos,
                                                                      currQuerySplit[j].info.frame,
                                                                      selectedHammingSum[k], 1, selectedHammings[k]};
                                } else {
                                    matchBuffer.buffer[posToWrite] = {currQuerySplit[j].info.sequenceID,
                                                                      taxIdList[targetInfoList2.data[selectedMatches[k]].sequenceID],
                                                                      spTaxIdList[targetInfoList2.data[selectedMatches[k]].sequenceID],
                                                                      genusTaxIdList[targetInfoList2.data[selectedMatches[k]].sequenceID],
                                                                      currQuerySplit[j].info.pos,
                                                                      currQuerySplit[j].info.frame,
                                                                      selectedHammingSum[k], 0, selectedHammings[k]};
                                }
                                posToWrite++;
                            }
                        }
                        continue;
                    }
                    candidateTargetKmers.clear();

                    // Get next query, and start to find
                    currentQuery = currQuerySplit[j].ADkmer;
                    currentQueryAA = AminoAcid(currentQuery);

                    // Skip target k-mers that are not matched in amino acid level
                    while (AminoAcid(currentQuery) > AminoAcid(currentTargetKmer) &&
                           (targetInfoIdx < numOfTargetKmer) && (diffIdxPos != numOfDiffIdx)) {
                        currentTargetKmer = getNextTargetKmer(currentTargetKmer, targetDiffIdxList2.data, diffIdxPos);
                        targetInfoIdx++;
                    }

                    if (AminoAcid(currentQuery) !=
                        AminoAcid(currentTargetKmer)) // Move to next query k-mer if there isn't any match.
                        continue;
                    else
                        startIdxOfAAmatch = targetInfoIdx;

                    // Load target k-mers that are matched in amino acid level
                    while (AminoAcid(currentQuery) == AminoAcid(currentTargetKmer) &&
                           (targetInfoIdx < numOfTargetKmer) && (diffIdxPos != numOfDiffIdx)) {
                        candidateTargetKmers.push_back(currentTargetKmer);
                        currentTargetKmer = getNextTargetKmer(currentTargetKmer, targetDiffIdxList2.data, diffIdxPos);
                        targetInfoIdx++;
                    }

                    // Compare the current query and the loaded target k-mers and select
                    compareDna(currentQuery, candidateTargetKmers, startIdxOfAAmatch, selectedMatches,
                               selectedHammingSum, selectedHammings);
                    posToWrite = matchBuffer.reserveMemory(selectedMatches.size());
                    if (posToWrite + selectedMatches.size() >= matchBuffer.bufferSize) {
                        hasOverflow = true;
                        querySplits[i].start = j;
#pragma omp atomic
                        matchBuffer.startIndexOfReserve -= selectedMatches.size();
                        break;
                    } else {
                        range = selectedMatches.size();
                        for (size_t k = 0; k < range; k++) {
                            if (targetInfoList2.data[selectedMatches[k]].redundancy) {
                                matchBuffer.buffer[posToWrite] = {currQuerySplit[j].info.sequenceID,
                                                                  spTaxIdList[targetInfoList2.data[selectedMatches[k]].sequenceID],
                                                                  spTaxIdList[targetInfoList2.data[selectedMatches[k]].sequenceID],
                                                                  genusTaxIdList[targetInfoList2.data[selectedMatches[k]].sequenceID],
                                                                  currQuerySplit[j].info.pos,
                                                                  currQuerySplit[j].info.frame,
                                                                  selectedHammingSum[k], 1, selectedHammings[k]};
                            } else {
                                matchBuffer.buffer[posToWrite] = {currQuerySplit[j].info.sequenceID,
                                                                  taxIdList[targetInfoList2.data[selectedMatches[k]].sequenceID],
                                                                  spTaxIdList[targetInfoList2.data[selectedMatches[k]].sequenceID],
                                                                  genusTaxIdList[targetInfoList2.data[selectedMatches[k]].sequenceID],
                                                                  currQuerySplit[j].info.pos,
                                                                  currQuerySplit[j].info.frame,
                                                                  selectedHammingSum[k], 0, selectedHammings[k]};
                            }
                            posToWrite++;
                        }
                    }
                }

                // Check whether current split is completed or not
                if (querySplits[i].start - 1 == querySplits[i].end) {
                    splitCheckList[i] = true;
#pragma omp atomic
                    completedSplitCnt++;
                }
            }

        }
        if (hasOverflow)
            writeMatches(matchBuffer, matchFile);
    }
    munmap(targetDiffIdxList2.data, targetDiffIdxList2.fileSize + 1);
    munmap(targetInfoList2.data, targetInfoList2.fileSize + 1);
    cout << "Time spent for linearSearch: " << double(time(nullptr) - beforeSearch) << endl;
    free(splitCheckList);
    queryKmerCnt = 0;
    cout << "end of linear seach parallel" << endl;
}

void Classifier::writeMatches(Buffer<Match> &matchBuffer, FILE *matchFile) {
    fwrite(matchBuffer.buffer, sizeof(Match), matchBuffer.startIndexOfReserve, matchFile);
    matchBuffer.startIndexOfReserve = 0;
}

// It compares query k-mers to target k-mers.
// If a query has matches, the matches with the smallest hamming distance will be selected
void Classifier::compareDna(uint64_t &query, vector<uint64_t> &targetKmersToCompare, const size_t &startIdx,
                            vector<size_t> &selectedMatches, vector<uint8_t> &selectedHammingSum,
                            vector<uint16_t> &selectedHammings) {

    size_t size = targetKmersToCompare.size();
    auto *hammingSums = new uint8_t[size + 1];
    uint8_t currentHammingSum;
    uint8_t minHammingSum = UINT8_MAX;

    // Calculate hamming distance
    for (size_t i = 0; i < size; i++) {
        currentHammingSum = getHammingDistanceSum(query, targetKmersToCompare[i]);
        if (currentHammingSum < minHammingSum)
            minHammingSum = currentHammingSum;
        hammingSums[i] = currentHammingSum;
    }

    // Select target k-mers that passed hamming criteria
    for (size_t h = 0; h < size; h++) {
        if (hammingSums[h] == minHammingSum) {
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

    for (int i = 0; i < seqNum; i++) {
        ++taxCounts[queryList[i].classification];
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
                 int(matchList[i].frame) << " " << matchList[i].position << " " << int(matchList[i].hamming) << endl;
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
                 << " " <<
                 int(matchesForLCA[i].frame) << " " << matchesForLCA[i].position << " " << int(matchesForLCA[i].hamming)
                 << endl;
        }
    }

    if (PRINT) {
        sort(matchesForLCA.begin(), matchesForLCA.end(), Classifier::sortByGenusAndSpecies);
    }

    //If there is no proper genus for current query, it is un-classified.
    if (res == 3) {
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
                cout << i << " " << int(matchesForLCA[i].frame) << " " << matchesForLCA[i].position << " " <<
                     matchesForLCA[i].taxID << " " << int(matchesForLCA[i].hamming) << " " << matchesForLCA[i].red
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

    // Classify at the genus rank if more than one species are selected.
    // Classify at the genus rank if the score at species level is not enough.
    if (species.size() > 1
        || (speciesScrCov.score < 0.9 && !ncbiTaxonomy.IsAncestor(par.virusTaxId, matchesForLCA[0].taxID))) {
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
        minStrainSpecificCnt = 5;
        if (queryLength > 5000) {
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
            cout << i << " " << int(matchesForLCA[i].frame) << " " << matchesForLCA[i].position << " " <<
                 matchesForLCA[i].taxID << " " << int(matchesForLCA[i].hamming) << " " << matchesForLCA[i].red << endl;
        }
        cout << "Score: " << speciesScrCov.score << "  " << selectedSpecies << " "
             << ncbiTaxonomy.taxonNode(selectedSpecies)->rank
             << endl;
    }
}

//TaxID Classifier::chooseBestTaxon_index(NcbiTaxonomy &ncbiTaxonomy, uint32_t currentQuery,
//                                  size_t offset, size_t end, Match * matchList, Query * queryList,
//                                  const LocalParameters & par ) {
//    int queryLength = queryList[currentQuery].queryLength;
//    TaxID selectedTaxon;
//    if(PRINT) {
//        cout<<"# "<<currentQuery<<endl;
//        for (size_t i = offset; i < end + 1; i++) {
//            cout << matchList[i].genusTaxID<<" "<<matchList[i].speciesTaxID << " " << matchList[i].taxID << " " <<
//                 int(matchList[i].frame)<< " " << matchList[i].position  << " " << int(matchList[i].hamming) << endl;
//        }
//    }
//
//    // Get the best genus for current query
//    vector<size_t> genusMatches;
//    genusMatches.reserve(end-offset+1);
//    float highRankScore;
//    int res;
//    if(par.seqMode == 2){
//        res = getMatchesOfTheBestGenus_paired(matchesForLCA,
//                                              matchList,
//                                              end,
//                                              offset,
//                                              queryList[currentQuery].queryLength,
//                                              queryList[currentQuery].queryLength2,
//                                              highRankScore);
//    } else {
//        res = getMatchesOfTheBestGenus_index(genusMatches, matchList, end, offset, queryLength, highRankScore);
//    }
//
//    if(PRINT) {
//        cout<<"# "<<currentQuery<<" filtered"<<endl;
//        for (size_t i = 0; i < genusMatches.size(); i++) {
//            cout << matchList[genusMatches[i]].genusTaxID<<" "<< matchList[genusMatches[i]].speciesTaxID << " " <<  matchList[genusMatches[i]].taxID << " " <<
//                 int( matchList[genusMatches[i]].frame)<< " " <<  matchList[genusMatches[i]].position  << " " << int( matchList[genusMatches[i]].hamming) << endl;
//        }
//    }
//
////    if(PRINT){
////        sort(matchesForLCA.begin(), matchesForLCA.end(), Classifier::sortByGenusAndSpecies);
////    }
//
//    //If there is no proper genus for current query, it is un-classified.
//    if(res == 3){
//        queryList[currentQuery].isClassified = false;
//        queryList[currentQuery].classification = 0;
//        queryList[currentQuery].score = 0;
//        queryList[currentQuery].newSpecies = false;
//        return 0;
//    }
//
//    for(size_t i = 0; i < genusMatches.size(); i++ ){
//        queryList[currentQuery].taxCnt[matchList[genusMatches[i]].taxID] ++;
//    }
//
//    // If there are two or more good genus level candidates, find the LCA.
//    if(res == 2){
//        vector<TaxID> taxIdList;
//        taxIdList.reserve(genusMatches.size());
//        for(size_t i = 0; i < genusMatches.size(); i++ ){
//            taxIdList.push_back(matchList[genusMatches[i]].taxID);
//        }
//        selectedTaxon = ncbiTaxonomy.LCA(taxIdList)->taxId;
//        queryList[currentQuery].isClassified = true;
//        queryList[currentQuery].classification = selectedTaxon;
//        queryList[currentQuery].score = highRankScore;
//        if(PRINT) {
//            cout << "# " << currentQuery << " " << res << endl;
//            for (size_t i = 0; i < taxIdList.size(); i++) {
//                cout << i << " " << int(matchList[genusMatches[i]].frame) << " " << matchList[genusMatches[i]].position<< " " <<
//                matchList[genusMatches[i]].taxID << " " << int(matchList[genusMatches[i]].hamming) <<" "<< matchList[genusMatches[i]].red << endl;
//            }
//            cout << "Score: " << highRankScore << " " << selectedTaxon << " "
//                 << ncbiTaxonomy.taxonNode(selectedTaxon)->rank << endl;
//        }
//        return selectedTaxon;
//    }
//
//    // Choose the species with the highest coverage.
//    TaxID selectedSpecies;
//    ScrCov speciesScrCov(0.f, 0.f);
//    vector<TaxID> species;
//    if(par.seqMode == 2) {
//        classifyFurther_paired(matchesForLCA,
//                               ncbiTaxonomy,
//                               queryLength,
//                               queryList[currentQuery].queryLength2,
//                               speciesScrCov,
//                               species);
//    } else {
//        chooseSpecies(matchesForLCA,
//                      ncbiTaxonomy,
//                      queryLength,
//                      (float) queryList[currentQuery].kmerCnt / 6,
//                      speciesScrCov,
//                      species);
//    }
//
//    // Classify at the genus rank if more than one species are selected.
//    // Classify at the genus rank if the score at species level is not enough.
//    if(species.size() > 1
//       || (speciesScrCov.score < 0.9 && !ncbiTaxonomy.IsAncestor(par.virusTaxId, matchList[genusMatches[0]].taxID))){
//        queryList[currentQuery].isClassified = true;
//        queryList[currentQuery].classification = ncbiTaxonomy.getTaxIdAtRank(matchList[genusMatches[0]].taxID, "genus");
//        queryList[currentQuery].score = highRankScore;
//        return queryList[currentQuery].classification;
//    }
//
//    selectedSpecies = species[0];
//
//    // Check if it can be classified at the subspecies rank.
//    int numOfstrains = 0;
//    TaxID strainID = 0;
//    int count = 1;
//    int minStrainSpecificCnt = 1;
//    if(par.seqMode == 1){
//        minStrainSpecificCnt = 1;
//    } else if (par.seqMode == 2){
//        minStrainSpecificCnt = 2;
//    } else if (par.seqMode == 3){
//        minStrainSpecificCnt = 5;
//        if(queryLength > 5000){
//            minStrainSpecificCnt = queryLength / 1000;
//        }
//    }
//    if(NcbiTaxonomy::findRankIndex(ncbiTaxonomy.taxonNode(selectedSpecies)->rank) == 4){
//        unordered_map<TaxID, int> strainMatchCnt;
//        for(size_t i = 0; i < genusMatches.size(); i++){
//            if(selectedSpecies != matchList[genusMatches[i]].taxID
//               && ncbiTaxonomy.IsAncestor(selectedSpecies, matchList[genusMatches[i]].taxID)){
//                strainMatchCnt[matchesForLCA[i].taxID]++;
//            }
//        }
//        for(auto strainIt = strainMatchCnt.begin(); strainIt != strainMatchCnt.end(); strainIt ++){
//            if(strainIt->second > minStrainSpecificCnt){
//                strainID = strainIt->first;
//                numOfstrains++;
//                count = strainIt->second;
//            }
//        }
//        if(numOfstrains == 1 && count > minStrainSpecificCnt + 1) {
//            selectedSpecies = strainID;
//        }
//    }
//
//    // Store classification results
//    queryList[currentQuery].isClassified = true;
//    queryList[currentQuery].classification = selectedSpecies;
//    queryList[currentQuery].score = speciesScrCov.score;
//    queryList[currentQuery].newSpecies = false;
//
//    if(PRINT) {
//        cout << "# " << currentQuery << endl;
//        for (size_t i = 0; i < matchesForLCA.size(); i++) {
//            cout << i << " " << int(matchesForLCA[i].frame) << " " << matchesForLCA[i].position<< " " <<
//                 matchesForLCA[i].taxID << " " << int(matchesForLCA[i].hamming) <<" "<< matchesForLCA[i].red << endl;
//        }
//        cout << "Score: " << speciesScrCov.score << "  " << selectedSpecies << " " << ncbiTaxonomy.taxonNode(selectedSpecies)->rank
//             << endl;
//    }
//
//    return selectedSpecies;
//}

int Classifier::getMatchesOfTheBestGenus_paired(vector<Match> &matchesForMajorityLCA, Match *matchList, size_t end,
                                                size_t offset, int readLength1, int readLength2, float &bestScore) {
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
    while (i < end + 1) {
        currentGenus = matchList[i].genusTaxID;
        // For current genus
        while (currentGenus == matchList[i].genusTaxID && (i < end + 1)) {
            currentSpecies = matchList[i].speciesTaxID;
            // For current species
            // Filter un-consecutive matches (probably random matches)
            while (currentSpecies == matchList[i].speciesTaxID && (i < end + 1)) {
                offsetIdx = i;
                i++;
                newOffset = true;
                hammingSum = 0;
                hammingMean = matchList[i - 1].hamming;
                while ((i < end + 1) && currentSpecies == matchList[i].speciesTaxID &&
                       matchList[i].position <= matchList[i - 1].position + 3) {
                    if (newOffset && ((float) matchList[i].hamming) <= hammingMean + 3) {
                        newOffset = false;
                        hammingSum = matchList[offsetIdx].hamming;
                        filteredMatches.push_back(matchList[offsetIdx]);
                        continue;
                    } else if (!newOffset && ((float) matchList[i].hamming) <= hammingMean + 3) {
                        filteredMatches.push_back(matchList[i]);
                        hammingSum += matchList[i].hamming;
                        hammingMean = float(hammingSum) / float(filteredMatches.size());
                    }
                    i++;
                }
            }
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
        return 3;
    }

    float maxScore = *max_element(scoreOfEachGenus.begin(), scoreOfEachGenus.end());
    if (maxScore < 0.3)
        return 3;
    vector<size_t> maxIdx;
    for (size_t g = 0; g < scoreOfEachGenus.size(); g++) {
        if (scoreOfEachGenus[g] > maxScore * 0.9f) {
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

int Classifier::getMatchesOfTheBestGenus_index(vector<size_t> &matchesForMajorityLCA,
                                               Match *matchList,
                                               size_t end,
                                               size_t offset,
                                               int queryLength,
                                               float &bestScore) {
    int conCnt;
    uint32_t hammingSum;
    float hammingMean;
    TaxID currentGenus;
    TaxID currentSpecies;

    vector<size_t> filteredMatches;
    vector<vector<size_t>> matchesForEachGenus;
    vector<bool> conservedWithinGenus;
    vector<float> scoreOfEachGenus;
    size_t i = offset;
    size_t offsetIdx;
    bool newOffset;
    while (i < end + 1) {
        currentGenus = matchList[i].genusTaxID;
        // For current genus
        while (currentGenus == matchList[i].genusTaxID && (i < end + 1)) {
            currentSpecies = matchList[i].speciesTaxID;
            // For current species
            // Filter un-consecutive matches (probably random matches)
            while (currentSpecies == matchList[i].speciesTaxID && (i < end + 1)) {
                offsetIdx = i;
                i++;
                newOffset = true;
                hammingSum = 0;
                hammingMean = matchList[i - 1].hamming;
                while ((i < end + 1) && currentSpecies == matchList[i].speciesTaxID &&
                       matchList[i].position <= matchList[i - 1].position + 3) {
                    if (newOffset && ((float) matchList[i].hamming) <= hammingMean + 3) {
                        newOffset = false;
                        hammingSum = matchList[offsetIdx].hamming;
                        filteredMatches.push_back(offsetIdx);
                        continue;
                    } else if (!newOffset && ((float) matchList[i].hamming) <= hammingMean + 3) {
                        filteredMatches.push_back(i);
                        hammingSum += matchList[i].hamming;
                        hammingMean = float(hammingSum) / float(filteredMatches.size());
                    }
                    i++;
                }
            }
        }

        // Construct a match combination using filtered matches of current genus
        // so that it can best cover the query, and score the combination
        if (!filteredMatches.empty()) {
            constructMatchCombination_index(filteredMatches,
                                            matchList,
                                            matchesForEachGenus,
                                            scoreOfEachGenus,
                                            queryLength);
        }
        filteredMatches.clear();
    }

    // If there are no meaningful genus
    if (scoreOfEachGenus.empty()) {
        bestScore = 0;
        return 3;
    }

    float maxScore = *max_element(scoreOfEachGenus.begin(), scoreOfEachGenus.end());
    if (maxScore < 0.3)
        return 3;
    vector<size_t> maxIdx;
    for (size_t g = 0; g < scoreOfEachGenus.size(); g++) {
        if (scoreOfEachGenus[g] > maxScore * 0.9f) {
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
    while (i < end + 1) {
        currentGenus = matchList[i].genusTaxID;
        // For current genus
        while ((i < end + 1) && currentGenus == matchList[i].genusTaxID) {
            currentSpecies = matchList[i].speciesTaxID;
            // For current species
            // Filter un-consecutive matches (probably random matches)
            while ((i < end + 1) && currentSpecies == matchList[i].speciesTaxID) {
                offsetIdx = i;
                i++;
                newOffset = true;
                hammingSum = 0;
                hammingMean = matchList[i - 1].hamming;
                while ((i < end + 1) && currentSpecies == matchList[i].speciesTaxID &&
                       matchList[i].position <= matchList[i - 1].position + 3) {
                    if (newOffset && ((float) matchList[i].hamming) <= hammingMean + 3) {
                        newOffset = false;
                        hammingSum = matchList[offsetIdx].hamming;
                        filteredMatches.push_back(matchList[offsetIdx]);
                        continue;
                    } else if (!newOffset && ((float) matchList[i].hamming) <= hammingMean + 3) {
                        filteredMatches.push_back(matchList[i]);
                        hammingSum += matchList[i].hamming;
                        hammingMean = float(hammingSum) / float(filteredMatches.size());
                    }
                    i++;
                }
            }
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
    if (maxScore < 0.3)
        return 3;
    vector<size_t> maxIdx;
    for (size_t g = 0; g < scoreOfEachGenus.size(); g++) {
        if (scoreOfEachGenus[g] > maxScore * 0.9f) {
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

void Classifier::constructMatchCombination_index(vector<size_t> &index,
                                                 Match *matchList,
                                                 vector<vector<size_t>> &matchesForEachGenus,
                                                 vector<float> &scoreOfEachGenus,
                                                 int queryLength) {
    // Do not allow overlaps between the same species
    vector<size_t> m;
    size_t walker = 0;
    size_t numOfFitMat = index.size();
    Match *currentMatch = matchList;
    size_t currIndex;
    while (walker < numOfFitMat) {
        TaxID currentSpecies = matchList[index[walker]].speciesTaxID;
        int currentPosition = matchList[index[walker]].position / 3;
        currentMatch = matchList + index[walker];
        currIndex = walker;
        // Look through overlaps within a species
        while (matchList[index[walker]].speciesTaxID == currentSpecies
               && matchList[index[walker]].position / 3 == currentPosition
               && walker < numOfFitMat) {
            // Take the match with lower hamming distance
            if (matchList[index[walker]].hamming < currentMatch->hamming) {
                currentMatch = matchList + index[walker];
                currIndex = walker;
            }
                // Overlapping with same hamming distance but different subspecies taxonomy ID -> species level
            else if (currentMatch->taxID != matchList[index[walker]].taxID &&
                     currentMatch->hamming == matchList[index[walker]].hamming) {
                currentMatch->taxID = currentMatch->speciesTaxID;
            }
            walker++;
        }
        m.push_back(currIndex);
    }

    // Calculate Hamming distance & covered length
    int coveredPosCnt = 0;
    uint16_t currHammings;
    int aminoAcidNum = (int) queryLength / 3;
    int currPos;
    size_t matchNum = m.size();
    size_t f = 0;

    // Get the largest hamming distance at each position of query
    auto *hammingsAtEachPos = new signed char[aminoAcidNum + 1];
    memset(hammingsAtEachPos, -1, (aminoAcidNum + 1));
    while (f < matchNum) {
        currPos = matchList[m[f]].position / 3;
        currHammings = matchList[m[f]].rightEndHamming;
        if (GET_2_BITS(currHammings) > hammingsAtEachPos[currPos])
            hammingsAtEachPos[currPos] = GET_2_BITS(currHammings);
        if (GET_2_BITS(currHammings >> 2) > hammingsAtEachPos[currPos + 1])
            hammingsAtEachPos[currPos + 1] = GET_2_BITS(currHammings >> 2);
        if (GET_2_BITS(currHammings >> 4) > hammingsAtEachPos[currPos + 2])
            hammingsAtEachPos[currPos + 2] = GET_2_BITS(currHammings >> 4);
        if (GET_2_BITS(currHammings >> 6) > hammingsAtEachPos[currPos + 3])
            hammingsAtEachPos[currPos + 3] = GET_2_BITS(currHammings >> 6);
        if (GET_2_BITS(currHammings >> 8) > hammingsAtEachPos[currPos + 4])
            hammingsAtEachPos[currPos + 4] = GET_2_BITS(currHammings >> 8);
        if (GET_2_BITS(currHammings >> 10) > hammingsAtEachPos[currPos + 5])
            hammingsAtEachPos[currPos + 5] = GET_2_BITS(currHammings >> 10);
        if (GET_2_BITS(currHammings >> 12) > hammingsAtEachPos[currPos + 6])
            hammingsAtEachPos[currPos + 6] = GET_2_BITS(currHammings >> 12);
        if (GET_2_BITS(currHammings >> 14) > hammingsAtEachPos[currPos + 7])
            hammingsAtEachPos[currPos + 7] = GET_2_BITS(currHammings >> 14);
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
    matchesForEachGenus.push_back(m);

    if (PRINT) {
        cout << matchList[index[0]].genusTaxID << " " << coveredLength << " " << hammingSum << " "
             << ((float) coveredLength - hammingSum) / (float) queryLength <<
             " " << m.size()
             << endl;
    }
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

//    // Get the smallest hamming distance at each position of query
//    auto * hammingsAtEachPos = new signed char[aminoAcidNum + 1];
//    memset(hammingsAtEachPos, 10, (aminoAcidNum + 1));
//    while(f < matchNum){
//        currPos = matches[f].position / 3;
//        currHammings = matches[f].rightEndHamming;
//        if(GET_2_BITS(currHammings) < hammingsAtEachPos[currPos]) hammingsAtEachPos[currPos] = GET_2_BITS(currHammings);
//        if(GET_2_BITS(currHammings>>2) < hammingsAtEachPos[currPos + 1]) hammingsAtEachPos[currPos + 1] = GET_2_BITS(currHammings>>2);
//        if(GET_2_BITS(currHammings>>4) < hammingsAtEachPos[currPos + 2]) hammingsAtEachPos[currPos + 2] = GET_2_BITS(currHammings>>4);
//        if(GET_2_BITS(currHammings>>6) < hammingsAtEachPos[currPos + 3]) hammingsAtEachPos[currPos + 3] = GET_2_BITS(currHammings>>6);
//        if(GET_2_BITS(currHammings>>8) < hammingsAtEachPos[currPos + 4]) hammingsAtEachPos[currPos + 4] = GET_2_BITS(currHammings>>8);
//        if(GET_2_BITS(currHammings>>10) < hammingsAtEachPos[currPos + 5]) hammingsAtEachPos[currPos + 5] = GET_2_BITS(currHammings>>10);
//        if(GET_2_BITS(currHammings>>12) < hammingsAtEachPos[currPos + 6]) hammingsAtEachPos[currPos + 6] = GET_2_BITS(currHammings>>12);
//        if(GET_2_BITS(currHammings>>14) < hammingsAtEachPos[currPos + 7]) hammingsAtEachPos[currPos + 7] = GET_2_BITS(currHammings>>14);
//        f++;
//    }
//
//    // Sum up hamming distances and count the number of position covered by the matches.
//    float hammingSum = 0;
//    for(int h = 0; h < aminoAcidNum; h++){
//        if(hammingsAtEachPos[h] == 0) { // Add 0 for 0 hamming dist.
//            coveredPosCnt++;
//        }else if(hammingsAtEachPos[h] != 10){ // Add 1.5, 2, 2.5 for 1, 2, 3 hamming dist. respectively
//            hammingSum += 1.0f + (0.5f * hammingsAtEachPos[h]);
//            coveredPosCnt ++;
//        }
//    }

    // Get the largest hamming distance at each position of query
    auto *hammingsAtEachPos = new signed char[aminoAcidNum + 1];
    memset(hammingsAtEachPos, -1, (aminoAcidNum + 1));
    while (f < matchNum) {
        currPos = matches[f].position / 3;
        currHammings = matches[f].rightEndHamming;
        if (GET_2_BITS(currHammings) > hammingsAtEachPos[currPos])
            hammingsAtEachPos[currPos] = GET_2_BITS(currHammings);
        if (GET_2_BITS(currHammings >> 2) > hammingsAtEachPos[currPos + 1])
            hammingsAtEachPos[currPos + 1] = GET_2_BITS(currHammings >> 2);
        if (GET_2_BITS(currHammings >> 4) > hammingsAtEachPos[currPos + 2])
            hammingsAtEachPos[currPos + 2] = GET_2_BITS(currHammings >> 4);
        if (GET_2_BITS(currHammings >> 6) > hammingsAtEachPos[currPos + 3])
            hammingsAtEachPos[currPos + 3] = GET_2_BITS(currHammings >> 6);
        if (GET_2_BITS(currHammings >> 8) > hammingsAtEachPos[currPos + 4])
            hammingsAtEachPos[currPos + 4] = GET_2_BITS(currHammings >> 8);
        if (GET_2_BITS(currHammings >> 10) > hammingsAtEachPos[currPos + 5])
            hammingsAtEachPos[currPos + 5] = GET_2_BITS(currHammings >> 10);
        if (GET_2_BITS(currHammings >> 12) > hammingsAtEachPos[currPos + 6])
            hammingsAtEachPos[currPos + 6] = GET_2_BITS(currHammings >> 12);
        if (GET_2_BITS(currHammings >> 14) > hammingsAtEachPos[currPos + 7])
            hammingsAtEachPos[currPos + 7] = GET_2_BITS(currHammings >> 14);
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
        while (filteredMatches[walker].speciesTaxID == currentSpecies
               && filteredMatches[walker].position / 3 == currentPosition
               && walker < numOfFitMat) {
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

//    // Get the smallest hamming distance at each position of query
//    auto * hammingsAtEachPos = new signed char[aminoAcidNum + 1];
//    memset(hammingsAtEachPos, 10, (aminoAcidNum + 1));
//    while(f < matchNum){
//        currPos = matches[f].position / 3;
//        currHammings = matches[f].rightEndHamming;
//        if(GET_2_BITS(currHammings) < hammingsAtEachPos[currPos]) hammingsAtEachPos[currPos] = GET_2_BITS(currHammings);
//        if(GET_2_BITS(currHammings>>2) < hammingsAtEachPos[currPos + 1]) hammingsAtEachPos[currPos + 1] = GET_2_BITS(currHammings>>2);
//        if(GET_2_BITS(currHammings>>4) < hammingsAtEachPos[currPos + 2]) hammingsAtEachPos[currPos + 2] = GET_2_BITS(currHammings>>4);
//        if(GET_2_BITS(currHammings>>6) < hammingsAtEachPos[currPos + 3]) hammingsAtEachPos[currPos + 3] = GET_2_BITS(currHammings>>6);
//        if(GET_2_BITS(currHammings>>8) < hammingsAtEachPos[currPos + 4]) hammingsAtEachPos[currPos + 4] = GET_2_BITS(currHammings>>8);
//        if(GET_2_BITS(currHammings>>10) < hammingsAtEachPos[currPos + 5]) hammingsAtEachPos[currPos + 5] = GET_2_BITS(currHammings>>10);
//        if(GET_2_BITS(currHammings>>12) < hammingsAtEachPos[currPos + 6]) hammingsAtEachPos[currPos + 6] = GET_2_BITS(currHammings>>12);
//        if(GET_2_BITS(currHammings>>14) < hammingsAtEachPos[currPos + 7]) hammingsAtEachPos[currPos + 7] = GET_2_BITS(currHammings>>14);
//        f++;
//    }
//
//    // Sum up hamming distances and count the number of position covered by the matches.
//    float hammingSum = 0;
//    for(int h = 0; h < aminoAcidNum; h++){
//        if(hammingsAtEachPos[h] == 0) { // Add 0 for 0 hamming dist.
//            coveredPosCnt++;
//        }else if(hammingsAtEachPos[h] != 10){ // Add 1.5, 2, 2.5 for 1, 2, 3 hamming dist. respectively
//            hammingSum += 1.0f + (0.5f * hammingsAtEachPos[h]);
//            coveredPosCnt ++;
//        }
//    }

    // Get the largest hamming distance at each position of query
    auto *hammingsAtEachPos = new signed char[aminoAcidNum_total + 1];
    memset(hammingsAtEachPos, -1, (aminoAcidNum_total + 1));
    while (f < matchNum) {
        currPos = matches[f].position / 3;
        currHammings = matches[f].rightEndHamming;
        if (GET_2_BITS(currHammings) > hammingsAtEachPos[currPos])
            hammingsAtEachPos[currPos] = GET_2_BITS(currHammings);
        if (GET_2_BITS(currHammings >> 2) > hammingsAtEachPos[currPos + 1])
            hammingsAtEachPos[currPos + 1] = GET_2_BITS(currHammings >> 2);
        if (GET_2_BITS(currHammings >> 4) > hammingsAtEachPos[currPos + 2])
            hammingsAtEachPos[currPos + 2] = GET_2_BITS(currHammings >> 4);
        if (GET_2_BITS(currHammings >> 6) > hammingsAtEachPos[currPos + 3])
            hammingsAtEachPos[currPos + 3] = GET_2_BITS(currHammings >> 6);
        if (GET_2_BITS(currHammings >> 8) > hammingsAtEachPos[currPos + 4])
            hammingsAtEachPos[currPos + 4] = GET_2_BITS(currHammings >> 8);
        if (GET_2_BITS(currHammings >> 10) > hammingsAtEachPos[currPos + 5])
            hammingsAtEachPos[currPos + 5] = GET_2_BITS(currHammings >> 10);
        if (GET_2_BITS(currHammings >> 12) > hammingsAtEachPos[currPos + 6])
            hammingsAtEachPos[currPos + 6] = GET_2_BITS(currHammings >> 12);
        if (GET_2_BITS(currHammings >> 14) > hammingsAtEachPos[currPos + 7])
            hammingsAtEachPos[currPos + 7] = GET_2_BITS(currHammings >> 14);
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


//TODO kmer count -> covered length
//TODO hamming
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
            hammingsAtEachPos[currPos] = GET_2_BITS(currHammings);
        if (GET_2_BITS(currHammings >> 2) > hammingsAtEachPos[currPos + 1])
            hammingsAtEachPos[currPos + 1] = GET_2_BITS(currHammings >> 2);
        if (GET_2_BITS(currHammings >> 4) > hammingsAtEachPos[currPos + 2])
            hammingsAtEachPos[currPos + 2] = GET_2_BITS(currHammings >> 4);
        if (GET_2_BITS(currHammings >> 6) > hammingsAtEachPos[currPos + 3])
            hammingsAtEachPos[currPos + 3] = GET_2_BITS(currHammings >> 6);
        if (GET_2_BITS(currHammings >> 8) > hammingsAtEachPos[currPos + 4])
            hammingsAtEachPos[currPos + 4] = GET_2_BITS(currHammings >> 8);
        if (GET_2_BITS(currHammings >> 10) > hammingsAtEachPos[currPos + 5])
            hammingsAtEachPos[currPos + 5] = GET_2_BITS(currHammings >> 10);
        if (GET_2_BITS(currHammings >> 12) > hammingsAtEachPos[currPos + 6])
            hammingsAtEachPos[currPos + 6] = GET_2_BITS(currHammings >> 12);
        if (GET_2_BITS(currHammings >> 14) > hammingsAtEachPos[currPos + 7])
            hammingsAtEachPos[currPos + 7] = GET_2_BITS(currHammings >> 14);
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
    auto *hammingsAtEachPos = new signed char[aminoAcidNum_total + 1];
    memset(hammingsAtEachPos, -1, (aminoAcidNum_total + 1));

    int currPos;
    size_t walker = begin;
    uint16_t currHammings;
    while (walker < end) {
        currPos = matches[walker].position / 3;
        currHammings = matches[walker].rightEndHamming;
        if (GET_2_BITS(currHammings) > hammingsAtEachPos[currPos])
            hammingsAtEachPos[currPos] = GET_2_BITS(currHammings);
        if (GET_2_BITS(currHammings >> 2) > hammingsAtEachPos[currPos + 1])
            hammingsAtEachPos[currPos + 1] = GET_2_BITS(currHammings >> 2);
        if (GET_2_BITS(currHammings >> 4) > hammingsAtEachPos[currPos + 2])
            hammingsAtEachPos[currPos + 2] = GET_2_BITS(currHammings >> 4);
        if (GET_2_BITS(currHammings >> 6) > hammingsAtEachPos[currPos + 3])
            hammingsAtEachPos[currPos + 3] = GET_2_BITS(currHammings >> 6);
        if (GET_2_BITS(currHammings >> 8) > hammingsAtEachPos[currPos + 4])
            hammingsAtEachPos[currPos + 4] = GET_2_BITS(currHammings >> 8);
        if (GET_2_BITS(currHammings >> 10) > hammingsAtEachPos[currPos + 5])
            hammingsAtEachPos[currPos + 5] = GET_2_BITS(currHammings >> 10);
        if (GET_2_BITS(currHammings >> 12) > hammingsAtEachPos[currPos + 6])
            hammingsAtEachPos[currPos + 6] = GET_2_BITS(currHammings >> 12);
        if (GET_2_BITS(currHammings >> 14) > hammingsAtEachPos[currPos + 7])
            hammingsAtEachPos[currPos + 7] = GET_2_BITS(currHammings >> 14);
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
//    if (coveredLength_read1 >= queryLength) coveredLength_read1 = queryLength;
//    if (coveredLength_read2 >= queryLength2) coveredLength_read2 = queryLength2;
    if (coveredLength_read1 >= queryLength) coveredLength_read1 = queryLength;
    if (coveredLength_read2 >= queryLength2) coveredLength_read2 = queryLength2;

    return {((float) coveredLength_read1 + coveredLength_read2 - hammingSum) / ((float) queryLength + queryLength2),
            ((float) coveredLength_read1 + coveredLength_read2) / ((float) queryLength + queryLength2)};
}

int Classifier::getNumOfSplits() const { return this->numOfSplit; }

bool Classifier::compareForLinearSearch(const QueryKmer &a, const QueryKmer &b) {
    if (a.ADkmer < b.ADkmer) {
        return true;
    } else if (a.ADkmer == b.ADkmer) {
        return (a.info.sequenceID < b.info.sequenceID);
    }
    return false;
}

bool Classifier::sortByGenusAndSpecies(const Match &a, const Match &b) {
    if (a.queryId < b.queryId) return true;
    else if (a.queryId == b.queryId) {
        if (a.genusTaxID < b.genusTaxID) return true;
        else if (a.genusTaxID == b.genusTaxID) {
            if (a.speciesTaxID < b.speciesTaxID) return true;
            else if (a.speciesTaxID == b.speciesTaxID) {
                if (a.position < b.position) return true;
                else if (a.position == b.position) {
                    if (a.frame < b.frame) return true;
                }
            }
        }
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
                               << queryList[i].classification << "\t" << queryList[i].queryLength << "\t"
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

void Classifier::performanceTest(NcbiTaxonomy &ncbiTaxonomy, Query *queryList, int numOfquery, vector<int> &wrongs) {

    ///Load the mapping file
    const char *mappingFile = "../../gtdb_taxdmp/assacc_to_taxid_gtdb.tsv";
    unordered_map<string, int> assacc2taxid;
    string key, value;
    ifstream map;
    map.open(mappingFile);
    if (map.is_open()) {
        while (getline(map, key, '\t')) {
            getline(map, value, '\n');
            assacc2taxid[key] = stoi(value);
        }
    } else {
        cout << "Cannot open file for mappig from assemlby accession to tax ID" << endl;
    }
    map.close();

    regex regex1("(GC[AF]_[0-9]*\\.[0-9]*)");
    smatch assacc;
    TaxID classificationResult;
    int rightAnswer;
    string queryName;


    counts.classificationCnt = 0;
    counts.correct = 0;
    counts.highRank = 0;

    //number of targets at each rank
    counts.subspeciesTargetNumber = 0;
    counts.speciesTargetNumber = 0;
    counts.genusTargetNumber = 0;
    counts.familyTargetNumber = 0;
    counts.orderTargetNumber = 0;
    counts.classTargetNumber = 0;
    counts.phylumTargetNumber = 0;
    counts.superkingdomTargetNumber = 0;

    //number of classification at each rank
    counts.subspeciesCnt_try = 0;
    counts.speciesCnt_try = 0;
    counts.genusCnt_try = 0;
    counts.familyCnt_try = 0;
    counts.orderCnt_try = 0;
    counts.classCnt_try = 0;
    counts.phylumCnt_try = 0;
    counts.superkingdomCnt_try = 0;


    //number of correct classifications at each rank
    counts.subspeciesCnt_correct = 0;
    counts.speciesCnt_correct = 0;
    counts.genusCnt_correct = 0;
    counts.familyCnt_correct = 0;
    counts.orderCnt_correct = 0;
    counts.classCnt_correct = 0;
    counts.phylumCnt_correct = 0;
    counts.superkingdomCnt_correct = 0;

    for (int i = 0; i < numOfquery; i++) {
        classificationResult = queryList[i].classification;
        if (classificationResult == 0) {
            cout << "here " << i << endl;
            continue;
        } else {
            classifiedCnt++;
            queryName = queryList[i].name;
            regex_search(queryName, assacc, regex1);
            if (assacc2taxid.count(assacc[0].str())) {
                rightAnswer = assacc2taxid[assacc[0].str()];
            } else {
                cout << assacc[0].str() << " is not in the mapping file" << endl;
                continue;
            }
            //cout<<"compareTaxon"<<" "<<i<<endl;
            cout << i << " ";
            compareTaxon(classificationResult, rightAnswer, ncbiTaxonomy, wrongs, i);

        }
    }

    cout << "Num of queries: " << queryInfos.size() << endl;
    cout << "Num of classifications: " << counts.classificationCnt << endl;
    cout << "Num of correct classifications: " << counts.correct << endl;
    cout << "Num of correct but too broad classifications: " << counts.highRank << endl;
    cout << "classified/total = " << float(counts.classificationCnt) / float(queryInfos.size()) << endl;
    cout << "correct   /total = " << float(counts.correct) / float(queryInfos.size()) << endl;
    cout << "correct   /classifications = " << float(counts.correct) / float(counts.classificationCnt) << endl;
    cout << "high rank /classifications = " << float(counts.highRank) / float(counts.classificationCnt) << endl << endl;

    cout << "Number of targets at each rank / correct classification / tries" << endl;
    cout << "Superkingdom: " << counts.superkingdomTargetNumber << " / " << counts.superkingdomCnt_correct << " / "
         << counts.superkingdomCnt_try << endl;
    cout << "Phylum      : " << counts.phylumTargetNumber << " / " << counts.phylumCnt_correct << " / "
         << counts.phylumCnt_try << endl;
    cout << "Class       : " << counts.classTargetNumber << " / " << counts.classCnt_correct << " / "
         << counts.classCnt_try << endl;
    cout << "Order       : " << counts.orderTargetNumber << " / " << counts.orderCnt_correct << " / "
         << counts.orderCnt_try << endl;
    cout << "Family      : " << counts.familyTargetNumber << " / " << counts.familyCnt_correct << " / "
         << counts.familyCnt_try << endl;
    cout << "Genus       : " << counts.genusTargetNumber << " / " << counts.genusCnt_correct << " / "
         << counts.genusCnt_try << endl;
    cout << "Species     : " << counts.speciesTargetNumber << " / " << counts.speciesCnt_correct << " / "
         << counts.speciesCnt_try << endl;
    cout << "Subspecies  : " << counts.subspeciesTargetNumber << " / " << counts.subspeciesCnt_correct << " / "
         << counts.subspeciesCnt_try << endl;
}

void Classifier::compareTaxon(TaxID shot, TaxID target, NcbiTaxonomy &ncbiTaxonomy, vector<int> &wrongs,
                              int i) { ///target: subspecies or species
    const TaxonNode *shotNode = ncbiTaxonomy.taxonNode(shot);
    const TaxonNode *targetNode = ncbiTaxonomy.taxonNode(target);
    string shotRank = shotNode->rank;
    string targetRank = targetNode->rank;
    cout << shot << " " << target << " " << shotRank << " " << targetRank << " ";

    if (shot == 0) {
        if (shotRank == "species" || shotRank == "subspecies")
            wrongs.push_back(i);
        cout << i << endl;
        cout << "X" << endl;
        return;
    } else {
        counts.classificationCnt++;
    }

    bool isCorrect = false;
    if (shot == target) {
        counts.correct++;
        isCorrect = true;
        cout << "O" << endl;
    } else if (NcbiTaxonomy::findRankIndex(shotRank) <=
               NcbiTaxonomy::findRankIndex(targetRank)) { //classified into wrong taxon or too specifically
        if (shotRank == "species" || shotRank == "subspecies")
            wrongs.push_back(i);
        cout << "X" << endl;
    } else { // classified at higher rank (too safe classification)
        if (shotRank == "superkingdom") {
            if (shotRank == "species" || shotRank == "subspecies")
                wrongs.push_back(i);
            cout << "X" << endl;
        } else if (shot == ncbiTaxonomy.getTaxIdAtRank(target, shotRank)) { //on right branch
            counts.correct++;
            cout << "O" << endl;
            isCorrect = true;
        } else { //on wrong branch
            if (shotRank == "species" || shotRank == "subspecies") wrongs.push_back(i);
            cout << "X" << endl;
        }
    }

    //count the number of classification at each rank
    if (shotRank == "subspecies") {
        counts.subspeciesCnt_try++;
    } else if (shotRank == "species") {
        counts.speciesCnt_try++;
    } else if (shotRank == "genus") {
        counts.genusCnt_try++;
    } else if (shotRank == "family") {
        counts.familyCnt_try++;
    } else if (shotRank == "order") {
        counts.orderCnt_try++;
    } else if (shotRank == "class") {
        counts.classCnt_try++;
    } else if (shotRank == "phylum") {
        counts.phylumCnt_try++;
    } else if (shotRank == "superkingdom") {
        counts.superkingdomCnt_try++;
    }

    if (!isCorrect) return;

    //count the number of correct classification at each rank
    if (shotRank == "subspecies") {
        counts.subspeciesCnt_correct++;
    } else if (shotRank == "species") {
        counts.speciesCnt_correct++;
    } else if (shotRank == "genus") {
        counts.genusCnt_correct++;
    } else if (shotRank == "family") {
        counts.familyCnt_correct++;
    } else if (shotRank == "order") {
        counts.orderCnt_correct++;
    } else if (shotRank == "class") {
        counts.classCnt_correct++;
    } else if (shotRank == "phylum") {
        counts.phylumCnt_correct++;
    } else if (shotRank == "superkingdom") {
        counts.superkingdomCnt_correct++;
    }
}


