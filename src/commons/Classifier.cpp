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

Classifier::~Classifier() { }//delete seqIterator; }

void Classifier::startClassify(const char * queryFileName,
                               const char * targetDiffIdxFileName,
                               const char * targetInfoFileName,
                               const char * diffIdxSplitFileName,
                               vector<int> & taxIdList,
                               const LocalParameters & par,
                               NcbiTaxonomy & taxonomy) {
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
    sprintf(matchFileName,"%s_match2", queryFileName);
    FILE * matchFile = fopen(matchFileName, "wb");

//    struct MmapedData<TargetKmerInfo> testInfo = mmapData<TargetKmerInfo>("/data3/jaebeom/onegenome/0716test_0_info");
//    int testNum = testInfo.fileSize/sizeof(TargetKmerInfo);
//    cout<<"hi231"<<endl;
//    for(int i = 0; i < testNum; i++ ){
//        cout<<testInfo.data[i].sequenceID<<endl;
//    }

    // Load database
    struct MmapedData<uint16_t> targetDiffIdxList = mmapData<uint16_t>(targetDiffIdxFileName);
    struct MmapedData<TargetKmerInfo> targetInfoList = mmapData<TargetKmerInfo>(targetInfoFileName);
    struct MmapedData<DiffIdxSplit> diffIdxSplits = mmapData<DiffIdxSplit>(diffIdxSplitFileName);

    // Allocate memory for buffers
    QueryKmerBuffer kmerBuffer(kmerBufSize);
    Buffer<Match> matchBuffer(kmerBufSize);

    // Load query file
    MmapedData<char> queryFile{};
    MmapedData<char> queryFile2{};
    vector<Sequence> sequences;
    vector<Sequence> sequences2;
    Query * queryList;
    size_t numOfSeq;
    if(par.seqMode == 1){
        queryFile = mmapData<char>(par.filenames[0].c_str());
        IndexCreator::getSeqSegmentsWithHead(sequences, queryFile);
        numOfSeq = sequences.size();
        queryList = new Query[numOfSeq];
    } else if (par.seqMode == 2){
        string queryFileName1 = par.filenames[0] + "_1";
        string queryFileName2 = par.filenames[0] + "_2";
        queryFile = mmapData<char>(queryFileName1.c_str());
        queryFile2 = mmapData<char>(queryFileName2.c_str());
        IndexCreator::getSeqSegmentsWithHead(sequences, queryFile);
        IndexCreator::getSeqSegmentsWithHead(sequences2, queryFile2);
        numOfSeq = sequences.size();
        size_t seqCnt2 = sequences2.size();
        if(numOfSeq > seqCnt2){
            queryList = new Query[numOfSeq];
        } else {
            numOfSeq = seqCnt2;
            queryList = new Query[numOfSeq];
        }
    }

    // Checker for multi-threading
    bool * processedSeqChecker = new bool[numOfSeq];
    fill_n(processedSeqChecker, numOfSeq, false);
    size_t processedSeqCnt = 0;


    // Timer
    time_t afterSearch, afterAnalyze, beforeIO;
    size_t numOfTatalQueryKmerCnt = 0;

    // Extract k-mers from query sequences and compare them to target k-mer DB
    omp_set_num_threads(par.threads);
    while(processedSeqCnt < numOfSeq){
        time_t beforeKmerExtraction = time(nullptr);
        if(par.seqMode == 1 || par.seqMode == 3) { // Single-end short-read sequence or long-read sequence
            fillQueryKmerBufferParallel(kmerBuffer, queryFile, sequences, processedSeqChecker, processedSeqCnt,
                                        queryList, par);
        } else if(par.seqMode == 2){
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
        cout<<"Time spent for k-mer extraction: " << double(time(nullptr) - beforeKmerExtraction) << endl;

        time_t beforeQueryKmerSort = time(nullptr);
        SORT_PARALLEL(kmerBuffer.buffer, kmerBuffer.buffer + kmerBuffer.startIndexOfReserve, Classifier::compareForLinearSearch);
        cout<<"Time spent for sorting query k-mer list: " << double(time(nullptr) - beforeQueryKmerSort) << endl;

        time_t beforeSearch = time(nullptr);
        linearSearchParallel(kmerBuffer.buffer, kmerBuffer.startIndexOfReserve, targetDiffIdxList, targetInfoList,
                             diffIdxSplits, matchBuffer, taxIdList, speciesTaxIdList, genusTaxIdList, matchFile, par);
        cout<<"Time spent for linearSearch: " << double(time(nullptr) - beforeSearch) << endl;
    }
    cout<<"Number of query k-mers: "<<numOfTatalQueryKmerCnt<<endl;

    fclose(matchFile);

    // Sort Matches
    time_t beforeSortMatches = time(nullptr);
    SORT_PARALLEL(matchBuffer.buffer, matchBuffer.buffer + matchBuffer.startIndexOfReserve, Classifier::sortByGenusAndSpecies2);
    cout << "Time spent for sorting matches" << double(time(nullptr) - beforeSortMatches) << endl;

    // Analyze matches
    cout<<"Analyse Result ... "<<endl;
    time_t beforeAnalyze = time(nullptr);
    analyseResultParallel(taxonomy, matchBuffer, numOfSeq, queryList, par);
    cout<<"Time spent for analyzing: "<<double(time(nullptr)-beforeAnalyze)<<endl;

    // Write read classification results.
    ofstream readClassificationFile;
    readClassificationFile.open(par.filenames[3]+"/"+par.filenames[4]+"_ReadClassification.tsv");
    writeReadClassification(queryList,numOfSeq,readClassificationFile);
    readClassificationFile.close();

    // Write a composition report file.
    writeReportFile(par.filenames[3]+"/"+par.filenames[4]+"_CompositionReport.tsv", taxonomy, numOfSeq);

//    // Below is for developing
//    vector<int> wrongClassifications;
//    sequences.clear();
//    IndexCreator::getSeqSegmentsWithHead(sequences, queryFile);
//    performanceTest(taxonomy, queryList, numOfSeq, wrongClassifications);
//    ofstream wr;
//    wr.open(par.filenames[0]+"_wrong");
//
//    for (size_t i = 0; i < wrongClassifications.size(); i++) {
//            kseq_buffer_t buffer(const_cast<char *>(&queryFile.data[sequences[wrongClassifications[i]].start]), sequences[wrongClassifications[i]].length);
//            kseq_t *seq = kseq_init(&buffer);
//            kseq_read(seq);
//            wr<<">"<<seq->name.s<<endl;
//            wr<<seq->seq.s<<endl;
//            kseq_destroy(seq);
//    }

//    wr.close();
    free(kmerBuffer.buffer);
    free(matchBuffer.buffer);
    munmap(queryFile.data, queryFile.fileSize + 1);
    munmap(targetDiffIdxList.data, targetDiffIdxList.fileSize + 1);
    munmap(targetInfoList.data, targetInfoList.fileSize + 1);
}

void Classifier::fillQueryKmerBufferParallel(QueryKmerBuffer & kmerBuffer, MmapedData<char> & seqFile,
                                             vector<Sequence> & seqs, bool * checker, size_t & processedSeqCnt,
                                             Query * queryList, const LocalParameters & par) {
    bool hasOverflow = false;
    omp_set_num_threads(*(int * )par.PARAM_THREADS.value);
#pragma omp parallel default(none), shared(checker, hasOverflow, processedSeqCnt, kmerBuffer, seqFile, seqs, cout, queryList)
    {
        SeqIterator seqIterator;
        size_t posToWrite;
#pragma omp for schedule(dynamic, 1)
        for (size_t i = 0; i < seqs.size(); i++) {
            if(checker[i] == false && !hasOverflow) {
                kseq_buffer_t buffer(const_cast<char *>(&seqFile.data[seqs[i].start]), seqs[i].length);
                kseq_t *seq = kseq_init(&buffer);
                kseq_read(seq);
                seqIterator.sixFrameTranslation(seq->seq.s);
                size_t kmerCnt = SeqIterator::kmerNumOfSixFrameTranslation(seq->seq.s);
                posToWrite = kmerBuffer.reserveMemory(kmerCnt);
                if (posToWrite + kmerCnt < kmerBuffer.bufferSize) {
                    seqIterator.fillQueryKmerBuffer(seq->seq.s, kmerBuffer, posToWrite, i);
                    checker[i] = true;
                    //seqs[i].length = strlen(seq->seq.s);
                    queryList[i].queryLength = getMaxCoveredLength((int) strlen(seq->seq.s));
                    queryList[i].queryId = i;
                    queryList[i].name = string(seq->name.s);
                            //buffer.entry.name.s;
#pragma omp atomic
                    processedSeqCnt ++;
                } else{
                    #pragma omp atomic
                    kmerBuffer.startIndexOfReserve -= kmerCnt;
                    hasOverflow = true;
                }
                kseq_destroy(seq);
            }

        }
    }
}


int Classifier::getMaxCoveredLength(int queryLength){
    if(queryLength % 3 == 2){
        return queryLength - 2; // 2
    } else if(queryLength % 3 == 1){
        return queryLength - 4; // 4
    } else{
        return queryLength - 3; // 3
    }
}

void Classifier::fillQueryKmerBufferParallel_paired(QueryKmerBuffer & kmerBuffer,
                                                    MmapedData<char> & seqFile1,
                                                    MmapedData<char> & seqFile2,
                                                    vector<Sequence> & seqs,
                                                    vector<Sequence> & seqs2,
                                                    bool * checker,
                                                    size_t & processedSeqCnt,
                                                    Query * queryList,
                                                    size_t numOfSeq,
                                                    const LocalParameters & par) {
    bool hasOverflow = false;

#pragma omp parallel default(none), shared(checker, hasOverflow, processedSeqCnt, kmerBuffer, seqFile1, seqFile2, seqs, seqs2, cout, queryList, numOfSeq)
    {
        SeqIterator seqIterator;
        SeqIterator seqIterator2;
        size_t posToWrite;
#pragma omp for schedule(dynamic, 1)
        for (size_t i = 0; i < numOfSeq; i++) {
            if(checker[i] == false && !hasOverflow) {
                // Read 1
                kseq_buffer_t buffer(const_cast<char *>(&seqFile1.data[seqs[i].start]), seqs[i].length);
                kseq_t *seq = kseq_init(&buffer);
                kseq_read(seq);
                size_t kmerCnt = SeqIterator::kmerNumOfSixFrameTranslation(seq->seq.s);

                // Read 2
                kseq_buffer_t buffer2(const_cast<char *>(&seqFile2.data[seqs[2].start]), seqs[2].length);
                kseq_t *seq2 = kseq_init(&buffer2);
                kseq_read(seq2);
                kmerCnt += SeqIterator::kmerNumOfSixFrameTranslation(seq2->seq.s);
                posToWrite = kmerBuffer.reserveMemory(kmerCnt);

                if (posToWrite + kmerCnt < kmerBuffer.bufferSize) {
                    checker[i] = true;
                    // Read 1
                    seqIterator.sixFrameTranslation(seq->seq.s);
                    seqIterator.fillQueryKmerBuffer(seq->seq.s, kmerBuffer, posToWrite, (int) i);

                    // Read 2
                    seqIterator2.sixFrameTranslation(seq2->seq.s);
                    seqIterator2.fillQueryKmerBuffer(seq2->seq.s, kmerBuffer, posToWrite, (int) i, strlen(seq->seq.s));

                    // Query Info
                    queryList[i].queryLength = getMaxCoveredLength((int) strlen(seq->seq.s)) +
                            getMaxCoveredLength((int) strlen(seq2->seq.s));
                    queryList[i].queryId = (int) i;
                    queryList[i].name = string(seq->name.s);
#pragma omp atomic
                    processedSeqCnt ++;
                } else{
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


void Classifier::linearSearchParallel(QueryKmer * queryKmerList, size_t & queryKmerCnt, const MmapedData<uint16_t> & targetDiffIdxList,
                                      const MmapedData<TargetKmerInfo> & targetInfoList, const MmapedData<DiffIdxSplit> & diffIdxSplits,
                                      Buffer<Match> & matchBuffer, const vector<int> & taxIdList, const vector<int> & spTaxIdList, const vector<TaxID> & genusTaxIdList,
                                      FILE * matchFile, const LocalParameters & par){
    cout<<"linearSearch start..."<<endl;
    ///Find the first index of garbage query k-mer (UINT64_MAX) and discard from there
    for(size_t checkN = queryKmerCnt - 1; checkN > 0; checkN--){
        if(queryKmerList[checkN].ADkmer != UINT64_MAX){
            queryKmerCnt = checkN + 1;
            break;
        }
    }

    ///Filter out meaningless target querySplits
    size_t numOfDiffIdxSplits = diffIdxSplits.fileSize / sizeof(DiffIdxSplit);
    size_t numOfDiffIdxSplits_use = numOfDiffIdxSplits;
    for(size_t i = 1; i < numOfDiffIdxSplits; i++){
        if(diffIdxSplits.data[i].ADkmer == 0 || diffIdxSplits.data[i].ADkmer == UINT64_MAX){
            diffIdxSplits.data[i] = {UINT64_MAX, UINT64_MAX, UINT64_MAX};
            numOfDiffIdxSplits_use--;
        }
    }

    cout<<"Filtering out meaningless target splits ... done"<<endl;

    //Devide query k-mer list into blocks for multi threading.
    //Each split has start and end points of query list + proper offset point of target k-mer list

    vector<QueryKmerSplit> querySplits;
    int threadNum = par.threads;

    //int * threadNum = (int *) par.PARAM_THREADS.value;
    uint64_t queryAA;
    if(threadNum == 1){ //Single thread
        querySplits.emplace_back(0, queryKmerCnt - 1, queryKmerCnt, diffIdxSplits.data[0]);
    } else if(threadNum == 2){ //Two threads
        size_t splitWidth = queryKmerCnt / 2;
        querySplits.emplace_back(0, splitWidth - 1, splitWidth, diffIdxSplits.data[0]);
        for(size_t tSplitCnt = 0; tSplitCnt < numOfDiffIdxSplits_use; tSplitCnt++){
            queryAA = AminoAcid(queryKmerList[splitWidth].ADkmer);
            if(queryAA <= AminoAcid(diffIdxSplits.data[tSplitCnt].ADkmer)){
                tSplitCnt = tSplitCnt - (tSplitCnt != 0);
                querySplits.emplace_back(splitWidth, queryKmerCnt - 1, queryKmerCnt - splitWidth, diffIdxSplits.data[tSplitCnt]);
                break;
            }
        }
    } else{ //More than two threads
        size_t splitWidth = queryKmerCnt / (threadNum - 1);
        querySplits.emplace_back(0, splitWidth - 1, splitWidth, diffIdxSplits.data[0]);
        for(int i = 1; i < threadNum; i ++){
            queryAA = AminoAcid(queryKmerList[splitWidth * i].ADkmer);
            bool needLastTargetBlock = true;
            for(size_t j = 0; j < numOfDiffIdxSplits_use; j++){
               if(queryAA <= AminoAcid(diffIdxSplits.data[j].ADkmer)){
                   j = j - (j!=0);
                   if(i != threadNum - 1)
                       querySplits.emplace_back(splitWidth * i, splitWidth * (i + 1) - 1, splitWidth, diffIdxSplits.data[j]);
                   else {
                       querySplits.emplace_back(splitWidth * i, queryKmerCnt - 1, queryKmerCnt - splitWidth * i,
                                                diffIdxSplits.data[j]);
                   }
                   needLastTargetBlock = false;
                   break;
               }
            }
            if(needLastTargetBlock){
                cout<<"needLastTargetBlock"<<endl;
                if(i != threadNum - 1)
                    querySplits.emplace_back(splitWidth * i, splitWidth * (i + 1) - 1, splitWidth, diffIdxSplits.data[numOfDiffIdxSplits_use - 1]);
                else {
                    querySplits.emplace_back(splitWidth * i, queryKmerCnt - 1, queryKmerCnt - splitWidth * i,
                                             diffIdxSplits.data[numOfDiffIdxSplits_use-1]);
                }
            }
        }
    }

    bool * splitCheckList = (bool *)malloc(sizeof(bool) * threadNum);
    fill_n(splitCheckList, threadNum, false);
    int completedSplitCnt = 0;
    cout<<"Deviding query k-mer list into blocks for multi threading... done"<<endl;

    size_t numOfTargetKmer = targetInfoList.fileSize / sizeof(TargetKmerInfo);
    size_t numOfDiffIdx = targetDiffIdxList.fileSize / sizeof(uint16_t);
    cout<<"The number of target k-mers: "<<numOfTargetKmer<<endl;

    cout<<"Hi"<<querySplits.size()<<endl;

    while(completedSplitCnt < threadNum) {
        bool hasOverflow = false;
#pragma omp parallel default(none), shared(numOfDiffIdx, completedSplitCnt, splitCheckList, numOfTargetKmer, hasOverflow, querySplits, queryKmerList, targetDiffIdxList, targetInfoList, matchBuffer, cout, genusTaxIdList, taxIdList, spTaxIdList)
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
#pragma omp for schedule(dynamic, 1)
            for (size_t i = 0; i < querySplits.size(); i++){
                if(hasOverflow || splitCheckList[i]) {
                    continue;
                }
                targetInfoIdx = querySplits[i].diffIdxSplit.infoIdxOffset - (i != 0);
                diffIdxPos = querySplits[i].diffIdxSplit.diffIdxOffset;
                currentTargetKmer = querySplits[i].diffIdxSplit.ADkmer;
                if(i == 0){
                    currentTargetKmer = getNextTargetKmer(currentTargetKmer, targetDiffIdxList.data, diffIdxPos);
                }
                currentQuery = UINT64_MAX;
                currentQueryAA = UINT64_MAX;

                for(size_t j = querySplits[i].start; j < querySplits[i].end + 1; j ++){
                    querySplits[i].start++;
                    ///Reuse the comparison data if queries are exactly identical
                    if(currentQuery == queryKmerList[j].ADkmer){
                        posToWrite = matchBuffer.reserveMemory(selectedMatches.size());
                        if(posToWrite + selectedMatches.size() >= matchBuffer.bufferSize){
                            hasOverflow = true;
                            querySplits[i].start = j; ///TODO why??
#pragma omp atomic
                            matchBuffer.startIndexOfReserve -= selectedMatches.size();
                            break;
                        } else{
                            range = selectedMatches.size();
                            for (size_t k = 0; k < range; k++) {
                                if(targetInfoList.data[selectedMatches[k]].redundancy){
                                    matchBuffer.buffer[posToWrite] = {queryKmerList[j].info.sequenceID,
                                                                      spTaxIdList[targetInfoList.data[selectedMatches[k]].sequenceID],
                                                                      spTaxIdList[targetInfoList.data[selectedMatches[k]].sequenceID],
                                                                      genusTaxIdList[targetInfoList.data[selectedMatches[k]].sequenceID],
                                                                      queryKmerList[j].info.pos, queryKmerList[j].info.frame,
                                                                      selectedHammingSum[k], 1, selectedHammings[k]};
                                } else{
                                    matchBuffer.buffer[posToWrite] = {queryKmerList[j].info.sequenceID,
                                                                      taxIdList[targetInfoList.data[selectedMatches[k]].sequenceID],
                                                                      spTaxIdList[targetInfoList.data[selectedMatches[k]].sequenceID],
                                                                      genusTaxIdList[targetInfoList.data[selectedMatches[k]].sequenceID],
                                                                      queryKmerList[j].info.pos, queryKmerList[j].info.frame,
                                                                      selectedHammingSum[k],0, selectedHammings[k]};
                                }
                                posToWrite ++;
                            }
                        }
                        continue;
                    }
                    selectedMatches.clear();
                    selectedHammingSum.clear();
                    selectedHammings.clear();

                    ///Reuse the candidate target k-mers to compare in DNA level if queries are the same at amino acid level but not at DNA level
                    if(currentQueryAA == AminoAcid(queryKmerList[j].ADkmer)){
                        compareDna(queryKmerList[j].ADkmer, candidateTargetKmers, startIdxOfAAmatch, selectedMatches, selectedHammingSum, selectedHammings);
                        posToWrite = matchBuffer.reserveMemory(selectedMatches.size());
                        if(posToWrite + selectedMatches.size() >= matchBuffer.bufferSize){
                            hasOverflow = true;
                            querySplits[i].start = j;
#pragma omp atomic
                            matchBuffer.startIndexOfReserve -= selectedMatches.size();
                            break;
                        } else{
                            range = selectedMatches.size();
                            for (size_t k = 0; k < range; k++) {
                                if(targetInfoList.data[selectedMatches[k]].redundancy){
                                    matchBuffer.buffer[posToWrite] = {queryKmerList[j].info.sequenceID,
                                                                      spTaxIdList[targetInfoList.data[selectedMatches[k]].sequenceID],
                                                                      spTaxIdList[targetInfoList.data[selectedMatches[k]].sequenceID],
                                                                      genusTaxIdList[targetInfoList.data[selectedMatches[k]].sequenceID],
                                                                      queryKmerList[j].info.pos, queryKmerList[j].info.frame,
                                                                      selectedHammingSum[k], 1, selectedHammings[k]};
                                } else{
                                    matchBuffer.buffer[posToWrite] = {queryKmerList[j].info.sequenceID,
                                                                      taxIdList[targetInfoList.data[selectedMatches[k]].sequenceID],
                                                                      spTaxIdList[targetInfoList.data[selectedMatches[k]].sequenceID],
                                                                      genusTaxIdList[targetInfoList.data[selectedMatches[k]].sequenceID],
                                                                      queryKmerList[j].info.pos, queryKmerList[j].info.frame,
                                                                      selectedHammingSum[k], 0, selectedHammings[k]};
                                }
                                posToWrite ++;
                            }
                        }
                        continue;
                    }
                    candidateTargetKmers.clear();

                    // Get next query, and start to find
                    currentQuery = queryKmerList[j].ADkmer;
                    currentQueryAA = AminoAcid(currentQuery);

                    // Skip target k-mers that are not matched in amino acid level
                    while (AminoAcid(currentQuery) > AminoAcid(currentTargetKmer) && (targetInfoIdx < numOfTargetKmer) && (diffIdxPos != numOfDiffIdx)) {
                        currentTargetKmer = getNextTargetKmer(currentTargetKmer, targetDiffIdxList.data, diffIdxPos);
                        targetInfoIdx++;
                    }

                    if(AminoAcid(currentQuery) != AminoAcid(currentTargetKmer)) ///Move to next query k-mer if there isn't any match.
                        continue;
                    else
                        startIdxOfAAmatch = targetInfoIdx;

                    ///Load target k-mers that are matched in amino acid level
                    while (AminoAcid(currentQuery) == AminoAcid(currentTargetKmer) && (targetInfoIdx < numOfTargetKmer) && (diffIdxPos != numOfDiffIdx)) {
                        candidateTargetKmers.push_back(currentTargetKmer);
                        currentTargetKmer = getNextTargetKmer(currentTargetKmer, targetDiffIdxList.data, diffIdxPos);
                        targetInfoIdx++;
                    }

                    ///Compare the current query and the loaded target k-mers and select
                    compareDna(currentQuery, candidateTargetKmers, startIdxOfAAmatch, selectedMatches, selectedHammingSum, selectedHammings);
                    posToWrite = matchBuffer.reserveMemory(selectedMatches.size());
                    if(posToWrite + selectedMatches.size() >= matchBuffer.bufferSize){
                        hasOverflow = true;
                        querySplits[i].start = j;
#pragma omp atomic
                        matchBuffer.startIndexOfReserve -= selectedMatches.size();
                        break;
                    } else{
                        range = selectedMatches.size();
                        for (size_t k = 0; k < range; k++) {
                            if(targetInfoList.data[selectedMatches[k]].redundancy){
                                matchBuffer.buffer[posToWrite] = {queryKmerList[j].info.sequenceID,
                                                                  spTaxIdList[targetInfoList.data[selectedMatches[k]].sequenceID],
                                                                  spTaxIdList[targetInfoList.data[selectedMatches[k]].sequenceID],
                                                                  genusTaxIdList[targetInfoList.data[selectedMatches[k]].sequenceID],
                                                                  queryKmerList[j].info.pos, queryKmerList[j].info.frame,
                                                                  selectedHammingSum[k], 1, selectedHammings[k]};
                            } else{
                                matchBuffer.buffer[posToWrite] = {queryKmerList[j].info.sequenceID,
                                                                  taxIdList[targetInfoList.data[selectedMatches[k]].sequenceID],
                                                                  spTaxIdList[targetInfoList.data[selectedMatches[k]].sequenceID],
                                                                  genusTaxIdList[targetInfoList.data[selectedMatches[k]].sequenceID],
                                                                  queryKmerList[j].info.pos, queryKmerList[j].info.frame,
                                                                  selectedHammingSum[k], 0, selectedHammings[k]};
                            }
                            posToWrite ++;
                        }
                    }
                }

               ///Check whether current split is completed or not
                if(querySplits[i].start - 1 == querySplits[i].end){
                    splitCheckList[i] = true;
                    #pragma omp atomic
                    completedSplitCnt ++;
                }
            }
        }

        if(hasOverflow)
            writeMatches(matchBuffer, matchFile);
    }
    free(splitCheckList);
    queryKmerCnt = 0;
    cout<<"end of linear seach parallel"<<endl;
}


void Classifier::writeMatches(Buffer<Match> & matchBuffer, FILE * matchFile){
    //SORT_PARALLEL(matchBuffer.buffer, matchBuffer.buffer + matchBuffer.startIndexOfReserve, Classifier::compareForWritingMatches);
    fwrite(matchBuffer.buffer, sizeof(Match), matchBuffer.startIndexOfReserve, matchFile);
    matchBuffer.startIndexOfReserve = 0;
}

bool Classifier::compareForWritingMatches(const Match & a, const Match & b){
    if (a.queryId < b.queryId) return true;
    else if (a.queryId == b.queryId) {
        if (a.frame < b.frame) return true;
        else if (a.frame == b.frame) {
            if (a.position < b.position) return true;
        }
    }
    return false;
}

// It compares query k-mers to target k-mers.
// If a query has matches, the matches with the smallest hamming distance will be selected
void Classifier::compareDna(uint64_t & query, vector<uint64_t> & targetKmersToCompare, const size_t & startIdx,
                            vector<size_t> & selectedMatches, vector<uint8_t> & selectedHammingSum,
                            vector<uint16_t> & selectedHammings) {

    size_t size = targetKmersToCompare.size();
    auto * hammingSums = new uint8_t[size + 1];
    uint8_t currentHammingSum;
    uint8_t minHammingSum = UINT8_MAX;

    // Calculate hamming distance
    for(size_t i = 0; i < size; i++){
        currentHammingSum = getHammingDistanceSum(query, targetKmersToCompare[i]);
        if(currentHammingSum < minHammingSum)
            minHammingSum = currentHammingSum;
        hammingSums[i] = currentHammingSum;
    }

    // Select target k-mers that passed hamming criteria
    for(size_t h = 0; h < size; h++){
        if(hammingSums[h] == minHammingSum){
            selectedMatches.push_back(startIdx + h);
            selectedHammingSum.push_back(hammingSums[h]);
            selectedHammings.push_back(getHammings(query, targetKmersToCompare[h]));
        }
    }

    delete[] hammingSums;
}

// It analyses the result of linear search.
void Classifier::analyseResultParallel(NcbiTaxonomy & ncbiTaxonomy, Buffer<Match> & matchBuffer,
                                       int seqNum, Query * queryList, const LocalParameters & par) {
    // Mmap the file of matches
//    struct MmapedData<Match> matchList = mmapData<Match>(matchFileName);
//    size_t numOfMatches = matchList.fileSize / sizeof(Match);


    size_t numOfMatches = matchBuffer.startIndexOfReserve;

    // Sort matches in order to analyze
    //SORT_PARALLEL(matchList.data, matchList.data + numOfMatches, Classifier::sortByGenusAndSpecies2);

    // Devide matches into blocks for multi threading
    MatchBlock *matchBlocks = new MatchBlock[seqNum];
    size_t matchIdx = 0;
    size_t blockIdx = 0;
    uint32_t currentQuery;
    while (matchIdx < numOfMatches) {
        currentQuery = matchBuffer.buffer[matchIdx].queryId;
        matchBlocks[blockIdx].id = currentQuery;
        matchBlocks[blockIdx].start = matchIdx;
        while ((currentQuery == matchBuffer.buffer[matchIdx].queryId) && (matchIdx < numOfMatches)) ++matchIdx;
        matchBlocks[blockIdx].end = matchIdx - 1;
        blockIdx++;
    }


    if (PRINT) {
        omp_set_num_threads(1);
    } else {
        omp_set_num_threads(par.threads);
    }

    // Process each block
#pragma omp parallel default(none), shared(cout, matchBlocks, matchBuffer, seqNum, ncbiTaxonomy, queryList, blockIdx, par)
    {
#pragma omp for schedule(dynamic, 1)
        for (size_t i = 0; i < blockIdx; ++i) {
           // chooseBestTaxon(ncbiTaxonomy, seqSegments[matchBlocks[i].id].length, matchBlocks[i].id,
            //                matchBlocks[i].start, matchBlocks[i].end, matchBuffer.buffer, queryList, par);

            chooseBestTaxon(ncbiTaxonomy,
                            matchBlocks[i].id,
                            matchBlocks[i].start,
                            matchBlocks[i].end,
                            matchBuffer.buffer,
                            queryList,
                            par);
        }
    }

    for (int i = 0; i < seqNum; i++) {
        ++taxCounts[queryList[i].classification];
    }
    delete[] matchBlocks;
    //munmap(matchList.data, matchList.fileSize + 1);
    cout << "End of analyseResultParallel" << endl;
}

TaxID Classifier::chooseBestTaxon(NcbiTaxonomy &ncbiTaxonomy, uint32_t currentQuery,
                                  size_t offset, size_t end, Match * matchList, Query * queryList,
                                  const LocalParameters & par ) {
    int queryLength = queryList[currentQuery].queryLength;
    TaxID selectedTaxon;
    if(PRINT) {
        cout<<"# "<<currentQuery<<endl;
        for (size_t i = offset; i < end + 1; i++) {
            cout << matchList[i].genusTaxID<<" "<<matchList[i].speciesTaxID << " " << matchList[i].taxID << " " <<
            int(matchList[i].frame)<< " " << matchList[i].position  << " " << int(matchList[i].hamming) << endl;
        }
    }

    // Get the best genus for current query
    vector<Match> matchesForLCA;
    matchesForLCA.reserve(end-offset+1);
    float maxScore;
    int res = getMatchesOfTheBestGenus(matchesForLCA, matchList, end, offset, queryLength, maxScore);

    if(PRINT) {
        cout<<"# "<<currentQuery<<" filtered"<<endl;
        for (size_t i = 0; i < matchesForLCA.size(); i++) {
            cout << matchesForLCA[i].genusTaxID<<" "<<matchesForLCA[i].speciesTaxID << " " << matchesForLCA[i].taxID << " " <<
                 int(matchesForLCA[i].frame)<< " " << matchesForLCA[i].position  << " " << int(matchesForLCA[i].hamming) << endl;
        }
    }

    if(PRINT){
        sort(matchesForLCA.begin(), matchesForLCA.end(), Classifier::sortByGenusAndSpecies);
    }

    //If there is no proper genus for current query, it is un-classified.
    if(res == 3){
        queryList[currentQuery].isClassified = false;
        queryList[currentQuery].classification = 0;
        queryList[currentQuery].score = 0;
        queryList[currentQuery].newSpecies = false;
        return 0;
    }

    float hammingSum;
    for(size_t i = 0; i < matchesForLCA.size(); i++ ){
        queryList[currentQuery].taxCnt[matchesForLCA[i].taxID] ++;
        hammingSum += matchesForLCA[i].hamming;
    }

    // If there are two or more good genus level candidates, find the LCA.
    if(res == 2){
        vector<TaxID> taxIdList;
        taxIdList.reserve(matchesForLCA.size());
        for(size_t i = 0; i < matchesForLCA.size(); i++ ){
            taxIdList.push_back(matchesForLCA[i].taxID);
        }
        selectedTaxon = ncbiTaxonomy.LCA(taxIdList)->taxId;
        queryList[currentQuery].isClassified = true;
        queryList[currentQuery].classification = selectedTaxon;
        queryList[currentQuery].score = maxScore;
        if(PRINT) {
            cout << "# " << currentQuery << " " << res << endl;
            for (size_t i = 0; i < taxIdList.size(); i++) {
                cout << i << " " << int(matchesForLCA[i].frame) << " " << matchesForLCA[i].position<< " " <<
                matchesForLCA[i].taxID << " " << int(matchesForLCA[i].hamming) <<" "<< matchesForLCA[i].red << endl;
            }
            cout << "Score: " << maxScore << " " << selectedTaxon << " "
                 << ncbiTaxonomy.taxonNode(selectedTaxon)->rank << endl;
        }
        return selectedTaxon;
    }

    queryList[currentQuery].score = maxScore;

    // Classify prokaryotes in genus level for highly diverged queries, not for virus
    if(maxScore < 0.8 && !ncbiTaxonomy.IsAncestor(par.virusTaxId, matchesForLCA[0].taxID)){
        selectedTaxon = ncbiTaxonomy.getTaxIdAtRank(matchesForLCA[0].taxID, "genus");
        queryList[currentQuery].isClassified = true;
        queryList[currentQuery].classification = selectedTaxon;
        queryList[currentQuery].newSpecies = true;
        if(PRINT) {
            cout << "# " << currentQuery << "HH" << endl;
            for (size_t i = 0; i < matchesForLCA.size(); i++) {
                cout << i << " " << int(matchesForLCA[i].frame) << " " << matchesForLCA[i].position<< " " <<
                     matchesForLCA[i].taxID << " " << int(matchesForLCA[i].hamming) <<" "<< matchesForLCA[i].red << endl;
            }
            cout << "Score: " << maxScore << "  " << selectedTaxon << " "
                 << ncbiTaxonomy.taxonNode(selectedTaxon)->rank << endl;
        }
        return selectedTaxon;
    }

    // Classify in species or lower level for queries that have close matches in reference DB.
    TaxID selectedLCA = classifyFurther(matchesForLCA, ncbiTaxonomy, queryLength);

    ///TODO optimize strain specific classification criteria
    // Strain classification only for high coverage with LCA of species level
    int numOfstrains = 0;
    TaxID strainID = 0;
    int count = 1;
    if(NcbiTaxonomy::findRankIndex(ncbiTaxonomy.taxonNode(selectedLCA)->rank) == 4){
        unordered_map<TaxID, int> strainMatchCnt;
        for(size_t i = 0; i < matchesForLCA.size(); i++){
            if(selectedLCA != matchesForLCA[i].taxID
                && ncbiTaxonomy.IsAncestor(selectedLCA, matchesForLCA[i].taxID)){
                strainMatchCnt[matchesForLCA[i].taxID]++;
            }
        }

        for(auto strainIt = strainMatchCnt.begin(); strainIt != strainMatchCnt.end(); strainIt ++){
            if(strainIt->second > 1){
                strainID = strainIt->first;
                numOfstrains++;
                count = strainIt->second;
            }
        }

        if(numOfstrains == 1 && count > 2) {
            selectedLCA = strainID;
        }
    }

    if(PRINT) {
        cout << "# " << currentQuery << endl;
        for (size_t i = 0; i < matchesForLCA.size(); i++) {
            cout << i << " " << int(matchesForLCA[i].frame) << " " << matchesForLCA[i].position<< " " <<
                 matchesForLCA[i].taxID << " " << int(matchesForLCA[i].hamming) <<" "<< matchesForLCA[i].red << endl;
        }
        cout << "Score: " << maxScore << "  " << selectedLCA << " " << ncbiTaxonomy.taxonNode(selectedLCA)->rank
             << endl;
    }

    if(PRINT && NcbiTaxonomy::findRankIndex(ncbiTaxonomy.taxonNode(selectedLCA)->rank) == 4){
        cout<<"sp\t"<<maxScore<<endl;
    } else if(PRINT && NcbiTaxonomy::findRankIndex(ncbiTaxonomy.taxonNode(selectedLCA)->rank) == 3){
        cout<<"sub\t"<<maxScore<<endl;
    } else if(PRINT && NcbiTaxonomy::findRankIndex(ncbiTaxonomy.taxonNode(selectedLCA)->rank) == 8){
        cout<<"genus\t"<<maxScore<<endl;
    }
    ///store classification results
    queryList[currentQuery].isClassified = true;
    queryList[currentQuery].classification = selectedLCA;
    queryList[currentQuery].score = maxScore;
    queryList[currentQuery].newSpecies = false;
    return selectedLCA;
}


int Classifier::getMatchesOfTheBestGenus(vector<Match> & matchesForMajorityLCA, Match * matchList, size_t end,
                                         size_t offset, int queryLength, float & bestScore){
    int conCnt;
    uint32_t hammingSum;
    float hammingMean;
    TaxID currentGenus;
    TaxID currentSpecies;

//    int maxCoveredLength;
//    if(queryLength % 3 == 2){
//        maxCoveredLength = queryLength - 2; // 2
//    } else if(queryLength % 3 == 1){
//        maxCoveredLength = queryLength - 4; // 4
//    } else{
//        maxCoveredLength = queryLength - 3; // 3
//    }

    vector<Match> filteredMatches;
    vector<vector<Match>> matchesForEachGenus;
    vector<bool> conservedWithinGenus;
    vector<float> scoreOfEachGenus;
    size_t i = offset;
    size_t offsetIdx;
    bool newOffset;
    while(i < end + 1) {
        currentGenus = matchList[i].genusTaxID;
        // For current genus
        while (currentGenus == matchList[i].genusTaxID && (i < end + 1)) {
            currentSpecies = matchList[i].speciesTaxID;
            // For current species
            // Filter un-consecutive matches (probably random matches)
            while(currentSpecies == matchList[i].speciesTaxID && (i < end + 1)){
                offsetIdx = i;
                i++;
                newOffset = true;
                hammingSum = 0;
                conCnt = 0;
                hammingMean = matchList[i-1].hamming;
                while((i < end + 1) && currentSpecies == matchList[i].speciesTaxID &&
                      matchList[i].position <= matchList[i-1].position + 3){
                    if(newOffset && ((float)matchList[i].hamming) <= hammingMean + 3){
                        newOffset = false;
                        hammingSum = matchList[offsetIdx].hamming;
                        filteredMatches.push_back(matchList[offsetIdx]);
                        conCnt++;
                    }
                    if(!newOffset) {
                        filteredMatches.push_back(matchList[i]);
                        conCnt++;
                        hammingSum += matchList[i].hamming;
                        hammingMean = float(hammingSum) / float(filteredMatches.size());
                    }
                    i++;
                }
                //TODO Should I remove the offset match after checking if the hamming was two high?
            }
        }

        // Construct a match combination using filtered matches of current genus
        // so that it can best cover the query, and score the combination
        if(!filteredMatches.empty()) {
           constructMatchCombination(filteredMatches, matchesForEachGenus, scoreOfEachGenus, queryLength);
        }
        filteredMatches.clear();
    }

    // If there are no meaningful genus
    if(scoreOfEachGenus.empty()){
        bestScore = 0;
        return 3;
    }

    float maxScore = *max_element(scoreOfEachGenus.begin(), scoreOfEachGenus.end());
    if(maxScore < 0.25)
        return 3;
    vector<size_t> maxIdx;
    for(size_t g = 0; g < scoreOfEachGenus.size(); g++){
        if(scoreOfEachGenus[g] > maxScore * 0.9f){
            maxIdx.push_back(g);
        }
    }
    bestScore = maxScore;

    for(size_t g = 0; g < maxIdx.size(); g++){
        matchesForMajorityLCA.insert(matchesForMajorityLCA.end(), matchesForEachGenus[maxIdx[g]].begin(),
                                     matchesForEachGenus[maxIdx[g]].end());
    }

    if(maxIdx.size() > 1){
        return 2;
    }
    return 1;

    //Three cases
    //1. one genus
    //2. more than one genus
    //3. no genus
}


// Assumption : There will be no frame shift of the same gene between species
// What about intergenic region?
void Classifier::constructMatchCombination(vector<Match> & filteredMatches,
                                           vector<vector<Match>> & matchesForEachGenus,
                                           vector<float> & scoreOfEachGenus,
                                           int queryLength){
    // Sort
    sort(filteredMatches.begin(), filteredMatches.end(), Classifier::sortMatchesByPos);

    // Do not allow overlaps between the same species
    size_t numOfFitMat = filteredMatches.size();
    vector<Match> matches;
    Match overlappedMatch;
    int overlapCnt;
    matches.reserve(numOfFitMat);
    uint8_t minHamming = 0;
    bool isTheLastOverlapped = false;
    size_t i = 0;
    TaxID overlapID;
    int uniqueIdCnt;

    // Check overlaps between subspecies
    while(i + 1 < numOfFitMat){
        overlapCnt = 0;
        overlapID = 0;
        uniqueIdCnt = 0;
        while(filteredMatches[i].speciesTaxID == filteredMatches[i+1].speciesTaxID &&
                filteredMatches[i].position/3 == filteredMatches[i+1].position/3 && (i + 1 < numOfFitMat)){
            if(overlapCnt == 0){
                overlapCnt++;
                overlapID = filteredMatches[i].taxID;
                uniqueIdCnt ++;
                overlappedMatch = filteredMatches[i];
                minHamming = filteredMatches[i].hamming;
            } else if (filteredMatches[i].hamming == minHamming && (filteredMatches[i].taxID != overlapID)){
                overlapCnt++;
            }
            ++i;
        }
        if(overlapCnt){
            if(filteredMatches[i].hamming == minHamming && (filteredMatches[i].taxID != overlapID)) overlapCnt++;
            if(overlapCnt == 1){ // Overlapping with match of higher hamming distance -> ignore
                matches.push_back(overlappedMatch);
            } else { // Overlapping with match of the same hamming distance -> species
                overlappedMatch.taxID = overlappedMatch.speciesTaxID;
                matches.push_back(overlappedMatch);
            }
            isTheLastOverlapped = (i == numOfFitMat - 1);
        } else{
            matches.push_back(filteredMatches[i]);
        }
        i++;
    }
    if(!isTheLastOverlapped) {
        matches.push_back(filteredMatches[numOfFitMat - 1]);
    }

    // Calculate Hamming distance & covered length
    int coveredPosCnt = 0;
    uint16_t currHammings;
    int aminoAcidNum = (int)queryLength / 3;
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
    auto * hammingsAtEachPos = new signed char[aminoAcidNum + 1];
    memset(hammingsAtEachPos, -1, (aminoAcidNum + 1));
    while(f < matchNum){
        currPos = matches[f].position / 3;
        currHammings = matches[f].rightEndHamming;
        if(GET_2_BITS(currHammings) > hammingsAtEachPos[currPos]) hammingsAtEachPos[currPos] = GET_2_BITS(currHammings);
        if(GET_2_BITS(currHammings>>2) > hammingsAtEachPos[currPos + 1]) hammingsAtEachPos[currPos + 1] = GET_2_BITS(currHammings>>2);
        if(GET_2_BITS(currHammings>>4) > hammingsAtEachPos[currPos + 2]) hammingsAtEachPos[currPos + 2] = GET_2_BITS(currHammings>>4);
        if(GET_2_BITS(currHammings>>6) > hammingsAtEachPos[currPos + 3]) hammingsAtEachPos[currPos + 3] = GET_2_BITS(currHammings>>6);
        if(GET_2_BITS(currHammings>>8) > hammingsAtEachPos[currPos + 4]) hammingsAtEachPos[currPos + 4] = GET_2_BITS(currHammings>>8);
        if(GET_2_BITS(currHammings>>10) > hammingsAtEachPos[currPos + 5]) hammingsAtEachPos[currPos + 5] = GET_2_BITS(currHammings>>10);
        if(GET_2_BITS(currHammings>>12) > hammingsAtEachPos[currPos + 6]) hammingsAtEachPos[currPos + 6] = GET_2_BITS(currHammings>>12);
        if(GET_2_BITS(currHammings>>14) > hammingsAtEachPos[currPos + 7]) hammingsAtEachPos[currPos + 7] = GET_2_BITS(currHammings>>14);
        f++;
    }

    // Sum up hamming distances and count the number of position covered by the matches.
    float hammingSum = 0;
    for(int h = 0; h < aminoAcidNum; h++){
        if(hammingsAtEachPos[h] == 0) { // Add 0 for 0 hamming dist.
            coveredPosCnt++;
        }else if(hammingsAtEachPos[h] != -1){ // Add 1.5, 2, 2.5 for 1, 2, 3 hamming dist. respectively
            hammingSum += 1.0f + (0.5f * hammingsAtEachPos[h]);
            coveredPosCnt ++;
        }
    }
    delete[] hammingsAtEachPos;

    // Ignore too few matches
    if(coveredPosCnt < 2) return;

    // Score current genus
    int coveredLength = coveredPosCnt * 3;
    if (coveredLength > queryLength) coveredLength = queryLength;
    float score = ((float)coveredLength - hammingSum) / (float)queryLength;
    scoreOfEachGenus.push_back(score);
    matchesForEachGenus.push_back(matches);

    if(PRINT) {
        cout << filteredMatches[0].genusTaxID << " " << coveredLength << " " << hammingSum << " " <<((float)coveredLength - hammingSum) / (float)queryLength <<
        " "<<matches.size()
             << endl;
    }
}

void Classifier::constructMatchCombination2(vector<Match> & filteredMatches, int maxCoveredLength,
                                           vector<vector<Match>> & matchesForEachGenus,
                                           vector<float> & scoreOfEachGenus, size_t queryLength){
    // Sort
    sort(filteredMatches.begin(), filteredMatches.end(), Classifier::sortMatchesByPos);

    // Do not allow overlaps between the same species
    size_t numOfFitMat = filteredMatches.size();
    vector<Match> matches;
    Match overlappedMatch;
    int overlapCnt;
    matches.reserve(numOfFitMat);
    uint8_t minHamming = 0;
    bool isTheLastOverlapped = false;
    size_t i = 0;

    // Check overlaps between subspecies
    while(i + 1 < numOfFitMat){
        overlapCnt = 0;
        while(filteredMatches[i].speciesTaxID == filteredMatches[i+1].speciesTaxID &&
              filteredMatches[i].position/3 == filteredMatches[i+1].position/3 && (i + 1 < numOfFitMat)){
            if(overlapCnt == 0){
                overlapCnt++;
                overlappedMatch = filteredMatches[i];
                minHamming = filteredMatches[i].hamming;
            } else if (filteredMatches[i].hamming == minHamming){
                overlapCnt++;
            }
            ++i;
        }
        if(overlapCnt){
            if(filteredMatches[i].hamming == minHamming) overlapCnt++;
            if(overlapCnt == 1){ // Overlapping with match of higher hamming distance -> ignore
                matches.push_back(overlappedMatch);
            } else { // Overlapping with match of the same hamming distance -> species
                overlappedMatch.taxID = overlappedMatch.speciesTaxID;
                matches.push_back(overlappedMatch);
            }
            isTheLastOverlapped = (i == numOfFitMat - 1);
        } else{
            matches.push_back(filteredMatches[i]);
        }
        i++;
    }
    if(!isTheLastOverlapped) {
        matches.push_back(filteredMatches[numOfFitMat - 1]);
    }

    // Calculate Hamming distance & covered length
    int coveredPosCnt = 0;
    uint16_t currHammings;
    int aminoAcidNum = (int)queryLength / 3;
    int currPos;
    size_t matchNum = matches.size();
    size_t f = 0;

    // Get the sum of hamming distance at each position
    auto * hammingsSumAtEachPos = new int[aminoAcidNum + 1];
    memset(hammingsSumAtEachPos, 0, sizeof(int) * (aminoAcidNum + 1));
    auto * coveredCnt = new int[aminoAcidNum + 1];
    memset(coveredCnt, 0, sizeof(int) * (aminoAcidNum + 1));

    while(f < matchNum){
        currPos = matches[f].position / 3;
        currHammings = matches[f].rightEndHamming;
        hammingsSumAtEachPos[currPos    ] += GET_2_BITS(currHammings);
        hammingsSumAtEachPos[currPos + 1] += GET_2_BITS(currHammings>>2);
        hammingsSumAtEachPos[currPos + 2] += GET_2_BITS(currHammings>>4);
        hammingsSumAtEachPos[currPos + 3] += GET_2_BITS(currHammings>>6);
        hammingsSumAtEachPos[currPos + 4] += GET_2_BITS(currHammings>>8);
        hammingsSumAtEachPos[currPos + 5] += GET_2_BITS(currHammings>>10);
        hammingsSumAtEachPos[currPos + 6] += GET_2_BITS(currHammings>>12);
        hammingsSumAtEachPos[currPos + 7] += GET_2_BITS(currHammings>>14);
        coveredCnt[currPos] ++;
        coveredCnt[currPos + 1] ++;
        coveredCnt[currPos + 2] ++;
        coveredCnt[currPos + 3] ++;
        coveredCnt[currPos + 4] ++;
        coveredCnt[currPos + 5] ++;
        coveredCnt[currPos + 6] ++;
        coveredCnt[currPos + 7] ++;
        f++;
    }

    // Get the average of hamming distance at each position
    float hammingSum = 0.0f;
    for(int h = 0; h < aminoAcidNum; h++){
        if (coveredCnt[h] != 0){
            hammingSum += (float) hammingsSumAtEachPos[h] / (float) coveredCnt[h];
            coveredPosCnt ++;
        }
    }

    delete[] hammingsSumAtEachPos;
    delete[] coveredCnt;

    // Score current genus
    int coveredLength = coveredPosCnt * 3;
    if(coveredLength > maxCoveredLength) coveredLength = maxCoveredLength;
    float score = ((float)coveredLength - hammingSum) / (float)maxCoveredLength;
    if((float)coveredLength <= (float)queryLength * 0.2f || (matchNum < 4) || score < 0.2) return; // Ignore genus with low coverage
    scoreOfEachGenus.push_back(score);
    matchesForEachGenus.push_back(matches);

    if(PRINT) {
        cout << filteredMatches[0].genusTaxID << " " << coveredLength << " " << hammingSum << " " <<((float)coveredLength - hammingSum) / (float)maxCoveredLength <<
             " "<<matches.size()
             << endl;
    }
}


//TODO kmer count -> covered length
//TODO hamming
TaxID Classifier::classifyFurther(const vector<Match> & matches, NcbiTaxonomy & taxonomy, uint32_t queryLength) {

    std::unordered_map<TaxID, int> taxIdCounts;
    float majorityCutoff = 0.8;
    float coverageThreshold = 0.8;
    float maxKmerCnt = queryLength/3.0f - kmerLength;

    for(Match match : matches){
        taxIdCounts[match.taxID] += 1;
        taxIdCounts[match.speciesTaxID] += (match.speciesTaxID != match.taxID); // subspecies
    }

    float currentCoverage;
    float currnetPercentage;
    bool haveMetCovThr = false;
    bool haveMetMajorityThr = false;
    bool tied = false;
    size_t matchNum = matches.size();
    int maxCnt;
    TaxID bestOne;
    TaxID currRank;
    vector<TaxID> ties;
    int minRank = INT_MAX;
    float selectedPercent = 0;
    TaxID selectedTaxon;
    for(auto it = taxIdCounts.begin(); it != taxIdCounts.end(); it++){
        if(it->second >= maxKmerCnt) {
            it->second = (int)maxKmerCnt - 1;
        }
        currentCoverage = (float)it->second / maxKmerCnt;
        currnetPercentage = (float)it->second / matchNum;
        currRank = NcbiTaxonomy::findRankIndex(taxonomy.taxonNode(it->first)->rank);

        if(currentCoverage > coverageThreshold && currRank <= 4){
            if(!haveMetCovThr){
                haveMetCovThr = true;
                maxCnt = it->second;
                bestOne = it->first;
                ties.push_back(it->first);
            } else if(it->second == maxCnt){
                tied = true;
                ties.push_back(it->first);
            } else if(it->second > maxCnt){
                ties.clear();
                ties.push_back(it->first);
                tied = false;
                bestOne = it->first;
                maxCnt = it->second;
            }
        } else if (currnetPercentage >= majorityCutoff && (!haveMetCovThr)) {
            // TaxID currParentTaxId = node->parentTaxId;
            haveMetMajorityThr = true;
            if ((currRank < minRank) || ((currRank == minRank) && (currnetPercentage > selectedPercent))) {
                selectedTaxon = it->first;
                minRank = currRank;
                selectedPercent = currnetPercentage;
            }
        }
    }

    if (haveMetCovThr) {
        if(tied){
            return taxonomy.LCA(ties)->taxId;
        } else{
            return bestOne;
        }
    } else if (haveMetMajorityThr) {
        return selectedTaxon;
    }

    return matches[0].genusTaxID;
}

int Classifier::getNumOfSplits() const { return this->numOfSplit; }

bool Classifier::compareForLinearSearch(const QueryKmer & a, const QueryKmer & b){
    if(a.ADkmer < b.ADkmer){
        return true;
    } else if(a.ADkmer == b.ADkmer){
        return (a.info.sequenceID < b.info.sequenceID);
    }
    return false;
}

bool Classifier::sortByGenusAndSpecies(const Match & a, const Match & b) {
    if (a.queryId < b.queryId) return true;
    else if (a.queryId == b.queryId) {
        if(a.genusTaxID < b.genusTaxID) return true;
        else if(a.genusTaxID == b.genusTaxID) {
            if(a.speciesTaxID < b.speciesTaxID) return true;
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

bool Classifier::sortByGenusAndSpecies2(const Match & a, const Match & b) {
    if (a.queryId < b.queryId) return true;
    else if (a.queryId == b.queryId) {
        if(a.genusTaxID < b.genusTaxID) return true;
        else if(a.genusTaxID == b.genusTaxID) {
            if(a.speciesTaxID < b.speciesTaxID) return true;
            else if (a.speciesTaxID == b.speciesTaxID) {
                if (a.position < b.position) return true;
            }
        }
    }
    return false;
}

bool Classifier::sortMatchesByPos(const Match & a, const Match & b) {
    if (a.position/3 < b.position/3) return true;
    else if (a.position/3 == b.position/3) {
        if(a.speciesTaxID < b.speciesTaxID) return true;
        else if(a.speciesTaxID == b.speciesTaxID){
            return a.hamming < b.hamming;
        }
    }
    return false;
}

void Classifier::writeReadClassification(Query * queryList, int queryNum, ofstream & readClassificationFile){
    for(int i = 0; i < queryNum; i++){
        readClassificationFile <<queryList[i].isClassified << "\t" << queryList[i].name << "\t" << queryList[i].classification << "\t" << queryList[i].queryLength << "\t" << queryList[i].score << "\t";
        for(auto it = queryList[i].taxCnt.begin(); it != queryList[i].taxCnt.end(); ++it){
            readClassificationFile<<it->first<<":"<<it->second<<" ";
        }
        readClassificationFile<<"\n";
    }
}

void Classifier::writeReportFile(const string & reportFileName, NcbiTaxonomy & ncbiTaxonomy, int numOfQuery){
    unordered_map<TaxID, TaxonCounts> cladeCounts = ncbiTaxonomy.getCladeCounts(taxCounts);
    FILE * fp;
    fp = fopen(reportFileName.c_str(), "w");
    writeReport(fp, ncbiTaxonomy, cladeCounts, numOfQuery);
    fclose(fp);
}

void Classifier::writeReport(FILE * fp, const NcbiTaxonomy & ncbiTaxonomy, const unordered_map<TaxID, TaxonCounts> & cladeCounts, unsigned long totalReads, TaxID taxID, int depth){
    auto it = cladeCounts.find(taxID);
    unsigned int cladeCount = (it == cladeCounts.end()? 0 : it->second.cladeCount);
    unsigned int taxCount = (it == cladeCounts.end()? 0 : it->second.taxCount);
    if(taxID == 0) {
        if (cladeCount > 0) {
            fprintf(fp, "%.2f\t%i\t%i\t0\tno rank\tunclassified\n", 100 * cladeCount / double(totalReads), cladeCount, taxCount);
        }
        writeReport(fp, ncbiTaxonomy, cladeCounts, totalReads, 1);
    } else{
        if (cladeCount == 0){
            return;
        }
        const TaxonNode * taxon = ncbiTaxonomy.taxonNode(taxID);
        fprintf(fp, "%.2f\t%i\t%i\t%i\t%s\t%s%s\n", 100 * cladeCount / double(totalReads), cladeCount, taxCount, taxID, taxon->rank.c_str(), string(2*depth, ' ').c_str(), taxon->name.c_str());
        vector<TaxID> children = it -> second.children;
        sort(children.begin(), children.end(), [&](int a, int b) { return cladeCountVal(cladeCounts, a) > cladeCountVal(cladeCounts,b); });
        for (TaxID childTaxId : children) {
            if (cladeCounts.count(childTaxId)) {writeReport(fp, ncbiTaxonomy, cladeCounts, totalReads, childTaxId, depth + 1);
            } else {
                break;
            }
        }
    }
}

unsigned int Classifier::cladeCountVal(const std::unordered_map<TaxID, TaxonCounts>& map, TaxID key) {
    typename std::unordered_map<TaxID, TaxonCounts>::const_iterator it = map.find(key);
    if (it == map.end()) {
        return 0;
    } else {
        return it->second.cladeCount;
    }
}

void Classifier::performanceTest(NcbiTaxonomy & ncbiTaxonomy, Query * queryList, int numOfquery, vector<int>& wrongs){

    ///Load the mapping file
    const char * mappingFile = "../../gtdb_taxdmp/assacc_to_taxid_gtdb.tsv";
    unordered_map<string, int> assacc2taxid;
    string key, value;
    ifstream map;
    map.open(mappingFile);
    if(map.is_open()){
        while(getline(map,key,'\t')){
            getline(map, value, '\n');
            assacc2taxid[key] = stoi(value);
        }
    } else{
        cout<<"Cannot open file for mappig from assemlby accession to tax ID"<<endl;
    }
    map.close();

    regex regex1("(GC[AF]_[0-9]*\\.[0-9]*)");
    smatch assacc;
    TaxID classificationResult;
    int rightAnswer;
    string queryName;


    counts.classificationCnt=0;
    counts.correct=0;
    counts.highRank=0;

    //number of targets at each rank
    counts.subspeciesTargetNumber=0;
    counts.speciesTargetNumber=0;
    counts.genusTargetNumber=0;
    counts.familyTargetNumber=0;
    counts.orderTargetNumber=0;
    counts.classTargetNumber=0;
    counts.phylumTargetNumber=0;
    counts.superkingdomTargetNumber=0;

    //number of classification at each rank
    counts.subspeciesCnt_try=0;
    counts.speciesCnt_try=0;
    counts.genusCnt_try=0;
    counts.familyCnt_try=0;
    counts.orderCnt_try=0;
    counts.classCnt_try=0;
    counts.phylumCnt_try=0;
    counts.superkingdomCnt_try=0;


    //number of correct classifications at each rank
    counts.subspeciesCnt_correct=0;
    counts.speciesCnt_correct=0;
    counts.genusCnt_correct=0;
    counts.familyCnt_correct=0;
    counts.orderCnt_correct=0;
    counts.classCnt_correct=0;
    counts.phylumCnt_correct=0;
    counts.superkingdomCnt_correct=0;

    for(int i = 0; i < numOfquery; i++) {
        classificationResult = queryList[i].classification;
        if (classificationResult == 0) {
            cout<<"here "<<i<<endl;
            continue;
        } else {
            classifiedCnt ++;
            queryName = queryList[i].name;
            regex_search(queryName, assacc, regex1);
            if (assacc2taxid.count(assacc[0].str())) {
                rightAnswer = assacc2taxid[assacc[0].str()];
            } else {
                cout << assacc[0].str() << " is not in the mapping file" << endl;
                continue;
            }
            //cout<<"compareTaxon"<<" "<<i<<endl;
            cout<<i<<" ";
            compareTaxon(classificationResult, rightAnswer, ncbiTaxonomy, wrongs, i);

        }
    }

    cout<<"Num of queries: " << queryInfos.size()  << endl;
    cout<<"Num of classifications: "<< counts.classificationCnt << endl;
    cout<<"Num of correct classifications: "<<counts.correct<<endl;
    cout<<"Num of correct but too broad classifications: "<<counts.highRank<<endl;
    cout<<"classified/total = " << float(counts.classificationCnt)/float(queryInfos.size()) << endl;
    cout<<"correct   /total = "<< float(counts.correct) / float(queryInfos.size())<<endl;
    cout<<"correct   /classifications = "<<float(counts.correct) / float(counts.classificationCnt) <<endl;
    cout<<"high rank /classifications = "<<float(counts.highRank) / float(counts.classificationCnt) <<endl << endl;

    cout<<"Number of targets at each rank / correct classification / tries"<<endl;
    cout<<"Superkingdom: " << counts.superkingdomTargetNumber << " / " << counts.superkingdomCnt_correct << " / "<<counts.superkingdomCnt_try<<endl;
    cout<<"Phylum      : " << counts.phylumTargetNumber << " / " << counts.phylumCnt_correct << " / "<<counts.phylumCnt_try<<endl;
    cout<<"Class       : " << counts.classTargetNumber << " / " << counts.classCnt_correct << " / "<<counts.classCnt_try<<endl;
    cout<<"Order       : " << counts.orderTargetNumber<<" / "<<counts.orderCnt_correct<<" / "<<counts.orderCnt_try<<endl;
    cout<<"Family      : " << counts.familyTargetNumber << " / " << counts.familyCnt_correct << " / "<<counts.familyCnt_try<<endl;
    cout<<"Genus       : " << counts.genusTargetNumber<<" / "<<counts.genusCnt_correct<<" / "<<counts.genusCnt_try<<endl;
    cout<<"Species     : " << counts.speciesTargetNumber<<" / "<<counts.speciesCnt_correct<<" / "<<counts.speciesCnt_try<<endl;
    cout<<"Subspecies  : " << counts.subspeciesTargetNumber<<" / "<<counts.subspeciesCnt_correct<<" / "<<counts.subspeciesCnt_try<<endl;
}

void Classifier::compareTaxon(TaxID shot, TaxID target, NcbiTaxonomy & ncbiTaxonomy, vector<int> & wrongs, int i) { ///target: subspecies or species
    const TaxonNode * shotNode = ncbiTaxonomy.taxonNode(shot);
    const TaxonNode * targetNode = ncbiTaxonomy.taxonNode(target);
    string shotRank = shotNode->rank;
    string targetRank = targetNode->rank;
    cout<<shot<<" "<<target<<" "<<shotRank<<" "<<targetRank<<" ";

    if(shot == 0){
        if(shotRank == "species" || shotRank == "subspecies")
            wrongs.push_back(i);
        cout<<i<<endl;
        cout<<"X"<<endl;
        return;
    } else{
        counts.classificationCnt++;
    }

    bool isCorrect = false;
    if(shot == target){
        counts.correct ++;
        isCorrect = true;
        cout<<"O"<<endl;
    } else if(NcbiTaxonomy::findRankIndex(shotRank) <= NcbiTaxonomy::findRankIndex(targetRank)){ //classified into wrong taxon or too specifically
        if(shotRank == "species" || shotRank == "subspecies")
            wrongs.push_back(i);
        cout<<"X"<<endl;
    } else { // classified at higher rank (too safe classification)
        if(shotRank == "superkingdom"){
            if(shotRank == "species" || shotRank == "subspecies")
                wrongs.push_back(i);
            cout<<"X"<<endl;
        } else if(shot == ncbiTaxonomy.getTaxIdAtRank(target, shotRank)){ //on right branch
            counts.correct ++;
            cout<<"O"<<endl;
            isCorrect = true;
        } else{ //on wrong branch
            if(shotRank == "species" || shotRank == "subspecies") wrongs.push_back(i);
            cout<<"X"<<endl;
        }
    }

    //count the number of classification at each rank
    if(shotRank == "subspecies") {
        counts.subspeciesCnt_try++;
    } else if(shotRank == "species") {
        counts.speciesCnt_try ++;
    } else if(shotRank == "genus"){
        counts.genusCnt_try ++;
    } else if(shotRank == "family"){
        counts.familyCnt_try++;
    } else if(shotRank == "order") {
        counts.orderCnt_try++;
    } else if(shotRank == "class") {
        counts.classCnt_try++;
    } else if(shotRank == "phylum") {
        counts.phylumCnt_try++;
    } else if(shotRank == "superkingdom"){
        counts.superkingdomCnt_try++;
    }

    if(!isCorrect) return;

    //count the number of correct classification at each rank
    if(shotRank == "subspecies"){
        counts.subspeciesCnt_correct ++;
    } else if(shotRank == "species") {
        counts.speciesCnt_correct ++;
    } else if(shotRank == "genus"){
        counts.genusCnt_correct ++;
    } else if(shotRank == "family"){
        counts.familyCnt_correct++;
    } else if(shotRank == "order") {
        counts.orderCnt_correct++;
    } else if(shotRank == "class") {
        counts.classCnt_correct++;
    } else if(shotRank == "phylum") {
        counts.phylumCnt_correct++;
    } else if(shotRank == "superkingdom"){
        counts.superkingdomCnt_correct++;
    }
}


