//
// Created by KJB on 01/09/2020.
//

#include "Classifier.h"
#include "LocalParameters.h"
#include <ctime>

Classifier::Classifier() {
    seqIterator = new SeqIterator();
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

Classifier::~Classifier() { delete seqIterator; }

void Classifier::startClassify(const char * queryFileName, const char * targetDiffIdxFileName, const char * targetInfoFileName,
                               const char * diffIdxSplitFileName, vector<int> & taxIdList, const LocalParameters & par,
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

    //Load query file & database
    struct MmapedData<char> queryFile = mmapData<char>(par.filenames[0].c_str());
    struct MmapedData<uint16_t> targetDiffIdxList = mmapData<uint16_t>(targetDiffIdxFileName);
    struct MmapedData<TargetKmerInfo> targetInfoList = mmapData<TargetKmerInfo>(targetInfoFileName);
    struct MmapedData<DiffIdxSplit> diffIdxSplits = mmapData<DiffIdxSplit>(diffIdxSplitFileName);

    //Query sequences
    vector<Sequence> sequences;
    IndexCreator::getSeqSegmentsWithHead(sequences, queryFile);
    size_t numOfSeq = sequences.size();
    Query * queryList = new Query[numOfSeq];

    //Checker for multi-threading
    bool * processedSeqChecker = new bool[numOfSeq];
    fill_n(processedSeqChecker, numOfSeq, false);

    //allocate memory for buffers
    QueryKmerBuffer kmerBuffer(kmerBufSize);
    Buffer<Match> matchBuffer(kmerBufSize);

    //Progress checker
    size_t processedSeqCnt = 0;
    size_t processedKmerCnt = 0;

    //Timer
    time_t beforeSearch, afterSearch, afterAnalyze;

    size_t numOfTatalQueryKmerCnt = 0;

    //extact k-mers from query sequences and compare them to target k-mer DB
    ///TODO measure time for extract & sort & search separately
    beforeSearch = time(NULL);
    while(processedSeqCnt < numOfSeq){
        fillQueryKmerBufferParallel(kmerBuffer, queryFile, sequences, processedSeqChecker, processedSeqCnt, queryList);
        numOfTatalQueryKmerCnt += kmerBuffer.startIndexOfReserve;
        omp_set_num_threads(ThreadNum);
        SORT_PARALLEL(kmerBuffer.buffer, kmerBuffer.buffer + kmerBuffer.startIndexOfReserve, Classifier::compareForLinearSearch);
        linearSearchParallel(kmerBuffer.buffer, kmerBuffer.startIndexOfReserve, targetDiffIdxList, targetInfoList, diffIdxSplits, matchBuffer, taxIdList, speciesTaxIdList, genusTaxIdList, matchFile);
    }
    cout<<"Number of query k-mers: "<<numOfTatalQueryKmerCnt<<endl;
    writeMatches(matchBuffer, matchFile);
    fclose(matchFile);
    afterSearch = time(NULL);
    cout<<"Time spent for searching: "<<double(afterSearch-beforeSearch)<<endl;

    //load matches and analyze
    cout<<"Analyse Result ... "<<endl;
    analyseResultParallel(taxonomy, sequences, matchFileName, numOfSeq, queryList);
    afterAnalyze = time(NULL);
    cout<<"Time spent for analyzing: "<<double(afterAnalyze-afterSearch)<<endl;

    //write report files
    ofstream readClassificationFile;
    readClassificationFile.open(par.filenames[0]+"_ReadClassification.tsv");
    writeReadClassification(queryList,numOfSeq,readClassificationFile);
    writeReportFile(par.filenames[0].c_str(), taxonomy, numOfSeq);

    //
    vector<int> wrongClassifications;
    sequences.clear();
    IndexCreator::getSeqSegmentsWithHead(sequences, queryFile);
    performanceTest(taxonomy, queryList, numOfSeq, wrongClassifications);
    ofstream wr;
    wr.open(par.filenames[0]+"_wrong");

    for (size_t i = 0; i < wrongClassifications.size(); i++) {
            kseq_buffer_t buffer(const_cast<char *>(&queryFile.data[sequences[wrongClassifications[i]].start]), sequences[wrongClassifications[i]].length);
            kseq_t *seq = kseq_init(&buffer);
            kseq_read(seq);
            wr<<">"<<seq->name.s<<endl;
            wr<<seq->seq.s<<endl;
            kseq_destroy(seq);
    }

    wr.close();
    free(kmerBuffer.buffer);
    free(matchBuffer.buffer);
    munmap(queryFile.data, queryFile.fileSize + 1);
    munmap(targetDiffIdxList.data, targetDiffIdxList.fileSize + 1);
    munmap(targetInfoList.data, targetInfoList.fileSize + 1);
}

void Classifier::fillQueryKmerBufferParallel(QueryKmerBuffer & kmerBuffer, MmapedData<char> & seqFile, vector<Sequence> & seqs, bool * checker, size_t & processedSeqCnt, Query * queryList) {
    bool hasOverflow = false;
    omp_set_num_threads(ThreadNum);
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
                    seqs[i].length = strlen(seq->seq.s);
                    queryList[i].queryLength = seqs[i].length;
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

void Classifier::linearSearchParallel(QueryKmer * queryKmerList, size_t & queryKmerCnt, const MmapedData<uint16_t> & targetDiffIdxList,
                                      const MmapedData<TargetKmerInfo> & targetInfoList, const MmapedData<DiffIdxSplit> & diffIdxSplits,
                                      Buffer<Match> & matchBuffer, const vector<int> & taxIdList, const vector<int> & taxIdListAtRank, const vector<TaxID> & genusTaxIdList,
                                      FILE * matchFile){
    cout<<"linearSearch start..."<<endl;
    ///Find the first index of garbage query k-mer (UINT64_MAX) and discard from there
    for(size_t checkN = queryKmerCnt - 1; checkN >= 0; checkN--){
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
    int threadNum = ThreadNum;
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

//    cout<<"Query"<<endl;
//    for(int i = 0 ; i < threadNum; i++){
//        cout << querySplits[i].diffIdxSplit.infoIdxOffset << " " << querySplits[i].diffIdxSplit.diffIdxOffset << endl;
//    }


    bool * splitCheckList = (bool *)malloc(sizeof(bool)*threadNum);
    fill_n(splitCheckList, threadNum, false);
    int completedSplitCnt = 0;
    cout<<"Deviding query k-mer list into blocks for multi threading... done"<<endl;

    ///taxonomical ID at the lowest rank? or at the rank of redundancy reduced
    vector<const vector<int> *> taxID;
    taxID.push_back(& taxIdList);
    taxID.push_back(& taxIdListAtRank);


    size_t numOfcall = 0;
    size_t queryIdx=0;
    size_t numOfTargetKmer = targetInfoList.fileSize / sizeof(TargetKmerInfo);
    size_t numOfDiffIdx = targetDiffIdxList.fileSize / sizeof(uint16_t);
    cout<<"The number of target k-mers: "<<numOfTargetKmer<<endl;
    size_t tempKmer = 0;
    size_t tempPos = 0;

//    for(size_t z = 0 ; z < numOfTargetKmer ; z++){
//        tempKmer = getNextTargetKmer(tempKmer, targetDiffIdxList.data, tempPos);
//        seqIterator->printKmerInDNAsequence(tempKmer);
//    }

    cout<<"Hi"<<querySplits.size()<<endl;
    omp_set_num_threads(ThreadNum);
    while(completedSplitCnt < threadNum) {
        bool hasOverflow = false;
#pragma omp parallel default(none), shared(numOfDiffIdx, queryIdx, completedSplitCnt, splitCheckList, numOfTargetKmer, hasOverflow, numOfcall, querySplits, queryKmerList, targetDiffIdxList, targetInfoList, matchBuffer, taxID, cout, genusTaxIdList)
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
                if(hasOverflow || splitCheckList[i])
                    continue;
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
                                                                      taxID[1]->at(targetInfoList.data[selectedMatches[k]].sequenceID),
                                                                      taxID[1]->at(targetInfoList.data[selectedMatches[k]].sequenceID),
                                                                      genusTaxIdList[targetInfoList.data[selectedMatches[k]].sequenceID],
                                                                      queryKmerList[j].info.pos, queryKmerList[j].info.frame,
                                                                      selectedHammingSum[k], 1, selectedHammings[k]};
                                } else{
                                    matchBuffer.buffer[posToWrite] = {queryKmerList[j].info.sequenceID,
                                                                      taxID[0]->at(targetInfoList.data[selectedMatches[k]].sequenceID),
                                                                      taxID[1]->at(targetInfoList.data[selectedMatches[k]].sequenceID),
                                                                      genusTaxIdList[targetInfoList.data[selectedMatches[k]].sequenceID],
                                                                      queryKmerList[j].info.pos, queryKmerList[j].info.frame, selectedHammingSum[k], 0,
                                                                      selectedHammings[k]};
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
                                    matchBuffer.buffer[posToWrite] = {queryKmerList[j].info.sequenceID, taxID[1]->at(targetInfoList.data[selectedMatches[k]].sequenceID),
                                                                      taxID[1]->at(targetInfoList.data[selectedMatches[k]].sequenceID), genusTaxIdList[targetInfoList.data[selectedMatches[k]].sequenceID],
                                                                      queryKmerList[j].info.pos, queryKmerList[j].info.frame, selectedHammingSum[k], 1, selectedHammings[k]};
                                } else{
                                    matchBuffer.buffer[posToWrite] = {queryKmerList[j].info.sequenceID, taxID[0]->at(targetInfoList.data[selectedMatches[k]].sequenceID), taxID[1]->at(targetInfoList.data[selectedMatches[k]].sequenceID),
                                                                      genusTaxIdList[targetInfoList.data[selectedMatches[k]].sequenceID], queryKmerList[j].info.pos, queryKmerList[j].info.frame, selectedHammingSum[k], 0, selectedHammings[k]};
                                }
                                posToWrite ++;
                            }
                        }
                        continue;
                    }
                    candidateTargetKmers.clear();

                    ///Get next query and start to find
                    currentQuery = queryKmerList[j].ADkmer;
                    currentQueryAA = AminoAcid(currentQuery);

                    ///Skip target k-mers that are not matched in amino acid level
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
                                matchBuffer.buffer[posToWrite] = {queryKmerList[j].info.sequenceID, taxID[1]->at(targetInfoList.data[selectedMatches[k]].sequenceID), taxID[1]->at(targetInfoList.data[selectedMatches[k]].sequenceID),
                                                                  genusTaxIdList[targetInfoList.data[selectedMatches[k]].sequenceID], queryKmerList[j].info.pos, queryKmerList[j].info.frame, selectedHammingSum[k], 1, selectedHammings[k]};
                            } else{
                                matchBuffer.buffer[posToWrite] = {queryKmerList[j].info.sequenceID, taxID[0]->at(targetInfoList.data[selectedMatches[k]].sequenceID), taxID[1]->at(targetInfoList.data[selectedMatches[k]].sequenceID),
                                                                  genusTaxIdList[targetInfoList.data[selectedMatches[k]].sequenceID], queryKmerList[j].info.pos, queryKmerList[j].info.frame, selectedHammingSum[k], 0, selectedHammings[k]};
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
    bool endCheck = false;
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



///It compares query k-mers to target k-mers. If a query has matches, the matches with the smallest difference are selected.
void Classifier::compareDna(uint64_t & query, vector<uint64_t> & targetKmersToCompare, const size_t & startIdx,
                            vector<size_t> & selectedMatches, vector<uint8_t> & selectedHammingSum,
                            vector<uint16_t> & selectedHammings) {

    size_t size = targetKmersToCompare.size();
    auto * hammingSums = new uint8_t[size + 1];
    uint8_t currentHammingSum;
    uint8_t minHammingSum = UINT8_MAX;
    ///Calculate hamming distance
    for(size_t i = 0; i < size; i++){
        currentHammingSum = getHammingDistanceSum(query, targetKmersToCompare[i]);
        if(currentHammingSum < minHammingSum)
            minHammingSum = currentHammingSum;
        hammingSums[i] = currentHammingSum;
    }

    ///Select target k-mers that passed hamming criteria
    for(size_t h = 0; h < size; h++){
        if(hammingSums[h] == minHammingSum){
            selectedMatches.push_back(startIdx + h);
            selectedHammingSum.push_back(hammingSums[h]);
            selectedHammings.push_back(getHammings(query, targetKmersToCompare[h]));
        }
    }

    delete[] hammingSums;
}

///It analyses the result of linear search.
void Classifier::analyseResultParallel(NcbiTaxonomy & ncbiTaxonomy, vector<Sequence> & seqSegments, char * matchFileName, int seqNum, Query * queryList) {
    //Mmap the file of matches
    struct MmapedData<Match> matchList = mmapData<Match>(matchFileName);
    cout << matchList.fileSize << "!!" << endl;
    cout << matchFileName << endl;
    size_t numOfMatches = matchList.fileSize / sizeof(Match);
    cout << "num of matches" << numOfMatches << endl;

    //Sort matches in order to analyze
    SORT_PARALLEL(matchList.data, matchList.data + numOfMatches, Classifier::sortByGenusAndSpecies2);
    //Devide matches into blocks for multi threading
    MatchBlock *matchBlocks = new MatchBlock[seqNum];
    cout << seqNum << endl;
    size_t matchIdx = 0;
    size_t blockIdx = 0;
    uint32_t currentQuery;
    while (matchIdx < numOfMatches) {
        currentQuery = matchList.data[matchIdx].queryId;
        matchBlocks[blockIdx].id = currentQuery;
        matchBlocks[blockIdx].start = matchIdx;
        while ((currentQuery == matchList.data[matchIdx].queryId) && (matchIdx < numOfMatches)) ++matchIdx;
        matchBlocks[blockIdx].end = matchIdx - 1;
        blockIdx++;
    }


    if (PRINT) {
        omp_set_num_threads(1);
    } else {
        omp_set_num_threads(ThreadNum);
    }

    //Process each block
#pragma omp parallel default(none), shared(cout, matchBlocks, matchList, seqSegments, seqNum, ncbiTaxonomy, queryList, blockIdx)
    {
#pragma omp for schedule(dynamic, 1)
        for (size_t i = 0; i < blockIdx; ++i) {
            TaxID selectedLCA = chooseBestTaxon(ncbiTaxonomy, seqSegments[matchBlocks[i].id].length, matchBlocks[i].id,
                                                 matchBlocks[i].start, matchBlocks[i].end, matchList.data, queryList);
        }
    }

    for (int i = 0; i < seqNum; i++) {
        ++taxCounts[queryList[i].classification];
    }
    delete[] matchBlocks;
    munmap(matchList.data, matchList.fileSize + 1);
    cout << "end of analyseResultParallel" << endl;
}

TaxID Classifier::chooseBestTaxon(NcbiTaxonomy &ncbiTaxonomy, const size_t &queryLength, const int &currentQuery,
                                  const size_t &offset, const size_t &end, Match *matchList, Query *queryList) {
    TaxID selectedTaxon;
    if(PRINT) {
        cout<<"# "<<currentQuery<<endl;
        for (int i = offset; i < end + 1; i++) {
            cout << matchList[i].genusTaxID<<" "<<matchList[i].speciesTaxID << " " << int(matchList[i].frame) << " " <<
            matchList[i].position << " " << matchList[i].taxID << " " << int(matchList[i].hamming) << endl;
        }
    }

    //get the best genus for current query
    vector<Match> matchesForLCA;
    matchesForLCA.reserve(end-offset+1);
    float maxScore;
    int res = getMatchesOfTheBestGenus(matchesForLCA, matchList, end, offset, queryLength, maxScore);

    if(PRINT){
        sort(matchesForLCA.begin(), matchesForLCA.end(), Classifier::sortByGenusAndSpecies);
    }
    float maxNum = queryLength / 3 - kmerLength + 1;
    //float normalizedScore = maxScore / maxNum;

    //If there is no proper genus for current query, it is un-classified.
    if(res == 3){
        queryList[currentQuery].isClassified = false;
        queryList[currentQuery].classification = 0;
        queryList[currentQuery].score = 0;
        queryList[currentQuery].newSpecies = false;
        return 0;
    }


    float hammingSum = 0.0f;
    for(size_t i = 0; i < matchesForLCA.size(); i++ ){
        queryList[currentQuery].taxCnt[matchesForLCA[i].taxID] ++;
        hammingSum += matchesForLCA[i].hamming;
    }

    //If there are two or more good genus level candidates, find the LCA.
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

    //Classify in genus level for highly diverged queries
    if(maxScore < 0.8){
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

    //Classify in species or lower level for queries that have close matches in reference DB.
    double selectedPercent = 0;
    //TaxID selectedLCA = match2LCA(matchesForLCA, ncbiTaxonomy, 0.8, selectedPercent, queryLength);

    TaxID selectedLCA = classifyFurther(matchesForLCA, ncbiTaxonomy, queryLength);

    ///TODO optimize strain specific classification criteria
    //Strain classification only for high coverage with LCA of species level
    if(NcbiTaxonomy::findRankIndex(ncbiTaxonomy.taxonNode(selectedLCA)->rank) == 4){ /// There are more strain level classifications with lower coverage threshold, but also with more false postives. 0.8~0.85 looks good.
        unordered_map<TaxID, int> strainMatchCnt;
        for(size_t i = 0; i < matchesForLCA.size(); i++ ){
            if(selectedLCA != matchesForLCA[i].taxID && ncbiTaxonomy.IsAncestor(selectedLCA, matchesForLCA[i].taxID)){
                strainMatchCnt[matchesForLCA[i].taxID]++;
            }
        }

        if(strainMatchCnt.size() == 1 && strainMatchCnt.begin()->second > 2) {
            selectedLCA = strainMatchCnt.begin()->first;
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

    int maxCoveredLength;
    if(queryLength % 3 == 2){
        maxCoveredLength = queryLength - 2; // 2
    } else if(queryLength % 3 == 1){
        maxCoveredLength = queryLength - 4; // 4
    } else{
        maxCoveredLength = queryLength - 3; // 3
    }

    vector<Match> filteredMatches;
    vector<vector<Match>> matchesForEachGenus;
    vector<bool> conservedWithinGenus;
    vector<float> scoreOfEachGenus;
    size_t i = offset;
    size_t offsetIdx;
    bool newOffset;
    while(i < end + 1) {
        currentGenus = matchList[i].genusTaxID;
        //For current genus
        while (currentGenus == matchList[i].genusTaxID && (i < end + 1)) {
            currentSpecies = matchList[i].speciesTaxID;
            //For current species
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

        // Construct a list of matches for scoring current genus
        if(!filteredMatches.empty()) {
           constructMatchCombination(filteredMatches, maxCoveredLength, matchesForEachGenus, scoreOfEachGenus,
                                       queryLength);
        }
        filteredMatches.clear();
    }

    // If there are no meaningful genus
    if(scoreOfEachGenus.empty()){
        bestScore = 0;
        return 3;
    }

    float maxScore = *max_element(scoreOfEachGenus.begin(), scoreOfEachGenus.end());
    vector<size_t> maxIdx;
    for(size_t g = 0; g < scoreOfEachGenus.size(); g++){
        if(scoreOfEachGenus[g] > maxScore * 0.95f){
            maxIdx.push_back(g);
        }
//        if(scoreOfEachGenus[g] > maxScore){
//            maxScore = scoreOfEachGenus[g];
//            maxIdx.clear();
//            maxIdx.push_back(g);
//        } else if(scoreOfEachGenus[g] == maxScore){
//
//        }
    }
    bestScore = maxScore;

//    // Choose the best genus
//    float maxScore = -100;
//    vector<size_t> maxIdx;
//    for(size_t g = 0; g < scoreOfEachGenus.size(); g++){
//        if(scoreOfEachGenus[g] > maxScore){
//            maxScore = scoreOfEachGenus[g];
//            maxIdx.clear();
//            maxIdx.push_back(g);
//        } else if(scoreOfEachGenus[g] == maxScore){
//            maxIdx.push_back(g);
//        }
//    }
//    bestScore = maxScore;

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


// TODO How about now allowing overlapping between species
// Assumption : There will be no frame shift of the same gene between species
// What about intergenic region?
// I think this function can be combined with the calling function
void Classifier::constructMatchCombination(vector<Match> & filteredMatches, int maxCoveredLength,
                                            vector<vector<Match>> & matchesForEachGenus,
                                            vector<float> & scoreOfEachGenus, size_t queryLength){
    //Maximum covered length

    vector<Match> overlaps;
    // Sort
    sort(filteredMatches.begin(), filteredMatches.end(), Classifier::sortMatchesByPos);

    // Do not allow overlaps between the same species
    size_t l = filteredMatches.size();
    vector<Match> matches;
    matches.reserve(l);
    bool overlapped;
    uint8_t minHamming = 0;
    bool isTheLastOverlapped = false;
    size_t i = 0;
    while(i + 1 < l){
        //check overlap
        overlapped = false;
        while(filteredMatches[i].speciesTaxID == filteredMatches[i+1].speciesTaxID &&
            filteredMatches[i].position/3 == filteredMatches[i+1].position/3 && (i + 1 < l)){
            if(!overlapped) {
                overlapped = true;
                overlaps.push_back(filteredMatches[i]);
                minHamming = filteredMatches[i].hamming;
            } else if(filteredMatches[i].hamming == minHamming){
                overlaps.push_back(filteredMatches[i]);
            }
            i++;
        }
        if(overlapped) {
            if(filteredMatches[i].hamming == minHamming) overlaps.push_back(filteredMatches[i]);
            if(overlaps.size() == 1){
                matches.push_back(overlaps[0]);
            } else {
                overlaps[0].taxID = overlaps[0].speciesTaxID;
                overlaps[0].red = 1;
                matches.push_back(overlaps[0]);
            }
            overlaps.clear();
            isTheLastOverlapped = (i == l - 1);
        } else{
            matches.push_back(filteredMatches[i]);
        }
        i++;
    }
    if(!isTheLastOverlapped) {
        matches.push_back(filteredMatches[l-1]);
    }

    // Hamming distance & covered length
    int coveredPosCnt = 0;
    uint16_t currHammings;
    int size = (int)queryLength/3;
    auto * hammingsAtEachPos = new signed char[size + 1];
    memset(hammingsAtEachPos, -1, (size + 1));
    int currPos;
    //TODO Using max hamming at each position -> random match가 문제..
    size_t matchNum = matches.size();
    size_t f = 0;
    while(f < matchNum){
        currPos = matches[f].position / 3;
        currHammings = matches[f].rightEndHamming;
        for(int i2 = 0; i2 < 8; i2++){
            if((signed char)GET_2_BITS(currHammings>>2*i2) > hammingsAtEachPos[currPos + i2]){
                hammingsAtEachPos[currPos + i2] = GET_2_BITS(currHammings>>2*i2);
            }
        }
        f++;
    }
    float hammingSum = 0;
    for(int h = 0; h < size; h++){
       // cout<<(int)hammingsAtEachPos[h]<<" ";
        if(hammingsAtEachPos[h] == 0) {
            //hammingSum += hammingsAtEachPos[h];
            coveredPosCnt ++;
        } else if(hammingsAtEachPos[h] == 1){
            hammingSum += 1.5f;
            coveredPosCnt ++;
        } else if(hammingsAtEachPos[h] == 2){
            hammingSum += 2.0f;
            coveredPosCnt ++;
        } else if(hammingsAtEachPos[h] == 3){
            hammingSum += 2.5f;
            coveredPosCnt ++;
        }
    }

    delete[] hammingsAtEachPos;

//    for(size_t m = 1; m < matches.size(); m++){
//        gap = matches[m].position - matches[m-1].position;
//        if(gap > 24){
//            coveredLength += 8;
//
//            hammingSum += matches[m].hamming;
//        }else{
//            curHammings = matches[m].rightEndHamming;
//            for(int i2 = 0; i2 < gap/3; i2++){
//                hammingSum += GET_2_BITS(curHammings);
//                curHammings = curHammings >> 0X2U;
//            }
//            coveredLength += gap;
//        }
//    }
    int coveredLength = coveredPosCnt * 3;
    if(coveredLength > maxCoveredLength) coveredLength = maxCoveredLength;

    if((float)coveredLength <= (float)queryLength * 0.2f) return;
    scoreOfEachGenus.push_back(((float)coveredLength - hammingSum) / (float)maxCoveredLength);
    matchesForEachGenus.push_back(matches);
    if(PRINT) {
        cout << filteredMatches[0].genusTaxID << " " << coveredPosCnt << " " << hammingSum << " " << matches.size()
             << endl;
    }
}

// TODO It can be silplified
TaxID Classifier::match2LCA(const std::vector<Match> & matchList, NcbiTaxonomy & taxonomy, float majorityCutoff,
                            double &selectedPercent, uint32_t queryLength) {

    std::map<TaxID, taxNode> ancTaxIdsCounts;

    selectedPercent = 0;
    double totalAssignedSeqsWeights = 0.0;

    for (size_t i = 0; i < matchList.size(); ++i) {
        TaxID currTaxId = matchList[i].taxID;
        double currWeight = 1;
        // ignore unassigned sequences
        if (currTaxId == 0) {
            continue;
        }
        TaxonNode const *node = taxonomy.taxonNode(currTaxId, false);
        if (node == NULL) {
            Debug(Debug::ERROR) << "taxonid: " << currTaxId << " does not match a legal taxonomy node.\n";
            EXIT(EXIT_FAILURE);
        }
        totalAssignedSeqsWeights += currWeight;

        // each start of a path due to an orf is a candidate
        if (ancTaxIdsCounts.find(currTaxId) != ancTaxIdsCounts.end()) { //원소가 있다면
            ancTaxIdsCounts[currTaxId].update(currWeight, 0);
        } else {
            taxNode currNode;
            currNode.set(currWeight, true, 0);
            ancTaxIdsCounts.insert(std::pair<TaxID, taxNode>(currTaxId, currNode));
        }

        // iterate all ancestors up to root (including). add currWeight and candidate status to each
        TaxID currParentTaxId = node->parentTaxId;
        while (currParentTaxId != currTaxId) {
            if (ancTaxIdsCounts.find(currParentTaxId) != ancTaxIdsCounts.end()) {
                ancTaxIdsCounts[currParentTaxId].update(currWeight, currTaxId);
            } else {
                taxNode currParentNode;
                currParentNode.set(currWeight, false, currTaxId);
                ancTaxIdsCounts.insert(std::pair<TaxID, taxNode>(currParentTaxId, currParentNode));
            }
            // move up:
            currTaxId = currParentTaxId;
            node = taxonomy.taxonNode(currParentTaxId, false);
            currParentTaxId = node->parentTaxId;
        }
    }

    // select the lowest ancestor that meets the cutoff
    int minRank = INT_MAX;
    TaxID selectedTaxon = 0;
    float coverageThreshold = 0.8;

    float curCoverage;
    int spFisrtMaxWeight = 0;
    TaxID first;
    int maximunPossibleKmerNum = queryLength / 3 - kmerLength;
    bool haveMetCovThr = false;
    bool tied = false;
    vector<TaxID> ties;


    for (auto it = ancTaxIdsCounts.begin(); it != ancTaxIdsCounts.end(); it++) {
        // consider only candidates
        if (!(it->second.isCandidate)) {
            continue;
        }
        double currPercent = float(it->second.weight) / totalAssignedSeqsWeights;
        if(it->second.weight >= maximunPossibleKmerNum) {
            it->second.weight = maximunPossibleKmerNum - 1;
        }
        curCoverage = float(it->second.weight) / float(maximunPossibleKmerNum);

        TaxID currTaxId = it->first;
        TaxonNode const *node = taxonomy.taxonNode(currTaxId, false);
        int currRankIdx = NcbiTaxonomy::findRankIndex(node->rank);

        if (curCoverage > coverageThreshold && currRankIdx <= 4) {
            if (!haveMetCovThr) {
                haveMetCovThr = true;
                spFisrtMaxWeight = it->second.weight;
                first = currTaxId;
                ties.push_back(currTaxId);
            } else if (it->second.weight == spFisrtMaxWeight) {
                tied = true;
                ties.push_back(currTaxId);
            } else if (it->second.weight > spFisrtMaxWeight) {
                ties.clear();
                ties.push_back(currTaxId);
                tied = false;
                first = currTaxId;
                spFisrtMaxWeight = it->second.weight;
            }
        }

        else if (currPercent >= majorityCutoff && (!haveMetCovThr)) {
            // TaxID currParentTaxId = node->parentTaxId;
            if ((currRankIdx < minRank) || ((currRankIdx == minRank) && (currPercent > selectedPercent))) {
                selectedTaxon = it->first;
                minRank = currRankIdx;
                selectedPercent = currPercent;
            }
        }
    }

    if (haveMetCovThr) {
        if(tied){
            return taxonomy.LCA(ties)->taxId;
        } else{
            return first;
        }
    } else {
        return selectedTaxon;
    }
}


//TODO kmer count -> covered length
//TODO hamming
TaxID Classifier::classifyFurther(const vector<Match> & matches, NcbiTaxonomy & taxonomy, uint32_t queryLength) {

    std::map<TaxID, int> taxIdCounts;
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
            it->second = maxKmerCnt-1;
        }
        currentCoverage = (float)it->second/maxKmerCnt;
        currnetPercentage = (float)it->second/matchNum;
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
                if (a.frame < b.frame) return true;
                else if (a.frame == b.frame) {
                    if (a.position < b.position) return true;
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

void Classifier::writeReportFile(const char * queryFileName, NcbiTaxonomy & ncbiTaxonomy, int numOfQuery){
    unordered_map<TaxID, TaxonCounts> cladeCounts = ncbiTaxonomy.getCladeCounts(taxCounts);
    string outputFile = string(queryFileName) + "_REPORTFILE.tsv";
    FILE * fp;
    fp = fopen(outputFile.c_str(), "w");
    writeReport(fp, ncbiTaxonomy, cladeCounts, numOfQuery);
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
        wrongs.push_back(i);
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
        wrongs.push_back(i);
        cout<<"X"<<endl;
    } else { // classified at higher rank (too safe classification)
        if(shotRank == "superkingdom"){
            wrongs.push_back(i);
            cout<<"X"<<endl;
        } else if(shot == ncbiTaxonomy.getTaxIdAtRank(target, shotRank)){ //on right branch
            counts.correct ++;
            cout<<"O"<<endl;
            isCorrect = true;
        } else{ //on wrong branch
            wrongs.push_back(i);
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


