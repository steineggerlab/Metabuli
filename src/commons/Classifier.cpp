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
    totalMatchCount = 0;
    multipleMatchCount = 0;

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

void Classifier::startClassify(const char * queryFileName, const char * targetDiffIdxFileName, const char * targetInfoFileName, const char * diffIdxSplitFileName, vector<int> & taxIdList, const LocalParameters & par) {
    string names, nodes, merged;
    if(par.gtdbOrNcbi == 1 || par.gtdbOrNcbi == 0 ){
        cout<<"Classifying query sequences based on taxonomy of GTDB"<<endl;
        names = "../../gtdb_taxdmp/names.dmp";
        nodes = "../../gtdb_taxdmp/nodes.dmp";
        merged = "../../gtdb_taxdmp/merged.dmp";
    } else if(par.gtdbOrNcbi == 2 ){
        cout<<"Classifying query sequences based on taxonomy of NCBI"<<endl;
        names = "../../ncbi_taxdmp/names.dmp";
        nodes = "../../ncbi_taxdmp/nodes.dmp";
        merged = "../../ncbi_taxdmp/merged.dmp";
    } else{
        cout<<"ERROR"<<endl;
        return;
    }

    //taxonomical ID
    NcbiTaxonomy ncbiTaxonomy(names, nodes, merged);

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
    ncbiTaxonomy.createTaxIdListAtRank(taxIdList, speciesTaxIdList, "species");
    ncbiTaxonomy.createTaxIdListAtRank(taxIdList, genusTaxIdList, "genus");
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

    //query & target
    struct MmapedData<char> queryFile = mmapData<char>(par.filenames[0].c_str());
    struct MmapedData<uint16_t> targetDiffIdxList = mmapData<uint16_t>(targetDiffIdxFileName);
    struct MmapedData<TargetKmerInfo> targetInfoList = mmapData<TargetKmerInfo>(targetInfoFileName);
    struct MmapedData<DiffIdxSplit> diffIdxSplits = mmapData<DiffIdxSplit>(diffIdxSplitFileName);

    int testNum2 = targetInfoList.fileSize/sizeof(TargetKmerInfo);
    for(int i = 0 ; i < 1000; i++){
        cout<<i<<" "<<targetInfoList.data[i].sequenceID<<endl;
    }
    //query sequences
    vector<Sequence> sequences;
    IndexCreator::getSeqSegmentsWithHead(sequences, queryFile);
    size_t numOfSeq = sequences.size();
    Query * queryList = new Query[numOfSeq];

    //check for multi-threading
    bool * processedSeqChecker = new bool[numOfSeq];
    fill_n(processedSeqChecker, numOfSeq, false);

    //allocate memory for buffers
    QueryKmerBuffer kmerBuffer(kmerBufSize);
    Buffer<Match> matchBuffer(kmerBufSize);

    size_t processedSeqCnt = 0;
    size_t processedKmerCnt = 0;

    time_t beforeSearch, afterSearch, afterAnalyze;

    size_t numOfTatalQueryKmerCnt = 0;

    //extact k-mers from query sequences and compare them to target k-mer DB
    cout<<"hi"<<endl;
    ///TODO measure time for extract & sort & search separately
    beforeSearch = time(NULL);
    while(processedSeqCnt < numOfSeq){
        fillQueryKmerBufferParallel(kmerBuffer, queryFile, sequences, processedSeqChecker, processedSeqCnt, queryList);
        numOfTatalQueryKmerCnt += kmerBuffer.startIndexOfReserve;
        cout<<"buffer overflowed"<<endl;
        omp_set_num_threads(ThreadNum);
        SORT_PARALLEL(kmerBuffer.buffer, kmerBuffer.buffer + kmerBuffer.startIndexOfReserve, Classifier::compareForLinearSearch);
        cout<<"buffer sorted"<<endl;
        linearSearchParallel(kmerBuffer.buffer, kmerBuffer.startIndexOfReserve, targetDiffIdxList, targetInfoList, diffIdxSplits, matchBuffer, taxIdList, speciesTaxIdList, genusTaxIdList, matchFile);
        cout<<"after linear search parallel"<<endl;
    }
    cout<<"total kmer count: "<<numOfTatalQueryKmerCnt<<endl;
    writeMatches(matchBuffer, matchFile);
    fclose(matchFile);
    afterSearch = time(NULL);
    cout<<"Time spent for searching: "<<double(afterSearch-beforeSearch)<<endl;

    //load matches and analyze
    cout<<"analyse Result"<<endl;
    analyseResultParallel(ncbiTaxonomy, sequences, matchFileName, numOfSeq, queryList);
    afterAnalyze = time(NULL);
    cout<<"Time spent for analyzing: "<<double(afterAnalyze-afterSearch)<<endl;

    //write report files
    ofstream readClassificationFile;
    readClassificationFile.open(par.filenames[0]+"_ReadClassification.tsv");
    writeReadClassification(queryList,numOfSeq,readClassificationFile);
    writeReportFile(par.filenames[0].c_str(), ncbiTaxonomy, numOfSeq);
    performanceTest(ncbiTaxonomy, queryList, numOfSeq);

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
                KSeqBuffer buffer(const_cast<char *>(&seqFile.data[seqs[i].start]), seqs[i].length);
                buffer.ReadEntry();
                seqIterator.sixFrameTranslation(buffer.entry.sequence.s);
                size_t kmerCnt = seqIterator.kmerNumOfSixFrameTranslation(buffer.entry.sequence.s);
                posToWrite = kmerBuffer.reserveMemory(kmerCnt);
                if (posToWrite + kmerCnt < kmerBuffer.bufferSize) {
                    seqIterator.fillQueryKmerBuffer(buffer.entry.sequence.s, kmerBuffer, posToWrite, i);
                    checker[i] = true;
                    seqs[i].length = strlen(buffer.entry.sequence.s);
                    queryList[i].queryLength = seqs[i].length;
                    queryList[i].queryId = i;
                    queryList[i].name = buffer.entry.name.s;
#pragma omp atomic
                    processedSeqCnt ++;
                } else{
                    #pragma omp atomic
                    kmerBuffer.startIndexOfReserve -= kmerCnt;
                    hasOverflow = true;
                }
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

//    for(int i = 0 ; i < numOfDiffIdxSplits; i++){
//        cout<<diffIdxSplits.data[i].infoIdxOffset<<" "<<diffIdxSplits.data[i].diffIdxOffset<<endl;
//    }

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
            if(queryAA < AminoAcid(diffIdxSplits.data[tSplitCnt].ADkmer)){
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

    omp_set_num_threads(threadNum);
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
            vector<uint8_t> selectedHammings;
            vector<size_t> selectedMatches;

            size_t startIdxOfAAmatch = 0;
            size_t posToWrite;
            size_t range;
#pragma omp for schedule(dynamic, 1)
            for (size_t i = 0; i < querySplits.size(); i++){
                if(hasOverflow || splitCheckList[i])
                    continue;
                targetInfoIdx = querySplits[i].diffIdxSplit.infoIdxOffset;
                diffIdxPos = querySplits[i].diffIdxSplit.diffIdxOffset + 1;
                currentTargetKmer = querySplits[i].diffIdxSplit.ADkmer;
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
                                    matchBuffer.buffer[posToWrite] = {queryKmerList[j].info.sequenceID, taxID[1]->at(targetInfoList.data[selectedMatches[k]].sequenceID),
                                                                      genusTaxIdList[targetInfoList.data[selectedMatches[k]].sequenceID], queryKmerList[j].info.pos, queryKmerList[j].info.frame, selectedHammings[k],1};
                                } else{
                                    matchBuffer.buffer[posToWrite] = {queryKmerList[j].info.sequenceID, taxID[0]->at(targetInfoList.data[selectedMatches[k]].sequenceID),
                                                                      genusTaxIdList[targetInfoList.data[selectedMatches[k]].sequenceID], queryKmerList[j].info.pos, queryKmerList[j].info.frame, selectedHammings[k],0};
                                }
                                posToWrite ++;
                            }
                        }
                        continue;
                    }
                    selectedMatches.clear();
                    selectedHammings.clear();

                    ///Reuse the candidate target k-mers to compare in DNA level if queries are the same at amino acid level but not at DNA level
                    if(currentQueryAA == AminoAcid(queryKmerList[j].ADkmer)){
                        compareDna(queryKmerList[j].ADkmer, candidateTargetKmers, startIdxOfAAmatch, selectedMatches, selectedHammings);
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
                                                                      genusTaxIdList[targetInfoList.data[selectedMatches[k]].sequenceID], queryKmerList[j].info.pos, queryKmerList[j].info.frame, selectedHammings[k],1};
                                } else{
                                    matchBuffer.buffer[posToWrite] = {queryKmerList[j].info.sequenceID, taxID[0]->at(targetInfoList.data[selectedMatches[k]].sequenceID),
                                                                      genusTaxIdList[targetInfoList.data[selectedMatches[k]].sequenceID], queryKmerList[j].info.pos, queryKmerList[j].info.frame, selectedHammings[k],0};
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
                    compareDna(currentQuery, candidateTargetKmers, startIdxOfAAmatch, selectedMatches, selectedHammings);
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
                                                                  genusTaxIdList[targetInfoList.data[selectedMatches[k]].sequenceID], queryKmerList[j].info.pos, queryKmerList[j].info.frame, selectedHammings[k],1};
                            } else{
                                matchBuffer.buffer[posToWrite] = {queryKmerList[j].info.sequenceID, taxID[0]->at(targetInfoList.data[selectedMatches[k]].sequenceID),
                                                                  genusTaxIdList[targetInfoList.data[selectedMatches[k]].sequenceID], queryKmerList[j].info.pos, queryKmerList[j].info.frame, selectedHammings[k],0};
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
    SORT_PARALLEL(matchBuffer.buffer, matchBuffer.buffer + matchBuffer.startIndexOfReserve, Classifier::compareForWritingMatches);
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
void Classifier::compareDna(uint64_t & query, vector<uint64_t> & targetKmersToCompare, const size_t & startIdx, vector<size_t> & selectedMatches, vector<uint8_t> & selectedHamming) {
    vector<uint8_t> hammings;
    uint8_t currentHamming;
    uint8_t minHamming = UINT8_MAX;
    ///Calculate hamming distance
    for(size_t i = 0; i < targetKmersToCompare.size(); i++){
        currentHamming = getHammingDistance(query, targetKmersToCompare[i]);
        if(currentHamming < minHamming)
            minHamming = currentHamming;
        hammings.push_back(currentHamming);
    }

    if(minHamming > 5) {
       return;
    }

    ///Select target k-mers that passed hamming criteria
    for(size_t h = 0; h < hammings.size(); h++){
        if(hammings[h] == minHamming){
            selectedMatches.push_back(startIdx + h);
            selectedHamming.push_back(hammings[h]);
        }
    }
}

///It analyses the result of linear search.
void Classifier::analyseResultParallel(NcbiTaxonomy & ncbiTaxonomy, vector<Sequence> & seqSegments, char * matchFileName, int seqNum, Query * queryList){
    ///Mmap the file of matches
    struct MmapedData<Match> matchList = mmapData<Match>(matchFileName);
    cout<<matchList.fileSize<<"!!"<<endl;
    cout<<matchFileName<<endl;
    size_t numOfMatches = matchList.fileSize / sizeof(Match);
    cout<<"num of matches"<<numOfMatches<<endl;

    ///Sort matches in order to analyze
    SORT_PARALLEL(matchList.data, matchList.data + numOfMatches, Classifier::sortByTaxId);


    ///Devide matches into blocks for multi threading
    typedef Sequence Block;
    Block * matchBlocks = new Block[seqNum];
    cout<<seqNum<<endl;
    size_t matchIdx = 0;
    size_t blockIdx = 0;
    uint32_t currentQuery;
    while(matchIdx < numOfMatches){
        currentQuery = matchList.data[matchIdx].queryId;
        matchBlocks[blockIdx].start = matchIdx;
        while((currentQuery == matchList.data[matchIdx].queryId) && (matchIdx < numOfMatches)) ++matchIdx;
        matchBlocks[blockIdx].end = matchIdx - 1;
        blockIdx++;
    }

    ///Process each blocks
    omp_set_num_threads(ThreadNum);
#pragma omp parallel default(none), shared(cout,matchBlocks, matchList, seqSegments, seqNum, ncbiTaxonomy, queryList)
    {
#pragma omp for schedule(dynamic, 1)
        for(size_t i = 0; i < seqNum; ++ i ){
            TaxID selectedLCA = chooseBestTaxon(ncbiTaxonomy, seqSegments[i].length, i, matchBlocks[i].start,
                                                matchBlocks[i].end, matchList.data, queryList);
        }
    }


    for(int i = 0 ; i < seqNum; i++){
        ++ taxCounts[queryList[i].classification];
    }
    delete[] matchBlocks;
    munmap(matchList.data, matchList.fileSize + 1);
    cout<<"end of analyseResultParallel"<<endl;
}

///For a query read, assign the best Taxon, using k-mer matches
///문제점 redundancy reduced reference k-mer 임을 고려해야 한다. block을 species level에서 해줘야하지 않나.. 그리고 오버랩도 좀 허용해줘야할껄?
TaxID Classifier::chooseBestTaxon(NcbiTaxonomy & ncbiTaxonomy, const size_t & queryLength, const int & currentQuery, const size_t & offset, const size_t & end, Match * matchList, Query * queryList){
    vector<ConsecutiveMatches> matchCombi;
    getBestGenusLevelMatchCombination(matchCombi, matchList, end, offset);

    //un-classified
    if(matchCombi.empty()){
        return 0;
    }


    ///Check a query coverage
    float coverageThr = 0.2;
    int maxNum = queryLength / 3 - kmerLength + 1;
    int coveredKmerCnt = 0;
    float coverage;
    for(size_t cm = 0 ; cm < matchCombi.size(); cm ++){
        coveredKmerCnt += matchCombi[cm].diffPosCnt;
    }
    coverage = float(coveredKmerCnt) / float(maxNum);
    queryList[currentQuery].coverage = coverage;


    vector<TaxID> taxIdList;
//    vector<uint32_t> pos;
//    vector<uint8_t> frame;
//    vector<uint8_t> ham;
//    vector<int> redun;
    TaxID temp;
    for(size_t cs = 0; cs < matchCombi.size(); cs++ ){
        for(size_t k = matchCombi[cs].beginIdx ; k < matchCombi[cs].endIdx + 1; k++ ){
            temp = matchList[k].taxID;
            taxIdList.push_back(temp);
//            pos.push_back(matchList[k].position);
//            frame.push_back(matchList[k].frame);
//            ham.push_back(matchList[k].hamming);
//            redun.push_back(matchList[k].red);
            queryList[currentQuery].taxCnt[temp] ++;
        }
    }

    ///No classification for low query coverage.
    if(coverage < coverageThr) {
       return 0;
    }

    size_t numAssignedSeqs = 0;
    size_t numUnassignedSeqs = 0;
    size_t numSeqsAgreeWithSelectedTaxon = 0;
    double selectedPercent = 0;

    TaxID selectedLCA = match2LCA(taxIdList, ncbiTaxonomy, 0.7, numAssignedSeqs,
                                  numUnassignedSeqs, numSeqsAgreeWithSelectedTaxon,
                                  selectedPercent, queryLength);

    ///TODO optimize strain specific classification criteria
    ///Strain classification only for high coverage with LCA of species level
    if(coverage > 0.75 && NcbiTaxonomy::findRankIndex(ncbiTaxonomy.taxonNode(selectedLCA)->rank) == 4){ /// There are more strain level classifications with lower coverage threshold, but also with more false postives. 0.8~0.85 looks good.
        int strainCnt = 0;
        unordered_map<TaxID, int> strainMatchCnt;
        TaxID strainTaxId;

        for(size_t cs = 0; cs < matchCombi.size(); cs++ ){
            for(size_t k = matchCombi[cs].beginIdx ; k < matchCombi[cs].endIdx + 1; k++ ){
                temp = matchList[k].taxID;
                if(selectedLCA != temp && ncbiTaxonomy.IsAncestor(selectedLCA, temp)){
                    if(strainMatchCnt.find(temp) == strainMatchCnt.end()){
                        strainCnt ++;
                        strainTaxId = temp;
                        strainMatchCnt.insert(pair<TaxID, int>(temp, 1));
                    } else {
                        strainMatchCnt[temp] ++;
                    }
                }
            }
        }
        if(strainCnt == 1){
            selectedLCA = strainTaxId;
        }
    }


//    cout<<"# "<<currentQuery<<endl;
//    for(size_t i = 0; i < taxIdList.size(); i++){
//        cout<<i<<" "<<int(frame[i])<<" "<<pos[i]<<" "<<taxIdList[i]<<" "<<int(ham[i])<<" "<<redun[i]<<endl;
//    }
//    cout<<"coverage: "<<coverage<<"  "<<selectedLCA<<" "<<ncbiTaxonomy.taxonNode(selectedLCA)->rank<<endl;
    ///store classification results
    queryList[currentQuery].isClassified = true;
    queryList[currentQuery].classification = selectedLCA;
    queryList[currentQuery].coverage = coverage;
    return selectedLCA;
}

void Classifier::getBestGenusLevelMatchCombination(vector<ConsecutiveMatches> & chosenMatchCombination, Match * matchList, size_t end, size_t offset){
    vector<ConsecutiveMatches> coMatches;
    vector<vector<ConsecutiveMatches>> genus;
    int conCnt = 0;
    uint32_t diffPosCnt = 0;
    uint32_t hammingSum = 0;
    float hammingMean = 0.0;
    size_t beginIdx = 0;
    size_t endIdx = 0;
    uint32_t conBegin = 0;
    uint32_t conEnd = 0;
    uint32_t currentPos;
    uint8_t currentFrame;
    TaxID currentTaxID;


    size_t i = offset;
    while(i < end + 1) {
        currentTaxID = matchList[i].genusTaxID;
        //For current genus
        while (currentTaxID == matchList[i].genusTaxID && (i < end + 1)) {
            currentFrame = matchList[i].frame;
            //For current frame
            while (currentFrame == matchList[i].frame && currentTaxID == matchList[i].genusTaxID && (i < end + 1)){
                currentPos = matchList[i].position;
                hammingSum = 0;
                hammingMean = 0.0;
                conCnt = 0;
                diffPosCnt = 1;
                conBegin = currentPos;
                beginIdx = i;
                //Find consecutive matches
                //TODO: this can be faster
                while(matchList[i].position <= currentPos + 3 &&
                    (conCnt == 0 || matchList[i].hamming <= hammingMean + 3) && ///TODO: Is it okay?
                    currentFrame == matchList[i].frame &&
                    currentTaxID == matchList[i].genusTaxID && (i < end + 1)){
                    if(matchList[i].position != currentPos) {
                        diffPosCnt++;
                        currentPos = matchList[i].position;
                    }
                    conCnt ++;
                    hammingSum += matchList[i].hamming;
                    hammingMean = float(hammingSum) / float(conCnt);
                    i++;
                }
                if(diffPosCnt > 1){
                    coMatches.emplace_back(conBegin, currentPos, conCnt, hammingSum, diffPosCnt, beginIdx, i-1, currentFrame);
                }
            }
        }
        //choose the best combination of consecutive matches for current genus
        if(!coMatches.empty()) getMatchCombinationForCurGenus2(coMatches, genus, matchList);
        coMatches.clear();
    }
    //choose the best combination of consecutive-match among genus for current query
    if(!genus.empty())
        getTheBestGenus(genus, chosenMatchCombination);
}
void Classifier::getMatchCombinationForCurGenus(vector<ConsecutiveMatches> & coMatches, vector<vector<ConsecutiveMatches>> & genus, Match * matchList){
    //cout<<"663"<<endl;
    //    for(int i3 = 0; i3 < coMatches.size(); i3++){
//        cout<< coMatches[i3].begin << " " << coMatches[i3].end << " "<< coMatches[i3].matchCnt;
//        cout<<" "<<coMatches[i3].hamming << " "<<int(coMatches[i3].frame)<<endl;
//        cout<<matchList[coMatches[i3].beginIdx].taxID<<endl;
//        cout<<matchList[coMatches[i3].endIdx].taxID<<endl;
//    }
    int numOfSubsets = int(pow(2, coMatches.size()) - 1);
    vector<int> subset;
    vector<vector<int>> subsetList;
    size_t bestSubset;
    float bestScore = -FLT_MAX;
    float tiedScore;
    size_t tiedSubset;
    float currentScore = 0;
    getSubsets(subset, subsetList, 0, coMatches.size() - 1);
    vector<ConsecutiveMatches> matchesToScore;
    for(size_t j = 0; j < subsetList.size(); j++){
        if(!subsetList[j].empty()) {
            for (size_t i = 0; i < subsetList[j].size(); i++) {
                matchesToScore.push_back(coMatches[subsetList[j][i]]);
            }
            currentScore = scoreSubset(matchesToScore);
            matchesToScore.clear();
            if(currentScore > bestScore){
                bestScore = currentScore;
                bestSubset = j;
            } else if(currentScore == bestScore){
                tiedScore = bestScore;
                tiedSubset = j;
            }
        }
    }
    vector<ConsecutiveMatches> alignedCoMatches;
    for(size_t i = 0; i < subsetList[bestSubset].size(); i++){
        alignedCoMatches.push_back(coMatches[subsetList[bestSubset][i]]);
    }

    if(tiedScore == bestScore){
        cout<<"tie "<<tiedScore<<" "<<bestScore<<" "<<subsetList[tiedSubset].size()<<endl;
        for(size_t i = 0; i < subsetList[tiedSubset].size(); i++){
            alignedCoMatches.push_back(coMatches[subsetList[tiedSubset][i]]);
        }
    }

    //genus.push_back(coMatches);
    genus.push_back(alignedCoMatches);
}

void Classifier::getMatchCombinationForCurGenus2(vector<ConsecutiveMatches> & coMatches, vector<vector<ConsecutiveMatches>> & genus, Match * matchList){
    //sort consecutive match blocks
    sort(coMatches.begin(), coMatches.end(), Classifier::compareConsecutiveMatches);

    //container to store match blocks to be used.
    vector<ConsecutiveMatches> alignedCoMatches;

    //store the best one anyway
    alignedCoMatches.push_back(coMatches[0]);
    bool overlap = false;

    //TODO: Fix here to accept slight overlaps
    for(size_t i = 1; i < coMatches.size(); i++){
        overlap = false;
        for(size_t j = 0; j < alignedCoMatches.size(); j++){
            if((alignedCoMatches[j].begin < coMatches[j].end) && (alignedCoMatches[j].end > coMatches[j].begin)){ ///TODO check this condition
                overlap = true;
                break;
            }
        }

        if(overlap) continue;
        else{
            alignedCoMatches.push_back(coMatches[i]);
        }
    }

    genus.push_back(alignedCoMatches);
}

void Classifier::getSubsets(vector<int> & subset, vector<vector<int>> & uniqueSubset, int k, int n){
    if(k == n + 1){
        uniqueSubset.push_back(subset);
    } else{
        subset.push_back(k);
        getSubsets(subset, uniqueSubset, k+1, n);
        subset.pop_back();
        getSubsets(subset, uniqueSubset, k+1, n);
    }
}

float Classifier::scoreSubset(vector<ConsecutiveMatches> & subset){
    float score;
    bool overlapped = false;
    //give bad score if overlapping
    for(size_t i = 0; i + 1 < subset.size(); i++){
        for(size_t j = i + 1; j < subset.size(); j++){
            if((subset[i].begin <= subset[j].end) && (subset[i].end >= subset[j].begin))
                overlapped = true;
        }
    }
    if(overlapped) return -FLT_MAX;

    for(size_t i = 0; i < subset.size(); i++)
        score += float(subset[i].diffPosCnt) - float(subset[i].hamming);
    return score;
}
void Classifier::getTheBestGenus(vector<vector<ConsecutiveMatches>> & genus, vector<ConsecutiveMatches> & chosen){
    int chosenGenusIdx = INT_MAX;
    vector<TaxID> selecetedGenusList;
    int totalDiffPosCnt;
    int totalMatchCnt;
    int totalHamming;
    float maxScore = -FLT_MAX;
    float currScore;

    for(size_t i = 0; i < genus.size(); i++){
        totalDiffPosCnt = 0;
        totalMatchCnt = 0;
        totalHamming = 0;
        for(size_t j = 0; j < genus[i].size(); j++){
            totalDiffPosCnt += genus[i][j].diffPosCnt;
            totalMatchCnt += genus[i][j].matchCnt;
            totalHamming += genus[i][j].hamming;
        }
        currScore = totalDiffPosCnt - float(totalHamming)/float(totalMatchCnt);
        if(currScore > maxScore){
            chosenGenusIdx = i;
            selecetedGenusList.clear();
            selecetedGenusList.push_back(i);
            maxScore = currScore;
        } else if (currScore == maxScore){
            selecetedGenusList.push_back(i);
        }
    }

    for(size_t g = 0; g < selecetedGenusList.size(); g++) {
        for (size_t i = 0; i < genus[selecetedGenusList[g]].size(); i++) {
            chosen.push_back(genus[selecetedGenusList[g]][i]);
        }
    }
}


void Classifier::findConsecutiveMatches(vector<ConsecutiveMatches> & coMatches, Match * matchList, size_t end, size_t offset){
    int conCnt = 0;
    uint32_t hammingSum = 0;
    size_t beginIdx = 0;
    size_t endIdx = 0;
    ///This routine is for getting consecutive matched k-mer
    ///gapThr decides the maximun gap
    uint8_t currentFrame;
    int gapThr = 0;
    TaxID currentTaxID;
    uint32_t conBegin = 0;
    uint32_t conEnd = 0;
    uint32_t gapCnt = 0;
    size_t i = offset;
    while(i < end){
        currentTaxID = matchList[i].genusTaxID;
        currentFrame = matchList[i].frame;
        while(currentFrame == matchList[i+1].frame && currentTaxID == matchList[i+1].genusTaxID && (i < end)) {
            if (matchList[i + 1].position <= matchList[i].position + (gapThr + 1) * 3) {
                if (conCnt == 0) {
                    conBegin = matchList[i].position;
                    beginIdx = i;
                }
                conCnt++;
                hammingSum += matchList[i].hamming;
                if (matchList[i + 1].position != matchList[i].position) {
                    gapCnt += (matchList[i + 1].position - matchList[i].position) / 3 - 1;
                }
            } else {
                if (conCnt > 0) {
                    conCnt++;
                    hammingSum += matchList[i].hamming;
                    conEnd = matchList[i].position;
                    endIdx = i;
                    if(conBegin != conEnd)
                        coMatches.emplace_back(conBegin, conEnd, conCnt, hammingSum, gapCnt, beginIdx, endIdx, currentFrame);
                    conCnt = 0;
                    gapCnt = 0;
                    hammingSum = 0;
                }
            }
            i++;
        }
        if (conCnt > 0) {
            conCnt++;
            hammingSum += matchList[i].hamming;
            conEnd = matchList[i].position;
            endIdx = i;
            if(conBegin != conEnd)
                coMatches.emplace_back(conBegin, conEnd, conCnt, hammingSum, gapCnt, beginIdx, endIdx, currentFrame);
            conCnt = 0;
            gapCnt = 0;
            hammingSum = 0;
        }
        i++;
    }
}

TaxID Classifier::match2LCA(const std::vector<int> & taxIdList, NcbiTaxonomy const & taxonomy, const float majorityCutoff,
                            size_t &numAssignedSeqs, size_t &numUnassignedSeqs, size_t &numSeqsAgreeWithSelectedTaxon, double &selectedPercent, uint32_t queryLength){
    std::map<TaxID,taxNode> ancTaxIdsCounts;

    numAssignedSeqs = 0;
    numUnassignedSeqs = 0;
    numSeqsAgreeWithSelectedTaxon = 0;
    selectedPercent = 0;
    double totalAssignedSeqsWeights = 0.0;

    for (size_t i = 0; i < taxIdList.size(); ++i) {
        TaxID currTaxId = taxIdList[i];
        double currWeight = 1;
        // ignore unassigned sequences
        if (currTaxId == 0) {
            numUnassignedSeqs++;
            continue;
        }
        TaxonNode const * node = taxonomy.taxonNode(currTaxId, false);
        if (node == NULL) {
            Debug(Debug::ERROR) << "taxonid: " << currTaxId << " does not match a legal taxonomy node.\n";
            EXIT(EXIT_FAILURE);
        }
        totalAssignedSeqsWeights += currWeight;
        numAssignedSeqs++;

        // each start of a path due to an orf is a candidate
        if (ancTaxIdsCounts.find(currTaxId) != ancTaxIdsCounts.end()) { //원소가 있다면
            ancTaxIdsCounts[currTaxId].update(currWeight, 0);
        } else {
            taxNode currNode;
            currNode.set(currWeight, true, 0);
            ancTaxIdsCounts.insert(std::pair<TaxID,taxNode>(currTaxId, currNode));
        }

        // iterate all ancestors up to root (including). add currWeight and candidate status to each
        TaxID currParentTaxId = node->parentTaxId;
        while (currParentTaxId != currTaxId) {
            if (ancTaxIdsCounts.find(currParentTaxId) != ancTaxIdsCounts.end()) {
                ancTaxIdsCounts[currParentTaxId].update(currWeight, currTaxId);
            } else {
                taxNode currParentNode;
                currParentNode.set(currWeight, false, currTaxId);
                ancTaxIdsCounts.insert(std::pair<TaxID,taxNode>(currParentTaxId, currParentNode));
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
    float maxCoverage = -FLT_MAX;
    float spMaxCoverage = -FLT_MAX;

    int maximumKmerNum = queryLength / 3 - kmerLength + 1;
    bool haveMetCovThr = false;

    //한 위치에 중복되는 매치가 있다! 잘 생각해봅시다...
    for (std::map<TaxID,taxNode>::iterator it = ancTaxIdsCounts.begin(); it != ancTaxIdsCounts.end(); it++) {
        // consider only candidates:
        if (!(it->second.isCandidate)) {
            continue;
        }

        double currPercent = float(it->second.weight) / totalAssignedSeqsWeights;
        curCoverage = float(it->second.weight) / float(maximumKmerNum);

        TaxID currTaxId = it->first;
        TaxonNode const * node = taxonomy.taxonNode(currTaxId, false);
        int currRankInd = NcbiTaxonomy::findRankIndex(node->rank);
        if(curCoverage > coverageThreshold && currRankInd <= 4){
            if(!haveMetCovThr){
                haveMetCovThr = true;
                minRank = currRankInd;
                spMaxCoverage = curCoverage;
                selectedTaxon = it->first;
            } else if(currRankInd < minRank || (currRankInd == minRank && curCoverage > spMaxCoverage)){
                minRank = currRankInd;
                spMaxCoverage = curCoverage;
                selectedTaxon = it->first;
            }
        } else if (currPercent >= majorityCutoff && (!haveMetCovThr)) {
            // iterate all ancestors to find lineage min rank (the candidate is a descendant of a node with this rank)

            // int currMinRank = INT_MAX;
            // TaxID currParentTaxId = node->parentTaxId;
//            while (currParentTaxId != currTaxId) {
//                if ((currRankInd > 0) && (currRankInd < currMinRank)) {
//                    currMinRank = currRankInd;
//                    // the rank can only go up on the way to the root, so we can break
//                    break;
//                }
//                // move up:
//                currTaxId = currParentTaxId;
//                node = taxonomy.taxonNode(currParentTaxId, false);
//                currParentTaxId = node->parentTaxId;
//            }
            if ((currRankInd < minRank) || ((currRankInd == minRank) && (currPercent > selectedPercent))) {
                selectedTaxon = it->first;
                minRank = currRankInd;
                selectedPercent = currPercent;
            }
        }
    }

    return selectedTaxon;
}

TaxID Classifier::match2LCA2(const std::vector<int> & taxIdList, NcbiTaxonomy const & taxonomy, const float majorityCutoff,
                            size_t &numAssignedSeqs, size_t &numUnassignedSeqs, size_t &numSeqsAgreeWithSelectedTaxon, double &selectedPercent){
    std::map<TaxID,taxNode> ancTaxIdsCounts;

    numAssignedSeqs = 0;
    numUnassignedSeqs = 0;
    numSeqsAgreeWithSelectedTaxon = 0;
    selectedPercent = 0;
    double totalAssignedSeqsWeights = 0.0;

    for (size_t i = 0; i < taxIdList.size(); ++i) {
        TaxID currTaxId = taxIdList[i];
        double currWeight = 1;
        // ignore unassigned sequences
        if (currTaxId == 0) {
            numUnassignedSeqs++;
            continue;
        }
        TaxonNode const * node = taxonomy.taxonNode(currTaxId, false);
        if (node == NULL) {
            Debug(Debug::ERROR) << "taxonid: " << currTaxId << " does not match a legal taxonomy node.\n";
            EXIT(EXIT_FAILURE);
        }
        totalAssignedSeqsWeights += currWeight;
        numAssignedSeqs++;

        // each start of a path due to an orf is a candidate
        if (ancTaxIdsCounts.find(currTaxId) != ancTaxIdsCounts.end()) { //원소가 있다면
            ancTaxIdsCounts[currTaxId].update(currWeight, 0);
        } else {
            taxNode currNode;
            currNode.set(currWeight, true, 0);
            ancTaxIdsCounts.insert(std::pair<TaxID,taxNode>(currTaxId, currNode));
        }

        // iterate all ancestors up to root (including). add currWeight and candidate status to each
        TaxID currParentTaxId = node->parentTaxId;
        while (currParentTaxId != currTaxId) {
            if (ancTaxIdsCounts.find(currParentTaxId) != ancTaxIdsCounts.end()) {
                ancTaxIdsCounts[currParentTaxId].update(currWeight, currTaxId);
            } else {
                taxNode currParentNode;
                currParentNode.set(currWeight, false, currTaxId);
                ancTaxIdsCounts.insert(std::pair<TaxID,taxNode>(currParentTaxId, currParentNode));
            }
            // move up:
            currTaxId = currParentTaxId;
            node = taxonomy.taxonNode(currParentTaxId, false);
            currParentTaxId = node->parentTaxId;
        }
    }

    // select the lowest ancestor that meets the cutoff
    int minRank = INT_MAX;
    int weightOfMinRank;
    int currRank;
    TaxID selectedTaxon = 0;

    for (auto it = ancTaxIdsCounts.begin(); it != ancTaxIdsCounts.end(); it++) {
        // consider only candidates:
        if (!(it->second.isCandidate)) {
            continue;
        }

        if (it->second.weight > 25){
            TaxID currTaxId = it->first;
            TaxonNode const * node = taxonomy.taxonNode(currTaxId, false);
            currRank = NcbiTaxonomy::findRankIndex(node->rank);
            if(currRank == -1) continue;
            if((currRank < minRank) || (currRank == minRank && it->second.weight > weightOfMinRank)){
                minRank = currRank;
                weightOfMinRank = it->second.weight;
                selectedTaxon = it->first;
            }
        }
        double currPercent = float(it->second.weight) / totalAssignedSeqsWeights;
    }
    return selectedTaxon;
}

TaxID Classifier::match2LCA3(const std::vector<int> & taxIdList, NcbiTaxonomy const & taxonomy, const float majorityCutoff,
                            size_t &numAssignedSeqs, size_t &numUnassignedSeqs, size_t &numSeqsAgreeWithSelectedTaxon, double &selectedPercent){
    std::map<TaxID,taxNode> ancTaxIdsCounts;

    numAssignedSeqs = 0;
    numUnassignedSeqs = 0;
    numSeqsAgreeWithSelectedTaxon = 0;
    selectedPercent = 0;
    double totalAssignedSeqsWeights = 0.0;

    for (size_t i = 0; i < taxIdList.size(); ++i) {
        TaxID currTaxId = taxIdList[i];
        double currWeight = 1;
        // ignore unassigned sequences
        if (currTaxId == 0) {
            numUnassignedSeqs++;
            continue;
        }
        TaxonNode const * node = taxonomy.taxonNode(currTaxId, false);
        if (node == NULL) {
            Debug(Debug::ERROR) << "taxonid: " << currTaxId << " does not match a legal taxonomy node.\n";
            EXIT(EXIT_FAILURE);
        }
        totalAssignedSeqsWeights += currWeight;
        numAssignedSeqs++;

        // each start of a path due to an orf is a candidate
        if (ancTaxIdsCounts.find(currTaxId) != ancTaxIdsCounts.end()) { //원소가 있다면
            ancTaxIdsCounts[currTaxId].update(currWeight, 0);
        } else {
            taxNode currNode;
            currNode.set(currWeight, true, 0);
            ancTaxIdsCounts.insert(std::pair<TaxID,taxNode>(currTaxId, currNode));
        }

        // iterate all ancestors up to root (including). add currWeight and candidate status to each
        TaxID currParentTaxId = node->parentTaxId;
        while (currParentTaxId != currTaxId) {
            if (ancTaxIdsCounts.find(currParentTaxId) != ancTaxIdsCounts.end()) {
                ancTaxIdsCounts[currParentTaxId].update(currWeight, currTaxId);
            } else {
                taxNode currParentNode;
                currParentNode.set(currWeight, false, currTaxId);
                ancTaxIdsCounts.insert(std::pair<TaxID,taxNode>(currParentTaxId, currParentNode));
            }
            // move up:
            currTaxId = currParentTaxId;
            node = taxonomy.taxonNode(currParentTaxId, false);
            currParentTaxId = node->parentTaxId;
        }
    }

    // select the lowest ancestor that meets the cutoff
    int minRank = INT_MAX;
    TaxID selctedTaxon = 0;

    for (std::map<TaxID,taxNode>::iterator it = ancTaxIdsCounts.begin(); it != ancTaxIdsCounts.end(); it++) {
        // consider only candidates:
        if (!(it->second.isCandidate)) {
            continue;
        }

        double currPercent = float(it->second.weight) / totalAssignedSeqsWeights;
        if (currPercent >= majorityCutoff) {
            // iterate all ancestors to find lineage min rank (the candidate is a descendant of a node with this rank)
            TaxID currTaxId = it->first;
            TaxonNode const * node = taxonomy.taxonNode(currTaxId, false);
            int currMinRank = INT_MAX;
            TaxID currParentTaxId = node->parentTaxId;
            while (currParentTaxId != currTaxId) {
                int currRankInd = NcbiTaxonomy::findRankIndex(node->rank);
                if ((currRankInd > 0) && (currRankInd < currMinRank)) {
                    currMinRank = currRankInd;
                    // the rank can only go up on the way to the root, so we can break
                    break;
                }
                // move up:
                currTaxId = currParentTaxId;
                node = taxonomy.taxonNode(currParentTaxId, false);
                currParentTaxId = node->parentTaxId;
            }

            if ((currMinRank < minRank) || ((currMinRank == minRank) && (currPercent > selectedPercent))) {
                selctedTaxon = it->first;
                minRank = currMinRank;
                selectedPercent = currPercent;
            }
        }
    }

    return selctedTaxon;
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

bool Classifier::compareConsecutiveMatches(const ConsecutiveMatches & a, const ConsecutiveMatches & b){
    if((a.end - a.begin) > (b.end - b.begin)){ // compare query coverage fisrt
        if((a.end - a.begin) == (b.end - b.begin + 1)){ //
            return (a.endIdx - a.beginIdx + 1) * 2 / ((a.hamming+1)*(a.diffPosCnt + 1)) > (b.endIdx - b.beginIdx + 1) * 2 / ((b.hamming + 1) * (b.diffPosCnt + 1));
        }
        return true;
    }else if((a.end - a.begin) == (b.end- b.begin)){
            return (a.endIdx - a.beginIdx + 1) * 2 / ((a.hamming+1)*(a.diffPosCnt + 1)) > (b.endIdx - b.beginIdx + 1) * 2 / ((b.hamming + 1) * (b.diffPosCnt + 1));
    }
    return false;
}

bool Classifier::compareConsecutiveMatches2(const ConsecutiveMatches & a, const ConsecutiveMatches & b){
    if(a.matchCnt - a.hamming > b.matchCnt - b.hamming) return true;
    return false;
}

bool Classifier::compareConsecutiveMatches3(const ConsecutiveMatches & a, const ConsecutiveMatches & b){
    if(a.diffPosCnt - a.hamming > b.diffPosCnt - b.hamming) return true;
    return false;
}

bool Classifier::sortByTaxId(const Match & a, const Match & b){
    if (a.queryId < b.queryId) return true;
    else if (a.queryId == b.queryId) {
        if(a.genusTaxID < b.genusTaxID) return true;
        else if(a.genusTaxID == b.genusTaxID) {
            if (a.frame < b.frame) return true;
            else if (a.frame == b.frame) {
                if (a.position < b.position) return true;
            }
        }
    }
    return false;
}

void Classifier::writeReadClassification(Query * queryList, int queryNum, ofstream & readClassificationFile){
    for(size_t i = 0; i < queryNum; i++){
        readClassificationFile <<i<< "\t" << queryList[i].isClassified << "\t" << queryList[i].name << "\t" << queryList[i].classification << "\t" << queryList[i].queryLength << "\t" << queryList[i].coverage << "\t";
        for(auto it = queryList[i].taxCnt.begin(); it != queryList[i].taxCnt.end(); ++it){
            readClassificationFile<<it->first<<":"<<it->second<<" ";
        }
        readClassificationFile<<endl;
    }
}

void Classifier::writeReportFile(const char * queryFileName, NcbiTaxonomy & ncbiTaxonomy, const int numOfQuery){
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

void Classifier::performanceTest(NcbiTaxonomy & ncbiTaxonomy, Query * queryList, int numOfquery){

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

    for(size_t i = 0; i < numOfquery; i++) {
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
            compareTaxon(classificationResult, rightAnswer, ncbiTaxonomy);

        }
    }

    cout<<"Number of classification: "<< classifiedCnt << endl;
    cout<<"classified / total =" << float(classifiedCnt)/float(queryInfos.size()) << endl;
    cout<<"Superkingdom: "<< superCnt <<endl;
    cout<<"Phylum: "<<phylumCnt<<endl;
    cout<<"Class: "<<classCnt<<endl;
    cout<<"Order: "<<orderCnt<<endl;
    cout<<"Family: "<<familyCnt<<endl;
    cout<<"Genus: "<< genusCnt << endl;
    cout<<"Species: "<<speciesCnt<<endl;
    cout<<"Subspecies: "<<subspCnt<<endl;
    cout<<"(subS + S + G) / all classification" << float(genusCnt + speciesCnt + subspCnt) / float(classifiedCnt) <<endl;
    cout<<"Num of queries: " << queryInfos.size() << endl;
}

void Classifier::compareTaxon(TaxID shot, TaxID target, NcbiTaxonomy & ncbiTaxonomy) { ///target: subspecies or species
    const TaxonNode * shotNode = ncbiTaxonomy.taxonNode(shot);
    string shotRank = shotNode->rank;
    cout<<shot<<" "<<target<<" "<<shotRank<<" ";
    if(NcbiTaxonomy::findRankIndex(shotRank) <= 3){
        //cout<<"subspecies"<<endl;
        if(shot == target){
            cout<<"O"<<endl;
            subspCnt ++;
        } else cout<<"X"<<endl;

    } else if(shotRank == "species") {
        //cout<<"species"<<endl;
        if(shot == ncbiTaxonomy.getTaxIdAtRank(target, "species")){
            speciesCnt ++;
            cout<<"O"<<endl;
        }else cout<<"X"<<endl;
    } else if(shotRank == "genus"){
        //cout<<"genus"<<endl;
        if(shot == ncbiTaxonomy.getTaxIdAtRank(target, "genus")){
            genusCnt ++;
            cout<<"O"<<endl;
        } else cout<<"X"<<endl;
    } else if(shotRank == "family"){
        //cout<<"family"<<endl;
        if(shot == ncbiTaxonomy.getTaxIdAtRank(target, "family")) {
            familyCnt++;
        }
    }else if(shotRank == "order") {
        //cout<<"order"<<endl;
        if(shot == ncbiTaxonomy.getTaxIdAtRank(target, "order")) {
            orderCnt++;
        }
    }else if(shotRank == "class") {
        //cout<<"class"<<endl;
        if(shot == ncbiTaxonomy.getTaxIdAtRank(target, "class")) {
            classCnt++;
        }
    } else if(shotRank == "phylum") {
        //cout<<"phylum"<<endl;
        if(shot == ncbiTaxonomy.getTaxIdAtRank(target, "phylum")) {
            phylumCnt++;
        }
    } else {
        return;
    }
}
