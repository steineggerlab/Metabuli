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
    if(par.gtdbOrNcbi == 1 || par.gtdbOrNcbi == 0){
        cout<<"Classifying query sequences based on taxonomy of GTDB"<<endl;
        names = "../../gtdb_taxdmp/names.dmp";
        nodes = "../../gtdb_taxdmp/nodes.dmp";
        merged = "../../gtdb_taxdmp/merged.dmp";
    } else if(par.gtdbOrNcbi == 2){
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
    vector<int> taxIdListAtRank;
    ncbiTaxonomy.createTaxIdListAtRank(taxIdList, taxIdListAtRank, "species");

    //output file
    char matchFileName[300];
    sprintf(matchFileName,"%s_match", queryFileName);
    FILE * matchFile = fopen(matchFileName, "wb");
    //query & target
    struct MmapedData<char> queryFile = mmapData<char>(par.filenames[0].c_str());
    cout<<"1"<<endl;

    struct MmapedData<uint16_t> targetDiffIdxList = mmapData<uint16_t>(targetDiffIdxFileName);
    cout<<"fileSize: "<<targetDiffIdxList.fileSize<<endl;
    //cout<<targetDiffIdxList.data[targetDiffIdxList.fileSize / sizeof(uint16_t) - 1]<<endl;
    //targetDiffIdxList.data[targetDiffIdxList.fileSize / sizeof(uint16_t)] = 32768; //1000000000000000

    struct MmapedData<TargetKmerInfo> targetInfoList = mmapData<TargetKmerInfo>(targetInfoFileName);
    cout<<"4"<<endl;

    struct MmapedData<DiffIdxSplit> diffIdxSplits = mmapData<DiffIdxSplit>(diffIdxSplitFileName);
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
    ofstream readClassificationFile;
    readClassificationFile.open(par.filenames[0]+"_ReadClassification_temp.tsv");
    time_t beforeSearch, afterSearch;
    cout<<"1"<<endl;
    size_t numOfTatalQueryKmerCnt = 0;
    while(processedSeqCnt < numOfSeq){
        fillQueryKmerBufferParallel(kmerBuffer, queryFile, sequences, processedSeqChecker, processedSeqCnt, queryList);
        numOfTatalQueryKmerCnt += kmerBuffer.startIndexOfReserve;
        cout<<"buffer overflowed"<<endl;
        beforeSearch = time(NULL);
        omp_set_num_threads(64);
        SORT_PARALLEL(kmerBuffer.buffer, kmerBuffer.buffer + kmerBuffer.startIndexOfReserve, Classifier::compareForLinearSearch);
        cout<<"buffer sorted"<<endl;
        linearSearchParallel(kmerBuffer.buffer, kmerBuffer.startIndexOfReserve, targetDiffIdxList, targetInfoList, diffIdxSplits, matchBuffer, taxIdList, taxIdListAtRank, matchFile);
    }
    cout<<"total kmer count: "<<numOfTatalQueryKmerCnt<<endl;
    writeMatches(matchBuffer, matchFile);
    fclose(matchFile);
    afterSearch = time(NULL);
    cout<<"Time spent for searching: "<<double(afterSearch-beforeSearch)<<endl;

    //load matches and analyze
    cout<<"analyse Result"<<endl;
    //analyseResult(ncbiTaxonomy, sequences, matchFileName, queryList);
    analyseResultParallel(ncbiTaxonomy, sequences, matchFileName, numOfSeq);

    writeReadClassification(queryInfos,readClassificationFile);

    ///TODO split count 고려할 것
    cout<<"Sorting the 'queryfile_ReadClassification.tsv' file"<<endl;
    string sortCall = "sort -t '\t' -k1 -n " + par.filenames[0] + "_ReadClassification_temp.tsv > "+par.filenames[0]+"_ReadClassification.tsv";
    string rmCall = "rm " +par.filenames[0]+"_ReadClassification_temp.tsv";
    system(sortCall.c_str());
    system(rmCall.c_str());
    readClassificationFile.close();

    ///TODO: Merge ReportFiles

    writeReportFile(par.filenames[0].c_str(), ncbiTaxonomy, numOfSeq);
    performanceTest(ncbiTaxonomy);

    cout<<"Number of query k-mer                : "<<queryCount<<endl;
    cout<<"Number of total match                : "<<totalMatchCount <<endl;
    cout<<"mutipleMatch in AA level             : "<<multipleMatchCount << endl;
    cout<<"matches in DNA level                 : "<<perfectMatchCount<<endl;

    free(kmerBuffer.buffer);
    free(matchBuffer.buffer);
    munmap(queryFile.data, queryFile.fileSize + 1);
    munmap(targetDiffIdxList.data, targetDiffIdxList.fileSize + 1);
    munmap(targetInfoList.data, targetInfoList.fileSize + 1);
}

void Classifier::fillQueryKmerBufferParallel(QueryKmerBuffer & kmerBuffer, MmapedData<char> & seqFile, vector<Sequence> & seqs, bool * checker, size_t & processedSeqCnt, Query * queryList) {
    bool hasOverflow = false;
    omp_set_num_threads(ThreadNum);
#pragma omp parallel default(none), shared(checker, hasOverflow, processedSeqCnt, kmerBuffer, seqFile, seqs, queryInfos, cout, queryList)
    {
        vector<QueryInfo> infos;
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
                    infos.emplace_back(int(i), false, buffer.entry.name.s, 0, 0, seqs[i].length);
                    seqs[i].length = strlen(buffer.entry.sequence.s);
                    queryList[i].queryLength = seqs[i].length;
                    queryList[i].queryId = i;
                    queryList[i].name = buffer.entry.sequence.s;
#pragma omp atomic
                    processedSeqCnt ++;
                } else{
                    #pragma omp atomic
                    kmerBuffer.startIndexOfReserve -= kmerCnt;
                    hasOverflow = true;
                }
            }

        }
#pragma omp critical
        {
            queryInfos.insert(queryInfos.end(), make_move_iterator(infos.begin()), make_move_iterator(infos.end()));
        }
    }
}

void Classifier::linearSearchParallel(QueryKmer * queryKmerList, size_t & queryKmerCnt, const MmapedData<uint16_t> & targetDiffIdxList,
                                      const MmapedData<TargetKmerInfo> & targetInfoList, const MmapedData<DiffIdxSplit> & diffIdxSplits,
                                      Buffer<Match> & matchBuffer, const vector<int> & taxIdList, const vector<int> & taxIdListAtRank,
                                      FILE * matchFile){
    ///Find the first index of garbage query k-mer (UINT64_MAX) and discard from there
    for(size_t checkN = queryKmerCnt - 1; checkN >= 0; checkN--){
        if(queryKmerList[checkN].ADkmer != UINT64_MAX){
            queryKmerCnt = checkN + 1;
            break;
        }
    }

    ///Filter out meaningless target splits
    size_t numOfDiffIdxSplits = diffIdxSplits.fileSize / sizeof(DiffIdxSplit);
    size_t numOfDiffIdxSplits_use = numOfDiffIdxSplits;
    for(size_t i = 1; i < numOfDiffIdxSplits; i++){
        if(diffIdxSplits.data[i].ADkmer == 0 || diffIdxSplits.data[i].ADkmer == UINT64_MAX){
            diffIdxSplits.data[i] = {UINT64_MAX, UINT64_MAX, UINT64_MAX};
            numOfDiffIdxSplits_use--;
        }
    }

    ///Devide query k-mer list into blocks for multi threading.
    vector<QueryKmerSplit> splits;
    int threadNum = ThreadNum;
    size_t querySplitSize = queryKmerCnt / (threadNum - 1);
    uint64_t queryKmerAA;
    bool splitCheck = false;
    splits.emplace_back(0, querySplitSize - 1, querySplitSize, 0, 0, 0);
    for(int i = 1; i < threadNum; i++) {
        queryKmerAA = AminoAcid(queryKmerList[querySplitSize * i].ADkmer);
        splitCheck = false;
        for(size_t j = 0; j < numOfDiffIdxSplits; j++){
            if(queryKmerAA < AminoAcid(diffIdxSplits.data[j].ADkmer)){
                if(i == threadNum - 1)
                    splits.emplace_back(querySplitSize * i, queryKmerCnt - 1, querySplitSize, diffIdxSplits.data[numOfDiffIdxSplits_use - 1].ADkmer,
                                        diffIdxSplits.data[numOfDiffIdxSplits_use - 1].diffIdxOffset, diffIdxSplits.data[numOfDiffIdxSplits_use - 1].infoIdxOffset);
                else
                    splits.emplace_back(querySplitSize * i, querySplitSize * (i + 1) - 1, querySplitSize, diffIdxSplits.data[j - 1].ADkmer,
                                        diffIdxSplits.data[j - 1].diffIdxOffset,diffIdxSplits.data[j - 1].infoIdxOffset);
                break;
            }
        }
    }
    bool * splitCheckList = (bool *)malloc(sizeof(bool)*threadNum);
    fill_n(splitCheckList, threadNum, false);
    int completedSplitCnt = 0;

    ///taxonomical ID at the lowest rank? or at the rank of redundancy reduced
    vector<const vector<int> *> taxID;
    taxID.push_back(& taxIdList);
    taxID.push_back(& taxIdListAtRank);


    size_t numOfcall = 0;
    size_t queryIdx=0;
    size_t numOfTargetKmer = targetInfoList.fileSize / sizeof(TargetKmerInfo);
    omp_set_num_threads(ThreadNum);
    while( completedSplitCnt < threadNum) {
        bool hasOverflow = false;
#pragma omp parallel default(none), shared(queryIdx, completedSplitCnt, splitCheckList, numOfTargetKmer, hasOverflow, numOfcall, splits, queryKmerList, targetDiffIdxList, targetInfoList, matchBuffer, taxID, cout)
        {
            ///query variables
            uint64_t currentQuery = UINT64_MAX;
            uint64_t currentQueryAA = UINT64_MAX;

            ///target variables
            size_t diffIdxPos = 0;
            size_t targetInfoIdx = 0;
            vector<uint64_t> targetKmerCache; //vector for candidate target k-mer, some of which are selected after based on hamming distance
            uint64_t currentTargetKmer;

            ///vectors for selected target k-mers
            vector<uint8_t> selectedHammings;
            vector<size_t> selectedMatches;

            size_t startIdxOfAAmatch = 0;
            size_t posToWrite;
            size_t range;
#pragma omp for schedule(dynamic, 1)
            for (size_t i = 0; i < splits.size(); i++){
                if(hasOverflow || splitCheckList[i])
                    continue;

                targetInfoIdx = splits[i].diffIdxSplit.infoIdxOffset;
                diffIdxPos = splits[i].diffIdxSplit.diffIdxOffset + 1; //overflow 되었을 당시의 값들을 저장하여, target split의 처음부터 다시 시작하는 걸 방지할 수 있음             targetInfoIdx = splits[i].diffIdxSplit.infoIdxOffset;
                currentTargetKmer = splits[i].diffIdxSplit.ADkmer;
                currentQuery = UINT64_MAX;
                currentQueryAA = UINT64_MAX;

                for(size_t j = splits[i].start; j < splits[i].end + 1; j ++){
                    splits[i].start++;
                    ///Reuse the comparison data if queries are exactly identical
                    if(currentQuery == queryKmerList[j].ADkmer){
                        posToWrite = matchBuffer.reserveMemory(selectedMatches.size());
                        if(posToWrite + selectedMatches.size() >= matchBuffer.bufferSize){
                            hasOverflow = true;
                            splits[i].start = j;
#pragma omp atomic
                            matchBuffer.startIndexOfReserve -= selectedMatches.size();
                            break;
                        } else{
                            range = selectedMatches.size();
                            for (size_t k = 0; k < range; k++) {
                                matchBuffer.buffer[posToWrite] = {queryKmerList[j].info.sequenceID, taxID[targetInfoList.data[selectedMatches[k]].redundancy]->at(targetInfoList.data[selectedMatches[k]].sequenceID),
                                                                  queryKmerList[j].info.pos, queryKmerList[j].info.frame, selectedHammings[k]};
                                posToWrite ++;
                            }
                        }
                        continue;
                    }
                    selectedMatches.clear();
                    selectedHammings.clear();

                    ///Reuse the candidate target k-mers to compare in DNA level if queries are the same only at amino acid level
                    if(currentQueryAA == AminoAcid(queryKmerList[j].ADkmer)){
                        compareDna(queryKmerList[j].ADkmer, targetKmerCache, startIdxOfAAmatch, selectedMatches, selectedHammings);
                        posToWrite = matchBuffer.reserveMemory(selectedMatches.size());
                        if(posToWrite + selectedMatches.size() >= matchBuffer.bufferSize){
                            hasOverflow = true;
                            splits[i].start = j;
#pragma omp atomic
                            matchBuffer.startIndexOfReserve -= selectedMatches.size();
                            break;
                        } else{
                            range = selectedMatches.size();
                            for (size_t k = 0; k < range; k++) {
                                matchBuffer.buffer[posToWrite] = {queryKmerList[j].info.sequenceID, taxID[targetInfoList.data[selectedMatches[k]].redundancy]->at(targetInfoList.data[selectedMatches[k]].sequenceID),
                                                                  queryKmerList[j].info.pos, queryKmerList[j].info.frame, selectedHammings[k]};
                                posToWrite ++;
                            }
                        }
                        continue;
                    }
                    targetKmerCache.clear();

                    ///Get next query and start to find
                    currentQuery = queryKmerList[j].ADkmer;
                    currentQueryAA = AminoAcid(currentQuery);

                    ///Skip target k-mers that are not matched in amino acid level
                    while (AminoAcid(currentQuery) > AminoAcid(currentTargetKmer) && (targetInfoIdx < numOfTargetKmer)) {
                        currentTargetKmer = getNextTargetKmer(currentTargetKmer, targetDiffIdxList.data, diffIdxPos);
                        targetInfoIdx++;
                    }

                    startIdxOfAAmatch = targetInfoIdx;
                    ///Load target k-mers thatare matched in amino acid level
                    while (AminoAcid(currentQuery) == AminoAcid(currentTargetKmer) && (targetInfoIdx < numOfTargetKmer)) {
                        targetKmerCache.push_back(currentTargetKmer);
                        currentTargetKmer = getNextTargetKmer(currentTargetKmer, targetDiffIdxList.data, diffIdxPos);
                        targetInfoIdx++;
                    }

                    ///Compare the current query and the loaded target k-mers and select
                    compareDna(currentQuery, targetKmerCache, startIdxOfAAmatch, selectedMatches, selectedHammings);
                    posToWrite = matchBuffer.reserveMemory(selectedMatches.size());
                    if(posToWrite + selectedMatches.size() >= matchBuffer.bufferSize){
                        hasOverflow = true;
                        splits[i].start = j;
#pragma omp atomic
                        matchBuffer.startIndexOfReserve -= selectedMatches.size();
                        break;
                    } else{
                        range = selectedMatches.size();
                        for (size_t k = 0; k < range; k++) {
                            matchBuffer.buffer[posToWrite] = {queryKmerList[j].info.sequenceID, taxID[targetInfoList.data[selectedMatches[k]].redundancy]->at(targetInfoList.data[selectedMatches[k]].sequenceID),
                                                              queryKmerList[j].info.pos, queryKmerList[j].info.frame, selectedHammings[k]};
                            posToWrite ++;
                        }
                    }

                }

               ///Check whether current split is completed or not
                if(splits[i].start - 1 == splits[i].end){
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

void Classifier::compareDna(uint64_t & query, vector<uint64_t> & targetList, const size_t & startIdx, vector<size_t> & selectedMatches, vector<uint8_t> & selectedHamming) {
    vector<uint8_t> hammings;
    uint8_t currentHamming;
    uint8_t minHamming = UINT8_MAX;
    ///Calculate hamming distance
    for(size_t i = 0; i < targetList.size(); i++){
        currentHamming = getHammingDistance(query, targetList[i]);
        if(currentHamming < minHamming)
            minHamming = currentHamming;
        hammings.push_back(currentHamming);
    }

    if(minHamming > 2) return;

    ///Select target k-mers that passed hamming criteria
    for(size_t h = 0; h < hammings.size(); h++){
        if(hammings[h] == minHamming){
            selectedMatches.push_back(startIdx + h);
            selectedHamming.push_back(hammings[h]);
        }
    }
}

///It analyses the result of linear search.
void Classifier::analyseResult(NcbiTaxonomy & ncbiTaxonomy, vector<Sequence> & seqSegments, char * matchFileName, Query * queryList){
    struct MmapedData<Match> matchList = mmapData<Match>(matchFileName);
    size_t numOfMatches = matchList.fileSize / sizeof(Match);
    SORT_PARALLEL(matchList.data, matchList.data + numOfMatches , Classifier::compareForWritingMatches);

    uint32_t currentQuery;
    size_t i = 0;
    size_t queryOffset;
    size_t queryEnd;
    while(i < numOfMatches) {
        currentQuery = matchList.data[i].queryId;
        queryOffset = i;
        while((currentQuery == matchList.data[i].queryId) && (i < numOfMatches)) i++;
        queryEnd = i - 1;
        TaxID selectedLCA = chooseBestTaxon(ncbiTaxonomy, seqSegments[currentQuery].length, currentQuery, queryOffset,
                                            queryEnd, matchList.data);
        ++taxCounts[selectedLCA];
    }
    munmap(matchList.data, matchList.fileSize + 1);
}
void Classifier::analyseResultParallel(NcbiTaxonomy & ncbiTaxonomy, vector<Sequence> & seqSegments, char * matchFileName, int seqNum){
    struct MmapedData<Match> matchList = mmapData<Match>(matchFileName);
    size_t numOfMatches = matchList.fileSize / sizeof(Match);
    SORT_PARALLEL(matchList.data, matchList.data + numOfMatches , Classifier::compareForWritingMatches);
    cout<<"num of matches"<<numOfMatches<<endl;
    ///Get match blocks for multi threading
    typedef Sequence Block;
    Block * matchBlocks = new Block[seqNum];
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

    cout<<"here"<<endl;
    omp_set_num_threads(ThreadNum);
#pragma omp parallel default(none), shared(matchBlocks, matchList, seqSegments, seqNum, ncbiTaxonomy)
{
#pragma omp for schedule(dynamic, 1)
    for(size_t i = 0; i < seqNum; ++ i ){
        TaxID selectedLCA = chooseBestTaxon(ncbiTaxonomy, seqSegments[i].length, i, matchBlocks[i].start,
                                            matchBlocks[i].end, matchList.data);
#pragma omp atomic
        ++taxCounts[selectedLCA];
    }
}
    delete[] matchBlocks;
    munmap(matchList.data, matchList.fileSize + 1);
}

///For a query read, assign the best Taxon, using k-mer matches
TaxID Classifier::chooseBestTaxon(NcbiTaxonomy & ncbiTaxonomy, const size_t & queryLength, const int & currentQuery, const size_t & offset, const size_t & end, Match * matchList){
    vector<ConsecutiveMatches> coMatches;

    float coverageThr = 0.3;
    int conCnt = 0;
    uint32_t gapCnt = 0;
    uint32_t hammingSum = 0;
    uint32_t conBegin = 0;
    uint32_t conEnd = 0;
    size_t beginIdx = 0;
    size_t endIdx = 0;

    size_t i = offset;

    ///This routine is for getting consecutive matched k-mer
    ///gapThr decides the maximun gap
    int currentFrame;
    int gapThr = 0;
    while(i < end) {
        currentFrame = matchList[i].frame;
        while ((matchList[i + 1].frame == matchList[i].frame) && (i < end)) {
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
                    coMatches.emplace_back(conBegin, conEnd, hammingSum, gapCnt, beginIdx, endIdx);
//                    cout << currentFrame << " " << conBegin << " " << conEnd << " " << conCnt << " " << gapCnt << " "
//                         << hammingSum << endl;
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
            coMatches.emplace_back(conBegin, conEnd, hammingSum, gapCnt, beginIdx, endIdx);
//            cout << currentFrame << " " << conBegin << " " << conEnd << " " << conCnt << " " << gapCnt << " "
//                 << hammingSum << endl;
            conCnt = 0;
            gapCnt = 0;
            hammingSum = 0;
        }
        i++;
    }

    if (coMatches.size() == 0) return 0;
    sort(coMatches.begin(), coMatches.end(), Classifier::compareConsecutiveMatches);

    ///Align consecutive matches back to query.
    vector<ConsecutiveMatches> alignedCoMatches;

    alignedCoMatches.push_back(coMatches[0]);
    auto alignedBegin = alignedCoMatches.begin();
    int isOverlaped= 0;
    int overlappedIdx = 0;
    for(size_t i2 = 1; i2 < coMatches.size(); i2++){
        isOverlaped = 0;
        overlappedIdx = 0;
        for(size_t j = 0; j < alignedCoMatches.size(); j++){
            if((alignedCoMatches[j].begin < coMatches[i2].end) && (alignedCoMatches[j].end > coMatches[i2].begin)){ ///TODO check this condition
                isOverlaped = 1;
                overlappedIdx = j;
                break;
            }
        }

        if(1 == isOverlaped){ ///TODO what to do when overlaps
            continue;
        } else{
            alignedCoMatches.push_back(coMatches[i2]);
        }
    }

    ///Check a query coverage
    int maxNum = queryLength / 3 - kmerLength + 1;
    int matchedNum = 0;
    int coveredLen = 0;
    float coverage;
    for(size_t cm = 0 ; cm < alignedCoMatches.size(); cm ++){
        matchedNum += (alignedCoMatches[cm].end - alignedCoMatches[cm].begin)/3 + 1;
        coveredLen += alignedCoMatches[cm].end - alignedCoMatches[cm].begin + 24;
    }
    coverage = float(matchedNum) / float(maxNum);



    ///TODO: how about considering hamming distance here?
    ///Get a lowest common ancestor, and check whether strain taxIDs are existing
    vector<TaxID> taxIdList;
    TaxID temp;
    auto currentInfo = find(queryInfos.begin(), queryInfos.end(), currentQuery);
    for(size_t cs = 0; cs < alignedCoMatches.size(); cs++ ){
        for(size_t k = alignedCoMatches[cs].beginIdx ; k < alignedCoMatches[cs].endIdx + 1; k++ ){
            temp = matchList[k].taxID;
            taxIdList.push_back(temp);
            currentInfo->taxCnt[temp] ++;
        }
    }

    ///No classification for low coverage.
    if(coverage < coverageThr){
        currentInfo->coverage = coverage;
        //queryInfo[currentQuery].coverage = coverage;
        return 0;
    }

    size_t numAssignedSeqs = 0;
    size_t numUnassignedSeqs = 0;
    size_t numSeqsAgreeWithSelectedTaxon = 0;
    double selectedPercent = 0;

    TaxID selectedLCA = match2LCA(taxIdList, ncbiTaxonomy, 0.8, numAssignedSeqs,
                                  numUnassignedSeqs, numSeqsAgreeWithSelectedTaxon,
                                  selectedPercent);
//    cout<<"#"<<currentQuery<<endl;
//    cout<<"coverage: "<<coverage<<endl;
//    for(size_t k = 0; k < taxIdList.size(); k++){
//        cout<<taxIdList[k]<<endl;
//    }

    ///TODO optimize strain specific classification criteria
    ///Strain classification only for high coverage with LCA of species level
    if(coverage > 0.90 && NcbiTaxonomy::findRankIndex(ncbiTaxonomy.taxonNode(selectedLCA)->rank) == 4){ /// There are more strain level classifications with lower coverage threshold, but also with more false postives. 0.8~0.85 looks good.
        int strainCnt = 0;
        unordered_map<TaxID, int> strainMatchCnt;
        TaxID strainTaxId;

        for(size_t cs = 0; cs < alignedCoMatches.size(); cs++ ){
            for(size_t k = alignedCoMatches[cs].beginIdx ; k < alignedCoMatches[cs].endIdx + 1; k++ ){
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
//    cout<<"label: "<<selectedLCA<<endl<<endl;
//    if(NcbiTaxonomy::findRankIndex(ncbiTaxonomy.taxonNode(selectedLCA)->rank) == 3){
//        cout<<"strain level classification: "<<selectedLCA<<endl;
//    }else {
//        cout<<selectedLCA<<" "<<selectedPercent<<endl;
//    }

    ///store classification results
    currentInfo->isClassified = true;
    currentInfo->taxId = selectedLCA;
    currentInfo->coverage = coverage;
    return selectedLCA;
}



TaxID Classifier::match2LCA(const std::vector<int> & taxIdList, NcbiTaxonomy const & taxonomy, const float majorityCutoff,
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

///It reads differential index and return "current + (next - current)", which is equal to next.
inline uint64_t Classifier::getNextTargetKmer2(uint64_t lookingTarget, const uint16_t* targetDiffIdxList, size_t & diffIdxPos){
    uint16_t fragment;
    uint64_t diffIn64bit = 0;
    //bit packing, SIMD,
    for(int i = 0; i < 5  ; i++){
        fragment = targetDiffIdxList[diffIdxPos];
        diffIdxPos++;
        if(fragment & (0x1u << 15))
        {
            fragment &= ~(1<<15u);
            diffIn64bit |= fragment;
            break;
        }
        diffIn64bit |= fragment;
        diffIn64bit <<= 15u;
    }
    return diffIn64bit + lookingTarget;
}

inline uint64_t Classifier::getNextTargetKmer(uint64_t lookingTarget, const uint16_t* targetDiffIdxList, size_t & diffIdxPos){
    uint16_t fragment;
    uint16_t check = (0x1u << 15u);
    uint64_t diffIn64bit = 0;

    fragment = targetDiffIdxList[diffIdxPos];
    diffIdxPos++;
    while (!(fragment & check)){
        diffIn64bit |= fragment;
        diffIn64bit <<= 15u;
        fragment = targetDiffIdxList[diffIdxPos];
        diffIdxPos++;
    }
    fragment &= ~check;
    diffIn64bit |= fragment;

    return diffIn64bit + lookingTarget;
}

///
inline uint8_t Classifier::getHammingDistance(uint64_t kmer1, uint64_t kmer2) {
    uint8_t hammingDist = 0;
    for(int i = 0; i < 8 ; i++){
        hammingDist += hammingLookup[GET_3_BITS(kmer1)][GET_3_BITS(kmer2)];
        kmer1 >>= 3U;
        kmer2 >>= 3U;
    }
    return hammingDist;
}



int Classifier::getNumOfSplits() const {
    return this->numOfSplit;
}

bool Classifier::compareForLinearSearch(const QueryKmer & a, const QueryKmer & b){
    if(a.ADkmer < b.ADkmer){
        return true;
    } else if(a.ADkmer == b.ADkmer){
        return (a.info.sequenceID < b.info.sequenceID);
    }
    return false;
}

bool Classifier::compareConsecutiveMatches(const ConsecutiveMatches & a, const ConsecutiveMatches & b){
    if((a.end - a.begin) > (b.end- b.begin)){
        if((a.end - a.begin) == (b.end - b.begin + 1)){
            return (a.endIdx - a.beginIdx + 1) * 2 / ((a.hamming+1)*(a.gapCnt+1)) > (b.endIdx - b.beginIdx +1) * 2 / ((b.hamming + 1) * (b.gapCnt + 1));
        }
        return true;
    }else if((a.end - a.begin) == (b.end- b.begin)){
            return (a.endIdx - a.beginIdx + 1) * 2 / ((a.hamming+1)*(a.gapCnt+1)) > (b.endIdx - b.beginIdx +1) * 2 / ((b.hamming + 1) * (b.gapCnt + 1));
    }
    return false;
}

void Classifier::writeReadClassification(vector<QueryInfo> & queryInfos, ofstream & readClassificationFile){
    for(size_t i = 0; i < queryInfos.size(); i++){
        readClassificationFile << queryInfos[i].queryId << "\t" << queryInfos[i].isClassified << "\t" << queryInfos[i].name << "\t" << queryInfos[i].taxId << "\t" << queryInfos[i].queryLength << "\t" << queryInfos[i].coverage << "\t";
        for(auto it = queryInfos[i].taxCnt.begin(); it != queryInfos[i].taxCnt.end(); ++it){
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

void Classifier::performanceTest(NcbiTaxonomy & ncbiTaxonomy){

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

    for(size_t i = 0; i < queryInfos.size(); i++) {
        classificationResult = queryInfos[i].taxId;
        if (classificationResult == 0) {
            continue;
        } else {
            classifiedCnt ++;
            queryName = queryInfos[i].name;
            regex_search(queryName, assacc, regex1);
            if (assacc2taxid.count(assacc[0].str())) {
                rightAnswer = assacc2taxid[assacc[0].str()];
            } else {
                cout << assacc[0].str() << " is not in the mapping file" << endl;
                continue;
            }
            //cout<<"compareTaxon"<<" "<<i<<endl;
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

    if(NcbiTaxonomy::findRankIndex(shotRank) <= 3){
        //cout<<"subspecies"<<endl;
        if(shot == target){
            subspCnt ++;
        }
    } else if(shotRank == "species") {
        //cout<<"species"<<endl;
        if(shot == ncbiTaxonomy.getTaxIdAtRank(target, "species")){
            speciesCnt ++;
        }
    } else if(shotRank == "genus"){
        //cout<<"genus"<<endl;
        if(shot == ncbiTaxonomy.getTaxIdAtRank(target, "genus")){
            genusCnt ++;
        }
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
