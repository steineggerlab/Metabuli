//
// Created by KJB on 01/09/2020.
//

#include "Classifier.h"

Classifier::Classifier() : queryCount(0), multipleMatchCount(0), totalMatchCount(0), perfectMatchCount(0) {
    seqIterator = new SeqIterator();
    numOfSplit = 0;
    closestCount = 0;
    queryCount = 0;
    ESP = {0 ,0};
}

Classifier::~Classifier() { delete seqIterator; }

void Classifier::startClassify(const char * queryFileName, const char * targetDiffIdxFileName, const char * targetInfoFileName, const vector<int> & taxIdList) {

    NcbiTaxonomy ncbiTaxonomy("/Users/jaebeomkim/Desktop/pjt/taxdmp/names.dmp",
                              "/Users/jaebeomkim/Desktop/pjt/taxdmp/nodes.dmp",
                              "/Users/jaebeomkim/Desktop/pjt/taxdmp/merged.dmp");

    string dnaBuffer;
    size_t bufferIdx = 0;

    struct MmapedData<char> queryFile = mmapData<char>(queryFileName);
    struct MmapedData<uint16_t> targetDiffIdxList = mmapData<uint16_t>(targetDiffIdxFileName);
    targetDiffIdxList.data[targetDiffIdxList.fileSize/sizeof(uint16_t)] = 32768; //1000000000000000
    struct MmapedData<KmerInfo> targetInfoList = mmapData<KmerInfo>(targetInfoFileName);

    vector<Sequence> sequences;
    seqIterator->getSeqSegmentsWithHead(sequences, queryFile);
    size_t numOfSeq = sequences.size();

    bool processedSeqChecker[numOfSeq];
    fill_n(processedSeqChecker, numOfSeq, false);

    KmerBuffer kmerBuffer(kmerBufSize);
    size_t processedSeqCnt;

    cout<<"numOfseq: "<<numOfSeq<<endl;
    while(processedSeqCnt < numOfSeq){ ///check this condition
        seqIterator->fillQueryKmerBufferParallel(kmerBuffer, queryFile, sequences, processedSeqChecker, processedSeqCnt);
        linearSearch(kmerBuffer.buffer, kmerBuffer.startIndexOfReserve, targetDiffIdxList, targetInfoList, taxIdList);
    }


    //compare the rest query k-mers with target k-mers
//    linearSearch(kmerBuffer.buffer, bufferIdx, targetDiffIdxList, targetInfoList, taxIdList);
//    writeResultFile(matchedKmerList,queryFileName);

    cout<<"query count                          : "<<queryCount<<endl;
    cout<<"Total match count                    : "<<totalMatchCount <<endl;
    cout<<"mutipleMatch in AA level             : "<<multipleMatchCount << endl;
    cout<<"matches in DNA level                 : "<<perfectMatchCount<<endl;
    cout<<"number of closest matches            : "<<closestCount<<endl;

    analyseResult2(ncbiTaxonomy, sequences);

    free(kmerBuffer.buffer);
    munmap(queryFile.data, queryFile.fileSize + 1);
    munmap(targetDiffIdxList.data, targetDiffIdxList.fileSize + 1);
    munmap(targetInfoList.data, targetInfoList.fileSize + 1);
}

void Classifier::linearSearch(Kmer * queryKmerList, size_t & numOfKmer, const MmapedData<uint16_t> & targetDiffIdxList, const MmapedData<KmerInfo> & targetInfoList, const vector<int> & taxIdList) {
    //initialize
    size_t diffIdxPos = 0;
    uint64_t lastFirstMatch = 0;
    long lastFirstDiffIdxPos = 0;
    int lastFirstTargetIdx = 0;

    size_t numOfTargetKmer = targetInfoList.fileSize / sizeof(KmerInfo);
    uint8_t lowestHamming;
    SORT_PARALLEL(queryKmerList, queryKmerList + numOfKmer , [=](Kmer x, Kmer y) { return x.ADkmer < y.ADkmer; });
    uint64_t nextTargetKmer = getNextTargetKmer(0, targetDiffIdxList.data, diffIdxPos);
    size_t tarIter = 0;

    uint64_t currentQuery = UINT64_MAX;
    uint64_t currentTargetKmer = UINT64_MAX;
    uint64_t currentQueryAA;

    for(size_t i = 0; i < numOfKmer; i++)
    {
        /// get next query
        if(AminoAcid(currentQuery) == AminoAcid(queryKmerList[i].ADkmer)){
            nextTargetKmer = lastFirstMatch;
            diffIdxPos = lastFirstDiffIdxPos;
            tarIter = lastFirstTargetIdx;
        }
        currentQuery = queryKmerList[i].ADkmer;
        isMatched = 0;
        lowestHamming = 100;
        queryCount ++;
        currentQueryAA = AminoAcid(currentQuery);
        while((tarIter < numOfTargetKmer) && (AminoAcid(nextTargetKmer) <= currentQueryAA)){
            currentTargetKmer = nextTargetKmer;
            currentTargetPos = diffIdxPos;
            nextTargetKmer = getNextTargetKmer(currentTargetKmer, targetDiffIdxList.data, diffIdxPos);
            perfectMatchCount += (currentQuery == currentTargetKmer);
            if(AminoAcid(currentTargetKmer) == AminoAcid(currentQuery)){
                if (isMatched == 0) {
                    lastFirstMatch = currentTargetKmer;
                    lastFirstDiffIdxPos = currentTargetPos;
                    lastFirstTargetIdx = tarIter;
                    isMatched = 1;
                }
                totalMatchCount++;

                currentHamming = getHammingDistance(currentQuery, currentTargetKmer);
                if(currentHamming > lowestHamming){
                    tarIter ++;
                    continue;
                }
                else if(currentHamming < lowestHamming){
                    closestKmers.clear();
                    lowestHamming = currentHamming;
                }
                if(currentHamming == lowestHamming) closestKmers.push_back(tarIter);

//                if(currentHamming > lowestHamming){
//                    tarIter ++;
//                    continue;
//                } else if(currentHamming < lowestHamming){
//                    closestKmersCnt = (currentHamming == lowestHamming) ? closestKmersCnt : 0;
//                    closestKmers[closestKmersCnt] = tarIter;
//                    lowestHamming = currentHamming;
//                    closestKmersCnt++;
//
//                }

            }
            tarIter ++;
        }

        for(size_t k = 0; k < closestKmers.size(); k++){
            matchedKmerList.emplace_back(queryKmerList[i].info.sequenceID, targetInfoList.data[closestKmers[k]].sequenceID, taxIdList[targetInfoList.data[closestKmers[k]].sequenceID],
                                         queryKmerList[i].info.pos, lowestHamming, targetInfoList.data[closestKmers[k]].redundancy);
            closestCount++;
        }
        closestKmers.clear();
    }
    numOfKmer = 0;
    cout<<matchedKmerList.size()<<endl;
}

uint64_t Classifier::getNextTargetKmer(uint64_t lookingTarget, const uint16_t* targetDiffIdxList, size_t & diffIdxPos){
    uint16_t fragment;
    uint64_t diffIn64bit = 0;
    for(int i = 0; i < 5; i++)
    {
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

uint8_t Classifier::getHammingDistance(uint64_t kmer1, uint64_t kmer2) {
    uint8_t hammingDist = 0;
    for(int i = 0; i < 8 ; i++)
    {
        hammingDist += hammingLookup[GET_3_BITS(kmer1)][GET_3_BITS(kmer2)];
        kmer1 >>= 3U;
        kmer2 >>= 3U;
    }
    return hammingDist;
}

void Classifier::writeResultFile(vector<MatchedKmer> & matchList, const char * queryFileName) {
    char suffixedResultFileName[1000];
    sprintf(suffixedResultFileName,"%s_result", queryFileName);
    numOfSplit++;
    sort(matchList.begin(), matchList.end(), [=](MatchedKmer x, MatchedKmer y) { return x.queryID < y.queryID; });
    cout<<suffixedResultFileName<<endl;
    FILE * fp = fopen(suffixedResultFileName,"wb");
    fwrite(&(matchList[0]), sizeof(MatchedKmer), matchList.size(), fp);
    fclose(fp);
    matchList.clear();
}

int Classifier::getNumOfSplits() const {
    return this->numOfSplit;
}

void Classifier::analyseResult(const char * queryFileName, NcbiTaxonomy & ncbiTaxonomy, vector<Sequence> & seqSegments) {
    char suffixedResultFileName[1000];
    sprintf(suffixedResultFileName,"%s_result", queryFileName);
    string reduceLevel = "species";
    struct MmapedData<MatchedKmer> resultFile = mmapData<MatchedKmer>(suffixedResultFileName);
    size_t numOfMatches = resultFile.fileSize / sizeof(MatchedKmer);
    SORT_PARALLEL(resultFile.data, resultFile.data + numOfMatches, [](const MatchedKmer & a, const MatchedKmer & b) {
        return a.queryID < b.queryID || (a.queryID == b.queryID && a.queryPos < b.queryPos);});
    size_t currentReadId;
    size_t currentKmer;
    vector<int> matchedKmers;
//    vector<TaxID> matchedLCAs;
    vector<TaxID>::iterator matchedLCAsIt;
    vector<int> matchedLCACounts;
    unordered_map<TaxID, int> matchedLCAs; //int for counting
    unordered_map<int, TaxID> assignedReads; //map<queryID, TaxID>

//    for(size_t i = 0; i < numOfMatches; i++){

    size_t numAssignedSeqs = 0;
    size_t numUnassignedSeqs = 0;
    size_t numSeqsAgreeWithSelectedTaxon = 0;
    double selectedPercent = 0;
    cout<<"numOfMatches: "<<numOfMatches<<endl;
    size_t i = 0;
    while(i < numOfMatches) {
        currentReadId = resultFile.data[i].queryID;
        while((currentReadId == (size_t) resultFile.data[i].queryID) && (i < numOfMatches)){
            currentKmer = resultFile.data[i].queryPos;
            while((currentKmer == (size_t) resultFile.data[i].queryPos) && (i < numOfMatches)) {
//                cout<<currentKmer<<" "<<resultFile.data[i].queryPos<<" "<<resultFile.data[i].tragetID<<endl;
//                cout<<resultFile.data[i].taxID<<endl;
                if(resultFile.data[i].redundancy){
                    matchedKmers.push_back(ncbiTaxonomy.taxIdAtRank(resultFile.data[i].taxID, reduceLevel));
                }else {
                    matchedKmers.push_back(resultFile.data[i].taxID);
                }
                i++;
            }

            TaxID selectedLCA = selectLcaFromTaxIdList(matchedKmers, ncbiTaxonomy, 0.8, numAssignedSeqs,
                                                       numUnassignedSeqs, numSeqsAgreeWithSelectedTaxon,
                                                       selectedPercent);
            cout<<"selectedLCA: "<<selectedLCA<<endl;
            if (matchedLCAs.find(selectedLCA) == matchedLCAs.end()) {
                matchedLCAs.insert(pair<TaxID, int>(selectedLCA, 1));
            } else {
                matchedLCAs[selectedLCA]++;
            }
            matchedKmers.clear();
        }
//        cout<<"selected leaf:"<<selectALeaf(matchedLCAs, ncbiTaxonomy, seqSegments[currentReadId+1].length )<<endl;
        assignedReads.insert(pair<int, TaxID>(currentReadId, selectALeaf(matchedLCAs, ncbiTaxonomy,seqSegments[currentReadId].length)));
        cout<<currentReadId<<" "<<seqSegments[currentReadId].length<<endl;
        matchedLCAs.clear();
    }

    for(auto it: assignedReads){
        cout<<it.first<<" "<<it.second<<endl;
    }
}

void Classifier::analyseResult2(NcbiTaxonomy & ncbiTaxonomy, vector<Sequence> & seqSegments) {
    string reduceLevel = "species";
    size_t numOfMatches = matchedKmerList.size();
    SORT_PARALLEL(matchedKmerList.begin(), matchedKmerList.end(), [](const MatchedKmer & a, const MatchedKmer & b) {
        return a.queryID < b.queryID || (a.queryID == b.queryID && a.queryPos < b.queryPos);});
    int currentQuery;
    uint32_t currentKmer;
    vector<int> matchedKmers;
    unordered_map<TaxID, int> matchedLCAs; //int for counting
    unordered_map<int, TaxID> assignedReads; //map<queryID, TaxID>

    size_t numAssignedSeqs = 0;
    size_t numUnassignedSeqs = 0;
    size_t numSeqsAgreeWithSelectedTaxon = 0;
    double selectedPercent = 0;

    size_t i = 0;
    while(i < numOfMatches) {
        currentQuery = matchedKmerList[i].queryID;
        while((currentQuery ==  matchedKmerList[i].queryID) && (i < numOfMatches)){
            currentKmer = matchedKmerList[i].queryPos;
            while((currentKmer == matchedKmerList[i].queryPos) && (i < numOfMatches)) {
                if(matchedKmerList[i].redundancy){
                    matchedKmers.push_back(ncbiTaxonomy.taxIdAtRank(matchedKmerList[i].taxID, reduceLevel));
                }else {
                    matchedKmers.push_back(matchedKmerList[i].taxID);
                }
                i++;

            }
            TaxID selectedLCA = selectLcaFromTaxIdList(matchedKmers, ncbiTaxonomy, 0.8, numAssignedSeqs,
                                                       numUnassignedSeqs, numSeqsAgreeWithSelectedTaxon,
                                                       selectedPercent);

            if (matchedLCAs.find(selectedLCA) == matchedLCAs.end()) {
                matchedLCAs.insert(pair<TaxID, int>(selectedLCA, 1));
            } else {
                matchedLCAs[selectedLCA]++;
            }
            matchedKmers.clear();
        }
//        cout<<"selected leaf:"<<selectALeaf(matchedLCAs, ncbiTaxonomy, seqSegments[currentQuery+1].length )<<endl;
        assignedReads.insert(pair<int, TaxID>(currentQuery, selectALeaf(matchedLCAs, ncbiTaxonomy, seqSegments[currentQuery].length)));
        matchedLCAs.clear();
    }

    for(auto it: assignedReads){
        cout<<it.first<<" "<<it.second<<endl;
    }
}

TaxID Classifier::selectALeaf(unordered_map<TaxID, int> & listOfLCAs, NcbiTaxonomy & ncbiTaxonomy, const size_t & length, float coverageThr){
    float numberOfPossibleKmersFromThisRead = length/3 - 7;

//    if(length % 3 == 0){
//        numberOfPossibleKmersFromThisRead = (6 * (length/3) - 46);
//    }else if(length % 3 == 1){
//        numberOfPossibleKmersFromThisRead = (6 * (length/3) - 44);
//    }else if(length % 3 == 2){
//        numberOfPossibleKmersFromThisRead = (6 * (length/3) - 42);
//    }
    unordered_map<TaxID, float> taxonNodeCount; //int = count
    double totalCount = 0;
    TaxID currTaxId;
    float currCount;
    int highestRankSatisfyingThr = 0;
    int maxCount = 0;
    int currRank;
    TaxID selectedLeaf = 0;
    coverageThr = 0.1;
    float maxCoverage = 0;
    float currCoverage;

    for(unordered_map<TaxID, int>::iterator it = listOfLCAs.begin(); it != listOfLCAs.end(); ++it) {
        currTaxId = it->first;
        currCount = it->second;
        for (std::pair<TaxID, int> it2 : listOfLCAs) {
            if ((it2.first != currTaxId) && ncbiTaxonomy.IsAncestor(it2.first, currTaxId)) {
                currCount += it2.second;
            }
        }
        cout<<"count of leaf: "<<currTaxId<< " "<<currCount<<" "<<((float)currCount)/numberOfPossibleKmersFromThisRead<<endl;
        taxonNodeCount.insert(pair<TaxID, float>(currTaxId, currCount));

        totalCount += it->second;
    }

    for(auto it = taxonNodeCount.begin(); it != taxonNodeCount.end(); ++it ) {
        currCoverage = it->second / numberOfPossibleKmersFromThisRead;
        if (currCoverage >= coverageThr){
            currRank = NcbiTaxonomy::findRankIndex(ncbiTaxonomy.taxonNode(it->first)->rank);
//            cout<<"rank: "<<ncbiTaxonomy.taxonNode(it->first)->rank<<endl;
//            cout<<"currRank: "<<currRank<<endl;
            if((currRank > highestRankSatisfyingThr) || (currRank == highestRankSatisfyingThr && currCoverage > maxCoverage )) {
                highestRankSatisfyingThr = currRank;
                maxCount = it->second;
                selectedLeaf = it->first;
                cout<<"highestRank: "<<highestRankSatisfyingThr<<" "<<"maxCount: "<<maxCount<<" "<<"selected leaf: "<<selectedLeaf<<endl;
            }
        }
    }
    cout<<"LEAF: "<<selectedLeaf<<endl<<endl;
    return selectedLeaf; // 0 -> unclassified
}

TaxID Classifier::selectLcaFromTaxIdList(const std::vector<int> & taxIdList, NcbiTaxonomy const & taxonomy, const float majorityCutoff,
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