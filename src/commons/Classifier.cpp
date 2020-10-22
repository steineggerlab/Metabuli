//
// Created by KJB on 01/09/2020.
//

#include "Classifier.h"
#include "FastSort.h"
Classifier::Classifier() : queryCount(0), multipleMatchCount(0), totalMatchCount(0), perfectMatchCount(0) {
    seqAlterator = new SeqAlterator();
    numOfSplit = 0;
    closestCount = 0;
    queryCount = 0;
    totalMatchCount = 0;
    perfectMatchCount = 0;
    ESP = {0 ,0};
}

Classifier::~Classifier() { delete seqAlterator; }

void Classifier::startClassify(const char * queryFileName, const char * targetDiffIdxFileName, const char * targetInfoFileName, const vector<int> & taxIdList) {
    string dnaBuffer;
    size_t bufferIdx = 0;

    struct MmapedData<char> queryFile = mmapData<char>(queryFileName);
    struct MmapedData<uint16_t> targetDiffIdxList = mmapData<uint16_t>(targetDiffIdxFileName);
    targetDiffIdxList.data[targetDiffIdxList.fileSize/sizeof(uint16_t)] = 32768; //1000000000000000
    struct MmapedData<KmerInfo> targetInfoList = mmapData<KmerInfo>(targetInfoFileName);

    vector<SeqSegment> seqSegments;
    seqAlterator->getSeqSegments(seqSegments, queryFile);
    Kmer * kmerBuffer = (Kmer *)malloc(sizeof(Kmer) * kmerBufSize);

    int seqID = 0;
    for(size_t i = 1 ; i < seqSegments.size(); i++)
    {
        seqAlterator-> dna2aa2(seqSegments[i], queryFile);
        ESP = seqAlterator->fillKmerBuffer2(seqSegments[i], queryFile, kmerBuffer, seqID, bufferIdx, ESP);
        while (ESP.startOfFrame + ESP.frame != 0)
        {
            linearSearch(kmerBuffer, bufferIdx, targetDiffIdxList, targetInfoList, taxIdList);
            writeResultFile(matchedKmerList,queryFileName);
            ESP = seqAlterator->fillKmerBuffer2(seqSegments[i], queryFile, kmerBuffer, seqID, bufferIdx, ESP);
        }
        seqID ++;
    }
    //compare the rest query k-mers with target k-mers
    linearSearch(kmerBuffer, bufferIdx, targetDiffIdxList, targetInfoList, taxIdList);
    writeResultFile(matchedKmerList,queryFileName);

    cout<<"query count                          : "<<queryCount<<endl;
    cout<<"Total match count                    : "<<totalMatchCount <<endl;
    cout<<"mutipleMatch in AA level             : "<<multipleMatchCount << endl;
    cout<<"matches in DNA level                 : "<<perfectMatchCount<<endl;
    cout<<"number of closest matches            : "<<closestCount<<endl;

    free(kmerBuffer);
    munmap(queryFile.data, queryFile.fileSize + 1);
    munmap(targetDiffIdxList.data, targetDiffIdxList.fileSize + 1);
    munmap(targetInfoList.data, targetInfoList.fileSize + 1);
}

void Classifier::linearSearch(Kmer * queryKmerList, size_t & bufferIdx, const MmapedData<uint16_t> & targetDiffIdxList, const MmapedData<KmerInfo> & targetInfoList, const vector<int> & taxIdList) {
    cout<<"compare started"<<endl;

    //initialize
    size_t diffIdxPos = 0;
    uint64_t lastFirstMatch = 0;
    long lastFirstDiffIdxPos = 0;
    int lastFirstTargetIdx = 0;

    size_t maxTargetSize = targetInfoList.fileSize / sizeof(KmerInfo);

    uint8_t lowestHamming;

    SORT_PARALLEL(queryKmerList, queryKmerList + bufferIdx , [=](Kmer x, Kmer y) { return x.ADkmer < y.ADkmer; });

    uint64_t nextTargetKmer = getNextTargetKmer(0, targetDiffIdxList.data, diffIdxPos);
    size_t tarIter = 0;

    uint64_t currentQuery = UINT64_MAX;
    uint64_t currentTargetKmer = UINT64_MAX;
    uint64_t currentQueryAA;
//    size_t closestKmersCnt = 0;
    for(size_t i = 0; i < bufferIdx; i++)
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
        while(tarIter < maxTargetSize  && AminoAcid(nextTargetKmer) <= currentQueryAA){
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
    bufferIdx = 0;
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
//    sort(matchList.begin(), matchList.end(), [=](MatchedKmer x, MatchedKmer y) { return x.queryID < y.queryID; });
//    cout<<suffixedResultFileName<<endl;
    FILE * fp = fopen(suffixedResultFileName,"ab");
    fwrite(&(matchList[0]), sizeof(MatchedKmer), matchList.size(), fp);
    fclose(fp);
    matchList.clear();
}

int Classifier::getNumOfSplits() const {
    return this->numOfSplit;
}

void Classifier::analyseResult(const char * queryFileName, NcbiTaxonomy & ncbiTaxonomy) {
    char suffixedResultFileName[1000];
    sprintf(suffixedResultFileName,"%s_result", queryFileName);
    struct MmapedData<MatchedKmer> resultFile = mmapData<MatchedKmer>(suffixedResultFileName);
    size_t numOfMatches = resultFile.fileSize / sizeof(MatchedKmer);
    SORT_PARALLEL(resultFile.data, resultFile.data + numOfMatches, [](const MatchedKmer & a, const MatchedKmer & b) {
        return a.queryID < b.queryID || (a.queryID == b.queryID && a.posInfo < b.posInfo);});
    size_t currentRead;
    size_t currentKmer;
    vector<int> matchedKmers;
    vector<int> matchedLCAs;
    vector<int> matchedTaxons;

//    for(size_t i = 0; i < numOfMatches; i++){

    size_t numAssignedSeqs = 0;
    size_t numUnassignedSeqs = 0;
    size_t numSeqsAgreeWithSelectedTaxon = 0;
    double selectedPercent = 0;

    size_t i = 0;
    while(i < numOfMatches){
        currentRead = resultFile.data[i].queryID;
        while(currentRead == (size_t)resultFile.data[i].queryID){
            currentKmer = resultFile.data[i].posInfo;
            while(currentKmer == (size_t)resultFile.data[i].posInfo){
                matchedKmers.push_back(i);
                i++;
            }
            TaxID selectedLCA = selectLcaFromTaxIdList(matchedKmers,ncbiTaxonomy, 0.8, numAssignedSeqs, numUnassignedSeqs, numSeqsAgreeWithSelectedTaxon, selectedPercent);
            matchedLCAs.push_back(selectedLCA);
            matchedKmers.clear();
        }
       ///TODO: giveTaxonToRead(matchedLCA)
       ///TODO: matchedLCA.clear()

    }

}

//TaxID Classifier::selectTaxForSet(const std::vector<int> &setTaxa, NcbiTaxonomy const *taxonomy, const float majorityCutoff,
//                       size_t &numAssignedSeqs, size_t &numUnassignedSeqs, size_t &numSeqsAgreeWithSelectedTaxon, double &selectedPercent) {
//    // count num occurences of each ancestor, possibly weighted
//    std::map<TaxID,taxNode> ancTaxIdsCounts;
//
//    // initialize counters and weights
//    numAssignedSeqs = 0;
//    numUnassignedSeqs = 0;
//    numSeqsAgreeWithSelectedTaxon = 0;
//    selectedPercent = 0;
//    double totalAssignedSeqsWeights = 0.0;
//
//    for (size_t i = 0; i < setTaxa.size(); ++i) {
//        TaxID currTaxId = setTaxa[i].taxon;
//        double currWeight = setTaxa[i].weight;
//        // ignore unassigned sequences
//        if (currTaxId == 0) {
//            numUnassignedSeqs++;
//            continue;
//        }
//        TaxonNode const * node = taxonomy->taxonNode(currTaxId, false);
//        if (node == NULL) {
//            Debug(Debug::ERROR) << "taxonid: " << currTaxId << " does not match a legal taxonomy node.\n";
//            EXIT(EXIT_FAILURE);
//        }
//        totalAssignedSeqsWeights += currWeight;
//        numAssignedSeqs++;
//
//        // each start of a path due to an orf is a candidate
//        if (ancTaxIdsCounts.find(currTaxId) != ancTaxIdsCounts.end()) { //원소가 있다면
//            ancTaxIdsCounts[currTaxId].update(currWeight, 0);
//        } else {
//            taxNode currNode;
//            currNode.set(currWeight, true, 0);
//            ancTaxIdsCounts.insert(std::pair<TaxID,taxNode>(currTaxId, currNode));
//        }
//
//        // iterate all ancestors up to root (including). add currWeight and candidate status to each
//        TaxID currParentTaxId = node->parentTaxId;
//        while (currParentTaxId != currTaxId) {
//            if (ancTaxIdsCounts.find(currParentTaxId) != ancTaxIdsCounts.end()) {
//                ancTaxIdsCounts[currParentTaxId].update(currWeight, currTaxId);
//            } else {
//                taxNode currParentNode;
//                currParentNode.set(currWeight, false, currTaxId);
//                ancTaxIdsCounts.insert(std::pair<TaxID,taxNode>(currParentTaxId, currParentNode));
//            }
//            // move up:
//            currTaxId = currParentTaxId;
//            node = taxonomy->taxonNode(currParentTaxId, false);
//            currParentTaxId = node->parentTaxId;
//        }
//    }
//
//    // select the lowest ancestor that meets the cutoff
//    int minRank = INT_MAX;
//    TaxID selctedTaxon = 0;
//
//    for (std::map<TaxID,taxNode>::iterator it = ancTaxIdsCounts.begin(); it != ancTaxIdsCounts.end(); it++) {
//        // consider only candidates:
//        if (!(it->second.isCandidate)) {
//            continue;
//        }
//
//        double currPercent = float(it->second.weight) / totalAssignedSeqsWeights;
//        if (currPercent >= majorityCutoff) {
//            // iterate all ancestors to find lineage min rank (the candidate is a descendant of a node with this rank)
//            TaxID currTaxId = it->first;
//            TaxonNode const * node = taxonomy->taxonNode(currTaxId, false);
//            int currMinRank = ROOT_RANK;
//            TaxID currParentTaxId = node->parentTaxId;
//            while (currParentTaxId != currTaxId) {
//                int currRankInd = NcbiTaxonomy::findRankIndex(node->rank);
//                if ((currRankInd > 0) && (currRankInd < currMinRank)) {
//                    currMinRank = currRankInd;
//                    // the rank can only go up on the way to the root, so we can break
//                    break;
//                }
//                // move up:
//                currTaxId = currParentTaxId;
//                node = taxonomy->taxonNode(currParentTaxId, false);
//                currParentTaxId = node->parentTaxId;
//            }
//
//            if ((currMinRank < minRank) || ((currMinRank == minRank) && (currPercent > selectedPercent))) {
//                selctedTaxon = it->first;
//                minRank = currMinRank;
//                selectedPercent = currPercent;
//            }
//        }
//    }
//
//    // count the number of seqs who have selectedTaxon in their ancestors (agree with selection):
//    if (selctedTaxon == ROOT_TAXID) {
//        // all agree with "root"
//        numSeqsAgreeWithSelectedTaxon = numAssignedSeqs;
//        return (selctedTaxon);
//    }
//    if (selctedTaxon == 0) {
//        // nothing informative
//        return (selctedTaxon);
//    }
//    // otherwise, iterate over all seqs
//    for (size_t i = 0; i < setTaxa.size(); ++i) {
//        TaxID currTaxId = setTaxa[i].taxon;
//        // ignore unassigned sequences
//        if (currTaxId == 0) {
//            continue;
//        }
//        TaxonNode const * node = taxonomy->taxonNode(currTaxId, false);
//
//        // iterate all ancestors up to the root
//        TaxID currParentTaxId = node->parentTaxId;
//        while (currParentTaxId != currTaxId) {
//            if (currTaxId == selctedTaxon) {
//                numSeqsAgreeWithSelectedTaxon++;
//                break;
//            }
//            currTaxId = currParentTaxId;
//            node = taxonomy->taxonNode(currParentTaxId, false);
//            currParentTaxId = node->parentTaxId;
//        }
//    }
//
//    return (selctedTaxon);
//}

