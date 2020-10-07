//
// Created by KJB on 01/09/2020.
//

#include "Classifier.h"
#include "FastSort.h"
Classifier::Classifier() : queryCount(0), multipleMatchCount(0), totalMatchCount(0), perfectMatchCount(0)
{
    seqAlterator = new SeqAlterator();
    numOfSplit = 0;
    closestCount = 0;
    queryCount = 0;
    totalMatchCount = 0;
    perfectMatchCount = 0;
    ESP = {0 ,0};
}

Classifier::~Classifier() { delete seqAlterator; }

//void Classifier::startSearch(char * queryFileName, char * targetDiffIdxFileName, char * targetInfoFileName)
//{
//    string dnaBuffer;
//    string reads[2];
//
//    size_t bufferIdx = 0;
//    ifstream queryFile(queryFileName);
//    struct MmapedData<uint16_t> targetDiffIdxList = mmapData<uint16_t>(targetDiffIdxFileName);
//    struct MmapedData<KmerInfo> targetInfoList = mmapData<KmerInfo>(targetInfoFileName);
//
//    char taxID[100];
//    vector<int> taxIdList;
//    FILE * taxIdFile = fopen("/Users/kjb/Desktop/ADclassifier/refseq/taxIDs", "r");
//    while(feof(taxIdFile) == 0)
//    {
//        fscanf(taxIdFile,"%s",taxID);
//        taxIdList.push_back(atoi(taxID));
//    }
//    fclose(taxIdFile);
//
//    Kmer * kmerBuffer = (Kmer *)malloc(sizeof(Kmer) * kmerBufSize);
//
//    getline(queryFile, dnaBuffer);
//    int seqID = 1;
//
//    while(queryFile)
//    {
//        getline(queryFile, dnaBuffer);
//        if(dnaBuffer[0] == '>'){
//            reads[1] = kmerExtractor->reverseCompliment(reads[0]);
//            kmerExtractor->dna2aa(reads[0], reads[1]);
//            ESP = kmerExtractor->fillKmerBuffer(reads, kmerBuffer, seqID, bufferIdx, ESP);
//            while (ESP.startOfFrame + ESP.frame != 0)
//            {
//                linearSearch(kmerBuffer, bufferIdx, targetDiffIdxList, targetInfoList, taxIdList);
//                writeResultFile(matchedKmerList,queryFileName);
//                ESP = kmerExtractor->fillKmerBuffer(reads, kmerBuffer, seqID, bufferIdx, ESP);
//            }
//            reads[0].clear();
//            seqID++;
//            continue;
//        }
//        reads[0].append(dnaBuffer);
//    }
//    // For last one
//    reads[1] = kmerExtractor->reverseCompliment(reads[0]);
//    kmerExtractor->dna2aa(reads[0], reads[1]);
//    ESP = kmerExtractor->fillKmerBuffer(reads, kmerBuffer, seqID, bufferIdx, ESP);
//    while (ESP.startOfFrame + ESP.frame != 0)
//    {
//        linearSearch(kmerBuffer, bufferIdx, targetDiffIdxList, targetInfoList, taxIdList);
//        writeResultFile(matchedKmerList,queryFileName);
//        ESP = kmerExtractor->fillKmerBuffer(reads, kmerBuffer, seqID, bufferIdx, ESP);
//    }
//
//    //compare the rest query k-mers with target k-mers
//    linearSearch(kmerBuffer, bufferIdx, targetDiffIdxList, targetInfoList, taxIdList);
//    writeResultFile(matchedKmerList,queryFileName);
//
//    cout<<"query count                          : "<<queryCount<<endl;
//    cout<<"Total match count                    : "<<totalMatchCount <<endl;
//    cout<<"mutipleMatch in AA level             : "<<multipleMatchCount << endl;
//    cout<<"matches in DNA level                 : "<<perfectMatchCount<<endl;
//    cout<<"number of closest matches            : "<<closestCount<<endl;
//
//    queryFile.close();
//    free(kmerBuffer);
//    munmap(targetDiffIdxList.data, targetDiffIdxList.fileSize + 1);
//    munmap(targetInfoList.data, targetInfoList.fileSize + 1);
//}

void Classifier::startClassify(const char * queryFileName, const char * targetDiffIdxFileName, const char * targetInfoFileName, const vector<int> & taxIdList)
{
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

void Classifier::linearSearch(Kmer * kmerBuffer, size_t & bufferIdx, const MmapedData<uint16_t> & targetDiffIdxList, const MmapedData<KmerInfo> & targetInfoList, const vector<int> & taxIdList) {
    cout<<"compare started"<<endl;

    //initialize
    size_t diffIdxPos = 0;
    uint64_t lastFirstMatch = 0;
    long lastFirstDiffIdxPos = 0;
    int lastFirstTargetIdx = 0;

    size_t maxTargetSize = targetInfoList.fileSize / sizeof(KmerInfo);

    uint8_t lowestHamming;

    SORT_PARALLEL(kmerBuffer, kmerBuffer + bufferIdx , [=](Kmer x, Kmer y) { return x.ADkmer < y.ADkmer; });

    uint64_t nextTargetKmer = getNextTargetKmer(0, targetDiffIdxList.data, diffIdxPos);
    size_t tarIter = 0;

    uint64_t currentQuery = UINT64_MAX;
    uint64_t currentTargetKmer = UINT64_MAX;
    uint64_t currentQueryAA;
//    size_t closestKmersCnt = 0;
    for(size_t i = 0; i < bufferIdx; i++)
    {
        /// get next query
        if(AminoAcid(currentQuery) == AminoAcid(kmerBuffer[i].ADkmer)){
            nextTargetKmer = lastFirstMatch;
            diffIdxPos = lastFirstDiffIdxPos;
            tarIter = lastFirstTargetIdx;
        }
        currentQuery = kmerBuffer[i].ADkmer;
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

        for(size_t k = 0; k < closestKmers.size(); k++)
        {
            matchedKmerList.emplace_back(kmerBuffer[i].info.sequenceID, targetInfoList.data[closestKmers[k]].sequenceID, taxIdList[targetInfoList.data[closestKmers[k]].sequenceID],
                                         kmerBuffer[i].info.pos - targetInfoList.data[closestKmers[k]].pos,
                                         lowestHamming, targetInfoList.data[closestKmers[k]].redundancy);
            closestCount++;
        }
        closestKmers.clear();
    }
    bufferIdx = 0;
}

uint64_t Classifier::getNextTargetKmer(uint64_t lookingTarget, const uint16_t* targetDiffIdxList, size_t & diffIdxPos)
{
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

uint8_t Classifier::getHammingDistance(uint64_t kmer1, uint64_t kmer2)
{
    uint8_t hammingDist = 0;
    for(int i = 0; i < 8 ; i++)
    {
        hammingDist += hammingLookup[GET_3_BITS(kmer1)][GET_3_BITS(kmer2)];
        kmer1 >>= 3U;
        kmer2 >>= 3U;
    }
    return hammingDist;
}
void Classifier::writeResultFile(vector<MatchedKmer> & matchList, const char * queryFileName)
{
    char suffixedResultFileName[1000];
    sprintf(suffixedResultFileName,"%s_result_%zu", queryFileName,numOfSplit);
    numOfSplit++;
    cout<<suffixedResultFileName<<endl;
    FILE * fp = fopen(suffixedResultFileName,"wb");
    cout<<matchList.size();
    fwrite(&(matchList[0]), sizeof(MatchedKmer), matchList.size(), fp);
    fclose(fp);
    matchList.clear();
}