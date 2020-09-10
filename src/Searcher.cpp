//
// Created by KJB on 01/09/2020.
//

#include "Searcher.h"

Searcher::Searcher(){ kmerExtractor = new KmerExtractor();}

void Searcher::startSearch(char * queryFileName, char * targetDiffIdxFileName, char * targetInfoFileName)
{
    string buffer;
    string forwardRead;
    string reverseComplimentRead;
    string reads[2];
    ExtractStartPoint ESP = {0 ,0};
    size_t bufferIdx = 0;
    ifstream queryFile(queryFileName);
    struct MmapedData<uint16_t> targetDiffIdxList = mmapData<uint16_t>(targetDiffIdxFileName);
    struct MmapedData<KmerInfo> targetInfoList = mmapData<KmerInfo>(targetInfoFileName);

    Kmer * kmerBuffer = (Kmer *)malloc(sizeof(Kmer) * kmerBufSize);

    getline(queryFile, buffer);
    int seqID = 1;

    /// mmap을 써보자
    while(queryFile)
    {
        getline(queryFile, buffer);
        if(buffer[0] == '>'){
            reverseComplimentRead = kmerExtractor->reverseCompliment(forwardRead);
            reads[0] = forwardRead; reads[1] = reverseComplimentRead;
            kmerExtractor->dna2aa(forwardRead, reverseComplimentRead);

            ESP = kmerExtractor->fillKmerBuffer(reads, kmerBuffer, seqID, bufferIdx, ESP);
            while (ESP.startOfFrame + ESP.frame != 0)
            {
                linearSearch(kmerBuffer, bufferIdx, targetDiffIdxList, targetInfoList);
                ESP = kmerExtractor->fillKmerBuffer(reads, kmerBuffer, seqID, bufferIdx, ESP);
            }
            forwardRead.clear();
            seqID++;
            continue;
        }
        forwardRead.append(buffer);
    }
    // For last one
    reverseComplimentRead = kmerExtractor->reverseCompliment(forwardRead);
    reads[0] = forwardRead; reads[1] = reverseComplimentRead;
    kmerExtractor->dna2aa(forwardRead, reverseComplimentRead);
    ESP = kmerExtractor->fillKmerBuffer(reads, kmerBuffer, seqID, bufferIdx, ESP);
    while (ESP.startOfFrame + ESP.frame != 0)
    {
        linearSearch(kmerBuffer, bufferIdx, targetDiffIdxList, targetInfoList);
        ESP = kmerExtractor->fillKmerBuffer(reads, kmerBuffer, seqID, bufferIdx, ESP);
    }
    //compare the rest query k-mers with target k-mers
    linearSearch(kmerBuffer, bufferIdx, targetDiffIdxList, targetInfoList);

    queryFile.close();
    free(kmerBuffer);
    munmap(targetDiffIdxList.data, targetDiffIdxList.fileSize + 1);
    munmap(targetInfoList.data, targetInfoList.fileSize + 1);
}
void Searcher::linearSearch(Kmer * kmerBuffer, size_t & bufferIdx, const MmapedData<uint16_t> & targetDiffIdxList, const MmapedData<KmerInfo> & targetInfoList) {

    cout<<"compare started"<<endl;
    uint64_t lookingQuery = 0;
    uint64_t lookingTarget = 0;
    size_t lookingTargetPos = 0;
    uint64_t nextTarget = 0;
    uint8_t hammingDistance = 0;
    int matchCount  = 0;
    int isMatched;
    int queryCount = 0;
    uint64_t marker= ~0 & ~16777215;
    uint64_t lastFirstMatch = 0;
    long lastFirstDiffIdxPos = 0;
    int multipleMatch = 0;
    size_t maxTarget = targetInfoList.fileSize / sizeof(KmerInfo);

    size_t diffIdxPos = 0;
    int lastFirstTargetIdx = 0;

    int perfect = 0;

    sort(kmerBuffer, kmerBuffer + bufferIdx , [=](Kmer x, Kmer y) { return x.ADkmer < y.ADkmer; });
//    cout<<"first after sort : "<<kmerBuffer[0].ADkmer<<endl;
//    cout<<"last  after sort : "<<kmerBuffer[bufferIdx-1].ADkmer<<endl;
    vector<matchedKmer> matchedKmerList;
    nextTarget = getNextTargetKmer(0, targetDiffIdxList.data, diffIdxPos);

    for(size_t i = 0; i < bufferIdx; i++)
    {
        /// get next query
        lookingQuery = kmerBuffer[i].ADkmer;
        isMatched=0;
        queryCount ++;

        for (int j = lastFirstTargetIdx; j < maxTarget - 1; j++)
        {

            lookingTarget = nextTarget;
            lookingTargetPos = diffIdxPos;
            nextTarget = getNextTargetKmer(lookingTarget, targetDiffIdxList.data, diffIdxPos);

            if((lookingTarget & marker) == (lookingQuery & marker)) {
                matchCount++;
                matchedKmer temp = {lookingTarget, lookingQuery, kmerBuffer[i].info.sequenceID,
                                    targetInfoList.data[j].sequenceID,
                                    kmerBuffer[i].info.pos - targetInfoList.data[j].pos,
                                    getHammingDistance(lookingQuery, lookingTarget)};
                matchedKmerList.push_back(temp);
//                cout << queryCount << endl;
//                cout << "query : " << lookingQuery << endl;
//                cout << "target: " << lookingTarget << endl;

                if (isMatched == 0) {
                    lastFirstMatch = lookingTarget;
                    lastFirstDiffIdxPos = lookingTargetPos;
                    lastFirstTargetIdx = j;
                    isMatched = 1;
                }
                if ((nextTarget & marker) != (lookingQuery & marker)) {
                    break;
                }
                if ((nextTarget & marker) == (lookingQuery & marker)) {
//                    cout<<"Multiple Match "<<queryCount<<endl;
//                    cout<<"query : "<<lookingQuery<<endl;
//                    cout<<"target: "<<lookingTarget<<endl;
//                    cout<<"next  : "<<nextTarget<<endl;
                    multipleMatch++;
                }
            }

            if((nextTarget & marker) > (lookingQuery & marker)) {
                break;
            }

        }
        if((nextTarget & marker) == (lookingQuery & marker))
        {
            if(nextTarget == lookingQuery) perfect ++;
            cout<<queryCount<<"last"<<endl;
            cout<<"query : "<<lookingQuery<<endl;
            cout<<"target: "<<nextTarget<<endl;
            matchCount++;
            matchedKmer temp = {nextTarget, lookingQuery, kmerBuffer[i].info.sequenceID, targetInfoList.data[maxTarget-1].sequenceID,
                                kmerBuffer[i].info.pos - targetInfoList.data[maxTarget-1].pos,
                                getHammingDistance(lookingQuery, nextTarget)};
            matchedKmerList.push_back(temp);
        }
        nextTarget = lastFirstMatch;
        diffIdxPos = lastFirstDiffIdxPos;
    }
    int asdd = matchedKmerList.size();
    cout<<"query count                          : "<<queryCount<<endl;
    cout<<"Total match count                    : "<< matchCount <<endl;
    cout<<"mutipleMatch in AA level             : "<< multipleMatch << endl;
    int pm = 0;
    for(int i = 0 ; i < asdd; i++)
    {
        if(matchedKmerList[i].hammingDistance == 0)
            pm++;
    }
    cout<<"matches in DNA level                 : "<<pm;

}
uint64_t Searcher::getNextTargetKmer(uint64_t lookingTarget, const uint16_t* targetDiffIdxList, size_t & diffIdxPos)
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

uint8_t Searcher::getHammingDistance(uint64_t kmer1, uint64_t kmer2)
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