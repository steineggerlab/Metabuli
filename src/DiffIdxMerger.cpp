//
// Created by KJB on 01/09/2020.
//

#include "DiffIdxMerger.h"
DiffIdxMerger::DiffIdxMerger(char* mergedDiffFileName, char * mergedInfoFileNmae)
:mergedDiffFileName(mergedDiffFileName), mergedInfoFileName(mergedInfoFileNmae)
{
    cre = new IndexCreator();
}

void DiffIdxMerger::mergeTargetFiles(std::vector<char*> diffIdxFileNames, std::vector<char*> infoFileNames) {

    FILE * mergedDiffFIle = fopen(mergedDiffFileName, "wb");
    FILE * mergedInfoFIle = fopen(mergedInfoFileName, "wb");

    size_t fileCnt = diffIdxFileNames.size();
    size_t leftFile = fileCnt;
    uint64_t lookingKmers[fileCnt];
    KmerInfo lookingInfo[fileCnt];
    size_t diffFileIdx[fileCnt]; for(size_t i = 0; i < fileCnt; i++) diffFileIdx[i] = 0;
    size_t idxOfMin;
    struct MmapedData<uint16_t> *diffFileList = new struct MmapedData<uint16_t>[fileCnt];
    struct MmapedData<KmerInfo> *infoFileList = new struct MmapedData<KmerInfo>[fileCnt];

    uint16_t * diffBuffer = (uint16_t *)malloc(sizeof(uint16_t) * 250000000); size_t diffBufferIdx = 0;
    KmerInfo * infoBuffer = (KmerInfo *)malloc(sizeof(KmerInfo) * 250000000); size_t infoBufferIdx = 0;
    uint64_t lastEntry = 0;
    size_t maxIdxOfEachFiles[fileCnt];

    for (size_t file = 0; file < fileCnt; file++) {
        diffFileList[file] = mmapData<uint16_t>(diffIdxFileNames[file]);
        infoFileList[file] = mmapData<KmerInfo>(infoFileNames[file]);
        maxIdxOfEachFiles[file] = diffFileList[file].fileSize / sizeof(uint16_t);
    }

    size_t totalKmerCnt = 0;
    for (size_t file = 0; file < fileCnt; file++) {
        totalKmerCnt += infoFileList[file].fileSize / sizeof(KmerInfo);
    }

    // get the first entry
    for(size_t file = 0; file < fileCnt; file++)
    {
        lookingKmers[file] = getNextKmer(0, diffFileList[file], diffFileIdx[file]);
        lookingInfo[file] = infoFileList[file].data[0];
    }
    idxOfMin = smallest(lookingKmers, fileCnt);
    lastEntry = lookingKmers[idxOfMin];

    cre->writeKmerDiff(0, lastEntry, mergedDiffFIle, diffBuffer, diffBufferIdx);
    cre->writeInfo(&lookingInfo[idxOfMin], mergedInfoFIle, infoBuffer, infoBufferIdx);

    for (size_t i = 0; i < totalKmerCnt - 1; i++)
    while(leftFile){
        lookingKmers[idxOfMin] = getNextKmer(lastEntry, diffFileList[idxOfMin], diffFileIdx[idxOfMin]);
        //cout<<diffFileIdx[idxOfMin]<<" "<<maxIdxOfEachFiles[idxOfMin]<<endl;
        if( diffFileIdx[idxOfMin] > maxIdxOfEachFiles[idxOfMin] )
        {
            lookingKmers[idxOfMin] = 18446744073709551615;
            leftFile--;
        }
        idxOfMin = smallest(lookingKmers, fileCnt);
        //fill buffer with min. k-mer
        cre->writeKmerDiff(lastEntry, lookingKmers[idxOfMin], mergedDiffFIle, diffBuffer, diffBufferIdx);
        cre->writeInfo(&lookingInfo[idxOfMin], mergedInfoFIle, infoBuffer, infoBufferIdx);
        //update last entry
        lastEntry = lookingKmers[idxOfMin];
    }

    cre->flushInfoBuf(infoBuffer, mergedInfoFIle, infoBufferIdx);
    cre->flushKmerBuf(diffBuffer, mergedDiffFIle, diffBufferIdx);

    free(diffBuffer);
    free(infoBuffer);
    for(size_t file = 0; file < fileCnt; file++){
        munmap(diffFileList[file].data, diffFileList[file].fileSize + 1);
        munmap(infoFileList[file].data, infoFileList[file].fileSize + 1);
    }
}
uint64_t DiffIdxMerger::getNextKmer(uint64_t lookingTarget, const struct MmapedData<uint16_t> diffList, size_t & idx)
{
    uint16_t fragment = 0;
    uint64_t diffIn64bit = 0;

    for(int i = 0; i < 5; i++)
    {
        fragment = diffList.data[idx];
        idx++;
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

size_t DiffIdxMerger::smallest(const uint64_t lookingKmers[], const size_t & fileCnt)
{
    size_t idxOfMin = 0;
    uint64_t min = lookingKmers[0];

    for(size_t i = 1; i < fileCnt; i++)
    {
        if(lookingKmers[i] < min)
        {
            min = lookingKmers[i];
            idxOfMin = i;
        }
    }
    return idxOfMin;
}