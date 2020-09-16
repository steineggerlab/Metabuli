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

    char taxID[100];
    FILE * taxIdFile = fopen("/Users/kjb/Desktop/ADclassifier/tengenome/taxIDs", "r");
    vector<int> taxIdList;

    while(feof(taxIdFile) == 0)
    {
        fscanf(taxIdFile,"%s",taxID);
        taxIdList.push_back(atoi(taxID));
    }

    size_t fileCnt = diffIdxFileNames.size();
    size_t leftFile = fileCnt;
    uint64_t lookingKmers[fileCnt];
    KmerInfo lookingInfo[fileCnt];
    size_t diffFileIdx[fileCnt]; for(size_t i = 0; i < fileCnt; i++) diffFileIdx[i] = 0;
    size_t infoFileIdx[fileCnt]; for(size_t i = 0; i < fileCnt; i++) infoFileIdx[i] = 0;
    size_t idxOfMin;
    struct MmapedData<uint16_t> *diffFileList = new struct MmapedData<uint16_t>[fileCnt];
    struct MmapedData<KmerInfo> *infoFileList = new struct MmapedData<KmerInfo>[fileCnt];

    uint16_t * diffBuffer = (uint16_t *)malloc(sizeof(uint16_t) * 250000000); size_t diffBufferIdx = 0;
    KmerInfo * infoBuffer = (KmerInfo *)malloc(sizeof(KmerInfo) * 250000000); size_t infoBufferIdx = 0;
    uint64_t lastWrittenKmer = 0;
    uint64_t lastKmer = 0;
    KmerInfo lastInfo;
    size_t maxIdxOfEachFiles[fileCnt];
    size_t rightnumber = 0;
    //mmap split files
    for (size_t file = 0; file < fileCnt; file++) {
        diffFileList[file] = mmapData<uint16_t>(diffIdxFileNames[file]);
        infoFileList[file] = mmapData<KmerInfo>(infoFileNames[file]);
        maxIdxOfEachFiles[file] = diffFileList[file].fileSize / sizeof(uint16_t);
        rightnumber += infoFileList[file].fileSize/sizeof(KmerInfo);
    }

    size_t totalKmerCnt = 0;
    size_t writtenKmerCnt = 0;
    // get the first entry
    for(size_t file = 0; file < fileCnt; file++)
    {
        lookingKmers[file] = getNextKmer(0, diffFileList[file], diffFileIdx[file]);
        lookingInfo[file] = infoFileList[file].data[0];
        infoFileIdx[file] ++;
        totalKmerCnt ++;
    }
    idxOfMin = smallest(lookingKmers, fileCnt);
    lastKmer = lookingKmers[idxOfMin];
    lastInfo = lookingInfo[idxOfMin];

    cre->writeKmerDiff(0, lastKmer, mergedDiffFIle, diffBuffer, diffBufferIdx);
    lastWrittenKmer = lastKmer;
    ///info바르게 기록되는지 확인할 것
    cre->writeInfo(&lastInfo, mergedInfoFIle, infoBuffer, infoBufferIdx);;
    writtenKmerCnt++;
    size_t intraSpe = 0;
    size_t interSpe = 0;

    int sharedInSpecies = 0;
    int sharedBwSpecies = 0;
    int endFlag = 0;
    int both = 0;
    int interloss = 0;
    int intraloss = 0;
    ///끝부분 잘 되는지 확인할 것
    while(1){
        //update looking K mer
        lookingKmers[idxOfMin] = getNextKmer(lastKmer, diffFileList[idxOfMin], diffFileIdx[idxOfMin]);
        lookingInfo[idxOfMin] = infoFileList[idxOfMin].data[infoFileIdx[idxOfMin]]; infoFileIdx[idxOfMin] ++;
        if( diffFileIdx[idxOfMin] > maxIdxOfEachFiles[idxOfMin] )
        {
            lookingKmers[idxOfMin] = 18446744073709551615;
            leftFile--;
            if(leftFile == 0) break;
        }
        totalKmerCnt ++;
        idxOfMin = smallest(lookingKmers, fileCnt);

        int asd = 0;
        while(lastKmer == lookingKmers[idxOfMin])
        {
            if(taxIdList[lastInfo.sequenceID] == taxIdList[lookingInfo[idxOfMin].sequenceID])
            {
                sharedInSpecies = 1;
                intraSpe ++;
            }
            else{
                sharedBwSpecies = 1;
                interSpe ++;
            }
            lookingKmers[idxOfMin] = getNextKmer(lookingKmers[idxOfMin],diffFileList[idxOfMin],diffFileIdx[idxOfMin]);
            totalKmerCnt ++;
            lookingInfo[idxOfMin] = infoFileList[idxOfMin].data[infoFileIdx[idxOfMin]]; infoFileIdx[idxOfMin] ++;

            if( diffFileIdx[idxOfMin] > maxIdxOfEachFiles[idxOfMin] )
            {
                lookingKmers[idxOfMin] = 18446744073709551615;
                leftFile--;
                if(leftFile == 0)
                {
                    endFlag = 1;
                    break;
                }
            }
            asd++;
            idxOfMin = smallest(lookingKmers, fileCnt);
        }
        if(sharedBwSpecies == 0)
        {
            cre->writeKmerDiff(lastWrittenKmer, lastKmer, mergedDiffFIle, diffBuffer, diffBufferIdx);
            cre->writeInfo(&lastInfo, mergedInfoFIle, infoBuffer, infoBufferIdx);
            writtenKmerCnt++;
        }
        if(endFlag == 1)
        {
            break;
        }
        sharedBwSpecies = 0;
        sharedInSpecies = 0;

        //update last k-mer
        lastKmer = lookingKmers[idxOfMin];
        lastInfo = lookingInfo[idxOfMin];
    }

    cre->flushInfoBuf(infoBuffer, mergedInfoFIle, infoBufferIdx);
    cre->flushKmerBuf(diffBuffer, mergedDiffFIle, diffBufferIdx);

    free(diffBuffer);
    free(infoBuffer);
    for(size_t file = 0; file < fileCnt; file++){
        munmap(diffFileList[file].data, diffFileList[file].fileSize + 1);
        munmap(infoFileList[file].data, infoFileList[file].fileSize + 1);
    }
    cout<<"please: "<<rightnumber<<endl;
    cout<<"Creating target DB is done"<<endl;
    cout<<"Total k-mer count    : "<<totalKmerCnt<<endl;
    cout<<"Written k-mer count  : "<< writtenKmerCnt << endl;
    cout<<"intra loss           : "<<intraloss<<endl;
    cout<<"inter loss           : "<<interloss<<endl;
    cout<<"within a species     : "<<intraSpe<<endl;
    cout<<"between species      : "<<interSpe<<endl;
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