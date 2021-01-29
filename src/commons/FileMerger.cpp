//
// Created by KJB on 01/09/2020.
//

#include "FileMerger.h"
FileMerger::FileMerger(char* mergedDiffFileName, char * mergedInfoFileNmae)
:mergedDiffFileName(mergedDiffFileName), mergedInfoFileName(mergedInfoFileNmae)
{
    cre = new IndexCreator();
}

///Merge differential index and k-mer information files, reducing redundancy
void FileMerger::mergeTargetFiles(std::vector<char*> diffIdxFileNames, std::vector<char*> infoFileNames, vector<int> & taxIdListAtRank, vector<int> & taxIdList) {
    size_t numOfKmerBeforeMerge = 0;
    size_t writtenKmerCnt = 0;

    ///Files to write on & buffers to fill them
    FILE * mergedDiffFile = fopen(mergedDiffFileName, "wb");
    FILE * mergedInfoFile = fopen(mergedInfoFileName, "wb");
    uint16_t * diffBuffer = (uint16_t *)malloc(sizeof(uint16_t) * 250000000); size_t diffBufferIdx = 0;
    TargetKmerInfo * infoBuffer = (TargetKmerInfo *)malloc(sizeof(TargetKmerInfo) * 250000000); size_t infoBufferIdx = 0;


    ///Prepare files to merge
    size_t numOfSplitFiles = diffIdxFileNames.size();
    size_t numOfincompletedFiles = numOfSplitFiles;
    uint64_t lookingKmers[numOfSplitFiles];
    TargetKmerInfo lookingInfos[numOfSplitFiles];
    size_t diffFileIdx[numOfSplitFiles];
    memset(diffFileIdx, 0, sizeof(diffFileIdx));
    size_t infoFileIdx[numOfSplitFiles];
    memset(infoFileIdx, 0, sizeof(infoFileIdx));
    size_t maxIdxOfEachFiles[numOfSplitFiles];
    struct MmapedData<uint16_t> *diffFileList = new struct MmapedData<uint16_t>[numOfSplitFiles];
    struct MmapedData<TargetKmerInfo> *infoFileList = new struct MmapedData<TargetKmerInfo>[numOfSplitFiles];
    for (size_t file = 0; file < numOfSplitFiles; file++) {
        diffFileList[file] = mmapData<uint16_t>(diffIdxFileNames[file]);
        infoFileList[file] = mmapData<TargetKmerInfo>(infoFileNames[file]);
        maxIdxOfEachFiles[file] = diffFileList[file].fileSize / sizeof(uint16_t);
    }

    /// get the first k-mer to write
    for(size_t file = 0; file < numOfSplitFiles; file++){
        lookingKmers[file] = getNextKmer(0, diffFileList[file], diffFileIdx[file]);
        lookingInfos[file] = infoFileList[file].data[0];
        infoFileIdx[file] ++;
    }

    size_t idxOfMin = smallest(lookingKmers, numOfSplitFiles);
    uint64_t lastWrittenKmer = 0;
    uint64_t lastKmer = lookingKmers[idxOfMin];
    TargetKmerInfo lastInfo = lookingInfos[idxOfMin];

    ///write first k-mer
    cre->getDiffIdx(0, lastKmer, mergedDiffFile, diffBuffer, diffBufferIdx);
    lastWrittenKmer = lastKmer;
    cre->writeInfo(&lastInfo, mergedInfoFile, infoBuffer, infoBufferIdx);;
    writtenKmerCnt++;

    int endFlag = 0;
    int unique = 0;
    ///끝부분 잘 되는지 확인할 것
    while(1){
        ///update looking k-mers
        lookingKmers[idxOfMin] = getNextKmer(lastKmer, diffFileList[idxOfMin], diffFileIdx[idxOfMin]);
        lookingInfos[idxOfMin] = infoFileList[idxOfMin].data[infoFileIdx[idxOfMin]];
        infoFileIdx[idxOfMin] ++;
        if( diffFileIdx[idxOfMin] > maxIdxOfEachFiles[idxOfMin] ){
            lookingKmers[idxOfMin] = UINT64_MAX;
            numOfincompletedFiles--;
            if(numOfincompletedFiles == 0) break;
        }
        numOfKmerBeforeMerge ++;
        idxOfMin = smallest(lookingKmers, numOfSplitFiles);

        int asd = 0;
        while(taxIdListAtRank[lastInfo.sequenceID] == taxIdListAtRank[lookingInfos[idxOfMin].sequenceID]){
            if(lastKmer != lookingKmers[idxOfMin]) break;

            lookingKmers[idxOfMin] = getNextKmer(lastKmer, diffFileList[idxOfMin], diffFileIdx[idxOfMin]);
            lookingInfos[idxOfMin] = infoFileList[idxOfMin].data[infoFileIdx[idxOfMin]];
            infoFileIdx[idxOfMin] ++;
            numOfKmerBeforeMerge ++;
            asd++;

            if( diffFileIdx[idxOfMin] > maxIdxOfEachFiles[idxOfMin] ){
                lookingKmers[idxOfMin] = UINT64_MAX;
                numOfincompletedFiles--;
                if(numOfincompletedFiles == 0){
                    endFlag = 1;
                    break;
                }
            }
            idxOfMin = smallest(lookingKmers, numOfSplitFiles);
        }

        if(asd == 0){
            lastInfo.redundancy = false;
        }

        cre->getDiffIdx(lastWrittenKmer, lastKmer, mergedDiffFile, diffBuffer, diffBufferIdx);
        lastWrittenKmer = lastKmer;
        cre->writeInfo(&lastInfo, mergedInfoFile, infoBuffer, infoBufferIdx);
        writtenKmerCnt++;

        if(endFlag == 1){
            break;
        }

        ///update last k-mer
        lastKmer = lookingKmers[idxOfMin];
        lastInfo = lookingInfos[idxOfMin];
    }

    cre->flushInfoBuf(infoBuffer, mergedInfoFile, infoBufferIdx);
    cre->flushKmerBuf(diffBuffer, mergedDiffFile, diffBufferIdx);

    free(diffBuffer);
    free(infoBuffer);
    fclose(mergedDiffFile);
    fclose(mergedInfoFile);
    for(size_t file = 0; file < numOfSplitFiles; file++){
        munmap(diffFileList[file].data, diffFileList[file].fileSize + 1);
        munmap(infoFileList[file].data, infoFileList[file].fileSize + 1);
    }

    cout<<"Creating target DB is done"<<endl;
    cout<<"Total k-mer count    : " << numOfKmerBeforeMerge << endl;
    cout<<"Written k-mer count  : " << writtenKmerCnt << endl;
}

///It updates target database. Most of this function is the same with 'mergeTargetFiles', only some additional tasks are added to handle seqID
///It is not tested yet 2020.01.28
void FileMerger::updateTargetDatabase(std::vector<char*> diffIdxFileNames, std::vector<char*> infoFileNames, vector<int> & taxListAtRank, vector<int> & taxIdList, const int & seqIdOffset) {
    size_t totalKmerCnt = 0;
    size_t writtenKmerCnt = 0;

    ///Files to write on & buffers to fill them
    FILE * mergedDiffFile = fopen(mergedDiffFileName, "wb");
    FILE * mergedInfoFile = fopen(mergedInfoFileName, "wb");
    uint16_t * diffBuffer = (uint16_t *)malloc(sizeof(uint16_t) * 250000000); size_t diffBufferIdx = 0;
    TargetKmerInfo * infoBuffer = (TargetKmerInfo *)malloc(sizeof(TargetKmerInfo) * 250000000); size_t infoBufferIdx = 0;


    ///Prepare files to merge
    size_t numOfSplitFiles = diffIdxFileNames.size();
    size_t numOfincompletedFiles = numOfSplitFiles;
    uint64_t lookingKmers[numOfSplitFiles];
    TargetKmerInfo lookingInfos[numOfSplitFiles];
    size_t diffFileIdx[numOfSplitFiles];
    memset(diffFileIdx, 0, sizeof(diffFileIdx));
    size_t infoFileIdx[numOfSplitFiles];
    memset(infoFileIdx, 0, sizeof(infoFileIdx));
    size_t maxIdxOfEachFiles[numOfSplitFiles];
    struct MmapedData<uint16_t> *diffFileList = new struct MmapedData<uint16_t>[numOfSplitFiles];
    struct MmapedData<TargetKmerInfo> *infoFileList = new struct MmapedData<TargetKmerInfo>[numOfSplitFiles];
    for (size_t file = 0; file < numOfSplitFiles; file++) {
        diffFileList[file] = mmapData<uint16_t>(diffIdxFileNames[file]);
        infoFileList[file] = mmapData<TargetKmerInfo>(infoFileNames[file]);
        maxIdxOfEachFiles[file] = diffFileList[file].fileSize / sizeof(uint16_t);
    }

    ///Fix sequence IDs of new k-mers
    for (int i = 1; i < numOfSplitFiles; i++){
        size_t size = (infoFileList[i].fileSize - 1) / sizeof(TargetKmerInfo);
        for(int j = 0; j < size + 1; j++){
            infoFileList[i].data[j].sequenceID += seqIdOffset;
        }
    }

    /// get the first k-mer to write
    for(size_t file = 0; file < numOfSplitFiles; file++){
        lookingKmers[file] = getNextKmer(0, diffFileList[file], diffFileIdx[file]);
        lookingInfos[file] = infoFileList[file].data[0];
        infoFileIdx[file] ++;
    }

    size_t idxOfMin = smallest(lookingKmers, numOfSplitFiles);
    uint64_t lastWrittenKmer = 0;
    uint64_t lastKmer = lookingKmers[idxOfMin];
    TargetKmerInfo lastInfo = lookingInfos[idxOfMin];

    ///write first k-mer
    cre->getDiffIdx(0, lastKmer, mergedDiffFile, diffBuffer, diffBufferIdx);
    lastWrittenKmer = lastKmer;
    cre->writeInfo(&lastInfo, mergedInfoFile, infoBuffer, infoBufferIdx);;
    writtenKmerCnt++;

    int endFlag = 0;
    int isNew = 0;
    int isRedundant = 0;
    int hasSeenOtherStrains = 0;
    ///끝부분 잘 되는지 확인할 것
    while(1){
        ///update looking k-mers
        lookingKmers[idxOfMin] = getNextKmer(lastKmer, diffFileList[idxOfMin], diffFileIdx[idxOfMin]);
        lookingInfos[idxOfMin] = infoFileList[idxOfMin].data[infoFileIdx[idxOfMin]];
        infoFileIdx[idxOfMin] ++;
        if( diffFileIdx[idxOfMin] > maxIdxOfEachFiles[idxOfMin] ){
            lookingKmers[idxOfMin] = UINT64_MAX;
            numOfincompletedFiles--;
            if(numOfincompletedFiles == 0) break;
        }
        totalKmerCnt ++;
        idxOfMin = smallest(lookingKmers, numOfSplitFiles);

        hasSeenOtherStrains = 0;
        isRedundant = 0;
        while(taxListAtRank[lastInfo.sequenceID] == taxListAtRank[lookingInfos[idxOfMin].sequenceID]){
            if(lastKmer != lookingKmers[idxOfMin]) break;

            if (taxIdList[lastInfo.sequenceID] != taxIdList[lookingInfos[idxOfMin].sequenceID]){
                hasSeenOtherStrains = 1;
            }

            lookingKmers[idxOfMin] = getNextKmer(lastKmer, diffFileList[idxOfMin], diffFileIdx[idxOfMin]);
            lookingInfos[idxOfMin] = infoFileList[idxOfMin].data[infoFileIdx[idxOfMin]];
            infoFileIdx[idxOfMin] ++;
            totalKmerCnt ++;
            isRedundant = 1;

            if( diffFileIdx[idxOfMin] > maxIdxOfEachFiles[idxOfMin] ){
                lookingKmers[idxOfMin] = UINT64_MAX;
                numOfincompletedFiles--;
                if(numOfincompletedFiles == 0){
                    endFlag = 1;
                    break;
                }
            }
            idxOfMin = smallest(lookingKmers, numOfSplitFiles);
        }

        if(hasSeenOtherStrains == 1){
            lastInfo.redundancy = true;
        }

        if(isNew == 1) lastInfo.sequenceID += seqIdOffset;

        cre->getDiffIdx(lastWrittenKmer, lastKmer, mergedDiffFile, diffBuffer, diffBufferIdx);
        lastWrittenKmer = lastKmer;

        cre->writeInfo(&lastInfo, mergedInfoFile, infoBuffer, infoBufferIdx);
        writtenKmerCnt++;

        if(endFlag == 1){
            break;
        }

        ///update last k-mer
        lastKmer = lookingKmers[idxOfMin];
        lastInfo = lookingInfos[idxOfMin];
        if(idxOfMin > 0){
            isNew = 1;
        }
    }

    cre->flushInfoBuf(infoBuffer, mergedInfoFile, infoBufferIdx);
    cre->flushKmerBuf(diffBuffer, mergedDiffFile, diffBufferIdx);

    free(diffBuffer);
    free(infoBuffer);
    fclose(mergedDiffFile);
    fclose(mergedInfoFile);
    for(size_t file = 0; file < numOfSplitFiles; file++){
        munmap(diffFileList[file].data, diffFileList[file].fileSize + 1);
        munmap(infoFileList[file].data, infoFileList[file].fileSize + 1);
    }

    cout<<"Creating target DB is done"<<endl;
    cout << "Total k-mer count    : " << totalKmerCnt << endl;
    cout<<"Written k-mer count  : " << writtenKmerCnt << endl;
}

uint64_t FileMerger::getNextKmer(uint64_t lookingTarget, const struct MmapedData<uint16_t> & diffList, size_t & idx)
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

size_t FileMerger::smallest(const uint64_t lookingKmers[], const size_t & fileCnt)
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