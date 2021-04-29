//
// Created by KJB on 01/09/2020.
//

#include "FileMerger.h"

FileMerger::FileMerger(char* mergedDiffFileName, char * mergedInfoFileNmae, char * diffIdxSplitFileName)
:mergedDiffFileName(mergedDiffFileName), mergedInfoFileName(mergedInfoFileNmae), diffIdxSplitFileName(diffIdxSplitFileName)
{
    cre = new IndexCreator();
}

///Merge differential index and k-mer information files, reducing redundancy
void FileMerger::mergeTargetFiles(std::vector<char*> diffIdxFileNames, std::vector<char*> infoFileNames, vector<int> & taxIdListAtRank, vector<int> & taxIdList) {
    size_t writtenKmerCnt = 0;

    ///Files to write on & buffers to fill them
    FILE * mergedDiffFile = fopen(mergedDiffFileName, "wb");
    FILE * mergedInfoFile = fopen(mergedInfoFileName, "wb");
    FILE * diffIdxSplitFile = fopen(diffIdxSplitFileName, "wb");
    uint16_t * diffBuffer = (uint16_t *)malloc(sizeof(uint16_t) * kmerBufSize);
    size_t diffBufferIdx = 0;
    size_t totalBufferIdx = 0;
    TargetKmerInfo * infoBuffer = (TargetKmerInfo *)malloc(sizeof(TargetKmerInfo) * kmerBufSize);
    size_t infoBufferIdx = 0;
    size_t totalInfoIdx = 0;


    ///Prepare files to merge
    size_t numOfSplitFiles = diffIdxFileNames.size();
    size_t numOfincompletedFiles = numOfSplitFiles;
    size_t numOfKmerBeforeMerge = 0;
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
        numOfKmerBeforeMerge += infoFileList[file].fileSize / sizeof(TargetKmerInfo);
    }

    ///To make differential index splits
    uint64_t AAofTempSplitOffset = UINT64_MAX;
    size_t sizeOfSplit = numOfKmerBeforeMerge / (SplitNum - 1);
    size_t offsetList[SplitNum];
    int offsetListIdx = 1;
    for(size_t os = 0; os < SplitNum; os++){
        offsetList[os] = os * sizeOfSplit;
    }

    DiffIdxSplit splitList[SplitNum];
    memset(splitList, 0, sizeof(splitList));
    //splitList[0] = {0, 0 ,0};
    int splitListIdx = 1;

    /// get the first k-mer to write
    for(size_t file = 0; file < numOfSplitFiles; file++){
        lookingKmers[file] = getNextKmer(0, diffFileList[file], diffFileIdx[file]);
        lookingInfos[file] = infoFileList[file].data[0];
        infoFileIdx[file] ++;
    }

    size_t idxOfMin = smallest(lookingKmers, lookingInfos, taxIdListAtRank, numOfSplitFiles);
    uint64_t lastWrittenKmer = 0;
    uint64_t entryKmer = lookingKmers[idxOfMin];
    TargetKmerInfo entryInfo = lookingInfos[idxOfMin];

    ///write first k-mer
    cre->getDiffIdx(0, entryKmer, mergedDiffFile, diffBuffer, diffBufferIdx);
    lastWrittenKmer = entryKmer;
    cre->writeInfo(&entryInfo, mergedInfoFile, infoBuffer, infoBufferIdx);;
    writtenKmerCnt++;
    int splitCheck = 0;
    int endFlag = 0;
    while(true){
        ///update entry k-mer
        entryKmer = lookingKmers[idxOfMin];
        entryInfo = lookingInfos[idxOfMin];

        ///update looking k-mers
        lookingKmers[idxOfMin] = getNextKmer(entryKmer, diffFileList[idxOfMin], diffFileIdx[idxOfMin]);
        lookingInfos[idxOfMin] = infoFileList[idxOfMin].data[infoFileIdx[idxOfMin]];
        infoFileIdx[idxOfMin] ++;
        if( diffFileIdx[idxOfMin] > maxIdxOfEachFiles[idxOfMin] ){
            lookingKmers[idxOfMin] = UINT64_MAX;
            numOfincompletedFiles--;
            if(numOfincompletedFiles == 0) break;
        }
        idxOfMin = smallest(lookingKmers, lookingInfos, taxIdListAtRank, numOfSplitFiles);

        int hasSeenOtherStrains = 0;
        while(taxIdListAtRank[entryInfo.sequenceID] == taxIdListAtRank[lookingInfos[idxOfMin].sequenceID]){
            if(entryKmer != lookingKmers[idxOfMin]) break;

            hasSeenOtherStrains += (taxIdList[entryInfo.sequenceID] != taxIdList[lookingInfos[idxOfMin].sequenceID]);

            lookingKmers[idxOfMin] = getNextKmer(entryKmer, diffFileList[idxOfMin], diffFileIdx[idxOfMin]);
            lookingInfos[idxOfMin] = infoFileList[idxOfMin].data[infoFileIdx[idxOfMin]];
            infoFileIdx[idxOfMin] ++;

            if(diffFileIdx[idxOfMin] > maxIdxOfEachFiles[idxOfMin] ){
                lookingKmers[idxOfMin] = UINT64_MAX;
                numOfincompletedFiles--;
                if(numOfincompletedFiles == 0){
                    endFlag = 1;
                    break;
                }
            }
            idxOfMin = smallest(lookingKmers, lookingInfos, taxIdListAtRank, numOfSplitFiles);
        }

        entryInfo.redundancy = (hasSeenOtherStrains > 0 || entryInfo.redundancy);
        getDiffIdx(lastWrittenKmer, entryKmer, mergedDiffFile, diffBuffer, diffBufferIdx, totalBufferIdx);
        lastWrittenKmer = entryKmer;
        writeInfo(&entryInfo, mergedInfoFile, infoBuffer, infoBufferIdx, totalInfoIdx);
        writtenKmerCnt++;

        if(AminoAcid(lastWrittenKmer) != AAofTempSplitOffset && splitCheck == 1){
            splitList[splitListIdx++] = {lastWrittenKmer, totalBufferIdx, totalInfoIdx};
            splitCheck = 0;
        }

        if(writtenKmerCnt == offsetList[offsetListIdx]){
            AAofTempSplitOffset = AminoAcid(lastWrittenKmer);
            splitCheck = 1;
            offsetListIdx++;
        }

        if(endFlag == 1) break;
    }

    cre->flushInfoBuf(infoBuffer, mergedInfoFile, infoBufferIdx);
    cre->flushKmerBuf(diffBuffer, mergedDiffFile, diffBufferIdx);
    fwrite(splitList, sizeof(DiffIdxSplit), SplitNum, diffIdxSplitFile);
    for(int i = 0; i < SplitNum; i++){
        cout<<splitList[i].ADkmer<< " "<<splitList[i].diffIdxOffset<< " "<<splitList[i].infoIdxOffset<<endl;
    }
    free(diffBuffer);
    free(infoBuffer);
    fclose(mergedDiffFile);
    fclose(mergedInfoFile);
    fclose(diffIdxSplitFile);

    for(size_t file = 0; file < numOfSplitFiles; file++){
        munmap(diffFileList[file].data, diffFileList[file].fileSize + 1);
        munmap(infoFileList[file].data, infoFileList[file].fileSize + 1);
    }
    cout<<"Creating target DB is done"<<endl;
    cout<<"Total k-mer count    : " << numOfKmerBeforeMerge <<endl;
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
    uint16_t * diffBuffer = (uint16_t *)malloc(sizeof(uint16_t) * kmerBufSize); size_t diffBufferIdx = 0;
    TargetKmerInfo * infoBuffer = (TargetKmerInfo *)malloc(sizeof(TargetKmerInfo) * kmerBufSize); size_t infoBufferIdx = 0;


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

    ///TODO check the for loop condition; size
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

    size_t idxOfMin = smallest(lookingKmers, lookingInfos, taxListAtRank, numOfSplitFiles);
    uint64_t lastWrittenKmer = 0;
    uint64_t lastKmer = lookingKmers[idxOfMin];
    TargetKmerInfo lastInfo = lookingInfos[idxOfMin];

    ///write first k-mer
    cre->getDiffIdx(0, lastKmer, mergedDiffFile, diffBuffer, diffBufferIdx);
    lastWrittenKmer = lastKmer;
    cre->writeInfo(&lastInfo, mergedInfoFile, infoBuffer, infoBufferIdx);;
    writtenKmerCnt++;

    int endFlag = 0;
    int hasSeenOtherStrains;
    ///끝부분 잘 되는지 확인할 것
    while(true){
        ///update last k-mer
        lastKmer = lookingKmers[idxOfMin];
        lastInfo = lookingInfos[idxOfMin];

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
        idxOfMin = smallest(lookingKmers, lookingInfos, taxListAtRank, numOfSplitFiles);

        hasSeenOtherStrains = 0;
        while(taxListAtRank[lastInfo.sequenceID] == taxListAtRank[lookingInfos[idxOfMin].sequenceID]){
            if(lastKmer != lookingKmers[idxOfMin]) break;

            hasSeenOtherStrains += (taxIdList[lastInfo.sequenceID] != taxIdList[lookingInfos[idxOfMin].sequenceID]);
            lookingKmers[idxOfMin] = getNextKmer(lastKmer, diffFileList[idxOfMin], diffFileIdx[idxOfMin]);
            lookingInfos[idxOfMin] = infoFileList[idxOfMin].data[infoFileIdx[idxOfMin]];
            infoFileIdx[idxOfMin] ++;
            totalKmerCnt ++;

            if(diffFileIdx[idxOfMin] > maxIdxOfEachFiles[idxOfMin]){
                lookingKmers[idxOfMin] = UINT64_MAX;
                numOfincompletedFiles--;
                if(numOfincompletedFiles == 0){
                    endFlag = 1;
                    break;
                }
            }
            idxOfMin = smallest(lookingKmers, lookingInfos, taxListAtRank, numOfSplitFiles);
        }

        lastInfo.redundancy = (hasSeenOtherStrains > 0 || lastInfo.redundancy);
        cre->getDiffIdx(lastWrittenKmer, lastKmer, mergedDiffFile, diffBuffer, diffBufferIdx);
        lastWrittenKmer = lastKmer;
        cre->writeInfo(&lastInfo, mergedInfoFile, infoBuffer, infoBufferIdx);
        writtenKmerCnt++;

        if(endFlag == 1) break;
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

size_t FileMerger::smallest(const uint64_t lookingKmers[], const TargetKmerInfo lookingInfos[], vector<int> & taxListAtRank, const size_t & fileCnt)
{
    size_t idxOfMin = 0;
    uint64_t min = lookingKmers[0];
    int minTaxIdAtRank = taxListAtRank[lookingInfos[0].sequenceID];
    for(size_t i = 1; i < fileCnt; i++)
    {
        if(lookingKmers[i] < min ||(lookingKmers[i] == min && taxListAtRank[lookingInfos[i].sequenceID] < minTaxIdAtRank)){
            min = lookingKmers[i];
            minTaxIdAtRank = taxListAtRank[lookingInfos[i].sequenceID];
            idxOfMin = i;
        }
    }
    return idxOfMin;
}

void FileMerger::getDiffIdx(const uint64_t & lastKmer, const uint64_t & entryToWrite, FILE* handleKmerTable, uint16_t *kmerBuf, size_t & localBufIdx, size_t & totalBufIdx){
    uint64_t kmerdiff = entryToWrite - lastKmer;
    uint16_t buffer[5];
    int idx = 3;
    buffer[4] = SET_END_FLAG(GET_15_BITS(kmerdiff));
    kmerdiff >>= 15U;
    while (kmerdiff) {
        uint16_t toWrite = GET_15_BITS(kmerdiff);
        kmerdiff >>= 15U;
        buffer[idx] = toWrite;
        idx--;
    }
    writeDiffIdx(kmerBuf, handleKmerTable, (buffer + idx + 1), (4 - idx), localBufIdx, totalBufIdx);
}
void FileMerger::writeDiffIdx(uint16_t *buffer, FILE* handleKmerTable, uint16_t *toWrite, size_t size, size_t & localBufIdx, size_t & totalBufIdx) {
    if (localBufIdx + size >= kmerBufSize) {
        flushKmerBuf(buffer, handleKmerTable, localBufIdx);
    }
    memcpy(buffer + localBufIdx, toWrite, sizeof(uint16_t) * size);
    localBufIdx += size;
    totalBufIdx += size;
}

void FileMerger::flushKmerBuf(uint16_t *buffer, FILE *handleKmerTable, size_t & localBufIdx ) {
    fwrite(buffer, sizeof(uint16_t), localBufIdx, handleKmerTable);
    localBufIdx = 0;
}

void FileMerger::writeInfo(TargetKmerInfo * entryToWrite, FILE * infoFile, TargetKmerInfo * infoBuffer, size_t & infoBufferIdx, size_t & totalInfoIdx)
{
    if (infoBufferIdx >= kmerBufSize) {
        flushInfoBuf(infoBuffer, infoFile, infoBufferIdx);
    }
    memcpy(infoBuffer + infoBufferIdx, entryToWrite, sizeof(TargetKmerInfo));
    infoBufferIdx++;
    totalInfoIdx ++;
}
void FileMerger::flushInfoBuf(TargetKmerInfo * buffer, FILE * infoFile, size_t & localBufIdx ) {
    fwrite(buffer, sizeof(TargetKmerInfo), localBufIdx, infoFile);
    localBufIdx = 0;
}