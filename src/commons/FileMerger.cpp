#include "FileMerger.h"

FileMerger::FileMerger(const LocalParameters & par) {
    splitNum = par.splitNum;
    if (par.reducedAA == 1){
        MARKER = 0Xffffffff;
        MARKER = ~MARKER;
    } else {
        MARKER = 16777215;
        MARKER = ~ MARKER;
    }
}

FileMerger::~FileMerger() {

}

//void FileMerger::mergeTargetFiles(std::vector<char*> diffIdxFileNames, std::vector<char*> infoFileNames, vector<int> & taxIdListAtRank, vector<int> & taxIdList) {
//    size_t writtenKmerCnt = 0;
//
//    ///Files to write on & buffers to fill them
//    FILE * mergedDiffFile = fopen(mergedDiffFileName, "wb");
//    FILE * mergedInfoFile = fopen(mergedInfoFileName, "wb");
//    FILE * diffIdxSplitFile = fopen(diffIdxSplitFileName, "wb");
//    uint16_t * diffBuffer = (uint16_t *)malloc(sizeof(uint16_t) * kmerBufSize);
//    size_t diffBufferIdx = 0;
//    size_t totalBufferIdx = 0;
//    TargetKmerInfo * infoBuffer = (TargetKmerInfo *)malloc(sizeof(TargetKmerInfo) * kmerBufSize);
//    size_t infoBufferIdx = 0;
//    size_t totalInfoIdx = 0;
//
//    ///Prepare files to merge
//    size_t numOfSplitFiles = diffIdxFileNames.size();
//    size_t numOfincompletedFiles = numOfSplitFiles;
//    size_t numOfKmerBeforeMerge = 0;
//    uint64_t * lookingKmers = new uint64_t[numOfSplitFiles];
////    uint64_t lookingKmers[numOfSplitFiles];
////    TargetKmerInfo lookingInfos[numOfSplitFiles];
//    auto * lookingInfos = new TargetKmerInfo[numOfSplitFiles];
//    //size_t diffFileIdx[numOfSplitFiles];
//    auto * diffFileIdx = new size_t[numOfSplitFiles];
//    memset(diffFileIdx, 0, numOfSplitFiles * sizeof(size_t));
//    auto * infoFileIdx = new size_t[numOfSplitFiles];
////    size_t infoFileIdx[numOfSplitFiles];
//    memset(infoFileIdx, 0, numOfSplitFiles * sizeof(size_t));
//    size_t maxIdxOfEachFiles[numOfSplitFiles];
//    struct MmapedData<uint16_t> *diffFileList = new struct MmapedData<uint16_t>[numOfSplitFiles];
//    struct MmapedData<TargetKmerInfo> *infoFileList = new struct MmapedData<TargetKmerInfo>[numOfSplitFiles];
//    for (size_t file = 0; file < numOfSplitFiles; file++) {
//        diffFileList[file] = mmapData<uint16_t>(diffIdxFileNames[file]);
//        infoFileList[file] = mmapData<TargetKmerInfo>(infoFileNames[file]);
//        maxIdxOfEachFiles[file] = diffFileList[file].fileSize / sizeof(uint16_t);
//        numOfKmerBeforeMerge += infoFileList[file].fileSize / sizeof(TargetKmerInfo);
//    }
//
//    ///To make differential index splits
//    uint64_t AAofTempSplitOffset = UINT64_MAX;
//    size_t sizeOfSplit = numOfKmerBeforeMerge / (SplitNum - 1);
//    size_t offsetList[SplitNum + 1];
//    int offsetListIdx = 1;
//    for(size_t os = 0; os < SplitNum; os++){
//        offsetList[os] = os * sizeOfSplit;
//    }
//    offsetList[SplitNum] = UINT64_MAX;
//
//    DiffIdxSplit splitList[SplitNum];
//    memset(splitList, 0, sizeof(DiffIdxSplit) * SplitNum);
//    int splitListIdx = 1;
//
//    /// get the first k-mer to write
//    for(size_t file = 0; file < numOfSplitFiles; file++){
//        lookingKmers[file] = getNextKmer(0, diffFileList[file], diffFileIdx[file]);
//        lookingInfos[file] = infoFileList[file].data[0];
//        infoFileIdx[file] ++;
//    }
//
//    size_t idxOfMin = smallest(lookingKmers, lookingInfos, taxIdListAtRank, numOfSplitFiles);
//    uint64_t lastWrittenKmer = 0;
//    uint64_t entryKmer = lookingKmers[idxOfMin];
//    TargetKmerInfo entryInfo = lookingInfos[idxOfMin];
//
//    // write first k-mer
//    getDiffIdx(lastWrittenKmer, entryKmer, mergedDiffFile, diffBuffer, diffBufferIdx, totalBufferIdx);
//    lastWrittenKmer = entryKmer;
//    writeInfo(&entryInfo, mergedInfoFile, infoBuffer, infoBufferIdx, totalInfoIdx);
//    writtenKmerCnt++;
//    int splitCheck = 0;
//    int endFlag = 0;
//
//    while(true){
//        // update entry k-mer
//        entryKmer = lookingKmers[idxOfMin];
//        entryInfo = lookingInfos[idxOfMin];
//
//        ///update looking k-mers
//        lookingKmers[idxOfMin] = getNextKmer(entryKmer, diffFileList[idxOfMin], diffFileIdx[idxOfMin]);
//        lookingInfos[idxOfMin] = infoFileList[idxOfMin].data[infoFileIdx[idxOfMin]];
//        infoFileIdx[idxOfMin] ++;
//        if( diffFileIdx[idxOfMin] > maxIdxOfEachFiles[idxOfMin] ){
//            lookingKmers[idxOfMin] = UINT64_MAX;
//            numOfincompletedFiles--;
//            if(numOfincompletedFiles == 0) break;
//        }
//        idxOfMin = smallest(lookingKmers, lookingInfos, taxIdListAtRank, numOfSplitFiles);
//
//        int hasSeenOtherStrains = 0;
//        while(taxIdListAtRank[entryInfo.sequenceID] == taxIdListAtRank[lookingInfos[idxOfMin].sequenceID]){
//            if(entryKmer != lookingKmers[idxOfMin]) break;
//
//            hasSeenOtherStrains += (taxIdList[entryInfo.sequenceID] != taxIdList[lookingInfos[idxOfMin].sequenceID]);
//
//            lookingKmers[idxOfMin] = getNextKmer(entryKmer, diffFileList[idxOfMin], diffFileIdx[idxOfMin]);
//            lookingInfos[idxOfMin] = infoFileList[idxOfMin].data[infoFileIdx[idxOfMin]];
//            infoFileIdx[idxOfMin] ++;
//
//            if(diffFileIdx[idxOfMin] > maxIdxOfEachFiles[idxOfMin] ){
//                lookingKmers[idxOfMin] = UINT64_MAX;
//                numOfincompletedFiles--;
//                if(numOfincompletedFiles == 0){
//                    endFlag = 1;
//                    break;
//                }
//            }
//            idxOfMin = smallest(lookingKmers, lookingInfos, taxIdListAtRank, numOfSplitFiles);
//        }
//
//        entryInfo.redundancy = (hasSeenOtherStrains > 0 || entryInfo.redundancy);
//        getDiffIdx(lastWrittenKmer, entryKmer, mergedDiffFile, diffBuffer, diffBufferIdx, totalBufferIdx);
//        lastWrittenKmer = entryKmer;
//        writeInfo(&entryInfo, mergedInfoFile, infoBuffer, infoBufferIdx, totalInfoIdx);
//        writtenKmerCnt++;
//
//        if(AminoAcid(lastWrittenKmer) != AAofTempSplitOffset && splitCheck == 1){
//            splitList[splitListIdx++] = {lastWrittenKmer, totalBufferIdx, totalInfoIdx};
//            splitCheck = 0;
//        }
//
//        if(writtenKmerCnt == offsetList[offsetListIdx]){
//            AAofTempSplitOffset = AminoAcid(lastWrittenKmer);
//            splitCheck = 1;
//            offsetListIdx++;
//        }
//
//        if(endFlag == 1) break;
//    }
//
//    cre->flushInfoBuf(infoBuffer, mergedInfoFile, infoBufferIdx);
//    cre->flushKmerBuf(diffBuffer, mergedDiffFile, diffBufferIdx);
//    fwrite(splitList, sizeof(DiffIdxSplit), SplitNum, diffIdxSplitFile);
//    for(int i = 0; i < SplitNum; i++){
//        cout<<splitList[i].ADkmer<< " "<<splitList[i].diffIdxOffset<< " "<<splitList[i].infoIdxOffset<<endl;
//    }
//    free(diffBuffer);
//    free(infoBuffer);
//    fclose(mergedDiffFile);
//    fclose(mergedInfoFile);
//    fclose(diffIdxSplitFile);
//
//    for(size_t file = 0; file < numOfSplitFiles; file++){
//        munmap(diffFileList[file].data, diffFileList[file].fileSize + 1);
//        munmap(infoFileList[file].data, infoFileList[file].fileSize + 1);
//    }
//    cout<<"Creating target DB is done"<<endl;
//    cout<<"Total k-mer count    : " << numOfKmerBeforeMerge <<endl;
//    cout<<"Written k-mer count  : " << writtenKmerCnt << endl;
//
//    delete[] diffFileList;
//    delete[] infoFileList;
//    delete[] lookingInfos;
//    delete[] lookingKmers;
//    delete[] diffFileIdx;
//    delete[] infoFileIdx;
//}


// Merge differential index and k-mer information files, reducing redundancy
void FileMerger::mergeTargetFiles(const LocalParameters & par, int numOfSplits) {
    size_t writtenKmerCnt = 0;
    const string dbDirectory = par.filenames[0];

    // Taxonomy
    NcbiTaxonomy taxonomy(par.taxonomyPath + "/names.dmp",
                          par.taxonomyPath + "/nodes.dmp",
                          par.taxonomyPath + "/merged.dmp");

    // Load taxonomy ids
    FILE * taxIdFile;
    if((taxIdFile = fopen((string(dbDirectory) + "/taxID_list").c_str(),"r")) == NULL){
        cout<<"Cannot open the taxID list file."<<endl;
        return;
    }
    char taxID[100];
    unordered_map<TaxID, TaxID> taxId2speciesId;
    while(feof(taxIdFile) == 0)
    {
        fscanf(taxIdFile,"%s",taxID);
        TaxID taxId = atol(taxID);
        TaxonNode const * taxon = taxonomy.taxonNode(taxId);
        if (taxId == taxon->taxId){
            TaxID speciesTaxID = taxonomy.getTaxIdAtRank(taxId, "species");
            while (taxon->taxId != speciesTaxID) {
                taxId2speciesId[taxon->taxId] = speciesTaxID;
                taxon = taxonomy.taxonNode(taxon->parentTaxId);
            }
            taxId2speciesId[speciesTaxID] = speciesTaxID;
        } else { // merged
            TaxID speciesTaxID = taxonomy.getTaxIdAtRank(taxId, "species");
            while (taxon->taxId != speciesTaxID) {
                taxId2speciesId[taxon->taxId] = speciesTaxID;
                taxon = taxonomy.taxonNode(taxon->parentTaxId);
            }
            taxId2speciesId[speciesTaxID] = speciesTaxID;
            taxId2speciesId[taxId] = speciesTaxID;
        }
    }
    fclose(taxIdFile);

    // File names for the final DB
    string mergedDiffFileName = dbDirectory + "/diffIdx";
    string mergedInfoFileName = dbDirectory + "/info";
    string diffIdxSplitFileName = dbDirectory + "/split";

    // Files to write
    FILE * mergedDiffFile = fopen(mergedDiffFileName.c_str(), "wb");
    FILE * mergedInfoFile = fopen(mergedInfoFileName.c_str(), "wb");
    FILE * diffIdxSplitFile = fopen(diffIdxSplitFileName.c_str(), "wb");

    // Buffers to fill
    uint16_t * diffBuffer = (uint16_t *)malloc(sizeof(uint16_t) * kmerBufSize);
    size_t diffBufferIdx = 0;
    size_t totalBufferIdx = 0;
    TargetKmerInfo * infoBuffer = (TargetKmerInfo *)malloc(sizeof(TargetKmerInfo) * kmerBufSize);
    size_t infoBufferIdx = 0;
    size_t totalInfoIdx = 0;

    // Prepare files to merge
    size_t numOfKmerBeforeMerge = 0;
    uint64_t * lookingKmers = new uint64_t[numOfSplits];
    auto * lookingInfos = new TargetKmerInfo[numOfSplits];
    auto * diffFileIdx = new size_t[numOfSplits];
    memset(diffFileIdx, 0, numOfSplits * sizeof(size_t));
    auto * infoFileIdx = new size_t[numOfSplits];
    memset(infoFileIdx, 0, numOfSplits * sizeof(size_t));
    auto * maxIdxOfEachFiles = new size_t[numOfSplits];
    struct MmapedData<uint16_t> *diffFileList = new struct MmapedData<uint16_t>[numOfSplits];
    struct MmapedData<TargetKmerInfo> *infoFileList = new struct MmapedData<TargetKmerInfo>[numOfSplits];
    for (int file = 0; file < numOfSplits; file++) {
        diffFileList[file] = mmapData<uint16_t>((dbDirectory + "/" + to_string(file) + "_diffIdx").c_str());
        infoFileList[file] = mmapData<TargetKmerInfo>((dbDirectory + "/" + to_string(file) + "_info").c_str());
        maxIdxOfEachFiles[file] = diffFileList[file].fileSize / sizeof(uint16_t);
        numOfKmerBeforeMerge += infoFileList[file].fileSize / sizeof(TargetKmerInfo);
    }

    // To make differential index splits
    uint64_t AAofTempSplitOffset = UINT64_MAX;
    size_t sizeOfSplit = numOfKmerBeforeMerge / (splitNum - 1);
    size_t offsetList[splitNum + 1];
    int offsetListIdx = 1;
    for(size_t os = 0; os < splitNum; os++){
        offsetList[os] = os * sizeOfSplit;
    }
    offsetList[splitNum] = UINT64_MAX;
    DiffIdxSplit splitList[splitNum];
    memset(splitList, 0, sizeof(DiffIdxSplit) * splitNum);
    int splitListIdx = 1;

    // get the first k-mer to write
    for(size_t file = 0; file < numOfSplits; file++){
        lookingKmers[file] = getNextKmer(0, diffFileList[file], diffFileIdx[file]);
        lookingInfos[file] = infoFileList[file].data[0];
        infoFileIdx[file] ++;
    }


    size_t idxOfMin = smallest(lookingKmers, lookingInfos, taxId2speciesId, numOfSplits);
    uint64_t lastWrittenKmer = 0;
    uint64_t entryKmer = lookingKmers[idxOfMin];
    TargetKmerInfo entryInfo = lookingInfos[idxOfMin];

    // Write first k-mer
    getDiffIdx(lastWrittenKmer, entryKmer, mergedDiffFile, diffBuffer, diffBufferIdx, totalBufferIdx);
    lastWrittenKmer = entryKmer;
    writeInfo(&entryInfo, mergedInfoFile, infoBuffer, infoBufferIdx, totalInfoIdx);
    writtenKmerCnt++;
    int splitCheck = 0;
    int endFlag = 0;

    size_t numOfincompletedFiles = numOfSplits;
    vector<TaxID> taxIds;
    while(true){
        // update entry k-mer
        entryKmer = lookingKmers[idxOfMin];
        entryInfo = lookingInfos[idxOfMin];

        // update looking k-mers
        lookingKmers[idxOfMin] = getNextKmer(entryKmer, diffFileList[idxOfMin], diffFileIdx[idxOfMin]);
        lookingInfos[idxOfMin] = infoFileList[idxOfMin].data[infoFileIdx[idxOfMin]];
        infoFileIdx[idxOfMin] ++;
        if( diffFileIdx[idxOfMin] >= maxIdxOfEachFiles[idxOfMin] ){
            lookingKmers[idxOfMin] = UINT64_MAX;
            numOfincompletedFiles--;
            if(numOfincompletedFiles == 0) break;
        }
        idxOfMin = smallest(lookingKmers, lookingInfos, taxId2speciesId, numOfSplits);

        int hasSeenOtherStrains = 0;
        taxIds.clear();
        taxIds.push_back(entryInfo.sequenceID); // Wrong
        while(taxId2speciesId[entryInfo.sequenceID] == taxId2speciesId[lookingInfos[idxOfMin].sequenceID]){
            if(entryKmer != lookingKmers[idxOfMin]) break;

            hasSeenOtherStrains += (entryInfo.sequenceID != lookingInfos[idxOfMin].sequenceID);
            taxIds.push_back(lookingInfos[idxOfMin].sequenceID);

            lookingKmers[idxOfMin] = getNextKmer(entryKmer, diffFileList[idxOfMin], diffFileIdx[idxOfMin]);
            lookingInfos[idxOfMin] = infoFileList[idxOfMin].data[infoFileIdx[idxOfMin]];
            infoFileIdx[idxOfMin] ++;

            if(diffFileIdx[idxOfMin] >= maxIdxOfEachFiles[idxOfMin] ){
                lookingKmers[idxOfMin] = UINT64_MAX;
                numOfincompletedFiles--;
                if(numOfincompletedFiles == 0){
                    endFlag = 1;
                    break;
                }
            }
            idxOfMin = smallest(lookingKmers, lookingInfos, taxId2speciesId, numOfSplits);
        }

        if (taxIds.size() > 1) {
            entryInfo.sequenceID = taxonomy.LCA(taxIds)->taxId;
        } else {
            entryInfo.sequenceID = taxIds[0];
        }

        entryInfo.redundancy = (hasSeenOtherStrains > 0 || entryInfo.redundancy);
        getDiffIdx(lastWrittenKmer, entryKmer, mergedDiffFile, diffBuffer, diffBufferIdx, totalBufferIdx);
        lastWrittenKmer = entryKmer;
        writeInfo(&entryInfo, mergedInfoFile, infoBuffer, infoBufferIdx, totalInfoIdx);
        writtenKmerCnt++;

        if(AminoAcidPart(lastWrittenKmer) != AAofTempSplitOffset && splitCheck == 1){
            splitList[splitListIdx++] = {lastWrittenKmer, totalBufferIdx, totalInfoIdx};
            cout<<splitListIdx<< "\t" <<lastWrittenKmer<<"\t"<<totalBufferIdx<<"\t"<<totalInfoIdx<<endl;
            splitCheck = 0;
        }

        if(writtenKmerCnt == offsetList[offsetListIdx]){
            AAofTempSplitOffset = AminoAcidPart(lastWrittenKmer);
            splitCheck = 1;
            offsetListIdx++;
        }

        if(endFlag == 1) break;
    }

    IndexCreator::flushInfoBuf(infoBuffer, mergedInfoFile, infoBufferIdx);
    IndexCreator::flushKmerBuf(diffBuffer, mergedDiffFile, diffBufferIdx);
    fwrite(splitList, sizeof(DiffIdxSplit), splitNum, diffIdxSplitFile);
    for(int i = 0; i < splitNum; i++){
        cout<<splitList[i].ADkmer<< " "<<splitList[i].diffIdxOffset<< " "<<splitList[i].infoIdxOffset<<endl;
    }
    free(diffBuffer);
    free(infoBuffer);
    fclose(mergedDiffFile);
    fclose(mergedInfoFile);
    fclose(diffIdxSplitFile);

    for(size_t file = 0; file < numOfSplits; file++){
        munmap(diffFileList[file].data, diffFileList[file].fileSize + 1);
        munmap(infoFileList[file].data, infoFileList[file].fileSize + 1);
    }
    cout<<"Creating target DB is done"<<endl;
    cout<<"Total k-mer count    : " << numOfKmerBeforeMerge <<endl;
    cout<<"Written k-mer count  : " << writtenKmerCnt << endl;

    delete[] diffFileList;
    delete[] infoFileList;
    delete[] lookingInfos;
    delete[] lookingKmers;
    delete[] diffFileIdx;
    delete[] infoFileIdx;
    delete[] maxIdxOfEachFiles;

    cout<<"done"<<endl;

    cout<<"Reference DB files you need are as below"<<endl;
    cout<<mergedDiffFileName<<endl;
    cout<<mergedInfoFileName<<endl;
    cout<<string(dbDirectory) + "/taxID_list"<<endl;
    cout<<diffIdxSplitFileName<<endl;
}



///It updates target database. Most of this function is the same with 'mergeTargetFiles', only some additional tasks are added to handle seqID
///It is not tested yet 2020.01.28
//void FileMerger::updateTargetDatabase(std::vector<char*> diffIdxFileNames, std::vector<char*> infoFileNames, vector<int> & taxListAtRank, vector<int> & taxIdList, const int & seqIdOffset) {
//    size_t totalKmerCnt = 0;
//    size_t writtenKmerCnt = 0;
//
//    ///Files to write on & buffers to fill them
//    FILE * mergedDiffFile = fopen(mergedDiffFileName, "wb");
//    FILE * mergedInfoFile = fopen(mergedInfoFileName, "wb");
//    uint16_t * diffBuffer = (uint16_t *)malloc(sizeof(uint16_t) * kmerBufSize); size_t diffBufferIdx = 0;
//    TargetKmerInfo * infoBuffer = (TargetKmerInfo *)malloc(sizeof(TargetKmerInfo) * kmerBufSize); size_t infoBufferIdx = 0;
//
//
//    ///Prepare files to merge
//    size_t numOfSplitFiles = diffIdxFileNames.size();
//    size_t numOfincompletedFiles = numOfSplitFiles;
//    uint64_t lookingKmers[numOfSplitFiles];
//    TargetKmerInfo lookingInfos[numOfSplitFiles];
//    size_t diffFileIdx[numOfSplitFiles];
//    memset(diffFileIdx, 0, sizeof(diffFileIdx));
//    size_t infoFileIdx[numOfSplitFiles];
//    memset(infoFileIdx, 0, sizeof(infoFileIdx));
//    size_t maxIdxOfEachFiles[numOfSplitFiles];
//    struct MmapedData<uint16_t> *diffFileList = new struct MmapedData<uint16_t>[numOfSplitFiles];
//    struct MmapedData<TargetKmerInfo> *infoFileList = new struct MmapedData<TargetKmerInfo>[numOfSplitFiles];
//    for (size_t file = 0; file < numOfSplitFiles; file++) {
//        diffFileList[file] = mmapData<uint16_t>(diffIdxFileNames[file]);
//        infoFileList[file] = mmapData<TargetKmerInfo>(infoFileNames[file]);
//        maxIdxOfEachFiles[file] = diffFileList[file].fileSize / sizeof(uint16_t);
//    }
//
//    ///TODO check the for loop condition; size
//    ///Fix sequence IDs of new k-mers
//    for (int i = 1; i < numOfSplitFiles; i++){
//        size_t size = (infoFileList[i].fileSize - 1) / sizeof(TargetKmerInfo);
//        for(int j = 0; j < size + 1; j++){
//            infoFileList[i].data[j].sequenceID += seqIdOffset;
//        }
//    }
//
//    /// get the first k-mer to write
//    for(size_t file = 0; file < numOfSplitFiles; file++){
//        lookingKmers[file] = getNextKmer(0, diffFileList[file], diffFileIdx[file]);
//        lookingInfos[file] = infoFileList[file].data[0];
//        infoFileIdx[file] ++;
//    }
//
//    size_t idxOfMin = smallest(lookingKmers, lookingInfos, taxListAtRank, numOfSplitFiles);
//    uint64_t lastWrittenKmer = 0;
//    uint64_t lastKmer = lookingKmers[idxOfMin];
//    TargetKmerInfo lastInfo = lookingInfos[idxOfMin];
//
//    ///write first k-mer
//    cre->getDiffIdx(0, lastKmer, mergedDiffFile, diffBuffer, diffBufferIdx);
//    lastWrittenKmer = lastKmer;
//    cre->writeInfo(&lastInfo, mergedInfoFile, infoBuffer, infoBufferIdx);;
//    writtenKmerCnt++;
//
//    int endFlag = 0;
//    int hasSeenOtherStrains;
//    int hasSeenOtherGenus;
//    ///끝부분 잘 되는지 확인할 것
//    while(true){
//        ///update last k-mer
//        lastKmer = lookingKmers[idxOfMin];
//        lastInfo = lookingInfos[idxOfMin];
//
//        ///update looking k-mers
//        lookingKmers[idxOfMin] = getNextKmer(lastKmer, diffFileList[idxOfMin], diffFileIdx[idxOfMin]);
//        lookingInfos[idxOfMin] = infoFileList[idxOfMin].data[infoFileIdx[idxOfMin]];
//        infoFileIdx[idxOfMin] ++;
//        if( diffFileIdx[idxOfMin] > maxIdxOfEachFiles[idxOfMin] ){
//            lookingKmers[idxOfMin] = UINT64_MAX;
//            numOfincompletedFiles--;
//            if(numOfincompletedFiles == 0) break;
//        }
//        totalKmerCnt ++;
//        idxOfMin = smallest(lookingKmers, lookingInfos, taxListAtRank, numOfSplitFiles);
//
//        hasSeenOtherStrains = 0;
//        hasSeenOtherGenus = 0;
//
//        while(taxListAtRank[lastInfo.sequenceID] == taxListAtRank[lookingInfos[idxOfMin].sequenceID]){
//            if(lastKmer != lookingKmers[idxOfMin]) break;
//
//            hasSeenOtherStrains += (taxIdList[lastInfo.sequenceID] != taxIdList[lookingInfos[idxOfMin].sequenceID]);
//
//            lookingKmers[idxOfMin] = getNextKmer(lastKmer, diffFileList[idxOfMin], diffFileIdx[idxOfMin]);
//            lookingInfos[idxOfMin] = infoFileList[idxOfMin].data[infoFileIdx[idxOfMin]];
//            infoFileIdx[idxOfMin] ++;
//            totalKmerCnt ++;
//
//            if(diffFileIdx[idxOfMin] > maxIdxOfEachFiles[idxOfMin]){
//                lookingKmers[idxOfMin] = UINT64_MAX;
//                numOfincompletedFiles--;
//                if(numOfincompletedFiles == 0){
//                    endFlag = 1;
//                    break;
//                }
//            }
//            idxOfMin = smallest(lookingKmers, lookingInfos, taxListAtRank, numOfSplitFiles);
//        }
//
//        lastInfo.redundancy = (hasSeenOtherStrains > 0 || lastInfo.redundancy);
//        cre->getDiffIdx(lastWrittenKmer, lastKmer, mergedDiffFile, diffBuffer, diffBufferIdx);
//        lastWrittenKmer = lastKmer;
//        cre->writeInfo(&lastInfo, mergedInfoFile, infoBuffer, infoBufferIdx);
//        writtenKmerCnt++;
//
//        if(endFlag == 1) break;
//    }
//
//    cre->flushInfoBuf(infoBuffer, mergedInfoFile, infoBufferIdx);
//    cre->flushKmerBuf(diffBuffer, mergedDiffFile, diffBufferIdx);
//
//    free(diffBuffer);
//    free(infoBuffer);
//    fclose(mergedDiffFile);
//    fclose(mergedInfoFile);
//    for(size_t file = 0; file < numOfSplitFiles; file++){
//        munmap(diffFileList[file].data, diffFileList[file].fileSize + 1);
//        munmap(infoFileList[file].data, infoFileList[file].fileSize + 1);
//    }
//
//    cout<<"Creating target DB is done"<<endl;
//    cout << "Total k-mer count    : " << totalKmerCnt << endl;
//    cout<<"Written k-mer count  : " << writtenKmerCnt << endl;
//}

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

size_t FileMerger::smallest(const uint64_t lookingKmers[],
                            const TargetKmerInfo lookingInfos[],
                            const unordered_map<TaxID, TaxID> & taxId2speciesId,
                            const size_t & fileCnt)
{
    size_t idxOfMin = 0;
    uint64_t min = lookingKmers[0];
    int minTaxIdAtRank = taxId2speciesId.at((int) lookingInfos[0].sequenceID);
    for(size_t i = 1; i < fileCnt; i++)
    {
        if(lookingKmers[i] < min ||
          (lookingKmers[i] == min && taxId2speciesId.at((int) lookingInfos[i].sequenceID) < minTaxIdAtRank)){
            min = lookingKmers[i];
            minTaxIdAtRank = taxId2speciesId.at((int) lookingInfos[i].sequenceID);
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