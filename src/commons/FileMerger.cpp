#include "FileMerger.h"
#include "common.h"

FileMerger::FileMerger(const LocalParameters & par, const TaxonomyWrapper * taxonomy) : par(par), taxonomy(taxonomy) {
    // Load parameters
    removeRedundancyInfo = false;
    dbDir = par.filenames[0];
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

void FileMerger::addFilesToMerge(string diffIdxFileName, string infoFileName) {
    diffIdxFileNames.push_back(diffIdxFileName);
    infoFileNames.push_back(infoFileName);
}

void FileMerger::setMergedFileNames(string diffFileName, string infoFileName, string splitFileName) {
    mergedDiffFileName = diffFileName;
    mergedInfoFileName = infoFileName;
    diffIdxSplitFileName = splitFileName;
}

void FileMerger::updateTaxId2SpeciesTaxId(const string & taxIdListFileName) {
    FILE * taxIdFile;
    if((taxIdFile = fopen(taxIdListFileName.c_str(),"r")) == NULL){
        cout << "Cannot open the taxID list file: " << taxIdListFileName << endl;
        return;
    }

    char taxID[100];
    while(fscanf(taxIdFile,"%s",taxID) == 1) {
        TaxID taxId = atol(taxID);
        TaxonNode const * taxon = taxonomy->taxonNode(taxId);
        if (taxId == taxon->taxId){
            TaxID speciesTaxID = taxonomy->getTaxIdAtRank(taxId, "species");
            while (taxon->taxId != speciesTaxID) {
                taxId2speciesId[taxon->taxId] = speciesTaxID;
                taxon = taxonomy->taxonNode(taxon->parentTaxId);
            }
            taxId2speciesId[speciesTaxID] = speciesTaxID;
        } else { // merged
            TaxID speciesTaxID = taxonomy->getTaxIdAtRank(taxId, "species");
            while (taxon->taxId != speciesTaxID) {
                taxId2speciesId[taxon->taxId] = speciesTaxID;
                taxon = taxonomy->taxonNode(taxon->parentTaxId);
            }
            taxId2speciesId[speciesTaxID] = speciesTaxID;
            taxId2speciesId[taxId] = speciesTaxID;
        }
    }
    // Print taxId2speciesId
    // cout << "TaxID to speciesID" << endl;
    // for (auto & taxid: taxId2speciesId) {
    //     cout << "IN " << taxid.first << " " << taxid.second << endl;
    //     cout << "OR " << taxonomy->getOriginalTaxID(taxid.first) << " " << taxonomy->getOriginalTaxID(taxid.second) << endl;
    // }
    fclose(taxIdFile);
    Debug(Debug::INFO) << "Species-level taxonomy IDs are prepared.\n";
}

void FileMerger::printTaxIdList(const string & taxIdListFileName) {
    FILE * taxidListFile = fopen(taxIdListFileName.c_str(), "w");
    for (auto & taxid: taxId2speciesId) {
        fprintf(taxidListFile, "%d\n", taxid.first);
    }
    fclose(taxidListFile);
}

void FileMerger::mergeTargetFiles() {
    size_t writtenKmerCnt = 0;
   
    // Files to write
    FILE * mergedDiffFile = fopen(mergedDiffFileName.c_str(), "wb");
    FILE * mergedInfoFile = fopen(mergedInfoFileName.c_str(), "wb");
    FILE * diffIdxSplitFile = fopen(diffIdxSplitFileName.c_str(), "wb");

    // Buffers to fill
    size_t bufferSize = 1024 * 1024 * 32;
    uint16_t * diffBuffer = (uint16_t *)malloc(sizeof(uint16_t) * bufferSize);
    size_t diffBufferIdx = 0;
    size_t totalBufferIdx = 0;
    TaxID * infoBuffer = (TaxID *)malloc(sizeof(TaxID) * bufferSize);
    size_t infoBufferIdx = 0;
    size_t totalInfoIdx = 0;

    // Prepare files to merge
    size_t numOfKmerBeforeMerge = 0;
    size_t numOfSplits = diffIdxFileNames.size();
    uint64_t * lookingKmers = new uint64_t[numOfSplits];
    TaxID * lookingInfos = new TaxID[numOfSplits];
    size_t * diffFileIdx = new size_t[numOfSplits];
    memset(diffFileIdx, 0, numOfSplits * sizeof(size_t));
    size_t * infoFileIdx = new size_t[numOfSplits];
    memset(infoFileIdx, 0, numOfSplits * sizeof(size_t));
    size_t * maxIdxOfEachFiles = new size_t[numOfSplits];
    struct MmapedData<uint16_t> *diffFileList = new struct MmapedData<uint16_t>[numOfSplits];
    struct MmapedData<TaxID> *infoFileList = new struct MmapedData<TaxID>[numOfSplits];
    for (size_t file = 0; file < numOfSplits; file++) {
        diffFileList[file] = mmapData<uint16_t>(diffIdxFileNames[file].c_str());
        infoFileList[file] = mmapData<TaxID>(infoFileNames[file].c_str());
        maxIdxOfEachFiles[file] = diffFileList[file].fileSize / sizeof(uint16_t);
        numOfKmerBeforeMerge += infoFileList[file].fileSize / sizeof(TaxID);
    }

    // To make differential index splits
    uint64_t AAofTempSplitOffset = UINT64_MAX;
    size_t sizeOfSplit = numOfKmerBeforeMerge / (splitNum - 1);
    size_t offsetList[splitNum + 1];
    int offsetListIdx = 1;
    for(int os = 0; os < splitNum; os++){
        offsetList[os] = os * sizeOfSplit;
        // cout << os * sizeOfSplit << endl;
    }
    offsetList[splitNum] = UINT64_MAX;
    DiffIdxSplit splitList[splitNum];
    memset(splitList, 0, sizeof(DiffIdxSplit) * splitNum);
    int splitListIdx = 1;

    // get the first k-mer to write
    size_t lastFileIdx = numOfSplits - 1;
    unsigned int mask = ~((static_cast<unsigned int>(removeRedundancyInfo) << 31));
    for(size_t file = 0; file < numOfSplits; file++){
        lookingKmers[file] = getNextKmer(0, diffFileList[file], diffFileIdx[file]);
        if (removeRedundancyInfo && file == lastFileIdx) {
            lookingInfos[file] = infoFileList[file].data[0] & mask;
        } else {
            lookingInfos[file] = infoFileList[file].data[0];
        }
        infoFileIdx[file] ++;
    }

    size_t idxOfMin = smallest(lookingKmers, lookingInfos, numOfSplits);
    uint64_t lastWrittenKmer = 0;
    uint64_t entryKmer;
    TaxID entryInfo;
    int splitCheck = 0;
    size_t numOfincompletedFiles = numOfSplits;
    vector<TaxID> taxIds;
    int endPoint = 0;
    while(true){
        // update entry k-mer
        entryKmer = lookingKmers[idxOfMin];
        entryInfo = lookingInfos[idxOfMin];

        // update looking k-mers
        if (diffFileIdx[idxOfMin] >= maxIdxOfEachFiles[idxOfMin] ){
            lookingKmers[idxOfMin] = UINT64_MAX;
            numOfincompletedFiles--;
            if (numOfincompletedFiles == 0) {
                endPoint = 1;
                break;
            }
        } else {
            lookingKmers[idxOfMin] = getNextKmer(entryKmer, diffFileList[idxOfMin], diffFileIdx[idxOfMin]);
            if (removeRedundancyInfo && idxOfMin == lastFileIdx) {
                lookingInfos[idxOfMin] = infoFileList[idxOfMin].data[infoFileIdx[idxOfMin]] & mask;
            } else {
                lookingInfos[idxOfMin] = infoFileList[idxOfMin].data[infoFileIdx[idxOfMin]];
            }
            infoFileIdx[idxOfMin] ++;
        }

        idxOfMin = smallest(lookingKmers, lookingInfos, numOfSplits);

        taxIds.clear();
        taxIds.push_back(entryInfo);
        
        // Scan redundant k-mers
        while(taxId2speciesId[entryInfo] == taxId2speciesId[lookingInfos[idxOfMin]] && entryKmer == lookingKmers[idxOfMin]){
            taxIds.push_back(lookingInfos[idxOfMin]);
            if (diffFileIdx[idxOfMin] >= maxIdxOfEachFiles[idxOfMin]) {
                lookingKmers[idxOfMin] = UINT64_MAX;
                numOfincompletedFiles--;
                if (numOfincompletedFiles == 0) {
                    endPoint = 2;
                    break;
                }
            } else {
                lookingKmers[idxOfMin] = getNextKmer(entryKmer, diffFileList[idxOfMin], diffFileIdx[idxOfMin]);
                if (removeRedundancyInfo && idxOfMin == lastFileIdx) {
                    lookingInfos[idxOfMin] = infoFileList[idxOfMin].data[infoFileIdx[idxOfMin]] & mask;
                } else {
                    lookingInfos[idxOfMin] = infoFileList[idxOfMin].data[infoFileIdx[idxOfMin]];
                }
                infoFileIdx[idxOfMin] ++;
            }
            idxOfMin = smallest(lookingKmers, lookingInfos, numOfSplits);
        }

        if (taxIds.size() > 1) {
            entryInfo = taxonomy->LCA(taxIds)->taxId;
        } else {
            entryInfo = taxIds[0];
        }

        getDiffIdx(lastWrittenKmer, entryKmer, mergedDiffFile, diffBuffer, diffBufferIdx, totalBufferIdx, bufferSize);
        lastWrittenKmer = entryKmer;
        writeInfo(entryInfo, mergedInfoFile, infoBuffer, infoBufferIdx, totalInfoIdx, bufferSize);
        writtenKmerCnt++;

        if(AminoAcidPart(lastWrittenKmer) != AAofTempSplitOffset && splitCheck == 1){
            splitList[splitListIdx++] = {lastWrittenKmer, totalBufferIdx, totalInfoIdx};
            // cout << splitListIdx << "\t" << lastWrittenKmer << "\t" << totalBufferIdx << "\t" << totalInfoIdx <<endl;
            splitCheck = 0;
        }

        if(writtenKmerCnt == offsetList[offsetListIdx]){
            AAofTempSplitOffset = AminoAcidPart(lastWrittenKmer);
            splitCheck = 1;
            offsetListIdx++;
        }

        if(endPoint == 2) break;
    }

    if (endPoint == 1) {
        getDiffIdx(lastWrittenKmer, entryKmer, mergedDiffFile, diffBuffer, diffBufferIdx, totalBufferIdx, bufferSize);
        writeInfo(entryInfo, mergedInfoFile, infoBuffer, infoBufferIdx, totalInfoIdx, bufferSize);
        writtenKmerCnt++;
    }

    IndexCreator::flushInfoBuf(infoBuffer, mergedInfoFile, infoBufferIdx);
    IndexCreator::flushKmerBuf(diffBuffer, mergedDiffFile, diffBufferIdx);
    fwrite(splitList, sizeof(DiffIdxSplit), splitNum, diffIdxSplitFile);
    // for(int i = 0; i < splitNum; i++) {
    //     cout<<splitList[i].ADkmer<< " "<<splitList[i].diffIdxOffset<< " "<<splitList[i].infoIdxOffset<<endl;
    // }
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
    cout<<string(dbDir) + "/taxID_list"<<endl;
    cout<<diffIdxSplitFileName<<endl;
}


void FileMerger::mergeTargetFiles4() {
    size_t writtenKmerCnt = 0;
   
    // Files to write
    FILE * mergedDiffFile = fopen(mergedDiffFileName.c_str(), "wb");
    FILE * mergedInfoFile = fopen(mergedInfoFileName.c_str(), "wb");
    FILE * diffIdxSplitFile = fopen(diffIdxSplitFileName.c_str(), "wb");

    // Buffers to fill
    size_t bufferSize = 1024 * 1024 * 512;
    uint16_t * diffBuffer = (uint16_t *)malloc(sizeof(uint16_t) * bufferSize);
    size_t diffBufferIdx = 0;
    size_t totalBufferIdx = 0;
    TaxID * infoBuffer = (TaxID *)malloc(sizeof(TaxID) * bufferSize);
    size_t infoBufferIdx = 0;
    size_t totalInfoIdx = 0;

    // Prepare files to merge
    size_t numOfKmerBeforeMerge = 0;
    size_t numOfSplits = diffIdxFileNames.size();
    DeltaIdxReader ** deltaIdxReaders = new DeltaIdxReader*[numOfSplits];
    size_t valueBufferSize = 1024 * 1024 * 32;
    for (size_t file = 0; file < numOfSplits; file++) {
        deltaIdxReaders[file] = new DeltaIdxReader(diffIdxFileNames[file],
                                                   infoFileNames[file],
                                                   valueBufferSize, 
                                                   1024 * 1024);
        numOfKmerBeforeMerge += deltaIdxReaders[file]->getTotalValueNum();
    }

    // To make differential index splits
    uint64_t AAofTempSplitOffset = UINT64_MAX;
    size_t sizeOfSplit = numOfKmerBeforeMerge / (splitNum - 1);
    size_t offsetList[splitNum + 1];
    int offsetListIdx = 1;
    for(int os = 0; os < splitNum; os++){
        offsetList[os] = os * sizeOfSplit;
    }
    offsetList[splitNum] = UINT64_MAX;
    DiffIdxSplit splitList[splitNum];
    memset(splitList, 0, sizeof(DiffIdxSplit) * splitNum);
    int splitListIdx = 1;

    // get the first k-mer to write
    size_t lastFileIdx = numOfSplits - 1;
    unsigned int mask = ~((static_cast<unsigned int>(removeRedundancyInfo) << 31));
    uint64_t kmer;

    uint64_t lastWrittenKmer = 0;
    int splitCheck = 0;
    vector<TaxID> taxIds;
    bool haveLast = false;

    cout << "Merging target files..." << endl;
    Buffer<TargetKmer> metamerBuffer(1024 * 1024 * 1024);
    std::atomic<int> hasOverflow{0};
    std::vector<std::atomic<bool>> completedSplits(numOfSplits);
    int remainingSplits = numOfSplits;
    vector<pair<size_t, size_t>> uniqKmerIdxRanges;
    size_t write = 0;
    size_t lastKmer = 0;
    size_t totalValueNum = 0;
    while (remainingSplits > 0) {
        memset(metamerBuffer.buffer, 0, metamerBuffer.bufferSize * sizeof(TargetKmer));
        time_t start = time(nullptr);
        while (true) {
            size_t posToWrite;
            uint64_t max = UINT64_MAX;
            posToWrite = metamerBuffer.reserveMemory(valueBufferSize * numOfSplits);
            if (posToWrite + valueBufferSize * numOfSplits > metamerBuffer.bufferSize) {
                cout << "Buffer overflow, flushing..." << endl;
                metamerBuffer.startIndexOfReserve -= valueBufferSize * numOfSplits;
                break;
            }
            for (size_t i = 0; i < numOfSplits; i++) {
                if (completedSplits[i].load(std::memory_order_acquire))
                    continue; 
                if (deltaIdxReaders[i]->getLastValue() < max) {
                    max = deltaIdxReaders[i]->getLastValue();
                }
            }
            for (size_t split = 0; split < numOfSplits; split++) {
                if (completedSplits[split].load(std::memory_order_acquire))
                    continue; 
                size_t valueNum = deltaIdxReaders[split]->getValues(metamerBuffer.buffer + posToWrite + split * valueBufferSize, max);
                totalValueNum += valueNum;
                if (valueNum == 0 && deltaIdxReaders[split]->isCompleted()) {
                    completedSplits[split].store(true, std::memory_order_release);
                    remainingSplits--;
                    continue;
                }
                for (size_t i = 0; i < valueNum; i++) {
                    metamerBuffer.buffer[posToWrite + split * valueBufferSize + i].spTaxId 
                        = taxId2speciesId[metamerBuffer.buffer[posToWrite + split * valueBufferSize + i].metamer.id & mask];
                }    
            }
        }
        time_t end = time(nullptr);
        cout << "Time spent for reading k-mers: " << (double) (end - start) << endl;

        time_t beforeSort = time(nullptr);
        SORT_PARALLEL(metamerBuffer.buffer, metamerBuffer.buffer + metamerBuffer.startIndexOfReserve,
                      IndexCreator::compareForDiffIdx);
        time_t afterSort = time(nullptr);
        cout << "Sort time: " << afterSort - beforeSort << endl;

        // Reduce redundancy
        auto * uniqKmerIdx = new size_t[metamerBuffer.startIndexOfReserve + 1];
        size_t uniqKmerCnt = 0;
        uniqKmerIdxRanges.clear();
        reduceRedundancy(metamerBuffer, uniqKmerIdx, uniqKmerCnt, uniqKmerIdxRanges, par);
        time_t reduction = time(nullptr);
        cout << "Time spent for reducing redundancy: " << (double) (reduction - afterSort) << endl;

        // Write        
        for (size_t i = 0; i < uniqKmerIdxRanges.size(); i ++) {
            for (size_t j = uniqKmerIdxRanges[i].first; j < uniqKmerIdxRanges[i].second; j ++) {
                writeInfo(metamerBuffer.buffer[uniqKmerIdx[j]].metamer.id, mergedInfoFile, infoBuffer, infoBufferIdx, totalInfoIdx, bufferSize);
                write++;
                getDiffIdx(lastKmer, metamerBuffer.buffer[uniqKmerIdx[j]].metamer.metamer, mergedDiffFile, diffBuffer, diffBufferIdx, totalBufferIdx, bufferSize);
                lastKmer = metamerBuffer.buffer[uniqKmerIdx[j]].metamer.metamer;

                if (AminoAcidPart(lastWrittenKmer) != AAofTempSplitOffset && splitCheck == 1) {
                    splitList[splitListIdx++] = {lastWrittenKmer, totalBufferIdx, totalInfoIdx};
                    splitCheck = 0;
                }
                if (write == offsetList[offsetListIdx]) {
                    AAofTempSplitOffset = AminoAcidPart(lastWrittenKmer);
                    splitCheck = 1;
                    offsetListIdx++;
                }
            }
        }
        time_t writeTime = time(nullptr);
        cout << "Time spent for writing k-mers: " << (double) (writeTime - reduction) << endl;
        cout << "Written k-mers: " << write << endl;
        metamerBuffer.startIndexOfReserve = 0; // Reset the buffer for the next round
    }

    IndexCreator::flushInfoBuf(infoBuffer, mergedInfoFile, infoBufferIdx);
    IndexCreator::flushKmerBuf(diffBuffer, mergedDiffFile, diffBufferIdx);
    fwrite(splitList, sizeof(DiffIdxSplit), splitNum, diffIdxSplitFile);
    // for(int i = 0; i < splitNum; i++) {
    //     cout<<splitList[i].ADkmer<< " "<<splitList[i].diffIdxOffset<< " "<<splitList[i].infoIdxOffset<<endl;
    // }
    free(diffBuffer);
    free(infoBuffer);
    fclose(mergedDiffFile);
    fclose(mergedInfoFile);
    fclose(diffIdxSplitFile);

    cout<<"Creating target DB is done"<<endl;
    cout<<"Total k-mer count    : " << write <<endl;
    cout<<"done"<<endl;

    cout<<"Reference DB files you need are as below"<<endl;
    cout<<mergedDiffFileName<<endl;
    cout<<mergedInfoFileName<<endl;
    cout<<string(dbDir) + "/taxID_list"<<endl;
    cout<<diffIdxSplitFileName<<endl;
}

void FileMerger::reduceRedundancy(Buffer<TargetKmer> & kmerBuffer,
                                    size_t * uniqKmerIdx,
                                    size_t & uniqueKmerCnt,
                                    vector<pair<size_t, size_t>> & uniqKmerIdxRanges,
                                    const LocalParameters & par) {
    // Find the first index of garbage k-mer (UINT64_MAX)
    for(size_t checkN = kmerBuffer.startIndexOfReserve - 1; checkN != 0; checkN--){
        if(kmerBuffer.buffer[checkN].metamer.metamer != UINT64_MAX){
            kmerBuffer.startIndexOfReserve = checkN + 1;
            break;
        }
    }

    // Find the first index of meaningful k-mer
    size_t startIdx = 0;
    for(size_t i = 0; i < kmerBuffer.startIndexOfReserve ; i++){
        if(kmerBuffer.buffer[i].spTaxId != 0){
            startIdx = i;
            break;
        }
    }

    cout << "Number of k-mers before redundancy reduction: " << kmerBuffer.startIndexOfReserve - startIdx << endl;

    // Make splits
    vector<Split2> splits;
    size_t splitWidth = (kmerBuffer.startIndexOfReserve - startIdx) / par.threads;
    for (int i = 0; i < par.threads - 1; i++) {
        for (size_t j = startIdx + splitWidth; j + 1 < kmerBuffer.startIndexOfReserve; j++) {
            if (AminoAcidPart(kmerBuffer.buffer[j].metamer.metamer) != AminoAcidPart(kmerBuffer.buffer[j + 1].metamer.metamer)) {
                splits.emplace_back(startIdx, j);
                uniqKmerIdxRanges.emplace_back(pair<size_t, size_t>(startIdx, 0));
                startIdx = j + 1;
                break;
            }
        }
    }
    splits.emplace_back(startIdx, kmerBuffer.startIndexOfReserve - 1);
    uniqKmerIdxRanges.emplace_back(pair<size_t, size_t>(startIdx, 0));

    size_t * cntOfEachSplit = new size_t[splits.size()];
    for(size_t i = 0; i < splits.size(); i++) {
        cntOfEachSplit[i] = 0;
    }
#pragma omp parallel default(none), shared(kmerBuffer, cntOfEachSplit, splits, par, cout, uniqueKmerCnt, uniqKmerIdx, uniqKmerIdxRanges)
    {
        TargetKmer * lookingKmer;
        size_t lookingIndex;
        int endFlag;
        vector<TaxID> taxIds;
        size_t * tempUniqKmerIdx = new size_t[16 * 1024 * 1024];
        size_t tempUniqKmerCnt = 0;
#pragma omp for schedule(static, 1)
        for(size_t split = 0; split < splits.size(); split ++) {
            lookingKmer = & kmerBuffer.buffer[splits[split].offset];
            lookingIndex = splits[split].offset;
            endFlag = 0;
            for(size_t i = 1 + splits[split].offset; i < splits[split].end + 1 ; i++) {
                taxIds.clear();
                taxIds.push_back(lookingKmer->metamer.id);
                // Scan redundancy
                while(lookingKmer->spTaxId == kmerBuffer.buffer[i].spTaxId &&
                      lookingKmer->metamer.metamer == kmerBuffer.buffer[i].metamer.metamer){
                    taxIds.push_back(kmerBuffer.buffer[i].metamer.id);
                    i++;
                    if(i == splits[split].end + 1){
                        endFlag = 1;
                        break;
                    }
                }
                if(taxIds.size() > 1){
                    lookingKmer->metamer.id = taxonomy->LCA(taxIds)->taxId;
                } else {
                    lookingKmer->metamer.id = taxIds[0];
                }

                if (tempUniqKmerCnt >= 16 * 1024 * 1024) {
                    memcpy(uniqKmerIdx + splits[split].offset, tempUniqKmerIdx, tempUniqKmerCnt * sizeof(size_t));
                    splits[split].offset += tempUniqKmerCnt;
                    tempUniqKmerCnt = 0;
                }
                tempUniqKmerIdx[tempUniqKmerCnt++] = lookingIndex;
                cntOfEachSplit[split]++;
                if(endFlag == 1) break;
                lookingKmer = & kmerBuffer.buffer[i];
                lookingIndex = i;
            }

            //For the end part
            if(!((kmerBuffer.buffer[splits[split].end - 1].metamer.metamer == kmerBuffer.buffer[splits[split].end].metamer.metamer) &&
                 (kmerBuffer.buffer[splits[split].end - 1].spTaxId == kmerBuffer.buffer[splits[split].end].spTaxId))){
                if (tempUniqKmerCnt >= 16 * 1024 * 1024) {
                    memcpy(uniqKmerIdx + splits[split].offset, tempUniqKmerIdx, tempUniqKmerCnt * sizeof(size_t));
                    splits[split].offset += tempUniqKmerCnt;
                    tempUniqKmerCnt = 0;
                }
                tempUniqKmerIdx[tempUniqKmerCnt++] = splits[split].end;
                cntOfEachSplit[split]++;
            }
            memcpy(uniqKmerIdx + splits[split].offset, tempUniqKmerIdx, tempUniqKmerCnt * sizeof(size_t));
            splits[split].offset += tempUniqKmerCnt;
            uniqKmerIdxRanges[split].second = splits[split].offset;
            __sync_fetch_and_add(&uniqueKmerCnt, cntOfEachSplit[split]); 
        }
        delete[] tempUniqKmerIdx;
    }
    delete[] cntOfEachSplit;
}

void FileMerger::mergeDeltaIdxFiles() {
    size_t writtenKmerCnt = 0;
   
    // Files to write
    FILE * mergedDeltaIdxFile = fopen(mergedDiffFileName.c_str(), "wb");
    FILE * deltaIdxOffsetFile = fopen(diffIdxSplitFileName.c_str(), "wb");

    // Buffers to fill
    size_t bufferSize = 1024 * 1024 * 32;
    uint16_t * deltaIdxBuffer = (uint16_t *)malloc(sizeof(uint16_t) * bufferSize);
    size_t deltaBufferIdx = 0;
    size_t totalBufferIdx = 0;
    
    // Prepare files to merge
    size_t totalDeltaNum = 0;
    size_t numOfSplits = diffIdxFileNames.size();
    Metamer * lookingMetamers = new Metamer[numOfSplits];
    size_t * deltaFileIdx = new size_t[numOfSplits];
    memset(deltaFileIdx, 0, numOfSplits * sizeof(size_t));
    size_t * maxIdxOfEachFiles = new size_t[numOfSplits];
    struct MmapedData<uint16_t> *deltaFileList = new struct MmapedData<uint16_t>[numOfSplits];
    size_t numOfKmerBeforeMerge = 0;
    for (size_t file = 0; file < numOfSplits; file++) {
        deltaFileList[file] = mmapData<uint16_t>(diffIdxFileNames[file].c_str());
        maxIdxOfEachFiles[file] = deltaFileList[file].fileSize / sizeof(uint16_t);
        totalDeltaNum += maxIdxOfEachFiles[file];
    }

    // To make differential index splits
    uint64_t AAofTempSplitOffset = UINT64_MAX;
    size_t sizeOfSplit = totalDeltaNum / (splitNum - 1);
    size_t offsetList[splitNum + 1];
    int offsetListIdx = 1;
    for(int os = 0; os < splitNum; os++){
        offsetList[os] = os * sizeOfSplit;
    }
    offsetList[splitNum] = UINT64_MAX;
    DeltaIdxOffset splitList[splitNum];
    memset(splitList, 0, sizeof(DeltaIdxOffset) * splitNum);
    int splitListIdx = 1;

    // get the first k-mer to write
    size_t lastFileIdx = numOfSplits - 1;
    size_t asdf = 0;
    // SeqIterator * seqIterator = new SeqIterator(par); // DEBUG
    for(size_t file = 0; file < numOfSplits; file++){
        Metamer metamer;
        lookingMetamers[file] = KmerMatcher::getNextTargetKmer(metamer, deltaFileList[file].data, deltaFileIdx[file], asdf);
        numOfKmerBeforeMerge ++;
    }

    size_t idxOfMin = smallest(lookingMetamers, numOfSplits);
    Metamer lastWrittenMetamer;
    Metamer entryMetamer;
    int splitCheck = 0;
    size_t numOfincompletedFiles = numOfSplits;
    vector<Metamer> identicalMetamers;
    vector<Metamer> filteredMetamers;
    unordered_map<TaxID, TaxID> observedSpecies2LCA;
    int endPoint = 0;
    size_t totalFilteredMetamers = 0;
    while(true) {
        entryMetamer = lookingMetamers[idxOfMin];

        // Update looking k-mers
        if (unlikely(deltaFileIdx[idxOfMin] >= maxIdxOfEachFiles[idxOfMin])) {
            lookingMetamers[idxOfMin] = Metamer(UINT64_MAX, INT32_MAX);
            numOfincompletedFiles--;
            if (numOfincompletedFiles == 0) {
                endPoint = 1;
                break;
            }
        } else {
            lookingMetamers[idxOfMin] = KmerMatcher::getNextTargetKmer(entryMetamer, deltaFileList[idxOfMin].data, deltaFileIdx[idxOfMin], asdf);
            numOfKmerBeforeMerge ++;
        }

        idxOfMin = smallest(lookingMetamers, numOfSplits);
        
        // Collect identical k-mers
        identicalMetamers.clear();
        identicalMetamers.push_back(entryMetamer);
        while (entryMetamer.metamer == lookingMetamers[idxOfMin].metamer) {
            identicalMetamers.push_back(lookingMetamers[idxOfMin]);
            if (unlikely(deltaFileIdx[idxOfMin] >= maxIdxOfEachFiles[idxOfMin])) {
                lookingMetamers[idxOfMin] = Metamer(UINT64_MAX, INT32_MAX);
                numOfincompletedFiles--;
                if (numOfincompletedFiles == 0) {
                    endPoint = 2;
                    break;
                }
            } else {
                lookingMetamers[idxOfMin] = KmerMatcher::getNextTargetKmer(lookingMetamers[idxOfMin], deltaFileList[idxOfMin].data, deltaFileIdx[idxOfMin], asdf);
                numOfKmerBeforeMerge ++;
            }
            idxOfMin = smallest(lookingMetamers, numOfSplits);
        }

        // Remove duplicates wihin the same species
        filteredMetamers.clear();
        observedSpecies2LCA.clear();
        for (size_t i = 0; i < identicalMetamers.size(); i++) {
            if (taxId2speciesId[identicalMetamers[i].id] == 0) {
                cout << "TaxID not found 3: " << identicalMetamers[i].id << endl;
                continue;
            }
            if (observedSpecies2LCA.find(taxId2speciesId[identicalMetamers[i].id]) == observedSpecies2LCA.end()) {
                filteredMetamers.push_back(identicalMetamers[i]);
                observedSpecies2LCA[taxId2speciesId[identicalMetamers[i].id]] = identicalMetamers[i].id;
            } else {
                observedSpecies2LCA[taxId2speciesId[identicalMetamers[i].id]] = 
                    taxonomy->LCA(observedSpecies2LCA[taxId2speciesId[identicalMetamers[i].id]], identicalMetamers[i].id);
            }
        }

        // Label the LCA to the filteredMetamers
        for (size_t k = 0; k < filteredMetamers.size(); k++) {
            filteredMetamers[k].id = observedSpecies2LCA[taxId2speciesId[filteredMetamers[k].id]];
        }

        // Sort filtered Metamers
        sort(filteredMetamers.begin(), filteredMetamers.end(), IndexCreator::compareMetamerID);

        for (size_t k = 0; k < filteredMetamers.size(); k++) {
            IndexCreator::getDeltaIdx(lastWrittenMetamer, filteredMetamers[k], mergedDeltaIdxFile, deltaIdxBuffer, bufferSize, deltaBufferIdx, totalBufferIdx);
            lastWrittenMetamer = filteredMetamers[k];
            writtenKmerCnt++;
            if (AminoAcidPart(lastWrittenMetamer.metamer) != AAofTempSplitOffset && splitCheck == 1) {
                splitList[splitListIdx++] = {lastWrittenMetamer, totalBufferIdx};
                splitCheck = 0;
            }
            if(totalBufferIdx >= offsetList[offsetListIdx]) { 
                AAofTempSplitOffset = AminoAcidPart(lastWrittenMetamer.metamer);
                splitCheck = 1;
                offsetListIdx++;
            }
        }
        totalFilteredMetamers += filteredMetamers.size();
        if(endPoint == 2) break;
    }

    if (endPoint == 1) {
        IndexCreator::getDeltaIdx(lastWrittenMetamer, entryMetamer, mergedDeltaIdxFile, deltaIdxBuffer, bufferSize, deltaBufferIdx, totalBufferIdx);
        writtenKmerCnt++;
    }

    IndexCreator::flushKmerBuf(deltaIdxBuffer, mergedDeltaIdxFile, deltaBufferIdx);
    fwrite(splitList, sizeof(DeltaIdxOffset), splitNum, deltaIdxOffsetFile);

    free(deltaIdxBuffer);
    fclose(mergedDeltaIdxFile);
    fclose(deltaIdxOffsetFile);

    for(size_t file = 0; file < numOfSplits; file++){
        munmap(deltaFileList[file].data, deltaFileList[file].fileSize + 1);
    }
    cout<<"Creating target DB is done"<<endl;
    cout<<"Total k-mer count    : " << numOfKmerBeforeMerge <<endl;
    cout<<"Written k-mer count  : " << writtenKmerCnt << endl;
    cout << totalFilteredMetamers << endl;

    delete[] deltaFileList;
    delete[] lookingMetamers;
    delete[] deltaFileIdx;
    delete[] maxIdxOfEachFiles;

    cout<<"done"<<endl;

    cout<<"Reference DB files are as below"<<endl;
    cout<<mergedDiffFileName<<endl;
    cout<<string(dbDir) + "/taxID_list"<<endl;
    cout<<diffIdxSplitFileName<<endl;
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

size_t FileMerger::smallest(const uint64_t * lookingKmers,
                            const TaxID * lookingInfos,
                            size_t fileCnt) {
    size_t idxOfMin = 0;
    uint64_t min = lookingKmers[0];
    int minTaxIdAtRank;
    if (min == UINT64_MAX) {
        minTaxIdAtRank = INT_MAX;
    } else {
        minTaxIdAtRank = taxId2speciesId.at(lookingInfos[0]);
    }
    for(size_t i = 1; i < fileCnt; i++) {
        if(lookingKmers[i] < min ||
          (lookingKmers[i] == min && taxId2speciesId.at(lookingInfos[i]) < minTaxIdAtRank)){
            min = lookingKmers[i];
            minTaxIdAtRank = taxId2speciesId.at(lookingInfos[i]);
            idxOfMin = i;
        }
    }
    return idxOfMin;
}

size_t FileMerger::smallest(const Metamer * lookingKmers,
                            size_t fileCnt) {
    size_t idxOfMin = 0;
    uint64_t min = lookingKmers[0].metamer;
    int minTaxIdAtRank;
    if (min == UINT64_MAX) {
        minTaxIdAtRank = INT32_MAX;
    } else {
        minTaxIdAtRank = taxId2speciesId.at(lookingKmers[0].id);
    }
    for(size_t i = 1; i < fileCnt; i++) {
        if (lookingKmers[i].metamer == UINT64_MAX || lookingKmers[i].id == INT32_MAX) {
            continue;
        }
        if(lookingKmers[i].metamer < min ||
          (lookingKmers[i].metamer == min && taxId2speciesId.at(lookingKmers[i].id) < minTaxIdAtRank)){
            min = lookingKmers[i].metamer;
            minTaxIdAtRank = taxId2speciesId.at(lookingKmers[i].id);
            idxOfMin = i;
        }
    }
    return idxOfMin;
}

void FileMerger::getDiffIdx(const uint64_t & lastKmer, const uint64_t & entryToWrite, FILE* handleKmerTable, uint16_t *kmerBuf, size_t & localBufIdx, size_t & totalBufIdx, size_t bufferSize){
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
    writeDiffIdx(kmerBuf, handleKmerTable, (buffer + idx + 1), (4 - idx), localBufIdx, totalBufIdx, bufferSize);
}
void FileMerger::writeDiffIdx(uint16_t *buffer, FILE* handleKmerTable, uint16_t *toWrite, size_t size, size_t & localBufIdx, size_t & totalBufIdx, size_t bufferSize) {
    if (localBufIdx + size >= bufferSize) {
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

void FileMerger::writeInfo(TaxID entryToWrite, FILE * infoFile, TaxID * infoBuffer, size_t & infoBufferIdx, size_t & totalInfoIdx, size_t bufferSize)
{
    if (infoBufferIdx >= bufferSize) {
        flushInfoBuf(infoBuffer, infoFile, infoBufferIdx);
    }
    infoBuffer[infoBufferIdx] = entryToWrite;
    infoBufferIdx++;
    totalInfoIdx ++;
}
void FileMerger::flushInfoBuf(TaxID * buffer, FILE * infoFile, size_t & localBufIdx ) {
    fwrite(buffer, sizeof(TaxID), localBufIdx, infoFile);
    localBufIdx = 0;
}

void FileMerger::printFilesToMerge() {
    cout << "Files to merge :" << endl;
    for (size_t i = 0; i < diffIdxFileNames.size(); i++) {
        cout << diffIdxFileNames[i] << " " << infoFileNames[i] << endl;
    }
    // cout << "Merged files :" << endl;
    // cout << mergedDiffFileName << endl;
    // cout << mergedInfoFileName << endl;
}