#include "FileMerger.h"
#include "common.h"

FileMerger::FileMerger(const LocalParameters & par) : par(par) {
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
    taxonomy = loadTaxonomy(dbDir, "");
}

FileMerger::~FileMerger() {
    delete taxonomy;
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
    // Load taxonomy ids
    FILE * taxIdFile;
    if((taxIdFile = fopen(taxIdListFileName.c_str(),"r")) == NULL){
        cout << "Cannot open the taxID list file: " << taxIdListFileName << endl;
        return;
    }

    char taxID[100];
    while(feof(taxIdFile) == 0) {
        fscanf(taxIdFile,"%s",taxID);
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
    fclose(taxIdFile);
}

void FileMerger::printTaxIdList(const string & taxIdListFileName) {
    FILE * taxidListFile = fopen(taxIdListFileName.c_str(), "w");
    for (auto & taxid: taxId2speciesId) {
        fprintf(taxidListFile, "%d\n", taxid.first);
    }
    

    FILE * taxIdFile;
    if((taxIdFile = fopen(taxIdListFileName.c_str(),"r")) == NULL){
        cout << "Cannot open the taxID list file: " << taxIdListFileName << endl;
        return;
    }

    char taxID[100];
    while(feof(taxIdFile) == 0) {
        fscanf(taxIdFile,"%s",taxID);
        cout << taxID << endl;
    }
    fclose(taxIdFile);
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
    for (int file = 0; file < numOfSplits; file++) {
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
    for(size_t os = 0; os < splitNum; os++){
        offsetList[os] = os * sizeOfSplit;
        // cout << os * sizeOfSplit << endl;
    }
    offsetList[splitNum] = UINT64_MAX;
    DiffIdxSplit splitList[splitNum];
    memset(splitList, 0, sizeof(DiffIdxSplit) * splitNum);
    int splitListIdx = 1;

    // get the first k-mer to write
    size_t lastFileIdx = numOfSplits - 1;
    unsigned int mask = ~((static_cast<unsigned int>(par.skipRedundancy == 0) << 31));
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
        writeInfo(&entryInfo, mergedInfoFile, infoBuffer, infoBufferIdx, totalInfoIdx, bufferSize);
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
        writeInfo(&entryInfo, mergedInfoFile, infoBuffer, infoBufferIdx, totalInfoIdx, bufferSize);
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
        // if (taxId2speciesId.find((int) lookingInfos[0].sequenceID) == taxId2speciesId.end()) {
        //     cerr << lookingKmers[0] << endl;
        //     cerr << "TaxID not found 1: " << lookingInfos[0].sequenceID << endl;
        //     exit(1);
        // }
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

void FileMerger::writeInfo(TaxID * entryToWrite, FILE * infoFile, TaxID * infoBuffer, size_t & infoBufferIdx, size_t & totalInfoIdx, size_t bufferSize)
{
    if (infoBufferIdx >= bufferSize) {
        flushInfoBuf(infoBuffer, infoFile, infoBufferIdx);
    }
    memcpy(infoBuffer + infoBufferIdx, entryToWrite, sizeof(TaxID));
    infoBufferIdx++;
    totalInfoIdx ++;
}
void FileMerger::flushInfoBuf(TaxID * buffer, FILE * infoFile, size_t & localBufIdx ) {
    fwrite(buffer, sizeof(TaxID), localBufIdx, infoFile);
    localBufIdx = 0;
}