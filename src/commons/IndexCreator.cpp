//
// Created by KJB on 01/09/2020.
//

#include "IndexCreator.h"

IndexCreator::IndexCreator()
{
    seqIterator = new SeqIterator();

}

IndexCreator::~IndexCreator() {
    delete seqIterator;
}

///It reads a reference sequence file and write a differential index file and target k-mer info file.
void IndexCreator::startIndexCreatingParallel(const char * seqFileName, const char * outputFileName, vector<int> & taxIdListAtRank)
{
    struct MmapedData<char> seqFile = mmapData<char>(seqFileName);

    vector<int> startsOfTaxIDs; ///start points where sequences of next species start
    startsOfTaxIDs.reserve(taxIdListAtRank.size());
    startsOfTaxIDs.push_back(0);
    vector<int> seqCntOfTaxIDs; ///the number of sequences for each species
    seqCntOfTaxIDs.reserve(taxIdListAtRank.size());
    int seqCnt = 0;
    int currentTaxID = taxIdListAtRank[0];
    for(size_t i = 0; i < taxIdListAtRank.size(); i++){
        if(taxIdListAtRank[i] == currentTaxID){
            seqCnt++;
        }
        else{
            seqCntOfTaxIDs.push_back(seqCnt);
            startsOfTaxIDs.push_back(i);
            currentTaxID = taxIdListAtRank[i];
            seqCnt = 1;
        }
    }
    seqCntOfTaxIDs.push_back(seqCnt);

    vector<Sequence> sequences;
    seqIterator->getSeqSegmentsWithHead(sequences, seqFile);

    size_t numOfTaxIDs = startsOfTaxIDs.size();
    bool processedSeqChecker[numOfTaxIDs];
    fill_n(processedSeqChecker, numOfTaxIDs, false);

    TargetKmerBuffer kmerBuffer(kmerBufSize);
    size_t processedTaxIDCnt = 0;

    while(processedTaxIDCnt < numOfTaxIDs){ ///check this condition
        seqIterator->fillTargetKmerBuffer2(kmerBuffer, seqFile, sequences, processedSeqChecker, processedTaxIDCnt, startsOfTaxIDs, seqCntOfTaxIDs);
        writeTargetFiles(kmerBuffer.buffer, kmerBuffer.startIndexOfReserve, outputFileName, taxIdListAtRank);
    }

    free(kmerBuffer.buffer);
    munmap(seqFile.data, seqFile.fileSize + 1);
}

///This function sort the TargetKmerBuffer, do redundancy reducing task, write the differential index of them
void IndexCreator::writeTargetFiles(TargetKmer * kmerBuffer, size_t & kmerNum, const char * outputFileName, const vector<int> & taxIdListAtRank)
{
    char suffixedDiffIdxFileName[100];
    char suffixedInfoFileName[100];
    sprintf(suffixedDiffIdxFileName,"%s_diffIdx_%zu", outputFileName,numOfFlush);
    sprintf(suffixedInfoFileName,"%s_info_%zu", outputFileName,numOfFlush);
    numOfFlush++;
    FILE * diffIdxFile = fopen(suffixedDiffIdxFileName, "wb");
    FILE * infoFile = fopen(suffixedInfoFileName, "wb");

    uint16_t *kmerLocalBuf = (uint16_t *)malloc(sizeof(uint16_t) * kmerBufSize);
    size_t localBufIdx = 0;

    if (diffIdxFile == NULL || infoFile == NULL){
        cout<<"Cannot open the file for writing target DB"<<endl;
        return;
    }
    uint64_t lastKmer = 0;

    SORT_PARALLEL(kmerBuffer, kmerBuffer + kmerNum, [](const TargetKmer & a, const TargetKmer & b) {
        return a.ADkmer < b.ADkmer || (a.ADkmer == b.ADkmer && a.info.sequenceID < b.info.sequenceID);
    });

    TargetKmer lookingKmer = kmerBuffer[0];
    size_t write = 0;
    int endFlag = 0;
    ///Redundancy reduction task
    for(size_t i = 1 ; i < kmerNum ; i++) {
        while(taxIdListAtRank[lookingKmer.info.sequenceID] == taxIdListAtRank[kmerBuffer[i].info.sequenceID]){
            if (lookingKmer.ADkmer != kmerBuffer[i].ADkmer) {
                break;
            }
            else{
                lookingKmer.info.redundancy = true;
            }
            i++;
            if(i == kmerNum){
                endFlag = 1;
                break;
            }
        }

        fwrite(&lookingKmer.info, sizeof(TargetKmerInfo), 1, infoFile);
        write++;
        writeKmerDiff(lastKmer, lookingKmer.ADkmer, diffIdxFile, kmerLocalBuf, localBufIdx);

        if(endFlag == 1) break;

        lastKmer = lookingKmer.ADkmer;
        lookingKmer = kmerBuffer[i];
    }
    //For the end part
    if(!((kmerBuffer[kmerNum - 2].ADkmer == kmerBuffer[kmerNum - 1].ADkmer) &&
         (kmerBuffer[kmerNum - 2].info.sequenceID == kmerBuffer[kmerNum - 1].info.sequenceID))){
        fwrite(&lookingKmer.info, sizeof(TargetKmerInfo), 1, infoFile);
        write++;
        writeKmerDiff(lastKmer, lookingKmer.ADkmer, diffIdxFile, kmerLocalBuf, localBufIdx);
    }

    cout << "total k-mer count: " << kmerNum << endl;
    cout<<"written k-mer count: "<<write<<endl;
    flushKmerBuf(kmerLocalBuf, diffIdxFile, localBufIdx);

    free(kmerLocalBuf);
    fclose(diffIdxFile);
    fclose(infoFile);
    kmerNum = 0;
}
void IndexCreator::writeKmerDiff(const uint64_t & lastKmer, const uint64_t & entryToWrite, FILE* handleKmerTable, uint16_t *kmerBuf, size_t & localBufIdx ){
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
    writeKmer(kmerBuf, handleKmerTable, (buffer + idx + 1), (4 - idx), localBufIdx);
}

void IndexCreator::flushKmerBuf(uint16_t *buffer, FILE *handleKmerTable, size_t & localBufIdx ) {
    fwrite(buffer, sizeof(uint16_t), localBufIdx, handleKmerTable);
    localBufIdx = 0;
}

void IndexCreator::writeKmer(uint16_t *buffer, FILE* handleKmerTable, uint16_t *toWrite, size_t size, size_t & localBufIdx ) {
    if (localBufIdx + size >= kmerBufSize) {
        flushKmerBuf(buffer, handleKmerTable, localBufIdx);
    }
    memcpy(buffer + localBufIdx, toWrite, sizeof(uint16_t) * size);
    localBufIdx += size;
}

void IndexCreator::writeInfo(TargetKmerInfo * entryToWrite, FILE * infoFile, TargetKmerInfo * infoBuffer, size_t & infoBufferIdx)
{
    if (infoBufferIdx >= kmerBufSize) {
        flushInfoBuf(infoBuffer, infoFile, infoBufferIdx);
    }
    memcpy(infoBuffer + infoBufferIdx, entryToWrite, sizeof(TargetKmerInfo));
    infoBufferIdx++;
}
void IndexCreator::flushInfoBuf(TargetKmerInfo * buffer, FILE * infoFile, size_t & localBufIdx ) {
    fwrite(buffer, sizeof(TargetKmerInfo), localBufIdx, infoFile);
    localBufIdx = 0;
}
int IndexCreator::getNumOfFlush()
{
    return numOfFlush;
}