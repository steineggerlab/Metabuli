//
// Created by KJB on 01/09/2020.
//

#include "IndexCreator.h"

IndexCreator::IndexCreator()
{
    seqAlterator = new SeqAlterator();

}

IndexCreator::~IndexCreator() {
    delete seqAlterator;
}

void IndexCreator::startIndexCreating2(const char * seqFileName, const char * outputFileName, vector<int> & taxIdListAtRank)
{
    ExtractStartPoint ESP = {0 ,0};
    size_t bufferIdx = 0;

    struct MmapedData<char> seqFile = mmapData<char>(seqFileName);
    size_t numOfChar = seqFile.fileSize / sizeof(char);


    vector<SeqSegment> seqSegments;
    seqAlterator -> getSeqSegments(seqSegments, seqFile);


    Kmer * kmerBuffer = (Kmer *)malloc(sizeof(Kmer) * kmerBufSize);

    int seqID = 0;
    for(size_t i = 1 ; i < seqSegments.size(); i++)
    {
        seqAlterator -> dna2aa2(seqSegments[i], seqFile);
        ESP = seqAlterator->fillKmerBuffer2(seqSegments[i], seqFile, kmerBuffer, seqID, bufferIdx, ESP);
        while (ESP.startOfFrame + ESP.frame != 0)
        {
            writeTargetFiles(kmerBuffer, bufferIdx, outputFileName, taxIdListAtRank);
            ESP = seqAlterator->fillKmerBuffer2(seqSegments[i], seqFile, kmerBuffer, seqID, bufferIdx, ESP);
        }
        seqID ++;
    }

    //flush last buffer
    writeTargetFiles(kmerBuffer, bufferIdx, outputFileName, taxIdListAtRank);

    free(kmerBuffer);
    munmap(seqFile.data, seqFile.fileSize + 1);
}


void IndexCreator::startIndexCreating(const char * seqFileName, const char * outputFileName, vector<int> & taxIdListAtRank)
{
    ExtractStartPoint ESP = {0 ,0};
    size_t bufferIdx = 0;

    struct MmapedData<char> seqFile = mmapData<char>(seqFileName);
    size_t numOfChar = seqFile.fileSize / sizeof(char);

    kseq_buffer_t buffer(const_cast<char*>(seqFile.data), seqFile.fileSize + 1);
    kseq_t *seq = kseq_init(&buffer);

    Kmer * kmerBuffer = (Kmer *)malloc(sizeof(Kmer) * kmerBufSize);

    int seqID = 0;
    while(kseq_read(seq) >= 0){
        seqAlterator->dna2aa(seq->seq.s);
        ESP = seqAlterator->fillKmerBuffer(seq->seq.s, kmerBuffer, seqID, bufferIdx, ESP);
        while (ESP.startOfFrame + ESP.frame != 0)
        {
            writeTargetFiles(kmerBuffer, bufferIdx, outputFileName, taxIdListAtRank);
            ESP = seqAlterator->fillKmerBuffer(seq->seq.s, kmerBuffer, seqID, bufferIdx, ESP);
        }
        seqID ++;
    }

    //flush last buffer
    writeTargetFiles(kmerBuffer, bufferIdx, outputFileName, taxIdListAtRank);

    free(kmerBuffer);
    munmap(seqFile.data, seqFile.fileSize + 1);
}

void IndexCreator::writeTargetFiles(Kmer *kmerBuffer, size_t & bufferIdx, const char * outputFileName, vector<int> & taxIdListAtRank)
{
    char suffixedDiffIdxFileName[100];
    char suffixedInfoFileName[100];
    sprintf(suffixedDiffIdxFileName,"%s_diffIdx_%zu", outputFileName,numOfFlush);
    sprintf(suffixedInfoFileName,"%s_info_%zu", outputFileName,numOfFlush);
    numOfFlush++;
    FILE * diffIdxFile = fopen(suffixedDiffIdxFileName, "wb");
    FILE * idAndPosFile = fopen(suffixedInfoFileName, "wb");


    uint16_t *kmerLocalBuf = (uint16_t *)malloc(sizeof(uint16_t) * kmerBufSize);
    size_t localBufIdx = 0;

    if (diffIdxFile == NULL || idAndPosFile == NULL){
        cout<<"Cannot open the file for writing target DB"<<endl;
        return;
    }
    uint64_t lastKmer = 0;

    sort(kmerBuffer, kmerBuffer + bufferIdx, [](const Kmer & a, const Kmer & b) {
        return a.ADkmer < b.ADkmer || (a.ADkmer == b.ADkmer && a.info.sequenceID < b.info.sequenceID);
    });

    Kmer lookingKmer = kmerBuffer[0];
    size_t write = 0;
    int endFlag = 0;

    for(size_t i = 1 ; i < bufferIdx ; i++) {
        while(taxIdListAtRank[lookingKmer.info.sequenceID] == taxIdListAtRank[kmerBuffer[i].info.sequenceID]){
            if (lookingKmer.ADkmer != kmerBuffer[i].ADkmer) {
                break;
            }
            else{
                lookingKmer.info.redundancy = true;
            }
            i++;
            if(i == bufferIdx){
                endFlag = 1;
                break;
            }
        }

        fwrite(&lookingKmer.info, sizeof(KmerInfo), 1, idAndPosFile);
        write++;
        writeKmerDiff(lastKmer, lookingKmer.ADkmer, diffIdxFile, kmerLocalBuf, localBufIdx);

        if(endFlag == 1) break;

        lastKmer = lookingKmer.ADkmer;
        lookingKmer = kmerBuffer[i];
    }
    if(!((kmerBuffer[bufferIdx - 2].ADkmer == kmerBuffer[bufferIdx - 1].ADkmer) &&
        (kmerBuffer[bufferIdx - 2].info.sequenceID == kmerBuffer[bufferIdx - 1].info.sequenceID))){
        fwrite(&lookingKmer.info, sizeof(KmerInfo), 1, idAndPosFile);
        write++;
        writeKmerDiff(lastKmer, lookingKmer.ADkmer, diffIdxFile, kmerLocalBuf, localBufIdx);
    }

    cout<<"total k-mer count: "<<bufferIdx<<endl;
    cout<<"written k-mer count: "<<write<<endl;
    flushKmerBuf(kmerLocalBuf, diffIdxFile, localBufIdx);

    free(kmerLocalBuf);
    fclose(diffIdxFile);
    fclose(idAndPosFile);
    bufferIdx = 0;
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

void IndexCreator::writeInfo(KmerInfo * entryToWrite, FILE * infoFile, KmerInfo * infoBuffer, size_t & infoBufferIdx)
{
    if (infoBufferIdx >= kmerBufSize) {
        flushInfoBuf(infoBuffer, infoFile, infoBufferIdx);
    }
    memcpy(infoBuffer + infoBufferIdx, entryToWrite, sizeof(KmerInfo));
    infoBufferIdx++;
}
void IndexCreator::flushInfoBuf(KmerInfo * buffer, FILE * infoFile, size_t & localBufIdx ) {
    fwrite(buffer, sizeof(KmerInfo), localBufIdx, infoFile);
    localBufIdx = 0;
}
int IndexCreator::getNumOfFlush()
{
    return numOfFlush;
}