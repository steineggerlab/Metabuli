//
// Created by KJB on 01/09/2020.
//

#include "IndexCreator.h"
//#include "NcbiTaxonomy.h"
IndexCreator::IndexCreator()
{
    kmerExtractor = new KmerExtractor();
}

void IndexCreator::takeThisSequenceFile(char * referenceFileName, char * outputFileName)
{
    string buffer;
    string forwardRead;
    string reverseComplimentRead;
    string reads[2];
    ExtractStartPoint ESP = {0 ,0};
    size_t bufferIdx = 0;
    ifstream targetFile(referenceFileName);


    Kmer * kmerBuffer = (Kmer *)malloc(sizeof(Kmer) * kmerBufSize);
    //KmerInfo * infoBuffer = (KmerInfo *)malloc(sizeof(KmerInfo) * kmerBufSize);

    getline(targetFile, buffer);
    int seqID = 1;

    /// mmap을 써보자
    while(targetFile)
    {
        getline(targetFile, buffer);
        if(buffer[0] == '>'){
            reverseComplimentRead = kmerExtractor->reverseCompliment(forwardRead);
            reads[0] = forwardRead; reads[1] = reverseComplimentRead;
            kmerExtractor->dna2aa(forwardRead, reverseComplimentRead);

            ESP = kmerExtractor->fillKmerBuffer(reads, kmerBuffer, seqID, bufferIdx, ESP);
            while (ESP.startOfFrame + ESP.frame != 0)
            {
                writeTargetFiles(kmerBuffer, bufferIdx, outputFileName);
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
        writeTargetFiles(kmerBuffer, bufferIdx, outputFileName);
        ESP = kmerExtractor->fillKmerBuffer(reads, kmerBuffer, seqID, bufferIdx, ESP);
    }
    //
//    for(int i = 0; i <bufferIdx; i++)
//    {
//        cout<<kmerBuffer[i].ADkmer<<endl;
//    }

    //flush last buffer
    writeTargetFiles(kmerBuffer, bufferIdx, outputFileName);
    targetFile.close();
    free(kmerBuffer);
}

void IndexCreator::writeTargetFiles(Kmer *kmerBuffer, size_t & bufferIdx, char * outputFileName)
{
    char suffixedDiffIdxFileName[100];
    char suffixedInfoFileName[100];

    sprintf(suffixedDiffIdxFileName,"%s_diffIdx_%zu.txt", outputFileName,numOfFlush);
    sprintf(suffixedInfoFileName,"%s_info_%zu.txt", outputFileName,numOfFlush);


    FILE * diffIdxFile = fopen(suffixedDiffIdxFileName, "wb");
    FILE * idAndPosFile = fopen(suffixedInfoFileName, "wb");
    sort(kmerBuffer, kmerBuffer + bufferIdx, [=](Kmer x, Kmer y) { return x.ADkmer < y.ADkmer; });

    //// make it faster
    for(size_t i = 0 ; i < bufferIdx ; i++)
    {
        fwrite(&kmerBuffer[i].info, sizeof(KmerInfo), 1, idAndPosFile);
    }

    writeDiffIndexFile(kmerBuffer, bufferIdx, diffIdxFile);


    fclose(idAndPosFile);
    bufferIdx = 0;
    numOfFlush++;
}

void IndexCreator::writeDiffIndexFile(Kmer * kmerBuffer, const size_t & bufferIdx, FILE * diffIndexFile)
{

    size_t distinctKmerCount = 0;
    size_t numOfKmer = 0;
    uint64_t lastKmer = 0;

    uint16_t *kmerLocalBuf = (uint16_t *)malloc(sizeof(uint16_t) * kmerBufSize);
    size_t localBufIdx = 0;

    uint64_t marker= ~0 & ~16777215;
    for(size_t i = 0; i < bufferIdx ; i++)
    {
        //if(posInTable->ADkmer != entryToWrite->ADkmer){
        writeKmerDiff(lastKmer, kmerBuffer[i].ADkmer, diffIndexFile, kmerLocalBuf, localBufIdx);
        if((lastKmer & marker) != (kmerBuffer[i].ADkmer & marker))
       // if(lastKmer != kmerBuffer[i].ADkmer)
        {
            ++distinctKmerCount;
        }
        lastKmer = kmerBuffer[i].ADkmer;
        //entryToWrite = posInTable; posInTable++

        // }
        numOfKmer++;
    }
    // writeKmerDiff(lastKmer, entryToWrite, diffIndexFile, handleIDTable, kmerLocalBuf, IDLocalBuf);
    cout<<"Total target kmer count    : "<<numOfKmer<<endl;
    cout<<"Distinct target kmer count : "<<distinctKmerCount<<endl;

    flushKmerBuf(kmerLocalBuf, diffIndexFile, localBufIdx);
    free(kmerLocalBuf);
    fclose(diffIndexFile);
}

void IndexCreator::writeKmerDiff(uint64_t lastKmer, uint64_t & entryToWrite, FILE* handleKmerTable, uint16_t *kmerBuf, size_t & localBufIdx ){
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