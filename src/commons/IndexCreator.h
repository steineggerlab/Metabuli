//
// Created by KJB on 01/09/2020.
//

#ifndef ADKMER4_INDEXCREATOR_H
#define ADKMER4_INDEXCREATOR_H
#include <iostream>
#include <string>
#include <vector>
#include "printBinary.h"
#include "Mmap.h"
#include <fstream>
#include "Kmer.h"
#include <time.h>
#include "KmerExtractor.h"
#include "BitManipulateMacros.h"
#include "common.h"
#define kmerLength 8



using namespace std;

class IndexCreator{
private:
    size_t availableMemory;
  //  char * outPath;
    int numOfFlush=0;
    KmerExtractor * kmerExtractor;


    void writeTargetFiles(Kmer * kmerBuffer, size_t & bufferIdx, char * outputFileName, vector<int> & taxIdList);
    void writeDiffIndexFile(Kmer * kmerBuffer, const size_t & bufferIdx, FILE * diffIdxFile);
    void writeKmer(uint16_t *buffer, FILE* handleKmerTable, uint16_t *toWrite, size_t size, size_t & localBufIdx );
//    void addSequence(Kmer *kmerBuffer, size_t & bufferIdx, int & seqID, string & forwardRead);

public:

    IndexCreator();
    ~IndexCreator();
    int getNumOfFlush();
    void startIndexCreating(ifstream & targetFile, char * outputFileName, vector<int> & taxIdList);
    void startIndexCreating2(ifstream & targetFile, char * outputFileName, vector<int> & taxIdList);
    void writeKmerDiff(uint64_t lastKmer, uint64_t & entryToWrite, FILE* handleKmerTable, uint16_t *kmerBuf, size_t & localBufIdx );
    void writeInfo(KmerInfo * entryToWrite, FILE * infoFile, KmerInfo * infoBuffer, size_t & infoBufferIdx);
    void flushKmerBuf(uint16_t *buffer, FILE *handleKmerTable, size_t & localBufIdx);
    void flushInfoBuf(KmerInfo * buffer, FILE * infoFile, size_t & localBufIdx );
};
#endif //ADKMER4_INDEXCREATOR_H
