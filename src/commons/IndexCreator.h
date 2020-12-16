//
// Created by KJB on 01/09/2020.
//

#ifndef ADKMER4_INDEXCREATOR_H
#define ADKMER4_INDEXCREATOR_H
#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include <time.h>
#include <fstream>
#include <omp.h>
#include "printBinary.h"
#include "Mmap.h"
#include "Kmer.h"
#include "SeqIterator.h"
#include "BitManipulateMacros.h"
#include "common.h"
#include "NcbiTaxonomy.h"

#include "FastSort.h"


#define kmerLength 8



using namespace std;

class IndexCreator{
private:
    size_t availableMemory;
  //  char * outPath;
    size_t numOfFlush=0;
    SeqIterator * seqIterator;

    void writeTargetFiles(TargetKmer * kmerBuffer, size_t & bufferIdx, const char * outputFileName, vector<int> & taxIdList);
    void writeKmer(uint16_t *buffer, FILE* handleKmerTable, uint16_t *toWrite, size_t size, size_t & localBufIdx );


public:

    IndexCreator();
    ~IndexCreator();
    int getNumOfFlush();
    void startIndexCreating(const char * seqFileName, const char * outputFileName, vector<int> & taxIdList);
    void startIndexCreatingParallel(const char * seqFileName, const char * outputFileName, vector<int> & taxIdList);
    void writeKmerDiff(const uint64_t & lastKmer, const uint64_t & entryToWrite, FILE* handleKmerTable, uint16_t *kmerBuf, size_t & localBufIdx );
    void writeInfo(KmerInfo * entryToWrite, FILE * infoFile, KmerInfo * infoBuffer, size_t & infoBufferIdx);
    void flushKmerBuf(uint16_t *buffer, FILE *handleKmerTable, size_t & localBufIdx);
    void flushInfoBuf(KmerInfo * buffer, FILE * infoFile, size_t & localBufIdx );
};
#endif //ADKMER4_INDEXCREATOR_H
