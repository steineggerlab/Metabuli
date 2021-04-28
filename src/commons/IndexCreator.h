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
#include "printBinary.h"
#include "Mmap.h"
#include "Kmer.h"
#include "SeqIterator.h"
#include "BitManipulateMacros.h"
#include "common.h"
#include "NcbiTaxonomy.h"
#include "FastSort.h"
#include "Classifier.h"

#ifdef OPENMP
#include <omp.h>
#endif


#define kmerLength 8



using namespace std;

class IndexCreator{
private:
   struct FastaSplit{
        FastaSplit(size_t training, uint32_t offset, uint32_t cnt): training(training), offset(offset), cnt(cnt) {}
        size_t training;
        uint32_t offset;
        uint32_t cnt;
    };

    size_t availableMemory;
    size_t numOfFlush=0;
    SeqIterator * seqIterator;
    void writeTargetFiles(TargetKmer * kmerBuffer, size_t & kmerNum, const char * outputFileName,const vector<int> & taxIdList);
    void writeTargetFiles2(TargetKmer * kmerBuffer, size_t & kmerNum, const char * outputFileName,const vector<int> & taxIdList, const vector<TaxID> & taxIdListAtGenus);
    void writeDiffIdx(uint16_t *buffer, FILE* handleKmerTable, uint16_t *toWrite, size_t size, size_t & localBufIdx );
    static bool compareForDiffIdx(const TargetKmer & a, const TargetKmer & b);
    size_t fillTargetKmerBuffer(TargetKmerBuffer & kmerBuffer, MmapedData<char> & seqFile, vector<Sequence> & seqs, bool * checker, size_t & processedTaxIdCnt, const vector<FastaSplit> & splits, const vector<int> & taxIdList);
    static void getSeqSegmentsWithoutHead(vector<Sequence> & seqSegments, MmapedData<char> seqFile);
    void getFastaSplits(const vector<int> & taxIdListAtRank, vector<FastaSplit> & fastaSplit, vector<Sequence> & seqSegments);
public:
    static void getSeqSegmentsWithHead(vector<Sequence> & seqSegments, MmapedData<char> seqFile);
    IndexCreator();
    ~IndexCreator();
    int getNumOfFlush();
    void startIndexCreatingParallel(const char * seqFileName, const char * outputFileName, const vector<int> & taxIdListAtRank, const vector<int> & taxIdList);
    void startIndexCreatingParallel2(const char * seqFileName, const char * outputFileName, const vector<int> & taxIdListAtSpecies, const vector<int> & taxIdListAtGenus, const vector<int> & taxIdList);

    void getDiffIdx(const uint64_t & lastKmer, const uint64_t & entryToWrite, FILE* handleKmerTable, uint16_t *kmerBuf, size_t & localBufIdx );
    void writeInfo(TargetKmerInfo * entryToWrite, FILE * infoFile, TargetKmerInfo * infoBuffer, size_t & infoBufferIdx);
    void flushKmerBuf(uint16_t *buffer, FILE *handleKmerTable, size_t & localBufIdx);
    void flushInfoBuf(TargetKmerInfo * buffer, FILE * infoFile, size_t & localBufIdx );

};
#endif //ADKMER4_INDEXCREATOR_H
