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
#include "LocalParameters.h"
#include <omp.h>
#ifdef OPENMP
#include <omp.h>
#endif


#define kmerLength 8

struct TaxId2Fasta{
    TaxID species;
    TaxID taxid;
    string fasta;
    TaxId2Fasta(TaxID sp, TaxID ssp, string fasta): species(sp), taxid(ssp), fasta(std::move(fasta)) {}
};

using namespace std;

class IndexCreator{
private:
   struct FastaSplit{
        FastaSplit(size_t training, uint32_t offset, uint32_t cnt): training(training), offset(offset), cnt(cnt) {}
        size_t training;
        uint32_t offset;
        uint32_t cnt;
    };

   struct Split{
       Split(size_t offset, size_t end) : offset(offset), end(end) {}
       size_t offset;
       size_t end;
   };

    size_t availableMemory;
    size_t numOfFlush=0;
    //SeqIterator * seqIterator;
    void writeTargetFiles(TargetKmer * kmerBuffer, size_t & kmerNum, const char * outputFileName,const vector<TaxId2Fasta> & taxid2fasta, const LocalParameters & par);
    void writeTargetFiles(TargetKmer * kmerBuffer, size_t & kmerNum, const char * outputFileName,const vector<int> & taxIdList);
    void writeTargetFiles(TargetKmer * kmerBuffer, size_t & kmerNum, const char * outputFileName, const vector<TaxId2Fasta> & taxid2fasta, size_t * uniqeKmerIdx, size_t & uniqKmerCnt);
    void writeDiffIdx(uint16_t *buffer, FILE* handleKmerTable, uint16_t *toWrite, size_t size, size_t & localBufIdx );
    static bool compareForDiffIdx(const TargetKmer & a, const TargetKmer & b);
    static size_t fillTargetKmerBuffer2(TargetKmerBuffer & kmerBuffer,
                                        bool * checker,
                                        size_t & processedTaxIdCnt,
                                        const vector<FastaSplit> & splits,
                                        const vector<TaxId2Fasta> & taxid2fasta,
                                        const LocalParameters & par);

    static size_t fillTargetKmerBuffer2(TargetKmerBuffer & kmerBuffer,
                                        MmapedData<char> & seqFile,
                                        vector<Sequence> & seqs,
                                        bool * checker,
                                        size_t & processedTaxIdCnt,
                                        const vector<FastaSplit> & splits,
                                        const vector<int> & taxIdList,
                                        const LocalParameters & par);

    static size_t fillTargetKmerBuffer3(TargetKmerBuffer & kmerBuffer,
                                        bool * checker,
                                        size_t & processedTaxIdCnt,
                                        const vector<FastaSplit> & splits,
                                        const vector<TaxId2Fasta> & taxid2fasta,
                                        const LocalParameters & par);

    static void extractKmerFromFasta(SeqIterator &seqIterator, MmapedData<char> &seqFile, priority_queue<uint64_t> &standardList,
                     size_t lengthOfTrainingSeq, const vector<Sequence> &sequences, ProdigalWrapper &prodigal,
                     vector<uint64_t> &intergenicKmerList, TargetKmerBuffer &kmerBuffer, size_t posToWrite,
                     uint32_t seqID, int taxIdAtRank, size_t startIdx);



    static void getSeqSegmentsWithoutHead(vector<Sequence> & seqSegments, MmapedData<char> seqFile);
    void getFastaSplits(const vector<int> & taxIdListAtRank, vector<FastaSplit> & fastaSplit, vector<Sequence> & seqSegments);
    void getFastaSplits2(const vector<TaxId2Fasta> & taxIdListAtRank, vector<FastaSplit> & fastaSplit);
    void mappingFromTaxIDtoFasta(const string & fastaList_fname,
                                 unordered_map<string, int> & assacc2taxid,
                                 vector<TaxId2Fasta> & taxid2fasta,
                                 NcbiTaxonomy & taxonomy);

    void unzipAndList(const string & folder, const string & fastaList_fname){
        system(("./../../util/unzip_and_list.sh " + folder + " " + fastaList_fname).c_str());
    }

    void load_assacc2taxid(const string & mappingFile, unordered_map<string, int> & assacc2taxid);

    static size_t estimateKmerNum(const vector<TaxId2Fasta> & taxid2fasta, const FastaSplit & split);
    void reduceRedundancy(TargetKmerBuffer & kmerBuffer, size_t * uniqeKmerIdx, size_t & uniqKmerCnt, const LocalParameters & par,
                          const vector<TaxId2Fasta> & taxid2fasta);
public:
    static void getSeqSegmentsWithHead(vector<Sequence> & seqSegments, MmapedData<char> seqFile);
    static void getSeqSegmentsWithHead2(vector<Sequence> & seqSegments, const char * seqFileName);
    IndexCreator();
    ~IndexCreator();
    int getNumOfFlush();
    void startIndexCreatingParallel(const LocalParameters & par);
    void startIndexCreatingParallel(const char * seqFileName, const char * outputFileName,
                                    const vector<int> & taxIdListAtRank, const vector<int> & taxIdList,
                                    const LocalParameters & par);
    void startIndexCreatingParallel2(const char * seqFileName, const char * outputFileName, const vector<int> & taxIdListAtSpecies, const vector<int> & taxIdListAtGenus, const vector<int> & taxIdList);

    void getDiffIdx(const uint64_t & lastKmer, const uint64_t & entryToWrite, FILE* handleKmerTable, uint16_t *kmerBuf, size_t & localBufIdx );
    void writeInfo(TargetKmerInfo * entryToWrite, FILE * infoFile, TargetKmerInfo * infoBuffer, size_t & infoBufferIdx);
    void flushKmerBuf(uint16_t *buffer, FILE *handleKmerTable, size_t & localBufIdx);
    void flushInfoBuf(TargetKmerInfo * buffer, FILE * infoFile, size_t & localBufIdx );

};
#endif //ADKMER4_INDEXCREATOR_H
