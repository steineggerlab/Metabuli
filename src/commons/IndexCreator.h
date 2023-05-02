#ifndef ADKMER4_INDEXCREATOR_H
#define ADKMER4_INDEXCREATOR_H
#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include <ctime>
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

// For masking
#include "NucleotideMatrix.h"
#include "SubstitutionMatrix.h"
#include "tantan.h"
//#include "DBReader.h"
//#include "DBWriter.h"
//#include "Debug.h"
//#include "Util.h"
//#include "FileUtil.h"

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
    uint64_t MARKER;
    string tinfo_path;
    string tinfo_list;
    vector<TaxID> trainedSpecies;
    unordered_map<TaxID, _training> trainingInfo;
    int threadNum;
    BaseMatrix *subMat;

    NcbiTaxonomy * taxonomy;
    string dbDir;
    string fnaListFileName;
    string taxonomyDir;
    string acc2taxidFileName;

    struct FASTA {
        string path;
        TaxID speciesID;
        size_t trainingSeqIdx;
        vector<SequenceBlock> sequences;
    };

    vector<FASTA> fastaList;
    vector<TaxID> taxIdList;
    vector<size_t> processedSeqCnt; // Index of this vector is the same as the index of fnaList


    struct FnaSplit{
        // species, file_idx, training, offset, cnt
        size_t training;
        size_t offset;
        size_t cnt;
        TaxID speciesID;
        int file_idx;
        FnaSplit(size_t training, size_t offset, size_t cnt, TaxID speciesID, int file_idx):
                training(training), offset(offset), cnt(cnt), speciesID(speciesID), file_idx(file_idx) {}
    };
    vector<FnaSplit> fnaSplits;

    struct FastaSplit{
        FastaSplit(size_t training, uint32_t offset, uint32_t cnt, TaxID taxid):
            training(training), offset(offset), cnt(cnt), taxid(taxid) {}
        size_t training;
        uint32_t offset;
        uint32_t cnt;
        TaxID taxid;
    };

   struct Split{
       Split(size_t offset, size_t end) : offset(offset), end(end) {}
       size_t offset;
       size_t end;
   };

    size_t numOfFlush=0;

    void trainProdigal();

//    void writeTargetFiles(TargetKmer * kmerBuffer, size_t & kmerNum, const char * outputFileName,const vector<int> & taxIdList);
    void writeTargetFiles(TargetKmer * kmerBuffer, size_t & kmerNum, const LocalParameters & par, const size_t * uniqeKmerIdx, size_t & uniqKmerCnt);

    void writeTargetFilesAndSplits(TargetKmer * kmerBuffer, size_t & kmerNum, const LocalParameters & par, const size_t * uniqeKmerIdx, size_t & uniqKmerCnt);
    void writeDiffIdx(uint16_t *buffer, FILE* handleKmerTable, uint16_t *toWrite, size_t size, size_t & localBufIdx );
    static bool compareForDiffIdx(const TargetKmer & a, const TargetKmer & b);

    void maskLowComplexityRegions(char * seq, char * maskedSeq, ProbabilityMatrix & probMat,
                                  const LocalParameters & par);
    size_t fillTargetKmerBuffer(TargetKmerBuffer &kmerBuffer,
                                bool *checker,
                                size_t &processedSplitCnt,
                                const LocalParameters &par);


    void makeBlocksForParallelProcessing();

    void splitFastaForProdigalTraining(int file_idx, TaxID speciesID);

    void unzipAndList(const string & folder, const string & fastaList_fname){
        system(("./../../util/unzip_and_list.sh " + folder + " " + fastaList_fname).c_str());
    }

    void load_assacc2taxid(const string & mappingFile, unordered_map<string, int> & assacc2taxid);
    static void load_accession2taxid(const string & mappingFile, unordered_map<string, int> & assacc2taxid);

    void reduceRedundancy(TargetKmerBuffer & kmerBuffer, size_t * uniqeKmerIdx, size_t & uniqKmerCnt,
                          const LocalParameters & par);
    size_t AminoAcidPart(size_t kmer) {
        return (kmer) & MARKER;
    }

    int getNumberOfLines(const string & filename){
        ifstream file(filename);
        int cnt = 0;
        string line;
        while (getline(file, line)) {
            cnt++;
        }
        file.close();
        return cnt;
    }

public:
    static void splitSequenceFile(vector<SequenceBlock> & seqSegments, MmapedData<char> seqFile);

    string getSeqSegmentsWithHead(vector<SequenceBlock> & seqSegments, const string & seqFileName,
                                  const unordered_map<string, TaxID> & acc2taxid,
                                  unordered_map<string, TaxID> & foundAcc2taxid);
    static void getSeqSegmentsWithHead(vector<SequenceBlock> & seqSegments, const char * seqFileName);
    IndexCreator(const LocalParameters & par);
    IndexCreator(const LocalParameters & par, string dbDir, string fnaListFileName, string acc2taxidFile);
    IndexCreator() {taxonomy = nullptr;}
    ~IndexCreator();
    int getNumOfFlush();
    void startIndexCreatingParallel(const LocalParameters & par);

    void createIndex(const LocalParameters & par);

    void updateIndex(const LocalParameters & par);

    void getDiffIdx(const uint64_t & lastKmer, const uint64_t & entryToWrite, FILE* handleKmerTable,
                    uint16_t *kmerBuf, size_t & localBufIdx );
    void getDiffIdx(const uint64_t & lastKmer, const uint64_t & entryToWrite, FILE* handleKmerTable,
                    uint16_t *kmerBuf, size_t & localBufIdx, size_t & totalBufferIdx);
    void writeInfo(TargetKmerInfo * entryToWrite, FILE * infoFile, TargetKmerInfo * infoBuffer, size_t & infoBufferIdx);
    static void flushKmerBuf(uint16_t *buffer, FILE *handleKmerTable, size_t & localBufIdx);
    static void flushInfoBuf(TargetKmerInfo * buffer, FILE * infoFile, size_t & localBufIdx );

};
#endif //ADKMER4_INDEXCREATOR_H
