#ifndef ADKMER4_KMEREXTRACTOR_H
#define ADKMER4_KMEREXTRACTOR_H

#include <iostream>
#include <vector>
#include <queue>
#include <cstdint>
#include <algorithm>
#include <functional>

#include "common.h"
#include "Mmap.h"
#include "xxhash.h"
#include "LocalParameters.h"
#include "NucleotideMatrix.h"
#include "SubstitutionMatrix.h"
#include "tantan.h"
#include "LocalUtil.h"
#include "SyncmerScanner.h"
#include "Kmer.h"
#include "printBinary.h"
#include "kseq.h"
#include "KSeqBufferReader.h"
#include "ProdigalWrapper.h"
#ifdef OPENMP
    #include <omp.h>
#endif

KSEQ_INIT(kseq_buffer_t*, kseq_buffer_reader)

#define nuc2int(x) (x & 14u)>>1u

using namespace std;

class SeqIterator {
private:
    uint64_t * powers;
    uint32_t * mask;
    int * mask_int;
    uint32_t spaceNum;
    int spaceNum_int;
    int bitsForCodon;
    int bitsFor8Codons;
    uint32_t smerMask;
    uint64_t dnaMask;
    int kmerLen = 12;
    int smerLen;
    
public:
    SeqIterator(const LocalParameters &par); 
    
    ~SeqIterator();
    
    string reverseComplement(string &read) const;
    
    void devideToCdsAndNonCds(const char *maskedSeq,
                              size_t seqLen,
                              const vector<CDSinfo> &cdsInfo, 
                              vector<string> &cds,
                              vector<string> &nonCds);

    char *reverseComplement(char *read, size_t length) const;

    void generateIntergenicKmerList(struct _gene *genes, struct _node *nodes, int numberOfGenes,
                                    vector<uint64_t> &intergenicKmerList, const char *seq);

    void getMinHashList(priority_queue<uint64_t> &sortedHashQue, const char *seq);

    bool compareMinHashList(priority_queue<uint64_t> list1, priority_queue<uint64_t> &list2, size_t length1, size_t length2);

    bool isSyncmer(const vector<int> &aaSeq, int startPos, int k, int s) {
        size_t min_smer_value = UINT64_MAX;
        int min_smer_pos = -1;
        size_t current_value = 0;
        for (int i = 0; i <= k - s; ++i) {
            if (i == 0) {
                for (int j = 0; j < s; ++j) {
                    current_value = (current_value << 5) | aaSeq[startPos + i + j];
                }
            } else {
                current_value = (current_value << 5) | aaSeq[startPos + i + s - 1];
            }
            current_value = current_value & ((1ULL << (5 * s)) - 1); // Mask to keep only the last s amino acids
            if (current_value < min_smer_value) {
                min_smer_value = current_value;
                min_smer_pos = i;
            }
        }
        return (min_smer_pos == 0 || min_smer_pos == (k - s));
    }

    static void maskLowComplexityRegions(const unsigned char * seq, unsigned char * maskedSeq, ProbabilityMatrix & probMat,
                                         float maskProb, const BaseMatrix * subMat);

    void printKmerInDNAsequence(uint64_t kmer);

    void printAAKmer(uint64_t kmer, int shits = 28);

    void printSyncmer(const uint64_t syncmer) {
        string aminoacid = "ARNDCQEGHILKMFPSTWYVX";

        uint64_t dnaPart = syncmer;
        uint64_t aaPart = syncmer >> bitsFor8Codons;
        string aaStr = "01234567";
        string dnaStr;
        vector<string> dna24mer(8);
        for (int i = 0; i < 8; i++) {
            aaStr[7 - i] = aminoacid[(aaPart >> (5 * i)) & 0x1F];
            uint64_t dnaInfo = (dnaPart >> (3 * i)) & 0x07;
            switch ((aaPart >> (5 * i)) & 0x1F) {
                case 0: //A
//                cout << "A";
                    if (dnaInfo == 0) {
                        dna24mer[7 - i] = "GCA";
                    } else if (dnaInfo == 1) {
                        dna24mer[7 - i] = "GCC";
                    } else if (dnaInfo == 2) {
                        dna24mer[7 - i] = "GCT";
                    } else if (dnaInfo == 3){
                        dna24mer[7 - i] = "GCG";
                    } else {
                        cout << "Error in " << aaStr[7 - i] << endl;
                    }
                    break;
                case 1: //R
//                cout << "R";
                    if (dnaInfo == 0) {
                        dna24mer[7 - i] = "CGA";
                    } else if (dnaInfo == 1) {
                        dna24mer[7 - i] = "CGC";
                    } else if (dnaInfo == 2) {
                        dna24mer[7 - i] = "CGT";
                    } else if (dnaInfo == 3) {
                        dna24mer[7 - i] = "CGG";
                    } else if (dnaInfo == 4) {
                        dna24mer[7 - i] = "AGG";
                    } else if (dnaInfo == 5) {
                        dna24mer[7 - i] = "AGA";
                    } else{
                        cout << "Error in " << aaStr[7 - i] << endl;
                    }
                    break;
                case 2: //N
//                cout << "N";
                    if (dnaInfo == 0) {
                        cout << "dddddd";
                    } else if (dnaInfo == 1) {
                        dna24mer[7 - i] = "AAC";
                    } else if (dnaInfo == 2) {
                        dna24mer[7 - i] = "AAT";
                    } else {
                        cout << "dddddd";
                    }
                    break;
                case 3: //D
//                cout << "D";
                    if (dnaInfo == 0) {
                        cout << "dddddd";
                    } else if (dnaInfo == 1) {
                        dna24mer[7 - i] = "GAC";
                    } else if (dnaInfo == 2) {
                        dna24mer[7 - i] = "GAT";
                    } else {
                        cout << "dddddd";
                    }
                    break;
                case 4: //C
//                cout << "C";
                    if (dnaInfo == 0) {
                        cout << "dddddd";
                    } else if (dnaInfo == 1) {
                        dna24mer[7 - i] = "TGC";
                    } else if (dnaInfo == 2) {
                        dna24mer[7 - i] = "TGT";
                    } else {
                        cout << "dddddd";
                    }
                    break;
                case 5: // Q
//                cout << "Q";
                    if (dnaInfo == 0) {
                        dna24mer[7 - i] = "CAA";
                    } else if (dnaInfo == 1) {
                        cout << "dddddd";
                    } else if (dnaInfo == 2) {
                        cout << "dddddd";
                    } else {
                        dna24mer[7 - i] = "CAG";
                    }
                    break;
                case 6: //E
//                cout << "E";
                    if (dnaInfo == 0) {
                        dna24mer[7 - i] = "GAA";
                    } else if (dnaInfo == 1) {
                        cout << "dddddd";
                    } else if (dnaInfo == 2) {
                        cout << "dddddd" ;
                    } else {
                        dna24mer[7 - i] = "GAG";
                    }
                    break;
                case 7: //G
//                cout << "G";
                    if (dnaInfo == 0) {
                        dna24mer[7 - i] = "GGA";
                    } else if (dnaInfo == 1) {
                        dna24mer[7 - i] = "GGC";
                    } else if (dnaInfo == 2) {
                        dna24mer[7 - i] = "GGT";
                    } else {
                        dna24mer[7 - i] = "GGG";
                    }
                    break;
                case 8: //H
//                cout << "H";
                    if (dnaInfo == 0) {
                        cout << "dddddd";
                    } else if (dnaInfo == 1) {
                        dna24mer[7 - i] = "CAC";
                    } else if (dnaInfo == 2) {
                        dna24mer[7 - i] = "CAT";
                    } else {
                        cout << "dddddd";
                    }
                    break;
                case 9: //I
//                cout << "I";
                    if (dnaInfo == 0) {
                        dna24mer[7 - i] = "ATA";
                    } else if (dnaInfo == 1) {
                        dna24mer[7 - i] = "ATC";
                    } else if (dnaInfo == 2) {
                        dna24mer[7 - i] = "ATT";
                    } else {
                        cout << "dddddd";
                    }
                    break;
                case 10: //L
//                cout << "L";
                    if (dnaInfo == 0) {
                        dna24mer[7 - i] = "CTA";
                    } else if (dnaInfo == 1) {
                        dna24mer[7 - i] = "CTC";
                    } else if (dnaInfo == 2) {
                        dna24mer[7 - i] = "CTT";
                    } else if (dnaInfo == 3) {
                        dna24mer[7 - i] = "CTG";
                    } else if (dnaInfo == 4) {
                        dna24mer[7 - i] = "TTG";
                    } else {
                        dna24mer[7 - i] = "TTA";
                    }
                    break;
                case 11: //K
//                cout << "K";
                    if (dnaInfo == 0) {
                        dna24mer[7 - i] = "AAA";
                    } else if (dnaInfo == 1) {
                        cout << "dddddd";
                    } else if (dnaInfo == 2) {
                        cout << "dddddd";
                    } else {
                        dna24mer[7 - i] = "AAG";
                    }
                    break;
                case 12: // M
//                cout << "M";
                    if (dnaInfo == 0) {
                        cout << "dddddd";
                    } else if (dnaInfo == 1) {
                        cout << "dddddd";
                    } else if (dnaInfo == 2) {
                        cout << "dddddd";
                    } else {
                        dna24mer[7 - i] = "ATG";
                    }
                    break;
                case 13://F
//                cout << "F";
                    if (dnaInfo == 0) {
                        cout << "dddddd" ;
                    } else if (dnaInfo == 1) {
                        dna24mer[7 - i] = "TTC";
                    } else if (dnaInfo == 2) {
                        dna24mer[7 - i] = "TTT";
                    } else {
                        cout << "dddddd";
                    }
                    break;
                case 14: //P
//                cout << "P";
                    if (dnaInfo == 0) {
                        dna24mer[7 - i] = "CCA";
                    } else if (dnaInfo == 1) {
                        dna24mer[7 - i] = "CCC";
                    } else if (dnaInfo == 2) {
                        dna24mer[7 - i] = "CCT";
                    } else {
                        dna24mer[7 - i] = "CCG";
                    }
                    break;
                case 15: //S
//                cout << "S";
                    if (dnaInfo == 0) {
                        dna24mer[7 - i] = "TCA";
                    } else if (dnaInfo == 1) {
                        dna24mer[7 - i] = "TCC";
                    } else if (dnaInfo == 2) {
                        dna24mer[7 - i] = "TCT";
                    } else if (dnaInfo == 3) {
                        dna24mer[7 - i] = "TCG";
                    } else if (dnaInfo == 4) {
                        cout << "dddddd";
                    } else if (dnaInfo == 5) {
                        cout << "dddddd";
                    } else if (dnaInfo == 6) {
                        dna24mer[7 - i] = "AGT";
                    } else {
                        dna24mer[7 - i] = "AGC";
                    }
                    break;
                case 16: //T
//                cout << "T";
                    if (dnaInfo == 0) {
                        dna24mer[7 - i] = "ACA";
                    } else if (dnaInfo == 1) {
                        dna24mer[7 - i] = "ACC";
                    } else if (dnaInfo == 2) {
                        dna24mer[7 - i] = "ACT";
                    } else {
                        dna24mer[7 - i] = "ACG";
                    }
                    break;
                case 17: //W
//                cout << "W";
                    if (dnaInfo == 0) {
                        cout << "dddddd";
                    } else if (dnaInfo == 1) {
                        cout << "dddddd";
                    } else if (dnaInfo == 2) {
                        cout << "dddddd";
                    } else {
                        dna24mer[7 - i] = "TGG";
                    }
                    break;
                case 18: //Y
//                cout << "Y";
                    if (dnaInfo == 0) {
                        cout << "dddddd";
                    } else if (dnaInfo == 1) {
                        dna24mer[7 - i] = "TAC";
                    } else if (dnaInfo == 2) {
                        dna24mer[7 - i] = "TAT";
                    } else {
                        cout << "dddddd";
                    }
                    break;
                case 19: //V
//                cout << "V";
                    if (dnaInfo == 0) {
                        dna24mer[7 - i] = "GTA";
                    } else if (dnaInfo == 1) {
                        dna24mer[7 - i] = "GTC";
                    } else if (dnaInfo == 2) {
                        dna24mer[7 - i] = "GTT";
                    } else {
                        dna24mer[7 - i] = "GTG";
                    }
                    break;
                case 20: //stop
//                cout << "Z";
                    if (dnaInfo == 0) {
                        dna24mer[7 - i] = "TAA";
                    } else if (dnaInfo == 1) {
                        cout << "dddddd";
                    } else if (dnaInfo == 2) {
                        cout << "dddddd";
                    } else if (dnaInfo == 3) {
                        dna24mer[7 - i] = "TAG";
                    } else if (dnaInfo == 4) {
                        cout << "dddddd";
                    } else {
                        dna24mer[7 - i] = "TGA";
                    }
                    break;
            }
        }
        dnaStr = "";
        for (int i = 0; i < 8; i++) {
            dnaStr += dna24mer[i];
        }
        cout << "Syncmer: " << aaStr << " " << dnaStr;
    }


};

#endif //ADKMER4_KMEREXTRACTOR_H

