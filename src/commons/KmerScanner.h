#ifndef METABULI_KMERSCANNER_H
#define METABULI_KMERSCANNER_H

#include <iostream>

#include "Kmer.h"
#include "GeneticCode.h"
#include "common.h"

struct Kmer {
    uint64_t value;
    uint32_t pos;
    void printAA(const GeneticCode & code) const {
        uint64_t aaPart = value >> 24;
        for (int i = 0; i < 8; ++i) {
            int aa = (aaPart >> (35 - 5 * i)) & 0x1F;
            std::cout << code.aminoacids[aa];
        }
    }

    void printDNA(const GeneticCode & code) const {
        uint64_t dnaPart = value & 0xFFFFFF;
        uint64_t aaPart = value >> 24;
        for (int i = 0; i < 8; ++i) {
            int aa = (aaPart >> (35 - 5 * i)) & 0x1F;
            int codon = (dnaPart >> (21 - 3 * i)) & 0x7;
            std::cout << code.aa2codon[aa][codon];
        }
    }
};

class KmerScanner {
protected:
    // Internal values
    const GeneticCode &geneticCode;
    uint64_t dnaMask;

    const char *seq;
    size_t seqStart;
    size_t seqEnd;
    size_t seqLen;
    size_t aaLen;

    uint64_t dnaPart;
    uint64_t aaPart;
    int loadedCharCnt;
    int prevPos;
    int posStart;

public:
    KmerScanner(const GeneticCode &geneticCode) : geneticCode(geneticCode) {
        this->dnaMask = (1ULL << 24) - 1;    
    }

    void initScanner(const char * seq, size_t seqStart, size_t seqEnd) {
        this->seq = seq;
        this->seqStart = seqStart;
        this->seqEnd = seqEnd;
        this->seqLen = seqEnd - seqStart + 1;
        this->aaLen = seqLen / 3;
        this->dnaPart = 0;
        this->aaPart = 0;
        this->loadedCharCnt = 0;
        this->prevPos = -8;
        this->posStart = 0;
    }


    Kmer next(bool forward) {
        int aa = 0;
        int codon = 0;
        while (posStart <= aaLen - 8) {
            bool sawN = false;
            loadedCharCnt -= (loadedCharCnt == 8);
            while (loadedCharCnt < 8) {
                int ci;
                if (forward) {
                    ci = seqStart + (posStart + loadedCharCnt) * 3;
                    aa = geneticCode.getAA(atcg[seq[ci]], atcg[seq[ci + 1]], atcg[seq[ci + 2]]);
                    codon = geneticCode.getCodon(atcg[seq[ci]], atcg[seq[ci + 1]], atcg[seq[ci + 2]]);
                } else {
                    ci = seqEnd - (posStart + loadedCharCnt) * 3;
                    aa = geneticCode.getAA(iRCT[atcg[seq[ci]]], iRCT[atcg[seq[ci - 1]]], iRCT[atcg[seq[ci - 2]]]);          
                    codon = geneticCode.getCodon(iRCT[atcg[seq[ci]]], iRCT[atcg[seq[ci - 1]]], iRCT[atcg[seq[ci - 2]]]);
                }
                if (aa < 0) { sawN = true; break; }
                dnaPart = (dnaPart << 3) | (uint64_t)codon;
                aaPart = (aaPart << 5) | (uint64_t)aa;
                loadedCharCnt++;
            }
            if (sawN) {
                posStart += loadedCharCnt + 1;
                prevPos = posStart - 8; 
                dnaPart = aaPart = 0;
                loadedCharCnt = 0;
                continue;
            }
            prevPos = posStart;
            posStart++;
            if (forward) {
                return { (aaPart << 24) | (dnaPart & dnaMask), seqStart + prevPos * 3 };
            } else {
                return { (aaPart << 24) | (dnaPart & dnaMask), seqEnd - (prevPos + 8) * 3 + 1 };
            }
        }
        return { UINT64_MAX, 0 }; // No more kmers found
    }

};

#endif // METABULI_KMERSCANNER_H