#ifndef METABULI_KMERSCANNER_H
#define METABULI_KMERSCANNER_H

#include <iostream>
#include <deque>

#include "Kmer.h"
#include "GeneticCode.h"
#include "common.h"

struct Kmer {
    uint64_t value;
    uint32_t pos;
    Kmer() : value(0), pos(0) {}
    Kmer(uint64_t value, uint32_t pos) : value(value), pos(pos) {}
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
    uint32_t seqStart;
    uint32_t seqEnd;
    uint32_t seqLen;
    int aaLen;
    bool isForward;

    uint64_t dnaPart;
    uint64_t aaPart;
    int loadedCharCnt;
    int posStart;

public:
    KmerScanner(const GeneticCode &geneticCode) : geneticCode(geneticCode) {
        // std::cout << "KmerScanner initialized." << std::endl;
        this->dnaMask = (1ULL << 24) - 1;    
    }

    virtual ~KmerScanner() {
        // std::cout << "KmerScanner destroyed." << std::endl;
    }

    const GeneticCode &getGeneticCode() const {
        return geneticCode;
    }

    virtual void initScanner(const char * seq, size_t seqStart, size_t seqEnd, bool isForward) {
        this->seq = seq;
        this->seqStart = seqStart;
        this->seqEnd = seqEnd;
        this->seqLen = seqEnd - seqStart + 1;
        this->aaLen = seqLen / 3;
        this->dnaPart = 0;
        this->aaPart = 0;
        this->loadedCharCnt = 0;
        this->posStart = 0;
        this->isForward = isForward;
    }


    virtual Kmer next() {
        int aa = 0;
        int codon = 0;
        while (posStart <= aaLen - 8) {
            bool sawN = false;
            loadedCharCnt -= (loadedCharCnt == 8);
            while (loadedCharCnt < 8) {
                int ci;
                if (isForward) {
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
                dnaPart = aaPart = 0;
                loadedCharCnt = 0;
                continue;
            }
            if (isForward) {
                return { (aaPart << 24) | (dnaPart & dnaMask), seqStart + (posStart++) * 3 };
            } else {
                return { (aaPart << 24) | (dnaPart & dnaMask), seqEnd - ((posStart++) + 8) * 3 + 1 };
            }
        }
        return { UINT64_MAX, 0 }; // No more kmers found
    }
};

// OldKmerScanner is made to support searching old-format databases
// Implementation is very efficient and puzzling, but it works
class OldKmerScanner : public KmerScanner {
    private: 
        std::deque<size_t> dq;
    public:
        OldKmerScanner(const GeneticCode &geneticCode) : KmerScanner(geneticCode) {}

        ~OldKmerScanner() {
            // std::cout << "OldKmerScanner destroyed." << std::endl;
        }
        
        void initScanner(const char * seq, size_t seqStart, size_t seqEnd, bool isForward) override {
            KmerScanner::initScanner(seq, seqStart, seqEnd, isForward);
            dq.clear();
        }

        Kmer next() override {
            int aa = 0;
            int codon = 0;
            while (posStart <= aaLen - 8) {
                bool sawN = false;
                loadedCharCnt -= (loadedCharCnt == 8);
                while (loadedCharCnt < 8) {
                    int ci;
                    if (isForward) {
                        ci = seqEnd - (posStart + loadedCharCnt) * 3;
                        aa    = geneticCode.getAA   (atcg[seq[ci - 2]], atcg[seq[ci - 1]], atcg[seq[ci]]);
                        codon = geneticCode.getCodon(atcg[seq[ci - 2]], atcg[seq[ci - 1]], atcg[seq[ci]]);
                    } else {
                        ci = seqStart + (posStart + loadedCharCnt) * 3;
                        aa    = geneticCode.getAA   (iRCT[atcg[seq[ci + 2]]], iRCT[atcg[seq[ci + 1]]], iRCT[atcg[seq[ci]]]);
                        codon = geneticCode.getCodon(iRCT[atcg[seq[ci + 2]]], iRCT[atcg[seq[ci + 1]]], iRCT[atcg[seq[ci]]]);
                    }
                    if (aa < 0) { sawN = true; break; }
                    if (dq.size() == 8) {
                        aaPart = aaPart - dq.back();
                        dq.pop_back();
                    }
                    for (auto &x : dq) {
                        x *= 21;
                    }
                    dq.emplace_front(aa);
                    aaPart = aaPart * 21 + aa;
                    dnaPart = (dnaPart << 3) | (uint64_t)codon;
                    loadedCharCnt++;
                }
                if (sawN) {
                    posStart += loadedCharCnt + 1;
                    dnaPart = aaPart = 0;
                    loadedCharCnt = 0;
                    dq.clear();
                    continue;
                }
                if (isForward) {
                    return { (aaPart << 24) | (dnaPart & dnaMask), seqEnd - ((posStart++) + 8) * 3 + 1};
                } else {
                    return { (aaPart << 24) | (dnaPart & dnaMask), seqStart + (posStart++) * 3 };
                }    
            }
            return { UINT64_MAX, 0 }; // No more kmers found
        }    
};

#endif // METABULI_KMERSCANNER_H