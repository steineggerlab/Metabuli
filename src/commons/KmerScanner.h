#ifndef METABULI_KMERSCANNER_H
#define METABULI_KMERSCANNER_H

#include <iostream>
#include <deque>

#include "Kmer.h"
#include "common.h"


class KmerScanner {
protected:
    const char * seq;
    uint32_t seqStart;
    uint32_t seqEnd;
    uint32_t seqLen;
    int loadedCharCnt;
    int posStart;
    int kmerSize;
    bool isForward;

public:
    KmerScanner(int kmerSize) 
        : seq(nullptr), seqStart(0), seqEnd(0), seqLen(0), loadedCharCnt(0), posStart(0), kmerSize(kmerSize) {
    }

    KmerScanner(const char * seq, uint32_t seqStart, uint32_t seqEnd, int kmerSize) 
        : seq(seq), seqStart(seqStart), seqEnd(seqEnd), seqLen(seqEnd - seqStart + 1), loadedCharCnt(0), posStart(0), 
          kmerSize(kmerSize) {
    } 

    virtual ~KmerScanner() {
        // std::cout << "KmerScanner destroyed." << std::endl;
    }

    virtual void initScanner(const char * seq, size_t seqStart, size_t seqEnd, bool isForward = true) {
        this->seq = seq;
        this->seqStart = seqStart;
        this->seqEnd = seqEnd;
        this->seqLen = seqEnd - seqStart + 1;
        this->loadedCharCnt = 0;
        this->posStart = 0;
        this->isForward = isForward;
    }

    virtual Kmer next() = 0; 
};

class MetamerScanner : public KmerScanner {
protected:
    // Internal values
    const GeneticCode &geneticCode;
    uint64_t dnaMask;

    int aaLen;
    uint64_t dnaPart;
    uint64_t aaPart;

public:
    MetamerScanner(const GeneticCode &geneticCode) 
        : KmerScanner(8), geneticCode(geneticCode) {
        // std::cout << "KmerScanner initialized." << std::endl;
        this->dnaMask = (1ULL << 24) - 1;    
    }

    virtual ~MetamerScanner() {
        // std::cout << "KmerScanner destroyed." << std::endl;
    }

    const GeneticCode &getGeneticCode() const {
        return geneticCode;
    }

    void initScanner(const char * seq, size_t seqStart, size_t seqEnd, bool isForward = true) override {
        KmerScanner::initScanner(seq, seqStart, seqEnd, isForward);
        this->aaLen = seqLen / 3;
        this->dnaPart = 0;
        this->aaPart = 0;
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
class OldMetamerScanner : public MetamerScanner {
    private: 
        std::deque<size_t> dq;
    public:
        OldMetamerScanner(const GeneticCode &geneticCode) : MetamerScanner(geneticCode) {}

        ~OldMetamerScanner() {
            // std::cout << "OldKmerScanner destroyed." << std::endl;
        }
        
        void initScanner(const char * seq, size_t seqStart, size_t seqEnd, bool isForward = true) override {
            MetamerScanner::initScanner(seq, seqStart, seqEnd, isForward);
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


class KmerScanner_dna2aa : public KmerScanner {
private:
    const GeneticCode &geneticCode;
    uint64_t mask;
    
    int aaLen;
    uint64_t aaPart;

public:
    KmerScanner_dna2aa(const GeneticCode &geneticCode, int kmerSize) 
        : KmerScanner(kmerSize), geneticCode(geneticCode) 
    {
        this->mask = (1ULL << (5 * kmerSize)) - 1;
        if (kmerSize > 12 || kmerSize < 1) {
            std::cerr << "Error: k must be between 1 and 12, inclusive." << std::endl;
            exit(EXIT_FAILURE);
        } 
    }

    ~KmerScanner_dna2aa() {
        // std::cout << "KmerScanner destroyed." << std::endl;
    }

    const GeneticCode &getGeneticCode() const {
        return geneticCode;
    }

    void initScanner(
        const char * seq, 
        size_t seqStart, 
        size_t seqEnd, 
        bool isForward) override 
    {
        KmerScanner::initScanner(seq, seqStart, seqEnd, isForward);
        this->aaLen = seqLen / 3;
        this->aaPart = 0;
    }


    Kmer next() override {
        int aa = 0;
        while (posStart <= aaLen - kmerSize) {
            bool sawN = false;
            loadedCharCnt -= (loadedCharCnt == kmerSize);
            while (loadedCharCnt < kmerSize) {
                int ci;
                if (isForward) {
                    ci = seqStart + (posStart + loadedCharCnt) * 3;
                    aa = geneticCode.getAA(atcg[seq[ci]], atcg[seq[ci + 1]], atcg[seq[ci + 2]]);
                } else {
                    ci = seqEnd - (posStart + loadedCharCnt) * 3;
                    aa = geneticCode.getAA(iRCT[atcg[seq[ci]]], iRCT[atcg[seq[ci - 1]]], iRCT[atcg[seq[ci - 2]]]);          
                }
                if (aa < 0) { sawN = true; break; }
                aaPart = (aaPart << 5) | (uint64_t)aa;
                loadedCharCnt++;
            }
            if (sawN) {
                posStart += loadedCharCnt + 1;
                aaPart = 0;
                loadedCharCnt = 0;
                continue;
            }

            // Kmer result = { aaPart & mask, 0 };
            // result.printAA(geneticCode, k);
            // std::cout << "\n";

            if (isForward) {
                return { aaPart & mask, seqStart + (posStart++) * 3 };
            } else {
                return { aaPart & mask, seqEnd - ((posStart++) + kmerSize) * 3 + 1 };
            }
        }
        return { UINT64_MAX, 0 }; // No more kmers found
    }
};


class KmerScanner_aa2aa : public KmerScanner {
protected:
    // Internal values
    std::vector<int8_t> aacids;
    uint64_t aaPart = 0;
    uint64_t mask;

public:
    KmerScanner_aa2aa(int k) : KmerScanner(k) {
        aacids.resize(256);
        for (int i = 0; i < 256; i++) {
            aacids[i] = 27; // Initialize all to -1
        }
        aacids['A'] = 0; // Alanine
        aacids['R'] = 1; // Arginine
        aacids['N'] = 2; // Asparagine
        aacids['D'] = 3; // Aspartic acid
        aacids['C'] = 4; // Cysteine
        aacids['Q'] = 5; // Glutamine
        aacids['E'] = 6; // Glutamic acid
        aacids['G'] = 7; // Glycine
        aacids['H'] = 8; // Histidine
        aacids['I'] = 9; // Isoleucine
        aacids['L'] = 10; // Leucine
        aacids['K'] = 11; // Lysine
        aacids['M'] = 12; // Methionine
        aacids['F'] = 13; // Phenylalanine
        aacids['P'] = 14; // Proline
        aacids['S'] = 15; // Serine
        aacids['T'] = 16; // Threonine
        aacids['W'] = 17; // Tryptophan
        aacids['Y'] = 18; // Tyrosine
        aacids['V'] = 19; // Valine
        aacids['B'] = 20; // Asparagine or Aspartic acid
        aacids['Z'] = 21; // Glutamine or Glutamic acid
        aacids['U'] = 22; // Selenocysteine
        aacids['O'] = 23; // Pyrrolysine
        aacids['*'] = 24; // Stop codon
        aacids['-'] = 25; // Gap or missing data
        aacids['.'] = 25; // Gap or missing data
        aacids['?'] = 25; // Gap or missing data
        aacids['X'] = 26; // Unknown amino acid


        if (kmerSize > 12 || kmerSize < 1) {
            std::cerr << "Error: k must be between 1 and 12, inclusive." << std::endl;
            exit(EXIT_FAILURE);
        } 
        this->mask = (1ULL << (5 * kmerSize)) - 1;
    };

    ~KmerScanner_aa2aa() {
        // std::cout << "KmerScanner destroyed." << std::endl;
    }

    void initScanner(
        const char * seq, 
        size_t seqStart, 
        size_t seqEnd, 
        bool isForward = true) override 
    {
        KmerScanner::initScanner(seq, seqStart, seqEnd, isForward);
        this->aaPart = 0;
    }

    virtual Kmer next() {
        int aa = 0;
        while (posStart <= seqLen - kmerSize) {
            bool sawN = false;
            loadedCharCnt -= (loadedCharCnt == kmerSize);
            while (loadedCharCnt < kmerSize) {
                aa = aacids[seq[seqStart + posStart + loadedCharCnt]];
                if (aa > 23) { sawN = true; break; }
                aaPart = (aaPart << 5) | (uint64_t)aa;
                loadedCharCnt++;
            }
            if (sawN) {
                posStart += loadedCharCnt + 1;
                aaPart = 0;
                loadedCharCnt = 0;
                continue;
            }
            return { aaPart & mask, seqStart + (posStart++) };
        }
        return { UINT64_MAX, 0 }; // No more kmers found
    }
};


#endif // METABULI_KMERSCANNER_H