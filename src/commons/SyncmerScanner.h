#ifndef METABULI_SYNCMER_SCANNER_H
#define METABULI_SYNCMER_SCANNER_H

#include <iostream>
#include <deque>

#include "KmerScanner.h"

class SyncmerScanner : public MetamerScanner {
protected:
    // Internal values
    int smerLen;
    uint64_t smerMask;

    // Variables for syncmer scanning
    std::deque<Kmer> dq;
    int smerCnt;
    uint64_t smer;
    int prevPos;

public:
    SyncmerScanner(int smerLen, const GeneticCode &geneticCode) : MetamerScanner(geneticCode) {
        // std::cout << "SyncmerScanner initialized with smerLen: " << smerLen << std::endl;
        this->smerLen = smerLen;
        this->smerMask = (1ULL << (5 * smerLen)) - 1;
    }

    void initScanner(const char * seq, size_t seqStart, size_t seqEnd, bool isForward) override {
        MetamerScanner::initScanner(seq, seqStart, seqEnd, isForward);
        this->dq.clear();
        this->smerCnt = 0;
        this->smer = 0;
        this->prevPos = -8;
    }

    Kmer next() override {
        bool syncmerFound = false;
        int aa = 0;
        while (posStart <= aaLen - 8 && !syncmerFound) {
            bool sawN = false;
            smerCnt -= (smerCnt > 0);
            while (smerCnt < 8 - smerLen + 1) {
                loadedCharCnt -= (loadedCharCnt == smerLen);
                while (loadedCharCnt < smerLen) {
                    if (isForward) {
                        int ci = seqStart + (posStart + smerCnt + loadedCharCnt) * 3;
                        aa = geneticCode.getAA(atcg[seq[ci]], atcg[seq[ci + 1]], atcg[seq[ci + 2]]);
                    } else {
                        int ci = seqEnd - (posStart + smerCnt + loadedCharCnt) * 3;
                        aa = geneticCode.getAA(iRCT[atcg[seq[ci]]], iRCT[atcg[seq[ci - 1]]], iRCT[atcg[seq[ci - 2]]]);
                    }
                    if (aa < 0) { sawN = true; break; }
                    smer = (smer << 5) | (uint64_t)aa;                    
                    loadedCharCnt++;
                }
                if (sawN) break;
                smer &= smerMask;
                while (!dq.empty() && dq.back().value > smer) dq.pop_back();
                dq.emplace_back(smer, posStart + smerCnt);
                smerCnt++;
            }
            if (sawN) {
                posStart += smerCnt + loadedCharCnt + 1;
                prevPos = posStart - 8; // Reset previous position ??
                dq.clear(); 
                smerCnt = loadedCharCnt = 0; 
                smer = 0;
                continue;
            }
            if (!dq.empty() && dq.front().pos < posStart) dq.pop_front();
            uint32_t anchor1 = static_cast<uint32_t>(posStart);
            uint32_t anchor2 = static_cast<uint32_t>(posStart + (kmerSize - smerLen));
            if (!dq.empty() && (dq.front().pos == anchor1 || dq.front().pos == anchor2)) {
                int shifts = posStart - prevPos;
                if (isForward) {
                    for (int i = 0; i < shifts; ++i) {
                        int ci = seqStart + (prevPos + 8 + i) * 3;
                        aaPart = (aaPart << 5) | (uint64_t)geneticCode.getAA(atcg[seq[ci]], atcg[seq[ci + 1]], atcg[seq[ci + 2]]);
                        dnaPart = (dnaPart << 3) | (uint64_t)geneticCode.getCodon(atcg[seq[ci]], atcg[seq[ci + 1]], atcg[seq[ci + 2]]);
                    }
                } else {
                    for (int i = 0; i < shifts; ++i) {
                        int ci = seqEnd - (prevPos + 8 + i) * 3;
                        aaPart = (aaPart << 5) | (uint64_t)geneticCode.getAA(iRCT[atcg[seq[ci]]], iRCT[atcg[seq[ci - 1]]], iRCT[atcg[seq[ci - 2]]]);
                        dnaPart = (dnaPart << 3) | (uint64_t)geneticCode.getCodon(iRCT[atcg[seq[ci]]], iRCT[atcg[seq[ci - 1]]], iRCT[atcg[seq[ci - 2]]]);
                    }
                }
                prevPos = posStart;
                syncmerFound = true;
            }
            ++posStart;
        }
        if (syncmerFound) {
            if (isForward) {
                return {(aaPart << 24) | (dnaPart & dnaMask), seqStart + prevPos * 3};
            } else {
                return {(aaPart << 24) | (dnaPart & dnaMask), seqEnd - (prevPos + 8) * 3 + 1};
            }
        } else {
            return {UINT64_MAX, 0}; // No more syncmers found
        }
    }
};

class SyncmerScanner_aa2aa : public KmerScanner_aa2aa {
protected:
    // Internal values
    int smerLen;
    uint64_t smerMask;

    // Variables for syncmer scanning
    std::deque<Kmer> dq;
    int smerCnt;
    uint64_t smer;
    int prevPos;

public:
    SyncmerScanner_aa2aa(int k, int s) : KmerScanner_aa2aa(k), smerLen(s) {
        smerMask = (1ULL << (5 * smerLen)) - 1;
    }

    ~SyncmerScanner_aa2aa() {}

    void initScanner(
        const char * seq, 
        size_t seqStart, 
        size_t seqEnd, 
        bool isForward = true) override 
    {
        KmerScanner_aa2aa::initScanner(seq, seqStart, seqEnd, isForward);
        this->dq.clear();
        this->smerCnt = 0;
        this->smer = 0;
        this->prevPos = -kmerSize;
    }

    Kmer next() override {
        bool syncymerFound = false;
        int aa = 0;
        while (posStart <= seqLen - kmerSize && !syncymerFound) {
            bool sawN = false;
            smerCnt -= (smerCnt > 0);
            // Fill the deque with s-mers
            while (smerCnt < kmerSize - smerLen + 1) {
                loadedCharCnt -= (loadedCharCnt == smerLen);
                while (loadedCharCnt < smerLen) {
                    aa = aacids[seq[seqStart + posStart + smerCnt + loadedCharCnt]];
                    if (aa > 23) { sawN = true; break; }
                    smer = (smer << 5) | (uint64_t)aa;
                    loadedCharCnt++;
                }
                if (sawN) break;
                smer &= smerMask;
                while (!dq.empty() && dq.back().value > smer) dq.pop_back();
                dq.emplace_back(smer, posStart + smerCnt);
                smerCnt++;
            }
            if (sawN) {
                posStart += smerCnt + loadedCharCnt + 1;
                prevPos = posStart - kmerSize;
                loadedCharCnt = 0;
                smerCnt = 0;
                smer = 0;
                dq.clear();
                continue;
            }

            // Remove s-mers that are out of the current k-mer window
            if (!dq.empty() && dq.front().pos < posStart) dq.pop_front();
            uint32_t anchor1 = static_cast<uint32_t>(posStart);
            uint32_t anchor2 = static_cast<uint32_t>(posStart + (kmerSize - smerLen));

            // Check if the minimum s-mer is at one of the anchor positions
            if (!dq.empty() && (dq.front().pos == anchor1 || dq.front().pos == anchor2)) {
                int shift = posStart - prevPos;
                for (int i = 0; i < shift; ++i) {
                    int aa = aacids[seq[seqStart + prevPos + kmerSize + i]];
                    aaPart = (aaPart << 5) | (uint64_t)aa;
                }
                prevPos = posStart;
                syncymerFound = true;
            }
            ++posStart;
        }
        if (syncymerFound) {
            return { aaPart & mask, seqStart + (posStart - 1) };
        }
        return { UINT64_MAX, 0 }; // No more kmers found
    }
};

class SyncmerScanner_dna2aa : public KmerScanner_dna2aa {
protected:
    // Internal values
    int smerLen;
    uint64_t smerMask;

    // Variables for syncmer scanning
    std::deque<Kmer> dq;
    int smerCnt;
    uint64_t smer;
    int prevPos;

public:
    SyncmerScanner_dna2aa(const GeneticCode &geneticCode, int k, int s) 
        : KmerScanner_dna2aa(geneticCode, k), smerLen(s) 
    {
        this->smerMask = (1ULL << (5 * smerLen)) - 1;
    }

    void initScanner(
        const char * seq, 
        size_t seqStart, 
        size_t seqEnd, 
        bool isForward = true) override 
    {
        KmerScanner_dna2aa::initScanner(seq, seqStart, seqEnd, isForward);
        this->dq.clear();
        this->smerCnt = 0;
        this->smer = 0;
        this->prevPos = -kmerSize;
    }

    Kmer next() override {
        bool syncymerFound = false;
        int aa = 0;
        while (posStart <= aaLen - kmerSize && !syncymerFound) {
            bool sawN = false;
            smerCnt -= (smerCnt > 0);
            // Fill the deque with s-mers
            while (smerCnt < kmerSize - smerLen + 1) {
                loadedCharCnt -= (loadedCharCnt == smerLen);
                while (loadedCharCnt < smerLen) {
                    if (isForward) {
                        int ci = seqStart + (posStart + smerCnt + loadedCharCnt) * 3;
                        aa = geneticCode.getAA(atcg[seq[ci]], atcg[seq[ci + 1]], atcg[seq[ci + 2]]);
                    } else {
                        int ci = seqEnd - (posStart + smerCnt + loadedCharCnt) * 3;
                        aa = geneticCode.getAA(iRCT[atcg[seq[ci]]], iRCT[atcg[seq[ci - 1]]], iRCT[atcg[seq[ci - 2]]]);            
                    }
                    if (aa < 0) { sawN = true; break; }
                    smer = (smer << 5) | (uint64_t)aa;
                    loadedCharCnt++;
                }
                if (sawN) break;
                smer &= smerMask;
                while (!dq.empty() && dq.back().value > smer) dq.pop_back();
                dq.emplace_back(smer, posStart + smerCnt);
                smerCnt++;
            }
            if (sawN) {
                posStart += smerCnt + loadedCharCnt + 1;
                prevPos = posStart - kmerSize;
                loadedCharCnt = 0;
                smerCnt = 0;
                smer = 0;
                dq.clear();
                continue;
            }

            // Remove s-mers that are out of the current k-mer window
            if (!dq.empty() && dq.front().pos < posStart) dq.pop_front();
            uint32_t anchor1 = static_cast<uint32_t>(posStart);
            uint32_t anchor2 = static_cast<uint32_t>(posStart + (kmerSize - smerLen));

            // Check if the minimum s-mer is at one of the anchor positions
            if (!dq.empty() && (dq.front().pos == anchor1 || dq.front().pos == anchor2)) {
                int shifts = posStart - prevPos;
                if (isForward) {
                    for (int i = 0; i < shifts; ++i) {
                        int ci = seqStart + (prevPos + kmerSize + i) * 3;
                        aaPart = (aaPart << 5) | (uint64_t)geneticCode.getAA(atcg[seq[ci]], atcg[seq[ci + 1]], atcg[seq[ci + 2]]);
                    }
                } else {
                    for (int i = 0; i < shifts; ++i) {
                        int ci = seqEnd - (prevPos + kmerSize + i) * 3;
                        aaPart = (aaPart << 5) | (uint64_t)geneticCode.getAA(iRCT[atcg[seq[ci]]], iRCT[atcg[seq[ci - 1]]], iRCT[atcg[seq[ci - 2]]]);
                    }
                }
                prevPos = posStart;
                syncymerFound = true;
            }
            ++posStart;
        }
        if (syncymerFound) {
            if (isForward) {
                return { aaPart & mask, seqStart + (prevPos * 3)};
            } else {
                return { aaPart & mask, seqEnd - (prevPos + kmerSize) * 3 + 1};
            }
        }
        return {UINT64_MAX, 0}; // No more kmers found
    }

};

#endif //METABULI_SYNCMER_SCANNER_H