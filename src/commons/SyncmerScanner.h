#ifndef METABULI_SYNCMER_SCANNER_H
#define METABULI_SYNCMER_SCANNER_H

#include <iostream>
#include <deque>

#include "Kmer.h"
#include "GeneticCode.h"
#include "common.h"

struct Syncmer {
    uint64_t value;
    uint32_t pos;
    void printAA(const GeneticCode & code) const {
        uint64_t aaPart = value >> 24;
        for (int i = 0; i < 8; ++i) {
            int aa = (aaPart >> (35 - 5 * i)) & 0x1F;
            cout << code.aminoacids[aa];
        }
    }

    void printDNA(const GeneticCode & code) const {
        uint64_t dnaPart = value & 0xFFFFFF;
        uint64_t aaPart = value >> 24;
        for (int i = 0; i < 8; ++i) {
            int aa = (aaPart >> (35 - 5 * i)) & 0x1F;
            int codon = (dnaPart >> (21 - 3 * i)) & 0x7;
            cout << code.aa2codon[aa][codon];
        }
    }
};

class SyncmerScanner {
private:
    // Internal values
    const GeneticCode &geneticCode;
    int smerLen;
    uint64_t smerMask;
    uint64_t dnaMask;

    // Variables for syncmer scanning
    const char *seq;
    size_t seqStart;
    size_t seqEnd;
    size_t seqLen;
    size_t aaLen;

    uint64_t dnaPart;
    uint64_t aaPart;
    std::deque<Smer> dq;
    int smerCnt, loadedCharCnt;
    uint64_t smer;
    int prevPos;
    int posStart;
    uint64_t syncmer;

public:
    SyncmerScanner(int smerLen, const GeneticCode &geneticCode) : geneticCode(geneticCode) {
        this->smerLen = smerLen;
        this->smerMask = (1ULL << (5 * smerLen)) - 1;
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
        this->dq.clear();
        this->smerCnt = 0;
        this->loadedCharCnt = 0;
        this->smer = 0;
        this->prevPos = -8;
        this->posStart = 0;
    }

    Syncmer getNextSyncmer(bool forward) {
        bool syncmerFound = false;
        syncmer = 0;
        int aa = 0;
        while (posStart <= aaLen - 8 && !syncmerFound) {
            bool sawN = false;
            smerCnt -= (smerCnt > 0);
            while (smerCnt < 8 - smerLen + 1) {
                loadedCharCnt -= (loadedCharCnt == smerLen);
                while (loadedCharCnt < smerLen) {
                    if (forward) {
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
                dq.clear(); smerCnt = loadedCharCnt = 0; smer = 0;
                continue;
            }
            if (!dq.empty() && dq.front().pos < posStart) dq.pop_front();
            int anchor1 = posStart;
            int anchor2 = posStart + (8 - smerLen);
            if (!dq.empty() && (dq.front().pos == anchor1 || dq.front().pos == anchor2)) {
                int shifts = posStart - prevPos;
                if (forward) {
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
                syncmer = (aaPart << 24) | (dnaPart & dnaMask);
                syncmerFound = true;
            }
            ++posStart;
        }
        if (syncmerFound) {
            if (forward) {
                return {syncmer, seqStart + prevPos * 3};
            } else {
                return {syncmer, seqEnd - (prevPos + 8) * 3 + 1};
            }
        } else {
            return {UINT64_MAX, 0}; // No more syncmers found
        }
    }


    Syncmer getNextSyncmer_forward() {
        bool syncmerFound = false;
        syncmer = 0;
        while (posStart <= aaLen - 8 && !syncmerFound) {
            bool sawN = false;
            smerCnt -= (smerCnt > 0);
            while (smerCnt < 8 - smerLen + 1) {
                loadedCharCnt -= (loadedCharCnt == smerLen);
                while (loadedCharCnt < smerLen) {
                    int ci = seqStart + (posStart + smerCnt + loadedCharCnt) * 3;
                    int aa = geneticCode.getAA(atcg[seq[ci]], atcg[seq[ci + 1]], atcg[seq[ci + 2]]);
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
                dq.clear(); smerCnt = loadedCharCnt = 0; smer = 0;
                continue;
            }
            if (!dq.empty() && dq.front().pos < posStart) dq.pop_front();
            int anchor1 = posStart;
            int anchor2 = posStart + (8 - smerLen);
            if (!dq.empty() && (dq.front().pos == anchor1 || dq.front().pos == anchor2)) {
                int shifts = posStart - prevPos;
                for (int i = 0; i < shifts; ++i) {
                    int ci = seqStart + (prevPos + 8 + i) * 3;
                    aaPart = (aaPart << 5) | (uint64_t)geneticCode.getAA(atcg[seq[ci]], atcg[seq[ci + 1]], atcg[seq[ci + 2]]);
                    dnaPart = (dnaPart << 3) | (uint64_t)geneticCode.getCodon(atcg[seq[ci]], atcg[seq[ci + 1]], atcg[seq[ci + 2]]);
                }
                prevPos = posStart;
                syncmer = (aaPart << 24) | (dnaPart & dnaMask);
                syncmerFound = true;
            }
            ++posStart;
        }
        if (syncmerFound) {
            return {syncmer, seqStart + prevPos * 3};
        } else {
            return {UINT64_MAX, 0}; // No more syncmers found
        }
    }


    Syncmer getNextSyncmer_reverse() {
        bool syncmerFound = false;
        syncmer = 0;
        while (posStart <= aaLen - 8 && !syncmerFound) {
            bool sawN = false;
            smerCnt -= (smerCnt > 0);
            while (smerCnt < 8 - smerLen + 1) {
                loadedCharCnt -= (loadedCharCnt == smerLen);
                while (loadedCharCnt < smerLen) {                    
                    int ci = seqEnd - (posStart + smerCnt + loadedCharCnt) * 3;
                    int aa = geneticCode.getAA(iRCT[atcg[seq[ci]]], iRCT[atcg[seq[ci - 1]]], iRCT[atcg[seq[ci - 2]]]);                    
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
                dq.clear(); smerCnt = loadedCharCnt = 0; smer = 0;
                continue;
            }
            if (!dq.empty() && dq.front().pos < posStart) dq.pop_front();
            int anchor1 = posStart;
            int anchor2 = posStart + (8 - smerLen);
            if (!dq.empty() && (dq.front().pos == anchor1 || dq.front().pos == anchor2)) {
                int shifts = posStart - prevPos;
                for (int i = 0; i < shifts; ++i) {
                    int ci = seqEnd - (prevPos + 8 + i) * 3;
                    aaPart = (aaPart << 5) | (uint64_t)geneticCode.getAA(iRCT[atcg[seq[ci]]], iRCT[atcg[seq[ci - 1]]], iRCT[atcg[seq[ci - 2]]]);
                    dnaPart = (dnaPart << 3) | (uint64_t)geneticCode.getCodon(iRCT[atcg[seq[ci]]], iRCT[atcg[seq[ci - 1]]], iRCT[atcg[seq[ci - 2]]]);
                }
                prevPos = posStart;
                syncmer = (aaPart << 24) | (dnaPart & dnaMask);
                syncmerFound = true;
            }
            ++posStart;
        }
        if (syncmerFound) {
            return {syncmer, seqEnd - (prevPos + 8) * 3 + 1};
        } else {
            return {UINT64_MAX, 0}; // No more syncmers found
        }
    }

};

#endif //METABULI_SYNCMER_SCANNER_H