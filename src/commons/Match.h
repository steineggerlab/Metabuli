#ifndef ADCLASSIFIER2_MATCH_H
#define ADCLASSIFIER2_MATCH_H

#include "Kmer.h"
#include <cstdint>
#include <iostream>
#include "BitManipulateMacros.h"

struct Match { // 24 byte
    Match(){}
    Match(QueryKmerInfo qInfo,
          int targetId,
          TaxID speciesId,
          uint32_t dnaEncoding,
          uint16_t eachHamming,
          uint8_t hamming):
          qInfo(qInfo), targetId(targetId), speciesId(speciesId), dnaEncoding(dnaEncoding),
          rightEndHamming(eachHamming), hamming(hamming) { }

    QueryKmerInfo qInfo;      // 8 // Query K-mer information
    TaxID targetId;           // 4 // axonomy id infact
    TaxID speciesId;          // 4 // Used to group matches by species
    uint32_t dnaEncoding;     // 4 // Used to check if two matches are consecutive
    uint16_t rightEndHamming; // 2 // Used to calculate score
    uint8_t hamming;          // 1 // Used to filter redundant matches

    void printMatch() const {
        std::cout << qInfo.sequenceID << " " << qInfo.pos << " " << qInfo.frame << " "
        << targetId << " " << speciesId << " " << rightEndHamming << " " << (int)hamming << " " << getScore() << "\n";
    }

    float getScore(float score = 0.0f, int cnt = 0) const { 
        int currentHamming = GET_2_BITS(rightEndHamming >> (cnt * 2));
        if (currentHamming == 0) {
            score += 3.0f;
        } else {
            score += 2.0f - 0.5f * currentHamming;
        }
        if (cnt == 7) {
            return score;
        } else {
        return getScore(score, cnt + 1);    
        }
    }

    virtual float getRightPartScore(const int range, float score = 0.0f, int cnt = 0) const {
        if (cnt == range) {
            return score;
        }
        int currentHamming = GET_2_BITS(rightEndHamming >> (cnt * 2));
        if (currentHamming == 0) {
            score += 3.0f;
        } else {
            score += 2.0f - 0.5f * currentHamming;
        }
        return getRightPartScore(range, score, cnt + 1);    
    }

    virtual float getLeftPartScore(const int range, float score = 0.0f, int cnt = 0) const {
        if (cnt == range) {
            return score;
        }
        int currentHamming = GET_2_BITS(rightEndHamming >> (14 - cnt * 2));
        if (currentHamming == 0) {
            score += 3.0f;
        } else {
            score += 2.0f - 0.5f * currentHamming;
        }
        return getLeftPartScore(range, score, cnt + 1);    
    }

    virtual int getRightPartHammingDist(const int range) const {
        int sum = 0;
        for (int i = 0; i < range; i++) {
            sum += GET_2_BITS(rightEndHamming >> (i * 2));
        }
        return sum;
    }

    virtual int getLeftPartHammingDist(const int range) const {
        int sum = 0;
        for (int i = 0; i < range; i++) {
            sum += GET_2_BITS(rightEndHamming >> (14 - i * 2));
        }
        return sum;
    }
};

struct Match_AA {
    uint32_t queryId;
    uint32_t targetId;
    uint32_t pos;     // For developing purpose only
    uint64_t kmer;    // For developing purpose only

    Match_AA(uint32_t queryId, uint32_t targetId) : queryId(queryId), targetId(targetId) { }

    Match_AA(uint32_t queryId, uint32_t targetId, uint64_t kmer) 
        : queryId(queryId), targetId(targetId), kmer(kmer) { }
    
    Match_AA(uint32_t queryId, uint32_t targetId, uint32_t pos, uint64_t kmer) 
        : queryId(queryId), targetId(targetId), pos(pos), kmer(kmer) { }

    static bool compare(const Match_AA &a, const Match_AA &b) {
        if (a.queryId != b.queryId)
            return a.queryId < b.queryId;
        if (a.pos != b.pos)
            return a.pos < b.pos;
        return a.targetId < b.targetId;
    }
};


struct MatchBlock {
    MatchBlock(size_t start, size_t end, int id) : start(start), end(end), id(id) {}
    MatchBlock() : start(0), end(0), id(0) {}
    size_t start;
    size_t end;
    uint32_t id;
};


#endif //ADCLASSIFIER2_MATCH_H