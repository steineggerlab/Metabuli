#ifndef ADCLASSIFIER2_MATCH_H
#define ADCLASSIFIER2_MATCH_H

#include "Kmer.h"
#include <cstdint>
#include <iostream>
#include "BitManipulateMacros.h"

struct Match { // 20 byte
    Match(){}
    Match(QueryKmerInfo qInfo,
          int targetId,
          TaxID speciesId,
          uint16_t eachHamming,
          uint8_t hamming,
          bool redundancy):
          qInfo(qInfo), targetId(targetId), speciesId(speciesId),
          rightEndHamming(eachHamming), hamming(hamming), redundancy(redundancy) { }

    QueryKmerInfo qInfo; // 8
    TaxID targetId; // 4 taxonomy id infact
    TaxID speciesId; // 4
    uint16_t rightEndHamming; // 2
    uint8_t hamming; // 1
    bool redundancy; // 1

    void printMatch() const {
        std::cout << qInfo.sequenceID << " " << qInfo.pos << " " << qInfo.frame << " "
        << targetId << " " << speciesId << " " << rightEndHamming << " " << (int)hamming << " " << getScore() << std::endl;
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

    float getRightPartScore(const int range, float score = 0.0f, int cnt = 0) const {
        if (cnt == range) {
            return score;
        }
        int currentHamming = GET_2_BITS(rightEndHamming >> (14 - cnt * 2));
        if (currentHamming == 0) {
            score += 3.0f;
        } else {
            score += 2.0f - 0.5f * currentHamming;
        }
        return getRightPartScore(range, score, cnt + 1);    
    }

    float getLeftPartScore(const int range, float score = 0.0f, int cnt = 0) const {
        if (cnt == range) {
            return score;
        }
        int currentHamming = GET_2_BITS(rightEndHamming >> (cnt * 2));
        if (currentHamming == 0) {
            score += 3.0f;
        } else {
            score += 2.0f - 0.5f * currentHamming;
        }
        return getLeftPartScore(range, score, cnt + 1);    
    }
};

#endif //ADCLASSIFIER2_MATCH_H