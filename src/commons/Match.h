#ifndef ADCLASSIFIER2_MATCH_H
#define ADCLASSIFIER2_MATCH_H

#include <iostream>

struct Match { // 16 byte
    Match(){}
    Match(uint32_t queryId, int targetId, int position, uint16_t eachHamming, uint8_t hamming, bool redundancy)
            : queryId(queryId), targetId(targetId), position(position), rightEndHamming(eachHamming),
              hamming(hamming), redundacny(redundancy) { }
    uint32_t queryId; // 4
    int targetId; // 4
    int position; // 4
    uint16_t rightEndHamming; // 2
    uint8_t hamming; // 1
    bool redundacny; // 1
};

#endif //ADCLASSIFIER2_MATCH_H
