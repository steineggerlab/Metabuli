//
// Created by 김재범 on 2022/03/18.
//

#ifndef ADCLASSIFIER2_MATCH_H
#define ADCLASSIFIER2_MATCH_H

#include <iostream>

struct Match{ // 24(23) byte
    Match(){}
    Match(uint32_t queryId, int taxID, int speciesTaxID, int genusTaxID, int position, uint16_t eachHamming,
          uint8_t hamming, uint8_t frame)
            : queryId(queryId), taxID(taxID), speciesTaxID(speciesTaxID), genusTaxID(genusTaxID), position(position),
              rightEndHamming(eachHamming), hamming(hamming), frame(frame)  { }
    uint32_t queryId; // 4
    int taxID; // 4
    int speciesTaxID; // 4
    int genusTaxID; // 4
    int position; // 4
    uint16_t rightEndHamming; // 2
    uint8_t hamming; // 1
    uint8_t frame;
};

#endif //ADCLASSIFIER2_MATCH_H
