#ifndef ADCLASSIFIER2_MATCH_H
#define ADCLASSIFIER2_MATCH_H

#include "Kmer.h"
#include <iostream>

//struct Match_qInfo {
//    explicit Match_qInfo(uint32_t position = 0, uint32_t queryId = 0, uint8_t frame = 0)
//            : position(position), queryId(queryId),  frame(frame) {}
//    uint64_t position : 32;
//    uint64_t queryId : 29;
//    uint64_t frame : 3; // 0-5
//};

struct Match { // 24 byte
    Match(){}
    Match(QueryKmerInfo qInfo,
          int targetId,
          TaxID genusId,
          TaxID speciesId,
          uint16_t eachHamming,
          uint8_t hamming,
          bool redundancy):
          qInfo(qInfo), targetId(targetId), genusId(genusId), speciesId(speciesId),
          rightEndHamming(eachHamming), hamming(hamming), redundancy(redundancy) { }

    QueryKmerInfo qInfo; // 8
    TaxID targetId; // 4
    TaxID genusId; // 4
    TaxID speciesId; // 4
    uint16_t rightEndHamming; // 2
    uint8_t hamming; // 1
    bool redundancy; // 1

    void printMatch() const {
        std::cout << qInfo.sequenceID << " " << qInfo.pos << " " << qInfo.frame << " "
        << targetId << " " << genusId << " " << speciesId << " " << rightEndHamming << " " << (int)hamming << std::endl;
    }
};

#endif //ADCLASSIFIER2_MATCH_H
