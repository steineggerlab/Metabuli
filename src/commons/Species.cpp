//
// Created by 김재범 on 2022/03/18.
//

#include "Species.h"

Species::Species(size_t queryLength, TaxID speciesID) : queryLength(queryLength), speciesId(speciesID) {
    aminoAcidNum = (int) queryLength / 3;
    hammingsAtEachPos = new signed char[aminoAcidNum + 1];
    std::memset(hammingsAtEachPos, 10, (aminoAcidNum + 1));
    posCheckList = new signed char[aminoAcidNum + 1];
    std::memset(posCheckList, -1, (aminoAcidNum + 1));
}

TaxID Species::getSpeciesID() {
    return this->speciesId;
}

void Species::takeMatch(Match &match) {
    int currPos = match.position / 3;

    keepMatch(match);
    if (match.taxID == match.speciesTaxID) {
        posCheckList[currPos] = 1;
    } else {
        subspeciesList.emplace_back(match.taxID, currPos);
    }
//
//
//    // First match at the position
//    if (posCheckList[currPos] == -1) {
//        keepMatch(match);
//        if(match.speciesTaxID == match.taxID) { // species level
//            posCheckList[currPos] = 1;
//        } else { // subspecies level
//            posCheckList[currPos] = 2;
//            subspeciesCnt[match.taxID] ++;
//        }
//    }
//
//    // Stain match was already there
//    if (posCheckList[currPos] == 2) {
//        keepMatch(match);
//        if(match.taxID == match.speciesTaxID) { // species level
//            subspeciesCnt[]
//        }
//    }


}

void Species::keepMatch(Match &match) {
    int currPos = match.position / 3;
    uint16_t currHammings = match.rightEndHamming;
    if (GET_2_BITS(currHammings) < hammingsAtEachPos[currPos])
        hammingsAtEachPos[currPos] = GET_2_BITS(currHammings);
    if (GET_2_BITS(currHammings >> 2) < hammingsAtEachPos[currPos + 1])
        hammingsAtEachPos[currPos + 1] = GET_2_BITS(currHammings >> 2);
    if (GET_2_BITS(currHammings >> 4) < hammingsAtEachPos[currPos + 2])
        hammingsAtEachPos[currPos + 2] = GET_2_BITS(currHammings >> 4);
    if (GET_2_BITS(currHammings >> 6) < hammingsAtEachPos[currPos + 3])
        hammingsAtEachPos[currPos + 3] = GET_2_BITS(currHammings >> 6);
    if (GET_2_BITS(currHammings >> 8) < hammingsAtEachPos[currPos + 4])
        hammingsAtEachPos[currPos + 4] = GET_2_BITS(currHammings >> 8);
    if (GET_2_BITS(currHammings >> 10) < hammingsAtEachPos[currPos + 5])
        hammingsAtEachPos[currPos + 5] = GET_2_BITS(currHammings >> 10);
    if (GET_2_BITS(currHammings >> 12) < hammingsAtEachPos[currPos + 6])
        hammingsAtEachPos[currPos + 6] = GET_2_BITS(currHammings >> 12);
    if (GET_2_BITS(currHammings >> 14) < hammingsAtEachPos[currPos + 7])
        hammingsAtEachPos[currPos + 7] = GET_2_BITS(currHammings >> 14);

}

void Species::scoreSpecies() {
    float hammingSum = 0;
    int coveredPosCnt = 0;
    for (int h = 0; h < aminoAcidNum; h++) {
        if (hammingsAtEachPos[h] == 0) { // Add 0 for 0 hamming dist.
            coveredPosCnt++;
        } else if (hammingsAtEachPos[h] != -1) { // Add 1.5, 2, 2.5 for 1, 2, 3 hamming dist. respectively
            hammingSum += 1.0f + (0.5f * hammingsAtEachPos[h]);
            coveredPosCnt++;
        }
    }

    // Score current genus
    int coveredLength = coveredPosCnt * 3;
    if (coveredLength > queryLength) coveredLength = queryLength;
    score = ((float) coveredLength - hammingSum) / (float) queryLength;
}

float Species::getScore() { return score; }

void Species::countSubspecies() {
    std::sort(subspeciesList.begin(), subspeciesList.end(), Species::sortSubspecies);

    bool overlap;
    size_t subMatchNum = subspeciesList.size();
    if(subMatchNum == 0) {
        return;
    } else if (subMatchNum == 1) {
        subspeciesCnt[subspeciesList[0].id]++;
        return;
    }

    for (size_t i = 0; i + 1 < subMatchNum; i++) {
        overlap = false;
        while(i + 1 < subMatchNum && subspeciesList[i].pos == subspeciesList[i+1].pos){
            overlap = true;
            i ++;
        }
        if(overlap){
            continue;
        } else {
            subspeciesCnt[subspeciesList[i].id]++;
        }
    }

    if(subspeciesList[subMatchNum - 1].pos != subspeciesList[subMatchNum - 2].pos){
        subspeciesCnt[subspeciesList[subMatchNum - 1].id]++;
    }
}

bool Species::sortSubspecies(const Subspecies &a, const Subspecies &b) {
    return a.pos < b.pos;
}

TaxID Species::chooseSubspecies() {
    countSubspecies();
    if (subspeciesCnt.empty()) {
        return speciesId;
    }

    int maxCnt = 0;
    int subCnt = 0;
    int minCnt = 1;
    TaxID strain;
    for (auto sub = subspeciesCnt.begin(); sub != subspeciesCnt.end(); sub++) {
        if (sub->second > minCnt) {
            subCnt++;
            if (sub->second > maxCnt) {
                maxCnt = sub->second;
                strain = sub->first;
            }
        }
    }

    if (subCnt == 1 && maxCnt > maxCnt + 1) {
        return strain;
    } else {
        return speciesId;
    }
}