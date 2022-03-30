//
// Created by 김재범 on 2022/03/18.
//

#include "Genus.h"

Genus::Genus(size_t queryLength, TaxID genusID) : genusID(genusID), queryLength(queryLength) {
    aminoAcidNum = (int) queryLength / 3;
    hammingsAtEachPos = new signed char[aminoAcidNum + 1];
    std::memset(hammingsAtEachPos, -1, (aminoAcidNum + 1));
    posCheckList = new bool[aminoAcidNum + 1];
    std::memset(posCheckList, false, sizeof(bool) * (aminoAcidNum + 1));
}

// TODO hamming distance
void Genus::takeMatch(Match &match) {
    int currPos = match.position / 3;

    if (currPos == 0) { // begin
        if (posCheckList[currPos + 1]) {
            // Keep match info for scoring
            keepMatch(match);
            return;
        }
    } else {
        if (posCheckList[currPos + 1] || posCheckList[currPos - 1]) {
            // Keep match info for scoring
            keepMatch(match);
            return;
        }
    }

    // Compare to the matches in the waiting list.
    bool haveNeighbor = false;
    for (auto m = waitingMatches.begin(); m != waitingMatches.end(); m++) {
        if (match.position - 3 <= m->position && m->position <= match.position + 3) {
            haveNeighbor = true;
            // Keep match info for scoring
            keepMatch(*m);
            waitingMatches.erase(m);
        }
    }

    if (haveNeighbor) {
        keepMatch(match);
    } else {
        waitingMatches.push_back(match);
    }
}

void Genus::keepMatch(Match &match) {
    int currPos = match.position / 3;
    uint16_t currHammings = match.rightEndHamming;
    if (GET_2_BITS(currHammings) > hammingsAtEachPos[currPos])
        __sync_lock_test_and_set(hammingsAtEachPos + currPos, GET_2_BITS(currHammings));
    if (GET_2_BITS(currHammings >> 2) > hammingsAtEachPos[currPos + 1])
        __sync_lock_test_and_set(hammingsAtEachPos + currPos + 1, GET_2_BITS(currHammings >> 2));
    if (GET_2_BITS(currHammings >> 4) > hammingsAtEachPos[currPos + 2])
        __sync_lock_test_and_set(hammingsAtEachPos + currPos + 2, GET_2_BITS(currHammings >> 4));
    if (GET_2_BITS(currHammings >> 6) > hammingsAtEachPos[currPos + 3])
        __sync_lock_test_and_set(hammingsAtEachPos + currPos + 3, GET_2_BITS(currHammings >> 6));
    if (GET_2_BITS(currHammings >> 8) > hammingsAtEachPos[currPos + 4])
        __sync_lock_test_and_set(hammingsAtEachPos + currPos + 4, GET_2_BITS(currHammings >> 8));
    if (GET_2_BITS(currHammings >> 10) > hammingsAtEachPos[currPos + 5])
        __sync_lock_test_and_set(hammingsAtEachPos + currPos + 5, GET_2_BITS(currHammings >> 10));
    if (GET_2_BITS(currHammings >> 12) > hammingsAtEachPos[currPos + 6])
        __sync_lock_test_and_set(hammingsAtEachPos + currPos + 6, GET_2_BITS(currHammings >> 12));
    if (GET_2_BITS(currHammings >> 14) > hammingsAtEachPos[currPos + 7])
        __sync_lock_test_and_set(hammingsAtEachPos + currPos + 7, GET_2_BITS(currHammings >> 14));


    // Update species
    for (Species sp: speciesList) {
        if (match.speciesTaxID == sp.getSpeciesID()) {
            sp.takeMatch(match);
            return;
        }
    }

    // Add species
    speciesList.emplace_back(Species(queryLength, match.speciesTaxID));
    speciesList.back().takeMatch(match);
}

void Genus::scoreGenus() {
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

Species &Genus::chooseSpecies() {
    float maxScore = 0;
    Species *bestSpecies;
    for (Species s: speciesList) {
        s.scoreSpecies();
        if (s.getScore() > maxScore) {
            maxScore = s.getScore();
            bestSpecies = &s;
        }
    }
    return *bestSpecies;
}