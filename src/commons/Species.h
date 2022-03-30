//
// Created by 김재범 on 2022/03/18.
//

#ifndef ADCLASSIFIER2_SPECIES_H
#define ADCLASSIFIER2_SPECIES_H

#include <iostream>
#include <iostream>
#include "NcbiTaxonomy.h"
#include "Match.h"
#include "BitManipulateMacros.h"
#include <cstring>
#include <unordered_map>
#include <algorithm>

class Species {
private:
    struct Subspecies {
        TaxID id;
        int pos;
        Subspecies(TaxID id ,int pos) : id(id), pos(pos) {};
    };
    int queryLength;
    TaxID speciesId;
    float score;
    signed char *hammingsAtEachPos;
    int aminoAcidNum;
    std::unordered_map<TaxID, int> subspeciesCnt;
    signed char *posCheckList;
    std::vector<Subspecies> subspeciesList;

    void countSubspecies();

    static bool sortSubspecies(const Subspecies & a, const Subspecies & b);
public:
    void scoreSpecies();

    TaxID getSpeciesID();

    void takeMatch(Match &match);

    void keepMatch(Match &match);

    Species(size_t queryLength, TaxID speciesID);

    float getScore();

    TaxID chooseSubspecies();


};


#endif //ADCLASSIFIER2_SPECIES_H
