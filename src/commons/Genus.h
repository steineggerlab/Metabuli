//
// Created by 김재범 on 2022/03/18.
//

#ifndef ADCLASSIFIER2_GENUS_H
#define ADCLASSIFIER2_GENUS_H

#include <iostream>
#include "Species.h"
#include <vector>
#include "Match.h"
#include "NcbiTaxonomy.h"
#include <cstring>
#include "BitManipulateMacros.h"
#include <list>
class Genus {
private:
    TaxID genusID;
    int aminoAcidNum;
    int queryLength;
    float score;
    std::vector<Species> speciesList;
    signed char * hammingsAtEachPos;
    bool * posCheckList;
    std::list<Match> waitingMatches;

    void keepMatch(Match & match);
public:
    TaxID getGenusID() const { return genusID;}
    void takeMatch(Match & match);
    float getScore() const { return score;}
    Species& chooseSpecies();
    void scoreGenus();
    Genus(size_t queryLength, TaxID genusID);
    Genus(){};
};


#endif //ADCLASSIFIER2_GENUS_H
