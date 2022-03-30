//
// Created by 김재범 on 2022/03/18.
//

#include "Query.h"

void Query::takeMatch(Match & match) {
    for(Genus genus : genusList){
        if(match.genusTaxID == genus.getGenusID()){
            genus.takeMatch(match);
            return;
        }
    }
    genusList.emplace_back(Genus(queryLength, match.genusTaxID));
    genusList.back().takeMatch(match);
}

void Query::chooseBestTaxon2() {

    //Choose genus
    Genus * bestGenus;
    float bestScore = 0;
    for(Genus genus : genusList){
        genus.scoreGenus();
        if(genus.getScore() > bestScore){
            bestGenus = & genus;
        }
    }

    if(bestScore < 0.8f){
        // Classify in genus
        return;
    }

    // Choose species
    Species & sp = bestGenus->chooseSpecies();
    taxId = sp.chooseSubspecies();
    // Choose subspecies

}