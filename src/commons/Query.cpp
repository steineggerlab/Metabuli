//
// Created by 김재범 on 2022/03/18.
//

#include "Query.h"


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