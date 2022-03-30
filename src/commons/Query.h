//
// Created by 김재범 on 2022/03/18.
//

#ifndef ADCLASSIFIER2_QUERY_H
#define ADCLASSIFIER2_QUERY_H

#include <iostream>
#include <string>
#include <NcbiTaxonomy.h>
#include <unordered_map>
#include <Genus.h>
#include <Match.h>
class Query {
private:
    int queryId;
    bool isClassified;
    std::string name;
    int taxId;
    float coverage;
    std::unordered_map<TaxID,int> taxCnt; ///how about using it for REPORTFILE? --> k-mer count
    size_t queryLength;
    std::vector<Genus> genusList;

    void getBestGenus();
public:
    void takeMatch(Match & match);
    void chooseBestTaxon2();
};


#endif //ADCLASSIFIER2_QUERY_H
