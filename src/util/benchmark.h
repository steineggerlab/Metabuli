//
// Created by 김재범 on 2022/04/11.
//

#ifndef ADCLASSIFIER2_BENCHMARK_H
#define ADCLASSIFIER2_BENCHMARK_H

#include <iostream>
#include <string>
#include <vector>
#include "NcbiTaxonomy.h"

struct Score2{
    Score2(int tf, std::string rank, float score) : tf(tf), rank(rank), score(score) { }
    int tf; // 1 = t, 2 = f
    std::string rank;
    float score;
};

struct Counts{
    int classificationCnt;
    int correct;
    int highRank;

    //number of targets at each rank
    int subspeciesTargetNumber;
    int speciesTargetNumber;
    int genusTargetNumber;
    int familyTargetNumber;
    int orderTargetNumber;
    int classTargetNumber;
    int phylumTargetNumber;
    int superkingdomTargetNumber;

    //number of classification at each rank
    int subspeciesCnt_try;
    int speciesCnt_try;
    int genusCnt_try;
    int familyCnt_try;
    int orderCnt_try;
    int classCnt_try;
    int phylumCnt_try;
    int superkingdomCnt_try;


    //number of correct classifications at each rank
    int subspeciesCnt_correct;
    int speciesCnt_correct;
    int genusCnt_correct;
    int familyCnt_correct;
    int orderCnt_correct;
    int classCnt_correct;
    int phylumCnt_correct;
    int superkingdomCnt_correct;

    //FP at each rank
    int fp_subspecies;
    int fp_species;
    int fp_genus;
    int fp_family;
    int fp_order;
    int fp_class;
    int fp_phylum;
    int fp_superkingdom;

    int queryCnt;
};

struct CountAtRank {
    int total;
    int FP;
    int TP;
    int FN;
    float precision;
    float sensitivity;
    float f1;
    void calculate() {
        precision = (float)TP / (float)(TP + FP);
        sensitivity = (float)TP / (float)(total);
        f1 = 2 * precision * sensitivity / (precision + sensitivity);
    }
};

void compareTaxon(TaxID shot, TaxID target, NcbiTaxonomy & ncbiTaxonomy, Counts & counts,
                  std::vector<Score2> & tpOrFp, float score);

void compareTaxonAtRank(TaxID shot, TaxID target, NcbiTaxonomy & ncbiTaxonomy, CountAtRank & count, const std::string & rank);

#endif //ADCLASSIFIER2_BENCHMARK_H
