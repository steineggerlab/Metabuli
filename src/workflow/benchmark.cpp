//
// Created by 김재범 on 2022/04/11.
//
#include "benchmark.h"

using namespace std;
void compareTaxonAtRank(TaxID shot, TaxID target, NcbiTaxonomy & ncbiTaxonomy, CountAtRank & count, const string & rank) {
    const TaxonNode * shotNode = ncbiTaxonomy.taxonNode(shot);
    const TaxonNode * targetNode = ncbiTaxonomy.taxonNode(target);
    TaxID shotTaxIdAtRank = shotNode->taxId;
    TaxID targetTaxIdAtRank = targetNode->taxId;
    if(shot == 1) return;

    // Classification at higher rank -> ignore
    if(NcbiTaxonomy::findRankIndex(shotNode->rank) > NcbiTaxonomy::findRankIndex(rank)){
        return;
    }

    // If classification is at the lower rank, climb up the tree to the rank.
    if(NcbiTaxonomy::findRankIndex(shotNode->rank) < NcbiTaxonomy::findRankIndex(rank)){
        shotTaxIdAtRank = ncbiTaxonomy.getTaxIdAtRank(shotNode->taxId, rank);
    }
    if(NcbiTaxonomy::findRankIndex(targetNode->rank) < NcbiTaxonomy::findRankIndex(rank)){
        targetTaxIdAtRank = ncbiTaxonomy.getTaxIdAtRank(targetNode->taxId, rank);
    }

    // Correct classification at the rank
    if(shotTaxIdAtRank == targetTaxIdAtRank){
        count.TP++;
    } else {
        count.FP++;
    }
    count.total++;

    return;
}

void compareTaxon(TaxID shot, TaxID target, NcbiTaxonomy & ncbiTaxonomy, Counts& counts, vector<Score2> & scores,
                  float score) { ///target: subspecies or species

    const TaxonNode * shotNode = ncbiTaxonomy.taxonNode(shot);
    const TaxonNode * targetNode = ncbiTaxonomy.taxonNode(target);
    string shotRank = shotNode->rank;
    string targetRank = targetNode->rank;
    cout<<shot<<" "<<target<<" "<<shotRank<<" "<<targetRank<<" ";

    if(shot == 0){
        cout<<"X"<<endl;
        return;
    } else{
        counts.classificationCnt++;
    }

    bool isCorrect = false;
    if(shot == target){
        counts.correct ++;
        scores.emplace_back(1, shotRank, score);
        isCorrect = true;
        cout<<"O"<<endl;
    } else if(NcbiTaxonomy::findRankIndex(shotRank) <= NcbiTaxonomy::findRankIndex(targetRank)){ //classified into wrong taxon or too specifically
        cout<<"X"<<endl;
        scores.emplace_back(2, shotRank, score);
    } else { // classified at higher rank (too safe classification)
        if(shotRank == "superkingdom"){
            cout<<"X"<<endl;
        } else if(shot == ncbiTaxonomy.getTaxIdAtRank(target, shotRank)){ //on right branch
            counts.correct ++;
            cout<<"0"<<endl;
            isCorrect = true;
            scores.emplace_back(1, shotRank, score);
        } else{ //on wrong branch
            cout<<"X"<<endl;
            scores.emplace_back(2, shotRank, score);
        }
    }

    //count the number of classification at each rank
    if(shotRank == "subspecies") {
        counts.subspeciesCnt_try++;
    } else if(shotRank == "species") {
        counts.speciesCnt_try ++;
    } else if(shotRank == "genus"){
        counts.genusCnt_try ++;
    } else if(shotRank == "family"){
        counts.familyCnt_try++;
    } else if(shotRank == "order") {
        counts.orderCnt_try++;
    } else if(shotRank == "class") {
        counts.classCnt_try++;
    } else if(shotRank == "phylum") {
        counts.phylumCnt_try++;
    } else if(shotRank == "superkingdom"){
        counts.superkingdomCnt_try++;
    }

    if(!isCorrect) {
        const TaxonNode * LCAnode = ncbiTaxonomy.taxonNode(ncbiTaxonomy.LCA(shot, target));
        string LCArank = LCAnode->rank;
        if(LCArank == "species") {
            counts.fp_subspecies ++;
        } else if(LCArank == "genus"){
            counts.fp_species ++;
        } else if(LCArank == "family"){
            counts.fp_genus++;
        } else if(LCArank == "order") {
            counts.fp_family++;
        } else if(LCArank == "class") {
            counts.fp_order++;
        } else if(LCArank == "phylum") {
            counts.fp_class++;
        } else if(LCArank == "kingdom"){
            counts.fp_phylum++;
        }

        return;
    }

    //count the number of correct classification at each rank
    if(shotRank == "subspecies"){
        counts.subspeciesCnt_correct ++;
    } else if(shotRank == "species") {
        counts.speciesCnt_correct ++;
    } else if(shotRank == "genus"){
        counts.genusCnt_correct ++;
    } else if(shotRank == "family"){
        counts.familyCnt_correct++;
    } else if(shotRank == "order") {
        counts.orderCnt_correct++;
    } else if(shotRank == "class") {
        counts.classCnt_correct++;
    } else if(shotRank == "phylum") {
        counts.phylumCnt_correct++;
    } else if(shotRank == "superkingdom"){
        counts.superkingdomCnt_correct++;
    }

    //    const TaxonNode * shotNode = ncbiTaxonomy.taxonNode(shot);
//    string shotRank = shotNode->rank;
//    cout<<shot<<" "<<target<<" "<<shotRank<<" ";
//    if(shot == 0){
//        cout<<"X"<<endl;
//    }
//    if(NcbiTaxonomy::findRankIndex(shotRank) <= 3){
//        //cout<<"subspecies"<<endl;
//        if(shot == target){
//            cout<<"O"<<endl;
//            counts.subspCnt ++;
//        } else cout<<"X"<<endl;
//
//    } else if(shotRank == "species") {
//        //cout<<"species"<<endl;
//        if(shot == ncbiTaxonomy.getTaxIdAtRank(target, "species")){
//            counts.spCnt ++;
//            cout<<"O"<<endl;
//        }else cout<<"X"<<endl;
//    } else if(shotRank == "genus"){
//        //cout<<"genus"<<endl;
//        if(shot == ncbiTaxonomy.getTaxIdAtRank(target, "genus")){
//            counts.genusCnt ++;
//            cout<<"O"<<endl;
//        } else cout<<"X"<<endl;
//    } else if(shotRank == "family"){
//        //cout<<"family"<<endl;
//        if(shot == ncbiTaxonomy.getTaxIdAtRank(target, "family")) {
//            counts.familyCnt_correct++;
//        }
//    }else if(shotRank == "order") {
//        //cout<<"order"<<endl;
//        if(shot == ncbiTaxonomy.getTaxIdAtRank(target, "order")) {
//            counts.orderCnt_correct++;
//        }
//    }else if(shotRank == "class") {
//        //cout<<"class"<<endl;
//        if(shot == ncbiTaxonomy.getTaxIdAtRank(target, "class")) {
//            counts.classCnt_correct++;
//        }
//    } else if(shotRank == "phylum") {
//        //cout<<"phylum"<<endl;
//        if(shot == ncbiTaxonomy.getTaxIdAtRank(target, "phylum")) {
//            counts.phylumCnt_correct++;
//        }
//    } else {
//        return;
//    }
}
