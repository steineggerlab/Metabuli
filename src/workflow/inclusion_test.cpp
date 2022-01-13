//
// Created by 김재범 on 2021/05/10.
//

//#include "Classifier.h"
#include "NcbiTaxonomy.h"
#include "Parameters.h"
#include "LocalParameters.h"

#include <string>
#include <sstream>
#include <fstream>
#include <iostream>
#include <regex>

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

};

struct CountAtRank {
    int total;
    int FP;
    int TP;
    float precision;
    float sensitivity;
};
using namespace std;

void compareTaxon(TaxID shot, TaxID target, NcbiTaxonomy & ncbiTaxonomy, Counts & counts);
void compareTaxonAtRank(TaxID shot, TaxID target, NcbiTaxonomy & ncbiTaxonomy, CountAtRank & count, const string & rank);

int inclusiontest(int argc, const char **argv, const Command &command){

    LocalParameters &par = LocalParameters::getLocalInstance();
    par.parseParameters(argc, argv, command, false, Parameters::PARSE_ALLOW_EMPTY, 0);

    const string readClassificationFileName = par.filenames[0];

    string names = "../../gtdb_taxdmp/names.dmp";
    string nodes = "../../gtdb_taxdmp/nodes.dmp";
    string merged = "../../gtdb_taxdmp/merged.dmp";
    NcbiTaxonomy ncbiTaxonomy(names, nodes, merged);

    unordered_map<TaxID, unsigned int> taxCnt;
    for(int i = 0 ; i < ncbiTaxonomy.taxonNodes.size() ; i++) {
        taxCnt[ncbiTaxonomy.taxonNodes[i].taxId] = 1;
    }

    unordered_map<TaxID, TaxonCounts> cladeCnt = ncbiTaxonomy.getCladeCounts(taxCnt);


//    ///Load taxDB of kraken
//    unordered_map<int, int> child2parent;
//    string childString, parentString, throwaway;
//    int childInt, parentInt;
//    ifstream taxDB;
//    taxDB.open(krakenTaxDB);
//    if(taxDB.is_open()){
//        while(getline(taxDB,childString,'\t')){
//            getline(taxDB, parentString, '\t');
//            getline(taxDB, throwaway,'\n');
//            childInt = stoi(childString);
//            parentInt = stoi(parentString);
//            if(childInt > 1000000000)
//                child2parent[childInt] = parentInt;
//        }
//    } else{
//        cout<<"Cannot open taxDB"<<endl;
//    }
//    taxDB.close();

    ///Load the mapping file (assacc to taxID)
    const char * mappingFile = "../../gtdb_taxdmp/assacc_to_taxid_gtdb.tsv";
    unordered_map<string, int> assacc2taxid;
    string key, value;
    ifstream map;
    map.open(mappingFile);
    if(map.is_open()){
        while(getline(map,key,'\t')){
            getline(map, value, '\n');
            assacc2taxid[key] = stoi(value);
        }
    } else{
        cout<<"Cannot open file for mappig from assemlby accession to tax ID"<<endl;
    }
    map.close();


    ///read classification
    vector<int> rightAnswers;
    vector<int> classList;

    string classString;
    ifstream readClassification;
    readClassification.open(readClassificationFileName);
    vector<string> fields;
    string field;
    int classInt;

    regex regex1("(GC[AF]_[0-9]*\\.[0-9]*)");
    smatch assacc;

    while(getline(readClassification,classString,'\n')){
        istringstream lineStream(classString);
        fields.clear();
        while(getline(lineStream, field, '\t')){
            fields.push_back(field);
        }
        classInt = stoi(fields[2]);
        classList.push_back(classInt);

        regex_search(fields[1], assacc, regex1);
        rightAnswers.push_back(assacc2taxid[assacc[0]]);
    }
    cout<<"hi"<<endl;
    cout<<"num of classification: "<< classList.size()<<endl;

    ///Load query file -> name
    //regex regex1("(GC[AF]_[0-9]*\\.[0-9]*)");
    //
//    string queryName;
//    ifstream query;
//    query.open(queryFileName);
//    string queryLine;
//    vector<string> queryNameList;
//    while(getline(query,queryLine,'\n')){
//        if(queryLine[0] == '>'){
//            regex_search(queryLine, assacc, regex1);
//            queryNameList.push_back(assacc[0]);
//        }else{
//            continue;
//        }
//    }



//    ///right answer list
//    vector<int> rightAnswers;
//    for(size_t i = 0; i < queryNameList.size(); i++){
//        if (assacc2taxid.count(queryNameList[i])) {
//            rightAnswers.push_back(assacc2taxid[queryNameList[i]]);
//        } else{
//            cout << queryNameList[i] << " is not in the mapping file" << endl;
//            rightAnswers.push_back(-1);
//            continue;
//        }
//    }

    Counts counts = {0,0,0,0,0,0,0,0};
    CountAtRank SS;
    CountAtRank S;
    CountAtRank G;
    CountAtRank F;
    ///score the classification
    for(size_t i = 0; i < classList.size(); i++){
        cout<<i<<" ";
        compareTaxon(classList[i], rightAnswers[i], ncbiTaxonomy, counts);
        compareTaxonAtRank(classList[i], rightAnswers[i], ncbiTaxonomy, SS, "subspecies");
        compareTaxonAtRank(classList[i], rightAnswers[i], ncbiTaxonomy, S, "species");
        compareTaxonAtRank(classList[i], rightAnswers[i], ncbiTaxonomy, G, "genus");
        compareTaxonAtRank(classList[i], rightAnswers[i], ncbiTaxonomy, F, "family");
    }
    SS.precision = (float)SS.TP / (float)SS.total;
    S.precision = (float)S.TP / (float)S.total;
    G.precision = (float)G.TP / (float)G.total;
    F.precision = (float)F.TP / (float)F.total;

    SS.sensitivity = (float)SS.TP / classList.size();
    S.sensitivity = (float)S.TP / classList.size();
    G.sensitivity = (float)G.TP / classList.size();
    F.sensitivity = (float)F.TP / classList.size();

    cout<<"Num of queries: " << classList.size() << endl;
    cout<<"Num of classifications: "<< counts.classificationCnt << endl;
    cout<<"Num of correct classifications: "<<counts.correct<<endl;
    cout<<"Num of correct but too broad classifications: "<<counts.highRank<<endl;
    cout<<"classified/total = " << float(counts.classificationCnt)/float(classList.size()) << endl;
    cout<<"correct   /total = "<< float(counts.correct) / float(classList.size())<<endl;
    cout<<"correct   /classifications = "<<float(counts.correct) / float(counts.classificationCnt) <<endl;
    cout<<"high rank /classifications = "<<float(counts.highRank) / float(counts.classificationCnt) <<endl << endl;

    cout<<"Number of targets at each rank / correct classification / tries"<<endl;
    cout<<"Superkingdom: " << counts.superkingdomTargetNumber << " / " << counts.superkingdomCnt_correct << " / "<<counts.superkingdomCnt_try<<endl;
    cout<<"Phylum      : " << counts.phylumTargetNumber << " / " << counts.phylumCnt_correct << " / "<<counts.phylumCnt_try<<endl;
    cout<<"Class       : " << counts.classTargetNumber << " / " << counts.classCnt_correct << " / "<<counts.classCnt_try<<endl;
    cout<<"Order       : " << counts.orderTargetNumber<<" / "<<counts.orderCnt_correct<<" / "<<counts.orderCnt_try<<endl;
    cout<<"Family      : " << counts.familyTargetNumber << " / " << counts.familyCnt_correct << " / "<<counts.familyCnt_try<<endl;
    cout<<"Genus       : " << counts.genusTargetNumber<<" / "<<counts.genusCnt_correct<<" / "<<counts.genusCnt_try<<endl;
    cout<<"Species     : " << counts.speciesTargetNumber<<" / "<<counts.speciesCnt_correct<<" / "<<counts.speciesCnt_try<<endl;
    cout<<"Subspecies  : " << counts.subspeciesTargetNumber<<" / "<<counts.subspeciesCnt_correct<<" / "<<counts.subspeciesCnt_try<<endl;

    cout<<"False positive at each rank"<<endl;
    cout<<counts.fp_phylum<<endl;
    cout<<counts.fp_class<<endl;
    cout<<counts.fp_order<<endl;
    cout<<counts.fp_family<<endl;
    cout<<counts.fp_genus<<endl;
    cout<<counts.fp_species<<endl;
    cout<<counts.fp_subspecies<<endl;

    cout<<"NEW"<<endl;
    cout<<"Family      : " << F.total << " / " << F.TP << " / "<< F.FP << " / " << F.precision << F.sensitivity << endl;
    cout<<"Genus       : " << G.total << " / " << G.TP << " / "<< G.FP << " / " << G.precision << G.sensitivity << endl;
    cout<<"Species     : " << S.total << " / " << S.TP << " / "<< S.FP << " / " << S.precision << S.sensitivity << endl;
    cout<<"Subspecies  : " << SS.total << " / " << SS.TP << " / "<< SS.FP << " / " << SS.precision << SS.sensitivity << endl;

}

void compareTaxonAtRank(TaxID shot, TaxID target, NcbiTaxonomy & ncbiTaxonomy, CountAtRank & count, const string & rank) {
    const TaxonNode * shotNode = ncbiTaxonomy.taxonNode(shot);
    const TaxonNode * targetNode = ncbiTaxonomy.taxonNode(target);
    TaxID shotTaxIdAtRank = shotNode->taxId;
    TaxID targetTaxIdAtRank = targetNode->taxId;

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

void compareTaxon(TaxID shot, TaxID target, NcbiTaxonomy & ncbiTaxonomy, Counts& counts) { ///target: subspecies or species

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
        isCorrect = true;
        cout<<"O"<<endl;
    } else if(NcbiTaxonomy::findRankIndex(shotRank) <= NcbiTaxonomy::findRankIndex(targetRank)){ //classified into wrong taxon or too specifically
        cout<<"X"<<endl;
    } else { // classified at higher rank (too safe classification)
        if(shotRank == "superkingdom"){
            cout<<"X"<<endl;
        } else if(shot == ncbiTaxonomy.getTaxIdAtRank(target, shotRank)){ //on right branch
            counts.correct ++;
            cout<<"0"<<endl;
            isCorrect = true;
        } else{ //on wrong branch
            cout<<"X"<<endl;
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
