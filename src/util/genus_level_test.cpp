#include "Parameters.h"
#include "LocalParameters.h"
#include "NcbiTaxonomy.h"
#include <fstream>
#include <iostream>
#include <regex>
#include "benchmark.h"
#include <string>
#include <sstream>

using namespace std;

void compareTaxon(TaxID shot, TaxID target, NcbiTaxonomy & ncbiTaxonomy, CountAtRank & counts);
void compareTaxon_exclusion(TaxID shot, TaxID target, NcbiTaxonomy & ncbiTaxonomy, Counts & counts, float score,
                   vector<Score2> & Scores);

int genus_level_test(int argc, const char **argv, const Command &command){

    LocalParameters &par = LocalParameters::getLocalInstance();
    par.parseParameters(argc, argv, command, false, Parameters::PARSE_ALLOW_EMPTY, 0);

    const string readClassificationFileName = par.filenames[0];
    const string mappingFile = par.filenames[1];
    const string taxonomy = par.filenames[2];
    string scoreFileName = readClassificationFileName + "_SC";

    string names = taxonomy + "/names.dmp";
    string nodes =  taxonomy + "/nodes.dmp";
    string merged =  taxonomy + "/merged.dmp";
    NcbiTaxonomy ncbiTaxonomy(names, nodes, merged);

    Counts counts{};

    cout<<"Load the mapping file"<<endl;
    // Load the mapping file (assacc to taxID)
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


    // read classification
    vector<int> rightAnswers;
    vector<int> classList;
    string classString;
    ifstream readClassification;
    readClassification.open(readClassificationFileName);
    vector<string> fields;
    string field;
    int classification;
    regex regex1("(GC[AF]_[0-9]*\\.[0-9]*)");
    smatch assacc;
    while(getline(readClassification,classString,'\n')){
        istringstream lineStream(classString);
        fields.clear();
        while(getline(lineStream, field, '\t')){
            fields.push_back(field);
        }
        // 3rd field -> classification
        classification = stoi(fields[2]);
        classList.push_back(classification);

        // 2nd field -> assacc
        regex_search(fields[1], assacc, regex1);
        rightAnswers.push_back(ncbiTaxonomy.getTaxIdAtRank(assacc2taxid[assacc[0]], "genus"));
    }

    cout<<"The number of queries: "<< classList.size()<<endl;

    CountAtRank G = {};
    for(size_t i = 0; i < classList.size(); i++){
        compareTaxon(classList[i], rightAnswers[i], ncbiTaxonomy, G);
    }

    G.sensitivity = (float) G.TP / (float) classList.size();
    G.precision = (float) G.TP / (float) G.total;
    float f1 = (2 * G.sensitivity * G.precision) / (G.sensitivity + G.precision);


//    cout<<"Num of queries: " << counts.queryCnt << endl;
//    cout<<"Num of classifications: "<< counts.classificationCnt << endl;
//    cout<<"Num of correct classifications: "<<counts.correct<<endl;
//    cout<<"Num of correct but too broad classifications: "<<counts.highRank<<endl;
//    cout<<"classified/total = " << float(counts.classificationCnt)/float(counts.queryCnt) << endl;
//    cout<<"correct   /total (Sensitivity) = "<< float(counts.correct) / float(counts.queryCnt)<<endl;
//    cout<<"correct   /classifications (PPV) = "<<float(counts.correct) / float(counts.classificationCnt) <<endl;
//    cout<<"high rank /classifications = "<<float(counts.highRank) / float(counts.classificationCnt) <<endl << endl;
//
//    cout<<"Number of targets at each rank / correct classification / tries"<<endl;
//    cout<<"Superkingdom: " << counts.superkingdomTargetNumber << " / " << counts.superkingdomCnt_correct << " / "<<counts.superkingdomCnt_try<<endl;
//    cout<<"Phylum      : " << counts.phylumTargetNumber << " / " << counts.phylumCnt_correct << " / "<<counts.phylumCnt_try<<endl;
//    cout<<"Class       : " << counts.classTargetNumber << " / " << counts.classCnt_correct << " / "<<counts.classCnt_try<<endl;
//    cout<<"Order       : " << counts.orderTargetNumber<<" / "<<counts.orderCnt_correct<<" / "<<counts.orderCnt_try<<endl;
//    cout<<"Family      : " << counts.familyTargetNumber << " / " << counts.familyCnt_correct << " / "<<counts.familyCnt_try<<endl;
//    cout<<"Genus       : " << counts.genusTargetNumber<<" / "<<counts.genusCnt_correct<<" / "<<counts.genusCnt_try<<endl;
//    cout<<"Species     : " << counts.speciesTargetNumber<<" / "<<counts.speciesCnt_correct<<" / "<<counts.speciesCnt_try<<endl;
//    cout<<"Subspecies  : " << counts.subspeciesTargetNumber<<" / "<<counts.subspeciesCnt_correct<<" / "<<counts.subspeciesCnt_try<<endl;

    cout<<"Genus     : " << G.total << " / " << G.TP << " / "<< G.FP << " / " << G.precision << " / "<< G.sensitivity << " / " << f1 << endl;
    return 0;
}

void compareTaxon(TaxID shot, TaxID target, NcbiTaxonomy & ncbiTaxonomy, Counts & counts, float score, vector<Score2> & Scores) {
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
        Scores.emplace_back(1, shotRank, score);
    } else if(NcbiTaxonomy::findRankIndex(shotRank) <= NcbiTaxonomy::findRankIndex(targetRank)){ //classified into wrong taxon or too specifically
        Scores.emplace_back(2, shotRank, score);
        cout<<"X"<<endl;
    } else { // classified at higher rank (too safe classification)
        if(shotRank == "superkingdom"){
            cout<<"X"<<endl;
        } else if(shot == ncbiTaxonomy.getTaxIdAtRank(target, shotRank)){ //on right branch
            counts.highRank ++;
            Scores.emplace_back(2, shotRank, score);
            cout<<"U"<<endl;
        } else{ //on wrong branch
            cout<<"X"<<endl;
            Scores.emplace_back(2, shotRank, score);
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

    if(!isCorrect) return;

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
}
void compareTaxon(TaxID shot, TaxID target, NcbiTaxonomy & ncbiTaxonomy, CountAtRank & counts) {
    if (shot == 0) {
        return;
    }

    if (shot == target) {
        counts.TP++;
        counts.total++;
    } else {
        counts.FP++;
        counts.total++;
    }
}
//    const TaxonNode * shotNode = ncbiTaxonomy.taxonNode(shot);
//    const TaxonNode * targetNode = ncbiTaxonomy.taxonNode(target);
//    string shotRank = shotNode->rank;
//    string targetRank = targetNode->rank;
////    cout<<shot<<" "<<target<<" "<<shotRank<<" "<<targetRank<<" ";
//    if(shot == 0){
//        cout<<"X"<<endl;
//        return;
//    } else{
//        counts.classificationCnt++;
//    }
//
//    bool isCorrect = false;
//    if(shot == target){
//        counts.correct ++;
//        isCorrect = true;
//        cout<<"O"<<endl;
//    } else if(NcbiTaxonomy::findRankIndex(shotRank) <= NcbiTaxonomy::findRankIndex(targetRank)){ // classified into wrong taxon or too specifically
//        cout<<"X"<<endl;
//    } else { // classified at higher rank (too safe classification)
//        if(shotRank == "superkingdom"){
//            cout<<"X"<<endl;
//        } else if(shot == ncbiTaxonomy.getTaxIdAtRank(target, shotRank)){ // on right branch
//            counts.highRank ++;
//            cout<<"U"<<endl;
//        } else{ //on wrong branch
//            cout<<"X"<<endl;
//        }
//    }
//
//    //count the number of classification at each rank
//    if(shotRank == "subspecies") {
//        counts.subspeciesCnt_try++;
//    } else if(shotRank == "species") {
//        counts.speciesCnt_try ++;
//    } else if(shotRank == "genus"){
//        counts.genusCnt_try ++;
//    } else if(shotRank == "family"){
//        counts.familyCnt_try++;
//    } else if(shotRank == "order") {
//        counts.orderCnt_try++;
//    } else if(shotRank == "class") {
//        counts.classCnt_try++;
//    } else if(shotRank == "phylum") {
//        counts.phylumCnt_try++;
//    } else if(shotRank == "superkingdom"){
//        counts.superkingdomCnt_try++;
//    }
//
//    if(!isCorrect) return;
//
//    //count the number of correct classification at each rank
//    if(shotRank == "subspecies"){
//        counts.subspeciesCnt_correct ++;
//    } else if(shotRank == "species") {
//        counts.speciesCnt_correct ++;
//    } else if(shotRank == "genus"){
//        counts.genusCnt_correct ++;
//    } else if(shotRank == "family"){
//        counts.familyCnt_correct++;
//    } else if(shotRank == "order") {
//        counts.orderCnt_correct++;
//    } else if(shotRank == "class") {
//        counts.classCnt_correct++;
//    } else if(shotRank == "phylum") {
//        counts.phylumCnt_correct++;
//    } else if(shotRank == "superkingdom"){
//        counts.superkingdomCnt_correct++;
//    }

