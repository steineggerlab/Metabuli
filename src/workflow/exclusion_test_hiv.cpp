//
// Created by 김재범 on 2021/7/12.
//

//#include "Classifier.h"
#include "Parameters.h"
#include "LocalParameters.h"
#include "NcbiTaxonomy.h"
#include <fstream>
#include <iostream>
#include <regex>
#include <sstream>

using namespace std;
struct Counts13{
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
};

void compareTaxon112(TaxID shot, TaxID target, NcbiTaxonomy & ncbiTaxonomy, Counts13 & counts);



int exclusiontest_hiv(int argc, const char **argv, const Command &command){

    LocalParameters &par = LocalParameters::getLocalInstance();
    par.parseParameters(argc, argv, command, false, Parameters::PARSE_ALLOW_EMPTY, 0);

    const string readClassificationFileName = par.filenames[0];
    const string mappingFile = par.filenames[1];
    const string taxonomy = par.filenames[2];

    string names = taxonomy + "/names.dmp";
    string nodes =  taxonomy + "/nodes.dmp";
    string merged =  taxonomy + "/merged.dmp";
    NcbiTaxonomy ncbiTaxonomy(names, nodes, merged);

    Counts13 counts = {};

    // 1) Load mapping file
    cout<<"Load mapping from accession ID to taxonomy ID"<<endl;
    unordered_map<string, int> acc2taxid;
    string eachLine;
    string eachItem;
    ifstream map;
    map.open(mappingFile);
    vector<string> items;
    if(map.is_open()){
        while(getline(map,eachLine,'\n')){
            istringstream ss(eachLine);
            while (getline(ss, eachItem, '\t')){
                items.push_back(eachItem);
            }
            acc2taxid[items[0]] = stoi(items[1]);
            items.clear();
        }
    } else{
        cout<<"Cannot open file for mapping from accession to tax ID"<<endl;
    }
    map.close();

    // Load read classification
    vector<int> rightAnswers;
    vector<int> classifications;
    vector<float> scores;
    string classificationString;
    vector<string> fields;
    string field;
    int taxID;
    int rightAnswer;
    string accession;
    smatch assacc;
    ifstream readClassification;
    readClassification.open(readClassificationFileName);
    while(getline(readClassification, classificationString, '\n')){
        istringstream lineStream(classificationString);
        fields.clear();
        // 3rd field -> classification
        while(getline(lineStream, field, '\t')){
            fields.push_back(field);
        }
        taxID = stoi(fields[2]);
        classifications.push_back(ncbiTaxonomy.getTaxIdAtRank(taxID, par.testRank));

        // 2nd field -> accession
        accession = fields[1];
        int pos = accession.find(".");
        accession = accession.substr(0,pos+2);

        // 5th field -> score
        scores.push_back(stof(fields[4]));

        // Accession to right answer
        rightAnswer = ncbiTaxonomy.getTaxIdAtRank(acc2taxid[accession], par.testRank);
        rightAnswers.push_back(rightAnswer);

    }
    cout << "Num of classification: " << classifications.size() << endl;

    // Score the classification
    ofstream fp, tp;
    fp.open(readClassificationFileName + "_FP");
    tp.open(readClassificationFileName + "_TP");
    for(size_t i = 0; i < classifications.size(); i++){
        if(classifications[i] == rightAnswers[i]){
            tp << scores[i] << "\n";
        } else {
            fp << scores[i] << "\n";
        }
    }

    cout << "Num of queries: " << classifications.size() << endl;
    cout<<"Num of classifications: "<< counts.classificationCnt << endl;
    cout<<"Num of correct classifications: "<<counts.correct<<endl;
    cout<<"Num of correct but too broad classifications: "<<counts.highRank<<endl;
    cout << "classified/total = " << float(counts.classificationCnt)/float(classifications.size()) << endl;
    cout << "correct   /total = " << float(counts.correct) / float(classifications.size()) << endl;
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
    return 0;
}

void compareTaxon112(TaxID shot, TaxID target, NcbiTaxonomy & ncbiTaxonomy, Counts13 & counts) { ///target: subspecies or species
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
            counts.highRank ++;
            cout<<"U"<<endl;
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

