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
};

struct CountAtRank2 {
    int total;
    int FP;
    int TP;
    float precision;
    float sensitivity;
};

using namespace std;
void compareTaxon_hiv(TaxID shot, TaxID target, NcbiTaxonomy & ncbiTaxonomy, Counts & counts);
void compareTaxonAtRank2(TaxID shot, TaxID target, NcbiTaxonomy & ncbiTaxonomy, CountAtRank2 & count, const string & rank);


int inclusiontest_hiv(int argc, const char **argv, const Command &command){

    LocalParameters &par = LocalParameters::getLocalInstance();
    par.parseParameters(argc, argv, command, false, Parameters::PARSE_ALLOW_EMPTY, 0);

    const string readClassificationFileName = par.filenames[0];
    const char * mappingFile = par.filenames[1].c_str();
    string taxonomyPath = par.filenames[2];
    string names = taxonomyPath + "/names.dmp";
    string nodes = taxonomyPath + "/nodes.dmp";
    string merged = taxonomyPath + "/merged.dmp";
    NcbiTaxonomy ncbiTaxonomy(names, nodes, merged);
    NcbiTaxonomy * taxonomy = & ncbiTaxonomy;

    // Load the taxonomical ID list
    cout << "Loading taxonomy ID list ... ";
    FILE * taxIdFile;
    if((taxIdFile = fopen(mappingFile,"r")) == NULL){
        cout<<"Cannot open the taxID list file."<<endl;
        return 0;
    }
    vector<int> taxIdList; char taxID[100];
    while(feof(taxIdFile) == 0)
    {
        fscanf(taxIdFile,"%s",taxID);
        taxIdList.push_back(atol(taxID));
    }
    fclose(taxIdFile);

    vector<TaxID> genusTaxIdList;
    ncbiTaxonomy.createTaxIdListAtRank(taxIdList, genusTaxIdList, "genus");
    // Print TaxIDs in genusTaxIdList
    for(int i = 0; i < genusTaxIdList.size(); i++){
        cout << taxIdList[i] << "\t" << genusTaxIdList[i] << endl;
    }



    unordered_map<TaxID, unsigned int> taxonCnt;
    for(TaxID x : taxIdList){
        taxonCnt[x] = 1;
    }
    unordered_map<TaxID, TaxonCounts> cladeCnt = taxonomy->getCladeCounts(taxonCnt);
    vector<TaxID> genusWithMultipleSpecies;
    for(auto it = cladeCnt.begin(); it != cladeCnt.end(); it ++){
        if(taxonomy->taxonNode(it->first)->rank == "genus" && it->second.children.size() > 1){
            genusWithMultipleSpecies.push_back(it->first);
            cout << taxonomy->taxonNode(it->first)-> name << "\t";
            cout << it->second.children.size() << "\t";
            cout << taxonomy->taxonNode(it->second.children[0])->name << "\t"; //species 1
            if (cladeCnt[it->second.children[0]].children.size() != 0){
                cout << cladeCnt[it->second.children[0]].children[0] << "\t"; //subspecies 1
                cout << taxonomy->taxonNode(cladeCnt[it->second.children[0]].children[0])->name << "\t"; //subspecies
            } else {
                cout << it->second.children[0] << "\t"; //subspecies 1
                cout << taxonomy->taxonNode(it->second.children[0])->name << "\t"; //subspecies
            }

            if (cladeCnt[it->second.children[1]].children.size() != 0){
                cout << cladeCnt[it->second.children[1]].children[0] << "\t"; //subspecies 2
                cout << taxonomy->taxonNode(cladeCnt[it->second.children[1]].children[0])->name << "\t"; //subspecies
            } else {
                cout << it->second.children[1] << "\t"; //subspecies 2
                cout << taxonomy->taxonNode(it->second.children[1])->name << "\t"; //subspecies
            }
            cout << endl;
//            cout << taxonomy->taxonNode(cladeCnt[it->second.children[0]].children[0])->name << endl;
//            for(size_t i = 0; i < it->second.children.size(); i++){
//
//                for(size_t j = 0; j < cladeCnt[it->second.children[i]].children.size(); j++){
//                    cout << taxonomy->taxonNode(cladeCnt[it->second.children[i]].children[j])->name << "\t";
//                }
//                cout << endl;
//            }

//            for(size_t i = 0; i < 2; i++){
//                cout << taxonomy->taxonNode(it->second.children[0])->name << "\t" << taxonomy->taxonNode(it->second.children[i])->name << endl;
//            }
            // Iterate species
//            size_t selectedSp = 10000000;
//            for(size_t i = 0; i < it->second.children.size(); i++){
//                if(cladeCnt[it->second.children[i]].children.size() > 1){
//                    cout << taxonomy->taxonNode(it->second.children[i])->name << "\t" <<
//                            taxonomy->taxonNode(cladeCnt[it->second.children[i]].children[0])->name << "\t" <<
//                            taxonomy->taxonNode(cladeCnt[it->second.children[i]].children[1])->name << "\t";
//                            selectedSp = i;
//                    break;
//                }
//            }
//            if(selectedSp != 10000000) {
//                for (size_t i = 0; i < it->second.children.size(); i++) {
//                    if (i != selectedSp) {
//                        cout << taxonomy->taxonNode(cladeCnt[it->second.children[i]].children[0])->name << endl;
//                        break;
//                    }
//                }
//            }
        }
    }

    return 0;
    // 1) Load mapping file
    cout<<"Load mapping from accession ID to taxonomy ID ...";
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
    cout << "Done" <<endl;

    // Read classification
    vector<int> rightAnswers;
    vector<int> classList;

    string classString;
    ifstream readClassification;
    readClassification.open(readClassificationFileName);
    vector<string> fields;
    string field;
    string seqID;
    int classInt;

    smatch assacc;
    int pos;
    while(getline(readClassification,classString,'\n')){
        istringstream lineStream(classString);
        fields.clear();
        while(getline(lineStream, field, '\t')){
            fields.push_back(field);
        }

        classInt = stoi(fields[2]);
        seqID = fields[1];
        pos = seqID.find("_");
        seqID = seqID.substr(0,pos);
        classList.push_back(classInt);

        //regex_search(fields[1], assacc, regex1);
        rightAnswers.push_back(acc2taxid[seqID]);
    }
    cout<<"hi"<<endl;
    cout<<"num of classification: "<< classList.size()<<endl;

    Counts counts = {0,0,0,0,0,0,0,0};
    CountAtRank2 SS = {0, 0, 0, 0, 0};
    CountAtRank2 S = {0, 0, 0, 0, 0};
    CountAtRank2 G = {0, 0, 0, 0, 0};
    CountAtRank2 F = {0, 0, 0, 0, 0};
    ///score the classification
    for(size_t i = 0; i < classList.size(); i++){
        cout<<i<<" ";
        compareTaxon_hiv(classList[i], rightAnswers[i], ncbiTaxonomy, counts);
        compareTaxonAtRank2(classList[i], rightAnswers[i], ncbiTaxonomy, SS, "subspecies");
        compareTaxonAtRank2(classList[i], rightAnswers[i], ncbiTaxonomy, S, "species");
        compareTaxonAtRank2(classList[i], rightAnswers[i], ncbiTaxonomy, G, "genus");
        compareTaxonAtRank2(classList[i], rightAnswers[i], ncbiTaxonomy, F, "family");
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


    cout<<endl<<"NEW"<<endl;
    cout<<"Family      : " << F.total << " / " << F.TP << " / "<< F.FP << " / " << F.precision << " / "<< F.sensitivity << endl;
    cout<<"Genus       : " << G.total << " / " << G.TP << " / "<< G.FP << " / " << G.precision << " / "<< G.sensitivity << endl;
    cout<<"Species     : " << S.total << " / " << S.TP << " / "<< S.FP << " / " << S.precision << " / "<< S.sensitivity << endl;
    cout<<"Subspecies  : " << SS.total << " / " << SS.TP << " / "<< SS.FP << " / " << SS.precision << " / "<< SS.sensitivity << endl;

    return 0;
}

void compareTaxon_hiv(TaxID shot, TaxID target, NcbiTaxonomy & ncbiTaxonomy, Counts& counts) { ///target: subspecies or species

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

void compareTaxonAtRank2(TaxID shot, TaxID target, NcbiTaxonomy & ncbiTaxonomy, CountAtRank2 & count, const string & rank) {
    const TaxonNode * shotNode = ncbiTaxonomy.taxonNode(shot);
    const TaxonNode * targetNode = ncbiTaxonomy.taxonNode(target);
    TaxID shotTaxIdAtRank = shotNode->taxId;
    TaxID targetTaxIdAtRank = targetNode->taxId;
    if(shot == 1 || shot == 0) return;

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
