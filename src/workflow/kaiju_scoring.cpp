//
// Created by 김재범 on 2021/7/12.
//

#include "Classifier.h"
#include "Parameters.h"
#include "LocalParameters.h"

struct Counts{
    int classificationCnt;
    int subspCnt;
    int spCnt;
    int genusCnt;
    int familyCnt;
    int orderCnt;
    int classCnt;
    int phylumCnt;
};

void compareTaxon(TaxID shot, TaxID target, NcbiTaxonomy & ncbiTaxonomy, Counts & counts);



int krakenuniq_test(int argc, const char **argv, const Command &command){

    LocalParameters &par = LocalParameters::getLocalInstance();
    par.parseParameters(argc, argv, command, false, Parameters::PARSE_ALLOW_EMPTY, 0);

    const string queryFileName = par.filenames[0];
    const string readClassificationFileName = par.filenames[1];
    const string krakenTaxDB = par.filenames[2];

    string names = "../../gtdb_taxdmp/names.dmp";
    string nodes = "../../gtdb_taxdmp/nodes.dmp";
    string merged = "../../gtdb_taxdmp/merged.dmp";
    NcbiTaxonomy ncbiTaxonomy(names, nodes, merged);
    for(int i = 0 ; i < 1000 ; i++) {
        cout<<ncbiTaxonomy.taxonNodes[i].id<<ncbiTaxonomy.taxonNodes[i].taxId<<ncbiTaxonomy.taxonNodes[i].rank<<endl;
    }

    ///Load taxDB of kraken
    unordered_map<int, int> child2parent;
    string childString, parentString, throwaway;
    int childInt, parentInt;
    ifstream taxDB;
    taxDB.open(krakenTaxDB);
    if(taxDB.is_open()){
        while(getline(taxDB,childString,'\t')){
            getline(taxDB, parentString, '\t');
            getline(taxDB, throwaway,'\n');
            childInt = stoi(childString);
            parentInt = stoi(parentString);
            if(childInt > 1000000000)
                child2parent[childInt] = parentInt;
        }
    } else{
        cout<<"Cannot open taxDB"<<endl;
    }
    taxDB.close();

    ///read classification
    string classString;
    ifstream readClassification;
    readClassification.open(readClassificationFileName);
    vector<int> classList;
    int classInt;
    while(getline(readClassification,classString,'\n')){
        classInt = stoi(classString);
        if(classInt > 1000000000){
            classList.push_back(child2parent[classInt]);
        } else{
            classList.push_back(classInt);
        }

    }
    cout<<"num of classification: "<< classList.size()<<endl;

    ///Load query file -> name
    regex regex1("(GC[AF]_[0-9]*\\.[0-9]*)");
    smatch assacc;
    string queryName;
    ifstream query;
    query.open(queryFileName);
    string queryLine;
    vector<string> queryNameList;
    while(getline(query,queryLine,'\n')){
        if(queryLine[0] == '>'){
            regex_search(queryLine, assacc, regex1);
            queryNameList.push_back(assacc[0]);
        }else{
            continue;
        }
    }

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

    ///right answer list
    vector<int> rightAnswers;
    for(size_t i = 0; i < queryNameList.size(); i++){
        if (assacc2taxid.count(queryNameList[i])) {
            rightAnswers.push_back(assacc2taxid[queryNameList[i]]);
        } else{
            cout << queryNameList[i] << " is not in the mapping file" << endl;
            rightAnswers.push_back(-1);
            continue;
        }
    }

    Counts counts = {0,0,0,0,0,0,0,0};
    ///score the classification
    for(size_t i = 0; i < queryNameList.size(); i++){
        counts.classificationCnt ++;
        compareTaxon(classList[i], rightAnswers[i], ncbiTaxonomy, counts);
    }

    cout<<"Number of classification: "<< counts.classificationCnt << endl;
    cout<<"classified / total =" << float(counts.classificationCnt)/float(queryNameList.size()) << endl;
    //cout<<"Superkingdom: "<< counts.superCnt <<endl;
    cout<<"Phylum: "<<counts.phylumCnt<<endl;
    cout<<"Class: "<<counts.classCnt<<endl;
    cout<<"Order: "<<counts.orderCnt<<endl;
    cout<<"Family: "<<counts.familyCnt<<endl;
    cout<<"Genus: "<< counts.genusCnt << endl;
    cout<<"Species: "<<counts.spCnt<<endl;
    cout<<"Subspecies: "<<counts.subspCnt<<endl;
    cout<<"(subS + S + G) / all classification" << float(counts.genusCnt + counts.spCnt + counts.subspCnt) / float(counts.classificationCnt) <<endl;
    cout<<"Num of queries: " << queryNameList.size() << endl;

}

void compareTaxon(TaxID shot, TaxID target, NcbiTaxonomy & ncbiTaxonomy, Counts& counts) { ///target: subspecies or species
    const TaxonNode * shotNode = ncbiTaxonomy.taxonNode(shot);
    string shotRank = shotNode->rank;
    cout<<shot<<" "<<target<<" "<<shotRank<<" ";
    if(shot == 0){
        cout<<"X"<<endl;
    }
    if(NcbiTaxonomy::findRankIndex(shotRank) <= 3){
        //cout<<"subspecies"<<endl;
        if(shot == target){
            cout<<"O"<<endl;
            counts.subspCnt ++;
        } else cout<<"X"<<endl;

    } else if(shotRank == "species") {
        //cout<<"species"<<endl;
        if(shot == ncbiTaxonomy.getTaxIdAtRank(target, "species")){
            counts.spCnt ++;
            cout<<"O"<<endl;
        }else cout<<"X"<<endl;
    } else if(shotRank == "genus"){
        //cout<<"genus"<<endl;
        if(shot == ncbiTaxonomy.getTaxIdAtRank(target, "genus")){
            counts.genusCnt ++;
            cout<<"O"<<endl;
        } else cout<<"X"<<endl;
    } else if(shotRank == "family"){
        //cout<<"family"<<endl;
        if(shot == ncbiTaxonomy.getTaxIdAtRank(target, "family")) {
            counts.familyCnt++;
        }
    }else if(shotRank == "order") {
        //cout<<"order"<<endl;
        if(shot == ncbiTaxonomy.getTaxIdAtRank(target, "order")) {
            counts.orderCnt++;
        }
    }else if(shotRank == "class") {
        //cout<<"class"<<endl;
        if(shot == ncbiTaxonomy.getTaxIdAtRank(target, "class")) {
            counts.classCnt++;
        }
    } else if(shotRank == "phylum") {
        //cout<<"phylum"<<endl;
        if(shot == ncbiTaxonomy.getTaxIdAtRank(target, "phylum")) {
            counts.phylumCnt++;
        }
    } else {
        return;
    }
}