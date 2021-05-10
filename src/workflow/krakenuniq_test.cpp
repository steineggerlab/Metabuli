//
// Created by 김재범 on 2021/05/10.
//

#include "Classifier.h"
#include "Parameters.h"
#include "LocalParameters.h"

int krakenuniq_test(int argc, const char **argv, const Command &command){

    LocalParameters &par = LocalParameters::getLocalInstance();
    par.parseParameters(argc, argv, command, false, Parameters::PARSE_ALLOW_EMPTY, 0);

    const string queryFileName = par.filenames[0];
    const string readClassificationFileName = par.filenames[1];
    const string krakenTaxDB = par.filenames[2];

    ///Load taxDB of kraken
    unordered_map<int, int> child2parent;
    string childString, parentString;
    int childInt, parentInt;
    ifstream taxDB;
    taxDB.open(krakenTaxDB);
    if(taxDB.is_open()){
        while(getline(taxDB,childString,'\t')){
            getline(taxDB, parentString, '\n');
            childInt = stoi(childString);
            parentInt = stoi(parentString);
            if(childInt > 1000000000)
                child2parent[childInt] = parentInt;
        }
    } else{
        cout<<"Cannot open file for mappig from assemlby accession to tax ID"<<endl;
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
    for(int i = 0 ; i<classList.size(); i++){
        cout<<i<< " "<<classList[i]<<endl;
    }

    ///Load query file
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

    for(int i = 0 ; i<queryNameList.size(); i++){
        cout<<i<< " "<<queryNameList[i]<<endl;
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

}