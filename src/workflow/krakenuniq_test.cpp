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
    //const string krakenTaxId = par.filenames[2];

    ///read classification
    string classString;
    ifstream readClassification;
    readClassification.open(readClassificationFileName);
    vector<int> classList;
    while(getline(readClassification,classString,'\n')){
        cout<<classString<<endl;
        classList.push_back(stoi(classString));
    }
    cout<<"num of classification: "<< classList.size()<<endl;
    for(int i = 0 ; i<classList.size(); i++){
        cout<<i<< " "<<classList[i]<<endl;
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