#include "NcbiTaxonomy.h"
#include "Parameters.h"
#include "LocalParameters.h"

#include <string>
#include <sstream>
#include <fstream>
#include <iostream>
#include <regex>
#include "benchmark.h"

using namespace std;

int grade(int argc, const char **argv, const Command &command){

    LocalParameters &par = LocalParameters::getLocalInstance();
    par.parseParameters(argc, argv, command, false, Parameters::PARSE_ALLOW_EMPTY, 0);

    const string readClassificationFileName = par.filenames[0];
    const string mappingFile = par.filenames[1];
    const string taxonomy = par.filenames[2];

    string names = taxonomy + "/names.dmp";
    string nodes =  taxonomy + "/nodes.dmp";
    string merged =  taxonomy + "/merged.dmp";
    NcbiTaxonomy ncbiTaxonomy(names, nodes, merged);

    // Load the mapping file (answer sheet) (assacc to taxID)
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

    // Load classification results
    vector<int> rightAnswers;
    vector<int> classList;
    string classString;
    ifstream readClassification;
    readClassification.open(readClassificationFileName);
    vector<string> fields;
    string field;
    int classInt;
    vector<float> scores;
    vector<Score2> tpOrFp;
    regex regex1("(GC[AF]_[0-9]*\\.[0-9]*)");
    smatch assacc;
    while(getline(readClassification,classString,'\n')){
        istringstream lineStream(classString);
        fields.clear();
        while(getline(lineStream, field, '\t')){
            fields.push_back(field);
        }
        // Read ID -> right answer
        string id = fields[1];
        if (par.testType == "gtdb") {
            regex_search(fields[1], assacc, regex1);
            rightAnswers.push_back(assacc2taxid[assacc[0]]);
        } else if (par.testType == "hiv"){
            size_t pos = id.find('_');
            id = id.substr(0,pos);
            rightAnswers.push_back(assacc2taxid[id]);
        } else if (par.testType == "cami"){
            size_t pos = id.find('/');
            id = id.substr(0,pos);
            rightAnswers.push_back(assacc2taxid[id]);
        }

        // Read classification
        classInt = stoi(fields[2]);
        classList.push_back(classInt);
    }
    cout<<"The number of reads: "<< classList.size()<<endl;

    CountAtRank SS = {0, 0, 0, 0, 0};
    CountAtRank S = {0, 0, 0, 0, 0};
    CountAtRank G = {0, 0, 0, 0, 0};
    CountAtRank F = {0, 0, 0, 0, 0};

    // Score the classification
    for(size_t i = 0; i < classList.size(); i++){
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

    cout<<"Family      : " << F.total << " / " << F.TP << " / "<< F.FP << " / " << F.precision << " / "<< F.sensitivity << endl;
    cout<<"Genus       : " << G.total << " / " << G.TP << " / "<< G.FP << " / " << G.precision << " / "<< G.sensitivity << endl;
    cout<<"Species     : " << S.total << " / " << S.TP << " / "<< S.FP << " / " << S.precision << " / "<< S.sensitivity << endl;
    cout<<"Subspecies  : " << SS.total << " / " << SS.TP << " / "<< SS.FP << " / " << SS.precision << " / "<< SS.sensitivity << endl;

    return 0;
}

