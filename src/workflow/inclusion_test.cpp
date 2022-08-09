//
// Created by 김재범 on 2021/05/10.
//

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

int inclusiontest(int argc, const char **argv, const Command &command){

    LocalParameters &par = LocalParameters::getLocalInstance();
    par.parseParameters(argc, argv, command, false, Parameters::PARSE_ALLOW_EMPTY, 0);

    const string readClassificationFileName = par.filenames[0];
    string scoreFileName = readClassificationFileName + "_SC";

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
    const char * mappingFile = "../../gtdb_taxdmp/assacc_to_taxid.tsv";
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
        classInt = stoi(fields[2]);
        classList.push_back(classInt);
//        scores.push_back(stof(fields[4]));
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
    CountAtRank SS = {0, 0, 0, 0, 0};
    CountAtRank S = {0, 0, 0, 0, 0};
    CountAtRank G = {0, 0, 0, 0, 0};
    CountAtRank F = {0, 0, 0, 0, 0};
    ///score the classification
    for(size_t i = 0; i < classList.size(); i++){
        cout<<i<<" ";
        compareTaxon(classList[i], rightAnswers[i], ncbiTaxonomy, counts, tpOrFp, 3.0f);
        compareTaxonAtRank(classList[i], rightAnswers[i], ncbiTaxonomy, SS, "subspecies");
        compareTaxonAtRank(classList[i], rightAnswers[i], ncbiTaxonomy, S, "species");
        compareTaxonAtRank(classList[i], rightAnswers[i], ncbiTaxonomy, G, "genus");
        compareTaxonAtRank(classList[i], rightAnswers[i], ncbiTaxonomy, F, "family");
    }

    ofstream scoreFile;
    scoreFile.open(scoreFileName);
    for(size_t i = 0; i < tpOrFp.size(); i ++){
        scoreFile << tpOrFp[i].tf << "\t" << tpOrFp[i].rank << "\t" <<endl;
//<< tpOrFp[0].score << "\n";
    }
    scoreFile.close();

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
    cout<<"Family      : " << F.total << " / " << F.TP << " / "<< F.FP << " / " << F.precision << " / "<< F.sensitivity << endl;
    cout<<"Genus       : " << G.total << " / " << G.TP << " / "<< G.FP << " / " << G.precision << " / "<< G.sensitivity << endl;
    cout<<"Species     : " << S.total << " / " << S.TP << " / "<< S.FP << " / " << S.precision << " / "<< S.sensitivity << endl;
    cout<<"Subspecies  : " << SS.total << " / " << SS.TP << " / "<< SS.FP << " / " << SS.precision << " / "<< SS.sensitivity << endl;

}

