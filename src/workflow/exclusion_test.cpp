//
// Created by 김재범 on 2021/7/12.
//

#include "Classifier.h"
#include "Parameters.h"
#include "LocalParameters.h"

struct Counts{
    int classificationCnt;
    int correct;
    int highRank;

    //number of targets at each rank
    int speciesTargetNumber;
    int genusTargetNumber;
    int familyTargetNumber;
    int orderTargetNumber;
    int classTargetNumber;
    int phylumTargetNumber;
    int superkingdomTargetNumber;

    //number of correct classifications at each rank
    int speciesCnt;
    int genusCnt;
    int familyCnt;
    int orderCnt;
    int classCnt;
    int phylumCnt;
    int superkingdomCnt;
};

void compareTaxon2(TaxID shot, TaxID target, NcbiTaxonomy & ncbiTaxonomy, Counts & counts);



int exclusiontest(int argc, const char **argv, const Command &command){

    LocalParameters &par = LocalParameters::getLocalInstance();
    par.parseParameters(argc, argv, command, false, Parameters::PARSE_ALLOW_EMPTY, 0);

    const string queryFileName = par.filenames[0];
    const string readClassificationFileName = par.filenames[1];
    const string krakenTaxDB = par.filenames[2];

    string names = "../../gtdb_taxdmp/names.dmp";
    string nodes = "../../gtdb_taxdmp/nodes.dmp";
    string merged = "../../gtdb_taxdmp/merged.dmp";
    NcbiTaxonomy ncbiTaxonomy(names, nodes, merged);

    unordered_map<TaxID, unsigned int> taxCnt;
    for(int i = 0 ; i < ncbiTaxonomy.taxonNodes.size() ; i++) {
        if(ncbiTaxonomy.taxonNode(ncbiTaxonomy.taxonNodes[i].taxId)->rank == "subspecies") {
            taxCnt[ncbiTaxonomy.taxonNodes[i].taxId] = 1;
        }
    }

    unordered_map<TaxID, TaxonCounts> cladeCnt = ncbiTaxonomy.getCladeCounts(taxCnt);
    for(auto it = cladeCnt.begin(); it != cladeCnt.end(); it ++){
        if(ncbiTaxonomy.taxonNode(it->first)->rank == "genus") {
            cout << it->second.cladeCount << " " << it->second.taxCount << endl;
        }
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

    Counts counts = {0,0,0,0,0,0,0,0};
    ///right answer list
    vector<int> rightAnswers;
    int taxid;
    int taxid_sp;
    unsigned int cladeCnt_sp;
    for(size_t i = 0; i < queryNameList.size(); i++){
        if (assacc2taxid.count(queryNameList[i])) {
            taxid = assacc2taxid[queryNameList[i]];
            taxid_sp = ncbiTaxonomy.getTaxIdAtRank(taxid, "species");
            cladeCnt_sp = cladeCnt[taxid_sp].cladeCount;
            const TaxonNode * ancestor = ncbiTaxonomy.taxonNode(ncbiTaxonomy.getTaxIdAtRank(taxid_sp, "genus"));
            while(cladeCnt_sp != cladeCnt[ancestor->taxId].cladeCount){
                ancestor = ncbiTaxonomy.taxonNode(ancestor->parentTaxId);
                if(ancestor->rank == "superkingdom"){
                    break;
                }
            }
            rightAnswers.push_back(ancestor->taxId);
            if(ancestor->rank == "superkingdom"){
                counts.superkingdomTargetNumber ++;
            } else if(ancestor->rank == "phylum"){
                counts.phylumTargetNumber ++;
            } else if(ancestor->rank == "order"){
                counts.orderTargetNumber ++;
            } else if(ancestor->rank == "class"){
                counts.classTargetNumber ++;
            } else if(ancestor->rank == "family"){
                counts.familyTargetNumber ++;
            } else if(ancestor->rank == "genus"){
                counts.genusTargetNumber ++;
            }
        } else{
            cout << queryNameList[i] << " is not in the mapping file" << endl;
            rightAnswers.push_back(-1);
            continue;
        }
    }

    ///score the classification
    for(size_t i = 0; i < queryNameList.size(); i++){

        compareTaxon2(classList[i], rightAnswers[i], ncbiTaxonomy, counts);
    }

    cout<<"Num of queries: " << queryNameList.size() << endl;
    cout<<"Number of classification: "<< counts.classificationCnt << endl;
    cout<<"classified / total = " << float(counts.classificationCnt)/float(queryNameList.size()) << endl;
    cout<<"correct / total = "<< float(counts.correct) / float(queryNameList.size())<<endl;
    cout<<"correct / classification = "<<float(counts.correct) / float(counts.classificationCnt) <<endl<<endl;

    cout<<"Number of targets at each rank"<<endl;
    cout<<"Superkingdom: "<< counts.superkingdomTargetNumber<<endl;
    cout<<"Phylum: "<<counts.phylumTargetNumber<<endl;
    cout<<"Class: "<<counts.classTargetNumber<<endl;
    cout<<"Order: "<<counts.orderTargetNumber<<endl;
    cout<<"Family: "<<counts.familyTargetNumber<<endl;
    cout<<"Genus: "<< counts.genusTargetNumber<< endl;
    cout<<"Species: "<<counts.speciesTargetNumber<<endl<<endl;


}

void compareTaxon2(TaxID shot, TaxID target, NcbiTaxonomy & ncbiTaxonomy, Counts& counts) { ///target: subspecies or species
    const TaxonNode * shotNode = ncbiTaxonomy.taxonNode(shot);
    const TaxonNode * targetNode = ncbiTaxonomy.taxonNode(target);
    string shotRank = shotNode->rank;
    string targetRank = shotNode->rank;
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
    } else if(NcbiTaxonomy::findRankIndex(shotRank) <= NcbiTaxonomy::findRankIndex(shotRank)){ //classified into wrong taxon or too specifically
        cout<<"X"<<endl;
    } else { // classified at higher rank (too safe classification)
        if(shot == ncbiTaxonomy.getTaxIdAtRank(target, shotRank)){ //on right branch
            counts.highRank ++;
            cout<<"U"<<endl;
        } else{ //on wrong branch
            cout<<"X"<<endl;
        }
    }


    if(!isCorrect) return;

    if(shotRank == "species") {
        counts.speciesCnt ++;
    } else if(shotRank == "genus"){
        counts.genusCnt ++;
    } else if(shotRank == "family"){
        counts.familyCnt++;
    } else if(shotRank == "order") {
        counts.orderCnt++;
    } else if(shotRank == "class") {
        counts.classCnt++;
    } else if(shotRank == "phylum") {
        counts.phylumCnt++;
    } else if(shotRank == "superkingdom"){
        counts.superkingdomCnt++;
    }
}
