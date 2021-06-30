//
// Created by KJB on 10/09/2020.
//
//#include "createTargetDB.h"
#include "IndexCreator.h"
#include "FileMerger.h"
#include "LocalParameters.h"
#include <Command.h>
#include <regex>
#include "Classifier.h"
#include "omp.h"
#include <random>

void prepareForCreatingTargetDB(const LocalParameters & par, unordered_map<int, int> & speciesCnt);
void makeDiffIdxLookup(char * diffIdxFileName, char * infoFileName);

int createTargetDB(int argc, const char **argv, const Command &command)
{
    LocalParameters &par = LocalParameters::getLocalInstance();
    par.parseParameters(argc, argv, command, false, Parameters::PARSE_ALLOW_EMPTY, 0);
    const char * folder = par.filenames[0].c_str();
    const char * outputFileName = par.filenames[2].c_str();

    string genome_fname;
    string taxIdList_fname;
    string names, nodes, merged;

    ///Remove here
    unordered_map<int,int> subspeciesCnt; //(TaxID, 1)

    if(par.gtdbOrNcbi == 1 || par.gtdbOrNcbi == 0){
        cout<<"Creating target database based on taxonomy of GTDB"<<endl;
        prepareForCreatingTargetDB(par, subspeciesCnt);
        genome_fname = string(folder) + "/concatenated_genome_GTDB";
        taxIdList_fname = string(outputFileName) +"_taxID_list_GTDB";
        names = "../../gtdb_taxdmp/names.dmp";
        nodes = "../../gtdb_taxdmp/nodes.dmp";
        merged = "../../gtdb_taxdmp/merged.dmp";
    } else if(par.gtdbOrNcbi == 2){
        cout<<"Creating target database based on taxonomy of NCBI"<<endl;
        prepareForCreatingTargetDB(par, subspeciesCnt);
        genome_fname = string(folder) + "/concatenated_genome_NCBI";
        taxIdList_fname = string(outputFileName) +"_taxID_list_NCBI";
        names = "../../ncbi_taxdmp/names.dmp";
        nodes = "../../ncbi_taxdmp/nodes.dmp";
        merged = "../../ncbi_taxdmp/merged.dmp";
    } else{
        cout<<"ERROR"<<endl;
        return 0;
    }

    vector<int> speciesList;
//    for(auto it = subspeciesCnt.begin(); it != subspeciesCnt.end(); it ++){
//
//    }


    NcbiTaxonomy ncbiTaxonomy(names, nodes, merged);
    IndexCreator idxCre;
    cout<<"hi"<<endl;
    cout<<subspeciesCnt.begin()->first<< " "<< subspeciesCnt.begin()->second<<" "<<subspeciesCnt.size()<<endl;
    unordered_map<int,int> speciesTaxIdCnt; //<TAXid, cnt>
    for(auto it = subspeciesCnt.begin(); it != subspeciesCnt.end(); it++){
        cout<<it->first<<" "<<it->second<<" "<<ncbiTaxonomy.getTaxIdAtRank(it->first, "species")<<endl;
        speciesTaxIdCnt[ncbiTaxonomy.getTaxIdAtRank(it->first, "species")] ++;
    }
    cout<<"number of species: "<< speciesTaxIdCnt.size()<<endl;

    int max = 0;
    map<int,int> cntFre;
    for(auto it = speciesTaxIdCnt.begin(); it != speciesTaxIdCnt.end(); it++){
        if(it->second > max){
            max = it->second;
        }
        speciesList.push_back(it->first);
        cout<<max<<" "<<it->second<<endl;
        cntFre[it->second] ++;
    }

    cout<<"start"<<endl;
    for(auto it = cntFre.begin(); it != cntFre.end(); it ++){
        cout<<it->first<<'\t'<<it->second<<endl;
    }
    cout<<"Max "<< max<<endl;

    unordered_map<TaxID, int> familyCnt;
    for(auto it = speciesTaxIdCnt.begin(); it != speciesTaxIdCnt.end(); it++){
        familyCnt[ncbiTaxonomy.getTaxIdAtRank(it->first, "family")] = 0;
    }

    vector<TaxID> speciesToBeExcluded;
    srand(time(NULL));
    int cnt = 0;
    int randomN;
    while(cnt < familyCnt.size()){
        randomN = rand() % speciesList.size();
        if(familyCnt[ncbiTaxonomy.getTaxIdAtRank(speciesList[randomN], "family")] == 0){
            familyCnt[ncbiTaxonomy.getTaxIdAtRank(speciesList[randomN], "family")] = 1;
            speciesToBeExcluded.push_back(speciesList[randomN]);
            cnt++;
        }
    }

    for(auto x : speciesToBeExcluded){
        cout<<ncbiTaxonomy.taxonNode(x)->name<<endl;
    }

    int iii = 0;
    cout<<"Archaea"<<endl;
    for(auto it = subspeciesCnt.begin(); it != subspeciesCnt.end(); it++){
        if(speciesToBeExcluded.end() != find(speciesToBeExcluded.begin(), speciesToBeExcluded.end(), ncbiTaxonomy.getTaxIdAtRank(it->first, "species"))){
           // if(ncbiTaxonomy.getTaxIdAtRank(it->first, "superkingdom") == 2) {
                cout<< ncbiTaxonomy.taxonNode(it->first)->name << endl;
            //}
        }
    }

//    cout<<"Bacteria"<<endl;
//    for(auto it = subspeciesCnt.begin(); it != subspeciesCnt.end(); it++){
//        if(speciesToBeExcluded.end() != find(speciesToBeExcluded.begin(), speciesToBeExcluded.end(), ncbiTaxonomy.getTaxIdAtRank(it->first, "species"))){
//            if(ncbiTaxonomy.getTaxIdAtRank(it->first, "superkingdom") == 8034) {
//                cout << iii++ << " " << ncbiTaxonomy.taxonNode(it->first)->name << endl;
//            }
//        }
//    }

//    map<int,int> cntFre2;
//    for(auto it = familyCnt.begin(); it != familyCnt.end(); it++){
//        cntFre2[it->second]++;
//        if(it->second == 195)
//            cout<<"taxid "<<it->first<<endl;
//    }
//
//    cout<<"# species"<<'\t'<<"count"<<endl;
//    for(auto it = cntFre2.begin(); it != cntFre2.end(); it ++){
//        cout<<it->first<<'\t'<<it->second<<endl;
//    }








    return 0;

    const char * seqFileName = genome_fname.c_str();
    const char * taxIdFileName = taxIdList_fname.c_str();

    ifstream seqFile;
    seqFile.open(seqFileName);

    if (!seqFile.is_open()){
        cout<<"Cannot open the sequence file."<<endl;
        return 0;
    }
    seqFile.close();

    ///Make mapping from sequence ID to taxID. Index of vector is sequence ID.
    FILE * taxIdFile;
    if((taxIdFile = fopen(taxIdFileName,"r")) == NULL){
        cout<<"Cannot open the taxID list file."<<endl;
        return 0;
    }
    vector<int> taxIdList; char taxID[100];
    while(feof(taxIdFile) == 0) {
        fscanf(taxIdFile,"%s",taxID);
        taxIdList.push_back(atol(taxID));
    }
    fclose(taxIdFile);
    taxIdList.pop_back();
    vector<int> taxIdListAtSpecies;
    vector<int> taxIdListAtGenus;

    ncbiTaxonomy.createTaxIdListAtRank(taxIdList, taxIdListAtSpecies, "species");
    ncbiTaxonomy.createTaxIdListAtRank(taxIdList, taxIdListAtGenus, "genus");


    ///Make files of differential indexing and information of k-mers
    idxCre.startIndexCreatingParallel(seqFileName,outputFileName, taxIdListAtSpecies, taxIdList);

    int numOfSplits = idxCre.getNumOfFlush();
    char suffixedDiffIdxFileName[numOfSplits][100];
    char suffixedInfoFileName[numOfSplits][100];

    if(numOfSplits == 1){
        sprintf(suffixedDiffIdxFileName[0], "%s_diffIdx", outputFileName);
        sprintf(suffixedInfoFileName[0], "%s_info", outputFileName);
        cout<<"k-mer DB in: "<<endl;
        cout<<suffixedDiffIdxFileName[0]<<"and"<<endl;
        cout<<suffixedInfoFileName[0]<<endl;
        return 0;
    }

    ///Merge files
    vector<char *> diffSplits;
    vector<char *> infoSplits;
    for(int split = 0; split < numOfSplits ; split++){
        sprintf(suffixedDiffIdxFileName[split], "%s_%d_diffIdx", outputFileName, split);
        sprintf(suffixedInfoFileName[split], "%s_%d_info", outputFileName, split);
        diffSplits.push_back(suffixedDiffIdxFileName[split]);
        infoSplits.push_back(suffixedInfoFileName[split]);
    }

    char mergedDiffFileName[100];
    char mergedInfoFileName[100];
    char diffIdxSplitFileName[100];
    sprintf(mergedDiffFileName, "%s_diffIdx", outputFileName);
    sprintf(mergedInfoFileName, "%s_info", outputFileName);
    sprintf(diffIdxSplitFileName, "%s_split", outputFileName);
    FileMerger merger(mergedDiffFileName, mergedInfoFileName, diffIdxSplitFileName);
    merger.mergeTargetFiles(diffSplits, infoSplits,taxIdListAtSpecies, taxIdList);


    cout<<"k-mer DB in: "<<endl;
    cout<<mergedDiffFileName<<" and"<<endl;
    cout<<mergedInfoFileName<<endl;
    return 0;
}

void makeDiffIdxLookup(char * diffIdxFileName, char * infoFileName){

}

void prepareForCreatingTargetDB(const LocalParameters & par, unordered_map<int, int> & speciesCnt){
    const char * folder = par.filenames[0].c_str();
    const char * mappingFile = par.filenames[1].c_str();
    const char * outputFileName = par.filenames[2].c_str();

    int mode = par.gtdbOrNcbi;
    string taxid_fname_fname;
    string taxid_fname_sorted_fname;
    string fastList_fname;
    string taxidList_fname;
    string genome_fname;

    if(mode == 0 || mode == 1){
        taxid_fname_fname = string(folder) + "/taxid_filename_GTDB";
        taxid_fname_sorted_fname = string(folder) + "/taxid_filename_sorted_GTDB";
        fastList_fname = string(folder) + "/fasta_list_GTDB";
        taxidList_fname = string(outputFileName) + "_taxID_list_GTDB";
        genome_fname = string(folder) + "/concatenated_genome_GTDB";
        system("echo \"\t|\t\t|\" > ../../gtdb_taxdmp/merged.dmp");
    }else{
        taxid_fname_fname = string(folder) + "/taxid_filename_NCBI";
        taxid_fname_sorted_fname = string(folder) + "/taxid_filename_sorted_NCBI";
        fastList_fname = string(folder) + "/fasta_list_NCBI";
        taxidList_fname = string(outputFileName) + "_taxID_list_NCBI";
        genome_fname = string(folder) + "/concatenated_genome_NCBI";
    }

    system(("./../../util/unzip_and_list.sh "+ string(folder)+" "+fastList_fname).c_str());

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

    ifstream fastaList;
    ofstream taxID_fname;
    taxID_fname.open(taxid_fname_fname);
    fastaList.open(fastList_fname);
    string fileName;
    smatch assacc;
    int taxId;
    regex regex1("(GC[AF]_[0-9]*\\.[0-9]*)");
    if(fastaList.is_open()){
        cout<<"Writing taxID to fileName mapping file"<<endl;
        while(getline(fastaList,fileName,'\n')) {
            regex_search(fileName, assacc, regex1);
            if (assacc2taxid.count(assacc[0].str())) {
                taxId = assacc2taxid[assacc[0].str()];
                speciesCnt[taxId] = 1;
                taxID_fname << taxId << "\t" << fileName << endl;
            } else{
                cout<<assacc[0].str()<<" is excluded in creating target DB because it is not mapped to taxonomical ID"<<endl;
            }
        }
    }
    taxID_fname.close();

    return;
    system(("sort -k 1 -g "+taxid_fname_fname+" > "+taxid_fname_sorted_fname).c_str());
    system("chmod +x ./../../util/make_taxIdList_and_concatenatedGenome.sh");
    system(("./../../util/make_taxIdList_and_concatenatedGenome.sh "+taxidList_fname+" "+taxid_fname_sorted_fname+" "+genome_fname).c_str());
}

