#include "IndexCreator.h"
#include "FileMerger.h"
#include "LocalParameters.h"
#include <Command.h>
#include <regex>
#include <random>

void prepareForCreatingTargetDB(const LocalParameters & par);

int build_dir(int argc, const char **argv, const Command &command)
{
    LocalParameters &par = LocalParameters::getLocalInstance();
    par.parseParameters(argc, argv, command, false, Parameters::PARSE_ALLOW_EMPTY, 0);
    const char * folder = par.filenames[0].c_str();
    const string taxonomyDirectory = par.filenames[1];
    const char * dbDirectory = par.filenames[2].c_str();

    //Taxonomy
    const string names = taxonomyDirectory + "/names.dmp";
    const string nodes = taxonomyDirectory + "/nodes.dmp";
    const string merged = taxonomyDirectory + "/merged.dmp";
    NcbiTaxonomy ncbiTaxonomy(names, nodes, merged);

    prepareForCreatingTargetDB(par);
    string taxIdList_fname = string(dbDirectory) + "/taxID_list";
    string genome_fname = string(folder) + "/concatenated_genome";


    ifstream seqFile;
    seqFile.open(genome_fname.c_str());

    if (!seqFile.is_open()){
        cout<<"Cannot open the sequence file."<<endl;
        return 0;
    }
    seqFile.close();

    //Make mapping from sequence ID to taxID. Index of vector is sequence ID.
    FILE * taxIdFile;
    if((taxIdFile = fopen(taxIdList_fname.c_str(),"r")) == NULL){
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

    //Create lists of species taxonomical IDs of each sequences.
    cout<<"Create taxonomical ID list at species rank ... ";
    ncbiTaxonomy.createTaxIdListAtRank(taxIdList, taxIdListAtSpecies, "species");
    cout<<"done"<<endl;

    //Create lists of genus taxonomical IDs of each sequences.
    cout<<"Create taxonomical ID list at genus rank ... ";
    ncbiTaxonomy.createTaxIdListAtRank(taxIdList, taxIdListAtGenus, "genus");
    cout<<"done"<<endl;

    //Make files of differential indexing and information of k-mers
    cout<<"Start to creat reference DB file(s) ... ";
    IndexCreator idxCre;
    idxCre.startIndexCreatingParallel(genome_fname.c_str(), dbDirectory,
                                      taxIdListAtSpecies, taxIdList, par);
    cout<<"done"<<endl;

    //Merge files
    cout<<"Merge reference DB files ... "<<endl;
    int numOfSplits = idxCre.getNumOfFlush();
    char diffIdxSplitFileName[300];
    vector<char *> diffSplits;
    vector<char *> infoSplits;
    for(int split = 0; split < numOfSplits ; split++){
        char * suffixedDiffIdxFileName = new char[300];
        char * suffixedInfoFileName = new char[300];
        sprintf(suffixedDiffIdxFileName, "%s/%d_diffIdx", dbDirectory, split);
        sprintf(suffixedInfoFileName, "%s/%d_info", dbDirectory, split);
        diffSplits.push_back(suffixedDiffIdxFileName);
        infoSplits.push_back(suffixedInfoFileName);
    }

    char mergedDiffFileName[300];
    char mergedInfoFileName[300];
    sprintf(mergedDiffFileName, "%s/diffIdx", dbDirectory);
    sprintf(mergedInfoFileName, "%s/info", dbDirectory);
    sprintf(diffIdxSplitFileName, "%s/split", dbDirectory);

    FileMerger merger(mergedDiffFileName, mergedInfoFileName, diffIdxSplitFileName, par);
    merger.mergeTargetFiles(diffSplits, infoSplits,taxIdListAtSpecies, taxIdList);

    for(int split = 0; split < numOfSplits ; split++){
        delete[] diffSplits[split];
        delete[] infoSplits[split];
    }
    cout<<"done"<<endl;

    cout<<"Reference DB files you need are as below"<<endl;
    cout<<mergedDiffFileName<<endl;
    cout<<mergedInfoFileName<<endl;
    cout<<taxIdList_fname<<endl;
    cout<<diffIdxSplitFileName<<endl;

    return 0;
}

void prepareForCreatingTargetDB(const LocalParameters & par){
    const string folder = par.filenames[0].c_str();
    const string mappingFile = par.filenames[1] + "/assacc_to_taxid.tsv";
    const char * outputFileName = par.filenames[2].c_str();

    string taxid_fname_fname = folder + "/taxid_filename";
    string taxid_fname_sorted_fname = folder + "/taxid_filename_sorted";
    string fastaList_fname = folder + "/fasta_list_GTDB";
    string taxidList_fname = string(outputFileName) + "/taxID_list";
    string genome_fname = folder + "/concatenated_genome";

    system(("./../../util/unzip_and_list.sh " + folder + " " + fastaList_fname).c_str());

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
    fastaList.open(fastaList_fname);
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
                taxID_fname << taxId << "\t" << fileName << endl;
            } else{
                cout<<assacc[0].str()<<" is excluded in creating target DB because it is not mapped to taxonomical ID"<<endl;
            }
        }
    }
    taxID_fname.close();

    system(("sort -k 1 -g "+taxid_fname_fname+" > "+taxid_fname_sorted_fname).c_str());
    system("chmod +x ./../../util/make_taxIdList_and_concatenatedGenome.sh");
    system(("./../../util/make_taxIdList_and_concatenatedGenome.sh "+taxidList_fname+" "+taxid_fname_sorted_fname+" "+genome_fname).c_str());
}

