#include "IndexCreator.h"
#include "FileMerger.h"
#include "LocalParameters.h"
#include <Command.h>
#include <sstream>


int build_fasta(int argc, const char **argv, const Command &command)
{
    LocalParameters &par = LocalParameters::getLocalInstance();
    par.parseParameters(argc, argv, command, false, Parameters::PARSE_ALLOW_EMPTY, 0);
    const char * fastaName = par.filenames[0].c_str();
    const char * acc2taxidFile = par.filenames[1].c_str();
    const char * dbDirectory = par.filenames[2].c_str();
    const string taxonomyDirectory = par.filenames[3];

    string genome_fname;
    string taxIdList_fname;

    // Taxonomy
    const string names = taxonomyDirectory + "/names.dmp";
    const string nodes = taxonomyDirectory + "/nodes.dmp";
    const string merged = taxonomyDirectory + "/merged.dmp";
    NcbiTaxonomy ncbiTaxonomy(names, nodes, merged);

    // Make a tax ID list using mapping file (acc2taxID)
    // 1) Load mapping file
    unordered_map<string, int> acc2taxid;
    string key, value;
    ifstream map;
    map.open(acc2taxidFile);
    if(map.is_open()){
        while(getline(map,key,'\t')){
            getline(map, value, '\n');
            acc2taxid[key] = stoi(value);
        }
    } else{
        cout<<"Cannot open file for mapping from accession to tax ID"<<endl;
    }
    map.close();

    // 2) Make a tax ID list
    ifstream seqFile;
    seqFile.open(fastaName);
    string eachLine;
    string accessionID;
    vector<TaxID> taxIDs;
    if (seqFile.is_open()){
        while(getline(seqFile, eachLine, '\n')){
            if(eachLine[0] == '>'){
                istringstream ss(eachLine);
                getline(ss, accessionID, ' ');
                taxIDs.push_back(acc2taxid[accessionID]);
            }
        }
    } else{
        cout<<"Cannot open the FASTA file."<<endl;
        return 0;
    }
    seqFile.close();

    // 3) Write the list into a file
    const string taxIdFileName = string(dbDirectory) + "/taxID_list";
    ofstream taxIdFile;
    taxIdFile.open(taxIdFileName);
    if(taxIdFile.is_open()){
        for(int taxID : taxIDs){
            taxIdFile << taxID << '\n';
        }
    } else{
        cout<<"Cannot open a file for writing taxonomy ID list."<<endl;
    }

    //Create lists of species taxonomical IDs of each sequences.
    cout<<"Create taxonomical ID list at species rank ... ";
    vector<int> taxIdListAtSpecies;
    ncbiTaxonomy.createTaxIdListAtRank(taxIDs, taxIdListAtSpecies, "species");
    cout<<"done"<<endl;

    //Create lists of genus taxonomical IDs of each sequences.
    cout<<"Create taxonomical ID list at genus rank ... ";
    vector<int> taxIdListAtGenus;
    ncbiTaxonomy.createTaxIdListAtRank(taxIDs, taxIdListAtGenus, "genus");
    cout<<"done"<<endl;

    //Make files of differential indexing and information of k-mers
    cout<<"Start to creat reference DB file(s) ... ";
    IndexCreator idxCre;
    idxCre.startIndexCreatingParallel(fastaName, dbDirectory, taxIdListAtSpecies, taxIDs);
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

    FileMerger merger(mergedDiffFileName, mergedInfoFileName, diffIdxSplitFileName);
    merger.mergeTargetFiles(diffSplits, infoSplits,taxIdListAtSpecies, taxIDs);

    for(int split = 0; split < numOfSplits ; split++){
        delete[] diffSplits[split];
        delete[] infoSplits[split];
    }
    cout<<"done"<<endl;

    cout<<"Reference DB files you need are as below"<<endl;
    cout<<mergedDiffFileName<<endl;
    cout<<mergedInfoFileName<<endl;
    cout<<taxIdFileName<<endl;
    cout<<diffIdxSplitFileName<<endl;

    return 0;
}
