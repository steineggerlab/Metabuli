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
    cout<<"Loading Taxonomy"<<endl;
    const string names = taxonomyDirectory + "/names.dmp";
    const string nodes = taxonomyDirectory + "/nodes.dmp";
    const string merged = taxonomyDirectory + "/merged.dmp";
    NcbiTaxonomy ncbiTaxonomy(names, nodes, merged);

    // Make a tax ID list using mapping file (acc2taxID)
    // 1) Load mapping file
    cout<<"Load mapping from accession ID to taxonomy ID"<<endl;
    unordered_map<string, int> acc2taxid;
    string eachLine;
    string eachItem;
    ifstream map;
    map.open(acc2taxidFile);
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

    // 2) Make a tax ID list
    cout<<"<Make a taxonomy ID list"<<endl;
    ifstream seqFile;
    seqFile.open(fastaName);
    string accessionID;
    string accessionID2;
    vector<TaxID> taxIDs;
    if (seqFile.is_open()){
        while(getline(seqFile, eachLine, '\n')){
            if(eachLine[0] == '>'){
                istringstream ss(eachLine);
                getline(ss, accessionID, ' ');
                accessionID2 = accessionID.substr(1);
                if(acc2taxid.find(accessionID2) != acc2taxid.end()){
                    taxIDs.push_back(acc2taxid[accessionID2]);
                } else {
                    cout<<accessionID2<<" is not in the mapping file"<<endl;
                    taxIDs.push_back(0);
                }
            }
        }
    } else{
        cout<<"Cannot open the FASTA file."<<endl;
        return 0;
    }
    seqFile.close();

    // 3) Write the list into a file
    cout<<"Write the taxonomy list into a file"<<endl;
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
    taxIdFile.close();

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
    idxCre.startIndexCreatingParallel(fastaName, dbDirectory, taxIdListAtSpecies, taxIDs, par);
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
    char parametersFileName[300];
    sprintf(mergedDiffFileName, "%s/diffIdx", dbDirectory);
    sprintf(mergedInfoFileName, "%s/info", dbDirectory);
    sprintf(parametersFileName, "%s/parameters.txt", dbDirectory);
    sprintf(diffIdxSplitFileName, "%s/split", dbDirectory);

    FileMerger merger(par);
    merger.mergeTargetFiles(diffSplits, infoSplits,taxIdListAtSpecies, taxIDs);

    for(int split = 0; split < numOfSplits ; split++){
        delete[] diffSplits[split];
        delete[] infoSplits[split];
    }
    cout<<"done"<<endl;

    // Write parameters used
    ofstream params;
    params.open(parametersFileName);
    params.write(("Mask for spaced k-mer: " + par.spaceMask).c_str(), (int)("Mask for spaced k-mer: " + par.spaceMask).length());
    params.write(string("Number of alphabets for encoding amino acids: " + to_string(par.reducedAA)).c_str(),
                 (int) ("Number of alphabets for encoding amino acids: " + to_string(par.reducedAA)).length());
    cout<<"Reference DB files you need are as below"<<endl;
    cout<<mergedDiffFileName<<endl;
    cout<<mergedInfoFileName<<endl;
    cout<<taxIdFileName<<endl;
    cout<<diffIdxSplitFileName<<endl;

    return 0;
}
