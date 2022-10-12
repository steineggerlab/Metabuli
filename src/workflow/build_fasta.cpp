#include "IndexCreator.h"
#include "FileMerger.h"
#include "LocalParameters.h"
#include <Command.h>
#include <sstream>

void setDefaults_build_fasta(LocalParameters & par){
    par.reducedAA = 0;
    par.spaceMask = "11111111";
}


int build_fasta(int argc, const char **argv, const Command &command) {
    LocalParameters &par = LocalParameters::getLocalInstance();
    setDefaults_build_fasta(par);
    par.parseParameters(argc, argv, command, true, Parameters::PARSE_ALLOW_EMPTY, 0);
    const char *fastaName = par.filenames[0].c_str();
    const char *acc2taxidFile = par.filenames[2].c_str();
    const string dbDirectory = par.filenames[1];
    const string taxonomyDirectory = dbDirectory + "/taxonomy";

    string genome_fname;
    string taxIdList_fname;

    // Taxonomy
    cout << "Loading Taxonomy" << endl;
    const string names = taxonomyDirectory + "/names.dmp";
    const string nodes = taxonomyDirectory + "/nodes.dmp";
    const string merged = taxonomyDirectory + "/merged.dmp";
    NcbiTaxonomy ncbiTaxonomy(names, nodes, merged);

//    // Make a tax ID list using mapping file (acc2taxID)
//    // 1) Load mapping file
//    cout << "Load mapping from accession ID to taxonomy ID" << endl;
//    unordered_map<string, int> acc2taxid;
//    cout << acc2taxid.max_size() << endl; //2^32
//    string eachLine;
//    string eachItem;
//    if (FILE * mappingFile = fopen(acc2taxidFile, "r")) {
//        char buffer[512];
//        int taxID;
//        fscanf(mappingFile, "%*s\t%*s\t%*s\t%*s");
//        while (fscanf(mappingFile, "%*s\t%s\t%d\t%*d", buffer, &taxID) == 2 ){
//            acc2taxid[string(buffer)] = taxID;
//        }
//    } else {
//        cerr << "Cannot open file for mapping from accession to tax ID" << endl;
//    }
//
//    // 2) Make a tax ID list
//    cout << "Make a taxonomy ID list" << endl;
//    ifstream seqFile;
//    seqFile.open(fastaName);
//    string accessionID;
//    string accessionID2;
//    vector<TaxID> taxIDs;
//    if (seqFile.is_open()) {
//        while (getline(seqFile, eachLine, '\n')) {
//            if (eachLine[0] == '>') {
//                istringstream ss(eachLine);
//                getline(ss, accessionID, ' ');
//                accessionID2 = accessionID.substr(1);
//                if (acc2taxid.find(accessionID2) != acc2taxid.end()) {
//                    taxIDs.push_back(acc2taxid[accessionID2]);
//                } else {
//                    cout << accessionID2 << " is not in the mapping file" << endl;
//                    taxIDs.push_back(0);
//                }
//            }
//        }
//    } else {
//        cerr << "Cannot open the FASTA file." << endl;
//        return 0;
//    }
//    seqFile.close();
//
//    // 3) Write the list into a file
//    cout << "Write the taxonomy list into a file ... ";
//    const string taxIdFileName = string(dbDirectory) + "/taxID_list";
//    ofstream taxIdFile;
//    taxIdFile.open(taxIdFileName);
//    if (taxIdFile.is_open()) {
//        for (int taxID: taxIDs) {
//            taxIdFile << taxID << '\n';
//        }
//    } else {
//        cerr << "Cannot open a file for writing taxonomy ID list." << endl;
//        return 0;
//    }
//    taxIdFile.close();
//    cout << "Done" << endl;

    // Load the taxonomical ID list
    cout << "Loading taxonomy ID list ... ";
    FILE * taxIdFile;
    if((taxIdFile = fopen(string(dbDirectory + "/taxID_list").c_str(),"r")) == NULL){
        cout<<"Cannot open the taxID list file."<<endl;
        return 0;
    }
    vector<int> taxIDs; char taxID[100];
    while(feof(taxIdFile) == 0)
    {
        fscanf(taxIdFile,"%s",taxID);
        taxIDs.push_back(atol(taxID));
    }
    taxIDs.pop_back();
    fclose(taxIdFile);
    cout<<"Done"<<endl;

    //Create lists of species taxonomical IDs of each sequences.
    vector<int> taxIdListAtSpecies;
    vector<int> taxIdListAtSuperkingdom;
    ncbiTaxonomy.createTaxIdListAtRank(taxIDs, taxIdListAtSpecies, "species");
    ncbiTaxonomy.getSuperKingdoms(taxIdListAtSpecies, taxIdListAtSuperkingdom);

    for(int i = 0; i < taxIdListAtSpecies.size(); i++){
        cout << taxIdListAtSpecies[i] << " " << taxIdListAtSuperkingdom[i] << endl;
    }

    return 0;
    
    //Make files of differential indexing and information of k-mers
    cout << "Start to create reference DB file(s) ... " << endl;
    IndexCreator idxCre(par);
    idxCre.startIndexCreatingParallel(fastaName, dbDirectory.c_str(), taxIdListAtSuperkingdom,
                                      taxIdListAtSpecies, taxIDs, par);
    cout << "done" << endl;

    //Merge files
    cout << "Merge reference DB files ... " << endl;
    FileMerger merger(par);
    merger.mergeTargetFiles(par, idxCre.getNumOfFlush());

    // Write parameters used
    ofstream params;
    params.open(string(dbDirectory) + "/parameters");
    params.write(("Mask for spaced k-mer: " + par.spaceMask).c_str(),
                 (int) ("Mask for spaced k-mer: " + par.spaceMask).length());
    params.write(string("Number of alphabets for encoding amino acids: " + to_string(par.reducedAA)).c_str(),
                 (int) ("Number of alphabets for encoding amino acids: " + to_string(par.reducedAA)).length());
    return 0;
}
