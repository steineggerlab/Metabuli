//
// Created by KJB on 23/09/2020.
//

#include "Classifier.h"
#include "Parameters.h"
#include "LocalParameters.h"
#include "NcbiTaxonomy.h"

void setClassifyDefaults(LocalParameters & par){
    par.virusTaxId = 10239;// Taxonomy ID of virus taxon in NCBI
    par.seqMode = 1;
    par.memoryMode = 1;
    par.reducedAA = 0;
}
int classify(int argc, const char **argv, const Command& command)
{
    LocalParameters & par = LocalParameters::getLocalInstance();
    setClassifyDefaults(par);
    par.parseParameters(argc, argv, command, false, Parameters::PARSE_ALLOW_EMPTY, 0);
    omp_set_num_threads(par.threads);
    cout << "Number of threads: " << par.threads << endl;
    const char * queryFileName = par.filenames[0].c_str();
    const string databaseDirectory = par.filenames[1];
    const string taxonomyDirectory = par.filenames[2];

    const string names = taxonomyDirectory + "/names.dmp";
    const string nodes = taxonomyDirectory + "/nodes.dmp";
    const string merged = taxonomyDirectory + "/merged.dmp";
    NcbiTaxonomy taxonomy(names, nodes, merged);

    const string targetDiffIdxFileName = databaseDirectory+"/diffIdx";
    const string targetInfoFileName = databaseDirectory+"/info";
    string taxIdFileName = databaseDirectory+"/taxID_list";;
    const string diffIdxSplitFileName = databaseDirectory+"/split";;

    // Load the taxonomical ID list
    cout << "Loading taxonomy ID list ... ";
    FILE * taxIdFile;
    if((taxIdFile = fopen(taxIdFileName.c_str(),"r")) == NULL){
        cout<<"Cannot open the taxID list file."<<endl;
        return 0;
    }
    vector<int> taxIdList; char taxID[100];
    while(feof(taxIdFile) == 0)
    {
        fscanf(taxIdFile,"%s",taxID);
        taxIdList.push_back(atol(taxID));
    }
    fclose(taxIdFile);
    cout<<"Done"<<endl;

    Classifier classifier(par);
    classifier.startClassify(queryFileName, targetDiffIdxFileName.c_str(), targetInfoFileName.c_str(),
                             diffIdxSplitFileName.c_str(), taxIdList, par, taxonomy);

    return 0;
}