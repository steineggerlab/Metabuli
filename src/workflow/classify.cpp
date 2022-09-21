#include "Classifier.h"
#include "ReducedClassifier.h"
#include "Parameters.h"
#include "LocalParameters.h"
#include "NcbiTaxonomy.h"

void setClassifyDefaults(LocalParameters & par){
    par.virusTaxId = 10239;// Taxonomy ID of virus taxon in NCBI
    par.seqMode = 2;
    par.memoryMode = 1;
    par.reducedAA = 0;
    par.minScore = 0.1;
    par.spaceMask = "11111111";
    par.minConsCnt = 4;
    par.hammingMargin = 0;
    par.minSpScore = 0.7;
}

int classify(int argc, const char **argv, const Command& command)
{
    LocalParameters & par = LocalParameters::getLocalInstance();
    setClassifyDefaults(par);
    par.parseParameters(argc, argv, command, true, Parameters::PARSE_ALLOW_EMPTY, 0);

#ifdef OPENMP
    omp_set_num_threads(par.threads);
#endif

    cout << "Number of threads: " << par.threads << endl;
    const char * queryFileName = par.filenames[0].c_str();
    const string databaseDirectory = par.filenames[1];

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

    Classifier * classifier;
    if(par.reducedAA == 1){
        classifier = new ReducedClassifier(par, taxIdList);
    } else {
        classifier = new Classifier(par, taxIdList);
    }

    classifier->startClassify(targetDiffIdxFileName.c_str(), targetInfoFileName.c_str(),
                             diffIdxSplitFileName.c_str(), par);
    delete classifier;
    return 0;
}