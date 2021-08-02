//
// Created by KJB on 23/09/2020.
//

#include "Classifier.h"
#include "Parameters.h"
#include "LocalParameters.h"

int classify(int argc, const char **argv, const Command& command)
{
    LocalParameters &par = LocalParameters::getLocalInstance();
    par.parseParameters(argc, argv, command, false, Parameters::PARSE_ALLOW_EMPTY, 0);

    Classifier classifier;
    const char * queryFileName = par.filenames[0].c_str();
    //const string dbName = par.filenames[1];
    //const string outputDirectory = par.filenames[2];
    const string databaseDirectory = par.filenames[1];

    const string targetDiffIdxFileName = databaseDirectory+"/diffIdx";
    const string targetInfoFileName = databaseDirectory+"/info";
    string taxIdFileName = databaseDirectory+"/taxID";;
    const string diffIdxSplitFileName = databaseDirectory+"/split";;

//    const string targetDiffIdxFileName = par.filenames[1];
//    const string targetInfoFileName = par.filenames[2];
//    string taxIdFileName = par.filenames[3];
//    const string diffIdxSplitFileName = par.filenames[4];



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

    classifier.startClassify(queryFileName, targetDiffIdxFileName.c_str(), targetInfoFileName.c_str(), diffIdxSplitFileName.c_str(), taxIdList, par);
   // classifier.startClassify(queryFileName, targetDiffIdxFileName.c_str(), targetInfoFileName.c_str(), diffIdxSplitFileName.c_str(), taxIdList, par);

    return 0;
}