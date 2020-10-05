//
// Created by KJB on 23/09/2020.
//

#include "Classifier.h"
#include "Parameters.h"

using namespace std;

int classify(int argc, const char **argv, const Command& command)
{
    Classifier classifier;

    const char * queryFileName = argv[0];
    const char * targetDiffIdxFileName = argv[1];
    const char * targetInfoFileName = argv[2];
    const char * taxIdFileName = argv[3];
   // const char * outputdirectory = argv[4];

    FILE * taxIdFile;
    if((taxIdFile = fopen(taxIdFileName,"r")) == NULL){
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

    classifier.startClassify(queryFileName, targetDiffIdxFileName, targetInfoFileName, taxIdList);
    return 0;
}