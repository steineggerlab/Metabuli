//
// Created by KJB on 10/09/2020.
//
//#include "createTargetDB.h"
#include "IndexCreator.h"
#include "DiffIdxMerger.h"
#include "LocalParameters.h"
#include <Command.h>
int createTargetDB(int argc, const char **argv, const Command &command)
{
    LocalParameters &par = LocalParameters::getLocalInstance();
//    par.parseParameters(argc, argv, command, false, Parameters::PARSE_ALLOW_EMPTY, 0);

    vector<char*> mergedFileNames;
    IndexCreator idxCre;
    string name, node, merged;
//    cout<<"Input the directory for name.dmp"<<endl;
//    cin>>name;
//    cout<<"Input the directory for node.dmp"<<endl;
//    cin>>node;
//    cout<<"Input the directory for merged.dmp"<<endl;
//    cin>>merged;
    NcbiTaxonomy ncbiTaxonomy("/Users/jaebeomkim/Desktop/ADclassifier/taxdmp/names.dmp",
            "/Users/jaebeomkim/Desktop/ADclassifier/taxdmp/nodes.dmp",
                 "/Users/jaebeomkim/Desktop/ADclassifier/taxdmp/merged.dmp");
    const char * seqFileName = argv[0];
    const char * taxIdFileName = argv[1];
    const char * outputFileName = argv[2];

    ifstream seqFile;
    seqFile.open(seqFileName);

    if (!seqFile.is_open()){
        cout<<"Cannot open the sequence file."<<endl;
        return 0;
    }
    seqFile.close();

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

    vector<int> taxIdListAtRank;
    ncbiTaxonomy.makeTaxIdListAtRank(taxIdList, taxIdListAtRank, "species");

    idxCre.startIndexCreating2(seqFileName,outputFileName,taxIdListAtRank);

    int numOfSplits = idxCre.getNumOfFlush();
    char suffixedDiffIdxFileName[numOfSplits][100];
    char suffixedInfoFileName[numOfSplits][100];

    if(numOfSplits == 1)
    {
        sprintf(suffixedDiffIdxFileName[0], "%s_diffIdx_0", outputFileName);
        sprintf(suffixedInfoFileName[0], "%s_info_0", outputFileName);
        cout<<"k-mer DB in: "<<endl;
        cout<<suffixedDiffIdxFileName<<"and"<<endl;
        cout<<suffixedInfoFileName<<endl;
        return 0;
    }
    vector<char *> diffSplits;
    vector<char *> infoSplits;

    for(int split = 0; split < numOfSplits ; split++)
    {
        sprintf(suffixedDiffIdxFileName[split], "%s_diffIdx_%d", outputFileName, split);
        sprintf(suffixedInfoFileName[split], "%s_info_%d", outputFileName, split);
        diffSplits.push_back(suffixedDiffIdxFileName[split]);
        infoSplits.push_back(suffixedInfoFileName[split]);
    }

    for(int i = 0; i < numOfSplits; i++)
    {
        cout<<diffSplits[i]<<endl;
    }
    char mergedDiffFileName[100];
    char mergedInfoFileName[100];
    sprintf(mergedDiffFileName, "%s_diffIdx", outputFileName);
    sprintf(mergedInfoFileName, "%s_info", outputFileName);
    DiffIdxMerger merger(mergedDiffFileName, mergedInfoFileName);
    merger.mergeTargetFiles(diffSplits, infoSplits,taxIdListAtRank);

    cout<<"k-mer DB in: "<<endl;
    cout<<mergedDiffFileName<<" and"<<endl;
    cout<<mergedInfoFileName<<endl;

    mergedFileNames.push_back(mergedDiffFileName);
    mergedFileNames.push_back(mergedInfoFileName);

    return 0;
}


