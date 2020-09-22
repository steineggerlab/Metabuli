//
// Created by KJB on 10/09/2020.
//
#include "createTargetDB.h"

vector<char*> createTargetDB()
{
    vector<char*> mergedFileNames;
    IndexCreator idxCre;
    string seqFileName;
    char kmerFileName[100];
    char taxIdFileName[300];
    FILE * taxIdFile;

    cout<<"Input the directory of a sequence file from which you want to make k-mer DB"<<endl;
    cin >> seqFileName;
    ifstream seqFile;
    seqFile.open(seqFileName);

    if (!seqFile.is_open())
    {
        cout<<"Cannot open the sequence file."<<endl;
        return mergedFileNames;
    }

    cout<<"Input the directory of corresponding taxID list"<<endl;
    cin>>taxIdFileName;

    if((taxIdFile = fopen(taxIdFileName,"r")) == NULL){
        cout<<"Cannot open the taxID list file."<<endl;
        return mergedFileNames;
    }
    vector<int> taxIdList; char taxID[100];
    while(feof(taxIdFile) == 0)
    {
        fscanf(taxIdFile,"%s",taxID);
        taxIdList.push_back(atol(taxID));
    }
    fclose(taxIdFile);


    cout<<"Input the directory where k-mer DB would be created"<<endl;
    cin>>kmerFileName;
    idxCre.startIndexCreating(seqFile,kmerFileName,taxIdList);


    int numOfSplits = idxCre.getNumOfFlush();
    char suffixedDiffIdxFileName[numOfSplits][100];
    char suffixedInfoFileName[numOfSplits][100];

    if(numOfSplits == 1)
    {
        sprintf(suffixedDiffIdxFileName[0],"%s_diffIdx_0", kmerFileName);
        sprintf(suffixedInfoFileName[0],"%s_info_0", kmerFileName);
        cout<<"k-mer DB in: "<<endl;
        cout<<suffixedDiffIdxFileName<<"and"<<endl;
        cout<<suffixedInfoFileName<<endl;
        return mergedFileNames;
    }
    vector<char *> diffSplits;
    vector<char *> infoSplits;

    for(int split = 0; split < numOfSplits ; split++)
    {
        sprintf(suffixedDiffIdxFileName[split],"%s_diffIdx_%d", kmerFileName,split);
        sprintf(suffixedInfoFileName[split],"%s_info_%d", kmerFileName, split);
        diffSplits.push_back(suffixedDiffIdxFileName[split]);
        infoSplits.push_back(suffixedInfoFileName[split]);
    }

    for(int i = 0; i < numOfSplits; i++)
    {
        cout<<diffSplits[i]<<endl;
    }
    char mergedDiffFileName[100];
    char mergedInfoFileName[100];
    sprintf(mergedDiffFileName,"%s_diffIdx", kmerFileName);
    sprintf(mergedInfoFileName,"%s_info", kmerFileName);
    DiffIdxMerger merger(mergedDiffFileName, mergedInfoFileName);
    merger.mergeTargetFiles(diffSplits, infoSplits);

    cout<<"k-mer DB in: "<<endl;
    cout<<mergedDiffFileName<<"and"<<endl;
    cout<<mergedInfoFileName<<endl;

    mergedFileNames.push_back(mergedDiffFileName);
    mergedFileNames.push_back(mergedInfoFileName);

    return mergedFileNames;
}


