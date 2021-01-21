//
// Created by 김재범 on 2021/01/20.
//

#include "IndexCreator.h"
#include "FileMerger.h"
#include "LocalParameters.h"
#include <Command.h>

int updateTargetDB(int argc, const char **argv, const Command &command){

    LocalParameters &par = LocalParameters::getLocalInstance();

    NcbiTaxonomy ncbiTaxonomy("/Users/jaebeomkim/Desktop/pjt/taxdmp/names.dmp",
                              "/Users/jaebeomkim/Desktop/pjt/taxdmp/nodes.dmp",
                              "/Users/jaebeomkim/Desktop/pjt/taxdmp/merged.dmp");

    const char * seqFileName = argv[0];
    const char * taxIdFileName = argv[1];
    const char * outdatedFileName = argv[2];
    const char * outdatedTaxIdList = argv[3];///before _diffIdx or _info
    const char * updatedFileName = argv[4];

    ifstream seqFile;
    seqFile.open(seqFileName);

    if (!seqFile.is_open()){
        cout<<"Cannot open the sequence file."<<endl;
        return 0;
    }
    seqFile.close();

    ///Make mapping from sequence ID to TaxID. Index of vector is sequence ID.
    FILE * newTaxIdFile;
    FILE * oldTaxIdFile;
    if((newTaxIdFile = fopen(taxIdFileName, "r")) == NULL){
        cout<<"Cannot open the new taxID list file."<<endl;
        return 0;
    }

    if((oldTaxIdFile = fopen(outdatedTaxIdList, "r")) == NULL){
        cout<<"Cannot open the old taxID list file."<<endl;
        return 0;
    }

    vector<int> newTaxIdList; char taxID[100];
    vector<int> oldTaxIdList;
    while(feof(newTaxIdFile) == 0) {
        fscanf(newTaxIdFile, "%s", taxID);
        newTaxIdList.push_back(atol(taxID));
    }
    fclose(newTaxIdFile);

    while(feof(oldTaxIdFile) == 0) {
        fscanf(oldTaxIdFile, "%s", taxID);
        newTaxIdList.push_back(atol(taxID));
    }
    fclose(oldTaxIdFile);

    vector<int> newTaxIdListAtRank;
    vector<int> oldTaxIdListAtRank;
    ncbiTaxonomy.makeTaxIdListAtRank(newTaxIdList, newTaxIdListAtRank, "species");
    ncbiTaxonomy.makeTaxIdListAtRank(oldTaxIdList, oldTaxIdListAtRank, "species");

    IndexCreator idxCre;
    idxCre.startIndexCreatingParallel(seqFileName, updatedFileName, newTaxIdListAtRank, newTaxIdList);


    /**Merge new k-mer data into outdated database.**/

    ///Add outdated file names to the list of files to be merged
    vector<char *> diffSplits; //list of diff. index files
    vector<char *> infoSplits; //list of info. files
    char suffixedOutdatedDiffIdx[300];
    char suffixedOutdatedInfo[300];
    sprintf(suffixedOutdatedDiffIdx, "%s_diffIdx", outdatedFileName);
    sprintf(suffixedOutdatedInfo, "%s_info", outdatedFileName);
    diffSplits.push_back(suffixedOutdatedDiffIdx);
    infoSplits.push_back(suffixedOutdatedInfo);

    ///Add newly made files to the list
    int numOfSplits = idxCre.getNumOfFlush();
    char suffixedDiffIdxFileName[numOfSplits][100];
    char suffixedInfoFileName[numOfSplits][100];
    for(int split = 0; split < numOfSplits ; split++)
    {
        sprintf(suffixedDiffIdxFileName[split], "%s_diffIdx_%d", updatedFileName, split);
        sprintf(suffixedInfoFileName[split], "%s_info_%d", updatedFileName, split);
        diffSplits.push_back(suffixedDiffIdxFileName[split]);
        infoSplits.push_back(suffixedInfoFileName[split]);
    }

    ///Merge the files in the list
    char mergedDiffFileName[300];
    char mergedInfoFileName[300];
    sprintf(mergedDiffFileName, "%s_diffIdx", updatedFileName);
    sprintf(mergedInfoFileName, "%s_info", updatedFileName);
    FileMerger merger(mergedDiffFileName, mergedInfoFileName);
    merger.updateTargetDatabase(diffSplits, infoSplits, oldTaxIdListAtRank, oldTaxIdList, newTaxIdListAtRank, newTaxIdList, newTaxIdList.size()); ///이거 고쳐야함, 둘 다 필

    cout<<"k-mer DB in: "<<endl;
    cout<<mergedDiffFileName<<" and"<<endl;
    cout<<mergedInfoFileName<<endl;

    ///Make an updated list of taxonomical IDs
    char mergedTaxIdListName[300];
    sprintf(mergedTaxIdListName, "%s_taxIDs", updatedFileName);

    FILE * newTaxIdListFile;
    FILE * oldTaxIdList;
    FILE * mergedTaxIdList;
    newTaxIdListFile = fopen(taxIdFileName, "r");
    oldTaxIdList = fopen(outdatedTaxIdList, "r");
    mergedTaxIdList = fopen(mergedDiffFileName, "w");

    ///알아서 복사 붙여 넣기 해라 귀찮다.

    fclose(newTaxIdListFile);
    fclose(oldTaxIdList);
    fclose(mergedTaxIdList);

    return 0;
}