//
// Created by 김재범 on 2022/07/07.
//
#include "common.h"

void loadTaxIdList(const char * fileName, std::vector<TaxID> & taxids){
    FILE * taxIdFile;
    if((taxIdFile = fopen(fileName,"r")) == NULL){
        std::cout<<"Cannot open the taxID list file."<<std::endl;
        return;
    }
    char taxID[100];
    while(feof(taxIdFile) == 0) {
        fscanf(taxIdFile,"%s",taxID);
        taxids.push_back(atol(taxID));
    }
    fclose(taxIdFile);
    taxids.pop_back();
}