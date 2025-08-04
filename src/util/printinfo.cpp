#include "LocalParameters.h"
#include <Command.h>
#include <cstddef>
#include <string>
#include <iostream>
#include "KmerMatcher.h"
#include "common.h"
#include "SeqIterator.h"

using namespace std;

int printInfo(int argc, const char **argv, const Command &command){
    LocalParameters &par = LocalParameters::getLocalInstance();
    par.parseParameters(argc, argv, command, false, Parameters::PARSE_ALLOW_EMPTY, 0);
    string infoFileName = par.filenames[0];
    
    MmapedData<TaxID> infoFile = mmapData<TaxID>(infoFileName.c_str());
    size_t infoNum = infoFile.fileSize / sizeof(TaxID);
   
    size_t begin = par.infoBegin;
    size_t end = par.infoEnd;
    if (end > infoNum) {
        end = infoNum;
    }
    for (size_t i = begin; i < end; i++) {
        if (infoFile.data[i] == 0) {
            cout << "No taxID for index " << i << "\n";
            continue;
        }

        if (infoFile.data[i] == 372305) {
            cout << "372305 " << i << "\n";
            continue;
        }

        // cout << infoFile.data[i] << "\n";
    }
    munmap(infoFile.data, infoFile.fileSize + 1);

    return 0;
}