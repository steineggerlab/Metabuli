#include "LocalParameters.h"
#include <Command.h>
#include <cstddef>
#include <string>
#include <iostream>
#include "KmerMatcher.h"
#include "common.h"
#include "InfoIndex.h"
#include "SeqIterator.h"

using namespace std;

void setPrintInfoDefault(LocalParameters &par) {
    par.infoBegin = 0;
    par.infoEnd = 0;
}

int printInfo(int argc, const char **argv, const Command &command){
    LocalParameters &par = LocalParameters::getLocalInstance();
    par.parseParameters(argc, argv, command, false, Parameters::PARSE_ALLOW_EMPTY, 0);
    string infoFileName = par.filenames[0];
    
    InfoIndexReader infoIdxFile(infoFileName, 1024 * 1024 * 16);
    size_t maxIdx = infoIdxFile.getTotalValueNum();
    size_t begin = par.infoBegin;
    size_t end = par.infoEnd;
    if (end > maxIdx) {
        end = maxIdx;
    }
    if (begin > end) {
        begin = end;
    }

    infoIdxFile.loadAt(begin);
    uint32_t value = 0;
    for (size_t idx = begin; idx < end && infoIdxFile.next(value); ++idx) {
        std::cout << value << std::endl;
    }
    return 0;
}
