#ifndef METABULI_LOCALUTIL_H
#define METABULI_LOCALUTIL_H

#include "Util.h"
#include <string>
#include "common.h"
#include "KSeqWrapper.h"

class LocalUtil : public Util {
public:
    LocalUtil() = default;

    static std::string getQueryBaseName(const std::string & queryPath);

    template<typename T>
    static T getQueryKmerNumber(T queryLength, int spaceNum);

    static void splitQueryFile(std::vector<SequenceBlock> & seqSegments, const std::string & queryPath);

};


template <typename T>
T LocalUtil::getQueryKmerNumber(T queryLength, int spaceNum) {
    return (getMaxCoveredLength(queryLength) / 3 - kmerLength - spaceNum + 1) * 6;
}


#endif //METABULI_LOCALUTIL_H
