#ifndef METABULI_LOCALUTIL_H
#define METABULI_LOCALUTIL_H

#include "Util.h"
#include <string>
#include "common.h"
#include "KSeqWrapper.h"
#include <cstdint>

class LocalUtil : public Util {
public:
    LocalUtil() = default;

    static std::string getQueryBaseName(const std::string & queryPath);

    template<typename T>
    static T getQueryKmerNumber(T queryLength, int spaceNum);

    template<typename T>
    static T getMaxCoveredLength(T queryLength);

    static int getFirstWhiteSpacePos(const std::string & str);

    // static std::string getAccessionFromHeader(const std::string & header);
};


template <typename T>
T LocalUtil::getQueryKmerNumber(T queryLength, int spaceNum) {
    return (getMaxCoveredLength(queryLength) / 3 - kmerLength - spaceNum + 1) * 6;
}

template<typename T>
T LocalUtil::getMaxCoveredLength(T queryLength) {
    if (queryLength % 3 == 2) {
        return queryLength - 2; 
    } else if (queryLength % 3 == 1) {
        return queryLength - 4; 
    } else {
        return queryLength - 3; 
    }
}

#endif //METABULI_LOCALUTIL_H
