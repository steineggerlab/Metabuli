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

    static bool isFasta(const std::string & queryPath) {
        // Allow .fna, .fasta, .fna.gz, .fasta.gz, .fa, .fa.gz
        return Util::endsWith(".fna", queryPath) || Util::endsWith(".fasta", queryPath) ||
               Util::endsWith(".fna.gz", queryPath) || Util::endsWith(".fasta.gz", queryPath) ||
               Util::endsWith(".fa", queryPath) || Util::endsWith(".fa.gz", queryPath);
    }

    static bool isFastq(const std::string & queryPath) {
        // Allow .fq, .fastq, .fq.gz, .fastq.gz
        return Util::endsWith(".fq", queryPath) || Util::endsWith(".fastq", queryPath) ||
               Util::endsWith(".fq.gz", queryPath) || Util::endsWith(".fastq.gz", queryPath); 
    }

    static bool isValidQueryFile(const std::string & queryPath) {
        return isFasta(queryPath) || isFastq(queryPath);
    }

    template<typename T>
    static T getQueryKmerNumber(T queryLength, int spaceNum, bool onlyAA = false);

    template<typename T>
    static T getMaxCoveredLength(T queryLength);

    static int getFirstWhiteSpacePos(const std::string & str);

    // static std::string getAccessionFromHeader(const std::string & header);
};


template <typename T>
T LocalUtil::getQueryKmerNumber(T queryLength, int spaceNum, bool onlyAA) {
    if (onlyAA)
        return (getMaxCoveredLength(queryLength) / 3 - kmerLengthAA - spaceNum + 1) * 6;
    else
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
