#include "QueryIndexer.h"

QueryIndexer::QueryIndexer(const LocalParameters & par) {
    seqMode = par.seqMode;
    if (seqMode == 1 || seqMode == 3) {
        queryPath_1 = par.filenames[0];
        queryPath_2 = "";
    } else {
        queryPath_1 = par.filenames[0];
        queryPath_2 = par.filenames[1];
    }

    matchPerKmer = par.matchPerKmer;
    maxRam = par.ramUsage;
    threads = par.threads;
    bytesPerKmer = sizeof(QueryKmer) + matchPerKmer * sizeof(Match);
    readNum_1 = 0;
    readNum_2 = 0;
    spaceNum = par.spaceMask.length() - kmerLength;
    totalReadLength = 0;

    setAvailableRam();
}

void QueryIndexer::setAvailableRam() {
    availableRam = ((size_t) maxRam * (size_t) 1024 * 1024 * 1024)
                         - ((size_t) 134217728 * (size_t) threads);
}

void QueryIndexer::indexQueryFile() {
    // Read 1
    KSeqWrapper* kseq;
    kseq = KSeqFactory(queryPath_1.c_str());
    size_t kmerCnt = 0;
    size_t seqCnt = 0;
    size_t start = 0;
    while (kseq->ReadEntry()) {
        readNum_1++;
        const KSeqWrapper::KSeqEntry &e = kseq->entry;
        totalReadLength += e.sequence.l;
        size_t currentKmerCnt = LocalUtil::getQueryKmerNumber<size_t>(e.sequence.l, spaceNum);
        kmerCnt += currentKmerCnt;
        seqCnt++;
        if (bytesPerKmer * kmerCnt + ((size_t) 200 * seqCnt) > availableRam) {
            querySplits.emplace_back(start, readNum_1, kmerCnt - currentKmerCnt);
            kmerCnt = currentKmerCnt;
            start = readNum_1;
            seqCnt = 1;
        }
    }
    querySplits.emplace_back(start, readNum_1, kmerCnt);
    delete kseq;

    // Read 2
    if (seqMode == 2) {
        kseq = KSeqFactory(queryPath_2.c_str());
        kmerCnt = 0;
        seqCnt = 0;
        start = 0;
        while (kseq->ReadEntry()) {
            readNum_2++;
            const KSeqWrapper::KSeqEntry &e = kseq->entry;
            totalReadLength += e.sequence.l;
            size_t currentKmerCnt = LocalUtil::getQueryKmerNumber<size_t>(e.sequence.l, spaceNum);
            kmerCnt += currentKmerCnt;
            seqCnt++;
            if (bytesPerKmer * kmerCnt + ((size_t) 200 * seqCnt) > availableRam) {
                querySplits.emplace_back(start, readNum_2, kmerCnt - currentKmerCnt);
                kmerCnt = currentKmerCnt;
                start = readNum_2;
                seqCnt = 1;
            }
        }
        querySplits.emplace_back(start, readNum_2, kmerCnt);
        delete kseq;

        // Check if the number of reads in the two files are equal
        if (readNum_1 != readNum_2) {
            Debug(Debug::ERROR) << "The number of reads in the two files are not equal." << "\n";
            EXIT(EXIT_FAILURE);
        }
    }
}

size_t QueryIndexer::getReadNum_1() const {
    return readNum_1;
}

size_t QueryIndexer::getReadNum_2() const {
    return readNum_2;
}

const std::vector<QuerySplit> & QueryIndexer::getQuerySplits() const {
    return querySplits;
}

std::size_t QueryIndexer::getTotalReadLength() const {
    return totalReadLength;
}

size_t QueryIndexer::getAvailableRam() const {
    return availableRam;
}