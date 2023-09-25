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
    if (seqMode == 1 || seqMode == 3) {
        KSeqWrapper* kseq = KSeqFactory(queryPath_1.c_str());
        size_t kmerCnt = 0;
        size_t seqCnt = 0;
        size_t start = 0;
        while (kseq->ReadEntry()) {
            readNum_1++;
            seqCnt++;
            totalReadLength += kseq->entry.sequence.l;
            size_t currentKmerCnt = LocalUtil::getQueryKmerNumber<size_t>(kseq->entry.sequence.l, spaceNum);
            kmerCnt += currentKmerCnt;
            // std::cout << "currentKmerCnt: " << kmerCnt << "\n";
        
            if (bytesPerKmer * kmerCnt + ((size_t) 200 * seqCnt) > availableRam) {
                querySplits.emplace_back(start, readNum_1 - 1, kmerCnt - currentKmerCnt);
                kmerCnt = currentKmerCnt;
                start = readNum_1 - 1;
                seqCnt = 1;
            }
        }
        querySplits.emplace_back(start, readNum_1, kmerCnt);
        // Print elements
        for (auto & querySplit : querySplits) {
            std::cout << "start: " << querySplit.start << "\t";
            std::cout << "end: " << querySplit.end << "\t";
            std::cout << "kmerCnt: " << querySplit.kmerCnt << "\n";
        }
        delete kseq;
    } else {
        KSeqWrapper* kseq_1 = KSeqFactory(queryPath_1.c_str());
        KSeqWrapper* kseq_2 = KSeqFactory(queryPath_2.c_str());
        size_t kmerCnt = 0;
        size_t seqCnt_1 = 0;
        size_t seqCnt_2 = 0;
        size_t start = 0;
        size_t currentKmerCnt;
        bool end = false;
        while(true) {
            if (kseq_1->ReadEntry()) {
                readNum_1++;
                seqCnt_1++;
                totalReadLength += kseq_1->entry.sequence.l;
                currentKmerCnt = LocalUtil::getQueryKmerNumber<size_t>(kseq_1->entry.sequence.l, spaceNum);
                kmerCnt += currentKmerCnt;
            } else {
                end = true;
            }

            if (kseq_2->ReadEntry()) {
                readNum_2++;
                seqCnt_2++;
                totalReadLength += kseq_2->entry.sequence.l;
                currentKmerCnt += LocalUtil::getQueryKmerNumber<size_t>(kseq_2->entry.sequence.l, spaceNum);
                kmerCnt += currentKmerCnt;
            } else {
                end = true;
            }

            if (readNum_1 != readNum_2) {
                Debug(Debug::ERROR) << "The number of reads in the two files are not equal." << "\n";
                EXIT(EXIT_FAILURE);
            }

            if (bytesPerKmer * kmerCnt + ((size_t) 200 * seqCnt_1) > availableRam) {
                querySplits.emplace_back(start, readNum_1 - 1, kmerCnt - currentKmerCnt);
                kmerCnt = currentKmerCnt;
                start = readNum_1 - 1;
                seqCnt_1 = 1;
            }

            if (end) {
                querySplits.emplace_back(start, readNum_1, kmerCnt);
                break;
            }
        }
        delete kseq_1;
        delete kseq_2;
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