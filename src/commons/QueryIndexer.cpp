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

    // matchPerKmer = par.matchPerKmer;
    maxRam = par.ramUsage;
    threads = par.threads;
    // std::cout << "bytesPerKmer: " << bytesPerKmer << "\n";
    spaceNum = 0; // par.spaceMask.length() - kmerLength;
    totalReadLength = 0;
    kmerLen = 8;

    setAvailableRam();
}

void QueryIndexer::setAvailableRam() {
    availableRam = ((size_t) maxRam * (size_t) 1024 * 1024 * 1024)
                         - ((size_t) 128 * 1024 * 1024 * (size_t) threads);
    // std::cout << "availableRam: " << availableRam << "\n";
}

void QueryIndexer::indexQueryFile(size_t processedQueryNum) {
    querySplits.clear();
    readNum_1 = 0;
    readNum_2 = 0;
    // Read 1
    if (seqMode == 1 || seqMode == 3) {
        KSeqWrapper* kseq = KSeqFactory(queryPath_1.c_str());
        size_t seqPos = 0;
        
        // Skip processed reads
        for (size_t i = 0; i < processedQueryNum; i++) {
            kseq->ReadEntry();
            seqPos++;
        }
        size_t start = 0;
        size_t kmerCnt = 0;
        size_t seqCnt = 0;
        while (kseq->ReadEntry()) {
            seqPos++;
            if (unlikely(kseq->entry.sequence.l == 0 || kseq->entry.name.l == 0)) {
                std::cout << seqPos << "th entry has no sequence or name." << std::endl;
                exit(1);
            }
            readNum_1++;
            seqCnt++;
            totalReadLength += kseq->entry.sequence.l;
            int kmerCnt_int = LocalUtil::getQueryKmerNumber<int>(kseq->entry.sequence.l, spaceNum, this->kmerLen);
            if (kmerCnt_int > 0) {
                kmerCnt += (size_t) kmerCnt_int;  
            } else {
                continue;
            }
            if (bytesPerKmer * kmerCnt + ((size_t) 200 * seqCnt) > availableRam) {
                querySplits.emplace_back(start, readNum_1 - 1, kmerCnt - (size_t) kmerCnt_int, seqCnt - 1);
                kmerCnt = (size_t) kmerCnt_int;
                start = readNum_1 - 1;
                seqCnt = 1;
            }
        }
        querySplits.emplace_back(start, readNum_1, kmerCnt, seqCnt);
        delete kseq;
    } else {
        KSeqWrapper* kseq_1 = KSeqFactory(queryPath_1.c_str());
        KSeqWrapper* kseq_2 = KSeqFactory(queryPath_2.c_str());
        size_t seqPos1 = 0;
        size_t seqPos2 = 0;
        // Skip processed reads
        for (size_t i = 0; i < processedQueryNum; i++) {
            kseq_1->ReadEntry();
            kseq_2->ReadEntry();
            seqPos1++;
            seqPos2++;
        }
        size_t kmerCnt = 0;
        size_t seqCnt_1 = 0;
        size_t seqCnt_2 = 0;
        size_t start = 0;
        bool end = false;
        int kmerCnt_int_1 = 0;
        int kmerCnt_int_2 = 0;
        while(true) {
            if (kseq_1->ReadEntry()) {
                seqPos1++;
                if (unlikely(kseq_1->entry.sequence.l == 0 || kseq_1->entry.name.l == 0)) {
                    std::cout << "In file " << queryPath_1 << ", " << std::endl << "\t";
                    std::cout << seqPos1 << "th entry has no sequence or name." << std::endl;
                    exit(1);
                }
                readNum_1++;
                seqCnt_1++;
                totalReadLength += kseq_1->entry.sequence.l;
                kmerCnt_int_1 = LocalUtil::getQueryKmerNumber<int>(kseq_1->entry.sequence.l, spaceNum, this->kmerLen);
            } else {
                end = true;
            }

            if (kseq_2->ReadEntry()) {
                seqPos2++;
                if (unlikely(kseq_2->entry.sequence.l == 0 || kseq_2->entry.name.l == 0)) {
                    std::cout << "In file " << queryPath_2 << ", " << std::endl << "\t";
                    std::cout << seqPos2 << "th entry has no sequence or name." << std::endl;
                    exit(1);
                }
                readNum_2++;
                seqCnt_2++;
                totalReadLength += kseq_2->entry.sequence.l;
                kmerCnt_int_2 = LocalUtil::getQueryKmerNumber<int>(kseq_2->entry.sequence.l, spaceNum, this->kmerLen);                
            } else {
                end = true;
            }

            if (readNum_1 != readNum_2) {
                Debug(Debug::ERROR) << "The number of reads in the two files are not equal." << "\n";
                EXIT(EXIT_FAILURE);
            }

            if (kmerCnt_int_1 > 0 && kmerCnt_int_2 > 0) {
                kmerCnt += (size_t) kmerCnt_int_1 + (size_t) kmerCnt_int_2;
            } else {
                continue;
            }

            if (bytesPerKmer * kmerCnt + ((size_t) 200 * seqCnt_1) > availableRam) {
                querySplits.emplace_back(start, readNum_1 - 1, kmerCnt - ((size_t) kmerCnt_int_1 + (size_t) kmerCnt_int_2), seqCnt_1 - 1);
                kmerCnt = (size_t) kmerCnt_int_1 + (size_t) kmerCnt_int_2;
                start = readNum_1 - 1;
                seqCnt_1 = 1;
            }

            if (end) {
                querySplits.emplace_back(start, readNum_1, kmerCnt, seqCnt_1);
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