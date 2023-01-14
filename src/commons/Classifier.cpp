#include "Classifier.h"
#include "LocalParameters.h"
#include <ctime>

Classifier::Classifier(LocalParameters & par) {
    if (par.seqMode == 2){
        queryPath_1 = par.filenames[0];
        queryPath_2 = par.filenames[1];
        dbDir = par.filenames[2];
        outDir = par.filenames[3];
        jobId = par.filenames[4];
        cout << "Query file 1: " << queryPath_1 << endl;
        cout << "Query file 2: " << queryPath_2 << endl;
        cout << "Database directory: " << dbDir << endl;
        cout << "Output directory: " << outDir << endl;
        cout << "Job ID: " << jobId << endl;

    } else {
        queryPath_1 = par.filenames[0];
        dbDir = par.filenames[1];
        outDir = par.filenames[2];
        jobId = par.filenames[3];
        cout << "Query file: " << queryPath_1 << endl;
        cout << "Database directory: " << dbDir << endl;
        cout << "Output directory: " << outDir << endl;
        cout << "Job ID: " << jobId << endl;
    }

    MARKER = 16777215;
    MARKER = ~ MARKER;
    bitsForCodon = 3;
    numOfSplit = 0;
    minConsCnt = par.minConsCnt;
    minSpScore = par.minSpScore;
    verbosity = par.verbosity;

    // Mask for spaced k-mer
    size_t maskLen = par.spaceMask.length();
    mask = new uint32_t[maskLen];
    spaceNum = 0;
    spaceNum_int = 0;
    for(size_t i = 0, j = 0; i < maskLen; i++){
        mask[i] = par.spaceMask[i] - 48;
        spaceNum += (mask[i] == 0);
        spaceNum_int += (mask[i] == 0);
        if(mask[i]==1){
            unmaskedPos[j] = (int) i;
            j++;
        }
    }

    // Hamming Dist. margin
    hammingMargin = (uint8_t) par.hammingMargin;

    // Taxonomy
    const string taxonomyDirectory = dbDir + "/taxonomy";
    const string names = taxonomyDirectory + "/names.dmp";
    const string nodes = taxonomyDirectory + "/nodes.dmp";
    const string merged = taxonomyDirectory + "/merged.dmp";
    taxonomy = new NcbiTaxonomy(names, nodes, merged);

    // Taxonomy ID list
    // Load the taxonomical ID list
    FILE * taxIdFile;
    if((taxIdFile = fopen((dbDir + "/taxID_list").c_str(),"r")) == NULL){
        cout<<"Cannot open the taxID list file."<<endl;
        return;
    }
    char taxID[100];
    while(feof(taxIdFile) == 0)
    {
        fscanf(taxIdFile,"%s",taxID);
        this->taxIdList.push_back(atol(taxID));
    }
    fclose(taxIdFile);
    taxonomy->createTaxIdListAtRank(this->taxIdList, speciesTaxIdList, "species");
    taxonomy->createTaxIdListAtRank(speciesTaxIdList, genusTaxIdList, "genus");
    spORssp.push_back(&this->taxIdList);
    spORssp.push_back(&this->speciesTaxIdList);
}

Classifier::~Classifier() {
    delete[] mask;
    delete taxonomy;
}

static inline bool compareForLinearSearch(const QueryKmer &a, const QueryKmer &b) {
    if (a.ADkmer < b.ADkmer) {
        return true;
    } else if (a.ADkmer == b.ADkmer) {
        return (a.info.sequenceID < b.info.sequenceID);
    }
    return false;
}

void Classifier::startClassify(const LocalParameters &par) {

    size_t totalNumOfQueryKmer = 0;

    // Check if the query file is in FASTA format
    bool isFasta = false;
    ifstream check(queryPath_1.c_str());
    if (check.is_open()) {
        string line;
        getline(check, line);
        if (line[0] == '>') {
            isFasta = true;
        }
    }
    check.close();

    // Load query file
    cout << "Indexing query file ...";
    size_t totalReadLength = 0;
    MmapedData<char> queryFile{};
    MmapedData<char> queryFile2{};
    vector<Sequence> sequences;
    vector<Sequence> sequences2;
    Query *queryList;
    size_t numOfSeq = 0;

    if (par.seqMode == 1 || par.seqMode == 3) {
        queryFile = mmapData<char>(queryPath_1.c_str());
        madvise(queryFile.data, queryFile.fileSize, MADV_SEQUENTIAL);
        Util::touchMemory(queryFile.data, queryFile.fileSize);

        // Get start and end positions of each read
        if (isFasta) {
            IndexCreator::splitSequenceFile(sequences, queryFile);
        } else {
            splitFASTQ(sequences, queryPath_1);
        }

        // Allocate memory for query list
        numOfSeq = sequences.size();
        queryList = new Query[numOfSeq];

        // Calculate the total read length
        for (size_t i = 0; i < numOfSeq; i++) {
            totalReadLength += sequences[i].length;
        }
    } else if (par.seqMode == 2) {
        queryFile = mmapData<char>(queryPath_1.c_str());
        queryFile2 = mmapData<char>(queryPath_2.c_str());
        madvise(queryFile.data, queryFile.fileSize, MADV_SEQUENTIAL);
        madvise(queryFile2.data, queryFile2.fileSize, MADV_SEQUENTIAL);
        Util::touchMemory(queryFile.data, queryFile.fileSize);
        Util::touchMemory(queryFile2.data, queryFile2.fileSize);

        // Get start and end positions of each read
        if (isFasta) {
            IndexCreator::splitSequenceFile(sequences, queryFile);
            IndexCreator::splitSequenceFile(sequences2, queryFile2);
        } else {
            splitFASTQ(sequences, queryPath_1);
            splitFASTQ(sequences2, queryPath_2);
        }

        if (sequences.size() != sequences2.size()) {
            Debug(Debug::ERROR) << "The number of reads in the two files are not equal." << "\n";
            EXIT(EXIT_FAILURE);
        }

        // Allocate memory for query list
        numOfSeq = sequences.size();
        queryList = new Query[numOfSeq];

        // Calculate the total read length
        for (size_t i = 0; i < numOfSeq; i++) {
            totalReadLength += sequences[i].length;
            totalReadLength += sequences2[i].length;
        }
    }
    cout << "Done" << endl;
    cout << "Total number of sequences: " << numOfSeq << endl;
    cout << "Total read length: " << totalReadLength <<  "nt" << endl;

    // Allocate memory for buffers
    // 1. Calculate estimated maximum RAM usage
    size_t estimatedNumOfKmer = 2 * totalReadLength - 84 * numOfSeq;
    size_t memoryForReads = 200 * numOfSeq; // 200 bytes per read
    size_t memoryForThreads = (size_t) 128'000'000 * (size_t) par.threads; // 128 MB per thread
    size_t memoryForQueryKmer = estimatedNumOfKmer * sizeof(QueryKmer);
    size_t memoryForKmerMatch = estimatedNumOfKmer * sizeof(Match) * 5;
    size_t estimatedMaxRamUsage = memoryForReads + memoryForThreads + memoryForQueryKmer + memoryForKmerMatch;

    // 2. Check if the estimated memory usage is larger than the available memory
    size_t maxCount = 0;
    if (estimatedMaxRamUsage / 1'000'000'000 >= (size_t) par.ramUsage ) {
        maxCount = ((size_t) par.ramUsage * 1'000'000'000 - memoryForReads - memoryForThreads) /
                (96 * (par.ramUsage + 32) / par.ramUsage);
        memoryForQueryKmer = maxCount * sizeof(QueryKmer);
        memoryForKmerMatch = maxCount * sizeof(Match) * 5;
    } else {
        maxCount = estimatedNumOfKmer;
    }

    cout << memoryForReads << " " << memoryForThreads << " " << memoryForQueryKmer << " " << memoryForKmerMatch << endl;
    cout << "TOTAL: " << memoryForReads + memoryForThreads + memoryForQueryKmer + memoryForKmerMatch << endl;
    cout << maxCount << endl;


    QueryKmerBuffer kmerBuffer(maxCount);
    Buffer<Match> matchBuffer(size_t(maxCount) * size_t(5));


    // Checker for multi-threading
    bool *processedSeqChecker = new bool[numOfSeq];
    fill_n(processedSeqChecker, numOfSeq, false);
    size_t processedSeqCnt = 0;

    //
    size_t numOfTatalQueryKmerCnt = 0;
    size_t totalMatchCnt = 0;
    // Extract k-mers from query sequences and compare them to target k-mer DB
#ifdef OPENMP
    omp_set_num_threads(par.threads);
#endif

    while (processedSeqCnt < numOfSeq) {
        time_t beforeKmerExtraction = time(nullptr);

        // Initialize query k-mer buffer and match buffer
        kmerBuffer.startIndexOfReserve = 0;
        matchBuffer.startIndexOfReserve = 0;

        // Extract query k-mer
        cout << "K-mer extraction ... " << endl;
        if (par.seqMode == 1 || par.seqMode == 3) { // Single-end short-read sequence or long-read sequence
            fillQueryKmerBufferParallel(kmerBuffer,queryFile,sequences,processedSeqChecker,processedSeqCnt,
                                        queryList, par);
        } else if (par.seqMode == 2) {
            fillQueryKmerBufferParallel_paired(kmerBuffer,queryFile,queryFile2,sequences,sequences2,
                                               processedSeqChecker, processedSeqCnt, queryList, numOfSeq, par);
        }
        numOfTatalQueryKmerCnt += kmerBuffer.startIndexOfReserve;
        cout << "Time spent for k-mer extraction: " << double(time(nullptr) - beforeKmerExtraction) << endl;

        // Sort query k-mer
        time_t beforeQueryKmerSort = time(nullptr);
        SORT_PARALLEL(kmerBuffer.buffer, kmerBuffer.buffer + kmerBuffer.startIndexOfReserve, compareForLinearSearch);
        cout << "Time spent for sorting query k-mer list: " << double(time(nullptr) - beforeQueryKmerSort) << endl;

        // Search matches between query and target k-mers
        linearSearchParallel(kmerBuffer.buffer, kmerBuffer.startIndexOfReserve, matchBuffer, par);

        // Sort matches
        time_t beforeSortMatches = time(nullptr);
        totalMatchCnt += matchBuffer.startIndexOfReserve;
        SORT_PARALLEL(matchBuffer.buffer, matchBuffer.buffer + matchBuffer.startIndexOfReserve,
                      sortMatch(this));
        cout << "Time spent for sorting matches: " << double(time(nullptr) - beforeSortMatches) << endl;

        // Classify queries based on the matches
        time_t beforeAnalyze = time(nullptr);
        analyseResultParallel(matchBuffer.buffer, matchBuffer.startIndexOfReserve, (int) numOfSeq, queryList, par);
        cout << "Time spent for analyzing: " << double(time(nullptr) - beforeAnalyze) << endl;
        cout << "The number of processed sequences: " << processedSeqCnt << " (" << (double) processedSeqCnt / (double) numOfSeq << ")" << endl;
    }
    cout << "Number of query k-mers: " << numOfTatalQueryKmerCnt << endl;
    cout << "The number of matches: " << totalMatchCnt << endl;

    // Write report files
    ofstream readClassificationFile;
    readClassificationFile.open(outDir + "/" + jobId + "_classifications.tsv");
    writeReadClassification(queryList, (int) numOfSeq, readClassificationFile);
    readClassificationFile.close();
    writeReportFile(outDir + "/" + jobId + "_report.tsv", numOfSeq, taxCounts);

    // Memory deallocation
    free(matchBuffer.buffer);
    delete[] queryList;
    delete[] processedSeqChecker;
    munmap(queryFile.data, queryFile.fileSize + 1);
    if (par.seqMode == 2) {
        munmap(queryFile2.data, queryFile2.fileSize + 1);
    }
}

void Classifier::fillQueryKmerBufferParallel(QueryKmerBuffer &kmerBuffer,
                                             MmapedData<char> &seqFile,
                                             vector<Sequence> &seqs,
                                             bool *checker,
                                             size_t &processedSeqCnt,
                                             Query *queryList,
                                             const LocalParameters &par) {
    bool hasOverflow = false;

#ifdef OPENMP
    omp_set_num_threads(par.threads);
#endif

#pragma omp parallel default(none), shared(par, checker, hasOverflow, processedSeqCnt, kmerBuffer, seqFile, seqs, cout, queryList)
    {
        SeqIterator seqIterator(par);
        size_t posToWrite;
#pragma omp for schedule(dynamic, 1)
        for (size_t i = 0; i < seqs.size(); i++) {
            if (!checker[i] && !hasOverflow) {
                kseq_buffer_t buffer(const_cast<char *>(&seqFile.data[seqs[i].start]), seqs[i].length);
                kseq_t *seq = kseq_init(&buffer);
                kseq_read(seq);
                cout << seq->seq.s << endl;
                int kmerCnt = getQueryKmerNumber((int) seq->seq.l);

                // Ignore short read
                if (kmerCnt < 1) {
                    __sync_fetch_and_add(&processedSeqCnt, 1);
                    checker[i] = true;
                    continue;
                }

                posToWrite = kmerBuffer.reserveMemory(kmerCnt);
                if (posToWrite + kmerCnt < kmerBuffer.bufferSize) {
                    checker[i] = true;

                    seqIterator.sixFrameTranslation(seq->seq.s);
                    seqIterator.fillQueryKmerBuffer(seq->seq.s, (int)seq->seq.l, kmerBuffer, posToWrite,
                                                    (int) i);

                    queryList[i].queryLength = getMaxCoveredLength((int) seq->seq.l);
                    queryList[i].queryId = (int) i;
                    queryList[i].name = string(seq->name.s);
                    queryList[i].kmerCnt = kmerCnt;

                    __sync_fetch_and_add(&processedSeqCnt, 1);
                } else {
                    __sync_fetch_and_add(&(kmerBuffer.startIndexOfReserve), -kmerCnt);
                    hasOverflow = true;
                }
                kseq_destroy(seq);
            }
        }
    }
}


int Classifier::getMaxCoveredLength(int queryLength) {
    if (queryLength % 3 == 2) {
        return queryLength - 2; // 2
    } else if (queryLength % 3 == 1) {
        return queryLength - 4; // 4
    } else {
        return queryLength - 3; // 3
    }
}

int Classifier::getQueryKmerNumber(int queryLength) {
    return (getMaxCoveredLength(queryLength) / 3 - kmerLength - spaceNum_int + 1) * 6;
}

void Classifier::fillQueryKmerBufferParallel_paired(QueryKmerBuffer &kmerBuffer,
                                                    MmapedData<char> &seqFile1,
                                                    MmapedData<char> &seqFile2,
                                                    vector<Sequence> &seqs,
                                                    vector<Sequence> &seqs2,
                                                    bool *checker,
                                                    size_t &processedSeqCnt,
                                                    Query *queryList,
                                                    size_t numOfSeq,
                                                    const LocalParameters &par) {
    bool hasOverflow = false;

#pragma omp parallel default(none), shared(par, checker, hasOverflow, processedSeqCnt, kmerBuffer, seqFile1, seqFile2, seqs, seqs2, cout, queryList, numOfSeq)
    {
        SeqIterator seqIterator(par);
        SeqIterator seqIterator2(par);
        size_t posToWrite;
#pragma omp for schedule(dynamic, 1)
        for (size_t i = 0; i < numOfSeq; i++) {
            if (!checker[i] && !hasOverflow) {
                // Read 1
                kseq_buffer_t buffer(const_cast<char *>(&seqFile1.data[seqs[i].start]), seqs[i].length);
                kseq_t *seq = kseq_init(&buffer);
                kseq_read(seq);
                int kmerCnt = getQueryKmerNumber((int) seq->seq.l);

                // Read 2
                kseq_buffer_t buffer2(const_cast<char *>(&seqFile2.data[seqs2[i].start]), seqs2[i].length);
                kseq_t *seq2 = kseq_init(&buffer2);
                kseq_read(seq2);
                int kmerCnt2 = getQueryKmerNumber((int) seq2->seq.l);

                // Ignore short read
                if (kmerCnt2 < 1 || kmerCnt < 1) {
                    __sync_fetch_and_add(&processedSeqCnt, 1);
                    checker[i] = true;
                    continue;
                }

                posToWrite = kmerBuffer.reserveMemory(kmerCnt + kmerCnt2);
                if (posToWrite + kmerCnt + kmerCnt2 < kmerBuffer.bufferSize) {
                    checker[i] = true;
                    // Read 1
                    seqIterator.sixFrameTranslation(seq->seq.s);
                    seqIterator.fillQueryKmerBuffer(seq->seq.s, (int)seq->seq.l, kmerBuffer, posToWrite, (int) i);
                    queryList[i].queryLength = getMaxCoveredLength((int) seq->seq.l);

                    // Read 2
                    seqIterator2.sixFrameTranslation(seq2->seq.s);
                    seqIterator2.fillQueryKmerBuffer(seq2->seq.s, (int)seq2->seq.l, kmerBuffer, posToWrite, (int) i,
                                                     queryList[i].queryLength);

                    // Query Info
                    queryList[i].queryLength2 = getMaxCoveredLength((int) seq2->seq.l);
                    queryList[i].queryId = (int) i;
                    queryList[i].name = string(seq->name.s);
                    queryList[i].kmerCnt = kmerCnt + kmerCnt2;

                    __sync_fetch_and_add(&processedSeqCnt, 1);
                } else {
                    __sync_fetch_and_add(&(kmerBuffer.startIndexOfReserve), -(kmerCnt + kmerCnt2));
                    hasOverflow = true;
                }
                kseq_destroy(seq);
                kseq_destroy(seq2);
            }
        }
    }
}

void Classifier::linearSearchParallel(QueryKmer *queryKmerList, size_t &queryKmerCnt,
                                      Buffer<Match> &matchBuffer, const LocalParameters &par) {
    int threadNum = par.threads;
    string targetDiffIdxFileName = dbDir + "/diffIdx";
    string targetInfoFileName = dbDir + "/info";
    string diffIdxSplitFileName = dbDir + "/split";;

    struct stat diffIdxFileSt{};
    stat(targetDiffIdxFileName.c_str(), &diffIdxFileSt);
    size_t numOfDiffIdx = diffIdxFileSt.st_size / sizeof(uint16_t);

    struct MmapedData<DiffIdxSplit> diffIdxSplits = mmapData<DiffIdxSplit>(diffIdxSplitFileName.c_str());

    cout << "linearSearch start..." << endl;
    SeqIterator seqIterator1(par);
    // Find the first index of garbage query k-mer (UINT64_MAX) and discard from there
    for (size_t checkN = queryKmerCnt - 1; checkN > 0; checkN--) {
        if (queryKmerList[checkN].ADkmer != UINT64_MAX) {
            queryKmerCnt = checkN + 1;
            break;
        }
    }

    // Filter out meaningless target splits
    size_t numOfDiffIdxSplits = diffIdxSplits.fileSize / sizeof(DiffIdxSplit);
    size_t numOfDiffIdxSplits_use = numOfDiffIdxSplits;
    for (size_t i = 1; i < numOfDiffIdxSplits; i++) {
        if (diffIdxSplits.data[i].ADkmer == 0 || diffIdxSplits.data[i].ADkmer == UINT64_MAX) {
            diffIdxSplits.data[i] = {UINT64_MAX, UINT64_MAX, UINT64_MAX};
            numOfDiffIdxSplits_use--;
        }
    }

    // Divide query k-mer list into blocks for multi threading.
    // Each split has start and end points of query list + proper offset point of target k-mer list
    vector<QueryKmerSplit> querySplits;
    uint64_t queryAA;
    if (threadNum == 1) { //Single thread
        querySplits.emplace_back(0, queryKmerCnt - 1, queryKmerCnt, diffIdxSplits.data[0]);
    } else if (threadNum == 2) { //Two threads
        size_t splitWidth = queryKmerCnt / 2;
        querySplits.emplace_back(0, splitWidth - 1, splitWidth, diffIdxSplits.data[0]);
        for (size_t tSplitCnt = 0; tSplitCnt < numOfDiffIdxSplits_use; tSplitCnt++) {
            queryAA = AminoAcidPart(queryKmerList[splitWidth].ADkmer);
            if (queryAA <= AminoAcidPart(diffIdxSplits.data[tSplitCnt].ADkmer)) {
                tSplitCnt = tSplitCnt - (tSplitCnt != 0);
                querySplits.emplace_back(splitWidth, queryKmerCnt - 1, queryKmerCnt - splitWidth,
                                         diffIdxSplits.data[tSplitCnt]);
                break;
            }
        }
    } else { //More than two threads
        size_t splitWidth = queryKmerCnt / (threadNum - 1);
        querySplits.emplace_back(0, splitWidth - 1, splitWidth, diffIdxSplits.data[0]);
        for (int i = 1; i < threadNum; i++) {
            queryAA = AminoAcidPart(queryKmerList[splitWidth * i].ADkmer);
            bool needLastTargetBlock = true;
            for (size_t j = 0; j < numOfDiffIdxSplits_use; j++) {
                if (queryAA <= AminoAcidPart(diffIdxSplits.data[j].ADkmer)) {
                    j = j - (j != 0);
                    if (i != threadNum - 1) {
                        querySplits.emplace_back(splitWidth * i, splitWidth * (i + 1) - 1, splitWidth,
                                                 diffIdxSplits.data[j]);
                    } else {
                        querySplits.emplace_back(splitWidth * i, queryKmerCnt - 1, queryKmerCnt - splitWidth * i,
                                                 diffIdxSplits.data[j]);
                    }
                    needLastTargetBlock = false;
                    break;
                }
            }
            if (needLastTargetBlock) {
                if (i != threadNum - 1) {
                    querySplits.emplace_back(splitWidth * i, splitWidth * (i + 1) - 1, splitWidth,
                                             diffIdxSplits.data[numOfDiffIdxSplits_use - 1]);
                } else {
                    querySplits.emplace_back(splitWidth * i, queryKmerCnt - 1, queryKmerCnt - splitWidth * i,
                                             diffIdxSplits.data[numOfDiffIdxSplits_use - 1]);
                }
            }
        }
    }

    bool *splitCheckList = (bool *) malloc(sizeof(bool) * threadNum);
    fill_n(splitCheckList, threadNum, false);
    int completedSplitCnt = 0;

    time_t beforeSearch = time(nullptr);

    while (completedSplitCnt < threadNum) {
        bool hasOverflow = false;
#pragma omp parallel default(none), shared(completedSplitCnt, splitCheckList, hasOverflow, \
querySplits, queryKmerList, matchBuffer, cout, par, targetDiffIdxFileName, numOfDiffIdx, targetInfoFileName)
        {
            // FILE
            FILE * diffIdxFp = fopen(targetDiffIdxFileName.c_str(), "rb");
            FILE * kmerInfoFp = fopen(targetInfoFileName.c_str(), "rb");

            // Target K-mer buffer
            uint16_t * diffIdxBuffer = (uint16_t *) malloc(sizeof(uint16_t) * (BufferSize + 1));
            TargetKmerInfo * kmerInfoBuffer = (TargetKmerInfo *) malloc(sizeof(TargetKmerInfo) * (BufferSize+1));
            size_t kmerInfoBufferIdx = 0;
            size_t diffIdxBufferIdx = 0;

            //query variables
            uint64_t currentQuery = UINT64_MAX;
            uint64_t currentQueryAA = UINT64_MAX;

            //target variables
            size_t diffIdxPos = 0;
            vector<uint64_t> candidateTargetKmers; //vector for candidate target k-mer, some of which are selected after based on hamming distance
            vector<TargetKmerInfo> candidateKmerInfos;
            uint64_t currentTargetKmer;

            //Match buffer for each thread
            int localBufferSize = 2'000'000; // 32 Mb
            auto *matches = new Match[localBufferSize];
            int matchCnt = 0;

            // For debug
            SeqIterator seqIterator(par);

            //vectors for selected target k-mers
            vector<uint8_t> selectedHammingSum;
            vector<size_t> selectedMatches;
            vector<uint16_t> selectedHammings;
            size_t posToWrite;

            int currMatchNum;
            size_t idx;
#pragma omp for schedule(dynamic, 1)
            for (size_t i = 0; i < querySplits.size(); i++) {
                if (hasOverflow || splitCheckList[i]) {
                    continue;
                }

                currentTargetKmer = querySplits[i].diffIdxSplit.ADkmer;
                diffIdxBufferIdx = querySplits[i].diffIdxSplit.diffIdxOffset;
                kmerInfoBufferIdx = querySplits[i].diffIdxSplit.infoIdxOffset - (i != 0);
                diffIdxPos = querySplits[i].diffIdxSplit.diffIdxOffset;

                fseek(kmerInfoFp, 4 * (long)(kmerInfoBufferIdx), SEEK_SET);
                loadBuffer(kmerInfoFp, kmerInfoBuffer, kmerInfoBufferIdx, BufferSize);
                fseek(diffIdxFp, 2 * (long) (diffIdxBufferIdx), SEEK_SET);
                loadBuffer(diffIdxFp, diffIdxBuffer, diffIdxBufferIdx, BufferSize);

                if (i == 0) {
                    currentTargetKmer = getNextTargetKmer(currentTargetKmer, diffIdxBuffer,
                                                          diffIdxBufferIdx, diffIdxPos, BufferSize, diffIdxFp);
                }
                currentQuery = UINT64_MAX;
                currentQueryAA = UINT64_MAX;

                size_t lastMovedQueryIdx;
                for (size_t j = querySplits[i].start; j < querySplits[i].end + 1; j++) {
                    querySplits[i].start++;
                    // Reuse the comparison data if queries are exactly identical
                    if (currentQuery == queryKmerList[j].ADkmer) {
                        currMatchNum = selectedMatches.size();
                        // If local buffer is full, copy them to the shared buffer.
                        if (matchCnt + currMatchNum > localBufferSize) {
                            // Check if the shared buffer is full.
                            posToWrite = matchBuffer.reserveMemory(matchCnt);
                            if (posToWrite + matchCnt >= matchBuffer.bufferSize) {
                                hasOverflow = true;
                                querySplits[i].start = lastMovedQueryIdx + 1;
                                __sync_fetch_and_sub(& matchBuffer.startIndexOfReserve, matchCnt);
                                break;
                            } else { // not full -> copy matches to the shared buffer
                                moveMatches(matchBuffer.buffer + posToWrite, matches, matchCnt);
                                lastMovedQueryIdx = j;
                            }
                        }
                        for (int k = 0; k < currMatchNum; k++) {
                            idx = selectedMatches[k];
                            matches[matchCnt] = {queryKmerList[j].info.sequenceID,
                                                 candidateKmerInfos[idx].sequenceID,
                                                 queryKmerList[j].info.pos,
                                                 selectedHammings[k],
                                                 selectedHammingSum[k],
                                                 (bool) candidateKmerInfos[idx].redundancy};
                            matchCnt++;
                        }

                        if (currMatchNum != 0) {
                        }
                        continue;
                    }
                    selectedMatches.clear();
                    selectedHammingSum.clear();
                    selectedHammings.clear();

                    // Reuse the candidate target k-mers to compare in DNA level if queries are the same at amino acid level but not at DNA level
                    if (currentQueryAA == AminoAcidPart(queryKmerList[j].ADkmer)) {
                        compareDna(queryKmerList[j].ADkmer, candidateTargetKmers, selectedMatches,
                                   selectedHammingSum, selectedHammings);
                        currMatchNum = selectedMatches.size();

                        // If local buffer is full, copy them to the shared buffer.
                        if (matchCnt + currMatchNum > localBufferSize) {
                            // Check if the shared buffer is full.
                            posToWrite = matchBuffer.reserveMemory(matchCnt);
                            if (posToWrite + matchCnt >= matchBuffer.bufferSize) {
                                hasOverflow = true;
                                querySplits[i].start = lastMovedQueryIdx + 1;
                                __sync_fetch_and_sub(& matchBuffer.startIndexOfReserve, matchCnt);
                                break;
                            } else { // not full -> copy matches to the shared buffer
                                moveMatches(matchBuffer.buffer + posToWrite, matches, matchCnt);
                                lastMovedQueryIdx = j;
                            }
                        }

                        for (int k = 0; k < currMatchNum; k++) {
                            idx = selectedMatches[k];
                            matches[matchCnt] = {queryKmerList[j].info.sequenceID,
                                                 candidateKmerInfos[idx].sequenceID,
                                                 queryKmerList[j].info.pos,
                                                 selectedHammings[k],
                                                 selectedHammingSum[k],
                                                 (bool) candidateKmerInfos[idx].redundancy};

                            matchCnt++;
                        }
                        continue;
                    }
                    candidateTargetKmers.clear();
                    candidateKmerInfos.clear();

                    // Get next query, and start to find
                    currentQuery = queryKmerList[j].ADkmer;
                    currentQueryAA = AminoAcidPart(currentQuery);

                    // Skip target k-mers that are not matched in amino acid level
                    while (diffIdxPos != numOfDiffIdx
                        && (AminoAcidPart(currentQuery) > AminoAcidPart(currentTargetKmer))) {
                        if (unlikely(BufferSize < diffIdxBufferIdx + 7)){
                            loadBuffer(diffIdxFp, diffIdxBuffer, diffIdxBufferIdx, BufferSize, ((int)(BufferSize - diffIdxBufferIdx)) * -1 );
                        }
                        currentTargetKmer = getNextTargetKmer(currentTargetKmer, diffIdxBuffer,
                                                              diffIdxBufferIdx, diffIdxPos,
                                                              BufferSize, diffIdxFp);
                        kmerInfoBufferIdx ++;
                    }

                    if (AminoAcidPart(currentQuery) != AminoAcidPart(currentTargetKmer)) // Move to next query k-mer if there isn't any match.
                        continue;
                    else

                    // Load target k-mers that are matched in amino acid level
                    while (diffIdxPos != numOfDiffIdx &&
                        AminoAcidPart(currentQuery) == AminoAcidPart(currentTargetKmer)) {
//                        // Print the target k-mer
//                        cout << queryKmerList[j].info.sequenceID << "\t" << queryKmerList[j].info.pos << "\t" << (int) queryKmerList[j].info.frame << endl;
//                        cout << "Query  k-mer: " ;
//                        print_binary64(64, currentQuery); cout << "\t";
//                        seqIterator.printKmerInDNAsequence(currentQuery); cout << endl;
//                        cout << "Target k-mer: " ;
//                        print_binary64(64, currentTargetKmer); cout << "\t";
//                        seqIterator.printKmerInDNAsequence(currentTargetKmer); cout << "\t" << taxIdList[kmerInfoBuffer[kmerInfoBufferIdx].sequenceID] << endl;
//                        cout << (int) getHammingDistanceSum(currentQuery, currentTargetKmer) << endl;
                        candidateTargetKmers.push_back(currentTargetKmer);
                        candidateKmerInfos.push_back(getKmerInfo(BufferSize, kmerInfoFp, kmerInfoBuffer, kmerInfoBufferIdx));

                        if (unlikely(BufferSize < diffIdxBufferIdx + 7)){
                            loadBuffer(diffIdxFp, diffIdxBuffer, diffIdxBufferIdx,
                                       BufferSize, ((int)(BufferSize - diffIdxBufferIdx)) * -1 );
                        }

                        currentTargetKmer = getNextTargetKmer(currentTargetKmer, diffIdxBuffer,
                                                              diffIdxBufferIdx, diffIdxPos, BufferSize, diffIdxFp);
                        kmerInfoBufferIdx ++;
                    }

                    // Compare the current query and the loaded target k-mers and select
                    compareDna(currentQuery, candidateTargetKmers, selectedMatches, selectedHammingSum, selectedHammings);

                    // If local buffer is full, copy them to the shared buffer.
                    currMatchNum = selectedMatches.size();
                    if (matchCnt + currMatchNum > localBufferSize) {
                        // Check if the shared buffer is full.
                        posToWrite = matchBuffer.reserveMemory(matchCnt);
                        if (posToWrite + matchCnt >= matchBuffer.bufferSize) { // full -> write matches to file first
                            hasOverflow = true;
                            querySplits[i].start = lastMovedQueryIdx + 1;
                            __sync_fetch_and_sub(& matchBuffer.startIndexOfReserve, matchCnt);
                            break;
                        } else { // not full -> copy matches to the shared buffer
                            moveMatches(matchBuffer.buffer + posToWrite, matches, matchCnt);
                            lastMovedQueryIdx = j;
                        }
                    }

                    for (int k = 0; k < currMatchNum; k++) {
                        idx = selectedMatches[k];
                        matches[matchCnt] = {queryKmerList[j].info.sequenceID,
                                             candidateKmerInfos[idx].sequenceID,
                                             queryKmerList[j].info.pos,
                                             selectedHammings[k],
                                             selectedHammingSum[k],
                                             (bool) candidateKmerInfos[idx].redundancy};
                        matchCnt++;
                    }
                } // End of one split

                // Move matches in the local buffer to the shared buffer
                posToWrite = matchBuffer.reserveMemory(matchCnt);
                if (posToWrite + matchCnt >= matchBuffer.bufferSize) {
                    hasOverflow = true;
                    querySplits[i].start = lastMovedQueryIdx + 1;
                    __sync_fetch_and_sub(& matchBuffer.startIndexOfReserve, matchCnt);
                } else {
                    moveMatches(matchBuffer.buffer + posToWrite, matches, matchCnt);
                }

                // Check whether current split is completed or not
                if (querySplits[i].start - 1 == querySplits[i].end) {
                    splitCheckList[i] = true;
                    __sync_fetch_and_add(&completedSplitCnt, 1);
                }
            } // End of omp for (Iterating for splits)
            delete[] matches;
            fclose(diffIdxFp);
            fclose(kmerInfoFp);
            free(diffIdxBuffer);
            free(kmerInfoBuffer);
        } // End of omp parallel
        if (hasOverflow) {
            cout << "overflow!!!" << endl;
            break;
        }
    } // end of while(completeSplitCnt < threadNum)
    cout << "Time spent for linearSearch: " << double(time(nullptr) - beforeSearch) << endl;
    munmap(diffIdxSplits.data, diffIdxSplits.fileSize + 1);
    free(splitCheckList);
    queryKmerCnt = 0;
}

void Classifier::moveMatches(Match *dest, Match *src, int &matchNum) {
    memcpy(dest, src, sizeof(Match) * matchNum);
    matchNum = 0;
}

// It compares query k-mers to target k-mers.
// If a query has matches, the matches with the smallest hamming distance will be selected
void Classifier::compareDna(uint64_t query, vector<uint64_t> &targetKmersToCompare,
                            vector<size_t> &selectedMatches, vector<uint8_t> &selectedHammingSum,
                            vector<uint16_t> &selectedHammings) {

    size_t size = targetKmersToCompare.size();
    uint8_t *hammingSums = new uint8_t[size + 1];
    uint8_t currentHammingSum;
    uint8_t minHammingSum = UINT8_MAX;

    // Calculate hamming distance
    for (size_t i = 0; i < size; i++) {
        currentHammingSum = getHammingDistanceSum(query, targetKmersToCompare[i]);
        if (currentHammingSum < minHammingSum) {
            minHammingSum = currentHammingSum;
        }
        hammingSums[i] = currentHammingSum;
    }

    // Select target k-mers that passed hamming criteria
    for (size_t h = 0; h < size; h++) {
        if (hammingSums[h] <= minHammingSum + hammingMargin) {
            selectedMatches.push_back(h);
            selectedHammingSum.push_back(hammingSums[h]);
            selectedHammings.push_back(getHammings(query, targetKmersToCompare[h]));
        }
    }
    delete[] hammingSums;
}

// It analyses the result of linear search.
void Classifier::analyseResultParallel(Match *matchList,
                                       size_t numOfMatches,
                                       int seqNum,
                                       Query *queryList,
                                       const LocalParameters &par) {

    // Devide matches into blocks for multi threading
    cout << "Devide matches into blocks for multi threading" << endl;
    MatchBlock *matchBlocks = new MatchBlock[seqNum];
    size_t matchIdx = 0;
    size_t blockIdx = 0;
    uint32_t currentQuery;
    while (matchIdx < numOfMatches) {
        currentQuery = matchList[matchIdx].queryId;
        matchBlocks[blockIdx].id = currentQuery;
        matchBlocks[blockIdx].start = matchIdx;
        while ((currentQuery == matchList[matchIdx].queryId) && (matchIdx < numOfMatches)) ++matchIdx;
        matchBlocks[blockIdx].end = matchIdx - 1;
        blockIdx++;
    }

#ifdef OPENMP
    omp_set_num_threads(par.threads);
#endif

    // Process each block
#pragma omp parallel default(none), shared(cout, matchBlocks, matchList, seqNum, queryList, blockIdx, par)
    {
#pragma omp for schedule(dynamic, 1)
        for (size_t i = 0; i < blockIdx; ++i) {
            chooseBestTaxon(matchBlocks[i].id,
                            matchBlocks[i].start,
                            matchBlocks[i].end,
                            matchList,
                            queryList,
                            par);
        }
    }


    for (size_t i = 0; i < blockIdx; i++) {
        ++taxCounts[queryList[matchBlocks[i].id].classification];
    }
    delete[] matchBlocks;
    cout << "End of analyseResultParallel" << endl;
}


void Classifier::chooseBestTaxon(uint32_t currentQuery,
                                 size_t offset, size_t end, Match *matchList, Query *queryList,
                                 const LocalParameters &par) {
//    int queryLength = queryList[currentQuery].queryLength;
    TaxID selectedTaxon;
//    if (par.verbosity == 4) {
//        cout << "# " << currentQuery << endl;
//        for (size_t i = offset; i < end + 1; i++) {
//            cout << genusTaxIdList[matchList[i].targetId] << " " << speciesTaxIdList[matchList[i].targetId] << " " <<
//            taxIdList[matchList[i].targetId] << " " << matchList[i].position << " " << int(matchList[i].hamming) << endl;
//        }
//    }

    // Get the best genus for current query
    vector<Match> genusMatches;
    genusMatches.reserve(end - offset + 1);

    int res;
    TaxonScore genusScore(0, 0, 0, 0);
    if (par.seqMode == 2) {
        genusScore = getBestGenusMatches(genusMatches, matchList, end, offset,
                                         queryList[currentQuery].queryLength,
                                         queryList[currentQuery].queryLength2);
    } else {
        genusScore = getBestGenusMatches(genusMatches, matchList, end, offset,
                                         queryList[currentQuery].queryLength);
    }

//    if (par.verbosity == 4) {
//        cout << "# " << currentQuery << " filtered\n";
//        for (size_t i = 0; i < genusMatches.size(); i++) {
//            cout << genusTaxIdList[genusMatches[i].targetId] << " " << speciesTaxIdList[genusMatches[i].targetId] << " " <<
//                 taxIdList[genusMatches[i].targetId] << " " << genusMatches[i].position << " " << int(genusMatches[i].hamming) << "\n";
//        }
//        cout << "Genus score: " << genusScore.score << "\n";
//    }

    // If there is no proper genus for current query, it is un-classified.
    if (genusScore.score == 0 || genusScore.coverage < par.minCoverage || genusScore.score < par.minScore) {
        queryList[currentQuery].isClassified = false;
        queryList[currentQuery].classification = 0;
        queryList[currentQuery].score = genusScore.score;
        queryList[currentQuery].coverage = genusScore.coverage;
        queryList[currentQuery].hammingDist = genusScore.hammingDist;
        queryList[currentQuery].newSpecies = false;
        return;
    }

    // If there are two or more good genus level candidates, find the LCA.
    if (genusScore.taxId == 0) {
        vector<TaxID> genusList;
        genusList.reserve(genusMatches.size());
        for (auto & genusMatch : genusMatches) {
            genusList.push_back(genusTaxIdList[genusMatch.targetId]);
        }
        selectedTaxon = taxonomy->LCA(genusList)->taxId;
        queryList[currentQuery].isClassified = true;
        queryList[currentQuery].classification = selectedTaxon;
        queryList[currentQuery].score = genusScore.score;
        queryList[currentQuery].coverage = genusScore.coverage;
        queryList[currentQuery].hammingDist = genusScore.hammingDist;
        for (auto & genusMatch : genusMatches) {
            queryList[currentQuery].taxCnt[spORssp[genusMatch.redundacny]->operator[](genusMatch.targetId)]++;
        }

//        if (par.verbosity == 4) {
//            cout << "# " << currentQuery << " " << res << endl;
//            for (size_t i = 0; i < genusMatches.size(); i++) {
//                cout << i << " " << genusMatches[i].position << " " <<
//                     taxIdList[genusMatches[i].targetId] << " " << int(genusMatches[i].hamming) << " "
//                     << endl;
//            }
//            cout << "Genus score: " << genusScore.score << " " << selectedTaxon << " "
//                 << taxonomy->taxonNode(selectedTaxon)->rank << endl;
//        }
        return;
    }

    // Choose the species with the highest coverage.
    TaxonScore speciesScore;
    vector<TaxID> species;
    if (par.seqMode == 2) {
        speciesScore = chooseSpecies(genusMatches,
                                     queryList[currentQuery].queryLength,
                                     queryList[currentQuery].queryLength2,
                                     species);
    } else {
        speciesScore = chooseSpecies(genusMatches, queryList[currentQuery].queryLength, species);
    }

    // Classify to LCA if more than one species are selected
    if (species.size() > 1) {
        queryList[currentQuery].isClassified = true;
        queryList[currentQuery].classification = taxonomy->LCA(species)->taxId;
        queryList[currentQuery].score = genusScore.score;
        queryList[currentQuery].coverage = genusScore.coverage;
        queryList[currentQuery].hammingDist = genusScore.hammingDist;
        for (auto & genusMatch : genusMatches) {
            queryList[currentQuery].taxCnt[spORssp[genusMatch.redundacny]->operator[](genusMatch.targetId)]++;
        }
        return;
    }

    // If score is not enough, classify to the parent of the selected species
    if (speciesScore.score < minSpScore) {
        queryList[currentQuery].isClassified = true;
        queryList[currentQuery].classification = taxonomy->taxonNode(
                taxonomy->getTaxIdAtRank(species[0], "species"))->parentTaxId;
        queryList[currentQuery].score = genusScore.score;
        queryList[currentQuery].coverage = genusScore.coverage;
        queryList[currentQuery].hammingDist = genusScore.hammingDist;
        for (auto & genusMatch : genusMatches) {
            if(speciesTaxIdList[genusMatch.targetId] == species[0]){
                queryList[currentQuery].taxCnt[spORssp[genusMatch.redundacny]->operator[](genusMatch.targetId)]++;
            }
        }
        return;
    }
    selectedTaxon = species[0];

    // Check if it can be classified at rank lower than species.
    int numOfstrains = 0;
    TaxID strainID = 0;
    int count = 1;
    int minStrainSpecificCnt = 1;
    if (par.seqMode == 1) {
        minStrainSpecificCnt = 1;
    } else if (par.seqMode == 2) {
        minStrainSpecificCnt = 2;
    } else if (par.seqMode == 3) {
        minStrainSpecificCnt = 3;
        if (queryList[currentQuery].queryLength > 3000) {
            minStrainSpecificCnt = queryList[currentQuery].queryLength / 1000;
        }
    }
    unordered_map<TaxID, int> strainMatchCnt;
    for (auto & genusMatch : genusMatches){
        if(speciesTaxIdList[genusMatch.targetId] == selectedTaxon){
            // Record matches of selected species
            queryList[currentQuery].taxCnt[spORssp[genusMatch.redundacny]->operator[](genusMatch.targetId)]++;
            // Count the number of strain-specific matches
            if (!genusMatch.redundacny){
                strainMatchCnt[taxIdList[genusMatch.targetId]]++;
            }
        }
    }

    // Count the number of strains with enough strain-specific matches
    for (auto strainIt = strainMatchCnt.begin(); strainIt != strainMatchCnt.end(); strainIt++) {
        if (strainIt->second > minStrainSpecificCnt) {
            strainID = strainIt->first;
            numOfstrains++;
            count = strainIt->second;
        }
    }

    // If there are multiple strains with enough strain-specific matches, classify to the LCA of the strains.
    if (numOfstrains > 1) {
        vector<TaxID> strainList;
        strainList.reserve(strainMatchCnt.size());
        for (auto strainIt = strainMatchCnt.begin(); strainIt != strainMatchCnt.end(); strainIt++) {
            if (strainIt->second > minStrainSpecificCnt) {
                strainList.push_back(strainIt->first);
            }
        }
        selectedTaxon = taxonomy->LCA(strainList)->taxId;
    }
    // If there is only one strain with enough strain-specific matches, classify to the strain.
    else if (numOfstrains == 1 && count > minStrainSpecificCnt + 1) {
        selectedTaxon = strainID;
    }

    // Store classification results
    queryList[currentQuery].isClassified = true;
    queryList[currentQuery].classification = selectedTaxon;
    queryList[currentQuery].score = speciesScore.score;
    queryList[currentQuery].coverage = speciesScore.coverage;
    queryList[currentQuery].hammingDist = speciesScore.hammingDist;

//    if (par.verbosity == 4) {
//        cout << "# " << currentQuery << endl;
//        for (size_t i = 0; i < genusMatches.size(); i++) {
//            cout << i << " " << genusMatches[i].position << " " <<
//            taxIdList[genusMatches[i].targetId] << " " << int(genusMatches[i].hamming) << endl;
//        }
//        cout << "Score: " << speciesScore.score << "  " << selectedSpecies << " "
//             << taxonomy->taxonNode(selectedSpecies)->rank
//             << endl;
//    }
}

TaxonScore Classifier::getBestGenusMatches(vector<Match> & genusMatches,
                                           Match *matchList,
                                           size_t end, size_t offset,
                                           int readLength1,
                                           int readLength2) {
    TaxID currentGenus;
    TaxID currentSpecies;
    TaxonScore bestScore;
    vector<Match> tempMatchContainer;
    vector<Match> filteredMatches;
    vector<vector<Match>> matchesForEachGenus;
    vector<bool> conservedWithinGenus;
    vector<TaxonScore> genusScores;
    size_t i = offset;
    bool lastIn;
    size_t speciesMatchCnt;
    size_t speciesDiffPosCnt;
    size_t consecutiveCnt;
    int lastPos;

    while (i < end + 1) {
        currentGenus = genusTaxIdList[matchList[i].targetId];
        // For current genus
        while ((i < end + 1) && currentGenus == genusTaxIdList[matchList[i].targetId]) {
            currentSpecies = speciesTaxIdList[matchList[i].targetId];
            // For current species
            // Filter un-consecutive matches (probably random matches)
            speciesMatchCnt = 0;
            speciesDiffPosCnt = 0;
            consecutiveCnt = 0;
            lastPos = -1;
            lastIn = false;
            while ((i < end + 1) && currentSpecies == speciesTaxIdList[matchList[i + 1].targetId]) {
                if (matchList[i].position + 3 >= matchList[i + 1].position) {
                    tempMatchContainer.push_back(matchList[i]);
                    speciesMatchCnt++;
                    if (matchList[i].position / 3 != lastPos) {
                        lastPos = matchList[i].position / 3;
                        speciesDiffPosCnt++;
                        consecutiveCnt++;
                    }
                    lastIn = true;
                } else if (lastIn) {
                    lastIn = false;
                    tempMatchContainer.push_back(matchList[i]);
                    speciesMatchCnt++;
                    if (matchList[i].position / 3 != lastPos) {
                        lastPos = matchList[i].position / 3;
                        speciesDiffPosCnt++;
                        consecutiveCnt++;
                    }
                    if (consecutiveCnt >= minConsCnt) {
                        filteredMatches.insert(filteredMatches.end(), tempMatchContainer.begin(),
                                               tempMatchContainer.end());
                    }
                    consecutiveCnt = 0;
                    speciesMatchCnt = 0;
                    tempMatchContainer.clear();
                }
                i++;
            }
            if (lastIn) {
                tempMatchContainer.push_back(matchList[i]);
                speciesMatchCnt++;
                if (matchList[i].position / 3 != lastPos) {
                    lastPos = matchList[i].position / 3;
                    speciesDiffPosCnt++;
                    consecutiveCnt++;
                }
                if (consecutiveCnt >= minConsCnt) {
                    filteredMatches.insert(filteredMatches.end(), tempMatchContainer.begin(),
                                           tempMatchContainer.end());
                }
                tempMatchContainer.clear();
            }
            i++;
        }
        // Construct a match combination using filtered matches of current genus
        // so that it can best cover the query, and score the combination
        if (!filteredMatches.empty()) {
            genusScores.push_back(scoreGenus(filteredMatches, matchesForEachGenus, readLength1, readLength2));
        }
        filteredMatches.clear();
    }

    // If there are no meaningful genus
    if (genusScores.empty()) {
        bestScore.score = 0;
        return bestScore;
    }

    TaxonScore maxScore = *max_element(genusScores.begin(), genusScores.end(),
                                       [](const TaxonScore & a, const TaxonScore & b) { return a.score < b.score; });

    vector<size_t> maxIdx;

    for (size_t g = 0; g < genusScores.size(); g++) {
        if (genusScores[g].score > maxScore.score * 0.95f) {
            maxIdx.push_back(g);
        }
    }
    bestScore = maxScore;

    for (unsigned long g : maxIdx) {
        genusMatches.insert(genusMatches.end(),
                            matchesForEachGenus[g].begin(),
                            matchesForEachGenus[g].end());
    }

    // More than one genus
    if (maxIdx.size() > 1) {
        bestScore.taxId = 0;
        return bestScore;
    }
    return bestScore;

    //Three cases
    //1. one genus
    //2. more than one genus
    //4. no genus
}

TaxonScore Classifier::getBestGenusMatches(vector<Match> &genusMatches, Match *matchList, size_t end,
                                           size_t offset, int queryLength) {
    TaxID currentGenus;
    TaxID currentSpecies;

    vector<Match> tempMatchContainer;
    vector<Match> filteredMatches;
    vector<vector<Match>> matchesForEachGenus;
    vector<bool> conservedWithinGenus;
    vector<TaxonScore> genusScores;
    TaxonScore bestScore;
    size_t i = offset;
    bool lastIn;
    size_t speciesMatchCnt;
    size_t speciesDiffPosCnt;
    size_t consecutiveCnt;
    int lastPos;
    while (i < end + 1) {
        currentGenus = genusTaxIdList[matchList[i].targetId];
        // For current genus
        while ((i < end + 1) && currentGenus == genusTaxIdList[matchList[i].targetId]) {
            currentSpecies = speciesTaxIdList[matchList[i].targetId];
            // For current species
            // Filter un-consecutive matches (probably random matches)
            speciesMatchCnt = 0;
            speciesDiffPosCnt = 0;
            consecutiveCnt = 0;
            lastPos = -1;
            lastIn = false;
            while ((i < end + 1) && currentSpecies == speciesTaxIdList[matchList[i + 1].targetId]) {
                if (matchList[i].position + 3 >= matchList[i + 1].position) {
                    tempMatchContainer.push_back(matchList[i]);
                    speciesMatchCnt++;
                    if (matchList[i].position / 3 != lastPos) {
                        lastPos = matchList[i].position / 3;
                        speciesDiffPosCnt++;
                        consecutiveCnt++;
                    }
                    lastIn = true;
                } else if (lastIn) {
                    lastIn = false;
                    tempMatchContainer.push_back(matchList[i]);
                    speciesMatchCnt++;
                    if (matchList[i].position / 3 != lastPos) {
                        lastPos = matchList[i].position / 3;
                        speciesDiffPosCnt++;
                        consecutiveCnt++;
                    }
                    if (consecutiveCnt >= minConsCnt) {
                        filteredMatches.insert(filteredMatches.end(), tempMatchContainer.begin(),
                                               tempMatchContainer.end());
                    }
                    consecutiveCnt = 0;
                    speciesMatchCnt = 0;
                    tempMatchContainer.clear();
                }
                i++;
            }
            if (lastIn) {
                tempMatchContainer.push_back(matchList[i]);
                speciesMatchCnt++;
                if (matchList[i].position / 3 != lastPos) {
                    lastPos = matchList[i].position / 3;
                    speciesDiffPosCnt++;
                    consecutiveCnt++;
                }
                if (consecutiveCnt >= minConsCnt) {
                    filteredMatches.insert(filteredMatches.end(), tempMatchContainer.begin(),
                                           tempMatchContainer.end());
                }
                tempMatchContainer.clear();
            }
            i++;
        }

        // Construct a match combination using filtered matches of current genus
        // so that it can best cover the query, and score the combination
        if (!filteredMatches.empty()) {
            genusScores.push_back(scoreGenus(filteredMatches, matchesForEachGenus, queryLength));
        }
        filteredMatches.clear();
    }

    // If there are no meaningful genus
    if (genusScores.empty()) {
        bestScore.score = 0;
        return bestScore;
    }

    TaxonScore maxScore = *max_element(genusScores.begin(), genusScores.end(),
                                       [](const TaxonScore & a, const TaxonScore & b) { return a.score < b.score; });

    vector<size_t> maxIdx;
    for (size_t g = 0; g < genusScores.size(); g++) {
        if (genusScores[g].score > maxScore.score * 0.95f) {
            maxIdx.push_back(g);
        }
    }
    bestScore = maxScore;

    for (unsigned long g : maxIdx) {
        genusMatches.insert(genusMatches.end(),
                            matchesForEachGenus[g].begin(),
                            matchesForEachGenus[g].end());
    }

    // More than one genus
    if (maxIdx.size() > 1) {
        bestScore.taxId = 0;
        return bestScore;
    }
    return bestScore;

    //Three cases
    //1. one genus
    //2. more than one genus
    //4. no genus
}

TaxonScore Classifier::scoreGenus(vector<Match> &filteredMatches,
                                  vector<vector<Match>> &matchesForEachGenus,
                                  int queryLength) {
    // Do not allow overlaps between the same species
    vector<Match> matches;
    size_t walker = 0;
    size_t numOfFitMat = filteredMatches.size();
    Match &currentMatch = filteredMatches[0];
    while (walker < numOfFitMat) {
        TaxID currentSpecies = speciesTaxIdList[filteredMatches[walker].targetId];
        int currentPosition = filteredMatches[walker].position / 3;
        currentMatch = filteredMatches[walker];
        // Look through overlaps within a species
        while (walker < numOfFitMat
               && speciesTaxIdList[filteredMatches[walker].targetId] == currentSpecies
               && filteredMatches[walker].position / 3 == currentPosition) {
            // Take the match with lower hamming distance
            if (filteredMatches[walker].hamming < currentMatch.hamming) {
                currentMatch = filteredMatches[walker];
            }
                // Overlapping with same hamming distance but different subspecies taxonomy ID -> species level
            else if (taxIdList[currentMatch.targetId] != taxIdList[filteredMatches[walker].targetId] &&
                     currentMatch.hamming == filteredMatches[walker].hamming) {
                currentMatch.redundacny = true;
            }
            walker++;
        }
        matches.push_back(currentMatch);
    }

    // Calculate Hamming distance & covered length
    int coveredPosCnt = 0;
    uint16_t currHammings;
    int aminoAcidNum = (int) queryLength / 3;
    int currPos;
    size_t matchNum = matches.size();
    size_t f = 0;

    // Get the largest hamming distance at each position of query
    auto *hammingsAtEachPos = new signed char[aminoAcidNum + 1];
    memset(hammingsAtEachPos, -1, (aminoAcidNum + 1));
    while (f < matchNum) {
        currPos = matches[f].position / 3;
        currHammings = matches[f].rightEndHamming;
        if (GET_2_BITS(currHammings) > hammingsAtEachPos[currPos])
            hammingsAtEachPos[currPos + unmaskedPos[0]] = GET_2_BITS(currHammings);
        if (GET_2_BITS(currHammings >> 2) > hammingsAtEachPos[currPos + 1])
            hammingsAtEachPos[currPos + unmaskedPos[1]] = GET_2_BITS(currHammings >> 2);
        if (GET_2_BITS(currHammings >> 4) > hammingsAtEachPos[currPos + 2])
            hammingsAtEachPos[currPos + unmaskedPos[2]] = GET_2_BITS(currHammings >> 4);
        if (GET_2_BITS(currHammings >> 6) > hammingsAtEachPos[currPos + 3])
            hammingsAtEachPos[currPos + unmaskedPos[3]] = GET_2_BITS(currHammings >> 6);
        if (GET_2_BITS(currHammings >> 8) > hammingsAtEachPos[currPos + 4])
            hammingsAtEachPos[currPos + unmaskedPos[4]] = GET_2_BITS(currHammings >> 8);
        if (GET_2_BITS(currHammings >> 10) > hammingsAtEachPos[currPos + 5])
            hammingsAtEachPos[currPos + unmaskedPos[5]] = GET_2_BITS(currHammings >> 10);
        if (GET_2_BITS(currHammings >> 12) > hammingsAtEachPos[currPos + 6])
            hammingsAtEachPos[currPos + unmaskedPos[6]] = GET_2_BITS(currHammings >> 12);
        if (GET_2_BITS(currHammings >> 14) > hammingsAtEachPos[currPos + 7])
            hammingsAtEachPos[currPos + unmaskedPos[7]] = GET_2_BITS(currHammings >> 14);
        f++;
    }

    // Sum up hamming distances and count the number of position covered by the matches.
    float hammingSum = 0;
    for (int h = 0; h < aminoAcidNum; h++) {
        if (hammingsAtEachPos[h] == 0) { // Add 0 for 0 hamming dist.
            coveredPosCnt++;
        } else if (hammingsAtEachPos[h] != -1) { // Add 1.5, 2, 2.5 for 1, 2, 3 hamming dist. respectively
            hammingSum += 1.0f + (0.5f * hammingsAtEachPos[h]);
            coveredPosCnt++;
        }
    }
    delete[] hammingsAtEachPos;

    // Score current genus
    int coveredLength = coveredPosCnt * 3;
    if (coveredLength > queryLength) coveredLength = queryLength;
    float score = ((float) coveredLength - hammingSum) / (float) queryLength;
    float coverage = (float) (coveredLength) / (float) (queryLength);


//    if (verbosity == 4) {
//        cout << genusTaxIdList[filteredMatches[0].targetId] << " " << coveredLength << " " << hammingSum << " "
//             << ((float) coveredLength - hammingSum) / (float) queryLength <<
//             " " << matches.size()
//             << endl;
//    }
    matchesForEachGenus.push_back(move(matches));
    return {genusTaxIdList[filteredMatches[0].targetId], score, coverage, (int) hammingSum};
}

TaxonScore Classifier::scoreGenus(vector<Match> &filteredMatches,
                                  vector<vector<Match>> &matchesForEachGenus,
                                  int readLength1,
                                  int readLength2) {
    // Do not allow overlaps between the same species
    vector<Match> matches;
    size_t walker = 0;
    size_t numOfFitMat = filteredMatches.size();
    Match &currentMatch = filteredMatches[0];
    while (walker < numOfFitMat) {
        TaxID currentSpecies = speciesTaxIdList[filteredMatches[walker].targetId];
        int currentPosition = filteredMatches[walker].position / 3;
        currentMatch = filteredMatches[walker];
        // Look through overlaps within a species
        while (walker < numOfFitMat && speciesTaxIdList[filteredMatches[walker].targetId] == currentSpecies
               && filteredMatches[walker].position / 3 == currentPosition
               ) {
            // Take the match with lower hamming distance
            if (filteredMatches[walker].hamming < currentMatch.hamming) {
                currentMatch = filteredMatches[walker];
            }
            // Overlapping with same hamming distance but different subspecies taxonomy ID -> species level
            else if (taxIdList[currentMatch.targetId] != taxIdList[filteredMatches[walker].targetId] &&
                     currentMatch.hamming == filteredMatches[walker].hamming) {
                currentMatch.redundacny = true;
            }
            walker++;
        }
        matches.push_back(currentMatch);
    }

    // Calculate Hamming distance & covered length
    uint16_t currHammings;
    int aminoAcidNum_total = ((int) readLength1 / 3) + ((int) readLength2 / 3);
    int aminoAcidNum_read1 = ((int) readLength1 / 3);
    int currPos;
    size_t matchNum = matches.size();
    size_t f = 0;

    // Get the largest hamming distance at each position of query
    auto *hammingsAtEachPos = new signed char[aminoAcidNum_total + 3];
    memset(hammingsAtEachPos, -1, (aminoAcidNum_total + 3));
    while (f < matchNum) {
        currPos = matches[f].position / 3;
        currHammings = matches[f].rightEndHamming;
        if (GET_2_BITS(currHammings) > hammingsAtEachPos[currPos])
            hammingsAtEachPos[currPos + unmaskedPos[0]] = GET_2_BITS(currHammings);
        if (GET_2_BITS(currHammings >> 2) > hammingsAtEachPos[currPos + 1])
            hammingsAtEachPos[currPos + unmaskedPos[1]] = GET_2_BITS(currHammings >> 2);
        if (GET_2_BITS(currHammings >> 4) > hammingsAtEachPos[currPos + 2])
            hammingsAtEachPos[currPos + unmaskedPos[2]] = GET_2_BITS(currHammings >> 4);
        if (GET_2_BITS(currHammings >> 6) > hammingsAtEachPos[currPos + 3])
            hammingsAtEachPos[currPos + unmaskedPos[3]] = GET_2_BITS(currHammings >> 6);
        if (GET_2_BITS(currHammings >> 8) > hammingsAtEachPos[currPos + 4])
            hammingsAtEachPos[currPos + unmaskedPos[4]] = GET_2_BITS(currHammings >> 8);
        if (GET_2_BITS(currHammings >> 10) > hammingsAtEachPos[currPos + 5])
            hammingsAtEachPos[currPos + unmaskedPos[5]] = GET_2_BITS(currHammings >> 10);
        if (GET_2_BITS(currHammings >> 12) > hammingsAtEachPos[currPos + 6])
            hammingsAtEachPos[currPos + unmaskedPos[6]] = GET_2_BITS(currHammings >> 12);
        if (GET_2_BITS(currHammings >> 14) > hammingsAtEachPos[currPos + 7])
            hammingsAtEachPos[currPos + unmaskedPos[7]] = GET_2_BITS(currHammings >> 14);
        f++;
    }

    // Sum up hamming distances and count the number of position covered by the matches.
    float hammingSum = 0;
    int coveredPosCnt_read1 = 0;
    int coveredPosCnt_read2 = 0;
    for (int h = 0; h < aminoAcidNum_total; h++) {
        // Read 1
        if (h < aminoAcidNum_read1) {
            if (hammingsAtEachPos[h] == 0) { // Add 0 for 0 hamming dist.
                coveredPosCnt_read1++;
            } else if (hammingsAtEachPos[h] != -1) { // Add 1.5, 2, 2.5 for 1, 2, 3 hamming dist. respectively
                hammingSum += 1.0f + (0.5f * (float) hammingsAtEachPos[h]);
                coveredPosCnt_read1++;
            }
        }
        // Read 2
        else {
            if (hammingsAtEachPos[h] == 0) { // Add 0 for 0 hamming dist.
                coveredPosCnt_read2++;
            } else if (hammingsAtEachPos[h] != -1) { // Add 1.5, 2, 2.5 for 1, 2, 3 hamming dist. respectively
                hammingSum += 1.0f + (0.5f * (float) hammingsAtEachPos[h]);
                coveredPosCnt_read2++;
            }
        }
    }
    delete[] hammingsAtEachPos;

    // Score current genus
    int coveredLength_read1 = coveredPosCnt_read1 * 3;
    int coveredLength_read2 = coveredPosCnt_read2 * 3;
    if (coveredLength_read1 > readLength1) coveredLength_read1 = readLength1;
    if (coveredLength_read2 > readLength2) coveredLength_read2 = readLength2;
    float score =
            ((float) (coveredLength_read1 + coveredLength_read2) - hammingSum) / (float) (readLength1 + readLength2);
    float coverage = (float) (coveredLength_read1 + coveredLength_read2) / (float) (readLength1 + readLength2);

//    if (verbosity == 4) {
//        cout << genusTaxIdList[filteredMatches[0].targetId] << " " << coveredLength_read1 + coveredLength_read2 << " " << hammingSum
//             << " " << score <<
//             " " << matches.size()
//             << endl;
//    }
    matchesForEachGenus.push_back(move(matches));
    return {genusTaxIdList[filteredMatches[0].targetId], score, coverage, (int) hammingSum};
}

TaxonScore Classifier::chooseSpecies(const vector<Match> &matches, int queryLength, vector<TaxID> &species) {
    // Score each species
    std::unordered_map<TaxID, TaxonScore> speciesScores;
    size_t i = 0;
    TaxID currentSpeices;
    size_t numOfMatch = matches.size();
    size_t speciesBegin, speciesEnd;
    while (i < numOfMatch) {
        currentSpeices = speciesTaxIdList[matches[i].targetId];
        speciesBegin = i;
        while ((i < numOfMatch) && currentSpeices == speciesTaxIdList[matches[i].targetId]) {
            i++;
        }
        speciesEnd = i;
        speciesScores[currentSpeices] = scoreTaxon(matches, speciesBegin, speciesEnd, queryLength);
        speciesScores[currentSpeices].taxId = currentSpeices;
    }

    // Get the best species
    TaxonScore bestScore;
    for (auto & sp : speciesScores) {
        if (sp.second.score > bestScore.score) {
            species.clear();
            species.push_back(sp.first);
            bestScore = sp.second;
        } else if (sp.second.coverage == bestScore.coverage) {
            species.push_back(sp.first);
        }
    }
    return bestScore;
}

TaxonScore Classifier::chooseSpecies(const vector<Match> &matches,
                                     int read1Length,
                                     int read2Length,
                                     vector<TaxID> &species) {
    // Score each species
    std::unordered_map<TaxID, TaxonScore> speciesScores;
    size_t i = 0;
    TaxID currentSpeices;
    size_t numOfMatch = matches.size();
    size_t speciesBegin, speciesEnd;
    while (i < numOfMatch) {
        currentSpeices = speciesTaxIdList[matches[i].targetId];
        speciesBegin = i;
        while ((i < numOfMatch) && currentSpeices == speciesTaxIdList[matches[i].targetId]) {
            i++;
        }
        speciesEnd = i;
        speciesScores[currentSpeices] = scoreTaxon(matches, speciesBegin, speciesEnd, read1Length, read2Length);
        speciesScores[currentSpeices].taxId = currentSpeices;
    }

    // Get the best species
    TaxonScore bestScore;
    for (auto & sp : speciesScores) {
        if (sp.second.score > bestScore.score) {
            species.clear();
            species.push_back(sp.first);
            bestScore = sp.second;
        } else if (sp.second.coverage == bestScore.coverage) {
            species.push_back(sp.first);
        }
    }
    return bestScore;
}

TaxonScore Classifier::scoreTaxon(const vector<Match> &matches,
                                  size_t begin,
                                  size_t end,
                                  int queryLength) {
    // Get the largest hamming distance at each position of query
    int aminoAcidNum = queryLength / 3;
    auto *hammingsAtEachPos = new signed char[aminoAcidNum + 1];
    memset(hammingsAtEachPos, -1, (aminoAcidNum + 1));
    int currPos;
    size_t walker = begin;
    uint16_t currHammings;
    while (walker < end) {
        currPos = matches[walker].position / 3;
        currHammings = matches[walker].rightEndHamming;
        if (GET_2_BITS(currHammings) > hammingsAtEachPos[currPos])
            hammingsAtEachPos[currPos + unmaskedPos[0]] = GET_2_BITS(currHammings);
        if (GET_2_BITS(currHammings >> 2) > hammingsAtEachPos[currPos + 1])
            hammingsAtEachPos[currPos + unmaskedPos[1]] = GET_2_BITS(currHammings >> 2);
        if (GET_2_BITS(currHammings >> 4) > hammingsAtEachPos[currPos + 2])
            hammingsAtEachPos[currPos + unmaskedPos[2]] = GET_2_BITS(currHammings >> 4);
        if (GET_2_BITS(currHammings >> 6) > hammingsAtEachPos[currPos + 3])
            hammingsAtEachPos[currPos + unmaskedPos[3]] = GET_2_BITS(currHammings >> 6);
        if (GET_2_BITS(currHammings >> 8) > hammingsAtEachPos[currPos + 4])
            hammingsAtEachPos[currPos + unmaskedPos[4]] = GET_2_BITS(currHammings >> 8);
        if (GET_2_BITS(currHammings >> 10) > hammingsAtEachPos[currPos + 5])
            hammingsAtEachPos[currPos + unmaskedPos[5]] = GET_2_BITS(currHammings >> 10);
        if (GET_2_BITS(currHammings >> 12) > hammingsAtEachPos[currPos + 6])
            hammingsAtEachPos[currPos + unmaskedPos[6]] = GET_2_BITS(currHammings >> 12);
        if (GET_2_BITS(currHammings >> 14) > hammingsAtEachPos[currPos + 7])
            hammingsAtEachPos[currPos + unmaskedPos[7]] = GET_2_BITS(currHammings >> 14);
        walker++;
    }

    // Sum up hamming distances and count the number of position covered by the matches.
    float hammingSum = 0;
    int hammingDist = 0;
    int coveredPosCnt = 0;
    for (int h = 0; h < aminoAcidNum; h++) {
        if (hammingsAtEachPos[h] == 0) { // Add 0 for 0 hamming dist.
            coveredPosCnt++;
        } else if (hammingsAtEachPos[h] != -1) { // Add 1.5, 2, 2.5 for 1, 2, 3 hamming dist. respectively
            hammingSum += 1.0f + (0.5f * (float) hammingsAtEachPos[h]);
            hammingDist += hammingsAtEachPos[h];
            coveredPosCnt++;
        }
    }
    delete[] hammingsAtEachPos;
    // Score
    int coveredLength = coveredPosCnt * 3;
    if (coveredLength >= queryLength) coveredLength = queryLength;

    float score = ((float)coveredLength - hammingSum) / (float) queryLength;
    float coverage = (float) coveredLength / (float) (queryLength);

    return {0, score, coverage, hammingDist};
}

TaxonScore Classifier::scoreTaxon(const vector<Match> &matches,
                                  size_t begin,
                                  size_t end,
                                  int queryLength,
                                  int queryLength2) {

    // Get the largest hamming distance at each position of query
    int aminoAcidNum_total = queryLength / 3 + queryLength2 / 3;
    int aminoAcidNum_read1 = queryLength / 3;
    auto *hammingsAtEachPos = new signed char[aminoAcidNum_total + 3];
    memset(hammingsAtEachPos, -1, (aminoAcidNum_total + 3));

    int currPos;
    size_t walker = begin;
    uint16_t currHammings;

    while (walker < end) {
        currPos = matches[walker].position / 3;
        currHammings = matches[walker].rightEndHamming;
        if (GET_2_BITS(currHammings) > hammingsAtEachPos[currPos])
            hammingsAtEachPos[currPos + unmaskedPos[0]] = GET_2_BITS(currHammings);
        if (GET_2_BITS(currHammings >> 2) > hammingsAtEachPos[currPos + 1])
            hammingsAtEachPos[currPos + unmaskedPos[1]] = GET_2_BITS(currHammings >> 2);
        if (GET_2_BITS(currHammings >> 4) > hammingsAtEachPos[currPos + 2])
            hammingsAtEachPos[currPos + unmaskedPos[2]] = GET_2_BITS(currHammings >> 4);
        if (GET_2_BITS(currHammings >> 6) > hammingsAtEachPos[currPos + 3])
            hammingsAtEachPos[currPos + unmaskedPos[3]] = GET_2_BITS(currHammings >> 6);
        if (GET_2_BITS(currHammings >> 8) > hammingsAtEachPos[currPos + 4])
            hammingsAtEachPos[currPos + unmaskedPos[4]] = GET_2_BITS(currHammings >> 8);
        if (GET_2_BITS(currHammings >> 10) > hammingsAtEachPos[currPos + 5])
            hammingsAtEachPos[currPos + unmaskedPos[5]] = GET_2_BITS(currHammings >> 10);
        if (GET_2_BITS(currHammings >> 12) > hammingsAtEachPos[currPos + 6])
            hammingsAtEachPos[currPos + unmaskedPos[6]] = GET_2_BITS(currHammings >> 12);
        if (GET_2_BITS(currHammings >> 14) > hammingsAtEachPos[currPos + 7])
            hammingsAtEachPos[currPos + unmaskedPos[7]] = GET_2_BITS(currHammings >> 14);
        walker++;
    }

    // Sum up hamming distances and count the number of position covered by the matches.
    float hammingSum = 0;
    int hammingDist = 0;
    int coveredPosCnt_read1 = 0;
    int coveredPosCnt_read2 = 0;
    for (int h = 0; h < aminoAcidNum_total; h++) {
        // Read 1
        if (h < aminoAcidNum_read1) {
            if (hammingsAtEachPos[h] == 0) { // Add 0 for 0 hamming dist.
                coveredPosCnt_read1++;
            } else if (hammingsAtEachPos[h] != -1) { // Add 1.5, 2, 2.5 for 1, 2, 3 hamming dist. respectively
                hammingSum += 1.0f + (0.5f * (float) hammingsAtEachPos[h]);
                hammingDist += hammingsAtEachPos[h];
                coveredPosCnt_read1++;
            }
        }
        // Read 2
        else {
            if (hammingsAtEachPos[h] == 0) { // Add 0 for 0 hamming dist.
                coveredPosCnt_read2++;
            } else if (hammingsAtEachPos[h] != -1) { // Add 1.5, 2, 2.5 for 1, 2, 3 hamming dist. respectively
                hammingSum += 1.0f + (0.5f * (float) hammingsAtEachPos[h]);
                hammingDist += hammingsAtEachPos[h];
                coveredPosCnt_read2++;
            }
        }
    }
    delete[] hammingsAtEachPos;

    // Score
    int coveredLength_read1 = coveredPosCnt_read1 * 3;
    int coveredLength_read2 = coveredPosCnt_read2 * 3;
    if (coveredLength_read1 >= queryLength) coveredLength_read1 = queryLength;
    if (coveredLength_read2 >= queryLength2) coveredLength_read2 = queryLength2;

    float score = ((float) (coveredLength_read1 + coveredLength_read2) - hammingSum) / (float) (queryLength + queryLength2);
    float coverage = (float) (coveredLength_read1 + coveredLength_read2) / (float) (queryLength + queryLength2);

    return {0, score, coverage, hammingDist};
}

void Classifier::writeReadClassification(Query *queryList, int queryNum, ofstream &readClassificationFile) {
    for (int i = 0; i < queryNum; i++) {
        readClassificationFile << queryList[i].isClassified << "\t" << queryList[i].name << "\t"
                               << queryList[i].classification << "\t"
                               << queryList[i].queryLength + queryList[i].queryLength2 << "\t"
                               << queryList[i].score << "\t"
                               << queryList[i].coverage << "\t"
                               << queryList[i].hammingDist << "\t"
                               << taxonomy->taxonNode(queryList[i].classification)->rank << "\t";
        for (auto it = queryList[i].taxCnt.begin(); it != queryList[i].taxCnt.end(); ++it) {
            readClassificationFile << it->first << ":" << it->second << " ";
        }
        readClassificationFile << "\n";
    }
}

void Classifier::writeReportFile(const string &reportFileName, int numOfQuery, unordered_map<TaxID, unsigned int> &taxCnt) {
    unordered_map<TaxID, TaxonCounts> cladeCounts = taxonomy->getCladeCounts(taxCnt);
    FILE *fp;
    fp = fopen(reportFileName.c_str(), "w");
    writeReport(fp, cladeCounts, numOfQuery);
    fclose(fp);
}

void Classifier::writeReport(FILE *fp, const unordered_map<TaxID, TaxonCounts> &cladeCounts, unsigned long totalReads,
                             TaxID taxID, int depth) {
    auto it = cladeCounts.find(taxID);
    unsigned int cladeCount = (it == cladeCounts.end() ? 0 : it->second.cladeCount);
    unsigned int taxCount = (it == cladeCounts.end() ? 0 : it->second.taxCount);
    if (taxID == 0) {
        if (cladeCount > 0) {
            fprintf(fp, "%.2f\t%i\t%i\t0\tno rank\tunclassified\n", 100 * cladeCount / double(totalReads), cladeCount,
                    taxCount);
        }
        writeReport(fp, cladeCounts, totalReads, 1);
    } else {
        if (cladeCount == 0) {
            return;
        }
        const TaxonNode *taxon = taxonomy->taxonNode(taxID);
        fprintf(fp, "%.2f\t%i\t%i\t%i\t%s\t%s%s\n", 100 * cladeCount / double(totalReads), cladeCount, taxCount, taxID,
                taxon->rank.c_str(), string(2 * depth, ' ').c_str(), taxon->name.c_str());
        vector<TaxID> children = it->second.children;
        sort(children.begin(), children.end(),
             [&](int a, int b) { return cladeCountVal(cladeCounts, a) > cladeCountVal(cladeCounts, b); });
        for (TaxID childTaxId: children) {
            if (cladeCounts.count(childTaxId)) {
                writeReport(fp,  cladeCounts, totalReads, childTaxId, depth + 1);
            } else {
                break;
            }
        }
    }
}

unsigned int Classifier::cladeCountVal(const std::unordered_map<TaxID, TaxonCounts> &map, TaxID key) {
    typename std::unordered_map<TaxID, TaxonCounts>::const_iterator it = map.find(key);
    if (it == map.end()) {
        return 0;
    } else {
        return it->second.cladeCount;
    }
}

// FASTQ
// First line: ID
// Second line: Sequence
// Third line: +
// Fourth line: Quality
// Repeat
// Store file pointer for the start of the first line and the end of the second line
// Because we don't need to store the quality line
void Classifier::splitFASTQ(vector<Sequence> & seqSegments, const string & queryPath) {
    ifstream fastq;
    fastq.open(queryPath);
    if (!fastq.is_open()) {
        cerr << "Error: Cannot open file " << queryPath << endl;
        exit(1);
    }

    string line;
    size_t lineCnt = 0;
    size_t start;
    size_t end;
    size_t pos;
    while (getline(fastq, line)) {
        if (lineCnt % 4 == 0){
            start = (size_t) fastq.tellg() - line.length() - 1;
        }
        if (lineCnt % 4 == 1){
            end = (size_t) fastq.tellg() - 1;
            seqSegments.emplace_back(start, end, end - start + 1);
        }
        lineCnt++;
    }
}
