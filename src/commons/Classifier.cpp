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
    if (par.taxonomyPath == "DBDIR/taxonomy/") par.taxonomyPath = dbDir + "/taxonomy/";

    MARKER = 16777215;
    MARKER = ~ MARKER;
    bitsForCodon = 3;
    numOfSplit = 0;
    minCoveredPos = par.minCoveredPos;
    minSpScore = par.minSpScore;
    verbosity = par.verbosity;
    maxGap = par.maxGap;

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
    const string names = par.taxonomyPath + "/names.dmp";
    const string nodes = par.taxonomyPath + "/nodes.dmp";
    const string merged = par.taxonomyPath + "/merged.dmp";
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
        TaxID taxId = atol(taxID);
        TaxonNode const * taxon = taxonomy->taxonNode(taxId);
        if (taxId == taxon->taxId) {
            TaxID speciesTaxID = taxonomy->getTaxIdAtRank(taxId, "species");
            TaxID genusTaxID = taxonomy->getTaxIdAtRank(taxId, "genus");
            while (taxon->taxId != speciesTaxID) {
                taxId2speciesId[taxon->taxId] = speciesTaxID;
                taxId2genusId[taxon->taxId] = genusTaxID;
                taxon = taxonomy->taxonNode(taxon->parentTaxId);
            }
            taxId2speciesId[speciesTaxID] = speciesTaxID;
            taxId2genusId[speciesTaxID] = genusTaxID;
        } else {
            TaxID speciesTaxID = taxonomy->getTaxIdAtRank(taxId, "species");
            TaxID genusTaxID = taxonomy->getTaxIdAtRank(taxId, "genus");
            while (taxon->taxId != speciesTaxID) {
                taxId2speciesId[taxon->taxId] = speciesTaxID;
                taxId2genusId[taxon->taxId] = genusTaxID;
                taxon = taxonomy->taxonNode(taxon->parentTaxId);
            }
            taxId2speciesId[speciesTaxID] = speciesTaxID;
            taxId2genusId[speciesTaxID] = genusTaxID;
            taxId2speciesId[taxId] = speciesTaxID;
            taxId2genusId[taxId] = genusTaxID;
        }
    }
    fclose(taxIdFile);

//    localIndexBufferSize =  16 * 1024 * 1024;
//    localMatchBufferSize = 2 * 1024 * 1024;
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

    // Calculate maximum number of k-mers for each iteration.
    size_t maxNumOfKmer = 0;
    size_t matchPerKmer = 4;
    size_t c = sizeof(QueryKmer) + matchPerKmer * sizeof(Match);
    size_t ram_threads = ((size_t) par.ramUsage * (size_t) 1024 * 1024 * 1024)
                            - ((size_t) 134217728 * (size_t) par.threads);
//    size_t ram_per_thread = (sizeof(TargetKmerInfo) + sizeof(uint16_t)) * localIndexBufferSize
//            + sizeof(Match) * localMatchBufferSize;
//    size_t ram_threads = ((size_t) par.ramUsage * (size_t) 1073741824) - ((size_t) ram_per_thread * (size_t) par.threads);
                           // N GB - 128 MB * M threads
//    cout << "RAM per thread: " << ram_per_thread << endl;
    cout << "The rest RAM: " << ram_threads << endl;
    // Load query file
    cout << "Indexing query file ...";
    size_t totalReadLength = 0;
    MmapedData<char> queryFile{};
    MmapedData<char> queryFile2{};
    vector<SequenceBlock> sequences;
    vector<SequenceBlock> sequences2;
    vector<pair<size_t, size_t>> queryReadSplit;
    vector<size_t> splitKmerCnt;
    size_t numOfSeq = 0;
    size_t start = 0;
    size_t kmerCnt = 0;
    size_t currentKmerCnt = 0;
    size_t seqCnt = 0;
    if (par.seqMode == 1 || par.seqMode == 3) {
        queryFile = mmapData<char>(queryPath_1.c_str());
        madvise(queryFile.data, queryFile.fileSize, MADV_SEQUENTIAL);

        // Get start and end positions of each read
        if (isFasta) {
            splitFASTA(sequences, queryPath_1);
        } else {
            splitFASTQ(sequences, queryPath_1);
        }

        // Make query read splits
        numOfSeq = sequences.size();
        for (size_t i = 0; i < numOfSeq; i++) {
            currentKmerCnt = getQueryKmerNumber<size_t>(sequences[i].seqLength);
            kmerCnt += currentKmerCnt;
            seqCnt++;
            if (c * kmerCnt + ((size_t) 200 * seqCnt) > ram_threads) {
                splitKmerCnt.push_back(kmerCnt - currentKmerCnt);
                queryReadSplit.emplace_back(start, i);
                kmerCnt = currentKmerCnt;
                start = i;
                seqCnt = 1;
            }
            totalReadLength += sequences[i].seqLength;
        }
        queryReadSplit.emplace_back(start, numOfSeq);
        splitKmerCnt.push_back(kmerCnt);
    } else if (par.seqMode == 2) {

        // Get start and end positions of each read
        if (isFasta) {
            splitFASTA(sequences, queryPath_1);
            splitFASTA(sequences2, queryPath_2);
        } else {
            splitFASTQ(sequences, queryPath_1);
            splitFASTQ(sequences2, queryPath_2);
        }

        if (sequences.size() != sequences2.size()) {
            Debug(Debug::ERROR) << "The number of reads in the two files are not equal." << "\n";
            EXIT(EXIT_FAILURE);
        }

        numOfSeq = sequences.size();
        // Make query read splits
        for (size_t i = 0; i < numOfSeq; i++) {
            totalReadLength += sequences[i].seqLength;
            totalReadLength += sequences2[i].seqLength;
            currentKmerCnt = getQueryKmerNumber<size_t>(sequences[i].seqLength) +
                    getQueryKmerNumber<size_t>(sequences2[i].seqLength);
            kmerCnt += currentKmerCnt;
            seqCnt ++;
            if (c * kmerCnt + ((size_t) 200 * seqCnt) > ram_threads) {
                splitKmerCnt.push_back(kmerCnt - currentKmerCnt);
                queryReadSplit.emplace_back(start, i);
                kmerCnt = currentKmerCnt;
                start = i;
                seqCnt = 1;
            }
        }
        queryReadSplit.emplace_back(start, numOfSeq);
        splitKmerCnt.push_back(kmerCnt);
    }
    cout << "Done" << endl;
    cout << "Total number of sequences: " << numOfSeq << endl;
    cout << "Total read length: " << totalReadLength <<  "nt" << endl;

    QueryKmerBuffer kmerBuffer;
    Buffer<Match> matchBuffer;
    vector<Query> queryList;

    size_t numOfTatalQueryKmerCnt = 0;
    size_t totalMatchCnt = 0;
    size_t processedSeqCnt = 0;

    ofstream readClassificationFile;
    readClassificationFile.open(outDir + "/" + jobId + "_classifications.tsv");
#ifdef OPENMP
    omp_set_num_threads(par.threads);
#endif
    // Extract k-mers from query sequences and compare them to target k-mer DB
    double vm, rss;
    for (size_t splitIdx = 0; splitIdx < queryReadSplit.size(); splitIdx++) {
        // Allocate memory for query list
        queryList.clear();
        queryList.resize(queryReadSplit[splitIdx].second - queryReadSplit[splitIdx].first + 1);

        // Allocate memory for query k-mer list and match list
        kmerBuffer.reallocateMemory(splitKmerCnt[splitIdx]);
        matchBuffer.reallocateMemory(splitKmerCnt[splitIdx] * matchPerKmer);

        // Initialize query k-mer buffer and match buffer
        kmerBuffer.startIndexOfReserve = 0;
        matchBuffer.startIndexOfReserve = 0;

        // Extract query k-mer
        time_t beforeKmerExtraction = time(nullptr);
        cout << "Extracting query metamers ... " << endl;
        if (par.seqMode == 1 || par.seqMode == 3) { // Single-end short-read sequence or long-read sequence
            fillQueryKmerBufferParallel(kmerBuffer, queryFile, sequences, queryList, queryReadSplit[splitIdx], par);
        } else if (par.seqMode == 2) {
            fillQueryKmerBufferParallel(kmerBuffer, sequences, sequences2, queryList, queryReadSplit[splitIdx], par);
        }
        numOfTatalQueryKmerCnt += kmerBuffer.startIndexOfReserve;
        cout << "Time spent for metamer extraction: " << double(time(nullptr) - beforeKmerExtraction) << endl;

        // Sort query k-mer
        time_t beforeQueryKmerSort = time(nullptr);
        cout << "Sorting query metamer list ..." << endl;
        SORT_PARALLEL(kmerBuffer.buffer, kmerBuffer.buffer + kmerBuffer.startIndexOfReserve, compareForLinearSearch);
        cout << "Time spent for sorting query metamer list: " << double(time(nullptr) - beforeQueryKmerSort) << endl;

//#ifdef OPENMP
//        if (par.printLog == 1) {
//            omp_set_num_threads(1);
//        } else {
//            omp_set_num_threads(par.threads);
//        }
//#endif
        // Search matches between query and target k-mers
        linearSearchParallel(kmerBuffer.buffer, kmerBuffer.startIndexOfReserve, matchBuffer, par);

#ifdef OPENMP
        omp_set_num_threads(par.threads);
#endif
        // Sort matches
        time_t beforeSortMatches = time(nullptr);
        totalMatchCnt += matchBuffer.startIndexOfReserve;
        cout << "Sorting matches ..." << endl;
        SORT_PARALLEL(matchBuffer.buffer, matchBuffer.buffer + matchBuffer.startIndexOfReserve,
                      sortMatch());
        cout << "Time spent for sorting matches: " << double(time(nullptr) - beforeSortMatches) << endl;

//        for (size_t i = 0; i < matchBuffer.startIndexOfReserve; i++) {
//            cout << matchBuffer.buffer[i].queryId << " " <<  matchBuffer.buffer[i].splitIdx << " " <<
//            matchBuffer.buffer[i].targetSplitIdx << " " << matchBuffer.buffer[i].targetId << " " <<
//            genusTaxIdList[matchBuffer.buffer[i].targetId] << " " << speciesTaxIdList[matchBuffer.buffer[i].targetId] << " "
//            << matchBuffer.buffer[i].position << " " << (int) matchBuffer.buffer[i].hamming << " " << taxIdList[matchBuffer.buffer[i].targetId] << endl;
//        }

        // Classify queries based on the matches
//#ifdef OPENMP
//        if (par.printLog == 1) {
//            omp_set_num_threads(1);
//        } else {
//            omp_set_num_threads(par.threads);
//        }
//#endif

        time_t beforeAnalyze = time(nullptr);
        cout << "Analyzing matches ..." << endl;
        fromMatchToClassification(matchBuffer.buffer, matchBuffer.startIndexOfReserve, queryList, par);
        cout << "Time spent for analyzing: " << double(time(nullptr) - beforeAnalyze) << endl;
        processedSeqCnt += queryReadSplit[splitIdx].second - queryReadSplit[splitIdx].first;
        cout << "The number of processed sequences: " << processedSeqCnt << " (" << (double) processedSeqCnt / (double) numOfSeq << ")" << endl;

        // Write classification results
        writeReadClassification(queryList,
                                (int) (queryReadSplit[splitIdx].second - queryReadSplit[splitIdx].first),
                                readClassificationFile);
    }
    cout << "Number of query k-mers: " << numOfTatalQueryKmerCnt << endl;
    cout << "The number of matches: " << totalMatchCnt << endl;
    readClassificationFile.close();

    // Write report files
    writeReportFile(outDir + "/" + jobId + "_report.tsv", numOfSeq, taxCounts);

    // Memory deallocation
    free(matchBuffer.buffer);

    munmap(queryFile.data, queryFile.fileSize + 1);
    if (par.seqMode == 2) {
        munmap(queryFile2.data, queryFile2.fileSize + 1);
    }
}

void Classifier::fillQueryKmerBufferParallel(QueryKmerBuffer &kmerBuffer,
                                             MmapedData<char> &seqFile,
                                             const vector<SequenceBlock> &seqs,
                                             vector<Query> & queryList,
                                             const pair<size_t, size_t> & currentSplit,
                                             const LocalParameters &par) {

#pragma omp parallel default(none), shared(par, kmerBuffer, seqFile, seqs, cout, queryList, currentSplit)
    {
        SeqIterator seqIterator(par);
        size_t posToWrite;
#pragma omp for schedule(dynamic, 1)
        for (size_t i = currentSplit.first; i < currentSplit.second; i++) {
            size_t queryIdx = i - currentSplit.first;
            kseq_buffer_t buffer(const_cast<char *>(&seqFile.data[seqs[i].start]), seqs[i].length);
            kseq_t *seq = kseq_init(&buffer);
            kseq_read(seq);
            auto kmerCnt = getQueryKmerNumber<size_t> (seq->seq.l);
            // Ignore short read
            if (kmerCnt < 1) {
                continue;
            }
            posToWrite = kmerBuffer.reserveMemory(kmerCnt);

            seqIterator.sixFrameTranslation(seq->seq.s);
            seqIterator.fillQueryKmerBuffer(seq->seq.s, (int)seq->seq.l, kmerBuffer, posToWrite,
                                            (int) queryIdx);

            queryList[queryIdx].queryLength = getMaxCoveredLength((int) seq->seq.l);
            queryList[queryIdx].queryId = (int) queryIdx;
            queryList[queryIdx].name = string(seq->name.s);
            queryList[queryIdx].kmerCnt = (int) kmerCnt;
            kseq_destroy(seq);
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

template <typename T>
T Classifier::getQueryKmerNumber(T queryLength) {
    return (getMaxCoveredLength(queryLength) / 3 - kmerLength - spaceNum_int + 1) * 6;
}

void Classifier::fillQueryKmerBufferParallel(QueryKmerBuffer &kmerBuffer,
                                             const vector<SequenceBlock> &seqs,
                                             const vector<SequenceBlock> &seqs2,
                                             vector<Query> & queryList,
                                             const pair<size_t, size_t> & currentSplit,
                                             const LocalParameters &par) {
    vector<pair<size_t, size_t>> queryReadSplit;
//    size_t numOfSeq = currentSplit.second - currentSplit.first;
    size_t numOfSeqPerSplit = size_t(1000);
    for (size_t i = currentSplit.first; i < currentSplit.second; i += numOfSeqPerSplit) {
        queryReadSplit.emplace_back(i, std::min(i + numOfSeqPerSplit - 1, currentSplit.second- 1));
    }

#pragma omp parallel default(none), shared(par, kmerBuffer, seqs, seqs2, cout, queryList, currentSplit, queryReadSplit, numOfSeqPerSplit)
    {
        FILE * query1 = fopen(par.filenames[0].c_str(), "r");
        FILE * query2 = fopen(par.filenames[1].c_str(), "r");
        char * readBuffer1 = (char *) malloc(3000 * numOfSeqPerSplit);
        char * readBuffer2 = (char *) malloc(3000 * numOfSeqPerSplit);
        size_t readBufferIdx1 = 0;
        size_t readBufferIdx2 = 0;
        SeqIterator seqIterator(par);
        SeqIterator seqIterator2(par);
        size_t posToWrite;

#pragma omp for schedule(dynamic, 1)
        for (size_t j = 0; j < queryReadSplit.size(); j ++) {
            // Load query reads of current split
            fseek(query1, (long) seqs[queryReadSplit[j].first].start, SEEK_SET);
            loadBuffer(query1, readBuffer1, readBufferIdx1, seqs[queryReadSplit[j].second].end - seqs[queryReadSplit[j].first].start + 1);
            fseek(query2, (long) seqs2[queryReadSplit[j].first].start, SEEK_SET);
            loadBuffer(query2, readBuffer2, readBufferIdx2, seqs2[queryReadSplit[j].second].end - seqs2[queryReadSplit[j].first].start + 1);
            for (size_t i = queryReadSplit[j].first; i <= queryReadSplit[j].second; i++) {
                size_t queryIdx = i - currentSplit.first;
                // Load Read 1
                kseq_buffer_t buffer1(readBuffer1 + seqs[i].start - seqs[queryReadSplit[j].first].start,
                                      seqs[i].length);
                kseq_t *seq1 = kseq_init(&buffer1);
                kseq_read(seq1);
                auto kmerCnt = getQueryKmerNumber<size_t>(seq1->seq.l);

                // Load Read 2
                kseq_buffer_t buffer2(readBuffer2 + seqs2[i].start - seqs2[queryReadSplit[j].first].start,
                                      seqs2[i].length);
                kseq_t *seq2 = kseq_init(&buffer2);
                kseq_read(seq2);
                auto kmerCnt2 = getQueryKmerNumber<size_t>(seq2->seq.l);

                // Ignore short read
                if (kmerCnt2 < 1 || kmerCnt < 1) {
                    continue;
                }

                posToWrite = kmerBuffer.reserveMemory(kmerCnt + kmerCnt2);

                // Process Read 1
                seqIterator.sixFrameTranslation(seq1->seq.s);
                seqIterator.fillQueryKmerBuffer(seq1->seq.s, (uint32_t) seq1->seq.l, kmerBuffer, posToWrite,
                                                (uint32_t) queryIdx);
                queryList[queryIdx].queryLength = getMaxCoveredLength((uint32_t) seq1->seq.l);

                // Process Read 2
                seqIterator2.sixFrameTranslation(seq2->seq.s);
                seqIterator2.fillQueryKmerBuffer(seq2->seq.s, (uint32_t) seq2->seq.l, kmerBuffer, posToWrite,
                                                 (uint32_t) queryIdx, queryList[queryIdx].queryLength);

                // Query Info
                queryList[queryIdx].queryLength2 = getMaxCoveredLength((int) seq2->seq.l);
                queryList[queryIdx].queryId = (int) queryIdx;
                queryList[queryIdx].name = string(seq1->name.s);
                queryList[queryIdx].kmerCnt = (int) (kmerCnt + kmerCnt2);

                kseq_destroy(seq1);
                kseq_destroy(seq2);
            }
        }
        free(readBuffer1);
        free(readBuffer2);
        fclose(query1);
        fclose(query2);
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

    cout << "Comparing qeury and reference metamers..." << endl;

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
    vector<int> targetSplitIdxs;
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
        // Devide query k-mers into blocks
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
                    targetSplitIdxs.emplace_back(j);
                    needLastTargetBlock = false;
                    break;
                }
            }
            if (needLastTargetBlock) {
                if (i != threadNum - 1) { // If it is not the last split
                    querySplits.emplace_back(splitWidth * i, splitWidth * (i + 1) - 1, splitWidth,
                                             diffIdxSplits.data[numOfDiffIdxSplits_use - 2]);
                    targetSplitIdxs.emplace_back(numOfDiffIdxSplits_use - 2);
                } else {
                    querySplits.emplace_back(splitWidth * i, queryKmerCnt - 1, queryKmerCnt - splitWidth * i,
                                             diffIdxSplits.data[numOfDiffIdxSplits_use - 2]);
                    targetSplitIdxs.emplace_back(numOfDiffIdxSplits_use - 2);
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
querySplits, queryKmerList, matchBuffer, cout, par, targetDiffIdxFileName, numOfDiffIdx, targetInfoFileName, targetSplitIdxs)
        {
            // FILE
            FILE * diffIdxFp = fopen(targetDiffIdxFileName.c_str(), "rb");
            FILE * kmerInfoFp = fopen(targetInfoFileName.c_str(), "rb");

            // Target K-mer buffer
            uint16_t * diffIdxBuffer = (uint16_t *) malloc(sizeof(uint16_t) * (BufferSize + 1)); // size = 32 Mb
            TargetKmerInfo * kmerInfoBuffer = (TargetKmerInfo *) malloc(sizeof(TargetKmerInfo) * (BufferSize+1)); // 64 Mb
            size_t kmerInfoBufferIdx = 0;
            size_t diffIdxBufferIdx = 0;

            //query variables
            uint64_t currentQuery = UINT64_MAX;
            uint64_t currentQueryAA = UINT64_MAX;
            QueryKmerInfo currentQueryInfo;

            //target variables
            size_t diffIdxPos = 0;
            vector<uint64_t> candidateTargetKmers; //vector for candidate target k-mer, some of which are selected after based on hamming distance
            vector<TargetKmerInfo> candidateKmerInfos;
            uint64_t currentTargetKmer;

            //Match buffer for each thread
            int localBufferSize = 2'000'000; // 32 Mb
            auto *matches = new Match[localBufferSize]; // 16 * 2'000'000 = 32 Mb
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
                kmerInfoBufferIdx = querySplits[i].diffIdxSplit.infoIdxOffset
                                    - (querySplits[i].diffIdxSplit.ADkmer != 0);
                diffIdxPos = querySplits[i].diffIdxSplit.diffIdxOffset;

                fseek(kmerInfoFp, 4 * (long)(kmerInfoBufferIdx), SEEK_SET);
                loadBuffer(kmerInfoFp, kmerInfoBuffer, kmerInfoBufferIdx, BufferSize);
                fseek(diffIdxFp, 2 * (long) (diffIdxBufferIdx), SEEK_SET);
                loadBuffer(diffIdxFp, diffIdxBuffer, diffIdxBufferIdx, BufferSize);

                if (i == 0) {
                    currentTargetKmer = getNextTargetKmer(currentTargetKmer, diffIdxBuffer,
                                                          diffIdxBufferIdx, diffIdxPos);
                }
                currentQuery = UINT64_MAX;
                currentQueryAA = UINT64_MAX;

                size_t lastMovedQueryIdx = 0;
                for (size_t j = querySplits[i].start; j < querySplits[i].end + 1; j++) {
                    querySplits[i].start++;

                    // Reuse the comparison data if queries are exactly identical
                    if (currentQuery == queryKmerList[j].ADkmer
                        && (currentQueryInfo.frame/3 == queryKmerList[j].info.frame/3)) {
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
                            matches[matchCnt] = {queryKmerList[j].info,
                                                 candidateKmerInfos[idx].sequenceID,
                                                 taxId2genusId[candidateKmerInfos[idx].sequenceID],
                                                 taxId2speciesId[candidateKmerInfos[idx].sequenceID],
                                                 selectedHammings[k],
                                                 selectedHammingSum[k],
                                                 (bool) candidateKmerInfos[idx].redundancy};
                            matchCnt++;
                        }
                        continue;
                    }
                    selectedMatches.clear();
                    selectedHammingSum.clear();
                    selectedHammings.clear();

                    // Reuse the candidate target k-mers to compare in DNA level if queries are the same at amino acid level but not at DNA level
                    if (currentQueryAA == AminoAcidPart(queryKmerList[j].ADkmer)) {
                        compareDna(queryKmerList[j].ADkmer, candidateTargetKmers, selectedMatches,
                                   selectedHammingSum, selectedHammings,queryKmerList[j].info.frame);
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
                            matches[matchCnt] = {queryKmerList[j].info,
                                                 candidateKmerInfos[idx].sequenceID,
                                                 taxId2genusId[candidateKmerInfos[idx].sequenceID],
                                                 taxId2speciesId[candidateKmerInfos[idx].sequenceID],
                                                 selectedHammings[k],
                                                 selectedHammingSum[k],
                                                 (bool) candidateKmerInfos[idx].redundancy};
                            matchCnt++;
                        }
                        currentQuery = queryKmerList[j].ADkmer;
                        currentQueryAA = AminoAcidPart(currentQuery);
                        currentQueryInfo = queryKmerList[j].info;
                        continue;
                    }
                    candidateTargetKmers.clear();
                    candidateKmerInfos.clear();

                    // Get next query, and start to find
                    currentQuery = queryKmerList[j].ADkmer;
                    currentQueryAA = AminoAcidPart(currentQuery);
                    currentQueryInfo = queryKmerList[j].info;

                    // Skip target k-mers that are not matched in amino acid level
                    while (diffIdxPos != numOfDiffIdx
                        && (currentQueryAA > AminoAcidPart(currentTargetKmer))) {
                        if (unlikely(BufferSize < diffIdxBufferIdx + 7)){
                            loadBuffer(diffIdxFp, diffIdxBuffer, diffIdxBufferIdx, BufferSize, ((int)(BufferSize - diffIdxBufferIdx)) * -1 );
                        }
                        currentTargetKmer = getNextTargetKmer(currentTargetKmer, diffIdxBuffer,
                                                              diffIdxBufferIdx, diffIdxPos);
                        kmerInfoBufferIdx ++;
                    }

                    if (currentQueryAA != AminoAcidPart(currentTargetKmer)) // Move to next query k-mer if there isn't any match.
                        continue;

                    // Load target k-mers that are matched in amino acid level
                    while (diffIdxPos != numOfDiffIdx &&
                    currentQueryAA == AminoAcidPart(currentTargetKmer)) {
                        candidateTargetKmers.push_back(currentTargetKmer);
                        candidateKmerInfos.push_back(getKmerInfo(BufferSize, kmerInfoFp, kmerInfoBuffer, kmerInfoBufferIdx));
                        // Print the target k-mer
//                        if (par.printLog == 1) {
//                            cout << queryKmerList[j].info.sequenceID << "\t" << queryKmerList[j].info.pos << "\t"
//                                 << (int) queryKmerList[j].info.frame << endl;
//                            cout << "Query  k-mer: ";
//                            print_binary64(64, currentQuery);
//                            cout << "\t";
//                            seqIterator.printKmerInDNAsequence(currentQuery);
//                            cout << endl;
//                            cout << "Target k-mer: ";
//                            print_binary64(64, currentTargetKmer);
//                            cout << "\t";
//                            seqIterator.printKmerInDNAsequence(currentTargetKmer);
//                            cout << "\t" << kmerInfoBuffer[kmerInfoBufferIdx].sequenceID
//                                 << "\t" << taxId2speciesId[kmerInfoBuffer[kmerInfoBufferIdx].sequenceID] << endl;
//                            cout << (int) getHammingDistanceSum(currentQuery, currentTargetKmer) << "\t";
//                            print_binary16(16, getHammings(currentQuery, currentTargetKmer)); cout << endl;
//                        }

                        if (unlikely(BufferSize < diffIdxBufferIdx + 7)){
                            loadBuffer(diffIdxFp, diffIdxBuffer, diffIdxBufferIdx,
                                       BufferSize, ((int)(BufferSize - diffIdxBufferIdx)) * -1 );
                        }

                        currentTargetKmer = getNextTargetKmer(currentTargetKmer, diffIdxBuffer,
                                                              diffIdxBufferIdx, diffIdxPos);
                        kmerInfoBufferIdx ++;
                    }

                    // Compare the current query and the loaded target k-mers and select
                    compareDna(currentQuery, candidateTargetKmers, selectedMatches, selectedHammingSum,
                               selectedHammings, queryKmerList[j].info.frame);

                    // If local buffer is full, copy them to the shared buffer.
                    currMatchNum = selectedMatches.size();
                    if (matchCnt + currMatchNum > localBufferSize) {
                        // Check if the shared buffer is full.
                        posToWrite = matchBuffer.reserveMemory(matchCnt);
                        if (posToWrite + matchCnt >= matchBuffer.bufferSize) { // full -> write matches to file first
                            hasOverflow = true;
                            querySplits[i].start = lastMovedQueryIdx + 1;
                            __sync_fetch_and_sub(&matchBuffer.startIndexOfReserve, matchCnt);
                            break;
                        } else { // not full -> copy matches to the shared buffer
                            moveMatches(matchBuffer.buffer + posToWrite, matches, matchCnt);
                            lastMovedQueryIdx = j;
                        }
                    }

                    for (int k = 0; k < currMatchNum; k++) {
                        idx = selectedMatches[k];
                        matches[matchCnt] = {queryKmerList[j].info,
                                             candidateKmerInfos[idx].sequenceID,
                                             taxId2genusId[candidateKmerInfos[idx].sequenceID],
                                             taxId2speciesId[candidateKmerInfos[idx].sequenceID],
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
    cout << "Time spent for the comparison: " << double(time(nullptr) - beforeSearch) << endl;
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
                            vector<uint16_t> &selectedHammings, uint8_t frame) {

    size_t size = targetKmersToCompare.size();
    auto *hammingSums = new uint8_t[size + 1];
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
            if (frame < 3) {
                selectedHammings.push_back(getHammings(query, targetKmersToCompare[h]));
            } else {
                selectedHammings.push_back(getHammings_reverse(query, targetKmersToCompare[h]));
            }
        }
    }
    delete[] hammingSums;
}

// It analyses the result of linear search.
void Classifier::fromMatchToClassification(const Match *matchList,
                                           size_t numOfMatches,
                                           vector<Query> & queryList,
                                           const LocalParameters &par) {

    // Devide matches into blocks for multi threading
    size_t seqNum = queryList.size();
    MatchBlock *matchBlocks = new MatchBlock[seqNum];
    size_t matchIdx = 0;
    size_t blockIdx = 0;
    uint32_t currentQuery;
    while (matchIdx < numOfMatches) {
        currentQuery = matchList[matchIdx].qInfo.sequenceID;
        matchBlocks[blockIdx].id = currentQuery;
        matchBlocks[blockIdx].start = matchIdx;
        while ((currentQuery == matchList[matchIdx].qInfo.sequenceID) && (matchIdx < numOfMatches)) ++matchIdx;
        matchBlocks[blockIdx].end = matchIdx - 1;
        blockIdx++;
    }

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

    for (size_t i = 0; i < seqNum; i++) {
        ++taxCounts[queryList[i].classification];
    }
    delete[] matchBlocks;
}


void Classifier::chooseBestTaxon(uint32_t currentQuery,
                                 size_t offset,
                                 size_t end,
                                 const Match *matchList,
                                 vector<Query> & queryList,
                                 const LocalParameters &par) {
    TaxID selectedTaxon;
//    if (par.printLog) {
//        cout << "# " << currentQuery << " " << queryList[currentQuery].name << endl;
//        for (size_t i = offset; i < end + 1; i++) {
//            cout << taxId2genusId[matchList[i].targetId] << " " << taxId2speciesId[matchList[i].targetId] <<
//            " "  << matchList[i].targetId << " " << matchList[i].qInfo.frame << " ";
//            print_binary16(16, matchList[i].rightEndHamming);
//            cout << " " << matchList[i].qInfo.position << " " << int(matchList[i].hamming) <<  " "  << int(matchList[i].redundancy) << endl;
//        }
//    }

    // Get the best genus for current query
    vector<Match> genusMatches;
    genusMatches.reserve(end - offset + 1);

    int res;
    TaxonScore genusScore(0, 0, 0, 0);
    if (par.seqMode == 2) {
        if (par.spaceMask != "11111111"){
            genusScore = getBestGenusMatches_spaced(genusMatches, matchList, end, offset,
                                                    queryList[currentQuery].queryLength,
                                                    queryList[currentQuery].queryLength2);
        } else {
            genusScore = getBestGenusMatches(genusMatches, matchList, end, offset,
                                              queryList[currentQuery].queryLength,
                                              queryList[currentQuery].queryLength2, par);
        }
    } else {
        if (par.spaceMask != "11111111") {
            genusScore = getBestGenusMatches_spaced(genusMatches, matchList, end, offset,
                                                    queryList[currentQuery].queryLength);
        } else {
            genusScore = getBestGenusMatches(genusMatches, matchList, end, offset,
                                              queryList[currentQuery].queryLength, par);
        }
    }

//    if (par.printLog) {
//        cout << "# " << currentQuery << " " << queryList[currentQuery].name << " filtered\n";
//        for (size_t i = 0; i < genusMatches.size(); i++) {
//            cout << taxId2genusId[genusMatches[i].targetId] << " " << taxId2speciesId[genusMatches[i].targetId] <<
//                 " "  << genusMatches[i].targetId << " " << genusMatches[i].qInfo.frame << " ";
//            print_binary16(16, genusMatches[i].rightEndHamming);
//            cout << " " << genusMatches[i].qInfo.position << " " << int(genusMatches[i].hamming) <<  " "  << int(genusMatches[i].redundancy) << endl;
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
            genusList.push_back(genusMatch.genusId);
        }
        selectedTaxon = taxonomy->LCA(genusList)->taxId;
        queryList[currentQuery].isClassified = true;
        queryList[currentQuery].classification = selectedTaxon;
        queryList[currentQuery].score = genusScore.score;
        queryList[currentQuery].coverage = genusScore.coverage;
        queryList[currentQuery].hammingDist = genusScore.hammingDist;
        for (auto & genusMatch : genusMatches) {
            queryList[currentQuery].taxCnt[genusMatch.targetId]++;
        }
        return;
    }

    // Choose the species with the highest coverage.
    TaxID selectedSpecies;
    TaxonScore speciesScore;
    vector<TaxID> species;
    unordered_map<TaxID, pair<int, int>> speciesMatchRange;
    if (par.seqMode == 2) {
        speciesScore = chooseSpecies(genusMatches,
                                     queryList[currentQuery].queryLength,
                                     queryList[currentQuery].queryLength2,
                                     species,
                                     speciesMatchRange);
    } else {
        speciesScore = chooseSpecies(genusMatches,
                                     queryList[currentQuery].queryLength,
                                     species,
                                     speciesMatchRange);
    }


    // Classify to LCA if more than one species are selected
    if (species.size() > 1) {
        queryList[currentQuery].isClassified = true;
        queryList[currentQuery].classification = taxonomy->LCA(species)->taxId;
        queryList[currentQuery].score = genusScore.score;
        queryList[currentQuery].coverage = genusScore.coverage;
        queryList[currentQuery].hammingDist = genusScore.hammingDist;
        for (auto & genusMatch : genusMatches) {
            queryList[currentQuery].taxCnt[genusMatch.targetId]++;
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
            if(genusMatch.speciesId == species[0]){
                queryList[currentQuery].taxCnt[genusMatch.targetId]++;
            }
        }
        return;
    }

    // Sort matches by the position of the query sequence
    selectedSpecies = species[0];
//    sort(genusMatches.begin() + speciesMatchRange[selectedSpecies].first,
//         genusMatches.begin() + speciesMatchRange[selectedSpecies].second,
//         [](const Match & a, const Match & b) {
//        if (a.qInfo.position / 3 == b.qInfo.position / 3)
//            return a.hamming < b.hamming;
//        else
//            return a.qInfo.position / 3 < b.qInfo.position / 3;
//    });

    sort(genusMatches.begin() + speciesMatchRange[selectedSpecies].first,
         genusMatches.begin() + speciesMatchRange[selectedSpecies].second,
         [](const Match & a, const Match & b) { return a.qInfo.pos > b.qInfo.pos; });


    TaxID result = lowerRankClassification(genusMatches, speciesMatchRange[selectedSpecies], selectedSpecies);

    // Record matches of selected species
    for (size_t i = speciesMatchRange[selectedSpecies].first; i < speciesMatchRange[selectedSpecies].second; i++) {
        queryList[currentQuery].taxCnt[genusMatches[i].targetId]++;
    }


    // Store classification results
    queryList[currentQuery].isClassified = true;
    queryList[currentQuery].classification = result;
    queryList[currentQuery].score = speciesScore.score;
    queryList[currentQuery].coverage = speciesScore.coverage;
    queryList[currentQuery].hammingDist = speciesScore.hammingDist;
    queryList[currentQuery].newSpecies = false;
//    if (par.printLog) {
//        cout << "# " << currentQuery << endl;
//        for (size_t i = 0; i < genusMatches.size(); i++) {
//            cout << i << " " << genusMatches[i].qInfo.position << " " <<
//            genusMatches[i].targetId << " " << int(genusMatches[i].hamming) << endl;
//        }
//        cout << "Score: " << speciesScore.score << "  " << selectedSpecies << " "
//             << taxonomy->getString(taxonomy->taxonNode(selectedSpecies)->rankIdx)
//
//             << endl;
//    }
}

TaxID Classifier::lowerRankClassification(vector<Match> &matches, pair<int, int> &matchRange, TaxID spTaxId) {
    int i = matchRange.second - 1;
    unordered_map<TaxID, unsigned int> taxCnt;

    while ( i >= matchRange.first ) {
        size_t currQuotient = matches[i].qInfo.pos / 3;
        uint8_t minHamming = matches[i].hamming;
        Match * minHammingMatch = & matches[i];
        TaxID minHammingTaxId = minHammingMatch->targetId;
        bool first = true;
        i --;
        while ( (i >= matchRange.first) && (currQuotient == matches[i].qInfo.pos / 3) ) {
            if (matches[i].hamming < minHamming) {
                minHamming = matches[i].hamming;
                minHammingMatch = & matches[i];
                minHammingTaxId = minHammingMatch->targetId;
            } else if (matches[i].hamming == minHamming) {
                minHammingTaxId = taxonomy->LCA(minHammingTaxId, matches[i].targetId);
                minHammingMatch->redundancy = true;
                matches[i].redundancy = true;
            }
            i--;
        }
        taxCnt[minHammingTaxId]++;
    }

    unordered_map<TaxID, TaxonCounts> cladeCnt;
    getSpeciesCladeCounts(taxCnt, cladeCnt, spTaxId);

    return BFS(cladeCnt, spTaxId);
}

void Classifier::getSpeciesCladeCounts(const unordered_map<TaxID, unsigned int> &taxCnt,
                                       unordered_map<TaxID, TaxonCounts> & cladeCount,
                                       TaxID speciesTaxID) {
    for (auto it = taxCnt.begin(); it != taxCnt.end(); ++it) {
//        cladeCount[it->first].taxCount = it->second;
//        cladeCount[it->first].cladeCount += it->second;
        TaxonNode const * taxon = taxonomy->taxonNode(it->first);
        cladeCount[taxon->taxId].taxCount = it->second;
        cladeCount[taxon->taxId].cladeCount += it->second;
        while (taxon->taxId != speciesTaxID) {
            if (find(cladeCount[taxon->parentTaxId].children.begin(),
                     cladeCount[taxon->parentTaxId].children.end(),
                     taxon->taxId) == cladeCount[taxon->parentTaxId].children.end()) {
                cladeCount[taxon->parentTaxId].children.push_back(taxon->taxId);
            }
            cladeCount[taxon->parentTaxId].cladeCount += it->second;
            taxon = taxonomy->taxonNode(taxon->parentTaxId);
        }
    }
}

TaxID Classifier::BFS(const unordered_map<TaxID, TaxonCounts> & cladeCnt, TaxID root) {
    if (cladeCnt.at(root).children.empty()) { // root is a leaf
        return root;
    }
    unsigned int maxCnt = 3;
    unsigned int currentCnt;
    vector<TaxID> bestChildren;
    for (auto it = cladeCnt.at(root).children.begin(); it != cladeCnt.at(root).children.end(); it++) {
        currentCnt = cladeCnt.at(*it).cladeCount;
        if (currentCnt > maxCnt) {
            bestChildren.clear();
            bestChildren.push_back(*it);
            maxCnt = currentCnt;
        } else if (currentCnt == maxCnt) {
            bestChildren.push_back(*it);
        }
    }
    if (bestChildren.size() == 1) {
        return BFS(cladeCnt, bestChildren[0]);
    } else {
        return root;
    }
}

TaxonScore Classifier::getBestGenusMatches(vector<Match> &genusMatches, const Match *matchList, size_t end,
                                           size_t offset, int readLength1, int readLength2, const LocalParameters & par) {
    TaxID currentGenus;
    TaxID currentSpecies;

    vector<const Match *> filteredMatches;
    vector<vector<const Match *>> matchesForEachGenus;
    vector<TaxonScore> genusScores;
    TaxonScore bestScore;
    size_t i = offset;
    uint8_t curFrame;
    vector<const Match *> curFrameMatches;
    while (i  < end + 1) {
//        currentGenus = taxId2genusId[matchList[i].targetId];
        currentGenus = matchList[i].genusId;
        // For current genus
        while ((i < end + 1) && currentGenus == matchList[i].genusId) {
//            currentSpecies = taxId2speciesId[matchList[i].targetId];
            currentSpecies = matchList[i].speciesId;
//            if (par.printLog) {
//                cout << currentGenus << " " << currentSpecies << endl;
//            }
            // For current species
            while ((i < end + 1) && currentSpecies == matchList[i].speciesId) {
                curFrame = matchList[i].qInfo.frame;
                curFrameMatches.clear();

                // For current frame
                while ((i < end + 1) && currentSpecies == matchList[i].speciesId
                    && curFrame == matchList[i].qInfo.frame) {
                    curFrameMatches.push_back(&matchList[i]);
                    i ++;
                }
                if (curFrameMatches.size() > 1) {
                    remainConsecutiveMatches(curFrameMatches, filteredMatches, currentGenus, par);
                }
            }
        }

        // Construct a match combination using filtered matches of current genus
        // so that it can best cover the query, and score the combination
        if (!filteredMatches.empty()) {
            matchesForEachGenus.push_back(filteredMatches);
            genusScores.push_back(scoreGenus(filteredMatches, readLength1, readLength2));
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
        for (const Match * m : matchesForEachGenus[g]) {
            genusMatches.push_back(*m);
        }
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



void Classifier::remainConsecutiveMatches(vector<const Match *> & curFrameMatches,
                                          vector<const Match *> & filteredMatches,
                                          TaxID genusId,
                                          const LocalParameters & par) {
    size_t i = 0;
    size_t end = curFrameMatches.size();
    vector<pair<const Match *, size_t>> curPosMatches; // <match, index>
    vector<pair<const Match *, size_t>> nextPosMatches;
    map<size_t, vector<size_t>> linkedMatches; // <index, linked indexes>

    size_t currPos = curFrameMatches[0]->qInfo.pos;
    while ( i < end && curFrameMatches[i]->qInfo.pos == currPos) {
        curPosMatches.emplace_back(curFrameMatches[i], i);
        i++;
    }
    while (i < end) {
        uint32_t nextPos = curFrameMatches[i]->qInfo.pos;
        while (i < end  && nextPos == curFrameMatches[i]->qInfo.pos) {
            nextPosMatches.emplace_back(curFrameMatches[i], i);
            ++ i;
        }
        // Check if current position and next position are consecutive
        if (currPos + 3 == nextPos) {
            // Compare curPosMatches and nextPosMatches
            for (auto &curPosMatch: curPosMatches) {
                for (auto &nextPosMatch: nextPosMatches) {
                    if (isConsecutive(curPosMatch.first, nextPosMatch.first)) {
                        linkedMatches[curPosMatch.second].push_back(nextPosMatch.second);
                    }
                }
            }

        }
        // Update curPosMatches and nextPosMatches
        curPosMatches = nextPosMatches;
        nextPosMatches.clear();
        currPos = nextPos;
    }
    // Print linkedMatches
//    if (par.printLog) {
//        cout << "linkedMatches: " << endl;
//        for (const auto &entry: linkedMatches) {
//            cout << entry.first << ": ";
//            for (auto &idx: entry.second) {
//                cout << idx << " ";
//            }
//            cout << endl;
//        }
//    }

    // Iterate linkedMatches to get filteredMatches
    int MIN_DEPTH = par.minConsCnt - 1;
    if (taxonomy->IsAncestor(par.eukaryotaTaxId, genusId)) {
        MIN_DEPTH = par.minConsCntEuk - 1;
    }
    unordered_set<size_t> used;
    vector<size_t> filteredMatchIdx;
    unordered_map<size_t, size_t> idx2depth;
    for (const auto& entry : linkedMatches) {
        if (!used.count(entry.first)) {
            used.insert(entry.first);
            vector<const Match*> curMatches;
            DFS(entry.first, linkedMatches, filteredMatchIdx, 0, MIN_DEPTH, used, idx2depth);
        }
    }

    // Print filteredMatchIdx
//    if (par.printLog) {
//        cout << "filteredMatchIdx: ";
//        for (auto &idx: filteredMatchIdx) {
//            cout << idx << " ";
//        }
//        cout << endl;
//    }
    for (auto &idx: filteredMatchIdx) {
        filteredMatches.push_back(curFrameMatches[idx]);
    }
}


size_t Classifier::DFS(size_t curMatchIdx, const map<size_t, vector<size_t>> & linkedMatches,
                       vector<size_t>& filteredMatches, size_t depth, size_t MIN_DEPTH, unordered_set<size_t>& used,
                       unordered_map<size_t, size_t> & idx2depth) {
    depth++;
    size_t maxDepth = 0;
    size_t returnDepth = 0;
    if (linkedMatches.find(curMatchIdx) == linkedMatches.end()) { //|| linkedMatches.at(curMatchIdx).empty()) {
        // reached a leaf node
        idx2depth[curMatchIdx] = depth;
        if (depth > MIN_DEPTH) {
            filteredMatches.push_back(curMatchIdx);
        }
        return depth;
    } else { // not a leaf node
        for (auto &nextMatchIdx: linkedMatches.at(curMatchIdx)) {
            used.insert(nextMatchIdx);
            if (idx2depth.find(nextMatchIdx) != idx2depth.end()) {
                returnDepth = idx2depth[nextMatchIdx];
                maxDepth = max(maxDepth, returnDepth);
                continue;
            }
            returnDepth = DFS(nextMatchIdx, linkedMatches, filteredMatches, depth, MIN_DEPTH, used, idx2depth);
            maxDepth = max(maxDepth, returnDepth);
        }
        if (maxDepth > MIN_DEPTH) {
            filteredMatches.push_back(curMatchIdx);
            idx2depth[curMatchIdx] = maxDepth;
        }
    }
    return maxDepth;
}

TaxonScore Classifier::getBestGenusMatches_spaced(vector<Match> &genusMatches, const Match *matchList, size_t end,
                                                  size_t offset, int readLength1, int readLength2) {
    TaxID currentGenus;
    TaxID currentSpecies;

    vector<const Match *> tempMatchContainer;
    vector<const Match *> filteredMatches;
    vector<vector<const Match *>> matchesForEachGenus;
    vector<bool> conservedWithinGenus;
    vector<TaxonScore> genusScores;
    TaxonScore bestScore;
    size_t i = offset;
    bool lastIn;
    while (i + 1 < end + 1) {
        currentGenus = matchList[i].genusId;
        // For current genus
        while ((i + 1 < end + 1) && currentGenus == matchList[i].genusId) {
//            currentSpecies = taxId2speciesId[matchList[i].targetId];
            currentSpecies = matchList[i].speciesId;
            // For current species
            // Filter un-consecutive matches (probably random matches)
            lastIn = false;
            int distance = 0;
            int diffPosCntOfCurrRange = 1;
            int dnaDist = 0;

            // For the same species
            while ((i + 1 < end + 1) && currentSpecies == matchList[i + 1].speciesId) {
                distance = matchList[i+1].qInfo.pos / 3 - matchList[i].qInfo.pos / 3;
                dnaDist = matchList[i+1].qInfo.pos - matchList[i].qInfo.pos;
                if (distance == 0) { // At the same position
                    tempMatchContainer.push_back(matchList + i);
                } else if (dnaDist < (8 + spaceNum_int + maxGap) * 3) { // Overlapping
                    lastIn = true;
                    tempMatchContainer.push_back(matchList + i);
                    diffPosCntOfCurrRange ++;
                } else { // Not consecutive --> End range
                    if (lastIn){
                        tempMatchContainer.push_back(matchList + i);
                        if (diffPosCntOfCurrRange >= minCoveredPos) {
                            filteredMatches.insert(filteredMatches.end(), tempMatchContainer.begin(),
                                                   tempMatchContainer.end());
                        }
                    }
                    lastIn = false;
                    // Initialize range info
                    tempMatchContainer.clear();
                    diffPosCntOfCurrRange = 1;
                }
                i++;
            }

            // Met next species
            if (lastIn) {
                tempMatchContainer.push_back(matchList + i);
                if (diffPosCntOfCurrRange >= minCoveredPos) {
                    filteredMatches.insert(filteredMatches.end(), tempMatchContainer.begin(),
                                           tempMatchContainer.end());
                }
            }
            tempMatchContainer.clear();
            i++;
        }

        // Construct a match combination using filtered matches of current genus
        // so that it can best cover the query, and score the combination
        if (!filteredMatches.empty()) {
            genusScores.push_back(scoreGenus(filteredMatches, readLength1, readLength2));
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
        for (const Match * m : matchesForEachGenus[g]) {
            genusMatches.push_back(*m);
        }
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

TaxonScore Classifier::getBestGenusMatches(vector<Match> &genusMatches, const Match *matchList, size_t end,
                                           size_t offset, int queryLength, const LocalParameters & par) {
    TaxID currentGenus;
    TaxID currentSpecies;

    vector<const Match *> filteredMatches;
    vector<vector<const Match *>> matchesForEachGenus;
    vector<TaxonScore> genusScores;
    TaxonScore bestScore;
    size_t i = offset;
    uint8_t curFrame;
    vector<const Match *> curFrameMatches;
    while (i  < end + 1) {
        currentGenus = matchList[i].genusId;
        // For current genus
        while ((i < end + 1) && currentGenus == matchList[i].genusId) {
            currentSpecies = matchList[i].speciesId;

            // For current species
            while ((i < end + 1) && currentSpecies == matchList[i].speciesId) {
                curFrame = matchList[i].qInfo.frame;
                curFrameMatches.clear();

                // For current frame
                while ((i < end + 1) && currentSpecies == matchList[i].speciesId
                       && curFrame == matchList[i].qInfo.frame) {
                    curFrameMatches.push_back(&matchList[i]);
                    i ++;
                }
                if (curFrameMatches.size() > 1) {
                    remainConsecutiveMatches(curFrameMatches, filteredMatches, currentGenus, par);
                }
            }
        }

        // Construct a match combination using filtered matches of current genus
        // so that it can best cover the query, and score the combination
        if (!filteredMatches.empty()) {
            matchesForEachGenus.push_back(filteredMatches);
            genusScores.push_back(scoreGenus(filteredMatches, queryLength));
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
        for (const Match * m : matchesForEachGenus[g]) {
            genusMatches.push_back(*m);
        }
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

TaxonScore Classifier::getBestGenusMatches_spaced(vector<Match> &genusMatches, const Match *matchList, size_t end,
                                                  size_t offset, int readLength) {
    TaxID currentGenus;
    TaxID currentSpecies;

    vector<const Match *> tempMatchContainer;
    vector<const Match *> filteredMatches;
    vector<vector<Match>> matchesForEachGenus;
    vector<bool> conservedWithinGenus;
    vector<TaxonScore> genusScores;
    TaxonScore bestScore;
    size_t i = offset;
    bool lastIn;
    size_t speciesMatchCnt;
    while (i + 1 < end + 1) {
        currentGenus = matchList[i].genusId;
        // For current genus
        while ((i + 1 < end + 1) && currentGenus == matchList[i].genusId) {
            currentSpecies = matchList[i].speciesId;
            // For current species
            // Filter un-consecutive matches (probably random matches)
            lastIn = false;
            int distance = 0;
            int diffPosCntOfCurrRange = 1;
            int dnaDist = 0;

            // For the same species
            while ((i + 1 < end + 1) && currentSpecies == matchList[i + 1].speciesId) {
                distance = matchList[i + 1].qInfo.pos / 3 - matchList[i].qInfo.pos / 3;
                dnaDist = matchList[i + 1].qInfo.pos - matchList[i].qInfo.pos;
                if (distance == 0) { // At the same position
                    tempMatchContainer.push_back(matchList + i);
                } else if (dnaDist < (8 + spaceNum_int + maxGap) * 3) { // Overlapping
                    lastIn = true;
                    tempMatchContainer.push_back(matchList + i);
                    diffPosCntOfCurrRange++;
                } else { // Not consecutive --> End range
                    if (lastIn) {
                        tempMatchContainer.push_back(matchList + i);
                        if (diffPosCntOfCurrRange >= minCoveredPos) {
                            filteredMatches.insert(filteredMatches.end(), tempMatchContainer.begin(),
                                                   tempMatchContainer.end());
                        }
                    }
                    lastIn = false;
                    // Initialize range info
                    tempMatchContainer.clear();
                    diffPosCntOfCurrRange = 1;
                }
                i++;
            }

            // Met next species
            if (lastIn) {
                tempMatchContainer.push_back(matchList + i);
                if (diffPosCntOfCurrRange >= minCoveredPos) {
                    filteredMatches.insert(filteredMatches.end(), tempMatchContainer.begin(),
                                           tempMatchContainer.end());
                }
            }
            tempMatchContainer.clear();
            i++;
        }

        // Construct a match combination using filtered matches of current genus
        // so that it can best cover the query, and score the combination
        if (!filteredMatches.empty()) {
            genusScores.push_back(scoreGenus(filteredMatches, readLength));
        }
        filteredMatches.clear();
    }

    // If there are no meaningful genus
    if (genusScores.empty()) {
        bestScore.score = 0;
        return bestScore;
    }

    TaxonScore maxScore = *max_element(genusScores.begin(), genusScores.end(),
                                       [](const TaxonScore &a, const TaxonScore &b) { return a.score < b.score; });

    vector<size_t> maxIdx;
    for (size_t g = 0; g < genusScores.size(); g++) {
        if (genusScores[g].score > maxScore.score * 0.95f) {
            maxIdx.push_back(g);
        }
    }
    bestScore = maxScore;

    for (unsigned long g: maxIdx) {
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

TaxonScore Classifier::scoreGenus(vector<const Match *> &filteredMatches,
                                  int queryLength) {
    // Calculate Hamming distance & covered length
    int coveredPosCnt = 0;
    uint16_t currHammings;
    int aminoAcidNum = (int) queryLength / 3;
    int currPos;
    size_t matchNum = filteredMatches.size();
    size_t f = 0;

    // Get the largest hamming distance at each position of query
    auto *hammingsAtEachPos = new signed char[aminoAcidNum + 1];
    memset(hammingsAtEachPos, -1, (aminoAcidNum + 1));
    while (f < matchNum) {
        currPos = filteredMatches[f]->qInfo.pos / 3;
        currHammings = filteredMatches[f]->rightEndHamming;
        if (GET_2_BITS(currHammings) > hammingsAtEachPos[currPos + unmaskedPos[0]])
            hammingsAtEachPos[currPos + unmaskedPos[0]] = GET_2_BITS(currHammings);
        if (GET_2_BITS(currHammings >> 2) > hammingsAtEachPos[currPos + unmaskedPos[1]])
            hammingsAtEachPos[currPos + unmaskedPos[1]] = GET_2_BITS(currHammings >> 2);
        if (GET_2_BITS(currHammings >> 4) > hammingsAtEachPos[currPos + unmaskedPos[2]])
            hammingsAtEachPos[currPos + unmaskedPos[2]] = GET_2_BITS(currHammings >> 4);
        if (GET_2_BITS(currHammings >> 6) > hammingsAtEachPos[currPos + unmaskedPos[3]])
            hammingsAtEachPos[currPos + unmaskedPos[3]] = GET_2_BITS(currHammings >> 6);
        if (GET_2_BITS(currHammings >> 8) > hammingsAtEachPos[currPos + unmaskedPos[4]])
            hammingsAtEachPos[currPos + unmaskedPos[4]] = GET_2_BITS(currHammings >> 8);
        if (GET_2_BITS(currHammings >> 10) > hammingsAtEachPos[currPos + unmaskedPos[5]])
            hammingsAtEachPos[currPos + unmaskedPos[5]] = GET_2_BITS(currHammings >> 10);
        if (GET_2_BITS(currHammings >> 12) > hammingsAtEachPos[currPos + unmaskedPos[6]])
            hammingsAtEachPos[currPos + unmaskedPos[6]] = GET_2_BITS(currHammings >> 12);
        if (GET_2_BITS(currHammings >> 14) > hammingsAtEachPos[currPos + unmaskedPos[7]])
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

    return {filteredMatches[0]->genusId, score, coverage, (int) hammingSum};
}

TaxonScore Classifier::scoreGenus(vector<const Match *> &filteredMatches,
                                  int readLength1,
                                  int readLength2) {

    // Calculate Hamming distance & covered length
    uint16_t currHammings;
    int aminoAcidNum_total = ((int) readLength1 / 3) + ((int) readLength2 / 3);
    int aminoAcidNum_read1 = ((int) readLength1 / 3);
    int currPos;
    size_t matchNum = filteredMatches.size();
    size_t f = 0;

    // Get the largest hamming distance at each position of query
    auto *hammingsAtEachPos = new signed char[aminoAcidNum_total + 3];
    memset(hammingsAtEachPos, -1, (aminoAcidNum_total + 3));
    while (f < matchNum) {
        currPos = (int) filteredMatches[f]->qInfo.pos / 3;
        currHammings = filteredMatches[f]->rightEndHamming;
        if (GET_2_BITS(currHammings) > hammingsAtEachPos[currPos + unmaskedPos[0]])
            hammingsAtEachPos[currPos + unmaskedPos[0]] = GET_2_BITS(currHammings);
        if (GET_2_BITS(currHammings >> 2) > hammingsAtEachPos[currPos + unmaskedPos[1]])
            hammingsAtEachPos[currPos + unmaskedPos[1]] = GET_2_BITS(currHammings >> 2);
        if (GET_2_BITS(currHammings >> 4) > hammingsAtEachPos[currPos + unmaskedPos[2]])
            hammingsAtEachPos[currPos + unmaskedPos[2]] = GET_2_BITS(currHammings >> 4);
        if (GET_2_BITS(currHammings >> 6) > hammingsAtEachPos[currPos + unmaskedPos[3]])
            hammingsAtEachPos[currPos + unmaskedPos[3]] = GET_2_BITS(currHammings >> 6);
        if (GET_2_BITS(currHammings >> 8) > hammingsAtEachPos[currPos + unmaskedPos[4]])
            hammingsAtEachPos[currPos + unmaskedPos[4]] = GET_2_BITS(currHammings >> 8);
        if (GET_2_BITS(currHammings >> 10) > hammingsAtEachPos[currPos + unmaskedPos[5]])
            hammingsAtEachPos[currPos + unmaskedPos[5]] = GET_2_BITS(currHammings >> 10);
        if (GET_2_BITS(currHammings >> 12) > hammingsAtEachPos[currPos + unmaskedPos[6]])
            hammingsAtEachPos[currPos + unmaskedPos[6]] = GET_2_BITS(currHammings >> 12);
        if (GET_2_BITS(currHammings >> 14) > hammingsAtEachPos[currPos + unmaskedPos[7]])
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

//    matchesForEachGenus.push_back(move(filteredMatches));
    return {filteredMatches[0]->genusId, score, coverage, (int) hammingSum};
}

TaxonScore Classifier::chooseSpecies(const vector<Match> &matches,
                                     int queryLength,
                                     vector<TaxID> &species,
                                     unordered_map<TaxID, pair<int, int>> & speciesMatchRange) {
    // Score each species
    std::unordered_map<TaxID, TaxonScore> speciesScores;
    size_t i = 0;
    TaxID currentSpeices;
    size_t numOfMatch = matches.size();
    size_t speciesBegin, speciesEnd;
    while (i < numOfMatch) {
        currentSpeices = matches[i].speciesId;
        speciesBegin = i;
        while ((i < numOfMatch) && currentSpeices == matches[i].speciesId) {
            i++;
        }
        speciesEnd = i;
        speciesScores[currentSpeices] = scoreSpecies(matches, speciesBegin, speciesEnd, queryLength);
        speciesMatchRange[currentSpeices] = {(int) speciesBegin, (int) speciesEnd};
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
                                     vector<TaxID> &species,
                                     unordered_map<TaxID, pair<int, int>> & speciesMatchRange) {
    // Score each species
    std::unordered_map<TaxID, TaxonScore> speciesScores;


    size_t i = 0;
    TaxID currentSpeices;
    size_t numOfMatch = matches.size();
    size_t speciesBegin, speciesEnd;
    while (i < numOfMatch) {
        currentSpeices = matches[i].speciesId;
        speciesBegin = i;
        while ((i < numOfMatch) && currentSpeices == matches[i].speciesId) {
            i++;
        }
        speciesEnd = i;
        speciesScores[currentSpeices] = scoreSpecies(matches, speciesBegin, speciesEnd, read1Length, read2Length);
        speciesMatchRange[currentSpeices] = {(int) speciesBegin, (int) speciesEnd};
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

TaxonScore Classifier::scoreSpecies(const vector<Match> &matches,
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
        currPos = matches[walker].qInfo.pos / 3;
        currHammings = matches[walker].rightEndHamming;
        if (GET_2_BITS(currHammings) > hammingsAtEachPos[currPos + unmaskedPos[0]])
            hammingsAtEachPos[currPos + unmaskedPos[0]] = GET_2_BITS(currHammings);
        if (GET_2_BITS(currHammings >> 2) > hammingsAtEachPos[currPos + unmaskedPos[1]])
            hammingsAtEachPos[currPos + unmaskedPos[1]] = GET_2_BITS(currHammings >> 2);
        if (GET_2_BITS(currHammings >> 4) > hammingsAtEachPos[currPos + unmaskedPos[2]])
            hammingsAtEachPos[currPos + unmaskedPos[2]] = GET_2_BITS(currHammings >> 4);
        if (GET_2_BITS(currHammings >> 6) > hammingsAtEachPos[currPos + unmaskedPos[3]])
            hammingsAtEachPos[currPos + unmaskedPos[3]] = GET_2_BITS(currHammings >> 6);
        if (GET_2_BITS(currHammings >> 8) > hammingsAtEachPos[currPos + unmaskedPos[4]])
            hammingsAtEachPos[currPos + unmaskedPos[4]] = GET_2_BITS(currHammings >> 8);
        if (GET_2_BITS(currHammings >> 10) > hammingsAtEachPos[currPos + unmaskedPos[5]])
            hammingsAtEachPos[currPos + unmaskedPos[5]] = GET_2_BITS(currHammings >> 10);
        if (GET_2_BITS(currHammings >> 12) > hammingsAtEachPos[currPos + unmaskedPos[6]])
            hammingsAtEachPos[currPos + unmaskedPos[6]] = GET_2_BITS(currHammings >> 12);
        if (GET_2_BITS(currHammings >> 14) > hammingsAtEachPos[currPos + unmaskedPos[7]])
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

TaxonScore Classifier::scoreSpecies(const vector<Match> &matches,
                                  size_t begin,
                                  size_t end,
                                  int queryLength,
                                  int queryLength2) {

    // Get the smallest hamming distance at each position of query
    int aminoAcidNum_total = queryLength / 3 + queryLength2 / 3;
    int aminoAcidNum_read1 = queryLength / 3;
    auto *hammingsAtEachPos = new signed char[aminoAcidNum_total + 3];
    memset(hammingsAtEachPos, -1, (aminoAcidNum_total + 3));

    int currPos;
    size_t walker = begin;
    uint16_t currHammings;

    while (walker < end) {
        currPos = matches[walker].qInfo.pos / 3;
        currHammings = matches[walker].rightEndHamming;
        if (GET_2_BITS(currHammings) > hammingsAtEachPos[currPos + unmaskedPos[0]])
            hammingsAtEachPos[currPos + unmaskedPos[0]] = GET_2_BITS(currHammings);
        if (GET_2_BITS(currHammings >> 2) > hammingsAtEachPos[currPos + unmaskedPos[1]])
            hammingsAtEachPos[currPos + unmaskedPos[1]] = GET_2_BITS(currHammings >> 2);
        if (GET_2_BITS(currHammings >> 4) > hammingsAtEachPos[currPos + unmaskedPos[2]])
            hammingsAtEachPos[currPos + unmaskedPos[2]] = GET_2_BITS(currHammings >> 4);
        if (GET_2_BITS(currHammings >> 6) > hammingsAtEachPos[currPos + unmaskedPos[3]])
            hammingsAtEachPos[currPos + unmaskedPos[3]] = GET_2_BITS(currHammings >> 6);
        if (GET_2_BITS(currHammings >> 8) > hammingsAtEachPos[currPos + unmaskedPos[4]])
            hammingsAtEachPos[currPos + unmaskedPos[4]] = GET_2_BITS(currHammings >> 8);
        if (GET_2_BITS(currHammings >> 10) > hammingsAtEachPos[currPos + unmaskedPos[5]])
            hammingsAtEachPos[currPos + unmaskedPos[5]] = GET_2_BITS(currHammings >> 10);
        if (GET_2_BITS(currHammings >> 12) > hammingsAtEachPos[currPos + unmaskedPos[6]])
            hammingsAtEachPos[currPos + unmaskedPos[6]] = GET_2_BITS(currHammings >> 12);
        if (GET_2_BITS(currHammings >> 14) > hammingsAtEachPos[currPos + unmaskedPos[7]])
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

void Classifier::writeReadClassification(const vector<Query> & queryList, int queryNum, ofstream &readClassificationFile) {
    for (int i = 0; i < queryNum; i++) {
        readClassificationFile << queryList[i].isClassified << "\t" << queryList[i].name << "\t"
                               << queryList[i].classification << "\t"
                               << queryList[i].queryLength + queryList[i].queryLength2 << "\t"
                               << queryList[i].score << "\t"
                               << queryList[i].coverage << "\t"
                               << queryList[i].hammingDist << "\t"
                               << taxonomy->getString(taxonomy->taxonNode(queryList[i].classification)->rankIdx) << "\t";
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

void Classifier::writeReport(FILE *FP, const std::unordered_map<TaxID, TaxonCounts> &cladeCounts,
                             unsigned long totalReads, TaxID taxID, int depth) {
    std::unordered_map<TaxID, TaxonCounts>::const_iterator it = cladeCounts.find(taxID);
    unsigned int cladeCount = it == cladeCounts.end() ? 0 : it->second.cladeCount;
    unsigned int taxCount = it == cladeCounts.end() ? 0 : it->second.taxCount;
    if (taxID == 0) {
        if (cladeCount > 0) {
            fprintf(FP, "%.4f\t%i\t%i\tno rank\t0\tunclassified\n",
                    100 * cladeCount / double(totalReads),
                    cladeCount, taxCount);
        }
        writeReport(FP, cladeCounts, totalReads, 1);
    } else {
        if (cladeCount == 0) {
            return;
        }
        const TaxonNode *taxon = taxonomy->taxonNode(taxID);
        fprintf(FP, "%.4f\t%i\t%i\t%s\t%i\t%s%s\n",
                100 * cladeCount / double(totalReads), cladeCount, taxCount,
                taxonomy->getString(taxon->rankIdx), taxID, std::string(2 * depth, ' ').c_str(), taxonomy->getString(taxon->nameIdx));
        std::vector<TaxID> children = it->second.children;
        SORT_SERIAL(children.begin(), children.end(), [&](int a, int b) { return cladeCountVal(cladeCounts, a) > cladeCountVal(cladeCounts, b); });
        for (size_t i = 0; i < children.size(); ++i) {
            TaxID childTaxId = children[i];
            if (cladeCounts.count(childTaxId)) {
                writeReport(FP, cladeCounts, totalReads, childTaxId, depth + 1);
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
void Classifier::splitFASTQ(vector<SequenceBlock> & seqSegments, const string & queryPath) {
    ifstream fastq;
    fastq.open(queryPath);
    if (!fastq.is_open()) {
        cerr << "Error: Cannot open file " << queryPath << endl;
        exit(1);
    }

    string line;
    size_t lineCnt = 0;
    size_t start = 0;
    size_t end;
    while (getline(fastq, line)) {
        if (lineCnt % 4 == 0){
            start = (size_t) fastq.tellg() - line.length() - 1;
        }
        if (lineCnt % 4 == 1){
            end = (size_t) fastq.tellg() - 1;
            seqSegments.emplace_back(start, end, end - start + 1, line.length());
        }
        lineCnt++;
    }
}

void Classifier::splitFASTA(vector<SequenceBlock> & seqSegments, const string & queryPath) {
    ifstream fasta;
    fasta.open(queryPath);
    if (!fasta.is_open()) {
        cerr << "Error: Cannot open file " << queryPath << endl;
        exit(1);
    }

    string line;
    size_t start = 0;
    size_t end;
    getline(fasta, line);
    size_t seqLength = 0;
    while (getline(fasta, line)) {
        if (line[0] == '>') {
            seqSegments.emplace_back(start, end, end - start + 1, seqLength);
            start = end + 1;
            seqLength = 0;
        } else {
            end = (size_t) fasta.tellg() - 1;
            seqLength += line.length();
        }
    }
    seqSegments.emplace_back(start, end, end - start + 1, seqLength);
}


bool Classifier::isConsecutive(const Match * match1, const Match * match2) {
//    uint16_t hamming1 = match1->rightEndHamming;
//    uint16_t hamming2 = match2->rightEndHamming;
//    // move bits to right by 2
//    hamming1 >>= 2;
//    // set most significant two bits to 0
//    hamming2 &= 0x3FFF;
    return (match1->rightEndHamming >> 2) == (match2->rightEndHamming & 0x3FFF);
}

bool Classifier::isConsecutive(const Match & match1, const Match & match2, const LocalParameters & par) {
    uint16_t hamming1 = match1.rightEndHamming;
    uint16_t hamming2 = match2.rightEndHamming;
    if (par.printLog) {
        print_binary16(16, hamming1); cout << endl;
        print_binary16(16, hamming2); cout << endl;
    }

    // set most significant two bits to 0
    hamming2 &= 0x3FFF; // 07654321
    // move bits to right by 2
    hamming1 >>= 2; // 07654321
    if (par.printLog) {
        print_binary16(16, hamming1); cout << endl;
        print_binary16(16, hamming2); cout << endl;
    }

    return hamming1 == hamming2;
}

