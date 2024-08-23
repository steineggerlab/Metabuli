#include "Classifier.h"
#include "FileUtil.h"
#include "QueryIndexer.h"
#include "common.h"

Classifier::Classifier(LocalParameters & par) {
    // Load parameters
    dbDir = par.filenames[1 + (par.seqMode == 2)];
    matchPerKmer = par.matchPerKmer;
    loadDbParameters(par);

    cout << "DB name: " << par.dbName << endl;
    cout << "DB creation date: " << par.dbDate << endl;
    
    // Taxonomy
    taxonomy = loadTaxonomy(dbDir, par.taxonomyPath);

    // Agents
    queryIndexer = new QueryIndexer(par);
    kmerExtractor = new KmerExtractor(par);
    if (par.reducedAA) {
        kmerMatcher = new ReducedKmerMatcher(par, taxonomy);
    } else {
        kmerMatcher = new KmerMatcher(par, taxonomy);
    }
    // taxonomer = new Taxonomer(par, taxonomy);
    reporter = new Reporter(par, taxonomy);
}

Classifier::~Classifier() {
    delete taxonomy;
    delete queryIndexer;
    delete kmerExtractor;
    delete kmerMatcher;
    // delete taxonomer;
    delete reporter;
}

void Classifier::startClassify(const LocalParameters &par) {
    QueryKmerBuffer queryKmerBuffer;
    Buffer<Match> matchBuffer;
    vector<Query> queryList;
    size_t numOfTatalQueryKmerCnt = 0;
    reporter->openReadClassificationFile();

    bool complete = false;
    size_t processedReadCnt = 0;
    size_t tries = 0;
    size_t totalSeqCnt = 0;
    
    // Extract k-mers from query sequences and compare them to target k-mer DB
    while (!complete) {
        tries++;

        // Get splits for remaining sequences
        if (tries == 1) {
                cout << "Indexing query file ...";
        }
        queryIndexer->setBytesPerKmer(matchPerKmer);
        queryIndexer->indexQueryFile(processedReadCnt);
        const vector<QuerySplit> & queryReadSplit = queryIndexer->getQuerySplits();

        if (tries == 1) {
            totalSeqCnt = queryIndexer->getReadNum_1();
            cout << "Done" << endl;
            cout << "Total number of sequences: " << queryIndexer->getReadNum_1() << endl;
            cout << "Total read length: " << queryIndexer->getTotalReadLength() <<  "nt" << endl;
        }

        // Set up kseq
        KSeqWrapper* kseq1 = KSeqFactory(par.filenames[0].c_str());
        KSeqWrapper* kseq2 = nullptr;
        if (par.seqMode == 2) { kseq2 = KSeqFactory(par.filenames[1].c_str()); }

        // Move kseq to unprocessed reads
        for (size_t i = 0; i < processedReadCnt; i++) {
            kseq1->ReadEntry();
            if (par.seqMode == 2) { kseq2->ReadEntry(); }
        }

        for (size_t splitIdx = 0; splitIdx < queryReadSplit.size(); splitIdx++) {
            // Allocate memory for query list
            queryList.clear();
            queryList.resize(queryReadSplit[splitIdx].end - queryReadSplit[splitIdx].start);

            // Allocate memory for query k-mer buffer
            queryKmerBuffer.reallocateMemory(queryReadSplit[splitIdx].kmerCnt);

            // Allocate memory for match buffer
            if (queryReadSplit.size() == 1) {
                size_t remain = queryIndexer->getAvailableRam() 
                                - (queryReadSplit[splitIdx].kmerCnt * sizeof(QueryKmer)) 
                                - (queryIndexer->getReadNum_1() * 200); // TODO: check it later
                matchBuffer.reallocateMemory(remain / sizeof(Match));
            } else {
                matchBuffer.reallocateMemory(queryReadSplit[splitIdx].kmerCnt * matchPerKmer);
            }

            // Initialize query k-mer buffer and match buffer
            queryKmerBuffer.startIndexOfReserve = 0;
            matchBuffer.startIndexOfReserve = 0;

            // Extract query k-mers
            kmerExtractor->extractQueryKmers(queryKmerBuffer,
                                             queryList,
                                             queryReadSplit[splitIdx],
                                             par,
                                             kseq1,
                                             kseq2); // sync kseq1 and kseq2
            
            // Search matches between query and target k-mers
            if (kmerMatcher->matchKmers(&queryKmerBuffer, &matchBuffer)) {
                kmerMatcher->sortMatches(&matchBuffer);
                
                // Classify queries based on the matches.
                assignTaxonomy(matchBuffer.buffer, matchBuffer.startIndexOfReserve, queryList, par);

                // Write classification results
                reporter->writeReadClassification(queryList);

                // Print progress
                processedReadCnt += queryReadSplit[splitIdx].readCnt;
                cout << "The number of processed sequences: " << processedReadCnt << " (" << (double) processedReadCnt / (double) totalSeqCnt << ")" << endl;

                numOfTatalQueryKmerCnt += queryKmerBuffer.startIndexOfReserve;
            } else { // search was incomplete
                // Increase matchPerKmer and try again 
                matchPerKmer += 4;
                // delete kseq1;
                // delete kseq2;
                cout << "--match-per-kmer was increased to " << matchPerKmer << " and searching again..." << endl;
                // cout << "The search was incomplete. Increasing --match-per-kmer to " << matchPerKmer << " and trying again..." << endl;
                break;
            }
         }
         
        delete kseq1;
        if (par.seqMode == 2) {
            delete kseq2;
        }
        if (processedReadCnt == totalSeqCnt) {
            complete = true;
        }
    }

    cout << "Number of query k-mers: " << numOfTatalQueryKmerCnt << endl;
    cout << "The number of matches: " << kmerMatcher->getTotalMatchCnt() << endl;
    reporter->closeReadClassificationFile();

    // Write report files
    reporter->writeReportFile(totalSeqCnt, taxCounts);

    // Memory deallocation
    free(matchBuffer.buffer);
}

void Classifier::assignTaxonomy(const Match *matchList,
                               size_t numOfMatches,
                               std::vector<Query> &queryList,
                               const LocalParameters &par) {
    time_t beforeAnalyze = time(nullptr);
    cout << "Analyzing matches ..." << endl;

    // Divide matches into blocks for multi threading
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
    cout << "Time spent for spliting matches: " << double(time(nullptr) - beforeAnalyze) << endl;
    // Process each block
#pragma omp parallel default(none), shared(cout, matchBlocks, matchList, seqNum, queryList, blockIdx, par)
    {
        Taxonomer taxonomer(par, taxonomy);
#pragma omp for schedule(dynamic, 1)
        for (size_t i = 0; i < blockIdx; ++i) {
            taxonomer.chooseBestTaxon(matchBlocks[i].id,
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
    cout << "Time spent for analyzing: " << double(time(nullptr) - beforeAnalyze) << endl;

}