#include "Classifier.h"
#include "FileUtil.h"
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
    taxonomer = new Taxonomer(par, taxonomy);
    reporter = new Reporter(par, taxonomy);
}

Classifier::~Classifier() {
    delete taxonomy;
    delete queryIndexer;
    delete kmerExtractor;
    delete kmerMatcher;
    delete taxonomer;
    delete reporter;
}

void Classifier::startClassify(const LocalParameters &par) {

    cout << "Indexing query file ...";
   
    size_t numOfSeq = queryIndexer->getReadNum_1();
    size_t totalReadLength = queryIndexer->getTotalReadLength();
    const vector<QuerySplit> & queryReadSplit = queryIndexer->getQuerySplits();
    cout << "Done" << endl;
    cout << "Total number of sequences: " << numOfSeq << endl;
    cout << "Total read length: " << totalReadLength <<  "nt" << endl;

    QueryKmerBuffer queryKmerBuffer;
    Buffer<Match> matchBuffer;
    vector<Query> queryList;

    size_t numOfTatalQueryKmerCnt = 0;

    reporter->openReadClassificationFile();
#ifdef OPENMP
    omp_set_num_threads(par.threads);
#endif
    // Extract k-mers from query sequences and compare them to target k-mer DB
    
    
    bool complete = false;
    size_t processedReadCnt = 0;

    while (!complete) {
        // Get splits for remaining sequences
        queryIndexer->setBytesPerKmer(matchPerKmer);
        queryIndexer->indexQueryFile(processedReadCnt);

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
                size_t remain = queryIndexer->getAvailableRam() - queryReadSplit[splitIdx].kmerCnt * sizeof(QueryKmer) - numOfSeq * 200; // TODO: check it later
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
            numOfTatalQueryKmerCnt += queryKmerBuffer.startIndexOfReserve;

            // Search matches between query and target k-mers
            if (kmerMatcher->matchKmers(&queryKmerBuffer, &matchBuffer)) {
                kmerMatcher->sortMatches(&matchBuffer);
                
                // Classify queries based on the matches.
                taxonomer->assignTaxonomy(matchBuffer.buffer, matchBuffer.startIndexOfReserve, queryList, par);
                processedReadCnt += queryReadSplit[splitIdx].readCnt;
                // processedSeqCnt += queryReadSplit[splitIdx].end - queryReadSplit[splitIdx].start;
                cout << "The number of processed sequences: " << processedReadCnt << " (" << (double) processedReadCnt / (double) numOfSeq << ")" << endl;

                // Write classification results
                reporter->writeReadClassification(queryList);
            } else { // search was incomplete
                matchPerKmer *= 2;
                delete kseq1;
                delete kseq2;
                break;
            }
        }
        delete kseq1;
        delete kseq2;
        complete = true;
    }

    cout << "Number of query k-mers: " << numOfTatalQueryKmerCnt << endl;
    cout << "The number of matches: " << kmerMatcher->getTotalMatchCnt() << endl;
    reporter->closeReadClassificationFile();

    // Write report files
    reporter->writeReportFile(numOfSeq, taxonomer->getTaxCounts());

    // Memory deallocation
    free(matchBuffer.buffer);
}
