#include "QueryFilter.h"
#include "common.h"
#include <unordered_map>

QueryFilter::QueryFilter(LocalParameters & par) {
    // Load parameters
    dbDir = par.filenames[1 + (par.seqMode == 2)];
    if(FileUtil::fileExists(string(dbDir + "/diffIdx").c_str())) {
        isNewDB = false;
    } else {
        isNewDB = true; 
    }

    matchPerKmer = par.matchPerKmer;
    printMode = par.printMode;
    seqMode = par.seqMode;
    contams = Util::split(par.contamList, ",");
    loadDbParameters(par, dbDir);
    
    // Taxonomy
    taxonomy = loadTaxonomy(dbDir, par.taxonomyPath);
 
    // Agents
    geneticCode = new GeneticCode(par.reducedAA == 1);
    queryIndexer = new QueryIndexer(par);
    kmerExtractor = new KmerExtractor(par, *geneticCode, isNewDB);
    if (par.reducedAA) { kmerMatcher = new ReducedKmerMatcher(par, taxonomy, isNewDB);} 
    else { kmerMatcher = new KmerMatcher(par, taxonomy, isNewDB);}
    taxonomer = new Taxonomer(par, taxonomy, isNewDB);
    reporter = new Reporter(par, taxonomy);
    setInputAndOutputFiles(par);
    reporter->setReadClassificationFileName(readClassificationFileName);
    reporter->setReportFileName(reportFileName);
    cout << "Filtered reads: " << f1 << endl;
    if (par.seqMode == 2) { cout << "Filtered reads: " << f2 << endl; }
    if (printMode == 2) {
        cout << "Removed reads: " << rm1 << endl;
        if (par.seqMode == 2) { cout << "Removed reads: " << rm2 << endl; }
    }

    filter_kseq1 = KSeqFactory(in1.c_str());
    if (par.seqMode == 2) { filter_kseq2 = KSeqFactory(in2.c_str()); }

    readCounter = 0;

    // Open output files
    f1_fp = fopen(f1.c_str(), "w");
    if (par.seqMode == 2) { f2_fp = fopen(f2.c_str(), "w"); }
    if (printMode == 2) {
        rm1_fp = fopen(rm1.c_str(), "w");
        if (par.seqMode == 2) { rm2_fp = fopen(rm2.c_str(), "w"); }
    }


}

QueryFilter::~QueryFilter() {
    delete taxonomy;
    delete queryIndexer;
    delete kmerExtractor;
    delete kmerMatcher;
    delete taxonomer;
    delete reporter;
    delete filter_kseq1;
    delete filter_kseq2;
    delete[] isFiltered;
    fclose(f1_fp);
    if (seqMode == 2) { fclose(f2_fp); }
    if (printMode == 2) {
        fclose(rm1_fp);
        if (seqMode == 2) { fclose(rm2_fp); }
    }
}

void QueryFilter::setInputAndOutputFiles(const LocalParameters & par) {
    cout << "Setting output file names" << endl;
    // Get the base name of in1
    in1 = par.filenames[0];
    string baseName = LocalUtil::getQueryBaseName(in1);

    // Set the output file names
    f1 = baseName + "_filtered.fna";
    rm1 = baseName + "_removed.fna";
    reportFileName = baseName + "_report.tsv";
    readClassificationFileName = baseName + "_classifications.tsv";

    // For paired-end reads
    if (seqMode == 2) {
        in2 = par.filenames[1];
        f2 = LocalUtil::getQueryBaseName(in2) + "_filtered.fna";
        rm2 = LocalUtil::getQueryBaseName(in2) + "_removed.fna";
    }
}

void QueryFilter::recordFilteredReads(const vector<Query> & queryList) {
    for (auto query : queryList) {
        isFiltered[readCounter++] = query.isClassified;
    }
}

void QueryFilter::printFilteredReads() {
    for (size_t i = 0; i < readCounter; i ++) {
        // Read query reads
        filter_kseq1->ReadEntry();
        if (seqMode == 2) { filter_kseq2->ReadEntry(); }

        // Print reads
        if (!isFiltered[i]) { // Print filtered reads
            fprintf(f1_fp, ">%s\n%s\n", filter_kseq1->entry.name.s, filter_kseq1->entry.sequence.s);
            if (seqMode == 2) { fprintf(f2_fp, ">%s\n%s\n", filter_kseq2->entry.name.s, filter_kseq2->entry.sequence.s); }
        } else if (printMode == 2) { // Print removed reads
            fprintf(rm1_fp, ">%s\n%s\n", filter_kseq1->entry.name.s, filter_kseq1->entry.sequence.s);
            if (seqMode == 2) { fprintf(rm2_fp, ">%s\n%s\n", filter_kseq2->entry.name.s, filter_kseq2->entry.sequence.s); }
        }
    }
}

void QueryFilter::filterReads(LocalParameters & par) {

    cout << "Indexing query file ...";
    queryIndexer->indexQueryFile(0);
    size_t numOfSeq = queryIndexer->getReadNum_1();
    size_t totalReadLength = queryIndexer->getTotalReadLength();
    const vector<QuerySplit> & queryReadSplit = queryIndexer->getQuerySplits();
    // print queryReadSplit
    // for (size_t i = 0; i < queryReadSplit.size(); i++) {
    //     cout << queryReadSplit[i].start << " " << queryReadSplit[i].end << " " << queryReadSplit[i].kmerCnt << endl;
    // }
    cout << "Done" << endl;
    cout << "Total number of sequences: " << numOfSeq << endl;
    cout << "Total read length: " << totalReadLength <<  "nt" << endl;

    isFiltered = new bool[queryIndexer->getReadNum_1()];
    memset(isFiltered, 0, sizeof(bool) * queryIndexer->getReadNum_1());
    Buffer<Kmer> kmerBuffer;
    Buffer<Match> matchBuffer;
    vector<Query> queryList;

    size_t numOfTatalQueryKmerCnt = 0;
    size_t processedSeqCnt = 0;
    reporter->openReadClassificationFile();

#ifdef OPENMP
    omp_set_num_threads(par.threads);
#endif

    KSeqWrapper* kseq1 = KSeqFactory(in1.c_str());
    KSeqWrapper* kseq2 = nullptr;
    if (par.seqMode == 2) { kseq2 = KSeqFactory(in2.c_str()); }

    for (size_t splitIdx = 0; splitIdx < queryReadSplit.size(); splitIdx++) {
        // Allocate memory for query list
        queryList.clear();
        queryList.resize(queryReadSplit[splitIdx].end - queryReadSplit[splitIdx].start);

        // Allocate memory for query k-mer list and match list
        kmerBuffer.reallocateMemory(queryReadSplit[splitIdx].kmerCnt);
        if (queryReadSplit.size() == 1) {
            size_t remain = queryIndexer->getAvailableRam() - queryReadSplit[splitIdx].kmerCnt * sizeof(Kmer) - numOfSeq * 200;
            matchBuffer.reallocateMemory(remain / sizeof(Match));
        } else {
            matchBuffer.reallocateMemory(queryReadSplit[splitIdx].kmerCnt * matchPerKmer);
        }

        // Initialize query k-mer buffer and match buffer
        kmerBuffer.startIndexOfReserve = 0;
        matchBuffer.startIndexOfReserve = 0;

        // Extract query k-mer
        kmerExtractor->extractQueryKmers(kmerBuffer,
                                         queryList,
                                         queryReadSplit[splitIdx],
                                         par,
                                         kseq1,
                                         kseq2);
        numOfTatalQueryKmerCnt += kmerBuffer.startIndexOfReserve;

        // new code
        // Search matches between query and target k-mers
        for (auto db : contams) {
            cout <<"";
        }
        kmerMatcher->sortMatches(&matchBuffer);

        // Classify queries based on the matches
        taxonomer->assignTaxonomy(matchBuffer.buffer, matchBuffer.startIndexOfReserve, queryList, par);
        processedSeqCnt += queryReadSplit[splitIdx].end - queryReadSplit[splitIdx].start;
        cout << "The number of processed sequences: " << processedSeqCnt << " (" << (double) processedSeqCnt / (double) numOfSeq << ")" << endl;
        
        // Write classification results
        reporter->writeReadClassification(queryList, true);

        recordFilteredReads(queryList);
    }

    cout << "Number of query k-mers: " << numOfTatalQueryKmerCnt << endl;
    cout << "The number of matches: " << kmerMatcher->getTotalMatchCnt() << endl;
    printFilteredReads();
    reporter->writeReportFile(numOfSeq, taxonomer->getTaxCounts(), ReportType::Default);
    reporter->closeReadClassificationFile();

    // Memory deallocation
    free(matchBuffer.buffer);
    delete kseq1;
    delete kseq2;
}

