#include "QueryFilter.h"

QueryFilter::QueryFilter(LocalParameters & par) {
    queryIndexer = new QueryIndexer(par);

    setInputAndOutputFiles(par);
}

QueryFilter::~QueryFilter() {
    delete queryIndexer;
}

void QueryFilter::setInputAndOutputFiles(const LocalParameters & par) {
    // Get the base name of in1
    in1 = par.filenames[0];
    string baseName = LocalUtil::getQueryBaseName(in1);

    // Set the output file names
    out1 = baseName + "_filtered.fna.gz";
    reportFileName = baseName + "_filter_report.tsv";

    // For paired-end reads
    if (par.seqMode == 2) {
        in2 = par.filenames[1];
        out2 = LocalUtil::getQueryBaseName(in2) + "_filtered.fna.gz";
    }
}

void QueryFilter::filterReads(LocalParameters & par) {

    cout << "Indexing query file ...";
    queryIndexer->indexQueryFile();
    size_t numOfSeq = queryIndexer->getReadNum_1();
    size_t totalReadLength = queryIndexer->getTotalReadLength();
    const vector<QuerySplit> & queryReadSplit = queryIndexer->getQuerySplits();
    cout << "Done" << endl;
    cout << "Total number of sequences: " << numOfSeq << endl;
    cout << "Total read length: " << totalReadLength <<  "nt" << endl;

    QueryKmerBuffer kmerBuffer;
    Buffer<Match> matchBuffer;
    vector<Query> queryList;


}

