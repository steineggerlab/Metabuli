#include "Classifier.h"
#include "FileUtil.h"
#include "QueryIndexer.h"
#include "common.h"

Classifier::Classifier(LocalParameters & par) {
    dbDir = par.filenames[1 + (par.seqMode == 2)];
    // if(FileUtil::fileExists(string(dbDir + "/diffIdx").c_str())) {
    //     isNewDB = false;
    // } else {
    //     isNewDB = true; 
    // }

    matchPerKmer = par.matchPerKmer;
    loadDbParameters(par, par.filenames[1 + (par.seqMode == 2)]);
    kmerFormat = par.kmerFormat;

    cout << "Database name : " << par.dbName << endl;
    cout << "Creation date : " << par.dbDate << endl;
    
    taxonomy = loadTaxonomy(dbDir, par.taxonomyPath);

    geneticCode = new GeneticCode(par.reducedAA == 1);
    queryIndexer = new QueryIndexer(par);
    kmerExtractor = new KmerExtractor(par, *geneticCode, kmerFormat);
    if (par.reducedAA) {
        kmerMatcher = new ReducedKmerMatcher(par, taxonomy, kmerFormat);
    } else {
        kmerMatcher = new KmerMatcher(par, taxonomy, kmerFormat);
    }
    reporter = new Reporter(par, taxonomy);
}

Classifier::~Classifier() {
    delete taxonomy;
    delete queryIndexer;
    delete kmerExtractor;
    delete kmerMatcher;
    delete reporter;
    delete geneticCode;
    if (mappingResList) delete[] mappingResList;
}

void Classifier::startClassify(const LocalParameters &par) {
    Buffer<Kmer> queryKmerBuffer;
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
        queryIndexer->setBytesPerKmer(matchPerKmer);
        queryIndexer->indexQueryFile(processedReadCnt);
        const vector<QuerySplit> & queryReadSplit = queryIndexer->getQuerySplits();
        if (tries == 1) {
            totalSeqCnt = queryIndexer->getReadNum_1();
            cout << "--------------------" << endl;
            cout << "Total read count : " << queryIndexer->getReadNum_1() << endl;
            cout << "Total read length: " << queryIndexer->getTotalReadLength() <<  "nt" << endl;
            cout << "--------------------" << endl;
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
            queryKmerBuffer.init();
            // memset(queryKmerBuffer.buffer, 0, queryReadSplit[splitIdx].kmerCnt * sizeof(Kmer));

            // Allocate memory for match buffer
            if (queryReadSplit.size() == 1) {
                size_t remain = queryIndexer->getAvailableRam() 
                                - (queryReadSplit[splitIdx].kmerCnt * sizeof(Kmer)) 
                                - (queryIndexer->getReadNum_1() * 200); // TODO: check it later
                matchBuffer.reallocateMemory(remain / sizeof(Match));
            } else {
                matchBuffer.reallocateMemory(queryReadSplit[splitIdx].kmerCnt * matchPerKmer);
            }

            // Initialize query k-mer buffer and match buffer
            matchBuffer.startIndexOfReserve = 0;

            // Extract query k-mers
            kmerExtractor->extractQueryKmers(queryKmerBuffer,
                                             queryList,
                                             queryReadSplit[splitIdx],
                                             par,
                                             kseq1,
                                             kseq2);
            
            // Search matches between query and target k-mers
            bool searchComplete = false;
            searchComplete = kmerMatcher->matchKmers(&queryKmerBuffer, &matchBuffer);
            if (searchComplete) {
                cout << "K-mer match count      : " << kmerMatcher->getTotalMatchCnt() << endl;
                kmerMatcher->sortMatches(&matchBuffer);
                assignTaxonomy(matchBuffer.buffer, matchBuffer.startIndexOfReserve, queryList, par);
                reporter->writeReadClassification(queryList);
                if (par.em) {
                    reporter->writeMappings(queryList, processedReadCnt);
                    getTopSpecies(queryList);
                }
                processedReadCnt += queryReadSplit[splitIdx].readCnt;
                cout << "Processed read count   : " << processedReadCnt << " (" << (double) processedReadCnt / (double) totalSeqCnt << ")" << endl;
                // numOfTatalQueryKmerCnt += queryKmerBuffer.startIndexOfReserve;
            } else { // search was incomplete
                matchPerKmer += 4;
                cout << "--match-per-kmer was increased to " << matchPerKmer << " and searching again..." << endl;
                break;
            }
            cout << "--------------------" << endl;
        }
         
        delete kseq1;
        if (par.seqMode == 2) {
            delete kseq2;
        }
        if (processedReadCnt == totalSeqCnt) {
            complete = true;
        }

        cout << "--------------------" << endl;
    }

    // Finalize original classification results
    cout << "Total k-mer match count: " << kmerMatcher->getTotalMatchCnt() << endl;    
    reporter->closeReadClassificationFile();
    reporter->writeReportFile(totalSeqCnt, taxCounts, ReportType::Default);
    cout << "Taxonomic classification completed." << endl;

    // Perform EM algorithm    
    if (par.em) {
        cout << "Running EM algorithm for taxonomic re-assignment..." << endl;
        reporter->freeMappingWriteBuffer();        
        loadOriginalResults(reporter->getClassificationFileName(), totalSeqCnt);   
        em(totalSeqCnt);
        reporter->writeReclassifyResults(emResults);
        reporter->writeReportFile(totalSeqCnt, emTaxCounts, ReportType::EM);
        reporter->writeReportFile(totalSeqCnt, reclassifyTaxCounts, ReportType::EM_RECLASSIFY);
    }
    
    return;
}

void Classifier::assignTaxonomy(const Match *matchList,
                               size_t numOfMatches,
                               std::vector<Query> &queryList,
                               const LocalParameters &par) {
    time_t beforeAnalyze = time(nullptr);
    std::cout << "K-mer match analysis   : " << std::flush;
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
    // Process each block
#pragma omp parallel default(none), shared(cout, matchBlocks, matchList, seqNum, queryList, blockIdx, par)
    {
        Taxonomer taxonomer(par, taxonomy, kmerFormat);
        #pragma omp for schedule(dynamic, 1)
        for (size_t i = 0; i < blockIdx; ++i) {
            taxonomer.chooseBestTaxon(matchBlocks[i].id - 1,
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
    cout << double(time(nullptr) - beforeAnalyze) << " s" << endl;

}

void Classifier::em(size_t totalQueryCnt) {
    loadMappings(reporter->getMappingFileName());
    int maxTaxId = taxonomy->getMaxTaxID();
    vector<uint32_t> sp2uniqKmerCnt(maxTaxId + 1, 0);
    countUniqueKmerPerSpecies(sp2uniqKmerCnt);
    vector<double> sp2lengthFactor(maxTaxId + 1, 0.0);
    for (size_t i = 0; i < sp2uniqKmerCnt.size(); i++) {
        if (sp2uniqKmerCnt[i] > 0) {
            sp2lengthFactor[i] = 1.0 / log(sp2uniqKmerCnt[i]);
        } else {
            sp2lengthFactor[i] = 0.0;
        }
    }
    
    time_t st = time(nullptr);
    cout << "EM algorithm for taxonomic assignment: " << std::endl;
    size_t idx = 0;
    std::vector<std::pair<size_t, size_t>> queryRanges;
    while (idx < mappingResListSize) {
        TaxID currentQuery = mappingResList[idx].queryId;
        size_t startIdx = idx;
        while (idx < mappingResListSize && mappingResList[idx].queryId == currentQuery) {
            idx++;
        }
        queryRanges.emplace_back(startIdx, idx);
    }

    cout << "Species count: " << topSpeciesSet.size() << std::endl;
    std::unordered_map<TaxID, double> Fnew;
    std::vector<TaxID> spList;
    for (const auto & sp : topSpeciesSet) {
        taxProbs[sp] = 1.0 / topSpeciesSet.size();
        Fnew[sp] = 0.0;
        spList.push_back(sp);
    }

    std::atomic<int> converged{0};
    size_t queryCount = 0;
    #pragma omp parallel default(none) shared(queryRanges, queryCount, converged, mappingResList, taxProbs, Fnew, sp2lengthFactor, spList, std::cout)
    {
    std::unordered_map<TaxID, double> localFnew;
    size_t localQueryCount = 0;

    for (size_t iter = 0; iter < 1000; ++iter) {
        if (converged.load(std::memory_order_acquire))
            continue;

        // Init local variables
        for (auto & sp : localFnew) { sp.second = 0.0; }
        localQueryCount = 0;

        // Init shared variables
        #pragma omp single
        {
            for (auto & sp : Fnew) { sp.second = 0.0; }
            queryCount = 0;
            cout << "Iteration " << iter + 1 << ": " << std::endl;
        }

        // Calculate new probabilities for each species
        #pragma omp for schedule(dynamic, 1)
        for (size_t i = 0; i < queryRanges.size(); ++i) {
            double denom = 0.0;
            for (size_t j = queryRanges[i].first; j < queryRanges[i].second; ++j) {
                denom += mappingResList[j].score * taxProbs[mappingResList[j].speciesId] * sp2lengthFactor[mappingResList[j].speciesId];
            }
            if (denom == 0.0) { continue; }
            localQueryCount++;
            for (size_t j = queryRanges[i].first; j < queryRanges[i].second; ++j) {
                localFnew[mappingResList[j].speciesId] += (mappingResList[j].score * taxProbs[mappingResList[j].speciesId] * sp2lengthFactor[mappingResList[j].speciesId]) / denom;
            }
        }
    
        // Accumulate results from all threads
        #pragma omp critical
        {
            for (const auto & sp : localFnew) {
                Fnew[sp.first] += sp.second; 
            }
            queryCount += localQueryCount;
        }
        #pragma omp barrier

        // Update probabilities and check convergence
        #pragma omp single
        {
            for (size_t i = 0; i < spList.size(); i++) {
                Fnew[spList[i]] /= queryCount;
            }
            double delta = 0.0;
            for (size_t i = 0; i < spList.size(); i++) {
                delta += fabs(Fnew[spList[i]] - taxProbs[spList[i]]);
                if (iter > 10 && Fnew[spList[i]] < 1e-5) {
                    Fnew[spList[i]] = 0.0;
                }
            }
            cout << "Delta: " << delta << endl;
            taxProbs.swap(Fnew);    
            if (delta < 1e-6) {
                converged.fetch_add(1, std::memory_order_relaxed);
            }
        }
    }
    } // end of parallel region

    size_t explainedQueryCount = 0;
    for (size_t i = 0; i < spList.size(); i++) {
        emTaxCounts[spList[i]] = taxProbs[spList[i]] * queryCount;
        explainedQueryCount += emTaxCounts[spList[i]];
    }
    emTaxCounts[0] = totalQueryCnt - explainedQueryCount; // Unclassified queries
    cout << double(time(nullptr) - st) << " s" << endl;
    reclassify(queryRanges, mappingResList, sp2lengthFactor, totalQueryCnt);
}


void Classifier::reclassify(
    const std::vector<std::pair<size_t, size_t>> & queryRanges,
    const MappingRes * mappingResList,
    const vector<double> & sp2lengthFactor,
    size_t totalQueryCnt) 
{
    #pragma omp parallel default(none) shared(queryRanges, mappingResList, sp2lengthFactor, std::cout)
    {
        std::vector<std::pair<TaxID, double>> sp2prob;
        std::vector<TaxID> candidateSpecies;
        std::unordered_map<TaxID, unsigned int> localTaxCounts;

        #pragma omp for schedule(dynamic, 1)
        for (size_t i = 0; i < queryRanges.size(); ++i) {
            double denom = 0.0;
            candidateSpecies.clear();
            sp2prob.clear();

            TaxID currentQuery = mappingResList[queryRanges[i].first].queryId;
            for (size_t j = queryRanges[i].first; j < queryRanges[i].second; ++j) {
                double score = taxProbs[mappingResList[j].speciesId] * mappingResList[j].score * sp2lengthFactor[mappingResList[j].speciesId];
                denom += score;
                sp2prob.emplace_back(mappingResList[j].speciesId, score);
            }
            if (denom == 0.0) {
                emResults[currentQuery].taxId = 0;
                emResults[currentQuery].score = 0.0;
                continue;
            }

            for (auto & sp : sp2prob) {
                sp.second /= denom; // Normalize scores
            }

            std::sort(sp2prob.begin(), sp2prob.end(),
                      [](const std::pair<TaxID, double> &a, const std::pair<TaxID, double> &b) {
                          return a.second > b.second;
                      });

            double sum = 0.0;
            for (size_t j = 0; j < sp2prob.size() && sum < 0.5; ++j) {
                sum += sp2prob[j].second;
                candidateSpecies.push_back(sp2prob[j].first);
            }
            emResults[currentQuery].taxId = taxonomy->LCA(candidateSpecies)->taxId;
            emResults[currentQuery].score = sum;
            localTaxCounts[emResults[currentQuery].taxId] += 1;
        }
        #pragma omp critical
        {
            for (const auto & taxCount : localTaxCounts) {
                reclassifyTaxCounts[taxCount.first] += taxCount.second;
            }
        }
    }
    size_t explainedQueryCnt = 0;
    for (const auto & taxCount : reclassifyTaxCounts) {
        if (taxCount.first != 0) { // Skip unclassified
            explainedQueryCnt += taxCount.second;
        }
    }
//    std::cout << double(time(nullptr) - st) << " s" << endl;
}

void Classifier::countUniqueKmerPerSpecies(
    vector<uint32_t> & sp2uniqKmerCnt) 
{
    string sp2uniqKmerCntFileName = dbDir + "/sp2uniqKmerCnt";
    if (FileUtil::fileExists(sp2uniqKmerCntFileName.c_str())) {
        cout << "Loading unique k-mer count per species from file: " << sp2uniqKmerCntFileName << endl;
        ifstream sp2uniqKmerCntFile(sp2uniqKmerCntFileName);
        if (!sp2uniqKmerCntFile.is_open()) {
            cerr << "Error: Could not open file " << sp2uniqKmerCntFileName << endl;
            return;
        }
        TaxID taxId;
        uint32_t count;
        while (sp2uniqKmerCntFile >> taxId >> count) {
            if (taxId < sp2uniqKmerCnt.size()) {
                sp2uniqKmerCnt[taxId] = count;
            } else {
                cerr << "Warning: TaxID " << taxId << " exceeds the size of sp2uniqKmerCnt vector." << endl;
            }
        }
        sp2uniqKmerCntFile.close();
        return;
    } else {
        string infoFileName = dbDir + "/info";
        ReadBuffer<TaxID> readBuffer(infoFileName, 1000000);
        TaxID taxId;
        time_t st = time(nullptr);
        unordered_map<TaxID, TaxID> & taxid2speciesId = kmerMatcher->getTaxId2SpeciesId();
        cout << "Count unique k-mers per species: " << std::flush;
        while((taxId = readBuffer.getNext()) != 0) {
            TaxID speciesTaxId = taxid2speciesId[taxId];
            if (speciesTaxId) {
                ++sp2uniqKmerCnt[speciesTaxId];
            }
        }
        cout << double(time(nullptr) - st) << " s" << endl;
    
        // Save the unique k-mer count per species to a file
        ofstream sp2uniqKmerCntFile(sp2uniqKmerCntFileName);
        if (!sp2uniqKmerCntFile.is_open()) {
            cerr << "Error: Could not open file " << sp2uniqKmerCntFileName << " for writing." << endl;
            return;
        }
        for (size_t i = 0; i < sp2uniqKmerCnt.size(); ++i) {
            if (sp2uniqKmerCnt[i] > 0) {
                sp2uniqKmerCntFile << i << " " << sp2uniqKmerCnt[i] << "\n";
            }
        }
        sp2uniqKmerCntFile.close();
    }
}

void Classifier::loadMappings(const string & mappingResFileName) {
    size_t fileSize = FileUtil::getFileSize(mappingResFileName);
    mappingResListSize = fileSize / sizeof(MappingRes);
    mappingResList = new MappingRes[mappingResListSize];
    FILE * mappingResFile = fopen(mappingResFileName.c_str(), "rb");
    if (mappingResFile == nullptr) {
        cerr << "Error: Could not open mapping results file " << mappingResFileName << endl;
        return;
    }
    size_t readCount = fread(mappingResList, sizeof(MappingRes), mappingResListSize, mappingResFile);
    if (readCount != mappingResListSize) {
        cerr << "Error: Could not read all mapping results from file " << mappingResFileName << endl;
        fclose(mappingResFile);
        return;
    }
    fclose(mappingResFile);
}

void Classifier::loadOriginalResults(
    const string & classificationFileName,
    size_t seqNum) 
{
    ifstream classificationFile(classificationFileName);
    if (!classificationFile.is_open()) {
        cerr << "Error: Could not open classification results file " << classificationFileName << endl;
        return;
    }
    emResults.reserve(seqNum);
    string line;
    while (getline(classificationFile, line)) {
        if (line.empty() || line[0] == '#') {
            continue; // Skip empty lines and comments
        }
        vector<string> columns = TaxonomyWrapper::splitByDelimiter(line, "\t", 4);
        if (columns.size() < 4) {
            cerr << "Error: Invalid line format in classification results file: " << line << endl;
            continue;
        }
        string queryName = columns[1];
        int length = atoi(columns[3].c_str());
        emResults.emplace_back(queryName, length);
    }
}