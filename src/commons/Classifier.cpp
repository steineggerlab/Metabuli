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
}

void Classifier::startClassify(const LocalParameters &par) {
    Buffer<QueryKmer> queryKmerBuffer;
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
        // if (tries == 1) {
        //         cout << "Deviding a query file ... " << std::flush;
        // }
        queryIndexer->setBytesPerKmer(matchPerKmer);
        queryIndexer->indexQueryFile(processedReadCnt);
        const vector<QuerySplit> & queryReadSplit = queryIndexer->getQuerySplits();

        if (tries == 1) {
            totalSeqCnt = queryIndexer->getReadNum_1();
            // cout << "Done" << endl;
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
            memset(queryKmerBuffer.buffer, 0, queryReadSplit[splitIdx].kmerCnt * sizeof(QueryKmer));

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
            bool searchComplete = false;
            searchComplete = kmerMatcher->matchKmers(&queryKmerBuffer, &matchBuffer);
            if (searchComplete) {
                cout << "K-mer match count      : " << kmerMatcher->getTotalMatchCnt() << endl;
                kmerMatcher->sortMatches(&matchBuffer);
                assignTaxonomy(matchBuffer.buffer, matchBuffer.startIndexOfReserve, queryList, par);
                reporter->writeReadClassification(queryList);
                
                 for (size_t i = 0; i < queryList.size(); i++) {
                    if (queryList[i].isClassified) {
                        mappingResList.emplace_back(queryList[i].name, processedReadCnt + i, queryList[i].taxScore);
                    }
                }
                processedReadCnt += queryReadSplit[splitIdx].readCnt;
                cout << "Processed read count   : " << processedReadCnt << " (" << (double) processedReadCnt / (double) totalSeqCnt << ")" << endl;
                // numOfTatalQueryKmerCnt += queryKmerBuffer.startIndexOfReserve;
            } else { // search was incomplete
                matchPerKmer += 4;
                cout << "--match-per-kmer was increased to " << matchPerKmer << " and searching again..." << endl;
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

        cout << "--------------------" << endl;
    }


    em();
    // cout << "Number of query k-mers: " << numOfTatalQueryKmerCnt << endl;
    cout << "Total k-mer match count: " << kmerMatcher->getTotalMatchCnt() << endl;
    reporter->closeReadClassificationFile();
    reporter->writeReportFile(totalSeqCnt, taxCounts);
    reporter->writeReportFile2(totalSeqCnt, taxCounts_EM);
    reporter->writeMappingRes(mappingResList);

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

void Classifier::em() {
    // unordered_map<TaxID, unsigned int> sp2uniqKmerCnt;
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
    std::unordered_set<TaxID> speciesSet;
    time_t st = time(nullptr);
    cout << "EM algorithm for taxonomic assignment: " << std::flush;
    for (size_t i = 0; i < mappingResList.size(); i++) {
        MappingRes &res = mappingResList[i];
        for (const auto & sp2score : res.sp2score) {
            if (sp2score.second > 0) {
                speciesSet.insert(sp2score.first);
            }
        }
    }

    cout << "Species count: " << speciesSet.size() << std::endl;

    std::unordered_map<TaxID, double> F, Fnew;
    for (const auto & sp : speciesSet) {
        F[sp] = 1.0 / speciesSet.size();
        Fnew[sp] = 0.0;
    }

    for (size_t iter = 0; iter < 1000; ++iter) {
        cout << "Iteration " << iter + 1 << ": " << std::endl;
        for (auto & sp : Fnew) {
            sp.second = 0.0; // Reset Fnew
        }
        for (size_t i = 0; i < mappingResList.size(); i++) {
            MappingRes &res = mappingResList[i];
            double denom = 0.0;
            for (const auto & sp2score : res.sp2score) {
                denom += sp2score.second * F[sp2score.first] * sp2lengthFactor[sp2score.first];
            }
            for (const auto & sp2score : res.sp2score) {
                Fnew[sp2score.first] += (sp2score.second * F[sp2score.first] * sp2lengthFactor[sp2score.first]) / denom;
            }
        }
        for (const auto & sp : speciesSet) {
            Fnew[sp] /= mappingResList.size(); // Normalize
        }
        // cout << "Total probability: " << total << endl;
        double delta = 0.0;
        for (const auto & sp : speciesSet) {
            delta += fabs(Fnew[sp] - F[sp]);
        }
        if (delta < 1e-6) {
            break; // Convergence
        }
        F.swap(Fnew); // Update F with Fnew
        cout << "Delta: " << delta << endl;
    }

    // Update mapping results with final probabilities
    for (size_t i = 0; i < mappingResList.size(); i++) {
        MappingRes &res = mappingResList[i];
        double bestScore = -1.0;
        TaxID bestTaxId = 0;
        for (const auto & sp2score : res.sp2score) {
            double score = F[sp2score.first] * sp2score.second * sp2lengthFactor[sp2score.first];
            if (score > bestScore) {
                bestScore = score;
                bestTaxId = sp2score.first;
            }
        }
        res.taxId = bestTaxId;
        taxCounts_EM[bestTaxId] += 1;
    }
    cout << double(time(nullptr) - st) << " s" << endl;
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