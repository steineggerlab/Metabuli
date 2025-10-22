#include "IndexCreator.h"
#include "FileUtil.h"
#include "LocalUtil.h"
#include "ProdigalWrapper.h"
#include <cstdint>
#include <cstdio>
#include <unordered_map>
#include <utility>
#include "NcbiTaxonomy.cpp"
#include "common.h"

extern const char *version;

IndexCreator::IndexCreator(
    const LocalParameters & par, 
    TaxonomyWrapper * taxonomy,
    int kmerFormat) 
    : par(par), taxonomy(taxonomy), kmerFormat(kmerFormat) 
{
    dbDir = par.filenames[0];
    fnaListFileName = par.filenames[1];
    if (par.filenames.size() >= 3) {
        acc2taxidFileName = par.filenames[2];
    }
    
    taxidListFileName = dbDir + "/taxID_list";
    taxonomyBinaryFileName = dbDir + "/taxonomyDB";
    versionFileName = dbDir + "/db.version";
    paramterFileName = dbDir + "/db.parameters";

    if (par.reducedAA == 1){
        MARKER = 0Xffffffff;
        MARKER = ~ MARKER;
    } else {
        MARKER = 16777215;
        MARKER = ~ MARKER;
    }
    geneticCode = new GeneticCode(par.reducedAA == 1);
    kmerExtractor = new KmerExtractor(par, *geneticCode, kmerFormat);
    isUpdating = false;
    subMat = new NucleotideMatrix(par.scoringMatrixFile.values.nucleotide().c_str(), 1.0, 0.0);
}

IndexCreator::IndexCreator(
    const LocalParameters & par, 
    UnirefTree * unirefTree,
    int kmerFormat) 
    : par(par), unirefTree(unirefTree), kmerFormat(kmerFormat) 
{
    dbDir = par.filenames[0];
    geneticCode = new GeneticCode(par.reducedAA == 1);
    kmerExtractor = new KmerExtractor(par, *geneticCode, kmerFormat);
    isUpdating = false;
}

IndexCreator::IndexCreator(
    const LocalParameters & par, 
    int kmerFormat) : par(par), kmerFormat(kmerFormat) 
{
    dbDir = par.filenames[0];
    geneticCode = new GeneticCode(par.reducedAA == 1);
    kmerExtractor = new KmerExtractor(par, *geneticCode, kmerFormat);
    isUpdating = false;
    subMat = new NucleotideMatrix(par.scoringMatrixFile.values.nucleotide().c_str(), 1.0, 0.0);
}


IndexCreator::~IndexCreator() {
    delete geneticCode;
    delete kmerExtractor;
    delete subMat;
}

void IndexCreator::createLcaKmerIndex() {
    Buffer<Kmer> kmerBuffer(Buffer<Kmer>::calculateBufferSize(par.ramUsage, par.threads, sizeof(Kmer) + sizeof(size_t)));
    Buffer<size_t> uniqKmerIdx(kmerBuffer.bufferSize);
    vector<pair<size_t, size_t>> uniqKmerIdxRanges;

    string fileName = par.filenames[1];
    KSeqWrapper * kseq = KSeqFactory(fileName.c_str());
    std::unordered_map<string, uint32_t> name2id;
    uint32_t idOffset = 0;

    std::cout << "Filling UniRef100 name to ID mapping ... " << std::endl;
    time_t start = time(nullptr);
    unirefTree->getName2Id(name2id);
    cout << "Mapping filled in " << time(nullptr) - start << " s." << endl;


    bool complete = false;
    SeqEntry savedSeq;
    uint32_t processedSeqCnt = 0;
    while (!complete) {
        // Extract k-mers
        time_t start = time(nullptr);
        cout << "K-mer extraction    : " << flush;
        bool moreData = kmerExtractor->extractUnirefKmers(kseq, kmerBuffer, name2id, processedSeqCnt, savedSeq);
        complete = !moreData;
        cout << double(time(nullptr) - start) << " s" << endl;
        cout << "Processed sequences : " << processedSeqCnt << endl;

        // Sort the k-mers
        start = time(nullptr);
        cout << "Sort k-mers         : " << flush;
        SORT_PARALLEL(kmerBuffer.buffer, kmerBuffer.buffer + kmerBuffer.startIndexOfReserve, Kmer::compareKmer);
        cout << double(time(nullptr) - start) << " s" << endl;

        // Filter k-mers
        start = time(nullptr);
        size_t selectedKmerCnt = 0;
        uniqKmerIdxRanges.clear();
        filterKmers<FilterMode::UNIREF_LCA>(kmerBuffer, uniqKmerIdx.buffer, selectedKmerCnt, uniqKmerIdxRanges);
        cout << "Filter k-mers       : " << time(nullptr) - start << " s" << endl; 
        cout << "Selected k-mers     : " << selectedKmerCnt << endl;

        // Write k-mers
        start = time(nullptr);
        if (complete && numOfFlush == 0 && !isUpdating) {
            writeTargetFilesAndSplits(kmerBuffer, uniqKmerIdx.buffer, selectedKmerCnt, uniqKmerIdxRanges, true);
        } else {
            writeTargetFiles(kmerBuffer, uniqKmerIdx.buffer, uniqKmerIdxRanges);
        }
        cout << "Write k-mers        : " << time(nullptr) - start << " s" << endl;

        if (!complete) {
            kmerBuffer.init();
            uniqKmerIdx.init();
        }
        cout << "--------" << endl;
    }
    
    // for (int i = 0; i < 47; i++) {
    //     addFilesToMerge(dbDir + "/" + to_string(i) + "_diffIdx",
    //                     dbDir + "/" + to_string(i) + "_info");
    // }

    if (numOfFlush == 1) {
        cout << "Index creation completed." << endl;
        return;
    }
    cout << "Merge reference DB files ... " << endl;
    printFilesToMerge();
    setMergedFileNames(
        par.filenames[0] + "/diffIdx",  
        par.filenames[0] + "/info", 
        par.filenames[0] + "/split");
    mergeTargetFiles<FilterMode::UNIREF_LCA>();

}

void IndexCreator::createUniqueKmerIndex() {
    Buffer<Kmer> kmerBuffer(Buffer<Kmer>::calculateBufferSize(par.ramUsage, par.threads, sizeof(Kmer) + sizeof(size_t)));
    Buffer<size_t> uniqKmerIdx(kmerBuffer.bufferSize);
    bool complete = false;
    string fileName = par.filenames[1];
    KSeqWrapper * kseq = KSeqFactory(fileName.c_str());
    std::unordered_map<string, uint32_t> accession2index;
    uint32_t idOffset = 0;
    SeqEntry savedSeq;
    vector<pair<size_t, size_t>> uniqKmerIdxRanges;
    while (!complete) {
        // Extract k-mers
        time_t start = time(nullptr);
        cout << "K-mer extraction    : " << flush;
        bool moreData = kmerExtractor->extractKmers(kseq, kmerBuffer, accession2index, idOffset, savedSeq);
        complete = !moreData;
        cout << double(time(nullptr) - start) << " s" << endl;
        cout << "Processed sequences : " << idOffset << endl;

        // Sort the k-mers
        start = time(nullptr);
        cout << "Sort k-mers         : " << flush;
        SORT_PARALLEL(kmerBuffer.buffer, kmerBuffer.buffer + kmerBuffer.startIndexOfReserve, Kmer::compareKmer);
        cout << double(time(nullptr) - start) << " s" << endl;

        // Filter k-mers
        
        start = time(nullptr);
        size_t selectedKmerCnt = 0;
        uniqKmerIdxRanges.clear();
        filterKmers<FilterMode::UNIQ_KMER>(kmerBuffer, uniqKmerIdx.buffer, selectedKmerCnt, uniqKmerIdxRanges);
        cout << "Filter k-mers       : " << time(nullptr) - start << " s" << endl; 
        cout << "Selected k-mers     : " << selectedKmerCnt << endl;

        // Write k-mers
        
        start = time(nullptr);
        if (complete && numOfFlush == 0 && !isUpdating) {
            writeTargetFilesAndSplits(kmerBuffer, uniqKmerIdx.buffer, selectedKmerCnt, uniqKmerIdxRanges, true);
        } else {
            writeTargetFiles(kmerBuffer, uniqKmerIdx.buffer, uniqKmerIdxRanges);
        }
        cout << "Write k-mers        : " << time(nullptr) - start << " s" << endl;

        // Write accession to index mapping
        FILE * accIndexFile = fopen((dbDir + "/accession2index").c_str(), "a");
        for (const auto & entry : accession2index) {
            fprintf(accIndexFile, "%s\t%u\n", entry.first.c_str(), entry.second);
        }
        fclose(accIndexFile);

        if (!complete) {
            accession2index.clear();
            kmerBuffer.init();
            uniqKmerIdx.init();
        }
        cout << "--------" << endl;
    }
    

    if (numOfFlush == 1) {
        cout << "Index creation completed." << endl;
        return;
    }
    cout << "Merge reference DB files ... " << endl;

    // for (int i = 0; i < 10; i++) {
    //     addFilesToMerge(dbDir + "/" + to_string(i) + "_diffIdx",
    //                     dbDir + "/" + to_string(i) + "_info");
    // }

    printFilesToMerge();
    setMergedFileNames(
        par.filenames[0] + "/diffIdx",  
        par.filenames[0] + "/info", 
        par.filenames[0] + "/split");
    mergeTargetFiles<FilterMode::UNIQ_KMER>();

}

void IndexCreator::createCommonKmerIndex() {
    size_t sizE = calculateBufferSize(par.ramUsage);
    indexReferenceSequences(sizE);
    Buffer<Kmer> kmerBuffer(sizE);
    Buffer<size_t> uniqKmerIdx(sizE);
    if (!par.cdsInfo.empty()) {
        cout << "Loading CDS info from: " << par.cdsInfo << endl;
        loadCdsInfo(par.cdsInfo);
    }

    // Write taxonomy ID list
    FILE * taxidListFile = fopen(taxidListFileName.c_str(), "w");
    for (auto & taxid: taxIdSet) {
        fprintf(taxidListFile, "%d\n", taxid);
    }
    fclose(taxidListFile);

    std::vector<std::atomic<bool>> batchChecker(accessionBatches.size());
    size_t processedBatchCnt = 0;
    
    vector<pair<size_t, size_t>> uniqKmerIdxRanges;
    while(processedBatchCnt < accessionBatches.size()) {
        // Extract target k-mers
        time_t start = time(nullptr);
        cout << "K-mer extraction : " << flush;
        fillTargetKmerBuffer(kmerBuffer, batchChecker, processedBatchCnt, par);
        cout << double(time(nullptr) - start) << " s" << endl;
        
        // Sort the k-mers
        start = time(nullptr);
        cout << "Sort k-mers      : " << flush;
        SORT_PARALLEL(kmerBuffer.buffer, kmerBuffer.buffer + kmerBuffer.startIndexOfReserve, Kmer::compareTargetKmer);
        cout << time(nullptr) - start << " s" << endl;

        // Filter k-mers
        start = time(nullptr);
        size_t selectedKmerCnt = 0;
        filterKmers<FilterMode::DB_CREATION>(kmerBuffer, uniqKmerIdx.buffer, selectedKmerCnt, uniqKmerIdxRanges);
        cout << "Filter k-mers    : " << time(nullptr) - start << " s" << endl; 
        cout << "Selected k-mers  : " << selectedKmerCnt << endl;

        // Write the target files
        start = time(nullptr);
        if(processedBatchCnt == accessionBatches.size() && numOfFlush == 0) {
            writeTargetFilesAndSplits(kmerBuffer, uniqKmerIdx.buffer, selectedKmerCnt, uniqKmerIdxRanges, true);
        } else {
            writeTargetFiles(kmerBuffer, uniqKmerIdx.buffer, uniqKmerIdxRanges);
        }
        cout << "Write k-mers     : " << time(nullptr) - start << " s" << endl;

        // Reset buffers
        if (processedBatchCnt < accessionBatches.size()) {
            kmerBuffer.init();
            uniqKmerIdx.init();
            uniqKmerIdxRanges.clear();
        }
    }

    taxonomy->writeTaxonomyDB(par.filenames[0] + "/taxonomyDB");
    writeDbParameters();
    
    if (numOfFlush == 1) {
        cout << "Index creation completed." << endl;
        return;
    }
    cout << "Merge reference DB files ... " << endl;

    // for (int i = 0; i < 243; i++) {
    //     addFilesToMerge(dbDir + "/" + to_string(i) + "_diffIdx",
    //                     dbDir + "/" + to_string(i) + "_info");
    // }
    // updateTaxId2SpeciesTaxId(dbDir + "/taxID_list");

    printFilesToMerge();
    setMergedFileNames(
        par.filenames[0] + "/diffIdx",  
        par.filenames[0] + "/info", 
        par.filenames[0] + "/split");
    mergeTargetFiles<FilterMode::COMMON_KMER>();
}

void IndexCreator::createIndex() {
    Buffer<Kmer> kmerBuffer(calculateBufferSize(par.ramUsage));
    cout << "Target metamer buffer size: " << kmerBuffer.bufferSize << endl;
    
    indexReferenceSequences(kmerBuffer.bufferSize);
    cout << "Made blocks for each thread" << endl;

    if (!par.cdsInfo.empty()) {
        cout << "Loading CDS info from: " << par.cdsInfo << endl;
        loadCdsInfo(par.cdsInfo);
    }

    // Write taxonomy id list
    FILE * taxidListFile = fopen(taxidListFileName.c_str(), "w");
    for (auto & taxid: taxIdSet) {
        fprintf(taxidListFile, "%d\n", taxid);
    }
    fclose(taxidListFile);

    // Process the splits until all are processed
    std::vector<std::atomic<bool>> batchChecker(accessionBatches.size());
    size_t processedBatchCnt = 0;
    
    vector<pair<size_t, size_t>> uniqKmerIdxRanges;
#ifdef OPENMP
    omp_set_num_threads(par.threads);
#endif
    while(processedBatchCnt < accessionBatches.size()) {
        kmerBuffer.init();
        cout << "Buffer initialized" << endl;

        // Extract target k-mers
        fillTargetKmerBuffer(kmerBuffer, batchChecker, processedBatchCnt, par);

        // Sort the k-mers
        time_t start = time(nullptr);
        SORT_PARALLEL(kmerBuffer.buffer, kmerBuffer.buffer + kmerBuffer.startIndexOfReserve,
                      Kmer::compareTargetKmer);
        time_t sort = time(nullptr);
        cout << "Reference k-mer sort : " << sort - start << endl;

        // Reduce redundancy
        auto * uniqKmerIdx = new size_t[kmerBuffer.startIndexOfReserve + 1];
        size_t uniqKmerCnt = 0;
        uniqKmerIdxRanges.clear();
        filterKmers<FilterMode::DB_CREATION>(kmerBuffer, uniqKmerIdx, uniqKmerCnt, uniqKmerIdxRanges);
        time_t reduction = time(nullptr);
        cout << "Unique k-mer count   : " << uniqKmerCnt << endl;
        cout << "Redundancy reduction : " << (double) (reduction - sort) << " s" << endl;

        // Write the target files
        if(processedBatchCnt == accessionBatches.size() && numOfFlush == 0 && !isUpdating) {
            writeTargetFilesAndSplits(kmerBuffer, uniqKmerIdx, uniqKmerCnt, uniqKmerIdxRanges);
        } else {
            writeTargetFiles(kmerBuffer, uniqKmerIdx, uniqKmerIdxRanges);
        }
        delete[] uniqKmerIdx;
    }

    writeDbParameters();
}


string IndexCreator::addToLibrary(const std::string & dbDir,
                                  const std::string & fnaListFileName,
                                  const std::string & acc2taxIdFileName) {
    // Make library directory
    time_t now = time(0);
    tm *ltm = localtime(&now);
    string timeStr = to_string(1900 + ltm->tm_year) + "-" + to_string(1 + ltm->tm_mon) + "-" + to_string(ltm->tm_mday) + "-" + to_string(ltm->tm_hour) + "-" + to_string(ltm->tm_min);
    string libraryDir = dbDir + "/" + timeStr;
    if (!FileUtil::directoryExists(libraryDir.c_str())) {
        FileUtil::makeDir(libraryDir.c_str());
    }

    unordered_map<std::string, TaxID> accession2taxid;
    vector<std::string> fileNames;
    vector<string> libraryFiles;
    unordered_set<TaxID> observedSpecies;
    getObservedAccessionList(fnaListFileName, fileNames, accession2taxid);
    fillAcc2TaxIdMap(accession2taxid, acc2taxIdFileName);

    unordered_map<TaxID, TaxID> original2internalTaxId;
    if (taxonomy->hasInternalTaxID()) {
      taxonomy->getOriginal2InternalTaxId(original2internalTaxId);
    } 
    
    vector<std::string> unmapped;
    std::string accession;
    for (size_t i = 0; i < fileNames.size(); ++i) {
      KSeqWrapper* kseq = KSeqFactory(fileNames[i].c_str());
      while (kseq->ReadEntry()) {
        const KSeqWrapper::KSeqEntry & e = kseq->entry;
        accession = string(e.name.s);
        size_t pos = accession.find('.');
        if (pos != std::string::npos) {
          accession = accession.substr(0, pos);
        }
        TaxID taxId = accession2taxid[accession];
        if (taxId == 0) {
          std::cout << "During processing " << fileNames[i] << ", accession " << e.name.s <<
               " is not found in the mapping file. It is skipped.\n";
          unmapped.push_back(e.name.s);
          continue;
        }
        if (taxonomy->hasInternalTaxID()) {
          if (original2internalTaxId.find(taxId) == original2internalTaxId.end()) {
            std::cout << "During processing " << fileNames[i] << ", accession " << e.name.s <<
            ", " << taxId << " is not included in the taxonomy. It is skipped.\n";
            unmapped.push_back(e.name.s);
            continue;
          }
          taxId = original2internalTaxId[taxId];
        }
        int speciesTaxID = taxonomy->getTaxIdAtRank(taxId, "species");
        if (speciesTaxID == 0) {
          cout << "During processing " << fileNames[i] << ", accession " << e.name.s <<
               " is not matched to any species. It is skipped.\n";
          unmapped.push_back(e.name.s);
          continue;
        }
        if (taxonomy->hasInternalTaxID()) {
            speciesTaxID = taxonomy->getOriginalTaxID(speciesTaxID);
        }
        string speciesFileName = libraryDir + "/" + to_string(speciesTaxID) + ".fna";
        if (observedSpecies.find(speciesTaxID) == observedSpecies.end()) {
            observedSpecies.insert(speciesTaxID);
            libraryFiles.push_back(speciesFileName);
        }
        FILE *file = fopen(speciesFileName.c_str(), "a");
        fprintf(file, ">%s %s\n", e.name.s, e.comment.s);
        fprintf(file, "%s\n", e.sequence.s);
        fclose(file);
      }
      delete kseq;
    }

    // Write unmapped accession to file
    if (unmapped.empty()) {
        cout << "All accessions are mapped to taxonomy" << endl;
    } else {
        FILE *file = fopen((libraryDir + "/unmapped.txt").c_str(), "w");
        for (const auto & i : unmapped) {
            fprintf(file, "%s\n", i.c_str());
        }
        fclose(file);
        cout << "Unmapped accessions are written to " << libraryDir + "/unmapped.txt" << endl;
    }  
    // Write the list of absolute paths of library files
    string libraryListFileName = libraryDir + "/library.list";
    FILE *libraryListFile = fopen(libraryListFileName.c_str(), "w");
    for (size_t i = 0; i < libraryFiles.size(); ++i) {
        fprintf(libraryListFile, "%s\n", libraryFiles[i].c_str());
    }
    fclose(libraryListFile);
    return libraryListFileName;
}


void IndexCreator::indexReferenceSequences(size_t bufferSize) {
    time_t start = time(nullptr);
    cout << "Make sequence batches for each thread : " << flush;
    
    vector<Accession> observedAccessionsVec;       // vector of observed accessions
    unordered_map<string, size_t> accession2index; // map accession to its index in the observedAccessionsVec 
    if (par.makeLibrary) {
        getObservedAccessions(addToLibrary(dbDir, fnaListFileName, acc2taxidFileName), observedAccessionsVec, accession2index);
    } else {
        getObservedAccessions(fnaListFileName, observedAccessionsVec, accession2index);
    }
    cout << "Number of observed accessions: " << observedAccessionsVec.size() << endl;
    getTaxonomyOfAccessions(observedAccessionsVec, accession2index, acc2taxidFileName);
    cout << "Taxonomy of accessions is obtained" << endl;
    vector<Accession> accessionsWithTaxonomy;
    getAccessionBatches(observedAccessionsVec, bufferSize);
    cout << "Number of accession batches: " << accessionBatches.size() << endl;
    sort(accessionBatches.begin(), accessionBatches.end(),
         [](const AccessionBatch &a, const AccessionBatch &b) {
             return a.totalLength < b.totalLength;
         });

    time_t end = time(nullptr);
    cout << end - start << " s" << endl;
}


void IndexCreator::getObservedAccessions(
    const string & fnaListFileName,
    vector<Accession> & observedAccessionsVec,
    unordered_map<string, size_t> & accession2index) 
{
    ifstream fileListFile(fnaListFileName);
    if (fileListFile.is_open()) {
        for (string eachLine; getline(fileListFile, eachLine);) {
            fastaPaths.push_back(eachLine);
        }
    } else {
        cout << "Cannot open file for file list" << endl;
    } 

    // Iterate through the fasta files to get observed accessions
    size_t accCnt = 0;
    size_t copyCount = 0;
    #pragma omp parallel default(none), shared(accession2index, cout, observedAccessionsVec, accCnt, copyCount)
    {
        vector<Accession> localObservedAccessionsVec;
        localObservedAccessionsVec.reserve(4096 * 4);
        unordered_set<string> duplicateCheck;
        
        #pragma omp for schedule(static, 1)
        for (size_t i = 0; i < fastaPaths.size(); ++i) {
            KSeqWrapper* kseq = KSeqFactory(fastaPaths[i].c_str());
            uint32_t order = 0;
            while (kseq->ReadEntry()) {
                const KSeqWrapper::KSeqEntry & e = kseq->entry;
                // Get the accession ID without version
                char* pos = strchr(e.name.s, '.'); 
                if (pos != nullptr) {
                    *pos = '\0';
                }
                if (duplicateCheck.find(e.name.s) != duplicateCheck.end()) {
                    continue;
                } else {
                    duplicateCheck.insert(e.name.s);
                } 
                localObservedAccessionsVec.emplace_back(string(e.name.s), i, order, e.sequence.l);
                order++; 
            }
            delete kseq;
        } 
        __sync_fetch_and_add(&accCnt, localObservedAccessionsVec.size()); 
        #pragma omp barrier
       
        #pragma omp critical
        {
            if (observedAccessionsVec.size() < accCnt) {
                observedAccessionsVec.resize(accCnt);
                accession2index.reserve(accCnt);
            }
        }   

        #pragma omp critical
        {
            size_t start = copyCount;
            for (size_t j = start ; j < start + localObservedAccessionsVec.size(); ++j) {
                observedAccessionsVec[j] = localObservedAccessionsVec[j - start];
                accession2index[observedAccessionsVec[j].accession] = j;
            }
            copyCount += localObservedAccessionsVec.size();
        }                     
    }
}

void IndexCreator::getTaxonomyOfAccessions(vector<Accession> & observedAccessionsVec,
                                           const unordered_map<string, size_t> & accession2index,
                                           const string & acc2taxidFileName) {
    unordered_map<TaxID, TaxID> old2merged;
    taxonomy->getMergedNodeMap(old2merged, true);
    
    vector<pair<string, pair<TaxID, TaxID>>> acc2accId;   
    int fd = open(acc2taxidFileName.c_str(), O_RDONLY);
    if (fd < 0) {
        cout << "Cannot open file for mapping from accession to tax ID" << endl;
        return;
    }

    // Get the size of the file
    struct stat sb;
    if (fstat(fd, &sb) == -1) {
        cout << "Cannot get the size of the file for mapping from accession to tax ID" << endl;
        close(fd);
        return;
    }

    size_t fileSize = sb.st_size;

    // Map the file to memory
    char* fileData = static_cast<char*>(mmap(nullptr, fileSize, PROT_READ, MAP_PRIVATE, fd, 0));
    if (fileData == MAP_FAILED) {
        cout << "mmap failed" << endl;
        close(fd);
        return;
    }
    close(fd);  // Close the file descriptor as it is no longer needed after mmap.

    // Parse the file
    char* current = fileData;
    char* end = fileData + fileSize;

    // Skip the header line
    while (current < end && *current != '\n') {
        ++current;
    }
    ++current;  // Move past the newline

    char accession[16384];
    TaxID taxID;
    std::unordered_set<TaxID> usedExternalTaxIDs;
    std::vector<NewTaxon> newTaxons;
    if (par.accessionLevel == 1) {
        taxonomy->getUsedExternalTaxIDs(usedExternalTaxIDs);
    }

    // First, label the accessions with external taxIDs
    while (current < end) {
        // Read a line
        char* lineStart = current;
        while (current < end && *current != '\n') {
            ++current;
        }
        std::string line(lineStart, current - lineStart);

        // Parse the line
        if (sscanf(line.c_str(), "%s\t%*s\t%d\t%*d", accession, &taxID) == 2) {
            // Get the accession ID without version
            // char* pos = strchr(accession, '.');
            // if (pos != nullptr) {
            //     *pos = '\0';
            // }
            auto it = accession2index.find(accession);
            if (it != accession2index.end()) {
                if (old2merged.count(taxID) > 0) {
                    taxID = old2merged[taxID];
                }
                if (par.accessionLevel == 1) {
                    TaxID accTaxId = taxonomy->getSmallestUnusedExternalTaxID(usedExternalTaxIDs);
                    acc2accId.emplace_back(accession, make_pair(taxID, accTaxId));
                    if (accTaxId == 0) {
                        cout << "accTaxId is 0 for accession " << accession << " " << taxID << endl;
                    }
                    observedAccessionsVec[it->second].taxID = accTaxId;
                    newTaxons.emplace_back(accTaxId, taxID, "accession", accession);
                } else {
                    observedAccessionsVec[it->second].taxID = taxID;
                }
            }
        }
        ++current;  // Move to the next line
    }

    if (munmap(fileData, fileSize) == -1) {
        cout << "munmap failed" << endl;
    }                 

    if (par.accessionLevel == 1) {
        TaxonomyWrapper * newTaxonomy = taxonomy->addNewTaxa(newTaxons);
        delete taxonomy;
        taxonomy = newTaxonomy;
    }

    // Second, convert external taxIDs to internal taxIDs
    cout << "Converting external taxIDs to internal taxIDs" << endl;
    vector<std::string> unmappedAccessions;
    std::unordered_map<TaxID, TaxID> external2internalTaxID;
    taxonomy->getExternal2internalTaxID(external2internalTaxID);
    std::unordered_set<TaxID> spIdSet;
    for (size_t i = 0; i < observedAccessionsVec.size(); ++i) {    
        auto it = external2internalTaxID.find(observedAccessionsVec[i].taxID);
        if (it == external2internalTaxID.end() || observedAccessionsVec[i].taxID == 0) {
            cout << "TaxID is not found for accession " << observedAccessionsVec[i].accession << " " << observedAccessionsVec[i].taxID << endl;
            unmappedAccessions.push_back(observedAccessionsVec[i].accession);
            continue;
        }
        observedAccessionsVec[i].taxID = it->second; // store the internal taxID
        observedAccessionsVec[i].speciesID = taxonomy->getTaxIdAtRank(it->second, "species");
        spIdSet.insert(observedAccessionsVec[i].speciesID);
        taxIdSet.insert(it->second);
        taxId2speciesId[it->second] = observedAccessionsVec[i].speciesID;
        if (observedAccessionsVec[i].speciesID == 0) {
            cout << "Species ID is not found for accession " << observedAccessionsVec[i].accession << " " << it->second << endl;
            unmappedAccessions.push_back(observedAccessionsVec[i].accession);
            exit(1);
        }
        taxId2speciesId[observedAccessionsVec[i].speciesID] = observedAccessionsVec[i].speciesID;
    }

    cout << "Number of unique taxIDs: " << taxIdSet.size() << endl;
    cout << "Number of unique speciesIDs: " << spIdSet.size() << endl;
    
    string mappingFileName = dbDir + "/acc2taxid.map";
    FILE * acc2taxidFile = nullptr;
    if (isUpdating) {
        acc2taxidFile = fopen(mappingFileName.c_str(), "a");
    } else {
        acc2taxidFile = fopen(mappingFileName.c_str(), "w");
    }
    if (par.accessionLevel == 1) {
        for (auto & acc2taxid: acc2accId) {
            fprintf(acc2taxidFile, "%s\t%d\t%d\n", acc2taxid.first.c_str(), acc2taxid.second.first, acc2taxid.second.second);
        }
    } else {
        for (size_t i = 0; i < observedAccessionsVec.size(); ++i) {
            if (observedAccessionsVec[i].taxID == 0) {
                continue;
            }
            fprintf(acc2taxidFile, "%s\t%d\n", observedAccessionsVec[i].accession.c_str(), taxonomy->getOriginalTaxID(observedAccessionsVec[i].taxID));
        }        
    }   

    if (unmappedAccessions.empty()) {
        cout << "All accessions are mapped to taxonomy" << endl;
    } else {
        FILE *file = fopen((dbDir + "/unmapped.txt").c_str(), "w");
        for (const auto & i : unmappedAccessions) {
            fprintf(file, "%s\n", i.c_str());
        }
        fclose(file);
        cout << "Unmapped accessions are written to " << dbDir + "/unmapped.txt" << endl;
    }

    fclose(acc2taxidFile);                   
}

void IndexCreator::getAccessionBatches(std::vector<Accession> & observedAccessionsVec, size_t bufferSize) {
    size_t accCnt = observedAccessionsVec.size();
    SORT_PARALLEL(observedAccessionsVec.begin(), observedAccessionsVec.end(), Accession::compare);
    vector<uint32_t> orders;
    vector<TaxID> taxIDs;
    vector<uint32_t> lengths;
    for (size_t i = 0; i < accCnt;) {
        if (observedAccessionsVec[i].speciesID == 0 || observedAccessionsVec[i].taxID == 0) {
            i++;
            continue;
        }
        TaxID currentSpeciesID = observedAccessionsVec[i].speciesID;
        uint32_t maxLength = 0, trainingFasta = 0, trainingSeq = 0;
        size_t firstBatchOfSpecies = accessionBatches.size();
        while (i < accCnt && currentSpeciesID == observedAccessionsVec[i].speciesID) {    
            uint32_t lengthSum = 0;
            size_t kmerCntSum = 0;
            orders.clear();
            lengths.clear();
            taxIDs.clear();
            uint32_t currentFasta = observedAccessionsVec[i].whichFasta;
            while (i < accCnt && currentSpeciesID == observedAccessionsVec[i].speciesID 
                   && currentFasta == observedAccessionsVec[i].whichFasta) {
                if (observedAccessionsVec[i].length > maxLength) {
                    maxLength = observedAccessionsVec[i].length;
                    trainingFasta = observedAccessionsVec[i].whichFasta;
                    trainingSeq = observedAccessionsVec[i].order;
                }
                lengthSum += observedAccessionsVec[i].length;
                kmerCntSum += static_cast<size_t>(observedAccessionsVec[i].length * 0.4);
                orders.push_back(observedAccessionsVec[i].order);
                lengths.push_back(observedAccessionsVec[i].length);
                taxIDs.push_back(observedAccessionsVec[i].taxID);
                i++;
                if (kmerCntSum > bufferSize || lengthSum > 100'000'000 || orders.size() > 300 || (orders.size() > 100 && lengthSum > 50'000'000)) {
                    break;
                }
            }
            // Add the batch
            accessionBatches.emplace_back(currentFasta, currentSpeciesID, 0, 0, lengthSum);
            accessionBatches.back().orders = orders;
            accessionBatches.back().lengths = lengths;
            accessionBatches.back().taxIDs = taxIDs;
        }
        // Update training sequence information
        for (size_t j = firstBatchOfSpecies; j < accessionBatches.size(); ++j) {
            accessionBatches[j].trainingSeqFasta = trainingFasta;
            accessionBatches[j].trainingSeqIdx = trainingSeq;
        }
    }
}


void IndexCreator::writeTargetFiles(
    Buffer<Kmer> & kmerBuffer, 
    const size_t * uniqKmerIdx,
    const vector<pair<size_t, size_t>> & uniqKmerIdxRanges) 
{
    string diffIdxFileName = dbDir + "/" + to_string(numOfFlush) + "_diffIdx";
    string infoFileName = dbDir + "/" + to_string(numOfFlush) + "_info";
    deltaIdxFileNames.push_back(diffIdxFileName);
    infoFileNames.push_back(infoFileName);

    numOfFlush++;

    size_t bufferSize = 1024 * 1024 * 32;
    WriteBuffer<uint16_t> diffBuffer(diffIdxFileName, bufferSize);
    WriteBuffer<uint32_t> infoBuffer(infoFileName, bufferSize); 
    uint64_t lastKmer = 0;

    for (size_t i = 0; i < uniqKmerIdxRanges.size(); i ++) {
        for (size_t j = uniqKmerIdxRanges[i].first; j < uniqKmerIdxRanges[i].second; j ++) {
            infoBuffer.write(&kmerBuffer.buffer[uniqKmerIdx[j]].id);
            getDiffIdx(lastKmer, kmerBuffer.buffer[uniqKmerIdx[j]].value, diffBuffer);
        }
    }
    infoBuffer.flush();
    diffBuffer.flush();

    kmerBuffer.startIndexOfReserve = 0; // Reset the buffer for the next batch
}

void IndexCreator::writeTargetFilesAndSplits(
    Buffer<Kmer> & kmerBuffer,
    const size_t * uniqKmerIdx,
    size_t & uniqKmerCnt,
    const vector<pair<size_t, size_t>> & uniqKmerIdxRanges,
    bool writeInfo)
{
    DiffIdxSplit * splitList = new DiffIdxSplit[par.splitNum];
    memset(splitList, 0, sizeof(DiffIdxSplit) * par.splitNum);    
    size_t splitWidth = uniqKmerCnt / par.splitNum;
    size_t remainder = uniqKmerCnt % par.splitNum;
    size_t splitCnt = 1;
    size_t start = 0;
    for (size_t i = 1; i < (size_t) par.splitNum; i++) {
        start = start + splitWidth;
        if (remainder > 0) {
            start++;
            remainder--;
        }
        size_t counter = 0;
        for (size_t j = 0; j < uniqKmerIdxRanges.size(); j++) {
            if (counter + (uniqKmerIdxRanges[j].second - uniqKmerIdxRanges[j].first) < start) {
                counter += uniqKmerIdxRanges[j].second - uniqKmerIdxRanges[j].first;
                continue;
            }
            bool found = false;
            for (size_t k = uniqKmerIdxRanges[j].first + start - counter; k + 1 < uniqKmerIdxRanges[j].second; k++) {
                if (AminoAcidPart(kmerBuffer.buffer[uniqKmerIdx[k]].value) 
                    != AminoAcidPart(kmerBuffer.buffer[uniqKmerIdx[k + 1]].value)) {
                    splitList[splitCnt++].ADkmer = kmerBuffer.buffer[uniqKmerIdx[k + 1]].value;
                    found = true;
                    break;
                }
            }
            if (!found) {
                splitList[splitCnt++].ADkmer = kmerBuffer.buffer[uniqKmerIdx[uniqKmerIdxRanges[j+1].first]].value;
            }
            break;
        }
    }

    numOfFlush++;
    size_t bufferSize = 1024 * 1024 * 32;
    uint64_t lastKmer = 0;
    size_t splitIdx = 1;
    WriteBuffer<uint16_t> diffBuffer(dbDir + "/diffIdx", bufferSize);
    
    if (writeInfo) {
        WriteBuffer<uint32_t> infoBuffer(dbDir + "/info", bufferSize); 
        for (size_t i = 0; i < uniqKmerIdxRanges.size(); i ++) {
            for (size_t j = uniqKmerIdxRanges[i].first; j < uniqKmerIdxRanges[i].second; j ++) {
                infoBuffer.write(&kmerBuffer.buffer[uniqKmerIdx[j]].id);
                getDiffIdx(lastKmer, kmerBuffer.buffer[uniqKmerIdx[j]].value, diffBuffer);
                if((splitIdx < splitCnt) && (lastKmer == splitList[splitIdx].ADkmer)){
                    splitList[splitIdx].diffIdxOffset = diffBuffer.writeCnt;
                    splitList[splitIdx].infoIdxOffset = infoBuffer.writeCnt;
                    splitIdx ++;
                }
            }
        }
        cout << "Written k-mer count : " << infoBuffer.writeCnt << endl;
    } else {
        size_t writeCnt = 0;
        for (size_t i = 0; i < uniqKmerIdxRanges.size(); i ++) {
            for (size_t j = uniqKmerIdxRanges[i].first; j < uniqKmerIdxRanges[i].second; j ++) {
                getDiffIdx(lastKmer, kmerBuffer.buffer[uniqKmerIdx[j]].value, diffBuffer);
                writeCnt++;
                if((splitIdx < splitCnt) && (lastKmer == splitList[splitIdx].ADkmer)){
                    splitList[splitIdx].diffIdxOffset = diffBuffer.writeCnt;
                    splitList[splitIdx].infoIdxOffset = 0;
                    splitIdx ++;
                }
            }
        }
        cout << "Written k-mer count : " << writeCnt << endl;
    }
    
    FILE * deltaIdxSplitFile = fopen((dbDir + "/split").c_str(), "wb");
    if (deltaIdxSplitFile == nullptr) {
        cout << "Cannot open the file for writing target DB" << endl;
        return;
    }
    fwrite(splitList, sizeof(DiffIdxSplit), par.splitNum, deltaIdxSplitFile);
    delete[] splitList;
    fclose(deltaIdxSplitFile);
    
    kmerBuffer.startIndexOfReserve = 0; // Reset the buffer for the next batch
}

void IndexCreator::getDiffIdx(
    uint64_t & lastKmer,
    uint64_t entryToWrite,
    WriteBuffer<uint16_t> & diffBuffer) 
{
    uint64_t kmerdiff = entryToWrite - lastKmer;
    uint16_t buffer[5];
    int idx = 3;
    buffer[4] = SET_END_FLAG(GET_15_BITS(kmerdiff));
    kmerdiff >>= 15U;
    while (kmerdiff) {
        uint16_t toWrite = GET_15_BITS(kmerdiff);
        kmerdiff >>= 15U;
        buffer[idx] = toWrite;
        idx--;
    }
    diffBuffer.write((buffer + idx + 1), (4 - idx));
    lastKmer = entryToWrite;
}


void IndexCreator::load_assacc2taxid(const string & mappingFile, unordered_map<string, int> & assacc2taxid){
    string key, value;
    ifstream map;
    map.open(mappingFile);
    if(map.is_open()){
        while(getline(map,key,'\t')){
            getline(map, value, '\n');
            assacc2taxid[key] = stoi(value);
            if(key[2] == 'F'){
                key[2] = 'A';
                assacc2taxid[key] = stoi(value);
            }
        }
    } else{
        cout<<"Cannot open file for mappig from assemlby accession to tax ID"<<endl;
    }
    map.close();
}

size_t IndexCreator::fillTargetKmerBuffer(Buffer<Kmer> &kmerBuffer,
                                          std::vector<std::atomic<bool>> & batchChecker,
                                          size_t &processedBatchCnt,
                                          const LocalParameters &par) {
    std::atomic<int> hasOverflow{0};
#pragma omp parallel default(none), shared(kmerBuffer, batchChecker, processedBatchCnt, hasOverflow, par, cout)
    {
        ProbabilityMatrix probMatrix(*subMat);
        SeqIterator seqIterator(par);
        size_t posToWrite;
        size_t orfNum;
        vector<SequenceBlock> extendedORFs;
        priority_queue<uint64_t> standardList;
        priority_queue<uint64_t> currentList;
        size_t lengthOfTrainingSeq = 0;
        char *reverseComplement;
        vector<uint64_t> intergenicKmers;
        vector<string> cds;
        vector<string> nonCds;
        bool trained = false;
        size_t estimatedKmerCnt = 0;
#pragma omp for schedule(dynamic, 1)
        for (size_t batchIdx = 0; batchIdx < accessionBatches.size(); batchIdx ++) {
            if (hasOverflow.load(std::memory_order_acquire))
                continue;
            
            if (batchChecker[batchIdx].exchange(true, std::memory_order_acq_rel))
                continue; 
            
            intergenicKmers.clear();
            standardList = priority_queue<uint64_t>();

            // Estimate the number of k-mers to be extracted from current split
            size_t totalLength = 0;
            for (size_t p = 0; p < accessionBatches[batchIdx].lengths.size(); p++) {
                totalLength += accessionBatches[batchIdx].lengths[p];
            }

            if (par.syncmer) {
                estimatedKmerCnt = static_cast<size_t>(
                    (totalLength * 1.3 / 3.0) / ((8 - par.smerLen + 1) / 2.0)
                );
            } else {
                estimatedKmerCnt = static_cast<size_t>(
                    (totalLength * 1.3) / 3.0
                );
            }
                
            ProdigalWrapper * prodigal = new ProdigalWrapper();
            trained = false;

            // Process current split if buffer has enough space.
            posToWrite = kmerBuffer.reserveMemory(estimatedKmerCnt);
            if (posToWrite + estimatedKmerCnt < kmerBuffer.bufferSize) {
                KSeqWrapper* kseq = KSeqFactory(fastaPaths[accessionBatches[batchIdx].whichFasta].c_str());
                size_t seqCnt = 0;
                size_t idx = 0;
                while (kseq->ReadEntry()) {
                    if (seqCnt == accessionBatches[batchIdx].orders[idx]) {
                        if (accessionBatches[batchIdx].taxIDs[idx] == 0) {
                            #pragma omp critical
                            {
                            accessionBatches[batchIdx].print();
                            exit(1);
                            }
                        }
                        const KSeqWrapper::KSeqEntry & e = kseq->entry;
                        // Mask low complexity regions
                        char *maskedSeq = nullptr;
                        if (par.maskMode) {
                            maskedSeq = new char[e.sequence.l + 1]; // TODO: reuse the buffer
                            SeqIterator::maskLowComplexityRegions((unsigned char *) e.sequence.s, (unsigned char *) maskedSeq, probMatrix, par.maskProb, subMat);
                            maskedSeq[e.sequence.l] = '\0';
                        } else {
                            maskedSeq = e.sequence.s;
                        }

                        orfNum = 0;
                        extendedORFs.clear();
                        int tempCheck = 0;
                        if (cdsInfoMap.find(string(e.name.s)) != cdsInfoMap.end()) {
                            // Get CDS and non-CDS
                            cds.clear();
                            nonCds.clear();
                            seqIterator.devideToCdsAndNonCds(maskedSeq,
                                                             e.sequence.l,
                                                             cdsInfoMap[string(e.name.s)],
                                                             cds,
                                                             nonCds);

                            for (size_t cdsCnt = 0; cdsCnt < cds.size(); cdsCnt ++) {
                                tempCheck = kmerExtractor->extractTargetKmers(
                                                cds[cdsCnt].c_str(),
                                                kmerBuffer,
                                                posToWrite,
                                                accessionBatches[batchIdx].taxIDs[idx],
                                                accessionBatches[batchIdx].speciesID,
                                                {0, (int) cds[cdsCnt].length() - 1, 1});
                                if (tempCheck == -1) {
                                    cout << "ERROR: Buffer overflow " << e.name.s << e.sequence.l << endl;
                                }
                            }
                            for (size_t nonCdsCnt = 0; nonCdsCnt < nonCds.size(); nonCdsCnt ++) {
                                tempCheck = kmerExtractor->extractTargetKmers(
                                                nonCds[nonCdsCnt].c_str(),
                                                kmerBuffer,
                                                posToWrite,
                                                accessionBatches[batchIdx].taxIDs[idx],
                                                accessionBatches[batchIdx].speciesID,
                                                {0, (int) cds[nonCdsCnt].length() - 1, 1});
                                if (tempCheck == -1) {
                                    cout << "ERROR: Buffer overflow " << e.name.s << e.sequence.l << endl;
                                }
                            }
                        } else {
                            // USE PRODIGAL
                            if (!trained) {
                                KSeqWrapper* training_seq = KSeqFactory(fastaPaths[accessionBatches[batchIdx].trainingSeqFasta].c_str()); 
                                size_t seqCnt = 0;
                                while (training_seq->ReadEntry()) {
                                    if (seqCnt == accessionBatches[batchIdx].trainingSeqIdx) {
                                        break;
                                    }
                                    seqCnt++;
                                }
                                lengthOfTrainingSeq = training_seq->entry.sequence.l;
                                prodigal->is_meta = 0;
                                if (lengthOfTrainingSeq < 100'000 || 
                                    ((taxonomy->getEukaryotaTaxID() != 0) && 
                                     (taxonomy->IsAncestor(accessionBatches[batchIdx].speciesID, taxonomy->getEukaryotaTaxID()))
                                    )) {
                                    prodigal->is_meta = 1;
                                    prodigal->trainMeta((unsigned char *) training_seq->entry.sequence.s, 
                                                        training_seq->entry.sequence.l);
                                } else {
                                    prodigal->trainASpecies((unsigned char *) training_seq->entry.sequence.s,
                                                            training_seq->entry.sequence.l);
                                }

                                // Generate intergenic 23-mer list. It is used to determine extension direction of intergenic sequences.
                                prodigal->getPredictedGenes((unsigned char *) training_seq->entry.sequence.s,
                                                            training_seq->entry.sequence.l);
                                seqIterator.generateIntergenicKmerList(prodigal->genes, prodigal->nodes,
                                                                       prodigal->getNumberOfPredictedGenes(),
                                                                       intergenicKmers, training_seq->entry.sequence.s);
                                // Get min k-mer hash list for determining strandness
                                seqIterator.getMinHashList(standardList, training_seq->entry.sequence.s);
                                delete training_seq;
                                trained = true;
                            }
                            currentList = priority_queue<uint64_t>();
                            seqIterator.getMinHashList(currentList, e.sequence.s);
                            if (seqIterator.compareMinHashList(standardList, currentList, lengthOfTrainingSeq, e.sequence.l)) {
                                // Get extended ORFs
                                prodigal->getPredictedGenes((unsigned char *) e.sequence.s, e.sequence.l);
                                prodigal->removeCompletelyOverlappingGenes();
                                prodigal->getExtendedORFs(prodigal->finalGenes, prodigal->nodes, extendedORFs,
                                                             prodigal->fng, e.sequence.l,
                                                        orfNum, intergenicKmers, e.sequence.s);
                                // Get k-mers from extended ORFs
                                for (size_t orfCnt = 0; orfCnt < orfNum; orfCnt++) {
                                    tempCheck = kmerExtractor->extractTargetKmers(
                                                    maskedSeq,
                                                    kmerBuffer,
                                                    posToWrite,
                                                    accessionBatches[batchIdx].taxIDs[idx],
                                                    accessionBatches[batchIdx].speciesID,
                                                    extendedORFs[orfCnt]);
                                    if (tempCheck == -1) {
                                        cout << "ERROR: Buffer overflow " << e.name.s << e.sequence.l << endl;
                                    }
                                }
                            } else { // Reverse complement
                                reverseComplement = seqIterator.reverseComplement(e.sequence.s, e.sequence.l);
                                // Get extended ORFs
                                prodigal->getPredictedGenes((unsigned char *) reverseComplement, e.sequence.l);
                                prodigal->removeCompletelyOverlappingGenes();
                                prodigal->getExtendedORFs(prodigal->finalGenes, prodigal->nodes, extendedORFs,
                                                                 prodigal->fng, e.sequence.l,
                                                            orfNum, intergenicKmers, reverseComplement);

                                // Get reverse masked sequence
                                if (par.maskMode) {
                                    delete[] maskedSeq;
                                    maskedSeq = new char[e.sequence.l + 1];
                                    SeqIterator::maskLowComplexityRegions((unsigned char *) reverseComplement, (unsigned char *) maskedSeq, probMatrix, par.maskProb, subMat);
                                    maskedSeq[e.sequence.l] = '\0';
                                } else {
                                    maskedSeq = reverseComplement;
                                }

                                for (size_t orfCnt = 0; orfCnt < orfNum; orfCnt++) {
                                    tempCheck = kmerExtractor->extractTargetKmers(
                                                    maskedSeq,
                                                    kmerBuffer,
                                                    posToWrite,
                                                    accessionBatches[batchIdx].taxIDs[idx],
                                                    accessionBatches[batchIdx].speciesID,
                                                    extendedORFs[orfCnt]);
                                    if (tempCheck == -1) {
                                        cout << "ERROR: Buffer overflow " << e.name.s << e.sequence.l << endl;
                                    }
                                }
                                free(reverseComplement);  
                            }                            
                        }
                        idx++;
                        if (par.maskMode) {
                            delete[] maskedSeq;
                        }
                        if (idx == accessionBatches[batchIdx].lengths.size()) {
                            break;
                        }
                    }
                    seqCnt++;
                }
                delete kseq;
                __sync_fetch_and_add(&processedBatchCnt, 1);
                #pragma omp critical
                {
                    cout << processedBatchCnt << " batches processed out of " << accessionBatches.size() << endl;
                        // cout << fastaPaths[accessionBatches[batchIdx].whichFasta] << " processed\n";
                }
            } else {
                batchChecker[batchIdx].store(false, std::memory_order_release);
                hasOverflow.fetch_add(1, std::memory_order_relaxed);
                __sync_fetch_and_sub(&kmerBuffer.startIndexOfReserve, estimatedKmerCnt);
            }
            delete prodigal;   
        }
    }

    // cout << "Before return: " << kmerBuffer.startIndexOfReserve << endl;
    return 0;
}


void IndexCreator::writeDbParameters() {
    FILE *handle = fopen(paramterFileName.c_str(), "w");
    if (handle == NULL) {
        Debug(Debug::ERROR) << "Could not open " << paramterFileName << " for writing\n";
        EXIT(EXIT_FAILURE);
    }
    fprintf(handle, "DB_name\t%s\n", par.dbName.c_str());
    fprintf(handle, "Creation_date\t%s\n", par.dbDate.c_str());
    fprintf(handle, "Metabuli commit used to create the DB\t%s\n", version);
    fprintf(handle, "Reduced_alphabet\t%d\n", par.reducedAA);
    // fprintf(handle, "Spaced_kmer_mask\t%s\n", spaceMask.c_str());
    fprintf(handle, "Accession_level\t%d\n", par.accessionLevel);
    fprintf(handle, "Mask_mode\t%d\n", par.maskMode);
    fprintf(handle, "Mask_prob\t%f\n", par.maskProb);
    fprintf(handle, "Skip_redundancy\t1\n");
    fprintf(handle, "Syncmer\t%d\n", par.syncmer);
    if (par.syncmer == 1) {
        fprintf(handle, "Syncmer_len\t%d\n", par.smerLen);
    }
    fprintf(handle, "Kmer_format\t%d\n", kmerFormat);
    fclose(handle);
}


void IndexCreator::loadCdsInfo(const string & cdsInfoFileList) {
    uint32_t prtId = 1;
    ifstream cdsInfoList(cdsInfoFileList);
    if (!cdsInfoList.is_open()) {
        Debug(Debug::ERROR) << "Could not open " << cdsInfoFileList << " for reading\n";
        EXIT(EXIT_FAILURE);
    }

    string cdsInfoFile;
    while (getline(cdsInfoList, cdsInfoFile)) {
        if (!FileUtil::fileExists(cdsInfoFile.c_str())) {
            Debug(Debug::ERROR) << "Could not open " << cdsInfoFile << " for reading\n";
            EXIT(EXIT_FAILURE);
        }
        KSeqWrapper* kseq = KSeqFactory(cdsInfoFile.c_str());
        while (kseq->ReadEntry()) {
            const KSeqWrapper::KSeqEntry & e = kseq->entry;
            string ename = e.name.s;
            size_t start = ename.find('|') + 1;
            size_t end = ename.find('.', start);
            string accession = ename.substr(start, end - start + 2);
            int frame = 1;
            string comment = e.comment.s;
            while (true) {
                start = comment.find('[', end) + 1;
                end = comment.find(']', start);
                if (start == string::npos) { break;}
                size_t equalPos = comment.find('=', start);
                string feature = comment.substr(start, equalPos - start);
                string value = comment.substr(equalPos + 1, end - equalPos - 1);
                if (feature == "pseudo") {
                    break;
                } else if (feature == "protein" && value == "hypothetical protein") {
                    break;
                } else if (feature == "frame") {
                    frame = stoi(value);
                } else if (feature == "protein_id") {
                    cdsInfoMap[accession].emplace_back(CDSinfo(prtId++, frame));
                } else if (feature == "location") {
                    size_t complementPos = value.find('c');
                    bool isComplement = (complementPos != string::npos);
                    if (isComplement) {
                        cdsInfoMap[accession].back().isComplement = true;
                        value = value.substr(complementPos + 11, value.size() - complementPos - 12);
                    } else {
                        cdsInfoMap[accession].back().isComplement = false;
                    }
                    size_t joinPos = value.find('j');
                    if (joinPos != string::npos) {
                        value = value.substr(joinPos + 5, value.size() - joinPos - 6);
                    }                
                    size_t commaPos = value.find(',');
                    size_t dotPos;
                    string locationBegin, locationEnd;
                    while (commaPos != string::npos) {
                        dotPos = value.find('.');
                        if (dotPos > commaPos) {
                            locationBegin = value.substr(0, commaPos);
                            locationEnd = value.substr(0, commaPos);
                        } else {
                            locationBegin = value.substr(0, dotPos);
                            locationEnd = value.substr(dotPos + 2, commaPos - dotPos - 2);
                        }
                        if (locationBegin[0] == '<') {
                            locationBegin = locationBegin.substr(1, locationBegin.size() - 1);
                        }
                        if (locationEnd[0] == '>') {
                            locationEnd = locationEnd.substr(1, locationEnd.size() - 1);
                        }
                        cdsInfoMap[accession].back().loc.emplace_back(stoi(locationBegin), stoi(locationEnd));
                        value = value.substr(commaPos + 1, value.size() - commaPos - 1);
                        commaPos = value.find(',');
                    }
                    dotPos = value.find('.');
                    if (dotPos == string::npos) {
                        locationBegin = value;
                        locationEnd = value;
                    } else {
                        locationBegin = value.substr(0, dotPos);
                        locationEnd = value.substr(dotPos + 2, commaPos - dotPos - 2);
                    }
                    // cout << value << endl;
                    // cout << locationBegin << endl;
                    // cout << locationEnd << endl;
                    if (locationBegin[0] == '<') {
                        locationBegin = locationBegin.substr(1, locationBegin.size() - 1);
                    }
                    if (locationEnd[0] == '>') {
                        locationEnd = locationEnd.substr(1, locationEnd.size() - 1);
                    }
                    cdsInfoMap[accession].back().loc.emplace_back(stoi(locationBegin), stoi(locationEnd));

                    if (frame != 1) {
                        if (!isComplement) {
                            cdsInfoMap[accession].back().loc[0].first += frame - 1;
                        } else {
                            cdsInfoMap[accession].back().loc.back().second -= frame - 1;
                        }
                    }
                    break;
                }
            }
        }
        delete kseq;
    }
    for (auto & cdsInfo : cdsInfoMap) {
        cout << "CDS info for " << cdsInfo.first << ": ";   
    }
    cdsInfoList.close();
}

void IndexCreator::loadMergedTaxIds(const std::string &mergedFile, unordered_map<TaxID, TaxID> & old2new) {
    std::ifstream ss(mergedFile);
    if (ss.fail()) {
        cout << "File " << mergedFile << " not found!\n";
        EXIT(EXIT_FAILURE);
    }

    std::string line;
    while (std::getline(ss, line)) {
        std::vector<std::string> result = splitByDelimiter(line, "\t|\t", 2);
        if (result.size() != 2) {
            Debug(Debug::ERROR) << "Invalid name entry!\n";
            EXIT(EXIT_FAILURE);
        }
        TaxID oldId = (TaxID) atoi(result[0].c_str());
        TaxID newId = (TaxID) atoi(result[1].c_str());
        old2new[oldId] = newId;
    }
}


void IndexCreator::addFilesToMerge(string diffIdxFileName, string infoFileName) {
    deltaIdxFileNames.push_back(diffIdxFileName);
    infoFileNames.push_back(infoFileName);
}

void IndexCreator::setMergedFileNames(string diffFileName, string infoFileName, string splitFileName) {
    mergedDeltaIdxFileName = diffFileName;
    mergedInfoFileName = infoFileName;
    deltaIdxSplitFileName = splitFileName;
}

void IndexCreator::updateTaxId2SpeciesTaxId(const string & taxIdListFileName) {
    FILE * taxIdFile;
    if((taxIdFile = fopen(taxIdListFileName.c_str(),"r")) == NULL){
        cout << "Cannot open the taxID list file: " << taxIdListFileName << endl;
        return;
    }

    char taxID[100];
    while(fscanf(taxIdFile,"%s",taxID) == 1) {
        TaxID taxId = atol(taxID);
        TaxonNode const * taxon = taxonomy->taxonNode(taxId);
        if (taxId == taxon->taxId){
            TaxID speciesTaxID = taxonomy->getTaxIdAtRank(taxId, "species");
            while (taxon->taxId != speciesTaxID) {
                taxId2speciesId[taxon->taxId] = speciesTaxID;
                taxon = taxonomy->taxonNode(taxon->parentTaxId);
            }
            taxId2speciesId[speciesTaxID] = speciesTaxID;
        } else { // merged
            TaxID speciesTaxID = taxonomy->getTaxIdAtRank(taxId, "species");
            while (taxon->taxId != speciesTaxID) {
                taxId2speciesId[taxon->taxId] = speciesTaxID;
                taxon = taxonomy->taxonNode(taxon->parentTaxId);
            }
            taxId2speciesId[speciesTaxID] = speciesTaxID;
            taxId2speciesId[taxId] = speciesTaxID;
        }
    }
    fclose(taxIdFile);
    Debug(Debug::INFO) << "Species-level taxonomy IDs are prepared.\n";
}



// void IndexCreator::editTaxonomyDumpFiles(const vector<pair<string, pair<TaxID, TaxID>>> & newAcc2taxid) {
//     // Load merged.dmp
//     string mergedFileName = taxonomyDir + "/merged.dmp";
//     std::ifstream ss(mergedFileName);
//     if (ss.fail()) {
//         Debug(Debug::ERROR) << "File " << mergedFileName << " not found!\n";
//         EXIT(EXIT_FAILURE);
//     }

//     std::string line;
//     unordered_map<int, int> mergedMap;
//     while (std::getline(ss, line)) {
//         std::vector<std::string> result = splitByDelimiter(line, "\t|\t", 2);
//         if (result.size() != 2) {
//             Debug(Debug::ERROR) << "Invalid name entry!\n";
//             EXIT(EXIT_FAILURE);
//         }
//         mergedMap[atoi(result[0].c_str())] = atoi(result[1].c_str());
//     }

//     // Edit names.dmp
//     string nameFileName = taxonomyDir + "/names.dmp";
//     string newNameFileName = taxonomyDir + "/names.dmp.new";
//     FileUtil::copyFile(nameFileName.c_str(), newNameFileName.c_str());
//     FILE *nameFile = fopen(newNameFileName.c_str(), "a");
//     if (nameFile == NULL) {
//         Debug(Debug::ERROR) << "Could not open " << newNameFileName << " for writing\n";
//         EXIT(EXIT_FAILURE);
//     }

//     for (size_t i = 0; i < newAcc2taxid.size() - 1; i++) {
//         fprintf(nameFile, "%d\t|\t%s\t|\t\t|\tscientific name\t|\n", newAcc2taxid[i].second.second, newAcc2taxid[i].first.c_str());
//     }
//     fprintf(nameFile, "%d\t|\t%s\t|\t\t|\tscientific name\t|", newAcc2taxid.back().second.second, newAcc2taxid.back().first.c_str());
//     fclose(nameFile);

//     // Edit nodes.dmp
//     string nodeFileName = taxonomyDir + "/nodes.dmp";
//     string newNodeFileName = taxonomyDir + "/nodes.dmp.new";
//     FileUtil::copyFile(nodeFileName.c_str(), newNodeFileName.c_str());
//     FILE *nodeFile = fopen(newNodeFileName.c_str(), "a");
//     if (nodeFile == NULL) {
//         Debug(Debug::ERROR) << "Could not open " << newNodeFileName << " for writing\n";
//         EXIT(EXIT_FAILURE);
//     }

//     for (size_t i = 0; i < newAcc2taxid.size() - 1; i++) {
//         // Check if the parent taxon is merged
//         if (mergedMap.find(newAcc2taxid[i].second.first) != mergedMap.end()) { // merged
//             fprintf(nodeFile, "%d\t|\t%d\t|\t\t|\tscientific name\t|\t\t|\t\t|\t\t|\t\t|\t\t|\t\t|\t\t|\n", newAcc2taxid[i].second.second, mergedMap[newAcc2taxid[i].second.first]);
//         } else {
//             fprintf(nodeFile, "%d\t|\t%d\t|\t\t|\tscientific name\t|\t\t|\t\t|\t\t|\t\t|\t\t|\t\t|\t\t|\n", newAcc2taxid[i].second.second, newAcc2taxid[i].second.first);
//         }
//         // fprintf(nodeFile, "%d\t|\t%d\t|\t\t|\tscientific name\t|\t\t|\t\t|\t\t|\t\t|\t\t|\t\t|\t\t|\n", newAcc2taxid[i].second.second, taxonomy->taxonNode(newAcc2taxid[i].second.first)->taxId);
//     }
//     // Check if the parent taxon is merged
//     if (mergedMap.find(newAcc2taxid.back().second.first) != mergedMap.end()) { // merged
//         fprintf(nodeFile, "%d\t|\t%d\t|\t\t|\tscientific name\t|\t\t|\t\t|\t\t|\t\t|\t\t|\t\t|\t\t|", newAcc2taxid.back().second.second, mergedMap[newAcc2taxid.back().second.first]);
//     } else {
//         fprintf(nodeFile, "%d\t|\t%d\t|\t\t|\tscientific name\t|\t\t|\t\t|\t\t|\t\t|\t\t|\t\t|\t\t|", newAcc2taxid.back().second.second, newAcc2taxid.back().second.first);
//     }
//     fclose(nodeFile);
// }


// TaxID IndexCreator::getMaxTaxID() {
//     // Check nodes.dmp
//     string nodeFileName = taxonomyDir + "/nodes.dmp";
//     std::ifstream ss(nodeFileName);
//     if (ss.fail()) {
//         Debug(Debug::ERROR) << "File " << nodeFileName << " not found!\n";
//         EXIT(EXIT_FAILURE);
//     }

//     std::string line;
//     TaxID maxTaxID = 0;
//     while (std::getline(ss, line)) {
//         std::vector<std::string> result = splitByDelimiter(line, "\t|\t", 2);
//         if (result.size() != 2) {
//             Debug(Debug::ERROR) << "Invalid name entry!\n";
//             EXIT(EXIT_FAILURE);
//         }
//         maxTaxID = std::max(maxTaxID, (TaxID) atoi(result[0].c_str()));
//     }
//     ss.close();

//     // Check names.dmp
//     string nameFileName = taxonomyDir + "/names.dmp";
//     ss = std::ifstream(nameFileName);
//     if (ss.fail()) {
//         Debug(Debug::ERROR) << "File " << nameFileName << " not found!\n";
//         EXIT(EXIT_FAILURE);
//     }

//     while (std::getline(ss, line)) {
//         std::vector<std::string> result = splitByDelimiter(line, "\t|\t", 2);
//         if (result.size() != 2) {
//             Debug(Debug::ERROR) << "Invalid name entry!\n";
//             EXIT(EXIT_FAILURE);
//         }
//         maxTaxID = std::max(maxTaxID, (TaxID) atoi(result[0].c_str()));
//     }
//     ss.close();

//     return maxTaxID;
// }
