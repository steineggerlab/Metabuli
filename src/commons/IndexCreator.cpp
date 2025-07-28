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
    : par(par), taxonomy(taxonomy), kmerFormat(kmerFormat) {
    dbDir = par.filenames[0];
    if (par.taxonomyPath.empty()) {
        taxonomyDir = dbDir + "/taxonomy/";
    } else {
        taxonomyDir = par.taxonomyPath + "/";
    }
    fnaListFileName = par.filenames[1];
    acc2taxidFileName = par.filenames[2];

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


IndexCreator::~IndexCreator() {
    delete geneticCode;
    delete kmerExtractor;
    delete subMat;
}

void IndexCreator::createIndex(const LocalParameters &par) {
    Buffer<TargetKmer> kmerBuffer(calculateBufferSize(par.ramUsage));
    cout << "Target metamer buffer size: " << kmerBuffer.bufferSize << endl;
    
    indexReferenceSequences(kmerBuffer.bufferSize);
    cout << "Made blocks for each thread" << endl;

    if (!par.cdsInfo.empty()) {
        loadCdsInfo(par.cdsInfo);
    }

    // Write taxonomy id list
    FILE * taxidListFile = fopen(taxidListFileName.c_str(), "w");
    for (auto & taxid: taxIdSet) {
        fprintf(taxidListFile, "%d\n", taxid);
    }
    fclose(taxidListFile);

    // Process the splits until all are processed
    size_t batchNum = accessionBatches.size();
    std::vector<std::atomic<bool>> batchChecker(batchNum);
    size_t processedBatchCnt = 0;
    
    vector<pair<size_t, size_t>> uniqKmerIdxRanges;
#ifdef OPENMP
    omp_set_num_threads(par.threads);
#endif
    while(processedBatchCnt < batchNum) {
        memset(kmerBuffer.buffer, 0, kmerBuffer.bufferSize * sizeof(TargetKmer));
        cout << "Buffer initialized" << endl;
        fillTargetKmerBuffer(kmerBuffer, batchChecker, processedBatchCnt, par);

        // Sort the k-mers
        time_t start = time(nullptr);
        SORT_PARALLEL(kmerBuffer.buffer, kmerBuffer.buffer + kmerBuffer.startIndexOfReserve,
                      IndexCreator::compareForDiffIdx);
        time_t sort = time(nullptr);
        cout << "Sort time: " << sort - start << endl;

        // Reduce redundancy
        auto * uniqKmerIdx = new size_t[kmerBuffer.startIndexOfReserve + 1];
        size_t uniqKmerCnt = 0;
        uniqKmerIdxRanges.clear();
        reduceRedundancy(kmerBuffer, uniqKmerIdx, uniqKmerCnt, uniqKmerIdxRanges);
        time_t reduction = time(nullptr);
        cout << "Time spent for reducing redundancy: " << (double) (reduction - sort) << endl;
        if(processedBatchCnt == batchNum && numOfFlush == 0 && !isUpdating) {
            writeTargetFilesAndSplits(kmerBuffer.buffer, kmerBuffer.startIndexOfReserve, uniqKmerIdx, uniqKmerCnt, uniqKmerIdxRanges);
        } else {
            writeTargetFiles(kmerBuffer.buffer, kmerBuffer.startIndexOfReserve, uniqKmerIdx, uniqKmerIdxRanges);
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
            char* pos = strchr(accession, '.');
            if (pos != nullptr) {
                *pos = '\0';
            }
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
    cout << "external2internalTaxID" << endl;
    for (size_t i = 0; i < observedAccessionsVec.size(); ++i) {    
        auto it = external2internalTaxID.find(observedAccessionsVec[i].taxID);
        if (it == external2internalTaxID.end() || observedAccessionsVec[i].taxID == 0) {
            // cout << "TaxID is not found for accession " << observedAccessionsVec[i].accession << " " << observedAccessionsVec[i].taxID << endl;
            unmappedAccessions.push_back(observedAccessionsVec[i].accession);
            continue;
        }
        observedAccessionsVec[i].taxID = it->second; // store the internal taxID
        observedAccessionsVec[i].speciesID = taxonomy->getTaxIdAtRank(it->second, "species");
        taxIdSet.insert(it->second);
        taxId2speciesId[it->second] = observedAccessionsVec[i].speciesID;
    }
    
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
    TargetKmer * kmerBuffer, 
    size_t & kmerNum,
    const size_t * uniqKmerIdx,
    const vector<pair<size_t, size_t>> & uniqKmerIdxRanges) 
{
    string diffIdxFileName;
    string infoFileName;
    FILE * diffIdxFile = nullptr;
    FILE * infoFile = nullptr;
    // if (isNewFormat) {
    //     diffIdxFileName = dbDir + "/" + to_string(numOfFlush) + "_deltaIdx.mtbl";
    //     diffIdxFile = fopen(diffIdxFileName.c_str(), "wb");
    //     deltaIdxFileNames.push_back(diffIdxFileName);
    // } else 
    // {
    diffIdxFileName = dbDir + "/" + to_string(numOfFlush) + "_diffIdx";
    infoFileName = dbDir + "/" + to_string(numOfFlush) + "_info";
    diffIdxFile = fopen(diffIdxFileName.c_str(), "wb");
    infoFile = fopen(infoFileName.c_str(), "wb");
    deltaIdxFileNames.push_back(diffIdxFileName);
    infoFileNames.push_back(infoFileName);
    // }
    
    if (diffIdxFile == nullptr || infoFile == nullptr) {
        cout<<"Cannot open the file for writing target DB"<<endl;
        return;
    }

    numOfFlush++;

    size_t bufferSize = 1024 * 1024 * 32;
    uint16_t *deltaIdxBuffer = (uint16_t *) malloc(sizeof(uint16_t) * bufferSize); // 64MB
    TaxID * infoBuffer = (TaxID *) malloc(sizeof(TaxID) * bufferSize); 
    size_t localBufIdx = 0;
    size_t localInfoBufIdx = 0;
    Metamer prevMetamer;
    size_t write = 0;

    vector<Metamer> identicalMetamers;
    uint64_t currentMetamer;
    for (size_t i = 0; i < uniqKmerIdxRanges.size(); i ++) {
        size_t j = uniqKmerIdxRanges[i].first;
        while (j < uniqKmerIdxRanges[i].second) {
            identicalMetamers.clear();
            currentMetamer = kmerBuffer[uniqKmerIdx[j]].metamer.metamer;
            while (j < uniqKmerIdxRanges[i].second && 
                   currentMetamer == kmerBuffer[uniqKmerIdx[j]].metamer.metamer) {
                identicalMetamers.push_back(kmerBuffer[uniqKmerIdx[j]].metamer);
                j++;
            }
            sort(identicalMetamers.begin(), identicalMetamers.end(), compareMetamerID);
            for (size_t k = 0; k < identicalMetamers.size(); k++) {
                writeInfo((int) identicalMetamers[k].id, infoFile, infoBuffer, bufferSize, localInfoBufIdx);
                getDiffIdx(prevMetamer.metamer, identicalMetamers[k].metamer, diffIdxFile, deltaIdxBuffer, bufferSize, localBufIdx);
                prevMetamer = identicalMetamers[k];
                write++;               
            }
        }
    }
    flushKmerBuf(deltaIdxBuffer, diffIdxFile, localBufIdx);
    free(deltaIdxBuffer);
    fclose(diffIdxFile);
    flushInfoBuf(infoBuffer, infoFile, localInfoBufIdx);
    free(infoBuffer);
    fclose(infoFile);
    cout<<"total k-mer count  : "<< kmerNum << endl;
    cout<<"Number of written k-mers: "<< write << endl;    
    kmerNum = 0;
}

void IndexCreator::writeTargetFilesAndSplits(
    TargetKmer * kmerBuffer,
    size_t & kmerNum,
    const size_t * uniqKmerIdx,
    size_t & uniqKmerCnt,
    const vector<pair<size_t, size_t>> & uniqKmerIdxRanges)
{
    string diffIdxFileName;
    string infoFileName;
    string splitFileName;
    FILE * diffIdxFile = nullptr;
    FILE * infoFile = nullptr;
    FILE * deltaIdxSplitFile = nullptr;

    diffIdxFileName = dbDir + "/diffIdx";
    splitFileName = dbDir + "/split";
    infoFileName = dbDir + "/info";
    infoFile = fopen(infoFileName.c_str(), "wb");
    diffIdxFile = fopen(diffIdxFileName.c_str(), "wb");
    deltaIdxSplitFile = fopen(splitFileName.c_str(), "wb");
    if (diffIdxFile == nullptr || deltaIdxSplitFile == nullptr || infoFile == nullptr) {
        cout << "Cannot open the file for writing target DB" << endl;
        return;
    }

    // DeltaIdxOffset * offsetList = nullptr;
    DiffIdxSplit * splitList = nullptr;
    // if (isNewFormat) {
    //     offsetList = new DeltaIdxOffset[par.splitNum];
    //     memset(offsetList, 0, sizeof(DeltaIdxOffset) * par.splitNum);
    // } else {
    splitList = new DiffIdxSplit[par.splitNum];
    memset(splitList, 0, sizeof(DiffIdxSplit) * par.splitNum);
    // }
    
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
                if (AminoAcidPart(kmerBuffer[uniqKmerIdx[k]].metamer.metamer) 
                    != AminoAcidPart(kmerBuffer[uniqKmerIdx[k + 1]].metamer.metamer)) {
                    // if (isNewFormat) {
                    //     offsetList[splitCnt++].metamer = kmerBuffer[uniqKmerIdx[k + 1]].metamer;
                    // } else {
                        splitList[splitCnt++].ADkmer = kmerBuffer[uniqKmerIdx[k + 1]].metamer.metamer;
                    // }
                    found = true;
                    break;
                }
            }
            if (!found) {
                // if (isNewFormat) {
                //     offsetList[splitCnt++].metamer = kmerBuffer[uniqKmerIdx[uniqKmerIdxRanges[j+1].first]].metamer;
                // } else {
                    splitList[splitCnt++].ADkmer = kmerBuffer[uniqKmerIdx[uniqKmerIdxRanges[j+1].first]].metamer.metamer;
                // }
                // cout << "Split " << splitCnt << " at " << offsetList[splitCnt-1].metamer.metamer << endl;
            }
            break;
        }
    }

    numOfFlush++;
    size_t bufferSize = 1024 * 1024 * 32;
    uint16_t *deltaIdxBuffer = (uint16_t *)malloc(sizeof(uint16_t) * bufferSize); // 64MB
    TaxID * infoBuffer = (TaxID *)malloc(sizeof(TaxID) * bufferSize);
    // if (!isNewFormat) {
    //     infoBuffer = (TaxID *)malloc(sizeof(TaxID) * bufferSize); // 128MB
    // }
    size_t localBufIdx = 0;
    size_t localInfoBufIdx = 0;
    Metamer prevMetamer;
    size_t splitIdx = 1;
    size_t totalDiffIdx = 0;
    size_t write = 0;

    vector<Metamer> identicalMetamers;
    uint64_t currentMetamer;

    cout << "Writing k-mers to disk" << endl;
    for (size_t i = 0; i < uniqKmerIdxRanges.size(); i ++) {
        size_t j = uniqKmerIdxRanges[i].first;
        while (j < uniqKmerIdxRanges[i].second) {
            identicalMetamers.clear();
            currentMetamer = kmerBuffer[uniqKmerIdx[j]].metamer.metamer;
            while (j < uniqKmerIdxRanges[i].second && 
                   currentMetamer == kmerBuffer[uniqKmerIdx[j]].metamer.metamer) {
                identicalMetamers.push_back(kmerBuffer[uniqKmerIdx[j]].metamer);
                j++;
            }
            // Sort the identicalMetamers using compareMetamer
            sort(identicalMetamers.begin(), identicalMetamers.end(), compareMetamerID);
            for (size_t k = 0; k < identicalMetamers.size(); k++) {
                // if (isNewFormat) {
                //     write++;
                //     getDeltaIdx(prevMetamer, identicalMetamers[k], diffIdxFile, deltaIdxBuffer, bufferSize, localBufIdx, totalDiffIdx);
                //     prevMetamer = identicalMetamers[k];
                //     if ((splitIdx < splitCnt) && (prevMetamer.metamer == offsetList[splitIdx].metamer.metamer)) {
                //         offsetList[splitIdx].metamer = prevMetamer;
                //         offsetList[splitIdx].offset = totalDiffIdx;
                //         splitIdx ++;
                //     }
                // } else {
                    writeInfo((int) identicalMetamers[k].id, infoFile, infoBuffer, bufferSize, localInfoBufIdx);
                    write++;
                    getDiffIdx(prevMetamer.metamer, identicalMetamers[k].metamer, diffIdxFile, deltaIdxBuffer, bufferSize, localBufIdx, totalDiffIdx);
                    prevMetamer = identicalMetamers[k];
                    if((splitIdx < splitCnt) && (prevMetamer.metamer == splitList[splitIdx].ADkmer)){
                        splitList[splitIdx].diffIdxOffset = totalDiffIdx;
                        splitList[splitIdx].infoIdxOffset = write;
                        splitIdx ++;
                    }
                // }
            }
        }
    }
    cout<<"total k-mer count  : "<< kmerNum << endl;
    cout<<"written k-mer count: "<< write << endl;

    flushKmerBuf(deltaIdxBuffer, diffIdxFile, localBufIdx);
    free(deltaIdxBuffer);
    // if (isNewFormat) {
    //     fwrite(offsetList, sizeof(DeltaIdxOffset), par.splitNum, deltaIdxSplitFile);
    //     delete[] offsetList;
    // } else {
        fwrite(splitList, sizeof(DiffIdxSplit), par.splitNum, deltaIdxSplitFile);
        delete[] splitList;
        flushInfoBuf(infoBuffer, infoFile, localInfoBufIdx);
        free(infoBuffer);
        fclose(infoFile);
    // }
    fclose(deltaIdxSplitFile);
    fclose(diffIdxFile);
    kmerNum = 0;
    cout << "Finished writing k-mers to disk" << endl;
}

void IndexCreator::reduceRedundancy(
    Buffer<TargetKmer> & kmerBuffer,
    size_t * uniqKmerIdx,
    size_t & uniqueKmerCnt,
    vector<pair<size_t, size_t>> & uniqKmerIdxRanges) 
{
    // Find the first index of garbage k-mer (UINT64_MAX)
    for(size_t checkN = kmerBuffer.startIndexOfReserve - 1; checkN != 0; checkN--){
        if(kmerBuffer.buffer[checkN].metamer.metamer != UINT64_MAX){
            kmerBuffer.startIndexOfReserve = checkN + 1;
            break;
        }
    }

    // Find the first index of meaningful k-mer
    size_t startIdx = 0;
    for (size_t i = 0; i < kmerBuffer.startIndexOfReserve ; i++) {
        if(kmerBuffer.buffer[i].spTaxId != 0){
            startIdx = i;
            break;
        }
    }

    for (size_t i = startIdx; i < kmerBuffer.startIndexOfReserve ; i++) {
        if (kmerBuffer.buffer[i].metamer.id == 0) {
            cout << "Error: k-mer with ID 0 found at index " << i << endl;
            exit(1);
        }

    }

    // cout << "Number of k-mers before redundancy reduction: " << kmerBuffer.startIndexOfReserve - startIdx << endl;

    // Make splits
    vector<Split> splits;
    size_t splitWidth = (kmerBuffer.startIndexOfReserve - startIdx) / par.threads;
    for (int i = 0; i < par.threads - 1; i++) {
        for (size_t j = startIdx + splitWidth; j + 1 < kmerBuffer.startIndexOfReserve; j++) {
            if (AminoAcidPart(kmerBuffer.buffer[j].metamer.metamer) != AminoAcidPart(kmerBuffer.buffer[j + 1].metamer.metamer)) {
                splits.emplace_back(startIdx, j);
                uniqKmerIdxRanges.emplace_back(pair<size_t, size_t>(startIdx, 0));
                startIdx = j + 1;
                break;
            }
        }
    }
    splits.emplace_back(startIdx, kmerBuffer.startIndexOfReserve - 1);
    uniqKmerIdxRanges.emplace_back(pair<size_t, size_t>(startIdx, 0));

    size_t * cntOfEachSplit = new size_t[splits.size()];
    for(size_t i = 0; i < splits.size(); i++) {
        cntOfEachSplit[i] = 0;
    }
#pragma omp parallel default(none), shared(kmerBuffer, cntOfEachSplit, splits, par, cout, uniqueKmerCnt, uniqKmerIdx, uniqKmerIdxRanges)
    {
        TargetKmer * lookingKmer;
        size_t lookingIndex;
        int endFlag;
        vector<TaxID> taxIds;
        size_t * tempUniqKmerIdx = new size_t[16 * 1024 * 1024];
        size_t tempUniqKmerCnt = 0;
#pragma omp for schedule(static, 1)
        for(size_t split = 0; split < splits.size(); split ++) {
            lookingKmer = & kmerBuffer.buffer[splits[split].offset];
            lookingIndex = splits[split].offset;
            endFlag = 0;
            for(size_t i = 1 + splits[split].offset; i < splits[split].end + 1 ; i++) {
                taxIds.clear();
                taxIds.push_back(lookingKmer->metamer.id);
                // Scan redundancy
                while(lookingKmer->spTaxId == kmerBuffer.buffer[i].spTaxId &&
                      lookingKmer->metamer.metamer == kmerBuffer.buffer[i].metamer.metamer){
                    taxIds.push_back(kmerBuffer.buffer[i].metamer.id);
                    i++;
                    if(i == splits[split].end + 1){
                        endFlag = 1;
                        break;
                    }
                }
                if(taxIds.size() > 1){
                    lookingKmer->metamer.id = taxonomy->LCA(taxIds)->taxId;
                } else {
                    lookingKmer->metamer.id = taxIds[0];
                }

                if (tempUniqKmerCnt >= 16 * 1024 * 1024) {
                    memcpy(uniqKmerIdx + splits[split].offset, tempUniqKmerIdx, tempUniqKmerCnt * sizeof(size_t));
                    splits[split].offset += tempUniqKmerCnt;
                    tempUniqKmerCnt = 0;
                }
                tempUniqKmerIdx[tempUniqKmerCnt++] = lookingIndex;
                cntOfEachSplit[split]++;
                if(endFlag == 1) break;
                lookingKmer = & kmerBuffer.buffer[i];
                lookingIndex = i;
            }

            //For the end part
            if(!((kmerBuffer.buffer[splits[split].end - 1].metamer.metamer == kmerBuffer.buffer[splits[split].end].metamer.metamer) &&
                 (kmerBuffer.buffer[splits[split].end - 1].spTaxId == kmerBuffer.buffer[splits[split].end].spTaxId))){
                if (tempUniqKmerCnt >= 16 * 1024 * 1024) {
                    memcpy(uniqKmerIdx + splits[split].offset, tempUniqKmerIdx, tempUniqKmerCnt * sizeof(size_t));
                    splits[split].offset += tempUniqKmerCnt;
                    tempUniqKmerCnt = 0;
                }
                tempUniqKmerIdx[tempUniqKmerCnt++] = splits[split].end;
                cntOfEachSplit[split]++;
            }
            memcpy(uniqKmerIdx + splits[split].offset, tempUniqKmerIdx, tempUniqKmerCnt * sizeof(size_t));
            splits[split].offset += tempUniqKmerCnt;
            uniqKmerIdxRanges[split].second = splits[split].offset;
            __sync_fetch_and_add(&uniqueKmerCnt, cntOfEachSplit[split]); 
        }
        delete[] tempUniqKmerIdx;
    }
    delete[] cntOfEachSplit;
}

void IndexCreator::getDeltaIdx(const Metamer & previousMetamer,
                               const Metamer & currentMetamer,
                               FILE* handleKmerTable,
                               uint16_t * deltaIndexBuffer,
                               size_t bufferSize,
                               size_t & localBufIdx,
                               size_t & totalBufferIdx) {
    bitset<96> diff = Metamer::substract(currentMetamer, previousMetamer);                       
    uint16_t buffer[7];
    int idx = 5;
    buffer[6] = SET_END_FLAG(static_cast<uint16_t>((diff & bitset<96>(0x7FFF)).to_ulong()));
    diff >>= 15U;
    while (diff.any()) {
        uint16_t toWrite = GET_15_BITS(static_cast<uint16_t>((diff & bitset<96>(0x7FFF)).to_ulong()));
        diff >>= 15U;
        buffer[idx] = toWrite;
        idx--;
    }
    totalBufferIdx += 6 - idx;
    writeDiffIdx(deltaIndexBuffer, bufferSize, handleKmerTable, (buffer + idx + 1), (6 - idx), localBufIdx);
}

void IndexCreator::getDeltaIdx(const Metamer & previousMetamer,
                               const Metamer & currentMetamer,
                               FILE* handleKmerTable,
                               uint16_t * deltaIndexBuffer,
                               size_t bufferSize,
                               size_t & localBufIdx) {
    bitset<96> diff = Metamer::substract(currentMetamer, previousMetamer);                               
    uint16_t buffer[7];
    int idx = 5;
    buffer[6] = SET_END_FLAG(static_cast<uint16_t>((diff & bitset<96>(0x7FFF)).to_ulong()));
    diff >>= 15U;
    while (diff.any()) {
        uint16_t toWrite = GET_15_BITS(static_cast<uint16_t>((diff & bitset<96>(0x7FFF)).to_ulong()));
        diff >>= 15U;
        buffer[idx] = toWrite;
        idx--;
    }
    writeDiffIdx(deltaIndexBuffer, bufferSize, handleKmerTable, (buffer + idx + 1), (6 - idx), localBufIdx);
}

void IndexCreator::getDiffIdx(const uint64_t & lastKmer,
                              const uint64_t & entryToWrite,
                              FILE* handleKmerTable,
                              uint16_t *kmerBuf, 
                              size_t bufferSize, 
                              size_t & localBufIdx, 
                              size_t & totalBufferIdx) {
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
    totalBufferIdx += 4 - idx;
    writeDiffIdx(kmerBuf, bufferSize, handleKmerTable, (buffer + idx + 1), (4 - idx), localBufIdx);
}

void IndexCreator::getDiffIdx(const uint64_t & lastKmer,
                              const uint64_t & entryToWrite, 
                              FILE* handleKmerTable, 
                              uint16_t *kmerBuf,
                              size_t bufferSize, 
                              size_t & localBufIdx) {
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
    writeDiffIdx(kmerBuf, bufferSize, handleKmerTable, (buffer + idx + 1), (4 - idx), localBufIdx);
}

void IndexCreator::flushKmerBuf(uint16_t *buffer, FILE *handleKmerTable, size_t & localBufIdx ) {
    fwrite(buffer, sizeof(uint16_t), localBufIdx, handleKmerTable);
    localBufIdx = 0;
}

void IndexCreator::writeDiffIdx(uint16_t *buffer,
                                size_t bufferMaxSize,
                                FILE* handleKmerTable,
                                uint16_t *toWrite,
                                size_t size,
                                size_t & localBufIdx ) {
    if (localBufIdx + size >= bufferMaxSize) {
        flushKmerBuf(buffer, handleKmerTable, localBufIdx);
    }
    memcpy(buffer + localBufIdx, toWrite, sizeof(uint16_t) * size);
    localBufIdx += size;
}

void IndexCreator::writeInfo(
    TaxID entryToWrite, 
    FILE * infoFile, 
    TaxID * infoBuffer, 
    size_t bufferSize, 
    size_t & infoBufferIdx) 
{
    if (infoBufferIdx >= bufferSize) 
    {
        flushInfoBuf(infoBuffer, infoFile, infoBufferIdx);
    }
    infoBuffer[infoBufferIdx++] = entryToWrite;
}

void IndexCreator::flushInfoBuf(TaxID * buffer, FILE * infoFile, size_t & localBufIdx) {
    fwrite(buffer, sizeof(TaxID), localBufIdx, infoFile);
    localBufIdx = 0;
}

int IndexCreator::getNumOfFlush() {return numOfFlush;}

inline bool IndexCreator::compareForDiffIdx(const TargetKmer & a, const TargetKmer & b){
    if (a.metamer.metamer != b.metamer.metamer) {
        return a.metamer.metamer < b.metamer.metamer;
    }

    return a.spTaxId < b.spTaxId;

    // if (a.spTaxId != b.spTaxId) {
    //     return a.spTaxId < b.spTaxId;
    // }

    // return a.metamer.id < b.metamer.id;
}

inline bool IndexCreator::compareMetamerID(const Metamer & a, const Metamer & b) {
    // if (a.metamer != b.metamer) {
    //     return a.metamer < b.metamer;
    // }
    return a.id < b.id;
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

size_t IndexCreator::fillTargetKmerBuffer(Buffer<TargetKmer> &kmerBuffer,
                                          std::vector<std::atomic<bool>> & batchChecker,
                                          size_t &processedBatchCnt,
                                          const LocalParameters &par) {
    std::atomic<int> hasOverflow{0};
#pragma omp parallel default(none), shared(kmerBuffer, batchChecker, processedBatchCnt, hasOverflow, par, cout)
    {
        ProbabilityMatrix probMatrix(*subMat);
        // ProdigalWrapper * prodigal = new ProdigalWrapper();
        SeqIterator seqIterator(par);
        size_t posToWrite;
        size_t orfNum;
        vector<SequenceBlock> extendedORFs;
        priority_queue<uint64_t> standardList;
        priority_queue<uint64_t> currentList;
        size_t lengthOfTrainingSeq;
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

    cout << "Before return: " << kmerBuffer.startIndexOfReserve << endl;
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

void IndexCreator::editTaxonomyDumpFiles(const vector<pair<string, pair<TaxID, TaxID>>> & newAcc2taxid) {
    // Load merged.dmp
    string mergedFileName = taxonomyDir + "/merged.dmp";
    std::ifstream ss(mergedFileName);
    if (ss.fail()) {
        Debug(Debug::ERROR) << "File " << mergedFileName << " not found!\n";
        EXIT(EXIT_FAILURE);
    }

    std::string line;
    unordered_map<int, int> mergedMap;
    while (std::getline(ss, line)) {
        std::vector<std::string> result = splitByDelimiter(line, "\t|\t", 2);
        if (result.size() != 2) {
            Debug(Debug::ERROR) << "Invalid name entry!\n";
            EXIT(EXIT_FAILURE);
        }
        mergedMap[atoi(result[0].c_str())] = atoi(result[1].c_str());
    }

    // Edit names.dmp
    string nameFileName = taxonomyDir + "/names.dmp";
    string newNameFileName = taxonomyDir + "/names.dmp.new";
    FileUtil::copyFile(nameFileName.c_str(), newNameFileName.c_str());
    FILE *nameFile = fopen(newNameFileName.c_str(), "a");
    if (nameFile == NULL) {
        Debug(Debug::ERROR) << "Could not open " << newNameFileName << " for writing\n";
        EXIT(EXIT_FAILURE);
    }

    for (size_t i = 0; i < newAcc2taxid.size() - 1; i++) {
        fprintf(nameFile, "%d\t|\t%s\t|\t\t|\tscientific name\t|\n", newAcc2taxid[i].second.second, newAcc2taxid[i].first.c_str());
    }
    fprintf(nameFile, "%d\t|\t%s\t|\t\t|\tscientific name\t|", newAcc2taxid.back().second.second, newAcc2taxid.back().first.c_str());
    fclose(nameFile);

    // Edit nodes.dmp
    string nodeFileName = taxonomyDir + "/nodes.dmp";
    string newNodeFileName = taxonomyDir + "/nodes.dmp.new";
    FileUtil::copyFile(nodeFileName.c_str(), newNodeFileName.c_str());
    FILE *nodeFile = fopen(newNodeFileName.c_str(), "a");
    if (nodeFile == NULL) {
        Debug(Debug::ERROR) << "Could not open " << newNodeFileName << " for writing\n";
        EXIT(EXIT_FAILURE);
    }

    for (size_t i = 0; i < newAcc2taxid.size() - 1; i++) {
        // Check if the parent taxon is merged
        if (mergedMap.find(newAcc2taxid[i].second.first) != mergedMap.end()) { // merged
            fprintf(nodeFile, "%d\t|\t%d\t|\t\t|\tscientific name\t|\t\t|\t\t|\t\t|\t\t|\t\t|\t\t|\t\t|\n", newAcc2taxid[i].second.second, mergedMap[newAcc2taxid[i].second.first]);
        } else {
            fprintf(nodeFile, "%d\t|\t%d\t|\t\t|\tscientific name\t|\t\t|\t\t|\t\t|\t\t|\t\t|\t\t|\t\t|\n", newAcc2taxid[i].second.second, newAcc2taxid[i].second.first);
        }
        // fprintf(nodeFile, "%d\t|\t%d\t|\t\t|\tscientific name\t|\t\t|\t\t|\t\t|\t\t|\t\t|\t\t|\t\t|\n", newAcc2taxid[i].second.second, taxonomy->taxonNode(newAcc2taxid[i].second.first)->taxId);
    }
    // Check if the parent taxon is merged
    if (mergedMap.find(newAcc2taxid.back().second.first) != mergedMap.end()) { // merged
        fprintf(nodeFile, "%d\t|\t%d\t|\t\t|\tscientific name\t|\t\t|\t\t|\t\t|\t\t|\t\t|\t\t|\t\t|", newAcc2taxid.back().second.second, mergedMap[newAcc2taxid.back().second.first]);
    } else {
        fprintf(nodeFile, "%d\t|\t%d\t|\t\t|\tscientific name\t|\t\t|\t\t|\t\t|\t\t|\t\t|\t\t|\t\t|", newAcc2taxid.back().second.second, newAcc2taxid.back().second.first);
    }
    fclose(nodeFile);
}

TaxID IndexCreator::getMaxTaxID() {
    // Check nodes.dmp
    string nodeFileName = taxonomyDir + "/nodes.dmp";
    std::ifstream ss(nodeFileName);
    if (ss.fail()) {
        Debug(Debug::ERROR) << "File " << nodeFileName << " not found!\n";
        EXIT(EXIT_FAILURE);
    }

    std::string line;
    TaxID maxTaxID = 0;
    while (std::getline(ss, line)) {
        std::vector<std::string> result = splitByDelimiter(line, "\t|\t", 2);
        if (result.size() != 2) {
            Debug(Debug::ERROR) << "Invalid name entry!\n";
            EXIT(EXIT_FAILURE);
        }
        maxTaxID = std::max(maxTaxID, (TaxID) atoi(result[0].c_str()));
    }
    ss.close();

    // Check names.dmp
    string nameFileName = taxonomyDir + "/names.dmp";
    ss = std::ifstream(nameFileName);
    if (ss.fail()) {
        Debug(Debug::ERROR) << "File " << nameFileName << " not found!\n";
        EXIT(EXIT_FAILURE);
    }

    while (std::getline(ss, line)) {
        std::vector<std::string> result = splitByDelimiter(line, "\t|\t", 2);
        if (result.size() != 2) {
            Debug(Debug::ERROR) << "Invalid name entry!\n";
            EXIT(EXIT_FAILURE);
        }
        maxTaxID = std::max(maxTaxID, (TaxID) atoi(result[0].c_str()));
    }
    ss.close();

    return maxTaxID;
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
        KSeqWrapper* kseq = KSeqFactory(cdsInfoFileList.c_str());
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

void IndexCreator::mergeTargetFiles() {
    // Files to write
    FILE * mergedDiffFile = fopen(mergedDeltaIdxFileName.c_str(), "wb");
    FILE * mergedInfoFile = fopen(mergedInfoFileName.c_str(), "wb");
    FILE * diffIdxSplitFile = fopen(deltaIdxSplitFileName.c_str(), "wb");

    // Buffers to fill
    size_t bufferSize = 1024 * 1024 * 512;
    uint16_t * diffBuffer = (uint16_t *)malloc(sizeof(uint16_t) * bufferSize);
    size_t diffBufferIdx = 0;
    size_t totalBufferIdx = 0;
    TaxID * infoBuffer = (TaxID *)malloc(sizeof(TaxID) * bufferSize);
    size_t infoBufferIdx = 0;

    // Prepare files to merge
    size_t numOfKmerBeforeMerge = 0;
    size_t numOfSplits = deltaIdxFileNames.size();
    DeltaIdxReader ** deltaIdxReaders = new DeltaIdxReader*[numOfSplits];
    size_t valueBufferSize = 1024 * 1024 * 32;
    for (size_t file = 0; file < numOfSplits; file++) {
        deltaIdxReaders[file] = new DeltaIdxReader(deltaIdxFileNames[file],
                                                   infoFileNames[file],
                                                   valueBufferSize, 
                                                   1024 * 1024 * 8);
        numOfKmerBeforeMerge += deltaIdxReaders[file]->getTotalValueNum();
    }

    // To make differential index splits
    uint64_t AAofTempSplitOffset = UINT64_MAX;
    size_t sizeOfSplit = numOfKmerBeforeMerge / (par.splitNum - 1);
    size_t offsetList[par.splitNum + 1];
    int offsetListIdx = 1;
    for(int os = 0; os < par.splitNum; os++){
        offsetList[os] = os * sizeOfSplit;
    }
    offsetList[par.splitNum] = UINT64_MAX;
    DiffIdxSplit splitList[par.splitNum];
    memset(splitList, 0, sizeof(DiffIdxSplit) * par.splitNum);
    int splitListIdx = 1;

    // get the first k-mer to write
    unsigned int mask = ~((static_cast<unsigned int>(par.skipRedundancy == 0) << 31));
    int splitCheck = 0;
    vector<TaxID> taxIds;

    cout << "Merging target files..." << endl;
    Buffer<TargetKmer> metamerBuffer(1024 * 1024 * 1024);
    std::atomic<int> hasOverflow{0};
    std::vector<std::atomic<bool>> completedSplits(numOfSplits);
    int remainingSplits = numOfSplits;
    vector<pair<size_t, size_t>> uniqKmerIdxRanges;
    size_t write = 0;
    size_t lastKmer = 0;
    size_t totalValueNum = 0;

    while (remainingSplits > 0) {
        memset(metamerBuffer.buffer, 0, metamerBuffer.bufferSize * sizeof(TargetKmer));
        time_t start = time(nullptr);
        while (true) {
            size_t posToWrite;
            uint64_t max = UINT64_MAX;
            posToWrite = metamerBuffer.reserveMemory(valueBufferSize * numOfSplits);
            if (posToWrite + valueBufferSize * numOfSplits > metamerBuffer.bufferSize) {
                metamerBuffer.startIndexOfReserve -= valueBufferSize * numOfSplits;
                break;
            }
            for (size_t i = 0; i < numOfSplits; i++) {
                if (completedSplits[i].load(std::memory_order_acquire))
                    continue; 
                if (deltaIdxReaders[i]->getLastValue() < max) {
                    max = deltaIdxReaders[i]->getLastValue();
                }
            }
#pragma omp parallel for default(none), shared(cout, metamerBuffer, deltaIdxReaders, completedSplits, posToWrite, max, numOfSplits, valueBufferSize, taxId2speciesId, mask, remainingSplits, totalValueNum)
            for (size_t split = 0; split < numOfSplits; split++) {
                if (completedSplits[split].load(std::memory_order_acquire))
                    continue; 
                size_t valueNum = deltaIdxReaders[split]->getValues(metamerBuffer.buffer + posToWrite + split * valueBufferSize, max);
                totalValueNum += valueNum;
                if (deltaIdxReaders[split]->isCompleted()) {
                    completedSplits[split].store(true, std::memory_order_release);
                    __sync_fetch_and_sub(&remainingSplits, 1);
                    // continue;
                }
                for (size_t i = 0; i < valueNum; i++) {
                    if (metamerBuffer.buffer[posToWrite + split * valueBufferSize + i].metamer.id == 0) {
                        cout << "valueNum: " << valueNum << endl;
                        cout << "split: " << split << endl;
                        cout << "posToWrite: " << posToWrite << endl;
                        cout << posToWrite + split * valueBufferSize + i << endl;
                        cout << "speciesId: " << taxId2speciesId[metamerBuffer.buffer[posToWrite + split * valueBufferSize + i].metamer.id] << endl;


                    }
                    metamerBuffer.buffer[posToWrite + split * valueBufferSize + i].spTaxId 
                        = taxId2speciesId[metamerBuffer.buffer[posToWrite + split * valueBufferSize + i].metamer.id & mask];
                }    
            }
        }
        time_t end = time(nullptr);
        cout << "K-mer loading       : " << (double) (end - start) << " s" << endl;

        time_t beforeSort = time(nullptr);
        SORT_PARALLEL(metamerBuffer.buffer, metamerBuffer.buffer + metamerBuffer.startIndexOfReserve,
                      IndexCreator::compareForDiffIdx);
        time_t afterSort = time(nullptr);
        cout << "Sorting k-mer list  : " << afterSort - beforeSort << " s" << endl;

        // Reduce redundancy
        auto * uniqKmerIdx = new size_t[metamerBuffer.startIndexOfReserve + 1];
        size_t uniqKmerCnt = 0;
        uniqKmerIdxRanges.clear();
        reduceRedundancy(metamerBuffer, uniqKmerIdx, uniqKmerCnt, uniqKmerIdxRanges);
        time_t reduction = time(nullptr);
        cout << "Redundancy reduction: " << (double) (reduction - afterSort) << " s" << endl;

        // Write       
        for (size_t i = 0; i < uniqKmerIdxRanges.size(); i ++) {
            for (size_t j = uniqKmerIdxRanges[i].first; j < uniqKmerIdxRanges[i].second; j ++) {
                writeInfo((int) metamerBuffer.buffer[uniqKmerIdx[j]].metamer.id,
                          mergedInfoFile, 
                          infoBuffer, 
                          bufferSize, 
                          infoBufferIdx);
                write++;
                getDiffIdx(lastKmer,
                           metamerBuffer.buffer[uniqKmerIdx[j]].metamer.metamer, 
                           mergedDiffFile,
                           diffBuffer,
                           bufferSize,
                           diffBufferIdx,
                           totalBufferIdx);
                lastKmer = metamerBuffer.buffer[uniqKmerIdx[j]].metamer.metamer;

                if (AminoAcidPart(lastKmer) != AAofTempSplitOffset && splitCheck == 1) {
                    splitList[splitListIdx++] = {lastKmer, totalBufferIdx, write};
                    splitCheck = 0;
                }
                if (write == offsetList[offsetListIdx]) {
                    AAofTempSplitOffset = AminoAcidPart(lastKmer);
                    splitCheck = 1;
                    offsetListIdx++;
                }
            }
        }   
        flushInfoBuf(infoBuffer, mergedInfoFile, infoBufferIdx);
        flushKmerBuf(diffBuffer, mergedDiffFile, diffBufferIdx);
        time_t writeTime = time(nullptr);
        cout << "Writing k-mers      : " << (double) (writeTime - reduction) << " s" << endl;
        cout << "Written k-mers      : " << write << " " << std::fixed << std::setprecision(2) << ((float) write / numOfKmerBeforeMerge) * 100 << "%" << endl;
        cout << "--------------------" << endl;
        cout << "remaining splits   : " << remainingSplits << endl;
        std::cout.unsetf(std::ios::fixed | std::ios::scientific);
        metamerBuffer.startIndexOfReserve = 0;
    }


    fwrite(splitList, sizeof(DiffIdxSplit), par.splitNum, diffIdxSplitFile);
    for(int i = 0; i < par.splitNum; i++) {
        cout<<splitList[i].ADkmer<< " "<<splitList[i].diffIdxOffset<< " "<<splitList[i].infoIdxOffset<<endl;
    }
    free(diffBuffer);
    free(infoBuffer);
    fclose(mergedDiffFile);
    fclose(mergedInfoFile);
    fclose(diffIdxSplitFile);

    cout<<"Creating target DB is done"<<endl;
    cout<<"Total k-mer count    : " << write <<endl;
    cout<<"done"<<endl;

    cout<<"Reference DB files you need are as below"<<endl;
    cout<<mergedDeltaIdxFileName<<endl;
    cout<<mergedInfoFileName<<endl;
    cout<<string(dbDir) + "/taxID_list"<<endl;
    cout<<deltaIdxSplitFileName<<endl;
}
