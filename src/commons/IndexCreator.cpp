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

IndexCreator::IndexCreator(const LocalParameters & par) : par(par) {
   // Parameters
    threadNum = par.threads;
    reducedAA = par.reducedAA;
    accessionLevel = par.accessionLevel;
    lowComplexityMasking = par.maskMode;
    lowComplexityMaskingThreshold = par.maskProb;
    dbName = par.dbName;
    dbDate = par.dbDate;
    

    // Input files
    dbDir = par.filenames[0];
    if (par.taxonomyPath.empty()) {
        taxonomyDir = dbDir + "/taxonomy/";
    } else {
        taxonomyDir = par.taxonomyPath + "/";
    }
    cout << "Taxonomy path: " << par.taxonomyPath << endl;
    fnaListFileName = par.filenames[1];
    acc2taxidFileName = par.filenames[2];

    // Output files
    taxidListFileName = dbDir + "/taxID_list";
    taxonomyBinaryFileName = dbDir + "/taxonomyDB";
    versionFileName = dbDir + "/db.version";
    paramterFileName = dbDir + "/db.parameters";

    taxonomy = new NcbiTaxonomy(taxonomyDir + "/names.dmp",
                                taxonomyDir + "/nodes.dmp",
                                taxonomyDir + "/merged.dmp");    

    if (par.reducedAA == 1){
        MARKER = 0Xffffffff;
        MARKER = ~ MARKER;
    } else {
        MARKER = 16777215;
        MARKER = ~ MARKER;
    }

    // For masking low complexity regions
    subMat = new NucleotideMatrix(par.scoringMatrixFile.values.nucleotide().c_str(), 1.0, 0.0);
}


IndexCreator::~IndexCreator() {
    delete taxonomy;
    delete subMat;
}

void IndexCreator::createIndex(const LocalParameters &par) {
    // Index library reference sequences
    indexReferenceSequences();
    
    if (!par.cdsInfo.empty()) {
        loadCdsInfo(par.cdsInfo);   
    }

    cout << "Made blocks for each thread" << endl;

    // Write taxonomy id list
    FILE * taxidListFile = fopen(taxidListFileName.c_str(), "w");
    for (auto & taxid: taxIdSet) {
        fprintf(taxidListFile, "%d\n", taxid);
    }
    fclose(taxidListFile);

    // Process the splits until all are processed
    size_t batchNum = accessionBatches.size();
    bool * batchChecker = new bool[batchNum];
    fill_n(batchChecker, batchNum, false);
    size_t processedBatchCnt = 0;
    TargetKmerBuffer kmerBuffer(par.bufferSize);
    vector<pair<size_t, size_t>> uniqKmerIdxRanges;
    cout << "Kmer buffer size: " << kmerBuffer.bufferSize << endl;
#ifdef OPENMP
    omp_set_num_threads(par.threads);
#endif
    while(processedBatchCnt < batchNum) { // Check this condition
        // Initialize the k-mer buffer
        memset(kmerBuffer.buffer, 0, kmerBuffer.bufferSize * sizeof(TargetKmer));

        // Extract Target k-mers
        fillTargetKmerBuffer(kmerBuffer, batchChecker, processedBatchCnt, par);

        // Sort the k-mers
        time_t start = time(nullptr);
        SORT_PARALLEL(kmerBuffer.buffer, kmerBuffer.buffer + kmerBuffer.startIndexOfReserve,
                      IndexCreator::compareForDiffIdx);
        time_t sort = time(nullptr);
        cout << "Sort time: " << sort - start << endl;
        auto * uniqKmerIdx = new size_t[kmerBuffer.startIndexOfReserve + 1];
        size_t uniqKmerCnt = 0;
        uniqKmerIdxRanges.clear();

        // Reduce redundancy
        reduceRedundancy(kmerBuffer, uniqKmerIdx, uniqKmerCnt, uniqKmerIdxRanges, par);
        time_t reduction = time(nullptr);
        cout << "Time spent for reducing redundancy: " << (double) (reduction - sort) << endl;
        if(processedBatchCnt == batchNum && numOfFlush == 0){
            writeTargetFilesAndSplits(kmerBuffer.buffer, kmerBuffer.startIndexOfReserve, par, uniqKmerIdx, uniqKmerCnt, uniqKmerIdxRanges);
        } else {
            writeTargetFiles(kmerBuffer.buffer, kmerBuffer.startIndexOfReserve, par, uniqKmerIdx, uniqKmerIdxRanges);
        }
        delete[] uniqKmerIdx;
    }
    delete[] batchChecker;
    writeDbParameters();
    writeTaxonomyDB();
}

void IndexCreator::updateIndex(const LocalParameters &par) {

}


// void IndexCreator::makeBlocksForParallelProcessing_accession_level() {
//     unordered_map<string, TaxID> acc2taxid;
//     TaxID maxTaxID = load_accession2taxid(fnaListFileName, acc2taxidFileName, acc2taxid);
//     newTaxID = std::max(getMaxTaxID() + 1, maxTaxID + 1);
    
//     vector<pair<string,pair<TaxID, TaxID>>> newAcc2taxid; // accession.version -> (parent, newTaxID)

//     // Make blocks of sequences that can be processed in parallel
//     int fileNum = getNumberOfLines(fnaListFileName);
//     fastaList.resize(fileNum);

//     ifstream fnaListFile;
//     fnaListFile.open(fnaListFileName);
//     if (!fnaListFile.is_open()) {
//         Debug(Debug::ERROR) << "Cannot open file for file list" << "\n";
//         EXIT(EXIT_FAILURE);
//     }
//     string eachFile;
//     string seqHeader;
//     string accession_version;
//     string accession;
//     vector<TaxID> tempTaxIDList;

//     for (int i = 0; i < fileNum; ++i) {
//         // Get start and end position of each sequence in the file
//         getline(fnaListFile, eachFile);
//         fastaList[i].path = eachFile;
//         processedSeqCnt.push_back(taxIdList.size());
//         seqHeader = getSeqSegmentsWithHead(fastaList[i].sequences, eachFile, acc2taxid, newAcc2taxid);
//         accession_version = seqHeader.substr(1, LocalUtil::getFirstWhiteSpacePos(seqHeader) - 1);
//         accession = accession_version.substr(0, accession_version.find('.'));
//         tempTaxIDList.push_back(acc2taxid.at(accession));   
//     }

//     // Edit taxonomy dump files
//     editTaxonomyDumpFiles(newAcc2taxid);

//     // Load taxonomy
//     taxonomy = new NcbiTaxonomy(taxonomyDir + "/names.dmp.new",
//                                 taxonomyDir + "/nodes.dmp.new",
//                                 taxonomyDir + "/merged.dmp");


//     for (int i = 0; i < fileNum; ++i) {
//         TaxID speciesTaxid = taxonomy->getTaxIdAtRank(tempTaxIDList[i], "species");
//         splitFastaForProdigalTraining(i, speciesTaxid);
//         fastaList[i].speciesID = speciesTaxid;
//     }
//     fnaListFile.close();

//     // Write accession to taxid map to file
//     string acc2taxidFileName2 = dbDir + "/acc2taxid.map";
//     FILE * acc2taxidFile = fopen(acc2taxidFileName2.c_str(), "w");
//     for (auto it : newAcc2taxid) {
//         fprintf(acc2taxidFile, "%s\t%d\t%d\n", it.first.c_str(), it.second.first, it.second.second);
//     }
//     fclose(acc2taxidFile);
// }

void IndexCreator::indexReferenceSequences() {
    vector<Accession> observedAccessionsVec;
    unordered_map<string, size_t> accession2index;
    getObservedAccessions(par.filenames[1], observedAccessionsVec, accession2index);
    cout << "Number of observed accessions: " << observedAccessionsVec.size() << endl;
    getTaxonomyOfAccessions(observedAccessionsVec, accession2index, par.filenames[2]);
    // cout << "Number of observed taxids: " << taxIdSet.size() << endl;
    getAccessionBatches(observedAccessionsVec);
    cout << "Number of accession batches: " << accessionBatches.size() << endl;
    // for (size_t i = 0; i < accessionBatches.size(); ++i) {
    //     cout << "Batch " << i << endl;
    //     accessionBatches[i].print();
    // }
}

void IndexCreator::getObservedAccessions(const string & fnaListFileName,
                                         vector<Accession> & observedAccessionsVec,
                                         unordered_map<string, size_t> & accession2index) {
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
        
        #pragma omp for schedule(dynamic, 1)
        for (size_t i = 0; i < fastaPaths.size(); ++i) {
            KSeqWrapper* kseq = KSeqFactory(fastaPaths[i].c_str());
            uint32_t order = 0;
            while (kseq->ReadEntry()) {
                const KSeqWrapper::KSeqEntry & e = kseq->entry;
                // Get the accession ID without version
                char* pos = strchr(e.name.s, '.'); 
                if (pos != nullptr) {
                    *pos = '\0';
                    localObservedAccessionsVec.emplace_back(string(e.name.s), i, order, e.sequence.l + e.name.l + e.comment.l);
                } else {
                    localObservedAccessionsVec.emplace_back(string(e.name.s), i, order, e.sequence.l + e.name.l + e.comment.l);
                }
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
            // __sync_fetch_and_add(&copyCount, localObservedAccessionsVec.size());            
        }                     
    }
}

void IndexCreator::getTaxonomyOfAccessions(vector<Accession> & observedAccessionsVec,
                                           const unordered_map<string, size_t> & accession2index,
                                           const string & acc2taxidFileName) {
    vector<pair<string, pair<TaxID, TaxID>>> acc2accId;                                       

    // Open the file
    int fd = open(acc2taxidFileName.c_str(), O_RDONLY);
    if (fd < 0) {
        cerr << "Cannot open file for mapping from accession to tax ID" << endl;
        return;
    }

    // Get the size of the file
    struct stat sb;
    if (fstat(fd, &sb) == -1) {
        cerr << "Cannot get the size of the file for mapping from accession to tax ID" << endl;
        close(fd);
        return;
    }

    size_t fileSize = sb.st_size;

    // Map the file to memory
    char* fileData = static_cast<char*>(mmap(nullptr, fileSize, PROT_READ, MAP_PRIVATE, fd, 0));
    if (fileData == MAP_FAILED) {
        cerr << "mmap failed" << endl;
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
    TaxID taxID, newTaxID = 0, maxTaxID = 0;
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
                maxTaxID = std::max(maxTaxID, taxID);
                // TODO: Accession level taxonomy
                if (par.accessionLevel == 1) {
                    acc2accId.emplace_back(accession, make_pair(taxID, newTaxID));
                    taxIdSet.insert(newTaxID);
                    observedAccessionsVec[it->second].taxID = newTaxID;
                    observedAccessionsVec[it->second].speciesID = taxonomy->getTaxIdAtRank(taxID, "species");
                    newTaxID++;
                } else {
                    taxIdSet.insert(taxID);
                    observedAccessionsVec[it->second].taxID = taxID;
                    observedAccessionsVec[it->second].speciesID = taxonomy->getTaxIdAtRank(taxID, "species");
                }
            }
        }
        ++current;  // Move to the next line
    }
    // Unmap the file
    if (munmap(fileData, fileSize) == -1) {
        cerr << "munmap failed" << endl;
    }                 

    if (par.accessionLevel == 1) {
        maxTaxID = std::max(getMaxTaxID() + 1, maxTaxID + 1);
        for (size_t i = 0; i < acc2accId.size(); ++i) {
            acc2accId[i].second.second += maxTaxID;
        }
        for (size_t i = 0; i < observedAccessionsVec.size(); ++i) {
            observedAccessionsVec[i].taxID += maxTaxID;
        }
        std::unordered_set<TaxID> updatedTaxIdSet;
        updatedTaxIdSet.reserve(taxIdSet.size());
        for (const auto& taxId : taxIdSet) {
            updatedTaxIdSet.insert(taxId + maxTaxID);
        }
        taxIdSet = std::move(updatedTaxIdSet);
        

        editTaxonomyDumpFiles(acc2accId);
        delete taxonomy;
        cout << "Accession level database is being created (--accession-level 1)" << endl;
        cout << "New taxonomic file containing accessions as children of taxa were created as follows:" << endl;
        cout << taxonomyDir + "/names.dmp.new" << endl;
        cout << taxonomyDir + "/nodes.dmp.new" << endl;
        cout << "Loading new taxonomy files" << endl;
        taxonomy = new NcbiTaxonomy(taxonomyDir + "/names.dmp.new",
                                    taxonomyDir + "/nodes.dmp.new",
                                    taxonomyDir + "/merged.dmp");
    }                       
}

void IndexCreator::getAccessionBatches(std::vector<Accession> & observedAccessionsVec) {
    size_t accCnt = observedAccessionsVec.size();
    SORT_PARALLEL(observedAccessionsVec.begin(), observedAccessionsVec.end(), Accession::compare);
    vector<uint32_t> orders;
    vector<TaxID> taxIDs;
    vector<uint32_t> lengths;
    for (size_t i = 0; i < accCnt;) {
        TaxID currentSpeciesID = observedAccessionsVec[i].speciesID;
        if (currentSpeciesID == 0) {
            continue;
        }
        uint32_t maxLength = 0, trainingFasta = 0, trainingSeq = 0;
        size_t firstBatchOfSpecies = accessionBatches.size();
        while (i < accCnt && currentSpeciesID == observedAccessionsVec[i].speciesID) {    
            uint32_t lengthSum = 0;
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
                orders.push_back(observedAccessionsVec[i].order);
                lengths.push_back(observedAccessionsVec[i].length);
                taxIDs.push_back(observedAccessionsVec[i].taxID);
                i++;
                if (lengthSum > 100'000'000 || orders.size() > 300 || (orders.size() > 100 && lengthSum > 50'000'000)) {
                    break;
                }
            }
            // Add the batch
            accessionBatches.emplace_back(currentFasta, currentSpeciesID, 0, 0);
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

// This function sort the TargetKmerBuffer, do redundancy reducing task, write the differential index of them
void IndexCreator::writeTargetFiles(TargetKmer * kmerBuffer, 
                                    size_t & kmerNum,
                                    const LocalParameters & par,
                                    const size_t * uniqKmerIdx,
                                    const vector<pair<size_t, size_t>> & uniqKmerIdxRanges) {
    string diffIdxFileName;
    string infoFileName;
    diffIdxFileName = dbDir + "/" + to_string(numOfFlush) + "_diffIdx";
    infoFileName = dbDir + "/" + to_string(numOfFlush) + "_info";

    FILE * diffIdxFile = fopen(diffIdxFileName.c_str(), "wb");
    FILE * infoFile = fopen(infoFileName.c_str(), "wb");
    if (diffIdxFile == nullptr || infoFile == nullptr){
        cout<<"Cannot open the file for writing target DB"<<endl;
        return;
    }
    numOfFlush++;
    size_t bufferSize = 1024 * 1024 * 32;
    uint16_t *diffIdxBuffer = (uint16_t *)malloc(sizeof(uint16_t) * bufferSize); // 64MB
    TaxID *infoBuffer = (TaxID *)malloc(sizeof(TaxID) * bufferSize);
    size_t localBufIdx = 0;
    size_t localInfoBufIdx = 0;
    uint64_t lastKmer = 0;
    size_t write = 0;

    for (size_t i = 0; i < uniqKmerIdxRanges.size(); i ++) {
        for (size_t j = uniqKmerIdxRanges[i].first; j < uniqKmerIdxRanges[i].second; j ++) {
            writeInfo(&kmerBuffer[uniqKmerIdx[j]].seqId, infoFile, infoBuffer, bufferSize, localInfoBufIdx);
            // fwrite(& kmerBuffer[uniqKmerIdx[j]].info, sizeof (TargetKmerInfo), 1, infoFile);
            write++;
            getDiffIdx(lastKmer, kmerBuffer[uniqKmerIdx[j]].ADkmer, diffIdxFile, diffIdxBuffer, bufferSize, localBufIdx);
            lastKmer = kmerBuffer[uniqKmerIdx[j]].ADkmer;
        }
    }
    flushKmerBuf(diffIdxBuffer, diffIdxFile, localBufIdx);
    flushInfoBuf(infoBuffer, infoFile, localInfoBufIdx);

    cout<<"total k-mer count  : "<< kmerNum << endl;
    cout<<"written k-mer count: "<< write << endl;
    free(infoBuffer);
    free(diffIdxBuffer);
    fclose(diffIdxFile);
    fclose(infoFile);
    kmerNum = 0;
}

void IndexCreator::writeTargetFilesAndSplits(TargetKmer * kmerBuffer,
                                             size_t & kmerNum,
                                             const LocalParameters & par,
                                             const size_t * uniqKmerIdx,
                                             size_t & uniqKmerCnt,
                                             const vector<pair<size_t, size_t>> & uniqKmerIdxRanges){
    string diffIdxFileName;
    string infoFileName;
    string splitFileName;

    diffIdxFileName = dbDir + "/diffIdx";
    infoFileName = dbDir + "/info";
    splitFileName = dbDir + "/split";

    // Make splits
    FILE * diffIdxSplitFile = fopen(splitFileName.c_str(), "wb");
    DiffIdxSplit splitList[par.splitNum];
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
                if (AminoAcidPart(kmerBuffer[uniqKmerIdx[k]].ADkmer) 
                    != AminoAcidPart(kmerBuffer[uniqKmerIdx[k + 1]].ADkmer)) {
                    splitList[splitCnt].ADkmer = kmerBuffer[uniqKmerIdx[k + 1]].ADkmer;
                    splitCnt++;
                    found = true;
                    break;
                }
            }
            if (!found) {
                splitList[splitCnt].ADkmer = kmerBuffer[uniqKmerIdx[uniqKmerIdxRanges[j+1].first]].ADkmer;
                cout << "Split " << splitCnt << " at " << splitList[splitCnt].ADkmer << endl;
                splitCnt++;
            }
            break;
        }
    }

    FILE * diffIdxFile = fopen(diffIdxFileName.c_str(), "wb");
    FILE * infoFile = fopen(infoFileName.c_str(), "wb");
    if (diffIdxFile == nullptr || infoFile == nullptr){
        cout<<"Cannot open the file for writing target DB"<<endl;
        return;
    }
    numOfFlush++;

    size_t bufferSize = 1024 * 1024 * 32;
    uint16_t *diffIdxBuffer = (uint16_t *)malloc(sizeof(uint16_t) * bufferSize); // 64MB
    TaxID *infoBuffer = (TaxID *)malloc(sizeof(TaxID) * bufferSize);
    size_t localBufIdx = 0;
    size_t localInfoBufIdx = 0;
    uint64_t lastKmer = 0;
    size_t write = 0;

    size_t splitIdx = 1;
    size_t totalDiffIdx = 0;
    cout << "Writing k-mers to disk" << endl;
    for (size_t i = 0; i < uniqKmerIdxRanges.size(); i ++) {
        for (size_t j = uniqKmerIdxRanges[i].first; j < uniqKmerIdxRanges[i].second; j ++) {
            writeInfo(&kmerBuffer[uniqKmerIdx[j]].seqId, infoFile, infoBuffer, bufferSize, localInfoBufIdx);
            // fwrite(& kmerBuffer[uniqKmerIdx[j]].info, sizeof (TargetKmerInfo), 1, infoFile);
            write++;
            getDiffIdx(lastKmer, kmerBuffer[uniqKmerIdx[j]].ADkmer, diffIdxFile, diffIdxBuffer, bufferSize, localBufIdx, totalDiffIdx);
            lastKmer = kmerBuffer[uniqKmerIdx[j]].ADkmer;
            if((splitIdx < splitCnt) && (lastKmer == splitList[splitIdx].ADkmer)){
                splitList[splitIdx].diffIdxOffset = totalDiffIdx;
                splitList[splitIdx].infoIdxOffset = write;
                splitIdx ++;
            }
        }
    }
    cout<<"total k-mer count  : "<< kmerNum << endl;
    cout<<"written k-mer count: "<< write << endl;

    flushKmerBuf(diffIdxBuffer, diffIdxFile, localBufIdx);
    flushInfoBuf(infoBuffer, infoFile, localInfoBufIdx);
    fwrite(splitList, sizeof(DiffIdxSplit), par.splitNum, diffIdxSplitFile);

    free(diffIdxBuffer);
    free(infoBuffer);
    fclose(diffIdxSplitFile);
    fclose(diffIdxFile);
    fclose(infoFile);
    kmerNum = 0;
    cout << "Finished writing k-mers to disk" << endl;
}

void IndexCreator::reduceRedundancy(TargetKmerBuffer & kmerBuffer,
                                    size_t * uniqKmerIdx,
                                    size_t & uniqueKmerCnt,
                                    vector<pair<size_t, size_t>> & uniqKmerIdxRanges,
                                    const LocalParameters & par) {
    // Find the first index of garbage k-mer (UINT64_MAX)
    for(size_t checkN = kmerBuffer.startIndexOfReserve - 1; checkN != 0; checkN--){
        if(kmerBuffer.buffer[checkN].ADkmer != UINT64_MAX){
            kmerBuffer.startIndexOfReserve = checkN + 1;
            break;
        }
    }

    // Find the first index of meaningful k-mer
    size_t startIdx = 0;
    for(size_t i = 0; i < kmerBuffer.startIndexOfReserve ; i++){
        if(kmerBuffer.buffer[i].taxIdAtRank != 0){
            startIdx = i;
            break;
        }
    }

    cout << "Starting Idx: " << startIdx << endl;
    cout << "Ending Idx: " << kmerBuffer.startIndexOfReserve << endl;

    // Make splits
    vector<Split> splits;
    size_t splitWidth = (kmerBuffer.startIndexOfReserve - startIdx) / par.threads;
    for (int i = 0; i < par.threads - 1; i++) {
        for (size_t j = startIdx + splitWidth; j + 1 < kmerBuffer.startIndexOfReserve; j++) {
            if (AminoAcidPart(kmerBuffer.buffer[j].ADkmer) != AminoAcidPart(kmerBuffer.buffer[j + 1].ADkmer)) {
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
#pragma omp for schedule(dynamic, 1)
        for(size_t split = 0; split < splits.size(); split ++) {
            lookingKmer = & kmerBuffer.buffer[splits[split].offset];
            lookingIndex = splits[split].offset;
            endFlag = 0;
            for(size_t i = 1 + splits[split].offset; i < splits[split].end + 1 ; i++) {
                taxIds.clear();
                taxIds.push_back(lookingKmer->seqId);
                // Scan redundancy
                while(lookingKmer->taxIdAtRank == kmerBuffer.buffer[i].taxIdAtRank &&
                      lookingKmer->ADkmer == kmerBuffer.buffer[i].ADkmer){
                    taxIds.push_back(kmerBuffer.buffer[i].seqId);
                    i++;
                    if(i == splits[split].end + 1){
                        endFlag = 1;
                        break;
                    }
                }
                if(taxIds.size() > 1){
                    lookingKmer->seqId = taxonomy->LCA(taxIds)->taxId;
                } else {
                    lookingKmer->seqId= taxIds[0];
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
            if(!((kmerBuffer.buffer[splits[split].end - 1].ADkmer == kmerBuffer.buffer[splits[split].end].ADkmer) &&
                 (kmerBuffer.buffer[splits[split].end - 1].taxIdAtRank == kmerBuffer.buffer[splits[split].end].taxIdAtRank))){
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

void IndexCreator::getDiffIdx(const uint64_t & lastKmer,
                              const uint64_t & entryToWrite, 
                              FILE* handleKmerTable, 
                              uint16_t *kmerBuf,
                              size_t bufferSize, 
                              size_t & localBufIdx){
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

void IndexCreator::getDiffIdx(const uint64_t & lastKmer, const uint64_t & entryToWrite, FILE* handleKmerTable,
                              uint16_t *kmerBuf, size_t bufferSize, size_t & localBufIdx, size_t & totalBufferIdx){
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

void IndexCreator::writeInfo(TaxID * entryToWrite, FILE * infoFile, TaxID * infoBuffer, size_t bufferSize, size_t & infoBufferIdx) {
    if (infoBufferIdx >= bufferSize) {
        flushInfoBuf(infoBuffer, infoFile, infoBufferIdx);
    }
    memcpy(infoBuffer + infoBufferIdx, entryToWrite, sizeof(TaxID));
    infoBufferIdx++;
}

void IndexCreator::flushInfoBuf(TaxID * buffer, FILE * infoFile, size_t & localBufIdx) {
    fwrite(buffer, sizeof(TaxID), localBufIdx, infoFile);
    localBufIdx = 0;
}

int IndexCreator::getNumOfFlush() {return numOfFlush;}

inline bool IndexCreator::compareForDiffIdx(const TargetKmer & a, const TargetKmer & b){
    if (a.ADkmer != b.ADkmer) {
        return a.ADkmer < b.ADkmer;
    }
    return a.taxIdAtRank < b.taxIdAtRank;
}

// string IndexCreator::getSeqSegmentsWithHead(vector<SequenceBlock> & seqSegments,
//                                             const string & seqFileName,
//                                             const unordered_map<string, TaxID> & acc2taxid,
//                                             vector<pair<string,pair<TaxID, TaxID>>> & newAcc2taxid) {
//     struct stat stat1{};
//     stat(seqFileName.c_str(), &stat1);
//     size_t numOfChar = stat1.st_size;
//     string firstLine;
//     ifstream seqFile;
//     seqFile.open(seqFileName);
//     string eachLine;
//     size_t start = 0;
//     size_t pos;
//     vector<SequenceBlock> seqSegmentsTmp;
//     string accession;
//     string accession_version;
//     int taxid;

//     if (seqFile.is_open()) {
//         getline(seqFile, firstLine, '\n');
//         accession_version = firstLine.substr(1, LocalUtil::getFirstWhiteSpacePos(firstLine) - 1);
//         accession = accession_version.substr(0, accession_version.find('.'));
//         taxid = acc2taxid.at(accession);
//         if (taxid == 0) {
//             cerr << "Cannot find accession: " << accession_version << endl;
//             cerr << "Please run 'add-to-library' first." << endl;
//             exit(1);
//         }
//         newAcc2taxid.emplace_back(accession_version, make_pair(taxid, newTaxID));
//         taxIdList.push_back(newTaxID++);

//         while (getline(seqFile, eachLine, '\n')) {
//             if (eachLine[0] == '>') {
//                 accession_version = eachLine.substr(1, LocalUtil::getFirstWhiteSpacePos(eachLine) - 1);
//                 accession = accession_version.substr(0, accession_version.find('.'));
//                 taxid = acc2taxid.at(accession);
//                 if (taxid == 0) {
//                     cerr << "Cannot find accession: " << accession_version << endl;
//                     cerr << "Please run 'add-to-library' first." << endl;
//                     exit(1);
//                 }
//                 newAcc2taxid.emplace_back(accession_version, make_pair(taxid, newTaxID));
//                 taxIdList.push_back(newTaxID++);
//                 pos = (size_t) seqFile.tellg();
//                 seqSegmentsTmp.emplace_back(start, pos - eachLine.length() - 3,pos - eachLine.length() - start - 2);
//                 start = pos - eachLine.length() - 1;
//             }
//         }
//         seqSegmentsTmp.emplace_back(start, numOfChar - 2, numOfChar - start - 1);
//     } else {
//         cerr << "Unable to open file: " << seqFileName << endl;
//     }
//     seqFile.close();
//     seqSegments = std::move(seqSegmentsTmp);
//     return firstLine;
// }

// string IndexCreator::getSeqSegmentsWithHead(vector<SequenceBlock> & seqSegments,
//                                             const string & seqFileName,
//                                             const unordered_map<string, TaxID> & acc2taxid) {
//     struct stat stat1{};
//     stat(seqFileName.c_str(), &stat1);
//     size_t numOfChar = stat1.st_size;
//     ifstream seqFile;
//     seqFile.open(seqFileName);
//     string firstLine, eachLine;
//     size_t start = 0;
//     size_t pos;
//     vector<SequenceBlock> seqSegmentsTmp;
//     vector<string> headers;
//     size_t seqCnt = taxIdList.size();
//     string accession_version, accession;
//     int taxid;

//     if (seqFile.is_open()) {
//         getline(seqFile, firstLine, '\n');
//         accession_version = firstLine.substr(1, LocalUtil::getFirstWhiteSpacePos(firstLine) - 1);
//         accession = accession_version.substr(0, accession_version.find('.'));
//         taxid = acc2taxid.at(accession);
//         if (taxid == 0) {
//             cerr << "Cannot find accession: " << accession_version << endl;
//             cerr << "Please run 'add-to-library' first." << endl;
//             exit(1);
//         }
//         taxIdList.push_back(taxid);

//         while (getline(seqFile, eachLine, '\n')) {
//             if (eachLine[0] == '>') {
//                 accession_version = eachLine.substr(1, LocalUtil::getFirstWhiteSpacePos(eachLine) - 1);
//                 accession = accession_version.substr(0, accession_version.find('.'));
//                 taxid = acc2taxid.at(accession);
//                 if (taxid == 0) {
//                     cerr << "Cannot find accession: " << accession_version << endl;
//                     cerr << "Please run 'add-to-library' first." << endl;
//                     exit(1);
//                 }
//                 taxIdList.push_back(taxid);
//                 pos = (size_t) seqFile.tellg();
//                 seqSegmentsTmp.emplace_back(start, pos - eachLine.length() - 3,pos - eachLine.length() - start - 2);
//                 start = pos - eachLine.length() - 1;
//             }
//         }
//         seqSegmentsTmp.emplace_back(start, numOfChar - 2, numOfChar - start - 1, seqCnt);
//     } else {
//         cerr << "Unable to open file: " << seqFileName << endl;
//     }
//     seqFile.close();
//     seqSegments = std::move(seqSegmentsTmp);
//     return firstLine;                            
// }

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
        cerr<<"Cannot open file for mappig from assemlby accession to tax ID"<<endl;
    }
    map.close();
}

size_t IndexCreator::fillTargetKmerBuffer(TargetKmerBuffer &kmerBuffer,
                                          bool *checker,
                                          size_t &processedBatchCnt,
                                          const LocalParameters &par) {
    int hasOverflow = 0;
#pragma omp parallel default(none), shared(kmerBuffer, checker, processedBatchCnt, hasOverflow, par, cout)
    {
        ProbabilityMatrix probMatrix(*subMat);
        // ProdigalWrapper prodigal;
        SeqIterator seqIterator(par);
        size_t posToWrite;
        size_t orfNum;
        vector<PredictedBlock> extendedORFs;
        priority_queue<uint64_t> standardList;
        priority_queue<uint64_t> currentList;
        size_t lengthOfTrainingSeq;
        char *reverseCompliment;
        vector<bool> strandness;
        vector<uint64_t> intergenicKmers;
        vector<int> aaSeq;
        vector<string> cds;
        vector<string> nonCds;
        bool trained = false;
#pragma omp for schedule(dynamic, 1)
        for (size_t batchIdx = 0; batchIdx < accessionBatches.size(); batchIdx ++) {
            if (!checker[batchIdx] && !hasOverflow) {
                checker[batchIdx] = true;
                intergenicKmers.clear();
                strandness.clear();
                standardList = priority_queue<uint64_t>();

                // Estimate the number of k-mers to be extracted from current split
                size_t totalLength = 0;
                for (size_t p = 0; p < accessionBatches[batchIdx].lengths.size(); p++) {
                    totalLength += accessionBatches[batchIdx].lengths[p];
                }
                size_t estimatedKmerCnt = (totalLength + totalLength / 5) / 3;

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
                            const KSeqWrapper::KSeqEntry & e = kseq->entry;
                            // #pragma omp critical
                            // {
                            //     cout << "Processing " << e.name.s << "\t" 
                            //          << e.sequence.l << "\t" 
                            //          << accessionBatches[batchIdx].speciesID << "\t" 
                            //          << accessionBatches[batchIdx].taxIDs[idx] << "\n";
                            // }

                            // Mask low complexity regions
                            char *maskedSeq = nullptr;
                            if (par.maskMode) {
                                maskedSeq = new char[e.sequence.l + 1]; // TODO: reuse the buffer
                                SeqIterator::maskLowComplexityRegions(e.sequence.s, maskedSeq, probMatrix, par.maskProb, subMat);
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
                                    seqIterator.translate(cds[cdsCnt], aaSeq);
                                    tempCheck = seqIterator.fillBufferWithKmerFromBlock(
                                                    cds[cdsCnt].c_str(),
                                                    kmerBuffer,
                                                    posToWrite,
                                                    accessionBatches[batchIdx].taxIDs[idx],
                                                    accessionBatches[batchIdx].speciesID,
                                                    aaSeq);
                                    if (tempCheck == -1) {
                                        cout << "ERROR: Buffer overflow " << e.name.s << e.sequence.l << endl;
                                    }
                                }

                                for (size_t nonCdsCnt = 0; nonCdsCnt < nonCds.size(); nonCdsCnt ++) {
                                    seqIterator.translate(nonCds[nonCdsCnt], aaSeq);
                                    tempCheck = seqIterator.fillBufferWithKmerFromBlock(
                                                    nonCds[nonCdsCnt].c_str(),
                                                    kmerBuffer,
                                                    posToWrite,
                                                    accessionBatches[batchIdx].taxIDs[idx],
                                                    accessionBatches[batchIdx].speciesID,
                                                    aaSeq);
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
                                    // #pragma omp critical
                                    // {
                                    //     cout << "Training " << accessionBatches[batchIdx].speciesID << "\t" << training_seq->entry.name.s << "\t" << training_seq->entry.sequence.l << endl;
                                    // }
                                    lengthOfTrainingSeq = training_seq->entry.sequence.l;
                                    prodigal->is_meta = 0;
                                    if (lengthOfTrainingSeq < 100'000) {
                                        prodigal->is_meta = 1;
                                        prodigal->trainMeta(training_seq->entry.sequence.s);
                                    } else {
                                        prodigal->trainASpecies(training_seq->entry.sequence.s);
                                    }

                                    // Generate intergenic 23-mer list. It is used to determine extension direction of intergenic sequences.
                                    prodigal->getPredictedGenes(training_seq->entry.sequence.s);
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
                                    prodigal->getPredictedGenes(e.sequence.s);
                                    prodigal->removeCompletelyOverlappingGenes();
                                    prodigal->getExtendedORFs(prodigal->finalGenes, prodigal->nodes, extendedORFs,
                                                                 prodigal->fng, e.sequence.l,
                                                            orfNum, intergenicKmers, e.sequence.s);

                                    // Get k-mers from extended ORFs
                                    for (size_t orfCnt = 0; orfCnt < orfNum; orfCnt++) {
                                        aaSeq.clear();
                                        seqIterator.translateBlock(maskedSeq, extendedORFs[orfCnt], aaSeq, e.sequence.l);
                                        tempCheck = seqIterator.fillBufferWithKmerFromBlock(
                                                extendedORFs[orfCnt],
                                                maskedSeq,
                                                kmerBuffer,
                                                posToWrite,
                                                accessionBatches[batchIdx].taxIDs[idx],
                                                accessionBatches[batchIdx].speciesID,
                                                aaSeq);
                                        if (tempCheck == -1) {
                                            cout << "ERROR: Buffer overflow " << e.name.s << e.sequence.l << endl;
                                        }
                                    }
                                } else { // Reverse complement
                                    reverseCompliment = seqIterator.reverseCompliment(e.sequence.s, e.sequence.l);

                                    // Get extended ORFs
                                    prodigal->getPredictedGenes(reverseCompliment);
                                    prodigal->removeCompletelyOverlappingGenes();
                                    prodigal->getExtendedORFs(prodigal->finalGenes, prodigal->nodes, extendedORFs,
                                                                     prodigal->fng, e.sequence.l,
                                                                orfNum, intergenicKmers, reverseCompliment);

                                    // Get reverse masked sequence
                                    if (par.maskMode) {
                                        delete[] maskedSeq;
                                        maskedSeq = new char[e.sequence.l + 1];
                                        SeqIterator::maskLowComplexityRegions(reverseCompliment, maskedSeq, probMatrix, par.maskProb, subMat);
                                        maskedSeq[e.sequence.l] = '\0';
                                    } else {
                                        maskedSeq = reverseCompliment;
                                    }

                                    for (size_t orfCnt = 0; orfCnt < orfNum; orfCnt++) {
                                        aaSeq.clear();
                                        seqIterator.translateBlock(maskedSeq, extendedORFs[orfCnt], aaSeq, e.sequence.l);
                                        tempCheck = seqIterator.fillBufferWithKmerFromBlock(
                                                extendedORFs[orfCnt],
                                                maskedSeq,
                                                kmerBuffer,
                                                posToWrite,
                                                accessionBatches[batchIdx].taxIDs[idx],
                                                accessionBatches[batchIdx].speciesID,
                                                aaSeq);
                                        if (tempCheck == -1) {
                                            cout << "ERROR: Buffer overflow " << e.name.s << e.sequence.l << endl;
                                        }
                                    }
                                    free(reverseCompliment);  
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
                        cout << processedBatchCnt << "batches processed out of " << accessionBatches.size() << endl;
                    }
                } else {
                    checker[batchIdx] = false;
                    __sync_fetch_and_add(&hasOverflow, 1);
                    __sync_fetch_and_sub(&kmerBuffer.startIndexOfReserve, estimatedKmerCnt);
                }
                delete prodigal;   
            }
        }
    }

    cout << "Before return: " << kmerBuffer.startIndexOfReserve << endl;
    return 0;
}

void IndexCreator::writeTaxonomyDB() {
    std::pair<char *, size_t> serialized = NcbiTaxonomy::serialize(*taxonomy);
    FILE *handle = fopen(taxonomyBinaryFileName.c_str(), "w");
    if (handle == NULL) {
        Debug(Debug::ERROR) << "Could not open " << taxonomyBinaryFileName << " for writing\n";
        EXIT(EXIT_FAILURE);
    }
    fwrite(serialized.first, serialized.second, sizeof(char), handle);
    fclose(handle);
    free(serialized.first);
}

void IndexCreator::writeDbParameters() {
    FILE *handle = fopen(paramterFileName.c_str(), "w");
    if (handle == NULL) {
        Debug(Debug::ERROR) << "Could not open " << paramterFileName << " for writing\n";
        EXIT(EXIT_FAILURE);
    }
    fprintf(handle, "DB_name\t%s\n", dbName.c_str());
    fprintf(handle, "Creation_date\t%s\n", dbDate.c_str());
    fprintf(handle, "Metabuli commit used to create the DB\t%s\n", version);
    fprintf(handle, "Reduced_alphabet\t%d\n", reducedAA);
    // fprintf(handle, "Spaced_kmer_mask\t%s\n", spaceMask.c_str());
    fprintf(handle, "Accession_level\t%d\n", accessionLevel);
    fprintf(handle, "Mask_mode\t%d\n", lowComplexityMasking);
    fprintf(handle, "Mask_prob\t%f\n", lowComplexityMaskingThreshold);
    fprintf(handle, "Skip_redundancy\t1\n");
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
    if (cdsInfoList.is_open()) {
        string cdsInfoFile;
        while (getline(cdsInfoList, cdsInfoFile)) { // Read each CDS info file
            ifstream cdsInfo(cdsInfoFile);
            if (cdsInfo.is_open()) {
                string line;
                while (getline(cdsInfo, line)) { // Read each line of the CDS info file
                    if (line[0] == '>') { // Check if the line starts with ">"
                        // Get the accession number between the '|' and '.'.
                        size_t start = line.find('|') + 1;
                        size_t end = line.find('.', start);
                        string accession = line.substr(start, end - start + 2);
                        // cout << "Accession: " << accession << endl;
                        int frame = 1;
                        while (true) {
                            start = line.find('[', end) + 1;
                            end = line.find(']', start);
                            if (start == string::npos) { break;}
                            size_t equalPos = line.find('=', start);
                            string feature = line.substr(start, equalPos - start);
                            string value = line.substr(equalPos + 1, end - equalPos - 1);
                            if (feature == "pseudo") {
                                break;
                            } else if (feature == "protein" && value == "hypothetical protein") {
                                break;
                            } else if (feature == "frame") {
                                frame = stoi(value);
                            } else if (feature == "protein_id") {
                                // cout << "Protein ID: " << value << "\t" << frame << endl;
                                cdsInfoMap[accession].emplace_back(CDSinfo(prtId++, frame));
                            } else if (feature == "location") {
                                // cout << "Location: " << value << endl;
                                // Check if the location is complement
                                size_t complementPos = value.find('c');
                                bool isComplement = (complementPos != string::npos);
                                if (isComplement) {
                                    cdsInfoMap[accession].back().isComplement = true;
                                    value = value.substr(complementPos + 11, value.size() - complementPos - 12);
                                } else {
                                    cdsInfoMap[accession].back().isComplement = false;
                                }

                                // Check if spliced
                                size_t joinPos = value.find('j');
                                if (joinPos != string::npos) {
                                    value = value.substr(joinPos + 5, value.size() - joinPos - 6);
                                }
                                
                                // Load the locations
                                size_t commaPos = value.find(',');
                                size_t dotPos;
                                string locationBegin, locationEnd;
                                while (commaPos != string::npos) {
                                    dotPos = value.find('.');
                                    if (dotPos > commaPos) {
                                        // AAA,BBB..CCC
                                        locationBegin = value.substr(0, commaPos);
                                        locationEnd = value.substr(0, commaPos);
                                    } else {
                                        // AAA..BBB,CCC..DDD
                                        locationBegin = value.substr(0, dotPos);
                                        locationEnd = value.substr(dotPos + 2, commaPos - dotPos - 2);
                                    }
                                    // Check < and > signs
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
                                cout << value << endl;
                                cout << locationBegin << endl;
                                cout << locationEnd << endl;
                                
                                // Check < and > signs
                                if (locationBegin[0] == '<') {
                                    locationBegin = locationBegin.substr(1, locationBegin.size() - 1);
                                }
                                if (locationEnd[0] == '>') {
                                    locationEnd = locationEnd.substr(1, locationEnd.size() - 1);
                                }
                                cdsInfoMap[accession].back().loc.emplace_back(stoi(locationBegin), stoi(locationEnd));

                                // Frame correction
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
                }
            } else {
                Debug(Debug::ERROR) << "Cannot open file " << cdsInfoFile << "\n";
                EXIT(EXIT_FAILURE);
            }
        }
    } else {
        Debug(Debug::ERROR) << "Cannot open file " << cdsInfoFileList << "\n";
        EXIT(EXIT_FAILURE);
    }
}