#include "KmerMatcher.h"
#include "BitManipulateMacros.h"
#include "IndexCreator.h"
#include "Kmer.h"
#include "Mmap.h"
#include <ostream>
#include <vector>

KmerMatcher::KmerMatcher(const LocalParameters & par,
                         NcbiTaxonomy * taxonomy) : par(par) {                        
    // Parameters
    threads = par.threads;
    dbDir = par.filenames[1 + (par.seqMode == 2)];
    hammingMargin = par.hammingMargin;
    
    MARKER = 16777215;
    MARKER = ~ MARKER;
    totalMatchCnt = 0;

    this->taxonomy = taxonomy;
    loadTaxIdList(par);
}


KmerMatcher::~KmerMatcher() {
}

void KmerMatcher::loadTaxIdList(const LocalParameters & par) {
    cout << "Loading the list for taxonomy IDs ... ";
    if (par.contamList != "") {
        vector<string> contams = Util::split(par.contamList, ",");
        for (auto &contam : contams) {
            FILE *taxIdFile;
            cout << dbDir + "/" + contam + "/taxID_list" << endl;
            if ((taxIdFile = fopen((dbDir + "/" + contam + "/taxID_list").c_str(), "r")) == NULL) {
                std::cout << "Cannot open the taxID list file." << std::endl;
                return;
            }
            char taxID[100];
            while (feof(taxIdFile) == 0) {
                fscanf(taxIdFile, "%s", taxID);
                TaxID taxId = atol(taxID);
                TaxonNode const *taxon = taxonomy->taxonNode(taxId);
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
        }
    } else {
        FILE *taxIdFile;
        if ((taxIdFile = fopen((dbDir + "/taxID_list").c_str(), "r")) == NULL) {
            std::cout << "Cannot open the taxID list file." << std::endl;
            return;
        }
        char taxID[100];
        while (feof(taxIdFile) == 0) {
            fscanf(taxIdFile, "%s", taxID);
            TaxID taxId = atol(taxID);
            TaxonNode const *taxon = taxonomy->taxonNode(taxId);
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
    }
    cout << "Done" << endl;
}


bool KmerMatcher::matchKmers(QueryKmerBuffer * queryKmerBuffer,
                             Buffer<Match> * matchBuffer,
                             const string & db){
    std::cout << "Comparing query and reference metamers..." << std::endl;
    
    // Set database files
    if (db.empty()) {
        targetDiffIdxFileName = dbDir + "/diffIdx";
        targetInfoFileName = dbDir + "/info";
        diffIdxSplitFileName = dbDir + "/split";
    } else { // for the case of multiple databases
        targetDiffIdxFileName = dbDir + "/" + db + "/diffIdx";
        targetInfoFileName = dbDir + "/" + db + "/info";
        diffIdxSplitFileName = dbDir + "/" + db + "/split";
    }
 
    MmapedData<DiffIdxSplit> diffIdxSplits = mmapData<DiffIdxSplit>(diffIdxSplitFileName.c_str(), 3);
    size_t numOfDiffIdx = FileUtil::getFileSize(targetDiffIdxFileName) / sizeof(uint16_t);
    size_t queryKmerNum = queryKmerBuffer->startIndexOfReserve;
    QueryKmer * queryKmerList = queryKmerBuffer->buffer;
    
    // Find the first index of garbage query k-mer (UINT64_MAX) and discard from there
    for (size_t checkN = queryKmerNum - 1; checkN > 0; checkN--) {
        if (queryKmerList[checkN].ADkmer != UINT64_MAX) {
            queryKmerNum = checkN + 1;
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
    std::vector<QueryKmerSplit> querySplits;
    uint64_t queryAA;
    size_t quotient = queryKmerNum / threads;
    size_t remainder = queryKmerNum % threads;
    size_t startIdx = 0;
    size_t endIdx = 0; // endIdx is inclusive
    for (size_t i = 0; i < threads; i++) {
        endIdx = startIdx + quotient - 1;
        if (remainder > 0) {
            endIdx++;
            remainder--;
        }
        bool needLastTargetBlock = true;
        queryAA = AMINO_ACID_PART(queryKmerList[startIdx].ADkmer);
        for (size_t j = 0; j < numOfDiffIdxSplits_use; j ++) {
            if (queryAA <= AMINO_ACID_PART(diffIdxSplits.data[j].ADkmer)) {
                j = j - (j != 0);
                querySplits.emplace_back(startIdx, endIdx, endIdx - startIdx + 1, diffIdxSplits.data[j]);
                needLastTargetBlock = false;
                break;
            }
        }
        if (needLastTargetBlock) {
            querySplits.emplace_back(startIdx, endIdx, endIdx - startIdx + 1, diffIdxSplits.data[numOfDiffIdxSplits_use - 2]);
        }
        startIdx = endIdx + 1;
    }

    if (querySplits.size() != threads) {
        threads = querySplits.size();
    }

    bool *splitCheckList = (bool *) malloc(sizeof(bool) * threads);
    std::fill_n(splitCheckList, threads, false);
    time_t beforeSearch = time(nullptr);
    size_t totalOverFlowCnt = 0;
#pragma omp parallel default(none), shared(splitCheckList, totalOverFlowCnt, \
querySplits, queryKmerList, matchBuffer, cout, targetDiffIdxFileName, numOfDiffIdx, targetInfoFileName)
{
    // FILE
    FILE * diffIdxFp = fopen(targetDiffIdxFileName.c_str(), "rb");
    FILE * kmerInfoFp = fopen(targetInfoFileName.c_str(), "rb");

    // Target K-mer buffer
    uint16_t * diffIdxBuffer = (uint16_t *) malloc(sizeof(uint16_t) * (BufferSize + 1)); // size = 32 Mb
    TaxID * kmerInfoBuffer = (TaxID *) malloc(sizeof(TaxID) * (BufferSize + 1)); // 64 Mb
    size_t kmerInfoBufferIdx = 0;
    size_t diffIdxBufferIdx = 0;
    
    // Query variables
    uint64_t currentQuery = UINT64_MAX;
    uint64_t currentQueryAA = UINT64_MAX;
    QueryKmerInfo currentQueryInfo;
        
    // Target variables
    size_t diffIdxPos = 0;
    std::vector<uint64_t> candidateTargetKmers; // vector for candidate target k-mer, some of which are selected after based on hamming distance
    std::vector<TaxID> candidateKmerInfos;
    std::vector<uint8_t> hammingDists;
    uint64_t currentTargetKmer;

    // Match buffer for each thread
    size_t localBufferSize = 2'000'000; // 32 Mb
    auto *matches = new Match[localBufferSize]; // 16 * 2'000'000 = 32 Mb
    size_t matchCnt = 0;

    // Vectors for selected target k-mers
    std::vector<uint8_t> selectedHammingSum;
    std::vector<size_t> selectedMatches;
    std::vector<uint16_t> selectedHammings;
    std::vector<uint32_t> selectedDnaEncodings;
    selectedHammingSum.resize(1024);
    selectedMatches.resize(1024);
    selectedHammings.resize(1024);
    selectedDnaEncodings.resize(1024);
    size_t selectedMatchCnt = 0;

    size_t posToWrite;
    size_t idx;
    bool hasOverflow = false;

#pragma omp for schedule(dynamic, 1)
    for (size_t i = 0; i < querySplits.size(); i++) {
        if (totalOverFlowCnt > 0 || splitCheckList[i]) {
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
                
        if (querySplits[i].diffIdxSplit.ADkmer == 0 && querySplits[i].diffIdxSplit.diffIdxOffset == 0 
            && querySplits[i].diffIdxSplit.infoIdxOffset == 0) {
            currentTargetKmer = getNextTargetKmer(currentTargetKmer, diffIdxBuffer,
                                                  diffIdxBufferIdx, diffIdxPos);
        }
        
        currentQuery = UINT64_MAX;
        currentQueryAA = UINT64_MAX;
        size_t lastMovedQueryIdx = 0;
        for (size_t j = querySplits[i].start; j < querySplits[i].end + 1; j++) {
            // Reuse the comparison data if queries are exactly identical
            if (currentQuery == queryKmerList[j].ADkmer
                && (currentQueryInfo.frame/3 == queryKmerList[j].info.frame/3)) {
                // If local buffer is full, copy them to the shared buffer.
                if (unlikely(matchCnt + selectedMatchCnt > localBufferSize)) {
                    // Check if the shared buffer is full.
                    posToWrite = matchBuffer->reserveMemory(matchCnt);
                    if (unlikely(posToWrite + matchCnt >= matchBuffer->bufferSize)) {
                        hasOverflow = true;
                        __sync_fetch_and_add(&totalOverFlowCnt, 1);
                        __sync_fetch_and_sub(& matchBuffer->startIndexOfReserve, matchCnt);
                        break;
                    } 
                    moveMatches(matchBuffer->buffer + posToWrite, matches, matchCnt);
                    lastMovedQueryIdx = j;
                }
                for (size_t k = 0; k < selectedMatchCnt; k++) {
                    idx = selectedMatches[k];
                    matches[matchCnt] = {queryKmerList[j].info,
                                         candidateKmerInfos[idx],
                                         taxId2speciesId[candidateKmerInfos[idx]],
                                         selectedDnaEncodings[k],
                                         selectedHammings[k],
                                         selectedHammingSum[k]};
                    matchCnt++;
                }
                continue;
            }
            selectedMatchCnt = 0;

            // Reuse the candidate target k-mers to compare in DNA level if queries are the same at amino acid level but not at DNA level
            if (currentQueryAA == AMINO_ACID_PART(queryKmerList[j].ADkmer)) {
                compareDna(queryKmerList[j].ADkmer, candidateTargetKmers, hammingDists, selectedMatches,
                           selectedHammingSum, selectedHammings, selectedDnaEncodings, selectedMatchCnt, queryKmerList[j].info.frame);
                // If local buffer is full, copy them to the shared buffer.
                if (unlikely(matchCnt + selectedMatchCnt > localBufferSize)) {
                    // Check if the shared buffer is full.
                    posToWrite = matchBuffer->reserveMemory(matchCnt);
                    if (unlikely(posToWrite + matchCnt >= matchBuffer->bufferSize)) {
                        hasOverflow = true;
                        __sync_fetch_and_add(&totalOverFlowCnt, 1);
                        __sync_fetch_and_sub(& matchBuffer->startIndexOfReserve, matchCnt);
                        break;
                    } 
                    moveMatches(matchBuffer->buffer + posToWrite, matches, matchCnt);
                    lastMovedQueryIdx = j;
                }
                for (size_t k = 0; k < selectedMatchCnt; k++) {
                    idx = selectedMatches[k];
                    matches[matchCnt] = {queryKmerList[j].info,
                                         candidateKmerInfos[idx],
                                         taxId2speciesId[candidateKmerInfos[idx]],
                                         selectedDnaEncodings[k],
                                         selectedHammings[k],
                                         selectedHammingSum[k]};
                    matchCnt++;
                }
                currentQuery = queryKmerList[j].ADkmer;
                currentQueryAA = AMINO_ACID_PART(currentQuery);
                currentQueryInfo = queryKmerList[j].info;
                continue;
            }
            candidateTargetKmers.clear();
            candidateKmerInfos.clear();

            // Get next query, and start to find
            currentQuery = queryKmerList[j].ADkmer;
            currentQueryAA = AMINO_ACID_PART(currentQuery);
            currentQueryInfo = queryKmerList[j].info;

            // Skip target k-mers that are not matched in amino acid level
            while (diffIdxPos != numOfDiffIdx
                   && (currentQueryAA > AMINO_ACID_PART(currentTargetKmer))) {  
                if (unlikely(BufferSize < diffIdxBufferIdx + 7)){
                    loadBuffer(diffIdxFp, diffIdxBuffer, diffIdxBufferIdx, BufferSize, ((int)(BufferSize - diffIdxBufferIdx)) * -1 );
                }
                currentTargetKmer = getNextTargetKmer(currentTargetKmer, diffIdxBuffer,
                                                      diffIdxBufferIdx, diffIdxPos);
                kmerInfoBufferIdx ++;
            }

            if (currentQueryAA != AMINO_ACID_PART(currentTargetKmer)) {
                continue;
            } 
                    
            // Load target k-mers that are matched in amino acid level
            while (diffIdxPos != numOfDiffIdx &&
                   currentQueryAA == AMINO_ACID_PART(currentTargetKmer)) {
                candidateTargetKmers.push_back(currentTargetKmer);
                candidateKmerInfos.push_back(getKmerInfo(BufferSize, kmerInfoFp, kmerInfoBuffer, kmerInfoBufferIdx));
                if (unlikely(BufferSize < diffIdxBufferIdx + 7)){
                    loadBuffer(diffIdxFp, diffIdxBuffer, diffIdxBufferIdx, BufferSize, ((int)(BufferSize - diffIdxBufferIdx)) * -1 );
                }
                currentTargetKmer = getNextTargetKmer(currentTargetKmer, diffIdxBuffer,
                                                      diffIdxBufferIdx, diffIdxPos);
                kmerInfoBufferIdx ++;
            }

            if (candidateTargetKmers.size() > selectedMatches.size()) {
                selectedMatches.resize(candidateTargetKmers.size());
                selectedHammingSum.resize(candidateTargetKmers.size());
                selectedHammings.resize(candidateTargetKmers.size());
                selectedDnaEncodings.resize(candidateTargetKmers.size());
            }

            // Compare the current query and the loaded target k-mers and select
            compareDna(currentQuery, candidateTargetKmers, hammingDists, selectedMatches, selectedHammingSum,
                       selectedHammings, selectedDnaEncodings, selectedMatchCnt, queryKmerList[j].info.frame);

            // If local buffer is full, copy them to the shared buffer.
            if (unlikely(matchCnt + selectedMatchCnt > localBufferSize)) {
                // Check if the shared buffer is full.
                posToWrite = matchBuffer->reserveMemory(matchCnt);
                if (unlikely(posToWrite + matchCnt >= matchBuffer->bufferSize)) {
                    hasOverflow = true;
                    __sync_fetch_and_add(&totalOverFlowCnt, 1);
                    __sync_fetch_and_sub(&matchBuffer->startIndexOfReserve, matchCnt);
                    break;
                } 
                moveMatches(matchBuffer->buffer + posToWrite, matches, matchCnt);
                lastMovedQueryIdx = j;
            }

            for (size_t k = 0; k < selectedMatchCnt; k++) {
                idx = selectedMatches[k];
                matches[matchCnt] = {queryKmerList[j].info,
                                     candidateKmerInfos[idx],
                                     taxId2speciesId[candidateKmerInfos[idx]],
                                     selectedDnaEncodings[k],
                                     selectedHammings[k],
                                     selectedHammingSum[k]};
                matchCnt++;
            }
        } // End of one split

        // Move matches in the local buffer to the shared buffer
        posToWrite = matchBuffer->reserveMemory(matchCnt);
        if (unlikely(posToWrite + matchCnt >= matchBuffer->bufferSize)) {
            hasOverflow = true;
            __sync_fetch_and_add(&totalOverFlowCnt, 1);
            __sync_fetch_and_sub(& matchBuffer->startIndexOfReserve, matchCnt);
        } else {
            moveMatches(matchBuffer->buffer + posToWrite, matches, matchCnt);
        }

        // Check whether current split is completed or not
        if (!hasOverflow) {
            splitCheckList[i] = true;
        }
    } // End of omp for (Iterating for splits)
    delete[] matches;
    fclose(diffIdxFp);
    fclose(kmerInfoFp);
    free(diffIdxBuffer);
    free(kmerInfoBuffer);
} // End of omp parallel
        
    if (totalOverFlowCnt > 0) {
        return false;
    }
    std::cout << "Time spent for the comparison: " << double(time(nullptr) - beforeSearch) << std::endl;
    free(splitCheckList);
    totalMatchCnt += matchBuffer->startIndexOfReserve;
    return true;
}

void KmerMatcher::sortMatches(Buffer<Match> * matchBuffer) {
    time_t beforeSortMatches = time(nullptr);
    std::cout << "Sorting matches ..." << std::endl;
    SORT_PARALLEL(matchBuffer->buffer,
                  matchBuffer->buffer + matchBuffer->startIndexOfReserve,
                  compareMatches);
    std::cout << "Time spent for sorting matches: " << double(time(nullptr) - beforeSortMatches) << std::endl;
}

void KmerMatcher::moveMatches(Match *dest, Match *src, size_t & matchNum) {
    memcpy(dest, src, sizeof(Match) * matchNum);
    matchNum = 0;
}

// It compares query k-mers to target k-mers.
// If a query has matches, the matches with the smallest hamming distance will be selected
void KmerMatcher::compareDna(uint64_t query,
                             std::vector<uint64_t> &targetKmersToCompare,
                             std::vector<uint8_t> & hammingDists,
                             std::vector<size_t> &selectedMatches,
                             std::vector<uint8_t> &selectedHammingSum,
                             std::vector<uint16_t> &selectedHammings,
                             std::vector<uint32_t> &selectedDnaEncodings,
                             size_t & selectedMatchIdx,
                             uint8_t frame) {
    hammingDists.resize(targetKmersToCompare.size());
    uint8_t minHammingSum = UINT8_MAX;

    // Calculate hamming distance
    for (size_t i = 0; i < targetKmersToCompare.size(); i++) {
        hammingDists[i] = getHammingDistanceSum(query, targetKmersToCompare[i]);
        minHammingSum = min(minHammingSum, hammingDists[i]);
    }

    // Select target k-mers that passed hamming criteria
    selectedMatchIdx = 0;
    uint8_t maxHamming = min(minHammingSum * 2, 7);
    for (size_t h = 0; h < targetKmersToCompare.size(); h++) {
        if (hammingDists[h] <= maxHamming) {
            selectedHammingSum[selectedMatchIdx] = hammingDists[h];
            selectedDnaEncodings[selectedMatchIdx] = GET_24_BITS_UINT(targetKmersToCompare[h]);
            selectedHammings[selectedMatchIdx] = (frame < 3)
                ? getHammings(query, targetKmersToCompare[h])
                : getHammings_reverse(query, targetKmersToCompare[h]);
            selectedMatches[selectedMatchIdx++] = h;
        }
    }
}



void KmerMatcher::compareDna2(uint64_t query,
                              const uint64_t * targetKmersToCompare,
                              size_t candidateCnt,
                              std::vector<uint8_t> & hammingDists,
                              std::vector<size_t> &selectedMatches,
                              std::vector<uint8_t> &selectedHammingSum,
                              std::vector<uint16_t> &selectedHammings,
                              std::vector<uint32_t> &selectedDnaEncodings,
                              size_t & selectedMatchIdx,
                              uint8_t frame) {
    hammingDists.resize(candidateCnt);
    uint8_t minHammingSum = UINT8_MAX;

    // Calculate hamming distance
    for (size_t i = 0; i < candidateCnt; i++) {
        hammingDists[i] = getHammingDistanceSum(query, targetKmersToCompare[i]);
        if (hammingDists[i] < minHammingSum) {
            minHammingSum = hammingDists[i];
        }
    }

    // Select target k-mers that passed hamming criteria
    selectedMatchIdx = 0;
    uint8_t maxHamming = min(minHammingSum * 2, 7);
    for (size_t h = 0; h < candidateCnt; h++) {
        if (hammingDists[h] <= maxHamming) {
            selectedMatches[selectedMatchIdx] = h;
            selectedHammingSum[selectedMatchIdx] = hammingDists[h];
            if (frame < 3) { // Frame of query k-mer
                selectedHammings[selectedMatchIdx] = getHammings(query, targetKmersToCompare[h]);
            } else {
                selectedHammings[selectedMatchIdx] = getHammings_reverse(query, targetKmersToCompare[h]);
            }
            // Store right 24 bits of the target k-mer in selectedDnaEncodings
            selectedDnaEncodings[selectedMatchIdx] = GET_24_BITS_UINT(targetKmersToCompare[h]);
            selectedMatchIdx++;
        }
    }
}


bool KmerMatcher::compareMatches(const Match& a, const Match& b) {
    if (a.qInfo.sequenceID != b.qInfo.sequenceID)
        return a.qInfo.sequenceID < b.qInfo.sequenceID;

    if (a.speciesId != b.speciesId)
        return a.speciesId < b.speciesId;

    if (a.qInfo.frame != b.qInfo.frame)
        return a.qInfo.frame < b.qInfo.frame;

    if (a.qInfo.pos != b.qInfo.pos)
        return a.qInfo.pos < b.qInfo.pos;

    if (a.hamming != b.hamming)
        return a.hamming < b.hamming;

    return a.dnaEncoding < b.dnaEncoding;
}

bool KmerMatcher::matchKmers_skipDecoding(QueryKmerBuffer * queryKmerBuffer,
                                          Buffer<Match> * matchBuffer,
                                          const string & db){
    // Set database files
    if (db.empty()) {
        targetDiffIdxFileName = dbDir + "/diffIdx";
        targetInfoFileName = dbDir + "/info";
        diffIdxSplitFileName = dbDir + "/split";
    } else { // for the case of multiple databases
        targetDiffIdxFileName = dbDir + "/" + db + "/diffIdx";
        targetInfoFileName = dbDir + "/" + db + "/info";
        diffIdxSplitFileName = dbDir + "/" + db + "/split";
    }
 
    MmapedData<DiffIdxSplit> diffIdxSplits = mmapData<DiffIdxSplit>(diffIdxSplitFileName.c_str(), 3);
    size_t numOfDiffIdx = FileUtil::getFileSize(targetDiffIdxFileName) / sizeof(uint16_t);

    // Load diffIdx count list
    string diffIdxAAFileName = targetDiffIdxFileName + ".aa";
    MmapedData<uint64_t> aminoacids = mmapData<uint64_t>(diffIdxAAFileName.c_str(), 3);
    size_t aaOffsetCnt = FileUtil::getFileSize(diffIdxAAFileName) / sizeof(uint64_t);
    // cout << "aaOffsetCnt: " << aaOffsetCnt << endl;

    // // Print target k-mer information
    // MmapedData<TargetKmerInfo> targetKmerInfo2 = mmapData<TargetKmerInfo>(targetInfoFileName.c_str(), 3);
    // size_t numOfTargetKmer = targetKmerInfo2.fileSize / sizeof(TargetKmerInfo);
    // for (size_t i = 0; i < numOfTargetKmer; i++) {
    //     cout << targetKmerInfo2.data[i].sequenceID << "\t" << (int) targetKmerInfo2.data[i].redundancy << endl;
    // }

    size_t queryKmerNum = queryKmerBuffer->startIndexOfReserve;
    QueryKmer *queryKmerList = queryKmerBuffer->buffer;
    
    std::cout << "Comparing query and reference metamers..." << std::endl;

    // Find the first index of garbage query k-mer (UINT64_MAX) and discard from there
    for (size_t checkN = queryKmerNum - 1; checkN > 0; checkN--) {
        if (queryKmerList[checkN].ADkmer != UINT64_MAX) {
            queryKmerNum = checkN + 1;
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
    std::vector<QueryKmerSplit> querySplits;
    uint64_t queryAA;
    size_t quotient = queryKmerNum / threads;
    size_t remainder = queryKmerNum % threads;
    size_t startIdx = 0;
    size_t endIdx = 0; // endIdx is inclusive
    for (size_t i = 0; i < threads; i++) {
        endIdx = startIdx + quotient - 1;
        if (remainder > 0) {
            endIdx++;
            remainder--;
        }
        bool needLastTargetBlock = true;
        queryAA = AMINO_ACID_PART(queryKmerList[startIdx].ADkmer);
        for (size_t j = 0; j < numOfDiffIdxSplits_use; j ++) {
            if (queryAA <= AMINO_ACID_PART(diffIdxSplits.data[j].ADkmer)) {
                j = j - (j != 0);
                querySplits.emplace_back(startIdx, endIdx, endIdx - startIdx + 1, diffIdxSplits.data[j]);
                needLastTargetBlock = false;
                break;
            }
        }
        if (needLastTargetBlock) {
            querySplits.emplace_back(startIdx, endIdx, endIdx - startIdx + 1, diffIdxSplits.data[numOfDiffIdxSplits_use - 2]);
        }
        startIdx = endIdx + 1;
    }
    munmap(diffIdxSplits.data, diffIdxSplits.fileSize + 1);

    if (querySplits.size() != threads) {
        threads = querySplits.size();
    }

    // Map querySplits to amino acid offset
    vector<uint64_t> offsets;
    size_t aaOffset = 0;
    bool isFound = false;
    for (size_t i = 0; i < querySplits.size(); i++) {
        for (; aaOffset < aaOffsetCnt; aaOffset++) {
            if (AMINO_ACID_PART(querySplits[i].diffIdxSplit.ADkmer) < aminoacids.data[aaOffset]) {
                offsets.push_back(aaOffset);
                isFound = true;
                break;
            }
        }
        if (!isFound) {
            offsets.push_back(offsets.back());
        }
        isFound = false;
        // cout << offsets.back() << endl;
    }
    munmap(aminoacids.data, aminoacids.fileSize + 1);


    bool *splitCheckList = (bool *) malloc(sizeof(bool) * threads);
    std::fill_n(splitCheckList, threads, false);
    size_t completedSplitCnt = 0;

    time_t beforeSearch = time(nullptr);

    // development
    size_t totalSkip = 0;

    while (completedSplitCnt < threads) {
        bool hasOverflow = false;
#pragma omp parallel default(none), shared(completedSplitCnt, splitCheckList, hasOverflow, \
querySplits, queryKmerList, matchBuffer, cout, targetDiffIdxFileName, numOfDiffIdx, targetInfoFileName, \
offsets, aaOffsetCnt, totalSkip)
        {
            // FILE
            FILE * diffIdxFp = fopen(targetDiffIdxFileName.c_str(), "rb");
            FILE * kmerInfoFp = fopen(targetInfoFileName.c_str(), "rb");
            FILE * aaFp = fopen((targetDiffIdxFileName + ".aa").c_str(), "rb");
            FILE * cntFp = fopen((targetDiffIdxFileName + ".deltaCnt").c_str(), "rb");
            FILE * kmerCntFp = fopen((targetDiffIdxFileName + ".kmerCnt").c_str(), "rb");
            FILE * kmerFp = fopen((targetDiffIdxFileName + ".kmers").c_str(), "rb");

            uint64_t * aaBuffer = (uint64_t *) malloc(sizeof(uint64_t) * (BufferSize + 1)); 
            uint64_t * nextKmers = (uint64_t *) malloc(sizeof(uint64_t) * (BufferSize + 1));
            uint32_t * cntBuffer = (uint32_t *) malloc(sizeof(uint32_t) * (BufferSize + 1));
            uint32_t * kmerCntBuffer = (uint32_t *) malloc(sizeof(uint32_t) * (BufferSize + 1));
            size_t aaOffsetIdx = 0;
            size_t aaOffsetIdx2 = 0;
            size_t totalOffsetIdx = 0;

            // Target K-mer buffer
            uint16_t * diffIdxBuffer = (uint16_t *) malloc(sizeof(uint16_t) * (BufferSize + 1)); // size = 32 Mb
            TaxID * kmerInfoBuffer = (TaxID *) malloc(sizeof(TaxID) * (BufferSize + 1)); // 64 Mb
            size_t kmerInfoBufferIdx = 0;
            size_t diffIdxBufferIdx = 0;

            // Query variables
            uint64_t currentQuery = UINT64_MAX;
            uint64_t currentQueryAA = UINT64_MAX;
            QueryKmerInfo currentQueryInfo;
            
            // Target variables
            size_t diffIdxPos = 0;
            std::vector<uint64_t> candidateTargetKmers; // vector for candidate target k-mer, some of which are selected after based on hamming distance
            std::vector<TaxID> candidateKmerInfos;
            std::vector<uint8_t> hammingDists;
            uint64_t currentTargetKmer;

            // Match buffer for each thread
            size_t localBufferSize = 2'000'000; // 32 Mb
            auto *matches = new Match[localBufferSize]; // 16 * 2'000'000 = 32 Mb
            size_t matchCnt = 0;

            // Vectors for selected target k-mers
            std::vector<uint8_t> selectedHammingSum;
            std::vector<size_t> selectedMatches;
            std::vector<uint16_t> selectedHammings;
            std::vector<uint32_t> selectedDnaEncodings;
            selectedHammingSum.resize(1024);
            selectedMatches.resize(1024);
            selectedHammings.resize(1024);
            selectedDnaEncodings.resize(1024);
            size_t selectedMatchCnt = 0;

            size_t posToWrite;
            size_t idx;

            size_t localSkip = 0;

#pragma omp for schedule(dynamic, 1)
            for (size_t i = 0; i < querySplits.size(); i++) {
                if (hasOverflow || splitCheckList[i]) {
                    continue;
                }
                aaOffsetIdx = offsets[i];
                aaOffsetIdx2 = offsets[i];
                totalOffsetIdx = aaOffsetIdx;
                currentTargetKmer = querySplits[i].diffIdxSplit.ADkmer;
                diffIdxBufferIdx = querySplits[i].diffIdxSplit.diffIdxOffset;
                kmerInfoBufferIdx = querySplits[i].diffIdxSplit.infoIdxOffset
                                    - (querySplits[i].diffIdxSplit.ADkmer != 0);
                diffIdxPos = querySplits[i].diffIdxSplit.diffIdxOffset;

                fseek(aaFp, 8 * (long) (aaOffsetIdx), SEEK_SET);
                loadBuffer(aaFp, aaBuffer, aaOffsetIdx, BufferSize);

                fseek(kmerFp, 8 * (long) (aaOffsetIdx2), SEEK_SET);
                loadBuffer(kmerFp, nextKmers, BufferSize);
                fseek(cntFp, 4 * (long) (aaOffsetIdx2), SEEK_SET);
                loadBuffer(cntFp, cntBuffer, BufferSize);
                fseek(kmerCntFp, 4 * (long) (aaOffsetIdx2), SEEK_SET);
                loadBuffer(kmerCntFp, kmerCntBuffer, aaOffsetIdx2, BufferSize);

                fseek(kmerInfoFp, 4 * (long)(kmerInfoBufferIdx), SEEK_SET);
                loadBuffer(kmerInfoFp, kmerInfoBuffer, kmerInfoBufferIdx, BufferSize);
                fseek(diffIdxFp, 2 * (long) (diffIdxBufferIdx), SEEK_SET);
                loadBuffer(diffIdxFp, diffIdxBuffer, diffIdxBufferIdx, BufferSize);
                
                if (querySplits[i].diffIdxSplit.ADkmer == 0 && querySplits[i].diffIdxSplit.diffIdxOffset == 0 
                    && querySplits[i].diffIdxSplit.infoIdxOffset == 0) {
                    currentTargetKmer = getNextTargetKmer(currentTargetKmer, diffIdxBuffer,
                                                          diffIdxBufferIdx, diffIdxPos);
                }
            
                currentQuery = UINT64_MAX;
                currentQueryAA = UINT64_MAX;

                size_t lastMovedQueryIdx = 0;
                for (size_t j = querySplits[i].start; j < querySplits[i].end + 1; j++) {
                    // Reuse the comparison data if queries are exactly identical
                    if (currentQuery == queryKmerList[j].ADkmer
                        && (currentQueryInfo.frame/3 == queryKmerList[j].info.frame/3)) {
                        // If local buffer is full, copy them to the shared buffer.
                        if (unlikely(matchCnt + selectedMatchCnt > localBufferSize)) {
                            // Check if the shared buffer is full.
                            posToWrite = matchBuffer->reserveMemory(matchCnt);
                            if (unlikely(posToWrite + matchCnt >= matchBuffer->bufferSize)) {
                                hasOverflow = true;
                                __sync_fetch_and_sub(& matchBuffer->startIndexOfReserve, matchCnt);
                                break;
                            } 
                            moveMatches(matchBuffer->buffer + posToWrite, matches, matchCnt);
                            lastMovedQueryIdx = j;
                        }
                        for (size_t k = 0; k < selectedMatchCnt; k++) {
                            idx = selectedMatches[k];
                            matches[matchCnt] = {queryKmerList[j].info,
                                                 candidateKmerInfos[idx],
                                                 taxId2speciesId[candidateKmerInfos[idx]],
                                                 selectedDnaEncodings[k],
                                                 selectedHammings[k],
                                                 selectedHammingSum[k]};
                            matchCnt++;
                        }
                        continue;
                    }
                    selectedMatchCnt = 0;

                    // Reuse the candidate target k-mers to compare in DNA level if queries are the same at amino acid level but not at DNA level
                    if (currentQueryAA == AMINO_ACID_PART(queryKmerList[j].ADkmer)) {
                        compareDna(queryKmerList[j].ADkmer, candidateTargetKmers, hammingDists, selectedMatches,
                                   selectedHammingSum, selectedHammings, selectedDnaEncodings, selectedMatchCnt, queryKmerList[j].info.frame);

                        // If local buffer is full, copy them to the shared buffer.
                        if (unlikely(matchCnt + selectedMatchCnt > localBufferSize)) {
                            // Check if the shared buffer is full.
                            posToWrite = matchBuffer->reserveMemory(matchCnt);
                            if (unlikely(posToWrite + matchCnt >= matchBuffer->bufferSize)) {
                                hasOverflow = true;
                                __sync_fetch_and_sub(& matchBuffer->startIndexOfReserve, matchCnt);
                                break;
                            } 
                            moveMatches(matchBuffer->buffer + posToWrite, matches, matchCnt);
                            lastMovedQueryIdx = j;
                        }
                        for (size_t k = 0; k < selectedMatchCnt; k++) {
                            idx = selectedMatches[k];
                            matches[matchCnt] = {queryKmerList[j].info,
                                                 candidateKmerInfos[idx],
                                                 taxId2speciesId[candidateKmerInfos[idx]],
                                                 selectedDnaEncodings[k],
                                                 selectedHammings[k],
                                                 selectedHammingSum[k]};
                            matchCnt++;
                        }
                        currentQuery = queryKmerList[j].ADkmer;
                        currentQueryAA = AMINO_ACID_PART(currentQuery);
                        currentQueryInfo = queryKmerList[j].info;
                        continue;
                    }
                    candidateTargetKmers.clear();
                    candidateKmerInfos.clear();

                    // Get next query, and start to find
                    currentQuery = queryKmerList[j].ADkmer;
                    currentQueryAA = AMINO_ACID_PART(currentQuery);
                    currentQueryInfo = queryKmerList[j].info;

                    // Skip target k-mers that are not matched in amino acid level
                    while (diffIdxPos != numOfDiffIdx
                           && (currentQueryAA > AMINO_ACID_PART(currentTargetKmer))) {
                        // seqIterator->printKmerInDNAsequence(currentTargetKmer); cout << " " << diffIdxPos << " " << diffIdxBufferIdx << " ";
                        // seqIterator->printAAKmer(AMINO_ACID_PART(aaBuffer[aaOffsetIdx]), 24); cout << "\n";
                        if (AMINO_ACID_PART(currentTargetKmer) == aaBuffer[aaOffsetIdx]) {
                            size_t temp = aaOffsetIdx2;
                            diffIdxBufferIdx += getElement(BufferSize, cntFp, cntBuffer, aaOffsetIdx2);
                            diffIdxPos += getElement(BufferSize, cntFp, cntBuffer, aaOffsetIdx2);
                            localSkip += getElement(BufferSize, cntFp, kmerCntBuffer, aaOffsetIdx2);
                            aaOffsetIdx2 = temp;
                            kmerInfoBufferIdx += getElement(BufferSize, kmerCntFp, kmerCntBuffer, aaOffsetIdx2);
                            aaOffsetIdx2 = temp;
                            currentTargetKmer = getElement(BufferSize, kmerFp, nextKmers, aaOffsetIdx2);
                            aaOffsetIdx++;
                            aaOffsetIdx2++;
                            totalOffsetIdx ++;
                            if (unlikely(aaOffsetIdx == BufferSize)) {
                                loadBuffer(aaFp, aaBuffer, aaOffsetIdx, BufferSize);
                            }
                            if (unlikely(BufferSize < diffIdxBufferIdx + 7)){
                                loadBuffer(diffIdxFp, diffIdxBuffer, diffIdxBufferIdx, BufferSize, ((int)(BufferSize - diffIdxBufferIdx)) * -1 );
                            }
                            continue;
                        } else {
                            while (totalOffsetIdx < aaOffsetCnt && AMINO_ACID_PART(currentTargetKmer) >= aaBuffer[aaOffsetIdx]) {
                                aaOffsetIdx++;
                                aaOffsetIdx2++;
                                totalOffsetIdx++;
                                if (unlikely(aaOffsetIdx == BufferSize)) {
                                    loadBuffer(aaFp, aaBuffer, aaOffsetIdx, BufferSize);
                                }
                            }                            
                        }  
                        if (unlikely(BufferSize < diffIdxBufferIdx + 7)){
                            loadBuffer(diffIdxFp, diffIdxBuffer, diffIdxBufferIdx, BufferSize, ((int)(BufferSize - diffIdxBufferIdx)) * -1 );
                        }
                        currentTargetKmer = getNextTargetKmer(currentTargetKmer, diffIdxBuffer,
                                                              diffIdxBufferIdx, diffIdxPos);
                        
                        kmerInfoBufferIdx ++;
                    }

                    if (currentQueryAA != AMINO_ACID_PART(currentTargetKmer)) {
                        continue;
                    } 
                    
                    // Load target k-mers that are matched in amino acid level
                    while (diffIdxPos != numOfDiffIdx &&
                           currentQueryAA == AMINO_ACID_PART(currentTargetKmer)) {
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

                    if (candidateTargetKmers.size() > selectedMatches.size()) {
                        selectedMatches.resize(candidateTargetKmers.size());
                        selectedHammingSum.resize(candidateTargetKmers.size());
                        selectedHammings.resize(candidateTargetKmers.size());
                        selectedDnaEncodings.resize(candidateTargetKmers.size());
                    }

                    // Compare the current query and the loaded target k-mers and select
                    compareDna(currentQuery, candidateTargetKmers, hammingDists, selectedMatches, selectedHammingSum,
                               selectedHammings, selectedDnaEncodings, selectedMatchCnt, queryKmerList[j].info.frame);

                    // If local buffer is full, copy them to the shared buffer.
                    if (unlikely(matchCnt + selectedMatchCnt > localBufferSize)) {
                        // Check if the shared buffer is full.
                        posToWrite = matchBuffer->reserveMemory(matchCnt);
                        if (unlikely(posToWrite + matchCnt >= matchBuffer->bufferSize)) { // full -> write matches to file first
                            hasOverflow = true;
                            __sync_fetch_and_sub(&matchBuffer->startIndexOfReserve, matchCnt);
                            break;
                        } 
                        moveMatches(matchBuffer->buffer + posToWrite, matches, matchCnt);
                        lastMovedQueryIdx = j;
                    }

                    for (size_t k = 0; k < selectedMatchCnt; k++) {
                        idx = selectedMatches[k];
                        matches[matchCnt] = {queryKmerList[j].info,
                                             candidateKmerInfos[idx],
                                             taxId2speciesId[candidateKmerInfos[idx]],
                                             selectedDnaEncodings[k],
                                             selectedHammings[k],
                                             selectedHammingSum[k]};
                        matchCnt++;
                    }
                } // End of one split

                // Move matches in the local buffer to the shared buffer
                posToWrite = matchBuffer->reserveMemory(matchCnt);
                if (unlikely(posToWrite + matchCnt >= matchBuffer->bufferSize)) {
                    hasOverflow = true;
                    __sync_fetch_and_sub(& matchBuffer->startIndexOfReserve, matchCnt);
                } else {
                    moveMatches(matchBuffer->buffer + posToWrite, matches, matchCnt);
                }

                // Check whether current split is completed or not
                if (!hasOverflow) {
                    splitCheckList[i] = true;
                    __sync_fetch_and_add(&completedSplitCnt, 1);
                }
            } // End of omp for (Iterating for splits)
            delete[] matches;
            fclose(diffIdxFp);
            fclose(kmerInfoFp);
            fclose(aaFp);
            fclose(cntFp);
            fclose(kmerCntFp);
            fclose(kmerFp);
            free(diffIdxBuffer);
            free(kmerInfoBuffer);
            free(aaBuffer);
            free(nextKmers);
            free(cntBuffer);
            free(kmerCntBuffer);
            __sync_fetch_and_add(&totalSkip, localSkip);
        } // End of omp parallel
        
        if (hasOverflow) {
            return false;
        }
    } // end of while(completeSplitCnt < threadNum)
    cout << "Total skipped diffIdx: " << totalSkip << endl;
    std::cout << "Time spent for the comparison: " << double(time(nullptr) - beforeSearch) << std::endl;
    free(splitCheckList);
    cout << "Match count: " << matchBuffer->startIndexOfReserve << endl;
    totalMatchCnt += matchBuffer->startIndexOfReserve;
    return true;
}