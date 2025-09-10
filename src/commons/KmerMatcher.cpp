#include "KmerMatcher.h"
#include "BitManipulateMacros.h"
#include "IndexCreator.h"
#include "Kmer.h"
#include "Mmap.h"
#include <ostream>
#include <vector>

constexpr uint16_t KmerMatcher::HAMMING_LUT0[64];
constexpr uint16_t KmerMatcher::HAMMING_LUT1[64];
constexpr uint16_t KmerMatcher::HAMMING_LUT2[64];
constexpr uint16_t KmerMatcher::HAMMING_LUT3[64];
constexpr uint16_t KmerMatcher::HAMMING_LUT4[64];
constexpr uint16_t KmerMatcher::HAMMING_LUT5[64];
constexpr uint16_t KmerMatcher::HAMMING_LUT6[64];
constexpr uint16_t KmerMatcher::HAMMING_LUT7[64];


KmerMatcher::KmerMatcher(
    const LocalParameters & par,
    TaxonomyWrapper * taxonomy,
    int kmerFormat) 
    : par(par), 
      kmerFormat(kmerFormat) 
{                        
    // Parameters
    threads = par.threads;
    dbDir = par.filenames[1 + (par.seqMode == 2)];
    hammingMargin = par.hammingMargin;
    DNA_MASK = 16777215;
    DNA_MASK = ~ DNA_MASK;
    AA_MASK = 0xffffffU;     
    totalMatchCnt = 0;
    this->taxonomy = taxonomy;
    geneticCode = new GeneticCode(par.reducedAA == 1);
    loadTaxIdList(par);
}

KmerMatcher::KmerMatcher(
    const LocalParameters &par,
    int kmerFormat) : par(par),
                      kmerFormat(kmerFormat)
{
    geneticCode = new GeneticCode(false);
    dbDir = par.filenames[1];
    targetDiffIdxFileName = dbDir + "/diffIdx";
    targetInfoFileName = dbDir + "/info";
    diffIdxSplitFileName = dbDir + "/split";
    totalMatchCnt = 0;
}

KmerMatcher::~KmerMatcher() {
    delete geneticCode;
}

void KmerMatcher::loadTaxIdList(const LocalParameters & par) {
    // cout << "Loading the list for taxonomy IDs ... " << std::flush;
    if (par.contamList != "") {
        vector<string> contams = Util::split(par.contamList, ",");
        for (auto &contam : contams) {
            string fileName = dbDir + "/" + contam + "/taxID_list";
            if (!FileUtil::fileExists(fileName.c_str())) {
                std::cout << "TaxID list file for " << contam << " does not exist." << std::endl;
                continue;
            }
            std::ifstream in{fileName};
            if (!in.is_open()) {
                std::cout << "Cannot open the taxID list file." << std::endl;
                return;
            }
            std::string line;
            while (std::getline(in, line)) {
                if (line.empty()) continue;                   
                TaxID taxId = static_cast<TaxID>(std::stoul(line));
                TaxID speciesTaxID = taxonomy->getTaxIdAtRank(taxId, "species");
                TaxID genusTaxID = taxonomy->getTaxIdAtRank(taxId, "genus");
                const TaxonNode* taxon = taxonomy->taxonNode(taxId);            
                if (taxId != taxon->taxId) {
                    taxId2speciesId[taxId] = speciesTaxID;
                    taxId2genusId[taxId] = genusTaxID;
                }
                while (taxon->taxId != speciesTaxID) {
                    taxId2speciesId[taxon->taxId] = speciesTaxID;
                    taxId2genusId[taxon->taxId] = genusTaxID;
                    taxon = taxonomy->taxonNode(taxon->parentTaxId);
                }
                taxId2speciesId[speciesTaxID] = speciesTaxID;
                taxId2genusId[speciesTaxID] = genusTaxID;
            }
            in.close();
        }
    } else {
        std::ifstream in{dbDir + "/taxID_list"};
        if (!in.is_open()) {
            std::cout << "Cannot open the taxID list file." << std::endl;
            return;
        }
        std::string line;
        while (std::getline(in, line)) {
            if (line.empty()) continue;                   
            TaxID taxId = static_cast<TaxID>(std::stoul(line));
            TaxID speciesTaxID = taxonomy->getTaxIdAtRank(taxId, "species");
            TaxID genusTaxID = taxonomy->getTaxIdAtRank(taxId, "genus");
            const TaxonNode* taxon = taxonomy->taxonNode(taxId);            
            if (taxId != taxon->taxId) {
                taxId2speciesId[taxId] = speciesTaxID;
                taxId2genusId[taxId] = genusTaxID;
            }
            while (taxon->taxId != speciesTaxID) {
                taxId2speciesId[taxon->taxId] = speciesTaxID;
                taxId2genusId[taxon->taxId] = genusTaxID;
                taxon = taxonomy->taxonNode(taxon->parentTaxId);
            }
            taxId2speciesId[speciesTaxID] = speciesTaxID;
            taxId2genusId[speciesTaxID] = genusTaxID;
        }
        in.close();
    }
    // cout << "Done" << endl;
}


bool KmerMatcher::matchKmers(Buffer<Kmer> * queryKmerBuffer,
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
    size_t queryKmerNum = queryKmerBuffer->startIndexOfReserve;
    Kmer * queryKmerList = queryKmerBuffer->buffer;
    

    size_t blankCnt = 0;
    for (size_t i = 0; i < queryKmerNum; i++) {
        if (queryKmerList[i].qInfo.sequenceID == 0) {
            blankCnt++;
        } else {
            break;
        }
    }
    queryKmerNum -= blankCnt;
    std::cout << "Query k-mer number     : " << queryKmerNum << endl;

    std::cout << "Reference k-mer search : " << flush;
    
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
    size_t startIdx = blankCnt;
    size_t endIdx = 0; // endIdx is inclusive
    for (size_t i = 0; i < threads; i++) {
        endIdx = startIdx + quotient - 1;
        if (remainder > 0) {
            endIdx++;
            remainder--;
        }
        bool needLastTargetBlock = true;
        queryAA = AMINO_ACID_PART(queryKmerList[startIdx].value);
        for (size_t j = 0; j < numOfDiffIdxSplits_use; j ++) {
            if (queryAA <= AMINO_ACID_PART(diffIdxSplits.data[j].ADkmer)) {
                j = j - (j != 0);
                querySplits.emplace_back(startIdx, endIdx, diffIdxSplits.data[j]);
                needLastTargetBlock = false;
                break;
            }
        }
        if (needLastTargetBlock) {
            querySplits.emplace_back(startIdx, endIdx, diffIdxSplits.data[numOfDiffIdxSplits_use - 2]);
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
    int redundancyStored = (par.skipRedundancy == 0);
    unsigned int mask = ~((static_cast<unsigned int>(redundancyStored) << 31));
    // 
#pragma omp parallel default(none), shared(splitCheckList, totalOverFlowCnt, \
querySplits, queryKmerList, matchBuffer, cout, mask, targetDiffIdxFileName, numOfDiffIdx, redundancyStored, targetInfoFileName)
{
    SeqIterator seqIter(par);
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
    selectedHammingSum.resize(1024);
    selectedMatches.resize(1024);
    selectedHammings.resize(1024);
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
        for (size_t j = querySplits[i].start; j < querySplits[i].end + 1; j++) {
            // Reuse the comparison data if queries are exactly identical
            if (currentQuery == queryKmerList[j].value
                && (currentQueryInfo.frame/3 == queryKmerList[j].qInfo.frame/3)) {
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
                }
                for (size_t k = 0; k < selectedMatchCnt; k++) {
                    if (candidateKmerInfos[idx] == 0 || taxId2speciesId[candidateKmerInfos[idx]] == 0) {
                        cout << candidateKmerInfos[idx] << "\t" << taxId2speciesId[candidateKmerInfos[idx]] << endl;
                        const TaxonNode * node = taxonomy->taxonNode(candidateKmerInfos[idx]);
                        cout << "TaxID: " << node->taxId << "\t"
                             << "Species: " << taxonomy->getOriginalTaxID(taxId2speciesId[candidateKmerInfos[idx]]) << endl;
                        cout << taxonomy->getString(node->nameIdx) << endl;
                        cout << taxonomy->getString(node->rankIdx) << endl;
                        exit(1);
                    }
                    idx = selectedMatches[k];
                    matches[matchCnt] = {queryKmerList[j].qInfo,
                                         candidateKmerInfos[idx],
                                         taxId2speciesId[candidateKmerInfos[idx]],
                                         (unsigned int) (candidateTargetKmers[idx] & AA_MASK),
                                         selectedHammings[k],
                                         selectedHammingSum[k]};
                    matchCnt++;
                }
                continue;
            }
            selectedMatchCnt = 0;

            // Reuse the candidate target k-mers to compare in DNA level if queries are the same at amino acid level but not at DNA level
            if (currentQueryAA == AMINO_ACID_PART(queryKmerList[j].value)) {
                compareDna(queryKmerList[j].value, candidateTargetKmers, hammingDists, selectedMatches,
                           selectedHammingSum, selectedHammings, selectedMatchCnt, queryKmerList[j].qInfo.frame);
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
                }
                for (size_t k = 0; k < selectedMatchCnt; k++) {
                    idx = selectedMatches[k];
                    if (candidateKmerInfos[idx] == 0 || taxId2speciesId[candidateKmerInfos[idx]] == 0) {
                        cout << candidateKmerInfos[idx] << "\t" << taxId2speciesId[candidateKmerInfos[idx]] << endl;
                        const TaxonNode * node = taxonomy->taxonNode(candidateKmerInfos[idx]);
                        cout << "TaxID: " << node->taxId << "\t"
                             << "Species: " << taxonomy->getOriginalTaxID(taxId2speciesId[candidateKmerInfos[idx]]) << endl;
                        cout << taxonomy->getString(node->nameIdx) << endl;
                        cout << taxonomy->getString(node->rankIdx) << endl;
                        exit(1);
                    }
                    matches[matchCnt] = {queryKmerList[j].qInfo,
                                         candidateKmerInfos[idx],
                                         taxId2speciesId[candidateKmerInfos[idx]],
                                         (unsigned int) (candidateTargetKmers[idx] & AA_MASK),
                                         selectedHammings[k],
                                         selectedHammingSum[k]};
                    matchCnt++;
                }
                currentQuery = queryKmerList[j].value;
                currentQueryAA = AMINO_ACID_PART(currentQuery);
                currentQueryInfo = queryKmerList[j].qInfo;
                continue;
            }
            candidateTargetKmers.clear();
            candidateKmerInfos.clear();

            // Get next query, and start to find
            currentQuery = queryKmerList[j].value;
            currentQueryAA = AMINO_ACID_PART(currentQuery);
            currentQueryInfo = queryKmerList[j].qInfo;

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
                candidateKmerInfos.push_back(getKmerInfo<TaxID>(BufferSize, kmerInfoFp, kmerInfoBuffer, kmerInfoBufferIdx) & mask);
                // if (par.printLog) {
                //     cout << "# "<< queryKmerList[j].info.sequenceID << "\t" << queryKmerList[j].info.pos << "\t"
                //          << (int) queryKmerList[j].info.frame << endl;
                //     cout << "Query  k-mer: ";
                //     print_binary64(64, currentQuery);
                //     cout << "\t";
                //     seqIter.printKmerInDNAsequence(currentQuery);
                //     cout << endl;
                //     cout << "Target k-mer: ";
                //     print_binary64(64, currentTargetKmer);
                //     cout << "\t";
                //     seqIter.printKmerInDNAsequence(currentTargetKmer);
                //     cout << endl;
                //     cout << "\t" << taxonomy->getOriginalTaxID(kmerInfoBuffer[kmerInfoBufferIdx])
                //          << "\t" << taxonomy->getOriginalTaxID(taxId2speciesId[kmerInfoBuffer[kmerInfoBufferIdx]])<< endl;
                //     cout << (int) getHammingDistanceSum(currentQuery, currentTargetKmer) << "\t";
                //     print_binary16(16, getHammings(currentQuery, currentTargetKmer)); cout << endl;
                // }
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
            }

            // Compare the current query and the loaded target k-mers and select
            compareDna(currentQuery, candidateTargetKmers, hammingDists, selectedMatches, selectedHammingSum,
                       selectedHammings, selectedMatchCnt, queryKmerList[j].qInfo.frame);

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
            }

            for (size_t k = 0; k < selectedMatchCnt; k++) {
                idx = selectedMatches[k];
                    if (candidateKmerInfos[idx] == 0 || taxId2speciesId[candidateKmerInfos[idx]] == 0) {
                        cout << candidateKmerInfos[idx] << "\t" << taxId2speciesId[candidateKmerInfos[idx]] << endl;
                        const TaxonNode * node = taxonomy->taxonNode(candidateKmerInfos[idx]);
                        cout << "TaxID: " << node->taxId << "\t"
                             << "Species: " << taxonomy->getOriginalTaxID(taxId2speciesId[candidateKmerInfos[idx]]) << endl;
                        cout << taxonomy->getString(node->nameIdx) << endl;
                        cout << taxonomy->getString(node->rankIdx) << endl;
                        exit(1);
                    }
                matches[matchCnt] = {queryKmerList[j].qInfo,
                                     candidateKmerInfos[idx],
                                     taxId2speciesId[candidateKmerInfos[idx]],
                                     (unsigned int) (candidateTargetKmers[idx] & AA_MASK),
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
    std::cout << double(time(nullptr) - beforeSearch) << " s" << std::endl;
    free(splitCheckList);
    totalMatchCnt += matchBuffer->startIndexOfReserve;
    return true;
}

std::vector<QueryKmerSplit> KmerMatcher::makeQueryKmerSplits(const Buffer<Kmer> * qKmers) 
{
    size_t blankCnt = std::find_if(qKmers->buffer,
                                   qKmers->buffer + qKmers->startIndexOfReserve, 
                                   [](const auto& kmer) { return kmer.id != 0;}
                                  ) - qKmers->buffer;
    size_t queryKmerNum = qKmers->startIndexOfReserve - blankCnt;
    std::cout << "Query k-mer number     : " << queryKmerNum << endl;

    // Filter out meaningless target splits
    MmapedData<DiffIdxSplit> diffIdxSplits = mmapData<DiffIdxSplit>(diffIdxSplitFileName.c_str(), 3);
    size_t numOfDiffIdxSplits = diffIdxSplits.fileSize / sizeof(DiffIdxSplit);
    size_t numOfDiffIdxSplits_use = numOfDiffIdxSplits;
    for (size_t i = 1; i < numOfDiffIdxSplits; i++) {
        if (diffIdxSplits.data[i].ADkmer == 0 || diffIdxSplits.data[i].ADkmer == UINT64_MAX) {
            numOfDiffIdxSplits_use--;
        }
    }

    // Divide query k-mer list into blocks for multi threading.
    std::vector<QueryKmerSplit> querySplits;
    size_t quotient = queryKmerNum / par.threads;
    size_t remainder = queryKmerNum % par.threads;
    size_t startIdx = blankCnt;
    size_t endIdx = 0; // endIdx is inclusive
    for (size_t i = 0; i < par.threads; i++) {
        endIdx = startIdx + quotient - 1;
        if (remainder > 0) {
            endIdx++;
            remainder--;
        }
        bool needLastTargetBlock = true;
        uint64_t queryAA = AminoAcidPart(qKmers->buffer[startIdx].value);
        for (size_t j = 0; j < numOfDiffIdxSplits_use; j ++) {
            if (queryAA <= AminoAcidPart(diffIdxSplits.data[j].ADkmer)) {
                // Kmer qkmer(queryAA, 0);
                // Kmer tkmer (diffIdxSplits.data[j].ADkmer, 0);
                // cout << j << "\t" << startIdx << "\t";
                // qkmer.printAA(*geneticCode, 12); cout << "\t";
                // tkmer.printAA(*geneticCode, 12); cout << endl;
                querySplits.emplace_back(startIdx, endIdx, diffIdxSplits.data[j - (j != 0)]);
                needLastTargetBlock = false;
                break;
            }
        }
        if (needLastTargetBlock) {
            querySplits.emplace_back(startIdx, endIdx, diffIdxSplits.data[numOfDiffIdxSplits_use - 2]);
        }
        startIdx = endIdx + 1;
    }

    for (size_t i = 0; i < querySplits.size(); i++) {
        std::cout << "Split " << i << ": Query k-mer index " << querySplits[i].start << " - " << querySplits[i].end
                  << ", Target k-mer offset " << querySplits[i].diffIdxSplit.ADkmer
                  << " (diffIdx: " << querySplits[i].diffIdxSplit.diffIdxOffset
                  << ", infoIdx: " << querySplits[i].diffIdxSplit.infoIdxOffset << ")" << std::endl;
    }

    munmap(diffIdxSplits.data, diffIdxSplits.fileSize);
    return querySplits;
}


bool KmerMatcher::matchKmers2(
    const Buffer<Kmer> * queryKmerBuffer,
    Buffer<Match> * totalMatches,
    const string & db)
{
    std::vector<QueryKmerSplit> querySplits = makeQueryKmerSplits(queryKmerBuffer);
    const Kmer * qKmers = queryKmerBuffer->buffer;
    bool *splitCheckList = (bool *) calloc(querySplits.size(), sizeof(bool));
    size_t totalOverFlowCnt = 0;
    unsigned int mask = ~((static_cast<unsigned int>(par.skipRedundancy == 0) << 31)); 

    std::cout << "Reference k-mer search : " << flush;
    time_t beforeSearch = time(nullptr);
#pragma omp parallel default(none), shared(splitCheckList, totalOverFlowCnt, \
querySplits, qKmers, totalMatches, cout, mask)
{
    Buffer<Match> localMatches(1024 * 1024 * 2);
    DeltaIdxReader * deltaIdxReaders = new DeltaIdxReader(targetDiffIdxFileName,
                                                    targetInfoFileName, 
                                                    1024 * 1024 * 32, 1024 * 1024);
    std::vector<Kmer> candidates;
    std::vector<Match> filteredMatches;
    bool hasOverflow = false;
#pragma omp for schedule(dynamic, 1)
    for (size_t i = 0; i < querySplits.size(); i++) {
        if (totalOverFlowCnt > 0 || splitCheckList[i]) { continue; } 

        deltaIdxReaders->setReadPosition(querySplits[i].diffIdxSplit);
        Kmer tKmer = deltaIdxReaders->next();
        uint64_t nextOffsetKmer = (i == querySplits.size() - 1)
                         ? UINT64_MAX
                         : querySplits[i + 1].diffIdxSplit.ADkmer;
        
        Kmer qKmer(UINT64_MAX, 0);
        uint64_t qKmerAA = UINT64_MAX;
        for (size_t j = querySplits[i].start; j < querySplits[i].end + 1; j++) {
            // Reuse the filtered matches if queries are exactly identical
            if (qKmer.value == qKmers[j].value && (qKmer.qInfo.frame/3 == qKmers[j].qInfo.frame/3)) {
                size_t filteredMatchCnt = filteredMatches.size();
                if (unlikely(!localMatches.afford(filteredMatchCnt))) {
                    if (!Buffer<Match>::moveSmallToLarge(&localMatches, totalMatches)) {
                        hasOverflow = true;
                        __sync_fetch_and_add(&totalOverFlowCnt, 1);
                        break;
                    }
                }
                size_t posToWrite = localMatches.reserveMemory(filteredMatchCnt);
                memcpy(localMatches.buffer + posToWrite, filteredMatches.data(),
                       sizeof(Match) * filteredMatchCnt);
                for (size_t k = 0; k < filteredMatchCnt; k++) {
                    localMatches.buffer[posToWrite + k].qInfo = qKmers[j].qInfo;
                }
                continue;
            }

            // Reuse the candidate target k-mers to compare in DNA level if queries are the same at amino acid level but not at DNA level
            if (qKmerAA == AminoAcidPart(qKmers[j].value)) {
                filterCandidates(qKmers[j], candidates, filteredMatches);
                size_t filteredMatchCnt = filteredMatches.size();
                if (unlikely(!localMatches.afford(filteredMatchCnt))) {
                    if (!Buffer<Match>::moveSmallToLarge(&localMatches, totalMatches)) {
                        hasOverflow = true;
                        __sync_fetch_and_add(&totalOverFlowCnt, 1);
                        break;
                    }
                }
                size_t posToWrite = localMatches.reserveMemory(filteredMatchCnt);
                memcpy(localMatches.buffer + posToWrite, filteredMatches.data(),
                       sizeof(Match) * filteredMatchCnt);
                for (size_t k = 0; k < filteredMatchCnt; k++) {
                    localMatches.buffer[posToWrite + k].qInfo = qKmers[j].qInfo;
                }
                qKmer = qKmers[j];
                qKmerAA = AminoAcidPart(qKmer.value);
                continue;
            }
            candidates.clear();

            // Get next query, and start to find
            qKmer = qKmers[j];
            qKmerAA = AminoAcidPart(qKmer.value);

            // Skip target k-mers lexiocographically smaller at amino acid level
            while (tKmer.value < nextOffsetKmer 
                    && !deltaIdxReaders->isCompleted() 
                    && qKmerAA > AminoAcidPart(tKmer.value)) {
                tKmer = deltaIdxReaders->next();
            }

            // No match found - skip to the next query
            if (qKmerAA != AminoAcidPart(tKmer.value)) {
                continue;
            } 
                    
            // Match found - load target k-mers matching at amino acid level
            while (tKmer.value < nextOffsetKmer 
                    && !deltaIdxReaders->isCompleted() 
                    && qKmerAA == AminoAcidPart(tKmer.value)) {
                candidates.emplace_back(tKmer.value, tKmer.id & mask);
                tKmer = deltaIdxReaders->next();                                      
            }

            filterCandidates(qKmer, candidates, filteredMatches);
            if (unlikely(!localMatches.afford(filteredMatches.size()))) {
                if (!Buffer<Match>::moveSmallToLarge(&localMatches, totalMatches)) {
                    hasOverflow = true;
                    __sync_fetch_and_add(&totalOverFlowCnt, 1);
                    break;
                }
            }

            size_t posToWrite = localMatches.reserveMemory(filteredMatches.size());
            memcpy(localMatches.buffer + posToWrite, filteredMatches.data(),
                   sizeof(Match) * filteredMatches.size());
        } // End of one split

        // Move matches in the local buffer to the shared buffer
        if (!Buffer<Match>::moveSmallToLarge(&localMatches, totalMatches)) {
            hasOverflow = true;
            __sync_fetch_and_add(&totalOverFlowCnt, 1);
        }

        // Check whether current split is completed or not
        if (!hasOverflow) {
            splitCheckList[i] = true;
        }
    } // End of omp for (Iterating for splits)
} // End of omp parallel
        
    if (totalOverFlowCnt > 0) {
        return false;
    }
    std::cout << double(time(nullptr) - beforeSearch) << " s" << std::endl;
    free(splitCheckList);
    totalMatchCnt += totalMatches->startIndexOfReserve;
    return true;
}


bool KmerMatcher::matchKmers_AA(
    const Buffer<Kmer> * queryKmerBuffer,
    Buffer<Match_AA> * totalMatches,
    const string & db)
{
    std::vector<QueryKmerSplit> querySplits = makeQueryKmerSplits(queryKmerBuffer);
    const Kmer * qKmers = queryKmerBuffer->buffer;
    std::atomic<int> totalOverFlowCnt{0};
    #pragma omp parallel default(none), shared(totalOverFlowCnt, querySplits, qKmers, totalMatches, cout)
    {
        Buffer<Match_AA> localMatches(1024 * 1024 * 2);  // 16 Mb
        DeltaIdxReader * deltaIdxReaders 
            = new DeltaIdxReader(targetDiffIdxFileName,
                                 targetInfoFileName, 
                                 1024 * 1024, 1024 * 1024);
        std::vector<Match_AA> tempMatches;  
        bool hasOverflow = false;
    
        #pragma omp for schedule(dynamic, 1)
        for (size_t i = 0; i < querySplits.size(); i++) {
            if (totalOverFlowCnt.load() > 0) { continue; }
            deltaIdxReaders->setReadPosition(querySplits[i].diffIdxSplit);
            Kmer tKmer = deltaIdxReaders->next();
            uint64_t nextOffsetKmer = (i == querySplits.size() - 1)
                             ? UINT64_MAX
                             : querySplits[i + 1].diffIdxSplit.ADkmer;

            Kmer qKmer(UINT64_MAX, 0);
            for (size_t j = querySplits[i].start; j < querySplits[i].end + 1; j++) {
                // Reuse the AA matches if queries are identical
                if (qKmer.value == qKmers[j].value) {
                    if (unlikely(!localMatches.afford(tempMatches.size()))) {
                        if (!Buffer<Match_AA>::moveSmallToLarge(&localMatches, totalMatches)) {
                            hasOverflow = true;
                            totalOverFlowCnt++;
                            break;
                        }
                    }
                    size_t posToWrite = localMatches.reserveMemory(tempMatches.size());
                    memcpy(localMatches.buffer + posToWrite, tempMatches.data(),
                           sizeof(Match) * tempMatches.size());
                    for (size_t k = 0; k < tempMatches.size(); k++) {
                        localMatches.buffer[posToWrite + k].queryId = qKmers[j].qInfo.sequenceID;
                        localMatches.buffer[posToWrite + k].pos = qKmers[j].qInfo.pos;
                    }
                    continue;
                }
                tempMatches.clear();
                // Get next query, and start to find
                qKmer = qKmers[j];

                // Skip target k-mers lexiocographically smaller
                while (!deltaIdxReaders->isCompleted() && qKmer.value > tKmer.value) {
                    tKmer = deltaIdxReaders->next();
                }

                // No match found - skip to the next query
                if (qKmer.value != tKmer.value) { continue; } 

                // Match found - load target k-mers matching at amino acid level
                while (!deltaIdxReaders->isCompleted() && qKmer.value == tKmer.value) {
                    tempMatches.emplace_back((uint32_t) qKmer.qInfo.sequenceID, tKmer.id, (uint32_t) qKmer.qInfo.pos, tKmer.value);
                    tKmer = deltaIdxReaders->next();                                      
                }

                if (unlikely(!localMatches.afford(tempMatches.size()))) {
                    if (!Buffer<Match_AA>::moveSmallToLarge(&localMatches, totalMatches)) {
                        hasOverflow = true;
                        totalOverFlowCnt++;
                        break;
                    }
                }

                size_t posToWrite = localMatches.reserveMemory(tempMatches.size());
                memcpy(localMatches.buffer + posToWrite, tempMatches.data(),
                       sizeof(Match_AA) * tempMatches.size());
            } // End of one split

            // Move matches in the local buffer to the shared buffer
            if (!Buffer<Match_AA>::moveSmallToLarge(&localMatches, totalMatches)) {
                hasOverflow = true;
                totalOverFlowCnt++;
            }
        } // End of omp for (Iterating for splits)
    } // End of omp parallel
        
    if (totalOverFlowCnt.load() > 0) {
        return false;
    }
    this->totalMatchCnt += totalMatches->startIndexOfReserve;
    return true;
}


bool KmerMatcher::matchMetamers(Buffer<Kmer> * queryKmerBuffer,
                                Buffer<Match> * matchBuffer,
                                const string & db){
    std::cout << "Comparing query and reference metamers..." << std::endl;
    targetDiffIdxFileName = dbDir + "/deltaIdx.mtbl";
    diffIdxSplitFileName = dbDir + "/deltaIdxSplits.mtbl";
    
    MmapedData<DeltaIdxOffset> diffIdxSplits = mmapData<DeltaIdxOffset>(diffIdxSplitFileName.c_str(), 3);
    size_t numOfDiffIdx = FileUtil::getFileSize(targetDiffIdxFileName) / sizeof(uint16_t);
    size_t queryKmerNum = queryKmerBuffer->startIndexOfReserve;
    Kmer * queryKmerList = queryKmerBuffer->buffer;

    size_t blankCnt = 0;
    for (size_t i = 0; i < queryKmerNum; i++) {
        if (queryKmerList[i].qInfo.sequenceID == 0) {
            blankCnt++;
        } else {
            break;
        }
    }
    queryKmerNum -= blankCnt;
    cout << "Total query k-mers: " << queryKmerNum << endl;
    // Filter out meaningless target splits
    size_t numOfDiffIdxSplits = diffIdxSplits.fileSize / sizeof(DiffIdxSplit);
    size_t numOfDiffIdxSplits_use = numOfDiffIdxSplits;
    for (size_t i = 1; i < numOfDiffIdxSplits; i++) {
        if (diffIdxSplits.data[i].metamer.metamer == 0 || diffIdxSplits.data[i].metamer.metamer == UINT64_MAX) {
            diffIdxSplits.data[i].metamer = {UINT64_MAX, UINT32_MAX};
            numOfDiffIdxSplits_use--;
        }
    }

    // Divide query k-mer list into blocks for multi threading.
    // Each split has start and end points of query list + proper offset point of target k-mer list
    std::vector<QueryKmerSplit2> querySplits;
    uint64_t queryAA;
    size_t quotient = queryKmerNum / threads;
    size_t remainder = queryKmerNum % threads;
    size_t startIdx = blankCnt;
    size_t endIdx = 0; // endIdx is inclusive
    for (size_t i = 0; i < threads; i++) {
        endIdx = startIdx + quotient - 1;
        if (remainder > 0) {
            endIdx++;
            remainder--;
        }
        bool needLastTargetBlock = true;
        queryAA = AMINO_ACID_PART(queryKmerList[startIdx].value);
        for (size_t j = 0; j < numOfDiffIdxSplits_use; j ++) {
            if (queryAA <= AMINO_ACID_PART(diffIdxSplits.data[j].metamer.metamer)) {
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
    int redundancyStored = par.skipRedundancy == 0;
    unsigned int mask = ~((static_cast<unsigned int>(par.skipRedundancy == 0) << 31));
#pragma omp parallel default(none), shared(splitCheckList, totalOverFlowCnt, \
querySplits, queryKmerList, matchBuffer, cout, mask, targetDiffIdxFileName, numOfDiffIdx, redundancyStored, targetInfoFileName)
{
    FILE * diffIdxFp = fopen(targetDiffIdxFileName.c_str(), "rb");

    // Target K-mer buffer
    uint16_t * diffIdxBuffer = (uint16_t *) malloc(sizeof(uint16_t) * (BufferSize + 1)); // size = 32 Mb
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
    Metamer currentTargetKmer;

    // Match buffer for each thread
    size_t localBufferSize = 2'000'000; // 32 Mb
    auto *matches = new Match[localBufferSize]; // 16 * 2'000'000 = 32 Mb
    size_t matchCnt = 0;

    // Vectors for selected target k-mers
    std::vector<uint8_t> selectedHammingSum;
    std::vector<size_t> selectedMatches;
    std::vector<uint16_t> selectedHammings;
    selectedHammingSum.resize(1024);
    selectedMatches.resize(1024);
    selectedHammings.resize(1024);
    size_t selectedMatchCnt = 0;

    size_t posToWrite;
    size_t idx;
    bool hasOverflow = false;

#pragma omp for schedule(dynamic, 1)
    for (size_t i = 0; i < querySplits.size(); i++) {
        if (totalOverFlowCnt > 0 || splitCheckList[i]) {
            continue;
        }
        currentTargetKmer = querySplits[i].deltaIdxOffset.metamer;
        diffIdxBufferIdx = querySplits[i].deltaIdxOffset.offset;
        diffIdxPos = querySplits[i].deltaIdxOffset.offset;

        fseek(diffIdxFp, 2 * (long) (diffIdxBufferIdx), SEEK_SET);
        loadBuffer(diffIdxFp, diffIdxBuffer, diffIdxBufferIdx, BufferSize);
                
        if (querySplits[i].deltaIdxOffset.metamer.metamer == 0 && querySplits[i].deltaIdxOffset.offset == 0) {
            currentTargetKmer = getNextTargetKmer(currentTargetKmer, diffIdxBuffer,
                                                  diffIdxBufferIdx, diffIdxPos);
        }
        currentQuery = UINT64_MAX;
        currentQueryAA = UINT64_MAX;
        for (size_t j = querySplits[i].start; j < querySplits[i].end + 1; j++) {
            // Reuse the comparison data if queries are exactly identical
            if (currentQuery == queryKmerList[j].value
                && (currentQueryInfo.frame/3 == queryKmerList[j].qInfo.frame/3)) {
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
                }
                for (size_t k = 0; k < selectedMatchCnt; k++) {
                    idx = selectedMatches[k];
                    matches[matchCnt] = {queryKmerList[j].qInfo,
                                         candidateKmerInfos[idx],
                                         taxId2speciesId[candidateKmerInfos[idx]],
                                         (unsigned int) (candidateTargetKmers[idx] & AA_MASK),
                                         selectedHammings[k],
                                         selectedHammingSum[k]};
                    matchCnt++;
                }
                continue;
            }
            selectedMatchCnt = 0;

            // Reuse the candidate target k-mers to compare in DNA level if queries are the same at amino acid level but not at DNA level
            if (currentQueryAA == AMINO_ACID_PART(queryKmerList[j].value)) {
                compareDna(queryKmerList[j].value, candidateTargetKmers, hammingDists, selectedMatches,
                           selectedHammingSum, selectedHammings, selectedMatchCnt, queryKmerList[j].qInfo.frame);
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
                }
                for (size_t k = 0; k < selectedMatchCnt; k++) {
                    idx = selectedMatches[k];
                    matches[matchCnt] = {queryKmerList[j].qInfo,
                                         candidateKmerInfos[idx],
                                         taxId2speciesId[candidateKmerInfos[idx]],
                                         (unsigned int) (candidateTargetKmers[idx] & AA_MASK),
                                         selectedHammings[k],
                                         selectedHammingSum[k]};
                    matchCnt++;
                }
                currentQuery = queryKmerList[j].value;
                currentQueryAA = AMINO_ACID_PART(currentQuery);
                currentQueryInfo = queryKmerList[j].qInfo;
                continue;
            }
            candidateTargetKmers.clear();
            candidateKmerInfos.clear();

            // Get next query, and start to find
            currentQuery = queryKmerList[j].value;
            currentQueryAA = AMINO_ACID_PART(currentQuery);
            currentQueryInfo = queryKmerList[j].qInfo;

            // Skip target k-mers that are not matched in amino acid level
            while (diffIdxPos != numOfDiffIdx
                   && (currentQueryAA > AMINO_ACID_PART(currentTargetKmer.metamer))) {  
                if (unlikely(BufferSize < diffIdxBufferIdx + 10)){
                    loadBuffer(diffIdxFp, diffIdxBuffer, diffIdxBufferIdx, BufferSize, ((int)(BufferSize - diffIdxBufferIdx)) * -1 );
                }
                currentTargetKmer = getNextTargetKmer(currentTargetKmer, diffIdxBuffer,
                                                      diffIdxBufferIdx, diffIdxPos);
            }

            if (currentQueryAA != AMINO_ACID_PART(currentTargetKmer.metamer)) {
                continue;
            } 
                    
            // Load target k-mers that are matched in amino acid level
            while (diffIdxPos != numOfDiffIdx &&
                   currentQueryAA == AMINO_ACID_PART(currentTargetKmer.metamer)) {
                candidateTargetKmers.push_back(currentTargetKmer.metamer);
                if (redundancyStored) {
                    candidateKmerInfos.push_back(currentTargetKmer.id & mask);
                } else {
                    candidateKmerInfos.push_back(currentTargetKmer.id);
                }
                if (unlikely(BufferSize < diffIdxBufferIdx + 10)){
                    loadBuffer(diffIdxFp, diffIdxBuffer, diffIdxBufferIdx, BufferSize, ((int)(BufferSize - diffIdxBufferIdx)) * -1 );
                }
                currentTargetKmer = getNextTargetKmer(currentTargetKmer, diffIdxBuffer,
                                                      diffIdxBufferIdx, diffIdxPos);
            }

            if (candidateTargetKmers.size() > selectedMatches.size()) {
                selectedMatches.resize(candidateTargetKmers.size());
                selectedHammingSum.resize(candidateTargetKmers.size());
                selectedHammings.resize(candidateTargetKmers.size());
            }

            // Compare the current query and the loaded target k-mers and select
            compareDna(currentQuery, candidateTargetKmers, hammingDists, selectedMatches, selectedHammingSum,
                       selectedHammings, selectedMatchCnt, queryKmerList[j].qInfo.frame);

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
            }
            for (size_t k = 0; k < selectedMatchCnt; k++) {
                idx = selectedMatches[k];
                matches[matchCnt] = {queryKmerList[j].qInfo,
                                     candidateKmerInfos[idx],
                                     taxId2speciesId[candidateKmerInfos[idx]],
                                     (unsigned int) (candidateTargetKmers[idx] & AA_MASK),
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
    free(diffIdxBuffer);
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
    std::cout << "K-mer match sorting    : " << std::flush;
    SORT_PARALLEL(matchBuffer->buffer,
                  matchBuffer->buffer + matchBuffer->startIndexOfReserve,
                  compareMatches);
    std::cout << double(time(nullptr) - beforeSortMatches) << " s" << std::endl;
}

void KmerMatcher::moveMatches(Match *dest, Match *src, size_t & matchNum) {
    memcpy(dest, src, sizeof(Match) * matchNum);
    matchNum = 0;
}

void KmerMatcher::filterCandidates(
    Kmer qKmer,
    const std::vector<Kmer> &candidates,
    std::vector<Match> &filteredMatches
) {
    filteredMatches.clear();
    std::vector<uint8_t> hammingDists(candidates.size());
    uint8_t hDistCutoff = UINT8_MAX;
    for (size_t i = 0; i < candidates.size(); i++) {
        hammingDists[i] = getHammingDistanceSum(qKmer.value, candidates[i].value);
        hDistCutoff = min(hDistCutoff, hammingDists[i]);
    }
    hDistCutoff = min(hDistCutoff * 2, 7);
    for (size_t h = 0; h < candidates.size(); h++) {
        if (hammingDists[h] <= hDistCutoff) {
            uint16_t hDists = !((qKmer.qInfo.frame < 3) ^ (kmerFormat == 2))
                ? getHammings(qKmer.value, candidates[h].value)
                : getHammings_reverse(qKmer.value, candidates[h].value);
            filteredMatches.emplace_back(
                qKmer.qInfo,                                    // query k-mer info
                candidates[h].id,                               // target TaxID
                taxId2speciesId[candidates[h].id],              // target species ID
                (unsigned int) (candidates[h].value & AA_MASK), // target DNA encoding
                hDists,                                         // Hamming dist. per codon                        
                hammingDists[h]                                 // Hamming dist. sum
            );
        }
    }
}

// It compares query k-mers to target k-mers.
// If a query has matches, the matches with the smallest hamming distance will be selected
void KmerMatcher::compareDna(uint64_t query,
                             std::vector<uint64_t> &targetKmersToCompare,
                             std::vector<uint8_t> & hammingDists,
                             std::vector<size_t> &selectedMatches,
                             std::vector<uint8_t> &selectedHammingSum,
                             std::vector<uint16_t> &selectedHammings,
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
            selectedHammings[selectedMatchIdx] = !((frame < 3) ^ (kmerFormat == 2))
                ? getHammings(query, targetKmersToCompare[h])
                : getHammings_reverse(query, targetKmersToCompare[h]);
            selectedMatches[selectedMatchIdx++] = h;
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
