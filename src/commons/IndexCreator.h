#ifndef ADKMER4_INDEXCREATOR_H
#define ADKMER4_INDEXCREATOR_H
#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include <ctime>
#include <fstream>
#include <unordered_set>
#include <atomic>
#include <cstdint>
#include <iomanip>
#ifdef OPENMP
    #include <omp.h>
#endif


#include "printBinary.h"
#include "Mmap.h"
#include "Kmer.h"
#include "SeqIterator.h"
#include "TaxonomyWrapper.h"
#include "BitManipulateMacros.h"
#include "common.h"
#include "FastSort.h"
#include "LocalParameters.h"
#include "NucleotideMatrix.h"
#include "SubstitutionMatrix.h"
#include "tantan.h"
#include "LocalUtil.h"
#include "TaxonomyWrapper.h"
#include "fasta_validate.h"
#include "GeneticCode.h"
#include "KmerExtractor.h"
#include "DeltaIdxReader.h"
#include "UnirefTree.h"


enum class FilterMode { DB_CREATION, COMMON_KMER, UNIQ_KMER, UNIREF_LCA};
struct Accession {
    Accession() = default;
    Accession(const string & accession, uint32_t whichFasta, uint32_t order, uint32_t length) 
        : accession(accession), whichFasta(whichFasta), order(order), length(length), speciesID(0), taxID(0) {}
    string accession;
    uint32_t whichFasta;
    uint32_t order;
    uint32_t length;
    TaxID speciesID;
    TaxID taxID;

    void print() {
        std::cout << accession << " " << whichFasta << " " << order << " " << length << " " << speciesID << " " << taxID << std::endl;
    }

    bool operator < (const Accession & a) const {
        if (speciesID != a.speciesID)
            return speciesID < a.speciesID;
        
        if (whichFasta != a.whichFasta)
            return whichFasta < a.whichFasta;

        if (order != a.order)
            return order < a.order;
        
        return false;
    }

    static bool compare(const Accession & a, const Accession & b) {
        if (a.speciesID != b.speciesID)
            return a.speciesID < b.speciesID;
        
        if (a.whichFasta != b.whichFasta)
            return a.whichFasta < b.whichFasta;

        if (a.order != b.order)
            return a.order < b.order;
            
        return false;
    }
};

struct AccessionBatch {
    uint32_t whichFasta;
    TaxID speciesID;
    uint32_t trainingSeqFasta;
    uint32_t trainingSeqIdx;
    uint64_t totalLength;
    vector<uint32_t> orders;
    vector<TaxID> taxIDs;
    vector<uint32_t> lengths;

    void print() const {
        std::cout << "whichFasta: " << whichFasta << " speciesID: " << speciesID << " trainingSeqFasta: " << trainingSeqFasta << " trainingSeqIdx: " << trainingSeqIdx << endl;
        for (size_t i = 0; i < orders.size(); ++i) {
            std::cout << "order: " << orders[i] << " taxID: " << taxIDs[i] << " length: " << lengths[i] << endl;
        }
    }

    AccessionBatch(uint32_t whichFasta, TaxID speciesID, uint32_t trainingSeqFasta, uint32_t trainingSeqIdx, uint64_t totalLength) 
        : whichFasta(whichFasta), speciesID(speciesID), trainingSeqFasta(trainingSeqFasta), trainingSeqIdx(trainingSeqIdx), totalLength(totalLength) {}
};

using namespace std;

class IndexCreator{
protected:
    // Parameters
    const LocalParameters & par;
    bool isUpdating;
    int kmerFormat;
    int kmerLen;

    uint64_t MARKER;
    BaseMatrix *subMat;
    bool removeRedundancyInfo;
    unordered_map<TaxID, TaxID> taxId2speciesId;

    // Inputs
    TaxonomyWrapper * taxonomy = nullptr;
    UnirefTree * unirefTree = nullptr;
    GeneticCode * geneticCode;
    KmerExtractor * kmerExtractor;

    bool externTaxonomy;
    string dbDir;
    string fnaListFileName;
    string acc2taxidFileName;

    // Outputs
    string taxidListFileName;
    string taxonomyBinaryFileName;
    string versionFileName;
    string paramterFileName;

    std::unordered_map<string, vector<CDSinfo>> cdsInfoMap;
    std::vector<AccessionBatch> accessionBatches;
    std::unordered_set<TaxID> taxIdSet;
    vector<string> fastaPaths;
    size_t numOfFlush=0;

    // Database splits
    std::vector<std::string> deltaIdxFileNames;
    std::vector<std::string> infoFileNames;
    std::string mergedDeltaIdxFileName;
    std::string mergedInfoFileName;
    std::string deltaIdxSplitFileName;
    struct Split{
        Split(size_t offset, size_t end) : offset(offset), end(end) {}
        size_t offset;
        size_t end;

        void print() {
            std::cout << "offset: " << offset << " end: " << end << std::endl;
        }
    };

    void writeTargetFiles(
        Buffer<Kmer> & kmerBuffer,
        const size_t * uniqeKmerIdx,
        const vector<pair<size_t, size_t>> & uniqKmerIdxRanges);

    void writeTargetFilesAndSplits(
        Buffer<Kmer> & kmerBuffer,
        const size_t * uniqeKmerIdx, 
        size_t & uniqKmerCnt, 
        const vector<pair<size_t, size_t>> & uniqKmerIdxRanges,
        bool writeInfo = true);

    void writeDbParameters();

    size_t fillTargetKmerBuffer(
        Buffer<Kmer> &kmerBuffer,                 
        std::vector<std::atomic<bool>> & batchChecker,
        size_t &processedSplitCnt,
        const LocalParameters &par);
    
    bool extractKmerFromSixFrames(
        Buffer<Kmer> & kmerBuffer,
        std::vector<std::atomic<bool>> & batchChecker,
        size_t &processedBatchCnt);

    void indexReferenceSequences(size_t bufferSize);

    void getAccessionBatches(std::vector<Accession> & observedAccessionsVec, size_t bufferSize);

    void getObservedAccessions(const string & fnaListFileName,
                               vector<Accession> & observedAccessionsVec,
                               unordered_map<string, size_t> & accession2index);

    void getTaxonomyOfAccessions(vector<Accession> & observedAccessionsVec,
                                 const unordered_map<string, size_t> & accession2index,
                                 const string & acc2taxidFileName);

    void unzipAndList(const string & folder, const string & fastaList_fname){
        system(("./../../util/unzip_and_list.sh " + folder + " " + fastaList_fname).c_str());
    }

    void load_assacc2taxid(const string & mappingFile, unordered_map<string, int> & assacc2taxid);

    // TaxID getMaxTaxID();

    // void editTaxonomyDumpFiles(const vector<pair<string, pair<TaxID, TaxID>>> & newAcc2taxid);

    template <FilterMode M>
    bool areKmersDuplicate(
        const Kmer & kmer1,
        const Kmer & kmer2
    );

    size_t AminoAcidPart(size_t kmer) {
        if (kmerFormat == 3 || kmerFormat == 4) {
            return kmer;
        }
        return (kmer) & MARKER;
    }

    void loadCdsInfo(const string & cdsInfoFileList);

    size_t calculateBufferSize(size_t maxRam) {
        constexpr double GIGABYTE = 1024.0 * 1024.0 * 1024.0;
        constexpr double MEGABYTE = 1024.0 * 1024.0;
        constexpr double MEMORY_PER_THREAD_MB = 50.0;
        
        float c = 0.7f;
        if (maxRam < 16) {
            c = 0.5f;
        } else if (maxRam <= 32) {
            c = 0.6f;
        }
        double totalRamBytes = maxRam * GIGABYTE * c;
        double memoryForThreads = par.threads * MEMORY_PER_THREAD_MB * MEGABYTE;
        double availableMemory = totalRamBytes - memoryForThreads;

        if (availableMemory <= 0.0) {
            cerr << "Not enough memory to create index" << endl;
            cerr << "Please increase the RAM usage or decrease the number of threads" << endl;
            exit(EXIT_FAILURE);
        }

        size_t size_per_item = sizeof(Kmer) + sizeof(size_t);

        return static_cast<size_t>(availableMemory / size_per_item);
    }

    size_t calculateBufferSizeForMerge(size_t maxRam, int fileCnt) {
        float c = 0.7;
        if (maxRam <= 32) {
            c = 0.6;
        } else if (maxRam < 16) {
            c = 0.5;
        }
        if ((maxRam * 1024.0 * 1024.0 * 1024.0 * c - (par.threads * 50.0 * 1024.0 * 1024.0)) <= 0.0) {
            cerr << "Not enough memory to create index" << endl;
            cerr << "Please increase the RAM usage or decrease the number of threads" << endl;
            exit(EXIT_FAILURE);
        }
        return static_cast<size_t>((maxRam * 1024.0 * 1024.0 * 1024.0 * c - (par.threads * 50.0 * 1024.0 * 1024.0))/ 
                                  (sizeof(Kmer) + sizeof(size_t)));
    }

    void loadMergedTaxIds(const std::string &mergedFile, unordered_map<TaxID, TaxID> & old2new);

    string addToLibrary(const std::string & dbDir,
                        const std::string & fileList,
                        const std::string & acc2taxIdFileName);



public:
    IndexCreator(const LocalParameters & par, TaxonomyWrapper * taxonomy, int kmerFormat);
    IndexCreator(const LocalParameters & par, UnirefTree * unirefTree, int kmerFormat);
    IndexCreator(const LocalParameters & par, int kmerFormat);
    ~IndexCreator();
    void createIndex();
    void createCommonKmerIndex();
    void createUniqueKmerIndex();
    void createLcaKmerIndex();

    template <FilterMode M>
    void mergeTargetFiles();

    template <FilterMode M>
    void filterKmers(
        Buffer<Kmer> & kmerBuffer,
        size_t * uniqeKmerIdx,
        size_t & uniqKmerCnt,
        vector<pair<size_t, size_t>> & uniqKmerIdxRanges);

    static void getDiffIdx(
        uint64_t & lastKmer,
        uint64_t entryToWrite,
        WriteBuffer<uint16_t> & diffBuffer);

    // Getters
    int getNumOfFlush() const { return numOfFlush; }
    TaxonomyWrapper* getTaxonomy() const { return taxonomy;}
    unordered_set<TaxID> getTaxIdSet() { return taxIdSet; }

    // Setters
    void setIsUpdating(bool isUpdating) { this->isUpdating = isUpdating; }
    void setIsNewFormat(int kmerFormat) { this->kmerFormat = kmerFormat; }
    void addFilesToMerge(string diffIdxFileName, string infoFileName);
    void updateTaxId2SpeciesTaxId(const string & taxIdListFileName);
    void setMergedFileNames(string diffFileName, string infoFileName, string splitFileName);

    void printFilesToMerge() {
        cout << "Files to merge :" << endl;
        for (size_t i = 0; i < deltaIdxFileNames.size(); i++) {
            cout << deltaIdxFileNames[i] << " " << infoFileNames[i] << endl;
        }
    }
    
    static void printIndexSplitList(DiffIdxSplit * splitList) {
        for (int i = 0; i < 4096; i++) {
            std::cout << splitList[i].infoIdxOffset << " " << 
                    splitList[i].diffIdxOffset << " " << 
                    splitList[i].ADkmer << std::endl;
        }
    }
};


template <FilterMode M>
void IndexCreator::mergeTargetFiles() {
    size_t bufferSize = 1024 * 1024 * 512;
    WriteBuffer<uint16_t> diffBuffer(mergedDeltaIdxFileName, bufferSize);
    WriteBuffer<uint32_t> infoBuffer(mergedInfoFileName, bufferSize);
    
    // Prepare files to merge
    size_t numOfKmerBeforeMerge = 0;
    size_t splitNum = deltaIdxFileNames.size();
    DeltaIdxReader ** deltaIdxReaders = new DeltaIdxReader*[splitNum];
    size_t valueBufferSize = 1024 * 1024 * 16;
    for (size_t file = 0; file < splitNum; file++) {
        deltaIdxReaders[file] = new DeltaIdxReader(deltaIdxFileNames[file],
                                                   infoFileNames[file],
                                                   valueBufferSize, 
                                                   1024 * 1024 * 4); 
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

    cout << "Merging target files..." << endl;
    size_t kmerBufferSize = 1024 * 1024 * 1024;
    if (kmerBufferSize < valueBufferSize * splitNum * 2) {
        kmerBufferSize = valueBufferSize * splitNum * 2; 
    }
    Buffer<Kmer> kmerBuffer(kmerBufferSize); // 64GB
    std::atomic<int> hasOverflow{0};
    std::vector<std::atomic<bool>> completedSplits(splitNum);
    int remainingSplits = splitNum;
    vector<pair<size_t, size_t>> uniqKmerIdxRanges;
    size_t lastKmer = 0;
    auto * uniqKmerIdx = new size_t[kmerBuffer.bufferSize];
    vector<size_t> splitToProcess;
    while (remainingSplits > 0) {
        kmerBuffer.init();
        memset(uniqKmerIdx, 0, kmerBuffer.bufferSize * sizeof(size_t));
        time_t start = time(nullptr);
        while (remainingSplits > 0
                && kmerBuffer.startIndexOfReserve + valueBufferSize * remainingSplits <= kmerBuffer.bufferSize) {
            size_t posToWrite = kmerBuffer.reserveMemory(valueBufferSize * remainingSplits);
            uint64_t max = UINT64_MAX;
            splitToProcess.clear();
            for (size_t i = 0; i < splitNum; i++) {
                if (completedSplits[i].load(std::memory_order_acquire))
                    continue; 
                if (deltaIdxReaders[i]->getLastValue() < max) {
                    max = deltaIdxReaders[i]->getLastValue();
                }
                splitToProcess.push_back(i);
            }
#pragma omp parallel for default(none), shared(cout, kmerBuffer, deltaIdxReaders, splitToProcess, completedSplits, posToWrite, max, splitNum, valueBufferSize, taxId2speciesId, mask, remainingSplits)
            for (size_t i = 0; i < splitToProcess.size(); i ++) {
                size_t split = splitToProcess[i];
                size_t offset = posToWrite + i * valueBufferSize;
                size_t valueNum = deltaIdxReaders[split]->getValues(kmerBuffer.buffer + offset, max);
                if (deltaIdxReaders[split]->isCompleted()) {
                    completedSplits[split].store(true, std::memory_order_release);
                    __sync_fetch_and_sub(&remainingSplits, 1);
                }
                if constexpr (M == FilterMode::DB_CREATION || M == FilterMode::COMMON_KMER) {            
                    for (size_t j = 0; j < valueNum; j++) {
                        kmerBuffer.buffer[offset + j].tInfo.speciesId 
                            = taxId2speciesId[kmerBuffer.buffer[offset + j].tInfo.taxId & mask];
                    }
                }
            }
        }
        time_t end = time(nullptr);
        cout << "K-mer loading       : " << (double) (end - start) << " s" << endl;

        time_t beforeSort = time(nullptr);

        if constexpr (M == FilterMode::DB_CREATION || M == FilterMode::COMMON_KMER) {
            SORT_PARALLEL(kmerBuffer.buffer, 
                          kmerBuffer.buffer + kmerBuffer.startIndexOfReserve,
                          Kmer::compareTargetKmer);
        } else if constexpr (M == FilterMode::UNIQ_KMER || M == FilterMode::UNIREF_LCA) {
            SORT_PARALLEL(kmerBuffer.buffer, 
                          kmerBuffer.buffer + kmerBuffer.startIndexOfReserve,
                          Kmer::compareKmer);
        }
        
        time_t afterSort = time(nullptr);
        cout << "Sorting k-mer list  : " << afterSort - beforeSort << " s" << endl;

        // Reduce redundancy
        size_t uniqKmerCnt = 0;
        uniqKmerIdxRanges.clear();
        filterKmers<M>(kmerBuffer, uniqKmerIdx, uniqKmerCnt, uniqKmerIdxRanges);
        time_t reduction = time(nullptr);
        cout << "Filtering k-mers    : " << (double) (reduction - afterSort) << " s" << endl;
        cout << "Selected count      : " << uniqKmerCnt << endl;

        // Write       
        for (size_t i = 0; i < uniqKmerIdxRanges.size(); i ++) {
            for (size_t j = uniqKmerIdxRanges[i].first; j < uniqKmerIdxRanges[i].second; j ++) {
                infoBuffer.write(&kmerBuffer.buffer[uniqKmerIdx[j]].id);
                getDiffIdx(lastKmer, kmerBuffer.buffer[uniqKmerIdx[j]].value, diffBuffer);

                // Write split info
                if (AminoAcidPart(lastKmer) != AAofTempSplitOffset && splitCheck == 1) {
                    splitList[splitListIdx++] = {lastKmer, diffBuffer.writeCnt, infoBuffer.writeCnt};
                    splitCheck = 0;
                }
                if (infoBuffer.writeCnt == offsetList[offsetListIdx]) {
                    AAofTempSplitOffset = AminoAcidPart(lastKmer);
                    splitCheck = 1;
                    offsetListIdx++;
                }
            }
        }  
        time_t writeTime = time(nullptr);
        cout << "Writing k-mers      : " << (double) (writeTime - reduction) << " s" << endl;
        cout << "Written k-mers      : " << infoBuffer.writeCnt << " " << std::fixed << std::setprecision(2) << ((float) infoBuffer.writeCnt / numOfKmerBeforeMerge) * 100 << "%" << endl;
        cout << "--------------------" << endl;
        cout << "remaining splits    : " << remainingSplits << endl;
        std::cout.unsetf(std::ios::fixed | std::ios::scientific);
        kmerBuffer.startIndexOfReserve = 0;
    }

    FILE * diffIdxSplitFile = fopen(deltaIdxSplitFileName.c_str(), "wb");
    fwrite(splitList, sizeof(DiffIdxSplit), par.splitNum, diffIdxSplitFile);
    fclose(diffIdxSplitFile);
    // for(int i = 0; i < par.splitNum; i++) {
    //     cout<<splitList[i].ADkmer<< " "<<splitList[i].diffIdxOffset<< " "<<splitList[i].infoIdxOffset<<endl;
    // }
    cout<<"DB creation completed"<<endl;
    cout<<"Total k-mer count   : " << infoBuffer.writeCnt <<endl;

    cout<<"DB files you need   : " << endl;
    cout<<mergedDeltaIdxFileName<<endl;
    cout<<mergedInfoFileName<<endl;
    cout<<string(dbDir) + "/taxID_list"<<endl;
    cout<<deltaIdxSplitFileName<<endl;
}


template <FilterMode M>
void IndexCreator::filterKmers(
    Buffer<Kmer> & kmerBuffer,
    size_t * selectedKmerIdx,
    size_t & selectedKmerCnt,
    vector<pair<size_t, size_t>> & selectedKmerIdxRanges) 
{
    // Find the first index of garbage k-mer (UINT64_MAX)
    for(size_t checkN = kmerBuffer.startIndexOfReserve - 1; checkN != 0; checkN--){
        if(kmerBuffer.buffer[checkN].value != UINT64_MAX){
            kmerBuffer.startIndexOfReserve = checkN + 1;
            break;
        }
    }

    // Find the first index of meaningful k-mer
    size_t startIdx = 0;
    for (size_t i = 0; i < kmerBuffer.startIndexOfReserve ; i++) {
        if(!kmerBuffer.buffer[i].isEmpty()){
            startIdx = i;
            break;
        }
    }
    cout << "Loaded k-mer count  : " << (kmerBuffer.startIndexOfReserve - startIdx) << endl;
    // Make splits
    vector<Split> splits;
    size_t splitWidth = (kmerBuffer.startIndexOfReserve - startIdx) / par.threads;
    for (int i = 0; i < par.threads - 1; i++) {
        for (size_t j = startIdx + splitWidth; j + 1 < kmerBuffer.startIndexOfReserve; j++) {
            if (AminoAcidPart(kmerBuffer.buffer[j].value) != AminoAcidPart(kmerBuffer.buffer[j + 1].value)) {
                splits.emplace_back(startIdx, j);
                selectedKmerIdxRanges.emplace_back(pair<size_t, size_t>(startIdx, 0));
                startIdx = j + 1;
                break;
            }
        }
    }
    splits.emplace_back(startIdx, kmerBuffer.startIndexOfReserve - 1);
    selectedKmerIdxRanges.emplace_back(pair<size_t, size_t>(startIdx, 0));

    size_t * cntOfEachSplit = new size_t[splits.size()];
    for(size_t i = 0; i < splits.size(); i++) {
        cntOfEachSplit[i] = 0;
    }
#pragma omp parallel default(none), shared(kmerBuffer, cntOfEachSplit, splits, par, cout, selectedKmerCnt, selectedKmerIdx, selectedKmerIdxRanges)
    {
        Kmer * lookingKmer;
        size_t lookingIndex;
        vector<TaxID> taxIds;
        vector<uint32_t> ids;
        size_t * localSelectedIdx = new size_t[16 * 1024 * 1024];
        size_t tempSelectedKmerCnt = 0;
#pragma omp for schedule(static, 1)
        for(size_t split = 0; split < splits.size(); split ++) {
            size_t i = splits[split].offset;
            i++;
            while (i < splits[split].end + 1) {
                lookingKmer = & kmerBuffer.buffer[i - 1];
                lookingIndex = i - 1;
                bool selected = false;
                if constexpr (M == FilterMode::DB_CREATION) {
                    taxIds.clear();
                    taxIds.push_back(lookingKmer->tInfo.taxId);
                } else if constexpr (M == FilterMode::COMMON_KMER) {
                    taxIds.clear();
                    taxIds.push_back(lookingKmer->tInfo.speciesId);
                } else if constexpr (M == FilterMode::UNIQ_KMER || M == FilterMode::UNIREF_LCA) {
                    ids.clear();
                    ids.push_back(lookingKmer->id);
                }

                while ((i < splits[split].end + 1) && areKmersDuplicate<M>(*lookingKmer, kmerBuffer.buffer[i])) {
                    if constexpr (M == FilterMode::DB_CREATION) {
                        taxIds.push_back(kmerBuffer.buffer[i].tInfo.taxId);
                    } else if constexpr (M == FilterMode::COMMON_KMER) {
                        taxIds.push_back(kmerBuffer.buffer[i].tInfo.speciesId);
                    } else if constexpr (M == FilterMode::UNIQ_KMER || M == FilterMode::UNIREF_LCA) {
                        ids.push_back(kmerBuffer.buffer[i].id);
                    }
                    i++;
                }
                
                if constexpr (M == FilterMode::DB_CREATION || M == FilterMode::UNIREF_LCA) {
                    selected = true;
                } else if constexpr (M == FilterMode::COMMON_KMER) {
                    for (size_t i = 0; i < taxIds.size(); i++) {
                        if (taxIds[i] != lookingKmer->tInfo.speciesId) {
                            selected = true;
                            break;
                        }
                    }
                } else if constexpr (M == FilterMode::UNIQ_KMER) {
                    selected = true;
                    for (size_t i = 1; i < ids.size(); i++) {
                        if (ids[i] != lookingKmer->id) {
                            selected = false;
                            break;
                        }
                    }
                }

                if (selected) {
                    if constexpr (M == FilterMode::DB_CREATION || M == FilterMode::COMMON_KMER) {
                        lookingKmer->tInfo.taxId = taxonomy->LCA(taxIds)->taxId;
                    } else if constexpr (M == FilterMode::UNIREF_LCA) {
                        lookingKmer->id = unirefTree->getLCA(ids);
                    }

                    if (tempSelectedKmerCnt >= 16 * 1024 * 1024) {
                        memcpy(selectedKmerIdx + splits[split].offset, localSelectedIdx, tempSelectedKmerCnt * sizeof(size_t));
                        splits[split].offset += tempSelectedKmerCnt;
                        tempSelectedKmerCnt = 0;
                    }
                    localSelectedIdx[tempSelectedKmerCnt++] = lookingIndex;
                    cntOfEachSplit[split]++;
                }

                i++;                
            }

            // Check the last k-mer
            if constexpr (M == FilterMode::DB_CREATION || M == FilterMode::UNIQ_KMER || M == FilterMode::UNIREF_LCA) {
                if(!areKmersDuplicate<M>(kmerBuffer.buffer[splits[split].end - 1], kmerBuffer.buffer[splits[split].end])){
                    if (tempSelectedKmerCnt >= 16 * 1024 * 1024) {
                        memcpy(selectedKmerIdx + splits[split].offset, localSelectedIdx, tempSelectedKmerCnt * sizeof(size_t));
                        splits[split].offset += tempSelectedKmerCnt;
                        tempSelectedKmerCnt = 0;
                    }
                    localSelectedIdx[tempSelectedKmerCnt++] = splits[split].end;
                    cntOfEachSplit[split]++;
                }
            }
            memcpy(selectedKmerIdx + splits[split].offset, localSelectedIdx, tempSelectedKmerCnt * sizeof(size_t));
            splits[split].offset += tempSelectedKmerCnt;
            selectedKmerIdxRanges[split].second = splits[split].offset;
            __sync_fetch_and_add(&selectedKmerCnt, cntOfEachSplit[split]); 
        }
        delete[] localSelectedIdx;
    }
    delete[] cntOfEachSplit;
}

template <FilterMode M>
bool IndexCreator::areKmersDuplicate(
    const Kmer & kmer1,
    const Kmer & kmer2) 
{
    if constexpr (M == FilterMode::DB_CREATION) {
        return kmer1.tInfo.speciesId == kmer2.tInfo.speciesId &&
               kmer1.value == kmer2.value;
    } else if constexpr (M == FilterMode::COMMON_KMER || M == FilterMode::UNIQ_KMER || M == FilterMode::UNIREF_LCA) {
        return kmer1.value == kmer2.value;
    }
    return false; // Default case, should not be reached
}

#endif