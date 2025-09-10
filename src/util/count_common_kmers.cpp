#ifndef METABULI_COUNT_COMMON_KMERS_H
#define METABULI_COUNT_COMMON_KMERS_H

#include <iostream>
#include <stdint.h>

#include "LocalParameters.h"
#include "common.h"
#include "DeltaIdxReader.h"
#include "TaxonomyWrapper.h"   
#include "Mmap.h" 

void set_count_common_kmers_default(LocalParameters & par) {
    
}

int count_common_kmers(int argc, const char **argv, const Command &command) {
    std::cout << "Counting common k-mers by taxonomic rank..." << std::endl;
    LocalParameters &par = LocalParameters::getLocalInstance();
    par.parseParameters(argc, argv, command, true, Parameters::PARSE_ALLOW_EMPTY, 0);
    std::cout << "Parsing parameters..." << std::endl;

    std::string dbDir = par.filenames[0];
    std::string diffIdxFileName = dbDir + "/diffIdx";
    std::string infoFileName = dbDir + "/info";
    std::string diffIdxSplitFileName = dbDir + "/split";
    std::string taxonomyFileName = dbDir + "/taxonomyDB";

    std::cout << "Loading taxonomy from: " << taxonomyFileName << std::endl;
    TaxonomyWrapper * taxonomy = loadTaxonomy(dbDir);

    MmapedData<DiffIdxSplit> diffIdxSplits = mmapData<DiffIdxSplit>(diffIdxSplitFileName.c_str(), 3);
    
    size_t numOfDiffIdxSplits = diffIdxSplits.fileSize / sizeof(DiffIdxSplit);
    size_t numOfDiffIdxSplits_use = numOfDiffIdxSplits;
    for (size_t i = 1; i < numOfDiffIdxSplits; i++) {
        if (diffIdxSplits.data[i].ADkmer == 0 || diffIdxSplits.data[i].ADkmer == UINT64_MAX) {
            diffIdxSplits.data[i] = {UINT64_MAX, UINT64_MAX, UINT64_MAX};
            numOfDiffIdxSplits_use--;
        }
    }

    std::unordered_map<TaxID, size_t> taxon2Count;
    std::unordered_map<TaxID, size_t> taxon2UniqKmerCount;
    size_t distinctKmerCount = 0;
    GeneticCode geneticCode(false);
    // geneticCode = new GeneticCode(false);
#pragma omp parallel default(none), shared(std::cout, geneticCode, diffIdxFileName, infoFileName, diffIdxSplits, taxonomy, numOfDiffIdxSplits_use, taxon2Count, taxon2UniqKmerCount, distinctKmerCount)
    {
        DeltaIdxReader * deltaIdxReaders = new DeltaIdxReader(diffIdxFileName, infoFileName, 1024 * 1024 * 32, 1024 * 1024);
        std::vector<TaxID> taxIds; taxIds.reserve(4096);
        std::unordered_map<TaxID, size_t> localTaxon2Count;
        std::unordered_map<TaxID, size_t> localTaxon2UniqKmerCount;
        uint64_t mask = 16777215;
        mask = ~ mask;
        unsigned int mask2 = ~((static_cast<unsigned int>(1) << 31));
        Kmer nextKmer;
        size_t localDistinctKmerCount = 0;
    #pragma omp for schedule(dynamic, 1)
        for (size_t i = 0; i < numOfDiffIdxSplits_use; ++i) {
            DiffIdxSplit & currOffset = diffIdxSplits.data[i];
            uint64_t nextOffsetKmer;
            if (i != numOfDiffIdxSplits_use - 1) {
                nextOffsetKmer = diffIdxSplits.data[i + 1].ADkmer;
            } else {
                nextOffsetKmer = UINT64_MAX;
            }

            deltaIdxReaders->setReadPosition(currOffset);
            Kmer kmer = deltaIdxReaders->next();
            while (!deltaIdxReaders->isCompleted() &&
                (kmer.value < nextOffsetKmer)) 
            {
                taxIds.push_back(kmer.tInfo.taxId & mask2);
                while (((nextKmer = deltaIdxReaders->next()).isEmpty() == false)
                    && ((nextKmer.value) == (kmer.value))) 
                {
                    taxIds.push_back(nextKmer.tInfo.taxId & mask2);
                }
                TaxID lcaTaxId = taxonomy->LCA(taxIds)->taxId;
                localTaxon2Count[lcaTaxId] += taxIds.size();
                localTaxon2UniqKmerCount[lcaTaxId] += 1; // Count distinct k-mers
                ++localDistinctKmerCount; // Count distinct k-mers
                taxIds.clear();
                kmer = nextKmer;
            }

            #pragma omp critical
            {
                std::cout << "Processed " << i << "th split out of " << numOfDiffIdxSplits_use << std::endl;
            }   
        }

        #pragma omp critical
        {
            for (const auto &pair : localTaxon2Count) {
                taxon2Count[pair.first] += pair.second;
            }
            for (const auto &pair : localTaxon2UniqKmerCount) {
                taxon2UniqKmerCount[pair.first] += pair.second;
            }
            distinctKmerCount += localDistinctKmerCount;
        }
        delete deltaIdxReaders;
    }

    std::vector<std::string> rankOrders = {
        "root", "superkingdom", "phylum", "class", "order", "family", "genus", "species", "subspecies", "no rank"
    };

    std::unordered_map<std::string, size_t> rank2Count;
    for (const auto &pair : taxon2Count) {
        TaxID taxId = pair.first;
        size_t count = pair.second;
        const TaxonNode * node = taxonomy->taxonNode(taxId);
        if (taxId == 1) {
            rank2Count["root"] += count;        
            continue;
        }
        if (node) {
            rank2Count[taxonomy->getString(node->rankIdx)] += count;
        } else {
            std::cerr << "TaxID " << taxId << " not found in the taxonomy." << std::endl;
        }
    }

    std::cout << "Common k-mers count by rank:" << std::endl;
    for (const auto &rank : rankOrders) {
        if (rank2Count.find(rank) != rank2Count.end()) {
            std::cout << rank << ": " << rank2Count[rank] << std::endl;
        } 
    }
    // for (const auto &pair : rank2Count) {
    //     std::cout << pair.first << ": " << pair.second << std::endl;
    // }

    std::cout << "Total distinct common k-mers count: " << distinctKmerCount << std::endl;
    std::unordered_map<std::string, size_t> rank2UniqKmerCount;
    for (const auto &pair : taxon2UniqKmerCount) {
        TaxID taxId = pair.first;
        size_t count = pair.second;
        const TaxonNode * node = taxonomy->taxonNode(taxId);
        if (taxId == 1) {
            rank2UniqKmerCount["root"] += count;        
            continue;
        }
        if (node) {
            rank2UniqKmerCount[taxonomy->getString(node->rankIdx)] += count;
        } else {
            std::cerr << "TaxID " << taxId << " not found in the taxonomy." << std::endl;
        }
    }
    std::cout << "Distinct common k-mers count by rank:" << std::endl;
    for (const auto &rank : rankOrders) {
        if (rank2UniqKmerCount.find(rank) != rank2UniqKmerCount.end()) {
            std::cout << rank << ": " << rank2UniqKmerCount[rank] << std::endl;
        } 
    }
    // for (const auto &pair : rank2UniqKmerCount) {
    //     std::cout << pair.first << ": " << pair.second << std::endl;
    // }

    delete taxonomy;

    return 0;
}
#endif //METABULI_COUNT_COMMON_KMERS_H