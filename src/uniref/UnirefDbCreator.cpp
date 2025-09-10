#include "UnirefDbCreator.h"


UnirefDbCreator::UnirefDbCreator(
    const LocalParameters &par,
    const std::string & dbDir) 
    : dbDir(dbDir), par(par)
{
 
}


int UnirefDbCreator::getLCA(
    const std::vector<int> & ids,
    const std::unordered_map<int, std::pair<int, int>> & uniref100to90and50
){
    if (ids.size() == 1) {
        return ids[0];
    }

    auto it = uniref100to90and50.find(ids[0]);
    if (it == uniref100to90and50.end()) {
        std::cerr << "Error: ID " << ids[0] << " not found in UniRef100 to UniRef90/50 mapping." << std::endl;
        exit(EXIT_FAILURE);
    }
    const int first100 = ids[0];
    const int first90 = it->second.first;
    const int first50 = it->second.second;

    bool allSame100 = true;
    bool allSame90 = true;
    bool allSame50 = true;

    for (size_t i = 1; i < ids.size(); ++i) {
        // we can stop early as the LCA must be the root.
        if (!allSame50) {
            break;
        }

        auto current_it = uniref100to90and50.find(ids[i]);
        if (current_it == uniref100to90and50.end()) {
            std::cerr << "Error: ID " << ids[i] << " not found in UniRef100 to UniRef90/50 mapping." << std::endl;
            exit(EXIT_FAILURE);
        }

        // Check each level and set the flag to false on the first mismatch.
        if (allSame100 && ids[i] != first100) {
            allSame100 = false;
        }
        if (allSame90 && current_it->second.first != first90) {
            allSame90 = false;
        }
        if (allSame50 && current_it->second.second != first50) {
            allSame50 = false;
        }
    }

    // Return the most specific (lowest number) level where all IDs matched.
    if (allSame100) return first100;
    if (allSame90) return first90;
    if (allSame50) return first50;
    
    return 1; // If nothing matched, the LCA is the root.
}



// void UnirefDbCreator::createUnirefDb() {
//     std::cout << "Creating UniRef database..." << std::endl;
    
//     Buffer<Kmer> kmerBuffer(Buffer<Kmer>::calculateBufferSize(par.ramUsage, par.threads, sizeof(Kmer) + sizeof(size_t)));
//     Buffer<size_t> uniqKmerIdx(kmerBuffer.bufferSize);
//     vector<pair<size_t, size_t>> uniqKmerIdxRanges;
    
//     KSeqWrapper * kseq = KSeqFactory(uniref100fasta.c_str());
//     std::unordered_map<string, uint32_t> uniref100toTaxId;
//     uint32_t idOffset = 0;
    
    
//     bool complete = false;
//     SeqEntry savedSeq;
//     size_t processedSeqCnt = 0;
//     while (!complete) {
//         // Extract k-mers
//         time_t start = time(nullptr);
//         cout << "K-mer extraction    : " << flush;
//         bool moreData = kmerExtractor->extractUnirefKmers(kseq, kmerBuffer, uniref100toTaxId, processedSeqCnt, savedSeq);
//         complete = !moreData;
//         cout << double(time(nullptr) - start) << " s" << endl;
//         cout << "Processed sequences : " << processedSeqCnt << endl;

//         // Sort the k-mers
//         start = time(nullptr);
//         cout << "Sort k-mers         : " << flush;
//         SORT_PARALLEL(kmerBuffer.buffer, kmerBuffer.buffer + kmerBuffer.startIndexOfReserve, Kmer::compareKmer);
//         cout << double(time(nullptr) - start) << " s" << endl;

//         // Filter k-mers
//         start = time(nullptr);
//         size_t selectedKmerCnt = 0;
//         uniqKmerIdxRanges.clear();
//         filterKmers<FilterMode::LCA>(kmerBuffer, uniqKmerIdx.buffer, selectedKmerCnt, uniqKmerIdxRanges);
//         cout << "Reduce k-mers       : " << time(nullptr) - start << " s" << endl; 
//         cout << "Selected k-mers     : " << selectedKmerCnt << endl;

//         // Write k-mers
//         start = time(nullptr);
//         if (complete && numOfFlush == 0 && !isUpdating) {
//             writeTargetFilesAndSplits(kmerBuffer, uniqKmerIdx.buffer, selectedKmerCnt, uniqKmerIdxRanges, true);
//         } else {
//             writeTargetFiles(kmerBuffer, uniqKmerIdx.buffer, uniqKmerIdxRanges);
//         }
//         cout << "Write k-mers        : " << time(nullptr) - start << " s" << endl;

//         // Write accession to index mapping
//         FILE * accIndexFile = fopen((dbDir + "/accession2index").c_str(), "a");
//         for (const auto & entry : accession2index) {
//             fprintf(accIndexFile, "%s\t%u\n", entry.first.c_str(), entry.second);
//         }
//         fclose(accIndexFile);

//         if (!complete) {
//             accession2index.clear();
//             kmerBuffer.init();
//             uniqKmerIdx.init();
//         }
//         cout << "--------" << endl;
//     }
    

//     if (numOfFlush == 1) {
//         cout << "Index creation completed." << endl;
//         return;
//     }
//     cout << "Merge reference DB files ... " << endl;

//     // for (int i = 0; i < 10; i++) {
//     //     addFilesToMerge(dbDir + "/" + to_string(i) + "_diffIdx",
//     //                     dbDir + "/" + to_string(i) + "_info");
//     // }

//     indexCreator->printFilesToMerge();
//     indexCreator->setMergedFileNames(
//         dbDir + "/diffIdx",  
//         dbDir + "/info", 
//         dbDir + "/split");
    
//     indexCreator->mergeTargetFiles<FilterMode::LCA>();


//     std::cout << "UniRef database creation completed." << std::endl;
// }