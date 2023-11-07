#include "Taxonomer.h"
#include "Match.h"
#include "NcbiTaxonomy.h"
#include <cstdint>
#include <sys/types.h>
#include <unordered_map>


Taxonomer::Taxonomer(const LocalParameters &par, NcbiTaxonomy *taxonomy) : taxonomy(taxonomy) {
    // Parameters
    auto mask = new uint32_t[par.spaceMask.length()];
    for(size_t i = 0, j = 0; i < par.spaceMask.length(); i++){
        mask[i] = par.spaceMask[i] - 48;
        spaceNum += (mask[i] == 0);
        if(mask[i] == 1){
            unmaskedPos[j] = (int) i;
            j++;
        }
    }
    delete[] mask;
    maxGap = par.maxGap;
    minCoveredPos = par.minCoveredPos;
    accessionLevel = par.accessionLevel;
    minSSMatch = par.minSSMatch;
    minConsCnt = par.minConsCnt;
    minConsCntEuk = par.minConsCntEuk;
    eukaryotaTaxId = par.eukaryotaTaxId;
}

Taxonomer::~Taxonomer() {

}

void Taxonomer::assignTaxonomy(const Match *matchList,
                               size_t numOfMatches,
                               std::vector<Query> &queryList,
                               const LocalParameters &par) {
    time_t beforeAnalyze = time(nullptr);
    cout << "Analyzing matches ..." << endl;

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
#pragma omp for schedule(dynamic, 1)
        for (size_t i = 0; i < blockIdx; ++i) {
            chooseBestTaxon2(matchBlocks[i].id,
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
    cout << "Time spent for analyzing: " << double(time(nullptr) - beforeAnalyze) << endl;

}

void Taxonomer::chooseBestTaxon2(uint32_t currentQuery,
                                 size_t offset,
                                 size_t end,
                                 const Match *matchList,
                                 vector<Query> & queryList,
                                 const LocalParameters &par) {

//    if (true) {
//        cout << "# " << currentQuery << " " << queryList[currentQuery].name << endl;
//        for (size_t i = offset; i < end + 1; i++) {
//            cout << matchList[i].targetId << " " << matchList[i].qInfo.frame << " " << matchList[i].qInfo.pos << " " << int(matchList[i].hamming) <<  " "  << int(matchList[i].redundancy) << endl;
//        }
//    }

    // Get the best species for current query
    vector<Match> speciesMatches;
    speciesMatches.reserve(end - offset + 1);
    TaxonScore speciesScore(0, 0, 0, 0, 0);
    if (par.seqMode == 2) {
        speciesScore = getBestSpeciesMatches(speciesMatches, matchList, end, offset,
                                             queryList[currentQuery].queryLength,
                                             queryList[currentQuery].queryLength2);
    } else {
        speciesScore = getBestSpeciesMatches(speciesMatches, matchList, end, offset,
                                             queryList[currentQuery].queryLength);
    }

//    if (true) {
//        cout << "# " << currentQuery << " " << queryList[currentQuery].name << " filtered\n";
//        for (size_t i = 0; i < genusMatches.size(); i++) {
//           cout << genusMatches[i].targetId << " " << genusMatches[i].qInfo.frame << " " << genusMatches[i].qInfo.pos << " " << int(genusMatches[i].hamming) <<  " "  << int(genusMatches[i].redundancy) << endl;
//         }
//        cout << "Genus score: " << genusScore.score << "\n";
//    }

    // If there is no proper species for current query, it is un-classified.
    if (speciesScore.score == 0 || speciesScore.coverage < par.minCoverage || speciesScore.score < par.minScore) {
        queryList[currentQuery].isClassified = false;
        queryList[currentQuery].classification = 0;
        queryList[currentQuery].score = speciesScore.score;
        queryList[currentQuery].coverage = speciesScore.coverage;
        queryList[currentQuery].hammingDist = speciesScore.hammingDist;
        queryList[currentQuery].newSpecies = false;
        return;
    }

    // If there are two or more good species level candidates, find the LCA.
    if (speciesScore.LCA) {
        // cout << "LCA" << endl;
        // vector<TaxID> genusList;
        // genusList.reserve(speciesMatches.size());
        // for (auto & genusMatch : speciesMatches) {
        //     genusList.push_back(genusMatch.genusId);
        // }
        // selectedTaxon = taxonomy->LCA(genusList)->taxId;
        queryList[currentQuery].isClassified = true;
        queryList[currentQuery].classification = speciesScore.taxId;
        queryList[currentQuery].score = speciesScore.score;
        queryList[currentQuery].coverage = speciesScore.coverage;
        queryList[currentQuery].hammingDist = speciesScore.hammingDist;
        // for (auto & spMatch : speciesMatches) {
        //     queryList[currentQuery].taxCnt[spMatch.targetId]++;
        // }
        return;
    }

    // If score is not enough, classify to the parent of the selected species
    if (speciesScore.score < par.minSpScore) {
        queryList[currentQuery].isClassified = true;
        queryList[currentQuery].classification = taxonomy->taxonNode(
                taxonomy->getTaxIdAtRank(speciesScore.taxId, "species"))->parentTaxId;
        queryList[currentQuery].score = speciesScore.score;
        queryList[currentQuery].coverage = speciesScore.coverage;
        queryList[currentQuery].hammingDist = speciesScore.hammingDist;
        for (auto & spMatch : speciesMatches) {
            queryList[currentQuery].taxCnt[spMatch.targetId]++;
        }
        return;
    }

    // Sort matches by the position of the query sequence
//    sort(genusMatches.begin() + speciesMatchRange[selectedSpecies].first,
//         genusMatches.begin() + speciesMatchRange[selectedSpecies].second,
//         [](const Match & a, const Match & b) {
//        if (a.qInfo.position / 3 == b.qInfo.position / 3)
//            return a.hamming < b.hamming;
//        else
//            return a.qInfo.position / 3 < b.qInfo.position / 3;
//    });

    // sort(speciesMatches.begin(), speciesMatches.end(),
    //      [](const Match & a, const Match & b) { return a.qInfo.pos < b.qInfo.pos; });

    cout << "7 " << currentQuery << endl;

    TaxID result = lowerRankClassification(speciesMatches, speciesScore.taxId);

    // Record matches of selected species
    for (auto & spMatch : speciesMatches) {
            queryList[currentQuery].taxCnt[spMatch.targetId]++;
    }

    // Store classification results
    queryList[currentQuery].isClassified = true;
    queryList[currentQuery].classification = result;
    queryList[currentQuery].score = speciesScore.score;
    queryList[currentQuery].coverage = speciesScore.coverage;
    queryList[currentQuery].hammingDist = speciesScore.hammingDist;
    queryList[currentQuery].newSpecies = false;
    cout << "8" << currentQuery << endl;
//    if (par.printLog) {
//        cout << "# " << currentQuery << endl;
//        for (size_t i = 0; i < genusMatches.size(); i++) {
//            cout << i << " " << genusMatches[i].qInfo.pos << " " <<
//            genusMatches[i].targetId << " " << int(genusMatches[i].hamming) << endl;
//        }
//        cout << "Score: " << speciesScore.score << "  " << selectedSpecies << " "
//             << taxonomy->getString(taxonomy->taxonNode(selectedSpecies)->rankIdx)
//
//             << endl;
//    }
}

// void Taxonomer::chooseBestTaxon(uint32_t currentQuery,
//                                  size_t offset,
//                                  size_t end,
//                                  const Match *matchList,
//                                  vector<Query> & queryList,
//                                  const LocalParameters &par) {
//     TaxID selectedTaxon;

// //    if (true) {
// //        cout << "# " << currentQuery << " " << queryList[currentQuery].name << endl;
// //        for (size_t i = offset; i < end + 1; i++) {
// //            cout << matchList[i].targetId << " " << matchList[i].qInfo.frame << " " << matchList[i].qInfo.pos << " " << int(matchList[i].hamming) <<  " "  << int(matchList[i].redundancy) << endl;
// //        }
// //    }

//     // Get the best genus for current query
//     vector<Match> genusMatches;
//     genusMatches.reserve(end - offset + 1);
//     TaxonScore genusScore(0, 0, 0, 0);
//     if (par.seqMode == 2) {
//         genusScore = getBestGenusMatches(genusMatches, matchList, end, offset,
//                                              queryList[currentQuery].queryLength,
//                                              queryList[currentQuery].queryLength2);
//     } else {
//         genusScore = getBestGenusMatches(genusMatches, matchList, end, offset,
//                                              queryList[currentQuery].queryLength);
//     }

// //    if (true) {
// //        cout << "# " << currentQuery << " " << queryList[currentQuery].name << " filtered\n";
// //        for (size_t i = 0; i < genusMatches.size(); i++) {
// //           cout << genusMatches[i].targetId << " " << genusMatches[i].qInfo.frame << " " << genusMatches[i].qInfo.pos << " " << int(genusMatches[i].hamming) <<  " "  << int(genusMatches[i].redundancy) << endl;
// //         }
// //        cout << "Genus score: " << genusScore.score << "\n";
// //    }

//     // If there is no proper genus for current query, it is un-classified.
//     if (genusScore.score == 0 || genusScore.coverage < par.minCoverage || genusScore.score < par.minScore) {
//         queryList[currentQuery].isClassified = false;
//         queryList[currentQuery].classification = 0;
//         queryList[currentQuery].score = genusScore.score;
//         queryList[currentQuery].coverage = genusScore.coverage;
//         queryList[currentQuery].hammingDist = genusScore.hammingDist;
//         queryList[currentQuery].newSpecies = false;
//         return;
//     }

//     // If there are two or more good genus level candidates, find the LCA.
//     if (genusScore.taxId == 0) {
//         vector<TaxID> genusList;
//         genusList.reserve(genusMatches.size());
//         for (auto & genusMatch : genusMatches) {
//             genusList.push_back(genusMatch.genusId);
//         }
//         selectedTaxon = taxonomy->LCA(genusList)->taxId;
//         queryList[currentQuery].isClassified = true;
//         queryList[currentQuery].classification = selectedTaxon;
//         queryList[currentQuery].score = genusScore.score;
//         queryList[currentQuery].coverage = genusScore.coverage;
//         queryList[currentQuery].hammingDist = genusScore.hammingDist;
//         for (auto & genusMatch : genusMatches) {
//             queryList[currentQuery].taxCnt[genusMatch.targetId]++;
//         }
//         return;
//     }

//     // Choose the species with the highest coverage.
//     TaxID selectedSpecies;
//     TaxonScore speciesScore;
//     vector<TaxID> species;
//     unordered_map<TaxID, pair<int, int>> speciesMatchRange;
//     if (par.seqMode == 2) {
//         speciesScore = chooseSpecies(genusMatches,
//                                      queryList[currentQuery].queryLength,
//                                      queryList[currentQuery].queryLength2,
//                                      species,
//                                      speciesMatchRange);
//     } else {
//         speciesScore = chooseSpecies(genusMatches,
//                                      queryList[currentQuery].queryLength,
//                                      species,
//                                      speciesMatchRange);
//     }


//     // Classify to LCA if more than one species are selected
//     if (species.size() > 1) {
//         queryList[currentQuery].isClassified = true;
//         queryList[currentQuery].classification = taxonomy->LCA(species)->taxId;
//         queryList[currentQuery].score = genusScore.score;
//         queryList[currentQuery].coverage = genusScore.coverage;
//         queryList[currentQuery].hammingDist = genusScore.hammingDist;
//         for (auto & genusMatch : genusMatches) {
//             queryList[currentQuery].taxCnt[genusMatch.targetId]++;
//         }
//         return;
//     }

//     // If score is not enough, classify to the parent of the selected species
//     if (speciesScore.score < par.minSpScore) {
//         queryList[currentQuery].isClassified = true;
//         queryList[currentQuery].classification = taxonomy->taxonNode(
//                 taxonomy->getTaxIdAtRank(species[0], "species"))->parentTaxId;
//         queryList[currentQuery].score = genusScore.score;
//         queryList[currentQuery].coverage = genusScore.coverage;
//         queryList[currentQuery].hammingDist = genusScore.hammingDist;
//         for (auto & genusMatch : genusMatches) {
//             if(genusMatch.speciesId == species[0]){
//                 queryList[currentQuery].taxCnt[genusMatch.targetId]++;
//             }
//         }
//         return;
//     }

//     // Sort matches by the position of the query sequence
//     selectedSpecies = species[0];
// //    sort(genusMatches.begin() + speciesMatchRange[selectedSpecies].first,
// //         genusMatches.begin() + speciesMatchRange[selectedSpecies].second,
// //         [](const Match & a, const Match & b) {
// //        if (a.qInfo.position / 3 == b.qInfo.position / 3)
// //            return a.hamming < b.hamming;
// //        else
// //            return a.qInfo.position / 3 < b.qInfo.position / 3;
// //    });

//     vector<Match>::const_iterator first = genusMatches.begin() + speciesMatchRange[selectedSpecies].first;
//     vector<Match>::const_iterator last = genusMatches.begin() + speciesMatchRange[selectedSpecies].second;
//     vector<Match> speciesMatches(first, last);


//     sort(speciesMatches.begin(), speciesMatches.end(),
//          [](const Match & a, const Match & b) { return a.qInfo.pos < b.qInfo.pos; });
         
//     TaxID result = lowerRankClassification(speciesMatches, selectedSpecies);

//     // Record matches of selected species
//     for (size_t i = speciesMatchRange[selectedSpecies].first; i < speciesMatchRange[selectedSpecies].second; i++) {
//         queryList[currentQuery].taxCnt[genusMatches[i].targetId]++;
//     }

//     // Store classification results
//     queryList[currentQuery].isClassified = true;
//     queryList[currentQuery].classification = result;
//     queryList[currentQuery].score = speciesScore.score;
//     queryList[currentQuery].coverage = speciesScore.coverage;
//     queryList[currentQuery].hammingDist = speciesScore.hammingDist;
//     queryList[currentQuery].newSpecies = false;
// //    if (par.printLog) {
// //        cout << "# " << currentQuery << endl;
// //        for (size_t i = 0; i < genusMatches.size(); i++) {
// //            cout << i << " " << genusMatches[i].qInfo.pos << " " <<
// //            genusMatches[i].targetId << " " << int(genusMatches[i].hamming) << endl;
// //        }
// //        cout << "Score: " << speciesScore.score << "  " << selectedSpecies << " "
// //             << taxonomy->getString(taxonomy->taxonNode(selectedSpecies)->rankIdx)
// //
// //             << endl;
// //    }
// }

TaxID Taxonomer::lowerRankClassification(vector<Match> &matches, TaxID spTaxId) {
    unordered_map<TaxID, unsigned int> taxCnt;
    size_t matchNum = matches.size();
    // cout << spTaxId << endl;
    

    for (size_t i = 0; i < matchNum; i++) {
        // cout << matches[i].targetId << endl;
        taxCnt[matches[i].targetId] ++;
        // size_t currQuotient = matches[i].qInfo.pos / 3;
        // uint8_t minHamming = 0; //matches[i].hamming;
        // Match * minHammingMatch = & matches[i];
        // TaxID minHammingTaxId = minHammingMatch->targetId;
        // while ((i < matchNum) && (currQuotient == matches[i].qInfo.pos / 3)) {
        //     if (matches[i].hamming < minHamming) {
        //         minHamming = matches[i].hamming;
        //         minHammingMatch = & matches[i];
        //         minHammingTaxId = minHammingMatch->targetId;
        //     } else if (matches[i].hamming == minHamming) {
        //         minHammingTaxId = taxonomy->LCA(minHammingTaxId, matches[i].targetId);
        //         minHammingMatch->redundancy = true;
        //         matches[i].redundancy = true;
        //     }
        //     i++;
        // }
        // taxCnt[minHammingTaxId]++;       
    }

    // int i = matchRange.second - 1;
    // while ( i >= matchRange.first ) {
    //     size_t currQuotient = matches[i].qInfo.pos / 3;
    //     uint8_t minHamming = matches[i].hamming;
    //     Match * minHammingMatch = & matches[i];
    //     TaxID minHammingTaxId = minHammingMatch->targetId;
    //     i --;
    //     while ( (i >= matchRange.first) && (currQuotient == matches[i].qInfo.pos / 3) ) {
    //         if (matches[i].hamming < minHamming) {
    //             minHamming = matches[i].hamming;
    //             minHammingMatch = & matches[i];
    //             minHammingTaxId = minHammingMatch->targetId;
    //         } else if (matches[i].hamming == minHamming) {
    //             minHammingTaxId = taxonomy->LCA(minHammingTaxId, matches[i].targetId);
    //             minHammingMatch->redundancy = true;
    //             matches[i].redundancy = true;
    //         }
    //         i--;
    //     }
    //     taxCnt[minHammingTaxId]++;
    // }

    unordered_map<TaxID, TaxonCounts> cladeCnt;
    // cout << "8" << endl;
    getSpeciesCladeCounts(taxCnt, cladeCnt, spTaxId);
    // // print cladeCnt
    // for (auto it = cladeCnt.begin(); it != cladeCnt.end(); it++) {
    //     cout << it->first << " " << it->second.taxCount << " " << it->second.cladeCount << endl;
    // }
    // cout << "9" << endl;
    if (accessionLevel == 2) { // Don't do accession-level classification
        // Remove leaf nodes
        // cout << "10" << endl;
        for (auto it = cladeCnt.begin(); it != cladeCnt.end(); it++) {
            TaxonNode const * taxon = taxonomy->taxonNode(it->first);
            if (strcmp(taxonomy->getString(taxon->rankIdx), "") == 0) {
                // Remove current node from its parent's children list
                cladeCnt[taxon->parentTaxId].children.erase(find(cladeCnt[taxon->parentTaxId].children.begin(),
                                                                 cladeCnt[taxon->parentTaxId].children.end(),
                                                                 it->first));
            } 
        }
        
        return BFS(cladeCnt, spTaxId);
    } else {
        // cout << "10-2" << endl;
        return BFS(cladeCnt, spTaxId);
    }
}

void Taxonomer::getSpeciesCladeCounts(const unordered_map<TaxID, unsigned int> &taxCnt,
                                       unordered_map<TaxID, TaxonCounts> & cladeCount,
                                       TaxID speciesTaxID) {
    for (auto it = taxCnt.begin(); it != taxCnt.end(); ++it) {
//        cladeCount[it->first].taxCount = it->second;
//        cladeCount[it->first].cladeCount += it->second;
        TaxonNode const * taxon = taxonomy->taxonNode(it->first);
        cladeCount[taxon->taxId].taxCount = it->second;
        cladeCount[taxon->taxId].cladeCount += it->second;
        while (taxon->taxId != speciesTaxID) {
            if (find(cladeCount[taxon->parentTaxId].children.begin(),
                     cladeCount[taxon->parentTaxId].children.end(),
                     taxon->taxId) == cladeCount[taxon->parentTaxId].children.end()) {
                cladeCount[taxon->parentTaxId].children.push_back(taxon->taxId);
            }
            cladeCount[taxon->parentTaxId].cladeCount += it->second;
            taxon = taxonomy->taxonNode(taxon->parentTaxId);
        }
    }
}

TaxID Taxonomer::BFS(const unordered_map<TaxID, TaxonCounts> & cladeCnt, TaxID root) {
    if (cladeCnt.at(root).children.empty()) { // root is a leaf
        return root;
    }
    unsigned int maxCnt = minSSMatch;
    unsigned int currentCnt;
    vector<TaxID> bestChildren;
    for (auto it = cladeCnt.at(root).children.begin(); it != cladeCnt.at(root).children.end(); it++) {
        currentCnt = cladeCnt.at(*it).cladeCount;
        if (currentCnt > maxCnt) {
            bestChildren.clear();
            bestChildren.push_back(*it);
            maxCnt = currentCnt;
        } else if (currentCnt == maxCnt) {
            bestChildren.push_back(*it);
        }
    }
    if (bestChildren.size() == 1) {
        return BFS(cladeCnt, bestChildren[0]);
    } else {
        return root;
    }
}

TaxonScore Taxonomer::getBestSpeciesMatches(vector<Match> &speciesMatches,
                                            const Match *matchList,
                                            size_t end,
                                            size_t offset,
                                            int queryLength) {
    TaxID currentSpecies;
    vector<const Match *> filteredMatches;
    vector<vector<const Match *>> matchesForEachSpecies;
    vector<TaxonScore> speciesScores;
    TaxonScore bestScore;
    size_t i = offset;
    uint8_t curFrame;
    vector<const Match *> curFrameMatches;
    vector<MatchPath> matchPaths;

     while (i  < end + 1) {
        currentSpecies = matchList[i].speciesId;
        // For current species
        while ((i < end + 1) && currentSpecies == matchList[i].speciesId) {
            curFrame = matchList[i].qInfo.frame;
            curFrameMatches.clear();
            // For current frame
            while ((i < end + 1) && currentSpecies == matchList[i].speciesId && curFrame == matchList[i].qInfo.frame) {
                curFrameMatches.push_back(&matchList[i]);
                i ++;
            }
            if (curFrameMatches.size() > 1) {
                remainConsecutiveMatches(curFrameMatches, matchPaths, currentSpecies);
            }
        }
        // Construct a match combination using filtered matches of current species
        // so that it can best cover the query, and score the combination
        if (!filteredMatches.empty()) {
            matchesForEachSpecies.push_back(filteredMatches);
            speciesScores.push_back(scoreTaxon(filteredMatches, currentSpecies, queryLength));
        }
        filteredMatches.clear();
    }
    
    // If there are no meaningful species
    if (speciesScores.empty()) {
        bestScore.score = 0;
        return bestScore;
    }

    TaxonScore maxScore = *max_element(speciesScores.begin(), speciesScores.end(),
                                       [](const TaxonScore & a, const TaxonScore & b) { return a.score < b.score; });

    vector<size_t> maxIdx;
    for (size_t g = 0; g < speciesScores.size(); g++) {
        if (speciesScores[g].score == maxScore.score) {
            maxIdx.push_back(g);
        }
    }
    bestScore = maxScore;

    for (unsigned long g : maxIdx) {
        for (const Match * m : matchesForEachSpecies[g]) {
            speciesMatches.push_back(*m);
        }
    }

    // More than one species
    if (maxIdx.size() > 1) {
        bestScore.taxId = 0;
    }

    return bestScore;                    
}

TaxonScore Taxonomer::getBestSpeciesMatches(vector<Match> &speciesMatches,
                                            const Match *matchList,
                                            size_t end,
                                            size_t offset,
                                            int readLength1,
                                            int readLength2) {
    TaxID currentSpecies;
    vector<const Match *> filteredMatches;
    vector<vector<const Match *>> matchesForEachSpecies;
    vector<TaxonScore> speciesScores;
    TaxonScore bestScore;
    size_t i = offset;
    uint8_t curFrame;
    vector<const Match *> curFrameMatches;
    vector<MatchPath> matchPaths;
    unordered_map<TaxID, float> species2score;
    unordered_map<TaxID, vector<MatchPath>> species2matchPaths;
    float bestSpScore = 0;

     while (i  < end + 1) {
        currentSpecies = matchList[i].speciesId;
        // For current species
        while ((i < end + 1) && currentSpecies == matchList[i].speciesId) {
            curFrame = matchList[i].qInfo.frame;
            curFrameMatches.clear();
            // For current frame
            while ((i < end + 1) && currentSpecies == matchList[i].speciesId && curFrame == matchList[i].qInfo.frame) {
                curFrameMatches.push_back(&matchList[i]);
                i ++;
            }
            if (curFrameMatches.size() > 1) {
                // cout << "1" << endl;
                remainConsecutiveMatches(curFrameMatches, matchPaths, currentSpecies);
            }
        }
        // Combine MatchPaths
        // so that it can best cover the query, and score the combination
        if (!matchPaths.empty()) {
            // Initialize species2matchPaths
            species2matchPaths[currentSpecies].emplace_back(0, 0, 0, 0);
            // cout << "2" << endl;
            cout << currentSpecies << endl;
            float score = combineMatchPaths(matchPaths, species2matchPaths[currentSpecies], readLength1 + readLength2);
            cout << endl;
            species2score[currentSpecies] = score;
            if (score > bestSpScore) {
                bestSpScore = score;
            }
        }
        matchPaths.clear();
    }
    
    // If there are no meaningful species
    if (species2score.empty()) {
        bestScore.score = 0;
        return bestScore;
    }
    // cout << "4" << endl;
    vector<TaxID> maxSpecies;
    for (auto & spScore : species2score) {
        cout << spScore.first << " " << spScore.second << endl;
        if (spScore.second == bestSpScore) {
            maxSpecies.push_back(spScore.first);
        }
    }
    // cout << "5" << endl;
    // More than one species --> LCA
    if (maxSpecies.size() > 1) {
        bestScore.LCA = true;
        bestScore.taxId = taxonomy->LCA(maxSpecies)->taxId;
        for (auto & sp : maxSpecies) {
            bestScore.score += species2score[sp];
        }
        bestScore.score /= maxSpecies.size();
        return bestScore;
    }
    
    // One species
    bestScore.taxId = maxSpecies[0];
    bestScore.score = species2score[maxSpecies[0]];
    float coveredLength = 0.f;
    int hammingDist = 0;
    for (auto & matchPath : species2matchPaths[maxSpecies[0]]) {
        // cout << "here" << endl;
        coveredLength += matchPath.end - matchPath.start + 1;
        hammingDist += matchPath.hammingDist;
        for (auto match : matchPath.matches) {
            // cout << match->targetId << endl;
            // match->printMatch();
            speciesMatches.push_back(*match);
            // speciesMatches.back().printMatch();
        }
    }
    bestScore.coverage = coveredLength / (readLength1 + readLength2);
    bestScore.hammingDist = hammingDist;

// cout << "6" << endl;
    return bestScore;                    
}

float Taxonomer::combineMatchPaths(vector<MatchPath> & matchPaths,
                                   vector<MatchPath> & combinedMatchPaths,
                                   int readLength) {
    combinedMatchPaths.clear();
    // Sort matchPaths by the their score
    sort(matchPaths.begin(), matchPaths.end(),
         [](const MatchPath & a, const MatchPath & b) { return a.score > b.score;});

    // Combine matchPaths
    // 1. Add the matchPath with the highest score to combinedMatchPaths
    // 2. Add the matchPath with the highest score that is not overlapped with the matchPath in combinedMatchPaths
    // 3. Repeat 2 until no matchPath can be added
    for (size_t i = 0; i < matchPaths.size(); i++) {
        cout << matchPaths[i].start << " " << matchPaths[i].end << " " << matchPaths[i].score << " " << matchPaths[i].matches.back()->targetId << " " << matchPaths[i].matches.back()->qInfo.frame <<endl;
        if (combinedMatchPaths.empty()) {
            combinedMatchPaths.push_back(matchPaths[i]);
            combinedMatchPaths.back().matches = matchPaths[i].matches;
            // cout << matchPaths[i].start << " " << matchPaths[i].end << " " << matchPaths[i].score << endl;
            // for (auto & match : matchPaths[i].matches) {
            //     match->printMatch();
            // }
        } else {
            bool isOverlapped = false;
            for (size_t j = 0; j < combinedMatchPaths.size(); j++) {
                if (!isMatchPathNotOverlapped(matchPaths[i], combinedMatchPaths[j])) {
                    isOverlapped = true;
                    break;
                } else {
                    // cout << matchPaths[i].start << " " << matchPaths[i].end << endl;
                    // cout << combinedMatchPaths[j].start << " " << combinedMatchPaths[j].end << endl << endl;;
                }
            }
            if (!isOverlapped) {
                combinedMatchPaths.push_back(matchPaths[i]);
                combinedMatchPaths.back().matches = matchPaths[i].matches;
                // cout << matchPaths[i].start << " " << matchPaths[i].end << " " << matchPaths[i].score << endl;
                // for (auto & match : matchPaths[i].matches) {
                //     match->printMatch();
                // }
            }
        }
    }
    // cout << endl;



    // Calculate the score of combinedMatchPaths
    float score = 0;
    for (auto & matchPath : combinedMatchPaths) {
        score += matchPath.score;
    }
    return score / readLength;
}

bool Taxonomer::isMatchPathNotOverlapped(const MatchPath & matchPath1,
                                         const MatchPath & matchPath2) {
    return (matchPath1.end < matchPath2.start) || (matchPath2.end < matchPath1.start);                                       
}

\

//     if (matchPath1.start > matchPath2.start) {
//         return isMatchPathOverlapped(matchPath2, matchPath1, readLength);
//     }
//     if (matchPath1.end < matchPath2.start) {
//         return false;
//     }
//     if (matchPath1.endPos >= matchPath2.startPos) {
//         if (matchPath1.endPos <= matchPath2.endPos) {
//             return true;
//         } else {
//             if (matchPath1.startPos + readLength - 1 >= matchPath2.startPos) {
//                 return true;
//             } else {
//                 return false;
//             }
//         }
//     }
//     return false;
// }


// TaxonScore Taxonomer::getBestGenusMatches(vector<Match> &genusMatches, const Match *matchList, size_t end,
//                                            size_t offset, int readLength1, int readLength2) {
//     TaxID currentGenus;
//     TaxID currentSpecies;

//     vector<const Match *> filteredMatches;
//     vector<vector<const Match *>> matchesForEachGenus;
//     vector<TaxonScore> genusScores;
//     TaxonScore bestScore;
//     size_t i = offset;
//     uint8_t curFrame;
//     vector<const Match *> curFrameMatches;
//     while (i  < end + 1) {
// //        currentGenus = taxId2genusId[matchList[i].targetId];
//         currentGenus = matchList[i].genusId;
//         // For current genus
//         while ((i < end + 1) && currentGenus == matchList[i].genusId) {
// //            currentSpecies = taxId2speciesId[matchList[i].targetId];
//             currentSpecies = matchList[i].speciesId;
// //            if (par.printLog) {
// //                cout << currentGenus << " " << currentSpecies << endl;
// //            }
//             // For current species
//             while ((i < end + 1) && currentSpecies == matchList[i].speciesId) {
//                 curFrame = matchList[i].qInfo.frame;
//                 curFrameMatches.clear();

//                 // For current frame
//                 while ((i < end + 1) && currentSpecies == matchList[i].speciesId
//                        && curFrame == matchList[i].qInfo.frame) {
//                     curFrameMatches.push_back(&matchList[i]);
//                     i ++;
//                 }
//                 if (curFrameMatches.size() > 1) {
//                     remainConsecutiveMatches(curFrameMatches, filteredMatches, currentGenus);
//                 }
//             }
//         }

//         // Construct a match combination using filtered matches of current genus
//         // so that it can best cover the query, and score the combination
//         if (!filteredMatches.empty()) {
//             matchesForEachGenus.push_back(filteredMatches);
//             genusScores.push_back(scoreTaxon(filteredMatches, currentGenus, readLength1, readLength2));
//         }
//         filteredMatches.clear();
//     }

//     // If there are no meaningful genus
//     if (genusScores.empty()) {
//         bestScore.score = 0;
//         return bestScore;
//     }

//     TaxonScore maxScore = *max_element(genusScores.begin(), genusScores.end(),
//                                        [](const TaxonScore & a, const TaxonScore & b) { return a.score < b.score; });

//     vector<size_t> maxIdx;
//     for (size_t g = 0; g < genusScores.size(); g++) {
//         if (genusScores[g].score > maxScore.score * 0.95f) {
//             maxIdx.push_back(g);
//         }
//     }
//     bestScore = maxScore;

//     for (unsigned long g : maxIdx) {
//         for (const Match * m : matchesForEachGenus[g]) {
//             genusMatches.push_back(*m);
//         }
//     }



//     // More than one genus
//     if (maxIdx.size() > 1) {
//         bestScore.taxId = 0;
//         return bestScore;
//     }

//     return bestScore;

//     //Three cases
//     //1. one genus
//     //2. more than one genus
//     //4. no genus
// }

void Taxonomer::remainConsecutiveMatches(const vector<const Match *> & curFrameMatches,
                                         vector<MatchPath> & matchPaths,
                                         TaxID genusId) {
    size_t i = 0;
    size_t end = curFrameMatches.size();
    vector<pair<const Match *, size_t>> curPosMatches; // <match, index>
    vector<pair<const Match *, size_t>> nextPosMatches;
    map<size_t, vector<size_t>> linkedMatches; // <index, linked indexes>

    size_t currPos = curFrameMatches[0]->qInfo.pos;
    while ( i < end && curFrameMatches[i]->qInfo.pos == currPos) {
        curPosMatches.emplace_back(curFrameMatches[i], i);
        i++;
    }
    while (i < end) {
        uint32_t nextPos = curFrameMatches[i]->qInfo.pos;
        while (i < end  && nextPos == curFrameMatches[i]->qInfo.pos) {
            nextPosMatches.emplace_back(curFrameMatches[i], i);
            ++ i;
        }
        // Check if current position and next position are consecutive
        if (currPos + 3 == nextPos) {
            // Compare curPosMatches and nextPosMatches
            for (auto &curPosMatch: curPosMatches) {
                for (auto &nextPosMatch: nextPosMatches) {
                    if (isConsecutive(curPosMatch.first, nextPosMatch.first)) {
                        linkedMatches[curPosMatch.second].push_back(nextPosMatch.second);
                    }
                }
            }

        }
        // Update curPosMatches and nextPosMatches
        curPosMatches = nextPosMatches;
        nextPosMatches.clear();
        currPos = nextPos;
    }
    // Print linkedMatches
//    if (par.printLog) {
//        cout << "linkedMatches: " << endl;
//        for (const auto &entry: linkedMatches) {
//            cout << entry.first << ": ";
//            for (auto &idx: entry.second) {
//                cout << idx << " ";
//            }
//            cout << endl;
//        }
//    }

    // Iterate linkedMatches to get filteredMatches 
    //(ignore matches not enoughly consecutive)
    size_t MIN_DEPTH = minConsCnt - 1;
    if (taxonomy->IsAncestor(eukaryotaTaxId, genusId)) {
        MIN_DEPTH = minConsCntEuk - 1;
    }
    unordered_set<size_t> used;
    unordered_map<size_t, depthScore> idx2depthScore;
    // unordered_map<size_t, size_t> edges;
    unordered_map<const Match *, const Match *> edges;

    for (const auto& entry : linkedMatches) {
        if (!used.count(entry.first)) {
            used.insert(entry.first);
            depthScore bestPath{};
            size_t bestNextIdx = 0;
            float curScore = curFrameMatches[entry.first]->getScore();
            // cout << curFrameMatches[entry.first]
            for (size_t j = 0; j < entry.second.size(); j++) {
                used.insert(entry.second[j]);
                depthScore curPath = DFS(curFrameMatches,
                                         entry.second[j],
                                         linkedMatches,
                                        1,
                                         MIN_DEPTH, 
                                         used, 
                                         idx2depthScore,
                                         edges, 
                                         curScore, 0);
                if (curPath.score > bestPath.score && curPath.depth > MIN_DEPTH) {
                    bestNextIdx = entry.second[j];
                    bestPath = curPath;
                }
            }
            // Store the best path
            if (bestPath.depth > MIN_DEPTH) {
                // cout << entry.first << endl;
                // curFrameMatches[entry.first]->printMatch();
                matchPaths.emplace_back(curFrameMatches[entry.first]->qInfo.pos, // start coordinate on query
                                        curFrameMatches[entry.first]->qInfo.pos + bestPath.depth * 3 + 20, // end coordinate on query
                                        bestPath.score, bestPath.hammingDist);
                const Match * curMatch = curFrameMatches[entry.first];
                edges[curMatch] = curFrameMatches[bestNextIdx];
                matchPaths.back().matches.push_back(curMatch);
                // curMatch = edges[curMatch];                        
                // edges2[curFrameMatches[entry.first]] = curFrameMatches[bestNextIdx];
                // Retrieve the best path
                // cout << bestPath.depth << endl;
                while (edges.find(curMatch) != edges.end()) {
                    // cout << curMatch << " ";
                    matchPaths.back().matches.push_back(edges[curMatch]);
                    curMatch = edges[curMatch];
                }
                // cout << endl;
            }
        }
    }

//    if (par.printLog) {
//        cout << "filteredMatchIdx: ";
//        for (auto &idx: filteredMatchIdx) {
//            cout << idx << " ";
//        }
//        cout << endl;
//    }
}

// size_t Taxonomer::DFS(size_t curMatchIdx, const map<size_t, vector<size_t>> & linkedMatches,
//                        vector<size_t>& filteredMatches, size_t depth, size_t MIN_DEPTH, unordered_set<size_t>& used,
//                        unordered_map<size_t, size_t> & idx2depth) {
//     depth++;
//     size_t maxDepth = 0;
//     size_t returnDepth = 0;
//     if (linkedMatches.find(curMatchIdx) == linkedMatches.end()) { 
//         // reached a leaf node
//         idx2depth[curMatchIdx] = depth;
//         if (depth > MIN_DEPTH) {
//             filteredMatches.push_back(curMatchIdx);
//         }
//         return depth;
//     } else { // not a leaf node
//         for (auto &nextMatchIdx: linkedMatches.at(curMatchIdx)) {
//             used.insert(nextMatchIdx);
//             if (idx2depth.find(nextMatchIdx) != idx2depth.end()) {
//                 returnDepth = idx2depth[nextMatchIdx];
//                 maxDepth = max(maxDepth, returnDepth);
//                 continue;
//             }
//             returnDepth = DFS(nextMatchIdx, linkedMatches, filteredMatches, depth, MIN_DEPTH, used, idx2depth);
//             maxDepth = max(maxDepth, returnDepth);
//         }
//         if (maxDepth > MIN_DEPTH) {
//             filteredMatches.push_back(curMatchIdx);
//             idx2depth[curMatchIdx] = maxDepth;
//         }
//     }
//     return maxDepth;
// }

// return: end
depthScore Taxonomer::DFS(const vector<const Match *> &matches,
                          size_t curMatchIdx,
                          const map<size_t, vector<size_t>> &linkedMatches,
                          size_t depth, size_t MIN_DEPTH,
                          unordered_set<size_t> &used,
                          unordered_map<size_t, depthScore> &idx2depthScore,
                          unordered_map<const Match *, const Match *> & edges, float score, int hammingDist) {
    depth++;
    depthScore bestDepthScore = depthScore(0, 0, 0);
    depthScore returnDepthScore;
    if (linkedMatches.find(curMatchIdx) == linkedMatches.end()) { // reached a leaf node
        uint8_t lastEndHamming = (matches[curMatchIdx]->rightEndHamming >> 14);
        if (lastEndHamming == 0) {
            score += 3.0f;
        } else {
            score += 2.0f - 0.5f * lastEndHamming;
        }
        idx2depthScore[curMatchIdx] = depthScore(depth, score, hammingDist + lastEndHamming);
        return depthScore(depth, score, hammingDist + lastEndHamming);
    } else { // not a leaf node
        uint8_t lastEndHamming = (matches[curMatchIdx]->rightEndHamming >> 14);
        if (lastEndHamming == 0) {
            score += 3.0f;
        } else {
            score += 2.0f - 0.5f * lastEndHamming;
        }
        for (auto &nextMatchIdx: linkedMatches.at(curMatchIdx)) {
            used.insert(nextMatchIdx);

            // Reuse the depth score of nextMatchIdx if it has been calculated
            if (idx2depthScore.find(nextMatchIdx) != idx2depthScore.end()) {
                returnDepthScore = idx2depthScore[nextMatchIdx];
                if (returnDepthScore.score > bestDepthScore.score
                    && returnDepthScore.depth > MIN_DEPTH) {
                    bestDepthScore = returnDepthScore;
                    edges[matches[curMatchIdx]] = matches[nextMatchIdx];
                }   
                continue;
            }
            returnDepthScore = DFS(matches, nextMatchIdx, linkedMatches, depth, MIN_DEPTH, used, idx2depthScore, edges, score, hammingDist + lastEndHamming);
            if (returnDepthScore.score > bestDepthScore.score
                && returnDepthScore.depth > MIN_DEPTH) {
                bestDepthScore = returnDepthScore;
                edges[matches[curMatchIdx]] = matches[nextMatchIdx];
            } 
        }    
        if (bestDepthScore.depth > MIN_DEPTH) {
            idx2depthScore[curMatchIdx] = bestDepthScore;
        }
    }
    return bestDepthScore;
}

// TaxonScore Taxonomer::getBestGenusMatches_spaced(vector<Match> &genusMatches, const Match *matchList, size_t end,
//                                                   size_t offset, int readLength1, int readLength2) {
//     TaxID currentGenus;
//     TaxID currentSpecies;

//     vector<const Match *> tempMatchContainer;
//     vector<const Match *> filteredMatches;
//     vector<vector<const Match *>> matchesForEachGenus;
//     vector<bool> conservedWithinGenus;
//     vector<TaxonScore> genusScores;
//     TaxonScore bestScore;
//     size_t i = offset;
//     bool lastIn;
//     while (i + 1 < end + 1) {
//         currentGenus = matchList[i].genusId;
//         // For current genus
//         while ((i + 1 < end + 1) && currentGenus == matchList[i].genusId) {
// //            currentSpecies = taxId2speciesId[matchList[i].targetId];
//             currentSpecies = matchList[i].speciesId;
//             // For current species
//             // Filter un-consecutive matches (probably random matches)
//             lastIn = false;
//             int distance = 0;
//             int diffPosCntOfCurrRange = 1;
//             int dnaDist = 0;

//             // For the same species
//             while ((i + 1 < end + 1) && currentSpecies == matchList[i + 1].speciesId) {
//                 distance = matchList[i+1].qInfo.pos / 3 - matchList[i].qInfo.pos / 3;
//                 dnaDist = matchList[i+1].qInfo.pos - matchList[i].qInfo.pos;
//                 if (distance == 0) { // At the same position
//                     tempMatchContainer.push_back(matchList + i);
//                 } else if (dnaDist < (8 + spaceNum + maxGap) * 3) { // Overlapping
//                     lastIn = true;
//                     tempMatchContainer.push_back(matchList + i);
//                     diffPosCntOfCurrRange ++;
//                 } else { // Not consecutive --> End range
//                     if (lastIn){
//                         tempMatchContainer.push_back(matchList + i);
//                         if (diffPosCntOfCurrRange >= minCoveredPos) {
//                             filteredMatches.insert(filteredMatches.end(), tempMatchContainer.begin(),
//                                                    tempMatchContainer.end());
//                         }
//                     }
//                     lastIn = false;
//                     // Initialize range info
//                     tempMatchContainer.clear();
//                     diffPosCntOfCurrRange = 1;
//                 }
//                 i++;
//             }

//             // Met next species
//             if (lastIn) {
//                 tempMatchContainer.push_back(matchList + i);
//                 if (diffPosCntOfCurrRange >= minCoveredPos) {
//                     filteredMatches.insert(filteredMatches.end(), tempMatchContainer.begin(),
//                                            tempMatchContainer.end());
//                 }
//             }
//             tempMatchContainer.clear();
//             i++;
//         }

//         // Construct a match combination using filtered matches of current genus
//         // so that it can best cover the query, and score the combination
//         if (!filteredMatches.empty()) {
//             genusScores.push_back(scoreTaxon(filteredMatches, readLength1, readLength2));
//         }
//         filteredMatches.clear();
//     }

//     // If there are no meaningful genus
//     if (genusScores.empty()) {
//         bestScore.score = 0;
//         return bestScore;
//     }

//     TaxonScore maxScore = *max_element(genusScores.begin(), genusScores.end(),
//                                        [](const TaxonScore & a, const TaxonScore & b) { return a.score < b.score; });

//     vector<size_t> maxIdx;
//     for (size_t g = 0; g < genusScores.size(); g++) {
//         if (genusScores[g].score > maxScore.score * 0.95f) {
//             maxIdx.push_back(g);
//         }
//     }
//     bestScore = maxScore;

//     for (unsigned long g : maxIdx) {
//         for (const Match * m : matchesForEachGenus[g]) {
//             genusMatches.push_back(*m);
//         }
//     }

//     // More than one genus
//     if (maxIdx.size() > 1) {
//         bestScore.taxId = 0;
//         return bestScore;
//     }
//     return bestScore;

//     //Three cases
//     //1. one genus
//     //2. more than one genus
//     //4. no genus
// }

// TaxonScore Taxonomer::getBestGenusMatches(vector<Match> &genusMatches, const Match *matchList, size_t end,
//                                            size_t offset, int queryLength) {
//     TaxID currentGenus;
//     TaxID currentSpecies;

//     vector<const Match *> filteredMatches;
//     vector<vector<const Match *>> matchesForEachGenus;
//     vector<TaxonScore> genusScores;
//     TaxonScore bestScore;
//     size_t i = offset;
//     uint8_t curFrame;
//     vector<const Match *> curFrameMatches;
//     while (i  < end + 1) {
//         currentGenus = matchList[i].genusId;
//         // For current genus
//         while ((i < end + 1) && currentGenus == matchList[i].genusId) {
//             currentSpecies = matchList[i].speciesId;

//             // For current species
//             while ((i < end + 1) && currentSpecies == matchList[i].speciesId) {
//                 curFrame = matchList[i].qInfo.frame;
//                 curFrameMatches.clear();

//                 // For current frame
//                 while ((i < end + 1) && currentSpecies == matchList[i].speciesId
//                        && curFrame == matchList[i].qInfo.frame) {
//                     curFrameMatches.push_back(&matchList[i]);
//                     i ++;
//                 }
//                 if (curFrameMatches.size() > 1) {
//                     remainConsecutiveMatches(curFrameMatches, filteredMatches, currentGenus);
//                 }
//             }
//         }

//         // Construct a match combination using filtered matches of current genus
//         // so that it can best cover the query, and score the combination

//         if (!filteredMatches.empty()) {
//             matchesForEachGenus.push_back(filteredMatches);
//             genusScores.push_back(scoreTaxon(filteredMatches, currentGenus, queryLength));
//         }
//         filteredMatches.clear();
//     }

//     // If there are no meaningful genus
//     if (genusScores.empty()) {
//         bestScore.score = 0;
//         return bestScore;
//     }

//     TaxonScore maxScore = *max_element(genusScores.begin(), genusScores.end(),
//                                        [](const TaxonScore & a, const TaxonScore & b) { return a.score < b.score; });

//     vector<size_t> maxIdx;
//     for (size_t g = 0; g < genusScores.size(); g++) {
//         if (genusScores[g].score > maxScore.score * 0.95f) {
//             maxIdx.push_back(g);
//         }
//     }

//     bestScore = maxScore;

//     for (unsigned long g : maxIdx) {
//         for (const Match * m : matchesForEachGenus[g]) {
//             genusMatches.push_back(*m);
//         }
//     }

//     // More than one genus
//     if (maxIdx.size() > 1) {
//         bestScore.taxId = 0;
//         return bestScore;
//     }
//     return bestScore;

//     //Three cases
//     //1. one genus
//     //2. more than one genus
//     //4. no genus
// }

// TaxonScore Taxonomer::getBestGenusMatches_spaced(vector<Match> &genusMatches, const Match *matchList, size_t end,
//                                                  size_t offset, int readLength) {
//     TaxID currentGenus;
//     TaxID currentSpecies;

//     vector<const Match *> tempMatchContainer;
//     vector<const Match *> filteredMatches;
//     vector<vector<Match>> matchesForEachGenus;
//     vector<bool> conservedWithinGenus;
//     vector<TaxonScore> genusScores;
//     TaxonScore bestScore;
//     size_t i = offset;
//     bool lastIn;
//     size_t speciesMatchCnt;
//     while (i + 1 < end + 1) {
//         currentGenus = matchList[i].genusId;
//         // For current genus
//         while ((i + 1 < end + 1) && currentGenus == matchList[i].genusId) {
//             currentSpecies = matchList[i].speciesId;
//             // For current species
//             // Filter un-consecutive matches (probably random matches)
//             lastIn = false;
//             int distance = 0;
//             int diffPosCntOfCurrRange = 1;
//             int dnaDist = 0;

//             // For the same species
//             while ((i + 1 < end + 1) && currentSpecies == matchList[i + 1].speciesId) {
//                 distance = matchList[i + 1].qInfo.pos / 3 - matchList[i].qInfo.pos / 3;
//                 dnaDist = matchList[i + 1].qInfo.pos - matchList[i].qInfo.pos;
//                 if (distance == 0) { // At the same position
//                     tempMatchContainer.push_back(matchList + i);
//                 } else if (dnaDist < (8 + spaceNum + maxGap) * 3) { // Overlapping
//                     lastIn = true;
//                     tempMatchContainer.push_back(matchList + i);
//                     diffPosCntOfCurrRange++;
//                 } else { // Not consecutive --> End range
//                     if (lastIn) {
//                         tempMatchContainer.push_back(matchList + i);
//                         if (diffPosCntOfCurrRange >= minCoveredPos) {
//                             filteredMatches.insert(filteredMatches.end(), tempMatchContainer.begin(),
//                                                    tempMatchContainer.end());
//                         }
//                     }
//                     lastIn = false;
//                     // Initialize range info
//                     tempMatchContainer.clear();
//                     diffPosCntOfCurrRange = 1;
//                 }
//                 i++;
//             }

//             // Met next species
//             if (lastIn) {
//                 tempMatchContainer.push_back(matchList + i);
//                 if (diffPosCntOfCurrRange >= minCoveredPos) {
//                     filteredMatches.insert(filteredMatches.end(), tempMatchContainer.begin(),
//                                            tempMatchContainer.end());
//                 }
//             }
//             tempMatchContainer.clear();
//             i++;
//         }

//         // Construct a match combination using filtered matches of current genus
//         // so that it can best cover the query, and score the combination
//         if (!filteredMatches.empty()) {
//             genusScores.push_back(scoreTaxon(filteredMatches, readLength));
//         }
//         filteredMatches.clear();
//     }

//     // If there are no meaningful genus
//     if (genusScores.empty()) {
//         bestScore.score = 0;
//         return bestScore;
//     }

//     TaxonScore maxScore = *max_element(genusScores.begin(), genusScores.end(),
//                                        [](const TaxonScore &a, const TaxonScore &b) { return a.score < b.score; });

//     vector<size_t> maxIdx;
//     for (size_t g = 0; g < genusScores.size(); g++) {
//         if (genusScores[g].score > maxScore.score * 0.95f) {
//             maxIdx.push_back(g);
//         }
//     }
//     bestScore = maxScore;

//     for (unsigned long g: maxIdx) {
//         genusMatches.insert(genusMatches.end(),
//                             matchesForEachGenus[g].begin(),
//                             matchesForEachGenus[g].end());
//     }

//     // More than one genus
//     if (maxIdx.size() > 1) {
//         bestScore.taxId = 0;
//         return bestScore;
//     }
//     return bestScore;

//     //Three cases
//     //1. one genus
//     //2. more than one genus
//     //4. no genus
// }

TaxonScore Taxonomer::scoreTaxon(vector<const Match *> &filteredMatches,
                                 TaxID taxId,
                                 int queryLength) {
    // Calculate Hamming distance & covered length
    int coveredPosCnt = 0;
    uint16_t currHammings;
    int aminoAcidNum = (int) queryLength / 3;
    int currPos;
    size_t matchNum = filteredMatches.size();
    size_t f = 0;

    // Get the smallest hamming distance at each position of query
    auto *hammingsAtEachPos = new signed char[aminoAcidNum + 1];
    memset(hammingsAtEachPos, 24, (aminoAcidNum + 1));
    while (f < matchNum) {
        currPos = filteredMatches[f]->qInfo.pos / 3;
        currHammings = filteredMatches[f]->rightEndHamming;
        if (GET_2_BITS(currHammings) < hammingsAtEachPos[currPos + unmaskedPos[0]])
            hammingsAtEachPos[currPos + unmaskedPos[0]] = GET_2_BITS(currHammings);
        if (GET_2_BITS(currHammings >> 2) < hammingsAtEachPos[currPos + unmaskedPos[1]])
            hammingsAtEachPos[currPos + unmaskedPos[1]] = GET_2_BITS(currHammings >> 2);
        if (GET_2_BITS(currHammings >> 4) < hammingsAtEachPos[currPos + unmaskedPos[2]])
            hammingsAtEachPos[currPos + unmaskedPos[2]] = GET_2_BITS(currHammings >> 4);
        if (GET_2_BITS(currHammings >> 6) < hammingsAtEachPos[currPos + unmaskedPos[3]])
            hammingsAtEachPos[currPos + unmaskedPos[3]] = GET_2_BITS(currHammings >> 6);
        if (GET_2_BITS(currHammings >> 8) < hammingsAtEachPos[currPos + unmaskedPos[4]])
            hammingsAtEachPos[currPos + unmaskedPos[4]] = GET_2_BITS(currHammings >> 8);
        if (GET_2_BITS(currHammings >> 10) < hammingsAtEachPos[currPos + unmaskedPos[5]])
            hammingsAtEachPos[currPos + unmaskedPos[5]] = GET_2_BITS(currHammings >> 10);
        if (GET_2_BITS(currHammings >> 12) < hammingsAtEachPos[currPos + unmaskedPos[6]])
            hammingsAtEachPos[currPos + unmaskedPos[6]] = GET_2_BITS(currHammings >> 12);
        if (GET_2_BITS(currHammings >> 14) < hammingsAtEachPos[currPos + unmaskedPos[7]])
            hammingsAtEachPos[currPos + unmaskedPos[7]] = GET_2_BITS(currHammings >> 14);
        f++;
    }

    // Sum up hamming distances and count the number of position covered by the matches.
    float hammingSum = 0;
    for (int h = 0; h < aminoAcidNum; h++) {
        if (hammingsAtEachPos[h] == 0) { // Add 0 for 0 hamming dist.
            coveredPosCnt++;
        } else if (hammingsAtEachPos[h] != 24) { // Add 1.5, 2, 2.5 for 1, 2, 3 hamming dist. respectively
            hammingSum += 1.0f + (0.5f * hammingsAtEachPos[h]);
            coveredPosCnt++;
        }
    }
    delete[] hammingsAtEachPos;

    // Score current genus
    int coveredLength = coveredPosCnt * 3;
    if (coveredLength > queryLength) coveredLength = queryLength;
    float score = ((float) coveredLength - hammingSum) / (float) queryLength;
    float coverage = (float) (coveredLength) / (float) (queryLength);

    return {taxId, score, coverage, (int) hammingSum, 0};
}

TaxonScore Taxonomer::scoreTaxon(vector<const Match *> &filteredMatches,
                                 TaxID taxId,       
                                 int readLength1,
                                 int readLength2) {

    // Calculate Hamming distance & covered length
    int aminoAcidNum_total = ((int) readLength1 / 3) + ((int) readLength2 / 3);
    int aminoAcidNum_read1 = ((int) readLength1 / 3);
    int currPos;
    size_t matchNum = filteredMatches.size();
    size_t f = 0;

    // Get the smallest hamming distance at each position of query
    auto *hammingsAtEachPos = new signed char[aminoAcidNum_total + 3];
    memset(hammingsAtEachPos, 24, (aminoAcidNum_total + 3));
    while (f < matchNum) {
        uint8_t minHammingDist = 24;
        uint16_t currHammings = 0;
        currPos = (int) filteredMatches[f]->qInfo.pos / 3;
        // Find the closest match at current position
        while ((f < matchNum) && currPos == (int) filteredMatches[f]->qInfo.pos / 3) {
            if (filteredMatches[f]->hamming < minHammingDist) {
                minHammingDist = filteredMatches[f]->hamming;
                currHammings = filteredMatches[f]->rightEndHamming;
            }
            f++;
        }
        // Update hamming distance at each position
        if (GET_2_BITS(currHammings) < hammingsAtEachPos[currPos + unmaskedPos[0]])
            hammingsAtEachPos[currPos + unmaskedPos[0]] = GET_2_BITS(currHammings);
        if (GET_2_BITS(currHammings >> 2) < hammingsAtEachPos[currPos + unmaskedPos[1]])
            hammingsAtEachPos[currPos + unmaskedPos[1]] = GET_2_BITS(currHammings >> 2);
        if (GET_2_BITS(currHammings >> 4) < hammingsAtEachPos[currPos + unmaskedPos[2]])
            hammingsAtEachPos[currPos + unmaskedPos[2]] = GET_2_BITS(currHammings >> 4);
        if (GET_2_BITS(currHammings >> 6) < hammingsAtEachPos[currPos + unmaskedPos[3]])
            hammingsAtEachPos[currPos + unmaskedPos[3]] = GET_2_BITS(currHammings >> 6);
        if (GET_2_BITS(currHammings >> 8) < hammingsAtEachPos[currPos + unmaskedPos[4]])
            hammingsAtEachPos[currPos + unmaskedPos[4]] = GET_2_BITS(currHammings >> 8);
        if (GET_2_BITS(currHammings >> 10) < hammingsAtEachPos[currPos + unmaskedPos[5]])
            hammingsAtEachPos[currPos + unmaskedPos[5]] = GET_2_BITS(currHammings >> 10);
        if (GET_2_BITS(currHammings >> 12) < hammingsAtEachPos[currPos + unmaskedPos[6]])
            hammingsAtEachPos[currPos + unmaskedPos[6]] = GET_2_BITS(currHammings >> 12);
        if (GET_2_BITS(currHammings >> 14) < hammingsAtEachPos[currPos + unmaskedPos[7]])
            hammingsAtEachPos[currPos + unmaskedPos[7]] = GET_2_BITS(currHammings >> 14);
        f++;
    }

    // Sum up hamming distances and count the number of position covered by the matches.
    float hammingSum = 0;
    int coveredPosCnt_read1 = 0;
    int coveredPosCnt_read2 = 0;
    for (int h = 0; h < aminoAcidNum_total; h++) {
        // Read 1
        if (h < aminoAcidNum_read1) {
            if (hammingsAtEachPos[h] == 0) { // Add 0 for 0 hamming dist.
                coveredPosCnt_read1++;
            } else if (hammingsAtEachPos[h] != 24) { // Add 1.5, 2, 2.5 for 1, 2, 3 hamming dist. respectively
                hammingSum += 1.0f + (0.5f * (float) hammingsAtEachPos[h]);
                coveredPosCnt_read1++;
            }
        }
            // Read 2
        else {
            if (hammingsAtEachPos[h] == 0) { // Add 0 for 0 hamming dist.
                coveredPosCnt_read2++;
            } else if (hammingsAtEachPos[h] != 24) { // Add 1.5, 2, 2.5 for 1, 2, 3 hamming dist. respectively
                hammingSum += 1.0f + (0.5f * (float) hammingsAtEachPos[h]);
                coveredPosCnt_read2++;
            }
        }
    }
    delete[] hammingsAtEachPos;

    // Score current genus
    int coveredLength_read1 = coveredPosCnt_read1 * 3;
    int coveredLength_read2 = coveredPosCnt_read2 * 3;
    if (coveredLength_read1 > readLength1) coveredLength_read1 = readLength1;
    if (coveredLength_read2 > readLength2) coveredLength_read2 = readLength2;
    float score =
            ((float) (coveredLength_read1 + coveredLength_read2) - hammingSum) / (float) (readLength1 + readLength2);
    float coverage = (float) (coveredLength_read1 + coveredLength_read2) / (float) (readLength1 + readLength2);

//    matchesForEachGenus.push_back(move(filteredMatches));
    return {taxId, score, coverage, (int) hammingSum, 0};
}

TaxonScore Taxonomer::chooseSpecies(const vector<Match> &matches,
                                     int queryLength,
                                     vector<TaxID> &species,
                                     unordered_map<TaxID, pair<int, int>> & speciesMatchRange) {
    // Score each species
    std::unordered_map<TaxID, TaxonScore> speciesScores;
    size_t i = 0;
    TaxID currentSpeices;
    size_t numOfMatch = matches.size();
    size_t speciesBegin, speciesEnd;
    while (i < numOfMatch) {
        currentSpeices = matches[i].speciesId;
        speciesBegin = i;
        while ((i < numOfMatch) && currentSpeices == matches[i].speciesId) {
            i++;
        }
        speciesEnd = i;
        speciesScores[currentSpeices] = scoreSpecies(matches, speciesBegin, speciesEnd, queryLength);
        speciesMatchRange[currentSpeices] = {(int) speciesBegin, (int) speciesEnd};
        speciesScores[currentSpeices].taxId = currentSpeices;
    }

    // Get the best species
    TaxonScore bestScore;
    for (auto & sp : speciesScores) {
        if (sp.second.score > bestScore.score) {
            species.clear();
            species.push_back(sp.first);
            bestScore = sp.second;
        } else if (sp.second.coverage == bestScore.coverage) {
            species.push_back(sp.first);
        }
    }
    return bestScore;
}

TaxonScore Taxonomer::chooseSpecies(const vector<Match> &matches,
                                     int read1Length,
                                     int read2Length,
                                     vector<TaxID> &species,
                                     unordered_map<TaxID, pair<int, int>> & speciesMatchRange) {
    // Score each species
    std::unordered_map<TaxID, TaxonScore> speciesScores;

    size_t i = 0;
    TaxID currentSpeices;
    size_t numOfMatch = matches.size();
    size_t speciesBegin, speciesEnd;
    while (i < numOfMatch) {
        currentSpeices = matches[i].speciesId;
        speciesBegin = i;
        while ((i < numOfMatch) && currentSpeices == matches[i].speciesId) {
            i++;
        }
        speciesEnd = i;
        speciesScores[currentSpeices] = scoreSpecies(matches, speciesBegin, speciesEnd, read1Length, read2Length);
        speciesMatchRange[currentSpeices] = {(int) speciesBegin, (int) speciesEnd};
        speciesScores[currentSpeices].taxId = currentSpeices;
    }

    // Get the best species
    TaxonScore bestScore;
    for (auto & sp : speciesScores) {
        if (sp.second.score > bestScore.score) {
            species.clear();
            species.push_back(sp.first);
            bestScore = sp.second;
        } else if (sp.second.coverage == bestScore.coverage) {
            species.push_back(sp.first);
        }
    }
    return bestScore;
}

TaxonScore Taxonomer::scoreSpecies(const vector<Match> &matches,
                                    size_t begin,
                                    size_t end,
                                    int queryLength) {

    // Get the largest hamming distance at each position of query
    int aminoAcidNum = queryLength / 3;
    auto *hammingsAtEachPos = new signed char[aminoAcidNum + 1];
    memset(hammingsAtEachPos, -1, (aminoAcidNum + 1));
    int currPos;
    size_t walker = begin;
    uint16_t currHammings;
    while (walker < end) {
        currPos = matches[walker].qInfo.pos / 3;
        currHammings = matches[walker].rightEndHamming;
        if (GET_2_BITS(currHammings) > hammingsAtEachPos[currPos + unmaskedPos[0]])
            hammingsAtEachPos[currPos + unmaskedPos[0]] = GET_2_BITS(currHammings);
        if (GET_2_BITS(currHammings >> 2) > hammingsAtEachPos[currPos + unmaskedPos[1]])
            hammingsAtEachPos[currPos + unmaskedPos[1]] = GET_2_BITS(currHammings >> 2);
        if (GET_2_BITS(currHammings >> 4) > hammingsAtEachPos[currPos + unmaskedPos[2]])
            hammingsAtEachPos[currPos + unmaskedPos[2]] = GET_2_BITS(currHammings >> 4);
        if (GET_2_BITS(currHammings >> 6) > hammingsAtEachPos[currPos + unmaskedPos[3]])
            hammingsAtEachPos[currPos + unmaskedPos[3]] = GET_2_BITS(currHammings >> 6);
        if (GET_2_BITS(currHammings >> 8) > hammingsAtEachPos[currPos + unmaskedPos[4]])
            hammingsAtEachPos[currPos + unmaskedPos[4]] = GET_2_BITS(currHammings >> 8);
        if (GET_2_BITS(currHammings >> 10) > hammingsAtEachPos[currPos + unmaskedPos[5]])
            hammingsAtEachPos[currPos + unmaskedPos[5]] = GET_2_BITS(currHammings >> 10);
        if (GET_2_BITS(currHammings >> 12) > hammingsAtEachPos[currPos + unmaskedPos[6]])
            hammingsAtEachPos[currPos + unmaskedPos[6]] = GET_2_BITS(currHammings >> 12);
        if (GET_2_BITS(currHammings >> 14) > hammingsAtEachPos[currPos + unmaskedPos[7]])
            hammingsAtEachPos[currPos + unmaskedPos[7]] = GET_2_BITS(currHammings >> 14);
        walker++;
    }

    // Sum up hamming distances and count the number of position covered by the matches.
    float hammingSum = 0;
    int hammingDist = 0;
    int coveredPosCnt = 0;
    for (int h = 0; h < aminoAcidNum; h++) {
        if (hammingsAtEachPos[h] == 0) { // Add 0 for 0 hamming dist.
            coveredPosCnt++;
        } else if (hammingsAtEachPos[h] != -1) { // Add 1.5, 2, 2.5 for 1, 2, 3 hamming dist. respectively
            hammingSum += 1.0f + (0.5f * (float) hammingsAtEachPos[h]);
            hammingDist += hammingsAtEachPos[h];
            coveredPosCnt++;
        }
    }
    delete[] hammingsAtEachPos;
    // Score
    int coveredLength = coveredPosCnt * 3;
    if (coveredLength >= queryLength) coveredLength = queryLength;

    float score = ((float)coveredLength - hammingSum) / (float) queryLength;
    float coverage = (float) coveredLength / (float) (queryLength);

    return {0, score, coverage, hammingDist, 0};
}

TaxonScore Taxonomer::scoreSpecies(const vector<Match> &matches,
                                    size_t begin,
                                    size_t end,
                                    int queryLength,
                                    int queryLength2) {

    // Get the largest hamming distance at each position of query
    int aminoAcidNum_total = queryLength / 3 + queryLength2 / 3;
    int aminoAcidNum_read1 = queryLength / 3;
    auto *hammingsAtEachPos = new signed char[aminoAcidNum_total + 3];
    memset(hammingsAtEachPos, -1, (aminoAcidNum_total + 3));

    int currPos;
    size_t walker = begin;
    uint16_t currHammings;

    while (walker < end) {
        currPos = matches[walker].qInfo.pos / 3;
        currHammings = matches[walker].rightEndHamming;
        if (GET_2_BITS(currHammings) > hammingsAtEachPos[currPos + unmaskedPos[0]])
            hammingsAtEachPos[currPos + unmaskedPos[0]] = GET_2_BITS(currHammings);
        if (GET_2_BITS(currHammings >> 2) > hammingsAtEachPos[currPos + unmaskedPos[1]])
            hammingsAtEachPos[currPos + unmaskedPos[1]] = GET_2_BITS(currHammings >> 2);
        if (GET_2_BITS(currHammings >> 4) > hammingsAtEachPos[currPos + unmaskedPos[2]])
            hammingsAtEachPos[currPos + unmaskedPos[2]] = GET_2_BITS(currHammings >> 4);
        if (GET_2_BITS(currHammings >> 6) > hammingsAtEachPos[currPos + unmaskedPos[3]])
            hammingsAtEachPos[currPos + unmaskedPos[3]] = GET_2_BITS(currHammings >> 6);
        if (GET_2_BITS(currHammings >> 8) > hammingsAtEachPos[currPos + unmaskedPos[4]])
            hammingsAtEachPos[currPos + unmaskedPos[4]] = GET_2_BITS(currHammings >> 8);
        if (GET_2_BITS(currHammings >> 10) > hammingsAtEachPos[currPos + unmaskedPos[5]])
            hammingsAtEachPos[currPos + unmaskedPos[5]] = GET_2_BITS(currHammings >> 10);
        if (GET_2_BITS(currHammings >> 12) > hammingsAtEachPos[currPos + unmaskedPos[6]])
            hammingsAtEachPos[currPos + unmaskedPos[6]] = GET_2_BITS(currHammings >> 12);
        if (GET_2_BITS(currHammings >> 14) > hammingsAtEachPos[currPos + unmaskedPos[7]])
            hammingsAtEachPos[currPos + unmaskedPos[7]] = GET_2_BITS(currHammings >> 14);
        walker++;
    }

    // Sum up hamming distances and count the number of position covered by the matches.
    float hammingSum = 0;
    int hammingDist = 0;
    int coveredPosCnt_read1 = 0;
    int coveredPosCnt_read2 = 0;
    for (int h = 0; h < aminoAcidNum_total; h++) {
        // Read 1
        if (h < aminoAcidNum_read1) {
            if (hammingsAtEachPos[h] == 0) { // Add 0 for 0 hamming dist.
                coveredPosCnt_read1++;
            } else if (hammingsAtEachPos[h] != -1) { // Add 1.5, 2, 2.5 for 1, 2, 3 hamming dist. respectively
                hammingSum += 1.0f + (0.5f * (float) hammingsAtEachPos[h]);
                hammingDist += hammingsAtEachPos[h];
                coveredPosCnt_read1++;
            }
        }
            // Read 2
        else {
            if (hammingsAtEachPos[h] == 0) { // Add 0 for 0 hamming dist.
                coveredPosCnt_read2++;
            } else if (hammingsAtEachPos[h] != -1) { // Add 1.5, 2, 2.5 for 1, 2, 3 hamming dist. respectively
                hammingSum += 1.0f + (0.5f * (float) hammingsAtEachPos[h]);
                hammingDist += hammingsAtEachPos[h];
                coveredPosCnt_read2++;
            }
        }
    }
    delete[] hammingsAtEachPos;

    // Score
    int coveredLength_read1 = coveredPosCnt_read1 * 3;
    int coveredLength_read2 = coveredPosCnt_read2 * 3;
    if (coveredLength_read1 >= queryLength) coveredLength_read1 = queryLength;
    if (coveredLength_read2 >= queryLength2) coveredLength_read2 = queryLength2;

    float score = ((float) (coveredLength_read1 + coveredLength_read2) - hammingSum) / (float) (queryLength + queryLength2);
    float coverage = (float) (coveredLength_read1 + coveredLength_read2) / (float) (queryLength + queryLength2);

    return {0, score, coverage, hammingDist, 0};
}

bool Taxonomer::isConsecutive(const Match * match1, const Match * match2) {
    // match1 87654321 -> 08765432
    // match2 98765432 -> 08765432
    return (match1->rightEndHamming >> 2) == (match2->rightEndHamming & 0x3FFF);
}