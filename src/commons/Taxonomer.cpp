#include "Taxonomer.h"
#include "BitManipulateMacros.h"
#include "Match.h"
#include "NcbiTaxonomy.h"
#include "common.h"
#include "printBinary.h"
#include <cstdint>
#include <sys/types.h>
#include <unordered_map>


Taxonomer::Taxonomer(const LocalParameters &par, TaxonomyWrapper *taxonomy) : taxonomy(taxonomy), par(par) {
    // Parameters
    string spaceMask = "11111111";
    auto mask = new uint32_t[spaceMask.length()];
    for(size_t i = 0, j = 0; i < spaceMask.length(); i++){
        mask[i] = spaceMask[i] - 48;
        spaceNum += (mask[i] == 0);
        if(mask[i] == 1){
            unmaskedPos[j] = (int) i;
            j++;
        }
    }
    delete[] mask;
    maxGap = par.maxGap;
    accessionLevel = par.accessionLevel;
    minSSMatch = par.minSSMatch;
    minConsCnt = (size_t) par.minConsCnt;
    minConsCntEuk = (size_t) par.minConsCntEuk;
    eukaryotaTaxId = taxonomy->getEukaryotaTaxID();
    tieRatio = par.tieRatio;

    if (par.seqMode == 1 || par.seqMode == 2) {
        denominator = 100;
    } else {
        denominator = 1000;
    }

    if (par.reducedAA == 1) {
        bitsPerCodon = 4;
        totalDnaBits = 32;
        lastCodonMask = 0x0FFFFFFF;
    } else {
        bitsPerCodon = 3;
        totalDnaBits = 24;
        lastCodonMask = 0x1FFFFF;
    }

    // chooseBestTaxon
    taxCnt.reserve(4096);

    // getBestSpeciesMatches
    matchPaths.reserve(4096);
    combinedMatchPaths.reserve(4096);
    maxSpecies.reserve(4096);
    speciesList.reserve(4096);
    speciesScores.reserve(4096);

    // remainConsecutiveMatches
    linkedMatchKeys.reserve(4096);
    linkedMatchValues.reserve(4096);
    linkedMatchValuesIdx.reserve(4096);
    match2depthScore.reserve(4096); 

    // lowerRankClassification
    cladeCnt.reserve(4096);

    // filterRedundantMatches
    arraySize_filterRedundantMatches = 4096;
    bestMatchForQuotient = new const Match*[arraySize_filterRedundantMatches]();
    bestMatchTaxIdForQuotient = new TaxID[arraySize_filterRedundantMatches];
    minHammingForQuotient = new uint8_t[arraySize_filterRedundantMatches];


    // Output
    taxCounts.reserve(4096);
}

Taxonomer::~Taxonomer() {
    delete[] bestMatchForQuotient;
    delete[] bestMatchTaxIdForQuotient;
    delete[] minHammingForQuotient;
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
    cout << "Time spent for spliting matches: " << double(time(nullptr) - beforeAnalyze) << endl;
#pragma omp parallel default(none), shared(cout, matchBlocks, matchList, seqNum, queryList, blockIdx, par)
    {
        #pragma omp for schedule(dynamic, 1)
        for (size_t i = 0; i < blockIdx; ++i) {
            chooseBestTaxon(matchBlocks[i].id,
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

void Taxonomer::chooseBestTaxon(uint32_t currentQuery,
                                size_t offset,
                                size_t end,
                                const Match *matchList,
                                vector<Query> & queryList,
                                const LocalParameters &par) {
    // cout << "Current query: " << currentQuery << endl;
    // for (size_t i = offset; i < end+1; i ++) {
    //     matchList[i].printMatch();
    // }
    TaxonScore speciesScore(0, 0, 0, 0, 0);
    std::pair<size_t, size_t> bestSpeciesRange;
    speciesScore = getBestSpeciesMatches(bestSpeciesRange,
                                         matchList,
                                         end,
                                         offset,                        
                                         queryList[currentQuery].queryLength + queryList[currentQuery].queryLength2);
    
    // If there is no proper species for current query, it is un-classified.
    if (speciesScore.score == 0 || speciesScore.score < par.minScore) {
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
        queryList[currentQuery].isClassified = true;
        queryList[currentQuery].classification = speciesScore.taxId;
        queryList[currentQuery].score = speciesScore.score;
        queryList[currentQuery].coverage = speciesScore.coverage;
        queryList[currentQuery].hammingDist = speciesScore.hammingDist;
        return;
    }

    // Filter redundant matches
    taxCnt.clear();
    filterRedundantMatches(matchList,
                           bestSpeciesRange,
                           taxCnt,
                           queryList[currentQuery].queryLength + queryList[currentQuery].queryLength2);
    for (auto & tax : taxCnt) {
      queryList[currentQuery].taxCnt[tax.first] = tax.second;    
    }
    
    // If score is not enough, classify to the parent of the selected species
    if (speciesScore.score < par.minSpScore) {
      queryList[currentQuery].isClassified = true;
      queryList[currentQuery].classification = taxonomy->taxonNode(
              taxonomy->getTaxIdAtRank(speciesScore.taxId, "species"))->parentTaxId;
      queryList[currentQuery].score = speciesScore.score;
      queryList[currentQuery].coverage = speciesScore.coverage;
      queryList[currentQuery].hammingDist = speciesScore.hammingDist;
      return;
    }

    // Lower rank classification
    TaxID result = lowerRankClassification(taxCnt,
                                           speciesScore.taxId,
                                           queryList[currentQuery].queryLength + queryList[currentQuery].queryLength2);
    
    // TaxID result = speciesScore.taxId;

    // Store classification results
    queryList[currentQuery].isClassified = true;
    queryList[currentQuery].classification = result;
    queryList[currentQuery].score = speciesScore.score;
    queryList[currentQuery].coverage = speciesScore.coverage;
    queryList[currentQuery].hammingDist = speciesScore.hammingDist;
    queryList[currentQuery].newSpecies = false;
}

void Taxonomer::chooseBestTaxon_syncmer(
    uint32_t currentQuery,
    size_t offset,
    size_t end,
    const Match *matchList,
    vector<Query> & queryList,
    const LocalParameters &par) 
{    
    TaxonScore speciesScore(0, 0, 0, 0, 0);
    vector<std::pair<size_t, bool>> optimalMatchIdx;
    speciesScore = getBestSpeciesMatches2(
                        optimalMatchIdx, 
                        matchList,
                        end,
                        offset,                        
                        queryList[currentQuery].queryLength + queryList[currentQuery].queryLength2);
    
    // If there is no proper species for current query, it is un-classified.
    if (speciesScore.score == 0 || speciesScore.score < par.minScore) {
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
        queryList[currentQuery].isClassified = true;
        queryList[currentQuery].classification = speciesScore.taxId;
        queryList[currentQuery].score = speciesScore.score;
        queryList[currentQuery].coverage = speciesScore.coverage;
        queryList[currentQuery].hammingDist = speciesScore.hammingDist;
        return;
    }


    cout << "Match on best path" << endl;
    for (size_t i = 0; i < optimalMatchIdx.size(); i++) {
        matchList[optimalMatchIdx[i].first].printMatch();
    }
        
    taxCnt.clear();
    for (size_t i = 0; i < optimalMatchIdx.size(); i++) {
        if (optimalMatchIdx[i].second) {
            taxCnt[matchList[optimalMatchIdx[i].first].speciesId]++;
        } else {
            taxCnt[matchList[optimalMatchIdx[i].first].targetId]++;
        }
    }

    for (auto & tax : taxCnt) {
      queryList[currentQuery].taxCnt[tax.first] = tax.second;    
    }
    
    // If score is not enough, classify to the parent of the selected species
    if (speciesScore.score < par.minSpScore) {
      queryList[currentQuery].isClassified = true;
      queryList[currentQuery].classification = taxonomy->taxonNode(
              taxonomy->getTaxIdAtRank(speciesScore.taxId, "species"))->parentTaxId;
      queryList[currentQuery].score = speciesScore.score;
      queryList[currentQuery].coverage = speciesScore.coverage;
      queryList[currentQuery].hammingDist = speciesScore.hammingDist;
      return;
    }
    // Lower rank classification
    TaxID result = lowerRankClassification(taxCnt,
                                           speciesScore.taxId,
                                           queryList[currentQuery].queryLength + queryList[currentQuery].queryLength2);

    // Store classification results
    queryList[currentQuery].isClassified = true;
    queryList[currentQuery].classification = result;
    queryList[currentQuery].score = speciesScore.score;
    queryList[currentQuery].coverage = speciesScore.coverage;
    queryList[currentQuery].hammingDist = speciesScore.hammingDist;
    queryList[currentQuery].newSpecies = false;
}

void Taxonomer::filterRedundantMatches(const Match *matchList,
                                       const std::pair<size_t, size_t> & bestSpeciesRange,
                                       unordered_map<TaxID, unsigned int> & taxCnt,
                                       int queryLength) {    
    // Determine the maximum quotient we need to handle
    size_t maxQuotient = (queryLength + 3) / 6;
    
    ensureArraySize(maxQuotient + 1);

    // std::fill_n(bestMatchTaxIdForQuotient, maxQuotient + 1, 0);

    for (size_t i = bestSpeciesRange.first; i < bestSpeciesRange.second; i ++) {
        size_t currQuotient = matchList[i].qInfo.pos / 6;
        uint8_t hamming = matchList[i].hamming;

        if (bestMatchForQuotient[currQuotient] == nullptr) {
            bestMatchForQuotient[currQuotient] = matchList + i;
            bestMatchTaxIdForQuotient[currQuotient] = matchList[i].targetId;
            minHammingForQuotient[currQuotient] = hamming;
        } else {
            if (hamming < minHammingForQuotient[currQuotient]) {
                bestMatchForQuotient[currQuotient] = matchList + i;
                bestMatchTaxIdForQuotient[currQuotient] = matchList[i].targetId;
                minHammingForQuotient[currQuotient] = hamming;
            } else if (hamming == minHammingForQuotient[currQuotient]) {
                bestMatchTaxIdForQuotient[currQuotient] = taxonomy->LCA(
                    bestMatchTaxIdForQuotient[currQuotient], matchList[i].targetId);
            }
        }
    }

    for (size_t i = 0; i <= maxQuotient; ++i) {
        if (bestMatchForQuotient[i] != nullptr) {
            taxCnt[bestMatchTaxIdForQuotient[i]]++;
        }
    }
}

void Taxonomer::printSpeciesMatches(
    const Match *matchList,
    const std::pair<size_t, size_t> & bestSpeciesRange) 
{
    for (size_t i = bestSpeciesRange.first; i < bestSpeciesRange.second; i ++) {
        matchList[i].printMatch();
    }
}

TaxID Taxonomer::lowerRankClassification(const unordered_map<TaxID, unsigned int> & taxCnt, TaxID spTaxId, int queryLength) {
    unsigned int maxCnt = (queryLength - 1)/denominator + 1;
    cladeCnt.clear();
    getSpeciesCladeCounts(taxCnt, cladeCnt, spTaxId);
    if (accessionLevel == 2) { // Don't do accession-level classification
        // Remove leaf nodes
        for (auto it = cladeCnt.begin(); it != cladeCnt.end(); it++) {
            TaxonNode const * taxon = taxonomy->taxonNode(it->first);
            if (strcmp(taxonomy->getString(taxon->rankIdx), "") == 0 || strcmp(taxonomy->getString(taxon->rankIdx), "accession") == 0) {
                // Remove current node from its parent's children list
                cladeCnt[taxon->parentTaxId].children.erase(find(cladeCnt[taxon->parentTaxId].children.begin(),
                                                                 cladeCnt[taxon->parentTaxId].children.end(),
                                                                 it->first));
            } 
        }
        return BFS(cladeCnt, spTaxId, maxCnt);
    } else {
        return BFS(cladeCnt, spTaxId, maxCnt);
    }
}

void Taxonomer::getSpeciesCladeCounts(const unordered_map<TaxID, unsigned int> &taxCnt,
                                      unordered_map<TaxID, TaxonCounts> & cladeCount,
                                      TaxID speciesTaxID) {
    for (auto it = taxCnt.begin(); it != taxCnt.end(); ++it) {
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

TaxID Taxonomer::BFS(const unordered_map<TaxID, TaxonCounts> & cladeCnt, TaxID root, unsigned int maxCnt) {
    unsigned int maxCnt2 = maxCnt;
    if (cladeCnt.at(root).children.empty()) { // root is a leaf
        return root;
    }
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
        return BFS(cladeCnt, bestChildren[0], maxCnt2);
    } else {
        return root;
    }
}

TaxonScore Taxonomer::getBestSpeciesMatches(std::pair<size_t, size_t> & bestSpeciesRange,
                                            const Match *matchList,
                                            size_t end,
                                            size_t offset,
                                            int queryLength) {
    matchPaths.clear();
    combinedMatchPaths.clear();
    speciesList.clear();
    speciesScores.clear();
    
    TaxonScore bestScore;
    float bestSpScore = 0;
    size_t i = offset;
    size_t meaningfulSpecies = 0;
    while (i  < end + 1) {
        TaxID currentSpecies = matchList[i].speciesId;
        size_t start = i;
        size_t previousPathSize = matchPaths.size();
        // For current species
        while ((i < end + 1) && currentSpecies == matchList[i].speciesId) {
            uint8_t curFrame = matchList[i].qInfo.frame;
            size_t start = i;
            // For current frame
            while ((i < end + 1) && currentSpecies == matchList[i].speciesId && curFrame == matchList[i].qInfo.frame) {
                i ++;
            }
            if (i - start >= minConsCnt) {
                getMatchPaths(matchList, start, i, matchPaths, currentSpecies);
            }
        }
        size_t pathSize = matchPaths.size();
        // Combine MatchPaths
        if (par.printLog) {
            cout << "Current species: " << taxonomy->getOriginalTaxID(currentSpecies) << " " << currentSpecies << endl;
            for (size_t kk = previousPathSize; kk < matchPaths.size(); kk++) {
                matchPaths[kk].printMatchPath();
            }
        }
        if (pathSize > previousPathSize) {
            speciesList.push_back(currentSpecies);
            float score = combineMatchPaths(matchPaths, previousPathSize, combinedMatchPaths, combinedMatchPaths.size(), queryLength);
            score = min(score, 1.0f);
            speciesScores.push_back(score);
            if (score > 0.f) {
                meaningfulSpecies++;
            }
            if (score > bestSpScore) {
                bestSpScore = score;
                bestSpeciesRange = make_pair(start, i);
            }
        }
    }
    
    // If there are no meaningful species
    if (meaningfulSpecies == 0) {
        bestScore.score = 0;
        return bestScore;
    }

    maxSpecies.clear();
    // float coveredLength = 0.f;
    for (size_t i = 0; i < speciesList.size(); i++) {
        if (speciesScores[i] >= bestSpScore * tieRatio) {
            maxSpecies.push_back(speciesList[i]);
            bestScore.score += speciesScores[i];   
        }
    }
    
    // More than one species --> LCA    
    if (maxSpecies.size() > 1) {
        bestScore.LCA = true;
        bestScore.taxId = taxonomy->LCA(maxSpecies)->taxId;
        bestScore.score /= maxSpecies.size();
        // bestScore.coverage /= maxSpecies.size();
        return bestScore;
    }
    
    // One species
    bestScore.taxId = maxSpecies[0];
    
    return bestScore;                                  
}

TaxonScore Taxonomer::getBestSpeciesMatches2(
    vector<std::pair<size_t, bool>> & optimalMatchIdx,
    const Match *matchList,
    size_t end,
    size_t offset,
    int queryLength) 
{
    matchPaths2.clear();
    speciesList.clear();
    speciesScores.clear();
    vector<size_t> selectedMatchPaths;
    vector<size_t> bestSpeciesMatchPaths;
    std::pair<size_t, size_t> bestSpeciesRange;
    TaxonScore bestScore;
    float bestSpScore = 0;
    size_t i = offset;
    size_t meaningfulSpecies = 0;
    while (i  < end + 1) {
        TaxID currentSpecies = matchList[i].speciesId;
        size_t start = i;
        size_t previousPathSize = matchPaths2.size();
        // For current species
        while ((i < end + 1) && currentSpecies == matchList[i].speciesId) {
            uint8_t curFrame = matchList[i].qInfo.frame;
            size_t start = i;
            // For current frame
            while ((i < end + 1) && currentSpecies == matchList[i].speciesId && curFrame == matchList[i].qInfo.frame) {
                i ++;
            }
            if (i - start >= minConsCnt) {
                getMatchPaths2(matchList, start, i, matchPaths2, currentSpecies);
            }
        }
        size_t pathSize = matchPaths2.size();
        // Combine MatchPaths
        if (par.printLog && pathSize > previousPathSize) {
            cout << "Current species: " << taxonomy->getOriginalTaxID(currentSpecies) << " " << currentSpecies << endl;
            for (size_t kk = previousPathSize; kk < matchPaths2.size(); kk++) {
                matchPaths2[kk].printMatchPath();
            }
        }
        if (pathSize > previousPathSize) {
            selectedMatchPaths.clear();
            speciesList.push_back(currentSpecies);
            float score = combineMatchPaths2(matchPaths2, selectedMatchPaths, previousPathSize, queryLength, matchList);
            score = min(score, 1.0f);
            speciesScores.push_back(score);
            if (score > 0.f) {
                meaningfulSpecies++;
            }
            if (score > bestSpScore) {
                bestSpeciesMatchPaths = selectedMatchPaths;
                bestSpScore = score;
                bestSpeciesRange = make_pair(start, i);
            }
        }
    }

    for (size_t i = 0; i < bestSpeciesMatchPaths.size(); i++) {
        for (size_t j = 0; j < matchPaths2[bestSpeciesMatchPaths[i]].matchIdxs.size(); j++) {
            matchList[matchPaths2[bestSpeciesMatchPaths[i]].matchIdxs[j]].printMatch();
            cout << matchPaths2[bestSpeciesMatchPaths[i]].cumulativeScore[j] << endl;
            optimalMatchIdx.emplace_back(matchPaths2[bestSpeciesMatchPaths[i]].matchIdxs[j], matchPaths2[bestSpeciesMatchPaths[i]].tied);
        }
    }

    printSpeciesMatches(matchList, bestSpeciesRange);
    
    // If there are no meaningful species
    if (meaningfulSpecies == 0) {
        bestScore.score = 0;
        return bestScore;
    }

    maxSpecies.clear();
    // float coveredLength = 0.f;
    for (size_t i = 0; i < speciesList.size(); i++) {
        if (speciesScores[i] >= bestSpScore * tieRatio) {
            maxSpecies.push_back(speciesList[i]);
            bestScore.score += speciesScores[i]; 
        }
    }
    
    // More than one species --> LCA    
    if (maxSpecies.size() > 1) {
        bestScore.LCA = true;
        bestScore.taxId = taxonomy->LCA(maxSpecies)->taxId;
        bestScore.score /= maxSpecies.size();
        // bestScore.coverage /= maxSpecies.size();
        return bestScore;
    }
    
    // One species
    bestScore.taxId = maxSpecies[0];
    
    return bestScore;                                  
}


float Taxonomer::combineMatchPaths2(
    vector<MatchPath2> & matchPaths,
    vector<size_t> & selectedMatchPaths,
    size_t matchPathStart,
    int readLength,
    const Match *matchList) 
{
    sort(matchPaths.begin() + matchPathStart, matchPaths.end(),
         [](const MatchPath2 &a, const MatchPath2 &b) {
           if (a.score != b.score) {
             return a.score > b.score;
           }
           if (a.hammingDist != b.hammingDist) {
             return a.hammingDist < b.hammingDist;
           }
           return a.start > b.start;
         });

    // Combine matchPaths
    // 1. Add the matchPath with the highest score to combinedMatchPaths
    // 2. Add the matchPath with the highest score that is not overlapped with the matchPath in combinedMatchPaths
    // 3. Repeat 2 until no matchPath can be added
    float score = 0;    
    for (size_t i = matchPathStart; i < matchPaths.size(); i++) {  
        if (selectedMatchPaths.empty()) {
            selectedMatchPaths.push_back(i);
            cout << "Path added ";
            matchPaths[i].printMatchPath();
            score += matchPaths[i].score;
        } else {
            bool isOverlapped = false;
            for (size_t j = 0; j < selectedMatchPaths.size(); j++) {
                if (isMatchPathOverlapped2(matchPaths[i], matchPaths[selectedMatchPaths[j]])) {
                    compareOverlapRegion(matchPaths[i], matchPaths[selectedMatchPaths[j]], matchList);
                    // Compare the score overlapped region 

                    if (MatchPath2::isTied(matchPaths[selectedMatchPaths[j]], matchPaths[i])) {
                        matchPaths[selectedMatchPaths[j]].tied = true;
                        isOverlapped = true;
                        break;
                    }
                    int overlappedLength = min(matchPaths[i].end, matchPaths[selectedMatchPaths[j]].end) 
                                            - max(matchPaths[i].start, matchPaths[selectedMatchPaths[j]].start) + 1;
                    if (overlappedLength == matchPaths[i].end - matchPaths[i].start + 1) { // Current path is completely overlapped
                        isOverlapped = true;
                        break;
                    }
                    if (overlappedLength < 24) {
                        trimMatchPath2(matchPaths[i], matchPaths[selectedMatchPaths[j]], overlappedLength);
                        cout << "Path trimmed "; matchPaths[i].printMatchPath();
                        continue;
                    } else {
                        isOverlapped = true;
                        break;
                    }
                    if (matchPaths[i].end - matchPaths[i].start + 1 < 24) { // Current path trimmed too much 
                        cout << "Path trimmed too much ";
                        matchPaths[i].printMatchPath();
                        isOverlapped = true;
                        break;
                    }  
                } 
            }
            if (!isOverlapped) {
                selectedMatchPaths.push_back(i);
                cout << "Path added "; matchPaths[i].printMatchPath();
                score += matchPaths[i].score;           
            }
        }
    }
    return score / readLength;
}

float Taxonomer::combineMatchPaths(vector<MatchPath> & matchPaths,
                                   size_t matchPathStart,
                                   vector<MatchPath> & combinedMatchPaths,
                                   size_t combMatchPathStart,
                                   int readLength) {
    // cout << "combineMatchPaths" << endl;
    // Sort matchPaths by the their score
    sort(matchPaths.begin() + matchPathStart, matchPaths.end(),
         [](const MatchPath &a, const MatchPath &b) {
           if (a.score != b.score) {
             return a.score > b.score;
           }
           if (a.hammingDist != b.hammingDist) {
             return a.hammingDist < b.hammingDist;
           }
           return a.start > b.start;
         });

    // Combine matchPaths
    // 1. Add the matchPath with the highest score to combinedMatchPaths
    // 2. Add the matchPath with the highest score that is not overlapped with the matchPath in combinedMatchPaths
    // 3. Repeat 2 until no matchPath can be added
    float score = 0;
    for (size_t i = matchPathStart; i < matchPaths.size(); i++) {  
        if (combMatchPathStart == combinedMatchPaths.size()) {
            combinedMatchPaths.push_back(matchPaths[i]);
            // cout << "Path added ";
            // matchPaths[i].printMatchPath();
            score += matchPaths[i].score;
        } else {
            bool isOverlapped = false;
            for (size_t j = combMatchPathStart; j < combinedMatchPaths.size(); j++) {
                if (isMatchPathOverlapped(matchPaths[i], combinedMatchPaths[j])) { // overlap!
                    int overlappedLength = min(matchPaths[i].end, combinedMatchPaths[j].end) 
                                            - max(matchPaths[i].start, combinedMatchPaths[j].start) + 1;
                    if (overlappedLength == matchPaths[i].end - matchPaths[i].start + 1) { // Current path is completely overlapped
                        isOverlapped = true;
                        break;
                    }
                    
                    if (overlappedLength < 24) {
                        trimMatchPath(matchPaths[i], combinedMatchPaths[j], overlappedLength);
                        // cout << "Path trimmed "; matchPaths[i].printMatchPath();
                        continue;
                    } else {
                        isOverlapped = true;
                        break;
                    }
                    if (matchPaths[i].end - matchPaths[i].start + 1 < 24) { // Current path trimmed too much 
                        isOverlapped = true;
                        break;
                    }  
                } 
            }
            if (!isOverlapped) {
                combinedMatchPaths.push_back(matchPaths[i]);
                // cout << "Path added "; matchPaths[i].printMatchPath();
                score += matchPaths[i].score;           
            }
        }
    }
    return score / readLength;
}

bool Taxonomer::isMatchPathOverlapped(const MatchPath & matchPath1,
                                      const MatchPath & matchPath2) {
    return !((matchPath1.end < matchPath2.start) || (matchPath2.end < matchPath1.start));                                       
}

bool Taxonomer::isMatchPathOverlapped2(const MatchPath2 & matchPath1,
                                      const MatchPath2 & matchPath2) {
    return !((matchPath1.end < matchPath2.start) || (matchPath2.end < matchPath1.start));                                       
}

void Taxonomer::compareOverlapRegion(
    MatchPath2 & matchPath1, 
    MatchPath2 & matchPath2,
    const Match *matchList) 
{   
    

    int overlapStart = max(matchPath1.start, matchPath2.start);
    int overlapEnd = min(matchPath1.end, matchPath2.end);
    int overlapLength = overlapEnd - overlapStart + 1;
    cout << "Overlap region: " << overlapStart << " " << overlapEnd << endl;
    matchPath1.printMatchPath();
    matchPath2.printMatchPath();
    float scoreOfOverlap1, scoreOfOverlap2;
    if (overlapLength <= 24) {
        if (overlapStart == matchPath1.start) {
            scoreOfOverlap1 = matchList[matchPath1.matchIdxs.front()].getScore();
            scoreOfOverlap2 = matchList[matchPath2.matchIdxs.back()].getScore();
        } else {
            scoreOfOverlap1 = matchList[matchPath1.matchIdxs.back()].getScore();
            scoreOfOverlap2 = matchList[matchPath2.matchIdxs.front()].getScore();
        }
    } else {

    

    bool endFound = false;
    bool startFound = false;
    size_t overlapEndIdx = 0;
    size_t overlapStartIdx = 0;
    
    cout << "Path 1" << endl;
    for (size_t i = 0; i < matchPath1.matchIdxs.size(); i++) {
        const Match * match = matchList + matchPath1.matchIdxs[i];
        if (!endFound && match->qInfo.pos + 23  <= overlapEnd) {
            endFound = true;
            overlapEndIdx = i;
            cout << "End found: " << match->qInfo.pos + 23 << endl;
            cout << "score: " << matchPath1.cumulativeScore[overlapEndIdx] << endl;
        }
        if (endFound && match->qInfo.pos <= overlapStart) {
            overlapStartIdx = i;
            cout << "Start found: " << match->qInfo.pos << endl;
            cout << "score: " << matchPath1.cumulativeScore[overlapStartIdx] << endl;
            startFound = true;
            break;
        }
    }
    if (startFound) {
        scoreOfOverlap1 = matchPath1.cumulativeScore[overlapEndIdx] - 
            (matchPath1.cumulativeScore[overlapStartIdx] - matchList[matchPath1.matchIdxs[overlapStartIdx]].getScore());
    } else {
        scoreOfOverlap1 = matchPath1.cumulativeScore[overlapEndIdx];
    }
    startFound = false;
    endFound = false;
    cout << "Path 2" << endl;
    for (size_t i = 0; i < matchPath2.matchIdxs.size(); i++) {
        const Match * match = matchList + matchPath2.matchIdxs[i];
        if (!endFound && match->qInfo.pos + 23 <= overlapEnd) {
            endFound = true;
            overlapEndIdx = i;
            cout << "End found: " << match->qInfo.pos + 23 << endl;
            cout << "score: " << matchPath2.cumulativeScore[overlapEndIdx] << endl;
        }
        if (endFound && match->qInfo.pos <= overlapStart) {
            overlapStartIdx = i;
            startFound = true;
            cout << "Start found: " << match->qInfo.pos << endl;
            cout << "score: " << matchPath2.cumulativeScore[overlapStartIdx] << endl;
            break;
        }
    }
    if (startFound) {
        scoreOfOverlap2 = matchPath2.cumulativeScore[overlapEndIdx] - 
            (matchPath2.cumulativeScore[overlapStartIdx] - matchList[matchPath2.matchIdxs[overlapStartIdx]].getScore());
    } else {
        scoreOfOverlap2 = matchPath2.cumulativeScore[overlapEndIdx];
    }
    // float scoreOfOverlap2 = matchPath2.cumulativeScore[overlapEndIdx] - matchPath2.cumulativeScore[overlapStartIdx];
    
    cout << "MatchPath 1 score: " << scoreOfOverlap1 << endl;
    cout << "MatchPath 2 score: " << scoreOfOverlap2 << endl;
}
}

void Taxonomer::trimMatchPath(MatchPath & path1, const MatchPath & path2, int overlapLength) {
    if (path1.start < path2.start) { 
        path1.end = path2.start - 1;
        path1.hammingDist = max(0, path1.hammingDist - path1.endMatch->getRightPartHammingDist(overlapLength/3));
        path1.score = path1.score - path1.endMatch->getRightPartScore(overlapLength/3) - (overlapLength % 3);
    } else {
        path1.start = path2.end + 1;
        path1.hammingDist = max(0, path1.hammingDist - path1.startMatch->getLeftPartHammingDist(overlapLength/3));
        path1.score = path1.score - path1.startMatch->getLeftPartScore(overlapLength/3) - (overlapLength % 3);
    }
}

void Taxonomer::trimMatchPath2(MatchPath2 & path1, const MatchPath2 & path2, int overlapLength) {
    if (path1.start < path2.start) { 
        path1.end = path2.start - 1;
        path1.hammingDist = max(0, path1.hammingDist - path1.endMatch->getRightPartHammingDist(overlapLength/3));
        path1.score = path1.score - path1.endMatch->getRightPartScore(overlapLength/3) - (overlapLength % 3);
    } else {
        path1.start = path2.end + 1;
        path1.hammingDist = max(0, path1.hammingDist - path1.startMatch->getLeftPartHammingDist(overlapLength/3));
        path1.score = path1.score - path1.startMatch->getLeftPartScore(overlapLength/3) - (overlapLength % 3);
    }
}

void Taxonomer::getMatchPaths2(
    const Match * matchList,
    size_t start,
    size_t end,
    vector<MatchPath2> & filteredMatchPaths,
    TaxID speciesId) 
{
    size_t i = start;
    size_t currPos = matchList[start].qInfo.pos;
    uint64_t frame = matchList[start].qInfo.frame;
    int MIN_DEPTH = (int) minConsCnt - 1;
    if (taxonomy->IsAncestor(eukaryotaTaxId, speciesId)) {
        MIN_DEPTH = (int) minConsCntEuk - 1;
    }
        
    connectedToNext.resize(end - start + 1);
    fill(connectedToNext.begin(), connectedToNext.end(), false);
    localMatchPaths2.clear();
    localMatchPaths2.resize(end - start + 1);

    if (frame < 3) { // Forward frame
        size_t curPosMatchStart = i;        
        while (i < end && matchList[i].qInfo.pos == currPos) {
            localMatchPaths2[i - start] = MatchPath2(matchList + i, i);
            ++ i;
        }
        size_t curPosMatchEnd = i; // exclusive

        while (i < end) {
            uint32_t nextPos = matchList[i].qInfo.pos;
            size_t nextPosMatchStart = i;
            while (i < end  && nextPos == matchList[i].qInfo.pos) {
                localMatchPaths2[i - start] = MatchPath2(matchList + i, i);
                ++ i;
            }
            size_t nextPosMatchEnd = i; // exclusive

            // Check if current position and next position are consecutive
            int shift = (nextPos - currPos) / 3;
            if (shift > 0 && shift <= par.maxShift) {
                // Compare curPosMatches and nextPosMatches.
                for (size_t nextIdx = nextPosMatchStart; nextIdx < nextPosMatchEnd; nextIdx++) {
                    float scoreIncrement = calScoreIncrement(matchList[nextIdx].rightEndHamming, shift);
                    const MatchPath2 * bestPath = nullptr;
                    float bestScore = 0;
                    for (size_t curIdx = curPosMatchStart; curIdx < curPosMatchEnd; ++curIdx) {
                        if (isConsecutive(matchList + curIdx, matchList + nextIdx, shift)) {
                            connectedToNext[curIdx - start] = true;
                            if (localMatchPaths2[curIdx - start].score > bestScore) {
                                bestPath = &localMatchPaths2[curIdx - start];
                                bestScore = localMatchPaths2[curIdx - start].score;
                            }
                        }
                    }
                    if (bestPath != nullptr) {
                        localMatchPaths2[nextIdx - start].start = bestPath->start;                        
                        localMatchPaths2[nextIdx - start].score = bestPath->score + scoreIncrement;
                        localMatchPaths2[nextIdx - start].hammingDist = bestPath->hammingDist + calHammingDistIncrement(matchList[nextIdx].rightEndHamming, shift);
                        localMatchPaths2[nextIdx - start].depth = bestPath->depth + shift;
                        localMatchPaths2[nextIdx - start].startMatch = bestPath->startMatch;
                        localMatchPaths2[nextIdx - start].prevMatchPath = bestPath;
                    }
                }
            } 
            for (size_t curIdx = curPosMatchStart; curIdx < curPosMatchEnd; ++curIdx) {
                if (!connectedToNext[curIdx - start] && localMatchPaths2[curIdx - start].depth > MIN_DEPTH) {
                    filteredMatchPaths.push_back(localMatchPaths2[curIdx - start]);
                    filteredMatchPaths.back().tied = false;
                    const MatchPath2 * p = &localMatchPaths2[curIdx - start];
                    do {
                        filteredMatchPaths.back().matchIdxs.push_back(p->matchIdx);
                        filteredMatchPaths.back().cumulativeScore.push_back(p->score);
                        p = p->prevMatchPath;
                    } while (p != nullptr);
                }
            }
            if (i == end) {
                for (size_t nextIdx = nextPosMatchStart; nextIdx < nextPosMatchEnd; ++nextIdx) {
                    if (localMatchPaths2[nextIdx - start].depth > MIN_DEPTH) {
                        filteredMatchPaths.push_back(localMatchPaths2[nextIdx - start]);
                        filteredMatchPaths.back().tied = false;
                        const MatchPath2 * p = &localMatchPaths2[nextIdx - start];
                        do {
                            filteredMatchPaths.back().matchIdxs.push_back(p->matchIdx);
                            filteredMatchPaths.back().cumulativeScore.push_back(p->score);
                            p = p->prevMatchPath;
                        } while (p != nullptr);
                    }
                }
            }
            curPosMatchStart = nextPosMatchStart;
            curPosMatchEnd = nextPosMatchEnd;
            currPos = nextPos;    
        }
    } else {
        size_t curPosMatchStart = i;        
        while (i < end && matchList[i].qInfo.pos == currPos) {
            localMatchPaths2[i - start] = MatchPath2(matchList + i, i);
            ++ i;
        }
        size_t curPosMatchEnd = i; // exclusive

        while (i < end) {
            uint32_t nextPos = matchList[i].qInfo.pos;
            size_t nextPosMatchStart = i;
            while (i < end  && nextPos == matchList[i].qInfo.pos) {
                localMatchPaths2[i - start] = MatchPath2(matchList + i, i);
                ++ i;
            }
            size_t nextPosMatchEnd = i; // exclusive

            // Check if current position and next position are consecutive
            int shift = (nextPos - currPos) / 3;
            if (shift > 0 && shift <= par.maxShift) {
                for (size_t nextIdx = nextPosMatchStart; nextIdx < nextPosMatchEnd; nextIdx++) {
                    float scoreIncrement = calScoreIncrement(matchList[nextIdx].rightEndHamming, shift);
                    const MatchPath2 * bestPath = nullptr;
                    float bestScore = 0;
                    for (size_t curIdx = curPosMatchStart; curIdx < curPosMatchEnd; ++curIdx) {
                        if (isConsecutive(matchList + nextIdx, matchList + curIdx, shift)) {
                            connectedToNext[curIdx - start] = true;
                            if (localMatchPaths2[curIdx - start].score > bestScore) {
                                bestPath = &localMatchPaths2[curIdx - start];
                                bestScore = localMatchPaths2[curIdx - start].score;
                            }
                        }
                    }
                    if (bestPath != nullptr) {
                        localMatchPaths2[nextIdx - start].start = bestPath->start;                        
                        localMatchPaths2[nextIdx - start].score = bestPath->score + scoreIncrement;
                        localMatchPaths2[nextIdx - start].hammingDist = bestPath->hammingDist + calHammingDistIncrement(matchList[nextIdx].rightEndHamming, shift);
                        localMatchPaths2[nextIdx - start].depth = bestPath->depth + shift;
                        localMatchPaths2[nextIdx - start].startMatch = bestPath->startMatch;
                        localMatchPaths2[nextIdx - start].prevMatchPath = bestPath;
                    }
                }
            } 
            for (size_t curIdx = curPosMatchStart; curIdx < curPosMatchEnd; ++curIdx) {
                if (!connectedToNext[curIdx - start] && localMatchPaths2[curIdx - start].depth > MIN_DEPTH) {
                    filteredMatchPaths.push_back(localMatchPaths2[curIdx - start]);
                    filteredMatchPaths.back().tied = false;
                    const MatchPath2 * p = &localMatchPaths2[curIdx - start];
                    do {
                        filteredMatchPaths.back().matchIdxs.push_back(p->matchIdx);
                        filteredMatchPaths.back().cumulativeScore.push_back(p->score);
                        p = p->prevMatchPath;
                    } while (p != nullptr);
                    
                }
            }
            if (i == end) {
                for (size_t nextIdx = nextPosMatchStart; nextIdx < nextPosMatchEnd; ++nextIdx) {
                    if (localMatchPaths2[nextIdx - start].depth > MIN_DEPTH) {
                        filteredMatchPaths.push_back(localMatchPaths2[nextIdx - start]);
                        filteredMatchPaths.back().tied = false;
                        const MatchPath2 * p = &localMatchPaths2[nextIdx - start];
                        do {
                            filteredMatchPaths.back().matchIdxs.push_back(p->matchIdx);
                            filteredMatchPaths.back().cumulativeScore.push_back(p->score);
                            p = p->prevMatchPath;
                        } while (p != nullptr);
                    }
                }
            }
            curPosMatchStart = nextPosMatchStart;
            curPosMatchEnd = nextPosMatchEnd;
            currPos = nextPos;    
        }
    }
}

void Taxonomer::getMatchPaths(
    const Match * matchList,
    size_t start,
    size_t end,
    vector<MatchPath> & filteredMatchPaths,
    TaxID speciesId) 
{
    size_t i = start;
    size_t currPos = matchList[start].qInfo.pos;
    uint64_t frame = matchList[start].qInfo.frame;
    int MIN_DEPTH = (int) minConsCnt - 1;
    if (taxonomy->IsAncestor(eukaryotaTaxId, speciesId)) {
        MIN_DEPTH = (int) minConsCntEuk - 1;
    }
        
    connectedToNext.resize(end - start + 1);
    fill(connectedToNext.begin(), connectedToNext.end(), false);
    localMatchPaths.clear();
    localMatchPaths.resize(end - start + 1);

    if (frame < 3) { // Forward frame
        size_t curPosMatchStart = i;        
        while (i < end && matchList[i].qInfo.pos == currPos) {
            localMatchPaths[i - start] = MatchPath(matchList + i);
            ++ i;
        }
        size_t curPosMatchEnd = i; // exclusive

        while (i < end) {
            uint32_t nextPos = matchList[i].qInfo.pos;
            size_t nextPosMatchStart = i;
            while (i < end  && nextPos == matchList[i].qInfo.pos) {
                localMatchPaths[i - start] = MatchPath(matchList + i);
                ++ i;
            }
            size_t nextPosMatchEnd = i; // exclusive

            // Check if current position and next position are consecutive
            int shift = (nextPos - currPos) / 3;
            if (shift > 0 && shift <= par.maxShift) {
                // Compare curPosMatches and nextPosMatches.
                for (size_t nextIdx = nextPosMatchStart; nextIdx < nextPosMatchEnd; nextIdx++) {
                    // uint8_t lastEndHamming = (matchList[nextIdx].rightEndHamming >> 14);
                    // float scoreIncrement = (lastEndHamming == 0) ? 3.0f : 2.0f - 0.5f * lastEndHamming;
                    float scoreIncrement = calScoreIncrement(matchList[nextIdx].rightEndHamming, shift);
                    const MatchPath * bestPath = nullptr;
                    float bestScore = 0;
                    for (size_t curIdx = curPosMatchStart; curIdx < curPosMatchEnd; ++curIdx) {
                        if (isConsecutive(matchList + curIdx, matchList + nextIdx, shift)) {
                            connectedToNext[curIdx - start] = true;
                            if (localMatchPaths[curIdx - start].score > bestScore) {
                                bestPath = &localMatchPaths[curIdx - start];
                                bestScore = localMatchPaths[curIdx - start].score;
                            }
                        }
                    }
                    if (bestPath != nullptr) {
                        localMatchPaths[nextIdx - start].start = bestPath->start;                        
                        localMatchPaths[nextIdx - start].score = bestPath->score + scoreIncrement;
                        localMatchPaths[nextIdx - start].hammingDist = bestPath->hammingDist + calHammingDistIncrement(matchList[nextIdx].rightEndHamming, shift);
                        localMatchPaths[nextIdx - start].depth = bestPath->depth + shift;
                        localMatchPaths[nextIdx - start].startMatch = bestPath->startMatch;
                    }
                }
            } 
            for (size_t curIdx = curPosMatchStart; curIdx < curPosMatchEnd; ++curIdx) {
                if (!connectedToNext[curIdx - start] && localMatchPaths[curIdx - start].depth > MIN_DEPTH) {
                    filteredMatchPaths.push_back(localMatchPaths[curIdx - start]);
                }
            }
            if (i == end) {
                for (size_t nextIdx = nextPosMatchStart; nextIdx < nextPosMatchEnd; ++nextIdx) {
                    if (localMatchPaths[nextIdx - start].depth > MIN_DEPTH) {
                        filteredMatchPaths.push_back(localMatchPaths[nextIdx - start]);
                    }
                }
            }
            curPosMatchStart = nextPosMatchStart;
            curPosMatchEnd = nextPosMatchEnd;
            currPos = nextPos;    
        }
    } else {
        size_t curPosMatchStart = i;        
        while (i < end && matchList[i].qInfo.pos == currPos) {
            localMatchPaths[i - start] = MatchPath(matchList + i);
            ++ i;
        }
        size_t curPosMatchEnd = i; // exclusive

        while (i < end) {
            uint32_t nextPos = matchList[i].qInfo.pos;
            size_t nextPosMatchStart = i;
            while (i < end  && nextPos == matchList[i].qInfo.pos) {
                localMatchPaths[i - start] = MatchPath(matchList + i);
                ++ i;
            }
            size_t nextPosMatchEnd = i; // exclusive

            // Check if current position and next position are consecutive
            int shift = (nextPos - currPos) / 3;
            if (shift > 0 && shift <= par.maxShift) {
                for (size_t nextIdx = nextPosMatchStart; nextIdx < nextPosMatchEnd; nextIdx++) {
                    float scoreIncrement = calScoreIncrement(matchList[nextIdx].rightEndHamming, shift);
                    const MatchPath * bestPath = nullptr;
                    float bestScore = 0;
                    for (size_t curIdx = curPosMatchStart; curIdx < curPosMatchEnd; ++curIdx) {
                        if (isConsecutive(matchList + nextIdx, matchList + curIdx, shift)) {
                            connectedToNext[curIdx - start] = true;
                            if (localMatchPaths[curIdx - start].score > bestScore) {
                                bestPath = &localMatchPaths[curIdx - start];
                                bestScore = localMatchPaths[curIdx - start].score;
                            }
                        }
                    }
                    if (bestPath != nullptr) {
                        localMatchPaths[nextIdx - start].start = bestPath->start;                        
                        localMatchPaths[nextIdx - start].score = bestPath->score + scoreIncrement;
                        localMatchPaths[nextIdx - start].hammingDist = bestPath->hammingDist + calHammingDistIncrement(matchList[nextIdx].rightEndHamming, shift);
                        localMatchPaths[nextIdx - start].depth = bestPath->depth + shift;
                        localMatchPaths[nextIdx - start].startMatch = bestPath->startMatch;
                    }
                }
            } 
            for (size_t curIdx = curPosMatchStart; curIdx < curPosMatchEnd; ++curIdx) {
                if (!connectedToNext[curIdx - start] && localMatchPaths[curIdx - start].depth > MIN_DEPTH) {
                    filteredMatchPaths.push_back(localMatchPaths[curIdx - start]);
                }
            }
            if (i == end) {
                for (size_t nextIdx = nextPosMatchStart; nextIdx < nextPosMatchEnd; ++nextIdx) {
                    if (localMatchPaths[nextIdx - start].depth > MIN_DEPTH) {
                        filteredMatchPaths.push_back(localMatchPaths[nextIdx - start]);
                    }
                }
            }
            curPosMatchStart = nextPosMatchStart;
            curPosMatchEnd = nextPosMatchEnd;
            currPos = nextPos;    
        }
    }
}

float Taxonomer::calScoreIncrement(uint16_t hammings, int shift) {
    float scoreIncrement = 0;
    for (int i = 0; i < shift; i++) {
        uint8_t hamming = (hammings >> (14 - i * 2)) & 0b11;
        if (hamming == 0) {
            scoreIncrement += 3.0f;
        } else {
            scoreIncrement += 2.0f - 0.5f * hamming;
        }
    }
    return scoreIncrement;
}

int Taxonomer::calHammingDistIncrement(uint16_t hammings, int shift) {
    int hammingDistIncrement = 0;
    for (int i = 0; i < shift; i++) {
        hammingDistIncrement += (hammings >> (14 - i * 2)) & 0b11;
    }
    return hammingDistIncrement;
}

void Taxonomer::remainConsecutiveMatches(const Match * matchList,
                                         size_t start,
                                         size_t end,
                                         vector<MatchPath> & matchPaths,
                                         TaxID speciesId) {
    size_t i = start;
    linkedMatchKeys.clear();
    linkedMatchValues.clear();
    linkedMatchValuesIdx.clear();
    
    linkedMatchKeys.resize(end - start);
    linkedMatchValuesIdx.resize(end - start);
    size_t linkedMatchKeysIdx = 0;

    size_t currPos = matchList[start].qInfo.pos;
    uint64_t frame = matchList[start].qInfo.frame;
    if (frame < 3) { // Forward frame
        size_t curPosMatchStart = i;        
        while (i < end && matchList[i].qInfo.pos == currPos) {
            ++ i;
        }
        size_t curPosMatchEnd = i; // exclusive

        while (i < end) {
            uint32_t nextPos = matchList[i].qInfo.pos;
            size_t nextPosMatchStart = i;
            while (i < end  && nextPos == matchList[i].qInfo.pos) {
                ++ i;
            }
            size_t nextPosMatchEnd = i; // exclusive

            // Check if current position and next position are consecutive
            if (currPos + 3 == nextPos) {
                // Compare curPosMatches and nextPosMatches
                for (size_t curIdx = curPosMatchStart; curIdx < curPosMatchEnd; ++curIdx) {
                    size_t startIdx = linkedMatchValues.size();
                    bool found = false;
                    for (size_t nextIdx = nextPosMatchStart; nextIdx < nextPosMatchEnd; nextIdx++) {
                        if (isConsecutive(matchList + curIdx, matchList + nextIdx)){
                            linkedMatchValues.push_back(matchList + nextIdx);
                            found = true;
                        }
                    }
                    if (found) {
                        linkedMatchKeys[linkedMatchKeysIdx] = matchList + curIdx;
                        linkedMatchValuesIdx[linkedMatchKeysIdx++] = startIdx;
                    }
                }
            }
            // Update curPosMatches and nextPosMatches
            curPosMatchStart = nextPosMatchStart;
            curPosMatchEnd = nextPosMatchEnd;
            currPos = nextPos;
        }
    } else {
        size_t curPosMatchStart = i;

        while ( i < end && matchList[i].qInfo.pos == currPos) {
            i++;
        }
        size_t curPosMatchEnd = i; // exclusive

        while (i < end) {
            uint32_t nextPos = matchList[i].qInfo.pos;
            size_t nextPosMatchStart = i;
            while (i < end  && nextPos == matchList[i].qInfo.pos) {
                ++ i;
            }
            size_t nextPosMatchEnd = i; // exclusive

            // Check if current position and next position are consecutive
            if (currPos + 3 == nextPos) {
                // Compare curPosMatches and nextPosMatches
                for (size_t curIdx = curPosMatchStart; curIdx < curPosMatchEnd; ++curIdx) {
                    size_t startIdx = linkedMatchValues.size();
                    bool found = false;
                    for (size_t nextIdx = nextPosMatchStart; nextIdx < nextPosMatchEnd; nextIdx++) {
                        if (isConsecutive(matchList + nextIdx, matchList + curIdx)){
                            linkedMatchValues.push_back(matchList + nextIdx);
                            found = true;
                        }
                    }
                    if (found) {
                        linkedMatchKeys[linkedMatchKeysIdx] = matchList + curIdx;
                        linkedMatchValuesIdx[linkedMatchKeysIdx++] = startIdx;
                    }
                }
            }
            // Update curPosMatches and nextPosMatches
            curPosMatchStart = nextPosMatchStart;
            curPosMatchEnd = nextPosMatchEnd;
            currPos = nextPos;
        }
    }

    // Iterate linkedMatches to get filteredMatches 
    // (ignore matches not enoughly consecutive)
    size_t MIN_DEPTH = minConsCnt - 1;
    if (taxonomy->IsAncestor(eukaryotaTaxId, speciesId)) {
        MIN_DEPTH = minConsCntEuk - 1;
    }

    match2depthScore.clear();

    for (size_t k = 0; k < linkedMatchKeysIdx; k++) {
        const Match * key = linkedMatchKeys[k];
        size_t startIdx = linkedMatchValuesIdx[k];
        size_t endIdx = (k + 1 < linkedMatchKeysIdx) ? linkedMatchValuesIdx[k + 1] : linkedMatchValues.size();

        if (match2depthScore.find(key) == match2depthScore.end()) {
            depthScore bestPath{};
            for (size_t j = startIdx; j < endIdx; j++) {
                depthScore curPath = DFS(linkedMatchValues[j],
                                         linkedMatchKeys,
                                         linkedMatchKeysIdx,
                                         linkedMatchValues,
                                         linkedMatchValuesIdx,
                                         1,
                                         MIN_DEPTH, 
                                         match2depthScore,
                                         key->getScore(),
                                         key->hamming);
                if (curPath.score > bestPath.score && curPath.depth > MIN_DEPTH) {
                    bestPath = curPath;
                }
            }
            match2depthScore[key] = bestPath;
            if (bestPath.depth > MIN_DEPTH) {
                matchPaths.emplace_back(key->qInfo.pos, // start coordinate on query
                                        key->qInfo.pos + bestPath.depth * 3 + 20, // end coordinate on query
                                        bestPath.score, bestPath.hammingDist, bestPath.depth,
                                        key,
                                        bestPath.endMatch);
            }
        }
    }
}

depthScore Taxonomer::DFS(
    const Match *curMatch,
    const vector<const Match *> &linkedMatchesKeys,
    size_t linkedMatchKeysIdx,
    const vector<const Match *> &linkedMatchesValues,
    const vector<size_t> &linkedMatchesIndices,
    size_t depth, size_t MIN_DEPTH,
    unordered_map<const Match *, depthScore> &match2depthScore,
    float score, int hammingDist) 
{
    depth++;
    depthScore bestDepthScore{};
    depthScore returnDepthScore;
    depthScore curDepthScore;
    float receivedScore = score;

    auto it = find(linkedMatchesKeys.begin(), linkedMatchesKeys.end(), curMatch);
    if (it == linkedMatchesKeys.end()) { // Reached a leaf node
        uint8_t lastEndHamming = (curMatch->rightEndHamming >> 14);
        if (lastEndHamming == 0) {
            score += 3.0f;
        } else {
            score += 2.0f - 0.5f * lastEndHamming;
        }
        match2depthScore[curMatch] = depthScore(1, score - receivedScore, lastEndHamming, curMatch);
        return depthScore(depth, score, hammingDist + lastEndHamming, curMatch);
    } else { // Not a leaf node
        size_t index = it - linkedMatchesKeys.begin();
        size_t startIdx = linkedMatchesIndices[index];
        size_t endIdx = (index + 1 < linkedMatchKeysIdx) ? linkedMatchesIndices[index + 1] : linkedMatchesValues.size();

        uint8_t lastEndHamming = (curMatch->rightEndHamming >> 14);
        if (lastEndHamming == 0) {
            score += 3.0f;
        } else {
            score += 2.0f - 0.5f * lastEndHamming;
        }
        for (size_t i = startIdx; i < endIdx; ++i) {
            const Match *nextMatch = linkedMatchesValues[i];
            if (match2depthScore.find(nextMatch) != match2depthScore.end()) {
                returnDepthScore = match2depthScore[nextMatch];
                curDepthScore = depthScore(returnDepthScore.depth + depth,
                                           returnDepthScore.score + score,
                                           returnDepthScore.hammingDist + hammingDist + lastEndHamming,
                                           returnDepthScore.endMatch);
            } else {
                curDepthScore = DFS(nextMatch, linkedMatchesKeys, linkedMatchKeysIdx, linkedMatchesValues, linkedMatchesIndices, depth, MIN_DEPTH, match2depthScore, score, hammingDist + lastEndHamming);
            }
            if (curDepthScore.score > bestDepthScore.score && curDepthScore.depth > MIN_DEPTH) {
                bestDepthScore = curDepthScore;
            }
        }
        if (bestDepthScore.depth > MIN_DEPTH) {
            match2depthScore[curMatch] = depthScore(bestDepthScore.depth - depth + 1,
                                                     bestDepthScore.score - receivedScore,
                                                     bestDepthScore.hammingDist - hammingDist,
                                                     bestDepthScore.endMatch);
        }
    }
    return bestDepthScore;
}

bool Taxonomer::isConsecutive(const Match * match1, const Match * match2) {
    // match1 87654321 -> 08765432
    // match2 98765432 -> 08765432
    return (match1->dnaEncoding >> bitsPerCodon) == (match2->dnaEncoding & lastCodonMask);
}

bool Taxonomer::isConsecutive(const Match * match1, const Match * match2, int shift) {
    // match1 ---76543210 -> ---**765432
    // match2 ---98765432 -> ---**765432
    // uint32_t dnaEncoding1 = match1->dnaEncoding >> (3 * shift);
    // uint32_t dnaEncoding2 = match2->dnaEncoding & ((1U << (24 - 3 * shift)) - 1);
    return (match1->dnaEncoding >> (bitsPerCodon * shift)) == (match2->dnaEncoding & ((1U << (totalDnaBits - bitsPerCodon * shift)) - 1));
}

// bool Taxonomer::isConsecutive_diffFrame(const Match * match1, const Match * match2) {
//     // int hamming1 = match1->hamming - GET_2_BITS(match1->rightEndHamming);
//     // int hamming2 = match2->hamming - GET_2_BITS(match2->rightEndHamming >> 14);
//     // match1 87654321 -> 08765432
//     // match2 98765432 -> 08765432
//     return (match1->hamming - GET_2_BITS(match1->rightEndHamming)) == (match2->hamming - GET_2_BITS(match2->rightEndHamming >> 14));
// }


void Taxonomer::ensureArraySize(size_t newSize) {
    if (newSize > arraySize_filterRedundantMatches) {
        delete[] bestMatchForQuotient;
        delete[] bestMatchTaxIdForQuotient;
        delete[] minHammingForQuotient;
        bestMatchForQuotient = new const Match*[newSize]();
        bestMatchTaxIdForQuotient = new TaxID[newSize]();
        minHammingForQuotient = new uint8_t[newSize];
        arraySize_filterRedundantMatches = newSize;
    }
    std::memset(bestMatchForQuotient, 0, newSize * sizeof(const Match*));
    std::memset(minHammingForQuotient, std::numeric_limits<uint8_t>::max(), newSize * sizeof(uint8_t));
}