#include "Taxonomer.h"
#include "BitManipulateMacros.h"
#include "Match.h"
#include "NcbiTaxonomy.h"
#include "common.h"
#include "printBinary.h"
#include <cstdint>
#include <sys/types.h>
#include <unordered_map>


Taxonomer::Taxonomer(const LocalParameters &par, NcbiTaxonomy *taxonomy) : taxonomy(taxonomy) {
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
    minCoveredPos = par.minCoveredPos;
    accessionLevel = par.accessionLevel;
    minSSMatch = par.minSSMatch;
    minConsCnt = par.minConsCnt;
    minConsCntEuk = par.minConsCntEuk;
    eukaryotaTaxId = par.eukaryotaTaxId;
    tieRatio = par.tieRatio;

    if (par.seqMode == 1 || par.seqMode == 2) {
        denominator = 100;
    } else {
        denominator = 1000;
    }


    // chooseBestTaxon
    taxCnt.reserve(4096);

    // getBestSpeciesMatches
    speciesMatches.reserve(4096);
    curFrameMatches.reserve(4096);
    matchPaths.reserve(4096);
    combinedMatchPaths.reserve(4096);
    maxSpecies.reserve(4096);
    speciesList.reserve(4096);
    speciesPathIdx.reserve(4096);
    speciesCombPathIdx.reserve(4096);
    speciesScores.reserve(4096);

    // remainConsecutiveMatches
    curPosMatches.reserve(4096);
    nextPosMatches.reserve(4096);
    linkedMatchKeys.reserve(4096);
    linkedMatchValues.reserve(4096);
    linkedMatchValuesIdx.reserve(4096);
    used.reserve(4096);
    idx2depthScore.reserve(4096); 

    // lowerRankClassification
    cladeCnt.reserve(4096);

    // Output
    taxCounts.reserve(4096);
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
    cout << "Time spent for spliting matches: " << double(time(nullptr) - beforeAnalyze) << endl;
    // Process each block
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
    // Get the best species and its matches for the current query
    if (speciesMatches.size() < end - offset + 1) {
        speciesMatches.reserve(end - offset + 1);
    }
    speciesMatches.clear();

    TaxonScore speciesScore(0, 0, 0, 0, 0);
    speciesScore = getBestSpeciesMatches(speciesMatches,
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
    filterRedundantMatches(speciesMatches, taxCnt);
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

void Taxonomer::filterRedundantMatches(vector<const Match *> & speciesMatches,
                                        unordered_map<TaxID, unsigned int> & taxCnt) {
    // Sort matches by the coordinate on the query
    sort(speciesMatches.begin(), speciesMatches.end(),
         [](const Match * a, const Match * b) { return a->qInfo.pos < b->qInfo.pos; });
    
    // Remove redundant matches
    size_t matchNum = speciesMatches.size();
    for (size_t i = 0; i < matchNum;) {
        size_t currQuotient = speciesMatches[i]->qInfo.pos / 3;
        uint8_t minHamming = speciesMatches[i]->hamming;
        Match minHammingMatch = *speciesMatches[i];
        TaxID minHammingTaxId = minHammingMatch.targetId;
        while ((i < matchNum) && (currQuotient == speciesMatches[i]->qInfo.pos / 3)) {
            if (speciesMatches[i]->hamming < minHamming) {
                minHamming = speciesMatches[i]->hamming;
                minHammingMatch = *speciesMatches[i];
                minHammingTaxId = minHammingMatch.targetId;
            } else if (speciesMatches[i]->hamming == minHamming) {
                minHammingTaxId = taxonomy->LCA(minHammingTaxId, speciesMatches[i]->targetId);
            }
            i++;
        }
        taxCnt[minHammingTaxId]++;
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
            if (strcmp(taxonomy->getString(taxon->rankIdx), "") == 0) {
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

TaxonScore Taxonomer::getBestSpeciesMatches(vector<const Match * > & speciesMatches,
                                            const Match *matchList,
                                            size_t end,
                                            size_t offset,
                                            int queryLength) {
    matchPaths.clear();
    combinedMatchPaths.clear();
    speciesList.clear();
    speciesPathIdx.clear();
    speciesCombPathIdx.clear();
    speciesScores.clear();
    
    pair<size_t, size_t> bestSpeciesRange;
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
        size_t pathSize = matchPaths.size();
        // Combine MatchPaths
        if (pathSize > previousPathSize) {
            speciesList.push_back(currentSpecies);
            speciesPathIdx.push_back(previousPathSize);
            speciesCombPathIdx.push_back(combinedMatchPaths.size());
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
        // matchPaths.clear();
    }
    
    // If there are no meaningful species
    if (meaningfulSpecies == 0) {
        bestScore.score = 0;
        return bestScore;
    }

    maxSpecies.clear();
    float coveredLength = 0.f;
    for (size_t i = 0; i < speciesList.size(); i++) {
        if (speciesScores[i] >= bestSpScore * tieRatio) {
            maxSpecies.push_back(speciesList[i]);
            bestScore.score += speciesScores[i];
            // Calculate coverage
            // size_t start = speciesCombPathIdx[i];
            // size_t end = (i + 1 < speciesCombPathIdx.size()) ? speciesCombPathIdx[i + 1] : combinedMatchPaths.size();
            // coveredLength = 0;
            // for (size_t j = start; j < end; j++) {
            //     coveredLength += combinedMatchPaths[j].end - combinedMatchPaths[j].start + 1;
            // }
            // bestScore.coverage += coveredLength / queryLength;    
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
    speciesMatches.reserve(bestSpeciesRange.second - bestSpeciesRange.first + 1);

    for (size_t j = bestSpeciesRange.first; j < bestSpeciesRange.second; j++) {
        speciesMatches.push_back(&matchList[j]);
    }
    
    return bestScore;                                  
}


float Taxonomer::combineMatchPaths(vector<MatchPath> & matchPaths,
                                   size_t matchPathStart,
                                   vector<MatchPath> & combinedMatchPaths,
                                   size_t combMatchPathStart,
                                   int readLength) {
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
            score += matchPaths[i].score;
        } else {
            bool isOverlapped = false;
            for (size_t j = combMatchPathStart; j < combinedMatchPaths.size(); j++) {
                if (isMatchPathOverlapped(matchPaths[i], combinedMatchPaths[j])) { // overlap!
                    int overlappedLength = min(matchPaths[i].end, combinedMatchPaths[j].end) 
                                            - max(matchPaths[i].start, combinedMatchPaths[j].start) + 1;
                    if (overlappedLength < 24) {
                        trimMatchPath(matchPaths[i], combinedMatchPaths[j], overlappedLength);
                       
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

void Taxonomer::remainConsecutiveMatches(const vector<const Match *> & curFrameMatches,
                                         vector<MatchPath> & matchPaths,
                                         TaxID speciesId) {
    size_t i = 0;
    size_t end = curFrameMatches.size();
    curPosMatches.clear();
    nextPosMatches.clear();
    linkedMatchKeys.clear();
    linkedMatchValues.clear();
    linkedMatchValuesIdx.clear();

    size_t currPos = curFrameMatches[0]->qInfo.pos;
    uint64_t frame = curFrameMatches[0]->qInfo.frame;
    if (frame < 3) { // Forward frame
        while ( i < end && curFrameMatches[i]->qInfo.pos == currPos) {
            curPosMatches.emplace_back(curFrameMatches[i]);
            i++;
        }
        while (i < end) {
            uint32_t nextPos = curFrameMatches[i]->qInfo.pos;
            while (i < end  && nextPos == curFrameMatches[i]->qInfo.pos) {
                nextPosMatches.emplace_back(curFrameMatches[i]);
                ++ i;
            }
            // Check if current position and next position are consecutive
            if (currPos + 3 == nextPos) {
                // Compare curPosMatches and nextPosMatches
                for (auto &curPosMatch: curPosMatches) {
                    size_t startIdx = linkedMatchValues.size();
                    bool found = false;
                    for (auto &nextPosMatch: nextPosMatches) {
                        if (isConsecutive(curPosMatch, nextPosMatch)){
                            linkedMatchValues.push_back(nextPosMatch);
                            found = true;
                            // linkedMatches[curPosMatch].push_back(nextPosMatch);
                        }
                    }
                    if (found) {
                        linkedMatchKeys.push_back(curPosMatch);
                        linkedMatchValuesIdx.push_back(startIdx);
                    }
                }
            }
            // Update curPosMatches and nextPosMatches
            curPosMatches = std::move(nextPosMatches);
            currPos = nextPos;
        }
    } else {
        while ( i < end && curFrameMatches[i]->qInfo.pos == currPos) {
            curPosMatches.emplace_back(curFrameMatches[i]);
            i++;
        }
        while (i < end) {
            uint32_t nextPos = curFrameMatches[i]->qInfo.pos;
            while (i < end  && nextPos == curFrameMatches[i]->qInfo.pos) {
                nextPosMatches.emplace_back(curFrameMatches[i]);
                ++ i;
            }
            // Check if current position and next position are consecutive
            if (currPos + 3 == nextPos) {
                // Compare curPosMatches and nextPosMatches
                for (auto &curPosMatch: curPosMatches) {
                    size_t startIdx = linkedMatchValues.size();
                    bool found = false;
                    for (auto &nextPosMatch: nextPosMatches) {
                        if (isConsecutive(nextPosMatch, curPosMatch)){
                            linkedMatchValues.push_back(nextPosMatch);
                            found = true;
                            // linkedMatches[curPosMatch].push_back(nextPosMatch);
                        }
                    }
                    if (found) {
                        linkedMatchKeys.push_back(curPosMatch);
                        linkedMatchValuesIdx.push_back(startIdx);
                    }
                }

            }
            // Update curPosMatches and nextPosMatches
            curPosMatches = std::move(nextPosMatches);
            currPos = nextPos;
        }
    }

    // Iterate linkedMatches to get filteredMatches 
    // (ignore matches not enoughly consecutive)
    size_t MIN_DEPTH = minConsCnt - 1;
    if (taxonomy->IsAncestor(eukaryotaTaxId, speciesId)) {
        MIN_DEPTH = minConsCntEuk - 1;
    }

    used.clear();
    idx2depthScore.clear();

    for (size_t k = 0; k < linkedMatchKeys.size(); k++) {
        const Match * key = linkedMatchKeys[k];
        size_t startIdx = linkedMatchValuesIdx[k];
        size_t endIdx = (k + 1 < linkedMatchKeys.size()) ? linkedMatchValuesIdx[k + 1] : linkedMatchValues.size();

        if (!used.count(key)) {
            used.insert(key);
            depthScore bestPath{};
            for (size_t j = startIdx; j < endIdx; j++) {
                used.insert(linkedMatchValues[j]);
                depthScore curPath = DFS(curFrameMatches,
                                         linkedMatchValues[j],
                                         linkedMatchKeys,
                                         linkedMatchValues,
                                         linkedMatchValuesIdx,
                                         1,
                                         MIN_DEPTH, 
                                         used, 
                                         idx2depthScore,
                                         key->getScore(),
                                         key->hamming);
                if (curPath.score > bestPath.score && curPath.depth > MIN_DEPTH) {
                    bestPath = curPath;
                }
            }
            // Store the best path
            if (bestPath.depth > MIN_DEPTH) {
                matchPaths.emplace_back(key->qInfo.pos, // start coordinate on query
                                        key->qInfo.pos + bestPath.depth * 3 + 20, // end coordinate on query
                                        bestPath.score, bestPath.hammingDist,
                                        key,
                                        bestPath.endMatch);
            }
        }
    }
}

depthScore Taxonomer::DFS(
    const vector<const Match *> &matches,
    const Match *curMatch,
    const vector<const Match *> &linkedMatchesKeys,
    const vector<const Match *> &linkedMatchesValues,
    const vector<size_t> &linkedMatchesIndices,
    size_t depth, size_t MIN_DEPTH,
    unordered_set<const Match *> &used,
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
        size_t endIdx = (index + 1 < linkedMatchesIndices.size()) ? linkedMatchesIndices[index + 1] : linkedMatchesValues.size();

        uint8_t lastEndHamming = (curMatch->rightEndHamming >> 14);
        if (lastEndHamming == 0) {
            score += 3.0f;
        } else {
            score += 2.0f - 0.5f * lastEndHamming;
        }
        for (size_t i = startIdx; i < endIdx; ++i) {
            const Match *nextMatch = linkedMatchesValues[i];
            used.insert(nextMatch);
            if (match2depthScore.find(nextMatch) != match2depthScore.end()) {
                returnDepthScore = match2depthScore[nextMatch];
                curDepthScore = depthScore(returnDepthScore.depth + depth,
                                           returnDepthScore.score + score,
                                           returnDepthScore.hammingDist + hammingDist + lastEndHamming,
                                           returnDepthScore.endMatch);
            } else {
                curDepthScore = DFS(matches, nextMatch, linkedMatchesKeys, linkedMatchesValues, linkedMatchesIndices, depth, MIN_DEPTH, used, match2depthScore, score, hammingDist + lastEndHamming);
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
    return (match1->dnaEncoding >> 3) == (match2->dnaEncoding & 0x1FFFFF);
}

bool Taxonomer::isConsecutive_diffFrame(const Match * match1, const Match * match2) {
    // int hamming1 = match1->hamming - GET_2_BITS(match1->rightEndHamming);
    // int hamming2 = match2->hamming - GET_2_BITS(match2->rightEndHamming >> 14);
    // match1 87654321 -> 08765432
    // match2 98765432 -> 08765432
    return (match1->hamming - GET_2_BITS(match1->rightEndHamming)) == (match2->hamming - GET_2_BITS(match2->rightEndHamming >> 14));
}