#include "Taxonomer.h"
#include "BitManipulateMacros.h"
#include "Match.h"
#include "NcbiTaxonomy.h"
#include "common.h"
#include "printBinary.h"
#include <cstdint>
#include <sys/types.h>
#include <unordered_map>


Taxonomer::Taxonomer(const LocalParameters &par, TaxonomyWrapper *taxonomy, const MetamerPattern *metamerPattern) :  
    par(par), 
    taxonomy(taxonomy), 
    metamerPattern(metamerPattern), 
    kmerLen(metamerPattern->kmerLen), 
    windowSize(metamerPattern->windowSize),
    windowMask((1U << windowSize) - 1)
{
    if (par.substitutionMatrices.empty()) {
        substitutionMatrix = nullptr;
    } else {
        substitutionMatrix = new SubstitutionMatrix(par.scoringMatrixFile.values.aminoacid().c_str(), 2.0, par.scoreBias);
    }
    // Parameters
    maxGap = par.maxGap;
    accessionLevel = par.accessionLevel;
    minSSMatch = par.minSSMatch;
    minConsCnt = (size_t) par.minConsCnt;
    minConsCntEuk = (size_t) par.minConsCntEuk;
    eukaryotaTaxId = taxonomy->getEukaryotaTaxID();
    tieRatio = par.tieRatio;


    maxCodonShift = par.maxShift;
    dnaShift = 3;
    if (par.syncmer) {
        dnaShift = (kmerLen - par.smerLen) * 3;
    }
    // if (par.syncmer) {
    //     dnaShift = (8 - par.smerLen) * 3;
    //     maxCodonShift = 8 - par.smerLen;
    // }

    if (par.seqMode == 1 || par.seqMode == 2) {
        denominator = 100;
    } else {
        denominator = 1000;
    }

    // chooseBestTaxon
    taxCnt.reserve(4096);

    // getBestSpeciesMatches
    matchPaths.reserve(4096);
    combinedMatchPaths.reserve(4096);
    maxSpecies.reserve(4096);

    // lowerRankClassification
    cladeCnt.reserve(4096);

    // filterRedundantMatches
    arraySize_filterRedundantMatches = 4096;
    bestMatchForQuotient = new const Match*[arraySize_filterRedundantMatches]();
    bestMatchTaxIdForQuotient = new TaxID[arraySize_filterRedundantMatches];
    minHammingForQuotient = new uint8_t[arraySize_filterRedundantMatches];

    // Output
    taxCounts.reserve(4096);

    if (par.maxEValue > 0) {
        logMaxEValue = std::log(par.maxEValue);
        useEvalueFilter = true;
    } else {
        useEvalueFilter = false;
    }
    dbSize = par.dbTotalLength == 0 ? readDbSize(par.filenames[1 + (par.seqMode == 2)]) : par.dbTotalLength;
    // size_t tempSize = 562762599 * 2; // for nr database
}

Taxonomer::~Taxonomer() {
    delete substitutionMatrix;
    delete[] bestMatchForQuotient;
    delete[] bestMatchTaxIdForQuotient;
    delete[] minHammingForQuotient;
    delete evaluer;
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
    TaxonScore speciesScore;
    std::pair<size_t, size_t> bestSpeciesRange;
    speciesScore = getBestSpeciesMatches(bestSpeciesRange,
                                         matchList,
                                         end,
                                         offset,                        
                                         queryList[currentQuery]);
    
    // If there is no proper species for current query, it is un-classified.
    if (speciesScore.score.idScore == 0 || speciesScore.score.idScore < par.minScore) {
        queryList[currentQuery].classification = 0;
        queryList[currentQuery].idScore = speciesScore.score.idScore;
        queryList[currentQuery].subScore = speciesScore.score.subScore;
        queryList[currentQuery].eValue = std::exp(speciesScore.score.logE);
        queryList[currentQuery].hammingDist = speciesScore.hammingDist;
        return;
    }

    // If there are two or more good species level candidates, find the LCA.
    if (speciesScore.LCA) {
        queryList[currentQuery].classification = speciesScore.taxId;
        queryList[currentQuery].idScore = speciesScore.score.idScore;
        queryList[currentQuery].subScore = speciesScore.score.subScore;
        queryList[currentQuery].eValue = std::exp(speciesScore.score.logE);
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
    if (speciesScore.score.idScore < par.minSpScore) {
      queryList[currentQuery].classification = taxonomy->taxonNode(
              taxonomy->getTaxIdAtRank(speciesScore.taxId, "species"))->parentTaxId;
      queryList[currentQuery].idScore = speciesScore.score.idScore;
      queryList[currentQuery].subScore = speciesScore.score.subScore;
      queryList[currentQuery].eValue = std::exp(speciesScore.score.logE);
      queryList[currentQuery].hammingDist = speciesScore.hammingDist;
      return;
    }

    // Store classification results
    queryList[currentQuery].idScore = speciesScore.score.idScore;
    queryList[currentQuery].subScore = speciesScore.score.subScore;
    queryList[currentQuery].eValue = std::exp(speciesScore.score.logE);
    queryList[currentQuery].hammingDist = speciesScore.hammingDist;

    if (!par.em) {
        queryList[currentQuery].classification
            = lowerRankClassification(
                taxCnt,
                speciesScore.taxId,
                queryList[currentQuery].queryLength + queryList[currentQuery].queryLength2);
    } else {
        queryList[currentQuery].classification = speciesScore.taxId;
    }
}


void Taxonomer::filterRedundantMatches(const Match *matchList,
                                       const std::pair<size_t, size_t> & bestSpeciesRange,
                                       unordered_map<TaxID, unsigned int> & taxCnt,
                                       int queryLength) {    
    // Determine the maximum quotient we need to handle
    size_t maxQuotient = (queryLength + 3) / dnaShift;
    
    ensureArraySize(maxQuotient + 1);

    // std::fill_n(bestMatchTaxIdForQuotient, maxQuotient + 1, 0);

    for (size_t i = bestSpeciesRange.first; i < bestSpeciesRange.second; i ++) {
        size_t currQuotient = matchList[i].qKmer.qInfo.pos / dnaShift;
        uint8_t hamming = metamerPattern->hammingDistSum(matchList[i].qKmer.value, matchList[i].tKmer.value, kmerLen, true);

        if (bestMatchForQuotient[currQuotient] == nullptr) {
            bestMatchForQuotient[currQuotient] = matchList + i;
            bestMatchTaxIdForQuotient[currQuotient] = matchList[i].tKmer.tInfo.taxId;
            minHammingForQuotient[currQuotient] = hamming;
        } else {
            if (hamming < minHammingForQuotient[currQuotient]) {
                bestMatchForQuotient[currQuotient] = matchList + i;
                bestMatchTaxIdForQuotient[currQuotient] = matchList[i].tKmer.tInfo.taxId;
                minHammingForQuotient[currQuotient] = hamming;
            } else if (hamming == minHammingForQuotient[currQuotient]) {
                bestMatchTaxIdForQuotient[currQuotient] = taxonomy->LCA(
                    bestMatchTaxIdForQuotient[currQuotient], matchList[i].tKmer.tInfo.taxId);
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
    unsigned int minSubSpeciesMatch = ((queryLength - 1)/denominator);
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
        return BFS(cladeCnt, spTaxId, minSubSpeciesMatch);
    } else {
        return BFS(cladeCnt, spTaxId, minSubSpeciesMatch);
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
                                            Query & query) {
    matchPaths.clear();
    combinedMatchPaths.clear();
    vector<pair<TaxID, MatchScore>> sp2score;
    int queryLength = query.queryLength + query.queryLength2;
    if (par.printLog==1) {
        cout << "## " << query.name << " ##" << endl;
    }
    TaxonScore bestScore;
    MatchScore bestSpScore;
    size_t i = offset;
    size_t meaningfulSpecies = 0;
    if (par.printLog==1) {
        cout << "## " << query.name << " ##" << endl;
        for (size_t j = offset; j < end + 1; j++) {
            metamerPattern->printDNA(matchList[j].qKmer.value);
            cout << " " ;
            metamerPattern->printAA(matchList[j].qKmer.value);
            cout << " | " ;
            metamerPattern->printDNA(matchList[j].tKmer.value);
            cout << " " ;
            metamerPattern->printAA(matchList[j].tKmer.value);
            cout << " | " << taxonomy->getOriginalTaxID(matchList[j].tKmer.tInfo.taxId) << " " << taxonomy->getOriginalTaxID(matchList[j].tKmer.tInfo.speciesId) << " ";
            cout << (int) matchList[j].qKmer.qInfo.frame << " " << (int) matchList[j].qKmer.qInfo.pos << endl;
        }
    }
    while (i  < end + 1) {
        TaxID currentSpecies = matchList[i].tKmer.tInfo.speciesId;
        size_t start = i;
        size_t previousPathSize = matchPaths.size();
        // For current species
        while ((i < end + 1) && currentSpecies == matchList[i].tKmer.tInfo.speciesId) {
            uint8_t curFrame = matchList[i].qKmer.qInfo.frame;
            size_t start = i;
            // For current frame
            while ((i < end + 1) && currentSpecies == matchList[i].tKmer.tInfo.speciesId && curFrame == matchList[i].qKmer.qInfo.frame) {
                i ++;
            }
            if (i - start > 1) {
                if (windowSize == kmerLen) {
                    getMatchPaths(matchList + start, i - start, matchPaths, currentSpecies);
                } else {
                    getMatchPaths2(matchList + start, i - start, matchPaths, currentSpecies);
                }
            }
        }
        size_t pathSize = matchPaths.size();
        // Combine MatchPaths
        if (par.printLog==1) {
            if (pathSize > previousPathSize) {
                cout << "Current species: " << taxonomy->getOriginalTaxID(currentSpecies) << " " << currentSpecies << endl;
                for (size_t kk = previousPathSize; kk < matchPaths.size(); kk++) {
                    matchPaths[kk].printMatchPath();
                }
            }
        }
        if (pathSize > previousPathSize) {
            MatchScore score = combineMatchPaths(matchPaths, previousPathSize, combinedMatchPaths, combinedMatchPaths.size(), queryLength);
            if (par.printLog==1) {   
                cout << "Combined score: " << score.idScore << " " << score.subScore << endl;
            }
            score.idScore = score.idScore / queryLength;
            score.idScore = min(score.idScore, 1.0f);
            if ((score.idScore < par.minScore) || (score.logE > logMaxEValue)) {
                continue;
            }

            sp2score.emplace_back(currentSpecies, score);
            if (score.idScore > 0.f) {
                meaningfulSpecies++;
            }
            if (score.isLargerThan(bestSpScore, par.scoreMode)) {
                bestSpScore = score;
                bestSpeciesRange = make_pair(start, i);
            }
        }
    }
    
    // If there are no meaningful species
    if (meaningfulSpecies == 0) {
        bestScore.score = {0.0f, 0.0f};
        return bestScore;
    }

    // if (par.em && !sp2score.empty()) {
    //     sort(sp2score.begin(), sp2score.end(),
    //          [](const pair<TaxID, float> &a, const pair<TaxID, float> &b) {
    //              return a.second > b.second;
    //          });
    //     query.topSpeciesId = sp2score[0].first;
    //     for (size_t i = 0; i < 10 && i < sp2score.size(); i++) {
    //         query.species2Score.emplace_back(sp2score[i].first, sp2score[i].second * sp2score[i].second);
    //     }
    // }

    maxSpecies.clear();
    for (size_t i = 0; i < sp2score.size(); i++) {
        if (sp2score[i].second.isLargerThan(bestSpScore * tieRatio, par.scoreMode)) {
            maxSpecies.push_back(sp2score[i].first);
            bestScore.score += sp2score[i].second;   
        }
    }
    
    // More than one species --> LCA    
    if (maxSpecies.size() > 1) {
        bestScore.LCA = true;
        bestScore.taxId = taxonomy->LCA(maxSpecies)->taxId;
        bestScore.score.idScore /= maxSpecies.size();
        bestScore.score.subScore /= maxSpecies.size();
        bestScore.score.logE = bestSpScore.logE;
        return bestScore;
    }
    
    // One species
    bestScore.taxId = maxSpecies[0];
    bestScore.score.logE = bestSpScore.logE;
    
    return bestScore;                                  
}

void Taxonomer::sortMatchPath(std::vector<MatchPath> & matchPaths, size_t i) {
    if (par.scoreMode == 0) {
        sort(matchPaths.begin() + i, matchPaths.end(),
         [](const MatchPath &a, const MatchPath &b) {
           if (a.score.idScore != b.score.idScore) {
             return a.score.idScore > b.score.idScore;
           }
           if (a.hammingDist != b.hammingDist) {
             return a.hammingDist < b.hammingDist;
           }
           return a.start > b.start;
         });
    } else if (par.scoreMode == 1) {
        sort(matchPaths.begin() + i, matchPaths.end(),
         [](const MatchPath &a, const MatchPath &b) {
           if (a.score.subScore != b.score.subScore) {
             return a.score.subScore > b.score.subScore;
           }
           if (a.hammingDist != b.hammingDist) {
             return a.hammingDist < b.hammingDist;
           }
           return a.start > b.start;
         });
    } else {
        sort(matchPaths.begin() + i, matchPaths.end(),
         [](const MatchPath &a, const MatchPath &b) {
           float aTotal = a.score.idScore + a.score.subScore;
           float bTotal = b.score.idScore + b.score.subScore;
           if (aTotal != bTotal) {
             return aTotal > bTotal;
           }
           if (a.hammingDist != b.hammingDist) {
             return a.hammingDist < b.hammingDist;
           }
           return a.start > b.start;
         });
    }
}

MatchScore Taxonomer::combineMatchPaths(
    vector<MatchPath> & matchPaths,
    size_t matchPathStart,
    vector<MatchPath> & combinedMatchPaths,
    size_t combMatchPathStart,
    int queryLength) 
{   
    // Combine matchPaths
    // 1. Add the matchPath with the highest score to combinedMatchPaths
    // 2. Add the matchPath with the highest score that is not overlapped with the matchPath in combinedMatchPaths
    // 3. Repeat 2 until no matchPath can be added
    sortMatchPath(matchPaths, matchPathStart);
    MatchScore score;
    int spanLength = 0;
    for (size_t i = matchPathStart; i < matchPaths.size(); i++) {  
        if (combMatchPathStart == combinedMatchPaths.size()) {
            combinedMatchPaths.push_back(matchPaths[i]);
            score += matchPaths[i].score;
            spanLength += matchPaths[i].end - matchPaths[i].start + 1;

            // score.logE = computeLEM_logE(
            //     queryLength,
            //     matchPaths[i].end - matchPaths[i].start + 1,
            //     dbSize,
            //     score.logP);
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
                        continue;
                    } else {
                        isOverlapped = true;
                        break;
                    }
                    if (matchPaths[i].end - matchPaths[i].start + 1 < windowSize * 3) { // Current path trimmed too much 
                        isOverlapped = true;
                        break;
                    }  
                } 
            }
            if (!isOverlapped) {
                combinedMatchPaths.push_back(matchPaths[i]);
                score += matchPaths[i].score;     
                spanLength += matchPaths[i].end - matchPaths[i].start + 1;      
            }
        }
    }

    score.logE = computeLEM_logE(
        queryLength,
        spanLength,
        dbSize,
        score.logP);

    return score;
}

bool Taxonomer::isMatchPathOverlapped(const MatchPath & matchPath1,
                                      const MatchPath & matchPath2) {
    return !((matchPath1.end < matchPath2.start) || (matchPath2.end < matchPath1.start));                                       
}

void Taxonomer::trimMatchPath(
    MatchPath & newPath, 
    const MatchPath & path2, 
    int overlapLength) 
{
    if (newPath.start < path2.start) { 
        newPath.end = path2.start - 1;
        if (newPath.endMatch->qKmer.qInfo.frame < 3) {
            newPath.hammingDist = max(0, newPath.hammingDist - metamerPattern->hammingDistSum(
                newPath.endMatch->qKmer.value,
                newPath.endMatch->tKmer.value,
                overlapLength/3, // (overlapLength + 2) / 3
                true));

            newPath.score -= metamerPattern->calMatchScore(
                newPath.endMatch->qKmer.value,
                newPath.endMatch->tKmer.value,
                overlapLength/3,
                *substitutionMatrix,
                true);

            newPath.score.idScore -= overlapLength % 3;
        } else {
            newPath.hammingDist = max(0, newPath.hammingDist - metamerPattern->hammingDistSum(
                newPath.endMatch->qKmer.value,
                newPath.endMatch->tKmer.value, 
                overlapLength/3,
                false));
            
            newPath.score -= metamerPattern->calMatchScore(
                newPath.endMatch->qKmer.value,
                newPath.endMatch->tKmer.value,
                overlapLength/3,
                *substitutionMatrix,
                false);

            newPath.score.idScore -= overlapLength % 3;
        }
    } else {
        newPath.start = path2.end + 1;
        if (newPath.startMatch->qKmer.qInfo.frame < 3) {
            newPath.hammingDist = max(0, newPath.hammingDist - metamerPattern->hammingDistSum(
                newPath.startMatch->qKmer.value,
                newPath.startMatch->tKmer.value, 
                overlapLength/3,
                false));
            
            newPath.score -= metamerPattern->calMatchScore(
                newPath.startMatch->qKmer.value,
                newPath.startMatch->tKmer.value,
                overlapLength/3,
                *substitutionMatrix,
                false);
            newPath.score.idScore -= overlapLength % 3;
        } else {
            newPath.hammingDist = max(0, newPath.hammingDist - metamerPattern->hammingDistSum(
                newPath.startMatch->qKmer.value,
                newPath.startMatch->tKmer.value, 
                overlapLength/3,
                true));

            newPath.score -= metamerPattern->calMatchScore(
                newPath.startMatch->qKmer.value,
                newPath.startMatch->tKmer.value,
                overlapLength/3,
                *substitutionMatrix,
                true);
            newPath.score.idScore -= overlapLength % 3;
        }
    }
}

// void Taxonomer::trimMatchPath2(
//     MatchPath & newPath, 
//     const MatchPath & path2, 
//     int overlapLength) 
// {
//     const int shift = overlapLength / 3;
//     if (newPath.start < path2.start) { 
//         newPath.end = path2.start - 1;
//         if (newPath.endMatch->qKmer.qInfo.frame < 3) {
//             newPath.hammingDist = max(0, newPath.hammingDist - metamerPattern->hammingDistSum(
//                 newPath.endMatch->qKmer.value,
//                 newPath.endMatch->tKmer.value,
//                 overlapLength/3,
//                 true));
            
//             newPath.score -= metamerPattern->calMatchScore(
//                 newPath.endMatch->qKmer.value,
//                 newPath.endMatch->tKmer.value,
//                 newPath.historyMask & ((1U << shift) - 1),
//                 *substitutionMatrix);    

//             newPath.score.idScore -= overlapLength % 3;
//         } else {
//             newPath.hammingDist = max(0, newPath.hammingDist - metamerPattern->hammingDistSum(
//                 newPath.endMatch->qKmer.value,
//                 newPath.endMatch->tKmer.value, 
//                 overlapLength/3,
//                 false));
            
//             newPath.score -= metamerPattern->calMatchScore(
//                 newPath.endMatch->qKmer.value,
//                 newPath.endMatch->tKmer.value,
//                 newPath.historyMask & (~0U << (32 - shift)),
//                 *substitutionMatrix);    

//             newPath.score.idScore -= overlapLength % 3;
//         }
//     } else {
//         newPath.start = path2.end + 1;
//         if (newPath.startMatch->qKmer.qInfo.frame < 3) {
//             newPath.hammingDist = max(0, newPath.hammingDist - metamerPattern->hammingDistSum(
//                 newPath.startMatch->qKmer.value,
//                 newPath.startMatch->tKmer.value, 
//                 overlapLength/3,
//                 false));
            
//             newPath.score -= metamerPattern->calMatchScore(
//                 newPath.startMatch->qKmer.value,
//                 newPath.startMatch->tKmer.value,
//                 newPath.firstHistoryMask & (~0U << (32 - shift)),
//                 *substitutionMatrix); 
//             newPath.score.idScore -= overlapLength % 3;
//         } else {
//             newPath.hammingDist = max(0, newPath.hammingDist - metamerPattern->hammingDistSum(
//                 newPath.startMatch->qKmer.value,
//                 newPath.startMatch->tKmer.value, 
//                 overlapLength/3,
//                 true));

//             newPath.score -= metamerPattern->calMatchScore(
//                 newPath.startMatch->qKmer.value,
//                 newPath.startMatch->tKmer.value,
//                 newPath.firstHistoryMask & ((1U << shift) - 1),
//                 *substitutionMatrix);
//             newPath.score.idScore -= overlapLength % 3;
//         }
//     }
// }


MatchPath Taxonomer::makeMatchPath(
    const Match * match) 
{
    return MatchPath(
        match, 
        metamerPattern->calMatchScore(
            match->qKmer.value, 
            match->tKmer.value, 
            metamerPattern->windowSize, 
            *substitutionMatrix,
            true),
        metamerPattern->hammingDistSum(
            match->qKmer.value, 
            match->tKmer.value),
        metamerPattern->windowSize * 3);
}

void Taxonomer::getMatchPaths(
    const Match * matchList,
    size_t matchNum,
    vector<MatchPath> & filteredMatchPaths,
    TaxID speciesId) 
{
    size_t i = 0;
    size_t currPos = matchList[0].qKmer.qInfo.pos;  
    uint64_t frame = matchList[0].qKmer.qInfo.frame;
    int MIN_DEPTH = (int) minConsCnt;
    if (taxonomy->IsAncestor(eukaryotaTaxId, speciesId)) {
        MIN_DEPTH = (int) minConsCntEuk;
    }
        
    connectedToNext.resize(matchNum);
    fill(connectedToNext.begin(), connectedToNext.end(), false);
    localMatchPaths.clear();
    localMatchPaths.resize(matchNum);

    if (frame < 3) { // Forward frame
        size_t curPosMatchStart = i;        
        while (i < matchNum && matchList[i].qKmer.qInfo.pos == currPos) {
            localMatchPaths[i] = makeMatchPath(matchList + i);
            ++ i;
        }
        size_t curPosMatchEnd = i; // exclusive

        while (i < matchNum) {
            uint32_t nextPos = matchList[i].qKmer.qInfo.pos;
            size_t nextPosMatchStart = i;
            while (i < matchNum  && nextPos == matchList[i].qKmer.qInfo.pos) {
                localMatchPaths[i] = makeMatchPath(matchList + i);
                ++ i;
            }
            size_t nextPosMatchEnd = i; // exclusive

            // Check if current position and next position are consecutive
            int shift = (nextPos - currPos) / 3;
            if (shift > 0 && shift <= maxCodonShift) {
                // Compare curPosMatches and nextPosMatches.
                for (size_t nextIdx = nextPosMatchStart; nextIdx < nextPosMatchEnd; nextIdx++) {
                    const MatchPath * bestPath = nullptr;
                    MatchScore bestScore;
                    for (size_t curIdx = curPosMatchStart; curIdx < curPosMatchEnd; ++curIdx) {
                        if (metamerPattern->checkOverlap(matchList[curIdx].tKmer.value, matchList[nextIdx].tKmer.value, shift)) {
                            connectedToNext[curIdx] = true;
                            if (localMatchPaths[curIdx].score.isLargerThan(bestScore, par.scoreMode)) {
                                bestPath = &localMatchPaths[curIdx];
                                bestScore = localMatchPaths[curIdx].score;
                            }
                        }
                    }
                    if (bestPath != nullptr) {
                        localMatchPaths[nextIdx].start = bestPath->start;                        
                        localMatchPaths[nextIdx].score = bestPath->score +
                            metamerPattern->calMatchScore(
                                matchList[nextIdx].qKmer.value, 
                                matchList[nextIdx].tKmer.value, 
                                shift,
                                *substitutionMatrix,
                                true);

                        localMatchPaths[nextIdx].hammingDist = bestPath->hammingDist + 
                            metamerPattern->hammingDistSum(
                                matchList[nextIdx].qKmer.value, 
                                matchList[nextIdx].tKmer.value, 
                                shift,
                                true);
                        localMatchPaths[nextIdx].depth = bestPath->depth + shift;
                        localMatchPaths[nextIdx].startMatch = bestPath->startMatch;
                    }
                }
            } 
            for (size_t curIdx = curPosMatchStart; curIdx < curPosMatchEnd; ++curIdx) {
                if (!connectedToNext[curIdx] && localMatchPaths[curIdx].depth >= MIN_DEPTH) {
                    filteredMatchPaths.push_back(localMatchPaths[curIdx]);
                }
            }
            if (i == matchNum) {
                for (size_t nextIdx = nextPosMatchStart; nextIdx < nextPosMatchEnd; ++nextIdx) {
                    if (localMatchPaths[nextIdx].depth >= MIN_DEPTH) {
                        filteredMatchPaths.push_back(localMatchPaths[nextIdx]);
                    }
                }
            }
            curPosMatchStart = nextPosMatchStart;
            curPosMatchEnd = nextPosMatchEnd;
            currPos = nextPos;    
        }
    } else {
        size_t curPosMatchStart = i;        
        while (i < matchNum && matchList[i].qKmer.qInfo.pos == currPos) {
            localMatchPaths[i] = makeMatchPath(matchList + i);
            ++ i;
        }
        size_t curPosMatchEnd = i; // exclusive

        while (i < matchNum) {
            uint32_t nextPos = matchList[i].qKmer.qInfo.pos;
            size_t nextPosMatchStart = i;
            while (i < matchNum  && nextPos == matchList[i].qKmer.qInfo.pos) {
                localMatchPaths[i] = makeMatchPath(matchList + i);
                ++ i;
            }
            size_t nextPosMatchEnd = i; // exclusive

            // Check if current position and next position are consecutive
            int shift = (nextPos - currPos) / 3;
            if (shift > 0 && shift <= maxCodonShift) {
                for (size_t nextIdx = nextPosMatchStart; nextIdx < nextPosMatchEnd; nextIdx++) {
                    const MatchPath * bestPath = nullptr;
                    MatchScore bestScore;
                    for (size_t curIdx = curPosMatchStart; curIdx < curPosMatchEnd; ++curIdx) {
                        if (metamerPattern->checkOverlap(matchList[nextIdx].tKmer.value, matchList[curIdx].tKmer.value, shift)) {
                            connectedToNext[curIdx] = true;
                            if (localMatchPaths[curIdx].score.isLargerThan(bestScore, par.scoreMode)) {
                                bestPath = &localMatchPaths[curIdx];
                                bestScore = localMatchPaths[curIdx].score;
                            }
                        }
                    }
                    if (bestPath != nullptr) {
                        localMatchPaths[nextIdx].start = bestPath->start;                        
                        localMatchPaths[nextIdx].score = bestPath->score +
                            metamerPattern->calMatchScore(
                                matchList[nextIdx].qKmer.value, 
                                matchList[nextIdx].tKmer.value, 
                                shift,
                                *substitutionMatrix,
                                false);

                        localMatchPaths[nextIdx].hammingDist = bestPath->hammingDist + 
                            metamerPattern->hammingDistSum(
                                matchList[nextIdx].qKmer.value, 
                                matchList[nextIdx].tKmer.value, 
                                shift,
                                false);
                        localMatchPaths[nextIdx].depth = bestPath->depth + shift;
                        localMatchPaths[nextIdx].startMatch = bestPath->startMatch;
                    }
                }
            } 
            for (size_t curIdx = curPosMatchStart; curIdx < curPosMatchEnd; ++curIdx) {
                if (!connectedToNext[curIdx] && localMatchPaths[curIdx].depth >= MIN_DEPTH) {
                    filteredMatchPaths.push_back(localMatchPaths[curIdx]);
                }
            }
            if (i == matchNum) {
                for (size_t nextIdx = nextPosMatchStart; nextIdx < nextPosMatchEnd; ++nextIdx) {
                    if (localMatchPaths[nextIdx].depth >= MIN_DEPTH) {
                        filteredMatchPaths.push_back(localMatchPaths[nextIdx]);
                    }
                }
            }
            curPosMatchStart = nextPosMatchStart;
            curPosMatchEnd = nextPosMatchEnd;
            currPos = nextPos;    
        }
    }
}

void Taxonomer::makeMatchPath(
    const Match * match,
    size_t index)
{
    localMatchPaths[index] = MatchPath(match, metamerPattern->windowSize * 3);
    localMatchPaths[index].score = metamerPattern->calMatchScore(
                                        match->qKmer.value, 
                                        match->tKmer.value,
                                        *substitutionMatrix);

    localMatchPaths[index].hammingDist = metamerPattern->hammingDistSum(
                                            match->qKmer.value, 
                                            match->tKmer.value);
    localMatchPaths[index].historyMask = metamerPattern->spaceMask;
    localMatchPaths[index].firstHistoryMask = localMatchPaths[index].historyMask;
}

void Taxonomer::getMatchPaths2(
    const Match * matchList,
    size_t matchNum,
    vector<MatchPath> & filteredMatchPaths,
    TaxID speciesId) 
{
    size_t i = 0;
    size_t currPos = matchList[0].qKmer.qInfo.pos;  
    uint64_t frame = matchList[0].qKmer.qInfo.frame;
    int MIN_DEPTH = (int) minConsCnt;
    if (taxonomy->IsAncestor(eukaryotaTaxId, speciesId)) {
        MIN_DEPTH = (int) minConsCntEuk;
    }
        
    connectedToNext.resize(matchNum);
    fill(connectedToNext.begin(), connectedToNext.end(), false);
    localMatchPaths.clear();
    localMatchPaths.resize(matchNum);

    if (frame < 3) { // Forward frame
        size_t curPosMatchStart = i;        
        while (i < matchNum && matchList[i].qKmer.qInfo.pos == currPos) {
            makeMatchPath(matchList + i, i);
            ++ i;
        }
        size_t curPosMatchEnd = i; // exclusive

        while (i < matchNum) {
            uint32_t nextPos = matchList[i].qKmer.qInfo.pos;
            size_t nextPosMatchStart = i;
            while (i < matchNum  && nextPos == matchList[i].qKmer.qInfo.pos) {
                makeMatchPath(matchList + i, i);
                ++ i;
            }
            size_t nextPosMatchEnd = i; // exclusive

            // Check if current position and next position are consecutive
            int shift = (nextPos - currPos) / 3;
            if (shift > 0 && shift <= maxCodonShift) {
                // Compare curPosMatches and nextPosMatches.
                for (size_t nextIdx = nextPosMatchStart; nextIdx < nextPosMatchEnd; nextIdx++) {
                    const MatchPath * bestPath = nullptr;
                    MatchScore bestScore;
                    for (size_t curIdx = curPosMatchStart; curIdx < curPosMatchEnd; ++curIdx) {                    
                        if (metamerPattern->checkOverlap(matchList[curIdx].tKmer.value, matchList[nextIdx].tKmer.value, shift)) {
                            connectedToNext[curIdx] = true;
                            const uint32_t shiftedHistoryMask = (localMatchPaths[curIdx].historyMask << shift) & windowMask;
                            const uint32_t validPosMask = metamerPattern->spaceMask & (~shiftedHistoryMask);
                            const MatchScore totalScore = localMatchPaths[curIdx].score +
                                                            metamerPattern->calMatchScore(
                                                                matchList[nextIdx].qKmer.value, 
                                                                matchList[nextIdx].tKmer.value, 
                                                                validPosMask,
                                                                *substitutionMatrix);

                            if (totalScore.isLargerThan(bestScore, par.scoreMode)) {
                                bestPath = &localMatchPaths[curIdx];
                                bestScore = totalScore;
                                localMatchPaths[nextIdx].historyMask = shiftedHistoryMask | metamerPattern->spaceMask;
                            }
                        }
                    }
                    if (bestPath != nullptr) {
                        localMatchPaths[nextIdx].start = bestPath->start;                        
                        localMatchPaths[nextIdx].score = bestScore;
                        localMatchPaths[nextIdx].hammingDist = bestPath->hammingDist + 
                            metamerPattern->hammingDistSum(
                                matchList[nextIdx].qKmer.value, 
                                matchList[nextIdx].tKmer.value, 
                                shift,
                                true);
                        localMatchPaths[nextIdx].depth = bestPath->depth + shift;
                        localMatchPaths[nextIdx].startMatch = bestPath->startMatch;
                        localMatchPaths[nextIdx].firstHistoryMask = 
                            bestPath->firstHistoryMask | safe_right_shift_32(metamerPattern->spaceMask, (localMatchPaths[nextIdx].depth - 1));
                    }
                }
            } 
            for (size_t curIdx = curPosMatchStart; curIdx < curPosMatchEnd; ++curIdx) {
                if (!connectedToNext[curIdx] && localMatchPaths[curIdx].depth >= MIN_DEPTH) {
                    filteredMatchPaths.push_back(localMatchPaths[curIdx]);
                }
            }
            if (i == matchNum) {
                for (size_t nextIdx = nextPosMatchStart; nextIdx < nextPosMatchEnd; ++nextIdx) {
                    if (localMatchPaths[nextIdx].depth >= MIN_DEPTH) {
                        filteredMatchPaths.push_back(localMatchPaths[nextIdx]);
                    }
                }
            }
            curPosMatchStart = nextPosMatchStart;
            curPosMatchEnd = nextPosMatchEnd;
            currPos = nextPos;    
        }
    } else {
        size_t curPosMatchStart = i;        
        while (i < matchNum && matchList[i].qKmer.qInfo.pos == currPos) {
            makeMatchPath(matchList + i, i);
            ++ i;
        }
        size_t curPosMatchEnd = i; // exclusive

        while (i < matchNum) {
            uint32_t nextPos = matchList[i].qKmer.qInfo.pos;
            size_t nextPosMatchStart = i;
            while (i < matchNum  && nextPos == matchList[i].qKmer.qInfo.pos) {
                makeMatchPath(matchList + i, i);
                ++ i;
            }
            size_t nextPosMatchEnd = i; // exclusive

            // Check if current position and next position are consecutive
            int shift = (nextPos - currPos) / 3;
            if (shift > 0 && shift <= maxCodonShift) {
                for (size_t nextIdx = nextPosMatchStart; nextIdx < nextPosMatchEnd; nextIdx++) {
                    const MatchPath * bestPath = nullptr;
                    MatchScore bestScore;
                    for (size_t curIdx = curPosMatchStart; curIdx < curPosMatchEnd; ++curIdx) {
                        if (metamerPattern->checkOverlap(matchList[nextIdx].tKmer.value, matchList[curIdx].tKmer.value, shift)) {
                            connectedToNext[curIdx] = true;
                            const uint32_t shiftedHistoryMask = (localMatchPaths[curIdx].historyMask >> shift); 
                            const uint32_t validPosMask = metamerPattern->spaceMask & (~shiftedHistoryMask);
                            const MatchScore totalScore = 
                                localMatchPaths[curIdx].score + 
                                metamerPattern->calMatchScore(
                                    matchList[nextIdx].qKmer.value, 
                                    matchList[nextIdx].tKmer.value, 
                                    validPosMask,
                                    *substitutionMatrix);

                            if (totalScore.isLargerThan(bestScore, par.scoreMode)) {
                                bestPath = &localMatchPaths[curIdx];
                                bestScore = totalScore;
                                localMatchPaths[nextIdx].historyMask = shiftedHistoryMask | metamerPattern->spaceMask;
                            }
                        }
                    }
                    if (bestPath != nullptr) {
                        localMatchPaths[nextIdx].start = bestPath->start;                        
                        localMatchPaths[nextIdx].score = bestScore;
                        localMatchPaths[nextIdx].hammingDist = bestPath->hammingDist + 
                            metamerPattern->hammingDistSum(
                                matchList[nextIdx].qKmer.value, 
                                matchList[nextIdx].tKmer.value, 
                                shift,
                                false);
                        localMatchPaths[nextIdx].depth = bestPath->depth + shift;
                        localMatchPaths[nextIdx].startMatch = bestPath->startMatch;
                        localMatchPaths[nextIdx].firstHistoryMask = 
                            bestPath->firstHistoryMask | safe_left_shift_32(metamerPattern->spaceMask, (localMatchPaths[nextIdx].depth - 1));
                        localMatchPaths[nextIdx].firstHistoryMask &= windowMask;
                    }
                }
            } 
            for (size_t curIdx = curPosMatchStart; curIdx < curPosMatchEnd; ++curIdx) {
                if (!connectedToNext[curIdx] && localMatchPaths[curIdx].depth >= MIN_DEPTH) {
                    filteredMatchPaths.push_back(localMatchPaths[curIdx]);
                }
            }
            if (i == matchNum) {
                for (size_t nextIdx = nextPosMatchStart; nextIdx < nextPosMatchEnd; ++nextIdx) {
                    if (localMatchPaths[nextIdx].depth >= MIN_DEPTH) {
                        filteredMatchPaths.push_back(localMatchPaths[nextIdx]);
                    }
                }
            }
            curPosMatchStart = nextPosMatchStart;
            curPosMatchEnd = nextPosMatchEnd;
            currPos = nextPos;    
        }
    }
}

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


void Taxonomer::getMatchPaths3(
    const Match * matchList,
    size_t matchNum,
    vector<MatchPath> & filteredMatchPaths,
    TaxID speciesId) 
{
    if (matchNum == 0) return;

    int MIN_DEPTH = (int) minConsCnt;
    if (taxonomy->IsAncestor(eukaryotaTaxId, speciesId)) {
        MIN_DEPTH = (int) minConsCntEuk;
    }
    
    // --- Step 1: Pre-calculate Groups (Start/End indices for each Position) ---
    // This replaces the complex 'while' loops with a clean structure we can jump around in.
    struct MatchGroup {
        uint32_t pos;
        size_t start;
        size_t end; // exclusive
    };
    std::vector<MatchGroup> groups;

    // Reset DP states
    connectedToNext.resize(matchNum);
    fill(connectedToNext.begin(), connectedToNext.end(), false);
    localMatchPaths.clear();
    localMatchPaths.resize(matchNum);

    // Scan input to build groups
    for (size_t i = 0; i < matchNum; ) {
        uint32_t currentPos = matchList[i].qKmer.qInfo.pos;
        size_t groupStart = i;

        // Process all matches with the same position
        while (i < matchNum && matchList[i].qKmer.qInfo.pos == currentPos) {
            makeMatchPath(matchList + i, i);
            i++; 
        }

        // i now represents the exclusive end of this group
        groups.push_back({currentPos, groupStart, i});
    }

    // Only process Forward Frame as per original code logic
    // (Note: You might want to handle reverse frame similarly if needed)
    if (matchList[0].qKmer.qInfo.frame < 3) { 
        
        // --- Step 2: Flexible Chaining Loop ---
        
        // Outer Loop: Iterate through every group acting as the "Destination"
        for (size_t destG = 0; destG < groups.size(); ++destG) {
            
            const MatchGroup& destGroup = groups[destG];

            // For every match in the Destination Group...
            for (size_t nextIdx = destGroup.start; nextIdx < destGroup.end; ++nextIdx) {
                
                const MatchPath * bestPath = nullptr;
                MatchScore bestScore; // Defaults to -inf or 0
                bool foundParent = false;
                int winningShift = 0;

                // Inner Loop: Look BACKWARDS at previous groups ("Sources")
                // We check as many groups as allowed by maxCodonShift
                for (int srcG = (int)destG - 1; srcG >= 0; --srcG) {
                    
                    const MatchGroup& srcGroup = groups[srcG];
                    
                    int shift = (destGroup.pos - srcGroup.pos) / 3;

                    // If we exceeded the max gap, stop looking further back
                    if (shift > maxCodonShift) break; 

                    bool connectedToThisGroup = false;

                    // Compare Destination Match vs All Candidates in this Source Group
                    for (size_t curIdx = srcGroup.start; curIdx < srcGroup.end; ++curIdx) {
                        
                        // Check geometrical overlap
                        if (metamerPattern->checkOverlap(matchList[curIdx].tKmer.value, matchList[nextIdx].tKmer.value, shift)) {
                            
                            // Mark the source as having a valid successor (it is not a tail)
                            connectedToNext[curIdx] = true; 
                            connectedToThisGroup = true;

                            const uint32_t shiftedHistoryMask = (localMatchPaths[curIdx].historyMask << shift) & windowMask;
                            const uint32_t validPosMask = metamerPattern->spaceMask & (~shiftedHistoryMask);
                            
                            const MatchScore totalScore = localMatchPaths[curIdx].score +
                                                            metamerPattern->calMatchScore(
                                                                matchList[nextIdx].qKmer.value, 
                                                                matchList[nextIdx].tKmer.value, 
                                                                validPosMask,
                                                                *substitutionMatrix);

                            // We want the BEST score across ALL checked source groups
                            if (!foundParent || totalScore.isLargerThan(bestScore, par.scoreMode)) {
                                bestPath = &localMatchPaths[curIdx];
                                bestScore = totalScore;
                                localMatchPaths[nextIdx].historyMask = shiftedHistoryMask | metamerPattern->spaceMask;
                                winningShift = shift;
                                foundParent = true;
                            }
                        }
                    }
                    if (connectedToThisGroup) {
                           break; 
                    }
                } // End Source Group Loop

                // If we found a valid parent among any of the previous allowable groups, update the state
                if (foundParent && bestPath != nullptr) {
                    localMatchPaths[nextIdx].start = bestPath->start;                        
                    localMatchPaths[nextIdx].score = bestScore;
                    localMatchPaths[nextIdx].hammingDist = bestPath->hammingDist + 
                        metamerPattern->hammingDistSum(
                            matchList[nextIdx].qKmer.value, 
                            matchList[nextIdx].tKmer.value, 
                            winningShift,
                            true);
                    localMatchPaths[nextIdx].depth = bestPath->depth + winningShift;
                    localMatchPaths[nextIdx].startMatch = bestPath->startMatch;
                    localMatchPaths[nextIdx].firstHistoryMask = 
                        bestPath->firstHistoryMask | safe_right_shift_32(metamerPattern->spaceMask, (localMatchPaths[nextIdx].depth - 1));
                 }
            }
        }

        // --- Step 3: Collect Results ---
        // Any match that didn't connect to a next one (is a "tail") and is deep enough is a valid path.
        for (size_t i = 0; i < matchNum; ++i) {
            if (!connectedToNext[i] && localMatchPaths[i].depth >= MIN_DEPTH) {
                filteredMatchPaths.push_back(localMatchPaths[i]);
            }
        }
    } else {
        for (size_t destG = 0; destG < groups.size(); ++destG) {
            
            const MatchGroup& destGroup = groups[destG];

            // For every match in the Destination Group...
            for (size_t nextIdx = destGroup.start; nextIdx < destGroup.end; ++nextIdx) {
                
                const MatchPath * bestPath = nullptr;
                MatchScore bestScore; // Defaults to -inf or 0
                bool foundParent = false;
                int winningShift = 0;

                // Inner Loop: Look BACKWARDS at previous groups ("Sources")
                // We check as many groups as allowed by maxCodonShift
                for (int srcG = (int)destG - 1; srcG >= 0; --srcG) {
                    
                    const MatchGroup& srcGroup = groups[srcG];
                    
                    int shift = (destGroup.pos - srcGroup.pos) / 3;

                    // If we exceeded the max gap, stop looking further back
                    if (shift > maxCodonShift) break; 

                    bool connectedToThisGroup = false;

                    // Compare Destination Match vs All Candidates in this Source Group
                    for (size_t curIdx = srcGroup.start; curIdx < srcGroup.end; ++curIdx) {
                        
                        // Check geometrical overlap
                        if (metamerPattern->checkOverlap(matchList[nextIdx].tKmer.value, matchList[curIdx].tKmer.value, shift)) {
                            
                            // Mark the source as having a valid successor (it is not a tail)
                            connectedToNext[curIdx] = true; 
                            connectedToThisGroup = true;

                            const uint32_t shiftedHistoryMask = (localMatchPaths[curIdx].historyMask >> shift); 
                            const uint32_t validPosMask = metamerPattern->spaceMask & (~shiftedHistoryMask);
                            
                            const MatchScore totalScore = localMatchPaths[curIdx].score +
                                                            metamerPattern->calMatchScore(
                                                                matchList[nextIdx].qKmer.value, 
                                                                matchList[nextIdx].tKmer.value, 
                                                                validPosMask,
                                                                *substitutionMatrix);

                            // We want the BEST score across ALL checked source groups
                            if (!foundParent || totalScore.isLargerThan(bestScore, par.scoreMode)) {
                                bestPath = &localMatchPaths[curIdx];
                                bestScore = totalScore;
                                localMatchPaths[nextIdx].historyMask = shiftedHistoryMask | metamerPattern->spaceMask;
                                winningShift = shift;
                                foundParent = true;
                            }
                        }
                    }

                    if (connectedToThisGroup) {
                       break; 
                    }
                } // End Source Group Loop

                // If we found a valid parent among any of the previous allowable groups, update the state
                if (foundParent && bestPath != nullptr) {
                    localMatchPaths[nextIdx].start = bestPath->start;                        
                    localMatchPaths[nextIdx].score = bestScore;
                    localMatchPaths[nextIdx].hammingDist = bestPath->hammingDist + 
                        metamerPattern->hammingDistSum(
                            matchList[nextIdx].qKmer.value, 
                            matchList[nextIdx].tKmer.value, 
                            winningShift,
                            false);
                    localMatchPaths[nextIdx].depth = bestPath->depth + winningShift;
                    localMatchPaths[nextIdx].startMatch = bestPath->startMatch;
                    localMatchPaths[nextIdx].firstHistoryMask = 
                        bestPath->firstHistoryMask | safe_left_shift_32(metamerPattern->spaceMask, (localMatchPaths[nextIdx].depth - 1));
                    localMatchPaths[nextIdx].firstHistoryMask &= windowMask;
                 }
            }
        }

        // --- Step 3: Collect Results ---
        // Any match that didn't connect to a next one (is a "tail") and is deep enough is a valid path.
        for (size_t i = 0; i < matchNum; ++i) {
            if (!connectedToNext[i] && localMatchPaths[i].depth >= MIN_DEPTH) {
                filteredMatchPaths.push_back(localMatchPaths[i]);
            }
        }

    }
}