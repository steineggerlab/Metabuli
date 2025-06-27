#include "SeqIterator.h"

SeqIterator::~SeqIterator() {
    delete[] mask;
    delete[] mask_int;
}

SeqIterator::SeqIterator(const LocalParameters &par, const GeneticCode &geneticCode) 
    : geneticCode(geneticCode) {

    // Mask for spaced k-mer
    size_t maskLen = 8; // par.spaceMask.length();
    mask = new uint32_t[maskLen+1];
    mask_int = new int[maskLen+1];
    spaceNum = 0;
    spaceNum_int = 0;
    string spaceMask = "11111111"; // par.spaceMask;
    smerLen = par.smerLen;
    smerMask = (1u << (5 * smerLen)) - 1;

    for(size_t i = 0; i < maskLen; i++){
        mask[i] = spaceMask[i] - 48;
        mask_int[i] = spaceMask[i] - 48;
        spaceNum += (mask[i] == 0);
        spaceNum_int += (mask[i] == 0);
    }

    if (par.syncmer) {
        kmerScanner = new SyncmerScanner(smerLen, geneticCode);
    } else {
        kmerScanner = new KmerScanner(geneticCode);
    }

    // powers
    size_t pow = 1;
    size_t numOfAlphabets = 0;
    if (par.reducedAA == 0) {
        numOfAlphabets = 21;
        bitsForCodon = 3;
        bitsFor8Codons = 24;
    } else if (par.reducedAA == 1) {
        numOfAlphabets = 16;
        bitsForCodon = 4;
        bitsFor8Codons = 32;
    }
    for (int i = 0; i < kmerLength; i++) {
        powers[i] = pow;
        pow *= numOfAlphabets;
    }
    dnaMask = (1ULL << bitsFor8Codons) - 1;
}

// It translates a DNA sequence to amino acid sequences with six frames
void SeqIterator::sixFrameTranslation(const char *seq, int len, vector<int> *aaFrames) {
    for (int i = 0; i < 6; i++) { aaFrames[i].clear(); }

    size_t end = len - 1;
    // Translate DNA to AA.
    if (len % 3 == 2) {
        for (int i = 0; i < len - 4; i = i + 3) {
            aaFrames[0].push_back(geneticCode.getAA(atcg[seq[i]],atcg[seq[i + 1]],atcg[seq[i + 2]]));
            aaFrames[1].push_back(
                    geneticCode.getAA(atcg[seq[i + 1]], atcg[seq[i + 2]], atcg[seq[i + 3]]));
            aaFrames[2].push_back(
                    geneticCode.getAA(atcg[seq[i + 2]], atcg[seq[i + 3]], atcg[seq[i + 4]]));
            aaFrames[3].push_back(
                    geneticCode.getAA(iRCT[atcg[seq[end - (i + 0)]]], iRCT[atcg[seq[end - (i + 1)]]],
                                      iRCT[atcg[seq[end - (i + 2)]]]));
            aaFrames[4].push_back(
                    geneticCode.getAA(iRCT[atcg[seq[end - (i + 1)]]], iRCT[atcg[seq[end - (i + 2)]]],
                                      iRCT[atcg[seq[end - (i + 3)]]]));
            aaFrames[5].push_back(
                    geneticCode.getAA(iRCT[atcg[seq[end - (i + 2)]]], iRCT[atcg[seq[end - (i + 3)]]],
                                      iRCT[atcg[seq[end - (i + 4)]]]));
        }
    } else if (len % 3 == 1) {
        for (int i = 0; i < len - 4; i = i + 3) {
            aaFrames[0].push_back(geneticCode.getAA(atcg[seq[i]],atcg[seq[i + 1]],atcg[seq[i + 2]]));
            aaFrames[1].push_back(
                    geneticCode.getAA(atcg[seq[i + 1]], atcg[seq[i + 2]], atcg[seq[i + 3]]));
            aaFrames[2].push_back(
                    geneticCode.getAA(atcg[seq[i + 2]], atcg[seq[i + 3]], atcg[seq[i + 4]]));
            aaFrames[3].push_back(
                    geneticCode.getAA(iRCT[atcg[seq[end - (i + 0) - 3]]], iRCT[atcg[seq[end - (i + 1) - 3]]],
                                      iRCT[atcg[seq[end - (i + 2) - 3]]]));
            aaFrames[4].push_back(
                    geneticCode.getAA(iRCT[atcg[seq[end - (i + 1) - 3]]], iRCT[atcg[seq[end - (i + 2) - 3]]],
                                      iRCT[atcg[seq[end - (i + 3) - 3]]]));
            aaFrames[5].push_back(
                    geneticCode.getAA(iRCT[atcg[seq[end - (i + 2)]]], iRCT[atcg[seq[end - (i + 3)]]],
                                      iRCT[atcg[seq[end - (i + 4)]]]));
        }
    } else {
        for (int i = 0; i < len - 4; i = i + 3) {
            aaFrames[0].push_back(geneticCode.getAA(atcg[seq[i]],atcg[seq[i + 1]],atcg[seq[i + 2]]));
            aaFrames[1].push_back(
                    geneticCode.getAA(atcg[seq[i + 1]], atcg[seq[i + 2]], atcg[seq[i + 3]]));
            aaFrames[2].push_back(
                    geneticCode.getAA(atcg[seq[i + 2]], atcg[seq[i + 3]], atcg[seq[i + 4]]));
            aaFrames[3].push_back(
                    geneticCode.getAA(iRCT[atcg[seq[end - (i + 0) - 3]]], iRCT[atcg[seq[end - (i + 1) - 3]]],
                                      iRCT[atcg[seq[end - (i + 2) - 3]]]));
            aaFrames[4].push_back(
                    geneticCode.getAA(iRCT[atcg[seq[end - (i + 1)]]], iRCT[atcg[seq[end - (i + 2)]]],
                                      iRCT[atcg[seq[end - (i + 3)]]]));
            aaFrames[5].push_back(
                    geneticCode.getAA(iRCT[atcg[seq[end - (i + 2)]]], iRCT[atcg[seq[end - (i + 3)]]],
                                      iRCT[atcg[seq[end - (i + 4)]]]));
        }
    }
}

bool SeqIterator::translateBlock(const char *seq, SequenceBlock block, vector<int> & aaSeq, size_t length) {
    aaSeq.clear();
    if(aaSeq.capacity() < length / 3 + 1) {
        aaSeq.reserve(length / 3 + 1);
    }
    if (block.strand == 1) {
        for (int i = block.start; i + 2 <= block.end; i = i + 3) {
            aaSeq.push_back(geneticCode.getAA(atcg[seq[i]], atcg[seq[i + 1]], atcg[seq[i + 2]]));
        }
    } else {
        for (int i = block.end; i >= (int) block.start + 2; i = i - 3) {
            aaSeq.push_back(geneticCode.getAA(iRCT[atcg[seq[i]]], iRCT[atcg[seq[i - 1]]], iRCT[atcg[seq[i - 2]]]));
        }
    }
    return true;
}

string SeqIterator::reverseComplement(string &read) const {
    int len = read.length();
    string out;
    for (int i = 0; i < len; i++) {
        out.push_back(iRCT[read[i]]);
    }
    reverse(out.begin(), out.end());
    return out;
}

char *SeqIterator::reverseComplement(char *read, size_t length) const {
    char *revCom = (char *) malloc(sizeof(char) * (length + 1));
    for (size_t i = 0; i < length; i++) {
        revCom[length - i - 1] = iRCT[read[i]];
    }
    revCom[length] = '\0';
    return revCom;
}

// It extracts kmers from amino acid sequence with DNA information and fill the kmerBuffer with them.
int SeqIterator::fillBufferWithKmerFromBlock(const char *seq,
                                             Buffer<TargetKmer> &kmerBuffer,
                                             size_t &posToWrite, 
                                             int seqID, 
                                             int taxIdAtRank, 
                                             const vector<int> & aaSeq,
                                             int blockStrand,
                                             int blockStart,
                                             int blockEnd) {
    uint64_t tempKmer = 0;
    int len = (int) aaSeq.size();
    int checkN;
    for (int kmerCnt = 0; kmerCnt < len - kmerLength - spaceNum_int + 1; kmerCnt++) {
        tempKmer = 0;
        checkN = 0;
        for (uint32_t i = 0, j = 0; i < kmerLength + spaceNum; i++, j += mask[i]) {
            if (-1 == aaSeq[kmerCnt + i]) {
                checkN = 1;
                break;
            }
            tempKmer += aaSeq[kmerCnt + i] * powers[j] * mask[i];
        }
        if (checkN == 1) {
            kmerBuffer.buffer[posToWrite] = {UINT64_MAX, -1, 0};
        } else {
            addDNAInfo_TargetKmer(tempKmer, seq, kmerCnt, blockStrand, blockStart, blockEnd);
            kmerBuffer.buffer[posToWrite] = {tempKmer, taxIdAtRank, seqID};
        }
        posToWrite++;
    }
    return 0;
}

int SeqIterator::extractTargetKmers(
    const char *seq,
    Buffer<TargetKmer> &kmerBuffer,
    size_t &posToWrite,
    int seqID,
    int taxIdAtRank,
    SequenceBlock block) 
{
    kmerScanner->initScanner(seq, block.start, block.end);
    bool isForward = block.strand > -1;
    while (true) {
        Kmer syncmer = kmerScanner->next(isForward);
        if (syncmer.value == UINT64_MAX) break; // No more syncmers
        kmerBuffer.buffer[posToWrite++] = {syncmer.value, taxIdAtRank, seqID};
    }
    return 0;
}

// Core syncmer extraction with lambda abstraction
// Simplified emission loop by merging AA and DNA append into one loop
// template<typename AAGetter, typename CodonGetter>
// int SeqIterator::extactTargetSyncmers(const char *seq,
//                                Buffer<TargetKmer> &buffer,
//                                size_t &pos,
//                                int seqID,
//                                int taxId,
//                                SequenceBlock block,
//                                AAGetter getAA,
//                                CodonGetter getCodonBits) {
//     int len = (block.end - block.start + 1) / 3;
//     uint64_t aaPart = 0, dnaPart = 0;
//     deque<Smer> dq;
//     int smerCnt = 0, loaded = 0;
//     uint64_t smer = 0;
//     int prevPos = -kmerLength;
//     int posStart = 0;

//     while (posStart <= len - kmerLength) {
//         bool sawN = false;
//         // Slide the sub-mer window
//         smerCnt -= (smerCnt > 0);
//         while (smerCnt < kmerLength - smerLen + 1) {
//             loaded -= (loaded == smerLen);
//             while (loaded < smerLen) {
//                 int aa = getAA(posStart + smerCnt + loaded);
//                 if (aa < 0) { sawN = true; break; }
//                 smer = (smer << 5) | (uint64_t)aa;
//                 loaded++;
//             }
//             if (sawN) break;
//             smer &= smerMask;
//             while (!dq.empty() && dq.back().value > smer) dq.pop_back();
//             dq.emplace_back(smer, posStart + smerCnt);
//             smerCnt++;
//         }
//         if (sawN) {
//             posStart += smerCnt + loaded + 1;
//             dq.clear(); smerCnt = loaded = 0; smer = 0;
//             continue;
//         }
//         if (!dq.empty() && dq.front().pos < posStart) dq.pop_front();

//         int anchor1 = posStart;
//         int anchor2 = posStart + (kmerLength - smerLen);
//         if (!dq.empty() && (dq.front().pos == anchor1 || dq.front().pos == anchor2)) {
//             int shifts = posStart - prevPos;
//             int aaStart = posStart + kmerLength - shifts;
//             // Single loop: append both AA and DNA bits
//             for (int i = 0; i < shifts; ++i) {
//                 aaPart = (aaPart << 5) | (uint64_t)getAA(aaStart + i); // posStart + kmerLength - shifts + i
//                 dnaPart = (dnaPart << bitsForCodon) | getCodonBits(aaStart + i);
//             }
//             prevPos = posStart;

//             uint64_t key = (aaPart << bitsFor8Codons) | (dnaPart & dnaMask);
//             buffer.buffer[pos++] = { key, taxId, seqID };
//         }
//         ++posStart;
//     }
//     return 0;
// }

// It adds DNA information to kmers referring the original DNA sequence.
// void SeqIterator::addDNAInfo_TargetKmer(uint64_t &kmer, const char *seq, int kmerCnt) {
//     kmer <<= bitsFor8Codons;
//     int start = kmerCnt * 3;
//     for (int i = 0, j = 0; i < kmerLength + spaceNum_int; i ++, j += mask_int[i]) {
//             kmer |= (nuc2num[nuc2int(atcg[seq[start + i*3]])][nuc2int(atcg[seq[start + i*3 + 1]])][
//                     nuc2int(atcg[seq[start + i*3 + 2]])] * mask[i]) << (j * bitsForCodon);
//     }
// }

// // It adds DNA information to kmers referring the original DNA sequence.
// void SeqIterator::addDNAInfo_TargetKmer(uint64_t &kmer, const char *seq, const PredictedBlock &block, int kmerCnt) {
//     kmer <<= bitsFor8Codons;
//     if (block.strand == 1) {
//         int start = block.start + (kmerCnt * 3);
//         for (int i = 0, j = 0; i < kmerLength + spaceNum_int; i ++, j += mask_int[i]) {
//             kmer |= (nuc2num[nuc2int(atcg[seq[start + i*3]])][nuc2int(atcg[seq[start + i*3 + 1]])][
//                     nuc2int(atcg[seq[start + i*3 + 2]])] * mask[i]) << (j * bitsForCodon);
//         }
//     } else {
//         int start = block.end - (kmerCnt * 3);
//         for (int i = 0, j = 0; i < kmerLength + spaceNum_int; i++, j += mask_int[i]) {
//             kmer |= (nuc2num[nuc2int(iRCT[atcg[seq[start - i*3]]])][nuc2int(iRCT[atcg[seq[start - i*3 - 1]]])][
//                     nuc2int(iRCT[atcg[seq[start - i*3 - 2]]])] * mask[i]) << (j * bitsForCodon);
//         }
//     }
// }

void SeqIterator::addDNAInfo_TargetKmer(uint64_t & kmer, const char * seq, int kmerCnt, int strand = 0, int start = 0, int end = 0) {
    kmer <<= bitsFor8Codons;
    if (strand > -1) {
        int startPos = start + (kmerCnt * 3);
        for (int i = 0, j = 0; i < kmerLength + spaceNum_int; i ++, j += mask_int[i]) {
            kmer |= (geneticCode.getCodon(atcg[seq[startPos + i*3]], atcg[seq[startPos + i*3 + 1]], atcg[seq[startPos + i*3 + 2]])
                    * mask[i]) << (j * bitsForCodon);
        }
    } else {
        int startPos = end - (kmerCnt * 3);
        for (int i = 0, j = 0; i < kmerLength + spaceNum_int; i++, j += mask_int[i]) {
            kmer |= (geneticCode.getCodon(iRCT[atcg[seq[startPos - i*3]]],
                    iRCT[atcg[seq[startPos - i*3 - 1]]],
                    iRCT[atcg[seq[startPos - i*3 - 2]]]) * mask[i]) << (j * bitsForCodon);
        }
    }

}

size_t SeqIterator::kmerNumOfSixFrameTranslation(const char *seq) {
    size_t len = strlen(seq);
    return (len / 3 - kmerLength) * 6;
}

size_t SeqIterator::getNumOfKmerForBlock(const SequenceBlock &block) {
    size_t len = block.end - block.start + 1;
    return len / 3 - 7;
}


bool SeqIterator::compareMinHashList(priority_queue <uint64_t> list1, priority_queue <uint64_t> &list2, size_t length1,
                                     size_t length2) {
    float lengthRatio = float(length2) / float(length1);
    float identicalCount = 0;
    float list1Size = list1.size();
    while (!list1.empty() && !list2.empty()) {
        if (list1.top() == list2.top()) {
            identicalCount++;
            list1.pop();
            list2.pop();
        } else if (list1.top() > list2.top()) {
            list1.pop();
        } else if (list1.top() < list2.top()) {
            list2.pop();
        }
    }
    if (identicalCount > list1Size * lengthRatio * 0.5) {
        return true;
    } else {
        return false;
    }
}

void SeqIterator::getMinHashList(priority_queue <uint64_t> &sortedHashQue, const char *seq) {
    size_t seqLength = strlen(seq);
    size_t kmerLegnth = 24;
    char *kmer = (char *) malloc(sizeof(char) * (kmerLegnth + 1));
    kmer[kmerLegnth] = '\0';
    size_t queLength = 0;
    size_t maxLength = 3000;
    size_t currHash;
    sortedHashQue.push(UINT64_MAX);

    for (size_t i = 0; i + kmerLegnth - 1 < seqLength; i++) {
        strncpy(kmer, seq + i, kmerLegnth);
        currHash = XXH64(kmer, kmerLegnth, 0);
        if (currHash < sortedHashQue.top()) {
            if (queLength < maxLength) {
                sortedHashQue.push(currHash);
                queLength++;
            } else {
                sortedHashQue.pop();
                sortedHashQue.push(currHash);
            }
        }
    }
    free(kmer);
}

void SeqIterator::generateIntergenicKmerList(_gene *genes, _node *nodes, int numberOfGenes,
                                             vector <uint64_t> &intergenicKmerList,
                                             const char *seq) {
    if (numberOfGenes == 0) return;

    int k = 23;
    char *kmer = (char *) malloc(sizeof(char) * (k + 1));
    char *reverseKmer = (char *) malloc(sizeof(char) * (k + 1));

    // Use the frame of the first gene for the first intergenic region
    int beginOfFisrtGene = genes[0].begin - 1;
    if (beginOfFisrtGene > k - 1) {
        strncpy(kmer, seq + beginOfFisrtGene - k, k);
        if (nodes[genes[0].start_ndx].strand == 1) {
            intergenicKmerList.push_back(XXH64(kmer, k, 0));
        } else {
            for (int j = k - 1; j >= 0; j--) {
                reverseKmer[k - j - 1] = iRCT[kmer[j]];
            }
            intergenicKmerList.push_back(XXH64(reverseKmer, k, 0));
        }
    }

    //
    for (int i = 0; i < numberOfGenes; i++) {
        strncpy(kmer, seq + genes[i].end, k);
        if (nodes[genes[i].start_ndx].strand == 1) {
            intergenicKmerList.push_back(XXH64(kmer, k, 0));
        } else {
            for (int j = k - 1; j >= 0; j--) {
                reverseKmer[k - j - 1] = iRCT[kmer[j]];
            }
            intergenicKmerList.push_back(XXH64(reverseKmer, k, 0));
        }
    }

    free(reverseKmer);
    free(kmer);
}

void SeqIterator::maskLowComplexityRegions(const char *seq, char *maskedSeq, ProbabilityMatrix & probMat,
                                           float maskProb, const BaseMatrix * subMat) {
    unsigned int seqLen = 0;
    while (seq[seqLen] != '\0') {
        maskedSeq[seqLen] = (char) subMat->aa2num[static_cast<int>(seq[seqLen])];
        seqLen++;
    }
    tantan::maskSequences(maskedSeq,
                          maskedSeq + seqLen,
                          50 /*options.maxCycleLength*/,
                          probMat.probMatrixPointers,
                          0.005 /*options.repeatProb*/,
                          0.05 /*options.repeatEndProb*/,
                          0.9 /*options.repeatOffsetProbDecay*/,
                          0, 0,
                          maskProb /*options.minMaskProb*/,
                          probMat.hardMaskTable);
    for (unsigned int pos = 0; pos < seqLen; pos++) {
        char nt = seq[pos];
        maskedSeq[pos] = (maskedSeq[pos] == probMat.hardMaskTable[0]) ? 'N' : nt;
    }
}




void SeqIterator::devideToCdsAndNonCds(const char *maskedSeq,
                                       size_t seqLen,
                                       const vector<CDSinfo> &cdsInfo, 
                                       vector<string> &cds,
                                       vector<string> &nonCds) {
    string tmpMasked;
    for (size_t i = 0; i < cdsInfo.size(); i++) {
        size_t locNum = cdsInfo[i].loc.size();
        int currStartCodonPos = 0;
        int extend = 0;
        for (size_t j = 0; j < locNum; j++) {
            // Extend 21 bases to both sides for k-mer from CDS boudaries
            size_t begin = cdsInfo[i].loc[j].first - 1;
            size_t end = cdsInfo[i].loc[j].second - 1;
            if (j == 0) {
                int k = 0;
                while (k < 7 && begin >= 3) {
                    begin -= 3;
                    currStartCodonPos += 3;
                    extend += 3;
                    k++;
                }
            }
            if (j == locNum - 1) {
                int k = 0;
                while (k < 7 && end + 3 < seqLen) {
                    end += 3;
                    extend += 3;
                    k++;
                }
            }    
            // cout << seqLen << " begin: " << begin << " end: " << end << endl;
            tmpMasked += string(maskedSeq + begin, end - begin + 1);    
        }
        // Reverse complement if needed
        if (cdsInfo[i].isComplement) {
            cds.emplace_back(reverseComplement(tmpMasked));
        } else {
            cds.emplace_back(tmpMasked);
        }
        tmpMasked.clear();
    }

    // Get non-CDS that are not in the CDS list
    bool * checker = new bool[seqLen];
    memset(checker, true, seqLen * sizeof(bool));
    for (size_t i = 0; i < cdsInfo.size(); i++) {
        for (size_t j = 0; j < cdsInfo[i].loc.size(); j++) {
            for (size_t k = cdsInfo[i].loc[j].first - 1; k < cdsInfo[i].loc[j].second; k++) {
                checker[k] = false;
            }
        }
    }

    size_t len;
    size_t i = 0;
    while (i < seqLen) {
        len = 0;
        while (i < seqLen && checker[i]) {
            i++;
            len++;
        }
        if (len > 32) {
            nonCds.emplace_back(string(maskedSeq + i - len, len));
        }
        i ++;
    }
}

void SeqIterator::printKmerInDNAsequence(uint64_t kmer) {
    if (bitsForCodon == 4) {
        uint64_t copy = kmer;
        kmer >>= bitsFor8Codons;
        int quotient;
        int dnaInfo;
        vector<int> aa8mer(8);
        vector<string> dna24mer(8);
        for (int i = 0; i < 8; i++) {
            quotient = kmer / powers[7 - i];
            kmer = kmer - (quotient * powers[7 - i]);
            aa8mer[7 - i] = quotient;
        }

        string aminoacid = "ARNDCQGHILKFPSTXVWYME";
        // Print Amino Acid 8 mer
        for (int i = 0; i < 8; i++) {
            cout << aminoacid[aa8mer[i]];
        }
        cout << " ";

        for (int i = 0; i < 8; i++) {
            if (bitsForCodon == 4) {
                dnaInfo = copy & 0Xfu;
            } else {
                dnaInfo = copy & 7u;
            }

            copy >>= bitsForCodon;
            switch (aa8mer[i]) {
                case 0: //A
//                cout << "A";
                    if (dnaInfo == 0) {
                        dna24mer[7 - i] = "GCA";
                    } else if (dnaInfo == 1) {
                        dna24mer[7 - i] = "GCC";
                    } else if (dnaInfo == 2) {
                        dna24mer[7 - i] = "GCT";
                    } else if (dnaInfo == 3) {
                        dna24mer[7 - i] = "GCG";
                    } else {
                        cout << "Error in A" << endl;
                    }
                    break;
                case 1: //R
//                cout << "R";
                    if (dnaInfo == 0) {
                        dna24mer[7 - i] = "CGA";
                    } else if (dnaInfo == 1) {
                        dna24mer[7 - i] = "CGC";
                    } else if (dnaInfo == 2) {
                        dna24mer[7 - i] = "CGT";
                    } else if (dnaInfo == 3) {
                        dna24mer[7 - i] = "CGG";
                    } else if (dnaInfo == 4) {
                        dna24mer[7 - i] = "AGA";
                    } else if (dnaInfo == 7) {
                        dna24mer[7 - i] = "AGG";
                    } else {
                        cout << "Error in R" << endl;
                    }
                    break;
                case 2: //N
//                cout << "N";
                    if (dnaInfo == 1) {
                        dna24mer[7 - i] = "AAC";
                    } else if (dnaInfo == 2) {
                        dna24mer[7 - i] = "AAT";
                    } else {
                        cout << "Error in " << aminoacid[aa8mer[i]] << endl;
                    }
                    break;
                case 3: //D
//                cout << "D";
                    if (dnaInfo == 1) {
                        dna24mer[7 - i] = "GAC";
                    } else if (dnaInfo == 2) {
                        dna24mer[7 - i] = "GAT";
                    } else {
                        cout << "Error in " << aminoacid[aa8mer[i]] << endl;
                    }
                    break;
                case 4: //C
//                cout << "C";
                    if (dnaInfo == 1) {
                        dna24mer[7 - i] = "TGC";
                    } else if (dnaInfo == 2) {
                        dna24mer[7 - i] = "TGT";
                    } else {
                        cout << "Error in " << aminoacid[aa8mer[i]] << endl;
                    }
                    break;
                case 5: // Q
//                cout << "Q";
                    if (dnaInfo == 0) {
                        dna24mer[7 - i] = "CAA";
                    } else if (dnaInfo == 3) {
                        dna24mer[7 - i] = "CAG";
                    } else if (dnaInfo == 4) {
                        dna24mer[7 - i] = "GAA";
                    } else if (dnaInfo == 7) {
                        dna24mer[7 - i] = "GAG";
                    } else {
                        cout << "Error in " << aminoacid[aa8mer[i]] << endl;
                    }
                    break;
                case 6: //G
                    if (dnaInfo == 0) {
                        dna24mer[7 - i] = "GGA";
                    } else if (dnaInfo == 1) {
                        dna24mer[7 - i] = "GGC";
                    } else if (dnaInfo == 2) {
                        dna24mer[7 - i] = "GGT";
                    } else if (dnaInfo == 3) {
                        dna24mer[7 - i] = "GGG";
                    } else {
                        cout << "Error in " << aminoacid[aa8mer[i]] << endl;
                    }
                    break;
                case 7: //H
                    if (dnaInfo == 0) {
                        cout << "Error in " << aminoacid[aa8mer[i]] << endl;
                    } else if (dnaInfo == 1) {
                        dna24mer[7 - i] = "CAC";
                    } else if (dnaInfo == 2) {
                        dna24mer[7 - i] = "CAT";
                    } else {
                        cout << "Error in " << aminoacid[aa8mer[i]] << endl;
                    }
                    break;
                case 8: //I
                    // IV

                    if (dnaInfo == 0) {
                        dna24mer[7 - i] = "GTA";
                    } else if (dnaInfo == 1) {
                        dna24mer[7 - i] = "GTC";
                    } else if (dnaInfo == 2) {
                        dna24mer[7 - i] = "GTT";
                    } else if (dnaInfo == 3) {
                        dna24mer[7 - i] = "GTG";
                    } else if (dnaInfo == 4) {
                        dna24mer[7 - i] = "ATA";
                    } else if (dnaInfo == 5) {
                        dna24mer[7 - i] = "ATC";
                    } else if (dnaInfo == 6) {
                        dna24mer[7 - i] = "ATT";
                    } else {
                        cout << "Error in " << aminoacid[aa8mer[i]] << endl;
                    }
                    break;
                case 9: //L
//                cout << "I";
                    if (dnaInfo == 0) {
                        dna24mer[7 - i] = "CTA";
                    } else if (dnaInfo == 1) {
                        dna24mer[7 - i] = "CTC";
                    } else if (dnaInfo == 2) {
                        dna24mer[7 - i] = "CTT";
                    } else if (dnaInfo == 3) {
                        dna24mer[7 - i] = "CTG";
                    } else if (dnaInfo == 4) {
                        dna24mer[7 - i] = "TTA";
                    } else if (dnaInfo == 7) {
                        dna24mer[7 - i] = "TTG";
                    } else if (dnaInfo == 8) {
                        dna24mer[7 - i] = "ATG";
                    } else {
                        cout << "Error in " << aminoacid[aa8mer[i]] << endl;
                    }
                    break;
                case 10: //K
//                cout << "L";
                    if (dnaInfo == 0) {
                        dna24mer[7 - i] = "AAA";
                    } else if (dnaInfo == 3) {
                        dna24mer[7 - i] = "AAG";
                    } else {

                    }
                    break;
                case 11: //FYW
                    if (dnaInfo == 1) {
                        dna24mer[7 - i] = "TTC";
                    } else if (dnaInfo == 2) {
                        dna24mer[7 - i] = "TTT";
                    } else if (dnaInfo == 5) {
                        dna24mer[7 - i] = "TAC";
                    } else if (dnaInfo == 6) {
                        dna24mer[7 - i] = "TAT";
                    } else if (dnaInfo == 7) {
                        dna24mer[7 - i] = "TGG";
                    } else {
                        cout << "Error in " << aminoacid[aa8mer[i]] << endl;
                    }
                    break;
                case 12: //P
                    if (dnaInfo == 0) {
                        dna24mer[7 - i] = "CCA";
                    } else if (dnaInfo == 1) {
                        dna24mer[7 - i] = "CCC";
                    } else if (dnaInfo == 2) {
                        dna24mer[7 - i] = "CCT";
                    } else if (dnaInfo == 3) {
                        dna24mer[7 - i] = "CCG";
                    } else {
                        cout << "Error in " << aminoacid[aa8mer[i]] << endl;
                    }
                    break;
                case 13://S
                    if (dnaInfo == 0) {
                        dna24mer[7 - i] = "TCA";
                    } else if (dnaInfo == 1) {
                        dna24mer[7 - i] = "TCC";
                    } else if (dnaInfo == 2) {
                        dna24mer[7 - i] = "TCT";
                    } else if (dnaInfo == 3) {
                        dna24mer[7 - i] = "TCG";
                    } else if (dnaInfo == 9) {
                        dna24mer[7 - i] = "AGC";
                    } else if (dnaInfo == 10) {
                        dna24mer[7 - i] = "AGT";
                    } else {
                        cout << "Error in " << aminoacid[aa8mer[i]] << endl;
                    }
                    break;
                case 14: //T
//                cout << "P";
                    if (dnaInfo == 0) {
                        dna24mer[7 - i] = "ACA";
                    } else if (dnaInfo == 1) {
                        dna24mer[7 - i] = "ACC";
                    } else if (dnaInfo == 2) {
                        dna24mer[7 - i] = "ACT";
                    } else if (dnaInfo == 3) {
                        dna24mer[7 - i] = "ACG";
                    } else {
                        cout << "Error in " << aminoacid[aa8mer[i]] << endl;
                    }
                    break;
                case 15: //STOP
                    if (dnaInfo == 0) {
                        dna24mer[7 - i] = "TAA";
                    } else if (dnaInfo == 4) {
                        dna24mer[7 - i] = "TGA";
                    } else if (dnaInfo == 3) {
                        dna24mer[7 - i] = "TAG";
                    } else {
                        cout << "Error in " << aminoacid[aa8mer[i]] << endl;
                    }
                    break;
                case 16: //V
                    if (dnaInfo == 0) {
                        dna24mer[7 - i] = "GTA";
                    } else if (dnaInfo == 1) {
                        dna24mer[7 - i] = "GTC";
                    } else if (dnaInfo == 2) {
                        dna24mer[7 - i] = "GTT";
                    } else if (dnaInfo == 3) {
                        dna24mer[7 - i] = "GTG";
                    } else {
                        cout << "Error in " << aminoacid[aa8mer[i]] << endl;
                    }
                    break;
                case 17: //W
                    if (dnaInfo == 3) {
                        dna24mer[7 - i] = "TGG";
                    } else {
                        cout << "Error in " << aminoacid[aa8mer[i]] << endl;
                    }
                    break;
                case 18: //Y
                    if (dnaInfo == 1) {
                        dna24mer[7 - i] = "TAC";
                    } else if (dnaInfo == 2) {
                        dna24mer[7 - i] = "TAT";
                    } else {
                        cout << "Error in " << aminoacid[aa8mer[i]] << endl;
                    }
                    break;
                case 19: //M
                    if (dnaInfo == 3) {
                        dna24mer[7 - i] = "ATG";
                    } else {
                        cout << "Error in " << aminoacid[aa8mer[i]] << endl;
                    }
                    break;
                case 20: //E
                    if (dnaInfo == 0) {
                        dna24mer[7 - i] = "GAA";
                    } else if (dnaInfo == 3) {
                        dna24mer[7 - i] = "GAG";
                    } else {
                        cout << "Error in " << aminoacid[aa8mer[i]] << endl;
                    }
                    break;
            }
        }

        for (int i = 0; i < 8; i++) {
            cout << dna24mer[7 - i];
        }
    }
    else {
        uint64_t copy = kmer;
        kmer >>= 24;
        int quotient;
        int dnaInfo;
        vector<int> aa8mer(8);
        vector<string> dna24mer(8);
        for (int i = 0; i < 8; i++) {
            quotient = kmer / powers[7 - i];
            kmer = kmer - (quotient * powers[7 - i]);
            aa8mer[7 - i] = quotient;
        }

//        ///Print Amino Acid 8 mer
//        for(int i  = 0; i<8; i++){
//            cout<<aa8mer[i]<<" ";
//        }
//        cout<<endl;

        string aminoacid = "ARNDCQEGHILKMFPSTWYVX";
        for (int i = 0; i < 8; i++) {
            cout << aminoacid[aa8mer[i]];
        }
        cout<<" ";
        for (int i = 0; i < 8; i++) {
            dnaInfo = copy & 7u;
            copy >>= 3;
            switch (aa8mer[i]) {
                case 0: //A
//                cout << "A";
                    if (dnaInfo == 0) {
                        dna24mer[7 - i] = "GCA";
                    } else if (dnaInfo == 1) {
                        dna24mer[7 - i] = "GCC";
                    } else if (dnaInfo == 2) {
                        dna24mer[7 - i] = "GCT";
                    } else if (dnaInfo == 3){
                        dna24mer[7 - i] = "GCG";
                    } else {
                        cout << "Error in " << aminoacid[aa8mer[i]] << endl;
                    }
                    break;
                case 1: //R
//                cout << "R";
                    if (dnaInfo == 0) {
                        dna24mer[7 - i] = "CGA";
                    } else if (dnaInfo == 1) {
                        dna24mer[7 - i] = "CGC";
                    } else if (dnaInfo == 2) {
                        dna24mer[7 - i] = "CGT";
                    } else if (dnaInfo == 3) {
                        dna24mer[7 - i] = "CGG";
                    } else if (dnaInfo == 4) {
                        dna24mer[7 - i] = "AGG";
                    } else if (dnaInfo == 5) {
                        dna24mer[7 - i] = "AGA";
                    } else{
                        cout << "Error in " << aminoacid[aa8mer[i]] << endl;
                    }
                    break;
                case 2: //N
//                cout << "N";
                    if (dnaInfo == 0) {
                        cout << "dddddd";
                    } else if (dnaInfo == 1) {
                        dna24mer[7 - i] = "AAC";
                    } else if (dnaInfo == 2) {
                        dna24mer[7 - i] = "AAT";
                    } else {
                        cout << "dddddd";
                    }
                    break;
                case 3: //D
//                cout << "D";
                    if (dnaInfo == 0) {
                        cout << "dddddd";
                    } else if (dnaInfo == 1) {
                        dna24mer[7 - i] = "GAC";
                    } else if (dnaInfo == 2) {
                        dna24mer[7 - i] = "GAT";
                    } else {
                        cout << "dddddd";
                    }
                    break;
                case 4: //C
//                cout << "C";
                    if (dnaInfo == 0) {
                        cout << "dddddd";
                    } else if (dnaInfo == 1) {
                        dna24mer[7 - i] = "TGC";
                    } else if (dnaInfo == 2) {
                        dna24mer[7 - i] = "TGT";
                    } else {
                        cout << "dddddd";
                    }
                    break;
                case 5: // Q
//                cout << "Q";
                    if (dnaInfo == 0) {
                        dna24mer[7 - i] = "CAA";
                    } else if (dnaInfo == 1) {
                        cout << "dddddd";
                    } else if (dnaInfo == 2) {
                        cout << "dddddd";
                    } else {
                        dna24mer[7 - i] = "CAG";
                    }
                    break;
                case 6: //E
//                cout << "E";
                    if (dnaInfo == 0) {
                        dna24mer[7 - i] = "GAA";
                    } else if (dnaInfo == 1) {
                        cout << "dddddd";
                    } else if (dnaInfo == 2) {
                        cout << "dddddd" ;
                    } else {
                        dna24mer[7 - i] = "GAG";
                    }
                    break;
                case 7: //G
//                cout << "G";
                    if (dnaInfo == 0) {
                        dna24mer[7 - i] = "GGA";
                    } else if (dnaInfo == 1) {
                        dna24mer[7 - i] = "GGC";
                    } else if (dnaInfo == 2) {
                        dna24mer[7 - i] = "GGT";
                    } else {
                        dna24mer[7 - i] = "GGG";
                    }
                    break;
                case 8: //H
//                cout << "H";
                    if (dnaInfo == 0) {
                        cout << "dddddd";
                    } else if (dnaInfo == 1) {
                        dna24mer[7 - i] = "CAC";
                    } else if (dnaInfo == 2) {
                        dna24mer[7 - i] = "CAT";
                    } else {
                        cout << "dddddd";
                    }
                    break;
                case 9: //I
//                cout << "I";
                    if (dnaInfo == 0) {
                        dna24mer[7 - i] = "ATA";
                    } else if (dnaInfo == 1) {
                        dna24mer[7 - i] = "ATC";
                    } else if (dnaInfo == 2) {
                        dna24mer[7 - i] = "ATT";
                    } else {
                        cout << "dddddd";
                    }
                    break;
                case 10: //L
//                cout << "L";
                    if (dnaInfo == 0) {
                        dna24mer[7 - i] = "CTA";
                    } else if (dnaInfo == 1) {
                        dna24mer[7 - i] = "CTC";
                    } else if (dnaInfo == 2) {
                        dna24mer[7 - i] = "CTT";
                    } else if (dnaInfo == 3) {
                        dna24mer[7 - i] = "CTG";
                    } else if (dnaInfo == 4) {
                        dna24mer[7 - i] = "TTG";
                    } else {
                        dna24mer[7 - i] = "TTA";
                    }
                    break;
                case 11: //K
//                cout << "K";
                    if (dnaInfo == 0) {
                        dna24mer[7 - i] = "AAA";
                    } else if (dnaInfo == 1) {
                        cout << "dddddd";
                    } else if (dnaInfo == 2) {
                        cout << "dddddd";
                    } else {
                        dna24mer[7 - i] = "AAG";
                    }
                    break;
                case 12: // M
//                cout << "M";
                    if (dnaInfo == 0) {
                        cout << "dddddd";
                    } else if (dnaInfo == 1) {
                        cout << "dddddd";
                    } else if (dnaInfo == 2) {
                        cout << "dddddd";
                    } else {
                        dna24mer[7 - i] = "ATG";
                    }
                    break;
                case 13://F
//                cout << "F";
                    if (dnaInfo == 0) {
                        cout << "dddddd" ;
                    } else if (dnaInfo == 1) {
                        dna24mer[7 - i] = "TTC";
                    } else if (dnaInfo == 2) {
                        dna24mer[7 - i] = "TTT";
                    } else {
                        cout << "dddddd";
                    }
                    break;
                case 14: //P
//                cout << "P";
                    if (dnaInfo == 0) {
                        dna24mer[7 - i] = "CCA";
                    } else if (dnaInfo == 1) {
                        dna24mer[7 - i] = "CCC";
                    } else if (dnaInfo == 2) {
                        dna24mer[7 - i] = "CCT";
                    } else {
                        dna24mer[7 - i] = "CCG";
                    }
                    break;
                case 15: //S
//                cout << "S";
                    if (dnaInfo == 0) {
                        dna24mer[7 - i] = "TCA";
                    } else if (dnaInfo == 1) {
                        dna24mer[7 - i] = "TCC";
                    } else if (dnaInfo == 2) {
                        dna24mer[7 - i] = "TCT";
                    } else if (dnaInfo == 3) {
                        dna24mer[7 - i] = "TCG";
                    } else if (dnaInfo == 4) {
                        cout << "dddddd";
                    } else if (dnaInfo == 5) {
                        cout << "dddddd";
                    } else if (dnaInfo == 6) {
                        dna24mer[7 - i] = "AGT";
                    } else {
                        dna24mer[7 - i] = "AGC";
                    }
                    break;
                case 16: //T
//                cout << "T";
                    if (dnaInfo == 0) {
                        dna24mer[7 - i] = "ACA";
                    } else if (dnaInfo == 1) {
                        dna24mer[7 - i] = "ACC";
                    } else if (dnaInfo == 2) {
                        dna24mer[7 - i] = "ACT";
                    } else {
                        dna24mer[7 - i] = "ACG";
                    }
                    break;
                case 17: //W
//                cout << "W";
                    if (dnaInfo == 0) {
                        cout << "dddddd";
                    } else if (dnaInfo == 1) {
                        cout << "dddddd";
                    } else if (dnaInfo == 2) {
                        cout << "dddddd";
                    } else {
                        dna24mer[7 - i] = "TGG";
                    }
                    break;
                case 18: //Y
//                cout << "Y";
                    if (dnaInfo == 0) {
                        cout << "dddddd";
                    } else if (dnaInfo == 1) {
                        dna24mer[7 - i] = "TAC";
                    } else if (dnaInfo == 2) {
                        dna24mer[7 - i] = "TAT";
                    } else {
                        cout << "dddddd";
                    }
                    break;
                case 19: //V
//                cout << "V";
                    if (dnaInfo == 0) {
                        dna24mer[7 - i] = "GTA";
                    } else if (dnaInfo == 1) {
                        dna24mer[7 - i] = "GTC";
                    } else if (dnaInfo == 2) {
                        dna24mer[7 - i] = "GTT";
                    } else {
                        dna24mer[7 - i] = "GTG";
                    }
                    break;
                case 20: //stop
//                cout << "Z";
                    if (dnaInfo == 0) {
                        dna24mer[7 - i] = "TAA";
                    } else if (dnaInfo == 1) {
                        cout << "dddddd";
                    } else if (dnaInfo == 2) {
                        cout << "dddddd";
                    } else if (dnaInfo == 3) {
                        dna24mer[7 - i] = "TAG";
                    } else if (dnaInfo == 4) {
                        cout << "dddddd";
                    } else {
                        dna24mer[7 - i] = "TGA";
                    }
                    break;
            }
            //cout<<dnaInfo<<" ";
        }

//    cout<<endl;
        for (int i = 0; i < 8; i++) {
            cout << dna24mer[7-i];
        }
    }
}


void SeqIterator::printAAKmer(uint64_t kmer, int shifts) {
    kmer >>= shifts;
    vector<int> aa8mer(8);
    for (int i = 0; i < 8; i++) {
        int quotient = kmer / powers[7 - i];
        kmer = kmer - (quotient * powers[7 - i]);
        aa8mer[7 - i] = quotient;
    }
    string aminoacid = "ARNDCQEGHILKMFPSTWYVX";
    for (int i = 0; i < 8; i++) {
        cout << aminoacid[aa8mer[i]];
    }
} 