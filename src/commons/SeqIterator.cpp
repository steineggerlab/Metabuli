//
// Created by KJB on 01/09/2020.
//

#include "SeqIterator.h"

const string SeqIterator::atcg = "................................................................"
                                 ".AGCG..GT..G.CN...ACTG.A.T......agcg..gt..g.cn...actg.a.t......."
                                 "................................................................"
                                 "................................................................";

const string SeqIterator::iRCT = "................................................................"
                                 ".TVGH..CD..M.KN...YSAABW.R.......tvgh..cd..m.kn...ysaabw.r......"
                                 "................................................................"
                                 "................................................................";



SeqIterator::~SeqIterator() {
    delete[] mask;
    delete[] mask_int;
}

SeqIterator::SeqIterator(const LocalParameters &par) {

    // Mask for spaced k-mer
    size_t maskLen = par.spaceMask.length();
    mask = new uint32_t[maskLen+1];
    mask_int = new int[maskLen+1];
    spaceNum = 0;
    spaceNum_int = 0;
    for(size_t i = 0; i < maskLen; i++){
        mask[i] = par.spaceMask[i] - 48;
        mask_int[i] = par.spaceMask[i] - 48;
        spaceNum += (mask[i] == 0);
        spaceNum_int += (mask[i] == 0);
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

    // Codon table
    if (par.reducedAA == 0) {
        // A
        nuc2aa[3][1][0] = 0; nuc2aa[3][1][1] = 0; nuc2aa[3][1][2] = 0 ; nuc2aa[3][1][3] = 0;
        // R
        nuc2aa[1][3][0] = 1; nuc2aa[1][3][1] = 1; nuc2aa[1][3][2] = 1 ; nuc2aa[1][3][3] = 1; nuc2aa[0][3][0] = 1; nuc2aa[0][3][3] = 1;
        // N
        nuc2aa[0][0][2] = 2; nuc2aa[0][0][1] = 2;
        // D
        nuc2aa[3][0][2] = 3; nuc2aa[3][0][1] = 3;
        // C
        nuc2aa[2][3][2] = 4; nuc2aa[2][3][1] = 4;
        // Q
        nuc2aa[1][0][0] = 5; nuc2aa[1][0][3] = 5;
        // E
        nuc2aa[3][0][0] = 6; nuc2aa[3][0][3] = 6;
        // G
        nuc2aa[3][3][0] = 7; nuc2aa[3][3][1] = 7; nuc2aa[3][3][2] = 7; nuc2aa[3][3][3] =7;
        // H
        nuc2aa[1][0][2] = 8; nuc2aa[1][0][1] = 8;
        // I
        nuc2aa[0][2][2] = 9; nuc2aa[0][2][1] = 9; nuc2aa[0][2][0] = 9;
        // L
        nuc2aa[2][2][0] = 10; nuc2aa[2][2][3] = 10; nuc2aa[1][2][0] = 10; nuc2aa[1][2][1] = 10;nuc2aa[1][2][2] = 10; nuc2aa[1][2][3] =10 ;
        // K
        nuc2aa[0][0][0] = 11; nuc2aa[0][0][3] = 11;
        // M
        nuc2aa[0][2][3] = 12;
        // F
        nuc2aa[2][2][2] = 13; nuc2aa[2][2][1] = 13;
        // P
        nuc2aa[1][1][0] = 14; nuc2aa[1][1][1] = 14; nuc2aa[1][1][2] = 14; nuc2aa[1][1][3] = 14;
        // S
        nuc2aa[2][1][0] = 15; nuc2aa[2][1][1] = 15; nuc2aa[2][1][2] = 15; nuc2aa[2][1][3] = 15; nuc2aa[0][3][2] = 15; nuc2aa[0][3][1] = 15;
        // T
        nuc2aa[0][1][0] = 16; nuc2aa[0][1][1] = 16; nuc2aa[0][1][2] = 16; nuc2aa[0][1][3] = 16;
        // W
        nuc2aa[2][3][3] = 17;
        // Y
        nuc2aa[2][0][2] = 18; nuc2aa[2][0][1] = 18;
        // V
        nuc2aa[3][2][0] = 19; nuc2aa[3][2][1] = 19; nuc2aa[3][2][2] = 19; nuc2aa[3][2][3] = 19;
        // Stop
        nuc2aa[2][0][0] = 20; nuc2aa[2][3][0] = 20; nuc2aa[2][0][3] = 20;

        // triplet code with N's
        for(int i = 0; i < 8; i++){
            for(int j = 0; j < 8; j++){
                nuc2aa[7][i][j] = -1;
                nuc2aa[i][7][j] = -1;
                nuc2aa[i][j][7] = -1;
            }
        }

        // For encoding DNA information in k-mer
        for(int i =0; i < 4 ; i++){
            for(int i2 = 0; i2 < 4 ; i2++){
                nuc2num[i][i2][0] = 0;
                nuc2num[i][i2][1] = 1;
                nuc2num[i][i2][2] = 2;
                nuc2num[i][i2][3] = 3;
            }
        }
        // for Arg
        nuc2num[0][3][3] = 4;
        nuc2num[0][3][0] = 5;
        // for Leu
        nuc2num[2][2][3] = 4;
        nuc2num[2][2][0] = 5;
        // for Ser
        nuc2num[0][3][2] = 6;
        nuc2num[0][3][1] = 7;
        // for stop codon
        nuc2num[2][3][0] = 5;
    }
        // Using 15 alphabets for amino acid, one for stop codon.
        // FYW, ML, IV, QE
    else if (par.reducedAA == 1) {
        // Codon table
        // A
        nuc2aa[3][1][0] = 0;
        nuc2aa[3][1][1] = 0;
        nuc2aa[3][1][2] = 0;
        nuc2aa[3][1][3] = 0;
        // R
        nuc2aa[1][3][0] = 1;
        nuc2aa[1][3][1] = 1;
        nuc2aa[1][3][2] = 1;
        nuc2aa[1][3][3] = 1;
        nuc2aa[0][3][0] = 1;
        nuc2aa[0][3][3] = 1;
        // N
        nuc2aa[0][0][2] = 2;
        nuc2aa[0][0][1] = 2;
        // D
        nuc2aa[3][0][2] = 3;
        nuc2aa[3][0][1] = 3;
        // C
        nuc2aa[2][3][2] = 4;
        nuc2aa[2][3][1] = 4;
        // QE
        nuc2aa[1][0][0] = 5;
        nuc2aa[1][0][3] = 5;
        nuc2aa[3][0][0] = 5;
        nuc2aa[3][0][3] = 5;
        // G
        nuc2aa[3][3][0] = 6;
        nuc2aa[3][3][1] = 6;
        nuc2aa[3][3][2] = 6;
        nuc2aa[3][3][3] = 6;
        // H
        nuc2aa[1][0][2] = 7;
        nuc2aa[1][0][1] = 7;
        // IV
        nuc2aa[0][2][2] = 8;
        nuc2aa[0][2][1] = 8;
        nuc2aa[0][2][0] = 8;
        nuc2aa[3][2][0] = 8;
        nuc2aa[3][2][1] = 8;
        nuc2aa[3][2][2] = 8;
        nuc2aa[3][2][3] = 8;
        // ML
        nuc2aa[2][2][0] = 9;
        nuc2aa[2][2][3] = 9;
        nuc2aa[1][2][0] = 9;
        nuc2aa[1][2][1] = 9;
        nuc2aa[1][2][2] = 9;
        nuc2aa[1][2][3] = 9;
        nuc2aa[0][2][3] = 9;
        // K
        nuc2aa[0][0][0] = 10;
        nuc2aa[0][0][3] = 10;
        // FYW
        nuc2aa[2][2][2] = 11;
        nuc2aa[2][2][1] = 11;
        nuc2aa[2][0][2] = 11;
        nuc2aa[2][0][1] = 11;
        nuc2aa[2][3][3] = 11;
        // P
        nuc2aa[1][1][0] = 12;
        nuc2aa[1][1][1] = 12;
        nuc2aa[1][1][2] = 12;
        nuc2aa[1][1][3] = 12;
        // S
        nuc2aa[2][1][0] = 13;
        nuc2aa[2][1][1] = 13;
        nuc2aa[2][1][2] = 13;
        nuc2aa[2][1][3] = 13;
        nuc2aa[0][3][2] = 13;
        nuc2aa[0][3][1] = 13;
        // T
        nuc2aa[0][1][0] = 14;
        nuc2aa[0][1][1] = 14;
        nuc2aa[0][1][2] = 14;
        nuc2aa[0][1][3] = 14;
        // Stop
        nuc2aa[2][0][0] = 15;
        nuc2aa[2][3][0] = 15;
        nuc2aa[2][0][3] = 15;
        // Triplet code with N's
        for (int i = 0; i < 8; i++) {
            for (int j = 0; j < 8; j++) {
                nuc2aa[7][i][j] = -1;
                nuc2aa[i][7][j] = -1;
                nuc2aa[i][j][7] = -1;
            }
        }

        // Codon encoding to store which DNA triplet code is used.
        // This encoding is used for indexing hamming distance look-up table.
        for (int i = 0; i < 4; i++) {
            for (int i2 = 0; i2 < 4; i2++) {
                nuc2num[i][i2][0] = 0;
                nuc2num[i][i2][1] = 1;
                nuc2num[i][i2][2] = 2;
                nuc2num[i][i2][3] = 3;
            }
        }
        // for Arg
        nuc2num[0][3][3] = 7; //AGG
        nuc2num[0][3][0] = 4; //AGA
        // for Leu
        nuc2num[2][2][3] = 7; //TTG
        nuc2num[2][2][0] = 4; //TTA
        nuc2num[0][2][3] = 8; //ATG
        // for Ser
        nuc2num[0][3][2] = 10; //AGT
        nuc2num[0][3][1] = 9; //AGC
        // for FYW
        nuc2num[2][0][1] = 5;
        nuc2num[2][0][2] = 6;
        nuc2num[2][3][3] = 7;
        // for IV
        nuc2num[0][2][0] = 4;
        nuc2num[0][2][1] = 5;
        nuc2num[0][2][2] = 6;
        // for QE
        nuc2num[3][0][0] = 4; //GAA
        nuc2num[3][0][3] = 7; //GAG
        // for stop codon
        nuc2num[2][3][0] = 4; //TGA
    }
}

// It translates a DNA sequence to amino acid sequences with six frames
void SeqIterator::sixFrameTranslation(const char *seq) {
    for (int i = 0; i < 6; i++) { aaFrames[i].clear(); }

    int len = strlen(seq);
    size_t end = len - 1;
    // Translate DNA to AA.
    if (len % 3 == 2) {
        for (int i = 0; i < len - 4; i = i + 3) {
            aaFrames[0].push_back(nuc2aa[nuc2int(atcg[seq[i]])][nuc2int(atcg[seq[i + 1]])][nuc2int(atcg[seq[i + 2]])]);
            aaFrames[1].push_back(
                    nuc2aa[nuc2int(atcg[seq[i + 1]])][nuc2int(atcg[seq[i + 2]])][nuc2int(atcg[seq[i + 3]])]);
            aaFrames[2].push_back(
                    nuc2aa[nuc2int(atcg[seq[i + 2]])][nuc2int(atcg[seq[i + 3]])][nuc2int(atcg[seq[i + 4]])]);
            aaFrames[3].push_back(
                    nuc2aa[nuc2int(iRCT[atcg[seq[end - (i + 0)]]])][nuc2int(iRCT[atcg[seq[end - (i + 1)]]])][
                            nuc2int(iRCT[atcg[seq[end - (i + 2)]]])]);
            aaFrames[4].push_back(
                    nuc2aa[nuc2int(iRCT[atcg[seq[end - (i + 1)]]])][nuc2int(iRCT[atcg[seq[end - (i + 2)]]])][
                            nuc2int(iRCT[atcg[seq[end - (i + 3)]]])]);
            aaFrames[5].push_back(
                    nuc2aa[nuc2int(iRCT[atcg[seq[end - (i + 2)]]])][nuc2int(iRCT[atcg[seq[end - (i + 3)]]])][
                            nuc2int(iRCT[atcg[seq[end - (i + 4)]]])]);
        }
    } else if (len % 3 == 1) {
        for (int i = 0; i < len - 4; i = i + 3) {
            aaFrames[0].push_back(nuc2aa[nuc2int(atcg[seq[i]])][nuc2int(atcg[seq[i + 1]])][nuc2int(atcg[seq[i + 2]])]);
            aaFrames[1].push_back(
                    nuc2aa[nuc2int(atcg[seq[i + 1]])][nuc2int(atcg[seq[i + 2]])][nuc2int(atcg[seq[i + 3]])]);
            aaFrames[2].push_back(
                    nuc2aa[nuc2int(atcg[seq[i + 2]])][nuc2int(atcg[seq[i + 3]])][nuc2int(atcg[seq[i + 4]])]);
            aaFrames[3].push_back(
                    nuc2aa[nuc2int(iRCT[atcg[seq[end - (i + 0) - 3]]])][nuc2int(iRCT[atcg[seq[end - (i + 1) - 3]]])][
                            nuc2int(iRCT[atcg[seq[end - (i + 2) - 3]]])]);
            aaFrames[4].push_back(
                    nuc2aa[nuc2int(iRCT[atcg[seq[end - (i + 1) - 3]]])][nuc2int(iRCT[atcg[seq[end - (i + 2) - 3]]])][
                            nuc2int(iRCT[atcg[seq[end - (i + 3) - 3]]])]);
            aaFrames[5].push_back(
                    nuc2aa[nuc2int(iRCT[atcg[seq[end - (i + 2)]]])][nuc2int(iRCT[atcg[seq[end - (i + 3)]]])][
                            nuc2int(iRCT[atcg[seq[end - (i + 4)]]])]);
        }
    } else {
        for (int i = 0; i < len - 4; i = i + 3) {
            aaFrames[0].push_back(nuc2aa[nuc2int(atcg[seq[i]])][nuc2int(atcg[seq[i + 1]])][nuc2int(atcg[seq[i + 2]])]);
            aaFrames[1].push_back(
                    nuc2aa[nuc2int(atcg[seq[i + 1]])][nuc2int(atcg[seq[i + 2]])][nuc2int(atcg[seq[i + 3]])]);
            aaFrames[2].push_back(
                    nuc2aa[nuc2int(atcg[seq[i + 2]])][nuc2int(atcg[seq[i + 3]])][nuc2int(atcg[seq[i + 4]])]);
            aaFrames[3].push_back(
                    nuc2aa[nuc2int(iRCT[atcg[seq[end - (i + 0) - 3]]])][nuc2int(iRCT[atcg[seq[end - (i + 1) - 3]]])][
                            nuc2int(iRCT[atcg[seq[end - (i + 2) - 3]]])]);
            aaFrames[4].push_back(
                    nuc2aa[nuc2int(iRCT[atcg[seq[end - (i + 1)]]])][nuc2int(iRCT[atcg[seq[end - (i + 2)]]])][
                            nuc2int(iRCT[atcg[seq[end - (i + 3)]]])]);
            aaFrames[5].push_back(
                    nuc2aa[nuc2int(iRCT[atcg[seq[end - (i + 2)]]])][nuc2int(iRCT[atcg[seq[end - (i + 3)]]])][
                            nuc2int(iRCT[atcg[seq[end - (i + 4)]]])]);
        }
    }
}

void SeqIterator::fillQueryKmerBuffer(const char *seq, int seqLen, QueryKmerBuffer &kmerBuffer, size_t &posToWrite, const int &seqID,
                                 uint32_t offset) {
    int forOrRev;
    uint64_t tempKmer = 0;
    int checkN;

    for (uint32_t frame = 0; frame < 6; frame++) {
        uint32_t len = aaFrames[frame].size();
        forOrRev = frame / 3;
        for (uint32_t kmerCnt = 0; kmerCnt < len - kmerLength - spaceNum + 1; kmerCnt++) {
            // Convert amino acid 8-mer to an integer
            tempKmer = 0;
            checkN = 0;
            for (uint32_t i = 0, j = 0; i < kmerLength + spaceNum; i++, j += mask[i]) { // j:0-7, i:0-11
                if (-1 == aaFrames[frame][kmerCnt + i]) {
                    checkN = 1;
                    break;
                }
                tempKmer += aaFrames[frame][kmerCnt + i] * powers[j] * mask[i];
            }

            // Add DNA info
            if (checkN == 1) {
                kmerBuffer.buffer[posToWrite] = {UINT64_MAX, 0, 0, frame};
            } else {
                addDNAInfo_QueryKmer(tempKmer, seq, forOrRev, kmerCnt, frame, seqLen);
                if (forOrRev == 0) {
                    kmerBuffer.buffer[posToWrite] = {tempKmer, seqID, (frame % 3) + (kmerCnt * 3) + offset, frame};
                } else {
                    if ((seqLen % 3 == 1 && (frame == 3 || frame == 4)) || (seqLen % 3 == 0 && frame == 3)) {
                        kmerBuffer.buffer[posToWrite] = {tempKmer, seqID,
                                                         seqLen - ((frame % 3) + (kmerCnt * 3)) - 24 + offset - 3,
                                                         frame};
                    } else {
                        kmerBuffer.buffer[posToWrite] = {tempKmer, seqID,
                                                         seqLen - ((frame % 3) + (kmerCnt * 3)) - 24 + offset, frame};
                    }
                }
            }
            posToWrite++;
        }
    }
}

void
SeqIterator::addDNAInfo_QueryKmer(uint64_t &kmer, const char *seq, int forOrRev, const int &kmerCnt, const int &frame,
                                  int seqLen) {
    int start = (frame % 3) + (kmerCnt * 3);
    kmer <<= bitsFor8Codons;
    size_t end = seqLen - 1;

    if (forOrRev == 0) { // Forward
        for (int i = 0, j = 0; i < kmerLength + spaceNum_int; i++, j += mask_int[i]) {
            kmer |= (nuc2num[nuc2int(atcg[seq[start + i * 3]])][nuc2int(atcg[seq[start + i * 3 + 1]])][
                    nuc2int(atcg[seq[start + i * 3 + 2]])] * mask[i]) << (j * bitsForCodon);
        }
    } else { // Reverse
        if ((seqLen % 3 == 1 && (frame == 3 || frame == 4)) || (seqLen % 3 == 0 && frame == 3)) {
            for (int i = 0, j = 0; i < kmerLength + spaceNum_int; i ++, j += mask_int[i]) {
                kmer |= (nuc2num[nuc2int(iRCT[atcg[seq[end - (start + i*3) - 3]]])][
                        nuc2int(iRCT[atcg[seq[end - (start + i*3 + 1) - 3]]])][
                        nuc2int(iRCT[atcg[seq[end - (start + i*3 + 2) - 3]]])] * mask[i]) << (j * bitsForCodon);
            }
        } else {
            for (int i = 0, j = 0; i < kmerLength + spaceNum_int; i ++, j += mask_int[i]) {
                kmer |= (nuc2num[nuc2int(iRCT[atcg[seq[end - (start + i*3)]]])][
                        nuc2int(iRCT[atcg[seq[end - (start + i*3 + 1)]]])][
                        nuc2int(iRCT[atcg[seq[end - (start + i*3 + 2)]]])] * mask[i]) << (j * bitsForCodon);
            }
        }
    }
}

bool SeqIterator::translateBlock(const char *seq, PredictedBlock block) {
    aaFrames[0].clear();
    if(aaFrames->capacity() < strlen(seq) / 3 + 1) {
        aaFrames->reserve(strlen(seq) / 3 + 1);
    }
    if (block.strand == 1) {
        for (int i = block.start; i + 2 <= block.end; i = i + 3) {
            aaFrames[0].push_back(nuc2aa[nuc2int(atcg[seq[i]])][nuc2int(atcg[seq[i + 1]])][nuc2int(atcg[seq[i + 2]])]);
        }
    } else {
        for (int i = block.end; i >= (int) block.start + 2; i = i - 3) {
            aaFrames[0].push_back(nuc2aa[nuc2int(iRCT[atcg[seq[i]]])][nuc2int(iRCT[atcg[seq[i - 1]]])][
                                          nuc2int(iRCT[atcg[seq[i - 2]]])]);
        }
    }
    return true;
}

string SeqIterator::reverseCompliment(string &read) const {
    int len = read.length();
    string out;
    for (int i = 0; i < len; i++) {
        out.push_back(iRCT[read[i]]);
    }
    reverse(out.begin(), out.end());
    return out;
}

char *SeqIterator::reverseCompliment(char *read, int length) const {
    char *revCom = (char *) malloc(sizeof(char) * (length + 1));
    for (int i = 0; i < length; i++) {
        revCom[length - i - 1] = iRCT[read[i]];
    }
    revCom[length] = '\0';
    return revCom;
}


// It extracts kmers from amino acid sequence with DNA information and fill the kmerBuffer with them.
int
SeqIterator::fillBufferWithKmerFromBlock(const PredictedBlock &block, const char *seq, TargetKmerBuffer &kmerBuffer,
                                         size_t &posToWrite, const uint32_t &seqID, int taxIdAtRank) {
    uint64_t tempKmer = 0;
    int len = (int) aaFrames[0].size();
    int checkN;
    for (int kmerCnt = 0; kmerCnt < len - kmerLength - spaceNum_int + 1; kmerCnt++) {
        tempKmer = 0;
        checkN = 0;
        for (uint32_t i = 0, j = 0; i < kmerLength + spaceNum; i++, j += mask[i]) {
            if (-1 == aaFrames[0][kmerCnt + i]) {
                checkN = 1;
                break;
            }
            tempKmer += aaFrames[0][kmerCnt + i] * powers[j] * mask[i];
        }
        if (checkN == 1) {
            kmerBuffer.buffer[posToWrite] = {UINT64_MAX, -1, 0, false};
        } else {
            addDNAInfo_TargetKmer(tempKmer, seq, block, kmerCnt);
            if(posToWrite >= kmerBuffer.bufferSize - 2) {
                cout << "HERE " << posToWrite << endl;
                return -1;
            }
            kmerBuffer.buffer[posToWrite] = {tempKmer, taxIdAtRank, seqID, false};
        }
        posToWrite++;
    }
    return 0;
}

// It adds DNA information to kmers referring the original DNA sequence.
void
SeqIterator::addDNAInfo_TargetKmer(uint64_t &kmer, const char *seq, const PredictedBlock &block, const int &kmerCnt) {
    kmer <<= bitsFor8Codons;
    if (block.strand == 1) {
        int start = block.start + (kmerCnt * 3);
        for (int i = 0, j = 0; i < kmerLength + spaceNum_int; i ++, j += mask_int[i]) {
            kmer |= (nuc2num[nuc2int(atcg[seq[start + i*3]])][nuc2int(atcg[seq[start + i*3 + 1]])][
                    nuc2int(atcg[seq[start + i*3 + 2]])] * mask[i]) << (j * bitsForCodon);
        }
    } else {
        int start = block.end - (kmerCnt * 3);
        for (int i = 0, j = 0; i < kmerLength + spaceNum_int; i++, j += mask_int[i]) {
            kmer |= (nuc2num[nuc2int(iRCT[atcg[seq[start - i*3]]])][nuc2int(iRCT[atcg[seq[start - i*3 - 1]]])][
                    nuc2int(iRCT[atcg[seq[start - i*3 - 2]]])] * mask[i]) << (j * bitsForCodon);
        }
    }
}

size_t SeqIterator::kmerNumOfSixFrameTranslation(const char *seq) {
    size_t len = strlen(seq);
    return (len / 3 - kmerLength) * 6;
}

size_t SeqIterator::getNumOfKmerForBlock(const PredictedBlock &block) {
    size_t len = block.end - block.start + 1;
    return len / 3 - 7;
}

// It makes the blocks for translation
// Each block has a predicted gene part and an intergenic region. When another gene shows up, new block starts.
void SeqIterator::getExtendedORFs(struct _gene *genes, struct _node *nodes, vector<PredictedBlock> &blocks,
                                       size_t numOfGene, size_t length,
                                       size_t &blockIdx, vector<uint64_t> &intergenicKmerList, const char *seq) {

    //Exceptional case 1: 0 prdicted gene
    if (numOfGene == 0) {
        blocks.emplace_back(0, length - 1, 1);
        blockIdx++;
        return;
    }

    //Exceptional case 2: Only 1 prdicted gene
    int frame;
    int rightEnd = 0;
    int leftEnd = 0;
    if (numOfGene == 1) {
        if (nodes[genes[0].start_ndx].strand == 1) { //forward
            frame = (genes[0].begin - 1) % 3;
            leftEnd = 0;
            while (leftEnd % 3 != frame) leftEnd++; // y - (y - x) % 3 .. which would be faster?
            blocks.emplace_back(leftEnd, length - 1, 1);
            blockIdx++;
        } else { //reverse
            frame = (genes[0].end - 1) % 3;
            rightEnd = length - 1;
            while (rightEnd % 3 != frame) rightEnd--;
            blocks.emplace_back(0, rightEnd, -1);
            blockIdx++;
        }
        return;
    }


    /* Main routine */

    bool hasBeenExtendedToLeft = false;
    int k = 23;
    char *newIntergenicKmer = (char *) malloc(sizeof(char) * (k + 1));
    char *leftKmer = (char *) malloc(sizeof(char) * (k + 1));
    char *rightKmer = (char *) malloc(sizeof(char) * (k + 1));
    char *leftKmerReverse = (char *) malloc(sizeof(char) * (k + 1));
    char *rightKmerReverse = (char *) malloc(sizeof(char) * (k + 1));
    bool isReverse = false;
    uint64_t leftKmerHash = 0, rightKmerHash = 0;


    // Extend the first gene to cover intergenic regions
    if (nodes[genes[0].start_ndx].strand == 1) { //forward
        frame = (genes[0].begin - 1) % 3;
        leftEnd = 0;
        while (leftEnd % 3 != frame) leftEnd++;
        blocks.emplace_back(leftEnd, genes[1].begin - 1 + 22, 1);
        blockIdx++;
    } else {
        frame = (genes[0].end - 1) % 3;
        rightEnd = genes[1].begin - 1 + 22;
        while (rightEnd % 3 != frame) rightEnd--;
        blocks.emplace_back(0, rightEnd, -1);
        blockIdx++;
    }

    // From the second gene to the second last gene
    for (size_t geneIdx = 1; geneIdx < numOfGene - 1; geneIdx++) {
        isReverse = false;

        // Make two k-mer hash; each from left and right of current gene. They are used for choosing extension direction.
        strncpy(leftKmer, seq + genes[geneIdx].begin - 1 - k, k);
        strncpy(rightKmer, seq + genes[geneIdx].end, k);
        if (nodes[genes[geneIdx].start_ndx].strand == 1) {
            leftKmerHash = XXH64(leftKmer, k, 0);
            rightKmerHash = XXH64(rightKmer, k, 0);
        } else {
            isReverse = true;
            for (int j = k - 1; j >= 0; j--) {
                leftKmerReverse[k - j - 1] = iRCT[leftKmer[j]];
                rightKmerReverse[k - j - 1] = iRCT[rightKmer[j]];
            }
            leftKmerHash = XXH64(leftKmerReverse, k, 0);
            rightKmerHash = XXH64(rightKmerReverse, k, 0);
        }

        // Extend genes to cover intergenic regions
        if (find(intergenicKmerList.begin(), intergenicKmerList.end(), leftKmerHash) !=
            intergenicKmerList.end()) { //Extension to left
            if (!hasBeenExtendedToLeft) {
                if (!isReverse) { // Forward
                    blocks.emplace_back(genes[geneIdx].begin - 1, genes[geneIdx].end - 1, 1);
                    blockIdx++;
                } else { // Reverse
                    blocks.emplace_back(genes[geneIdx].begin - 1, genes[geneIdx].end - 1, -1);
                    blockIdx++;
                }
            } else {
                if (!isReverse) { //forward
                    frame = (genes[geneIdx].begin - 1) % 3;
                    leftEnd = genes[geneIdx - 1].end - 1 - 22;
                    while (leftEnd % 3 != frame) leftEnd++;
                    blocks.emplace_back(leftEnd, genes[geneIdx].end - 1, 1);
                    blockIdx++;
                } else { // reverse
                    blocks.emplace_back(genes[geneIdx - 1].end - 22 - 1, genes[geneIdx].end - 1, -1);
                    blockIdx++;
                }
            }
            hasBeenExtendedToLeft = true;
        } else { // Extension to right
            if (hasBeenExtendedToLeft) {
                if (!isReverse) { //forward
                    frame = (genes[geneIdx].begin - 1) % 3;
                    leftEnd = genes[geneIdx - 1].end - 1 - 22;
                    while (leftEnd % 3 != frame) leftEnd++;
                    blocks.emplace_back(leftEnd, genes[geneIdx + 1].begin - 1 + 22, 1);
                    blockIdx++;
                } else {
                    frame = (genes[geneIdx].end - 1) % 3;
                    rightEnd = genes[geneIdx + 1].begin - 1 + 22;
                    while (rightEnd % 3 != frame) rightEnd--;
                    blocks.emplace_back(genes[geneIdx - 1].end - 1 - 22, rightEnd, -1);
                    blockIdx++;
                }
            } else {
                if (!isReverse) { //forward
                    blocks.emplace_back(genes[geneIdx].begin - 1, genes[geneIdx + 1].begin - 1 + 22, 1);
                    blockIdx++;
                } else {
                    frame = (genes[geneIdx].end - 1) % 3;
                    rightEnd = genes[geneIdx + 1].begin - 1 + 22;
                    while (rightEnd % 3 != frame) rightEnd--;
                    blocks.emplace_back(genes[geneIdx].begin - 1, rightEnd, -1);
                    blockIdx++;
                }
            }
            hasBeenExtendedToLeft = false;

            if (find(intergenicKmerList.begin(), intergenicKmerList.end(), rightKmerHash) == intergenicKmerList.end()) {
                intergenicKmerList.push_back(rightKmerHash);
            }
        }
    }

    //For the last gene
    if (find(intergenicKmerList.begin(), intergenicKmerList.end(), leftKmerHash) !=
        intergenicKmerList.end()) { //extension to left
        if (!isReverse) { //forward
            frame = (genes[numOfGene - 1].begin - 1) % 3;
            leftEnd = genes[numOfGene - 2].end - 1 - 22;
            while (leftEnd % 3 != frame) leftEnd++;
            blocks.emplace_back(leftEnd, length - 1, 1);
            blockIdx++;
        } else { // reverse
            frame = (genes[numOfGene - 1].end - 1) % 3;
            rightEnd = length - 1;
            while (rightEnd % 3 != frame) rightEnd--;
            blocks.emplace_back(genes[numOfGene - 2].end - 22 - 1, rightEnd, -1);
            blockIdx++;
        }
    } else { //extension to right
        if (hasBeenExtendedToLeft) {
            if (!isReverse) { //forward
                frame = (genes[numOfGene - 1].begin - 1) % 3;
                leftEnd = genes[numOfGene - 2].end - 1 - 22;
                while (leftEnd % 3 != frame) leftEnd++;
                blocks.emplace_back(leftEnd, length - 1, 1);
                blockIdx++;
            } else {
                frame = (genes[numOfGene - 1].end - 1) % 3;
                rightEnd = length - 1;
                while (rightEnd % 3 != frame) rightEnd--;
                blocks.emplace_back(genes[numOfGene - 2].end - 22 - 1, rightEnd, -1);
                blockIdx++;
            }
        } else {
            if (!isReverse) {
                blocks.emplace_back(genes[numOfGene - 1].begin, length - 1, 1);
                blockIdx++;
            } else {
                frame = (genes[numOfGene - 1].end - 1) % 3;
                rightEnd = length - 1;
                while (rightEnd % 3 != frame) rightEnd--;
                blocks.emplace_back(genes[numOfGene - 1].begin - 1, rightEnd, -1);
                blockIdx++;
            }
        }

        //If current intergenic sequences is new, update intergenicKmerList.
        if (find(intergenicKmerList.begin(), intergenicKmerList.end(), rightKmerHash) == intergenicKmerList.end()) {
            intergenicKmerList.push_back(rightKmerHash);
        }
    }

    free(newIntergenicKmer);
    free(leftKmer);
    free(rightKmer);
    free(leftKmerReverse);
    free(rightKmerReverse);
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
        cout << "\t";

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
        cout<<"\t";
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
