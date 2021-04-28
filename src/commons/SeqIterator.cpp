//
// Created by KJB on 01/09/2020.
//

#include "SeqIterator.h"

SeqIterator::SeqIterator() {

    ///powers
    size_t pow = 1;
    for(int i = 0 ; i < kmerLength; i++) {
        powers[i] = pow;
        pow *= 21;
    }

    ///Reverse compliment table
    atcg =
            "................................................................"
            ".AGCG..GT..G.CN...ACTG.A.T......agcg..gt..g.cn...actg.a.t......."
            "................................................................"
            "................................................................";

    iRCT =
            "................................................................"
            ".TVGH..CD..M.KN...YSAABW.R.......tvgh..cd..m.kn...ysaabw.r......"
            "................................................................"
            "................................................................";

    ///Codon table
    ///A
    nuc2aa[3][1][0] = 0; nuc2aa[3][1][1] = 0; nuc2aa[3][1][2] = 0 ; nuc2aa[3][1][3] = 0;
    ///R
    nuc2aa[1][3][0] = 1; nuc2aa[1][3][1] = 1; nuc2aa[1][3][2] = 1 ; nuc2aa[1][3][3] = 1; nuc2aa[0][3][0] = 1; nuc2aa[0][3][3] = 1;
    ///N
    nuc2aa[0][0][2] = 2; nuc2aa[0][0][1] = 2;
    ///D
    nuc2aa[3][0][2] = 3; nuc2aa[3][0][1] = 3;
    ///C
    nuc2aa[2][3][2] = 4; nuc2aa[2][3][1] = 4;
    ///Q
    nuc2aa[1][0][0] = 5; nuc2aa[1][0][3] = 5;
    ///E
    nuc2aa[3][0][0] = 6; nuc2aa[3][0][3] = 6;
    ///G
    nuc2aa[3][3][0] = 7; nuc2aa[3][3][1] = 7; nuc2aa[3][3][2] = 7; nuc2aa[3][3][3] =7;
    ///H
    nuc2aa[1][0][2] = 8; nuc2aa[1][0][1] = 8;
    ///I
    nuc2aa[0][2][2] = 9; nuc2aa[0][2][1] = 9; nuc2aa[0][2][0] = 9;
    ///L
    nuc2aa[2][2][0] = 10; nuc2aa[2][2][3] = 10; nuc2aa[1][2][0] = 10; nuc2aa[1][2][1] = 10;nuc2aa[1][2][2] = 10; nuc2aa[1][2][3] =10 ;
    ///K
    nuc2aa[0][0][0] = 11; nuc2aa[0][0][3] = 11;
    ///M
    nuc2aa[0][2][3] = 12;
    ///F
    nuc2aa[2][2][2] = 13; nuc2aa[2][2][1] = 13;
    ///P
    nuc2aa[1][1][0] = 14; nuc2aa[1][1][1] = 14; nuc2aa[1][1][2] = 14; nuc2aa[1][1][3] = 14;
    ///S
    nuc2aa[2][1][0] = 15; nuc2aa[2][1][1] = 15; nuc2aa[2][1][2] = 15; nuc2aa[2][1][3] = 15; nuc2aa[0][3][2] = 15; nuc2aa[0][3][1] = 15;
    ///T
    nuc2aa[0][1][0] = 16; nuc2aa[0][1][1] = 16; nuc2aa[0][1][2] = 16; nuc2aa[0][1][3] = 16;
    ///W
    nuc2aa[2][3][3] = 17;
    ///Y
    nuc2aa[2][0][2] = 18; nuc2aa[2][0][1] = 18;
    ///V
    nuc2aa[3][2][0] = 19; nuc2aa[3][2][1] = 19; nuc2aa[3][2][2] = 19; nuc2aa[3][2][3] = 19;
    ///Stop
    nuc2aa[2][0][0] = 20; nuc2aa[2][3][0] = 20; nuc2aa[2][0][3] = 20;

    ///triplet code with N's
    for(int i = 0; i < 8; i++){
        for(int j = 0; j < 8; j++){
            nuc2aa[7][i][j] = -1;
            nuc2aa[i][7][j] = -1;
            nuc2aa[i][j][7] = -1;
        }
    }

    ///For encoding DNA information in k-mer
    for(int i =0; i < 4 ; i++)
    {
        for(int i2 = 0; i2 < 4 ; i2++)
        {
            nuc2num[i][i2][0] = 0;
            nuc2num[i][i2][1] = 1;
            nuc2num[i][i2][2] = 2;
            nuc2num[i][i2][3] = 3;
        }
    }
    ///for Arg
    nuc2num[0][3][3] = 4;
    nuc2num[0][3][0] = 5;
    ///for Leu
    nuc2num[2][2][3] = 4;
    nuc2num[2][2][0] = 5;
    ///for Ser
    nuc2num[0][3][2] = 6;
    nuc2num[0][3][1] = 7;
    ///for stop codon
    nuc2num[2][3][0] = 5;

}

///It translates a DNA sequence to amino acid sequences with six frames
void SeqIterator::sixFrameTranslation(const char * seq){
    for(int i = 0 ; i < 6 ; i++){ aaFrames[i].clear(); }

    int len = strlen(seq);
    size_t end = len - 1;
    ///translation from DNA to AA. in each frame
    for(int i = 0; i < len - 4; i = i+3 )
    {
        aaFrames[0].push_back(nuc2aa[nuc2int(atcg[seq[i    ]])][nuc2int(atcg[seq[i + 1]])][nuc2int(atcg[seq[i + 2]])]);
        aaFrames[1].push_back(nuc2aa[nuc2int(atcg[seq[i + 1]])][nuc2int(atcg[seq[i + 2]])][nuc2int(atcg[seq[i + 3]])]);
        aaFrames[2].push_back(nuc2aa[nuc2int(atcg[seq[i + 2]])][nuc2int(atcg[seq[i + 3]])][nuc2int(atcg[seq[i + 4]])]);

        aaFrames[3].push_back(nuc2aa[nuc2int(iRCT[atcg[seq[end - (i + 0)]]])][nuc2int(iRCT[atcg[seq[end - (i + 1)]]])][nuc2int(iRCT[atcg[seq[end - (i + 2)]]])]);
        aaFrames[4].push_back(nuc2aa[nuc2int(iRCT[atcg[seq[end - (i + 1)]]])][nuc2int(iRCT[atcg[seq[end - (i + 2)]]])][nuc2int(iRCT[atcg[seq[end - (i + 3)]]])]);
        aaFrames[5].push_back(nuc2aa[nuc2int(iRCT[atcg[seq[end - (i + 2)]]])][nuc2int(iRCT[atcg[seq[end - (i + 3)]]])][nuc2int(iRCT[atcg[seq[end - (i + 4)]]])]);
    }
    if(len % 3 == 0){
        aaFrames[0].push_back(nuc2aa[nuc2int(atcg[seq[end - 2]])][nuc2int(atcg[seq[end - 1]])][nuc2int(atcg[seq[end]])]);
        aaFrames[3].push_back(nuc2aa[nuc2int(iRCT[atcg[seq[2]]])][nuc2int(iRCT[atcg[seq[1]]])][nuc2int(iRCT[atcg[seq[0]]])]);
    }
    if(len % 3 == 1 ){
        aaFrames[0].push_back(nuc2aa[nuc2int(atcg[seq[end - 3]])][nuc2int(atcg[seq[end - 2]])][nuc2int(atcg[seq[end - 1]])]);
        aaFrames[1].push_back(nuc2aa[nuc2int(atcg[seq[end - 2]])][nuc2int(atcg[seq[end - 1]])][nuc2int(atcg[seq[end]])]);
        aaFrames[3].push_back(nuc2aa[nuc2int(iRCT[atcg[seq[3]]])][nuc2int(iRCT[atcg[seq[2]]])][nuc2int(iRCT[atcg[seq[1]]])]);
        aaFrames[4].push_back(nuc2aa[nuc2int(iRCT[atcg[seq[2]]])][nuc2int(iRCT[atcg[seq[1]]])][nuc2int(iRCT[atcg[seq[0]]])]);
    }
}


void SeqIterator::fillQueryKmerBuffer(const char * seq, QueryKmerBuffer & kmerBuffer, size_t & posToWrite, const int & seqID) {

    uint32_t forOrRev;
    uint64_t tempKmer = 0;
    uint32_t seqLen = strlen(seq);
    int checkN;
    for(uint32_t frame = 0 ; frame < 6 ; frame++){
        uint32_t len = aaFrames[frame].size();
        forOrRev = frame / 3;
        for (uint32_t kmerCnt = 0 ; kmerCnt < len - kmerLength + 1 ; kmerCnt++) {
            ///Amino acid 2 number
            tempKmer = 0;
            checkN = 0;
            for (size_t i = 0; i < kmerLength; i++){
                if(-1 == aaFrames[frame][kmerCnt + i]){
                    checkN = 1;
                    break;
                }
                tempKmer += aaFrames[frame][kmerCnt + i] * powers[i];
            }

            if(checkN == 1){
                kmerBuffer.buffer[posToWrite] = {UINT64_MAX, 0, 0, frame};
            }else{
                addDNAInfo_QueryKmer(tempKmer, seq, forOrRev, kmerCnt, frame);
                if(forOrRev == 0) {
                    kmerBuffer.buffer[posToWrite] = {tempKmer, seqID, (frame % 3) + (kmerCnt * 3), frame};
                } else{
                    kmerBuffer.buffer[posToWrite] = {tempKmer, seqID, seqLen - ((frame%3) + (kmerCnt*3)) - 24 , frame};
                }
            }
            posToWrite++;
        }
    }
}

void SeqIterator::addDNAInfo_QueryKmer(uint64_t & kmer, const char * seq, int forOrRev, const int & kmerCnt, const int & frame){
    int start = (frame % 3) + (kmerCnt * 3);
    kmer <<= 25;
    size_t end = strlen(seq) - 1;
    if(forOrRev == 0){
        for( int i = 0; i < kmerLength * 3; i += 3) {
            kmer |= nuc2num[nuc2int(atcg[seq[start + i]])][nuc2int(atcg[seq[start + i + 1]])][nuc2int(atcg[seq[start + i + 2]])] << i;
        }
    } else{
        for( int i = 0; i < kmerLength * 3; i += 3) {
            kmer |= nuc2num[nuc2int(iRCT[atcg[seq[end - (start + i)]]])][nuc2int(iRCT[atcg[seq[end - (start + i + 1)]]])][nuc2int(iRCT[atcg[seq[end - (start + i + 2)]]])] << i;
        }
    }
}

void SeqIterator::translateBlock(const char * seq, PredictedBlock & block){
    aaFrames[0].clear();
    size_t blockLength = block.end - block.start + 1;
    aaFrames[0].reserve(blockLength / 3 + 1);
    if(block.strand == 1){
        for(size_t i = 0; i < blockLength - 2; i = i + 3){
            aaFrames[0].push_back(nuc2aa[nuc2int(atcg[seq[block.start + i]])][nuc2int(atcg[seq[block.start+i+1]])][nuc2int(atcg[seq[block.start+i+2]])]);
        }
    }else{
        for(size_t i = 0; i < blockLength - 2; i = i + 3){
            aaFrames[0].push_back(nuc2aa[nuc2int(iRCT[atcg[seq[block.end - i]]])][nuc2int(iRCT[atcg[seq[block.end-i-1]]])][nuc2int(iRCT[atcg[seq[block.end-i-2]]])]);
        }
    }
}

string SeqIterator::reverseCompliment(string & read) const {
    int len = read.length();
    string out;
    for(int i = 0; i < len; i++)
    {
        out.push_back(iRCT[read[i]]);
    }
    reverse(out.begin(),out.end());
    return out;
}



///It extracts kmers from amino acid sequence with DNA information and fill the kmerBuffer with them.
void SeqIterator::fillBufferWithKmerFromBlock(const PredictedBlock & block, const char * seq, TargetKmerBuffer & kmerBuffer, size_t & posToWrite, const uint32_t & seqID, int taxIdAtRank) {
    uint64_t tempKmer = 0;
    uint32_t len = aaFrames[0].size();
    int checkN;
    for (uint32_t startOfKmer = 0 ; startOfKmer < len - kmerLength + 1 ; startOfKmer++){
        ///Amino acid 2 number
        tempKmer = 0;
        checkN = 0;
        for (size_t i = 0; i < kmerLength; i++){
            if(-1 == aaFrames[0][startOfKmer + i]){
                checkN = 1;
                break;
            }
            tempKmer += aaFrames[0][startOfKmer + i] * powers[i];
        }
        if(checkN == 1){
            kmerBuffer.buffer[posToWrite] = {UINT64_MAX, -1, 0,false};
        }else{
            addDNAInfo_TargetKmer(tempKmer, seq, block, startOfKmer);
            kmerBuffer.buffer[posToWrite] = {tempKmer, taxIdAtRank, seqID, false};
        }
        posToWrite++;
    }
}

///It adds DNA information to kmers referring the original DNA sequence.
void SeqIterator::addDNAInfo_TargetKmer(uint64_t & kmer, const char * seq, const PredictedBlock& block, const int & startOfKmer) {
    kmer <<= 25;
    if(block.strand == 1){
        int start = block.start + (startOfKmer * 3);
        for( int i = 0; i < kmerLength * 3; i += 3) {
            kmer |= nuc2num[nuc2int(atcg[seq[start + i]])][nuc2int(atcg[seq[start + i + 1]])][nuc2int(atcg[seq[start + i + 2]])] << i;
        }
    } else{
        int start = block.end - (startOfKmer * 3);
        for( int i = 0; i < kmerLength * 3; i += 3) {
            kmer |= nuc2num[nuc2int(iRCT[atcg[seq[start - i]]])][nuc2int(iRCT[atcg[seq[start - i - 1]]])][nuc2int(iRCT[atcg[seq[start - i - 2]]])] << i;
        }
    }
}

size_t SeqIterator::kmerNumOfSixFrameTranslation(const string & seq){
    size_t len = seq.size();
    if(len % 3 == 0){
        return (6 * (len/3) - 46);
    }else if(len % 3 == 1){
        return (6 * (len/3) - 44);
    }else if(len % 3 == 2){
        return (6 * (len/3) - 42);
    }
    return 0;
}

size_t SeqIterator::getNumOfKmerForBlock(const PredictedBlock & block){
    size_t len = block.end - block.start + 1;
    return len/3 - 7;
}

///It makes the block for translation from DNA to AA
///Each block has a predicted gene part and intergenic region. When another gene shows up, new block starts.
void SeqIterator::getTranslationBlocks(struct _gene * genes, struct _node * nodes, PredictedBlock * blocks, size_t numOfGene, size_t length, size_t & blockIdx){

    if(numOfGene == 0){
        blocks[blockIdx].start = 0;
        blocks[blockIdx].end = length - 1;
        blocks[blockIdx].strand = 1;
        blockIdx++;
        return;
    }


    //for the first frame
    if(genes[0].begin > 23) {
        blocks[blockIdx].start = 0;
        blocks[blockIdx].end = genes[0].begin - 1;
        blocks[blockIdx].strand = 1;
        blockIdx++;
    }

    int frame;
    int rightEnd = 0;

    for (size_t geneIdx = 0 ; geneIdx < numOfGene - 1; geneIdx++) {
        if(genes[geneIdx].end > genes[geneIdx + 1].end){ //one gene completely includes another gene
            blocks[blockIdx].start = genes[geneIdx].begin;
            blocks[blockIdx].end = genes[geneIdx].end + 21;
            blocks[blockIdx].strand = nodes[genes[geneIdx].start_ndx].strand;
            geneIdx++;
            blockIdx++;
            continue;
        }

        if(nodes[genes[geneIdx].start_ndx].strand == 1){ //forward
            blocks[blockIdx].start = genes[geneIdx].begin;
            blocks[blockIdx].end = genes[geneIdx + 1].begin + 20;
            blocks[blockIdx].strand = nodes[genes[geneIdx].start_ndx].strand;
            blockIdx ++;
        } else { // reverse
            frame = genes[geneIdx].end % 3;
            rightEnd = genes[geneIdx+1].begin - 1;
            while(rightEnd%3 != frame){
                rightEnd--;
            }
            blocks[blockIdx].start = genes[geneIdx].begin;
            blocks[blockIdx].end = rightEnd + 21;
            blocks[blockIdx].strand = nodes[genes[geneIdx].start_ndx].strand;
            blockIdx++;
        }
    }

    //For the last block
    if(nodes[genes[numOfGene - 1].start_ndx].strand == 1){ //forward
        blocks[blockIdx].start = genes[numOfGene - 1].begin;
        blocks[blockIdx].end = length -1;
        blocks[blockIdx].strand = nodes[genes[numOfGene-1].start_ndx].strand;
        blockIdx ++;
    }else{ // reverse
        frame = genes[numOfGene-1].end % 3;
        rightEnd = length - 1;
        while(rightEnd%3 != frame){
            rightEnd--;
        }
        blocks[blockIdx].start = genes[numOfGene - 1].begin;
        blocks[blockIdx].end = rightEnd;
        blocks[blockIdx].strand = nodes[genes[numOfGene-1].start_ndx].strand;
        blockIdx ++;
    }
}

void SeqIterator::getTranslationBlocks2(struct _gene * genes, struct _node * nodes, vector<PredictedBlock> & blocks, size_t numOfGene, size_t length, size_t & blockIdx){

    if(numOfGene == 0){
        blocks.emplace_back(0, length - 1, 1);
        blockIdx++;
        return;
    }


    //for the first frame
    if(genes[0].begin > 23) {
        blocks.emplace_back(0, genes[0].begin - 1, 1);
        blockIdx++;
    }

    int frame;
    int rightEnd = 0;

    for (size_t geneIdx = 0 ; geneIdx < numOfGene - 1; geneIdx++) {
        if(genes[geneIdx].end > genes[geneIdx + 1].end){ //one gene completely includes another gene
            blocks.emplace_back(genes[geneIdx].begin, genes[geneIdx].end + 21, nodes[genes[geneIdx].start_ndx].strand);
            geneIdx++;
            blockIdx++;
            continue;
        }

        if(nodes[genes[geneIdx].start_ndx].strand == 1){ //forward
            blocks.emplace_back(genes[geneIdx].begin, genes[geneIdx + 1].begin + 20, nodes[genes[geneIdx].start_ndx].strand);
            blockIdx ++;
        } else { // reverse
            frame = genes[geneIdx].end % 3;
            rightEnd = genes[geneIdx+1].begin - 1;
            while(rightEnd%3 != frame){
                rightEnd--;
            }
            blocks.emplace_back(genes[geneIdx].begin, rightEnd + 21, nodes[genes[geneIdx].start_ndx].strand);
            blockIdx++;
        }
    }

    //For the last block
    if(nodes[genes[numOfGene - 1].start_ndx].strand == 1){ //forward
        blocks.emplace_back(genes[numOfGene - 1].begin, length -1, nodes[genes[numOfGene-1].start_ndx].strand);
        blockIdx ++;
    }else{ // reverse
        frame = genes[numOfGene-1].end % 3;
        rightEnd = length - 1;
        while(rightEnd%3 != frame){
            rightEnd--;
        }
        blocks.emplace_back(genes[numOfGene - 1].begin, rightEnd, nodes[genes[numOfGene-1].start_ndx].strand);
        blockIdx ++;
    }
}


void SeqIterator::printKmerInDNAsequence(uint64_t kmer) {
    uint64_t copy = kmer;
    kmer >>= 25;
    int quotient;
    int dnaInfo;
    vector<int> aa8mer(8);
    vector<string> dna24mer(8);
    for (int i = 0; i < 8; i++) {
        quotient = kmer / powers[7 - i];
        kmer = kmer - (quotient * powers[7 - i]);
        aa8mer[7 - i] = quotient;
    }

    ///Print Amino Acid 8 mer
//    for(int i  = 0; i<8; i++){
//        cout<<aa8mer[i]<<" ";
//    }
//    cout<<endl;


    for (int i = 0; i < 8; i++) {
        dnaInfo = copy & 7u;
        copy >>= 3;
        switch (aa8mer[i]) {
            case 0: //A
                cout << "A";
                if (dnaInfo == 0) {
                    dna24mer[7 - i] = "GCA";
                } else if (dnaInfo == 1) {
                    dna24mer[7 - i] = "GCC";
                } else if (dnaInfo == 2) {
                    dna24mer[7 - i] = "GCT";
                } else {
                    dna24mer[7 - i] = "GCG";
                }
                break;
            case 1: //R
                cout << "R";
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
                } else {
                    dna24mer[7 - i] = "AGA";
                }
                break;
            case 2: //N
                cout << "N";
                if (dnaInfo == 0) {
                    cout << "FUCK";
                } else if (dnaInfo == 1) {
                    dna24mer[7 - i] = "AAC";
                } else if (dnaInfo == 2) {
                    dna24mer[7 - i] = "AAT";
                } else {
                    cout << "FUCK";
                }
                break;
            case 3: //D
                cout << "D";
                if (dnaInfo == 0) {
                    cout << "FUCK";
                } else if (dnaInfo == 1) {
                    dna24mer[7 - i] = "GAC";
                } else if (dnaInfo == 2) {
                    dna24mer[7 - i] = "GAT";
                } else {
                    cout << "FUCK";
                }
                break;
            case 4: //C
                cout << "C";
                if (dnaInfo == 0) {
                    cout << "FUCK";
                } else if (dnaInfo == 1) {
                    dna24mer[7 - i] = "TGC";
                } else if (dnaInfo == 2) {
                    dna24mer[7 - i] = "TGT";
                } else {
                    cout << "FUCK";
                }
                break;
            case 5: // Q
                cout << "Q";
                if (dnaInfo == 0) {
                    dna24mer[7 - i] = "CAA";
                } else if (dnaInfo == 1) {
                    cout << "FUCK";
                } else if (dnaInfo == 2) {
                    cout << "FUCK";
                } else {
                    dna24mer[7 - i] = "CAG";
                }
                break;
            case 6: //E
                cout << "E";
                if (dnaInfo == 0) {
                    dna24mer[7 - i] = "GAA";
                } else if (dnaInfo == 1) {
                    cout << "FUCK";
                } else if (dnaInfo == 2) {
                    cout << "FUCK" ;
                } else {
                    dna24mer[7 - i] = "GAG";
                }
                break;
            case 7: //G
                cout << "G";
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
                cout << "H";
                if (dnaInfo == 0) {
                    cout << "FUCK";
                } else if (dnaInfo == 1) {
                    dna24mer[7 - i] = "CAC";
                } else if (dnaInfo == 2) {
                    dna24mer[7 - i] = "CAT";
                } else {
                    cout << "FUCK";
                }
                break;
                case 9: //I
                cout << "I";
                if (dnaInfo == 0) {
                    dna24mer[7 - i] = "ATA";
                } else if (dnaInfo == 1) {
                    dna24mer[7 - i] = "ATC";
                } else if (dnaInfo == 2) {
                    dna24mer[7 - i] = "ATT";
                } else {
                    cout << "FUCK";
                }
                break;
                case 10: //L
                cout << "L";
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
                cout << "K";
                if (dnaInfo == 0) {
                    dna24mer[7 - i] = "AAA";
                } else if (dnaInfo == 1) {
                    cout << "FUCK";
                } else if (dnaInfo == 2) {
                    cout << "FUCK";
                } else {
                    dna24mer[7 - i] = "AAG";
                }
                break;
                case 12: // M
                cout << "M";
                if (dnaInfo == 0) {
                    cout << "FUCK";
                } else if (dnaInfo == 1) {
                    cout << "FUCK";
                } else if (dnaInfo == 2) {
                    cout << "FUCK";
                } else {
                    dna24mer[7 - i] = "ATG";
                }
                break;
                case 13://F
                cout << "F";
                if (dnaInfo == 0) {
                    cout << "FUCK" ;
                } else if (dnaInfo == 1) {
                    dna24mer[7 - i] = "TTC";
                } else if (dnaInfo == 2) {
                    dna24mer[7 - i] = "TTT";
                } else {
                    cout << "FUCK";
                }
                break;
                case 14: //P
                cout << "P";
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
                cout << "S";
                if (dnaInfo == 0) {
                    dna24mer[7 - i] = "TCA";
                } else if (dnaInfo == 1) {
                    dna24mer[7 - i] = "TCC";
                } else if (dnaInfo == 2) {
                    dna24mer[7 - i] = "TCT";
                } else if (dnaInfo == 3) {
                    dna24mer[7 - i] = "TCG";
                } else if (dnaInfo == 4) {
                    cout << "FUCK";
                } else if (dnaInfo == 5) {
                    cout << "FUCK";
                } else if (dnaInfo == 6) {
                    dna24mer[7 - i] = "AGT";
                } else {
                    dna24mer[7 - i] = "AGC";
                }
                break;
                case 16: //T
                cout << "T";
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
                cout << "W";
                if (dnaInfo == 0) {
                    cout << "FUCK";
                } else if (dnaInfo == 1) {
                    cout << "FUCK";
                } else if (dnaInfo == 2) {
                    cout << "FUCK";
                } else {
                    dna24mer[7 - i] = "TGG";
                }
                break;
                case 18: //Y
                cout << "Y";
                if (dnaInfo == 0) {
                    cout << "FUCK";
                } else if (dnaInfo == 1) {
                    dna24mer[7 - i] = "TAC";
                } else if (dnaInfo == 2) {
                    dna24mer[7 - i] = "TAT";
                } else {
                    cout << "FUCK";
                }
                break;
                case 19: //V
                cout << "V";
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
                cout << "Z";
                if (dnaInfo == 0) {
                    dna24mer[7 - i] = "TAA";
                } else if (dnaInfo == 1) {
                    cout << "FUCK";
                } else if (dnaInfo == 2) {
                    cout << "FUCK";
                } else if (dnaInfo == 3) {
                    dna24mer[7 - i] = "TAG";
                } else if (dnaInfo == 4) {
                    cout << "FUCK";
                } else {
                    dna24mer[7 - i] = "TGA";
                }
                break;
        }
        //cout<<dnaInfo<<" ";
    }

    cout<<endl;
    for (int i = 0; i < 8; i++) {
        cout << dna24mer[7-i];
    }
    cout << endl;
    return;
}
