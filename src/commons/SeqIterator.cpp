//
// Created by KJB on 01/09/2020.
//

#include "SeqIterator.h"

SeqIterator::SeqIterator() {

    ///powers
    int pow = 1;
    for(int i = 0 ; i < kmerLength; i++) {
        powers[i] = pow;
        pow *= 21;
    }

    ///Reverse compliment table
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
    ///for stop codon
    nuc2num[2][3][0] = 5;

}

void SeqIterator::dna2aa(const string & seq){
    for(int i = 0 ; i < 6 ; i++){ aaFrames[i].clear(); }

    int len = seq.length();
    size_t end = len - 1;
    ///translation from DNA to AA. in each frame
    for(int i = 0; i < len - 4; i = i+3 )
    {
        aaFrames[0].push_back(nuc2aa[nuc2int(seq[i    ])][nuc2int(seq[i + 1])][nuc2int(seq[i + 2])]);
        aaFrames[1].push_back(nuc2aa[nuc2int(seq[i + 1])][nuc2int(seq[i + 2])][nuc2int(seq[i + 3])]);
        aaFrames[2].push_back(nuc2aa[nuc2int(seq[i + 2])][nuc2int(seq[i + 3])][nuc2int(seq[i + 4])]);

        aaFrames[3].push_back(nuc2aa[nuc2int(iRCT[seq[end - (i + 0)]])][nuc2int(iRCT[seq[end - (i + 1)]])][nuc2int(iRCT[seq[end - (i + 2)]])]);
        aaFrames[4].push_back(nuc2aa[nuc2int(iRCT[seq[end - (i + 1)]])][nuc2int(iRCT[seq[end - (i + 2)]])][nuc2int(iRCT[seq[end - (i + 3)]])]);
        aaFrames[5].push_back(nuc2aa[nuc2int(iRCT[seq[end - (i + 2)]])][nuc2int(iRCT[seq[end - (i + 3)]])][nuc2int(iRCT[seq[end - (i + 4)]])]);
    }
    if(len % 3 == 0){
        aaFrames[0].push_back(nuc2aa[nuc2int(seq[end - 2])][nuc2int(seq[end - 1])][nuc2int(seq[end])]);
        aaFrames[3].push_back(nuc2aa[nuc2int(iRCT[seq[2]])][nuc2int(iRCT[seq[1]])][nuc2int(iRCT[seq[0]])]);
    }
    if(len % 3 == 1 ){
        aaFrames[0].push_back(nuc2aa[nuc2int(seq[end - 3])][nuc2int(seq[end - 2])][nuc2int(seq[end - 1])]);
        aaFrames[1].push_back(nuc2aa[nuc2int(seq[end - 2])][nuc2int(seq[end - 1])][nuc2int(seq[end])]);
        aaFrames[3].push_back(nuc2aa[nuc2int(iRCT[seq[3]])][nuc2int(iRCT[seq[2]])][nuc2int(iRCT[seq[1]])]);
        aaFrames[4].push_back(nuc2aa[nuc2int(iRCT[seq[2]])][nuc2int(iRCT[seq[1]])][nuc2int(iRCT[seq[0]])]);
    }
}

void SeqIterator::dna2aa2(const Sequence & seq, const MmapedData<char> & seqFile){
    const size_t & start = seq.start;
    const size_t & end = seq.end;
    size_t len = end - start + 1;
    for(int i = 0 ; i < 6 ; i++){
        aaFrames[i].clear();
        aaFrames[i].reserve(len / 3 + 1);
    }

    ///translation from DNA to AA. in each frame
    for(size_t i = 0; i < len - 4; i = i+3 )
    {
        aaFrames[0].push_back(nuc2aa[nuc2int(seqFile.data[i + start    ])][nuc2int(seqFile.data[i + start + 1])][nuc2int(seqFile.data[i + start + 2])]);
        aaFrames[1].push_back(nuc2aa[nuc2int(seqFile.data[i + start + 1])][nuc2int(seqFile.data[i + start + 2])][nuc2int(seqFile.data[i + start + 3])]);
        aaFrames[2].push_back(nuc2aa[nuc2int(seqFile.data[i + start + 2])][nuc2int(seqFile.data[i + start + 3])][nuc2int(seqFile.data[i + start + 4])]);

        aaFrames[3].push_back(nuc2aa[nuc2int(iRCT[seqFile.data[end - (i + 0)]])][nuc2int(iRCT[seqFile.data[end - (i + 1)]])][nuc2int(iRCT[seqFile.data[end - (i + 2)]])]);
        aaFrames[4].push_back(nuc2aa[nuc2int(iRCT[seqFile.data[end - (i + 1)]])][nuc2int(iRCT[seqFile.data[end - (i + 2)]])][nuc2int(iRCT[seqFile.data[end - (i + 3)]])]);
        aaFrames[5].push_back(nuc2aa[nuc2int(iRCT[seqFile.data[end - (i + 2)]])][nuc2int(iRCT[seqFile.data[end - (i + 3)]])][nuc2int(iRCT[seqFile.data[end - (i + 4)]])]);

    }
    if(len % 3 == 0){
        aaFrames[0].push_back(nuc2aa[nuc2int(seqFile.data[len - 3])][nuc2int(seqFile.data[len - 2])][nuc2int(seqFile.data[len - 1])]);
        aaFrames[3].push_back(nuc2aa[nuc2int(iRCT[seqFile.data[start + 2]])][nuc2int(iRCT[seqFile.data[start + 1]])][nuc2int(iRCT[seqFile.data[start]])]);
    }
    if(len % 3 == 1 ){
        aaFrames[0].push_back(nuc2aa[nuc2int(seqFile.data[len - 4])][nuc2int(seqFile.data[len - 3])][nuc2int(seqFile.data[len - 2])]);
        aaFrames[1].push_back(nuc2aa[nuc2int(seqFile.data[len - 3])][nuc2int(seqFile.data[len - 2])][nuc2int(seqFile.data[len - 1])]);
        aaFrames[3].push_back(nuc2aa[nuc2int(iRCT[seqFile.data[start + 3]])][nuc2int(iRCT[seqFile.data[start + 2]])][nuc2int(iRCT[seqFile.data[start + 1]])]);
        aaFrames[4].push_back(nuc2aa[nuc2int(iRCT[seqFile.data[start + 2]])][nuc2int(iRCT[seqFile.data[start + 1]])][nuc2int(iRCT[seqFile.data[start + 0]])]);
    }
}

ExtractStartPoint SeqIterator::fillKmerBuffer(const string & seq, Kmer * kmerList, int seqID, size_t & kmerBufferIdx, ExtractStartPoint ESP){
    ExtractStartPoint defaultStartPoint = { 0, 0};
    if(ESP.startOfFrame + ESP.startOfFrame != 0){
        defaultStartPoint = {ESP.frame , ESP.startOfFrame};
    }
    uint64_t tempKmer = 0;
    uint32_t startOfKmer = defaultStartPoint.startOfFrame;
    uint32_t frame = defaultStartPoint.frame;
    int forOrRev;

    for( frame ; frame < 6 ; frame++)
    {
        int len = aaFrames[frame].size();
        forOrRev = frame / 3;
        for (startOfKmer ; startOfKmer < len - kmerLength + 1 ; startOfKmer++)
        {
            ///Amino acid 2 number
            for (size_t i = 0; i < kmerLength; i++)
            {
                tempKmer += aaFrames[frame][startOfKmer + i] * powers[i];
            }
            addDNAInfo(tempKmer, seq, forOrRev, startOfKmer, frame);

            ///memcpy를 써보자
            kmerList[kmerBufferIdx] = {tempKmer, seqID, startOfKmer*(frame+1), 0};
            kmerBufferIdx++;

            if(kmerBufferIdx == kmerBufSize)
            {
                ExtractStartPoint ESP = {frame, startOfKmer + 1};
                return ESP;
            }

            tempKmer = 0;
        }
        startOfKmer = 0;
    }

    ExtractStartPoint zero{0,0};
    return zero;


}

uint64_t SeqIterator::addDNAInfo(uint64_t & kmer, const string& seq, int forOrRev, const int startOfKmer, const int frame)
{
    int start = (frame % 3) + (startOfKmer * 3);
    kmer <<= 25;
    size_t end = seq.size() - 1;
    if(forOrRev == 0){
        for( int i = 0; i < kmerLength * 3; i += 3) {
            kmer |= nuc2num[nuc2int(seq[start + i])][nuc2int(seq[start + i + 1])][nuc2int(seq[start + i + 2])] << i;
        }
    } else{
        for( int i = 0; i < kmerLength * 3; i += 3) {
            kmer |= nuc2num[nuc2int(iRCT[seq[end - (start + i)]])][nuc2int(iRCT[seq[end - (start + i + 1)]])][nuc2int(iRCT[seq[end - (start + i + 2)]])] << i;
        }
    }
}

ExtractStartPoint SeqIterator::fillKmerBuffer2(Sequence seq, MmapedData<char> & seqFile, Kmer * kmerList, int seqID, size_t & kmerBufferIdx, ExtractStartPoint ESP)
{
    ExtractStartPoint defaultStartPoint = { 0, 0};
    if(ESP.startOfFrame + ESP.startOfFrame != 0){
        defaultStartPoint = {ESP.frame , ESP.startOfFrame};
    }
    uint64_t tempKmer = 0;
    uint32_t startOfKmer = defaultStartPoint.startOfFrame;
    uint32_t frame = defaultStartPoint.frame;
    int forOrRev;

    for( frame ; frame < 6 ; frame++)
    {
        int len = aaFrames[frame].size();
        forOrRev = frame / 3;
        for (startOfKmer ; startOfKmer < len - kmerLength + 1 ; startOfKmer++)
        {
            ///Amino acid 2 number
            for (size_t i = 0; i < kmerLength; i++)
            {
                tempKmer += aaFrames[frame][startOfKmer + i] * powers[i];
            }
            addDNAInfo2(tempKmer, seq, seqFile, forOrRev, startOfKmer, frame);

            ///memcpy를 써보자
            kmerList[kmerBufferIdx] = {tempKmer, seqID, startOfKmer*(frame+1), 0};
            kmerBufferIdx++;

            if(kmerBufferIdx == kmerBufSize)
            {
                ExtractStartPoint ESP = {frame, startOfKmer + 1};
                return ESP;
            }

            tempKmer = 0;
        }
        startOfKmer = 0;
    }

    ExtractStartPoint zero{0,0};
    return zero;
}

void SeqIterator::addDNAInfo2(uint64_t & kmer, Sequence & seq, MmapedData<char> & seqFile, const int & forOrRev, const int & startOfKmer, const int & frame)
{
    int start = (frame % 3) + (startOfKmer * 3);
    kmer <<= 25;

    if(forOrRev == 0){
        for( int i = 0; i < kmerLength * 3; i += 3) {
            kmer |= nuc2num[nuc2int(seqFile.data[seq.start + (start + i)])][nuc2int(seqFile.data[seq.start + (start + i + 1)])][nuc2int(seqFile.data[seq.start + (start + i + 2)])] << i;
        }
    } else{
        for( int i = 0; i < kmerLength * 3; i += 3) {
            kmer |= nuc2num[nuc2int(iRCT[seqFile.data[seq.end - (start + i)]])][nuc2int(iRCT[seqFile.data[seq.end - (start + i + 1)]])][nuc2int(iRCT[seqFile.data[seq.end - (start + i + 2)]])] << i;
        }
    }
    return;
}

string SeqIterator::reverseCompliment(string & read) const
{
    int len = read.length();
    string out;
    for(int i = 0; i < len; i++)
    {
        out.push_back(iRCT[read[i]]);
    }
    reverse(out.begin(),out.end());
    return out;
}

void SeqIterator::getSeqSegmentsWithoutHead(vector<Sequence> & seqSegments, MmapedData<char> seqFile)
{
    size_t start = 0;
    size_t numOfChar = seqFile.fileSize / sizeof(char);

    //Sequence temp(0,0);
    for(size_t i = 0; i < numOfChar; i++)
    {
        if(seqFile.data[i] == '>')
        {
            seqSegments.emplace_back(start, i-2, i - start - 1);// the first push_back is a garbage.
            while(seqFile.data[i] != '\n')
            {
                i++;
            }
            start = i + 1;
        }
    }
    seqSegments.emplace_back(start, numOfChar - 2, numOfChar - start - 1);
}

void SeqIterator::getSeqSegmentsWithHead(vector<Sequence> & seqSegments, MmapedData<char> seqFile)
{
    size_t start = 0;
    size_t numOfChar = seqFile.fileSize / sizeof(char);

    for(size_t i = 1; i < numOfChar; i++)
    {
        if(seqFile.data[i] == '>')
        {
            seqSegments.emplace_back(start, i-2, i - start - 1);// the first push_back is a garbage.
            start = i;
        }
    }
    seqSegments.emplace_back(start, numOfChar - 2, numOfChar - start - 1);
}

size_t SeqIterator::whatNameWouldBeGood(KmerBuffer & kmerBuffer, MmapedData<char> & seqFile, vector<Sequence> & seqs, bool * checker, size_t & processedSeqCnt) {
#pragma omp parallel
{
    SeqIterator seqIterator;
    size_t posToWrite;
    bool hasOverflow = false;
#pragma omp for schedule(dynamic, 1)
    for (size_t i = 0; i < seqs.size(); i++) {
        if(checker[i] == false && !hasOverflow) {
            kseq_buffer_t buffer(const_cast<char *>(&seqFile.data[seqs[i].start]), seqs[i].length);
            kseq_t *seq = kseq_init(&buffer);
            kseq_read(seq);
            seqs[i].length = strlen(seq->seq.s);
            seqIterator.dna2aa(seq->seq.s);
            size_t kmerCnt = getNumOfKmerForSeq(seq->seq.s);
            posToWrite = kmerBuffer.reserveMemory(kmerCnt);
            if (posToWrite + kmerCnt < kmerBufSize) {
                seqIterator.fillKmerBuffer3(seq->seq.s, kmerBuffer, posToWrite, i);
                checker[i] = true;
                processedSeqCnt ++;
            } else{
                    hasOverflow = true;
            }
        }

    }
}
}
size_t SeqIterator::whatNameWouldBeGoodWithFramePrediction(KmerBuffer & kmerBuffer, MmapedData<char> & seqFile, vector<Sequence> & seqs, bool * checker, size_t & processedSeqCnt) {
int z = 0;
int y = 0;
    omp_set_num_threads(1);
#pragma omp parallel
{
        ProdigalWrapper prodigal;
        SeqIterator seqIterator;
        size_t posToWrite;
        bool hasOverflow = false;
        PredictedBlock * frames;
        #pragma omp for schedule(dynamic, 1)
        for (size_t i = 0; i < seqs.size(); i++) {
            if(checker[i] == false && !hasOverflow) {
                kseq_buffer_t buffer(const_cast<char *>(&seqFile.data[seqs[i].start]), seqs[i].length);
                kseq_t *seq = kseq_init(&buffer);
                kseq_read(seq);
                cout<<"before prodigal"<<endl;
                prodigal.trainASpecies(seq->seq.s);
                prodigal.getPredictedFrames(seq->seq.s);
                cout<<"after"<<endl;
                int overlap = 0;

                vector<PredictedBlock> reverse;
                frames = (PredictedBlock*)malloc((prodigal.getNumberOfPredictedGenes() + 1) * sizeof(PredictedBlock));
                getTranslationBlocks(prodigal.genes, prodigal.nodes, frames, prodigal.getNumberOfPredictedGenes(),
                                     strlen(seq->seq.s));

                for(size_t i = 0; i < 1000 ; i++) {


                    if((prodigal.genes[i].end - prodigal.genes[i].begin + 1) % 3 != 0){
                        for(size_t j = prodigal.genes[i].begin ; j < prodigal.genes[i].end + 1 ; j++  ){
                            cout<<seq->seq.s[j];
                        } cout<<endl;
                        cout << prodigal.genes[i].begin << " " << prodigal.genes[i].end <<" "<< prodigal.nodes[prodigal.genes[i].start_ndx].type <<" "<< prodigal.nodes[prodigal.genes[i].stop_ndx].type<<endl;
                        if(prodigal.nodes[prodigal.genes[i].start_ndx].strand == 1){
                            z++;
                        } else{
                            y++;
                        }
                    }
                }
//                cout<<forward.size()<<endl;
//                for(size_t i = 0; i < forward.size() - 1; i ++) {
//                    if (forward[i].end > forward[i + 1].start) {
//                        overlap++;
//                    }
//                }
//                cout<<"overlapl: "<<overlap<<endl;
//                for(size_t i = 0; i < reverse.size() - 1; i ++) {
//                    if (reverse[i].end > reverse[i + 1].start) {
//                        overlap++;
//                    }
//                }
                cout<<"overlap: "<<overlap<<endl;

                        //cout << prodigal.genes[i].begin << " " << prodigal.genes[i].end <<" "<< prodigal.nodes[prodigal.genes[i].start_ndx].strand <<" "<< prodigal.nodes[prodigal.genes[i].stop_ndx].strand<<endl;





                cout<<"ng: "<<prodigal.getNumberOfPredictedGenes()<<endl;
                seqs[i].length = strlen(seq->seq.s);
                seqIterator.dna2aa(seq->seq.s);
                size_t kmerCnt = getNumOfKmerForSeq(seq->seq.s);
                posToWrite = kmerBuffer.reserveMemory(kmerCnt);
                if (posToWrite + kmerCnt < kmerBufSize) {
                    seqIterator.fillKmerBuffer3(seq->seq.s, kmerBuffer, posToWrite, i);
                    checker[i] = true;
                    processedSeqCnt ++;
                } else{
                    hasOverflow = true;
                }
            }
        }

}
cout<< z<<" " <<y<<endl;
}


void SeqIterator::fillKmerBuffer3(const string & seq,  KmerBuffer & kmerBuffer, size_t & posToWrite, const int & seqID)
{
    int tmp = 0;
    int forOrRev;
    uint64_t tempKmer = 0;
    uint32_t maxLength = 0;
    for(int i = 0; i < 6; i++){
        if(aaFrames[i].size() > maxLength)
            maxLength = aaFrames->size();
    }
    maxLength = maxLength - kmerLength + 1;
    for(uint32_t frame = 0 ; frame < 6 ; frame++)
    {
        int len = aaFrames[frame].size();
        forOrRev = frame / 3;
        for (uint32_t startOfKmer = 0 ; startOfKmer < len - kmerLength + 1 ; startOfKmer++)
        {
            ///Amino acid 2 number
            for (size_t i = 0; i < kmerLength; i++)
            {
                tempKmer += aaFrames[frame][startOfKmer + i] * powers[i];
            }
            addDNAInfo(tempKmer, seq, forOrRev, startOfKmer, frame);

            ///memcpy를 써보자
            kmerBuffer.buffer[posToWrite] = {tempKmer, seqID, startOfKmer + (frame * maxLength) , 0};
            tmp++;
            posToWrite++;
            tempKmer = 0;
        }

    }
}

size_t SeqIterator::getNumOfKmerForSeq(const string & seq){
    size_t len = seq.size();
    if(len % 3 == 0)
    {
        return (6 * (len/3) - 46);
    }
    if(len % 3 == 1)
    {
        return (6 * (len/3) - 44);
    }
    if(len % 3 == 2)
    {
        return (6 * (len/3) - 42);
    }
}

void SeqIterator::getTranslationBlocks(struct _gene * genes, struct _node * nodes, PredictedBlock * blocks, size_t numOfGene, size_t length){

    size_t blockIdx = 0;
    int cutNext = 0;
    //for the first frame
    blocks[0].start = 0;
    blocks[0].end = genes[0].begin - 1;
    blocks[0].strand = 1;
    blockIdx++;
    int frame;
    int rightEnd = 0;

    for(size_t geneIdx = 0 ; geneIdx < numOfGene - 1; geneIdx++){

        if(genes[geneIdx].end > genes[geneIdx + 1].end){ //one gene completely includes another gene
            blocks[blockIdx].start = genes[geneIdx].begin;
            blocks[blockIdx].end = genes[geneIdx].end;
            blocks[blockIdx].strand = nodes[genes[geneIdx].begin].strand;
            geneIdx++;
            blockIdx++;
            continue;
        }

        if(nodes[genes[geneIdx].start_ndx].strand == 1){ //forward
            blocks[blockIdx].start = genes[geneIdx].begin;
            blocks[blockIdx].end = genes[geneIdx + 1].begin - 1;
            blocks[blockIdx].strand = nodes[genes[geneIdx].start_ndx].strand;
            blockIdx ++;
        }else{ // reverse
            frame = genes[geneIdx].end % 3;
            rightEnd = genes[geneIdx+1].begin - 1;
            while(rightEnd%3 != frame){
                rightEnd--;
            }
            blocks[blockIdx].start = genes[geneIdx].begin;
            blocks[blockIdx].end = rightEnd;
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
    cout << "frameIdx: " << blockIdx << endl;
}