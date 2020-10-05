//
// Created by KJB on 01/09/2020.
//

#include "SeqAlterator.h"

SeqAlterator::SeqAlterator() {

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

void SeqAlterator::dna2aa(const string & forward, const string & reverse){

    for(int i = 0 ; i < 6 ; i++){ aaFrames[i].clear(); }

    int len = forward.length();

    ///translation from DNA to AA. in each frame
    for(int i = 0; i < len - 4; i = i+3 )
    {
        aaFrames[0].push_back(nuc2aa[nuc2int(forward[i ])][nuc2int(forward[i + 1])][nuc2int(forward[i + 2])]);
        aaFrames[1].push_back(nuc2aa[nuc2int(forward[i + 1])][nuc2int(forward[i + 2])][nuc2int(forward[i + 3])]);
        aaFrames[2].push_back(nuc2aa[nuc2int(forward[i + 2])][nuc2int(forward[i + 3])][nuc2int(forward[i + 4])]);
        aaFrames[3].push_back(nuc2aa[nuc2int(reverse[i ])][nuc2int(reverse[i + 1])][nuc2int(reverse[i + 2])]);
        aaFrames[4].push_back(nuc2aa[nuc2int(reverse[i + 1])][nuc2int(reverse[i + 2])][nuc2int(reverse[i + 3])]);
        aaFrames[5].push_back(nuc2aa[nuc2int(reverse[i + 2])][nuc2int(reverse[i + 3])][nuc2int(reverse[i + 4])]);
    }
    if(len % 3 == 0){
        aaFrames[0].push_back(nuc2aa[nuc2int(forward[len - 3])][nuc2int(forward[len - 2])][nuc2int(forward[len - 1])]);
        aaFrames[3].push_back(nuc2aa[nuc2int(reverse[len - 3])][nuc2int(reverse[len - 2])][nuc2int(reverse[len - 1])]);
    }
    if(len % 3 == 1 ){
        aaFrames[0].push_back(nuc2aa[nuc2int(forward[len - 4])][nuc2int(forward[len - 3])][nuc2int(forward[len - 2])]);
        aaFrames[1].push_back(nuc2aa[nuc2int(forward[len - 3])][nuc2int(forward[len - 2])][nuc2int(forward[len - 1])]);
        aaFrames[3].push_back(nuc2aa[nuc2int(reverse[len - 4])][nuc2int(reverse[len - 3])][nuc2int(reverse[len - 2])]);
        aaFrames[4].push_back(nuc2aa[nuc2int(reverse[len - 3])][nuc2int(reverse[len - 2])][nuc2int(reverse[len - 1])]);
    }
}

void SeqAlterator::dna2aa2(const SeqSegment & seq, const MmapedData<char> & seqFile){
    const size_t & start = seq.start;
    const size_t & end = seq.end;
    for(int i = 0 ; i < 6 ; i++){ aaFrames[i].clear(); }

    size_t len = end - start + 1;
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

ExtractStartPoint SeqAlterator::fillKmerBuffer2(SeqSegment seq, MmapedData<char> & seqFile, Kmer * kmerList, int seqID, size_t & kmerBufferIdx, ExtractStartPoint ESP)
{
    ExtractStartPoint defaultStartPoint = { 0, 0};
    if(ESP.startOfFrame + ESP.startOfFrame != 0){
        defaultStartPoint = {ESP.frame , ESP.startOfFrame};
    }
    uint64_t tempKmer = 0;
    size_t startOfKmer = defaultStartPoint.startOfFrame;
    size_t frame = defaultStartPoint.frame;
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
            kmerList[kmerBufferIdx].ADkmer = tempKmer;
            kmerList[kmerBufferIdx].info.sequenceID = seqID;
            kmerList[kmerBufferIdx].info.pos = startOfKmer;
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

ExtractStartPoint SeqAlterator::fillKmerBuffer(const string * dnaSeq, Kmer * kmerList, int seqID, size_t & kmerBufferIdx, ExtractStartPoint ESP)
{
    ExtractStartPoint defaultStartPoint = { 0, 0};
    if(ESP.startOfFrame + ESP.startOfFrame != 0){
        defaultStartPoint = {ESP.frame , ESP.startOfFrame};
    }
    uint64_t tempKmer = 0;
    size_t startOfKmer = defaultStartPoint.startOfFrame;
    size_t frame = defaultStartPoint.frame;
    int forOrRev;

    for( frame ; frame < 6 ; frame++)
    {
        int len = aaFrames[frame].size();
        forOrRev = frame / 3;
        for (startOfKmer ; startOfKmer < len - kmerLength ; startOfKmer++)
        {
            ///Amino acid 2 number
            for (size_t i = 0; i < kmerLength; i++)
            {
                tempKmer += aaFrames[frame][startOfKmer + i] * powers[i];
            }
            tempKmer = addDNAInfo(tempKmer, dnaSeq[forOrRev], startOfKmer, frame);

            ///memcpy를 써보자
            kmerList[kmerBufferIdx].ADkmer = tempKmer;
            kmerList[kmerBufferIdx].info.sequenceID = seqID;
            kmerList[kmerBufferIdx].info.pos = startOfKmer;
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

uint64_t SeqAlterator::addDNAInfo(uint64_t kmer, const string& read, const int startOfKmer, const int frame)
{
    int start = frame % 3 + (startOfKmer * 3);
    kmer <<= 25;

    for( int i = 0; i < kmerLength*3; i += 3)
    {
        kmer |= nuc2num[nuc2int(read[start + i])][nuc2int(read[start + i + 1])][nuc2int(read[start + i + 2])] << i;
    }
    return kmer;
}

void SeqAlterator::addDNAInfo2(uint64_t & kmer, SeqSegment & seq, MmapedData<char> & seqFile, const int & forOrRev, const int & startOfKmer, const int & frame)
{
    int start = (frame % 3) + (startOfKmer * 3);
    uint64_t temp = kmer;
    kmer <<= 25;

    if(forOrRev == 0){
        for( int i = 0; i < kmerLength * 3; i += 3) {
            kmer |= nuc2num[nuc2int(seqFile.data[seq.start + (start + i)])][nuc2int(seqFile.data[seq.start + (start + i + 1)])][nuc2int(seqFile.data[seq.start + (start + i + 2)])] << i;
        }
    }
    if(forOrRev == 1){
        for( int i = 0; i < kmerLength * 3; i += 3) {
            kmer |= nuc2num[nuc2int(iRCT[seqFile.data[seq.end - (start + i)]])][nuc2int(iRCT[seqFile.data[seq.end - (start + i + 1)]])][nuc2int(iRCT[seqFile.data[seq.end - (start + i + 2)]])] << i;
        }
    }
    return;
}

string SeqAlterator::reverseCompliment(string & read) const
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

void SeqAlterator::getSeqSegments(vector<SeqSegment> & seqSegments, MmapedData<char> seqFile)
{
    size_t start = 0;
    size_t end = 0;
    size_t numOfChar = seqFile.fileSize / sizeof(char);

    SeqSegment temp;
    for(size_t i = 0; i < numOfChar; i++)
    {
        if(seqFile.data[i] == '>')
        {
            end = i - 2;
            temp = {start, end};
            seqSegments.push_back(temp); // the first push_back is a garbage.
            while(seqFile.data[i] != '\n')
            {
                cout<<seqFile.data[i];
                i++;
            }
            cout<<endl;
            start = i + 1;
        }
    }
    temp = {start, numOfChar - 2};
    seqSegments.push_back(temp);
}

