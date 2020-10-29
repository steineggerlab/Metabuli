//
// Created by KJB on 24/09/2020.
//
#include "Mmap.h"
#include "Kmer.h"
#include "common.h"
#include <iostream>
#include <vector>
using namespace std;
void count8mers(const char * fileName) {
    struct MmapedData<char> seqFile = mmapData<char>(fileName);
    size_t maxNuc = seqFile.fileSize/sizeof(char);
    vector<SeqSegment> seqSegments;
    size_t start = 0;
    size_t end = 0;
    //SeqSegment temp;
    for(size_t i = 0; i < maxNuc; i++)
    {
        if(seqFile.data[i] == '>')
        {
            seqSegments.emplace_back(start, i - 2);
            while(seqFile.data[i] != '\n')
            {
                cout<<seqFile.data[i];
                i++;
            }
            cout<<endl;
            start = i + 1;
        }
    }
    seqSegments.emplace_back(start, maxNuc - 2);
    size_t len = 0;
    size_t sum = 0;
    for(int i = 1 ; i<seqSegments.size(); i++)
    {
        len = (seqSegments[i].end - seqSegments[i].start + 1);
        if(len % 3 == 0)
        {
            sum += 6 * ((seqSegments[i].end - seqSegments[i].start + 1) / 3) - 46;
        }
        if(len % 3 == 1)
        {
            sum += 6 * ((seqSegments[i].end - seqSegments[i].start + 1) / 3) - 44;
        }
        if(len % 3 == 2)
        {
            sum += 6 * ((seqSegments[i].end - seqSegments[i].start + 1) / 3) - 42;
        }
        cout<<seqFile.data[seqSegments[i].start]<<" "<<seqFile.data[seqSegments[i].end]<<" "<<seqSegments[i].end - seqSegments[i].start + 1<<" "<<((seqSegments[i].end - seqSegments[i].start + 1) / 3) * 6<<endl;
//        cout<<" "<<(seqSegments[i].end - seqSegments[i].start + 1) % 3<<endl;
    }
    cout<<sum<<endl;
    munmap(seqFile.data, seqFile.fileSize+1);
}