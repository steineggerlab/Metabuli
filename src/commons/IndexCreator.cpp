//
// Created by KJB on 01/09/2020.
//

#include "IndexCreator.h"

IndexCreator::IndexCreator()
{
    seqIterator = new SeqIterator();

}

IndexCreator::~IndexCreator() {
    delete seqIterator;
}

///It reads a reference sequence file and write a differential index file and target k-mer info file.
void IndexCreator::startIndexCreatingParallel(const char * seqFileName, const char * outputFileName, const vector<int> & taxIdListAtRank, const vector<int> & taxIdList)
{
    struct MmapedData<char> seqFile = mmapData<char>(seqFileName);
	cout<<"after mmap"<<endl;
    vector<int> startsOfTaxIDs; ///start points where sequences of next species start

    startsOfTaxIDs.push_back(0);
    vector<int> seqCntOfTaxIDs; ///the number of sequences for each species
    vector<FastaSplit> splits;

    getFastaSplits(taxIdListAtRank, splits);

    cout<<"here"<<endl;

    int seqCnt = 0;
    int currentTaxID = taxIdListAtRank[0];
    for(size_t i = 0; i < taxIdListAtRank.size(); i++){
        if(taxIdListAtRank[i] == currentTaxID){
            seqCnt++;
        }
        else{
            seqCntOfTaxIDs.push_back(seqCnt);
            startsOfTaxIDs.push_back(i);
            currentTaxID = taxIdListAtRank[i];
            seqCnt = 1;
        }
    }
    seqCntOfTaxIDs.push_back(seqCnt);
    vector<Sequence> sequences;
    SeqIterator::getSeqSegmentsWithHead(sequences, seqFile);
    size_t numOfSplits = splits.size();
    bool splitChecker[numOfSplits];
    fill_n(splitChecker, numOfSplits, false);

    size_t numOfTaxIDs = startsOfTaxIDs.size();
    bool processedSeqChecker[numOfTaxIDs];
    fill_n(processedSeqChecker, numOfTaxIDs, false);

    TargetKmerBuffer kmerBuffer(kmerBufSize);
    size_t processedTaxIDCnt = 0;
    size_t processedSplitCnt = 0;

    while(processedSplitCnt < numOfSplits){ ///check this condition
      //  fillTargetKmerBuffer(kmerBuffer, seqFile, sequences, processedSeqChecker, processedTaxIDCnt, startsOfTaxIDs, seqCntOfTaxIDs);
       	fillTargetKmerBuffer2(kmerBuffer, seqFile, sequences, splitChecker, processedSplitCnt, splits);
        cout<<"before writing"<<endl;
       	writeTargetFiles(kmerBuffer.buffer, kmerBuffer.startIndexOfReserve, outputFileName, taxIdListAtRank, taxIdList);
        cout<<"after writing "<<kmerBuffer.startIndexOfReserve<<endl;
    }

    free(kmerBuffer.buffer);
    munmap(seqFile.data, seqFile.fileSize + 1);
}

size_t IndexCreator::fillTargetKmerBuffer(TargetKmerBuffer & kmerBuffer, MmapedData<char> & seqFile, vector<Sequence> & seqs, bool * checker, size_t & processedTaxIdCnt, const vector<int> & startsOfTaxIDs, const vector<int> & seqCntOfTaxIDs) {
#ifdef OPENMP
   //omp_set_num_threads(1);
#endif
//    for(int i = 0; i < startsOfTaxIDs.size();i++){
//        cout<<startsOfTaxIDs[i]<<endl;
//        cout<<seqCntOfTaxIDs[i]<<endl<<endl;
//    }
    bool hasOverflow = false;
#pragma omp parallel default(none), shared(checker, hasOverflow, seqCntOfTaxIDs, seqFile, seqs, kmerBuffer, processedTaxIdCnt, startsOfTaxIDs, cout)
    {
        ProdigalWrapper prodigal;
        SeqIterator seqIterator;
        size_t posToWrite;

        PredictedBlock * blocks;
        size_t numOfBlocks;
        size_t totalSeqLengthForOneTaxID;
        size_t totalKmerCntForOneTaxID = 0;

#pragma omp for schedule(dynamic, 1)
        for (size_t i = 0; i < startsOfTaxIDs.size() ; i++) {
            if((checker[i] == false) && (!hasOverflow)) {
               // size_t numOfBlocksList[seqCntOfTaxIDs[i]];
	
                size_t * numOfBlocksList = (size_t*)malloc(seqCntOfTaxIDs[i] * sizeof(size_t));

                kseq_buffer_t buffer(const_cast<char *>(&seqFile.data[seqs[startsOfTaxIDs[i]].start]), seqs[startsOfTaxIDs[i]].length);
                kseq_t *seq = kseq_init(&buffer);
                kseq_read(seq);

                prodigal.is_meta = 0;
                if(strlen(seq->seq.s) < 20000){
                    prodigal.is_meta = 1;
                    prodigal.trainMeta(seq->seq.s);
                }else{
                    prodigal.trainASpecies(seq->seq.s);
                }

                prodigal.getPredictedFrames(seq->seq.s);
                ///fix here

                blocks = (PredictedBlock*)malloc(((prodigal.getNumberOfPredictedGenes() + 1) * (seqCntOfTaxIDs[i] + 1)) * 5 * sizeof(PredictedBlock));
                numOfBlocks = 0;

                seqIterator.getTranslationBlocks(prodigal.genes, prodigal.nodes, blocks, prodigal.getNumberOfPredictedGenes(), strlen(seq->seq.s), numOfBlocks);
                numOfBlocksList[0] = numOfBlocks;
                for(size_t p = 1; p < seqCntOfTaxIDs[i]; p++ ) {
                    kseq_buffer_t buffer2(const_cast<char *>(&seqFile.data[seqs[startsOfTaxIDs[i] + p].start]),
                                         seqs[startsOfTaxIDs[i] + p].length);
                    kseq_t *seq2 = kseq_init(&buffer2);
                    kseq_read(seq2);
                    prodigal.getPredictedFrames(seq2->seq.s);

                    seqIterator.getTranslationBlocks(prodigal.genes, prodigal.nodes, blocks,
                                                     prodigal.getNumberOfPredictedGenes(), strlen(seq->seq.s),
                                                     numOfBlocks);
                    cout<<"upper loop: "<<startsOfTaxIDs[i] + p<<endl;

                    numOfBlocksList[p] = numOfBlocks;

                }

                totalKmerCntForOneTaxID = 0;
                for(size_t block = 0; block < numOfBlocks; block++){
                    totalKmerCntForOneTaxID += seqIterator.getNumOfKmerForBlock(blocks[block]); ///이 함수 정확한지 확인해야
                }
                posToWrite = kmerBuffer.reserveMemory(totalKmerCntForOneTaxID);

                if(posToWrite + totalKmerCntForOneTaxID < kmerBuffer.bufferSize){
                    size_t start = 0;
                    for(size_t seqCnt = 0; seqCnt < seqCntOfTaxIDs[i]; seqCnt++){
                        kseq_buffer_t buffer(const_cast<char *>(&seqFile.data[seqs[startsOfTaxIDs[i]+seqCnt].start]), seqs[startsOfTaxIDs[i]+seqCnt].length);
                        kseq_t *seq = kseq_init(&buffer);
                        kseq_read(seq);
                        size_t end = numOfBlocksList[seqCnt];
                        cout<<"lower loop: "<<startsOfTaxIDs[i] + seqCnt<<endl;
                        for(size_t bl = start; bl < end ; bl++){
                            seqIterator.translateBlock(seq->seq.s,blocks[bl]);
                            seqIterator.fillBufferWithKmerFromBlock(blocks[bl], seq->seq.s, kmerBuffer, posToWrite, startsOfTaxIDs[i]+seqCnt);
                        }
                        cout<<"after filling: "<<startsOfTaxIDs[i] + seqCnt<<" "<<kmerBuffer.startIndexOfReserve<<endl;
                        start = numOfBlocksList[seqCnt];
                    }
                    checker[i] = true;
                    processedTaxIdCnt ++;
                }else {
                    kmerBuffer.startIndexOfReserve -= totalKmerCntForOneTaxID;
                    hasOverflow = true;
                }
                free(blocks);
                free(numOfBlocksList);
            }
        }
    }
    cout<<"before return: "<<kmerBuffer.startIndexOfReserve<<endl;
    return 0;
}

size_t IndexCreator::fillTargetKmerBuffer2(TargetKmerBuffer & kmerBuffer, MmapedData<char> & seqFile, vector<Sequence> & seqs, bool * checker, size_t & processedSplitCnt, const vector<FastaSplit> & splits) {
#ifdef OPENMP
  // omp_set_num_threads(1);
#endif
//    for(int i = 0; i < startsOfTaxIDs.size();i++){
//        cout<<startsOfTaxIDs[i]<<endl;
//        cout<<seqCntOfTaxIDs[i]<<endl<<endl;
//    }
    bool hasOverflow = false;
    size_t numOfGene =0;
#pragma omp parallel default(none), shared(checker, hasOverflow, splits, seqFile, seqs, kmerBuffer, processedSplitCnt, cout, numOfGene)
    {
        ProdigalWrapper prodigal;
        SeqIterator seqIterator;
        size_t posToWrite;

        PredictedBlock * blocks = (PredictedBlock*)malloc(10);
        size_t * numOfBlocksList =(size_t*)malloc(10);
        size_t numOfBlocks;

        size_t totalKmerCntForOneTaxID = 0;

#pragma omp for schedule(dynamic, 1)
        for (size_t i = 0; i < splits.size() ; i++) {
            if((checker[i] == false) && (!hasOverflow)) {
//                size_t * numOfBlocksList = (size_t*)malloc(splits[i].cnt * sizeof(size_t));
                realloc(numOfBlocksList, (splits[i].cnt + 1) * sizeof(size_t));
                kseq_buffer_t buffer(const_cast<char *>(&seqFile.data[seqs[splits[i].start].start]), seqs[splits[i].start].length);
                kseq_t *seq = kseq_init(&buffer);
                kseq_read(seq);
                prodigal.is_meta = 0;
                if(strlen(seq->seq.s) < 20000){
                    prodigal.is_meta = 1;
                    prodigal.trainMeta(seq->seq.s);
                }else{
                    prodigal.trainASpecies(seq->seq.s);
                }

                prodigal.getPredictedFrames(seq->seq.s);
                ///fix here
                realloc(blocks,(prodigal.getNumberOfPredictedGenes() + 1) * (splits[i].cnt + 1) * 5 * sizeof(PredictedBlock));
               // blocks = (PredictedBlock*)realloc(((prodigal.getNumberOfPredictedGenes() + 1) * (splits[i].cnt + 1)) * 5 * sizeof(PredictedBlock));
                numOfBlocks = 0;

                for(size_t p = 0; p < splits[i].cnt; p++ ) {
                    buffer = {const_cast<char *>(&seqFile.data[seqs[splits[i].start + splits[i].offset + p].start]), seqs[splits[i].start + splits[i].offset + p].length};
                    seq = kseq_init(&buffer);
                    kseq_read(seq);

                    prodigal.getPredictedFrames(seq->seq.s);
                    numOfGene += prodigal.getNumberOfPredictedGenes();
                    seqIterator.getTranslationBlocks(prodigal.genes, prodigal.nodes, blocks,
                                                     prodigal.getNumberOfPredictedGenes(), strlen(seq->seq.s),
                                                     numOfBlocks);
                    cout<<"upper loop: "<<splits[i].start + splits[i].offset + p<<endl;
                    numOfBlocksList[p] = numOfBlocks;
                }

                totalKmerCntForOneTaxID = 0;
                for(size_t block = 0; block < numOfBlocks; block++){
                    totalKmerCntForOneTaxID += seqIterator.getNumOfKmerForBlock(blocks[block]); ///이 함수 정확한지 확인해야
                }
                posToWrite = kmerBuffer.reserveMemory(totalKmerCntForOneTaxID);

                if(posToWrite + totalKmerCntForOneTaxID < kmerBuffer.bufferSize){
                    size_t start = 0;
                    for(size_t seqIdx = 0; seqIdx < splits[i].cnt; seqIdx++){
                        buffer = {const_cast<char *>(&seqFile.data[seqs[splits[i].start + splits[i].offset + seqIdx].start]), seqs[splits[i].start + splits[i].offset + seqIdx].length};
                        seq = kseq_init(&buffer);
                        kseq_read(seq);

                        size_t end = numOfBlocksList[seqIdx];
                        cout<<"before writing "<<kmerBuffer.startIndexOfReserve<<" "<<posToWrite<<endl;
                        for(size_t bl = start; bl < end ; bl++){
                            seqIterator.translateBlock(seq->seq.s,blocks[bl]);
                            seqIterator.fillBufferWithKmerFromBlock(blocks[bl], seq->seq.s, kmerBuffer, posToWrite, splits[i].start + splits[i].offset + seqIdx);
                        }
                        cout<<"after writing "<<kmerBuffer.startIndexOfReserve<<" "<<posToWrite<<endl;
                        start = numOfBlocksList[seqIdx];
                    }
                    checker[i] = true;
                    processedSplitCnt ++;
                }else {
                    kmerBuffer.startIndexOfReserve -= totalKmerCntForOneTaxID;
                    hasOverflow = true;
                }
            }
        }
        free(numOfBlocksList);
        free(blocks);
    }
    cout<<"before return: "<<kmerBuffer.startIndexOfReserve<<" "<<double(numOfGene)/double(4136)<<endl;
    return 0;
}
///This function sort the TargetKmerBuffer, do redundancy reducing task, write the differential index of them
void IndexCreator::writeTargetFiles(TargetKmer * kmerBuffer, size_t & kmerNum, const char * outputFileName, const vector<int> & taxIdListAtRank, const vector<int> & taxIdList)
{
    ///open a split file, which will be merged later.
    char suffixedDiffIdxFileName[100];
    char suffixedInfoFileName[100];
    sprintf(suffixedDiffIdxFileName,"%s_%zu_diffIdx", outputFileName,numOfFlush);
    sprintf(suffixedInfoFileName,"%s_%zu_info", outputFileName,numOfFlush);
    FILE * diffIdxFile = fopen(suffixedDiffIdxFileName, "wb");
    FILE * infoFile = fopen(suffixedInfoFileName, "wb");
    if (diffIdxFile == NULL || infoFile == NULL){
        cout<<"Cannot open the file for writing target DB"<<endl;
        return;
    }
    numOfFlush++;

    uint16_t *kmerLocalBuf = (uint16_t *)malloc(sizeof(uint16_t) * kmerBufSize);
    size_t localBufIdx = 0;
    uint64_t lastKmer = 0;

    ///Sort for differential indexing
    SORT_PARALLEL(kmerBuffer, kmerBuffer + kmerNum, IndexCreator::compareForDiffIdx);

    ///Redundancy reduction task
    TargetKmer lookingKmer = kmerBuffer[0];
    size_t write = 0;
    int endFlag = 0;
    int hasSeenOtherStrains;
    int redundant;
    size_t asd= 0;

    for(size_t i = 1 ; i < kmerNum ; i++) {
        hasSeenOtherStrains = 0;
        redundant = 0;
        while(taxIdListAtRank[lookingKmer.info.sequenceID] == taxIdListAtRank[kmerBuffer[i].info.sequenceID]){
            if (lookingKmer.ADkmer != kmerBuffer[i].ADkmer) {
                break;
            }
            if (taxIdList[lookingKmer.info.sequenceID] != taxIdList[kmerBuffer[i].info.sequenceID]){
                hasSeenOtherStrains = 1;
            }
            redundant = 1;
            i++;
            if(i == kmerNum){
                endFlag = 1;
                break;
            }
        }

        if(hasSeenOtherStrains == 1){
            lookingKmer.info.redundancy = true;
        }

        if((redundant == 1) && (hasSeenOtherStrains == 1)){
            asd ++;
        }


        fwrite(&lookingKmer.info, sizeof(TargetKmerInfo), 1, infoFile);
        write++;
        getDiffIdx(lastKmer, lookingKmer.ADkmer, diffIdxFile, kmerLocalBuf, localBufIdx);

        if(endFlag == 1) break;
        lastKmer = lookingKmer.ADkmer;
        lookingKmer = kmerBuffer[i];
    }

    //For the end part
    if(!((kmerBuffer[kmerNum - 2].ADkmer == kmerBuffer[kmerNum - 1].ADkmer) &&
        (taxIdListAtRank[kmerBuffer[kmerNum - 2].info.sequenceID] == taxIdListAtRank[kmerBuffer[kmerNum - 1].info.sequenceID]))){
        fwrite(&lookingKmer.info, sizeof(TargetKmerInfo), 1, infoFile);
        write++;
        getDiffIdx(lastKmer, lookingKmer.ADkmer, diffIdxFile, kmerLocalBuf, localBufIdx);
    }
    cout<<asd<<endl;
    cout<<"total k-mer count  : "<< kmerNum << endl;
    cout<<"written k-mer count: "<<write<<endl;

    flushKmerBuf(kmerLocalBuf, diffIdxFile, localBufIdx);
    free(kmerLocalBuf);
    fclose(diffIdxFile);
    fclose(infoFile);
    kmerNum = 0;
}


void IndexCreator::getDiffIdx(const uint64_t & lastKmer, const uint64_t & entryToWrite, FILE* handleKmerTable, uint16_t *kmerBuf, size_t & localBufIdx ){
    uint64_t kmerdiff = entryToWrite - lastKmer;
    uint16_t buffer[5];
    int idx = 3;
    buffer[4] = SET_END_FLAG(GET_15_BITS(kmerdiff));
    kmerdiff >>= 15U;
    while (kmerdiff) {
        uint16_t toWrite = GET_15_BITS(kmerdiff);
        kmerdiff >>= 15U;
        buffer[idx] = toWrite;
        idx--;
    }
    writeDiffIdx(kmerBuf, handleKmerTable, (buffer + idx + 1), (4 - idx), localBufIdx);
}

void IndexCreator::flushKmerBuf(uint16_t *buffer, FILE *handleKmerTable, size_t & localBufIdx ) {
    fwrite(buffer, sizeof(uint16_t), localBufIdx, handleKmerTable);
    localBufIdx = 0;
}

void IndexCreator::writeDiffIdx(uint16_t *buffer, FILE* handleKmerTable, uint16_t *toWrite, size_t size, size_t & localBufIdx ) {
    if (localBufIdx + size >= kmerBufSize) {
        flushKmerBuf(buffer, handleKmerTable, localBufIdx);
    }
    memcpy(buffer + localBufIdx, toWrite, sizeof(uint16_t) * size);
    localBufIdx += size;
}

void IndexCreator::writeInfo(TargetKmerInfo * entryToWrite, FILE * infoFile, TargetKmerInfo * infoBuffer, size_t & infoBufferIdx)
{
    if (infoBufferIdx >= kmerBufSize) {
        flushInfoBuf(infoBuffer, infoFile, infoBufferIdx);
    }
    memcpy(infoBuffer + infoBufferIdx, entryToWrite, sizeof(TargetKmerInfo));
    infoBufferIdx++;
}
void IndexCreator::flushInfoBuf(TargetKmerInfo * buffer, FILE * infoFile, size_t & localBufIdx ) {
    fwrite(buffer, sizeof(TargetKmerInfo), localBufIdx, infoFile);
    localBufIdx = 0;
}
int IndexCreator::getNumOfFlush()
{
    return numOfFlush;
}

//void IndexCreator::printOutTargetDB(char * diff, char * info){
//    struct MmapedData<uint16_t> diffFile = mmapData<uint16_t>(diff);
//    struct MmapedData<TargetKmerInfo> infoFile = mmapData<TargetKmerInfo>(info);
//    Classifier classifier;
//    size_t numOfkmer = infoFile.fileSize/sizeof(TargetKmerInfo);
//    size_t diffIdxPos = 0;
//    uint64_t nextTargetKmer = classifier.getNextTargetKmer(0, diffFile.data, diffIdxPos);
//    for( size_t i = 0; i < numOfkmer; i++){
//        cout<<infoFile.data[i].sequenceID<<endl;
//    }
//}

bool IndexCreator::compareForDiffIdx(const TargetKmer & a, const TargetKmer & b){
    return a.ADkmer < b.ADkmer || (a.ADkmer == b.ADkmer && a.info.sequenceID < b.info.sequenceID);
}

void IndexCreator::getFastaSplits(const vector<int> & taxIdListAtRank, vector<FastaSplit> & fastaSplit){
    size_t start = 0;
    size_t idx = 0;
    uint32_t offset = 0;
    uint32_t cnt = 0;
    int currentTaxId = taxIdListAtRank[0];
    while(idx < taxIdListAtRank.size()){
        offset = 0;
        cnt = 0;
        currentTaxId = taxIdListAtRank[start];
        while(currentTaxId == taxIdListAtRank[idx] && idx < taxIdListAtRank.size()){
            cnt ++;
            idx ++;
            if(cnt > 50){
                fastaSplit.emplace_back(start, offset, cnt-1);
                offset += cnt - 1;
                cnt = 1;
            }
        }
        fastaSplit.emplace_back(start, offset, cnt);
        start = idx;
    }
}