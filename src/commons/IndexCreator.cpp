#include "IndexCreator.h"

IndexCreator::IndexCreator()
{

}

IndexCreator::~IndexCreator() {
}

///It reads a reference sequence file and write a differential index file and target k-mer info file.
void IndexCreator::startIndexCreatingParallel(const char * seqFileName, const char * outputFileName,
                                              const vector<int> & taxIdListAtRank, const vector<int> & taxIdList,
                                              const LocalParameters & par)
{
    ///Mmap the input fasta file
    struct MmapedData<char> seqFile = mmapData<char>(seqFileName);

    ///Getting start and end position of each sequence
    vector<Sequence> sequences;
    getSeqSegmentsWithHead(sequences, seqFile);

    // Sequences in the same split share the sequence to be used for training the prodigal.
    vector<FastaSplit> splits;
    getFastaSplits(taxIdListAtRank, splits, sequences);
    size_t numOfSplits = splits.size();

    bool * splitChecker = new bool[numOfSplits];
    fill_n(splitChecker, numOfSplits, false);
    TargetKmerBuffer kmerBuffer(10'000'000'000);
    size_t processedSplitCnt = 0;
    while(processedSplitCnt < numOfSplits){ ///check this condition
        fillTargetKmerBuffer2(kmerBuffer, seqFile, sequences, splitChecker,processedSplitCnt, splits, taxIdListAtRank, par);
        writeTargetFiles(kmerBuffer.buffer, kmerBuffer.startIndexOfReserve, outputFileName, taxIdList);
    }

    //free(kmerBuffer.buffer);
    delete[] splitChecker;
    munmap(seqFile.data, seqFile.fileSize + 1);
}

size_t IndexCreator::fillTargetKmerBuffer(TargetKmerBuffer & kmerBuffer,
                                          MmapedData<char> & seqFile,
                                          vector<Sequence> & seqs,
                                          bool * checker, size_t & processedSplitCnt,
                                          const vector<FastaSplit> & splits,
                                          const vector<int> & taxIdListAtRank,
                                          const LocalParameters & par) {
#ifdef OPENMP
   omp_set_num_threads(*(int * )par.PARAM_THREADS.value);
#endif
    bool hasOverflow = false;

    #pragma omp parallel default(none), shared(par, checker, hasOverflow, splits, seqFile, seqs, kmerBuffer, processedSplitCnt, cout, taxIdListAtRank)
    {
        ProdigalWrapper prodigal;
        SeqIterator seqIterator(par);
        size_t posToWrite;
        size_t numOfBlocks;
        size_t totalKmerCntForOneTaxID;

#pragma omp for schedule(dynamic, 1)
        for (size_t i = 0; i < splits.size() ; i++) {
            if((checker[i] == false) && (!hasOverflow)) {
                size_t * numOfBlocksList = (size_t*)malloc(splits[i].cnt * sizeof(size_t));

                // Train Prodigal with a training sequence of i th split
                kseq_buffer_t buffer(const_cast<char *>(&seqFile.data[seqs[splits[i].training].start]), seqs[splits[i].training].length);
                kseq_t *seq = kseq_init(&buffer);
                kseq_read(seq);
                size_t lengthOfTrainingSeq = strlen(seq->seq.s);
                prodigal.is_meta = 0;
                string tmp = seq->name.s;
                if(strlen(seq->seq.s) < 20000){
                    prodigal.is_meta = 1;
                    cout<<"train Meta: "<<splits[i].training<<" "<<seqs[splits[i].training].start<<" "<<i<<seq->headerOffset<<" "<<splits[i].offset<<" "<<splits[i].cnt<<endl;
                    cout<<tmp<<endl;
                    prodigal.trainMeta(seq->seq.s);
                    cout<<"after train meta"<<endl;
                }else{
                    prodigal.trainASpecies(seq->seq.s);
                }

                ///Get min k-mer hash list for determining strandness
                priority_queue<uint64_t> standardList;
                seqIterator.getMinHashList(standardList, seq->seq.s);

                //prodigal.getPredictedGenes(seq->seq.s);
                vector<PredictedBlock> blocks;

                numOfBlocks = 0;

                ///Getting all the sequence blocks of current split. Each block will be translated later separately.
                priority_queue<uint64_t> currentList;
                for(size_t p = 0; p < splits[i].cnt; p++ ) {
                    buffer = {const_cast<char *>(&seqFile.data[seqs[splits[i].offset + p].start]), static_cast<size_t>(seqs[splits[i].offset + p].length)};
                    seq = kseq_init(&buffer);
                    kseq_read(seq);
                    seqIterator.getMinHashList(currentList, seq->seq.s);
                    prodigal.getPredictedGenes(seq->seq.s);
                    prodigal.removeCompletelyOverlappingGenes();
                    prodigal.printGenes();
                    if(seqIterator.compareMinHashList(standardList, currentList, lengthOfTrainingSeq, strlen(seq->seq.s))){
                        seqIterator.getTranslationBlocks(prodigal.finalGenes, prodigal.nodes, blocks,
                                                         prodigal.fng, strlen(seq->seq.s),
                                                         numOfBlocks);
                    } else {
                        SeqIterator::getTranslationBlocksReverse(prodigal.finalGenes, prodigal.nodes, blocks,
                                                                prodigal.fng, strlen(seq->seq.s),
                                                         numOfBlocks);
                    }
                    numOfBlocksList[p] = numOfBlocks;
                    currentList = priority_queue<uint64_t>();
                }

                /// Calculate the number of k-mers to reserve memory of k-mer buffer
                totalKmerCntForOneTaxID = 0;
                for(size_t block = 0; block < numOfBlocks; block++){
                    totalKmerCntForOneTaxID += seqIterator.getNumOfKmerForBlock(blocks[block]);
                }

                /// Fill k-mer buffer with k-mers of current split if the buffer has enough space
                posToWrite = kmerBuffer.reserveMemory(totalKmerCntForOneTaxID);
                if(posToWrite + totalKmerCntForOneTaxID < kmerBuffer.bufferSize){
                    size_t start = 0;
                    for(size_t seqIdx = 0; seqIdx < splits[i].cnt; seqIdx++){
                        buffer = {const_cast<char *>(&seqFile.data[seqs[splits[i].offset + seqIdx].start]), static_cast<size_t>(seqs[splits[i].offset + seqIdx].length)};
                        seq = kseq_init(&buffer);
                        kseq_read(seq);
                        size_t end = numOfBlocksList[seqIdx];
                        for(size_t bl = start; bl < end ; bl++){
//                            if(!(seqIterator.translateBlock(seq->seq.s,blocks[bl]))){
//                                cout<<seq->name.s<<endl;
//                            }
                            seqIterator.translateBlock(seq->seq.s,blocks[bl]);
                            seqIterator.fillBufferWithKmerFromBlock(blocks[bl], seq->seq.s, kmerBuffer, posToWrite, splits[i].offset + seqIdx, taxIdListAtRank[splits[i].offset + seqIdx]); //splits[i].offset + seqIdx
                        }
                        start = numOfBlocksList[seqIdx];
                    }
                    checker[i] = true;
                    #pragma omp atomic
                    processedSplitCnt ++;
                }else {
                    ///Withdraw the reservation if the buffer is full.
                    #pragma omp atomic
                    kmerBuffer.startIndexOfReserve -= totalKmerCntForOneTaxID;
                    cout<<"buffer is full"<<endl;
                    hasOverflow = true;
                }
                kseq_destroy(seq);
                free(numOfBlocksList);
                blocks.clear();
            }
        }
    }
    cout<<"before return: "<<kmerBuffer.startIndexOfReserve<<endl;
    return 0;
}


size_t IndexCreator::fillTargetKmerBuffer2(TargetKmerBuffer & kmerBuffer,
                                           MmapedData<char> & seqFile,
                                           vector<Sequence> & seqs,
                                           bool * checker,
                                           size_t & processedSplitCnt,
                                           const vector<FastaSplit> & splits,
                                           const vector<int> & taxIdListAtRank,
                                           const LocalParameters & par) {
#ifdef OPENMP
    omp_set_num_threads(*(int * )par.PARAM_THREADS.value);
#endif
    bool hasOverflow = false;

#pragma omp parallel default(none), shared(par, checker, hasOverflow, splits, seqFile, seqs, kmerBuffer, processedSplitCnt, cout, taxIdListAtRank)
    {
        ProdigalWrapper prodigal;
        SeqIterator seqIterator(par);
        size_t posToWrite;
        size_t numOfBlocks;
        size_t totalKmerCntForOneTaxID;
        vector<uint64_t> intergenicKmerList;
        vector<PredictedBlock> blocks;
        priority_queue<uint64_t> standardList;
        priority_queue<uint64_t> currentList;
        size_t lengthOfTrainingSeq;
        char * reverseCompliment;
        vector<bool> strandness;

#pragma omp for schedule(dynamic, 1)
        for (size_t i = 0; i < splits.size() ; i++) {
            if((checker[i] == false) && (!hasOverflow)) {
                size_t * numOfBlocksList = (size_t*)malloc(splits[i].cnt * sizeof(size_t));
                intergenicKmerList.clear();
                strandness.clear();
                standardList = priority_queue<uint64_t>();

                //Train Prodigal with a training sequence of i th split
                kseq_buffer_t buffer(const_cast<char *>(&seqFile.data[seqs[splits[i].training].start]), seqs[splits[i].training].length);
                kseq_t *seq = kseq_init(&buffer);
                kseq_read(seq);
                lengthOfTrainingSeq = strlen(seq->seq.s);
                prodigal.is_meta = 0;

                if(strlen(seq->seq.s) < 100000){
                    prodigal.is_meta = 1;
                    prodigal.trainMeta(seq->seq.s);
                }else{
                    prodigal.trainASpecies(seq->seq.s);
                }

                // Generate intergenic 23-mer list
                prodigal.getPredictedGenes(seq->seq.s);
                seqIterator.generateIntergenicKmerList(prodigal.genes, prodigal.nodes, prodigal.getNumberOfPredictedGenes(), intergenicKmerList, seq->seq.s);

                // Get min k-mer hash list for determining strandness
                seqIterator.getMinHashList(standardList, seq->seq.s);

                // Getting all the sequence blocks of current split. Each block will be translated later separately.
                numOfBlocks = 0;
                for(size_t p = 0; p < splits[i].cnt; p++ ) {
                    buffer = {const_cast<char *>(&seqFile.data[seqs[splits[i].offset + p].start]), static_cast<size_t>(seqs[splits[i].offset + p].length)};
                    kseq_destroy(seq);
                    seq = kseq_init(&buffer);
                    kseq_read(seq);
                    seqIterator.getMinHashList(currentList, seq->seq.s);
                    if(seqIterator.compareMinHashList(standardList, currentList, lengthOfTrainingSeq, strlen(seq->seq.s))){
                        prodigal.getPredictedGenes(seq->seq.s);
                        prodigal.removeCompletelyOverlappingGenes();
                        seqIterator.getTranslationBlocks2(prodigal.finalGenes, prodigal.nodes, blocks,
                                                          prodigal.fng, strlen(seq->seq.s),
                                                          numOfBlocks, intergenicKmerList, seq->seq.s);
                        strandness.push_back(true);
                    } else{
                        reverseCompliment = seqIterator.reverseCompliment(seq->seq.s, strlen(seq->seq.s));
                        prodigal.getPredictedGenes(reverseCompliment);
                        prodigal.removeCompletelyOverlappingGenes();
                        seqIterator.getTranslationBlocks2(prodigal.finalGenes, prodigal.nodes, blocks,
                                                          prodigal.fng, strlen(reverseCompliment),
                                                          numOfBlocks, intergenicKmerList, reverseCompliment);
                        free(reverseCompliment);
                        strandness.push_back(false);
                    }
                    numOfBlocksList[p] = numOfBlocks;
                    currentList = priority_queue<uint64_t>();
                    if(hasOverflow){
                        break;
                    }
                }
                if(hasOverflow){
                    kseq_destroy(seq);
                    free(numOfBlocksList);
                    blocks.clear();
                    break;
                }
                // Calculate the number of k-mers to reserve memory of k-mer buffer
                totalKmerCntForOneTaxID = 0;
                for(size_t block = 0; block < numOfBlocks; block++){
                    totalKmerCntForOneTaxID += seqIterator.getNumOfKmerForBlock(blocks[block]);
                }

                // Fill k-mer buffer with k-mers of current split if the buffer has enough space
                posToWrite = kmerBuffer.reserveMemory(totalKmerCntForOneTaxID);
                if(posToWrite + totalKmerCntForOneTaxID < kmerBuffer.bufferSize){
                    size_t start = 0;
                    for(size_t seqIdx = 0; seqIdx < splits[i].cnt; seqIdx++){
                        buffer = {const_cast<char *>(&seqFile.data[seqs[splits[i].offset + seqIdx].start]), static_cast<size_t>(seqs[splits[i].offset + seqIdx].length)};
                        kseq_destroy(seq);
                        seq = kseq_init(&buffer);
                        kseq_read(seq);
                        size_t end = numOfBlocksList[seqIdx];
                        if(strandness[seqIdx]){
                            for(size_t bl = start; bl < end ; bl++){
                                seqIterator.translateBlock(seq->seq.s,blocks[bl]);
                                seqIterator.fillBufferWithKmerFromBlock(blocks[bl], seq->seq.s, kmerBuffer, posToWrite, splits[i].offset + seqIdx, taxIdListAtRank[splits[i].offset + seqIdx]); //splits[i].offset + seqIdx
                            }
                        } else{
                            reverseCompliment = seqIterator.reverseCompliment(seq->seq.s, strlen(seq->seq.s));
                            for(size_t bl = start; bl < end ; bl++){
                                seqIterator.translateBlock(reverseCompliment,blocks[bl]);
                                seqIterator.fillBufferWithKmerFromBlock(blocks[bl], reverseCompliment, kmerBuffer, posToWrite, splits[i].offset + seqIdx, taxIdListAtRank[splits[i].offset + seqIdx]); //splits[i].offset + seqIdx
                            }
                            free(reverseCompliment);
                        }
                        start = numOfBlocksList[seqIdx];
                    }
                    checker[i] = true;
                    __sync_fetch_and_add(&processedSplitCnt, 1);
                }else {
                    // Withdraw the reservation if the buffer is full.
                   __sync_fetch_and_add(&kmerBuffer.startIndexOfReserve, -1);
                    cout<<"buffer is full"<<endl;
                    hasOverflow = true;
                }
                kseq_destroy(seq);
                free(numOfBlocksList);
                blocks.clear();
            }
        }
    }
    cout<<"before return: "<<kmerBuffer.startIndexOfReserve<<endl;
    return 0;
}

///This function sort the TargetKmerBuffer, do redundancy reducing task, write the differential index of them
void IndexCreator::writeTargetFiles(TargetKmer * kmerBuffer, size_t & kmerNum, const char * outputFileName, const vector<int> & taxIdList)
{
    // Write a split file, which will be merged later.
    char suffixedDiffIdxFileName[300];
    char suffixedInfoFileName[300];
    sprintf(suffixedDiffIdxFileName, "%s/%zu_diffIdx", outputFileName, numOfFlush);
    sprintf(suffixedInfoFileName, "%s/%zu_info", outputFileName, numOfFlush);

    FILE * diffIdxFile = fopen(suffixedDiffIdxFileName, "wb");
    FILE * infoFile = fopen(suffixedInfoFileName, "wb");
    if (diffIdxFile == nullptr || infoFile == nullptr){
        cout<<"Cannot open the file for writing target DB"<<endl;
        return;
    }
    numOfFlush++;

    uint16_t *diffIdxBuffer = (uint16_t *)malloc(sizeof(uint16_t) * kmerBufSize);
    size_t localBufIdx = 0;
    uint64_t lastKmer = 0;

    ///Sort for differential indexing
    SORT_PARALLEL(kmerBuffer, kmerBuffer + kmerNum, IndexCreator::compareForDiffIdx);

    ///Find the first index of garbage k-mer (UINT64_MAX)
    for(size_t checkN = kmerNum - 1; checkN >= 0; checkN--){
        if(kmerBuffer[checkN].ADkmer != UINT64_MAX){
            kmerNum = checkN + 1;
            break;
        }
    }

    ///Redundancy reduction task
    TargetKmer lookingKmer = kmerBuffer[0];
    size_t write = 0;
    int endFlag = 0;
    int hasSeenOtherStrains;

    for(size_t i = 1 ; i < kmerNum ; i++) {
        hasSeenOtherStrains = 0;
        while(lookingKmer.taxIdAtRank == kmerBuffer[i].taxIdAtRank){
            if (lookingKmer.ADkmer != kmerBuffer[i].ADkmer) {
                break;
            }
            hasSeenOtherStrains += (taxIdList[lookingKmer.info.sequenceID] != taxIdList[kmerBuffer[i].info.sequenceID]);
            i++;
            if(i == kmerNum){
                endFlag = 1;
                break;
            }
        }

        lookingKmer.info.redundancy = (hasSeenOtherStrains > 0);

        fwrite(&lookingKmer.info, sizeof(TargetKmerInfo), 1, infoFile);
        write++;
        getDiffIdx(lastKmer, lookingKmer.ADkmer, diffIdxFile, diffIdxBuffer, localBufIdx);

        if(endFlag == 1) break;
        lastKmer = lookingKmer.ADkmer;
        lookingKmer = kmerBuffer[i];
    }

    //For the end part
    if(!((kmerBuffer[kmerNum - 2].ADkmer == kmerBuffer[kmerNum - 1].ADkmer) &&
         (kmerBuffer[kmerNum - 2].taxIdAtRank == kmerBuffer[kmerNum - 1].taxIdAtRank))){
        fwrite(&lookingKmer.info, sizeof(TargetKmerInfo), 1, infoFile);
        write++;
        getDiffIdx(lastKmer, lookingKmer.ADkmer, diffIdxFile, diffIdxBuffer, localBufIdx);
    }
    cout<<"total k-mer count  : "<< kmerNum << endl;
    cout<<"written k-mer count: "<<write<<endl;

    flushKmerBuf(diffIdxBuffer, diffIdxFile, localBufIdx);
    free(diffIdxBuffer);
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

inline bool IndexCreator::compareForDiffIdx(const TargetKmer & a, const TargetKmer & b){
    return a.ADkmer < b.ADkmer || (a.ADkmer == b.ADkmer && a.taxIdAtRank < b.taxIdAtRank);
}

void IndexCreator::getSeqSegmentsWithoutHead(vector<Sequence> & seqSegments, MmapedData<char> seqFile) {
    size_t start = 0;
    size_t numOfChar = seqFile.fileSize / sizeof(char);
    for(size_t i = 0; i < numOfChar; i++){
        if(seqFile.data[i] == '>'){
            seqSegments.emplace_back(start, i-2, i - start - 1);// the first push_back is a garbage.
            while(seqFile.data[i] != '\n'){
                i++;
            }
            start = i + 1;
        }
    }
    seqSegments.emplace_back(start, numOfChar - 2, numOfChar - start - 1);
}

void IndexCreator::getSeqSegmentsWithHead(vector<Sequence> & seqSegments, MmapedData<char> seqFile) {
    size_t start = 0;
    size_t numOfChar = seqFile.fileSize / sizeof(char);
    for(size_t i = 1; i < numOfChar; i++){
        if(seqFile.data[i] == '>'){
            seqSegments.emplace_back(start, i-2, i - start - 1);
            start = i;
        }
    }
    seqSegments.emplace_back(start, numOfChar - 2, numOfChar - start - 1);
}

void IndexCreator::getFastaSplits(const vector<int> & taxIdListAtRank, vector<FastaSplit> & fastaSplit, vector<Sequence> & seqs){
    size_t training = 0;
    size_t idx = 0;
    uint32_t offset = 0;
    uint32_t cnt = 0;
    int theLargest;
    int isLeftover;

    int currentTaxId;

    while(idx < taxIdListAtRank.size()){
        offset = idx;
        training = idx;
        cnt = 0;
        isLeftover = 0;
        currentTaxId = taxIdListAtRank[idx];
        while(currentTaxId == taxIdListAtRank[idx] && idx < taxIdListAtRank.size()){
            cnt ++;
            idx ++;
            if(cnt > 100){ //The largest number of consecutive plasmid is the smallest. The smaller, the more training and the less time of single threading
                theLargest = 0;
                for(uint32_t i = 0; i < cnt - 1; i++){
                    if(seqs[offset + i].length > theLargest){
                        training = offset + i;
                        theLargest = seqs[offset + i].length;
                    }
                }
                fastaSplit.emplace_back(training, offset, cnt - 1);
                offset += cnt - 1;
                cnt = 1;
                isLeftover = 1;
            }
        }

        if(isLeftover == 1){
            fastaSplit.emplace_back(training, offset, cnt);
        }else {
            theLargest = 0;
            for (uint32_t i = 0; i < cnt; i++) {
                if (seqs[offset + i].length > theLargest) {
                    training = offset + i;
                    theLargest = seqs[offset + i].length;
                }
            }
            fastaSplit.emplace_back(training, offset, cnt);
        }
    }
}
