#include "IndexCreator.h"

#include <utility>

IndexCreator::IndexCreator(const LocalParameters & par)
{
    dbDir = par.filenames[0];
    fnaListFileName = par.filenames[1];
    taxonomyDir = par.filenames[0] + "/taxonomy";


    // Load taxonomy
    taxonomy = new NcbiTaxonomy(taxonomyDir + "/names.dmp",
                                taxonomyDir + "/nodes.dmp",
                                taxonomyDir + "/merged.dmp");


    // ==== To re-use prodigal training data if possible ==== //
    if (par.tinfoPath != ""){
        tinfo_path = par.tinfoPath;
        tinfo_list = tinfo_path + "/tinfo.list";
    } else {
        tinfo_path = par.filenames[1] + "/tinfo/";
        tinfo_list = tinfo_path + "tinfo.list";
    }
    // Load trained species list
    if (FILE * tinfoList = fopen(tinfo_list.c_str(), "r")){
        char buffer[512];
        TaxID taxid;
        while(feof(tinfoList) == 0)
        {
            fscanf(tinfoList,"%s",buffer);
            trainedSpecies.push_back(atol(buffer));
        }
        trainedSpecies.pop_back();
        fclose(tinfoList);
    }
    // ======================================================= //

    if (par.reducedAA == 1){
        MARKER = 0Xffffffff;
        MARKER = ~ MARKER;
    } else {
        MARKER = 16777215;
        MARKER = ~ MARKER;
    }
}

IndexCreator::IndexCreator(const LocalParameters &par, string dbDir, string fnaListFileName,
                           string taxonomyDir, string acc2taxidFile)
        : dbDir(std::move(dbDir)), fnaListFileName(move(fnaListFileName)),
          taxonomyDir(move(taxonomyDir)), acc2taxidFileName(std::move(acc2taxidFile))
{
    // Load taxonomy
    taxonomy = new NcbiTaxonomy(this->taxonomyDir + "/names.dmp",
                                this->taxonomyDir + "/nodes.dmp",
                                this->taxonomyDir + "/merged.dmp");

    if (par.reducedAA == 1){
        MARKER = 0Xffffffff;
        MARKER = ~ MARKER;
    } else {
        MARKER = 16777215;
        MARKER = ~ MARKER;
    }
}


IndexCreator::~IndexCreator() {
    delete taxonomy;
}

void IndexCreator::createIndex(const LocalParameters &par) {

    makeBlocksForParallelProcessing();

    size_t numOfSplits = fnaSplits.size();
    bool * splitChecker = new bool[numOfSplits];
    fill_n(splitChecker, numOfSplits, false);
    size_t processedSplitCnt = 0;

    TargetKmerBuffer kmerBuffer(kmerBufSize);
    while(processedSplitCnt < numOfSplits){ // Check this condition
        fillTargetKmerBuffer(kmerBuffer, splitChecker, processedSplitCnt, par);
        time_t start = time(nullptr);
        SORT_PARALLEL(kmerBuffer.buffer, kmerBuffer.buffer + kmerBuffer.startIndexOfReserve, IndexCreator::compareForDiffIdx);
        time_t sort = time(nullptr);
        cout << "Sort time: " << sort - start << endl;
        auto * uniqKmerIdx = new size_t[kmerBuffer.startIndexOfReserve + 1];
        size_t uniqKmerCnt = 0;
        reduceRedundancy(kmerBuffer, uniqKmerIdx, uniqKmerCnt, par);
        time_t reduction = time(nullptr);
        cout<<"Time spent for reducing redundancy: "<<(double) (reduction - sort) << endl;
        if(processedSplitCnt == numOfSplits && numOfFlush == 0){
            writeTargetFilesAndSplits(kmerBuffer.buffer, kmerBuffer.startIndexOfReserve, par, uniqKmerIdx, uniqKmerCnt);
        } else {
            writeTargetFiles(kmerBuffer.buffer, kmerBuffer.startIndexOfReserve, par,uniqKmerIdx, uniqKmerCnt);
        }
    }

    delete[] splitChecker;

}

void IndexCreator::makeBlocksForParallelProcessing(){
    unordered_map<string, TaxID> acc2taxid;
    load_accession2taxid(acc2taxidFileName, acc2taxid);

    // Make blocks of sequences that can be processed in parallel
    int fileNum = getNumberOfLines(fnaListFileName);
    ifstream fnaListFile;
    fnaListFile.open(fnaListFileName);
    string eachFile;
    string seqHeader;
    sequenceOfFastas.resize(fileNum);

    if (!fnaListFile.is_open()) {
        Debug(Debug::ERROR) << "Cannot open file for file list" << "\n";
        EXIT(EXIT_FAILURE);
    }
    for (int i = 0; i < fileNum; ++i) {
        // Get start and end position of each sequence in the file
        getline(fnaListFile, eachFile);
        cout << eachFile << endl;
        fnaList.push_back(eachFile);
        seqHeader = getSeqSegmentsWithHead(sequenceOfFastas[i], eachFile); //TODO : get the accession and taxid here
        seqHeader = seqHeader.substr(0, seqHeader.find(' '));
        TaxID taxid = acc2taxid[seqHeader];
        TaxID speciesTaxid = taxonomy->getTaxIdAtRank(taxid, "species");
        taxIdList.push_back(taxid);
        cout << "sp " << speciesTaxid << endl;
        // Split current file into blocks for parallel processing
        splitFasta(i, speciesTaxid);
    }
    for (int i = 0; i < fileNum; ++i) {
        for (int j = 0; j < sequenceOfFastas[i].size(); ++j) {
            cout << i << sequenceOfFastas[i][j].start << " " << sequenceOfFastas[i][j].end << " " << sequenceOfFastas[i][j].length << endl;
        }
    }

    for (int i = 0; i < fnaSplits.size(); ++i) {
        cout << fnaSplits[i].file_idx << " " << fnaSplits[i].offset << " " << fnaSplits[i].cnt << " " << fnaSplits[i].training << endl;
    }
    fnaListFile.close();
}

void IndexCreator::splitFasta(int fnaIdx, TaxID speciesTaxid) {
    uint32_t offset = 0;
    uint32_t cnt = 0;
    int maxLength = 0;
    size_t seqForTraining = 0;
    vector<FnaSplit> tempSplits;
    size_t seqIdx = 0;
    while(seqIdx < sequenceOfFastas[fnaIdx].size()){
        if(speciesTaxid == 0) {
            seqIdx++;
            continue;
        }
        if (sequenceOfFastas[fnaIdx][seqIdx].length > maxLength){
            maxLength = sequenceOfFastas[fnaIdx][seqIdx].length;
            seqForTraining = seqIdx;
        }
        cnt ++;
        if(cnt > 30){
            tempSplits.emplace_back(0, offset, cnt - 1, speciesTaxid, fnaIdx);
            offset += cnt - 1;
            cnt = 1;
        }
        seqIdx ++;
    }
    tempSplits.emplace_back(0, offset, cnt, speciesTaxid, fnaIdx);
    // Update the training sequence
    for(auto & x : tempSplits){
        x.training = seqForTraining;
        fnaSplits.push_back(x);
    }
}

void IndexCreator::load_accession2taxid(const string & mappingFileName, unordered_map<string, int> & acc2taxid) {
    cerr << "Load mapping from accession ID to taxonomy ID" << endl;
    string eachLine;
    string eachItem;
    if (FILE * mappingFile = fopen(mappingFileName.c_str(), "r")) {
        char buffer[512];
        int taxID;
        fscanf(mappingFile, "%*s\t%*s\t%*s\t%*s");
        while (fscanf(mappingFile, "%*s\t%s\t%d\t%*d", buffer, &taxID) == 2 ){
            acc2taxid[string(buffer)] = taxID;
        }
    } else {
        cerr << "Cannot open file for mapping from accession to tax ID" << endl;
    }
}

void IndexCreator::startIndexCreatingParallel(const char * seqFileName, const char * outputFileName,
                                              const vector<TaxID> & superkingdom,
                                              const vector<int> & speciesTaxIDs, const vector<int> & taxIdList,
                                              const LocalParameters & par){ //build_fasta

    // Getting start and end position of each sequence
    cerr<< "Get start and end position of each sequence" << endl;
    vector<Sequence> sequences;
    getSeqSegmentsWithHead(sequences, seqFileName);

    // Sequences in the same split share the sequence to be used for training the prodigal.
    cerr<< "Split the FASTA into blocks for prodigal" << endl;
    vector<FastaSplit> splits;
    splitAFastaFile(speciesTaxIDs, splits, sequences);
    size_t numOfSplits = splits.size();
    for (auto x : splits){
        cout << x.training << " " << x.offset << " " << x.cnt << " " << x.taxid << endl;
    }


    bool * splitChecker = new bool[numOfSplits];
    fill_n(splitChecker, numOfSplits, false);
    TargetKmerBuffer kmerBuffer(kmerBufSize);
    size_t processedSplitCnt = 0;

    // Mmap the input fasta file
    struct MmapedData<char> seqFile = mmapData<char>(seqFileName);

    // Training prodigal
    trainProdigal(splits, sequences, seqFile, par);

    // Load the trained prodigal model
    loadTrainingInfo();

    while(processedSplitCnt < numOfSplits) {
        fillTargetKmerBuffer(kmerBuffer, seqFile, sequences, splitChecker, processedSplitCnt, splits, speciesTaxIDs,
                             superkingdom, par);
        writeTargetFiles(kmerBuffer.buffer, kmerBuffer.startIndexOfReserve, outputFileName, taxIdList);
    }

    delete[] splitChecker;
    munmap(seqFile.data, seqFile.fileSize + 1);
}

// It reads a reference sequence file and write a differential index file and target k-mer info file.
void IndexCreator::startIndexCreatingParallel(const LocalParameters & par) //build_dir
{
    const string folder = par.filenames[0];
    const string dbDirectory = par.filenames[1];
    const string taxonomyDirectory = dbDirectory + "/taxonomy";

    // Taxonomy
    const string names = taxonomyDirectory + "/names.dmp";
    const string nodes = taxonomyDirectory + "/nodes.dmp";
    const string merged = taxonomyDirectory + "/merged.dmp";
    NcbiTaxonomy taxonomy(names, nodes, merged);

    // Unzip genomes and make a list of them
    unzipAndList(folder, folder + "/fasta_list_GTDB");

    // Load mapping from assembly accession to taxonomy ID
    unordered_map<string, int> assacc2taxid;
    load_assacc2taxid( taxonomyDirectory + "/assacc_to_taxid.tsv", assacc2taxid);

    // Make mapping from tax id to FASTA file
    vector<TaxId2Fasta> taxid2fasta;
    mappingFromTaxIDtoFasta(folder + "/fasta_list_GTDB", assacc2taxid, taxid2fasta, taxonomy);

    // Write a file of tax ids
    string taxIdList_fname = dbDirectory + "/taxID_list";
    ofstream taxIdList;
    taxIdList.open(taxIdList_fname);
    for(auto & cnt : taxid2fasta){
        taxIdList<<cnt.taxid<<'\n';
    }
    taxIdList.close();

    // Divide FASTA files
    // Sequences in the same split share the sequence to be used for training the prodigal.
    vector<FastaSplit> splits;
    groupFastaFiles(taxid2fasta, splits);
    size_t numOfSplits = splits.size();
    bool * splitChecker = new bool[numOfSplits];
    fill_n(splitChecker, numOfSplits, false);
    size_t processedSplitCnt = 0;

    TargetKmerBuffer kmerBuffer(kmerBufSize);
    while(processedSplitCnt < numOfSplits){ // Check this condition
        fillTargetKmerBuffer(kmerBuffer, splitChecker, processedSplitCnt, splits, taxid2fasta, par);
        time_t start = time(nullptr);
        SORT_PARALLEL(kmerBuffer.buffer, kmerBuffer.buffer + kmerBuffer.startIndexOfReserve, IndexCreator::compareForDiffIdx);
        time_t sort = time(nullptr);
        cout<<"Time spent for sort kmers: "<<(double) (sort - start) << endl;
        size_t * uniqKmerIdx = new size_t[kmerBuffer.startIndexOfReserve + 1];
        size_t uniqKmerCnt = 0;
        reduceRedundancy(kmerBuffer, uniqKmerIdx, uniqKmerCnt, par, taxid2fasta);
        time_t reduction = time(nullptr);
        cout<<"Time spent for reducing redundancy: "<<(double) (reduction - sort) << endl;
        if(processedSplitCnt == numOfSplits && numOfFlush == 0){
            writeTargetFilesAndSplits(kmerBuffer.buffer, kmerBuffer.startIndexOfReserve, par, uniqKmerIdx, uniqKmerCnt);
        } else {
            writeTargetFiles(kmerBuffer.buffer, kmerBuffer.startIndexOfReserve, par,uniqKmerIdx, uniqKmerCnt);
        }
    }

    delete[] splitChecker;
}

// build_dir
size_t IndexCreator::fillTargetKmerBuffer(TargetKmerBuffer & kmerBuffer,
                                           bool * checker,
                                           size_t & processedSplitCnt,
                                           const vector<FastaSplit> & splits,
                                           const vector<TaxId2Fasta> & taxid2fasta,
                                           const LocalParameters & par) {
#ifdef OPENMP
    omp_set_num_threads(par.threads);
#endif
    bool hasOverflow = false;

#pragma omp parallel default(none), shared(par, checker, hasOverflow, splits, taxid2fasta, kmerBuffer, processedSplitCnt, cout)
    {
        ProdigalWrapper prodigal;
        SeqIterator seqIterator(par);
        size_t posToWrite;
        size_t numOfBlocks = 0;
        size_t kmerCntOfCurrSplit;
        vector<uint64_t> intergenicKmerList;
        vector<PredictedBlock> blocks;
        priority_queue<uint64_t> standardList;
        size_t lengthOfTrainingSeq;
        vector<Sequence> sequences;

#pragma omp for schedule(dynamic, 1)
        for (size_t i = 0; i < splits.size(); i++) {
            if (!checker[i] && (!hasOverflow)) {
                // Init
                intergenicKmerList.clear();
                standardList = priority_queue<uint64_t>();
                sequences.clear();

                // Estimate the number of k-mers to be extracted from fasta files of current splits
                kmerCntOfCurrSplit = estimateKmerNum(taxid2fasta, splits[i]);
                cout<<i<<" "<<taxid2fasta[splits[i].offset].fasta<<" "<<kmerCntOfCurrSplit<<endl;

                posToWrite = kmerBuffer.reserveMemory(kmerCntOfCurrSplit);
                if (posToWrite + kmerCntOfCurrSplit < kmerBuffer.bufferSize) {
                    // Load FASTA file for training
                    struct MmapedData<char> fastaForTraining = mmapData<char>(
                            taxid2fasta[splits[i].training].fasta.c_str());
                    getSeqSegmentsWithHead(sequences, fastaForTraining);
                    sort(sequences.begin(), sequences.end(),
                         [](const Sequence &a, const Sequence &b) { return a.length > b.length; });

                    // Train Prodigal with a training sequence of i th split
                    kseq_buffer_t buffer(const_cast<char *>(&fastaForTraining.data[sequences[0].start]),
                                         sequences[0].length);
                    kseq_t *seq = kseq_init(&buffer);
                    kseq_read(seq);
                    lengthOfTrainingSeq = strlen(seq->seq.s);
                    prodigal.is_meta = 0;
                    if (strlen(seq->seq.s) < 100000) {
                        prodigal.is_meta = 1;
                        prodigal.trainMeta(seq->seq.s);
                    } else {
                        prodigal.trainASpecies(seq->seq.s);
                    }

                    // Generate intergenic 23-mer list. It is used to strandness of intergenic sequences.
                    prodigal.getPredictedGenes(seq->seq.s);
                    seqIterator.generateIntergenicKmerList(prodigal.genes, prodigal.nodes,
                                                           prodigal.getNumberOfPredictedGenes(), intergenicKmerList,
                                                           seq->seq.s);

                    // Get min k-mer hash list for determining strandness
                    seqIterator.getMinHashList(standardList, seq->seq.s);

                    // Extract k-mer from the training sequence.
                    prodigal.removeCompletelyOverlappingGenes();
                    blocks.clear();
                    numOfBlocks = 0;
                    seqIterator.getExtendedORFs(prodigal.finalGenes, prodigal.nodes, blocks,
                                                      prodigal.fng, strlen(seq->seq.s),
                                                      numOfBlocks, intergenicKmerList, seq->seq.s);
                    for (size_t bl = 0; bl < numOfBlocks; bl++) {
                        seqIterator.translateBlock(seq->seq.s, blocks[bl]);
                        seqIterator.fillBufferWithKmerFromBlock(blocks[bl], seq->seq.s, kmerBuffer,
                                                                posToWrite, splits[i].training,
                                                                taxid2fasta[splits[i].training].species);
                    }
                    kseq_destroy(seq);

                    // Extract k-mers from the sequences which are not used for training but in the same file.
                    // Ex) plasmid
                    if(sequences.size() > 1) {
                        extractKmerFromFasta(seqIterator, fastaForTraining, standardList, lengthOfTrainingSeq,
                                             sequences, prodigal, intergenicKmerList, kmerBuffer, posToWrite,
                                             splits[i].training,taxid2fasta[splits[i].training].species, 1);

                    }
                    munmap(fastaForTraining.data, fastaForTraining.fileSize + 1);
                    sequences.clear();

                    // For other FASTA files...
                    for (size_t fastaCnt = 1; fastaCnt < splits[i].cnt; fastaCnt++) {
                        struct MmapedData<char> seqFile = mmapData<char>(taxid2fasta[splits[i].offset + fastaCnt].fasta.c_str());
                        getSeqSegmentsWithHead(sequences, seqFile);
                        sort(sequences.begin(), sequences.end(),
                             [](const Sequence &a, const Sequence &b) { return a.length > b.length; });
                        extractKmerFromFasta(seqIterator, seqFile, standardList, lengthOfTrainingSeq, sequences,
                                    prodigal, intergenicKmerList, kmerBuffer, posToWrite, splits[i].offset + fastaCnt,
                                    taxid2fasta[splits[i].training].species, 0);
                        munmap(seqFile.data, seqFile.fileSize + 1);
                        sequences.clear();
                    }
                    checker[i] = true;
                    __sync_fetch_and_add(&processedSplitCnt, 1);
                }else {
                    // Withdraw the reservation if the buffer is full.
                    hasOverflow = true;
                    __sync_fetch_and_sub(&kmerBuffer.startIndexOfReserve, kmerCntOfCurrSplit);
                    cout << "buffer is full" << endl;
                }
            }
        }
    }
    cout << "before return: " << kmerBuffer.startIndexOfReserve << endl;
    return 0;
}

void IndexCreator::extractKmerFromFasta(SeqIterator & seqIterator, MmapedData<char> & seqFile, priority_queue<uint64_t> & standardList,
                 size_t lengthOfTrainingSeq, const vector<Sequence> & sequences, ProdigalWrapper & prodigal,
                 vector<uint64_t> & intergenicKmerList, TargetKmerBuffer & kmerBuffer, size_t & posToWrite,
                 uint32_t seqID, int taxIdAtRank, size_t startIdx){
    priority_queue<uint64_t> currentList;
    vector<PredictedBlock> blocks;
    kseq_buffer_t buffer;
    kseq_t *seq;
    char *reverseCompliment;
    // For each assembly
    for (size_t assembly = startIdx; assembly < sequences.size(); assembly++) {
        buffer = {const_cast<char *>(&seqFile.data[sequences[assembly].start]),
                  static_cast<size_t>(sequences[assembly].length)};
        seq = kseq_init(&buffer);
        kseq_read(seq);
        size_t numOfBlocks = 0;
        blocks.clear();
        // Get extended ORFs
        seqIterator.getMinHashList(currentList, seq->seq.s);
        if (seqIterator.compareMinHashList(standardList, currentList, lengthOfTrainingSeq,
                                           strlen(seq->seq.s))) {
            prodigal.getPredictedGenes(seq->seq.s);
            prodigal.removeCompletelyOverlappingGenes();
            seqIterator.getExtendedORFs(prodigal.finalGenes, prodigal.nodes, blocks,
                                              prodigal.fng, strlen(seq->seq.s),
                                              numOfBlocks, intergenicKmerList, seq->seq.s);
            for (size_t bl = 0; bl < numOfBlocks; bl++) {
                seqIterator.translateBlock(seq->seq.s, blocks[bl]);
                seqIterator.fillBufferWithKmerFromBlock(blocks[bl], seq->seq.s, kmerBuffer,
                                                        posToWrite, seqID,
                                                        taxIdAtRank);
            }
        } else {
            reverseCompliment = seqIterator.reverseCompliment(seq->seq.s, strlen(seq->seq.s));
            prodigal.getPredictedGenes(reverseCompliment);
            prodigal.removeCompletelyOverlappingGenes();
            seqIterator.getExtendedORFs(prodigal.finalGenes, prodigal.nodes, blocks,
                                              prodigal.fng, strlen(reverseCompliment),
                                              numOfBlocks, intergenicKmerList, reverseCompliment);
            for (size_t bl = 0; bl < numOfBlocks; bl++) {
                seqIterator.translateBlock(reverseCompliment, blocks[bl]);
                seqIterator.fillBufferWithKmerFromBlock(blocks[bl], reverseCompliment, kmerBuffer,
                                                        posToWrite, seqID,
                                                        taxIdAtRank); //splits[i].offset + seqIdx
            }
            free(reverseCompliment);
        }
        currentList = priority_queue<uint64_t>();
        kseq_destroy(seq);
    }

}

// This function sort the TargetKmerBuffer, do redundancy reducing task, write the differential index of them
void IndexCreator::writeTargetFiles(TargetKmer * kmerBuffer, size_t & kmerNum, const LocalParameters & par,
                                    const size_t * uniqeKmerIdx, size_t & uniqKmerCnt){
    string diffIdxFileName;
    string infoFileName;
    diffIdxFileName = par.filenames[1] + "/" + to_string(numOfFlush) + "_diffIdx";
    infoFileName = par.filenames[1] + "/" + to_string(numOfFlush) + "_info";

    FILE * diffIdxFile = fopen(diffIdxFileName.c_str(), "wb");
    FILE * infoFile = fopen(infoFileName.c_str(), "wb");
    if (diffIdxFile == nullptr || infoFile == nullptr){
        cout<<"Cannot open the file for writing target DB"<<endl;
        return;
    }
    numOfFlush++;

    uint16_t *diffIdxBuffer = (uint16_t *)malloc(sizeof(uint16_t) * 10'000'000'000);
    size_t localBufIdx = 0;
    uint64_t lastKmer = 0;
    size_t write = 0;

    for(size_t i = 0; i < uniqKmerCnt ; i++) {
        fwrite(& kmerBuffer[uniqeKmerIdx[i]].info, sizeof (TargetKmerInfo), 1, infoFile);
        write++;
        getDiffIdx(lastKmer, kmerBuffer[uniqeKmerIdx[i]].ADkmer, diffIdxFile, diffIdxBuffer, localBufIdx);
        lastKmer = kmerBuffer[uniqeKmerIdx[i]].ADkmer;
    }

    cout<<"total k-mer count  : "<< kmerNum << endl;
    cout<<"written k-mer count: "<<write<<endl;

    flushKmerBuf(diffIdxBuffer, diffIdxFile, localBufIdx);
    free(diffIdxBuffer);
    fclose(diffIdxFile);
    fclose(infoFile);
    kmerNum = 0;
}

void IndexCreator::writeTargetFilesAndSplits(TargetKmer * kmerBuffer, size_t & kmerNum, const LocalParameters & par,
                                    const size_t * uniqeKmerIdx, size_t & uniqKmerCnt){
    string diffIdxFileName;
    string infoFileName;
    string splitFileName;

    diffIdxFileName = par.filenames[1] + "/diffIdx";
    infoFileName = par.filenames[1] + "/info";
    splitFileName = par.filenames[1] + "/split";

    // Make splits
    FILE * diffIdxSplitFile = fopen(splitFileName.c_str(), "wb");
    DiffIdxSplit splitList[SplitNum];
    memset(splitList, 0, sizeof(DiffIdxSplit) * SplitNum);
    size_t splitWidth = uniqKmerCnt / par.threads;
    for (size_t i = 1; i < (size_t) par.threads; i++) {
        for (size_t j = uniqeKmerIdx[0] + splitWidth * i; j + 1 < uniqKmerCnt; j++) {
            if (AminoAcidPart(kmerBuffer[j].ADkmer) != AminoAcidPart(kmerBuffer[j + 1].ADkmer)) {
                splitList[i].ADkmer = kmerBuffer[j].ADkmer;
                break;
            }
        }
    }

    FILE * diffIdxFile = fopen(diffIdxFileName.c_str(), "wb");
    FILE * infoFile = fopen(infoFileName.c_str(), "wb");
    if (diffIdxFile == nullptr || infoFile == nullptr){
        cout<<"Cannot open the file for writing target DB"<<endl;
        return;
    }
    numOfFlush++;


    uint16_t *diffIdxBuffer = (uint16_t *)malloc(sizeof(uint16_t) * size_t(kmerBufSize));
    size_t localBufIdx = 0;
    uint64_t lastKmer = 0;
    size_t write = 0;

    size_t splitIdx = 1;
    size_t totalDiffIdx = 0;
    for(size_t i = 0; i < uniqKmerCnt ; i++) {
        fwrite(& kmerBuffer[uniqeKmerIdx[i]].info, sizeof (TargetKmerInfo), 1, infoFile);
        write++;
        getDiffIdx(lastKmer, kmerBuffer[uniqeKmerIdx[i]].ADkmer, diffIdxFile,
                   diffIdxBuffer, localBufIdx, totalDiffIdx);
        lastKmer = kmerBuffer[uniqeKmerIdx[i]].ADkmer;
        if(lastKmer == splitList[splitIdx].ADkmer){
            splitList[splitIdx].diffIdxOffset = totalDiffIdx;
            splitList[splitIdx].infoIdxOffset = write;
            splitIdx ++;
        }
    }

    cout<<"total k-mer count  : "<< kmerNum << endl;
    cout<<"written k-mer count: "<<write<<endl;

    flushKmerBuf(diffIdxBuffer, diffIdxFile, localBufIdx);
    fwrite(splitList, sizeof(DiffIdxSplit), SplitNum, diffIdxSplitFile);

    free(diffIdxBuffer);
    fclose(diffIdxSplitFile);
    fclose(diffIdxFile);
    fclose(infoFile);
    kmerNum = 0;
}

void IndexCreator::reduceRedundancy(TargetKmerBuffer & kmerBuffer, size_t * uniqeKmerIdx, size_t & uniqueKmerCnt,
                                    const LocalParameters & par) {
    // Find the first index of garbage k-mer (UINT64_MAX)
    for(size_t checkN = kmerBuffer.startIndexOfReserve - 1; checkN != 0; checkN--){
        if(kmerBuffer.buffer[checkN].ADkmer != UINT64_MAX){
            kmerBuffer.startIndexOfReserve = checkN + 1;
            break;
        }
    }

    // Find the first index of meaningful k-mer
    size_t startIdx = 0;
    for(size_t i = 0; i < kmerBuffer.startIndexOfReserve ; i++){
        if(kmerBuffer.buffer[i].taxIdAtRank != 0){
            startIdx = i;
            break;
        }
    }

    // Make splits
    vector<Split> splits;
    size_t splitWidth = (kmerBuffer.startIndexOfReserve - startIdx) / par.threads;
    for (size_t i = 0; i < par.threads - 1; i++) {
        for (size_t j = startIdx + splitWidth; j + 1 < kmerBuffer.startIndexOfReserve; j++) {
            if (kmerBuffer.buffer[j].taxIdAtRank != kmerBuffer.buffer[j + 1].taxIdAtRank) {
                splits.emplace_back(startIdx, j);
                startIdx = j + 1;
                break;
            }
        }
    }
    splits.emplace_back(startIdx, kmerBuffer.startIndexOfReserve - 1);

    //
    size_t ** idxOfEachSplit = new size_t * [splits.size()];
    size_t * cntOfEachSplit = new size_t[splits.size()];
    for(size_t i = 0; i < splits.size(); i++){
        idxOfEachSplit[i] = new size_t[splits[i].end - splits[i].offset + 2];
        cntOfEachSplit[i] = 0;
    }
#pragma omp parallel default(none), shared(kmerBuffer, idxOfEachSplit, cntOfEachSplit, splits)
    {
        TargetKmer * lookingKmer;
        size_t lookingIndex;
        int endFlag;
        int hasSeenOtherStrains;
#pragma omp for schedule(dynamic, 1)
        for(size_t split = 0; split < splits.size(); split ++){
            lookingKmer = & kmerBuffer.buffer[splits[split].offset];
            lookingIndex = splits[split].offset;
            endFlag = 0;
            for(size_t i = 1 + splits[split].offset; i < splits[split].end + 1 ; i++) {
                hasSeenOtherStrains = 0;
                while(lookingKmer->taxIdAtRank == kmerBuffer.buffer[i].taxIdAtRank){
                    if (lookingKmer->ADkmer != kmerBuffer.buffer[i].ADkmer) {
                        break;
                    }
                    hasSeenOtherStrains += (taxIdList[lookingKmer->info.sequenceID] != taxIdList[kmerBuffer.buffer[i].info.sequenceID]);
                    i++;
                    if(i == splits[split].end + 1){
                        endFlag = 1;
                        break;
                    }
                }

                lookingKmer->info.redundancy = (hasSeenOtherStrains > 0);
                idxOfEachSplit[split][cntOfEachSplit[split]] = lookingIndex;
                cntOfEachSplit[split] ++;
                if(endFlag == 1) break;
                lookingKmer = & kmerBuffer.buffer[i];
                lookingIndex = i;
            }

            //For the end part
            if(!((kmerBuffer.buffer[splits[split].end - 1].ADkmer == kmerBuffer.buffer[splits[split].end].ADkmer) &&
                 (kmerBuffer.buffer[splits[split].end - 1].taxIdAtRank == kmerBuffer.buffer[splits[split].end].taxIdAtRank))){
                idxOfEachSplit[split][cntOfEachSplit[split]] = splits[split].end;
                cntOfEachSplit[split] ++;
            }
        }
    }

    // Merge
    for(size_t i = 0; i < splits.size(); i++){
        memcpy(uniqeKmerIdx + uniqueKmerCnt, idxOfEachSplit[i], cntOfEachSplit[i] * sizeof(size_t));
        uniqueKmerCnt += cntOfEachSplit[i];
    }

    for(size_t i = 0; i < splits.size(); i++){
        delete[] idxOfEachSplit[i];
    }
    delete[] cntOfEachSplit;
}

void IndexCreator::reduceRedundancy(TargetKmerBuffer & kmerBuffer, size_t * uniqeKmerIdx, size_t & uniqueKmerCnt, const LocalParameters & par,
                                    const vector<TaxId2Fasta> & taxid2fasta){

    // Find the first index of garbage k-mer (UINT64_MAX)
    for(size_t checkN = kmerBuffer.startIndexOfReserve - 1; checkN >= 0; checkN--){
        if(kmerBuffer.buffer[checkN].ADkmer != UINT64_MAX){
            kmerBuffer.startIndexOfReserve = checkN + 1;
            break;
        }
    }

    // Find the first index of meaningful k-mer
    size_t startIdx = 0;
    for(size_t i = 0; i < kmerBuffer.startIndexOfReserve ; i++){
        if(kmerBuffer.buffer[i].taxIdAtRank != 0){
            startIdx = i;
            break;
        }
    }

    // Make splits
    vector<Split> splits;
    size_t splitWidth = (kmerBuffer.startIndexOfReserve - startIdx) / par.threads;
    for (size_t i = 0; i < par.threads - 1; i++) {
        for (size_t j = startIdx + splitWidth; j + 1 < kmerBuffer.startIndexOfReserve; j++) {
            if (kmerBuffer.buffer[j].taxIdAtRank != kmerBuffer.buffer[j + 1].taxIdAtRank) {
                splits.emplace_back(startIdx, j);
                startIdx = j + 1;
                break;
            }
        }
    }
    splits.emplace_back(startIdx, kmerBuffer.startIndexOfReserve - 1);

    //
    size_t ** idxOfEachSplit = new size_t * [splits.size()];
    size_t * cntOfEachSplit = new size_t[splits.size()];
    for(size_t i = 0; i < splits.size(); i++){
        idxOfEachSplit[i] = new size_t[splits[i].end - splits[i].offset + 2];
        cntOfEachSplit[i] = 0;
    }
#pragma omp parallel default(none), shared(kmerBuffer, taxid2fasta, idxOfEachSplit, cntOfEachSplit, splits)
    {
        TargetKmer * lookingKmer;
        size_t lookingIndex;
        int endFlag;
        int hasSeenOtherStrains;
#pragma omp for schedule(dynamic, 1)
        for(size_t split = 0; split < splits.size(); split ++){
            lookingKmer = & kmerBuffer.buffer[splits[split].offset];
            lookingIndex = splits[split].offset;
            endFlag = 0;
            for(size_t i = 1 + splits[split].offset; i < splits[split].end + 1 ; i++) {
                hasSeenOtherStrains = 0;
                while(lookingKmer->taxIdAtRank == kmerBuffer.buffer[i].taxIdAtRank){
                    if (lookingKmer->ADkmer != kmerBuffer.buffer[i].ADkmer) {
                        break;
                    }
                    hasSeenOtherStrains += (taxid2fasta[lookingKmer->info.sequenceID].taxid
                            != taxid2fasta[kmerBuffer.buffer[i].info.sequenceID].taxid);
                    i++;
                    if(i == splits[split].end + 1){
                        endFlag = 1;
                        break;
                    }
                }

                lookingKmer->info.redundancy = (hasSeenOtherStrains > 0);
                idxOfEachSplit[split][cntOfEachSplit[split]] = lookingIndex;
                cntOfEachSplit[split] ++;
                if(endFlag == 1) break;
                lookingKmer = & kmerBuffer.buffer[i];
                lookingIndex = i;
            }

            //For the end part
            if(!((kmerBuffer.buffer[splits[split].end - 1].ADkmer == kmerBuffer.buffer[splits[split].end].ADkmer) &&
                 (kmerBuffer.buffer[splits[split].end - 1].taxIdAtRank == kmerBuffer.buffer[splits[split].end].taxIdAtRank))){
                idxOfEachSplit[split][cntOfEachSplit[split]] = splits[split].end;
                cntOfEachSplit[split] ++;
            }
        }
    }

    // Merge
    for(size_t i = 0; i < splits.size(); i++){
        memcpy(uniqeKmerIdx + uniqueKmerCnt, idxOfEachSplit[i], cntOfEachSplit[i] * sizeof(size_t));
        uniqueKmerCnt += cntOfEachSplit[i];
    }

    for(size_t i = 0; i < splits.size(); i++){
        delete[] idxOfEachSplit[i];
    }
    delete[] cntOfEachSplit;
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

void IndexCreator::getDiffIdx(const uint64_t & lastKmer, const uint64_t & entryToWrite, FILE* handleKmerTable,
                              uint16_t *kmerBuf, size_t & localBufIdx, size_t & totalBufferIdx){
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
    totalBufferIdx += 4 - idx;
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

string IndexCreator::getSeqSegmentsWithHead(vector<Sequence> & seqSegments, const string & seqFileName) {
    struct stat stat1{};
    stat(seqFileName.c_str(), &stat1);
    size_t numOfChar = stat1.st_size;
    cout << seqFileName << ": " << numOfChar << endl;
    string firstLine;
    ifstream seqFile;
    seqFile.open(seqFileName);
    string eachLine;
    size_t start = 0;
    size_t pos;
    vector<Sequence> seqSegmentsTmp;
    if (seqFile.is_open()) {
        getline(seqFile, firstLine, '\n');
        while (getline(seqFile, eachLine, '\n')) {
            if (eachLine[0] == '>') {
                pos = (size_t) seqFile.tellg();
                seqSegmentsTmp.emplace_back(start, pos - eachLine.length() - 3,pos - eachLine.length() - start - 2);
                start = pos - eachLine.length() - 1;
            }
        }
        seqSegmentsTmp.emplace_back(start, numOfChar - 2, numOfChar - start - 1);
    } else {
        cerr << "Unable to open file: " << seqFileName << endl;
//        EXIT(EXIT_FAILURE);
    }
    seqFile.close();
    seqSegments = move(seqSegmentsTmp);
    // TODO SORT?
    return firstLine;
}

void IndexCreator::groupFastaFiles(const vector<TaxId2Fasta> & taxIdListAtRank, vector<FastaSplit> & fastaSplit){
    size_t i = 0;
    size_t training = 0;
    uint32_t offset = 0;
    uint32_t splitSize = 0;
    while(i < taxIdListAtRank.size()){
        TaxID currentSpecies = taxIdListAtRank[i].species;
        training = i;
        offset = i;
        splitSize = 0;
        while(currentSpecies == taxIdListAtRank[i].species && i < taxIdListAtRank.size() && splitSize < 30){
            splitSize ++;
            i++;
        }
        fastaSplit.emplace_back(training, offset, splitSize, currentSpecies);
    }
    for(auto x : fastaSplit){
        cout<<x.training<<" "<<x.offset<<" "<<x.cnt<<endl;
    }
}
void IndexCreator::splitAFastaFile(const vector<int> & speciesTaxIDs,
                                   vector<FastaSplit> & fastaSplit, // Each split is processed by a thread
                                   vector<Sequence> & seqSegments // Start and end position of each sequence
                                   ){
    size_t training = 0;
    size_t idx = 0;
    uint32_t offset = 0;
    uint32_t cnt = 0;
    int currSpecies;
    unordered_map<TaxID, size_t> species2training;
    unordered_map<TaxID, int> maxLenOfSpecies;

    while(idx < speciesTaxIDs.size()){
        offset = idx;
        cnt = 0;
        currSpecies = speciesTaxIDs[idx];
        if(currSpecies == 0) {
            idx++;
            continue;
        }
        if(maxLenOfSpecies.find(currSpecies) == maxLenOfSpecies.end()){ // First time to see this species
            maxLenOfSpecies[currSpecies] = 0;
        }
        while(idx < speciesTaxIDs.size() && currSpecies == speciesTaxIDs[idx] ){
            if(seqSegments[idx].length > maxLenOfSpecies[currSpecies]){ // Find the longest sequence of this species to be the training sequence
                maxLenOfSpecies[currSpecies] = seqSegments[idx].length;
                species2training[currSpecies] = idx;
            }
            cnt ++;
            if(cnt > 30){
                fastaSplit.emplace_back(0, offset, cnt - 1, currSpecies);
                offset += cnt - 1;
                cnt = 1;
            }
            idx ++;
        }
        fastaSplit.emplace_back(0, offset, cnt, currSpecies);
    }

    // Update the training sequence
    for(auto & x : fastaSplit){
        x.training = species2training[x.taxid];
    }
}

void IndexCreator::mappingFromTaxIDtoFasta(const string & fastaList_fname,
                             unordered_map<string, int> & assacc2taxid,
                             vector<TaxId2Fasta> & taxid2fasta,
                             NcbiTaxonomy & taxonomy){
    ifstream fastaList;
    fastaList.open(fastaList_fname);
    string fileName;
    smatch assacc;
    int taxId;
    regex regex1("(GC[AF]_[0-9]*\\.[0-9]*)");
    if(fastaList.is_open()){
        cout<<"Writing taxID to fileName mapping file"<<endl;
        while(getline(fastaList,fileName,'\n')) {
            regex_search(fileName, assacc, regex1);
            if (assacc2taxid.count(assacc[0].str())) {
                taxId = assacc2taxid[assacc[0].str()];
                taxid2fasta.emplace_back(taxonomy.getTaxIdAtRank(taxId,"species"), taxId, fileName);
            } else{
                cout<<assacc[0].str()<<" is excluded in creating target DB because it is not mapped to taxonomical ID"<<endl;
            }
        }
    }

    // Sort by species tax id
    sort(taxid2fasta.begin(), taxid2fasta.end(), [](TaxId2Fasta & a, TaxId2Fasta & b){return a.species < b.species;});
}

void IndexCreator::load_assacc2taxid(const string & mappingFile, unordered_map<string, int> & assacc2taxid){
    string key, value;
    ifstream map;
    map.open(mappingFile);
    if(map.is_open()){
        while(getline(map,key,'\t')){
            getline(map, value, '\n');
            assacc2taxid[key] = stoi(value);
            if(key[2] == 'F'){
                key[2] = 'A';
                assacc2taxid[key] = stoi(value);
            }
        }
    } else{
        cerr<<"Cannot open file for mappig from assemlby accession to tax ID"<<endl;
    }
    map.close();
}

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

    // Sort for differential indexing
    SORT_PARALLEL(kmerBuffer, kmerBuffer + kmerNum, IndexCreator::compareForDiffIdx);

    // Find the first index of garbage k-mer (UINT64_MAX)
    for(size_t checkN = kmerNum - 1; checkN >= 0; checkN--){
        if(kmerBuffer[checkN].ADkmer != UINT64_MAX){
            kmerNum = checkN + 1;
            break;
        }
    }

    // Find the first index of meaningful k-mer
    size_t startIdx = 0;
    for(size_t i = 0; i < kmerNum ; i++){
        if(kmerBuffer[i].taxIdAtRank != 0){
            startIdx = i;
            break;
        }
    }

    // Redundancy reduction task
    TargetKmer lookingKmer = kmerBuffer[startIdx];
    size_t write = 0;
    int endFlag = 0;
    int hasSeenOtherStrains;

    for(size_t i = startIdx + 1; i < kmerNum; i++) {
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

size_t IndexCreator::fillTargetKmerBuffer(TargetKmerBuffer &kmerBuffer,
                                          bool *checker,
                                          size_t &processedSplitCnt,
                                          const LocalParameters &par) {
#ifdef OPENMP
    omp_set_num_threads(par.threads);
#endif
    bool hasOverflow = false;

#pragma omp parallel default(none), shared(kmerBuffer, checker, processedSplitCnt, hasOverflow, par, cout)
    {
        ProdigalWrapper prodigal;
        SeqIterator seqIterator(par);
        size_t posToWrite;
        size_t orfNum;
        vector<uint64_t> intergenicKmerList;
        vector<PredictedBlock> extendedORFs;
        priority_queue<uint64_t> standardList;
        priority_queue<uint64_t> currentList;
        size_t lengthOfTrainingSeq;
        char *reverseCompliment;
        vector<bool> strandness;
        kseq_buffer_t buffer;
        kseq_t *seq;
#pragma omp for schedule(dynamic, 1)
        for (size_t i = 0; i < fnaSplits.size(); i++) {
            if (!checker[i] && !hasOverflow) {
                intergenicKmerList.clear();
                strandness.clear();
                standardList = priority_queue<uint64_t>();

                // Estimate the number of k-mers to be extracted from current split
                size_t totalLength = 0;
                for (size_t p = 0; p < fnaSplits[i].cnt; p++) {
                    totalLength += sequenceOfFastas[fnaSplits[i].file_idx][fnaSplits[i].offset + p].length;
                }
                size_t estimatedKmerCnt = totalLength / 3;

                // Process current split if buffer has enough space.
                posToWrite = kmerBuffer.reserveMemory(estimatedKmerCnt);
                if (posToWrite + estimatedKmerCnt < kmerBuffer.bufferSize) {
                    // MMap FASTA file of current split
                    struct MmapedData<char> fastaFile = mmapData<char>(fnaList[fnaSplits[i].file_idx].c_str());

                    // Train Prodigal with a training sequence of i th split
                    buffer = {const_cast<char *>(&fastaFile.data[sequenceOfFastas[fnaSplits[i].file_idx][fnaSplits[i].training].start]),
                              static_cast<size_t>(sequenceOfFastas[fnaSplits[i].file_idx][fnaSplits[i].training].length)};
                    kseq_t *seq = kseq_init(&buffer);
                    kseq_read(seq);
                    lengthOfTrainingSeq = strlen(seq->seq.s);
                    prodigal.is_meta = 0;
                    if (strlen(seq->seq.s) < 100000) {
                        prodigal.is_meta = 1;
                        prodigal.trainMeta(seq->seq.s);
                    } else {
                        prodigal.trainASpecies(seq->seq.s);
                    }

                    // Generate intergenic 23-mer list. It is used to strandness of intergenic sequences.
                    prodigal.getPredictedGenes(seq->seq.s);
                    seqIterator.generateIntergenicKmerList(prodigal.genes, prodigal.nodes,
                                                           prodigal.getNumberOfPredictedGenes(), intergenicKmerList,
                                                           seq->seq.s);

                    // Get min k-mer hash list for determining strandness
                    seqIterator.getMinHashList(standardList, seq->seq.s);
                    kseq_destroy(seq);

                    // Extract k-mer from the rest of the sequences of current split
                    for (size_t s_cnt = 0; s_cnt < fnaSplits[i].cnt; ++s_cnt) {
                        buffer = {const_cast<char *>(&fastaFile.data[sequenceOfFastas[fnaSplits[i].file_idx][
                                fnaSplits[i].offset + s_cnt].start]),
                                  static_cast<size_t>(sequenceOfFastas[fnaSplits[i].file_idx][fnaSplits[i].offset +
                                                                                              s_cnt].length)};
                        seq = kseq_init(&buffer);
                        kseq_read(seq);

                        seqIterator.getMinHashList(currentList, seq->seq.s);
                        orfNum = 0;
                        extendedORFs.clear();
                        if (seqIterator.compareMinHashList(standardList, currentList, lengthOfTrainingSeq,
                                                           strlen(seq->seq.s))) {
                            prodigal.getPredictedGenes(seq->seq.s);
                            prodigal.removeCompletelyOverlappingGenes();
                            seqIterator.getExtendedORFs(prodigal.finalGenes, prodigal.nodes, extendedORFs,
                                                             prodigal.fng, strlen(seq->seq.s),
                                                        orfNum, intergenicKmerList, seq->seq.s);
                            for (size_t orfCnt = 0; orfCnt < orfNum; orfCnt++) {
                                seqIterator.translateBlock(seq->seq.s, extendedORFs[orfCnt]);
                                seqIterator.fillBufferWithKmerFromBlock(extendedORFs[orfCnt], seq->seq.s, kmerBuffer, posToWrite,
                                                                        fnaSplits[i].file_idx, fnaSplits[i].speciesID);
                            }
                        } else {
                            reverseCompliment = seqIterator.reverseCompliment(seq->seq.s, strlen(seq->seq.s));
                            prodigal.getPredictedGenes(reverseCompliment);
                            prodigal.removeCompletelyOverlappingGenes();
                            seqIterator.getExtendedORFs(prodigal.finalGenes, prodigal.nodes, extendedORFs,
                                                             prodigal.fng, strlen(reverseCompliment),
                                                        orfNum, intergenicKmerList, reverseCompliment);
                            for (size_t orfCnt = 0; orfCnt < orfNum; orfCnt++) {
                                seqIterator.translateBlock(reverseCompliment, extendedORFs[orfCnt]);
                                seqIterator.fillBufferWithKmerFromBlock(extendedORFs[orfCnt], reverseCompliment, kmerBuffer,
                                                                        posToWrite, fnaSplits[i].file_idx,fnaSplits[i].speciesID);
                            }
                            free(reverseCompliment);
                        }
                        currentList = priority_queue<uint64_t>();
                        kseq_destroy(seq);
                    }
                    checker[i] = true;
                    __sync_fetch_and_add(&processedSplitCnt, 1);
                    munmap(fastaFile.data, fastaFile.fileSize + 1);
                }else {
                    // Withdraw the reservation if the buffer is full.
                    hasOverflow = true;
                    __sync_fetch_and_sub(&kmerBuffer.startIndexOfReserve, estimatedKmerCnt);
                }
            }
        }
    }

    cout << "Before return: " << kmerBuffer.startIndexOfReserve << endl;
    return 0;
}

size_t IndexCreator::fillTargetKmerBuffer(TargetKmerBuffer &kmerBuffer,
                                          MmapedData<char> &seqFile,
                                          vector<Sequence> &seqs,
                                          bool *checker,
                                          size_t &processedSplitCnt,
                                          const vector<FastaSplit> &splits,
                                          const vector<int> &taxIdListAtRank,
                                          const vector<TaxID> & superkingdoms,
                                          const LocalParameters &par) {
#ifdef OPENMP
    omp_set_num_threads(par.threads);
#endif
    bool hasOverflow = false;

#pragma omp parallel default(none), shared(checker, hasOverflow,par, splits, seqFile, seqs, kmerBuffer, \
processedSplitCnt, cout, taxIdListAtRank, superkingdoms)
    {
        ProdigalWrapper prodigal;
        SeqIterator seqIterator(par);
        size_t posToWrite;
        size_t numOfBlocks;
        vector<uint64_t> intergenicKmerList;
        vector<PredictedBlock> blocks;
        priority_queue<uint64_t> standardList;
        priority_queue<uint64_t> currentList;
        size_t lengthOfTrainingSeq;
        char * reverseCompliment;
        vector<bool> strandness;
        kseq_buffer_t buffer;
        kseq_t * seq;
#pragma omp for schedule(dynamic, 1)
        for (size_t i = 0; i < splits.size() ; i++) {
            if(!checker[i] && !hasOverflow) {
                intergenicKmerList.clear();
                strandness.clear();
                standardList = priority_queue<uint64_t>();

                // Estimate the number of k-mers to be extracted from current split
                size_t totalLength = 0;
                for(size_t p = 0; p < splits[i].cnt; p++) {
                    totalLength += seqs[splits[i].offset + p].length;
                }
                size_t estimatedKmerCnt = totalLength / 3;

                // Process current split if buffer has enough space.
                posToWrite = kmerBuffer.reserveMemory(estimatedKmerCnt);
                if (posToWrite + estimatedKmerCnt < kmerBuffer.bufferSize){
                    // Update prodigal model
                    prodigal.updateTrainingInfo(trainingInfo[splits[i].taxid]);

                    buffer = {const_cast<char *>(&seqFile.data[seqs[splits[i].training].start]),
                              static_cast<size_t>(seqs[splits[i].training].length)};
                    seq = kseq_init(&buffer);
                    kseq_read(seq);
                    lengthOfTrainingSeq = strlen(seq->seq.s);

                    // Generate intergenic 23-mer list for current split
                    prodigal.getPredictedGenes(seq->seq.s);
                    seqIterator.generateIntergenicKmerList(prodigal.genes, prodigal.nodes,
                                                           prodigal.getNumberOfPredictedGenes(),
                                                           intergenicKmerList, seq->seq.s);

                    // Get minimum k-mer hash list for determining strandness
                    seqIterator.getMinHashList(standardList, seq->seq.s);

//                    // Extract k-mer from the training sequence.
//                    prodigal.removeCompletelyOverlappingGenes();
//                    blocks.clear();
//                    numOfBlocks = 0;
//                    seqIterator.getTranslationBlocks(prodigal.finalGenes, prodigal.nodes, blocks,
//                                                     prodigal.fng, strlen(seq->seq.s),
//                                                     numOfBlocks, intergenicKmerList, seq->seq.s);
//                    for (size_t bl = 0; bl < numOfBlocks; bl++) {
//                        seqIterator.translateBlock(seq->seq.s, blocks[bl]);
//                        seqIterator.fillBufferWithKmerFromBlock(blocks[bl], seq->seq.s, kmerBuffer,
//                                                                posToWrite, splits[i].training,
//                                                                taxIdListAtRank[splits[i].training]);
//
//
//                    }
                    kseq_destroy(seq);

                    // Extract k-mer from the rest of the split
                    for (size_t p = 0; p < splits[i].cnt; p++) {
                        if (splits[i].offset + p == splits[i].training) continue;
                        buffer = {const_cast<char *>(&seqFile.data[seqs[splits[i].offset + p].start]),
                                  static_cast<size_t>(seqs[splits[i].offset + p].length)};
                        seq = kseq_init(&buffer);
                        kseq_read(seq);
                        seqIterator.getMinHashList(currentList, seq->seq.s);
                        numOfBlocks = 0;
                        blocks.clear();
                        if (seqIterator.compareMinHashList(standardList, currentList, lengthOfTrainingSeq,
                                                           strlen(seq->seq.s))) {
                            prodigal.getPredictedGenes(seq->seq.s);
                            prodigal.removeCompletelyOverlappingGenes();
                            seqIterator.getExtendedORFs(prodigal.finalGenes, prodigal.nodes, blocks,
                                                             prodigal.fng, strlen(seq->seq.s),
                                                             numOfBlocks, intergenicKmerList, seq->seq.s);
                            for (size_t bl = 0; bl < numOfBlocks; bl++) {
                                seqIterator.translateBlock(seq->seq.s, blocks[bl]);
                                seqIterator.fillBufferWithKmerFromBlock(blocks[bl], seq->seq.s, kmerBuffer, posToWrite,
                                                                        splits[i].offset + p,
                                                                        taxIdListAtRank[splits[i].offset +
                                                                                        p]); //splits[i].offset + seqIdx
                            }

                        } else {
                            reverseCompliment = seqIterator.reverseCompliment(seq->seq.s, strlen(seq->seq.s));
                            prodigal.getPredictedGenes(reverseCompliment);
                            prodigal.removeCompletelyOverlappingGenes();
                            seqIterator.getExtendedORFs(prodigal.finalGenes, prodigal.nodes, blocks,
                                                             prodigal.fng, strlen(reverseCompliment),
                                                             numOfBlocks, intergenicKmerList, reverseCompliment);
                            for (size_t bl = 0; bl < numOfBlocks; bl++) {
                                seqIterator.translateBlock(reverseCompliment, blocks[bl]);
                                seqIterator.fillBufferWithKmerFromBlock(blocks[bl], reverseCompliment, kmerBuffer,
                                                                        posToWrite, splits[i].offset + p,
                                                                        taxIdListAtRank[splits[i].offset +
                                                                                        p]); //splits[i].offset + seqIdx
                            }
                            free(reverseCompliment);
                        }
                        currentList = priority_queue<uint64_t>();
                        kseq_destroy(seq);
                    }
                    checker[i] = true;
                    __sync_fetch_and_add(&processedSplitCnt, 1);
                } else {
                    // Withdraw the reservation if the buffer is full.
                    hasOverflow = true;
                    __sync_fetch_and_sub(&kmerBuffer.startIndexOfReserve, estimatedKmerCnt);
                    cout << "buffer is full" << endl;
                }
            }
        }
    }
    cout<<"before return: "<<kmerBuffer.startIndexOfReserve<<endl;
    return 0;
}

size_t IndexCreator::estimateKmerNum(const vector<TaxId2Fasta> & taxid2fasta, const FastaSplit & split){
    struct stat stat1{};
    size_t kmerNum = 0;
    for(size_t i = split.offset ; i < split.offset + split.cnt; i ++){
        stat(taxid2fasta[i].fasta.c_str(), &stat1);
        size_t charNum = stat1.st_size;
        kmerNum += charNum/3 - kmerLength + 1;
    }
    return kmerNum;
}

void IndexCreator::trainProdigal(const vector<FastaSplit> &splits, const vector<Sequence> & seqs,
                                 MmapedData<char> &seqFile, const LocalParameters &par) {
    // Train prodigal for species of each split.
    ProdigalWrapper prodigal;
    kseq_buffer_t buffer;
    kseq_t * seq;
    size_t lengthOfTrainingSeq;
    for (size_t i = 0; i < splits.size(); i++){
        TaxID currentSpecies = splits[i].taxid;
        // Check if the species was trained before.
        if (std::find(trainedSpecies.begin(), trainedSpecies.end(), currentSpecies) != trainedSpecies.end()) {
            continue;
        }

        // Load sequence for species.
        buffer = {const_cast<char *>(&seqFile.data[seqs[splits[i].training].start]),
                  static_cast<size_t>(seqs[splits[i].training].length)};
        seq = kseq_init(&buffer);
        kseq_read(seq);

        // Train prodigal.
        prodigal.is_meta = 0;
        lengthOfTrainingSeq = strlen(seq->seq.s);
        if(lengthOfTrainingSeq < 100000){
            prodigal.is_meta = 1;
            prodigal.trainMeta(seq->seq.s);
        }else{
            prodigal.trainASpecies(seq->seq.s);
        }

        // Write training result into a file.
        string fileName = to_string(currentSpecies)+ ".tinfo";
        _training * trainingInfo = prodigal.getTrainingInfo();
        write_training_file( const_cast<char *>(fileName.c_str()), trainingInfo);

        // Add species to trainedSpecies.
        trainedSpecies.push_back(currentSpecies);
    }
    // Write trained species into a file.
    FILE * fp = fopen(tinfo_path.c_str(), "w");
    for (int trainedSpecie : trainedSpecies){
        fprintf(fp, "%d", trainedSpecie);
    }
    fclose(fp);
}

void IndexCreator::loadTrainingInfo() {
    // Load training info of each species.
    for (TaxID taxId : trainedSpecies){
        string fileName = to_string(taxId)+ ".tinfo";
        _training tinfo{};
        read_training_file(const_cast<char *>(fileName.c_str()), &tinfo);
        trainingInfo[taxId] = tinfo;
    }
}



