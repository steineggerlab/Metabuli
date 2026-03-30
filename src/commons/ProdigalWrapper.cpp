#include "ProdigalWrapper.h"
#include "SeqIterator.h"
#include <iostream>
ProdigalWrapper::~ProdigalWrapper() {
    for(size_t i = 0; i < NUM_META; i++){
        delete meta[i].tinf;
    }
    free(meta);
    free(seq);
    free(rseq);
    free(useq);
    free(nodes);
    free(genes);
    free(finalGenes);
}
ProdigalWrapper::ProdigalWrapper() {
    seq = (unsigned char *)malloc((MAX_SEQ/4 + 1)*sizeof(unsigned char)); // 8 Mb
    rseq = (unsigned char *)malloc((MAX_SEQ/4 + 1)*sizeof(unsigned char)); // 8 Mb
    useq = (unsigned char *)malloc((MAX_SEQ/8 + 1)*sizeof(unsigned char)); // 4 Mb
    nodes = (struct _node *)malloc(STT_NOD*sizeof(struct _node)); // 13.6 Mb
    genes = (struct _gene *)malloc(MAX_GENES*sizeof(struct _gene)); // 30 Mb
    finalGenes = (struct _gene *)malloc(MAX_GENES*sizeof(struct _gene));
    if(seq == NULL || rseq == NULL || nodes == NULL || genes == NULL) {
        fprintf(stderr, "\nError: Malloc failed on sequence/orfs\n\n"); exit(1);
    }

    memset(seq, 0, (MAX_SEQ/4 + 1)*sizeof(unsigned char));
    memset(rseq, 0, (MAX_SEQ/4 + 1)*sizeof(unsigned char));
    memset(useq, 0, (MAX_SEQ/8 + 1)*sizeof(unsigned char));
    memset(nodes, 0, STT_NOD*sizeof(struct _node));
    memset(genes, 0, MAX_GENES*sizeof(struct _gene));
    memset(&tinf, 0, sizeof(struct _training));

    nn = 0; slen = 0; ipath = 0; ng = 0; nmask = 0; fng = 0;
    user_tt = 0; is_meta = 0; num_seq = 0; quiet = 0;
    max_phase = 0; max_score = -100.0;
    train_file = NULL;
    start_file = NULL; trans_file = NULL; nuc_file = NULL;
    start_ptr = stdout; trans_ptr = stdout; nuc_ptr = stdout;
    input_file = NULL; output_file = NULL;
    max_slen = 0;
    output = 0; closed = 1; do_mask = 0; force_nonsd = 0;
    is_first_meta = 1;

    tinf.st_wt = 4.35;
    tinf.trans_table = 11;

    meta = (struct _metagenomic_bin *)malloc(NUM_META * sizeof(struct _metagenomic_bin));
    for(size_t i = 0; i < NUM_META; i++){
        meta[i].tinf = new _training();
    }
}

void ProdigalWrapper::trainASpecies(unsigned char * genome, size_t seqLength) {
    // Initialize memories to reuse them
    memset(seq, 0, (slen / 4 + 1) * sizeof(unsigned char));
    memset(rseq, 0, (slen / 4 + 1) * sizeof(unsigned char));
    memset(useq, 0, (slen / 8 + 1) * sizeof(unsigned char));
    memset(nodes, 0, nn * sizeof(struct _node));
    nn = 0; slen = 0; ipath = 0; nmask = 0;

    // Initialize training information
    memset(mlist, 0, MAX_MASKS*sizeof(mask));
    memset(&tinf, 0, sizeof(struct _training));
    tinf.st_wt = 4.35;
    tinf.trans_table = 11;

    slen = getNextSeq(genome, 1, seqLength);

    rcom_seq(seq, rseq, useq, slen);

    /***********************************************************************
      Find all the potential starts and stops, sort them, and create a
      comprehensive list of nodes for dynamic programming.
    ***********************************************************************/
    if(slen > max_slen && slen > STT_NOD*8) {
        nodes = (struct _node *)realloc(nodes, (int)(slen/8)*sizeof(struct _node));
        if(nodes == NULL) {
            fprintf(stderr, "Realloc failed on nodes\n\n");
            exit(11);
        }
        memset(nodes, 0, (int)(slen/8)*sizeof(struct _node));
        max_slen = slen;
    }
    nn = add_nodes(seq, rseq, slen, nodes, closed, mlist, nmask, &tinf);
    qsort(nodes, nn, sizeof(struct _node), &compare_nodes);

    /***********************************************************************
      Scan all the ORFS looking for a potential GC bias in a particular
      codon position.  This information will be used to acquire a good
      initial set of genes.
    ***********************************************************************/
    gc_frame = calc_most_gc_frame(seq, slen);
//    if(gc_frame == NULL) {
//        fprintf(stderr, "Malloc failed on gc frame plot\n\n");
//        exit(11);
//    }
    record_gc_bias(gc_frame, nodes, nn, &tinf);
    free(gc_frame);

    /***********************************************************************
      Do an initial dynamic programming routine with just the GC frame
      bias used as a scoring function.  This will get an initial set of
      genes to train on.
    ***********************************************************************/
    record_overlapping_starts(nodes, nn, &tinf, 0);

    ipath = dprog(nodes, nn, &tinf, 0);

    /***********************************************************************
      Gather dicodon statistics for the training set.  Score the entire set
      of nodes.
    ***********************************************************************/
    calc_dicodon_gene(&tinf, seq, rseq, slen, nodes, ipath);
    raw_coding_score(seq, rseq, slen, nodes, nn, &tinf);

    /***********************************************************************
      Determine if this organism uses Shine-Dalgarno or not and score the
      nodes appropriately.
    ***********************************************************************/

    rbs_score(seq, rseq, slen, nodes, nn, &tinf);
    train_starts_sd(seq, rseq, slen, nodes, nn, &tinf);
    determine_sd_usage(&tinf);
    if(force_nonsd == 1) tinf.uses_sd = 0;
    if(tinf.uses_sd == 0) train_starts_nonsd(seq, rseq, slen, nodes, nn, &tinf); 
}

void ProdigalWrapper::trainMeta(unsigned char *genome, size_t seqLength) {
    // Initialize memories to reuse them
    memset(seq, 0, (slen / 4 + 1) * sizeof(unsigned char));
    memset(rseq, 0, (slen / 4 + 1) * sizeof(unsigned char));
    memset(useq, 0, (slen / 8 + 1) * sizeof(unsigned char));
    memset(nodes, 0, nn * sizeof(struct _node));
    nn = 0; slen = 0; ipath = 0; nmask = 0;

    // Initialize training information
    memset(&tinf, 0, sizeof(struct _training));
    tinf.st_wt = 4.35;
    tinf.trans_table = 11;
    
    if (is_first_meta == 1) {
        initialize_metagenomic_bins(meta);
        is_first_meta = 0;
    }
    
    slen = getNextSeq(genome, 1, seqLength);

    rcom_seq(seq, rseq, useq, slen);

    if(slen > max_slen && slen > STT_NOD*8) {
        nodes = (struct _node *)realloc(nodes, (int)(slen/8)*sizeof(struct _node));
        if(nodes == NULL) {
            fprintf(stderr, "Realloc failed on nodes\n\n");
            exit(11);
        }
        memset(nodes, 0, (int)(slen/8)*sizeof(struct _node));
        max_slen = slen;
    }

    low = 0.88495*tinf.gc - 0.0102337;
    if(low > 0.65) low = 0.65;
    high = 0.86596*tinf.gc + .1131991;
    if(high < 0.35) high = 0.35;
    max_score = -100.0;
    for(int i = 0; i < NUM_META; i++) {
        if (i == 0 || meta[i].tinf->trans_table !=
                      meta[i - 1].tinf->trans_table) {
            memset(nodes, 0, nn * sizeof(struct _node));
            nn = add_nodes(seq, rseq, slen, nodes, closed, mlist, nmask,
                           meta[i].tinf);
            qsort(nodes, nn, sizeof(struct _node), &compare_nodes);
        }
        if (meta[i].tinf->gc < low || meta[i].tinf->gc > high) continue;
        reset_node_scores(nodes, nn);
        score_nodes(seq, rseq, slen, nodes, nn, meta[i].tinf, closed, is_meta);
        record_overlapping_starts(nodes, nn, meta[i].tinf, 1);
        ipath = dprog(nodes, nn, meta[i].tinf, 1);
        if(ipath == -1) continue;

        if (nodes[ipath].score > max_score) {
            max_phase = i;
            max_score = nodes[ipath].score;
        }
    }
}

void ProdigalWrapper::getPredictedGenes(unsigned char * genome, size_t seqLength) {
    // Initialize memories to reuse them
    // Initialization should be done here not at the end of the function
    memset(seq, 0, (slen / 4 + 1) * sizeof(unsigned char));
    memset(rseq, 0, (slen / 4 + 1) * sizeof(unsigned char));
    memset(useq, 0, (slen / 8 + 1) * sizeof(unsigned char));
    memset(nodes, 0, nn*sizeof(struct _node));
    nn = 0; slen = 0; nmask = 0; ipath=0;

    /* Initialize structure */
    slen = getNextSeq(genome, 0, seqLength);
    rcom_seq(seq, rseq, useq, slen);
    if(slen == 0) {
        fprintf(stderr, "\nSequence read failed (file must be Fasta, ");
        fprintf(stderr, "Genbank, or EMBL format).\n\n");
        exit(14);
    }

    /* Reallocate memory if this is the biggest sequence we've seen */
    if(slen > max_slen && slen > STT_NOD*8) {
        nodes = (struct _node *)realloc(nodes, (int)(slen/8)*sizeof(struct _node));
        if(nodes == NULL) {
            fprintf(stderr, "Realloc failed on nodes\n\n");
            exit(11);
        }
        memset(nodes, 0, (int)(slen/8)*sizeof(struct _node));
        max_slen = slen;
    }

    if(is_meta == 0) {
        /***********************************************************************
         Find all the potential starts and stops, sort them, and create
         comprehensive list of nodes for dynamic programming.
         ***********************************************************************/
        nn = add_nodes(seq, rseq, slen, nodes, closed, mlist, nmask, &tinf);
        qsort(nodes, nn, sizeof(struct _node), &compare_nodes);

        /***********************************************************************
        Second dynamic programming, using the dicodon statistics as the
        scoring function.
        ***********************************************************************/
        score_nodes(seq, rseq, slen, nodes, nn, &tinf, closed, is_meta);
        record_overlapping_starts(nodes, nn, &tinf, 1);

        ipath = dprog(nodes, nn, &tinf, 1);
        eliminate_bad_genes(nodes, ipath, &tinf);
        ng = add_genes(genes, nodes, ipath);

        tweak_final_starts(genes, ng, nodes, nn, &tinf);
        record_gene_data(genes, ng, nodes, &tinf, num_seq);
    }
    else{

        /// Metagenomic version
    
        nn = add_nodes(seq, rseq, slen, nodes, closed, mlist, nmask,
                       meta[max_phase].tinf);
        qsort(nodes, nn, sizeof(struct _node), &compare_nodes);
        score_nodes(seq, rseq, slen, nodes, nn, meta[max_phase].tinf, closed,
                    is_meta);
        record_overlapping_starts(nodes, nn, meta[max_phase].tinf, 1);
        ipath = dprog(nodes, nn, meta[max_phase].tinf, 1);
        eliminate_bad_genes(nodes, ipath, meta[max_phase].tinf);
        ng = add_genes(genes, nodes, ipath);
        tweak_final_starts(genes, ng, nodes, nn, meta[max_phase].tinf);
        record_gene_data(genes, ng, nodes, meta[max_phase].tinf, num_seq);
    }
}

int ProdigalWrapper::getNextSeq(unsigned char * line, int training, size_t length) {
    int bctr = 0, len = 0;
    int gc_cont = 0, mask_beg = -1;
    // size_t lengthOfLine = strlen(line);

    for(size_t i = 0; i < length; i++) {
        if(line[i] < 'A' || line[i] > 'z') {
            continue;
        }
        if(do_mask == 1 && mask_beg != -1 && line[i] != 'N' && line[i] != 'n') {
            if(len - mask_beg >= MASK_SIZE) {
                if(nmask == MAX_MASKS) {
                    fprintf(stderr, "Error: saw too many regions of 'N''s in the ");
                    fprintf(stderr, "sequence.\n");
                    exit(52);
                }
                mlist[nmask].begin = mask_beg;
                mlist[nmask].end = len-1;
                (nmask)++;
            }
            mask_beg = -1;
        }
        if(do_mask == 1 && mask_beg == -1 && (line[i] == 'N' || line[i] == 'n'))
            mask_beg = len;
        if(line[i] == 'g' || line[i] == 'G'){
            set(seq, bctr); gc_cont++;
        }else if(line[i] == 't' || line[i] == 'T') {
            set(seq, bctr);
            set(seq, bctr+1);
        }else if(line[i] == 'c' || line[i] == 'C') {
            set(seq, bctr+1);
            gc_cont++;
        }else if(line[i] != 'a' && line[i] != 'A') {
            set(seq, bctr+1);
            set(useq, len);
        }

        bctr+=2; len++;

        if(len >= MAX_SEQ) {
            // fprintf(stderr, "\n\nWarning:  Sequence is long (max %d).\n", MAX_SEQ);
            // fprintf(stderr, "Use the first %d bases.\n\n", MAX_SEQ);
            break;
        }
    }

    if(training == 1){
        tinf.gc = ((double)gc_cont / (double)len);
    }else{
        gc = ((double)gc_cont / (double)len);
    }

    return len;
}

int ProdigalWrapper::getNumberOfPredictedGenes(){ return ng; }

void ProdigalWrapper::printGenes() {
    for(int i = 0 ; i < ng; i++){
        cout<<i<<" "<<genes[i].begin<<" "<<genes[i].end<<" "<<nodes[genes[i].start_ndx].strand<<endl;
    }
}

void ProdigalWrapper::removeCompletelyOverlappingGenes() {
    fng = 0;
    if(ng == 0) return;
    for(int i = 0; i < ng - 1; i++){
        if(genes[i].begin >= genes[i+1].begin) {
            continue;
        }
        finalGenes[fng++] = genes[i];
    }
    finalGenes[fng++] = genes[ng-1];
}

_training * ProdigalWrapper::getTrainingInfo() { return & tinf; }

void ProdigalWrapper::setTrainingInfo(_training &tinf) {
    this->tinf = tinf;
}

void ProdigalWrapper::updateTrainingInfo(_training &tinf2) {
    this->tinf = tinf2;
}

// It makes the blocks for translation
// Each block has a predicted gene part and an intergenic region. When another gene shows up, new block starts.
void ProdigalWrapper::getExtendedORFs(struct _gene *genes, struct _node *nodes, vector<SequenceBlock> &blocks,
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
                    leftEnd = genes[geneIdx - 1].end -1 -22;
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

    // // For the last gene
    // // Extend to the end of the genome
    // isReverse = !(nodes[genes[numOfGene - 1].start_ndx].strand == 1);    
    // rightEnd = length - 1;
    // if (isReverse) {
    //     frame = (genes[numOfGene - 1].end - 1) % 3;
    //     while (rightEnd % 3 != frame) rightEnd--;
    // }
    // // If left region is not covered, cover it.
    // leftEnd = genes[numOfGene - 1].begin - 1;
    // if (hasBeenExtendedToLeft) {
    //     leftEnd = genes[numOfGene - 2].end - 1 - 22;
    //     if (!isReverse) {
    //         frame = (genes[numOfGene - 1].begin - 1) % 3;
    //         while (leftEnd % 3 != frame) leftEnd++;
    //     }
    // }
    // blocks.emplace_back(leftEnd, rightEnd, isReverse ? -1 : 1);
    // if (find(intergenicKmerList.begin(), intergenicKmerList.end(), rightKmerHash) == intergenicKmerList.end()) {
    //         intergenicKmerList.push_back(rightKmerHash);
    // }

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
