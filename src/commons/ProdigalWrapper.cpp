//
// Created by 김재범 on 2020/11/10.
//
#include "ProdigalWrapper.h"
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
}
ProdigalWrapper::ProdigalWrapper() {
    seq = (unsigned char *)malloc(MAX_SEQ/4*sizeof(unsigned char)); // 8 Mb
    rseq = (unsigned char *)malloc(MAX_SEQ/4*sizeof(unsigned char)); // 8 Mb
    useq = (unsigned char *)malloc(MAX_SEQ/8*sizeof(unsigned char)); // 4 Mb
    nodes = (struct _node *)malloc(STT_NOD*sizeof(struct _node)); // 13.6 Mb
    genes = (struct _gene *)malloc(MAX_GENES*sizeof(struct _gene)); // 30 Mb
    if(seq == NULL || rseq == NULL || nodes == NULL || genes == NULL) {
        fprintf(stderr, "\nError: Malloc failed on sequence/orfs\n\n"); exit(1);
    }

    memset(seq, 0, MAX_SEQ/4*sizeof(unsigned char));
    memset(rseq, 0, MAX_SEQ/4*sizeof(unsigned char));
    memset(useq, 0, MAX_SEQ/8*sizeof(unsigned char));
    memset(nodes, 0, STT_NOD*sizeof(struct _node));
    memset(genes, 0, MAX_GENES*sizeof(struct _gene));
    memset(&tinf, 0, sizeof(struct _training));

    nn = 0; slen = 0; ipath = 0; ng = 0; nmask = 0;
    user_tt = 0; is_meta = 0; num_seq = 0; quiet = 0;
    max_phase = 0; max_score = -100.0;
    train_file = NULL;
    start_file = NULL; trans_file = NULL; nuc_file = NULL;
    start_ptr = stdout; trans_ptr = stdout; nuc_ptr = stdout;
    input_file = NULL; output_file = NULL;
    max_slen = 0;
    output = 0; closed = 1; do_mask = 0; force_nonsd = 0;

    tinf.st_wt = 4.35;
    tinf.trans_table = 11;

    meta = (struct _metagenomic_bin *)malloc(NUM_META * sizeof(struct _metagenomic_bin));
    for(size_t i = 0; i < NUM_META; i++){
        meta[i].tinf = new _training();
    }

}

void ProdigalWrapper::
trainASpecies(char * genome){

    memset(seq, 0, (slen/4+1)*sizeof(unsigned char));
    memset(rseq, 0, (slen/4+1)*sizeof(unsigned char));
    memset(useq, 0, (slen/8+1)*sizeof(unsigned char));
    memset(nodes, 0, nn*sizeof(struct _node));
    memset(&tinf, 0, sizeof(struct _training));
    nn = 0; slen = 0; ipath = 0; nmask = 0;
    tinf.st_wt = 4.35;
    tinf.trans_table = 11;

    fprintf(stderr, "Request:  Single Genome, Phase:  Training\n");
    fprintf(stderr, "Reading in the sequence(s) to train...\n");

    slen = getNextSeq(genome, 1);
    if(slen == 0) {
        fprintf(stderr, "\n\nSequence read failed (file must be Fasta, ");
        fprintf(stderr, "Genbank, or EMBL format).\n\n");
        exit(9);
    }
//    if(slen < MIN_SINGLE_GENOME) {
//        fprintf(stderr, "\n\nError:  Sequence must be %d", MIN_SINGLE_GENOME);
//        fprintf(stderr, " characters (only %d read).\n(Consider", slen);
//        fprintf(stderr, " running with the -p meta option or finding");
//        fprintf(stderr, " more contigs from the same genome.)\n\n");
//        exit(10);
//    }
    if(slen < IDEAL_SINGLE_GENOME) {
        fprintf(stderr, "\n\nWarning:  ideally Prodigal should be given at");
        fprintf(stderr, " least %d bases for ", IDEAL_SINGLE_GENOME);
        fprintf(stderr, "training.\nYou may get better results with the ");
        fprintf(stderr, "-p meta option.\n\n");
    }
    rcom_seq(seq, rseq, useq, slen);
//    if(quiet == 0) {
//        fprintf(stderr, "%d bp seq created, %.2f pct GC\n", slen, tinf.gc*100.0);
//    }

    /***********************************************************************
      Find all the potential starts and stops, sort them, and create a
      comprehensive list of nodes for dynamic programming.
    ***********************************************************************/
//    if(quiet == 0) {
//        fprintf(stderr, "Locating all potential starts and stops...");
//    }
    if(slen > max_slen && slen > STT_NOD*8) {
        nodes = (struct _node *)realloc(nodes, (int)(slen/8)*sizeof(struct _node));
        if(nodes == NULL) {
            fprintf(stderr, "Realloc failed on nodes\n\n");
            exit(11);
        }
        max_slen = slen;
    }
    nn = add_nodes(seq, rseq, slen, nodes, closed, mlist, nmask, &tinf);
    qsort(nodes, nn, sizeof(struct _node), &compare_nodes);
//    if(quiet == 0) {
//        fprintf(stderr, "%d nodes\n", nn);
//    }

    /***********************************************************************
      Scan all the ORFS looking for a potential GC bias in a particular
      codon position.  This information will be used to acquire a good
      initial set of genes.
    ***********************************************************************/
//    if(quiet == 0) {
//        fprintf(stderr, "Looking for GC bias in different frames...");
//    }
    gc_frame = calc_most_gc_frame(seq, slen);
    if(gc_frame == NULL) {
        fprintf(stderr, "Malloc failed on gc frame plot\n\n");
        exit(11);
    }
    record_gc_bias(gc_frame, nodes, nn, &tinf);
//    if(quiet == 0) {
//        fprintf(stderr, "frame bias scores: %.2f %.2f %.2f\n", tinf.bias[0],
//                tinf.bias[1], tinf.bias[2]);
//    }
    free(gc_frame);

    /***********************************************************************
      Do an initial dynamic programming routine with just the GC frame
      bias used as a scoring function.  This will get an initial set of
      genes to train on.
    ***********************************************************************/
//    if(quiet == 0) {
//        fprintf(stderr, "Building initial set of genes to train from...");
//    }
    record_overlapping_starts(nodes, nn, &tinf, 0);

    ipath = dprog(nodes, nn, &tinf, 0);
//    if(quiet == 0) {
//        fprintf(stderr, "done!\n");
//    }

    /***********************************************************************
      Gather dicodon statistics for the training set.  Score the entire set
      of nodes.
    ***********************************************************************/
//    if(quiet == 0) {
//        fprintf(stderr, "Creating coding model and scoring nodes...");
//    }
    calc_dicodon_gene(&tinf, seq, rseq, slen, nodes, ipath);
    raw_coding_score(seq, rseq, slen, nodes, nn, &tinf);
    if(quiet == 0) {
        fprintf(stderr, "done!\n");
    }

    /***********************************************************************
      Determine if this organism uses Shine-Dalgarno or not and score the
      nodes appropriately.
    ***********************************************************************/
//    if(quiet == 0) {
//        fprintf(stderr, "Examining upstream regions and training starts...");
//    }
    rbs_score(seq, rseq, slen, nodes, nn, &tinf);
    train_starts_sd(seq, rseq, slen, nodes, nn, &tinf);
    determine_sd_usage(&tinf);
    if(force_nonsd == 1) tinf.uses_sd = 0;
    if(tinf.uses_sd == 0) train_starts_nonsd(seq, rseq, slen, nodes, nn, &tinf);
//    if(quiet == 0) {
//        fprintf(stderr, "done!\n");
//    }
}

void ProdigalWrapper::trainMeta(char *genome) {

    memset(seq, 0, (slen/4+1)*sizeof(unsigned char));
    memset(rseq, 0, (slen/4+1)*sizeof(unsigned char));
    memset(useq, 0, (slen/8+1)*sizeof(unsigned char));
    memset(nodes, 0, nn*sizeof(struct _node));
    nn = 0; slen = 0; ipath = 0; nmask = 0;

    if(1) {
        fprintf(stderr, "Request:  Metagenomic, Phase:  Training\n");
        fprintf(stderr, "Initializing training files...");
    }

    initialize_metagenomic_bins(meta);

    if(1) {
        fprintf(stderr, "done!\n");
        fprintf(stderr, "-------------------------------------\n");
    }

    slen = getNextSeq(genome, 1);

    cout<<"gc"<<gc<<" "<<tinf.gc<<endl;

    rcom_seq(seq, rseq, useq, slen);
    if(slen == 0) {
        fprintf(stderr, "\nSequence read failed (file must be Fasta, ");
        fprintf(stderr, "Genbank, or EMBL format).\n\n");
        exit(14);
    }
    if(slen > max_slen && slen > STT_NOD*8) {
        nodes = (struct _node *)realloc(nodes, (int)(slen/8)*sizeof(struct _node));
        if(nodes == NULL) {
            fprintf(stderr, "Realloc failed on nodes\n\n");
            exit(11);
        }
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
        if (nodes[ipath].score > max_score) {
            max_phase = i;
            max_score = nodes[ipath].score;
//            eliminate_bad_genes(nodes, ipath, meta[i].tinf);
//            ng = add_genes(genes, nodes, ipath);
//            tweak_final_starts(genes, ng, nodes, nn, meta[i].tinf);
//            record_gene_data(genes, ng, nodes, meta[i].tinf, num_seq);
        }
    }
}
void ProdigalWrapper::getPredictedFrames(char * genome){

    memset(seq, 0, (slen/4+1)*sizeof(unsigned char));
    memset(rseq, 0, (slen/4+1)*sizeof(unsigned char));
    memset(useq, 0, (slen/8+1)*sizeof(unsigned char));
    memset(nodes, 0, nn*sizeof(struct _node));
    nn = 0; slen = 0; nmask = 0; ipath=0;

    /* Initialize structure */
    slen = getNextSeq(genome, 0);
    rcom_seq(seq, rseq, useq, slen);
    if(slen == 0) {
        fprintf(stderr, "\nSequence read failed (file must be Fasta, ");
        fprintf(stderr, "Genbank, or EMBL format).\n\n");
        exit(14);
    }

//    if(1) {
//        fprintf(stderr, "Finding genes in sequence #%d (%d bp)...\n", num_seq, slen);
//    }

    /* Reallocate memory if this is the biggest sequence we've seen */
    if(slen > max_slen && slen > STT_NOD*8) {
        nodes = (struct _node *)realloc(nodes, (int)(slen/8)*sizeof(struct _node));
        if(nodes == NULL) {
            fprintf(stderr, "Realloc failed on nodes\n\n");
            exit(11);
        }
        max_slen = slen;
    }

    if(is_meta == 0) {
        ipath = 0;
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

    ///metagenomic version
        //trainMeta(genome);
        fprintf(stderr, "Request:  Metagenomic, Phase:  Gene Finding\n");

        nn = add_nodes(seq, rseq, slen, nodes, closed, mlist, nmask,
                       meta[max_phase].tinf);
        qsort(nodes, nn, sizeof(struct _node), &compare_nodes);
        score_nodes(seq, rseq, slen, nodes, nn, meta[max_phase].tinf, closed,
                    is_meta);
        record_overlapping_starts(nodes, nn, meta[max_phase].tinf, 1);
        ipath = dprog(nodes, nn, meta[max_phase].tinf, 1);
        eliminate_bad_genes(nodes, ipath, meta[max_phase].tinf);
        ng = add_genes(genes, nodes, ipath);
        fprintf(stderr, "done! gene count: %d (%d bp)\n", ng, slen);
        tweak_final_starts(genes, ng, nodes, nn, meta[max_phase].tinf);
        record_gene_data(genes, ng, nodes, meta[max_phase].tinf, num_seq);
    }
//     if(1) {
//         fprintf(stderr, "done! gene count: %d (%d bp)\n", ng, slen);
//     }
}

int ProdigalWrapper::getNextSeq(char * line, int training) {
    int bctr = 0, len = 0;
    int gc_cont = 0, mask_beg = -1;
    size_t lengthOfLine = strlen(line);

    for(size_t i = 0; i < lengthOfLine; i++) {
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
            fprintf(stderr, "\n\nWarning:  Sequence is long (max %d for training).\n",
                    MAX_SEQ);
            fprintf(stderr, "Training on the first %d bases.\n\n", MAX_SEQ);
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




