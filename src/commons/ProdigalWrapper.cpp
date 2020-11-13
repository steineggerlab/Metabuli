//
// Created by 김재범 on 2020/11/10.
//
#include "ProdigalWrapper.h"

ProdigalWrapper::ProdigalWrapper() {
    seq = (unsigned char *)malloc(MAX_SEQ/4*sizeof(unsigned char));
    rseq = (unsigned char *)malloc(MAX_SEQ/4*sizeof(unsigned char));
    useq = (unsigned char *)malloc(MAX_SEQ/8*sizeof(unsigned char));
    nodes = (struct _node *)malloc(STT_NOD*sizeof(struct _node));
    genes = (struct _gene *)malloc(MAX_GENES*sizeof(struct _gene));
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
    train_file = NULL; do_training = 0;
    start_file = NULL; trans_file = NULL; nuc_file = NULL;
    start_ptr = stdout; trans_ptr = stdout; nuc_ptr = stdout;
    input_file = NULL; output_file = NULL; piped = 0;
    output_ptr = stdout; max_slen = 0;
    output = 0; closed = 0; do_mask = 0; force_nonsd = 0;

}
void ProdigalWrapper::trainASpecies(char * genome){
    slen = read_seq_training(input_ptr, seq, useq, &(tinf.gc), do_mask, mlist,
                             &nmask);
    if(slen == 0) {
        fprintf(stderr, "\n\nSequence read failed (file must be Fasta, ");
        fprintf(stderr, "Genbank, or EMBL format).\n\n");
        exit(9);
    }
    if(slen < MIN_SINGLE_GENOME) {
        fprintf(stderr, "\n\nError:  Sequence must be %d", MIN_SINGLE_GENOME);
        fprintf(stderr, " characters (only %d read).\n(Consider", slen);
        fprintf(stderr, " running with the -p meta option or finding");
        fprintf(stderr, " more contigs from the same genome.)\n\n");
        exit(10);
    }
    if(slen < IDEAL_SINGLE_GENOME) {
        fprintf(stderr, "\n\nWarning:  ideally Prodigal should be given at");
        fprintf(stderr, " least %d bases for ", IDEAL_SINGLE_GENOME);
        fprintf(stderr, "training.\nYou may get better results with the ");
        fprintf(stderr, "-p meta option.\n\n");
    }
    rcom_seq(seq, rseq, useq, slen);
    if(quiet == 0) {
        fprintf(stderr, "%d bp seq created, %.2f pct GC\n", slen, tinf.gc*100.0);
    }

    /***********************************************************************
      Find all the potential starts and stops, sort them, and create a
      comprehensive list of nodes for dynamic programming.
    ***********************************************************************/
    if(quiet == 0) {
        fprintf(stderr, "Locating all potential starts and stops...");
    }
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
    if(quiet == 0) {
        fprintf(stderr, "%d nodes\n", nn);
    }

    /***********************************************************************
      Scan all the ORFS looking for a potential GC bias in a particular
      codon position.  This information will be used to acquire a good
      initial set of genes.
    ***********************************************************************/
    if(quiet == 0) {
        fprintf(stderr, "Looking for GC bias in different frames...");
    }
    gc_frame = calc_most_gc_frame(seq, slen);
    if(gc_frame == NULL) {
        fprintf(stderr, "Malloc failed on gc frame plot\n\n");
        exit(11);
    }
    record_gc_bias(gc_frame, nodes, nn, &tinf);
    if(quiet == 0) {
        fprintf(stderr, "frame bias scores: %.2f %.2f %.2f\n", tinf.bias[0],
                tinf.bias[1], tinf.bias[2]);
    }
    free(gc_frame);

    /***********************************************************************
      Do an initial dynamic programming routine with just the GC frame
      bias used as a scoring function.  This will get an initial set of
      genes to train on.
    ***********************************************************************/
    if(quiet == 0) {
        fprintf(stderr, "Building initial set of genes to train from...");
    }
    record_overlapping_starts(nodes, nn, &tinf, 0);
    ipath = dprog(nodes, nn, &tinf, 0);
    if(quiet == 0) {
        fprintf(stderr, "done!\n");
    }

    /***********************************************************************
      Gather dicodon statistics for the training set.  Score the entire set
      of nodes.
    ***********************************************************************/
    if(quiet == 0) {
        fprintf(stderr, "Creating coding model and scoring nodes...");
    }
    calc_dicodon_gene(&tinf, seq, rseq, slen, nodes, ipath);
    raw_coding_score(seq, rseq, slen, nodes, nn, &tinf);
    if(quiet == 0) {
        fprintf(stderr, "done!\n");
    }

    /***********************************************************************
      Determine if this organism uses Shine-Dalgarno or not and score the
      nodes appropriately.
    ***********************************************************************/
    if(quiet == 0) {
        fprintf(stderr, "Examining upstream regions and training starts...");
    }
    rbs_score(seq, rseq, slen, nodes, nn, &tinf);
    train_starts_sd(seq, rseq, slen, nodes, nn, &tinf);
    determine_sd_usage(&tinf);
    if(force_nonsd == 1) tinf.uses_sd = 0;
    if(tinf.uses_sd == 0) train_starts_nonsd(seq, rseq, slen, nodes, nn, &tinf);
    if(quiet == 0) {
        fprintf(stderr, "done!\n");
    }

    memset(seq, 0, (slen/4+1)*sizeof(unsigned char));
    memset(rseq, 0, (slen/4+1)*sizeof(unsigned char));
    memset(useq, 0, (slen/8+1)*sizeof(unsigned char));
    memset(nodes, 0, nn*sizeof(struct _node));
    nn = 0; slen = 0; ipath = 0; nmask = 0;

}

void ProdigalWrapper::getPredictedFrames(char * genome){
    rcom_seq(seq, rseq, useq, slen);
    if(slen == 0) {
        fprintf(stderr, "\nSequence read failed (file must be Fasta, ");
        fprintf(stderr, "Genbank, or EMBL format).\n\n");
        exit(14);
    }

    if(quiet == 0) {
        fprintf(stderr, "Finding genes in sequence #%d (%d bp)...", num_seq, slen);
    }

    /* Reallocate memory if this is the biggest sequence we've seen */
    if(slen > max_slen && slen > STT_NOD*8) {
        nodes = (struct _node *)realloc(nodes, (int)(slen/8)*sizeof(struct _node));
        if(nodes == NULL) {
            fprintf(stderr, "Realloc failed on nodes\n\n");
            exit(11);
        }
        max_slen = slen;
    }

    /* Calculate short header for this sequence */
    calc_short_header(cur_header, short_header, num_seq);

    if(is_meta == 0) { /* Single Genome Version */

        /***********************************************************************
          Find all the potential starts and stops, sort them, and create a
          comprehensive list of nodes for dynamic programming.
        ***********************************************************************/
        nn = add_nodes(seq, rseq, slen, nodes, closed, mlist, nmask, &tinf);
        qsort(nodes, nn, sizeof(struct _node), &compare_nodes);

        /***********************************************************************
          Second dynamic programming, using the dicodon statistics as the
          scoring function.
        ***********************************************************************/
        score_nodes(seq, rseq, slen, nodes, nn, &tinf, closed, is_meta);
        if(start_ptr != stdout)
            write_start_file(start_ptr, nodes, nn, &tinf, num_seq, slen, 0, NULL,
                             VERSION, cur_header);
        record_overlapping_starts(nodes, nn, &tinf, 1);
        ipath = dprog(nodes, nn, &tinf, 1);
        eliminate_bad_genes(nodes, ipath, &tinf);
        ng = add_genes(genes, nodes, ipath);
        tweak_final_starts(genes, ng, nodes, nn, &tinf);
        record_gene_data(genes, ng, nodes, &tinf, num_seq);
        if(quiet == 0) {
            fprintf(stderr, "done!\n");
        }


    }
}

int ProdigalWrapper::getNextSeq(char * line) {
    int hdr = 0, fhdr = 0, bctr = 0, len = 0, wrn = 0;
    int gc_cont = 0, mask_beg = -1;
    unsigned int i, gapsize = 0;

    if(hdr == 0 && line[strlen(line)-1] != '\n' && wrn == 0) {
        wrn = 1;
        fprintf(stderr, "\n\nWarning: saw non-sequence line longer than ");
        fprintf(stderr, "%d chars, sequence might not be read ", MAX_LINE);
        fprintf(stderr, "correctly.\n\n");
    }
    if(line[0] == '>' || (line[0] == 'S' && line[1] == 'Q') ||
       (strlen(line) > 6 && strncmp(line, "ORIGIN", 6) == 0)) {
        hdr = 1;
        if(fhdr > 0) {
            for(i = 0; i < 12; i++) {
                if(i%4 == 0 || i%4 == 1) { set(seq, bctr); set(seq, bctr+1); }
                bctr+=2; len++;
            }
        }
        fhdr++;
    }
    else if(hdr == 1 && (line[0] == '/' && line[1] == '/')) hdr = 0;
    else if(hdr == 1) {
        if(strstr(line, "Expand") != NULL && strstr(line, "gap") != NULL) {
            sscanf(strstr(line, "gap")+4, "%u", &gapsize);
            if(gapsize < 1 || gapsize > MAX_LINE) {
                fprintf(stderr, "Error: gap size in gbk file can't exceed line");
                fprintf(stderr, " size.\n");
                exit(51);
            }
            for(i = 0; i < gapsize; i++) line[i] = 'n';
            line[i] = '\0';
        }
        for(i = 0; i < strlen(line); i++) {
            if(line[i] < 'A' || line[i] > 'z') continue;
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
            if(line[i] == 'g' || line[i] == 'G') { set(seq, bctr); gc_cont++; }
            else if(line[i] == 't' || line[i] == 'T') {
                set(seq, bctr);
                set(seq, bctr+1);
            }
            else if(line[i] == 'c' || line[i] == 'C') {
                set(seq, bctr+1);
                gc_cont++;
            }
            else if(line[i] != 'a' && line[i] != 'A') {
                set(seq, bctr+1);
                set(useq, len);
            }
            bctr+=2; len++;
        }
    }
    if(len+MAX_LINE >= MAX_SEQ) {
        fprintf(stderr, "\n\nWarning:  Sequence is long (max %d for training).\n",
                MAX_SEQ);
        fprintf(stderr, "Training on the first %d bases.\n\n", MAX_SEQ);
        break;
    }

    if(fhdr > 1) {
        for(i = 0; i < 12; i++) {
            if(i%4 == 0 || i%4 == 1) { set(seq, bctr); set(seq, bctr+1); }
            bctr+=2; len++;
        }
    }
    *gc = ((double)gc_cont / (double)len);
    return len;
}

int read_seq_training(fptr fp, unsigned char *seq, unsigned char *useq,
                      double *gc, int do_mask, mask *mlist, int *nm) {
    char line[MAX_LINE+1];
    int hdr = 0, fhdr = 0, bctr = 0, len = 0, wrn = 0;
    int gc_cont = 0, mask_beg = -1;
    unsigned int i, gapsize = 0;

    line[MAX_LINE] = '\0';
    while(INPUT_GETS(line, MAX_LINE, fp) != NULL) {
        if(hdr == 0 && line[strlen(line)-1] != '\n' && wrn == 0) {
            wrn = 1;
            fprintf(stderr, "\n\nWarning: saw non-sequence line longer than ");
            fprintf(stderr, "%d chars, sequence might not be read ", MAX_LINE);
            fprintf(stderr, "correctly.\n\n");
        }
        if(line[0] == '>' || (line[0] == 'S' && line[1] == 'Q') ||
           (strlen(line) > 6 && strncmp(line, "ORIGIN", 6) == 0)) {
            hdr = 1;
            if(fhdr > 0) {
                for(i = 0; i < 12; i++) {
                    if(i%4 == 0 || i%4 == 1) { set(seq, bctr); set(seq, bctr+1); }
                    bctr+=2; len++;
                }
            }
            fhdr++;
        }
        else if(hdr == 1 && (line[0] == '/' && line[1] == '/')) hdr = 0;
        else if(hdr == 1) {
            if(strstr(line, "Expand") != NULL && strstr(line, "gap") != NULL) {
                sscanf(strstr(line, "gap")+4, "%u", &gapsize);
                if(gapsize < 1 || gapsize > MAX_LINE) {
                    fprintf(stderr, "Error: gap size in gbk file can't exceed line");
                    fprintf(stderr, " size.\n");
                    exit(51);
                }
                for(i = 0; i < gapsize; i++) line[i] = 'n';
                line[i] = '\0';
            }
            for(i = 0; i < strlen(line); i++) {
                if(line[i] < 'A' || line[i] > 'z') continue;
                if(do_mask == 1 && mask_beg != -1 && line[i] != 'N' && line[i] != 'n') {
                    if(len - mask_beg >= MASK_SIZE) {
                        if(*nm == MAX_MASKS) {
                            fprintf(stderr, "Error: saw too many regions of 'N''s in the ");
                            fprintf(stderr, "sequence.\n");
                            exit(52);
                        }
                        mlist[*nm].begin = mask_beg;
                        mlist[*nm].end = len-1;
                        (*nm)++;
                    }
                    mask_beg = -1;
                }
                if(do_mask == 1 && mask_beg == -1 && (line[i] == 'N' || line[i] == 'n'))
                    mask_beg = len;
                if(line[i] == 'g' || line[i] == 'G') { set(seq, bctr); gc_cont++; }
                else if(line[i] == 't' || line[i] == 'T') {
                    set(seq, bctr);
                    set(seq, bctr+1);
                }
                else if(line[i] == 'c' || line[i] == 'C') {
                    set(seq, bctr+1);
                    gc_cont++;
                }
                else if(line[i] != 'a' && line[i] != 'A') {
                    set(seq, bctr+1);
                    set(useq, len);
                }
                bctr+=2; len++;
            }
        }
        if(len+MAX_LINE >= MAX_SEQ) {
            fprintf(stderr, "\n\nWarning:  Sequence is long (max %d for training).\n",
                    MAX_SEQ);
            fprintf(stderr, "Training on the first %d bases.\n\n", MAX_SEQ);
            break;
        }
    }
    if(fhdr > 1) {
        for(i = 0; i < 12; i++) {
            if(i%4 == 0 || i%4 == 1) { set(seq, bctr); set(seq, bctr+1); }
            bctr+=2; len++;
        }
    }
    *gc = ((double)gc_cont / (double)len);
    return len;
}
