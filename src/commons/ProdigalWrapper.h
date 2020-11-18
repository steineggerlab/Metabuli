//
// Created by 김재범 on 2020/11/10.
//

#ifndef ADCLASSIFIER2_PRODIGALWRAPPER_H
#define ADCLASSIFIER2_PRODIGALWRAPPER_H

#include <sys/stat.h>
#include <unistd.h>
#include "bitmap.h"
#include "dprog.h"
#include "gene.h"
#include "fptr.h"
#include "metagenomic.h"
#include "node.h"
#include "prodigalsequence.h"
#include "training.h"

#define MIN_SINGLE_GENOME 20000
#define IDEAL_SINGLE_GENOME 100000
#define VERSION "2.6.3"
#define DATE "February, 2016"

class ProdigalWrapper{
private:
    int rv, slen, nn, ng, ipath, *gc_frame, do_training, output, max_phase;
    int closed, do_mask, nmask, force_nonsd, user_tt, is_meta, num_seq, quiet;
    int piped, max_slen, fnum;
    double max_score, gc, low, high;
    unsigned char *seq, *rseq, *useq;
    char *train_file, *start_file, *trans_file, *nuc_file;
    char *input_file, *output_file, input_copy[MAX_LINE];
    char cur_header[MAX_LINE], new_header[MAX_LINE], short_header[MAX_LINE];
    FILE *output_ptr, *start_ptr, *trans_ptr, *nuc_ptr;
    fptr input_ptr = NULL;
    struct stat fbuf;
    pid_t pid;


    struct _training tinf;
    struct _metagenomic_bin meta[NUM_META];
    mask mlist[MAX_MASKS];

    /* Allocate memory and initialize variables */


public:
    struct _node *nodes;
    struct _gene *genes;
    int getNumberOfPredictedGenes();
    void getPredictedFrames(char * genome);
    void trainASpecies(char * genome);
    int getNextSeq(char * seq);
    ProdigalWrapper();
    ~ProdigalWrapper();

};
#endif //ADCLASSIFIER2_PRODIGALWRAPPER_H
