#include "LocalParameters.h"

LocalParameters::LocalParameters() :
        Parameters(),
        VIRUS_TAX_ID(VIRUS_TAX_ID_ID,
                     "--virus-taxid",
                     "Taxonomy ID of virus taxon",
                     "NCBI: 10239 [Default]\nCUSTOM: Check names.dmp file ",
                     typeid(int),
                     (void *) &virusTaxId,
                     "[^[1-9]\\d*$]"),
        SEQ_MODE(SEQ_MODE_ID,
                 "--seq-mode",
                 "Sequencing type",
                 "Single-end: 1 \nPaired-end: 2\nLong read: 3",
                 typeid(int),
                 (void *) &seqMode,
                 "[1-3]"),
        MEMORY_MODE(MEMORY_MODE_ID,
                    "--memory-mode",
                    "Keeping k-mer matches in the RAM or writing into a file",
                    "Writing: 1 \nRAM:  2",
                    typeid(int),
                    (void *) &memoryMode,
                    "[1-2]"),
        REDUCED_AA(REDUCED_AA_ID,
                   "--reduced-aa",
                   "Using reduced 15 alphabets to encode amino acids. It increases sensitivity",
                   "Using 20 alphabets: 0 \nUsing 15 alphabets: 1",
                   typeid(int),
                   (void *) &reducedAA,
                   "[0-1]"),
        MIN_SCORE(MIN_SCORE_ID,
                  "--min-score",
                  "The minimum score for classification",
                  "You can set a value from 0.0 to 1.0",
                  typeid(float),
                  (void *) &minScore,
                  "^0(\\.[0-9]+)?|1(\\.0+)?$"),
        SPACED(SPACED_ID,
               "--spacing-mask",
               "Binary patterned mask for spaced k-mer.\nThe same mask must be used for DB creation and classification",
               "Binary patterned mask for spaced k-mer. The same mask must be used for DB creation and classification.\n"
               "A mask should contain at least eight '1's, and '0' means skip.",
               typeid(std::string),
               (void *) &spaceMask,
               ""),
        MIN_CONSECUTIVE(MIN_CONSECUTIVE_ID,
                        "--min-consecutive",
                        "Matches that are not consecutive for the specified number of times are ignored.",
                        "Matched k-mers from the same genus are pulled and aligned to query.\n"
                        "Matches that are not consecutive for the specified number of times are ignored.",
                        typeid(int),
                        (void *) &minConsCnt,
                        ""),
        HAMMING_MARGIN(HAMMING_MARGIN_ID,
                       "--hamming-margin",
                       "If a query k-mer has multiple matches, the matches with hamming distance lower than sum of \n"
                       "the minimum distance and this margin are selected for later steps.",
                       "If a query k-mer has multiple matches, the matches with hamming distance lower than sum of \n"
                       "the minimum distance and this margin are selected for later steps.",
                       typeid(int),
                       (void *) &hammingMargin,
                       ""),
        MIN_SP_SCORE(MIN_SP_SCORE_ID,
                       "--min-sp-score",
                       "Minimum score to be classified at species or lower rank.",
                       "Minimum score to be classified at the species level.",
                       typeid(float),
                       (void *) &minSpScore,
                     "^0(\\.[0-9]+)?|1(\\.0+)?$"),
        TEST_RANK(TEST_RANK_ID,
                  "--test-rank",
                  ".",
                  "Test Rank",
                  typeid(std::string),
                  (void *) &testRank,
                  ""),
        TEST_TYPE(TEST_TYPE_ID,
                  "--test-type",
                  ".",
                  "Test Type",
                  typeid(std::string),
                    (void *) &testType,
                    ""){
    //build_dir
    build_dir.push_back(&PARAM_THREADS);
    build_dir.push_back(&REDUCED_AA);
    build_dir.push_back(&SPACED);

    //build_fasta
    build_fasta.push_back(&PARAM_THREADS);
    build_fasta.push_back(&REDUCED_AA);
    build_fasta.push_back(&SPACED);

    //classify
    classify.push_back(&PARAM_THREADS);
    classify.push_back(&SEQ_MODE);
    classify.push_back(&VIRUS_TAX_ID);
    classify.push_back(&MEMORY_MODE);
    classify.push_back(&REDUCED_AA);
    classify.push_back(&MIN_SCORE);
    classify.push_back(&SPACED);
    classify.push_back(&MIN_CONSECUTIVE);
    classify.push_back(&HAMMING_MARGIN);
    classify.push_back(&MIN_SP_SCORE);

    //updateTargetDB
    exclusiontest_hiv.push_back(&TEST_RANK);
    inclusiontest.push_back(&TEST_TYPE);

}
