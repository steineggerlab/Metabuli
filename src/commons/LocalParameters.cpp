#include "LocalParameters.h"
#include "Parameters.h"

LocalParameters::LocalParameters() :
        Parameters(),
        VIRUS_TAX_ID(VIRUS_TAX_ID_ID,
                     "--virus-taxid",
                     "Taxonomy ID of virus taxon",
                     "NCBI: 10239 [Default]\nCUSTOM: Check names.dmp file ",
                     typeid(int),
                     (void *) &virusTaxId,
                     "^[0-9]+$"),
        BACTERIA_TAX_ID(BACTERIA_TAX_ID_ID,
                        "--bacteria-taxid",
                        "Taxonomy ID of bacteria taxon",
                        "NCBI: 2 [Default]\nCUSTOM: Check names.dmp file ",
                        typeid(int),
                        (void *) &bacteriaTaxId,
                        "^[0-9]+$"),
        ARCHAEA_TAX_ID(ARCHAEA_TAX_ID_ID,
                       "--archaea-taxid",
                       "Taxonomy ID of archaea taxon",
                       "NCBI: 2157 [Default]\nCUSTOM: Check names.dmp file ",
                       typeid(int),
                       (void *) &archaeaTaxId,
                       "^[0-9]+$"),
        EUKARYOTA_TAX_ID(EUKARYOTA_TAX_ID_ID,
                         "--eukaryota-taxid",
                         "Taxonomy ID of eukaryota taxon",
                         "NCBI: 2759 [Default]\nCUSTOM: Check names.dmp file ",
                         typeid(int),
                         (void *) &eukaryotaTaxId,
                         "^[0-9]+$"),
        SEQ_MODE(SEQ_MODE_ID,
                 "--seq-mode",
                 "Sequencing type",
                 "Single-end: 1 \nPaired-end: 2\nLong read: 3",
                 typeid(int),
                 (void *) &seqMode,
                 "[1-3]"),
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
        MIN_COVERAGE(MIN_COVERAGE_ID,
                     "--min-cov",
                     "The minimum coverage for classification",
                     "You can set a value from 0.0 to 1.0",
                     typeid(float),
                     (void *) &minCoverage,
                     "^0(\\.[0-9]+)?|1(\\.0+)?$"),
        SPACED(SPACED_ID,
               "--spacing-mask",
               "Binary patterned mask for spaced k-mer.\nThe same mask must be used for DB creation and classification",
               "Binary patterned mask for spaced k-mer. The same mask must be used for DB creation and classification.\n"
               "A mask should contain at least eight '1's, and '0' means skip.",
               typeid(std::string),
               (void *) &spaceMask,
               "^.*$"),
        MIN_COVERED_POS(MIN_COVERED_POS_ID,
                        "--min-covered-pos",
                        "Minimum number of covered positions of a range",
                        "Minimum number of covered positions of a range",
                        typeid(int),
                        (void *) &minCoveredPos,
                        "^[0-9]+$"),
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
        TINFO_PATH(TINFO_PATH_ID,
                   "--tinfo-path",
                   "Path to prodigal training information files",
                   "Path to prodigal training information files",
                   typeid(std::string),
                   (void *) &tinfoPath,
                   "^.*$"),
        RAM_USAGE(RAM_USAGE_ID,
                  "--max-ram",
                  "RAM usage in GiB",
                  "RAM usage in GiB",
                  typeid(int),
                  (void *) &ramUsage,
                  "^[0-9]+$"),
        PRINT_LOG(PRINT_LOG_ID,
                  "--print-log",
                  "Print logs to debug",
                  "Print logs to debug",
                  typeid(int),
                  (void *) &printLog,
                  "^[0-9]+$"),
        MAX_GAP(MAX_GAP_ID,
                "--max-gap",
                "Maximum gap between two consecutive k-mers (used only with spaced k-mer)",
                "Maximum gap between two consecutive k-mers (used only with spaced k-mer)",
                typeid(int),
                (void *) &maxGap,
                "^[0-9]+$"),
        MIN_CONS_CNT(MIN_CONS_CNT_ID,
                     "--min-cons-cnt",
                     "Minimum number of consecutive metamer matches to be used for prokaryote/virus classification",
                     "Minimum number of consecutive metamer matches to be used for prokaryote/virus classification",
                     typeid(int),
                     (void *) &minConsCnt,
                     "^[0-9]+$"),
        MIN_CONS_CNT_EUK(MIN_CONS_CNT_EUK_ID,
                         "--min-cons-cnt-euk",
                         "Minimum number of consecutive metamer matches to be used for eukaryote classification",
                         "Minimum number of consecutive metamer matches to be used for eukaryote classification",
                         typeid(int),
                         (void *) &minConsCntEuk,
                         "^[0-9]+$"),
        MATCH_PER_KMER(MATCH_PER_KMER_ID,
                       "--match-per-kmer",
                       "Number of matches per query k-mer",
                       "Number of matches per query k-mer. Larger values assign more memory for storing k-mer matches.",
                       typeid(int),
                       (void *) &matchPerKmer,
                       "^[0-9]+$"),
        LIBRARY_PATH(LIBRARY_PATH_ID,
                     "--library-path",
                     "Path to library where the FASTA files are stored",
                     "Path to library where the FASTA files are stored",
                     typeid(std::string),
                     (void *) &libraryPath,
                     "^.*$"),
        TAXONOMY_PATH(TAXONOMY_PATH_ID,
                      "--taxonomy-path",
                      "Directory where the taxonomy dump files are stored",
                      "Directory where the taxonomy dump files are stored",
                      typeid(std::string),
                      (void *) &taxonomyPath,
                      "^.*$"),
        IS_ASSEMBLY(IS_ASSEMBLY_ID,
                    "--assembly",
                    "Input is an assembly",
                    "Input is an assembly",
                    typeid(bool),
                    (void *) &assembly,
                    ""),
        SPLIT_NUM(SPLIT_NUM_ID,
                  "--split-num",
                  "A database is divided to N splits (offsets). During classification, unnecessary splits are skipped",
                  "A database is divided to N splits (offsets). During classification, unnecessary splits are skipped",
                  typeid(int),
                  (void *) &splitNum,
                  "^[0-9]+$"),
        BUFFER_SIZE(BUFFER_SIZE_ID,
                    "--buffer-size",
                    "Buffer size (the number of k-mers)",
                    "Buffer size (the number of k-mers)",
                    typeid(size_t),
                    (void *) &bufferSize,
                    "^[0-9]+$"),
        TEST_RANK(TEST_RANK_ID,
                  "--test-rank",
                  ".",
                  "csv of ranks to be tested",
                  typeid(std::string),
                  (void *) &testRank,
                  "^.*$"),
        TEST_TYPE(TEST_TYPE_ID,
                  "--test-type",
                  ".",
                  "Test Type",
                  typeid(std::string),
                  (void *) &testType,
                  "^.*$"),
        READID_COL(READID_COL_ID,
                   "--readid-col",
                   "Column number of accession in classification result",
                   "Column number of accession in classification result",
                   typeid(int),
                   (void *) &readIdCol,
                   "^[0-9]+$"),
        TAXID_COL(TAXID_COL_ID,
                  "--taxid-col",
                  "Column number of taxonomy ID in classification result",
                  "Column number of taxonomy ID in classification result",
                  typeid(int),
                  (void *) &taxidCol,
                  "^[0-9]+$"),
        SCORE_COL(SCORE_COL_ID,
                  "--score-col",
                  "Column number of score in classification result",
                  "Column number of score in classification result",
                  typeid(int),
                  (void *) &scoreCol,
                  "^[0-9]+$"),
        COVERAGE_COL(COVERAGE_COL_ID,
                     "--coverage-col",
                     "Column number of coverage in classification result",
                     "Column number of coverage in classification result",
                     typeid(int),
                     (void *) &coverageCol,
                     "^[0-9]+$"),
        PRINT_COLUMNS(PRINT_COLUMNS_ID,
                      "--print-columns",
                      "CSV of column numbers to be printed",
                      "CSV of column numbers to be printed",
                      typeid(std::string),
                      (void *) &printColumns,
                      "^.*$"),
        PRINT_MODE(PRINT_MODE_ID,
                        "--print-mode",
                       "[1] Only filtered reads [2] Both filtered and removed reads",
                       "[1] Only filtered reads [2] Both filtered and removed reads",
                       typeid(int),
                       (void *) &printMode,
                       "[1-2]"),
        CONTAM_LIST(CONTAM_LIST_ID, 
                   "--contam-list",
                   "List of contaminants to be filtered",
                   "List of taxids to be filtered",
                     typeid(std::string),
                        (void *) &contamList,
                        "^.*$") 
  {
    //add_to_library

    // build
    build.push_back(&PARAM_THREADS);
    build.push_back(&REDUCED_AA);
    build.push_back(&SPACED);
    build.push_back(&TAXONOMY_PATH);
    build.push_back(&SPLIT_NUM);
    build.push_back(&PARAM_MASK_PROBABILTY);
    build.push_back(&PARAM_MASK_RESIDUES);
    build.push_back(&BUFFER_SIZE);

    //classify
    classify.push_back(&PARAM_THREADS);
    classify.push_back(&SEQ_MODE);
    classify.push_back(&VIRUS_TAX_ID);
    classify.push_back(&REDUCED_AA);
    classify.push_back(&MIN_SCORE);
    classify.push_back(&MIN_COVERAGE);
    classify.push_back(&SPACED);
    classify.push_back(&HAMMING_MARGIN);
    classify.push_back(&MIN_SP_SCORE);
    classify.push_back(&PARAM_V);
    classify.push_back(&RAM_USAGE);
    classify.push_back(&MIN_COVERED_POS);
    classify.push_back(&PRINT_LOG);
    classify.push_back(&MAX_GAP);
    classify.push_back(&TAXONOMY_PATH);
    classify.push_back(&MIN_CONS_CNT);
    classify.push_back(&MIN_CONS_CNT_EUK);
    classify.push_back(&PARAM_MASK_RESIDUES);
    classify.push_back(&PARAM_MASK_PROBABILTY);
    classify.push_back(&MATCH_PER_KMER);

    // filter 
    filter.push_back(&PARAM_THREADS);
    filter.push_back(&SEQ_MODE);
    filter.push_back(&VIRUS_TAX_ID);
    filter.push_back(&REDUCED_AA);
    filter.push_back(&MIN_SCORE);
    filter.push_back(&MIN_COVERAGE);
    filter.push_back(&SPACED);
    filter.push_back(&HAMMING_MARGIN);
    filter.push_back(&MIN_SP_SCORE);
    filter.push_back(&PARAM_V);
    filter.push_back(&RAM_USAGE);
    filter.push_back(&MIN_COVERED_POS);
    filter.push_back(&PRINT_LOG);
    filter.push_back(&MAX_GAP);
    filter.push_back(&TAXONOMY_PATH);
    filter.push_back(&MIN_CONS_CNT);
    filter.push_back(&MIN_CONS_CNT_EUK);
    filter.push_back(&PARAM_MASK_RESIDUES);
    filter.push_back(&PARAM_MASK_PROBABILTY);
    filter.push_back(&MATCH_PER_KMER);
    filter.push_back(&PRINT_MODE);
    filter.push_back(&CONTAM_LIST);

    //updateTargetDB
    exclusiontest_hiv.push_back(&TEST_RANK);

    // grade
    grade.push_back(&TEST_RANK);
    grade.push_back(&TEST_TYPE);
    grade.push_back(&PARAM_THREADS);
    grade.push_back(&READID_COL);
    grade.push_back(&TAXID_COL);
    grade.push_back(&PARAM_V);
    grade.push_back(&SCORE_COL);
    grade.push_back(&COVERAGE_COL);
    grade.push_back(&PRINT_COLUMNS);

    // Apply thresholds
    applyThreshold.push_back(&MIN_SP_SCORE);
    applyThreshold.push_back(&MIN_SCORE);
    applyThreshold.push_back(&MIN_COVERAGE);

    // Binning to report
    binning2report.push_back(&READID_COL);
    binning2report.push_back(&TAXID_COL);

    // add to library
    addToLibrary.push_back(&IS_ASSEMBLY);
    addToLibrary.push_back(&TAXONOMY_PATH);
    addToLibrary.push_back(&LIBRARY_PATH);

    // db report
    databaseReport.push_back(&TAXONOMY_PATH);
}

