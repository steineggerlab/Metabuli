#include "LocalParameters.h"
#include "Parameters.h"
#include "Debug.h"
#include "CommandCaller.h"
#include "Parameters.cpp"
#include "ByteParser.h"
#include <iomanip>
#include "DistanceCalculator.h"

extern const char *version;

LocalParameters::LocalParameters() :
        Parameters(),
        VIRUS_TAX_ID(VIRUS_TAX_ID_ID,
                     "--virus-taxid",
                     "Taxonomy ID of virus taxon",
                     "NCBI: 10239 \nCUSTOM: Check names.dmp file ",
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
        SKIP_REDUNDANCY(SKIP_REDUNDANCY_ID,
                        "--skip-redundancy",
                        "Not storing k-mer's redundancy.",
                        "Not storing k-mer's redundancy.",
                        typeid(int),
                        (void *) &skipRedundancy,
                        "[0-1]"),
        SYNCMER(SYNCMER_ID,
                "--syncmer",
                "Using syncmer for k-mer selection",
                "Using syncmer for k-mer selection",
                typeid(int),
                (void *) &syncmer,
                "[0-1]"),
        SMER_LEN(SMER_LEN_ID,
                 "--smer-len",
                 "s-mer length for syncmer selection",
                 "s-mer length for syncmer selection",
                 typeid(int),
                 (void *) &smerLen,
                 "[0-6]"),
        SEQ_MODE(SEQ_MODE_ID,
                 "--seq-mode",
                 "Sequencing type",
                 "Single-end: 1, Paired-end: 2, Long read: 3",
                 typeid(int),
                 (void *) &seqMode,
                 "[1-3]"),
        REDUCED_AA(REDUCED_AA_ID,
                   "--reduced-aa",
                   "Using 15 alphabets to encode AAs for sensitivity",
                   "Set as 0 to use 15 alphabets to encode AAs for sensitivity",
                   typeid(int),
                   (void *) &reducedAA,
                   "[0-1]"),
        MIN_SCORE(MIN_SCORE_ID,
                  "--min-score",
                  "Min. sequence similarity score",
                  "Min. sequence similarity score (0.0-1.0)",
                  typeid(float),
                  (void *) &minScore,
                  "^0(\\.[0-9]+)?|1(\\.0+)?$"),
        MIN_COVERAGE(MIN_COVERAGE_ID,
                     "--min-cov",
                     "Min. query coverage",
                     "Min. query coverage (0.0-1.0)",
                     typeid(float),
                     (void *) &minCoverage,
                     "^0(\\.[0-9]+)?|1(\\.0+)?$"),
        // SPACED(SPACED_ID,
        //        "--spacing-mask",
        //        "Binary patterned mask for spaced k-mer.\nThe same mask must be used for DB creation and classification",
        //        "Binary patterned mask for spaced k-mer. The same mask must be used for DB creation and classification.\n"
        //        "A mask should contain at least eight '1's, and '0' means skip.",
        //        typeid(std::string),
        //        (void *) &spaceMask,
        //        "^.*$"),
        MIN_COVERED_POS(MIN_COVERED_POS_ID,
                        "--min-covered-pos",
                        "Minimum number of covered positions of a range",
                        "Minimum number of covered positions of a range",
                        typeid(int),
                        (void *) &minCoveredPos,
                        "^[0-9]+$"),
        HAMMING_MARGIN(HAMMING_MARGIN_ID,
                       "--hamming-margin",
                       "Allowed extra Hamming distance", 
                       "It allows extra Hamming distance than the minimum distance.",
                       typeid(int),
                       (void *) &hammingMargin,
                       ""),
        MIN_SP_SCORE(MIN_SP_SCORE_ID,
                     "--min-sp-score",
                     "Min. score for species- or lower-level classification.",
                     "Min. score for species- or lower-level classification.",
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
                     "Min. num. of cons. matches for non-euk. classification",
                     "Min. number of consecutive matches for prokaryote/virus classification",
                     typeid(int),
                     (void *) &minConsCnt,
                     "^[0-9]+$"),
        MIN_CONS_CNT_EUK(MIN_CONS_CNT_EUK_ID,
                         "--min-cons-cnt-euk",
                         "Min. num. of cons. matches for euk. classification",
                         "Min. number of consecutive matches for eukaryote classification",
                         typeid(int),
                         (void *) &minConsCntEuk,
                         "^[0-9]+$"),
        MATCH_PER_KMER(MATCH_PER_KMER_ID,
                       "--match-per-kmer",
                       "Number of matches per query k-mer. ",
                       "Num. of matches per query k-mer. Larger values assign more memory for storing k-mer matches. ",
                       typeid(int),
                       (void *) &matchPerKmer,
                       "^[0-9]+$"),
        MIN_SS_MATCH(MIN_SS_MATCH_ID,
                    "--min-ss-match",
                    "Min. num. of ssp.-specific matches for ssp. classification",
                    "Min. number of ssp.-specific matches for ssp. classification",
                    typeid(int),
                    (void *) &minSSMatch,
                    "^[0-9]+$"),
        TIE_RATIO(TIE_RATIO_ID,
                  "--tie-ratio",
                  "Best * --tie-ratio is considered as a tie",
                  "Best * --tie-ratio is considered as a tie",
                  typeid(float),
                  (void *) &tieRatio,
                  "^0(\\.[0-9]+)?|1(\\.0+)?$"),
        PRINT_LINEAGE(PRINT_LINEAGE_ID,
                      "--lineage",
                      "Print lineage information",
                      "Print lineage information",
                      typeid(int),
                      (void *) &printLineage,
                      "[0-1]"),
        MAX_SHIFT(MAX_SHIFT_ID,
                    "--max-shift",
                    "Max shift between two consecutive k-mers",
                    "Max shift between two consecutive k-mers",
                    typeid(int),
                    (void *) &maxShift,
                    "[1-7]"),
        TARGET_TAX_ID(TARGET_TAX_ID_ID,
               "--tax-id",
               "Tax. ID of clade to be extracted",
               "Tax. ID of clade to be extracted",
               typeid(int),
               (void *) &targetTaxId,
               "^[0-9]+$"),
        EXTRACT_MODE(EXTRACT_MODE_ID,
                     "--extract-format",
                     "0: original format, 1: FASTA, 2: FASTQ",
                     "0: original format, 1: FASTA, 2: FASTQ",
                     typeid(int),
                     (void *) &extractMode,
                     "[0-2]"),
        THR_K(THR_K_ID,
                    "--thr-k",
                    "Min. num. of shared kmer for read grouping",
                    "Min. num. of shared kmer for read grouping",
                    typeid(float),
                    (void *) &thresholdK,
                    "^(-?(10(\\.0+)?|[0-9](\\.[0-9]+)?))$"),
        GROUP_SCORE_THR(GROUP_SCORE_THR_ID,
                    "--group-score-thr",
                    "Min. score for read grouping",
                    "Min. score for read grouping",
                    typeid(float),
                    (void *) &groupScoreThr,
                    "^0(\\.[0-9]+)?|1(\\.0+)?$"),
        VOTE_MODE(VOTE_MODE_ID,
                    "--vote-mode",
                    "Vote mode of majority weighted LCA",
                    "Vote mode of majority weighted LCA",
                    typeid(int),
                    (void *) &voteMode,
                    "^(0|1|2)$"),
        MAJORITY_THR(MAJORITY_THR_ID,
                    "--majority-thr",
                    "Threshold for majority weighted LCA",
                    "Threshold for majority weighted LCA",
                    typeid(float),
                    (void *) &majorityThr,
                    "^0(\\.[0-9]+)?|1(\\.0+)?$"),
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
        ACCESSION_LEVEL(ACCESSION_LEVEL_ID,
                        "--accession-level",
                        "Accession-level DB build/search",
                        "Build or search a database for accession-level classification",
                        typeid(int),
                        (void *) &accessionLevel,
                        "[0-1]"),
        DB_NAME(DB_NAME_ID,
                "--db-name",
                "Name of the database (a random number as default)",
                "Name of the database",
                typeid(std::string),
                (void *) &dbName,
                "^.*$"),
        DB_DATE(DB_DATE_ID,
                "--db-date",
                "Date of the database creation (current date as default)",
                "Date of the database creation",
                typeid(std::string),
                (void *) &dbDate,
                "^.*$"),
        CDS_INFO(CDS_INFO_ID,
                 "--cds-info",
                 "List of CDS files",
                 "List of CDS files",
                 typeid(std::string),
                 (void *) &cdsInfo,
                 "^.*$"),
        MAKE_LIBRARY(MAKE_LIBRARY_ID,
                     "--make-library",
                     "Make library",
                     "Make library",
                     typeid(int),
                     (void *) &makeLibrary,
                     "[0-1]"),
        GTDB(GTDB_ID,
                "--gtdb",
                "GTDB-based database creation",
                "GTDB-based database creation",
                typeid(int),
                (void *) &gtdb,
                "[0-1]"),
        NEW_TAXA(NEW_TAXA_ID,
                "--new-taxa",
                "TSV file of new taxa to be added",
                "TSV file of new taxa to be added",
                typeid(std::string),
                (void *) &newTaxa,
                "^.*$"),
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
        CLADE_RANK(CLADE_RANK_ID,
                     "--clade-rank",
                     "Rank of clade to be tested",
                     "Rank of clade to be tested",
                     typeid(std::string),
                     (void *) &cladeRank,
                     "^.*$"),
        SKIP_SECONDARY(SKIP_SECONDARY_ID,
                          "--skip-secondary",
                          "Skip the results of already observed reads. (0: No, 1: Yes)",
                          "Skip secondary classification",
                          typeid(int),
                          (void *) &skipSecondary,
                          "[0-1]"),
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
                        "^.*$"),
        INFO_BEGIN(INFO_BEGIN_ID,
                "--info-begin",
                "Begin of the info to print",
                "Begin of the info to print",
                typeid(size_t),
                (void *) &infoBegin,
                "^[0-9]+$"), 
        INFO_END(INFO_END_ID,
                "--info-end",
                "End of the info to print",
                "End of the info to print",
                typeid(size_t),
                (void *) &infoEnd,
                "^[0-9]+$"),
        KMER_BEGIN(KMER_BEGIN_ID,
                "--kmer-begin",
                "First k-mer to print",
                "First k-mer to print",
                typeid(size_t),
                (void *) &kmerBegin,
                "^[0-9]+$"),
        KMER_END(KMER_END_ID,
                "--kmer-end",
                "Last k-mer to print",
                "Last k-mer to print",
                typeid(size_t),
                (void *) &kmerEnd,
                "^[0-9]+$")
  {
    // Initialize the parameters
    // Superkingdom taxonomy id
    virusTaxId = 10239;
    bacteriaTaxId = 2;
    archaeaTaxId = 2157;


    // Classify
    seqMode = 2;
    reducedAA = 0;
    minScore = 0;
    minConsCnt = 4;
    hammingMargin = 0;
    minSpScore = 0;
    minCoverage = 0;
    ramUsage = 0;
    minCoveredPos = 0;
    printLog = 0;
    maxGap = 0;
    minConsCntEuk = 0;
    matchPerKmer = 0;
    minSSMatch = 0;
    tieRatio = 0;

    // Group generation
    thresholdK = 0.5;
    voteMode = 0;
    majorityThr = 0.0;

    // Database creation
    tinfoPath = "";
    libraryPath = "";
    taxonomyPath = "";
    dbName = "";
    dbDate = "";
    splitNum = 0;
    bufferSize = 0;
    accessionLevel = 0;

    // Test parameters
    testRank = "";
    testType = "";
    printColumns = "";
    readIdCol = 0;
    taxidCol = 0;
    scoreCol = 0;
    coverageCol = 0;
    cladeRank = "";
    skipSecondary = 0;

    // Add to library
    assembly = false;

    // Filter
    printMode = 0;
    contamList = "";


    // build
    build.push_back(&PARAM_THREADS);
    build.push_back(&TAXONOMY_PATH);
    build.push_back(&SPLIT_NUM);
    build.push_back(&PARAM_MASK_PROBABILTY);
    build.push_back(&PARAM_MASK_RESIDUES);
    build.push_back(&ACCESSION_LEVEL);
    build.push_back(&DB_NAME);
    build.push_back(&DB_DATE);
    build.push_back(&CDS_INFO);
    build.push_back(&RAM_USAGE);
    build.push_back(&MAKE_LIBRARY);
    build.push_back(&GTDB);
    build.push_back(&SYNCMER);
    build.push_back(&SMER_LEN);

    // updateDB
    updateDB.push_back(&PARAM_THREADS);
    updateDB.push_back(&SPLIT_NUM);
    updateDB.push_back(&PARAM_MASK_PROBABILTY);
    updateDB.push_back(&PARAM_MASK_RESIDUES);
    updateDB.push_back(&ACCESSION_LEVEL);
    updateDB.push_back(&DB_NAME);
    updateDB.push_back(&DB_DATE);
    updateDB.push_back(&CDS_INFO);
    updateDB.push_back(&RAM_USAGE);
    updateDB.push_back(&NEW_TAXA);
    updateDB.push_back(&MAKE_LIBRARY);
    updateDB.push_back(&SYNCMER);

    //classify
    classify.push_back(&PARAM_THREADS);
    classify.push_back(&SEQ_MODE);
    classify.push_back(&MIN_SCORE);
    classify.push_back(&MIN_COVERAGE);
    classify.push_back(&MIN_CONS_CNT);
    classify.push_back(&MIN_CONS_CNT_EUK);
    classify.push_back(&MIN_SP_SCORE);
    classify.push_back(&HAMMING_MARGIN);
    classify.push_back(&TAXONOMY_PATH);
    classify.push_back(&PARAM_MASK_RESIDUES);
    classify.push_back(&PARAM_MASK_PROBABILTY);
    classify.push_back(&RAM_USAGE);
    classify.push_back(&MATCH_PER_KMER);
    classify.push_back(&ACCESSION_LEVEL);
    classify.push_back(&TIE_RATIO);
    classify.push_back(&SKIP_REDUNDANCY);
    classify.push_back(&PRINT_LINEAGE);
    classify.push_back(&SYNCMER);
    classify.push_back(&SMER_LEN);
    classify.push_back(&PRINT_LOG);
    classify.push_back(&MAX_SHIFT);


    // extract
    extract.push_back(&TAXONOMY_PATH);
    extract.push_back(&SEQ_MODE);
    extract.push_back(&TARGET_TAX_ID);
    extract.push_back(&EXTRACT_MODE);

    //groupGeneration
    groupGeneration.push_back(&PARAM_THREADS);
    groupGeneration.push_back(&SEQ_MODE);
    groupGeneration.push_back(&TAXONOMY_PATH);
    groupGeneration.push_back(&RAM_USAGE);
    groupGeneration.push_back(&MATCH_PER_KMER);
    groupGeneration.push_back(&THR_K);
    groupGeneration.push_back(&GROUP_SCORE_THR);
    groupGeneration.push_back(&VOTE_MODE);
    groupGeneration.push_back(&MAJORITY_THR);

    
    groupGeneration.push_back(&MIN_SCORE);
    groupGeneration.push_back(&MIN_COVERAGE);
    groupGeneration.push_back(&MIN_CONS_CNT);
    groupGeneration.push_back(&MIN_CONS_CNT_EUK);
    groupGeneration.push_back(&MIN_SP_SCORE);
    groupGeneration.push_back(&HAMMING_MARGIN);    
    groupGeneration.push_back(&PARAM_MASK_RESIDUES);
    groupGeneration.push_back(&PARAM_MASK_PROBABILTY);
    groupGeneration.push_back(&ACCESSION_LEVEL);
    groupGeneration.push_back(&TIE_RATIO);

    // filter 
    filter.push_back(&PARAM_THREADS);
    filter.push_back(&SEQ_MODE);
    filter.push_back(&VIRUS_TAX_ID);
    filter.push_back(&REDUCED_AA);
    filter.push_back(&MIN_SCORE);
    filter.push_back(&MIN_COVERAGE);
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
    filter.push_back(&ACCESSION_LEVEL);
    
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
    grade.push_back(&CLADE_RANK);
    grade.push_back(&SKIP_SECONDARY);

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

    // printInfo
    printInfo.push_back(&INFO_BEGIN);
    printInfo.push_back(&INFO_END);

    // expand_diffidx
    expand_diffidx.push_back(&KMER_BEGIN);
    expand_diffidx.push_back(&KMER_END);

    query2reference.push_back(&TEST_RANK);
}

void LocalParameters::printParameters(const std::string &module, int argc, const char* pargv[],
                                 const std::vector<MMseqsParameter*> &par){
    if (Debug::debugLevel < Debug::INFO) {
        return;
    }

    Debug(Debug::INFO) << module << " ";
    for (int i = 0; i < argc; i++) {
        // don't expose users to the interal b64 masking of whitespace characters
        if (strncmp("b64:", pargv[i], 4) == 0) {
            Debug(Debug::INFO) << "'" << base64_decode(pargv[i] + 4, strlen(pargv[i]) - 4) << "' ";
        } else {
            Debug(Debug::INFO) << pargv[i] << " ";
        }
    }
    Debug(Debug::INFO) << "\n\n";

    if (CommandCaller::getCallDepth() > 0) {
        return;
    }

    size_t maxWidth = 0;
    for(size_t i = 0; i < par.size(); i++) {
        maxWidth = std::max(strlen(par[i]->display), maxWidth);
    }

    std::stringstream ss;
    ss << std::boolalpha;

    ss << std::setw(maxWidth) << std::left  << "Metabuli Version (commit):" << "\t" << version << "\n";

    for (size_t i = 0; i < par.size(); i++) {
        if (par[i]->category & MMseqsParameter::COMMAND_HIDDEN) {
            continue;
        }
        ss << std::setw(maxWidth) << std::left << par[i]->display << "\t";
        if (typeid(int) == par[i]->type ) {
            ss << *((int *)par[i]->value);
        } else if(typeid(size_t) == par[i]->type ){
            ss << *((size_t *)par[i]->value);
        } else if(typeid(ByteParser) == par[i]->type) {
            ss << ByteParser::format(*((size_t *)par[i]->value), 'a', 'h');
        } else if(PARAM_SUB_MAT.uniqid == par[i]->uniqid || PARAM_SEED_SUB_MAT.uniqid == par[i]->uniqid) {
            MultiParam<NuclAA<std::string>>* param = ((MultiParam<NuclAA<std::string>>*) par[i]->value);
            MultiParam<NuclAA<std::string>> tmpPar(NuclAA<std::string>(
                BaseMatrix::unserializeName(param->values.aminoacid().c_str()),
                BaseMatrix::unserializeName(param->values.nucleotide().c_str())
            ));
            ss << MultiParam<NuclAA<std::string>>::format(tmpPar);
        } else if(typeid(MultiParam<NuclAA<std::string>>) == par[i]->type) {
            ss << MultiParam<NuclAA<std::string>>::format(*((MultiParam<NuclAA<std::string>> *)par[i]->value));
        } else if(typeid(MultiParam<NuclAA<int>>) == par[i]->type) {
            ss << MultiParam<NuclAA<int>>::format(*((MultiParam<NuclAA<int>> *)par[i]->value));
        } else if(typeid(MultiParam<NuclAA<float>>) == par[i]->type) {
            ss << MultiParam<NuclAA<float>>::format(*((MultiParam<NuclAA<float>> *)par[i]->value));
        } else if(typeid(MultiParam<SeqProf<int>>) == par[i]->type) {
            ss << MultiParam<SeqProf<int>>::format(*((MultiParam<SeqProf<int>> *)par[i]->value));
        } else if(typeid(MultiParam<PseudoCounts>) == par[i]->type) {
            ss << MultiParam<PseudoCounts>::format(*((MultiParam<PseudoCounts> *)par[i]->value));
        } else if(typeid(float) == par[i]->type) {
            ss << *((float *)par[i]->value);
        } else if(typeid(double) == par[i]->type) {
            ss << *((double *)par[i]->value);
        } else if(typeid(std::string) == par[i]->type) {
            ss << *((std::string *) par[i]->value);
        } else if (typeid(bool) == par[i]->type) {
            ss << *((bool *)par[i]->value);
        }
        ss << "\n";
    }

    Debug(Debug::INFO) << ss.str() << "\n";
}

void LocalParameters::parseParameters(int argc, const char *pargv[], const Command &command, bool printPar, int parseFlags,
                                 int outputFlags) {
    filenames.clear();
    std::vector<MMseqsParameter*> & par = *command.params;

    bool canHandleHelp = false;
    for (size_t parIdx = 0; parIdx < par.size(); parIdx++) {
        if (par[parIdx]->uniqid == PARAM_HELP_ID || par[parIdx]->uniqid == PARAM_HELP_LONG_ID) {
            canHandleHelp = true;
        }
    }

    size_t parametersFound = 0;
    for (int argIdx = 0; argIdx < argc; argIdx++) {
        // it is a parameter if it starts with - or --
        const bool longParameter = (pargv[argIdx][0] == '-' && pargv[argIdx][1] == '-');
        if (longParameter || (pargv[argIdx][0] == '-')) {
            if ((parseFlags & PARSE_REST) && longParameter && pargv[argIdx][2] == '\0') {
                restArgv = pargv + argIdx + 1;
                restArgc = argc - (argIdx + 1);
                break;
            }
            std::string parameter(pargv[argIdx]);
            if (canHandleHelp == false && (parameter.compare("-h") == 0 || parameter.compare("--help") == 0)) {
                printUsageMessage(command, 0xFFFFFFFF);
                EXIT(EXIT_SUCCESS);
            }

            bool hasUnrecognizedParameter = true;
            for (size_t parIdx = 0; parIdx < par.size(); parIdx++) {
                if(parameter.compare(par[parIdx]->name) == 0) {
                    if (typeid(bool) != par[parIdx]->type && argIdx + 1 == argc) {
                        printUsageMessage(command, outputFlags);
                        Debug(Debug::ERROR) << "Missing argument " << par[parIdx]->name << "\n";
                        EXIT(EXIT_FAILURE);
                    }

                    if (par[parIdx]->wasSet) {
                        printUsageMessage(command, outputFlags);
                        Debug(Debug::ERROR) << "Duplicate parameter " << par[parIdx]->name << "\n";
                        EXIT(EXIT_FAILURE);
                    }

                    if (typeid(int) == par[parIdx]->type) {
                        regex_t regex;
                        compileRegex(&regex, par[parIdx]->regex);
                        int nomatch = regexec(&regex, pargv[argIdx+1], 0, NULL, 0);
                        regfree(&regex);
                        // if no match found or two matches found (we want exactly one match)
                        if (nomatch){
                            printUsageMessage(command, 0xFFFFFFFF);
                            Debug(Debug::ERROR) << "Error in argument " << par[parIdx]->name << "\n";
                            EXIT(EXIT_FAILURE);
                        }else{
                            *((int *) par[parIdx]->value) = atoi(pargv[argIdx+1]);
                            par[parIdx]->wasSet = true;
                        }
                        argIdx++;
                    } else if (typeid(size_t) == par[parIdx]->type) {
                        regex_t regex;
                        compileRegex(&regex, par[parIdx]->regex);
                        int nomatch = regexec(&regex, pargv[argIdx+1], 0, NULL, 0);
                        regfree(&regex);
                        // if no match found or two matches found (we want exactly one match)
                        if (nomatch){
                            printUsageMessage(command, 0xFFFFFFFF);
                            Debug(Debug::ERROR) << "Error in argument " << par[parIdx]->name << "\n";
                            EXIT(EXIT_FAILURE);
                        }else{
                            *((size_t *) par[parIdx]->value) = strtoull(pargv[argIdx+1], nullptr, 10); //atoi(pargv[argIdx+1]);
                            // std::cout << "Value: " << *((size_t *) par[parIdx]->value) << std::endl;
                            par[parIdx]->wasSet = true;
                        }
                        argIdx++;
                    } else if (typeid(ByteParser) == par[parIdx]->type) {
                        regex_t regex;
                        compileRegex(&regex, par[parIdx]->regex);
                        int nomatch = regexec(&regex, pargv[argIdx+1], 0, NULL, 0);
                        regfree(&regex);

                        // if no match found or two matches found (we want exactly one match)
                        if (nomatch){
                            printUsageMessage(command, 0xFFFFFFFF);
                            Debug(Debug::ERROR) << "Error in argument regex " << par[parIdx]->name << "\n";
                            EXIT(EXIT_FAILURE);
                        } else {
                            size_t value = ByteParser::parse(pargv[argIdx+1]);
                            if (value == ByteParser::INVALID_SIZE) {
                                printUsageMessage(command, 0xFFFFFFFF);
                                Debug(Debug::ERROR) << "Error in value parsing " << par[parIdx]->name << "\n";
                                EXIT(EXIT_FAILURE);
                            } else {
                                *((size_t *) par[parIdx]->value) = value;
                                par[parIdx]->wasSet = true;
                            }
                        }
                        argIdx++;
                    } else if (typeid(MultiParam<NuclAA<std::string>>) == par[parIdx]->type) {
                        std::string val(pargv[argIdx+1]);
                        if (Util::startWith("b64:", val)) {
                            val = base64_decode(val.c_str() + 4, val.size() - 4);
                        }
                        NuclAA<std::string> value = MultiParam<NuclAA<std::string>>(val.c_str()).values;
                        if (value.first == "INVALID" || value.second == "INVALID") {
                            printUsageMessage(command, 0xFFFFFFFF);
                            Debug(Debug::ERROR) << "Error in value parsing " << par[parIdx]->name << "\n";
                            EXIT(EXIT_FAILURE);
                        } else {
                            *((MultiParam<NuclAA<std::string>> *) par[parIdx]->value) = value;
                            par[parIdx]->wasSet = true;
                        }
                        argIdx++;
                    }else if (typeid(MultiParam<NuclAA<int>>) == par[parIdx]->type) {
                        NuclAA<int> value = MultiParam<NuclAA<int>>(pargv[argIdx+1]).values;
                        if (value.first == INT_MAX || value.second == INT_MAX) {
                            printUsageMessage(command, 0xFFFFFFFF);
                            Debug(Debug::ERROR) << "Error in value parsing " << par[parIdx]->name << "\n";
                            EXIT(EXIT_FAILURE);
                        } else {
                            *((MultiParam<NuclAA<int>> *) par[parIdx]->value) = value;
                            par[parIdx]->wasSet = true;
                        }
                        argIdx++;
                    }else if (typeid(MultiParam<NuclAA<float>>) == par[parIdx]->type) {
                        NuclAA<float> value = MultiParam<NuclAA<float>>(pargv[argIdx + 1]).values;
                        if (value.first == FLT_MAX || value.second == FLT_MAX) {
                            printUsageMessage(command, 0xFFFFFFFF);
                            Debug(Debug::ERROR) << "Error in value parsing " << par[parIdx]->name << "\n";
                            EXIT(EXIT_FAILURE);
                        } else {
                            *((MultiParam<NuclAA<float>> *) par[parIdx]->value) = value;
                            par[parIdx]->wasSet = true;
                        }
                        argIdx++;
                    }else if (typeid(MultiParam<SeqProf<int>>) == par[parIdx]->type) {
                        SeqProf<int> value = MultiParam<SeqProf<int>>(pargv[argIdx+1]).values;
                        *((MultiParam<SeqProf<int>> *) par[parIdx]->value) = value;
                        par[parIdx]->wasSet = true;
                        argIdx++;
                    }else if (typeid(MultiParam<PseudoCounts>) == par[parIdx]->type) {
                        PseudoCounts value = MultiParam<PseudoCounts>(pargv[argIdx + 1]).values;
                        if (value.first == FLT_MAX || value.second == FLT_MAX) {
                            printUsageMessage(command, 0xFFFFFFFF);
                            Debug(Debug::ERROR) << "Error in value parsing " << par[parIdx]->name << "\n";
                            EXIT(EXIT_FAILURE);
                        } else {
                            *((MultiParam<PseudoCounts> *) par[parIdx]->value) = value;
                            par[parIdx]->wasSet = true;
                        }
                        argIdx++;
                    }else if (typeid(float) == par[parIdx]->type) {
                        regex_t regex;
                        compileRegex(&regex, par[parIdx]->regex);
                        int nomatch = regexec(&regex, pargv[argIdx+1], 0, NULL, 0);
                        regfree(&regex);
                        if (nomatch){
                            printUsageMessage(command, 0xFFFFFFFF);
                            Debug(Debug::ERROR) << "Error in argument " << par[parIdx]->name << "\n";
                            EXIT(EXIT_FAILURE);
                        }else{
                            double input = strtod(pargv[argIdx+1], NULL);
                            *((float *) par[parIdx]->value) = static_cast<float>(input);
                            par[parIdx]->wasSet = true;
                        }
                        argIdx++;
                    } else if (typeid(double) == par[parIdx]->type) {
                        regex_t regex;
                        compileRegex(&regex, par[parIdx]->regex);
                        int nomatch = regexec(&regex, pargv[argIdx+1], 0, NULL, 0);
                        regfree(&regex);
                        if (nomatch){
                            printUsageMessage(command, 0xFFFFFFFF);
                            Debug(Debug::ERROR) << "Error in argument " << par[parIdx]->name << "\n";
                            EXIT(EXIT_FAILURE);
                        }else{
                            *((double *) par[parIdx]->value) = strtod(pargv[argIdx+1], NULL);
                            par[parIdx]->wasSet = true;
                        }
                        argIdx++;
                    } else if (typeid(std::string) == par[parIdx]->type) {
                        std::string val(pargv[argIdx+1]);
                        if (Util::startWith("b64:", val)) {
                            val = base64_decode(val.c_str() + 4, val.size() - 4);
                        }
                        std::string* currVal = (std::string*)par[parIdx]->value;
                        currVal->assign(val);
                        par[parIdx]->wasSet = true;
                        argIdx++;
                    } else if (typeid(bool) == par[parIdx]->type) {
                        bool *value = (bool *) par[parIdx]->value;
                        if (argIdx + 1 == argc || pargv[argIdx+1][0] == '-') {
                            *value = !*value;
                        } else {
                            *value = parseBool(pargv[argIdx+1]);
                            argIdx++;
                        }
                        par[parIdx]->wasSet = true;
                    } else {
                        Debug(Debug::ERROR) << "Wrong parameter type in parseParameters. Please inform the developers\n";
                        EXIT(EXIT_FAILURE);
                    }

                    hasUnrecognizedParameter = false;
                    continue;
                }
            }

            if (hasUnrecognizedParameter) {
                printUsageMessage(command, 0xFFFFFFFF);

                // Suggest some parameter that the user might have meant
                std::vector<MMseqsParameter *>::const_iterator index = par.end();
                int maxDistance = 0;
                for (std::vector<MMseqsParameter *>::const_iterator it = par.begin(); it != par.end(); ++it) {
                    int distance = DistanceCalculator::localLevenshteinDistance(parameter, (*it)->name);
                    if (distance > maxDistance) {
                        maxDistance = distance;
                        index = it;
                    }
                }

                Debug(Debug::ERROR) << "Unrecognized parameter \"" << parameter << "\"";
                if (index != par.end()) {
                    Debug(Debug::ERROR) << ". Did you mean \"" << (*index)->name << "\" (" << (*index)->display << ")?\n";
                } else {
                    Debug(Debug::ERROR) << "\n";
                }

                EXIT(EXIT_FAILURE);
            }

            parametersFound++;
        } else {
            // parameter is actually a filename
#ifdef __CYGWIN__
            // normalize windows paths to cygwin unix paths
            const char *path = pargv[argIdx];
            ssize_t size = cygwin_conv_path(CCP_WIN_A_TO_POSIX | CCP_RELATIVE, path, NULL, 0);
            if (size < 0) {
                Debug(Debug::ERROR) << "Could not convert cygwin path!\n";
                EXIT(EXIT_FAILURE);
            } else {
                char *posix = new char[size];
                if (cygwin_conv_path(CCP_WIN_A_TO_POSIX | CCP_RELATIVE, path, posix, size)) {
                    Debug(Debug::ERROR) << "Could not convert cygwin path!\n";
                    EXIT(EXIT_FAILURE);
                }
                filenames.emplace_back(posix);
                delete posix;
            }
#else
            filenames.emplace_back(pargv[argIdx]);
#endif
        }
    }

    if (MMseqsMPI::isMaster()) {
        Debug::setDebugLevel(verbosity);
    }

#ifdef OPENMP
    omp_set_num_threads(threads);
#endif
#ifndef OPENMP
    threads = 1;
#endif


    bool ignorePathCountChecks = command.databases.empty() == false && command.databases[0].specialType & DbType::ZERO_OR_ALL && filenames.size() == 0;
    const size_t MAX_DB_PARAMETER = 6;
    if (ignorePathCountChecks == false && command.databases.size() > MAX_DB_PARAMETER) {
        Debug(Debug::ERROR) << "Use argv if you need more than " << MAX_DB_PARAMETER << " db parameters" << "\n";
        EXIT(EXIT_FAILURE);
    }

    if (ignorePathCountChecks == false && filenames.size() < command.databases.size()){
        printUsageMessage(command, outputFlags);
        Debug(Debug::ERROR) << "Not enough input paths provided. ";
        if (command.databases.size() == 1) {
            Debug(Debug::ERROR) << "1 path is required.\n";
        } else {
            Debug(Debug::ERROR) << command.databases.size() << " paths are required.\n";
        }
        EXIT(EXIT_FAILURE);
    }

    bool isVar = false;
    bool isStartVar = false;
    bool isMiddleVar = false;
    bool isEndVar = false;
    if(command.databases.empty() == false && command.databases[0].validator != NULL) {
        if (command.databases.size() >= 2) {
            for(size_t i = 0; i < command.databases.size();i++){
                if(i == 0){
                    isStartVar |= (command.databases[i].specialType & DbType::VARIADIC);
                } else if(i == command.databases.size() - 1){
                    isEndVar |= (command.databases[i].specialType & DbType::VARIADIC);
                } else {
                    isMiddleVar |= (command.databases[i].specialType & DbType::VARIADIC);
                }

            }
            isVar = isStartVar | isMiddleVar | isEndVar;
        }
        if (ignorePathCountChecks == false && isVar == false && filenames.size() > command.databases.size()) {
            printUsageMessage(command, outputFlags);
            Debug(Debug::ERROR) << "Too many input paths provided. Only " << SSTR(command.databases.size()) << " are allowed\n";
            EXIT(EXIT_FAILURE);
        }
    }
    switch (std::min(filenames.size(), MAX_DB_PARAMETER)) {
        case 6:
            db6 = filenames[5];
            db6Index = db6;
            db6Index.append(".index");
            db6dbtype = db6;
            db6dbtype.append(".dbtype");
            hdr6 = db6;
            hdr6.append("_h");
            hdr6Index = hdr6;
            hdr6Index.append(".index");
            hdr6dbtype = hdr6;
            hdr6dbtype.append(".dbtype");
            // FALLTHROUGH
        case 5:
            db5 = filenames[4];
            db5Index = db5;
            db5Index.append(".index");
            db5dbtype = db5;
            db5dbtype.append(".dbtype");
            hdr5 = db5;
            hdr5.append("_h");
            hdr5Index = hdr5;
            hdr5Index.append(".index");
            hdr5dbtype = hdr5;
            hdr5dbtype.append(".dbtype");
            // FALLTHROUGH
        case 4:
            db4 = filenames[3];
            db4Index = db4;
            db4Index.append(".index");
            db4dbtype = db4;
            db4dbtype.append(".dbtype");
            hdr4 = db4;
            hdr4.append("_h");
            hdr4Index = hdr4;
            hdr4Index.append(".index");
            hdr4dbtype = hdr4;
            hdr4dbtype.append(".dbtype");
            // FALLTHROUGH
        case 3:
            db3 = filenames[2];
            db3Index = db3;
            db3Index.append(".index");
            db3dbtype = db3;
            db3dbtype.append(".dbtype");
            hdr3 = db3;
            hdr3.append("_h");
            hdr3Index = hdr3;
            hdr3Index.append(".index");
            hdr3dbtype = hdr3;
            hdr3dbtype.append(".dbtype");
            // FALLTHROUGH
        case 2:
            db2 = filenames[1];
            db2Index = db2;
            db2Index.append(".index");
            db2dbtype = db2;
            db2dbtype.append(".dbtype");
            hdr2 = db2;
            hdr2.append("_h");
            hdr2Index = hdr2;
            hdr2Index.append(".index");
            hdr2dbtype = hdr2;
            hdr2dbtype.append(".dbtype");
            // FALLTHROUGH
        case 1:
            db1 = filenames[0];
            db1Index = db1;
            db1Index.append(".index");
            db1dbtype = db1;
            db1dbtype.append(".dbtype");
            hdr1 = db1;
            hdr1.append("_h");
            hdr1Index = hdr1;
            hdr1Index.append(".index");
            hdr1dbtype = hdr1;
            hdr1dbtype.append(".dbtype");
            break;
        default:
            // Do not abort execution if we expect a variable amount of parameters
            if (parseFlags & PARSE_VARIADIC)
                break;
            // FALLTHROUGH
        case 0:
            if (parseFlags & PARSE_ALLOW_EMPTY)
                break;
            printUsageMessage(command, outputFlags);
            printParameters(command.cmd, argc, pargv, par);
            Debug(Debug::ERROR) << "Unrecognized parameters!" << "\n";
            EXIT(EXIT_FAILURE);
    }

    initMatrices();

    if (ignorePathCountChecks == false) {
        checkIfDatabaseIsValid(command, argc, pargv, isStartVar, isMiddleVar, isEndVar);
    }

    if (printPar == true) {
        printParameters(command.cmd, argc, pargv, par);
    }
}