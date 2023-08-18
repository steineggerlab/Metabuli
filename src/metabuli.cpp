#include "Command.h"
#include "LocalCommandDeclarations.h"
#include "LocalParameters.h"
#include "DownloadDatabase.h"
#include "CommandDeclarations.h"
#include "Prefiltering.h"

#include "metabulidatabases.sh.h"

const char* binary_name = "metabuli";
const char* tool_name = "metabuli";
const char* tool_introduction = "Metabuli is a taxonomical classifier that jointly analyzes amino acid and DNA sequences.";
const char* main_author = "Jaebeom Kim <jbeom0731@gmail.com> ";
const char* show_extended_help = "1";
const char* show_bash_info = nullptr;
bool hide_base_commands = true;
extern const char* MMSEQS_CURRENT_INDEX_VERSION;
const char* index_version_compatible = MMSEQS_CURRENT_INDEX_VERSION;
bool hide_base_downloads = true;
void (*validatorUpdate)(void) = 0;




LocalParameters& localPar = LocalParameters::getLocalInstance();
std::vector<Command> commands = {
        {"databases",            databases,            &localPar.databases,            COMMAND_DATABASE_CREATION,
                "List and download databases",
                nullptr,
                "Milot Mirdita <milot@mirdita.de>",
                "<name> <o:sequenceDB> <tmpDir>",
                CITATION_SPACEPHARER, {{"selection", 0, DbType::ZERO_OR_ALL, &DbValidator::empty },
                                          {"sequenceDB", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::sequenceDb },
                                          {"tmpDir",     DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::directory }}},
       {"build", build, &localPar.build, COMMAND_DATABASE_CREATION,
                "Build database based on the list of FASTA files.",
               nullptr,
                "Jaebeom Kim <jbeom0731@gmail.com>",
                "<DB dir> <FASTA list> <Accesssion2taxid>",
                CITATION_SPACEPHARER,
                {{"Directory where the DB will be generated", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::empty},
                 {"A list of FASTA files", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::flatfile},
                 {"Mapping file (accession to tax ID)", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::flatfile}}},
        {"database-report", databaseReport, &localPar.databaseReport, COMMAND_DATABASE_CREATION,
                "It generates a report of taxa in a database.",
                nullptr,
                "Jaebeom Kim <jbeom0731@gmail.com>",
                "<i: DBDIR> ",
                CITATION_SPACEPHARER,
                {{"Directory where the DB will be generated", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::directory}}},
       {"updateDB", build, &localPar.build, COMMAND_DB,
               "Update database based on the list of FASTA files.",
                nullptr,
               "Jaebeom Kim <jbeom0731@gmail.com>",
               "<DB dir> <FASTA list> <Accesssion2taxid>",
               CITATION_SPACEPHARER,
               {{"DB directory to be updated", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::empty},
                       {"A list of FASTA files", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::flatfile},
                       {"Mapping file (accession to tax ID)", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::flatfile}}},
        {"classify", classify, &localPar.classify, COMMAND_TAXONOMY,
         "Assigning taxonomy label to query reads",
                nullptr,
         "Jaebeom Kim <jbeom0731@gmail.com>",
         "<i:QUERY> <i:DB dir> <o:out dir> <job ID> ",
         CITATION_SPACEPHARER,
         {{"FASTA", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA | DbType::VARIADIC, &DbValidator::flatfile},
          {"DB dir", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::directory},
          {"out dir", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::directory},
          {"job ID", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::flatfile}}},
        {"filter", classify, &localPar.classify, COMMAND_TAXONOMY,
                "Filtering reads based on the classification result",
                nullptr,
                "Jaebeom Kim <jbeom0731@gmail.com>",
                "<i:READ FILE> <i:FILTER DB>",
                CITATION_SPACEPHARER,
                {{"READ FILE", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA | DbType::VARIADIC, &DbValidator::flatfile},
                        {"FILTER DB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::directory}}},
        {"grade", grade, &localPar.grade, COMMAND_EXPERT,
                "Grade the classification result (only for benchmarking)",
                nullptr,
                "Jaebeom Kim <jbeom0731@gmail.com>",
                "<i:Result list> <i:Answer sheet list> <i:Taxonomy Dir>",
                CITATION_SPACEPHARER,
                {{"read-classification", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::flatfile},
                        {"Mapping file (Query ID to tax ID)", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::flatfile},
                        {"taxonomy dir", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA | DbType::VARIADIC, &DbValidator::directory}}},
        {"seqHeader2TaxId", seqHeader2TaxId, &localPar.seqHeader2TaxId, COMMAND_EXPERT,
                "It extracts k-mers from query sequences, and compares them to the target database",
                nullptr,
                "Jaebeom Kim <jbeom0731@gmail.com>",
                "<i:read-classification> <i:mapping>",
                CITATION_SPACEPHARER,
                {{"read-classification", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::flatfile},
                        {"Mapping file (accession to tax ID)", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::flatfile}}},
        {"add-to-library", addToLibrary, &localPar.addToLibrary, COMMAND_DATABASE_CREATION,
                    "It bins sequences into distinct files according to their species referring their accession number.\n "
                    "It requires a mapping file (accession to tax ID) and NCBI style tax dump files in a taxonomy directory.",
                nullptr,
                    "Jaebeom Kim <jbeom0731@gmail.com>",
                    "<I: FASTA list> <I: accession2taxid> <I: DB DIR>",
                    CITATION_SPACEPHARER,
                    {{"List of absolute paths of files to be added. One path per line.", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::flatfile},
                     {"NCBI style accession2taxid file. It should be consistent to tax dump files.", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::flatfile},
                     {"DB directory", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::directory}}},
        {"apply-threshold", applyThreshold, &localPar.applyThreshold, COMMAND_MAIN,
                "Assigning taxonomy label to query reads",
                nullptr,
                "Jaebeom Kim <jbeom0731@gmail.com>",
                "<i:Old readclassification> <o:OUT DIR> <o:JOB ID> <i: TAXONOMY DIR> ",
                CITATION_SPACEPHARER,
                {{"Old Results", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::flatfile},
                        {"OUT DIR", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::directory},
                        {"JOB ID", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::flatfile},
                        {"TAXONOMY DIR", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::directory}}},
        {"binning2report", binning2report, &localPar.binning2report, COMMAND_FORMAT_CONVERSION,
                "It generates Kraken style report file from binning results",
                nullptr,
                "Jaebeom Kim <jbeom0731@gmail.com>",
                "<i:Binning Result> <o:OUT DIR> <o:JOB ID> <i: TAXONOMY DIR> ",
                CITATION_SPACEPHARER,
                {{"Binning Result", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::flatfile},
                        {"OUT DIR", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::directory},
                        {"JOB ID", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::flatfile},
                        {"TAXONOMY DIR", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::directory}}},
        {"filter-by-genus", filterByGenus, &localPar.filterByGenus, COMMAND_EXPERT,
                    "It filters out reads classified as a specific genus",
                nullptr,
                    "Jaebeom Kim <jbeom0731@gmail.com>",
                    "<i:Binning Result> <i:Genus list> <i: TAXONOMY DIR> ",
                    CITATION_SPACEPHARER,
                    {{"Binning Result", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::flatfile},
                     {"Genus list", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::flatfile},
                     {"TAXONOMY DIR", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::directory}}}

};
std::vector<KmerThreshold> externalThreshold = {};

std::vector<DatabaseDownload> externalDownloads = {
        {
                "RefSeq",
                "Database built with NCBI RefSeq assemblies (Complete/Chromosome level only, Prokaryote & Virus) and Human genome (GRCh38.p14)",
                "O'Leary et al. Reference sequence (RefSeq) database at NCBI: current status, taxonomic expansion, and functional annotation. Nucleic Acids Res. (2016)",
                "https://www.ncbi.nlm.nih.gov/refseq/",
                true, LocalParameters::DBTYPE_INDEX_DB, metabulidatabases_sh, metabulidatabases_sh_len,
                {}
        },
        {
                "RefSeq217",
                "Database built with genomes of NCBI release 217 (Prokaryote & Virus) and Human genome (GRCh38.p14)",
                "O'Leary et al. Reference sequence (RefSeq) database at NCBI: current status, taxonomic expansion, and functional annotation. Nucleic Acids Res. (2016)",
                "https://www.ncbi.nlm.nih.gov/refseq/",
                true, LocalParameters::DBTYPE_INDEX_DB, metabulidatabases_sh, metabulidatabases_sh_len,
                {}
        },
        {
                "GTDB207",
                "Database built with genomes and taxonomy of GTDB207 (Complete/Chromosome level only, CheckM Completeness > 90, CheckM Contamination < 5) and Human genome (GRCh38.p14)",
                "Donovan et al. GTDB: an ongoing census of bacterial and archaeal diversity through a phylogenetically consistent, rank normalized and complete genome-based taxonomy. Nucleic Acids Res. (2022)",
                "https://gtdb.ecogenomic.org/",
                true, LocalParameters::DBTYPE_INDEX_DB, metabulidatabases_sh, metabulidatabases_sh_len,
                {}
        },
        {
                "RefSeq_virus",
                "Database built with NCBI RefSeq virus assemblies (Complete/Chromosome level only)",
                "O'Leary et al. Reference sequence (RefSeq) database at NCBI: current status, taxonomic expansion, and functional annotation. Nucleic Acids Res. (2016)",
                "https://www.ncbi.nlm.nih.gov/refseq/",
                true, LocalParameters::DBTYPE_INDEX_DB, metabulidatabases_sh, metabulidatabases_sh_len,
                {}
        }
};


