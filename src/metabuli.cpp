#include "Command.h"
#include "LocalCommandDeclarations.h"
#include "LocalParameters.h"

const char* binary_name = "metabuli";
const char* tool_name = "metabuli";
// TODO Write one full sentence
const char* tool_introduction = "Taxonomical classifier using 8-mers which have information of both of amino acid and DNA sequences";
const char* main_author = "Jaebeom Kim <jbeom0731@gmail.com> ";
const char* show_extended_help = "1";
const char* show_bash_info = NULL;
bool hide_base_commands = true;

LocalParameters& localPar = LocalParameters::getLocalInstance();
std::vector<Command> commands = {
        {"build_dir", build_dir, &localPar.build_dir, COMMAND_MAIN,
         "Building DB from multiple assemblies in input directory.",
         NULL,
         "Jaebeom Kim <jbeom0731@gmail.com>",
         "<i:directory of FASTA files> <o:DB directory>",
         CITATION_SPACEPHARER,
         {{"lowest directory including FASTA files", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA | DbType::VARIADIC, &DbValidator::directory},
          {"Directory where the DB will be generated", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::empty}}},

          {"build_fasta", build_fasta, &localPar.build_fasta, COMMAND_MAIN,
                "Building DB from a single FASTA file. ",
                NULL,
                "Jaebeom Kim <jbeom0731@gmail.com>",
                "<i:FASTA file> <i:mapping> <DB directory>",
                CITATION_SPACEPHARER,
                {{"A FASTA file", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA | DbType::VARIADIC, &DbValidator::flatfile},
                        {"Mapping file (accession to tax ID)", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::flatfile},
                        {"Directory where the DB will be generated", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::empty}}},
        {"classify", classify, &localPar.classify, COMMAND_MAIN,
         "Assigning taxonomy label to query reads",
         NULL,
         "Jaebeom Kim <jbeom0731@gmail.com>",
         "<i:FASTA> <i:DB dir> <o:out dir> <job ID> ",
         CITATION_SPACEPHARER,
         {{"FASTA", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA | DbType::VARIADIC, &DbValidator::flatfile},
          {"DB dir", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA | DbType::VARIADIC, &DbValidator::directory},
          {"out dir", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::directory},
          {"job ID", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::flatfile}}},

        {"inclusiontest", inclusiontest, &localPar.classify, COMMAND_MAIN,
                "It extracts k-mers from query sequences, and compares them to the target database",
                NULL,
                "Jaebeom Kim <jbeom0731@gmail.com>",
                "<i:read-classification>",
                CITATION_SPACEPHARER,
                {{"read-classification", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::flatfile}}},
        {"inclusiontest_hiv", inclusiontest_hiv, &localPar.classify, COMMAND_MAIN,
                "It extracts k-mers from query sequences, and compares them to the target database",
                NULL,
                "Jaebeom Kim <jbeom0731@gmail.com>",
                "<i:read-classification> <i:mapping> <i:taxonomy>",
                CITATION_SPACEPHARER,
                {{"read-classification", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::flatfile},
                        {"Mapping file (accession to tax ID)", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::flatfile},
                        {"taxonomy dir", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA | DbType::VARIADIC, &DbValidator::directory}}},

        {"exclusiontest", exclusiontest, &localPar.classify, COMMAND_MAIN,
                "It extracts k-mers from query sequences, and compares them to the target database",
                NULL,
                "Jaebeom Kim <jbeom0731@gmail.com>",
                "<i:read-classification> <i:included genome list>",
                CITATION_SPACEPHARER,
                {{"read-classification", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::flatfile},
                        {"included genome list", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::flatfile}}},
        {"exclusiontest_hiv", exclusiontest_hiv, &localPar.classify, COMMAND_MAIN,
                "It extracts k-mers from query sequences, and compares them to the target database",
                NULL,
                "Jaebeom Kim <jbeom0731@gmail.com>",
                "<i:read-classification> <i:mapping>",
                CITATION_SPACEPHARER,
                {{"read-classification", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::flatfile},
                        {"Mapping file (accession to tax ID)", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::flatfile}}},
        {"seqHeader2TaxId", seqHeader2TaxId, &localPar.seqHeader2TaxId, COMMAND_MAIN,
                "It extracts k-mers from query sequences, and compares them to the target database",
                NULL,
                "Jaebeom Kim <jbeom0731@gmail.com>",
                "<i:read-classification> <i:mapping>",
                CITATION_SPACEPHARER,
                {{"read-classification", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::flatfile},
                        {"Mapping file (accession to tax ID)", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::flatfile}}}

};

