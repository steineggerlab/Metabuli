#ifndef ADCLASSIFIER2_LOCALCOMMANDDECLARATIONS_H
#define ADCLASSIFIER2_LOCALCOMMANDDECLARATIONS_H
#include "Command.h"

extern int build(int argc, const char **argv, const Command& command);
extern int updateDB(int argc, const char **argv, const Command& command);
extern int classify(int argc, const char **argv, const Command& command);
extern int groupGeneration(int argc, const char **argv, const Command& command);
extern int filter(int argc, const char **argv, const Command& command);
extern int grade(int argc, const char **argv, const Command& command);
extern int gradeByCladeSize(int argc, const char **argv, const Command& command);
extern int seqHeader2TaxId(int argc, const char **argv, const Command& command);
extern int addToLibrary(int argc, const char **argv, const Command& command);
extern int filterByGenus(int argc, const char **argv, const Command& command);
extern int databaseReport(int argc, const char **argv, const Command& command);
extern int mapping2taxon(int argc, const char **argv, const Command& command);
extern int expand_diffidx(int argc, const char **argv, const Command& command);
extern int makeAAoffset(int argc, const char **argv, const Command& command);
extern int extract(int argc, const char **argv, const Command& command);
extern int printInfo(int argc, const char **argv, const Command& command);
extern int query2reference(int argc, const char **argv, const Command& command);
extern int ictvFormat(int argc, const char **argv, const Command& command);
extern int taxdump(int argc, const char **argv, const Command& command);
extern int accession2taxid(int argc, const char **argv, const Command& command);
extern int editNames(int argc, const char **argv, const Command& command);
extern int createnewtaxalist(int argc, const char **argv, const Command& command);
extern int classifiedRefiner(int argc, const char **argv, const Command& command);
extern int validateDatabase(int argc, const char **argv, const Command& command);
extern int printDeltaIdx(int argc, const char **argv, const Command& command); 
extern int makeBenchmarkSet(int argc, const char **argv, const Command &command);
extern int makeQuerySet(int argc, const char **argv, const Command &command);
extern int makeVirusBenchmarkSet(int argc, const char **argv, const Command &command);
extern int count_common_kmers(int argc, const char **argv, const Command &command);
extern int create_common_kmer_list(int argc, const char **argv, const Command &command);
extern int create_unique_kmer_list(int argc, const char **argv, const Command &command);
extern int create_unirefdb(int argc, const char **argv, const Command &command);
extern int create_uniref_tree(int argc, const char **argv, const Command &command);
extern int assign_uniref(int argc, const char **argv, const Command &command);

#endif //ADCLASSIFIER2_LOCALCOMMANDDECLARATIONS_H
