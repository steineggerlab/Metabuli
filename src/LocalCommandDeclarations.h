//
// Created by KJB on 25/09/2020.
//

#ifndef ADCLASSIFIER2_LOCALCOMMANDDECLARATIONS_H
#define ADCLASSIFIER2_LOCALCOMMANDDECLARATIONS_H
#include "Command.h"

//extern int download_databases(int argc, const char **argv, const Command& command);
extern int build(int argc, const char **argv, const Command& command);
extern int updataDB(int argc, const char **argv, const Command& command);
extern int classify(int argc, const char **argv, const Command& command);
extern int grade(int argc, const char **argv, const Command& command);
extern int seqHeader2TaxId(int argc, const char **argv, const Command& command);
extern int addToLibrary(int argc, const char **argv, const Command& command);
extern int applyThreshold(int argc, const char **argv, const Command& command);
extern int binning2report(int argc, const char **argv, const Command& command);
extern int filterByGenus(int argc, const char **argv, const Command& command);
#endif //ADCLASSIFIER2_LOCALCOMMANDDECLARATIONS_H
