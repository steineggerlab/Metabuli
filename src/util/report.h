#ifndef METABULI_REPORT_H
#define METABULI_REPORT_H

#include <iostream>
#include <NcbiTaxonomy.h>

using namespace std;
void write_report(FILE *fp, const unordered_map<TaxID, TaxonCounts> &cladeCounts, unsigned long totalReads,
                  NcbiTaxonomy & taxonomy, TaxID taxID = 0, int depth = 0);

void write_report_file(const string &reportFileName, int numOfQuery, unordered_map<TaxID, unsigned int> &taxCnt,
                       NcbiTaxonomy & taxonomy);

unsigned int clade_count_val(const std::unordered_map<TaxID, TaxonCounts> &map, TaxID key);



#endif //METABULI_REPORT_H
