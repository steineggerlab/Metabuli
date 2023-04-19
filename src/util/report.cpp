#include "report.h"

void write_report_file(const string &reportFileName, int numOfQuery, unordered_map<TaxID, unsigned int> &taxCnt,
                       NcbiTaxonomy & taxonomy) {
    unordered_map<TaxID, TaxonCounts> cladeCounts = taxonomy.getCladeCounts(taxCnt);
    FILE *fp;
    fp = fopen(reportFileName.c_str(), "w");
    write_report(fp, cladeCounts, numOfQuery, taxonomy);
    fclose(fp);
}

void write_report(FILE *FP, const unordered_map<TaxID, TaxonCounts> &cladeCounts, unsigned long totalReads,
                  NcbiTaxonomy & taxonomy, TaxID taxID, int depth) {
    std::unordered_map<TaxID, TaxonCounts>::const_iterator it = cladeCounts.find(taxID);
    unsigned int cladeCount = it == cladeCounts.end() ? 0 : it->second.cladeCount;
    unsigned int taxCount = it == cladeCounts.end() ? 0 : it->second.taxCount;
    if (taxID == 0) {
        if (cladeCount > 0) {
            fprintf(FP, "%.4f\t%i\t%i\tno rank\t0\tunclassified\n",
                    100 * cladeCount / double(totalReads),
                    cladeCount, taxCount);
        }
        write_report(FP, cladeCounts, totalReads, taxonomy, 1);
    } else {
        if (cladeCount == 0) {
            return;
        }
        const TaxonNode *taxon = taxonomy.taxonNode(taxID);
        fprintf(FP, "%.4f\t%i\t%i\t%s\t%i\t%s%s\n",
                100 * cladeCount / double(totalReads), cladeCount, taxCount,
                taxonomy.getString(taxon->rankIdx), taxID, std::string(2 * depth, ' ').c_str(), taxonomy.getString(taxon->nameIdx));
        std::vector<TaxID> children = it->second.children;
        SORT_SERIAL(children.begin(), children.end(), [&](int a, int b) { return clade_count_val(cladeCounts, a) > clade_count_val(cladeCounts, b); });
        for (size_t i = 0; i < children.size(); ++i) {
            TaxID childTaxId = children[i];
            if (cladeCounts.count(childTaxId)) {
                write_report(FP, cladeCounts, totalReads, taxonomy, childTaxId, depth + 1);
            } else {
                break;
            }
        }
    }


//    auto it = cladeCounts.find(taxID);
//    unsigned int cladeCount = (it == cladeCounts.end() ? 0 : it->second.cladeCount);
//    unsigned int taxCount = (it == cladeCounts.end() ? 0 : it->second.taxCount);
//    if (taxID == 0) {
//        if (cladeCount > 0) {
//            fprintf(fp, "%.2f\t%i\t%i\t0\tno rank\tunclassified\n", 100 * cladeCount / double(totalReads), cladeCount,
//                    taxCount);
//        }
//        write_report(fp, cladeCounts, totalReads, taxonomy, 1);
//    } else {
//        if (cladeCount == 0) {
//            return;
//        }
//        const TaxonNode *taxon = taxonomy.taxonNode(taxID);
//        fprintf(fp, "%.2f\t%i\t%i\t%i\t%s\t%s%s\n", 100 * cladeCount / double(totalReads), cladeCount, taxCount, taxID,
//                taxon->rank.c_str(), string(2 * depth, ' ').c_str(), taxon->name.c_str());
//        vector<TaxID> children = it->second.children;
//        sort(children.begin(), children.end(),
//             [&](int a, int b) { return clade_count_val(cladeCounts, a) > clade_count_val(cladeCounts, b); });
//        for (TaxID childTaxId: children) {
//            if (cladeCounts.count(childTaxId)) {
//                write_report(fp,  cladeCounts, totalReads, taxonomy, childTaxId, depth + 1);
//            } else {
//                break;
//            }
//        }
//    }
}

unsigned int clade_count_val(const std::unordered_map<TaxID, TaxonCounts> &map, TaxID key) {
    typename std::unordered_map<TaxID, TaxonCounts>::const_iterator it = map.find(key);
    if (it == map.end()) {
        return 0;
    } else {
        return it->second.cladeCount;
    }
}
