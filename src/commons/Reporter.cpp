#include "Reporter.h"
#include "taxonomyreport.cpp"

Reporter::Reporter(const LocalParameters &par, NcbiTaxonomy *taxonomy) : taxonomy(taxonomy) {
    if (par.contamList == "") { // classify module
        if (par.seqMode == 2) {
            outDir = par.filenames[3];
            jobId = par.filenames[4];
        } else {
            outDir = par.filenames[2];
            jobId = par.filenames[3];
        }
        // Output file names
        reportFileName = outDir + + "/" + jobId + "_report.tsv";
        readClassificationFileName = outDir + "/" + jobId + "_classifications.tsv";
    }

    
    
}

void Reporter::openReadClassificationFile() {
    readClassificationFile.open(readClassificationFileName);
}

void Reporter::writeReadClassification(const vector<Query> & queryList, bool classifiedOnly) {
    for (size_t i = 0; i < queryList.size(); i++) {
        if (classifiedOnly && !queryList[i].isClassified) {
            continue;
        }
        readClassificationFile << queryList[i].isClassified << "\t" << queryList[i].name << "\t"
                               << queryList[i].classification << "\t"
                               << queryList[i].queryLength + queryList[i].queryLength2 << "\t"
                               << queryList[i].score << "\t"
                               << taxonomy->getString(taxonomy->taxonNode(queryList[i].classification)->rankIdx) << "\t";
        for (auto it = queryList[i].taxCnt.begin(); it != queryList[i].taxCnt.end(); ++it) {
            readClassificationFile << it->first << ":" << it->second << " ";
        }
        readClassificationFile << "\n";
    }
}

void Reporter::closeReadClassificationFile() {
    readClassificationFile.close();
}

void Reporter::writeReportFile(int numOfQuery, unordered_map<TaxID, unsigned int> &taxCnt, bool krona) {
    unordered_map<TaxID, TaxonCounts> cladeCounts = taxonomy->getCladeCounts(taxCnt);
    FILE *fp;
    fp = fopen((reportFileName).c_str(), "w");
    writeReport(fp, cladeCounts, numOfQuery);
    fclose(fp);

    // Write Krona chart
    if (krona) {
        FILE *kronaFile = fopen((outDir + "/" + jobId + "_krona.html").c_str(), "w");
        fwrite(krona_prelude_html, krona_prelude_html_len, sizeof(char), kronaFile);
        fprintf(kronaFile, "<node name=\"all\"><magnitude><val>%zu</val></magnitude>", numOfQuery);
        kronaReport(kronaFile, *taxonomy, cladeCounts, numOfQuery);
        fprintf(kronaFile, "</node></krona></div></body></html>");
        fclose(kronaFile);
    }
}

void Reporter::writeReport(FILE *FP, const std::unordered_map<TaxID, TaxonCounts> &cladeCounts,
                             unsigned long totalReads, TaxID taxID, int depth) {
    std::unordered_map<TaxID, TaxonCounts>::const_iterator it = cladeCounts.find(taxID);
    unsigned int cladeCount = it == cladeCounts.end() ? 0 : it->second.cladeCount;
    unsigned int taxCount = it == cladeCounts.end() ? 0 : it->second.taxCount;
    if (taxID == 0) {
        if (cladeCount > 0) {
            fprintf(FP, "%.4f\t%i\t%i\tno rank\t0\tunclassified\n",
                    100 * cladeCount / double(totalReads),
                    cladeCount, taxCount);
        }
        writeReport(FP, cladeCounts, totalReads, 1);
    } else {
        if (cladeCount == 0) {
            return;
        }
        const TaxonNode *taxon = taxonomy->taxonNode(taxID);
        fprintf(FP, "%.4f\t%i\t%i\t%s\t%i\t%s%s\n",
                100 * cladeCount / double(totalReads), cladeCount, taxCount,
                taxonomy->getString(taxon->rankIdx), taxID, std::string(2 * depth, ' ').c_str(), taxonomy->getString(taxon->nameIdx));
        std::vector<TaxID> children = it->second.children;
        SORT_SERIAL(children.begin(), children.end(), [&](int a, int b) { return cladeCountVal(cladeCounts, a) > cladeCountVal(cladeCounts, b); });
        for (size_t i = 0; i < children.size(); ++i) {
            TaxID childTaxId = children[i];
            if (cladeCounts.count(childTaxId)) {
                writeReport(FP, cladeCounts, totalReads, childTaxId, depth + 1);
            } else {
                break;
            }
        }
    }
}

unsigned int Reporter::cladeCountVal(const std::unordered_map<TaxID, TaxonCounts> &map, TaxID key) {
    typename std::unordered_map<TaxID, TaxonCounts>::const_iterator it = map.find(key);
    if (it == map.end()) {
        return 0;
    } else {
        return it->second.cladeCount;
    }
}