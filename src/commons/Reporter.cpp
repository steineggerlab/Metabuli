#include "Reporter.h"
#include "taxonomyreport.cpp"

Reporter::Reporter(const LocalParameters &par, TaxonomyWrapper *taxonomy, const std::string &customReportFileName) : par(par), taxonomy(taxonomy) {
    if (!customReportFileName.empty()){
        reportFileName = customReportFileName;
    } else {
        if (par.targetTaxId != 0) {return;}
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
            if (par.em) {
                reportFileName_em            = outDir + "/" + jobId + "_EM_report.tsv";
                reportFileName_em_reclassify = outDir + "/" + jobId + "_EM+reclassify_report.tsv";
                reclassifyFileName           = outDir + "/" + jobId + "_EM+reclassify_results.tsv";
                mappingResFileName           = outDir + "/" + jobId + "_mapping_results.txt";
                mappingResBuffer = new WriteBuffer<MappingRes>(mappingResFileName, 1000000);
            }
        }
    }    
}

void Reporter::openReadClassificationFile() {
    readClassificationFile.open(readClassificationFileName);
}

void Reporter::writeReadClassification(const vector<Query> & queryList, bool classifiedOnly) {
    if (isFirstTime) {
        readClassificationFile << "#is_classified\tname\ttaxID\tquery_length\tscore\trank";
        if (par.printLineage) {
            readClassificationFile << "\tlineage";
        }
        readClassificationFile << "\ttaxID:match_count\n";
        isFirstTime = false;
    }
    for (size_t i = 0; i < queryList.size(); i++) {
        if (classifiedOnly && !queryList[i].isClassified) {
            continue;
        }
        if (queryList[i].isClassified != 0) {
            readClassificationFile 
                << queryList[i].isClassified << "\t" 
                << queryList[i].name << "\t"
                << taxonomy->getOriginalTaxID(queryList[i].classification) << "\t"
                << queryList[i].queryLength + queryList[i].queryLength2 << "\t"
                << queryList[i].score << "\t"
                << taxonomy->getString(taxonomy->taxonNode(queryList[i].classification)->rankIdx) << "\t";
            
            if (par.printLineage) {
                readClassificationFile << taxonomy->taxLineage2(taxonomy->taxonNode(queryList[i].classification)) << "\t";
            }
            
            for (auto it = queryList[i].taxCnt.begin(); it != queryList[i].taxCnt.end(); ++it) {
                readClassificationFile << taxonomy->getOriginalTaxID(it->first) << ":" << it->second << " ";
            }
            readClassificationFile << "\n";
        } else {
            readClassificationFile 
                << queryList[i].isClassified << "\t" 
                << queryList[i].name << "\t"
                << taxonomy->getOriginalTaxID(queryList[i].classification) << "\t"
                << queryList[i].queryLength + queryList[i].queryLength2 << "\t"
                << queryList[i].score << "\t"
                << "-" << "\t";
            
            if (par.printLineage) {
                readClassificationFile << "-\t";
            }
            readClassificationFile << "-\t\n";
        }
    }
}

void Reporter::closeReadClassificationFile() {
    readClassificationFile.close();
}

void Reporter::kronaReport(FILE *FP, const TaxonomyWrapper &taxDB, const std::unordered_map<TaxID, TaxonCounts> &cladeCounts, unsigned long totalReads, TaxID taxID, int depth) {
    std::unordered_map<TaxID, TaxonCounts>::const_iterator it = cladeCounts.find(taxID);
    unsigned int cladeCount = it == cladeCounts.end() ? 0 : it->second.cladeCount;
    if (taxID == 0) {
        if (cladeCount > 0) {
            fprintf(FP, "<node name=\"unclassified\"><magnitude><val>%d</val></magnitude></node>", cladeCount);
        }
        kronaReport(FP, taxDB, cladeCounts, totalReads, 1);
    } else {
        if (cladeCount == 0) {
            return;
        }
        const TaxonNode *taxon = taxDB.taxonNode(taxID);
        std::string escapedName = escapeAttribute(taxDB.getString(taxon->nameIdx));
        fprintf(FP, "<node name=\"%s\"><magnitude><val>%d</val></magnitude>", escapedName.c_str(), cladeCount);
        std::vector<TaxID> children = it->second.children;
        SORT_SERIAL(children.begin(), children.end(), [&](int a, int b) { return cladeCountVal(cladeCounts, a) > cladeCountVal(cladeCounts, b); });
        for (size_t i = 0; i < children.size(); ++i) {
            TaxID childTaxId = children[i];
            if (cladeCounts.count(childTaxId)) {
                kronaReport(FP, taxDB, cladeCounts, totalReads, childTaxId, depth + 1);
            } else {
                break;
            }
        }
        fprintf(FP, "</node>");
    }
}

void Reporter::writeReportFile(
    int numOfQuery, 
    unordered_map<TaxID, unsigned int> &taxCnt, 
    ReportType reportType,
    string kronaFileName) 
{
    std::unordered_map<TaxID, std::vector<TaxID>> parentToChildren = taxonomy->getParentToChildren();
    unordered_map<TaxID, TaxonCounts> cladeCounts = taxonomy->getCladeCounts(taxCnt, parentToChildren);
    FILE *fp = nullptr;
    if (reportType == ReportType::Default) {
        fp = fopen(reportFileName.c_str(), "w");
    } else if (reportType == ReportType::EM) {
        fp = fopen(reportFileName_em.c_str(), "w");
    } else if (reportType == ReportType::EM_RECLASSIFY) {
        fp = fopen(reportFileName_em_reclassify.c_str(), "w");
    }
    fprintf(fp, "#clade_proportion\tclade_count\ttaxon_count\trank\ttaxID\tname\n");
    writeReport(fp, cladeCounts, numOfQuery);
    fclose(fp);

    // Write Krona chart
    if (jobId.empty()) { return; }
    
    FILE *kronaFile = nullptr;

    if (reportType == ReportType::Default) {
        if (!kronaFileName.empty()) {
            kronaFile = fopen(kronaFileName.c_str(), "w");
        } else {
            kronaFile = fopen((outDir + "/" + jobId + "_krona.html").c_str(), "w");
        }
    } else if (reportType == ReportType::EM) {
        kronaFile = fopen((outDir + "/" + jobId + "_EM_krona.html").c_str(), "w");
    } else if (reportType == ReportType::EM_RECLASSIFY) {
        kronaFile = fopen((outDir + "/" + jobId + "_EM+reclassify_krona.html").c_str(), "w");
    }
    if (kronaFile == nullptr) {
        Debug(Debug::ERROR) << "Could not open Krona file for writing: " << kronaFileName << "\n";
        EXIT(EXIT_FAILURE);
    }
    fwrite(krona_prelude_html, krona_prelude_html_len, sizeof(char), kronaFile);
    fprintf(kronaFile, "<node name=\"all\"><magnitude><val>%zu</val></magnitude>", (size_t) numOfQuery);
    kronaReport(kronaFile, *taxonomy, cladeCounts, numOfQuery);
    fprintf(kronaFile, "</node></krona></div></body></html>");
    fclose(kronaFile);
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
                taxonomy->getString(taxon->rankIdx), taxonomy->getOriginalTaxID(taxID), std::string(2 * depth, ' ').c_str(), taxonomy->getString(taxon->nameIdx));
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

// void Reporter::writeEMreportFile(
//     unordered_map<TaxID, double> &taxProbs) 
// {
//     std::unordered_map<TaxID, std::vector<TaxID>> parentToChildren = taxonomy->getParentToChildren();
//     unordered_map<TaxID, TaxonCounts> cladeCounts = taxonomy->getCladeProbs(taxCnt, parentToChildren);
//     FILE *fp = fopen(reportFileName_em.c_str(), "w");
//     fprintf(fp, "#clade_proportion\tclade_count\ttaxon_count\trank\ttaxID\tname\n");
//     writeReport(fp, cladeCounts, numOfQuery);
//     fclose(fp);

//     // Write Krona chart
//     if (jobId.empty()) {
//         return;
//     }
//     FILE *kronaFile = nullptr;
//     if (!em) { 
//         if (!kronaFileName.empty()) {
//             kronaFile = fopen(kronaFileName.c_str(), "w");
//         } else {
//             kronaFile = fopen((outDir + "/" + jobId + "_krona.html").c_str(), "w");
//         }
//     } else {
//         kronaFile = fopen((outDir + "/" + jobId + "_EM_krona.html").c_str(), "w");
//     }
//     fwrite(krona_prelude_html, krona_prelude_html_len, sizeof(char), kronaFile);
//     fprintf(kronaFile, "<node name=\"all\"><magnitude><val>%zu</val></magnitude>", (size_t) numOfQuery);
//     kronaReport(kronaFile, *taxonomy, cladeCounts, numOfQuery);
//     fprintf(kronaFile, "</node></krona></div></body></html>");
//     fclose(kronaFile);
// }


unsigned int Reporter::cladeCountVal(const std::unordered_map<TaxID, TaxonCounts> &map, TaxID key) {
    typename std::unordered_map<TaxID, TaxonCounts>::const_iterator it = map.find(key);
    if (it == map.end()) {
        return 0;
    } else {
        return it->second.cladeCount;
    }
}

void Reporter::getReadsClassifiedToClade(TaxID cladeId,
                                         const string &readClassificationFileName,
                                         vector<size_t> &readIdxs) {
    FILE *results = fopen(readClassificationFileName.c_str(), "r");
    if (!results) {
        perror("Failed to open read-by-read classification file");
        return;
    }
    char line[4096];
    size_t idx = 0;
    if (taxonomy->hasInternalTaxID()) {
        unordered_map<TaxID, TaxID> extern2intern;
        taxonomy->getExternal2internalTaxID(extern2intern);
        while (fgets(line, sizeof(line), results)) {
            if (line[0] == '#') {
                continue;
            }
            int taxId;
            if (sscanf(line, "%*s %*s %d", &taxId) == 1) {            
                if (taxonomy->IsAncestor(cladeId, extern2intern[taxId])) {
                    readIdxs.push_back(idx);
                }
            }
            idx++;
        }
    } else {
        while (fgets(line, sizeof(line), results)) {
            if (line[0] == '#') {
                continue;
            }
            int taxId;
            if (sscanf(line, "%*s %*s %d", &taxId) == 1) {            
                if (taxonomy->IsAncestor(cladeId, taxId)) {
                    readIdxs.push_back(idx);
                }
            }
            idx++;
        }
    }
    fclose(results);
}

void Reporter::printSpecifiedReads(const vector<size_t> & readIdxs,
                                   const string & readFileName,
                                   string & outFileName) {
    // Check FASTA or FASTQ
    KSeqWrapper* tempKseq = KSeqFactory(readFileName.c_str());
    tempKseq->ReadEntry();
    bool isFasta = tempKseq->entry.qual.l == 0;

    if (isFasta && par.extractMode == 2) {
        Debug(Debug::ERROR) << "Cannot convert FASTA to FASTQ\n";
        EXIT(EXIT_FAILURE);
    }

    delete tempKseq;

    bool printFasta;
    if (isFasta || par.extractMode == 1) {
        printFasta = true;
        outFileName += ".fna";
    } else {
        printFasta = false;
        outFileName += ".fq";
    }
    
    KSeqWrapper* kseq = KSeqFactory(readFileName.c_str());
    FILE *outFile = fopen(outFileName.c_str(), "w");
    if (!outFile) {
        perror("Failed to open file");
        return;
    }

    size_t readCnt = 0;
    size_t idx = 0;

    if (printFasta) {
        while (kseq->ReadEntry()) {
            if (readCnt == readIdxs[idx]) {
                fprintf(outFile, ">%s", kseq->entry.name.s);
                if (kseq->entry.comment.l > 0) {
                    fprintf(outFile, " %s\n", kseq->entry.comment.s);
                } else {
                    fprintf(outFile, "\n");
                }
                fprintf(outFile, "%s\n", kseq->entry.sequence.s);
                idx++;
                if (idx == readIdxs.size()) {
                    break;
                }
            }
            readCnt++;
        }
    } else {
        while (kseq->ReadEntry()) {
            if (readCnt == readIdxs[idx]) {
                fprintf(outFile, "@%s", kseq->entry.name.s);
                if (kseq->entry.comment.l > 0) {
                    fprintf(outFile, " %s\n", kseq->entry.comment.s);
                } else {
                    fprintf(outFile, "\n");
                }
                fprintf(outFile, "%s\n", kseq->entry.sequence.s);
                fprintf(outFile, "+%s", kseq->entry.name.s);
                if (kseq->entry.comment.l > 0) {
                    fprintf(outFile, " %s\n", kseq->entry.comment.s);
                } else {
                    fprintf(outFile, "\n");
                }
                fprintf(outFile, "%s\n", kseq->entry.qual.s);
                idx++;
                if (idx == readIdxs.size()) {
                    break;
                }
            }
            readCnt++;
        }
    }
    delete kseq;
}

void Reporter::writeReclassifyResults(const std::vector<Classification> & results)
{   
    ofstream emResultFile(reclassifyFileName, std::ios::out | std::ios::trunc);
    if (!emResultFile.is_open()) {
        cerr << "Error: Could not open EM results file " << reclassifyFileName << endl;
        return;
    }

    emResultFile << "#is_classified\tname\ttaxID\tquery_length\tscore\trank";
    if (par.printLineage) {
        emResultFile << "\tlineage";
    }
    emResultFile << endl;

    for (size_t i = 0; i < results.size(); ++i) {
        const Classification &result = results[i];
        if (result.taxId != 0) {
            emResultFile << (result.taxId != 0) << "\t" 
                         << result.name << "\t"
                         << taxonomy->getOriginalTaxID(result.taxId) << "\t"
                         << result.length << "\t"
                         << result.score << "\t"
                         << taxonomy->getString(taxonomy->taxonNode(result.taxId)->rankIdx);
            if (par.printLineage) {
                emResultFile << "\t" << taxonomy->taxLineage2(taxonomy->taxonNode(result.taxId));
            }
        } else {
            emResultFile << (result.taxId != 0) << "\t" 
                         << result.name << "\t"
                         << taxonomy->getOriginalTaxID(result.taxId) << "\t"
                         << result.length << "\t"
                         << result.score << "\t"
                         << "-";
            if (par.printLineage) {
                emResultFile << "\t-";
            }        
        }
        emResultFile << "\n";
    }

    emResultFile.close();
    cout << "EM results written to " << reclassifyFileName << endl;
}