#include "LocalParameters.h"
#include <iostream>
#include <fstream>

using namespace std;

struct ICTV_element {
    string sequenceId;
    string realm;
    float realmScore;
    string subRealm;
    float subRealmScore;
    string kingdom;
    float kingdomScore;
    string subKingdom;
    float subKingdomScore;
    string phylum;
    float phylumScore;
    string subPhylum;
    float subPhylumScore;
    string class_;
    float classScore;
    string subClass;
    float subClassScore;
    string order;
    float orderScore;
    string subOrder;
    float subOrderScore;
    string family;
    float familyScore;
    string subFamily;
    float subFamilyScore;
    string genus;
    float genusScore;
    string subGenus;
    float subGenusScore;
    string species;
    float speciesScore;

    ICTV_element() {
        sequenceId = "";
        realm = "";
        realmScore = 0;
        subRealm = "";
        subRealmScore = 0;
        kingdom = "";
        kingdomScore = 0;
        subKingdom = "";
        subKingdomScore = 0;
        phylum = "";
        phylumScore = 0;
        subPhylum = "";
        subPhylumScore = 0;
        class_ = "";
        classScore = 0;
        subClass = "";
        subClassScore = 0;
        order = "";
        orderScore = 0;
        subOrder = "";
        subOrderScore = 0;
        family = "";
        familyScore = 0;
        subFamily = "";
        subFamilyScore = 0;
        genus = "";
        genusScore = 0;
        subGenus = "";
        subGenusScore = 0;
        species = "";
        speciesScore = 0;
    }

    void printLine(ofstream &out) {
        out << sequenceId << ",";
        out << realm << "," << realmScore << ",";
        out << subRealm << "," << subRealmScore << ",";
        out << kingdom << "," << kingdomScore << ",";
        out << subKingdom << "," << subKingdomScore << ",";
        out << phylum << "," << phylumScore << ",";
        out << subPhylum << "," << subPhylumScore << ",";
        out << class_ << "," << classScore << ",";
        out << subClass << "," << subClassScore << ",";
        out << order << "," << orderScore << ",";
        out << subOrder << "," << subOrderScore << ",";
        out << family << "," << familyScore << ",";
        out << subFamily << "," << subFamilyScore << ",";
        out << genus << "," << genusScore << ",";
        out << subGenus << "," << subGenusScore << ",";
        out << species << "," << speciesScore << endl;
    }

    // print score only when it is not 0
    void printLine2(ofstream &out) {
        out << sequenceId << ",";
        if (realmScore != 0) {
            out << realm << "," << realmScore << ",";
        } else {
            out << realm << ",,";
        }
        if (subRealmScore != 0) {
            out << subRealm << "," << subRealmScore << ",";
        } else {
            out << subRealm << ",,";
        }
        if (kingdomScore != 0) {
            out << kingdom << "," << kingdomScore << ",";
        } else {
            out << kingdom << ",,";
        }
        if (subKingdomScore != 0) {
            out << subKingdom << "," << subKingdomScore << ",";
        } else {
            out << subKingdom << ",,";
        }
        if (phylumScore != 0) {
            out << phylum << "," << phylumScore << ",";
        } else {
            out << phylum << ",,";
        }
        if (subPhylumScore != 0) {
            out << subPhylum << "," << subPhylumScore << ",";
        } else {
            out << subPhylum << ",,";
        }
        if (classScore != 0) {
            out << class_ << "," << classScore << ",";
        } else {
            out << class_ << ",,";
        }
        if (subClassScore != 0) {
            out << subClass << "," << subClassScore << ",";
        } else {
            out << subClass << ",,";
        }
        if (orderScore != 0) {
            out << order << "," << orderScore << ",";
        } else {
            out << order << ",,";
        }
        if (subOrderScore != 0) {
            out << subOrder << "," << subOrderScore << ",";
        } else {
            out << subOrder << ",,";
        }
        if (familyScore != 0) {
            out << family << "," << familyScore << ",";
        } else {
            out << family << ",,";
        }
        if (subFamilyScore != 0) {
            out << subFamily << "," << subFamilyScore << ",";
        } else {
            out << subFamily << ",,";
        }
        if (genusScore != 0) {
            out << genus << "," << genusScore << ",";
        } else {
            out << genus << ",,";
        }
        if (subGenusScore != 0) {
            out << subGenus << "," << subGenusScore << ",";
        } else {
            out << subGenus << ",,";
        }
        if (speciesScore != 0) {
            out << species << "," << speciesScore << endl;
        } else {
            out << species << endl;
        }
    }
};

int ictvFormat(int argc, const char **argv, const Command &command) {
    LocalParameters &par = LocalParameters::getLocalInstance();
    par.parseParameters(argc, argv, command, false, Parameters::PARSE_ALLOW_EMPTY, 0);
    cout << par.filenames[0] << endl;
    std::string resultFile = par.filenames[0];
    // Remove the extension
    std::string outputFile;
    size_t lastDotPos = resultFile.find_last_of('.');
    if (lastDotPos != std::string::npos) {
        outputFile = resultFile.substr(0, lastDotPos);
    }
    outputFile += "_ictv.csv";

    std::ofstream out(outputFile);
    if (!out.is_open()) {
        std::cerr << "Error: Cannot open file " << outputFile << std::endl;
        exit(1);
    }

    // Print header
    out << "SequenceID,Realm (-viria),Realm_score,Subrealm (-vira),Subrealm_score,Kingdom (-virae),Kingdom_score,Subkingdom (-virites),Subkingdom_score,Phylum (-viricota),Phylum_score,Subphylum (-viricotina),Subphylum_score,Class (-viricetes),Class_score,Subclass (-viricetidae),Subclass_score,Order (-virales),Order_score,Suborder (-virineae),Suborder_score,Family (-viridae),Family_score,Subfamily (-virinae),Subfamily_score,Genus (-virus),Genus_score,Subgenus (-virus),Subgenus_score,Species (binomial),Species_score" << endl;

    std::ifstream in(resultFile);
    if (!in.is_open()) {
        std::cerr << "Error: Cannot open file " << resultFile << std::endl;
        exit(1);
    }
    string line;
    while (getline(in, line)) {
        // split line by tab
        vector<string> fields = Util::split(line, "\t");
        string sequenceId = fields[1];
        float score = stof(fields[4]);
        string lineage = fields[6];

        ICTV_element ictv;
        ictv.sequenceId = sequenceId;

        // split lineage by ;
        vector<string> ranks = Util::split(lineage, ";");
        
        // iterate through ranks
        for (size_t i = 0; i < ranks.size(); i++) {
            // split a rank by '_'
            vector<string> rank = Util::split(ranks[i], "_");
            
            if (rank[0] == "r") {
                ictv.realm = rank[1];
                ictv.realmScore = score;
            } else if (rank[0] == "sr") {
                ictv.subRealm = rank[1];
                ictv.subRealmScore = score;
            } else if (rank[0] == "k") {
                ictv.kingdom = rank[1];
                ictv.kingdomScore = score;
            } else if (rank[0] == "sk") {
                ictv.subKingdom = rank[1];
                ictv.subKingdomScore = score;
            } else if (rank[0] == "p") {
                ictv.phylum = rank[1];
                ictv.phylumScore = score;
            } else if (rank[0] == "sp") {
                ictv.subPhylum = rank[1];
                ictv.subPhylumScore = score;
            } else if (rank[0] == "c") {
                ictv.class_ = rank[1];
                ictv.classScore = score;
            } else if (rank[0] == "sc") {
                ictv.subClass = rank[1];
                ictv.subClassScore = score;
            } else if (rank[0] == "o") {
                ictv.order = rank[1];
                ictv.orderScore = score;
            } else if (rank[0] == "so") {
                ictv.subOrder = rank[1];
                ictv.subOrderScore = score;
            } else if (rank[0] == "f") {
                ictv.family = rank[1];
                ictv.familyScore = score;
            } else if (rank[0] == "sf") {
                ictv.subFamily = rank[1];
                ictv.subFamilyScore = score;
            } else if (rank[0] == "g") {
                ictv.genus = rank[1];
                ictv.genusScore = score;
            } else if (rank[0] == "sg") {
                ictv.subGenus = rank[1];
                ictv.subGenusScore = score;
            } else if (rank[0] == "s") {
                ictv.species = rank[1];
                ictv.speciesScore = score;
            }
        }
        ictv.printLine2(out);
    }
    in.close();
    out.close();
    return 0;
}