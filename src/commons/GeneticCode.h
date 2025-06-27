#ifndef METABULI_GENETIC_CODE_H
#define METABULI_GENETIC_CODE_H

#include <iostream>

#define nuc2int(x) (x & 14u)>>1u
class GeneticCode {
    public:
        const std::string atcg = "................................................................"
                                 ".AGCG..GT..G.CN...ACTG.A.T.......agcg..gt..g.cn...actg.a.t......"
                                 "................................................................"
                                 "................................................................";

        const std::string iRCT = "................................................................"
                                 ".TVGH..CD..M.KN...YSAABW.R.......tvgh..cd..m.kn...ysaabw.r......"
                                 "................................................................"
                                 "................................................................";
        int nuc2aa[8][8][8];
        int nuc2num[8][8][8];
        std::string aminoacids;
        std::vector<std::vector<std::string>> aa2codon;

        int getAA(const char nuc1, const char nuc2, const char nuc3) const {
            return nuc2aa[nuc2int(nuc1)][nuc2int(nuc2)][nuc2int(nuc3)];
        }

        int getCodon(const char nuc1, const char nuc2, const char nuc3) const {
            return nuc2num[nuc2int(nuc1)][nuc2int(nuc2)][nuc2int(nuc3)];
        }
                

        explicit GeneticCode(bool reducedAA) {
            if (!reducedAA) {
                aminoacids = "ARNDCQEGHILKMFPSTWYVX";
                // A
                nuc2aa[3][1][0] = 0;
                nuc2aa[3][1][1] = 0;
                nuc2aa[3][1][2] = 0; 
                nuc2aa[3][1][3] = 0;
                aa2codon.push_back({"GCA", "GCC", "GCT", "GCG"});
                
                // R
                nuc2aa[1][3][0] = 1;
                nuc2aa[1][3][1] = 1;
                nuc2aa[1][3][2] = 1;
                nuc2aa[1][3][3] = 1; 
                nuc2aa[0][3][0] = 1; 
                nuc2aa[0][3][3] = 1;
                aa2codon.push_back({"CGA", "CGC", "CGT", "CGG", "AGG", "AGA"});

                // N
                nuc2aa[0][0][2] = 2;
                nuc2aa[0][0][1] = 2;
                aa2codon.push_back({"ERR", "AAC", "AAT"});

                // D
                nuc2aa[3][0][2] = 3;
                nuc2aa[3][0][1] = 3;
                aa2codon.push_back({"ERR", "GAC", "GAT"});
                
                // C
                nuc2aa[2][3][2] = 4;
                nuc2aa[2][3][1] = 4;
                aa2codon.push_back({"ERR", "TGC", "TGT"});

                // Q
                nuc2aa[1][0][0] = 5;
                nuc2aa[1][0][3] = 5;
                aa2codon.push_back({"CAA", "ERR", "ERR", "CAG"});

                // E
                nuc2aa[3][0][0] = 6;
                nuc2aa[3][0][3] = 6;
                aa2codon.push_back({"GAA", "ERR", "ERR", "GAG"});

                // G
                nuc2aa[3][3][0] = 7;
                nuc2aa[3][3][1] = 7;
                nuc2aa[3][3][2] = 7;
                nuc2aa[3][3][3] = 7;
                aa2codon.push_back({"GGA", "GGC", "GGT", "GGG"});

                // H
                nuc2aa[1][0][2] = 8;
                nuc2aa[1][0][1] = 8;
                aa2codon.push_back({"ERR", "CAC", "CAT"});

                // I
                nuc2aa[0][2][2] = 9;
                nuc2aa[0][2][1] = 9;
                nuc2aa[0][2][0] = 9;
                aa2codon.push_back({"ATA", "ATC", "ATT"});
                
                // L
                nuc2aa[2][2][0] = 10;
                nuc2aa[2][2][3] = 10;
                nuc2aa[1][2][0] = 10;
                nuc2aa[1][2][1] = 10;
                nuc2aa[1][2][2] = 10;
                nuc2aa[1][2][3] = 10;
                aa2codon.push_back({"CTA", "CTC", "CTT", "CTG", "TTG", "TTA"});
                
                // K
                nuc2aa[0][0][0] = 11; 
                nuc2aa[0][0][3] = 11;
                aa2codon.push_back({"AAA", "ERR", "ERR", "AAG"});

                // M
                nuc2aa[0][2][3] = 12;
                aa2codon.push_back({"ERR", "ERR", "ERR", "ATG"});
                
                // F
                nuc2aa[2][2][2] = 13;
                nuc2aa[2][2][1] = 13;
                aa2codon.push_back({"ERR", "TTC", "TTT"});

                // P
                nuc2aa[1][1][0] = 14; 
                nuc2aa[1][1][1] = 14;
                nuc2aa[1][1][2] = 14;
                nuc2aa[1][1][3] = 14;
                aa2codon.push_back({"CCA", "CCC", "CCT", "CCG"});

                // S
                nuc2aa[2][1][0] = 15;
                nuc2aa[2][1][1] = 15;
                nuc2aa[2][1][2] = 15;
                nuc2aa[2][1][3] = 15;
                nuc2aa[0][3][2] = 15;
                nuc2aa[0][3][1] = 15;
                aa2codon.push_back({"TCA", "TCC", "TCT", "TCG", "ERR", "ERR", "AGT", "AGC"});
                
                // T
                nuc2aa[0][1][0] = 16;
                nuc2aa[0][1][1] = 16;
                nuc2aa[0][1][2] = 16;
                nuc2aa[0][1][3] = 16;
                aa2codon.push_back({"ACA", "ACC", "ACT", "ACG"});

                // W
                nuc2aa[2][3][3] = 17;
                aa2codon.push_back({"ERR", "ERR", "ERR", "TGG"});
                
                // Y
                nuc2aa[2][0][2] = 18;
                nuc2aa[2][0][1] = 18;
                aa2codon.push_back({"ERR", "TAC", "TAT"});

                // V
                nuc2aa[3][2][0] = 19;
                nuc2aa[3][2][1] = 19;
                nuc2aa[3][2][2] = 19;
                nuc2aa[3][2][3] = 19;
                aa2codon.push_back({"GTA", "GTC", "GTT", "GTG"});

                // Stop
                nuc2aa[2][0][0] = 20; 
                nuc2aa[2][3][0] = 20;
                nuc2aa[2][0][3] = 20;
                aa2codon.push_back({"TAA", "ERR", "ERR", "TAG", "ERR", "TGA"});

                // triplet code with N's
                for(int i = 0; i < 8; i++){
                    for(int j = 0; j < 8; j++){
                        nuc2aa[7][i][j] = -1;
                        nuc2aa[i][7][j] = -1;
                        nuc2aa[i][j][7] = -1;
                        nuc2num[7][i][j] = -1;
                        nuc2num[i][7][j] = -1;
                        nuc2num[i][j][7] = -1;
                    }
                }

                // For encoding DNA information in k-mer
                for(int i =0; i < 4 ; i++){
                    for(int i2 = 0; i2 < 4 ; i2++){
                        nuc2num[i][i2][0] = 0;
                        nuc2num[i][i2][1] = 1;
                        nuc2num[i][i2][2] = 2;
                        nuc2num[i][i2][3] = 3;
                    }
                }
                // for Arg
                nuc2num[0][3][3] = 4;
                nuc2num[0][3][0] = 5;
                // for Leu
                nuc2num[2][2][3] = 4;
                nuc2num[2][2][0] = 5;
                // for Ser
                nuc2num[0][3][2] = 6;
                nuc2num[0][3][1] = 7;
                // for stop codon
                nuc2num[2][3][0] = 5;

            } else {
                aminoacids = "ARNDCQGHILKFPSTX";
                // Codon table
                // A
                nuc2aa[3][1][0] = 0;
                nuc2aa[3][1][1] = 0;
                nuc2aa[3][1][2] = 0;
                nuc2aa[3][1][3] = 0;
                // R
                nuc2aa[1][3][0] = 1;
                nuc2aa[1][3][1] = 1;
                nuc2aa[1][3][2] = 1;
                nuc2aa[1][3][3] = 1;
                nuc2aa[0][3][0] = 1;
                nuc2aa[0][3][3] = 1;
                // N
                nuc2aa[0][0][2] = 2;
                nuc2aa[0][0][1] = 2;
                // D
                nuc2aa[3][0][2] = 3;
                nuc2aa[3][0][1] = 3;
                // C
                nuc2aa[2][3][2] = 4;
                nuc2aa[2][3][1] = 4;
                // QE
                nuc2aa[1][0][0] = 5;
                nuc2aa[1][0][3] = 5;
                nuc2aa[3][0][0] = 5;
                nuc2aa[3][0][3] = 5;
                // G
                nuc2aa[3][3][0] = 6;
                nuc2aa[3][3][1] = 6;
                nuc2aa[3][3][2] = 6;
                nuc2aa[3][3][3] = 6;
                // H
                nuc2aa[1][0][2] = 7;
                nuc2aa[1][0][1] = 7;
                // IV
                nuc2aa[0][2][2] = 8;
                nuc2aa[0][2][1] = 8;
                nuc2aa[0][2][0] = 8;
                nuc2aa[3][2][0] = 8;
                nuc2aa[3][2][1] = 8;
                nuc2aa[3][2][2] = 8;
                nuc2aa[3][2][3] = 8;
                // ML
                nuc2aa[2][2][0] = 9;
                nuc2aa[2][2][3] = 9;
                nuc2aa[1][2][0] = 9;
                nuc2aa[1][2][1] = 9;
                nuc2aa[1][2][2] = 9;
                nuc2aa[1][2][3] = 9;
                nuc2aa[0][2][3] = 9;
                // K
                nuc2aa[0][0][0] = 10;
                nuc2aa[0][0][3] = 10;
                // FYW
                nuc2aa[2][2][2] = 11;
                nuc2aa[2][2][1] = 11;
                nuc2aa[2][0][2] = 11;
                nuc2aa[2][0][1] = 11;
                nuc2aa[2][3][3] = 11;
                // P
                nuc2aa[1][1][0] = 12;
                nuc2aa[1][1][1] = 12;
                nuc2aa[1][1][2] = 12;
                nuc2aa[1][1][3] = 12;
                // S
                nuc2aa[2][1][0] = 13;
                nuc2aa[2][1][1] = 13;
                nuc2aa[2][1][2] = 13;
                nuc2aa[2][1][3] = 13;
                nuc2aa[0][3][2] = 13;
                nuc2aa[0][3][1] = 13;
                // T
                nuc2aa[0][1][0] = 14;
                nuc2aa[0][1][1] = 14;
                nuc2aa[0][1][2] = 14;
                nuc2aa[0][1][3] = 14;
                // Stop
                nuc2aa[2][0][0] = 15;
                nuc2aa[2][3][0] = 15;
                nuc2aa[2][0][3] = 15;
                // Triplet code with N's
                for (int i = 0; i < 8; i++) {
                    for (int j = 0; j < 8; j++) {
                        nuc2aa[7][i][j] = -1;
                        nuc2aa[i][7][j] = -1;
                        nuc2aa[i][j][7] = -1;
                    }
                }

                // Codon encoding to store which DNA triplet code is used.
                // This encoding is used for indexing hamming distance look-up table.
                for (int i = 0; i < 4; i++) {
                    for (int i2 = 0; i2 < 4; i2++) {
                        nuc2num[i][i2][0] = 0;
                        nuc2num[i][i2][1] = 1;
                        nuc2num[i][i2][2] = 2;
                        nuc2num[i][i2][3] = 3;
                    }
                }
                // for Arg
                nuc2num[0][3][3] = 7; //AGG
                nuc2num[0][3][0] = 4; //AGA
                // for Leu
                nuc2num[2][2][3] = 7; //TTG
                nuc2num[2][2][0] = 4; //TTA
                nuc2num[0][2][3] = 8; //ATG
                // for Ser
                nuc2num[0][3][2] = 10; //AGT
                nuc2num[0][3][1] = 9; //AGC
                // for FYW
                nuc2num[2][0][1] = 5;
                nuc2num[2][0][2] = 6;
                nuc2num[2][3][3] = 7;
                // for IV
                nuc2num[0][2][0] = 4;
                nuc2num[0][2][1] = 5;
                nuc2num[0][2][2] = 6;
                // for QE
                nuc2num[3][0][0] = 4; //GAA
                nuc2num[3][0][3] = 7; //GAG
                // for stop codon
                nuc2num[2][3][0] = 4; //TGA
            }
        }  


};

#endif //METABULI_GENETIC_CODE_H