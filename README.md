# Metabuli
Metabuli taxonomically classifies metagenomic reads using both DNA and amino acid (AA) information.
It achieved specificity of DNA-based method and sensitivity of AA-method at the same time.


## Installation
### Precompiled binaries
```
# Linux AVX2 build (fast, recommended for most Linux system) (check using: cat /proc/cpuinfo | grep avx2)
wget https://mmseqs.com/metabuli/metabuli-linux-avx2.tar.gz; tar xvzf metabuli-linux-avx2.tar.gz; export PATH=$(pwd)/metabuli/bin/:$PATH

# Linux SSE2 build (slow, for old systems)
wget https://mmseqs.com/metabuli/metabuli-linux-sse2.tar.gz; tar xvzf metabuli-linux-sse2.tar.gz; export PATH=$(pwd)/metabuli/bin/:$PATH

# MacOS
wget https://mmseqs.com/metabuli/metabuli-osx-universal.tar.gz; tar xvzf metabuli-osx-universal.tar.gz; export PATH=$(pwd)/metabuli/bin/:$PATH
```
Metabuli also works on ARM64 systems. Please check [https://mmseqs.com/metabuli/](https://mmseqs.com/metabuli/)
### Compile from source code
Installation from Github source code.
```
git clone GITHUB_LINK
cd Metabuli
mkdir build
cd build
cmake ..
make -j 16
```
The built binary can be found in ./build/src

## Pre-built databases
You can download pre-built databases using `databases` command.
```
# RefSeq Complete/Chromosome
# - Complete Genome or Chromosome level assemblies of virus and prokaryotes in RefSeq (2023-04-04) and human genome (GRCh38.p14)
metabuli databases RefSeq refseq tmp

# RefSeq Releases 217
# - Viral and prokaryotic genomes of RefSeq release 217 and human genome (GRCh38.p14)
metabuli databases RefSeq217 refseq217 tmp

# GTDB 207
# - Complete Genome or Chromosome level assemblies in GTDB207 (CheckM Completeness > 90, CheckM Contamination < 5) with GTDB taxonomy.
metabuli databases GTDB207 gtdb tmp 
```


## Classification
```
metabuli classify <i:FASTA> <i:DBDIR> <o:OUTDIR> <Job ID> [options]
- FASTA : A FASTA file of reads you want to classify.
- DBDIR : The directory of reference DB. 
- OUTDIR : The directory where the result files will be generated.
- Job ID: It will be the prefix of result files.  
  
# Paired-end
metabuli classify read_1.fna read_2.fna dbdir outdir jobid

# Single-end
metabuli classify --seq-mode 1 read.fna dbdir outdir jobid

  * Options
   --threads : The number of CPU-cores used (all by default)
   --max-ram : The maximum RAM usage.
   --min-score : The minimum score to be classified (0.15 for precision mode)
   --min-sp-score : The minimum score to be classified at or below species rank. (0.5 for precision mode)
   --taxonomy-path: Directory where the taxonomy dump files are stored. (DBDIR/taxonomy by default)
   --reduced-aa : 0. Use 20 alphabets or 1. Use 15 alphabets to encode amino acids. Give the same value used for DB creation.
   --spacing-mask : Binary patterend mask for spaced k-mer. The same mask must be used for DB creation and classification. A mask should contain at least eight '1's, and '0' means skip.
```
It will generate two result files: 'Job ID_classifications.tsv' and 'Job ID_report.tsv'
#### Job ID_classifications.tsv
1. Classified or not
2. Read ID
3. Taxonomy identifier
4. Effective read length
5. DNA level identitiy score
6. Amino-acid level identity score
7. Total Hamming distance
8. Classification Rank
9. List of "taxID : k-mer match count"

```
#Example
1       read_1     2688    294     0.627551        0.806122        35      subspecies      2688:65
1       read_2     2688    294     0.816327        1       36      subspecies      2688:78
0       read_3     0       294     0       0       0       no rank
```

#### Job ID_report.tsv
Proportion of reads that are assigned to each taxon.
```
#Example
33.73   77571   77571   0       no rank unclassified
66.27   152429  132     1       no rank root
64.05   147319  2021    8034    superkingdom      d__Bacteria
22.22   51102   3       22784   phylum      p__Firmicutes
22.07   50752   361     22785   class         c__Bacilli
17.12   39382   57      123658  order           o__Bacillales
15.81   36359   3       126766  family            f__Bacillaceae
15.79   36312   26613   126767  genus               g__Bacillus
2.47    5677    4115    170517  species               s__Bacillus amyloliquefaciens
0.38    883     883     170531  subspecies                      RS_GCF_001705195.1
0.16    360     360     170523  subspecies                      RS_GCF_003868675.1
0.11    248     248     170525  subspecies                      RS_GCF_002209305.1
0.02    42      42      170529  subspecies                      RS_GCF_002173635.1
0.01    24      24      170539  subspecies                      RS_GCF_000204275.1
```

## Custom database
To build your custom database, you need three things.
1. **FASTA files** : Each sequence of your FASTA files must be separated by '>accession.version' like '>CP001849.1'
2. **accession2taxid** : Mapping from acession to taxonomy identifier. Sequences whose accessions are not listed in this file will be skipped.
3. **NCBI-style taxonomy dump** : 'names.dmp' , 'nodes.dmp', and 'merged.dmp' are required. Sequences whose taxid are not included here will be skipped.

Here, steps for creating a database based on a taxonomy of NCBI or GTDB are described.

#### 1. Prepare taxonomy and accession2taxid
  ##### NCBI taxonomy
  
  * accession2taxid can be downloaded from
  https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/
  
  * Taxonomy dump files can be downloaded from 
  https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/new_taxdump/
  
  ##### GTDB taxonomy
  
  Follow two steps below to generate GTDB taxonomy and accession2taxid file.
  * Requirements: You need assembly FASTA files whose file name (or path) include the assembly accession.  
    If you downloaded assemblies using [ncbi-genome-download](https://github.com/kblin/ncbi-genome-download), you probably don't have to care about it.  
    The regular experssion of assembly accession is (GC[AF]_[0-9].[0-9])
    
  ```
  # 1. 
  In 'util' directory
  ./prepare_gtdb_taxonomy.sh <DBDIR>
    - DBDIR : Result files are stored in 'DBDIR/taxonomy'. 
      Make sure that 'DBDIR/taxonomy' is exist and empty. 
      The same path should be used in step 1.
  ```
  It will generate taxonomy dump files and 'assacc_to_taxid.tsv' with other files.
    
  ```
  # 2. 
  ./metabuli add-to-library <FASTA list> <accession2taxid> <DBDIR> --assembly true
    - FASTA list : A list of absolute paths of each assembly files.
      Each absolute path must include assembly accession. 
    - accession2taxid : 'assacc_to_taxid.tsv' from the previous step
    - DBDIR : The same DBDIR from the previous step.
  ```
  It will add your FASTA files to 'DBDIR/library' according to their species taxonomy ID and generate 'my.accession2taxid' 

  
#### 2. Add to libarary (optional)
```
./metabuli add-to-library <FASTA list> <accession2taxid> <DBDIR>
  - FASTA list: A list of absolute paths of each FASTA files.
  - accession2taxid: A path to NCBI-style accession2taxid
  - DBDIR: The same DBDIR from the previous step.
```
This command groups your FASTA files of the same species and add stores them in separate files to DBDIR/library.  
You can skip this step in the case of
1. You have already used this command during the preparation for GTDB taxonomy.
2. Your FASTA list includes only one FASTA file per species.


#### 3. Build

```
metabuli build <DBDIR> <FASTA list> <accession2taxid> [options]
- DBDIR: The same DBDIR from the previous step. 
- FASTA list: A list of absolute paths to your FASTA files (in DBDIR/library)
- accession2taxid : accession2taxid file
  
  * Options
   --threads : The number of CPU-cores used (all by default)
   --tinfo-path : Path to prodigal training information files. (DBDIR/prodigal by default)
   --taxonomy-path: Directory where the taxonomy dump files are stored. (DBDIR/taxonomy by default)
   --reduced-aa : 0. Use 20 alphabets or 1. Use 15 alphabets to encode amino acids.
   --spacing-mask : Binary patterend mask for spaced k-mer. The same mask must be used for DB creation and classification. A mask should contain at least eight '1's, and '0' means skip.
```
It will generate **diffIdx**, **info**, **split**, and **taxID_list** and some other files. You can delete '\*\_diffIdx' and '\*\_info' if generated.
