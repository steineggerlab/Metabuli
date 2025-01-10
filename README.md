[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat)](http://bioconda.github.io/recipes/metabuli/README.html) 
![Platform](https://img.shields.io/badge/platform-Mac%20%7C%20Windows%20%7C%20Linux-brightgreen)
# Metabuli
***Metabuli*** classifies metagenomic reads by comparing them to reference genomes. You can use Metabuli to profile the taxonomic composition of your samples or to detect specific (pathogenic) species. 

***Sensitive and Specific.*** Metabuli uses a novel k-mer structure, called *metamer*, to analyze both amino acid (AA) and DNA sequences. It leverages AA conservation for sensitive homology detection and DNA mutations for specific differentiation between closely related taxa.

***A laptop is enough.*** Metabuli operates within user-specified RAM limits, allowing it to search any database that fits in storage. A PC with 8 GiB of RAM is sufficient for most analyses.

***A few clicks are enough.*** Metabuli App is now available [here](https://github.com/steineggerlab/Metabuli-App). With just a few clicks, you can run Metabuli and browse the results with Sankey and Krona plots on your PC.

***Short reads, long reads, and contigs.*** Metabuli can classify all types of sequences.


---

  
For more details, please see
[Nature Methods](https://www.nature.com/articles/s41592-024-02273-y), 
[PDF](https://www.nature.com/articles/s41592-024-02273-y.epdf?sharing_token=je_2D5Su0-xVOSjuKSAXF9RgN0jAjWel9jnR3ZoTv0M7gE7NDF_xi_3sW8QdRiwfSJNwqaXItSoeCvr7cvcoQxKLt0oROgWc6urmki9tP80cXEuHPN0D7b4y9y3i8Yv7sZw8MxxhAj7W6p9eZE2zaK3eozdOkXvwADVfso9cXIM%3D), 
[bioRxiv](https://www.biorxiv.org/content/10.1101/2023.05.31.543018v2), or [ISMB 2023 talk](https://www.youtube.com/watch?v=vz2fuRcVwyk).

Please cite: [Kim J, Steinegger M. Metabuli: sensitive and specific metagenomic classification via joint analysis of amino acid and DNA. Nature Methods (2024).](https://doi.org/10.1038/s41592-024-02273-y)

<p align="center"><img src="https://raw.githubusercontent.com/steineggerlab/Metabuli/master/.github/marv_metabuli_small.png" height="350" /></p>

---
### 🖥️  [Metabuli App](https://github.com/steineggerlab/Metabuli-App) for Windows, MacOS, and Linux are now available!
> Run taxonomic profiling in just a few clicks and explore results with Sankey and Krona plots.

> Download the app for your OS [here](https://github.com/steineggerlab/Metabuli-App/releases)—no separate Metabuli installation needed.
<p align="center"><img src="https://raw.githubusercontent.com/jaebeom-kim/Metabuli/master/.github/metabuli.jpg" height="500" /></p>


---
### Update in v1.0.9
- DB creation process improved
  - `updateDB` module to add new sequences to an existing database.
  - Users can provide CDS information to skip Prodigal's gene prediction.
  - `--max-ram` parameter added to `build` module.
  - Compatibility with taxdump files generated using [taxonkit](https://bioinf.shenwei.me/taxonkit/).
  - Please check release note for details.
  
### Update in v1.0.8
- Added `extract` module to extract reads classified into a certain taxon.

### Update in v1.0.7
- **Metabuli became faster 🚀**
  - Windows: *8.3* times faster
  - MacOS: *1.7* times faster
  - Linux: *1.3* times faster
  - Test details are in release note.
- Fixed a bug in score calculation that could affect classification results.

### Update in v1.0.6
- Windows OS is supported.
> Metabuli v1.0.6 is too slow on Windows OS. Please use v1.0.7 or later.

### Update in v1.0.4
- Fixed a minor reproducibility issue.
- Fixed a performance-harming bug occurring with sequences containing lowercased bases.
- Auto adjustment of `--match-per-kmer` parameter. Issue #20 solved.
- Record version info. in `db.parameter`

---

## Table of Contents
- [Installation](#installation)
  - [Precompiled binaries](#precompiled-binaries)
  - [Compile from source code](#compile-from-source-code)
- [Pre-built databases](#pre-built-databases)
- [Classification](#classification)
- [Extract](#extract)
- [Custom database](#custom-database)
  - [NCBI taxonomy based database](#ncbi-taxonomy-based-database)
  - [GTDB based database](#gtdb-based-database)
- [Update database](#update-database)
- [Example](#example)
  
## Installation
### Precompiled binaries
```
# install via conda
conda install -c conda-forge -c bioconda metabuli

# Linux AVX2 build (fast, recommended for most Linux system
# check using: cat /proc/cpuinfo | grep avx2)
wget https://mmseqs.com/metabuli/metabuli-linux-avx2.tar.gz; tar xvzf metabuli-linux-avx2.tar.gz; export PATH=$(pwd)/metabuli/bin/:$PATH

# Linux SSE2 build (slower, for old systems)
wget https://mmseqs.com/metabuli/metabuli-linux-sse2.tar.gz; tar xvzf metabuli-linux-sse2.tar.gz; export PATH=$(pwd)/metabuli/bin/:$PATH

# MacOS (Universal, works on Apple Silicon and Intel Macs)
wget https://mmseqs.com/metabuli/metabuli-osx-universal.tar.gz; tar xvzf metabuli-osx-universal.tar.gz; export PATH=$(pwd)/metabuli/bin/:$PATH

```
Metabuli also works on Linux ARM64 and Windows systems. Please check [https://mmseqs.com/metabuli](https://mmseqs.com/metabuli) for static builds for other architectures.

### Compile from source code
To compile Metabuli from source code use the following commands:
```
git clone https://github.com/steineggerlab/Metabuli.git
cd Metabuli
mkdir build && cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
make -j 16
```
The built binary can be found in `./build/src`.

---

## Pre-built databases
You can download [pre-built databases](https://metabuli.steineggerlab.workers.dev/) using `databases` workflow.

NOTE: The `databases` workflow may not work if you don't use the latest version of Metabuli.
In that case, please manually download databases from this [link](https://metabuli.steineggerlab.workers.dev/).


```
Usage:
metabuli databases DB_NAME OUTDIR tmp

# NOTE
- A human genome (T2T-CHM13v2.0) is included in all databases below.

1. RefSeq Virus (8.1 GiB)
- NCBI RefSeq release 223 virus genomes
- Database will be in OUT_DIR/refseq_virus
metabuli databases RefSeq_virus OUT_DIR tmp

2. RefSeq Prokaryote and Virus (115.6 GiB)
 - RefSeq prokaryote genomes (Complete Genome/Chromosome, 2024-03-26) + RefSeq Virus above.
 - Database will be in OUT_DIR/refseq_prokaryote_virus
metabuli databases RefSeq OUTDIR tmp

3. GTDB (101 GiB)
- GTDB 214.1 (Complete Genome/Chromosome, CheckM completeness > 90 and contamination < 5).
- Database will be in OUT_DIR/gtdb 
metabuli databases GTDB OUTDIR tmp

4. RefSeq Releases 224 (619 GiB)
- Viral and prokaryotic genomes of RefSeq release 224.
metabuli databases RefSeq_release OUTDIR tmp
```
Downloaded files are stored in `OUTDIR/DB_NAME` directory, which can be provided for `classify` module as `DBDIR`.

---

## Classification
```
metabuli classify <i:FASTA/Q> <i:DBDIR> <o:OUTDIR> <Job ID> [options]
- INPUT : FASTA/Q file of reads you want to classify. (gzip supported)
- DBDIR : The directory of reference DB. 
- OUTDIR : The directory where the result files will be generated.
- Job ID: It will be the prefix of result files.  
  
# Paired-end
metabuli classify read_1.fna read_2.fna dbdir outdir jobid

# Single-end
metabuli classify --seq-mode 1 read.fna dbdir outdir jobid

# Long-read 
metabuli classify --seq-mode 3 read.fna dbdir outdir jobid

  * Important parameters:
   --threads : The number of threads used (all by default)
   --max-ram : The maximum RAM usage. (128 GiB by default)
   --min-score : The minimum score to be classified 
   --min-sp-score : The minimum score to be classified at or below species rank. 
   --taxonomy-path: Directory where the taxonomy dump files are stored. (DBDIR/taxonomy by default)
   --accession-level : Set 1 to use accession level classification (0 by default). 
                       It is available when the DB is also built with accession level taxonomy.
```
- Paratemers for precision mode (Metabuli-P)
  - Illumina short reads: `--min-score 0.15 --min-sp-score 0.5`
  - PacBio HiFi reads: `--min-score 0.07 --min-sp-score 0.3`
  - PacBio Sequel II reads: `--min-score 0.005`
  - ONT reads: `--min-score 0.008`

This will generate three result files: `JobID_classifications.tsv`, `JobID_report.tsv`, and `JobID_krona.html`.
> Sankey diagram is available in the [GUI app](https://github.com/steineggerlab/Metabuli-App).

#### JobID_classifications.tsv
1. Classified or not
2. Read ID
3. Taxonomy identifier
4. Effective read length
5. DNA level identity score
6. Classification Rank
7. List of "taxID : k-mer match count"

```
1 read_1  2688  294     0.627551 subspecies  2688:65
1 read_2  2688  294     0.816327 subspecies  2688:78
0 read_3  0     294     0        no rank
```

#### JobID_report.tsv
The proportion of reads that are assigned to each taxon.
```
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

#### JobID_krona.html
It is for an interactive taxonomy report (Krona). You can use any modern web browser to open `JobID_krona.html`.
<p align="left"><img src="https://raw.githubusercontent.com/steineggerlab/Metabuli/master/.github/image.png" height="350" /></p>


### Resource requirements
Metabuli can classify reads against a database of any size as long as the database is fits in the hard disk, regardless of the machine's RAM size.
We tested it with a MacBook Air (2020, M1, 8 GiB), where we classified about 15 M paired-end 150 bp reads (~5 GiB in size) against a database built with ~23K prokaryotic genomes (~69 GiB in size).

---
## Extract 
After running the `classify` command, you can extract reads that are classified under a specific taxon.
This requires the FASTA/Q files used in the `classify` step and the `JobID_classifications.tsv` file, which is generated as one of the output files.

```
metabuli extract <i:FASTA/Q> <i:read-by-read classification> <i:DBDIR> --tax-id TAX_ID

- FASTA/Q : The FASTA/Q file(s) used during the `classify` step.
- read-by-read classification : The JobID_classifications.tsv file generated by the `classify` step.
- DBDIR : The same DBDIR used in the `classify` step.
- TAX_ID : The taxonomy ID of the taxon at any rank (e.g., species, genus) from which you want to extract the reads.


# Paired-end
metabuli extract read_1.fna read_2.fna JobID_classifications.tsv dbdir --tax-id TAX_ID

# Single-end
metabuli extract --seq-mode 1 read.fna JobID_classifications.tsv dbdir --tax-id TAX_ID

# Long-read 
metabuli extract --seq-mode 3 read.fna JobID_classifications.tsv dbdir --tax-id TAX_ID

```
#### Output
- For paired-end samples: `read_1_TAX_ID.fna` and `read_2_TAX_ID.fna`
- For single-end or long-read samples: `read_TAX_ID.fna`
  
---

## Custom database
To build a custom database, you need three things:
1. **FASTA files** : Each sequence of your FASTA files must be separated by '>accession.version' like '>CP001849.1'. The accession doesn't have to follow the NCBI format, but it must be unique and included in the accession2taxid file. 
2. **accession2taxid** : Mapping from accession to taxonomy ID. The sequences whose accessions are not listed here will be skipped.
3. **NCBI-style taxonomy dump** : 'names.dmp' , 'nodes.dmp', and 'merged.dmp' are required. The sequences whose taxonomy IDs are not included here will be skipped.
   
The procedure to build a database with NCBI or GTDB taxonomy are described below.

### User-provided CDS information (optional)
The `--cds-info` option in the `build` and `updateDB` command can be used to provide a list of absolute paths to CDS files. For the accessions included in the files, the provided CDS information will be used, and Prodigal's gene prediction will be skipped. Currently, only GenBank or RefSeq CDS files like below are supported.
```
Example:
GCA_000839185.1_ViralProj14174_cds_from_genomic.fna.gz 
```

### NCBI taxonomy based database
#### 1. Prepare taxonomy and accession2taxid
  * Download `accession2taxid` from [here](https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/).
  * Download `taxdump` files from [here](https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/new_taxdump/).
  
#### 2. Build

```
metabuli build <DBDIR> <FASTA_LIST> <accession2taxid> [options]
- DBDIR: directory where the database will be generated.
- FASTA_LIST: A file containing absolute paths to FASTA files.
- accession2taxid : A path to NCBI-style accession2taxid.
  
  * Options
   --threads : The number of threads used (all by default)
   --max-ram : The maximum RAM usage. (128 GiB by default)
   --taxonomy-path : Directory where the taxonomy dump files are stored. (DBDIR/taxonomy by default)
   --accession-level : Set 1 to creat a DB for accession level classification (0 by default).
   --cds-info : List of absolute paths to CDS files.
   --make-library : Make species library for faster execution (1 by default).
```
This will generate **diffIdx**, **info**, **split**, and **taxID_list** and some other files. You can delete '\*\_diffIdx' and '\*\_info' if generated.

### GTDB based database
>*Requirements*: 
You need assembly FASTA files whose file name (or path) includes the assembly accession like `GCF_028750015.1`.
If you downloaded assemblies using `ncbi-genome-download`, you probably don't have to care about it.
The regular expression of assembly accessions is (GC[AF]_[0-9].[0-9])
#### 1. Download GTDB taxdump files [here](https://github.com/shenwei356/gtdb-taxdump/releases).
#### 2. Build

```
GTDB_TAXDUMP: the directory where you downloaded the GTDB taxdump files.
FASTA_LIST: A file containing absolute paths of each assembly file.
DBDIR: The directory where the database will be generated.

# 1. Prepare accession2taxid and taxdump files for Metabuli
metabuli editNames <GTDB_TAXDUMP/names.dmp> <GTDB_TAXDUMP/taxid.map>
metabuli accession2taxid <FASTA_LIST> <GTDB_TAXDUMP/taxid.map>

# 2. Build a database
metabuli build <DBDIR> <FASTA_LIST> <GTDB_TAXDUMP/taxid.accession2taxid> --taxonomy-path GTDB_TAXDUMP [options]
  
  * Options are already described above. 
```
This will generate **diffIdx**, **info**, **split**, and **taxID_list** and some other files. You can delete '\*\_diffIdx' and '\*\_info' if generated.

---

## Update database 
You can add new sequences to an existing database.
The taxonomy of the previous database will be used.
You can add new taxa if the previous taxonomy does not include them (see "Add new taxa" below). 
Two files you need:
1. **FASTA files** : Each sequence must be separated by '>accession'.
2. **accession2taxid** : Mapping from accession to taxonomy ID. The sequences whose accessions are not listed here will be skipped.

```
metabuli updateDB <new DB directory> <FASTA list> <accession2taxid> <old DB directory> [options]
  - FASTA list: A file of paths to the FASTA file to be added.
  - accession2taxid : A path to NCBI-style accession2taxid.

* Options
  --threads : The number of threads used (all by default)
  --max-ram : The maximum RAM usage. (128 GiB by default)
  --accession-level : Set 1 to creat a DB for accession level classification (0 by default).
  --make-library : Make species library for faster execution (1 by default).
  --new-taxa : List of new taxa to be added.
```

### Add new taxa (optional)

#### 1. Check if you need to add new taxa
Check nodes.dmp, names.dmp, and merged.dmp files used in the previous DB creation.
You can retrieve them from `taxonomyDB` file in the previous DB directory.
```
metabuli taxdump <old DB directory/taxonomyDB>
```

#### 2. Prepare a list of new taxa
This list will be given to the `updateDB` command via the `--new-taxa` option.
The format is TSV as follows:
```
taxID parentID rank name
```
The added taxon must be linked to an existing taxon in the taxonomy tree.

## Example
> The example here was detecting SARS-CoV-2 variant-specific reads, but has changed since the pre-built DB no longer contains the variant genomes.

Classifying RNA-seq reads from a COVID-19 patient.
The whole process must take less than 10 mins using a personal machine.

#### 1. Download RefSeq Virus DB (1.5 GiB)
```
metabuli databases RefSeq_virus OUTDIR tmp
```

#### 2. Download an RNA-seq result (SRR14484345)
```
fasterq-dump --split-files SRR14484345
```
> Download SRA Toolkit containing `fasterq-dump` [here](https://github.com/ncbi/sra-tools/wiki/02.-Installing-SRA-Toolkit)
#### 3. Classify the reads using metabuli
   ```
   metabuli classify SRR14484345_1.fq SRR14484345_2.fq OUTDIR/refseq_virus RESULT_DIR JOB_ID --max-ram RAM_SIZE
   ```
#### 4. Check RESULT_DIR/JOB_ID_report.tsv
  Find a section like the example below
  ```
  92.2331 510490  442     species 694009  Severe acute respiratory syndrome-related coronavirus
  92.1433 509993  509993  no rank 2697049 Severe acute respiratory syndrome coronavirus 2
  ```
