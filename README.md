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

---
### 🖥️  [Metabuli App](https://github.com/steineggerlab/Metabuli-App) for Windows, MacOS, and Linux are now available!
> Run taxonomic profiling in just a few clicks and explore results with Sankey and Krona plots.

> Download the app for your OS [here](https://github.com/steineggerlab/Metabuli-App/releases)—no separate Metabuli installation needed.
<p align="center"><img src="https://raw.githubusercontent.com/jaebeom-kim/Metabuli/master/.github/metabuli.jpg" style="max-height: 500px; width: auto;" /></p>


---
### Update in v1.1.0
- Fix errors in v1.0.9
- Custom DB creation became easier
- Improve `updateDB` command
  
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


---

## Table of Contents
- [Installation](#installation)
  - [Precompiled binaries](#precompiled-binaries)
  - [Compile from source code](#compile-from-source-code)
- [Pre-built databases](#pre-built-databases)
- [Classification](#classification)
- [Extract](#extract)
- [GTDB-based custom database](#gtdb-based-custom-database)
  - [Creat a new database](#creat-a-new-database)
  - [Add new sequences to an existing database](#add-new-sequences-to-an-existing-database)
- [NCBI or custom taxonomy based database](#ncbi-or-custom-taxonomy-based-database)
  - [Creat a new database](#creat-a-new-database-1)
  - [Add new sequences to an existing database](#add-new-sequences-to-an-existing-database-1)
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
```
git clone --recurse-submodules https://github.com/steineggerlab/Metabuli.git
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
You can use `--lineage 1` option in `classify` module to print the full lineage next to `rank` column. 
1. `is_classified`: Classified or not
2. `name`: Read ID
3. `taxID`: Tax. ID in the tax. dump files used in database creation
4. `query_length`: Effective read length
5. `score`: DNA level identity score
6. `rank`: Taxonomic rank of the taxon
7. `taxID:match_count`: List of "taxID : k-mer match count"

```
#is_classified name  taxID query_length score rank taxID:match_count
1 read_1  2688  294     0.627551 subspecies  2688:65
1 read_2  2688  294     0.816327 subspecies  2688:78
0 read_3  0     294     0        no rank
```

#### JobID_report.tsv
It follows Kraken2's report format. The first line is a header, and the rest of the lines are tab-separated values. The columns are as follows:

1. `clade_proportion`: Percentage of reads classified to the clade rooted at this taxon
2. `clade_count`: Number of reads classified to the clade rooted at this taxon
3. `taxon_count`: Number of reads classified directly to this taxon
4. `rank`: Taxonomic rank of the taxon
5. `taxID`: Tax ID according to the taxonomy dump files used in the database creation
6. `name`: Taxonomic name of the taxon

```
#clade_proportion  clade_count  taxon_count rank  taxID name  
33.73   77571   77571   no rank 0        unclassified
66.27   152429  132     no rank 1        root
64.05   147319  2021    superkingdom  8034          d__Bacteria
22.22   51102   3       phylum        22784         p__Firmicutes
22.07   50752   361     class         22785            c__Bacilli
17.12   39382   57      order         123658             o__Bacillales
15.81   36359   3       family        126766              f__Bacillaceae
15.79   36312   26613   genus         126767                 g__Bacillus
2.47    5677    4115    species       170517                 s__Bacillus amyloliquefaciens
0.38    883     883     subspecies    170531                        RS_GCF_001705195.1
0.16    360     360     subspecies    170523                        RS_GCF_003868675.1

```

#### JobID_krona.html
It is for an interactive taxonomy report (Krona). You can use any modern web browser to open `JobID_krona.html`.
<p align="left"><img src="https://raw.githubusercontent.com/steineggerlab/Metabuli/master/.github/image.png" style="max-height: 350px; width: auto;" /></p>


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

## GTDB-based custom database
> ***User-provided CDS information (optional):***
Use `--cds-info` to provide absolute paths to CDS files. For included accessions, the provided CDS is used, and Prodigal's gene prediction is skipped. Only GenBank/RefSeq CDS files are supported (e.g., GCA_000839185.1_ViralProj14174_cds_from_genomic.fna.gz).

### <u>***Creat a new database***</u>
>[!IMPORTANT] 
>***Requirements***: Reference FASTA file name (or path) must include the assembly accession (e.g., `GCF_028750015.1`, regex`GC[AF]_[0-9]+\.[0-9]+`). Files from RefSeq or GenBank meet this requirement.

#### 1. Download taxonkit-generated GTDB taxdump files [here](https://github.com/shenwei356/gtdb-taxdump/releases).
#### 2. Build

```
# GTDB_TAXDUMP: Directory with downloaded GTDB taxdump files.
# FASTA_LIST: File of reference genome absolute paths.
# DBDIR: Directory where the database will be generated.

metabuli build --gtdb 1 <DBDIR> <FASTA_LIST> <GTDB_TAXDUMP/taxid.map> --taxonomy-path <GTDB_TAXDUMP>  [options]

* Options
   --threads : The number of threads to utilize (all by default)
   --max-ram : The maximum RAM usage. (128 GiB by default)
   --accession-level : Set 1 to creat a DB for accession level classification (0 by default).
   --cds-info : List of absolute paths to CDS files.
  
```
This will generate **diffIdx**, **info**, **split**, and **taxID_list** and some other files. You can delete `*_diffIdx` and `*_info` files.

### <u>***Add new sequences to an existing database***</u>
> [!NOTE] 
> If you want to use new GTDB release, please build a new database from scratch.

You can add new sequences to a GTDB-based database. Expanding the taxonomy for virus or eukaryote is also possible.

#### \<Add GTDB genomes>
```
# GTDB_TAXDUMP: Directory with downloaded GTDB taxdump files.
# FASTA_LIST: File of absolute paths to new sequences.
# NEW DBDIR: Updated database is generated here.
# OLD DBDIR: Directory of an existing database.

metabuli updateDB --gtdb 1 <NEW DBDIR> <FASTA_LIST> <GTDB_TAXDUMP/taxid.map> <OLD DBDIR> [options]

* Options 
  --make-library: When many species are in the same FASTA, enable it for faster execution (0 by default).
  --new-taxa: List of new taxa to be added.
  --threads: The number of threads to utilize (all by default)
  --max-ram: The maximum RAM usage. (128 GiB by default)
  --accession-level: Set 1 to add new sequences for accession level classification (0 by default).
  --cds-info: List of absolute paths to CDS files.
```

#### \<Add sequences of new taxa>
> [!WARNING] 
> Mixing taxonomies within the same domain is not recommended. For example, adding prokaryotes to a GTDB database using NCBI taxonomy will cause issues, but adding eukaryotes or viruses using NCBI taxonomy is fine since GTDB does not cover them.

1\. **Check taxdump** files to see if you need to add new taxa. `taxdump` command retrieves taxdump files of an existing database.

2-1\. **Create a new taxa list** 
  
If you have **accession2taxid** and **taxonomy dump** files of the new sequences, you can use `createnewtaxalist` to create an input for `--new-taxa` option. If not, you have to prepare the input manually (see below).
```
metabuli createnewtaxalist <OLD DBDIR> <FASTA_LIST> <new taxonomy dump> <accession2taxid> <OUTDIR>
```
  It generates `newtaxa.tsv` for `--new-taxa` option and `newtaxa.accession2taxid`.

##### Example
Suppose you're adding eukaryotes to a GTDB database. As GTDB doesn't include eukaryotes, you may want to use NCBI taxonomy for eukaryotes.
You can download `taxdump` files from [here](https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/new_taxdump/) and `accession2taxid` from [here](https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/).
```
metabuli createnewtaxalist <GTDB dir> <new seq list> <NCBI taxdump dir> <NCBI accession2taxid> <out dir>
metabuli updateDB <new db dir> <new seq list> <out dir/newtaxa.accession2taxid> <GTDB dir> --new-taxa <out dir/newtaxa.tsv>
```

2-2\. **Manually prepare a new taxa list**

For the `--new-tax` option, provide a four-column TSV file in the following format.
```
taxID parentID rank name
```
The new taxon must be linked to a taxon in the existing database's taxonomy.

##### Example
Suppose you want to add *Saccharomyces cerevisiae* to a GTDB database. 
After inspecting taxonomy with `taxdump`, you find that the taxonomy lacks the Fungi kingdom and only includes one eukaryote (*Homo sapiens*). In this scenario, your new taxa list and accession2taxid should be as follows.
```
# New taxa list
## taxid  parentTaxID rank  name // Don't put this header in your actual file.
10000013	10000012	species	Saccharomyces cerevisiae
10000012	10000011	genus	Saccharomyces
10000011	10000010	family	Saccharomycetaceae
10000010	10000009	order	Saccharomycetales
10000009	10000008	class	Saccharomycetes
10000008	10000007	phylum	Ascomycota
10000007	10000000	kingdom	Fungi // 10000000 is Eukaroyte taxID of the pre-built DB.

# accession2taxid
accession accession.version taxid gi
newseq1 newseq1 10000013  0
newseq2 newseq2 10000013  0
```

## NCBI or custom taxonomy based database
### <u>***Creat a new database***</u>

>[!IMPORTANT] 
Three requirements:
> 1. **FASTA files** : Each sequence must have a unique `>accession.version` or `>accesion` header (e.g., `>CP001849.1` or `>CP001849`).
> 2. **NCBI-style accession2taxid** : Sequences with accessions absent here are skipped, and versions are ignored.
> 3. **NCBI-style taxonomy dump** : `names.dmp`, `nodes.dmp`, and `merged.dmp`. Sequences with tax. IDs absent here are skipped.

#### 1. Prepare NCBI-format taxonomy dump files and accession2taxid
* Download `accession2taxid` from [here](https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/).
* Download `taxdump` files from [here](https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/new_taxdump/).
* For custom sequences, edit `accession2taxid` and `taxdump` files as follows.
  * `accession2taxid`
    * For a sequence whose header is `>custom`, add `custom[tab]custom[tab]taxid[tab]anynumber`.
    * As above, version number is not necessary.
    * `taxid` must be included in the `nodes.dmp` and `names.dmp`.
    * Put any number for the last column. It is not used in Metabuli.
  * `taxdump`
    * Edit `nodes.dmp` and `names.dmp` if you introduced a new `taxid` in `accession2taxid`.

#### 2. Build
```
# DBDIR: Directory where the database will be generated.
# FASTA_LIST: A file containing absolute paths to FASTA files.
# accession2taxid : NCBI-style accession2taxid file.
# TAXDUMP: Directory with taxonomy dump files.

metabuli build <DBDIR> <FASTA_LIST> <accession2taxid> --taxonomy-path <TAXDUMP> [options]

* Options
  --make-library: When many species are in the same FASTA, enable it for faster execution (0 by default).
  --threads: The number of threads to use (all by default)
  --max-ram: The maximum RAM usage. (128 GiB by default)
  --accession-level: Set 1 to creat a DB for accession level classification (0 by default).
  --cds-info: List of absolute paths to CDS files.
```
This will generate **diffIdx**, **info**, **split**, and **taxID_list** and some other files. You can delete `*_diffIdx` and `*_info` files and `DATE-TIME` folder (e.g., `2025-1-24-10-32`) if generated.

### <u>***Add new sequences to an existing database***</u>
You can add new sequences to an existing database, of which taxonomy will be used. You can add new taxa if the previous taxonomy does not include them (see "Add sequences of new taxa" below).

#### \<Add sequences of existing taxa>

1\. **Prepare two files**
  - **New FASTA file list** : Each sequence must have a unique `>accession.version` or `>accesion` header.
  - **NCBI-style accession2taxid** : Sequences with accessions absent here are skipped. Put any number in the GI column. Version number is ignored.
    ```
    accession  accession.version  taxID  gi 
    SequenceA  SequenceA.1 960611  0
    SequenceB  SequenceB.1 960612  0
    NoVersionOkay  NoVersionOkay 960613  0
    ```

2\. **Update database**
```
# NEW DBDIR: Directory where the updated database will be generated.
# FASTA_LIST: A file of paths to new FASTA files.
# accession2taxid : NCBI-style accession2taxid file.
# OLD DBDIR: Directory of an existing database.

metabuli updateDB <NEW DBDIR> <FASTA_LIST> <accession2taxid> <OLD DBDIR> [options]
  - FASTA list: A file of paths to the FASTA file to be added.
  - accession2taxid : A path to NCBI-style accession2taxid.

* Options
  --threads : The number of threads used (all by default)
  --max-ram : The maximum RAM usage. (128 GiB by default)
  --accession-level : Set 1 to create a DB for accession level classification (0 by default).
  --make-library : Make species library for faster execution (1 by default).
  --new-taxa : List of new taxa to be added.
```

#### \<Add sequences of new taxa> - Please refer [this section](#add-sequences-of-new-taxa).


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

## Reference
Shen, W., Ren, H., TaxonKit: a practical and efficient NCBI Taxonomy toolkit, Journal of Genetics and Genomics, https://doi.org/10.1016/j.jgg.2021.03.006