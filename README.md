# Metabuli
Metabuli is a taxonomical classifier using amino acid and DNA at the same time.
It is developed to achieve specificity of DNA based method and sensitivity of amino acid based method at the same time.

## Requirements (WIP)

## Installation
Installation from Github source code.
```
git clone --recursive GITHUB_LINK
cd Metabuli
mkdir build
cd build
cmake ..
make -j 16
```
The built binary can be found in ./build/src

## Database building
You can build DB from a directory of FASTA files or from a sinlge FASTA file with taxonomy of NCBI or GTDB
### 1. Prepare taxonomy and accession2taxid.map
  #### NCBI taxonomy
  
  Downlaod
  
  #### GTDB taxonomy
  
  Please follow two steps below to generate NCBI style taxonomy dump and accession2taxid.map file.
  * Requirements: The FASTA file name must include the assembly accession.  
    If you downloaded assemblies using "ncbi-genome-download", you probably don't have to care about it.  
    The regular experssion is (GC[AF]_[0-9].[0-9])
    
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
  It will add your FASTA files to 'DBDIR/library' according to their species taxonomy ID and generate 'accession2taxid.map' 

  
### 2. Add to libarary (optional)
```
./metabuli add-to-library <FASTA list> <accession2taxid> <DBDIR>
  - FASTA list: A list of absolute paths of each FASTA files.
  - accession2taxid: A path to NCBI-style accession2taxid.map
  - DBDIR: The same DBDIR from the previous step.
```
This command groups your FASTA files of the same species and add stores them in separate files to DBDIR/library.  
You can skip this step in the case of
1. You have already used this command to generate 'accession2taxid.map'.
2. Your FASTA list includes only one FASTA file per species.


### 3. Build

```
./metabuli build <DBDIR> <FASTA list> <accession2taxid> [options]
  - DBDIR: The same DBDIR from the previous step.
  - FASTA list: A list of absolute paths to your FASTA files (in DBDIR/library)
  - accession2taxid : accession2taxid.map
  * Option
    - --threads : The number of CPU-cores used (all by default)
    - --tinfo-path : Path to prodigal training information files.
    - --taxonomy-path
    - --reduced-aa : 0. Use 20 alphabets or 1. Use 15 alphabets to encode amino acids.
    - --spacing-mask : Binary patterend mask for spaced k-mer. The same mask must be used for DB creation and classification. A mask should contain at least eight '1's, and '0' means skip.
    
    

```



## Classification
```
./metabuli classify <i:FASTA> <i:DB dir> <o:out dir> <job ID> [options]
  
  - FASTA : A FASTA file of reads you want to classify.
  - DB dir : The directory where you bulit the reference DB. 
  - out dir : The directory where the report files will be generated.
  - job ID: For the result files.
  * Option
    - --threads : number of CPU-cores used (all by default)
```

### Output format
#### Read Classification
- Classified or Not
- Sequence name
- Taxonomical ID
- Length of read
- Classification score (0.0-1)
- List of "taxID : k-mer match count"

```
#Example
1       read1       1337    150     1       5334:16 1337:42
1       read2       139585  150     0.70    172034:5 167827:2 148906:3 148931:2 286486:5 139586:7
0       read3       0       150     0
```
#### Composition Report
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
