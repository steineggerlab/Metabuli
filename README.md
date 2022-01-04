# Metabuli
Metabuli is a taxonomical classifier using amino acid and DNA at the same time.
It is developed to achieve specificity of DNA based method and sensitivity of amino acid based method at the same time.

## Requirements (WIP)

## Installation
Installation from Github source code.
```
git clone https://github.com/JaebeomKim0731/ADclassifier
cd ADclassifier
mkdir build
cd build
cmake ..
make -j 16
```
The built binary can be found in ./build/src

## Database building
You can build DB from a directory of FASTA files or from a sinlge FASTA file with taxonomy of NCBI or GTDB
### 0. Generate taxonomy dump files & a mapping from assebmly accession to tax ID
You can choose between NCBI and GTDB
  - NCBI (WIP)
  ```
  # In 'util' directory
  ./make_assacc_to_ncbi_taxid.sh <o:outdir>
  
    - outdir : A directory where tax dump files will be generated. Make sure that the directory is exist and empty.
  ```
  
  - GTDB
  ```
  # In 'util' directory
  ./make_assacc_to_gtdb_taxid.sh <o:outdir>
  
    - outdir : A directory where tax dump files will be generated. Make sure that the directory is exist and empty.
  ```

### 1. Build DB from the directory of genome assemblies
- Requirements: The FASTA file name must include the assembly accession.
  If you downloaded assemblies using "ncbi-genome-download", you probably don't have to care about it.
  The regular experssion is (GC[AF]_[0-9]*\.[0-9]*)
```
./adclassifier build_dir <i:directory> <i:taxonomy dir> <o:output> <tmpDir> [options]

  - directory: A directory that contains all FASTA files from which you want to make a database. 
               It may contain subdirectories of FASTA file.
  - taxonomy dir: The directory where taxdump files & ass2taxID exist. 
                  The final DB will follow the taxonomy of this directory.
  - output : A directory where DB will be generated
  - tmpDir : Temporary directory
  * Option
    - --threads : number of CPU-cores used (all by default)

```

### 2. Build DB from a FASTA file. (WIP)
- Requirements


## Classification
```
./adclassifier classify <i:FASTA> <i:DB dir> <i:taxonomy dir> <o:out dir> <job ID> <tmpDir> [options]
  
  - FASTA : A FASTA file of reads you want to classify.
  - DB dir : The directory where you bulit the reference DB.
  - taxonomy dir : The directory of taxdump files. 
    !! Make sure that you are using the same taxonomy directory for DB build and classification !! 
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
