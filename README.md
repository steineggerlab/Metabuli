# Metabuli
Metabuli is a taxonomical classifier using amino acid and DNA at the same time.

## Requirements

## Installation

## Database building
You can build DB from a directory of FASTA files or from a sinlge FASTA file with taxonomy of NCBI or GTDB
### 0. Generate taxonomy dump files & a mapping from assebmly accession to tax ID(ass2taxID)
You can choose between NCBI and GTDB
  - NCBI
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
```
./adclassifier build_dir <i:directory> <i:taxonomy dir> <o:output> <tmpDir> [options]
  - directory: A directory that contains all FASTA files from which you want to make a database. 
                               It may contain subdirectories of FASTA file.
  - taxonomy dir: The directory where taxdump files and ass2taxID exist. 
                   The final DB will follow the taxonomy of this directory.
  - output : A directory where DB will be generated
  - tmpDir : Temporary directory
  * Option
    - --threads : number of CPU-cores used (all by default)
    - --tax-mode
```

### 2. Build DB from a FASTA file. (WIP)


## Classification
```
./adclassifier classify <i:FASTA> <i:DB dir> <i:taxonomy dir> <out dir> <tmpDir>[options]
  - FASTA : A FASTA file of reads you want to classify.
  - DB dir : The directory where you bulit the reference DB.
  - taxonomy dir : The directory of taxdump files. 
    !! Make sure that you are using the same taxonomy directory for DB build and classification !! 
  - out dir : The directory where the report files will be generated.
```

## Output format

