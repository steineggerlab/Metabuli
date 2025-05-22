/*
# =========================================================
# Copyright 2012-2021,  Nuno A. Fonseca (nuno dot fonseca at gmail dot com)
#
# This file is part of fastq_utils.
#
# This is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# if not, see <http://www.gnu.org/licenses/>.
#
#
# =========================================================
*/
#define VERSION "0.25.2"

#define DEFAULT  0
#define CASAVA18 1
#define INTEGERNAME 2
#define NOP 2

#ifndef MAX_READ_LENGTH
// 500000 bases should cover most of the cases for now :)
#define MAX_READ_LENGTH 2500000
#endif

#ifndef MAX_LABEL_LENGTH
#define MAX_LABEL_LENGTH 1000
#endif

#ifndef MAX_FILENAME_LENGTH
#define MAX_FILENAME_LENGTH 5000
#endif

#define MAX_BARCODE_LENGTH 50
#define MIN_READ_LENGTH 1
#define UNDEF -1
#define MAX_PHRED_QUAL 126

typedef enum  { COLORSPACE=1, SEQSPACE=0, UNDEFSPACE=-1 } READ_SPACE;
typedef enum  { TRUE=1, FALSE=0 } FASTQ_BOOLEAN;
#ifndef TRUE
#define TRUE 1
#endif
#ifndef FALSE
#define FALSE 0
#endif
typedef long FASTQ_READ_OFFSET;

#define DEFAULT_HASHSIZE 39000001

#include "hash.h"
#include <zlib.h> 
#include <stdio.h>


#define min(a,b) (a<b?a:b)
#define max(a,b) (a>b?a:b)


#define PRINT_INFO(s...) {fprintf(stderr,"INFO:"); fprintf(stderr,##s );fprintf(stderr,"\n");}
#define PRINT_ERROR(s...) {fprintf(stderr,"\nERROR: "); fprintf(stderr,##s );fprintf(stderr,"\n");}
#define FATAL_ERROR(e,s...) {fprintf(stderr,"\nERROR: "); fprintf(stderr,##s );fprintf(stderr,"\n");exit(e);}

#ifdef DEBUG
#define PRINT_DEBUG(s...) { fprintf(stderr,"DEBUG: "); fprintf(stderr,##s ); }
#else
#define PRINT_DEBUG(s...) 
#endif

#define PARAMS_ERROR_EXIT_STATUS 1
#define SYS_INT_ERROR_EXIT_STATUS 2
#define FASTQ_FORMAT_ERROR_EXIT_STATUS 3

#define PRINT_READS_PROCESSED(c,n) { if (c%n==0) { fprintf(stderr,"\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b%lu",c);fflush(stderr); }}

extern unsigned long index_mem;
extern char* encodings[];

struct index_entry {  
  // file offset: start entry
  // file offset: end entry
  // chat hdr(40)
  char *hdr;
  off_t entry_start;
  //unsigned  int  nbytes;
};
typedef struct index_entry INDEX_ENTRY;

struct fastq_entry {  
  // file offset: start entry
  // file offset: end entry
  // chat hdr(40)
  char hdr1[MAX_LABEL_LENGTH];
  char hdr2[MAX_LABEL_LENGTH];
  char seq[MAX_READ_LENGTH];
  char qual[MAX_READ_LENGTH];
  unsigned long read_len;
  long long offset;
};
typedef struct fastq_entry FASTQ_ENTRY;

typedef enum { FASTQ_GZIP, FASTQ_PLAIN } FASTQ_IO_MODE;

struct fastq_file {
  union {
    gzFile gz;
    FILE *fp;
  } fd;
  FASTQ_IO_MODE mode;

  long long cur_offset;
  unsigned long cline;
  char filename[MAX_FILENAME_LENGTH];

  unsigned long max_rl; // maximum read length
  unsigned long last_rl; // read length of the last read
  unsigned long min_rl;  // minimum read length
  unsigned long min_qual; // minimum quality
  unsigned long max_qual;   // maximum quality
  unsigned long num_rds; // number of reads
  unsigned long rdlen_ctr[MAX_READ_LENGTH]; // keep a tally on how many reads we observed per length
  
  int fix_dot;
  int fixed_dot;
  int is_pe;
  int readname_format;
  int is_casava_18;
  READ_SPACE space;
};
typedef struct fastq_file  FASTQ_FILE;

void fastq_print_version();
FASTQ_ENTRY* fastq_new_entry(void);
void fastq_write_entry(FASTQ_FILE* fd,FASTQ_ENTRY *e);

unsigned long get_elength(FASTQ_ENTRY*);
void fastq_index_delete(char *rname,hashtable index);
INDEX_ENTRY* fastq_index_lookup_header(hashtable sn_index,char *hdr);
char* fastq_get_readname(FASTQ_FILE*, FASTQ_ENTRY *,char* rn,unsigned long*,int is_header1);
int fastq_read_entry(FASTQ_FILE* fd,FASTQ_ENTRY *e);
void fastq_new_entry_stats(FASTQ_FILE *, FASTQ_ENTRY* );
int fastq_validate_entry(FASTQ_FILE *fd,FASTQ_ENTRY *e);
int fastq_read_next_entry(FASTQ_FILE* fd,FASTQ_ENTRY *e);

FASTQ_FILE* fastq_new(const char* filename, const int fix_dot, const char *mode);
void fastq_destroy(FASTQ_FILE*);
void fastq_is_pe(FASTQ_FILE* fd);
void fastq_index_readnames(FASTQ_FILE *,hashtable,long long,int);
void fastq_write_entry(FASTQ_FILE* fd,FASTQ_ENTRY *e);
void fastq_write_entry2stdout(FASTQ_ENTRY *e);
// void fastq_seek_copy_read(long offset,FASTQ_FILE* from,FASTQ_FILE *to);
char* fastq_qualRange2enc(unsigned int min_qual,unsigned int max_qual);
// void fastq_rewind(FASTQ_FILE* fd);

// void fastq_quick_copy_entry(long offset,FASTQ_FILE* from,FASTQ_FILE* to);
gzFile fastq_open(const char* filename,const char *mode);
void GZ_WRITE(gzFile fd,char *s);

static inline char* PLAIN_READ(FILE *fp, char *s, long max);
static inline char* FASTQ_READ(FASTQ_FILE *f, char *s, long max);
static inline void FASTQ_WRITE(FASTQ_FILE *f, const char *s);

static inline int fastq_eof(FASTQ_FILE *f) {
  return (f->mode == FASTQ_GZIP) ? gzeof(f->fd.gz) : feof(f->fd.fp);
};
