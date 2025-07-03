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
#include <stdio.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <string.h>
#include <stdlib.h>
#include <regex.h> 
#include <zlib.h> 
#include <inttypes.h>


#include "hash.h"
#include "fastq.h"



// approx. median read length
static inline unsigned int median_rl(FASTQ_FILE* fd1,FASTQ_FILE* fd2) {
  unsigned long long ctr=0;
  unsigned int crl=1;
  unsigned long nreads=fd1->num_rds;
  
  if ( fd1->num_rds==1 && fd2==NULL) return(fd1->min_rl);
  if ( fd2!=NULL) nreads+=fd2->num_rds;
  while ( crl < MAX_READ_LENGTH ) {    
    ctr+=fd1->rdlen_ctr[crl];
    if (fd2!=NULL) ctr+=fd2->rdlen_ctr[crl];
    //printf("%d-%lu\n",crl,rdlen_ctr[crl]);
    if ( fd1->num_rds>1 && ctr>nreads/2) break;
    ++crl;
  }

  return(crl);
}

FASTQ_FILE* validate_interleaved(char *f) {
  //unsigned long cline=1;

  fprintf(stdout,"Paired-end interleaved\n");

  FASTQ_FILE* fd1=fastq_new(f,FALSE,"r");
  fastq_is_pe(fd1);
  FASTQ_ENTRY *m1=fastq_new_entry(),
    *m2=fastq_new_entry();

  char rname1[MAX_LABEL_LENGTH];
  char rname2[MAX_LABEL_LENGTH];
  unsigned long nreads1=0;
  unsigned long len=0;    

  while (!fastq_eof(fd1)) {
    // read 1
    if (fastq_read_entry(fd1,m1)==0) break;
    // read 2
    if (fastq_read_entry(fd1,m2)==0) {
      PRINT_ERROR("Error in file %s: line %lu: file truncated?",f,fd1->cline);
      exit(FASTQ_FORMAT_ERROR_EXIT_STATUS);
    }
    // match
    char *readname1=fastq_get_readname(fd1,m1,&rname1[0],&len,TRUE);
    char *readname2=fastq_get_readname(fd1,m2,&rname2[0],&len,TRUE);

    // TODO
    // replace_dots(start_pos,seq1,hdr1,hdr1_2,qual1,fdf);    
    // replace_dots(start_pos,seq2,hdr2,hdr2_2,qual2,fdf);    

    if ( strcmp(readname1,readname2) ) {
      PRINT_ERROR("Error in file %s: line %lu: unpaired read - %s",f,fd1->cline,readname1);
      exit(FASTQ_FORMAT_ERROR_EXIT_STATUS);
    } 

    if (fastq_validate_entry(fd1,m1)) {
      exit(FASTQ_FORMAT_ERROR_EXIT_STATUS);
    }
    if (fastq_validate_entry(fd1,m2)) {
      exit(FASTQ_FORMAT_ERROR_EXIT_STATUS);
    }
    PRINT_READS_PROCESSED(fd1->cline/4,100000);
    nreads1+=2;
  }
  printf("\n");
  //close_fastq(fdf); ???
  //fastq_destroy(fd1);
  return(fd1);
}

FASTQ_FILE* validate_paired_sorted_fastq_file(char *f1,char *f2) {

  FASTQ_FILE* fd1=fastq_new(f1,FALSE,"r");
  FASTQ_FILE* fd2=fastq_new(f2,FALSE,"r");
  fastq_is_pe(fd1);
  fastq_is_pe(fd2);
  FASTQ_ENTRY *m1=fastq_new_entry();
  FASTQ_ENTRY *m2=fastq_new_entry();
  char rname1[MAX_LABEL_LENGTH];
  char rname2[MAX_LABEL_LENGTH];
    
  unsigned long nreads1=0;
  unsigned long len1,len2;
  while (!fastq_eof(fd1)) {
    // read 1
    if (fastq_read_entry(fd1,m1)==0) break;

    if (fastq_validate_entry(fd1,m1)) {
      exit(FASTQ_FORMAT_ERROR_EXIT_STATUS);
    }
    
    if (fastq_read_entry(fd2,m2)==0) break;
    if (fastq_validate_entry(fd2,m2)) {
      exit(FASTQ_FORMAT_ERROR_EXIT_STATUS);
    }
    fastq_get_readname(fd1,m1,&rname1[0],&len1,TRUE);
    fastq_get_readname(fd2,m2,&rname2[0],&len2,TRUE);
    if (strcmp(rname1,rname2) ) {
      PRINT_ERROR("Readnames do not match across files (read #%ld)",fd1->cline/4+1);
      exit(FASTQ_FORMAT_ERROR_EXIT_STATUS);
    }
    PRINT_READS_PROCESSED(fd1->cline/2,100000);
    nreads1+=1;
  }
  if ( fastq_read_entry(fd1,m1)!=0) {
      PRINT_ERROR("Premature end of file2");
      exit(FASTQ_FORMAT_ERROR_EXIT_STATUS);
  }
  if ( fastq_read_entry(fd2,m2)!=0) {
      PRINT_ERROR("Premature end of file1");
      exit(FASTQ_FORMAT_ERROR_EXIT_STATUS);
  }
  printf("\n");
  return(fd1);
}


FASTQ_FILE* validate_single_fastq_file(const char *f) {

  FASTQ_FILE* fd1=fastq_new(f,FALSE,"r");
  fastq_is_pe(fd1);
  FASTQ_ENTRY *m1=fastq_new_entry();

  unsigned long nreads1=0;

  while (!fastq_eof(fd1)) {
    // read 1
    if (fastq_read_entry(fd1,m1)==0) break;

    if (fastq_validate_entry(fd1,m1)) {
      exit(FASTQ_FORMAT_ERROR_EXIT_STATUS);
    }
    PRINT_READS_PROCESSED(fd1->cline/4,100000);
    nreads1+=1;
  }
  printf("\n");
  //fastq_destroy(fd1);
  return(fd1);
}

void print_usage(int verbose_usage) {

  printf("Usage: fastq_info [-r -e -s -q -h] fastq1 [fastq2 file|pe]\n");
  if ( verbose_usage ) {
    printf(" -h  : print this help message\n");
    printf(" -s  : the reads in the two fastq files have the same ordering\n");
    printf(" -e  : do not fail with empty files\n");
    printf(" -q  : do not fail if quality encoding cannot be determined\n");
    printf(" -r  : skip check for duplicated readnames\n");
  }
}

int fastq_info_main(int argc, char **argv ) {
  //long paired=0;
  unsigned long num_reads1=0,
    num_reads2=0;

  unsigned long max_rl; // maximum read length
  unsigned long min_rl;  // minimum read length
  unsigned long min_qual; // minimum quality
  unsigned long max_qual;   // maximum quality
 
  int is_paired_data=FALSE;
  int is_interleaved=FALSE;
  int is_sorted=FALSE;
  int empty_ok=FALSE;
  int no_encoding_ok=FALSE;
  int skip_readname_check=FALSE;
  //int fix_dot=FALSE;
  
  int nopt=0;
  int c;
  opterr = 0;

  fastq_print_version();
  
  while ((c = getopt (argc, argv, "esfrhq")) != -1)
    switch (c)
      {
      case 'q':
	no_encoding_ok=TRUE;
	++nopt;
	break;
      case 'e':
	empty_ok=TRUE;
	++nopt;
	break;
      case 's':
	is_sorted=TRUE;
	++nopt;
	break;
      case 'r':
	skip_readname_check=TRUE;
	++nopt;
	break;
      case 'h':
	print_usage(TRUE);
	exit(0);
	break;
      case 'f':
        //fix_dot = TRUE;
	fprintf(stdout,"Fixing (-f) enabled: Replacing . by N (creating .fix.gz files)\n");
	PRINT_ERROR("-f option is no longer valid.");
	exit(PARAMS_ERROR_EXIT_STATUS);
	++nopt;
        break;
      default:
	++nopt;
        PRINT_ERROR("Option -%c invalid",optopt);
	exit(PARAMS_ERROR_EXIT_STATUS);
      }
  
  if (argc-nopt<2 || argc-nopt>3) {
    PRINT_ERROR("Invalid number of arguments");
    print_usage(FALSE);
    //fprintf(stderr,"%d",argc);
    exit(PARAMS_ERROR_EXIT_STATUS);
  }


  if (argc-nopt ==3) {
    is_paired_data=TRUE;
    //fprintf(stderr,"%d %d %d %s\n",argc,nopt,argc-nopt,argv[2+nopt]);
    if ( strncmp(argv[2+nopt],"pe",2) ) {
      is_interleaved=FALSE;
    } else {
      //fprintf(stderr,"Expecting interleaved reads...\n");
      is_interleaved=TRUE;
    }
  }
  
  FASTQ_FILE* fd1=NULL;
  FASTQ_FILE* fd2=NULL;
  hashtable index=NULL;
  // ************************************************************
  if ( is_interleaved ) {
    // interleaved    
    fd1=validate_interleaved(argv[1+nopt]);
    num_reads1=fd1->num_rds;
  } else if ( is_paired_data && is_sorted && skip_readname_check ) {
    fprintf(stdout,"-s option used: assuming that reads have the same ordering in both files\n");
    fd1=validate_paired_sorted_fastq_file(argv[1+nopt],argv[2+nopt]);
    num_reads1=fd1->num_rds;
    
  } else if ( !is_paired_data && skip_readname_check) {
    // SE & skip readname check
    fprintf(stdout,"Skipping check for duplicated read names\n");
    fd1=validate_single_fastq_file(argv[1+nopt]);
    num_reads1=fd1->num_rds;
  } else {
    // single or pair of fastq file(s)
    fd1=fastq_new(argv[1+nopt],FALSE,"r");
    if ( is_paired_data) fastq_is_pe(fd1);   
    fprintf(stdout,"DEFAULT_HASHSIZE=%lu\n",(long unsigned int)DEFAULT_HASHSIZE);
    index=new_hashtable(DEFAULT_HASHSIZE);
    index_mem+=sizeof(hashtable);
    fprintf(stdout,"Scanning and indexing all reads from %s\n",fd1->filename);
    fastq_index_readnames(fd1,index,0,FALSE);
    fprintf(stdout,"Scanning complete.\n");    
    num_reads1=index->n_entries;
    fprintf(stdout,"\n");
    // print some info
    fprintf(stdout, "Reads processed: %" PRIu64 "\n", index->n_entries);
    // fprintf(stdout,"Reads processed: %llu\n",index->n_entries);    
    fprintf(stdout,"Memory used in indexing: ~%ld MB\n",index_mem/1024/1024);
  }
  
  if (num_reads1 == 0 ) {
    if ( empty_ok ) {
      fprintf(stdout,"Number of reads: %lu\n",0L);
      fprintf(stdout,"Quality encoding range: %lu %lu\n",0L,0L);
      fprintf(stdout,"Quality encoding: %s\n","");
      fprintf(stdout,"Read length: %lu %lu %u\n",0L,0L,0);
      exit(0);
    }
    PRINT_ERROR("No reads found in %s.",argv[1+nopt]);
    exit(FASTQ_FORMAT_ERROR_EXIT_STATUS);
  }

  min_rl=fd1->min_rl;
  max_rl=fd1->max_rl;
  min_qual=fd1->min_qual;
  max_qual=fd1->max_qual;

  // pair-end
  if (argc-nopt ==3 && !is_interleaved && ! is_sorted ) {
    fprintf(stdout,"File %s processed\n",argv[1+nopt]);  
    fprintf(stdout,"Next file %s\n",argv[2+nopt]);  
    // validate the second file and check if all reads are paired
    fd2=fastq_new(argv[2+nopt],FALSE,"r");
    fastq_is_pe(fd2);
    
    unsigned long len;
    FASTQ_ENTRY *m2=fastq_new_entry();
    char rname[MAX_LABEL_LENGTH];
    // 
    while (!fastq_eof(fd1)) {
      // read entry
      if (fastq_read_entry(fd2,m2)==0) break;
      char *readname=fastq_get_readname(fd2,m2,&rname[0],&len,TRUE);
      INDEX_ENTRY* e=fastq_index_lookup_header(index,readname);
      if (e==NULL) {
	// complain and exit if not found
	PRINT_ERROR("Error in file %s: line %lu: unpaired read - %s",argv[2+nopt],fd2->cline,readname);
	exit(FASTQ_FORMAT_ERROR_EXIT_STATUS);
      }
      fastq_index_delete(readname,index);
      //
      if (fastq_validate_entry(fd1,m2)) {
	exit(FASTQ_FORMAT_ERROR_EXIT_STATUS);
      }
      //replace_dots(start_pos,seq,hdr,hdr2,qual,fdf);
      PRINT_READS_PROCESSED(fd2->cline/4,100000);
    }
    printf("\n");
    //fastq_destroy(fdf);//???
    if (index->n_entries>0 ) {
      PRINT_ERROR("Error in file %s: found %" PRIu64 " unpaired reads",argv[1+nopt],index->n_entries);
      // PRINT_ERROR("Error in file %s: found %llu unpaired reads",argv[1+nopt],index->n_entries);
      exit(FASTQ_FORMAT_ERROR_EXIT_STATUS);
    }
    // stats
    min_rl=min(fd2->min_rl,min_rl);
    max_rl=max(fd2->max_rl,max_rl);
    min_qual=min(fd2->min_qual,min_qual);
    max_qual=max(fd2->max_qual,max_qual);    
  }

  // stats
  // min qual/max qual/read len
  FILE* out;  
  out=stdout;

  fprintf(out,"------------------------------------\n");
  if ( num_reads2>0 ) {
    fprintf(out,"Number of reads: %lu %lu\n",num_reads1,num_reads2);
  } else {
    fprintf(out,"Number of reads: %lu\n",num_reads1);
  }

  char *enc=fastq_qualRange2enc(min_qual,max_qual);
  if ( enc == NULL && no_encoding_ok==FALSE ) {
    if (max_qual>MAX_PHRED_QUAL) {
      PRINT_ERROR("Unable to determine quality encoding - unknown range [%lu,>%u]",min_qual,MAX_PHRED_QUAL);
    } else {
      PRINT_ERROR("Unable to determine quality encoding - unknown range [%lu,%lu]",min_qual,max_qual);
    }

    exit(FASTQ_FORMAT_ERROR_EXIT_STATUS);
  }
  fprintf(out,"Quality encoding range: %lu %lu\n",min_qual,max_qual);

  if ( enc==NULL && no_encoding_ok ) {
    fprintf(out,"Quality encoding: NA\n");
  } else {
    fprintf(out,"Quality encoding: %s\n",enc);
  }
  fprintf(out,"Read length: %lu %lu %u\n",min_rl-1,max_rl-1,median_rl(fd1,fd2)-1);
  fprintf(out,"OK\n"); 
  exit(0);
}

