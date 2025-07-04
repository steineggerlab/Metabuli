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
#include "fastq.h"
#include <stdio.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <string.h>
#include <stdlib.h>
#include <regex.h> 
#include <zlib.h> 

// Macros
//static char read_buffer[MAX_READ_LENGTH+1];

// public
unsigned long index_mem=0;
char* encodings[]={"33","64","solexa","33 *","sanger"};

#define READ_LINE(fd) gzgets(fd,&read_buffer[0],MAX_READ_LENGTH)


FASTQ_FILE* new_fastq_file(const char* filename,const int);
// functions
static int is_casava_1_8_readname(const char *s);
static int is_int_readname(const char *s);
static int is_nosuffix_readname(const char *s);
static READ_SPACE is_color_space(char *seq,FASTQ_FILE* f);

void free_indexentry(INDEX_ENTRY *e);
INDEX_ENTRY* new_indexentry(hashtable ht,char*hdr,int len,long start_pos);

static inline int compare_headers(const char *hdr1,const char *hdr2); //?


static inline char* GZ_READ(gzFile fd,char *s,long max);
//void GZ_WRITE(gzFile fd,char *s);

//
gzFile fastq_open(const char* filename,const char *mode);
static void fastq_close(gzFile fd);

void fastq_print_version() {
  fprintf(stdout,"fastq_utils %s\n",FASTQ_UTIL_VERSION);
}


/* ******************************************************************************* */
FASTQ_ENTRY* tmp_entry=NULL;
static FASTQ_ENTRY* get_tmp_entry() {
  if (tmp_entry!=NULL) return tmp_entry;
  tmp_entry=fastq_new_entry();
  return tmp_entry;
}

// void fastq_rewind(FASTQ_FILE* fd) {
//   fd->cline=1;
//   gzrewind(fd->fd);
// }
void fastq_write_entry2stdout(FASTQ_ENTRY *e) {
  fprintf(stdout,"%s",e->hdr1);
  fprintf(stdout,"%s",e->seq);
  fprintf(stdout,"%s",e->hdr2);
  fprintf(stdout,"%s",e->qual);
}

void fastq_is_pe(FASTQ_FILE* fd) {
  fd->is_pe=TRUE;
}

unsigned long get_elength(FASTQ_ENTRY* m) {
  // -2 = \n\0
  return(m->read_len-2);
}

void fastq_new_entry_stats(FASTQ_FILE *fd, FASTQ_ENTRY* entry) {

  unsigned long slen=entry->read_len;
  if (slen<fd->min_rl) {
    fd->min_rl=slen;
  }
  if (slen>fd->max_rl) {
    fd->max_rl=slen;
  }
  ++fd->num_rds;
  fd->last_rl=slen;
  fd->rdlen_ctr[slen]++;
  // update min/max quality
}

FASTQ_ENTRY* fastq_new_entry(void) {

  FASTQ_ENTRY* newEntry= (FASTQ_ENTRY*)malloc(sizeof(FASTQ_ENTRY));
  if (newEntry==NULL) {
    PRINT_ERROR("unable to allocate %ld bytes of memory",sizeof(FASTQ_ENTRY));
    exit(SYS_INT_ERROR_EXIT_STATUS);
  }
  newEntry->read_len=0;
  newEntry->offset=0;
  return newEntry;
}
unsigned long ctr_seek=0,ctr_noseek=0;
// void fastq_quick_copy_entry(long offset,FASTQ_FILE* from,FASTQ_FILE* to) {

//   if ( gztell(from->fd)!=offset ) {
//     //fprintf(stderr,"miss %lu / %lu\n",offset, gztell(from->fd));
//     // we need to seek
//     if (gzseek(from->fd,offset,SEEK_SET)<0) {
//       PRINT_ERROR("Error in file %s: line %lu: gzseek failed",from->filename,from->cline);
//       exit(SYS_INT_ERROR_EXIT_STATUS);
//     }
//     ++ctr_seek;
//   } else     ++ctr_noseek;
//   fprintf(stderr,"%lu / %lu\n",ctr_seek, ctr_noseek);
//   FASTQ_ENTRY* e=get_tmp_entry();
//   if( gzeof(from->fd)) {
//     PRINT_ERROR("Error in file %s: line %lu: premature eof",from->filename,from->cline);
//     exit(FASTQ_FORMAT_ERROR_EXIT_STATUS);
//   }
//   GZ_READ(from->fd,&e->hdr1[0],MAX_LABEL_LENGTH);
//   if ( e->hdr1[0]=='\0') {
//     PRINT_ERROR("Error in file %s: line %lu: file truncated",from->filename,from->cline);
//     exit(FASTQ_FORMAT_ERROR_EXIT_STATUS);
//   }
//   GZ_READ(from->fd,&e->seq[0],MAX_READ_LENGTH);
//   GZ_READ(from->fd,&e->hdr2[0],MAX_LABEL_LENGTH);
//   GZ_READ(from->fd,&e->qual[0],MAX_READ_LENGTH);
//   if (e->seq[0]=='\0' || e->hdr2[0]=='\0' || e->qual[0]=='\0' ) {
//     PRINT_ERROR("Error in file %s: line %lu: file truncated",from->filename,from->cline);

//     exit(FASTQ_FORMAT_ERROR_EXIT_STATUS);
//   }
//   GZ_WRITE(to->fd,&e->hdr1[0]);
//   GZ_WRITE(to->fd,&e->seq[0]);
//   GZ_WRITE(to->fd,&e->hdr2[0]);
//   GZ_WRITE(to->fd,&e->qual[0]);
//   from->cur_offset=gztell(from->fd);
// }
/*
 *
 */
FASTQ_FILE* fastq_new(const char* filename, const int fix_dot,const char *mode) {

  FASTQ_FILE* f = (FASTQ_FILE*)malloc(sizeof(FASTQ_FILE));
  if (!f) FATAL_ERROR(SYS_INT_ERROR_EXIT_STATUS, "fastq_new: malloc failed");
  memset(f, 0, sizeof(*f));

  strncpy(f->filename, filename, MAX_FILENAME_LENGTH - 1);
  f->filename[MAX_FILENAME_LENGTH - 1] = '\0';

  f->fix_dot = fix_dot;
  f->min_rl = MAX_READ_LENGTH;
  f->min_qual = MAX_PHRED_QUAL;
  f->space = UNDEFSPACE;
  if (strlen(filename) >= 3 && strcmp(filename + strlen(filename) - 3, ".gz") == 0) {
    f->mode = FASTQ_GZIP;
    f->fd.gz = gzopen(filename, mode);
    if (!f->fd.gz) FATAL_ERROR(PARAMS_ERROR_EXIT_STATUS, "Unable to open %s", filename);
    gzbuffer(f->fd.gz, 128000);
  } else {
    f->mode = FASTQ_PLAIN;
    f->fd.fp = fopen(filename, mode);
    if (!f->fd.fp) FATAL_ERROR(PARAMS_ERROR_EXIT_STATUS, "Unable to open %s", filename);
  }

  return f;
  
//   FASTQ_FILE* new=(FASTQ_FILE*)malloc(sizeof(FASTQ_FILE));
//   if (new==NULL) {
//     PRINT_ERROR("fastq_new: Error while processing file %s: unable to allocate %ld bytes of memory",filename,sizeof(FASTQ_FILE));
//     exit(SYS_INT_ERROR_EXIT_STATUS);
//   }
//   new->cur_offset=0L;
//   new->max_rl=0L;
//   new->last_rl=0L;
//   new->min_rl=MAX_READ_LENGTH;
//   new->min_qual=MAX_PHRED_QUAL;
//   new->max_qual=0;
//   new->num_rds=0;
//   new->fix_dot=fix_dot;
//   new->fixed_dot=FALSE;
//   new->is_pe=FALSE;
//   new->readname_format=UNDEF;
//   new->is_casava_18=UNDEF;
//   new->space=UNDEFSPACE;
//   strncpy(new->filename,filename,MAX_FILENAME_LENGTH-1);
//   new->filename[MAX_FILENAME_LENGTH-1]='\0';
//   new->fd=fastq_open(filename,mode);
//   memset(new->rdlen_ctr,0,sizeof(long)*MAX_READ_LENGTH);
//   return(new);
}


// void fastq_seek_copy_read(long offset,FASTQ_FILE* from,FASTQ_FILE* to) {
//   if (gzseek(from->fd,offset,SEEK_SET)<0) {
//     PRINT_ERROR("Error in file %s: line %lu: gzseek failed",from->filename,from->cline);    
//     exit(SYS_INT_ERROR_EXIT_STATUS);
//   }
//   FASTQ_ENTRY* e=get_tmp_entry();
//   fastq_read_entry(from,e);
//   fastq_write_entry(to,e);
// }


static inline char* GZ_READ(gzFile fd,char *s,long max) {
  char *s2=gzgets(fd,s,max);
  if ( s2!=Z_NULL ) return(s2);
  s[0]='\0';
  return(NULL);
  //fprintf(stderr,"ERROR: read from file\n");
  //exit(1);
}

void GZ_WRITE(gzFile fd,char *s) {
  int n=gzputs(fd,s);
  if ( n>0 ) return;
  if ( *s=='\0' ) return;
  const char *errmsg=gzerror(fd,&n);
  PRINT_ERROR("%s.\n",errmsg);
  /* switch(n) { */
  /* case Z_ERRNO: */
  /*   fprintf(stderr,"Error: Internal error (%d).\n",errno()); */
  /*   break; */
  /* case Z_STREAM_ERROR: */
  /*   fprintf(stderr,"Error: The stream is invalid, is not open for writing, or is in an invalid state.\n"); */
  /*   break; */
  /* case Z_BUF_ERROR: */
  /*   fprintf(stderr,"Error: internal error.\n"); */
  /*   break; */
  /* case Z_MEM_ERROR: */
  /*   fprint(stderr,"Error: Insufficient memory available to compress.\n"); */
  /*   break; */
  /* } */
  exit(SYS_INT_ERROR_EXIT_STATUS);
}
/*
 * Reads one entry from the fastq file fd and places it in e
 * Returns 0 on failure, 1 on success
 */
int fastq_read_next_entry(FASTQ_FILE* fd,FASTQ_ENTRY *e) {

  int r=fastq_read_entry(fd,e);
  if ( r <=0 ) return r;
  fastq_new_entry_stats(fd,e);
  return(1);
}
/* read the next entry e from the fastq stream fd */
int fastq_read_entry(FASTQ_FILE* fd,FASTQ_ENTRY *e) {
  if (fd->mode == FASTQ_GZIP && gzeof(fd->fd.gz)) return 0;
  if (fd->mode == FASTQ_PLAIN && feof(fd->fd.fp)) return 0;

  if (!FASTQ_READ(fd, e->hdr1, MAX_LABEL_LENGTH)) return 0;
  if (!FASTQ_READ(fd, e->seq, MAX_READ_LENGTH)) return 0;
  if (!FASTQ_READ(fd, e->hdr2, MAX_LABEL_LENGTH)) return 0;
  if (!FASTQ_READ(fd, e->qual, MAX_READ_LENGTH)) return 0;

  fd->cline += 4;
  e->read_len = strlen(e->seq);
  return 1;
}


/* read the next entry e from the fastq stream fd */
void fastq_write_entry(FASTQ_FILE* fd,FASTQ_ENTRY *e) {
  FASTQ_WRITE(fd, e->hdr1);
  FASTQ_WRITE(fd, e->seq);
  FASTQ_WRITE(fd, e->hdr2);
  FASTQ_WRITE(fd, e->qual); 
}


char* fastq_qualRange2enc(unsigned int min_qual,unsigned int max_qual) {
  int enc=0;
  if ( min_qual>=33 && min_qual <59 && max_qual>=90 ) {
    enc=4; // sanger: used by ONT and possibly by pacbio
  } else if ( min_qual >=33 && max_qual <=73 ) {
    enc=0; // 33
  } else if ( min_qual <59 ) {
    enc=0; // 33
  } else if ( min_qual >=64 && max_qual>74 ) {
    enc=1; // 64
  } else if (min_qual >=59 && max_qual>74 ) { // min_qual<64
    enc=2; // solexa
  } else {
    enc=3; // 33 is the default value (* means that the default value was used)
  }
  if ( max_qual > MAX_PHRED_QUAL )
    return(NULL);
  // raw reads should not have a value greater than min_qual+60
  // higher scores are possible in assemblies or read maps (http://en.wikipedia.org/wiki/FASTQ_format)
  if ( enc!=4 && max_qual > min_qual+60  ) {
    return(NULL);
  }
  return encodings[enc];
}

// return 0 on sucess, 1 otherwise
inline int fastq_validate_entry(FASTQ_FILE* fd,FASTQ_ENTRY *e) {
//(char *hdr,char *hdr2,char *seq,char *qual,unsigned long linenum,const char* filename) {
  char rname1[MAX_LABEL_LENGTH];
  char rname2[MAX_LABEL_LENGTH];

  // Sequence identifier
  if ( e->hdr1[0]!='@' ) {
    PRINT_ERROR("Error in file %s: line %lu: sequence identifier should start with an @ - %s",fd->filename,fd->cline,e->hdr1);
    return 1;
  }  
  if ( e->hdr1[1]=='\0' || e->hdr1[1]=='\n' || e->hdr1[1]=='\r') {
    PRINT_ERROR("Error in file %s: line %lu: sequence identifier should be longer than 1",fd->filename,fd->cline);	   
    return 1;
  }
  // sequence
  unsigned long slen=0;
  short found_T=FALSE, found_U=FALSE;
  while ( e->seq[slen]!='\0' && e->seq[slen]!='\n' && e->seq[slen]!='\r' ) {
    // check content: ACGT acgt nN 0123....include the .?
    if ( e->seq[slen]!='A' && e->seq[slen]!='C' && e->seq[slen]!='G' && e->seq[slen]!='T' && e->seq[slen]!='U' &&
	 e->seq[slen]!='a' && e->seq[slen]!='c' && e->seq[slen]!='g' && e->seq[slen]!='t' && e->seq[slen]!='u' &&
	 e->seq[slen]!='0' && e->seq[slen]!='1' && e->seq[slen]!='2' && e->seq[slen]!='3' &&
	 e->seq[slen]!='n' && e->seq[slen]!='N' && e->seq[slen]!='.' ) {      
      PRINT_ERROR("Error in file %s: line %lu: invalid character '%c' (hex. code:'%x'), expected ACGTUacgtu0123nN.",fd->filename,fd->cline+1,e->seq[slen],e->seq[slen]);
      return 1;
    }
    // soft check - this should probably be enforced on all reads in the file
    if ( e->seq[slen]=='U' || e->seq[slen]=='u') {
      found_U=TRUE;
      if (found_T) {
	      PRINT_ERROR("Error in file %s: line %lu: read contains both U and T bases",fd->filename,fd->cline-2);
	      return 1;
      }
    } else {
      if ( e->seq[slen]=='T' || e->seq[slen]=='t') {
	      found_T=TRUE;
	      if (found_U) {
	        PRINT_ERROR("Error in file %s: line %lu: read contains both U and T bases",fd->filename,fd->cline-2);
	        return 1;
	      }
      }	
    }
    slen++;
  }  
  fastq_new_entry_stats(fd,e);  
  // check len
  if (slen < MIN_READ_LENGTH ) {
    PRINT_ERROR("Error in file %s: line %lu: read length too small - %lu",fd->filename,fd->cline+1,slen);
    return 1;
  }
  // be tolerant
  //if (hdr2[1]!='\0' && hdr2[1]!='\n' && hdr2[1]!='\r') {
  //  fprintf(stderr,"Error in file %s, line %lu:  header2 wrong. The line should contain only '+' followed by a newline.\n",filename,linenum+2);
  //  return 1;
  //}  
  if (e->hdr2[0]!='+') {
    PRINT_ERROR("Error in file %s: line %lu:  header2 wrong. The line should contain only '+' followed by a newline or read name (header1).",fd->filename,fd->cline+2);
    return 1;
  }
  // length of hdr2 should be 1 or be the same has the hdr1
  // ignore the + sign
  //get_readname(&hdr2[1]);
  unsigned long len;
  if (e->hdr2[0]!='\0' && e->hdr2[0]!='\r' ) {    
    char *rn1=fastq_get_readname(fd,e,&rname1[0],&len,TRUE);
    char *rn2=fastq_get_readname(fd,e,&rname2[0],&len,FALSE);
    if ( !compare_headers(rn1,rn2) ) {
      PRINT_ERROR("Error in file %s: line %lu:  header2 differs from header1\nheader 1 \"%s\"\nheader 2 \"%s\"",fd->filename,fd->cline,e->hdr1,e->hdr2);
      return 1;
    }
  }
  // qual length==slen
  unsigned long qlen=0;
  while ( e->qual[qlen]!='\0' && e->qual[qlen]!='\n' && e->qual[qlen]!='\r') {
    unsigned int x=(unsigned int)e->qual[qlen];
    if (x<fd->min_qual) { fd->min_qual=x; }
    if (x>fd->max_qual) { fd->max_qual=x; }
    qlen++;    
  }  

  if ( fd->space==SEQSPACE && qlen!=slen ) {
    PRINT_ERROR("Error in file %s: line %lu: sequence and quality don't have the same length %lu!=%lu",fd->filename,fd->cline,slen,qlen);
    return 1;
  }
  
  if ( fd->space==COLORSPACE &&  qlen==(slen-1) ) return(0);
  if ( fd->space==COLORSPACE &&  qlen==slen ) return(0);
  if ( fd->space==COLORSPACE ) {
    PRINT_ERROR("Error in file %s: line %lu: sequence and quality length don't match %lu!=%lu",fd->filename,fd->cline,slen,qlen);
    return 1;
  }
  return 0;
}


// add option to replace dots
void fastq_index_readnames(FASTQ_FILE* fd1,hashtable index,long long start_offset,int replace_dots) {

  // replace dots not used anymore
  fd1->fix_dot=replace_dots;
  FASTQ_ENTRY *m1=fastq_new_entry();
  char rname[MAX_LABEL_LENGTH];

  if ((fd1->mode == FASTQ_GZIP && fd1->fd.gz == NULL) ||
      (fd1->mode == FASTQ_PLAIN && fd1->fd.fp == NULL)) {
    PRINT_ERROR("Unable to open %s",fd1->filename);
    exit(PARAMS_ERROR_EXIT_STATUS);
  }
  // move to the right position
  if(start_offset>0) {
    PRINT_ERROR(" Not implemented");
    exit(SYS_INT_ERROR_EXIT_STATUS);
  }
  unsigned long len;
  // index creation could be done in parallel...
  while (!fastq_eof(fd1)) {
    if ( fastq_read_next_entry(fd1,m1)==0) break;

    char* readname=fastq_get_readname(fd1,m1,&rname[0],&len,TRUE);
    //fprintf(stderr,"4---%s\n",readname);fflush(stderr);
    // TODO: replace dots() -> needs a new file
    //    replace_dots(start_pos,seq,hdr,hdr2,qual,fdf);    
    // check for duplicates
    if ( fastq_index_lookup_header(index,readname)!=NULL ) {
      PRINT_ERROR("Error in file %s: line %lu: duplicated sequence %s",fd1->filename,fd1->cline,readname);
      exit(FASTQ_FORMAT_ERROR_EXIT_STATUS);
    }

    if ( new_indexentry(index,readname,len,m1->offset)==NULL) {
      PRINT_ERROR("Error in file %s: line %lu: malloc failed?",fd1->filename,fd1->cline-4);
      exit(SYS_INT_ERROR_EXIT_STATUS);
    }
    // TODO validate option
    if (fastq_validate_entry(fd1,m1)!=0) {
      exit(FASTQ_FORMAT_ERROR_EXIT_STATUS);
    }
    PRINT_READS_PROCESSED(fd1->cline/4,100000);
  }  
  //fastq_close(fd1->fd);
  return;
}

			  
char* fastq_get_readname(FASTQ_FILE* fd, FASTQ_ENTRY* e,char* rn,unsigned long *len_p,int is_header1) {
  unsigned long len=0;
  char *hdr;
  if ( is_header1) hdr=e->hdr1;
  else  hdr=e->hdr2;

  if ( is_header1==TRUE && hdr[0]!='@' ) {
    PRINT_ERROR("Error in file %s: line %lu: wrong header %s",fd->filename,fd->cline,hdr);
    exit(FASTQ_FORMAT_ERROR_EXIT_STATUS);
  }

  // rn=&rn[1];// ignore/discard @
  if ( strncpy(rn,&hdr[1],MAX_LABEL_LENGTH-1)==NULL ) {
    PRINT_ERROR("Error in strcpy");
    exit(SYS_INT_ERROR_EXIT_STATUS);
  }
  // executed only once
  if ( fd->readname_format == UNDEF ) {
      fd->is_casava_18=is_casava_1_8_readname(rn);
      if (fd->is_casava_18) {
        fprintf(stdout,"CASAVA=1.8\n");
        fd->readname_format=CASAVA18;
      } else {
	int is_int_name=is_int_readname(rn);
	if ( is_int_name ) {
	  fprintf(stdout,"Read name provided as an integer\n");
	  fd->readname_format=INTEGERNAME;
	} else {
	  int no_suffix=is_nosuffix_readname(rn);
	  if ( no_suffix ) {
	    fprintf(stdout,"Read name provided with no suffix\n");
	    fd->readname_format=NOP;
	  } else 
	    fd->readname_format=DEFAULT;
	}
      }
  }
  
  if ( fd->space==UNDEFSPACE ) {
    fd->space=is_color_space(e->seq,fd);
    if ( fd->space==COLORSPACE ) {
      fprintf(stdout,"Color space\n");
    }
  }
  
  
  switch(fd->readname_format) {
  case DEFAULT:
    // discard last character if PE && not casava 1.8
    len=strlen(rn);
    if (fd->is_pe) 
      len--;
    rn[len-1]='\0';
    break;
    
  case INTEGERNAME: // == NOP
    // keep the sequence unchanged
    len=strlen(rn);
    rn[len-1]='\0';
    break;
  case CASAVA18:
    len=0;
    while (rn[len]!=' ' && rn[len]!='\0') ++len;
    rn[len]='\0';
    if  ( rn[len-2] == '/' ) {
      // discard /[12]
      rn[len-2]='\0';
      len=len-2;
    }
    break;
  }
  *len_p=len;
  //fprintf(stderr,"read=%s=\n",s);
  return(rn);
}


//unsigned long long qual_vals[126]; // distribution of quality values
/*
sdbm
this algorithm was created for sdbm (a public-domain reimplementation of ndbm) database library.
it was found to do well in scrambling bits, causing better distribution of the keys and fewer splits.
it also happens to be a good general hashing function with good distribution.
the actual function is hash(i) = hash(i - 1) * 65599 + str[i]; what is included below is the faster version used in gawk.
[there is even a faster, duff-device version] the magic constant 65599 was picked out of thin air while experimenting with
different constants, and turns out to be a prime. this is one of the algorithms used in berkeley db (see sleepycat) and elsewhere.
*/
static ulong hashit(char *str) {

  ulong hash = 0;
  int c;
  
  while ((c = *str++))
    hash = c + (hash << 6) + (hash << 16) - hash;
  
  return(hash);
}



// return 1 if the headers are the same...0 otherwise
static inline int compare_headers(const char *hdr1,const char *hdr2) {

  unsigned int slen=0;
  //fprintf(stderr,">%s<\n>%s<\n",hdr1,hdr2);
  // no readname in header2
  if ( hdr2[0]=='\n' || hdr2[0]=='\r' || hdr2[0]=='\0' ) {
    return 1;
  }
  while ( hdr1[slen]!='\0' && hdr2[slen]!='\0' ) {
    if ( hdr1[slen]!=hdr2[slen] ) break;
    slen++;
  }
  // ignore white spaces
  unsigned int slen2=slen;
  while ( hdr1[slen]!='\0' ) {
    if ( hdr1[slen]!='\r'  &&  hdr1[slen]!='\n' ) return 0;
    ++slen;
  }
  while ( hdr2[slen2]!='\0' ) {
    if ( hdr2[slen2]!='\r'  &&  hdr2[slen2]!='\n' ) return 0;
    ++slen2;
  }
  return 1;
}

void fastq_index_delete(char *rname,hashtable index) {
  unsigned long key=hashit(rname);
  INDEX_ENTRY* e=fastq_index_lookup_header(index,rname);
  if (hash_delete(index,key,e)!=e) {
    PRINT_ERROR("Unable to delete entry from index");
    exit(SYS_INT_ERROR_EXIT_STATUS);
  }
  free_indexentry(e);  
}
INDEX_ENTRY* fastq_index_lookup_header(hashtable sn_index,char *hdr) {
  // lookup hdr in sn_index
  ulong key=hashit(hdr);
  //printf("looking for %s: key=%lu\n",hdr,key);
  INDEX_ENTRY* e=(INDEX_ENTRY*)get_object(sn_index,key);
  while (e!=NULL) {      // confirm that hdr are equal
    if ( !strcmp(hdr,e->hdr)) break;
    e=(INDEX_ENTRY*)get_next_object(sn_index,key);
  }
  return e;
}

//long collisions[HASHSIZE+1];
INDEX_ENTRY* new_indexentry(hashtable ht,char*hdr,int len,long start_pos) {
  
  // Memory chunck: |[index_entry]len bytes+1|
  char *mem_block=(char*)malloc(sizeof(INDEX_ENTRY)+len+1);
  if (mem_block==NULL) { return(NULL);}  
  INDEX_ENTRY *e=(INDEX_ENTRY*)&mem_block[0];

  e->hdr=(char*)&mem_block[sizeof(INDEX_ENTRY)];
  e->entry_start=start_pos;
  
  strncpy(e->hdr,hdr,len);
  e->hdr[len]='\0';
  // add to hash table
  ulong key=hashit(e->hdr);
  //collisions[key%HASHSIZE]++;
  if(insere(ht,key,e)<0) {
    PRINT_ERROR("Error while adding %s to index",hdr);
    return(NULL);
  }
  index_mem+=sizeof(INDEX_ENTRY)+len+1+sizeof(hashnode);
  return(e);
}

void free_indexentry(INDEX_ENTRY *e) {
  free(e);
  // remove entry from hash table
  return;
}

void fastq_destroy(FASTQ_FILE* fd) {
  if (!fd) return;
  if (fd->mode == FASTQ_GZIP)
    gzclose(fd->fd.gz);
  else
    fclose(fd->fd.fp);
  free(fd);
}

static inline void fastq_close(gzFile fd) {
  //if (fd==NULL) { return; }
  if (gzclose(fd)!=Z_OK) {
    PRINT_ERROR("unable to close file descriptor");
    exit(SYS_INT_ERROR_EXIT_STATUS);
  }
}

gzFile fastq_open(const char* filename,const char *mode) {
  gzFile fd1;

  if ( filename[0]=='-' && filename[1]=='\0' ) {
    if (mode[0]=='r') {
      //SET_BINARY_MODE(stdin);
      fd1 = gzdopen(fileno(stdin), "rb");
      if (fd1 == NULL) {
	PRINT_ERROR("Unable to gzdopen stdin");
	exit(PARAMS_ERROR_EXIT_STATUS);
      }
    } else {
      //SET_BINARY_MODE(stdout);
      fd1 = gzdopen(fileno(stdout), "wb");
      if (fd1 == NULL) {
	PRINT_ERROR("Unable to gzdopen stdout");
	exit(PARAMS_ERROR_EXIT_STATUS);
      }
    }
  } else {
    fd1=gzopen(filename,mode);
    if (fd1==NULL) {
      PRINT_ERROR("Unable to open %s",filename);
      exit(PARAMS_ERROR_EXIT_STATUS);
    }
  }
  //gzbuffer(fd1,sizeof(FASTQ_ENTRY)*2);
  // too large value slows down seek
  gzbuffer(fd1,128000);
  return(fd1);
}


//  http://support.illumina.com/help/SequencingAnalysisWorkflow/Content/Vault/Informatics/Sequencing_Analysis/CASAVA/swSEQ_mCA_FASTQFiles.htm
// check if the read name format was generated by casava 1.8
int is_casava_1_8_readname(const char *s) {
  regex_t regex;
  int reti;
  int is_casava_1_8=FALSE;
  //reti = regcomp(&regex,"[A-Z0-9:]* [1234]:[YN]:[0-9]*:.*",0);
  // relaxed format
  reti = regcomp(&regex,"[A-Z0-9:]* [1234]:[YN]:[0-9]*.*",0);  
  if ( reti ) { 
    PRINT_ERROR("Internal error: Could not compile regex"); 
    exit(SYS_INT_ERROR_EXIT_STATUS); 
  }
  /* Execute regular expression */
  //fprintf(stderr,"%s\n",hdr);
  reti = regexec(&regex, s, 0, NULL, 0);
  if ( !reti ) {    // match
    is_casava_1_8=TRUE;
  } 
  regfree(&regex);
  return is_casava_1_8;
}

// is the read name a plain integer
// or it doesn't contain a #/ as the second last letter
static int is_int_readname(const char *s) {
  regex_t regex;
  int reti;
  int is_int_name=FALSE;
  // @ was alread removed
  reti = regcomp(&regex,"^[0-9]+[\n\r]?$",REG_EXTENDED);  
  if ( reti ) { 
    PRINT_ERROR("Internal error: Could not compile regex"); 
    exit(SYS_INT_ERROR_EXIT_STATUS); 
  }
  /* Execute regular expression */
  //fprintf(stderr,">%s<\n",s);
  reti = regexec(&regex, s, 0, NULL, 0);
  if ( !reti ) {    // match
    is_int_name=TRUE;
  } 
  regfree(&regex);
  return is_int_name;
}

static int is_nosuffix_readname(const char *s) {
  regex_t regex;
  int reti;
  int is_nosuffix_name=TRUE;
  // @ was alread removed
  reti = regcomp(&regex,"[# \t/:][0-9abAB][\n\r]?$",REG_EXTENDED);  
  if ( reti ) { 
    PRINT_ERROR("Internal error: Could not compile regex"); 
    exit(SYS_INT_ERROR_EXIT_STATUS); 
  }
  /* Execute regular expression */
  //fprintf(stderr,">%s<\n",s);
  reti = regexec(&regex, s, 0, NULL, 0);
  if ( !reti ) {    // match
    is_nosuffix_name=FALSE;
  } 
  regfree(&regex);
  return is_nosuffix_name;
}


/* check if the read is in colour space or sequence space */
READ_SPACE is_color_space(char *seq,FASTQ_FILE* f) {
  regex_t regex;
  int reti;
  if (f->space!=UNDEFSPACE)
    return(f->space);
  //fprintf(stderr,">>%s<<\n",seq);
  reti = regcomp(&regex,"^[GT]?[0123n\\.NtT]+\n?$",REG_EXTENDED);  
  if ( reti ) { 
    PRINT_ERROR("Internal error: Could not compile regex"); 
    exit(SYS_INT_ERROR_EXIT_STATUS); 
  }
  /* Execute regular expression */
  //fprintf(stderr,">%s<\n",s);
  reti = regexec(&regex, seq, 0, NULL, 0);
  regfree(&regex);
  READ_SPACE rs;
  if ( !reti ) {    // match
    rs=COLORSPACE;
  } else {
    rs=SEQSPACE;
  }
  f->space=rs;
  return(rs);  
}

static inline char* PLAIN_READ(FILE *fp, char *s, long max) {
  if (fgets(s, max, fp) != NULL) return s;
  s[0] = '\0';
  return NULL;
}

static inline char* FASTQ_READ(FASTQ_FILE *f, char *s, long max) {
  if (f->mode == FASTQ_GZIP)
    return gzgets(f->fd.gz, s, max);
  else
    return PLAIN_READ(f->fd.fp, s, max);
}

static inline void FASTQ_WRITE(FASTQ_FILE *f, const char *s) {
  if (f->mode == FASTQ_GZIP) {
    if (gzputs(f->fd.gz, s) < 0) {
      PRINT_ERROR("gzputs failed");
      exit(SYS_INT_ERROR_EXIT_STATUS);
    }
  } else {
    if (fputs(s, f->fd.fp) == EOF) {
      PRINT_ERROR("fputs failed");
      exit(SYS_INT_ERROR_EXIT_STATUS);
    }
  }
}

/* inline long replace_dot_by_N(char* seq) { */
/*   long n=0;   */
/*   long replaced=0; */
/*   while (seq[n]!='\0') { */
/*     if (seq[n]=='.') { */
/*       ++replaced; */
/*       seq[n]='N'; */
/*     } */
/*     ++n; */
/*   } */
/*   return replaced; */
/* } */

/* inline long replace_dots(long long start,char* seq, char *hdr1, char *hdr2, char *qual,gzFile fd) { */
/*   // FIX . */
/*   long replaced=0; */
/*   //if ( fix_dot ) { */
/*     replaced+=replace_dot_by_N(seq); */
/*     gzprintf(fd,"%s%s%s%s",hdr1,seq,hdr2,qual); */
/*     //} */
/*   return(replaced); */
/* } */
