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
#include <stdlib.h>
#include <string.h>

#include "hash.h"

#define BUCKET(table,i) table->buckets[i]
#define LAST_ENTRY(table,i) table->buckets_last[i]
#ifndef HASHSIZE
#define HASHSIZE(t) t->size
#endif
static ulong mhash(hashtable,ulong);
static hashnode* hash_lookup(hashtable,ulong);


static hashnode* hash_lookup(hashtable table,ulong  key){
  
  table->last_node = BUCKET(table,mhash(table,key)); /* set a pointer to the first bucket */
  while ( table->last_node != NULL ) {
    if( table->last_node->value==key) return table->last_node;
    table->last_node = table->last_node->next;
  }
  return NULL;
}
__ptr_t get_next_object(hashtable table,ulong key)
{
  if(table->last_node==NULL)
    return NULL; 
  table->last_node = table->last_node->next;
  while ( table->last_node != NULL ) {
     if( table->last_node->value==key) return table->last_node->obj;
     table->last_node = table->last_node->next;
  }
  return NULL;
}


/* removes the element with key 'key' and returns the object stored on him */
__ptr_t hash_delete(hashtable table,ulong key,__ptr_t obj)
{
  hashnode *b,*prev=NULL;
  ulong c=mhash(table,key);
  b=BUCKET(table,c); /* set a pointer to the first bucket */
  while( b!=NULL) {
    if( b->value==key && b->obj==obj){
      if(prev==NULL) /* first element */
	BUCKET(table,c)=b->next;
      else 
	prev->next=b->next;
      if (b->next==NULL) // last
	LAST_ENTRY(table,c)=prev;
      free(b);
      table->n_entries--;
      return obj;
    }
    prev = b;
    b = b->next;
  };
  return NULL;
}

/* __ptr_t replace_object(hashtable table,ulong  key,__ptr_t newobj) */
/* { */
/*   __ptr_t old; */
/*   hashnode *b=hash_lookup(table,key);  */
/*   if(b==NULL)return NULL; */
/*   old=b->obj; */
/*   b->obj=newobj; */
/*   return old;   */
/* } */

/* looks a 'bucket' in the hashing table whith 'key' and return the
 pointer to the object stored in that bucket or NULL if no bucket is found */ 
__ptr_t get_object(hashtable table,ulong key){
  
   hashnode *b=hash_lookup(table,key); 
   if(b==NULL)
       return NULL;

   return b->obj; 
}

/* Allocates space to a new hash table */
hashtable new_hashtable(ulong hashsize) {
  hashtable new_hash;

  if( (new_hash = (hashtable)malloc(sizeof(struct hashtable_s)))==NULL) return NULL;
  
  if( (new_hash->buckets = (hashnode**)malloc(sizeof(hashnode*)*hashsize))==NULL) 
     return NULL;
  if( (new_hash->buckets_last = (hashnode**)malloc(sizeof(hashnode*)*hashsize))==NULL) 
     return NULL;
  memset(new_hash->buckets,0,sizeof(hashnode*)*hashsize);
  memset(new_hash->buckets_last,0,sizeof(hashnode*)*hashsize);
  new_hash->size=hashsize;
  new_hash->last_bucket=0;
  new_hash->last_node=NULL;
  new_hash->n_entries=0;
  //for(i=0;i<hashsize;++i) { BUCKET(new,i) = NULL; LAST_ENTRY(new,i) = NULL; }
  return new_hash;
}

void hashtable_stats(hashtable table) {
  ulong zbuckets=0;
  ulong collisions=0;
  ulong max_col=0;
  ulong i,ctr;
  for(i=0;i<HASHSIZE(table);++i) {
    hashnode *b=BUCKET(table,i);
    if ( b==NULL ) ++zbuckets;
    else {
      ctr=0;
      while (b!=NULL) {
	b=b->next;
	++ctr;
      }
      if ( ctr > 1 ) {
	collisions+=ctr;
	if ( ctr > max_col ) max_col=ctr;
	//fprintf(stderr,"%llu %llu\n",i,ctr);
      }
    }
  }
  fprintf(stderr,"size: %llu\n",HASHSIZE(table));
  fprintf(stderr,"max. col: %llu\n",max_col);
  fprintf(stderr,"zbuckets: %llu\n",zbuckets);
  fprintf(stderr,"%% zbuckets: %.2f\n",zbuckets*1.0/HASHSIZE(table));
  fprintf(stderr,"collisions: %llu\n",collisions);
  fprintf(stderr,"avg. collisions: %.2f\n",collisions*1.0/(HASHSIZE(table)-zbuckets));
}

/* A very simple hashing function */
static ulong mhash(hashtable table,ulong key)
{
  return (ulong)(key%HASHSIZE(table));
}

/* inserts a new element in the hash table*/
int insere(hashtable table,ulong key,__ptr_t obj)
{
   ulong ind;
   hashnode *new_node;
   if((new_node=(hashnode *)malloc(sizeof(hashnode)))==NULL) return -1;
   ind=mhash(table,key);
   // add to the end of the list
   new_node->value=key;
   new_node->obj=obj;
   if ( LAST_ENTRY(table,ind)==NULL ) {
     // first element
     BUCKET(table,ind)=new_node;
   } else {
     LAST_ENTRY(table,ind)->next=new_node;
   }
   new_node->next=NULL;
   LAST_ENTRY(table,ind)=new_node;
   // add the new entry to the front of the list
   //new->next = BUCKET(table,ind);
   //BUCKET(table,ind)=new;   
   //
   table->n_entries++;
   return 1;
}

void free_hashtable(hashtable table)
{
   ulong i;
   hashnode *n,*tmp;
   //fprintf(stderr,"free_hashtable\n");fflush(stderr);
   if (table==NULL) return;
   for(i=0;i<HASHSIZE(table);++i) {
      n=BUCKET(table,i);
      while(n!=NULL) {
         tmp=n;
         n=n->next;
         free(tmp);
      }      
   }
   free(table->buckets);
   free(table->buckets_last);
   free(table);
}
void reset_hashtable(hashtable table)
{
   ulong i;
   hashnode *n,*tmp;
   //fprintf(stderr,"free_hashtable\n");fflush(stderr);
   if (table==NULL) return;
   for(i=0;i<HASHSIZE(table);++i) {
      n=BUCKET(table,i);
      while(n!=NULL) {
         tmp=n;
         n=n->next;
         free(tmp);
      }      
   }
}
/*********************************************************************************/
/*
 * Returns all objects stored in a basket by making successive calls
 */
void init_hash_traversal(hashtable table) {
  table->last_bucket=0;
  table->last_node=NULL;
}
/*
 * Returns all objects stored in a basket by making successive calls
 */
__ptr_t next_hash_object(hashtable table)
{
  // first time....
  if( table->last_bucket>=HASHSIZE(table)) 
    return NULL;
    
  if( table->last_node==NULL ) {
    // find bucket
    // find next bucket
    while ( table->last_node == NULL && table->last_bucket+1<HASHSIZE(table)) {
      ++table->last_bucket;
      table->last_node = BUCKET(table,table->last_bucket);
    }
    if (table->last_node==NULL)
      return NULL;
    return table->last_node->obj;
  } 
  // Next in bucket
  table->last_node=table->last_node->next;
  if (table->last_node==NULL) return next_hash_object(table);
  return table->last_node->obj;
}

/*
 * Returns all hash nodes stored in a basket by making successive calls
 */
__ptr_t next_hashnode(hashtable table)
{
  // first time....
  if( table->last_bucket>=HASHSIZE(table)) 
    return NULL;
    
  if( table->last_node==NULL ) {
    // find bucket
    // find next bucket
    while ( table->last_node == NULL && table->last_bucket+1<HASHSIZE(table)) {
      ++table->last_bucket;
      table->last_node = BUCKET(table,table->last_bucket);
    }

    if (table->last_node==NULL)
      return NULL;
    return table->last_node;
  } 
  // Next in bucket
  table->last_node=table->last_node->next;
  if (table->last_node==NULL) return next_hashnode(table);
  return table->last_node;
}
/*
 * Returns all objects stored in a basket by making successive calls
 * the returned objects are deleted from the hash table
 */
/* __ptr_t next_delete_hash_object(hashtable table) */
/* { */
/*   // first time.... */
/*   if( table->last_bucket>=HASHSIZE(table))  */
/*     return NULL; */
    
/*   if( table->last_node==NULL ) { */
/*     // find bucket */
/*     // find next bucket */
/*     while ( table->last_node == NULL && table->last_bucket+1<HASHSIZE(table)) { */
/*       ++table->last_bucket; */
/*       table->last_node = BUCKET(table,table->last_bucket); */
/*     } */
/*     if (table->last_node==NULL) */
/*       return NULL; */
/*     // delete */
/*     void *obj=table->last_node->obj; */
/*     free(table->last_node); */
/*     table->last_node=NULL; */
/*     return obj; */
/*   }  */
/*   // Next in bucket */
/*   table->last_node=table->last_node->next; */
/*   if (table->last_node==NULL) return next_hash_object(table); */
/*   void *obj=table->last_node->obj; */
/*   free(table->last_node); */
/*   table->last_node=NULL; */
/*   return obj; */
/* } */


