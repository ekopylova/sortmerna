/*
 * SortMeRNA - next-generation reads filter for metatranscriptomic or total RNA
 * Copyright (C) 2014 Bonsai Bioinformatics Research Group
 *
 * This file is part of SortMeRNA.
 *
 * SortMeRNA is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * SortMeRNA is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with Foobar.  If not, see <http://www.gnu.org/licenses/>.
 *
 * file: indexdb.cpp
 * contact: jenya.kopylov@gmail.com, laurent.noe@lifl.fr, helene.touzet@lifl.fr
 *
 */

/** @file */ 

#include <time.h>
#include <string.h>
#include <dirent.h>
#include <errno.h>
#include <sys/types.h>
#include <unistd.h>
#include "../include/indexdb.hpp"
#include "../cmph/cmph.h"
#include <map>
#include <iomanip>
#include <sstream>
#include <vector>
#include <sys/stat.h> //for creating tmp dir
#include <deque>



//! burst trie nucleotide map
/*! the trie nodes consist of an array holding four 
    NodeElement structs, they are traversed by their
    index ranging from 0-4. 

    Use ascii decimal value of a letter as an offset in this array such that:
	A/a -> 0, C/c -> 1, G/g -> 2, T/t -> 3, U/u -> 3
	R/r -> 0, Y/y -> 1, S/s -> 2, W/w -> 1, K/k -> 2
	M/m -> 0, B/b -> 1, D/d -> 0, H/h -> 0, V/v -> 0
	N/n -> 0	*/
char map_nt[122] = {
	// 0,   1,   2,   3,   4,   5,   6,   7,   8,   9   
    0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
	// 10,  11,  12,  13,  14,  15,  16,  17,  18,  19  
    0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
	// 20,  21,  22,  23,  24,  25,  26,  27,  28,  29  
    0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
	// 30,  31,  32,  33,  34,  35,  36,  37,  38,  39 
    0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
	// 40,  41,  42,  43,  44,  45,  46,  47,  48,  49  
    0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
	// 50,  51,  52,  53,  54,  55,  56,  57,  58,  59  
    0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
	// 60,  61,  62,  63,  64,  65,  66,  67,  68,  69  
    0,   0,   0,   0,   0,   0,   1,   1,   1,   0,
	// 70,  71,  72,  73,  74,  75,  76,  77,  78,  79 
    0,   2,   0,   0,   0,   2,   0,   0,   0,   0,
	// 80,  81,  82,  83,  84,  85,  86,  87,  88,  89 
    0,   0,   0,   2,   3,   3,   0,   1,   2,   1,
	// 90,  91,  92,  93,  94,  95,  96,  97,  98,  99 
    0,   0,   0,   0,   0,   0,   0,   0,   1,   1,
	// 100, 101, 102, 103, 104, 105, 106, 107, 108, 109
    0,   0,   0,   2,   0,   0,   0,   2,   0,   0,
	// 110, 111, 112, 113, 114, 115, 116, 117, 118, 119 
    0,   0,   0,   0,   0,   2,   3,   3,   0,   1,
	// 120, 121 
    2,   1};

char map_nt_mask_N[122] = {
  /* 0,   1,   2,   3,   4,   5,   6,   7,   8,   9   */
    5,   5,   5,   5,   5,   5,   5,   5,   5,   5,
  /* 10,  11,  12,  13,  14,  15,  16,  17,  18,  19  */
    5,   5,   5,   5,   5,   5,   5,   5,   5,   5,
  /* 20,  21,  22,  23,  24,  25,  26,  27,  28,  29  */
    5,   5,   5,   5,   5,   5,   5,   5,   5,   5,
  /* 30,  31,  32,  33,  34,  35,  36,  37,  38,  39  */
    5,   5,   5,   5,   5,   5,   5,   5,   5,   5,
  /* 40,  41,  42,  43,  44,  45,  46,  47,  48,  49  */
    5,   5,   5,   5,   5,   5,   5,   5,   5,   5,
  /* 50,  51,  52,  53,  54,  55,  56,  57,  58,  59  */
    5,   5,   5,   5,   5,   5,   5,   5,   5,   5,
  /* 60,  61,  62,  63,  64,  65,  66,  67,  68,  69  */
    5,   5,   5,   5,   5,   0,   5,   1,   5,   5,
  /* 70,  71,  72,  73,  74,  75,  76,  77,  78,  79  */
    5,   2,   5,   5,   5,   5,   5,   5,   5,   5,
  /* 80,  81,  82,  83,  84,  85,  86,  87,  88,  89  */
    5,   5,   5,   5,   3,   3,   5,   5,   5,   5,
  /* 90,  91,  92,  93,  94,  95,  96,  97,  98,  99  */
    5,   5,   5,   5,   5,   5,   5,   0,   5,   1,
  /* 100, 101, 102, 103, 104, 105, 106, 107, 108, 109 */
    5,   5,   5,   2,   5,   5,   5,   5,   5,   5,
  /* 110, 111, 112, 113, 114, 115, 116, 117, 118, 119 */
    5,   5,   5,   5,   5,   5,   3,   3,   5,   5,
  /* 120, 121 */
    5,   5};


// length of the sliding window parameters
uint32_t lnwin_gv = 0;
uint32_t pread_gv = 0;
uint32_t partialwin_gv = 0;

uint32_t numseq_gv = 0;

// bit masking during sliding of window by 1 character
uint32_t mask32 = 0;
uint64_t mask64 = 0;


uint32_t total_num_trie_nodes = 0;
uint32_t size_of_all_buckets = 0;
uint32_t sizeoftrie = 0;

// STATISTICS
uint32_t total_num_buckets = 0;
uint32_t largest_bucket_size = 0;
uint32_t high_num_elem_in_bucket = 0;
uint32_t low_num_elem_in_bucket = 1000;
uint32_t longest_elem_in_bucket = 0;
uint32_t shortest_elem_in_bucket = 1000;
uint32_t all_lengths_elements_in_buckets = 0;
uint32_t all_elem_in_buckets = 0;
uint32_t avg_len[11] = {0};
uint32_t num_elem[100] = {0};

bool verbose = false;

// change version number here
char version_num[] = "2.0, 16/07/2014";

/*! @function mean()
    @brief compute the mean of an integer vector
    @return double, mean
*/
inline double mean(vector<uint32_t>& counts)
{
  uint32_t sum = 0;
  uint32_t size = counts.size();
  for (uint32_t i = 0; i < size; i++)
    sum+=counts[i];
  double mean = static_cast<double>(sum)/static_cast<double>(size);
  return mean;
}


/*! @function variance()
    @brief compute the variance using a vector of integers
    and their mean
    @return double, variance
*/
inline double variance(vector<uint32_t>& counts, double mean)
{
  uint32_t size = counts.size();
  if (size == 0)
  {
    fprintf(stderr, "  ERROR: size of k-mer count vector is 0.\n");
    exit(EXIT_FAILURE);
  }
  double sum = 0;
  double inverse = 1.0 / static_cast<double>(size);
  for(uint32_t i = 0; i < size; i++)
  {
    sum += pow(static_cast<double>(counts[i]) - mean, 2);
  }
  //double stddev = sqrt(sum * inverse);
  double variance = sum * inverse;
  return variance;
}

/*! @function location_param()
    @brief compute the location parameter (mu) using the normal
    mean and variance values
    @return double, mu
*/
inline double location_param(double mean, double variance)
{
  double mean_sqr = pow(mean,2);
  double mu = log(mean_sqr/sqrt(variance+mean_sqr));
  return mu;
}

/*! @function scale_param()
    @brief compute the scale parameter (mu) using the normal
    mean and variance values
    @return double, sigma
*/
inline double scale_param(double mean, double variance)
{
  double sigma_sqr = log(1+(variance/pow(mean,2)));
  return sigma_sqr;
}

/*! @function lognormal_mean()
    @brief compute the lognormal mean using mu and sigma^2
    @return double, l_mean
*/
inline double lognormal_mean(double mu, double sigma_sqr)
{
  double l_mean = exp(mu + sigma_sqr/2);
  return l_mean;
}

/*! @function lognormal_variance()
    @brief compute the lognormal variance using mu and sigma^2
    @return double, l_variance
*/
inline double lognormal_variance(double mu, double sigma_sqr)
{
  double l_variance = (exp(sigma_sqr)-1)*exp(2*mu + sigma_sqr);
  return l_variance;
}



/*
 *
 * @function insert_prefix() add an element to the mini-burst trie,
 * the part of the 19-mer not represented as trie nodes, will be
 * recorded in a bucket using 2 bits per nt
 * @param NodeElement* trie_node
 * @param char* str
 * @param uint32_t sig
 * @return void
 * @version 1.0 Dec 11, 2012
 *
 *******************************************************************/
inline void insert_prefix( NodeElement* trie_node,
                          unsigned char *prefix)
{
	uint32_t depth = 0;
	uint32_t sig = 0;
    
	// get the trie node from which to start traversal (A,C,G or T)
	NodeElement *node_elem = (NodeElement*)(trie_node + *prefix++);
	depth++;
    
  // find the terminal trie node
	while ( node_elem->flag == 1 )
	{
		trie_node = node_elem->whichnode.trie;
		node_elem = (NodeElement*)(trie_node + *prefix++);
		depth++;
	}
        
	// a bucket does not exist, create one
	if ( node_elem->flag == 0 )
	{
	 	// set the baseptr to the beginning of the new bucket
		node_elem->whichnode.bucket = (void*)malloc(ENTRYSIZE);
		if ( node_elem->whichnode.bucket == NULL )
		{
      fprintf(stderr,"  %sERROR%s: could not allocate memory for bucket "
                     "(insert_prefix() in indexdb.cpp)\n","\033[0;31m","\033[0m");
			exit(EXIT_FAILURE);
		}
		// initialize bucket memory to 0
		memset(node_elem->whichnode.bucket, 0, ENTRYSIZE);
		// set flag to signify the existence of a bucket
		node_elem->flag = 2;
		node_elem->size = 0;
	}
    
	// a bucket does exist, allocate memory for 1 more entry
	else
	{
		uint32_t node_elem_size = node_elem->size;
		void* src = node_elem->whichnode.bucket;
        
		node_elem->whichnode.bucket = (void*)malloc(node_elem_size+ENTRYSIZE);
		if ( node_elem->whichnode.bucket == NULL )
		{
			fprintf(stderr,"  %sERROR%s: could not allocate memory for bucket "
                     "resize (insert_prefix() in indexdb.cpp): %s\n","\033[0;31m",
                     "\033[0m",strerror(errno));
			exit(EXIT_FAILURE);
		}
        
		memset(node_elem->whichnode.bucket, 0, node_elem_size+ENTRYSIZE);
		memcpy((void*)node_elem->whichnode.bucket, (void*)src, node_elem_size);
        
		free(src);
	}
    
	// add tail to bucket
	uint32_t* entry = (uint32_t*)((unsigned char*)node_elem->whichnode.bucket + node_elem->size);

	// length of entry to add to bucket
	int s = partialwin_gv+1-depth;
    
	// add the tail to the bucket
	uint32_t encode = 0;
	for ( int i = 0; i < s; i++ )
	{
		encode |= ((uint32_t)*prefix++)<<(2*i);
	}
    
	*entry++ = encode;
    
	// add the signature of the prefix following the tail
	*entry = sig;
    
	// record the new size of bucket
	(node_elem->size)+=ENTRYSIZE;
	
#define BURST
#ifdef BURST
	// (-3 = -2*k-1) smallest bucket must have at least 3-character strings
	if ( depth < (pread_gv-partialwin_gv-3) )
	{
		// burst if next string will exceed bucket limit
		if ( node_elem->size > THRESHOLD )
		{
			// create a new trie node
			NodeElement* child_node = (NodeElement*)malloc(4*sizeof(NodeElement));
			if ( child_node == NULL )
			{
        fprintf(stderr,"  %sERROR%s: could not allocate memory for child_node "
                       "(insert_prefix())\n","\033[0;31m","\033[0m");
				exit(EXIT_FAILURE);
			}
			memset(child_node, 0, 4*sizeof(NodeElement));
            
			unsigned char* start_bucket_to_burst = (unsigned char*)node_elem->whichnode.bucket;
			unsigned char* end_bucket_to_burst = (unsigned char*)node_elem->whichnode.bucket + node_elem->size;
            
			// read every entry in the bucket to burst
			while ( start_bucket_to_burst != end_bucket_to_burst )
			{
        NodeElement* node_elem_child = (NodeElement*)(child_node + ((*start_bucket_to_burst)&3));
        
        // create a new child bucket
        if ( node_elem_child->flag == 0 )
        {
          node_elem_child->whichnode.bucket = (void*)malloc(ENTRYSIZE);
          if (node_elem_child->whichnode.bucket == NULL)
          {
              fprintf(stderr,"  %sERROR%s: could not allocate memory for child bucket "
                             "(insert_prefix())\n","\033[0;31m","\033[0m");
              exit(EXIT_FAILURE);
          }
          memset(node_elem_child->whichnode.bucket, 0, ENTRYSIZE);
          node_elem_child->flag = 2;
          node_elem_child->size = 0;
        }
        // flag == 2, resize the existing bucket to add an extra element
        else
        {
          uint32_t child_bucket_size = node_elem_child->size;
          void* src = node_elem_child->whichnode.bucket;
          node_elem_child->whichnode.bucket = (void*)malloc(child_bucket_size+ENTRYSIZE);
          if ( node_elem_child->whichnode.bucket == NULL )
          {
              fprintf(stderr,"  %sERROR%s: could not allocate memory for child bucket resize "
                             "(insert_prefix() in indexdb.cpp)\n","\033[0;31m","\033[0m");
              exit(EXIT_FAILURE);
          }
          memset(node_elem_child->whichnode.bucket, 0, child_bucket_size+ENTRYSIZE);
          memcpy(node_elem_child->whichnode.bucket, src, child_bucket_size);
          free(src);
        }
        
        // shift all bits in an entry to make up for the removed nucleotide
        *((uint32_t*)start_bucket_to_burst)>>=2;
        
        // add an element to the bucket
        memcpy( (unsigned char*)node_elem_child->whichnode.bucket + node_elem_child->size, start_bucket_to_burst, ENTRYSIZE);
        
        (node_elem_child->size)+=ENTRYSIZE;
        
        // each entry is 2 uint32_ts
        start_bucket_to_burst+=ENTRYSIZE;
			}
			
			// reset the flag and size of the parent node
			free(node_elem->whichnode.bucket);
			node_elem->whichnode.trie = child_node;
      // now points to the child node
			node_elem->flag = 1;
      // size of bucket = 0
			node_elem->size = 0;

		}//~burst bucket
	}//~depth < pread_gv-2k-1
#endif
    
	return ;
}//~insert_prefix()



/*
 *
 * @function add a k-mer occurrence to the table of occurrences
 * @param kmer_origin* positions_tbl
 * @param uint32_t seq: the sequence number
 * @param uint32_t pos: the position on sequence
 * @return void
 * @version 1.0 Jan 16, 2013
 *
 *******************************************************************/
inline void add_kmer_to_table( kmer_origin* positions_tbl,
                               uint32_t seq,
                               uint32_t pos,
                               bool store)
{
	int32_t size = positions_tbl->size;

  // do not store positions for this abundant k-mer
  if (!store)
  {
    positions_tbl->size = -1;
    return;
  }
    
	// create an entry for the new unique k-mer occurrence
	if ( size == 0 )
	{
		positions_tbl->arr = (seq_pos*)malloc(sizeof(seq_pos));
		memset(positions_tbl->arr,0,sizeof(seq_pos));
	}
	// add a position for an existing k-mer occurrence
	else
	{        
		seq_pos* src = positions_tbl->arr;
		positions_tbl->arr = (seq_pos*)malloc((size+1)*sizeof(seq_pos));
		memset(positions_tbl->arr,0,(size+1)*sizeof(seq_pos));
		memcpy(positions_tbl->arr, src, size*sizeof(seq_pos));
        
		free(src);
	}
    
	positions_tbl->arr[size].seq = seq;
	positions_tbl->arr[size].pos = pos;
	positions_tbl->size++;
    
}//~add_kmer_to_table


/*
 *
 * @function free the mini-burst trie
 * @param NodeElement* trie_node
 * @return void
 * @version 1.0 Dec 11, 2012
 *
 *******************************************************************/
void freebursttrie( NodeElement* trie_node )
{
	/* traverse through the node elements in a trie node */
	for ( int i = 0; i < 4; i++ )
	{
		unsigned char value = trie_node->flag;
        
		/* the node element holds a pointer to another trie node */
		if ( value == 1 )
		{
			freebursttrie ( (NodeElement*)(trie_node->whichnode.trie) );
			free(trie_node++->whichnode.trie);
		}
        
		/* the node element points to a bucket */
		else if ( value == 2 )
		{
			free( trie_node++->whichnode.bucket );
		}
        
		/* the node element is empty, go to next node element */
		else if ( value == 0 ) { 	trie_node++; }
        
	}//~trie nodes
    
	return;
    
}//~freebursttrie()



/*
 *
 * @function search_burst_trie: test whether a 19-mer is in the trie
 * @param NodeElement* trie_node: pointer to mini-burst trie
 * @param unsigned char* kmer_short_key: pointer to second half of
 * 19-mer window
 * @param bool &new_position: true if 18-mer prefix of 19-mer exists in
 * the burst trie, false otherwise
 * @return bool: true if 19-mer is in the burst trie, false otherwise
 * @version 1.0 Mar 4, 2013
 *
 *******************************************************************/
bool search_burst_trie( NodeElement* trie_node, unsigned char* kmer_short_key, bool &new_position )
{
	uint32_t depth = 0;
    
	// find a terminal trie node
	NodeElement *node_elem = (NodeElement*)(trie_node + *kmer_short_key++);
	depth++;
    
	while ( node_elem->flag == 1 )
	{
		trie_node = node_elem->whichnode.trie;
		node_elem = (NodeElement*)(trie_node + *kmer_short_key++);
		depth++;
	}
    
	if ( node_elem->flag == 0 ) return false;
    
	// encode the remaining part of kmer_short_key using 4 nt per byte
	int s = partialwin_gv+1-depth;
    
	uint32_t encode = 0;
	for ( int i = 0; i < s; i++ )
	{
		encode |= ((uint32_t)*kmer_short_key++)<<(2*i);
	}
    
	// size of entry in a bucket
	unsigned char* start_bucket = (unsigned char*)node_elem->whichnode.bucket;
	unsigned char* end_bucket = start_bucket + node_elem->size;
    
	// to mask the last nucleotide of 19-mer
	uint32_t msk = (1<<(2*(s-1)))-1;
    
	// compare the 1 int representation of kmer_id_short_F with all elements in the bucket
	while ( start_bucket != end_bucket )
	{
		// 18-mer found
		if ( (encode&msk) == (*((uint32_t*)start_bucket)&msk ) )
		{
			new_position = false;
			// 19-mer found
			if ( encode == *((uint32_t*)start_bucket) ) return true;
		}
		start_bucket+=ENTRYSIZE;
	}
    
	// end of bucket reached, 19-mer not found
	return false;
	
}//~search_burst_trie()






/*
 *
 * @function add_id_to_burst_trie: initially all id's in the burst trie are set
 * to 0, here we set them to their proper MPHF values, forward 19-mer and reverse 19-mer must have the same id
 * @param NodeElement* trie_node: pointer to mini-burst trie
 * @param unsigned char* kmer_id_short_F_ptr: pointer to second half of
 * 19-mer window
 * @param uint32_t id: the true id of 19-mer
 * @return void
 * @version 1.0 Jan 11, 2013
 *
 *******************************************************************/
void add_id_to_burst_trie( NodeElement* trie_node, unsigned char* kmer_id_short_F_ptr, uint32_t id )
{
	uint32_t depth = 0;
    
	// find a terminal trie node
	NodeElement *node_elem = (NodeElement*)(trie_node + *kmer_id_short_F_ptr++);
	depth++;
    
	while ( node_elem->flag == 1 )
	{
		trie_node = node_elem->whichnode.trie;
		node_elem = (NodeElement*)(trie_node + *kmer_id_short_F_ptr++);
		depth++;
	}
    
	// encode the remaining part of kmer_id_short_F using 4 nt per byte
	int s = partialwin_gv+1-depth;
    
	uint32_t encode = 0;
	for ( int i = 0; i < s; i++ )
	{
		encode |= ((uint32_t)*kmer_id_short_F_ptr++)<<(2*i);
	}
    
	unsigned char* start_bucket = (unsigned char*)node_elem->whichnode.bucket;
	unsigned char* end_bucket = start_bucket + node_elem->size;
    
	// compare the 1 byte representation of kmer_id_short_F with all elements in the bucket
	while ( start_bucket != end_bucket )
	{
		// set the id
    if ( encode == *((uint32_t*)start_bucket ) )
		{
			*((uint32_t*)(start_bucket+sizeof(uint32_t))) = id;
		}
		start_bucket+=ENTRYSIZE;
	}
	
	return ;
}//~add_id_to_burst_trie()


void search_for_id( NodeElement* trie_node, unsigned char* kmer_id_short_F_ptr, uint32_t &id )
{
	uint32_t depth = 0;
    
	// find a terminal trie node
	NodeElement *node_elem = (NodeElement*)(trie_node + *kmer_id_short_F_ptr++);
	depth++;
    
	while ( node_elem->flag == 1 )
	{
		trie_node = node_elem->whichnode.trie;
		node_elem = (NodeElement*)(trie_node + *kmer_id_short_F_ptr++);
		depth++;
	}
    
	// encode the remaining part of kmer_id_short_F using 4 nt per byte
	int s = partialwin_gv+1-depth;
    
	uint32_t encode = 0;
	for ( int i = 0; i < s; i++ )
	{
		encode |= ((uint32_t)*kmer_id_short_F_ptr++)<<(2*i);
	}
    
	unsigned char* start_bucket = (unsigned char*)node_elem->whichnode.bucket;
	unsigned char* end_bucket = start_bucket + node_elem->size;
    
	// compare the 1 byte representation of kmer_id_short_F with all elements in the bucket
	while ( start_bucket != end_bucket )
	{
		// set the id
    if ( encode == *((uint32_t*)start_bucket ) )
		{
      id = *((uint32_t*)(start_bucket+sizeof(uint32_t)));
		}
		start_bucket+=ENTRYSIZE;
	}
	
	return ;
}//~search_for_id()




/*
 *
 * @function traversetrie: collect statistics on the mini-burst trie,
 * it's size, number of trie nodes vs. buckets
 * @param NodeElement* trie_node
 * @return void
 * @version 1.0 Jan 14, 2013
 *
 *******************************************************************/
void traversetrie( NodeElement* trie_node, uint32_t depth )
{
	total_num_trie_nodes++;
    
	// traverse through the node elements in a trie node
	for ( int i = 0; i < 4; i++ )
	{
		unsigned char value = trie_node->flag;
        
		// the node element holds a pointer to another trie node
		if ( value == 1 )
		{
			traversetrie( trie_node->whichnode.trie, ++depth );
			--depth;
		}
        
		// the node element points to a bucket
		else if ( value == 2 )
		{
			// pad to alignment length (16-byte line)
			// int padding = 16-((trie_node->size)%16);
			// size_of_all_buckets+=(trie_node->size + padding);
			size_of_all_buckets+=trie_node->size;
            
      if ( largest_bucket_size < trie_node->size ) largest_bucket_size = trie_node->size;
            
			// for STATISTICS
			total_num_buckets++;
			uint32_t s = partialwin_gv-depth;
            
			if ( s > longest_elem_in_bucket ) longest_elem_in_bucket = s;
			else if ( s < shortest_elem_in_bucket ) shortest_elem_in_bucket = s;
            
			all_lengths_elements_in_buckets+=s;
            
			avg_len[s]++;
            
			uint32_t numelem = (trie_node->size)/ENTRYSIZE;
            
			if ( numelem > high_num_elem_in_bucket ) high_num_elem_in_bucket = numelem;
			else if ( numelem < low_num_elem_in_bucket ) low_num_elem_in_bucket = numelem;
            
			all_elem_in_buckets+=numelem;
            
			num_elem[numelem]++;
		}
        
		// the node element is empty, go to next node element
		else if ( value == 0 )
    {
      ;
    }
        
		trie_node++;
        
	}//~trie nodes
    
	return;
    
}//~taversetrie()



#ifdef see_binary_output
/*
 *
 * @function traversetrie: collect statistics on the mini-burst trie,
 * it's size, number of trie nodes vs. buckets
 * @param NodeElement* trie_node
 * @return void
 * @version 1.0 Jan 14, 2013
 *
 *******************************************************************/
void traversetrie_debug( NodeElement* trie_node, uint32_t depth, uint32_t &total_entries, string &kmer_keep )
{
  char get_char[4] = {'A','C','G','T'};
    
	// traverse through the node elements in a trie node
	for ( int i = 0; i < 4; i++ )
	{
		unsigned char value = trie_node->flag;
        
		// the node element holds a pointer to another trie node
		if ( value == 1 )
		{
      kmer_keep.push_back((char)get_char[i]); //TESTING
			traversetrie_debug( trie_node->whichnode.trie, ++depth, total_entries, kmer_keep );
      kmer_keep.pop_back();
			--depth;
		}
        
		// the node element points to a bucket
		else if ( value == 2 )
		{
      kmer_keep.push_back((char)get_char[i]); //TESTING
      
      unsigned char* start_bucket = (unsigned char*)trie_node->whichnode.bucket;
      if ( start_bucket == NULL )
      {
        fprintf(stderr, "  ERROR: pointer start_bucket == NULL (paralleltraversal.cpp)\n");
        exit(EXIT_FAILURE);
      }
      unsigned char* end_bucket = start_bucket + trie_node->size;
      if ( end_bucket == NULL )
      {
        fprintf(stderr, "  ERROR: pointer end_bucket == NULL (paralleltraversal.cpp)\n");
        exit(EXIT_FAILURE);
      }
      
      // traverse the bucket
      while ( start_bucket != end_bucket )
      {
        uint32_t entry_str = *((uint32_t*)start_bucket);
        uint32_t s = partialwin_gv-depth;
        total_entries++;
        
        // for each nt in the string
        for ( uint32_t j = 0; j < s; j++ )
        {
            kmer_keep.push_back((char)get_char[entry_str&3]); //TESTING
            entry_str>>=2;
        }
        
        cout << kmer_keep << endl; //TESTING
        
        for ( uint32_t j = 0; j < s; j++ ) kmer_keep.pop_back(); //TESTING
        
        // next entry
        start_bucket+=ENTRYSIZE;
      }//~for each entry
      
      kmer_keep.pop_back(); //TESTING
            
		}
        
		/// the node element is empty, go to next node element
		else if ( value == 0 ) ;
        
		trie_node++;
        
	}//~trie nodes
    
	return;
    
}//~taversetrie_debug()
#endif



/*
 *
 * @function load index: write to binary file the 9-mer look-up
 * tables and the mini-burst tries
 * @param string root: the file name of the index
 * @param kmer* lookup_table: pointer to the 9-mer lookup table
 * @return void
 * @version 1.0 Jan 16, 2013
 *
 *******************************************************************/
void load_index( kmer* lookup_table, char* outfile )
{
	// output the mini-burst tries
	ofstream btrie ( outfile, ofstream::binary );
    
	uint32_t sizeoftries[2] = {0};
    
	// loop through all 9-mers
	for ( uint32_t i = 0; i < (uint32_t)(1<<lnwin_gv); i++ )
	{
		NodeElement* trienode = NULL;
        
#ifdef see_binary_output
    cout << "9-mer = " << i; //TESTING
#endif
        
		// 1. output size for the two mini-burst tries for each 9-mer
		for ( int j = 0; j < 2; j++ )
		{
			total_num_trie_nodes = 0;
			size_of_all_buckets = 0;
			if ( j == 0 ) trienode = lookup_table[i].trie_F;
			else trienode = lookup_table[i].trie_R;

      if ( trienode != NULL ) traversetrie (trienode, 0 );
            
			sizeoftrie = total_num_trie_nodes*sizeof(NodeElement)*4 + size_of_all_buckets*sizeof(char);
            
			sizeoftries[j] = sizeoftrie;
            
			btrie.write(reinterpret_cast<const char*>(&sizeoftrie), sizeof(uint32_t));
            
#ifdef see_binary_output
      if ( j == 0 ) cout << "\tsizeoftrie f = " << sizeoftrie; //TESTING
      else cout << "\tsizeoftrie r = " << sizeoftrie; //TESTING
#endif
		}
        
#ifdef see_binary_output
    cout << "\tlookup_tbl[i].count = " << lookup_table[i].count << endl; //TESTING
#endif
        
		// 2. output both mini-burst tries into binary file
		for ( int j = 0; j < 2; j++ )
		{
			// the mini-burst trie exists, load into memory
			if ( sizeoftries[j] != 0 )
			{
				if ( j == 0 )
        {
          trienode = lookup_table[i].trie_F;
#ifdef see_binary_output
          cout << "forward burst-trie \n"; //TESTING
#endif
        }
				else
        {
          trienode = lookup_table[i].trie_R;
#ifdef see_binary_output
          cout << "reverse burst-trie \n"; //TESTING
#endif
        }
                
				// queue of node elements for breadth-first traversal
				deque<NodeElement*> nodes;
                
				// load first set of NodeElements into the queue & write to file
				for ( int i = 0; i < 4; i++ )
				{
					nodes.push_back( trienode );
					btrie.write(reinterpret_cast<const char*>(&(trienode->flag)), sizeof(char));
#ifdef see_binary_output
          cout << " " << (int)trienode->flag; //TESTING
#endif
					trienode++;
				}
                
				int depth = 0;
				int poplimit = 4;
				int numpops = 0;
				int topop = 0;
                
				while ( !nodes.empty() )
				{
					// increment depth of burst trie
					if ( numpops == poplimit )
					{
						depth++;
						poplimit = topop;
						numpops = 0;
						topop = 0;
					}
                    
					trienode = nodes.front();
                    
					switch ( trienode->flag )
					{
            // empty node
						case 0:
						{
							;
						}
            break;
            // trie node, add child trie node to queue
						case 1:
						{
							NodeElement *child = trienode->whichnode.trie;
							for ( int i = 0; i < 4; i++ )
							{
								nodes.push_back( child );
								btrie.write(reinterpret_cast<const char*>(&(child->flag)), sizeof(char));
#ifdef see_binary_output
                cout << " " << (int)child->flag; //TESTING
#endif
								child++;
							}
							topop+=4;
						}
            break;
            // bucket node, add bucket to output file
						case 2:
						{
							char* bucket = (char*)(trienode->whichnode.bucket);
                            
							// bucket information
							uint32_t sizeofbucket = trienode->size;
                            
#ifdef see_binary_output
              cout << "\tsizeofbucket = " << sizeofbucket; //TESTING
#endif                 
							btrie.write(reinterpret_cast<const char*>(&sizeofbucket), sizeof(uint32_t));
                            
							// bucket content
							char* start = (char*)bucket;
                            
							btrie.write(reinterpret_cast<const char*>(start), sizeofbucket);
						}
            break;
            // ?
						default:
						{
							fprintf(stderr, "  %sERROR%s: flag is set to %d (load_index)\n","\033[0;31m","\033[0m",trienode->flag);
							exit(EXIT_FAILURE);
						}
            break;
					}
                    
					nodes.pop_front();
					numpops++;
                    
#ifdef see_binary_output
          if ( numpops%4 == 0) cout << "\n"; //TESTING
#endif             
				}//~while the queue is not empty
			}//~if mini-burst trie exists
		}//~for each mini-burst trie in the 9-mer
	}//~for each 9-mer
    
	btrie.close();
    
}//~load_index()



/*
 *
 * FUNCTION 	: welcome()
 * PARAMETERS	: none
 * PURPOSE	: program name, copyright information and contact
 *
 **************************************************************************************************************/
void welcome()
{
	printf("\n  Program:     SortMeRNA version %s\n",version_num );
	printf("  Copyright:   2012-2014 Bonsai Bioinformatics Research Group:\n");
	printf("               LIFL, University Lille 1, CNRS UMR 8022, INRIA Nord-Europe\n" );
  printf("               OTU-picking extensions developed in the Knight Lab,\n");
  printf("               BioFrontiers Institute, University of Colorado at Boulder\n");
  printf("  Disclaimer:  SortMeRNA comes with ABSOLUTELY NO WARRANTY; without even the\n");
  printf("               implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.\n");
  printf("               See the GNU Lesser General Public License for more details.\n");
  printf("  Contact:     Evguenia Kopylova, jenya.kopylov@gmail.com \n");
	printf("               Laurent Noé, laurent.noe@lifl.fr\n");
	printf("               Hélène Touzet, helene.touzet@lifl.fr\n\n");
}



/*
 *
 * @function printlist: output man page for SortMeRNA
 * @return void
 * @version 1.0 Jan 14, 2013
 *
 *******************************************************************/
void printlist()
{
	printf("\n  usage:   ./indexdb_rna --ref db.fasta,db.idx [OPTIONS]:\n\n");
  printf("  --------------------------------------------------------------------------------------------------------\n");
  printf("  | parameter        value           description                                                 default |\n");
  printf("  --------------------------------------------------------------------------------------------------------\n");
  printf("     %s--ref%s           %sSTRING,STRING%s   FASTA reference file, index file                            %smandatory%s\n","\033[1m","\033[0m","\033[4m","\033[0m","\033[0;32m","\033[0m");
  printf("                                      (ex. --ref /path/to/file1.fasta,/path/to/index1)\n");
  printf("                                       If passing multiple reference sequence files, separate\n");
  printf("                                       them by ':',\n");
  printf("                                      (ex. --ref /path/to/file1.fasta,/path/to/index1:/path/to/file2.fasta,path/to/index2)\n");
  printf("   [OPTIONS]:\n");
  printf("     %s--fast%s          %sFLAG%s            suggested option for aligning ~99%% related species          %soff%s\n","\033[1m","\033[0m","\033[4m","\033[0m","\033[4m","\033[0m");
  printf("     %s--sensitive%s     %sFLAG%s            suggested option for aligning ~75-98%% related species       %son%s\n","\033[1m","\033[0m","\033[4m","\033[0m","\033[4m","\033[0m");
  printf("     %s--tmpdir%s        %sSTRING%s          directory where to write temporary files\n","\033[1m","\033[0m","\033[4m","\033[0m");
  printf("     %s-m%s              %sINT%s             the amount of memory (in Mbytes) for building the index     %s3072%s \n","\033[1m","\033[0m","\033[4m","\033[0m","\033[4m","\033[0m");
  printf("     %s-L%s              %sINT%s             seed length                                                 %s18%s\n","\033[1m","\033[0m","\033[4m","\033[0m","\033[4m","\033[0m");
  #ifdef interval
  printf("     %s--interval%s      %sINT%s             index every INT L-mer in the reference database             %s1%s\n","\033[1m","\033[0m","\033[4m","\033[0m","\033[4m","\033[0m");
  #endif
  printf("     %s-S%s              %sINT%s             the scale factor for threshold = mean + S*stddev            %s50%s\n","\033[1m","\033[0m","\033[4m","\033[0m","\033[4m","\033[0m");
  printf("                                     this threshold sets the maximum number of occurrences\n");
  printf("                                     of a seed that can be stored.\n");
  printf("     %s--mask_N%s        %sFLAG%s            do not include k-mers having an ambiguous nucleotide        %soff%s\n","\033[1m","\033[0m","\033[4m","\033[0m","\033[4m","\033[0m");
  printf("                                             in the index\n");
	printf("     %s-v%s              %sFLAG%s            verbose\n","\033[1m","\033[0m","\033[4m","\033[0m");
	printf("     %s-h%s              %sFLAG%s            help	\n\n","\033[1m","\033[0m","\033[4m","\033[0m");
	exit(EXIT_FAILURE);
}//~printlist()




/*
 *
 * @function parse the fasta file, build the burst tries
 * @param ./indexdb --db [file name]
 * @return void
 * @version 1.0 Dec 10, 2012
 *
 **************************************************************************************************************/
int main (int argc, char** argv)
{
  int narg = 1;
  // time
  double	s = 0.0;
  double  f = 0.0;
  // memory of index
  double mem = 0;
  bool mem_is_set = false;
  bool fast_set = false;
	bool sensitive_set = false;
  bool lnwin_set = false;
  bool interval_set = false;
  // vector of (FASTA file, index name) pairs for constructing index
  vector< pair<string,string> > myfiles;
  // pointer to temporary directory
  char* ptr_tmpdir = NULL;
  uint32_t interval = 0;
  bool mask_N = false;
  bool mask_N_set = false;
  char *map_ptr = NULL;
  // threshold = mean + s_offset*stddev
  bool s_offset_set = false;
  uint32_t s_offset = 0;
    
	timeval t;
    
	if ( argc == 1 )
	{
		verbose = true;
		welcome();
		fprintf(stderr,"  For help or more information on usage, type `./indexdb_rna %s-h%s'\n\n","\033[1m","\033[0m");
		exit(EXIT_FAILURE);
	}
    
	while ( narg < argc )
	{
		switch( argv[narg][1] )
		{
			case '-':
			{
				char* myoption = argv[narg];
				// skip the '--'
				myoption+=2;
                
				// path to reference sequence file
				if ( strcmp ( myoption, "ref") == 0 )
				{
          // no files are given
					if ( argv[narg+1] == NULL )
					{
						fprintf(stderr,"\n  %sERROR%s: --ref must be followed by at least one entry (ex. --ref /path/to/file1.fasta,/path/to/index1)\n\n","\033[0;31m","\033[0m");
						exit(EXIT_FAILURE);
					}
          else
          {
            char *ptr = argv[narg+1];
            while ( *ptr != '\0' )
            {
              // get the FASTA file path + name
              char fastafile[2000];
              char *ptr_fastafile = fastafile;

              // the reference database FASTA file
              while ( *ptr != ',' && *ptr != '\0' )
              {
                *ptr_fastafile++ = *ptr++;
              }
              *ptr_fastafile = '\0';
              ptr++; //skip the ',' delimiter
                
              // check FASTA file exists & is not empty
              if ( FILE *file = fopen(fastafile, "r") )
             	{ 
             		// get file size
             		fseek(file, 0, SEEK_END);
             		size_t filesize = ftell(file);
             		// file size is 0, create empty index file and exit the program
             		if ( !filesize ) 
             		{
             			// get the index filepath
             			char indexfile[2000];
             			char *ptr_indexfile = indexfile;
             			while ( *ptr != ':' && *ptr != '\0' )
             			{
             				*ptr_indexfile++ = *ptr++;
             			}
             			*ptr_indexfile = '\0';

             			char bases[4][50] = {".bursttrie_0.dat", ".pos_0.dat", ".kmer_0.dat", ".stats"};

             			// output empty index files
             			for ( int file = 0; file < 4; file++ )
             			{
             				char str[2000];
             				strcpy (str, indexfile);
             				strcat (str, bases[file]);
             				FILE *t = fopen(str, "w");
             				fclose(t);
             			}
             			// exit
             			fprintf(stdout, "  The input file is empty, an index was not built.\n");
             			exit(EXIT_SUCCESS);
             		}
             		// file size > 0, reset file pointer to start of file
             		fseek(file, 0, SEEK_SET);
             		fclose(file);
             	}
              else
              {
                fprintf(stderr, "\n  %sERROR%s: the file %s could not be opened: %s.\n\n","\033[0;31m","\033[0m",fastafile,strerror(errno));
                exit(EXIT_FAILURE);
              }
                             
              // get the index path + name
              char indexfile[2000];
              char *ptr_indexfile = indexfile;
              // the reference database index name
              while ( *ptr != ':' && *ptr != '\0')
              {
                *ptr_indexfile++ = *ptr++;
              }
              *ptr_indexfile = '\0';
              if ( *ptr != '\0' ) ptr++; //skip the ':' delimiter
              
              // check the directory where to write the index exists
              char dir[500];
              char *ptr_end = strrchr( indexfile, '/');
              if ( ptr_end != NULL )
              {
                memcpy( dir, indexfile, (ptr_end-indexfile) );
                dir[(int)(ptr_end-indexfile)] = '\0';
              }
              else
              {
                strcpy( dir, "./" );
              }
              
              if ( DIR *dir_p = opendir(dir) ) closedir(dir_p);
              else
              {
                if ( ptr_end != NULL )
                  fprintf(stderr,"\n  %sERROR%s: the directory %s for writing "
                                 "index '%s' could not be opened. The full directory "
                                 "path must be provided (ex. no '~'). \n\n",
                                 "\033[0;31m","\033[0m",dir,ptr_end+1);
                else
                  fprintf(stderr,"\n  %sERROR%s: the directory %s for writing index "
                                 "'%s' could not be opened. The full directory path must "
                                 "be provided (ex. no '~'). \n\n","\033[0;31m",
                                  "\033[0m",dir,indexfile);
                
                exit(EXIT_FAILURE);
              }
                
              // check index file names are distinct
              for ( int i = 0; i < (int)myfiles.size(); i++ )
              {
                if ( (myfiles[i].first).compare(fastafile) == 0 )
                {
                  fprintf(stderr, "\n  %sWARNING%s: the FASTA file %s has "
                                  "been entered twice in the list. It will "
                                  "be indexed twice. \n\n","\033[0;33m",
                                  "\033[0m",fastafile);
                }
                else if ( (myfiles[i].second).compare(indexfile) == 0 )
                {
                  fprintf(stderr, "\n  %sERROR%s: the index name %s has "
                                  "been entered twice in the list. Index names "
                                  "must be unique.\n\n","\033[0;31m",
                                  "\033[0m",indexfile);
                  exit(EXIT_FAILURE);
                }
              }
              myfiles.push_back(pair<string,string>(fastafile,indexfile));    
            }
            narg+=2;
          }
				}
        // the tmpdir
        else if ( strcmp ( myoption, "tmpdir" ) == 0 )
        {
          if ( argv[narg+1] == NULL )
          {
            fprintf(stderr,"\n  %sERROR%s: a directory path must "
                           "follow the option --tmpdir (ex. /path/to/dir ) \n",
                           "\033[0;31m","\033[0m");
            exit(EXIT_FAILURE);
          }
          else
          {
              ptr_tmpdir = argv[narg+1];
              narg+=2;
          }
        }
        // deprecated option (not used)
        else if ( strcmp ( myoption, "max_pos" ) == 0 )
        {
          ;
        }
				else if ( strcmp ( myoption, "fast") == 0 )
				{
					if ( lnwin_gv > 0 )
					{
						fprintf(stderr,"\n  %sERROR%s: option -L INT cannot be used with "
                           "--fast or --sensitive (these options set default values "
                           "for L).\n\n","\033[0;31m","\033[0m");
						exit(EXIT_FAILURE);
					}
					else if ( !fast_set )
					{
						lnwin_gv = 24;
					}
					else
					{
            fprintf(stderr,"\n  %sERROR%s: --fast has already been set "
                           "once.\n\n","\033[0;31m","\033[0m");
            exit(EXIT_FAILURE);
					}
					fast_set = true;
					narg++;
				}
				else if ( strcmp ( myoption, "sensitive") == 0 )
				{
					if ( lnwin_gv > 0 )
					{
						fprintf(stderr,"\n  %sERROR%s: option -L INT cannot be used "
                           "with --sensitive or --fast (these options set "
                           "default values for L).\n\n","\033[0;31m","\033[0m");
						exit(EXIT_FAILURE);
					}
					else if ( !sensitive_set )
						lnwin_gv = 18;
					else
					{
            fprintf(stderr,"\n  %sERROR%s: --sensitive has already been "
                           "set once.\n\n","\033[0;31m","\033[0m");
            exit(EXIT_FAILURE);
					}
					sensitive_set = true;
					narg++;
				}
        else if ( strcmp ( myoption, "mask_N") == 0 )
        {
          if ( mask_N_set )
          {
            fprintf(stderr,"\n  %sERROR%s: --mask_N has already been "
                           "set once.\n\n","\033[0;31m","\033[0m");
            exit(EXIT_FAILURE);
          }
          mask_N_set = true;
          mask_N = true;
          narg++;
        }
        // Interval for constructing index on every INT words
        else if ( strcmp ( myoption, "interval" ) == 0 )
        {
          if (argv[narg+1] == NULL)
          {
            fprintf(stderr,"\n  %sERROR%s: --interval requires a positive "
                           "integer as input (ex. --interval 2).\n\n",
                           "\033[0;31m","\033[0m");
            exit(EXIT_FAILURE);
          }
          // set interval
          if ( !interval_set )
          {
            if ( argv[narg+1][0] == '-' )
            {
              fprintf(stderr,"\n  %sERROR%s: --interval requires a "
                             "positive integer as input (ex. --interval 2).\n\n",
                             "\033[0;31m","\033[0m");
              exit(EXIT_FAILURE);
            }
            else if (isdigit(argv[narg+1][0]))
            {
              interval = atoi(argv[narg+1]);
              narg+=2;
              interval_set = true;
            }
            else
            {
              fprintf(stderr,"\n  %sERROR%s: --interval requires a positive "
                             "integer as input (ex. --interval 2).\n\n",
                             "\033[0;31m","\033[0m");
              exit(EXIT_FAILURE);
            }
          }
          else
          {
            fprintf(stderr,"\n  %sERROR%s: --interval has been set "
                           "twice, please verify your choice\n\n",
                           "\033[0;31m","\033[0m");
            printlist();
          }
        }
				else
				{
					fprintf(stderr,"\n  %sERROR%s: unknown option --%s.\n\n",
                         "\033[0;31m","\033[0m",myoption);
					printlist();
					exit(EXIT_FAILURE);
				}
			}
      break;
      case 'S':
      {
        if ( !s_offset_set )
        {
          char *pEnd = NULL;
          s_offset = strtol( argv[narg+1], &pEnd, 10 );
          if ( !s_offset )
          {
            fprintf(stderr, "\n  %sERROR%s: -S [INT] must be a positive integer "
                            "value.\n\n","\033[0;31m","\033[0m");
            exit(EXIT_FAILURE);
          }
          narg+=2;
          s_offset_set = true;
        }
        else
        {
          fprintf(stderr,"\n  %sERROR%s: option -S can only be set once,\n\n",
                         "\033[0;31m","\033[0m");
          exit(EXIT_FAILURE);
        } 
      }
      break;      
			case 'L':
			{
				if ( fast_set || sensitive_set )
				{
					fprintf(stderr,"\n  %sERROR%s: -L INT cannot not be set when --fast "
                         "or --sensitive options are used (these options set default "
                         "values for L).\n\n","\033[0;31m","\033[0m");
					exit(EXIT_FAILURE);
				}
                
				if ( lnwin_set )
				{
					fprintf(stderr,"\n  %sERROR%s: option -L can only be set once.\n\n",
                         "\033[0;31m","\033[0m");
					exit(EXIT_FAILURE);
				}
                
				int lnwin_t = atoi(argv[narg+1]);
        lnwin_set = true;
                
				if ( lnwin_t <= 0 )
				{
					fprintf(stderr,"\n  %sERROR%s: -L must be a positive integer "
                         "(10, 12, 14, .. , 20).\n\n","\033[0;31m","\033[0m");
					exit(EXIT_FAILURE);
				}
				else if ( lnwin_t%2 == 1 )
				{
					fprintf(stderr,"\n  %sERROR%s: -L must be an even integer (10, 12, "
                         "14, .. , 20).\n\n","\033[0;31m","\033[0m");
					exit(EXIT_FAILURE);
				}
				else if ( (lnwin_t < 8) || (lnwin_t > 26) )
				{
					fprintf(stderr,"\n  %sERROR%s: -L must be between 8 and 26, inclusive.\n\n",
                         "\033[0;31m","\033[0m");
					exit(EXIT_FAILURE);
				}
				else
				{
          lnwin_gv = lnwin_t;
					pread_gv = lnwin_gv+1;
					partialwin_gv = lnwin_gv/2;
					narg+=2;
				}
			}
      break;          
			case 'm':
			{
				// set memory for index (in Mbytes)
				if ( !mem_is_set )
				{
          // RAM limit for mmap'ing reads in megabytes
          char *pEnd = NULL;
          mem = strtod( argv[narg+1], &pEnd );
          if ( !mem )
          {
            fprintf(stderr, "\n  %sERROR%s: -m [INT] must be a positive integer "
                            "value (in Mbyte).\n\n","\033[0;31m","\033[0m");
            exit(EXIT_FAILURE);
          }
					narg+=2;
					mem_is_set = true;
				}
				else
				{
					fprintf(stderr,"\n  %sERROR%s: option -m can only be set once,\n\n",
                         "\033[0;31m","\033[0m");
					exit(EXIT_FAILURE);
				}
			}
      break;      
			case 'v':
			{
				// verbose
				if ( !verbose )
				{
					verbose = true;
					narg++;
				}
			}
      break;
			case 'h':
			{
				// help
        welcome();
				printlist();
			}
      break;
			default :
			{
				fprintf(stderr,"\n  %sERROR%s: '%c' is not one of the options\n\n",
                       "\033[0;31m","\033[0m",argv[narg][1]);
				exit(EXIT_FAILURE);
			}
		}//~switch
	}//~while ( narg < argc )

	// check that the database file has been provided
  if ( myfiles.empty() )
  {
		fprintf(stderr,"\n  %sERROR%s: a FASTA reference database & index name "
                   "(--ref /path/to/file1.fasta,/path/to/index1) is mandatory "
                   "input.\n\n","\033[0;31m","\033[0m");
		exit(EXIT_FAILURE);
	}
    
	// set the default value for seed length
	if ( !lnwin_set ) lnwin_gv = 18;
  // set the default interval length
  if ( !interval_set ) interval = 1;
  // set the nucleotide to integer encode matrix
  if ( mask_N ) map_ptr = map_nt_mask_N;
  else map_ptr = map_nt;
  // default threshold for storing k-mer positions
  if ( !s_offset_set ) s_offset = 50;
    
	pread_gv = lnwin_gv+1;
	partialwin_gv = lnwin_gv/2;
    
	// default memory for building index (3072 Mbytes)
	if ( !mem_is_set ) mem = 3072;
    
	mask32 = (1<<lnwin_gv)-1;
	mask64 = (2ULL<<((pread_gv*2)-1))-1;
    
  // temporary file to store keys for building the CMPH
  int32_t pid = getpid();
  char pidStr[4000];
  sprintf(pidStr,"%d",pid);
  char keys_str[2000] = "";
  char* tmpdir_env = NULL;
    
  // tmpdir provided
  if ( ptr_tmpdir != NULL )
  {
    strcat(keys_str,ptr_tmpdir);
    char* ptr_tmpdir_t = ptr_tmpdir;
    while (*ptr_tmpdir_t++ != '\0');
    if (*(ptr_tmpdir_t-2) != '/') strcat(keys_str,"/");
    
    char try_str[4000] = "";
    strcat(try_str,keys_str);
    strcat(try_str,"test.txt");
    
    FILE *tmp = fopen(try_str, "w+");
    if ( tmp == NULL )
    {
      fprintf(stderr,"\n  %sERROR%s: cannot access directory %s: "
                     "%s\n\n","\033[0;31m","\033[0m",ptr_tmpdir,
                     strerror(errno));
      exit(EXIT_FAILURE);
    }
  }
  // tmpdir not provided, try $TMPDIR, /tmp and local directories
  else
  {
    bool try_further = true;
    
    // try TMPDIR
    tmpdir_env = getenv("TMPDIR");
    if ( (tmpdir_env != NULL) && (strcmp(tmpdir_env,"")!=0) )
    {
      strcat(keys_str,tmpdir_env);
      char* ptr_tmpdir_t = tmpdir_env;
      while (*ptr_tmpdir_t++ != '\0');
      if (*(ptr_tmpdir_t-2) != '/') strcat(keys_str,"/");
      
      char try_str[4000] = "";
      strcat(try_str,keys_str);
      strcat(try_str,"test.txt");
      
      FILE *tmp = fopen(try_str, "w+");
      if ( tmp == NULL )
      {
        fprintf(stderr,"\n  %sWARNING%s: no write permissions in "
                       "directory %s: %s\n","\033[0;33m","\033[0m",
                       tmpdir_env,strerror(errno));
        fprintf(stderr,"  will try /tmp/.\n\n");
        try_further = true;
      }
      else try_further = false;
    }
    // try "/tmp" directory
    if ( try_further )
    {
      FILE *tmp = fopen("/tmp/test.txt", "w+");
      if ( tmp == NULL )
      {
        fprintf(stderr,"\n  %sWARNING%s: no write permissions in "
                       "directory /tmp/: %s\n","\033[0;33m","\033[0m",
                       strerror(errno));
        fprintf(stderr,"  will try local directory.\n\n");
        try_further = true;
      }
      else
      {
        strcat(keys_str,"/tmp/");
        try_further = false;
      } 
    }
    // try the local directory
    if ( try_further )
    {
      FILE *tmp = fopen("./test.txt", "w+");
      if ( tmp == NULL )
      {
        fprintf(stderr,"\n  %sERROR%s: no write permissions in current "
                       "directory: %s\n","\033[0;31m","\033[0m",
                       strerror(errno));
        fprintf(stderr,"  Please set --tmpdir to a writable directory, "
                       "or change the write permissions in $TMPDIR, /tmp/ "
                       "or current directory.\n\n");
        exit(EXIT_FAILURE);
      }
      else
      {
        strcat(keys_str,"./");
        try_further = false;
      }
    }
  }
    
  strcat(keys_str, "sortmerna_keys_");
  strcat(keys_str,pidStr);
  strcat(keys_str,".txt");
  
  // the list of arguments is correct, welcome the user!
  if ( verbose ) welcome();
  
  eprintf("\n  Parameters summary: \n");
  eprintf("    K-mer size: %d\n",lnwin_gv+1);
  eprintf("    K-mer interval: %d\n",interval);
  if (mask_N)
    eprintf("    Do not include k-mers with ambiguous nucleotides in index\n");
  else
    eprintf("    Include k-mers with ambiguous nucleotides in index\n");
  
  eprintf("\n  Total number of databases to index: %d\n", (int)myfiles.size());
    
  // build index for each pair in --ref list
  for ( int newindex = 0; newindex < (int)myfiles.size(); newindex++ )
  {
    vector< pair<string,uint32_t> > sam_sq_header;
    // vector of structs storing information on which sequences from 
    // the original FASTA file were added to each index part
    vector<index_parts_stats> index_parts_stats_vec;
    
    // FASTA reference sequences input file
    FILE *fp = fopen((char*)(myfiles[newindex].first).c_str(),"r");
    if ( fp == NULL )
    {
      fprintf(stderr,"  %sERROR%s: could not open file %s\n",
                     "\033[0;31m","\033[0m",(char*)(myfiles[newindex].first).c_str());
      exit(EXIT_FAILURE);
    }
    
    eprintf("\n  Begin indexing file %s%s%s under index name %s%s%s: \n",
            "\033[0;34m",(char*)(myfiles[newindex].first).c_str(),
            "\033[0m","\033[0;34m",(char*)(myfiles[newindex].second).c_str(),
            "\033[0m");
    
    // get full file size
    fseek(fp,0L,SEEK_END);
    size_t filesize = ftell(fp);
    fseek(fp,0L,SEEK_SET);   
    
    // STEP 1 ************************************************************
    // For file part_0,
    // (a) compute the nucleotide background frequencies;
    // (b) length of total reference sequences (for Gumbel parameters);
    // (c) number of total sequences (for SW alignment)
    // (d) k-mer statistics
    
    // total number of reference sequences
    size_t strs = 0;
    // length of longest single sequence
    uint32_t maxlen = 0;
    // length of single sequence
    int32_t len = 0;
    // the percentage of each A/C/G/T in the database
    double background_freq[4] = {0.0,0.0,0.0,0.0};
    // total length of reference sequences
    int64_t full_len = 0;
    int nt = 0;
    
    eprintf("  Collecting sequence distribution statistics ..");
      
    TIME(s);
    do
    {
      nt = fgetc(fp);
        
      // name of sequence for SAM format @SQ
      char read_header[2000];
      char *pt_h = read_header;
      
      // start of read header
      if ( nt == '>' ) strs+=2;
      else
      {
        fprintf(stderr,"\n%sERROR%s: each read header of the database fasta "
                       "file must begin with '>';","\033[0;31m","\033[0m");
        fprintf(stderr,"\n  check sequence # %ld\n\n", strs);
        exit(EXIT_FAILURE);
      }
      
      // scan to end of header name
      bool stop = false;
      while ( nt != '\n' )
      {
        nt = fgetc(fp);
        if (nt != '\n' && nt != ' ' && nt != '\t' && !stop)
          *pt_h++ = nt;
        else
          stop = true;
      }
        
      *pt_h = '\0';          
      len = 0;
        
      // scan through the sequence, count its length
      nt = fgetc(fp);
      while ( nt != '>' && nt != EOF )
      {
        // skip line feed, carriage return or empty space in the sequence
        if ( nt != '\n' && nt != ' ' )
        {
          len++;
          if ( nt == 'A' || nt == 'a' ||
               nt == 'C' || nt == 'c' ||
               nt == 'G' || nt == 'g' ||
               nt == 'T' || nt == 't' ||
               nt == 'U' || nt == 'u') background_freq[(int)*(map_ptr+nt)]++;
        }
        nt = fgetc(fp);
      }
          
      // add sequence name and length to sam_header_
      string s(read_header);
      sam_sq_header.push_back(pair<string,uint32_t>(s,len));
      
      if ( nt != EOF ) ungetc(nt,fp);
      
      full_len+=len;
      
      // if ( len > maxlen ) then ( maxlen = rrnalen ) else ( do nothing )
      len > maxlen ? maxlen = len : maxlen;
        
    } while ( nt != EOF ); // read until end of file
    TIME(f);
    
    // set file pointer back to the beginning of file
    rewind(fp);
    
    eprintf("  done  [%f sec]\n", (f-s));
        
        
    /* STEP 1 END **********************************************/
    

    /* STEP 2 **************************************************/
    /* For every part of total index,
       (a) build the burst trie and set all 19-mer ids to 0
       (b) count the number of unique 19-mers in the database
       (c) output the unique 19-mers into a file for MPHF */
    
    // number of the index part
    uint16_t part = 0;
    // starting position given by ftell() where to
    // begin reading the reference sequences
    unsigned long int start_part = 0;
    // number of bytes of reference sequences to read
    unsigned long int seq_part_size = 0;
    // total size of index so far in bytes
    double index_size = 0;

    // for each unique 18-mer, its number of occurrences in the reference database
    map<unsigned long long int,uint32_t> kmer_count;
    
    // for each index part of the reference sequences
    do
    {
      // number of sequences in part size
      uint32_t numseq_part = 0;
      
      // set the file pointer to the beginning of the current part
      start_part = ftell(fp);
      
      // output file storing all s-mer (19-mer) words in the reference
      // sequences, required for CMPH to build minimal perfect
      // hash functions
      FILE *keys = fopen(keys_str, "w+");
      if ( keys == NULL )
      {
        fprintf(stderr,"  %sERROR%s: could not open %s file for writing\n",
                       "\033[0;31m","\033[0m",keys_str);
        exit(EXIT_FAILURE);
      }
            
      // count of unique 18-mers in database
      uint32_t number_elements = 0;
      
      // table storing occurrence of each 9-mer and pointers to
      // the forward and reverse burst tries
      kmer *lookup_table = (kmer*)malloc((1<<lnwin_gv)*sizeof(kmer));
      if ( lookup_table == NULL )
      {
        fprintf(stderr,"  %sERROR%s: could not allocate memory for "
                       "9-mer look-up table (indexdb.cpp)\n",
                       "\033[0;31m","\033[0m");
        exit(EXIT_FAILURE);
      }
        
      memset(lookup_table, 0, (1<<lnwin_gv)*sizeof(kmer));
      
      // bool vector to keep track which L/2-mers have been counted for
      // by the forward sliding L/2-mer
      vector<bool> incremented_by_forward((1<<lnwin_gv));
      
      // total size of index so far in bytes
      index_size = 0;
      
      eprintf("\n  start index part # %d: \n",part);
      
      eprintf("    (1/3) building burst tries ..");
                      
      TIME(s);
      // for the number of sequences for which the index is less than maximum (set by -m)
      // indexing the 19-mers will be done in the following manner:
      //
      //  000000000011111111112222222222333333333344444
      //  012345678901234567890123456789012345678901234
      //  ACTACTATCTAGTGTGCTAGCTAGTCATCGCTAGCTAGCTAGTCG
      //
      //  --------- ---------
      // we store all unique 18-mer positions (not 19-mer) because if an 18-mer on a read matches exactly to the prefix
      // or suffix of a 19-mer in the mini-burst trie, we need to recover all of the 18-mer occurrences in the database
      do
      {
        // start of current sequence in file
        long int start_seq = ftell(fp);
        nt = fgetc(fp);
        
        // scan to end of header name
        char header[1000] = " "; //tmp
        char* h_pt = header;
        while ( nt != '\n' )
        { 
          nt = fgetc(fp);
          *h_pt++ = nt; //tmp
        }
        *h_pt = '\n'; //tmp

        //cout << header; //tmp
        
        unsigned char* myseq = new unsigned char[maxlen];
        unsigned char* myseqr = new unsigned char[maxlen];
        uint32_t _j = 0;
        len = 0;
        
        nt = fgetc(fp);
        // encode each sequence using integer alphabet {0,1,2,3}
        while ( nt != '>' && nt != EOF )
        {
          // skip line feed, carriage return or empty space in the sequence
          if ( nt != '\n' && nt != ' ' )
          {
            len++;
            // exact character
            myseq[_j++] = *(map_ptr+nt);
          }
          nt = fgetc(fp);
        }

        // end of current sequence in file
        if ( nt != EOF ) ungetc(nt,fp);
        
        long int end_seq = ftell(fp);
        
        // check the addition of this sequence will not overflow the
        // maximum memory (estimated memory 10 bytes per L-mer)
        double estimated_seq_mem = (len-pread_gv+1)*9.5e-6;
          
        // the sequence alone is too large, it will not fit into maximum
        // memory, skip it
        if ( estimated_seq_mem > mem )
        {
          fseek(fp,start_seq,SEEK_SET);
          fprintf(stderr,"\n  %sWARNING%s: the index for sequence `","\033[0;33m","\033[0m");
          int c = 0;
          do
          {
            c = fgetc(fp);
            fprintf(stderr,"%c",(char)c);
          } while ( c != '\n' );
          
          fprintf(stderr,"` will not fit into %e Mbytes memory, it will be skipped.",mem);
          fprintf(stderr,"  If memory can be increased, please try `-m %e` Mbytes.",estimated_seq_mem);
          fseek(fp,end_seq,SEEK_SET);
          continue;
        }
        // the additional sequence will overflow the maximum index memory,
        // write existing index to disk and start a new index
        else if ( index_size+estimated_seq_mem > mem )
        {
          // set the character to something other than EOF
          if ( nt == EOF ) nt = 'A';
          
          // reset the file pointer to the beginning of current sequence
          // (which will be added to the next index part)
          fseek(fp,start_seq,SEEK_SET);
          
          break;
        }
        // add the additional sequence to the index
        else
        {
          index_size+=estimated_seq_mem;
          
          // record the number of bytes of raw reference sequences added to this part
          seq_part_size = ftell(fp) - start_part;
          // record the number of sequences in this part
          numseq_part++;
        }
                         
        // create a reverse sequence using the forward
        unsigned char* ptr = &myseq[len-1];
        
        for ( _j = 0; _j < len; _j++ ) myseqr[_j] = *ptr--;
        // 9-mer prefix of 19-mer
        uint32_t kmer_key_short_f = 0;
        // 9-mer suffix of 19-mer
        uint32_t kmer_key_short_r = 0;
        // pointer to next letter to add to 9-mer prefix
        unsigned char* kmer_key_short_f_p = &myseq[0];
        // pointer to next letter to add to 9-mer suffix
        unsigned char* kmer_key_short_r_p = &myseq[partialwin_gv+1];
        // pointer to 10-mer of reverse 19-mer to insert
        // into the mini-burst trie
        unsigned char* kmer_key_short_r_rp = &myseqr[len-partialwin_gv-1];
        // 19-mer
        unsigned long long int kmer_key = 0;
        // pointer to 19-mer
        unsigned char* kmer_key_ptr = &myseq[0];

        bool skip_kmer = false; //new
        uint32_t skip_kmer_pos = 0; //new
          
        // initialize the prefix and suffix 9-mers
        for ( uint32_t j = 0; j < partialwin_gv; j++ )
        {
            (kmer_key_short_f <<= 2) |= (int)*kmer_key_short_f_p++;
            (kmer_key_short_r <<= 2) |= (int)*kmer_key_short_r_p++;
        }
            
        // initialize the 19-mer
        for ( uint32_t j = 0; j < pread_gv; j++ )
        {
          // skip k-mers with ambiguous nucleotides
          if ( mask_N )
          {
            if ( (int)*kmer_key_ptr == 5 )
            {
              skip_kmer = true;
              skip_kmer_pos = j;
            }
          }
          (kmer_key <<= 2) |= (int)*kmer_key_ptr++;
        } 
        
        uint32_t numwin = (len-pread_gv+interval)/interval; //TESTING
        uint32_t index_pos = 0;

        // for all 19-mers on the sequence
        for ( uint32_t j = 0; j < numwin; j++ ) //TESTING
        {
          // no ambiguous nucleotides in the k-mer (will only be true if mask_N is set)
          if ( !skip_kmer )
          {
            lookup_table[kmer_key_short_f].count++;
            incremented_by_forward[kmer_key_short_f] = true;
            // increment 9-mer count only if it wasn't already
            // incremented by kmer_key_short_f before
            if ( !incremented_by_forward[kmer_key_short_r] ) lookup_table[kmer_key_short_r].count++;
            
            // ****** add the forward 19-mer
            
            // new position for 18-mer in positions_tbl
            bool new_position = true;
            
            // forward 19-mer does not exist in the burst trie (duplicates not allowed)
            if ( lookup_table[kmer_key_short_f].trie_F == NULL ||
                ( (lookup_table[kmer_key_short_f].trie_F != NULL) && !search_burst_trie( lookup_table[kmer_key_short_f].trie_F, kmer_key_short_f_p, new_position ) ) )
            {
              // create a trie node if it doesn't exist
              if ( lookup_table[kmer_key_short_f].trie_F == NULL )
              {
                lookup_table[kmer_key_short_f].trie_F = (NodeElement*)malloc(4*sizeof(NodeElement));
                if ( lookup_table[kmer_key_short_f].trie_F == NULL )
                {
                  fprintf(stderr,"  %sERROR%s: could not allocate memory for trie_node in indexdb.cpp\n","\033[0;31m","\033[0m");
                  exit(EXIT_FAILURE);
                }
                memset(lookup_table[kmer_key_short_f].trie_F, 0, 4*sizeof(NodeElement));
              }
              insert_prefix( lookup_table[kmer_key_short_f].trie_F, kmer_key_short_f_p );
            }

            // 18-mer doesn't exist in the burst trie, add it to keys file
            if ( new_position )
            {
              // increment number of unique 18-mers
              number_elements++;
              fprintf(keys,"%llu\n",(kmer_key>>2));
              // add new 18-mer to kmer_count
              kmer_count.insert(pair<unsigned long long int, uint32_t>(kmer_key>>2, 1));

              //if ((kmer_key>>2) == 66091079) cout << header << "\tpos = " << j << endl; //tmp
            }
            else
            {
              // incremember 18-mer count
              kmer_count[(kmer_key>>2)]++;

              //if ((kmer_key>>2) == 66091079) cout << header << "\tpos = " << j << endl; //tmp
            }
    
            // ****** add the reverse 19-mer
            new_position = true;
            
            // reverse 19-mer does not exist in the burst trie
            if ( lookup_table[kmer_key_short_r].trie_R == NULL ||
                ( (lookup_table[kmer_key_short_r].trie_R != NULL) && !search_burst_trie( lookup_table[kmer_key_short_r].trie_R, kmer_key_short_r_rp, new_position ) ) )
            {
              // create a trie node if it doesn't exist
              if ( lookup_table[kmer_key_short_r].trie_R == NULL )
              {
                lookup_table[kmer_key_short_r].trie_R = (NodeElement*)malloc(4*sizeof(NodeElement));
                if ( lookup_table[kmer_key_short_r].trie_R == NULL )
                {
                  fprintf(stderr,"  %sERROR%s: could not allocate memory for trie_node in indexdb.cpp\n","\033[0;31m","\033[0m");
                  exit(EXIT_FAILURE);
                }
                memset(lookup_table[kmer_key_short_r].trie_R, 0, 4*sizeof(NodeElement));
              }
              
              insert_prefix( lookup_table[kmer_key_short_r].trie_R, kmer_key_short_r_rp );
            }
          } 
          // shift 19-mer window and both 9-mers
          if ( j != numwin-1 )
          {
            for ( int shift = 0; shift < interval; shift++ )
            {
              (( kmer_key_short_f <<= 2 ) &= mask32 ) |= (int)*kmer_key_short_f_p++;
              (( kmer_key_short_r <<= 2 ) &= mask32 ) |= (int)*kmer_key_short_r_p++;

              // skip k-mers with ambiguous nucleotides
              if ( mask_N )
              {
                if ((int)*kmer_key_ptr == 5 )
                {
                  skip_kmer = true;
                  skip_kmer_pos = j+pread_gv;
                }
              }

              (( kmer_key <<= 2 ) &= mask64 ) |= (int)*kmer_key_ptr++;
              kmer_key_short_r_rp--;
              index_pos++;
            }

            // skip k-mers with ambiguous nucleotides
            if ( mask_N )
            {
              if ( skip_kmer && (j+1 > skip_kmer_pos) ) skip_kmer = false;
            }
          }
                  
        }//~for all 19-mers on the sequence
            
        delete [] myseq;
        delete [] myseqr;
            
      } while ( nt != EOF ); // all file

      // tmp
      map<unsigned long long int, uint32_t>::iterator it;
      vector<uint32_t> counts;
      //map<uint32_t,uint32_t> stats; //tmp

      //cout << "k-mer\t# occurrences\n";
      for (it = kmer_count.begin(); it != kmer_count.end(); it++)
      {
                        //0    1    2    3
        //char decode[] = {'A', 'C', 'G', 'T'};
        //unsigned long long mask_t = 18446744073709551612ULL;
        //char forward[18] = " ";
        //char *forward_p = forward;
        //unsigned long long int kkmer = it->first;

        //for (int i = 0; i < lnwin_gv; i++ )
        //{
        //  *forward_p++ = (char)decode[(kkmer & 3)];
        //  kkmer>>=2;
        //}

        //forward_p--;

        //cout << ">" << it->second << "\n";
        //for (int i = 0; i < lnwin_gv; i++ )
        //{
        //  cout << (char)*forward_p--;
        //}

        //cout << "\n";

        //cout << it->first << "\t" << it->second << endl; //tmp
        counts.push_back(it->second);
        //if (stats.find(it->second) == stats.end()) //tmp
        //  stats.insert(pair<uint32_t,uint32_t>(it->second,1)); //tmp
        //else
        //  stats[it->second]++; //tmp
      }

      double mean_c = mean(counts);
      double variance_c = variance(counts, mean_c);
      double stddev_c = sqrt(variance_c);
      uint32_t count_threshold = mean_c+(s_offset*stddev_c);

      cout << "\nmean = " << mean_c << endl;
      cout << "stddev = " << sqrt(variance_c) << endl;
      cout << "threshold = " << count_threshold << endl;

      //double mu = location_param(mean_c, variance_c);
      //cout << "mu = " << mu << endl;
      //double sigma_sqr = scale_param(mean_c, variance_c);
      //cout << "sigma_sqr = " << sigma_sqr << endl;
      //double mean_log = lognormal_mean(mu, sigma_sqr);
      //double variance_log = lognormal_variance(mu, sigma_sqr);
      //double stddev_log = sqrt(variance_log);
      //cout << "lognormal mean = " << mean_log << endl;
      //cout << "lognormal stddev = " << stddev_log << endl;
      //cout << "lognormal threshold = " << mean_log+(3*stddev_log) << endl;

      //map<uint32_t,uint32_t>::iterator it2; //tmp
      //cout << "# occurrences\t# k-mers\n"; //tmp
      //for (it2 = stats.begin(); it2 != stats.end(); it2++ ) //tmp
      //  cout << it2->first << "\t" << it2->second << endl; //tmp
    
      TIME(f);
      
      // no index can be created, all reference sequences are too large to
      // fit alone into maximum memory
      if ( index_size == 0 )
      {
        eprintf("\n  %sERROR%s: no index was created, all of your sequences are "
                "too large to be indexed with the current memory limit of %e Mbytes.\n",
                "\033[0;31m","\033[0m",mem);
        break;
      }
      // continue to build hash and positions tables
      else index_size = 0;
      
      rewind(keys);
      
      eprintf(" done  [%f sec]\n", (f-s));
      
      // 4. build MPHF on the unique 18-mers
      eprintf("    (2/3) building CMPH hash ..");
      TIME(s);
      cmph_t *hash = NULL;
      
      FILE * keys_fd = keys;
      if (keys_fd == NULL)
      {
        fprintf(stderr, "File \"%s\" not found\n",keys_str);
        exit(EXIT_FAILURE);
      }
      cmph_io_adapter_t *source = cmph_io_nlfile_adapter(keys_fd);
        
      cmph_config_t *config = cmph_config_new(source);
      cmph_config_set_algo(config, CMPH_CHM);
      hash = cmph_new(config);
      cmph_config_destroy(config);
        
      // Destroy file adapter
      cmph_io_nlfile_adapter_destroy(source);
      fclose(keys_fd);
      
      TIME(f);
      
      eprintf(" done  [%f sec]\n", (f-s));
      
      int ret = remove(keys_str);
      if ( ret != 0 )
      {
        fprintf(stderr, "  %sWARNING%s: could not delete temporary file %s\n",
                        "\033[0;33m",keys_str,"\033[0m");
      }
        
      // 5. add ids to burst trie
      // 6. build the positions lookup table using MPHF
      
      eprintf("    (3/3) building position lookup tables ..");
      
      // positions_tbl[kmer_id] will return a pointer to an array of pairs, each pair
      // stores the sequence number and index on the sequence of the kmer_id 19-mer
      kmer_origin* positions_tbl = NULL;
      
      positions_tbl = (kmer_origin*)malloc(number_elements*sizeof(kmer_origin));
      if ( positions_tbl == NULL )
      {
          fprintf(stderr,"  %sERROR%s: could not allocate memory for positions_tbl (main(), indexdb.cpp)\n","\033[0;31m","\033[0m");
          exit(EXIT_FAILURE);
      }
            
      memset(positions_tbl, 0, number_elements*sizeof(kmer_origin));     
        
      // sequence number
      uint32_t i = 0;
      
      // reset the file pointer to the beginning of the current part
      fseek(fp,start_part,SEEK_SET);   

      uint32_t max_run = 0; //tmp
      uint32_t min_run = 1000000; //tmp 
      uint32_t total_skips = 0; //tmp 
      // 18-mer, <min range, max range>
      map<unsigned long long int,pair<uint32_t,uint32_t> > kmer_count_map; //tmp   
        
      TIME(s);
      do
      {
        vector<uint32_t> skipped_kmers; //tmp
        uint32_t length_void = 0; //tmp 

        long int start_seq = ftell(fp);
        nt = fgetc(fp);

        // scan to end of header name
        while ( nt != '\n') nt = fgetc(fp);
        
        unsigned char* myseq = new unsigned char[maxlen];
        unsigned char* myseqr = new unsigned char[maxlen];
        uint32_t _j = 0;
        len = 0;
        
        // encode each sequence using integer alphabet {0,1,2,3}
        nt = fgetc(fp);
        while ( nt != '>' && nt != EOF )
        {
          // skip line feed, carriage return or empty space in the sequence
          if ( nt != '\n' && nt != ' ' )
          {
            len++;
            // exact character
            myseq[_j++] = *(map_ptr+nt);
          }
          nt = fgetc(fp);
        }
                
        // put back the '>'
        if ( nt != EOF ) ungetc(nt,fp);
        
        // check the addition of this sequence will not overflow the maximum memory
        double estimated_seq_mem = (len-pread_gv+1)*9.5e-6;
        
        // the sequence alone is too large, it will not fit into maximum memory, skip it
        if ( estimated_seq_mem > mem ) continue;
        // the additional sequence will overflow the maximum index memory,
        // write existing index to disk and start a new index
        else if ( index_size+estimated_seq_mem > mem )
        {
          // set the character to something other than EOF
          if ( nt == EOF ) nt = 'A';
          
          // scan back to start of sequence for next index part
          fseek(fp,start_seq,SEEK_SET);
          break;
        }
        // add the additional sequence to the index
        else
        {
          index_size+=estimated_seq_mem;
        }
              
        // create a reverse sequence using the forward
        unsigned char* ptr = &myseq[len-1];
        
        for ( _j = 0; _j < len; _j++ ) myseqr[_j] = *ptr--;
        
        uint32_t kmer_key_short_f = 0;
        uint32_t kmer_key_short_r = 0;
        unsigned char* kmer_key_short_f_p = &myseq[0];
        unsigned char* kmer_key_short_r_p = &myseq[partialwin_gv+1];
        unsigned char* kmer_key_short_r_rp = &myseqr[len-partialwin_gv-1];
        
        unsigned long long int kmer_key = 0;
        unsigned char* kmer_key_ptr = &myseq[0];

        bool skip_kmer = false;
        uint32_t skip_kmer_pos = 0;
        
        // initialize the 9-mers
        for ( uint32_t j = 0; j < partialwin_gv; j++ )
        {
          (kmer_key_short_f <<= 2) |= (int)*kmer_key_short_f_p++;
          (kmer_key_short_r <<= 2) |= (int)*kmer_key_short_r_p++;
        }

        // initialize the 19-mer
        for ( uint32_t j = 0; j < pread_gv; j++ )
        {
          // skip k-mers with ambiguous nucleotides
          if ( mask_N )
          {
            if ( (int)*kmer_key_ptr == 5 )
            {
              skip_kmer = true;
              skip_kmer_pos = j;
            }
          }
          (kmer_key <<= 2) |= (int)*kmer_key_ptr++;
        } 
        
        uint32_t numwin = (len-pread_gv+interval)/interval; //TESTING
        uint32_t id = 0;
        
        uint32_t index_pos = 0; //TESTING
              
        // for all 19-mers on the sequence
        for ( uint32_t j = 0; j < numwin; j++ ) //TESTING
        {
          // no ambiguous nucleotides in the k-mer (will only be true if mask_N is set)
          if ( !skip_kmer )
          {
            // character array to hold an unsigned long long integer for CMPH
            char a[38] = {0};
            sprintf(a,"%llu",(kmer_key>>2));
            const char *key = a;
            id = cmph_search(hash, key, (cmph_uint32)strlen(key));
            
            add_id_to_burst_trie( lookup_table[kmer_key_short_f].trie_F, kmer_key_short_f_p, id );
            add_id_to_burst_trie( lookup_table[kmer_key_short_r].trie_R, kmer_key_short_r_rp, id );
            
            // if # of k-mer occurrences exceeds threshold, do not record its positions
            bool store = true;
            if (kmer_count[(kmer_key>>2)] > count_threshold)
            { 
              store = false;
              //cout << "-"; //tmp
              length_void++; //tmp
              total_skips++; //tmp

              // k-mer seen first time, add to map
              if ( kmer_count_map.find(kmer_key>>2) == kmer_count_map.end() )
              {
                pair<uint32_t,uint32_t> tmp (j, j);
                kmer_count_map[kmer_key>>2] = tmp;
              }
              // k-mer seen more than once, update range
              else
              {
                if ( j < kmer_count_map[kmer_key>>2].first ) kmer_count_map[kmer_key>>2].first = j;
                else if ( j > kmer_count_map[kmer_key>>2].second ) kmer_count_map[kmer_key>>2].second = j;
              }
            }
            else
            { 
              //cout << "X"; //tmp
              skipped_kmers.push_back(length_void);
              length_void = 0;
            }

            add_kmer_to_table( positions_tbl+id, i, index_pos, store);
          }
          
          // shift the 19-mer and 9-mers
          if ( j != numwin-1 )
          {
            for ( int shift = 0; shift < interval; shift++ )
            {
              (( kmer_key_short_f <<= 2 ) &= mask32 ) |= (int)*kmer_key_short_f_p++;
              (( kmer_key_short_r <<= 2 ) &= mask32 ) |= (int)*kmer_key_short_r_p++;

              // skip k-mers with ambiguous nucleotides
              if ( mask_N )
              {
                if ((int)*kmer_key_ptr == 5 )
                {
                  skip_kmer = true;
                  skip_kmer_pos = j+pread_gv;
                }
              }

              (( kmer_key <<= 2 ) &= mask64 ) |= (int)*kmer_key_ptr++;
              kmer_key_short_r_rp--;
              index_pos++;
            }

            // skip k-mers with ambiguous nucleotides
            if ( mask_N )
            {
              if ( skip_kmer && (j+1 > skip_kmer_pos) ) skip_kmer = false;
            }
          }
        }

        //cout << endl;

        for ( int p = 0; p < skipped_kmers.size(); p++ )
        {
          // max
          skipped_kmers[p] > max_run ? max_run = skipped_kmers[p] : max_run;
          // min
          skipped_kmers[p] < min_run ? min_run = skipped_kmers[p] : min_run;
        }

        delete [] myseq;
        delete [] myseqr;
        
        // next sequence
        i++;
            
      } while ( nt != EOF ); // for all file

      //map<unsigned long long int,pair<uint32_t,uint32_t> >::iterator it3; //tmp
      //for ( it3 = kmer_count_map.begin(); it3 != kmer_count_map.end(); it3++ )
      //  cout << it3->first << "\t[" << it3->second.first << ", " << it3->second.second << "]\n"; 

      cout << "maximum run of skipped k-mers = " << max_run << endl;
      cout << "minimum run of skipped k-mers = " << min_run << endl;

      for (int f = 0; f < number_elements; f++)
      {
        if (positions_tbl[f].size > count_threshold)
          cout << "positions_tbl[" << f << "].size = " << positions_tbl[f].size << endl;
      }

      TIME(f);
      eprintf(" done [%f sec]\n", (f-s));
      
      eprintf("    total number of sequences in this part = %d\n",i);
      
      // Destroy hash
      cmph_destroy(hash);
            
      // *********** Check ID's in Burst trie are correct *****
#ifdef see_binary_output
      uint32_t total_entries_f = 0;
      uint32_t total_entries_r = 0;
      char get_char[4] = {'A','C','G','T'};
      
      for ( int p = 0; p < (1<<lnwin_gv); p++ )
      {
        uint32_t short_kmer = p;
        char kmer_out[9] = "";
        string kmer_keep = "";
        for ( int s = 0; s < partialwin_gv; s++ )
        {
            kmer_keep.push_back((char)get_char[short_kmer&3]);
            short_kmer>>=2;
        }
        string kmer_keep_rev = "";
        for ( std::string::reverse_iterator rit=kmer_keep.rbegin();  rit != kmer_keep.rend(); ++rit )
            kmer_keep_rev.push_back(*rit);
        
        if (lookup_table[p].trie_F != NULL) traversetrie_debug( lookup_table[p].trie_F, 0, total_entries_f, kmer_keep_rev );
        //if (lookup_table[p].trie_R != NULL) traversetrie_debug( lookup_table[p].trie_R, 0, total_entries_r );
      }
#endif            
      // Load constructed index part to binary file
      
      // covert part number into a string
      stringstream prt_str;
      prt_str << part;
      string part_str = prt_str.str();

      // 1. load the kmer 'count' variable /index/kmer.dat
      ofstream oskmer ( (char*)(myfiles[newindex].second + ".kmer_" + part_str + ".dat").c_str(), ios::binary );
      eprintf("      writing kmer data to %s\n",(myfiles[newindex].second + ".kmer_" + part_str + ".dat").c_str());
      
      index_parts_stats thispart;
      thispart.start_part = start_part;
      thispart.seq_part_size = seq_part_size;
      thispart.numseq_part = numseq_part;
      index_parts_stats_vec.push_back(thispart);
      
      // the 9-mer look up tables
      for ( uint32_t j = 0; j < (uint32_t)(1<<lnwin_gv); j++ )
      {
        oskmer.write(reinterpret_cast<const char*>(&(lookup_table[j].count)),
                     sizeof(uint32_t));
      }
      
      oskmer.close();
      
      // 2. mini-burst tries
      
      // load 9-mer look-up table and mini-burst tries to /index/bursttrief.dat
      eprintf("      writing burst tries to %s\n",
              (myfiles[newindex].second + ".bursttrie_" + part_str + ".dat").c_str());
      load_index( lookup_table, (char*)(myfiles[newindex].second + ".bursttrie_" + 
                  part_str + ".dat").c_str() );
      
      // 3. 19-mer position look up tables
      ofstream ospos ( (char*)(myfiles[newindex].second + ".pos_" +
                       part_str + ".dat").c_str(), ios::binary );
      eprintf("      writing position lookup table to %s\n",
              (myfiles[newindex].second + ".pos_" + part_str + ".dat").c_str());
            
      // number of unique 19-mers
      ospos.write(reinterpret_cast<const char*>(&number_elements), sizeof(uint32_t));
      
      // the positions
      for ( uint32_t j = 0; j < number_elements; j++ )
      {
        int32_t size = positions_tbl[j].size;
        ospos.write(reinterpret_cast<const char*>(&size), sizeof(int32_t));
        if (size > -1)
          ospos.write(reinterpret_cast<const char*>(positions_tbl[j].arr), sizeof(seq_pos)*size);
      }
      ospos.close();
      
      
      // Free malloc'd memory
      // Table of unique 19-mer positions
      for ( uint32_t z = 0; z < number_elements; z++ )
      {
        if (positions_tbl[z].size > -1) free(positions_tbl[z].arr);
      }
      free(positions_tbl);

      // 9-mer look-up table and mini-burst tries
      for ( uint32_t z = 0; z < (uint32_t)(1<<lnwin_gv); z++ )
      {
        if (lookup_table[z].trie_F != NULL )
        {
          freebursttrie(lookup_table[z].trie_F);
          free(lookup_table[z].trie_F);
        }
        if (lookup_table[z].trie_R != NULL )
        {
          freebursttrie(lookup_table[z].trie_R);
          free(lookup_table[z].trie_R);
        }
      }
      
      free(lookup_table);
      
      part++;
      
  } while ( nt != EOF ); // for all index parts
        
  if ( index_size != 0 )
  {
      eprintf("      writing nucleotide distribution statistics to %s\n",(myfiles[newindex].second + ".stats").c_str());
      ofstream stats ( (char*)(myfiles[newindex].second + ".stats").c_str(), ios::binary );
      if ( !stats.good() )
      {
        fprintf(stderr,"\n  %sERROR%s: The file '%s' cannot be created: %s\n\n",
                       "\033[0;31m","\033[0m",(char*)(myfiles[newindex].second + ".stats").c_str(),
                       strerror(errno));
        exit(EXIT_FAILURE);
      }
      
      // file size for file used to build the index
      stats.write(reinterpret_cast<const char*>(&filesize), sizeof(size_t));
      
      // length of fasta file name (incl. path)
      uint32_t fasta_len = (myfiles[newindex].first).length()+1;
      stats.write(reinterpret_cast<const char*>(&fasta_len), sizeof(uint32_t));
      
      // the fasta file name (incl. path)
      stats.write(reinterpret_cast<const char*>((myfiles[newindex].first).c_str()), sizeof(char)*fasta_len);
      
      // number of sequences in the reference file
      uint32_t num_sq = sam_sq_header.size();
      
      // background frequencies, size of the database, window length
      // and number of reference sequences in a separate stats file
      double total_nt = background_freq[0] + background_freq[1] + background_freq[2] + background_freq[3];
      background_freq[0] = background_freq[0]/total_nt;
      background_freq[1] = background_freq[1]/total_nt;
      background_freq[2] = background_freq[2]/total_nt;
      background_freq[3] = background_freq[3]/total_nt;
      
      // the A/C/G/T percentage distribution
      stats.write(reinterpret_cast<const char*>(&background_freq), sizeof(double)*4);
      
      // the length of all sequences in the database
      stats.write(reinterpret_cast<const char*>(&full_len), sizeof(uint64_t));
            
      // sliding window length
      stats.write(reinterpret_cast<const char*>(&lnwin_gv), sizeof(uint32_t));
      
      numseq_gv = (int)strs/2;
      
      // number of reference sequences in the database
      stats.write(reinterpret_cast<const char*>(&numseq_gv), sizeof(uint32_t));
      
      // number of index parts
      stats.write(reinterpret_cast<const char*>(&part), sizeof(uint16_t));
      
      // information on the location and size of sequences used to build each index part
      for ( uint16_t j = 0; j < part; j++ )
      {
        stats.write(reinterpret_cast<const char*>(&index_parts_stats_vec[j]), sizeof(index_parts_stats));
      }
      
      stats.write(reinterpret_cast<const char*>(&num_sq), sizeof(uint32_t));
      
      for ( uint32_t j = 0; j < sam_sq_header.size(); j++ )
      {
        // length of the sequence id
        uint32_t len_id = sam_sq_header[j].first.length();
        stats.write(reinterpret_cast<const char*>(&len_id), sizeof(uint32_t));
        
        // the sequence id
        stats.write(reinterpret_cast<const char*>(&(sam_sq_header[j].first[0])), sizeof(char)*len_id);
        
        // the length of the sequence itself
        stats.write(reinterpret_cast<const char*>(&(sam_sq_header[j].second)), sizeof(uint32_t));
      }
      
      stats.close();
      
      eprintf("    done.\n\n");
    }

    // Free map'd memory
    fclose(fp);
        
  } //for every FASTA file, index name pair listed after --ref option
    
  return 0;
    
}//~main