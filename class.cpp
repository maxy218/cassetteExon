/*
 * =====================================================================================
 *
 *       Filename:  class.cpp
 *
 *    Description:  the definitions of some classes.
 *
 *        Version:  1.0
 *        Created:  02/19/2014 10:21:42 PM
 *       Revision:  none
 *       Compiler:  g++
 *
 *         Author:  Xinyun Ma
 *   Organization:  Tsinghua University
 *
 * =====================================================================================
 */


#include <algorithm>
#include <cstdlib>
#include <list>
#include <string>
#include <vector>

#include <boost/unordered_set.hpp>

#include "common.h"
#include "const.h"
#include "class.h"
using namespace std;
using namespace boost;

//gene_info: constructor functions
gene_info::gene_info(){} // do nothing. Don't caculate based on the object initialized by default constructor.

//assignment member-by-member
//  when we are assigning the container-type members, maybe smart-pointer is a better choice. 
gene_info::gene_info(const gene_info& g){
  iso_num = g.iso_num;

  gene_name = g.gene_name;
  iso_name = g.iso_name; // smart pointer maybe is a better choice. 
  chrom = g.chrom;
  iso_chrom = g.iso_chrom; // smart pointer maybe is a better choice. 
  strand = g.strand;
  iso_strand = g.iso_strand; 
  exon_starts = g.exon_starts;
  exon_ends = g.exon_ends;

  is_anno_valid = g.is_anno_valid;
}

//default constructor
isoform_anno::isoform_anno():
  gene_name (""),
  name (""),
  chrom (""),
  strand (""),
  tx_start (0),
  tx_end (0),
  cds_start (0),
  cds_end (0),
  exon_cnt (0),
  exon_starts ( vector<_chr_coor>(0) ),
  exon_ends ( vector<_chr_coor>(0) )
{}

gene_info::gene_info(const vector<isoform_anno>& iso_vec){
  vector<isoform_anno>::const_iterator iter_gene = iso_vec.begin();

  gene_name = (*iter_gene).gene_name;

  iter_gene = iso_vec.begin();
  while(iter_gene != iso_vec.end()){
    iso_name.push_back((*iter_gene).name);
    chrom = (*iter_gene).chrom;
    iso_chrom.push_back(chrom);
    strand = (*iter_gene).strand;
    iso_strand.push_back(strand);
    exon_starts.push_back((*iter_gene).exon_starts);
    exon_ends.push_back((*iter_gene).exon_ends);

    iter_gene++;
  }
  iso_num = iso_name.size();
}

// return true if ovelap. false if NOT overlap
static bool exon_overlap(const vector<vector<int> > & starts, const vector<vector<int> > & ends){
  // here we use the property of map that the items are ordered by key.
  // so here we use map, instead of unordered_map of boost
  map<_chr_coor, _chr_coor> map_start_end;
  for(size_t i = 0; i < starts.size(); ++i){
    for(size_t j = 0; j < starts[i].size(); ++j){
      // if the exon hasn't existed, insert it
      if(map_start_end.find(starts[i][j]) == map_start_end.end()){
        map_start_end[starts[i][j]] = ends[i][j];
      }
      // if the start of this exon has existed but the end is not this exon's end, it's invalid.
      else if(map_start_end[starts[i][j]] != ends[i][j]){
        return true;
      }
    }
  }

  // check whether previous ends larger than the next starts.
  map<_chr_coor, _chr_coor>::iterator iter_start_end = map_start_end.begin();
  int pre_end = 0;
  for(; iter_start_end != map_start_end.end(); ++iter_start_end){
    if(iter_start_end->first < pre_end){
      return true;
    }
    pre_end = iter_start_end -> second;
  }
  return false;
}


//the following work is to change the return type from bool to int, different return value represent different error type
bool if_gene_anno_valid(const gene_info& gene){
  //  if there're duplicate isoforms
  //  sorted list may be not efficient. Maybe hash is faster!
  //  hash: just judge whether every isoform name is new, having no duplicate
  if(gene.iso_num > 1){
    unordered_set<string> set_iso;
    unordered_set<string>::iterator iter_set_iso;
    for(size_t i = 0; i < gene.iso_num; i++){
      if((iter_set_iso = set_iso.find(gene.iso_name[i])) != set_iso.end()){
        return false;
      }
      else{
        set_iso.insert(gene.iso_name[i]);
      }
    }
  }

  //if having different orient
  if(gene.iso_num > 1){
    for(int i = 0; i < gene.iso_num - 1; i++){
      if(gene.iso_strand[i] != gene.iso_strand[i+1]){
        return false;
      }
    }
  }

  //if from different chrome
  if(gene.iso_num > 1){
    for(int i = 0; i < gene.iso_num - 1; i++){
      if(gene.iso_chrom[i] != gene.iso_chrom[i+1]){
        return false;
      }
    }
  }

  //if the exons are overlapped. it's new for cassette exon project.
  //nurd doesn't filter this case. NURD use "exon splitting instead"
  if(exon_overlap(gene.exon_starts, gene.exon_ends)){
    // if the exon_overlap return true, it means it overlaped, so it's invalid.
    return false;
  }

  return true;
}
