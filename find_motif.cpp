/*
 * =====================================================================================
 *
 *       Filename:  find_motif.cpp
 *
 *    Description:  map the motif to 5 regions
 *
 *        Version:  1.0
 *        Created:  03/15/2014 11:03:31 PM
 *       Revision:  none
 *       Compiler:  g++
 *
 *         Author:  Xinyun Ma
 *   Organization:  Tsinghua University
 *
 * =====================================================================================
 */

#include <fstream>
#include <string>
#include <vector>

#include <boost/unordered_map.hpp>
#include <boost/unordered_set.hpp>

#include "common.h"
#include "fasta.h"

using namespace std;
using namespace boost;

const static size_t MAX_MOTIF_LEN = 15;

void get_map_sf_motif(ifstream & in_motif_seq, unordered_map<string, string> & map_sf_motif,
    unordered_map<string, string> & map_motif_sf){
  string line;
  const static unsigned int col_num_motif_seq = 3;
  vector<string> fields = vector<string>(col_num_motif_seq);

  while(getline(in_motif_seq, line)){
    delimiter_ret_ref(line, '\t', col_num_motif_seq, fields);
    string & sf = fields[0];
    string & motif_seq = fields[1];

    //filter out the motif that longer than MAX_MOTIF_LEN
    if(motif_seq.length() > MAX_MOTIF_LEN){
      continue;
    }

    // deal with map_sf_motif
    if(map_sf_motif.find(sf) == map_sf_motif.end()){
      map_sf_motif[sf] = motif_seq;
    }
    else{
      map_sf_motif[sf] += ":" + motif_seq;
    }

    // deal with map_motif_sf
    if(map_motif_sf.find(motif_seq) == map_motif_sf.end()){
      map_motif_sf[motif_seq] = sf;
    }
    else{
      map_motif_sf[motif_seq] += ":" + sf;
    }
  }
/*
  unordered_map<string, string>::iterator iter;
  for(iter = map_sf_motif.begin(); iter != map_sf_motif.end(); ++iter){
    cout << iter -> first << "\t" << iter -> second << endl;
  }

  cout << endl << endl;
  for(iter = map_motif_sf.begin(); iter != map_motif_sf.end(); ++iter){
    cout << iter -> first << "\t" << iter -> second << endl;
  }
*/
} 

static char table_DNA_2_RNA[256];
void reverse_cmpl_DNA_2_RNA(string & seq){

  table_DNA_2_RNA['A'] = 'U';
  table_DNA_2_RNA['T'] = 'A';
  table_DNA_2_RNA['G'] = 'C';
  table_DNA_2_RNA['C'] = 'G';

  reverse(seq);
  size_t size = seq.size();
  for(size_t i = 0; i < size; ++i){
    seq[i] = table_DNA_2_RNA[seq[i]];
  }
}

// get the hit results.
// return by bitmap. from the highest to lowest is : UU, UD, exon, DU, DD
unsigned int get_hit_for_motif(unordered_map<string, string> & reg_RNA_seq,
    const string & exon_name, const string & motif_seq){
  unsigned int hits = 0;

  if(reg_RNA_seq[exon_name + ":UU"].find(motif_seq) != string::npos){
    hits |= (1 << 4);
  }
  if(reg_RNA_seq[exon_name + ":UD"].find(motif_seq) != string::npos){
    hits |= (1 << 3);
  }
  if(reg_RNA_seq[exon_name + ":exon"].find(motif_seq) != string::npos){
    hits |= (1 << 2);
  }
  if(reg_RNA_seq[exon_name + ":DU"].find(motif_seq) != string::npos){
    hits |= (1 << 1);
  }
  if(reg_RNA_seq[exon_name + ":DD"].find(motif_seq) != string::npos){
    hits |= 1;
  }
  return hits;
}

string trans_hits_to_string(unsigned hits){
  string hits_str = "0:0:0:0:0";
  for(int bit_pos = 4; bit_pos >= 0; --bit_pos){
    if((hits & (1 << bit_pos)) != 0){
      hits_str[8 - 2*bit_pos] = '1';
    }
  }
  return hits_str;
}

void get_motif_hit_5_region(ifstream & in_motif_seq, ifstream & in_5_region_seq,
    ofstream & out_exon_motif_hit, ofstream & out_exon_sf_hit){
  unordered_map<string, string>  map_sf_motif, map_motif_sf;
  get_map_sf_motif(in_motif_seq, map_sf_motif, map_motif_sf);

  unordered_map<string, string> reg_RNA_seq;
  get_seq_from_fasta(reg_RNA_seq, in_5_region_seq);
  // since reg_RNA_seq is RNA sequence, we need to transfer to RNA sequence
  //   using reverse_cmpl_DNA_2_RNA
  unordered_map<string, string>::iterator iter_reg_RNA_seq = reg_RNA_seq.begin();
  for(; iter_reg_RNA_seq != reg_RNA_seq.end(); ++iter_reg_RNA_seq){
    reverse_cmpl_DNA_2_RNA(iter_reg_RNA_seq -> second);
  }

  //do the motif hits. exon by exon.
  //at first, we get the exon name list.
  unordered_set<string> exon_name_set;
  iter_reg_RNA_seq = reg_RNA_seq.begin();
  for(; iter_reg_RNA_seq != reg_RNA_seq.end(); ++iter_reg_RNA_seq){
    // find the last colon
    const string & key = iter_reg_RNA_seq -> first;
    unsigned found = key.find_last_of(':');
    exon_name_set.insert(key.substr(0, found));
  }

  unordered_map<string, string>::iterator iter_sf_motif;
  unordered_set<string>::iterator iter_exon_name = exon_name_set.begin();

  out_exon_sf_hit << "exon" << "\t";
  out_exon_motif_hit << "exon" << "\t";
  for(iter_sf_motif = map_sf_motif.begin();
        iter_sf_motif != map_sf_motif.end(); ++iter_sf_motif){
    out_exon_sf_hit << iter_sf_motif -> first << "\t";
    vector<string> sf_motif_seq = delimiter(iter_sf_motif -> second, ':');
    for(size_t idx = 0; idx < sf_motif_seq.size(); ++idx){
      out_exon_motif_hit << sf_motif_seq[idx] << "\t";
    }
  }
  out_exon_sf_hit << endl;
  out_exon_motif_hit << endl;

  for(; iter_exon_name != exon_name_set.end(); ++iter_exon_name){
    const string & exon_name = *iter_exon_name;
    out_exon_motif_hit << exon_name << "\t";
    out_exon_sf_hit << exon_name << "\t";

    for(iter_sf_motif = map_sf_motif.begin();
        iter_sf_motif != map_sf_motif.end(); ++iter_sf_motif){
      unsigned int hits_for_sf = 0;
      vector<string> sf_motif_seq = delimiter(iter_sf_motif -> second, ':');
      for(size_t idx = 0; idx < sf_motif_seq.size(); ++idx){
        unsigned int hits = get_hit_for_motif(reg_RNA_seq, exon_name, sf_motif_seq[idx]);
        out_exon_motif_hit << trans_hits_to_string(hits) << "\t";
        hits_for_sf |= hits;
      }
      out_exon_sf_hit << trans_hits_to_string(hits_for_sf) << "\t";
    }
    out_exon_motif_hit << endl;
    out_exon_sf_hit << endl;
  }
}


void usage(ostream& out){
  out << "./find_motif in_motif_seq in_5_region_seq.fa out_filtered_motif out_exon_motif_hit out_exon_sf_hit." << endl;
}

int main(int argc, char** argv){
  if(argc != 6){
    cerr << argv[0] << ": " << "ERROR: ";
    cerr << "invalid parameter!" << endl;
    usage(cerr);
    exit(1);
  }

  // input: motif sequence information. 
  ifstream in_motif_seq(argv[1]);
  if( !in_motif_seq.is_open() ){
    cerr << argv[0] << ": " << "ERROR: ";
    cerr << "cannot open file in_motif_seq: " << argv[1] << endl;
    exit(1);
  }

  // input: 5 region sequences. 
  ifstream in_5_region_seq(argv[2]);
  if( !in_5_region_seq.is_open() ){
    cerr << argv[0] << ": " << "ERROR: ";
    cerr << "cannot open file in_5_region_seq: " << argv[2] << endl;
    exit(1);
  }
  
  // output: filtered motif. 
  ofstream out_filtered_motif(argv[3]);
  if( !out_filtered_motif.is_open() ){
    cerr << argv[0] << ": " << "ERROR: ";
    cerr << "cannot open file out_filtered_motif: " << argv[3] << endl;
    exit(1);
  }

  // output: exon -> motif his. 
  ofstream out_exon_motif_hit(argv[4]);
  if( !out_exon_motif_hit.is_open() ){
    cerr << argv[0] << ": " << "ERROR: ";
    cerr << "cannot open file out_exon_motif_hit: " << argv[4] << endl;
    exit(1);
  }

  // output: exon -> splicing factor his. 
  ofstream out_exon_sf_hit(argv[5]);
  if( !out_exon_sf_hit.is_open() ){
    cerr << argv[0] << ": " << "ERROR: ";
    cerr << "cannot open file out_exon_sf_hit: " << argv[5] << endl;
    exit(1);
  }

  get_motif_hit_5_region(in_motif_seq, in_5_region_seq, out_exon_motif_hit, out_exon_sf_hit);

  return 0;
}

