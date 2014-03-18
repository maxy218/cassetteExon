/*
 * =====================================================================================
 *
 *       Filename:  get_seq.cpp
 *
 *    Description:  get the sequence of 5 region.
 *
 *        Version:  1.0
 *        Created:  02/25/2014 23:48:09 PM
 *       Revision:  none
 *       Compiler:  g++
 *
 *         Author:  Xinyun Ma
 *   Organization:  Tsinghua University
 *
 * =====================================================================================
 */

#include <cstdlib>
#include <cctype>
#include <iostream>
#include <string>
#include <vector>

#include <dirent.h>

#include <boost/unordered_map.hpp>

#include "common.h"

using namespace std;
using namespace boost;

void to_upper(string & seq){
  size_t size = seq.size();
  for(int i = 0; i < size; ++i){
    seq[i] = toupper(seq[i]);
  }
}

/*
template<typename T>
void swap(T & a, T & b){
  T & c = a;
  a = b;
  b = c;
}
*/

void reverse(string & seq){
  size_t size = seq.size();
  if(size == 0){
    return;
  }
  size_t i = 0, j = size - 1;
  while(i < j){
    swap(seq[i++], seq[j--]);
  }
}

static char table[256];

void reverse_cmpl(string & seq){

  table['A'] = 'T';
  table['T'] = 'A';
  table['G'] = 'C';
  table['C'] = 'G';

  reverse(seq);
  size_t size = seq.size();
  for(size_t i = 0; i < size; ++i){
    seq[i] = table[seq[i]];
  }
}

void get_chr_seq(unordered_map<string, string> & map_chr_seq, ifstream & fasta_file){
   string line;
   string chr_name;
   while(getline(fasta_file, line)){
     if(line.size() == 0){
       continue;
     }
     if(line[0] == '>'){
       chr_name = line.substr(1, line.size() - 1);
       map_chr_seq[chr_name] = "";
     }
     else{
       map_chr_seq[chr_name] += line;
     }
   }
}

// region: tab delimited. 
// 4 columns in total. 1: chr, 2: strand, 3: left, 4: right
string get_one_region(const unordered_map<string, string> & map_chr_seq, const string & region){
  const static unsigned int col_num_reg = 4;
  vector<string> fields = vector<string>(col_num_reg);
  delimiter_ret_ref(region, '\t', col_num_reg, fields);

  unordered_map<string, string>::const_iterator iter_chr_seq = map_chr_seq.find(fields[0]);
  if(iter_chr_seq == map_chr_seq.end()){
    return "";
  }

  int left = atoi(fields[2].c_str());
  int right = atoi(fields[3].c_str());
  string seq = (iter_chr_seq -> second).substr(left, right - left);
  to_upper(seq);

  if(fields[1] == "-"){
    // reverse complementary 
    reverse_cmpl(seq);
  }
  return seq;
}

void get_regs_seq(const unordered_map<string, string> & map_chr_seq, ifstream & in_reg_bndr,
    ofstream & out_reg_seq){
  string line;
  const static unsigned int col_num_reg_bndr = 5;
  vector<string> fields = vector<string>(col_num_reg_bndr);
  unordered_map<string, string>::const_iterator iter_chr_seq;
  while(getline(in_reg_bndr, line)){
    delimiter_ret_ref(line, '\t', col_num_reg_bndr, fields);
    iter_chr_seq = map_chr_seq.find(fields[1]);
    if(iter_chr_seq == map_chr_seq.end()){
      continue;
    }

    int left = atoi(fields[3].c_str());
    int right = atoi(fields[4].c_str());
    string seq = (iter_chr_seq -> second).substr(left, right - left);
    to_upper(seq);

    if(fields[2] == "-"){
      // reverse complementary
      reverse_cmpl(seq);
    }
    out_reg_seq << ">" << fields[0] << endl;
    out_reg_seq << seq << endl;
  }
}

void output_map(unordered_map<string, string> & map_data){
  unordered_map<string, string>::iterator iter = map_data.begin();
  for(; iter != map_data.end(); ++iter){
    cout << iter -> first << endl;
    cout << iter -> second << endl;
  }
}

void usage(ostream& out){
  out << "./get_seq chr_file_directory in_reg_bndr out_region_seq" << endl;
}

int main(int argc, char** argv){
  if(argc != 4){
    cerr << "ERROR: " << "invalid parameter!" << endl;
    usage(cerr);
    exit(1);
  }

  string chr_seq_dir = argv[1];
  DIR * dir_ptr; // the directory
  struct dirent * direntp; // each entry

  if( (dir_ptr = opendir(chr_seq_dir.c_str())) == NULL ){
    cerr << "cannot open directory: " << argv[1] << endl;
  }
  unordered_map<string, string> map_chr_seq;
  while( (direntp = readdir(dir_ptr)) != NULL ){
    string filename = chr_seq_dir + "//" + string(direntp -> d_name);
    if(filename.size() <= 3 || filename.substr(filename.size() - 3, 3) != ".fa"){
      continue;
    }
    ifstream fasta_file(filename.c_str());
    if( !fasta_file.is_open() ){
      cerr << "ERROR: " << "cannot open file gtf_anno_file: " << filename << endl;
      continue;
    }
    get_chr_seq(map_chr_seq, fasta_file);
  }
  closedir(dir_ptr);

  ifstream in_reg_bndr(argv[2]);
  if( !in_reg_bndr.is_open() ){
    cerr << "ERROR: " << "cannot open file in_reg_bndr: " << argv[2] << endl;
    exit(1);
  }
  ofstream out_reg_seq(argv[3]);
  if( !out_reg_seq.is_open() ){
    cerr << "ERROR: " << "cannot open file out_reg_seq: " << argv[3] << endl;
    exit(1);
  }
  get_regs_seq(map_chr_seq, in_reg_bndr, out_reg_seq);  

  return 0;
}

