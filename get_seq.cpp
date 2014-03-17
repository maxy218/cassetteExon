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

// region: tab delimited. 
// 4 columns in total. 1: chr, 2: strand, 3: left, 4: right
string get_region(const unordered_map<string, string> & map_chr_seq, const string & region){
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

void output_map(unordered_map<string, string> & map_data){
  unordered_map<string, string>::iterator iter = map_data.begin();
  for(; iter != map_data.end(); ++iter){
    cout << iter -> first << endl;
    cout << iter -> second << endl;
  }
}

void usage(ostream& out){
  out << "./get_seq chr_file_directory exon_boundary" << endl;
}

int main(int argc, char** argv){
  if(argc != 8){
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
/*
  output_map(map_chr_seq);
*/

  

  return 0;
}









