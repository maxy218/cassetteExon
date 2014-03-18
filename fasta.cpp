/*
 * =====================================================================================
 *
 *       Filename:  fasta.cpp
 *
 *    Description:  some interface to deal with fasta file.
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

#include <boost/unordered_map.hpp>

#include "common.h"
#include "fasta.h"

using namespace std;
using namespace boost;

void get_seq_from_fasta(unordered_map<string, string> & map_fasta_seq, ifstream & fasta_file){
  string line;
  string seq_name;
  while(getline(fasta_file, line)){
    if(line.size() == 0){
      continue;
    }
    if(line[0] == '>'){
      seq_name = line.substr(1, line.size() - 1);
      map_fasta_seq[seq_name] = "";
    }
    else{
      map_fasta_seq[seq_name] += line;
    }
  }
}

