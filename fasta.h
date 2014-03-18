/*
 * =====================================================================================
 *
 *       Filename:  fasta.h
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

#ifndef GTF__H_INCLUDED
#define GTF__H_INCLUDED

#include <fstream>
#include <string>

#include <boost/unordered_map.hpp>

using namespace std;
using namespace boost;

void get_seq_from_fasta(unordered_map<string, string> & map_fasta_seq, ifstream & fasta_file);

#endif
