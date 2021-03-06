/*
 * =====================================================================================
 *
 *       Filename:  anno_preproc.cpp
 *
 *    Description:  some algorithms
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
#include <iostream>
#include <string>

#include <boost/unordered_map.hpp>

#include "class.h"
#include "common.h"
#include "gtf.h"

using namespace std;
using namespace boost;

void usage(ostream& out){
  out << "anno_preproc origin_gtf_anno filtered_gtf_anno exon_info_file" << endl;
}

int main(int argc, char** argv){
  if(argc != 4){
    cerr << argv[0] << ": " << "ERROR: ";
    cerr << "invalid parameter!" << endl;
    usage(cerr);
    exit(1);
  }

  ifstream in_gtf_file(argv[1]);
  if( !in_gtf_file.is_open() ){
    cerr << argv[0] << ": " << "ERROR: ";
    cerr << "cannot open file gtf_anno_file: " << argv[1] << endl;
    exit(1);
  }

  ofstream out_filtered_gtf(argv[2]);
  if( !out_filtered_gtf.is_open() ){
    cerr << argv[0] << ": " << "ERROR: ";
    cerr << "cannot open file filtered_gtf_anno_file: " << argv[2] << endl;
    exit(1);
  }

  ofstream out_gene_exon(argv[3]);
  if( !out_gene_exon.is_open() ){
    cerr << argv[0] << ": " << "ERROR: ";
    cerr << "cannot open file exon_info_file: " << argv[3] << endl;
    exit(1);
  }

  unordered_map<string, gene_info> map_g_anno;
  get_anno_GTF_filtered(in_gtf_file, map_g_anno);
  output_anno_GTF_format(map_g_anno, out_filtered_gtf);
  output_gene_exon_information(map_g_anno, out_gene_exon);

  return 0;
}

