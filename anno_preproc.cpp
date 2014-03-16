#include <fstream>
#include <iostream>
#include <string>

#include <boost/unordered_map.hpp>

#include "class.h"
#include "common.h"
#include "parsing_gtf.h"

using namespace std;
using namespace boost;

void usage(ostream& out){
  out << "gtf_filter origin_gtf_anno filtered_gtf_anno exon_info_file" << endl;
}

int main(int argc, char** argv){
  if(argc != 4){
    cerr << "ERROR: " << "invalid parameter!" << endl;
    usage(cerr);
    exit(1);
  }

  ifstream gtf_file(argv[1]);
  if( !gtf_file.is_open() ){
    cerr << "ERROR: " << "cannot open file gtf_anno_file: " << argv[1] << endl;
    exit(1);
  }

  ofstream out_filtered_gtf(argv[2]);
  if( !gtf_file.is_open() ){
    cerr << "ERROR: " << "cannot open file filtered_gtf_anno_file: " << argv[2] << endl;
    exit(1);
  }

  ofstream out_exon(argv[3]);
  if( !gtf_file.is_open() ){
    cerr << "ERROR: " << "cannot open file exon_info_file: " << argv[3] << endl;
    exit(1);
  }

  unordered_map<string, gene_info> map_g_anno;
  get_anno_GTF_filtered(gtf_file, map_g_anno);
  output_anno_GTF_format(map_g_anno, out_filtered_gtf);
  output_gene_exon_information(map_g_anno, out_exon);

  return 0;
}

