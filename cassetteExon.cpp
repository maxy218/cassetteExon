/*
 * =====================================================================================
 *
 *       Filename:  cassetteExon.cpp
 *
 *    Description:  deal with the cassette exon
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


#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>
#include <sstream>
#include <vector>

#include <boost/unordered_map.hpp>
#include <boost/unordered_set.hpp>

#include <dirent.h>

#include "class.h"
#include "common.h"
#include "const.h"
#include "gtf.h"

using namespace std;
using namespace boost;

// return -1 if not found.
// the return value is not size_t
template<typename T>
static int find_in_vector(const vector<T> & vec, const T & x){
  size_t size = vec.size();
  for(size_t idx = 0; idx < size; ++idx){
    if(vec[idx] == x){
      return idx;
    }
  }
  return -1; // if run here, there's no match result.
}

//map_iso_anno: key: position information. (chr)_(pos/neg)_(start)_(end)
//              value: isoform name
void deal_with_cassetteExon(ifstream & in_exon_anno,
    const unordered_map<string, string>& exon_gene_map,
    unordered_map<string, string>& cas_exon_gene_map){
  string line;

  const static unsigned int col_num_cassetteExon = 7;
  vector<string> fields = vector<string>(col_num_cassetteExon);

  unordered_map<string, string>::const_iterator iter_map_ex_g;
  while(getline(in_exon_anno, line)){
    delimiter_ret_ref(line, '\t', col_num_cassetteExon, fields);

    if(fields[4] != "cassetteExon"){
      continue;
    }

    string key = fields[1] + KEY_DLM + fields[6] + KEY_DLM + fields[2] + KEY_DLM + fields[3];
    if( (iter_map_ex_g = exon_gene_map.find(key)) != exon_gene_map.end()){
      cas_exon_gene_map[key] = iter_map_ex_g -> second;
    }
  }
}

// get map of gene and isoform expression level
//     key: gene name
//     value: map of isoform and expression level.
// the second hash map is: a map of gene and total expression level.
// expr_type_choice: 1 -> NURD, 2 -> cufflinks.
void get_gene_expr_map_NURD(
    ifstream & in_expr_file,
    unordered_map<string, unordered_map<string, double> > & map_gene_expr,
    unordered_map<string, double> & map_gene_expr_tot){
  string line;
  const static unsigned int col_num_NURD = 6;
  vector<string> fields = vector<string>(col_num_NURD);

  while(getline(in_expr_file, line)){
    delimiter_ret_ref(line, '\t', col_num_NURD, fields);
    string gene_id = fields[0];
    map_gene_expr[gene_id] = unordered_map<string, double>();

    size_t iso_num = atoi(fields[1].c_str());
    vector<string> trans_id = vector<string>(iso_num);
    vector<string> trans_expr = vector<string>(iso_num);
    delimiter_ret_ref(fields[3], ',', iso_num, trans_id);
    delimiter_ret_ref(fields[4], ',', iso_num, trans_expr);

    double tot_expr = 0.0;
    for(size_t i = 0; i < iso_num; ++i){
      double single_trans_expr = atof(trans_expr[i].c_str());
      map_gene_expr[gene_id][ trans_id[i] ] = single_trans_expr;
      tot_expr += single_trans_expr;
    }
    map_gene_expr_tot[gene_id] = tot_expr;
  }
}

// get map of gene and isoform expression level
//     key: gene name
//     value: map of isoform and expression level.
// the second hash map is: a map of gene and total expression level.
// expr_type_choice: 1 -> NURD, 2 -> cufflinks.
void get_gene_expr_map_cufflinks(
    ifstream & in_expr_file,
    unordered_map<string, unordered_map<string, double> > & map_gene_expr,
    unordered_map<string, double> & map_gene_expr_tot){
  string line;
  const static unsigned int col_num_cufflinks = 10;
  vector<string> fields = vector<string>(col_num_cufflinks);

  // remove the first line of cufflinks's output
  getline(in_expr_file, line);

  while(getline(in_expr_file, line)){
    delimiter_ret_ref(line, '\t', col_num_cufflinks, fields);
    string gene_id = fields[3];
    string trans_id = fields[0];
    double trans_expr = atof(fields[9].c_str());
    if(map_gene_expr.find(gene_id) == map_gene_expr.end()){
      map_gene_expr[gene_id] = unordered_map<string, double>();
      map_gene_expr[gene_id][trans_id] = trans_expr;
    }
    else{
      map_gene_expr[gene_id][trans_id] = trans_expr;
    }

    if(map_gene_expr_tot.find(gene_id) == map_gene_expr_tot.end()){
      map_gene_expr_tot[gene_id] = trans_expr;
    }
    else{
      map_gene_expr_tot[gene_id] += trans_expr;
    }
  }
}

// get map of gene and isoform expression level
//     key: gene name
//     value: map of isoform and expression level.
// the second hash map is: a map of gene and total expression level.
// expr_type_choice: 1 -> NURD, 2 -> cufflinks.
void get_gene_expr_map(
    ifstream & in_expr_file, int expr_type_choice,
    unordered_map<string, unordered_map<string, double> > & map_gene_expr,
    unordered_map<string, double> & map_gene_expr_tot){
  if(expr_type_choice == 1){
    get_gene_expr_map_NURD(in_expr_file, map_gene_expr, map_gene_expr_tot);
  }
  else if(expr_type_choice == 2){
    get_gene_expr_map_cufflinks(in_expr_file, map_gene_expr, map_gene_expr_tot);
  }
}

void get_inclusion_level(ifstream & in_expr_esti, 
    const unordered_map<string, string>& cas_exon_gene_map,
    const int expr_esti_choice, ofstream & out_inclusion_level){
  unordered_map<string, unordered_map<string, double> > map_gene_expr;
  unordered_map<string, double> map_gene_expr_tot;

  get_gene_expr_map(in_expr_esti, expr_esti_choice, map_gene_expr, map_gene_expr_tot);

  // get the cas exon's inclusion level.
  vector<string> id = vector<string>(2);
  unordered_map<string, string>::const_iterator iter_exon_gene;
  for(iter_exon_gene = cas_exon_gene_map.begin(); 
      iter_exon_gene != cas_exon_gene_map.end(); ++iter_exon_gene){
    vector<string> pos = delimiter(iter_exon_gene -> first, ':');
    for(size_t idx = 0; idx < pos.size(); ++idx){
      out_inclusion_level << pos[idx] << "\t";
    }

    double incl_expr = 0.0;
    vector<string> id_pairs = delimiter(iter_exon_gene -> second, '\t');
    for(size_t i = 0; i < id_pairs.size(); ++i){
      delimiter_ret_ref(id_pairs[i], ':', 2, id);
      incl_expr += map_gene_expr[id[0]][id[1]];
    }

    //  output the gene name
    out_inclusion_level << id[0] << "\t";

    if(fabs(map_gene_expr_tot[id[0]]) < EPSILON){
      out_inclusion_level << 0.0 << endl;
    }
    else{
      out_inclusion_level << incl_expr / map_gene_expr_tot[id[0]] << endl;
    }
  }
}

void get_inclusion_level_map(ifstream & in_expr_esti, 
    const unordered_map<string, string>& cas_exon_gene_map, int expr_esti_choice,
    unordered_map<string, vector<double> >& map_exon_incls){
  unordered_map<string, unordered_map<string, double> > map_gene_expr;
  unordered_map<string, double> map_gene_expr_tot;

  get_gene_expr_map(in_expr_esti, expr_esti_choice, map_gene_expr, map_gene_expr_tot);

  // get the cas exon's inclusion level.
  vector<string> id = vector<string>(2);
  unordered_map<string, string>::const_iterator iter_exon_gene;
  for(iter_exon_gene = cas_exon_gene_map.begin(); 
      iter_exon_gene != cas_exon_gene_map.end(); ++iter_exon_gene){
    double incl_expr = 0.0;
    double incl_lev = 0.0;
    vector<string> id_pairs = delimiter(iter_exon_gene -> second, '\t');
    for(size_t i = 0; i < id_pairs.size(); ++i){
      delimiter_ret_ref(id_pairs[i], ':', 2, id);
      incl_expr += map_gene_expr[id[0]][id[1]];
    }

    //  add the gene name to the key
    string key = id[0] + KEY_DLM + iter_exon_gene -> first;

    if(fabs(map_gene_expr_tot[id[0]]) < EPSILON){
      incl_lev = 0.0;
    }
    else{
      incl_lev = incl_expr / map_gene_expr_tot[id[0]];
    }

    if(map_exon_incls.find(key) == map_exon_incls.end()){
      map_exon_incls[key] = vector<double>(1, incl_lev);
    }
    else{
      map_exon_incls[key].push_back(incl_lev);
    }
  }
}

void output_exon_incl_level(
    const unordered_map<string, vector<double> >& map_exon_incls,
    const vector<string> & sample_names, ofstream & out_inclusion_level){
  out_inclusion_level << "exons" << "\t";
  size_t sample_size = sample_names.size();
  for(size_t idx = 0; idx < sample_size; ++idx){
    out_inclusion_level << sample_names[idx] << "\t";
  }
  out_inclusion_level << endl;

  unordered_map<string, vector<double> >::const_iterator iter_map;
  for(iter_map = map_exon_incls.begin(); iter_map != map_exon_incls.end(); ++iter_map){
    out_inclusion_level << iter_map -> first << "\t";
    const vector<double> & vec_exon_incl = iter_map -> second;
    for(size_t idx = 0; idx < vec_exon_incl.size(); ++idx){
      out_inclusion_level << vec_exon_incl[idx] << "\t";
    }
    out_inclusion_level << endl;
  }
}

// this function will output the 5 regions of cassette exon to file 
void get_5_regions(ifstream & in_gene_exons_bndr,
   const unordered_map<string, string> & cas_exon_gene_map, ofstream & out_seq_regions){
  string line;

  unordered_map<string, vector<string> > map_gene_exons;

  const static unsigned int col_num_bndr = 4;
  vector<string> fields = vector<string>(col_num_bndr);
  while(getline(in_gene_exons_bndr, line)){
    delimiter_ret_ref(line, '\t', col_num_bndr, fields);
    map_gene_exons[fields[0]] = delimiter(fields[3], ',');
  }

  const static unsigned int col_num_cas_exon = 4;
  fields = vector<string>(col_num_cas_exon);
  const static unsigned int col_num_boundary = 2;
  vector<string> boundary_vec = vector<string>(col_num_boundary);
  int bound_idx = -1;
  unordered_map<string, string>::const_iterator iter_cas_exon_gene = cas_exon_gene_map.begin();
  for(; iter_cas_exon_gene != cas_exon_gene_map.end(); ++iter_cas_exon_gene){
    delimiter_ret_ref(iter_cas_exon_gene -> first, KEY_DLM[0], col_num_cas_exon, fields);
    string boundary = fields[2] + ":" + fields[3];

    size_t pos_of_colon = (iter_cas_exon_gene -> second).find(":");
    string gene = (iter_cas_exon_gene -> second).substr(0, pos_of_colon);
    const vector<string> & exon_boundaries = map_gene_exons[gene];
    bound_idx = find_in_vector(exon_boundaries, boundary);
    if(bound_idx <= 0 || bound_idx >= exon_boundaries.size() - 1){
      continue;
    }
    string key = gene + ":" + iter_cas_exon_gene -> first;
    if(fields[1] == "+"){
      // get the upper stream exon's right boundary
      delimiter_ret_ref(exon_boundaries[bound_idx - 1], ':', col_num_boundary, boundary_vec);
      _chr_coor upper_right = atoi(boundary_vec[1].c_str());
      out_seq_regions << key << ":" << "UU" << "\t";
      out_seq_regions << fields[0] << "\t" << fields[1] << "\t";
      out_seq_regions << upper_right << "\t" << upper_right + REGION_SIZE << endl;
  
      // get the middle three regions.
      _chr_coor left_boundary = atoi(fields[2].c_str());
      _chr_coor right_boundary = atoi(fields[3].c_str());
      out_seq_regions << key << ":" << "UD" << "\t";
      out_seq_regions << fields[0] << "\t" << fields[1] << "\t";
      out_seq_regions << left_boundary - REGION_SIZE << "\t" << left_boundary << endl;
      out_seq_regions << key << ":" << "exon" << "\t";
      out_seq_regions << fields[0] << "\t" << fields[1] << "\t";
      out_seq_regions << left_boundary << "\t" << right_boundary << endl;
      out_seq_regions << key << ":" << "DU" << "\t";
      out_seq_regions << fields[0] << "\t" << fields[1] << "\t";
      out_seq_regions << right_boundary << "\t" << right_boundary + REGION_SIZE<< endl;
  
      // get the down streams exons' left boundary.
      delimiter_ret_ref(exon_boundaries[bound_idx + 1], ':', col_num_boundary, boundary_vec);
      _chr_coor down_left = atoi(boundary_vec[0].c_str());
      out_seq_regions << key << ":" << "DD" << "\t";
      out_seq_regions << fields[0] << "\t" << fields[1] << "\t";
      out_seq_regions << down_left - REGION_SIZE << "\t" << down_left << endl;
    }
    else{
      // get the upper stream exon's right boundary
      delimiter_ret_ref(exon_boundaries[bound_idx - 1], ':', col_num_boundary, boundary_vec);
      _chr_coor upper_right = atoi(boundary_vec[1].c_str());
      out_seq_regions << key << ":" << "DD" << "\t";
      out_seq_regions << fields[0] << "\t" << fields[1] << "\t";
      out_seq_regions << upper_right << "\t" << upper_right + REGION_SIZE << endl;
  
      // get the middle three regions.
      _chr_coor left_boundary = atoi(fields[2].c_str());
      _chr_coor right_boundary = atoi(fields[3].c_str());
      out_seq_regions << key << ":" << "DU" << "\t";
      out_seq_regions << fields[0] << "\t" << fields[1] << "\t";
      out_seq_regions << left_boundary - REGION_SIZE << "\t" << left_boundary << endl;
      out_seq_regions << key << ":" << "exon" << "\t";
      out_seq_regions << fields[0] << "\t" << fields[1] << "\t";
      out_seq_regions << left_boundary << "\t" << right_boundary << endl;
      out_seq_regions << key << ":" << "UD" << "\t";
      out_seq_regions << fields[0] << "\t" << fields[1] << "\t";
      out_seq_regions << right_boundary << "\t" << right_boundary + REGION_SIZE<< endl;
  
      // get the down streams exons' left boundary.
      delimiter_ret_ref(exon_boundaries[bound_idx + 1], ':', col_num_boundary, boundary_vec);
      _chr_coor down_left = atoi(boundary_vec[0].c_str());
      out_seq_regions << key << ":" << "UU" << "\t";
      out_seq_regions << fields[0] << "\t" << fields[1] << "\t";
      out_seq_regions << down_left - REGION_SIZE << "\t" << down_left << endl;
    }
  }
}

void usage(ostream& out){
  out << "./cassetteExon in_gtf_anno.gtf in_cassetteExon_annotation 1|2 in_expression_estimation_file out_inclusion_level in_gene_exon_boundary out_5_regions." << endl;
  out << "Third parameter: 1 for NURD, 2 for Cufflinks." << endl;
}

int main(int argc, char** argv){

  if(argc != 8){
    cerr << "ERROR: " << "invalid parameter!" << endl;
    usage(cerr);
    exit(1);
  }

  // 1: NURD, 2: cufflinks
  int expr_esti_choice = atoi(argv[3]);
  if(expr_esti_choice <= 0 || expr_esti_choice >= 3){
    cerr << argv[0] << ": " << "ERROR: ";
    cerr << "invalid choice of expression estimation format." << endl;
    usage(cerr);
    exit(1);
  }
  
  // get exon-gene map from the gtf file.
  ifstream in_gtf_anno(argv[1]);
  if( !in_gtf_anno.is_open() ){
    cerr << argv[0] << ": " << "ERROR: ";
    cerr << "cannot open file in_gtf_anno: " << argv[1] << endl;
    exit(1);
  }
  unordered_map<string, string> exon_gene_map;
  get_exon_gene_map_gtf(in_gtf_anno, exon_gene_map);

  // get the cassette exon - gene map from the annotation of cassette exon
  ifstream in_exon_anno(argv[2]);
  if( !in_exon_anno.is_open() ){
    cerr << argv[0] << ": " << "ERROR: ";
    cerr << "cannot open file in_exon_anno: " << argv[2] << endl;
    exit(1);
  }
  unordered_map<string, string> cas_exon_gene_map;
  deal_with_cassetteExon(in_exon_anno, exon_gene_map, cas_exon_gene_map);

  string dir_expr_esti = argv[4];
  DIR * dir_ptr; // the directory
  struct dirent * direntp; // each entry

  if( (dir_ptr = opendir(dir_expr_esti.c_str())) == NULL ){
    cerr << argv[0] << ": " << "ERROR: ";
    cerr << "cannot open directory: " << argv[4] << endl;
  }
  unordered_map<string, vector<double> > map_exon_incls;
  vector<string> sample_names;
  while( (direntp = readdir(dir_ptr)) != NULL ){
    string filename_no_dir = string(direntp -> d_name);
    if(filename_no_dir == "." || filename_no_dir == ".."){
      continue;
    }
    sample_names.push_back(filename_no_dir);
    string filename = dir_expr_esti+ "//" + filename_no_dir;
    ifstream in_single_expr_esti(filename.c_str());
    if( !in_single_expr_esti.is_open() ){
      cerr << argv[0] << ": " << "ERROR: ";
      cerr << "cannot open file expression file: " << filename << endl;
      continue;
    }
    get_inclusion_level_map(in_single_expr_esti, cas_exon_gene_map,
        expr_esti_choice, map_exon_incls);
  }
  closedir(dir_ptr);

  ofstream out_inclusion_level(argv[5]);
  if( !out_inclusion_level.is_open() ){
    cerr << argv[0] << ": " << "ERROR: ";
    cerr << "cannot open file out_inclusion_level: " << argv[5] << endl;
    exit(1);
  }
  output_exon_incl_level(map_exon_incls, sample_names, out_inclusion_level);

  // get the 5 regions for each cassette exon
  ifstream in_gene_exons_bndr(argv[6]);
  if( !in_gene_exons_bndr.is_open() ){
    cerr << argv[0] << ": " << "ERROR: ";
    cerr << "cannot open file in_gene_exons_bndr: " << argv[6] << endl;
    exit(1);
  }
  ofstream out_seq_regions(argv[7]);
  if( !out_seq_regions.is_open() ){
    cerr << argv[0] << ": " << "ERROR: ";
    cerr << "cannot open file out_seq_regions: " << argv[7] << endl;
    exit(1);
  }
  get_5_regions(in_gene_exons_bndr, cas_exon_gene_map, out_seq_regions);
}

