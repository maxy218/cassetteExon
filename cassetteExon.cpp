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

#include "class.h"
#include "common.h"
#include "const.h"
#include "gtf.h"

using namespace std;
using namespace boost;

//map_iso_anno: key: position information. (chr)_(pos/neg)_(start)_(end)
//              value: isoform name
void deal_with_cassetteExon(ifstream & exon_anno_file,
    const unordered_map<string, string>& exon_gene_map,
    unordered_map<string, string>& cas_exon_gene_map){
  string line;

  const static unsigned int col_num_cassetteExon = 7;
  vector<string> fields = vector<string>(col_num_cassetteExon);

  unordered_map<string, string>::const_iterator iter_map_ex_g;
  while(getline(exon_anno_file, line)){
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

void get_inclusion_level(ifstream & expr_esti_file, 
    const unordered_map<string, string>& cas_exon_gene_map,
    const int expr_esti_choice, ofstream & inclusion_level_file){
  unordered_map<string, unordered_map<string, double> > map_gene_expr;
  unordered_map<string, double> map_gene_expr_tot;

  get_gene_expr_map(expr_esti_file, expr_esti_choice, map_gene_expr, map_gene_expr_tot);

  // get the cas exon's inclusion level.
  vector<string> id = vector<string>(2);
  unordered_map<string, string>::const_iterator iter_exon_gene;
  for(iter_exon_gene = cas_exon_gene_map.begin(); 
      iter_exon_gene != cas_exon_gene_map.end(); ++iter_exon_gene){
    vector<string> pos = delimiter(iter_exon_gene -> first, ':');
    for(size_t idx = 0; idx < pos.size(); ++idx){
      inclusion_level_file << pos[idx] << "\t";
    }

    double incl_expr = 0.0;
    vector<string> id_pairs = delimiter(iter_exon_gene -> second, '\t');
    for(size_t i = 0; i < id_pairs.size(); ++i){
      delimiter_ret_ref(id_pairs[i], ':', 2, id);
      incl_expr += map_gene_expr[id[0]][id[1]];
    }

    //  output the gene name
    inclusion_level_file << id[0] << "\t";

    if(fabs(map_gene_expr_tot[id[0]]) < EPSILON){
      inclusion_level_file << 0.0 << endl;
    }
    else{
      inclusion_level_file << incl_expr / map_gene_expr_tot[id[0]] << endl;
    }
  }
}

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

/*
template<typename T>
static int output_vector(const vector<T> & vec){
  for(size_t idx = 0; idx < vec.size(); ++idx){
    cout << vec[idx] << "\t";
  }
  cout << endl;
}

template<typename T1, typename T2>
static int output_map_key(const unordered_map<T1, T2> & data){
  typename unordered_map<T1, T2>::const_iterator iter = data.begin();
  for(; iter != data.end(); ++iter){
    cout << iter -> first << endl; 
  }
  cout << endl;
}
*/

// this function will output the 5 regions of cassette exon to file 
void get_5_regions(ifstream & gene_exons_bndr,
   const unordered_map<string, string> & cas_exon_gene_map, ofstream & seq_regions_file){
  string line;

  unordered_map<string, vector<string> > map_gene_exons;

  const static unsigned int col_num_bndr = 4;
  vector<string> fields = vector<string>(col_num_bndr);
  while(getline(gene_exons_bndr, line)){
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
    // get the upper stream exon's right boundary
    delimiter_ret_ref(exon_boundaries[bound_idx - 1], ':', col_num_boundary, boundary_vec);
    _chr_coor upper_right = atoi(boundary_vec[1].c_str());
    seq_regions_file << iter_cas_exon_gene -> first << ":" << 1 << "\t";
    seq_regions_file << upper_right << "\t" << upper_right + REGION_SIZE << endl;

    // get the middle three regions.
    _chr_coor left_boundary = atoi(fields[2].c_str());
    _chr_coor right_boundary = atoi(fields[3].c_str());
    seq_regions_file << iter_cas_exon_gene -> first << ":" << 2 << "\t";
    seq_regions_file << left_boundary - REGION_SIZE << "\t" << left_boundary << endl;
    seq_regions_file << iter_cas_exon_gene -> first << ":" << 3 << "\t";
    seq_regions_file << left_boundary << "\t" << right_boundary << endl;
    seq_regions_file << iter_cas_exon_gene -> first << ":" << 4 << "\t";
    seq_regions_file << right_boundary << "\t" << right_boundary + REGION_SIZE<< endl;

    // get the down streams exons' left boundary.
    delimiter_ret_ref(exon_boundaries[bound_idx + 1], ':', col_num_boundary, boundary_vec);
    _chr_coor down_left = atoi(boundary_vec[0].c_str());
    seq_regions_file << iter_cas_exon_gene -> first << ":" << 5 << "\t";
    seq_regions_file << down_left - REGION_SIZE << "\t" << down_left << endl;
  }
}

void usage(ostream& out){
  out << "./cassetteExon gtf_anno.gtf cassetteExon_annotation 1|2 expression_estimation_file output_inclusion_level gene_exon_boundary output_5_regions." << endl;
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
    cerr << "ERROR: " << "invalid choice of expression estimation format." << endl;
    usage(cerr);
    exit(1);
  }
  
  // get exon-gene map from the gtf file.
  ifstream gtf_anno_file(argv[1]);
  if( !gtf_anno_file.is_open() ){
    cerr << "ERROR: " << "cannot open file gtf_anno_file: " << argv[1] << endl;
    exit(1);
  }
  unordered_map<string, string> exon_gene_map;
  get_exon_gene_map_gtf(gtf_anno_file, exon_gene_map);

  // get the cassette exon - gene map from the annotation of cassette exon
  ifstream exon_anno_file(argv[2]);
  if( !exon_anno_file.is_open() ){
    cerr << "ERROR: " << "cannot open file exon_anno_file: " << argv[2] << endl;
    exit(1);
  }
  unordered_map<string, string> cas_exon_gene_map;
  deal_with_cassetteExon(exon_anno_file, exon_gene_map, cas_exon_gene_map);

  //get the expression estimation file
  // get the cassette exon - gene map from the annotation of cassette exon
  ifstream expr_esti_file(argv[4]);
  if( !expr_esti_file.is_open() ){
    cerr << "ERROR: " << "cannot open file expr_esti_file: " << argv[4] << endl;
    exit(1);
  }
  ofstream inclusion_level_file(argv[5]);
  if( !inclusion_level_file.is_open() ){
    cerr << "ERROR: " << "cannot open file inclusion_level_file: " << argv[5] << endl;
    exit(1);
  }
  get_inclusion_level(expr_esti_file, cas_exon_gene_map, expr_esti_choice, inclusion_level_file);

  
  // get the 5 regions for each cassette exon
  ifstream gene_exons_bndr(argv[6]);
  if( !gene_exons_bndr.is_open() ){
    cerr << "ERROR: " << "cannot open file expr_esti_file: " << argv[6] << endl;
    exit(1);
  }
  ofstream seq_regions_file(argv[7]);
  if( !seq_regions_file.is_open() ){
    cerr << "ERROR: " << "cannot open file seq_regions_file: " << argv[7] << endl;
    exit(1);
  }
  get_5_regions(gene_exons_bndr, cas_exon_gene_map, seq_regions_file);
}

