#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>
#include <sstream>
#include <vector>

#include <boost/unordered_map.hpp>

#include "common.h"

using namespace std;
using namespace boost;

typedef int _chr_coor;
const string KEY_DLM = ":";
const string VAL_DLM = ":";
const double EPSILON = 1.0e-8;

string itoa(int i){
  stringstream ss (stringstream::in | stringstream::out);
  ss.str("");
  ss << i;
  return ss.str();
}

// get gene or transcript id from the last field of gtf.
string get_id_gtf(const string& field, const string& id_type){
  size_t idx = field.find(id_type);
  size_t left = idx + id_type.size();
  size_t len = field.size();
  while(left < len && field[left] != '\"'){
    ++left;
  }

  size_t right = left + 1;
  while(right < len && field[right] != '\"'){
    ++right;
  }
 
  return field.substr(left + 1, right - left - 1); 
}

//This function deal with gtf file.
//exon_gene_map: key: position information. (chr)_(pos/neg)_(start)_(end)
//              value: gene name
void get_exon_gene_map_gtf(ifstream & gtf_anno_file, unordered_map<string, string>& exon_gene_map){
  string line;

  _chr_coor start_pos, end_pos;

  // the column number of GTF is <= 9
  const static unsigned int col_num_GTF = 9;
  vector<string> fields = vector<string>(col_num_GTF);

  while(getline(gtf_anno_file, line)){
    delimiter_ret_ref(line, '\t', col_num_GTF, fields);
    // only deal with exon.
    if(fields[2] != "exon"){
      continue;
    }

    string gene_id = get_id_gtf(fields[8], "gene_id");
    string trans_id = get_id_gtf(fields[8], "transcript_id");
    
    start_pos = atoi(fields[3].c_str()) - 1; // "-1" is because the GTF starts from 1, not 0 (which is refFlat style.)
    end_pos = atoi(fields[4].c_str()); // the end position is the same with refflat.

    string key = fields[0] + KEY_DLM + fields[6] + KEY_DLM + itoa(start_pos) + KEY_DLM + itoa(end_pos);

    if(exon_gene_map.find(key) == exon_gene_map.end()){
      exon_gene_map[key] = gene_id + VAL_DLM + trans_id;
    }
    else{
      exon_gene_map[key] += "\t" + gene_id + VAL_DLM + trans_id;
    }
  }
}

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
      //cout << iter_map_ex_g -> first << "\t" << iter_map_ex_g -> second << endl;
      cas_exon_gene_map[key] = iter_map_ex_g -> second;
    }
  }
}

void usage(ostream& out){
  out << "CaExon gtf_anno.gtf cassetteExon.txt" << endl;
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
    const int expr_esti_choice){
  unordered_map<string, unordered_map<string, double> > map_gene_expr;
  unordered_map<string, double> map_gene_expr_tot;

  get_gene_expr_map(expr_esti_file, expr_esti_choice, map_gene_expr, map_gene_expr_tot);

  // get the cas exon's inclusion level.
  vector<string> id = vector<string>(2);
  unordered_map<string, string>::const_iterator iter_exon_gene;
  for(iter_exon_gene = cas_exon_gene_map.begin(); 
      iter_exon_gene != cas_exon_gene_map.end(); ++iter_exon_gene){
    cout << iter_exon_gene -> first << "\t";

    double incl_expr = 0.0;
    vector<string> id_pairs = delimiter(iter_exon_gene -> second, '\t');
    for(size_t i = 0; i < id_pairs.size(); ++i){
      delimiter_ret_ref(id_pairs[i], ':', 2, id);
      incl_expr += map_gene_expr[id[0]][id[1]];
    }
    if(fabs(map_gene_expr_tot[id[0]]) < EPSILON){
      cout << 0.0 << endl;
    }
    else{
      cout << incl_expr / map_gene_expr_tot[id[0]] << endl;
    }
  }


}

int main(int argc, char** argv){
  if(argc != 5){
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
    cerr << "ERROR: " << "cannot open file expr_esti_file: " << argv[2] << endl;
    exit(1);
  }
  get_inclusion_level(expr_esti_file, cas_exon_gene_map, expr_esti_choice);

}










































