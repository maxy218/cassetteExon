/*
 * =====================================================================================
 *
 *       Filename:  read_gtf.cpp
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
#include <map>
#include <string>

#include <boost/unordered_map.hpp>
#include <boost/unordered_set.hpp>

#include "class.h"
#include "common.h"
#include "const.h"
#include "gtf.h"

using namespace std;
using namespace boost;

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

//only deal with the exon annotation. CDS and start/end_codon is ignored.
// 0 -> chr; 1 -> data source; 2 -> function; 3 -> start; 4 -> end; 5 -> score; 6 -> strand; 7 -> phase; 8 -> gene id and trans id
static void get_anno_GTF(ifstream& in_anno,
    unordered_map<string, vector<isoform_anno> >& g_iso_anno_map){
  unordered_map<string, gene_info>::iterator iter_map_g_anno;

  string line;

  // key: transcript name, value: transcript annotation.
  unordered_map<string, isoform_anno> iso_anno_map;
  unordered_map<string, isoform_anno>::iterator iter_map_iso_anno;

  // the column number of GTF is <= 9
  const static unsigned int col_num_GTF = 9;
  vector<string> str_vec = vector<string>(col_num_GTF);

  while(getline(in_anno, line)){
    delimiter_ret_ref(line, '\t', col_num_GTF, str_vec);

    string gene_id = get_id_gtf(str_vec[8], "gene_id");
    string trans_id = get_id_gtf(str_vec[8], "transcript_id");

    _chr_coor start_pos, end_pos;

    // only deal with exon.
    if(str_vec[2] != "exon"){
      continue;
    }
    start_pos = atoi(str_vec[3].c_str()) - 1; // "-1" is because the GTF starts from 1, not 0 (which is refFlat style.)
    end_pos = atoi(str_vec[4].c_str()); // the end position is the same with refflat.

    // trans name and chr name are needed to identify a transcript.
    // some trans may come from different chromosome.
    string iso_chr_combined = trans_id + "\t" + str_vec[0];

    // if has been dealt before.
    if( (iter_map_iso_anno = iso_anno_map.find(iso_chr_combined)) != iso_anno_map.end()){
      iter_map_iso_anno -> second.exon_starts.push_back(start_pos);
      iter_map_iso_anno -> second.exon_ends.push_back(end_pos);
    }
    else{ // this trans has not been dealt.
      iso_anno_map[ iso_chr_combined ] = isoform_anno();

      isoform_anno & iso_anno = iso_anno_map[ iso_chr_combined ];
      iso_anno.gene_name = gene_id;
      iso_anno.name = trans_id;
      iso_anno.chrom = str_vec[0];
      iso_anno.strand = str_vec[6];

      iso_anno.exon_starts.push_back(start_pos);
      iso_anno.exon_ends.push_back(end_pos);
    }
  }

  for(iter_map_iso_anno = iso_anno_map.begin(); iter_map_iso_anno != iso_anno_map.end(); iter_map_iso_anno++){
    isoform_anno & iso_anno = iter_map_iso_anno->second;
    sort(iso_anno.exon_starts.begin(), iso_anno.exon_starts.end());
    sort(iso_anno.exon_ends.begin(), iso_anno.exon_ends.end());
    unique(iso_anno.exon_starts.begin(), iso_anno.exon_starts.end());
    unique(iso_anno.exon_ends.begin(), iso_anno.exon_ends.end());

    iso_anno.exon_cnt = iso_anno.exon_starts.size();
    iso_anno.tx_start = iso_anno.exon_starts[0];
    iso_anno.tx_end = iso_anno.exon_ends[ iso_anno.exon_cnt - 1 ];

    g_iso_anno_map[iso_anno.gene_name].push_back( iso_anno );
  }
}

void get_anno_GTF_filtered(ifstream & in_anno,
    unordered_map<string, gene_info> & map_g_anno){
  unordered_map<string, gene_info>::iterator iter_map_g_anno;
  unordered_map<string, vector<isoform_anno> > g_iso_anno_map;
  unordered_map<string, vector<isoform_anno> >::iterator iter_g_iso;

  get_anno_GTF(in_anno, g_iso_anno_map);

  for(iter_g_iso = g_iso_anno_map.begin(); iter_g_iso != g_iso_anno_map.end(); iter_g_iso++){
    //map_g_anno[(*iter_g_iso).first] = gene_info( g_iso_anno_map[(*iter_g_iso).first] );
    map_g_anno[(*iter_g_iso).first] = gene_info( iter_g_iso -> second );
  }

  iter_map_g_anno = map_g_anno.begin();
  while(iter_map_g_anno != map_g_anno.end()){
    if(!if_gene_anno_valid( (*iter_map_g_anno).second) ){
      map_g_anno.erase(iter_map_g_anno++);
    }
    else{
      iter_map_g_anno++;
    }
  }
}

void output_anno_GTF_format(const unordered_map<string, gene_info> & map_g_anno,
    ofstream & out_anno){
  unordered_map<string, gene_info>::const_iterator iter_g_map = map_g_anno.begin();
  for(; iter_g_map != map_g_anno.end(); ++iter_g_map){
    // get the reference of the genn annotation
    const gene_info & g = iter_g_map -> second;

    // output the annotation information to file isoform-by-isoform
    for(size_t i = 0; i < g.iso_num; ++i){
      for(size_t j = 0; j < g.exon_starts[i].size(); ++j){
        out_anno << g.iso_chrom[i] << "\t";
        out_anno << "unknown" << "\t";
        out_anno << "exon" << "\t";
        // remenber that we minus one on start position
        out_anno << g.exon_starts[i][j] + 1<< "\t";
        out_anno << g.exon_ends[i][j] << "\t";
        out_anno << "." << "\t";
        out_anno << g.iso_strand[i] << "\t";
        out_anno << "." << "\t"; // for exon, this field is always "."
        // the following code is to output the id parts.
        out_anno << "gene_id \"" << g.gene_name << "\";";
        out_anno << " ";
        out_anno << "gene_name \"" << g.gene_name << "\";";
        out_anno << " ";
        out_anno << "transcript_id \"" << g.iso_name[i] << "\";";
        out_anno << endl;
      }
    }
  }
}

// output the position information of all the exons of all the genes.
// we only care about the gene's information, here we will not output the information of isoforms.
// the position will be sorted by the start point.
// since there's no overlap, it's also sorted by the end point
void output_gene_exon_information(const unordered_map<string, gene_info> & map_g_anno,
    ofstream & out_anno){
    unordered_map<string, gene_info>::const_iterator iter_g_map = map_g_anno.begin();
  for(; iter_g_map != map_g_anno.end(); ++iter_g_map){
    // get the reference of the genn annotation
    const gene_info & g = iter_g_map -> second;

    const vector<vector<_chr_coor> > & starts = g.exon_starts;
    const vector<vector<_chr_coor> > & ends = g.exon_ends;
    // here we use the property of map that the items are ordered by key.
    // so here we use map, instead of unordered_map of boost
    map<_chr_coor, _chr_coor> map_start_end;
    for(size_t i = 0; i < starts.size(); ++i){
      for(size_t j = 0; j < starts[i].size(); ++j){
        map_start_end[starts[i][j]] = ends[i][j];
      }
    }

    // the following code is to output the exon information for each gene
    out_anno << g.gene_name << "\t";
    out_anno << g.chrom << "\t";
    out_anno << g.strand << "\t";
    map<_chr_coor, _chr_coor>::iterator iter_start_end = map_start_end.begin();
    for(; iter_start_end != map_start_end.end(); ++iter_start_end){
      out_anno << iter_start_end -> first << ":" << iter_start_end -> second << ",";
    }
    out_anno << endl;
  }
}

