/*
 * =====================================================================================
 *
 *       Filename:  gtf.h
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


#ifndef GTF__H_INCLUDED
#define GTF__H_INCLUDED

#include <fstream>
#include <string>

#include <boost/unordered_map.hpp>

using namespace std;
using namespace boost;

string itoa(int i);

// get gene or transcript id from the last field of gtf.
string get_id_gtf(const string& field, const string& id_type);

static void get_anno_GTF(ifstream& in_anno,
    unordered_map<string, vector<isoform_anno> >& g_iso_anno_map);

void get_anno_GTF_filtered(ifstream & in_anno,
    unordered_map<string, gene_info> & map_g_anno);

void output_anno_GTF_format(const unordered_map<string, gene_info> & map_g_anno,
    ofstream & out_anno);

// output the position information of all the exons of all the genes.
// we only care about the gene's information, here we will not output the information of isoforms.
// the position will be sorted by the start point.
// since there's no overlap, it's also sorted by the end point
void output_gene_exon_information(const unordered_map<string, gene_info> & map_g_anno,
    ofstream & out_anno);

//This function deal with gtf file.
//exon_gene_map: key: position information. (chr)_(pos/neg)_(start)_(end)
//              value: gene name
void get_exon_gene_map_gtf(ifstream & gtf_anno_file, unordered_map<string, string>& exon_gene_map);

#endif
