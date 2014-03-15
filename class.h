/*
 * =====================================================================================
 *
 *     Filename:  class.h
 *
 *  Description:  the definitions of some classes.
 *
 *    Version:  1.0
 *    Created:  02/19/2014 10:21:42 PM
 *     Revision:  none
 *     Compiler:  g++
 *
 *     Author:  Xinyun Ma
 *   Organization:  Tsinghua University
 *
 * =====================================================================================
 */


#ifndef CLASS_H_INCLUDED
#define CLASS_H_INCLUDED

#include <string>
#include <vector>

#include "const.h"
using namespace std;

struct isoform_anno{
public:
  string gene_name;
  string name;
  string chrom;
  string strand;
  _chr_coor tx_start;
  _chr_coor tx_end;
  _chr_coor cds_start;
  _chr_coor cds_end;
  _chr_coor exon_cnt;
  vector<_chr_coor> exon_starts;
  vector<_chr_coor> exon_ends;

  isoform_anno();
};

struct gene_info{
public:
  int iso_num;

  //here, one elem in vector represent one isoform
  string gene_name;
  vector<string> iso_name; //maybe there are multiple isoforms from one gene
  string chrom;
  vector<string> iso_chrom; //which chromosomes are the isoforms from. 
                            //maybe different isoform comes from different chromosomes, which is invalid here.
  string strand;
  vector<string> iso_strand; //which chromosomes are the isoforms from.
                             //maybe different isoform comes from different chromosomes, which is invalid here.

  vector<vector<_chr_coor> > exon_starts; // exon boundaries.
  vector<vector<_chr_coor> > exon_ends; // exon boundaries.

  int is_anno_valid;

  //following are some constructor functions.
  gene_info();
  gene_info(const gene_info& g);
  gene_info(const vector<isoform_anno>& iso_vec);
};

//the following work is to change the return type from bool to int, different return value represent different error type
bool if_gene_anno_valid(const gene_info& gene);

#endif // CLASS_H_INCLUDED
