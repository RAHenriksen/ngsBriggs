/*
  For splitting reads between ancient and damage
  we have a likelihood function that needs to loops over read summary statistics
  return is double[3][nreads]={l_err,l_anc,readlen};
 */
#ifndef READ_ALL_READS_H
#define READ_ALL_READS_H
#include <htslib/sam.h>
#include <htslib/faidx.h>

double **read_all_reads(samFile *in,sam_hdr_t *hdr,faidx_t *seq_ref,int len_limit,double lambda,double delta,double delta_s,double nv,double Tol,int &ndim,int model);
#endif
