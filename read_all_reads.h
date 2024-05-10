#ifndef READ_ALL_READS_H
#define READ_ALL_READS_H
#include <vector>
#include <htslib/sam.h>
#include <htslib/faidx.h>

sam_hdr_t * read_all_reads(const char *htsname,const char *refName,std::vector<bam1_t*> &ret);

typedef struct{
  char *read;
  char *ref;
  uint8_t *qual;
  int len;
}asite;

sam_hdr_t * read_all_reads(const char *htsname,const char *refName, faidx_t *seq_ref,std::vector<asite> &ret,int len_limit);

double **read_all_reads(samFile *in,sam_hdr_t *hdr,faidx_t *seq_ref,int len_limit,double lambda,double delta,double delta_s,double nv,double Tol,int &ndim,int model);
#endif
