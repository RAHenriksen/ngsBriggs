#ifndef RECALIBRATION_H
#define RECALIBRATION_H

#include <htslib/hts.h>
#include <htslib/sam.h>
#include <htslib/faidx.h>

#include "read_all_reads.h"

double tsk_loglike_recalibration(const double *x,double **mat,int from,int to,sam_hdr_t *hdr, faidx_t *seq_ref,int len_limit, int len_min, char * model, double eps, double lambda, double delta, double delta_s, double nv,int threadid);

double all_loglike_recalibration(const double *x, const void *);
double tsk_all_loglike_recalibration(const double *x, const void *dats);
void *tsk_all_loglike_recalibration_slave(void *dats);
#if 0 
void tsk_loglike_recalibration_grad(const double *x, double *y, std::vector<bam1_t*> *reads,int from, int to,sam_hdr_t *hdr, faidx_t *seq_ref,int len_limit, int len_min, char * model, double eps, double lambda, double delta, double delta_s, double nv);
#endif
void tsk_loglike_recalibration_grad(const double *x, double *y, std::vector<asite> &reads,int from, int to,sam_hdr_t *hdr, faidx_t *seq_ref,int len_limit, int len_min, char * model, double eps, double lambda, double delta, double delta_s, double nv);
void all_loglike_recalibration_grad(const double *x, double *y,const void *);
void tsk_all_loglike_recalibration_grad(const double *x, double *y,const void *dats);
void *tsk_all_loglike_recalibration_grad_slave(void *dats);
#if 0 
void tsk_loglike_recalibration_hess(const double *x, double **y, std::vector<bam1_t*> *reads,int from, int to,sam_hdr_t *hdr, faidx_t *seq_ref,int len_limit, int len_min, char * model, double eps, double lambda, double delta, double delta_s, double nv);
#endif
void tsk_loglike_recalibration_hess(const double *x, double **y, std::vector<asite> &reads,int from, int to,sam_hdr_t *hdr, faidx_t *seq_ref,int len_limit, int len_min, char * model, double eps, double lambda, double delta, double delta_s, double nv);
void *tsk_all_loglike_recalibration_hess_slave(void *dats);

#endif
