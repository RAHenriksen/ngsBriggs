#ifndef RECALIBRATION_H
#define RECALIBRATION_H
#include <iostream>
#include <string> // Add this line at the top of Recalibration.h
#include <cstring>
#include <cstdlib>
#include <cstdio>
#include <htslib/hts.h>
#include <htslib/sam.h>
#include <htslib/kstring.h>
#include <htslib/faidx.h>
#include <zlib.h>
#include <cmath>
#include <iomanip>

#include <ctime>
using namespace std ;

//The log-likelihood for recalibration the ancient prob
double loglike_recalibration(const double *x, char *refName,char *fname, const char* chromname, const char* bedname,int mapped_only,int se_only, int mapq, faidx_t *seq_ref,int len_limit, int len_min, char * model, double eps, double lambda, double delta, double delta_s, double nv, std::string s);

//The log-likelihood for recalibration the ancient probs
double tsk_loglike_recalibration(const double *x, std::vector<bam1_t *> *reads,int from,int to,sam_hdr_t *hdr, faidx_t *seq_ref,int len_limit, int len_min, char * model, double eps, double lambda, double delta, double delta_s, double nv,int threadid);

double All_loglike_recalibration(const double *x, const void *);

double tsk_All_loglike_recalibration(const double *x, const void *dats);

void *tsk_All_loglike_recalibration_slave(void *dats);

void loglike_recalibration_grad(const double *x, double *y, char *refName,char *fname, const char* chromname, const char* bedname,int mapped_only,int se_only, int mapq, faidx_t *seq_ref,int len_limit, int len_min, char * model, double eps, double lambda, double delta, double delta_s, double nv, std::string s);

void tsk_loglike_recalibration_grad(const double *x, double *y, std::vector<bam1_t*> *reads,int from, int to,sam_hdr_t *hdr, faidx_t *seq_ref,int len_limit, int len_min, char * model, double eps, double lambda, double delta, double delta_s, double nv);

void All_loglike_recalibration_grad(const double *x, double *y,const void *);

void tsk_All_loglike_recalibration_grad(const double *x, double *y,const void *dats);

void *tsk_All_loglike_recalibration_grad_slave(void *dats);

void tsk_loglike_recalibration_hess(const double *x, double **y, std::vector<bam1_t*> *reads,int from, int to,sam_hdr_t *hdr, faidx_t *seq_ref,int len_limit, int len_min, char * model, double eps, double lambda, double delta, double delta_s, double nv);

void *tsk_All_loglike_recalibration_hess_slave(void *dats);

#endif
