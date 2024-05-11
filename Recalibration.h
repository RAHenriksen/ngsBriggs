#ifndef RECALIBRATION_H
#define RECALIBRATION_H

#include <htslib/hts.h>
#include <htslib/sam.h>
#include <htslib/faidx.h>

void *tsk_all_loglike_recalibration_slave(void *dats);
void *tsk_all_loglike_recalibration_grad_slave(void *dats);
void *tsk_all_loglike_recalibration_hess_slave(void *dats);
double tsk_all_loglike_recalibration(const double *x, const void *dats);
void tsk_all_loglike_recalibration_grad(const double *x, double *y,const void *dats);
/*
double tsk_loglike_recalibration(const double *x,double **mat,int from,int to,sam_hdr_t *hdr, faidx_t *seq_ref,int len_limit, int len_min, char * model, double eps, double lambda, double delta, double delta_s, double nv,int threadid);

double all_loglike_recalibration(const double *x, const void *);



void all_loglike_recalibration_grad(const double *x, double *y,const void *);



*/
#endif
