#ifndef LIKELIHOOD_H
#define LIKELIHOOD_H
#include <vector>
#include <htslib/sam.h>
#include <htslib/faidx.h>

#include "read_all_reads.h"

extern double PhredError[255];
extern double PhredErrorAThird[255];
extern int tsk_nthreads;


typedef struct{
  #if 0
  std::vector<bam1_t*> *reads;
  #endif
  std::vector<asite> *reads;
    sam_hdr_t *hdr;
    faidx_t *seq_ref;
    int len_limit;
    int len_min;
    char *model;
    double eps;
    double lambda;
    double delta;
    double delta_s;
    double nv;
    int from;
    int to;
    double llh_result;
    double *x;
    double *llh_result_grad; //Add the gradient feature
    double **llh_result_hess; //Add the hessian matrix feature
    int threadid;
    double Tol;
}tsk_struct;


typedef struct{
  double *freqCT;
  double *freqGA;
  double *scaleCT;
  double *scaleGA;
  double *seqError;
  int BinNum;
  double *Bin_Frag_len;
  double *Bin_Frag_freq;
  double Contam_eps;
}wrapOne;

//below are naive versions with out the nick. 
void loglike_grad(const double *x,double *y, double * freqCT, double * freqGA, double * scaleCT, double * scaleGA);
double loglike(const double *x, double * freqCT, double * freqGA, double * scaleCT, double * scaleGA);
double b_loglike(const double *x, const void *);



void loglike_hessian(const double *x, double ** z, double * freqCT, double * freqGA, double * scaleCT, double * scaleGA);

void b_loglike_grad(const double *x,double *y,const void*);

double loglike_complex3_full_b(const double *x, double * freqCT, double * freqGA, double * scaleCT, double * scaleGA, double * seqError, int BinNum, double* LEN, double* freqLEN, double eps);

double loglike_complex3_full_nb(const double *x, double * freqCT, double * freqGA, double * scaleCT, double * scaleGA, double * seqError, int BinNum, double* LEN, double* freqLEN, double eps);

double b_loglike_complex3_full(const double *x, const void *);

double nb_loglike_complex3_full(const double *x, const void *);

void loglike_complex3_grad_full_b(const double *x,double *y, double * freqCT, double * freqGA, double * scaleCT, double * scaleGA, double * seqError, int BinNum, double* LEN, double* freqLEN, double eps);

void loglike_complex3_grad_full_nb(const double *x,double *y, double * freqCT, double * freqGA, double * scaleCT, double * scaleGA, double * seqError, int BinNum, double* LEN, double* freqLEN, double eps);

void b_loglike_complex3_grad_full(const double *x,double *y,const void *);

void nb_loglike_complex3_grad_full(const double *x,double *y,const void *);

void loglike_complex3_hessian_full_b(const double *x, double ** z, double * freqCT, double * freqGA, double * scaleCT, double * scaleGA, double * seqError, int BinNum, double * LEN, double * freqLEN, double eps);

void loglike_complex3_hessian_full_nb(const double *x, double ** z, double * freqCT, double * freqGA, double * scaleCT, double * scaleGA, double * seqError, int BinNum, double * LEN, double * freqLEN, double eps);

double like_master(const double *xs,const void *);

void like_grad_master(const double *xs,double *y,const void *);

void like_hess_master(const double *xs,double **y);

double ErrorLik(char reffrag[], char frag[], int L, uint8_t seqError[]);

double PMDLik_b(char reffrag[], char frag[], int L, double lambda, double delta, double delta_s, double nv, uint8_t seqError[],double Tol);

double PMDLik_nb(char reffrag[], char frag[], int L, double lambda, double delta, double delta_s, double nv, uint8_t  seqError[],double Tol);

#endif
