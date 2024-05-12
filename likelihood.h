#ifndef LIKELIHOOD_H
#define LIKELIHOOD_H
#include <vector>
#include <htslib/sam.h>
#include <htslib/faidx.h>

#include "read_all_reads.h"

extern double PhredError[255];
extern double PhredErrorAThird[255];


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
  int counter[2];//<- ncall,ncall_gradient
  int ncycle;
}wrapOne;

//below are naive versions with out the nick. 
double naive_b_loglike(const double *x, const void *);
void naive_b_loglike_grad(const double *x,double *y,const void*);

//advanced versions that consider nicks
double b_loglike_complex3_full(const double *x, const void *);
double nb_loglike_complex3_full(const double *x, const void *);
void b_loglike_complex3_grad_full(const double *x,double *y,const void *);
void nb_loglike_complex3_grad_full(const double *x,double *y,const void *);

void loglike_complex3_hessian_full_b(const double *x, double ** z, double * freqCT, double * freqGA, double * scaleCT, double * scaleGA, double * seqError, int BinNum, double * LEN, double * freqLEN, double eps,int ncycle);
void loglike_complex3_hessian_full_nb(const double *x, double ** z, double * freqCT, double * freqGA, double * scaleCT, double * scaleGA, double * seqError, int BinNum, double * LEN, double * freqLEN, double eps,int ncycle);

double ErrorLik(char reffrag[], char frag[], int L, uint8_t seqError[],int l_check);
double PMDLik_b(char reffrag[], char frag[], int L, double lambda, double delta, double delta_s, double nv, uint8_t seqError[],double Tol,int l_check);
double PMDLik_nb(char reffrag[], char frag[], int L, double lambda, double delta, double delta_s, double nv, uint8_t  seqError[],double Tol,int l_check);

#endif
