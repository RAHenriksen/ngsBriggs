#ifndef LIKELIHOOD_H
#define LIKELIHOOD_H

double loglike(const double *x, double * freqCT, double * freqGA, double * scaleCT, double * scaleGA);

double b_loglike(const double *x, const void *);

void loglike_grad(const double *x,double *y, double * freqCT, double * freqGA, double * scaleCT, double * scaleGA);

void loglike_hessian(const double *x, double ** z, double * freqCT, double * freqGA, double * scaleCT, double * scaleGA);

void b_loglike_grad(const double *x,double *y,const void*);

double loglike_complex3_full_b(const double *x, double * freqCT, double * freqGA, double * scaleCT, double * scaleGA, double * seqError, int& BinNum, double*& LEN, double*& freqLEN, double eps);

double loglike_complex3_full_nb(const double *x, double * freqCT, double * freqGA, double * scaleCT, double * scaleGA, double * seqError, int& BinNum, double*& LEN, double*& freqLEN, double eps);

double b_loglike_complex3_full(const double *x, const void *);

double nb_loglike_complex3_full(const double *x, const void *);

void loglike_complex3_grad_full_b(const double *x,double *y, double * freqCT, double * freqGA, double * scaleCT, double * scaleGA, double * seqError, int& BinNum, double*& LEN, double*& freqLEN, double eps);

void loglike_complex3_grad_full_nb(const double *x,double *y, double * freqCT, double * freqGA, double * scaleCT, double * scaleGA, double * seqError, int& BinNum, double*& LEN, double*& freqLEN, double eps);

void b_loglike_complex3_grad_full(const double *x,double *y,const void *);

void nb_loglike_complex3_grad_full(const double *x,double *y,const void *);

void loglike_complex3_hessian_full_b(const double *x, double ** z, double * freqCT, double * freqGA, double * scaleCT, double * scaleGA, double * seqError, int BinNum, double * LEN, double * freqLEN, double eps);

void loglike_complex3_hessian_full_nb(const double *x, double ** z, double * freqCT, double * freqGA, double * scaleCT, double * scaleGA, double * seqError, int BinNum, double * LEN, double * freqLEN, double eps);

double like_master(const double *xs,const void *);

void like_grad_master(const double *xs,double *y,const void *);

void like_hess_master(const double *xs,double **y);

double ErrorLik(char reffrag[], char frag[], int L, int seqError[]);

double PMDLik_b(char reffrag[], char frag[], int L, double lambda, double delta, double delta_s, double nv, int seqError[]);

double PMDLik_nb(char reffrag[], char frag[], int L, double lambda, double delta, double delta_s, double nv, int seqError[]);

#endif