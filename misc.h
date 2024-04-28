#ifndef MISC_H
#define MISC_H

#include <cmath>
#include <string>
#include <vector>
#include <cstdlib>
#include <ctime>
#include <htslib/hts.h>
#include <htslib/sam.h>
#include <htslib/kstring.h>
#include <htslib/faidx.h>
#include <zlib.h>

double NormalCDF(double x);
double NormalPDF(double x);
double NormalPDF_grad(double x);
double NormalINC(double y, double x, double x_max, double x_min);
double NormalINC_grad_mu(double y, double x, double x_max, double x_min, double mu, double sigma);
double NormalINC_grad_si(double y, double x, double x_max, double x_min,double mu, double sigma);
double NormalINC_hess_mu2(double y, double x, double x_max, double x_min,double mu, double sigma);
double NormalINC_hess_si2(double y, double x, double x_max, double x_min,double mu, double sigma);
double NormalINC_hess_mu_si(double y, double x, double x_max, double x_min,double mu, double sigma);
void FragArrayBin(int number, int &BinNum, int*& Length, double *& Freq, double*& BinLength, double*& BinFreq);
void FragArrayReader(int len_limit, int& number, int*& Length, double *& Freq, const char* filename,int STRLENS);
int tabreader(char *tabname,int STRLENS,double** mm5p, double **mm3p);
void parse_tabledata(const char* filename,double** Table,int STRLENS);
void parse_tabledata2(const char* filename,double** Table);
int bamreader(char *fname, const char* chromname,faidx_t * seq_ref, int len_limit, int &len_min,int *Frag_len, double *Frag_freq,int &number,double **mm5p,double **mm3p);
void wrapperwithref(const bam1_t   * b,const bam_hdr_t  *hdr, char myread[512], char myref[512],faidx_t *seq_ref);
void CaldeamRate_b(double lambda, double delta, double delta_s, double nu, int len_limit,double **deamRateCT,double **deamRateGA);
void CaldeamRate_nb(double lambda, double delta, double delta_s, double nu, int len_limit,double **deamRateCT,double **deamRateGA);
void parse_sequencingdata1(char *refName,char *fname,const char* chromname,int mapped_only,int se_only,int mapq, faidx_t *seq_ref,int len_limit, int & len_min,int *Frag_len, double *Frag_freq,int &number,double** mm5p, double **mm3p);
#endif
