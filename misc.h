#ifndef MISC_H
#define MISC_H

#include <iostream>
#include <cmath>
#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Eigenvalues>
#include <string>
#include <sstream>
#include <fstream>
#include <vector>
#include <cstdlib>
#include <ctime>
#include <htslib/hts.h>
#include <htslib/sam.h>
#include <htslib/kstring.h>
#include <htslib/faidx.h>
#include <zlib.h>

typedef unsigned char uchar;

typedef struct{
    double nuclik[4];
    size_t pos;
    uchar chrid;
    char* chr;
}nuclist;

double NormalCDF(double x);

double NormalPDF(double x);

double NormalPDF_grad(double x);

double NormalINC(double y, double x, double x_max, double x_min);

double NormalINC_grad_mu(double y, double x, double x_max, double x_min, double mu, double sigma);

double NormalINC_grad_si(double y, double x, double x_max, double x_min,double mu, double sigma);

double NormalINC_hess_mu2(double y, double x, double x_max, double x_min,double mu, double sigma);

double NormalINC_hess_si2(double y, double x, double x_max, double x_min,double mu, double sigma);

double NormalINC_hess_mu_si(double y, double x, double x_max, double x_min,double mu, double sigma);

bool cmp(nuclist x,nuclist y);

void merge(std::vector<nuclist> &x_sort);

void FragArrayBin(int number, int &BinNum, int*& Length, double *& Freq, double*& BinLength, double*& BinFreq);

void FragArrayReader(int len_limit, int& number, int*& Length, double *& Freq, const char* filename);

int tabreader(char *tabname);

void parse_tabledata(const char* filename,double** Table);

void parse_tabledata2(const char* filename,double** Table);

int bamreader(char *fname, const char* chromname,const char* bedname, faidx_t * seq_ref, int len_limit, int &len_min);

void wrapperwithref(const bam1_t   * b,const bam_hdr_t  *hdr, char myread[512], char myref[512],faidx_t *seq_ref);

void CaldeamRate_b(double lambda, double delta, double delta_s, double nu, int len_limit);

void CaldeamRate_nb(double lambda, double delta, double delta_s, double nu, int len_limit);

void parse_sequencingdata1(char *refName,char *fname,const char* chromname, const char* bedname, int mapped_only,int se_only,int mapq, faidx_t *seq_ref,int len_limit, int & len_min);

void Calnuclik(char myread[], kstring_t *kstr, char* chromname, uchar chrid, bam1_t *b, double PostProb);

#endif