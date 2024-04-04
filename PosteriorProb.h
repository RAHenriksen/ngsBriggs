#ifndef POSTERPROB_H
#define POSTERPROB_H
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
#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Eigenvalues>
#include <ctime>

double NoPMDGivenAnc_b(char reffrag[], char frag[], int L, double lambda, double delta, double delta_s, double nv);

double NoPMDGivenAnc_nb(char reffrag[], char frag[], int L, double lambda, double delta, double delta_s, double nv);

double AncProb(char reffrag[], char frag[], int L, double lambda, double delta, double delta_s, double nv, int seqError[], char* model, double eps, double anc_mu, double anc_si, double mod_mu, double mod_si, int isrecal, int len_limit, int len_min);

double PMDProb(char reffrag[], char frag[], int L, double lambda, double delta, double delta_s, double nv, int seqError[], char* model);

bam_hdr_t* CalPostPMDProb(char *refName,char *fname, const char* chromname, const char* bedname, char* ofname, char* olik, int mapped_only,int se_only, int mapq, faidx_t *seq_ref, int len_limit, int len_min, char * model, double eps, double lambda, double delta, double delta_s, double nv, double anc_mu, double anc_si, double mod_mu, double mod_si, int isrecal, std::string s);

#endif