#ifndef POSTERPROB_H
#define POSTERPROB_H
double NoPMDGivenAnc_b(char reffrag[], char frag[], int L, double lambda, double delta, double delta_s, double nv,double Tol);
double NoPMDGivenAnc_nb(char reffrag[], char frag[], int L, double lambda, double delta, double delta_s, double nv,double Tol);
double AncProb(char reffrag[], char frag[], int L, double lambda, double delta, double delta_s, double nv, uint8_t seqError[], char* model, double eps, double anc_mu, double anc_si, double mod_mu, double mod_si, int isrecal, int len_limit, int len_min,double Tol);
double PMDProb(char reffrag[], char frag[], int L, double lambda, double delta, double delta_s, double nv, uint8_t seqError[], char* model,double Tol);
bam_hdr_t* calc_pp_pmd_prob(char *refName,char *fname, char* ofname, int mapped_only,int se_only, int mapq, faidx_t *seq_ref, int len_limit, int len_min, char * model, double eps, double lambda, double delta, double delta_s, double nv, double anc_mu, double anc_si, double mod_mu, double mod_si, int isrecal, kstring_t *str_cli,double **deamRateCT,double **deamRateGA,double Tol);
#endif
