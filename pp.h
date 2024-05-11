#ifndef POSTERPROB_H
#define POSTERPROB_H

bam_hdr_t* calc_pp_pmd_prob(char *refName,char *fname, char* ofname, int mapped_only,int se_only, int mapq, faidx_t *seq_ref, int len_limit, int len_min, int model, double eps, double lambda, double delta, double delta_s, double nv, double anc_mu, double anc_si, double mod_mu, double mod_si, int isrecal, kstring_t *str_cli,double **deamRateCT,double **deamRateGA,double Tol,int l_check);
#endif
