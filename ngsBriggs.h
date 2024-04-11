#ifndef NGSBRIGGS_H
#define NGSBRIGGS_H


extern int VERBOSE;
extern int SIG_COND;
extern int tsk_nthreads;

typedef struct{
    std::vector<bam1_t*> *reads;
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
}tsk_struct;

extern tsk_struct *my_tsk_struct;


//MAXLENGTH = 5;
typedef unsigned char uchar;
extern int nreads_per_pos;
typedef struct{
    char bp;
    int offset;
} mdField;

extern int number;

extern int STRLENS;

extern int* Frag_len;
extern double* Frag_freq;

extern double Contam_eps; //Contamination rate defined as the proportion of the number of the contaminated reads
extern double MAX0;
extern double MIN0;

extern int nproc1;//number of reads processed not used
extern char refToChar[256];
extern char toIndex[4][4];
extern char com[5];
extern unsigned char toDIndex[16][16];
extern int MAXLENGTH;
extern double** mm5p, **mm3p, **dmm5p, **dmm3p;

extern double **deamRateCT;
extern double **deamRateGA;

extern double PhredError[255];

extern double PhredErrorAThird[255]; 

extern double Tol; // Tolerance

static char DUMMYCHAR='#';

extern double l_check;

// likelihood global
//extern double *freqCT, *freqGA, *scaleCT, *scaleGA, *seqError;
extern int ncalls;
extern int ncalls_grad;

// I finally decide to use the accurate but the most complex model to check my inference
extern int MAXORDER; //This can be adjusted for different tolerance of errors.
//Check the code again

extern int BinNum;
extern double* Bin_Frag_len, *Bin_Frag_freq;

// OTHER DEFINES VALUES

extern char *refName2, *fname2, *model2;
extern const char* chromname2, * bedname2;
extern int mapped_only2, se_only2, mapq2, len_limit2, len_min2;
extern double eps2, lambda2, delta2, delta_s2, nv2;
extern std::string s2;
extern faidx_t *seq_ref2;

extern std::vector<nuclist> comp_nuc_llik;

int main(int argc, char **argv);

#endif
