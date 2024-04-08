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
#include <getopt.h>
#include <iostream>
#include "profile.h"
#include "bfgs.h"
#include "htslib/bgzf.h"
#include "briggs_writer.h"
#include "read_all_reads.h"
#include <array>  //Rasmus

#include "misc.h"
#include "Recalibration.h"
#include "Likelihood.h"
#include "PosteriorProb.h"
#include "ngsBriggs_cli.h"
#include "ngsBriggs.h"

// definining all global variables used across multiple scripts
int VERBOSE = 1;
int SIG_COND = 1;

int STRLENS = 4096;
int* Frag_len = new int[STRLENS];
double* Frag_freq = new double[STRLENS];

double MAX0 = 1-1e-8;
double MIN0 = 1e-8;

int nproc1 = 0;//number of reads processed not used
int MAXORDER = MAXLENGTH;

char *nuc = strdup("ACGTN");

double PhredError[255] = {1.00000000E+00,7.94328235E-01,6.30957344E-01,5.01187234E-01,3.98107171E-01,3.16227766E-01,2.51188643E-01,1.99526231E-01,1.58489319E-01,1.25892541E-01,1.00000000E-01,7.94328235E-02,6.30957344E-02,5.01187234E-02,3.98107171E-02,3.16227766E-02,2.51188643E-02,1.99526231E-02,1.58489319E-02,1.25892541E-02,1.00000000E-02,7.94328235E-03,6.30957344E-03,5.01187234E-03,3.98107171E-03,3.16227766E-03,2.51188643E-03,1.99526231E-03,1.58489319E-03,1.25892541E-03,1.00000000E-03,7.94328235E-04,6.30957344E-04,5.01187234E-04,3.98107171E-04,3.16227766E-04,2.51188643E-04,1.99526231E-04,1.58489319E-04,1.25892541E-04,1.00000000E-04,7.94328235E-05,6.30957344E-05,5.01187234E-05,3.98107171E-05,3.16227766E-05,2.51188643E-05,1.99526231E-05,1.58489319E-05,1.25892541E-05,1.00000000E-05,7.94328235E-06,6.30957344E-06,5.01187234E-06,3.98107171E-06,3.16227766E-06,2.51188643E-06,1.99526231E-06,1.58489319E-06,1.25892541E-06,1.00000000E-06,7.94328235E-07,6.30957344E-07,5.01187234E-07,3.98107171E-07,3.16227766E-07,2.51188643E-07,1.99526231E-07,1.58489319E-07,1.25892541E-07,1.00000000E-07,7.94328235E-08,6.30957344E-08,5.01187234E-08,3.98107171E-08,3.16227766E-08,2.51188643E-08,1.99526231E-08,1.58489319E-08,1.25892541E-08,1.00000000E-08,7.94328235E-09,6.30957344E-09,5.01187234E-09,3.98107171E-09,3.16227766E-09,2.51188643E-09,1.99526231E-09,1.58489319E-09,1.25892541E-09,1.00000000E-09,7.94328235E-10,6.30957344E-10,5.01187234E-10,3.98107171E-10,3.16227766E-10,2.51188643E-10,1.99526231E-10,1.58489319E-10,1.25892541E-10,1.00000000E-10,7.94328235E-11,6.30957344E-11,5.01187234E-11,3.98107171E-11,3.16227766E-11,2.51188643E-11,1.99526231E-11,1.58489319E-11,1.25892541E-11,1.00000000E-11,7.94328235E-12,6.30957344E-12,5.01187234E-12,3.98107171E-12,3.16227766E-12,2.51188643E-12,1.99526231E-12,1.58489319E-12,1.25892541E-12,1.00000000E-12,7.94328235E-13,6.30957344E-13,5.01187234E-13,3.98107171E-13,3.16227766E-13,2.51188643E-13,1.99526231E-13,1.58489319E-13,1.25892541E-13,1.00000000E-13,7.94328235E-14,6.30957344E-14,5.01187234E-14,3.98107171E-14,3.16227766E-14,2.51188643E-14,1.99526231E-14,1.58489319E-14,1.25892541E-14,1.00000000E-14,7.94328235E-15,6.30957344E-15,5.01187234E-15,3.98107171E-15,3.16227766E-15,2.51188643E-15,1.99526231E-15,1.58489319E-15,1.25892541E-15,1.00000000E-15,7.94328235E-16,6.30957344E-16,5.01187234E-16,3.98107171E-16,3.16227766E-16,2.51188643E-16,1.99526231E-16,1.58489319E-16,1.25892541E-16,1.00000000E-16,7.94328235E-17,6.30957344E-17,5.01187234E-17,3.98107171E-17,3.16227766E-17,2.51188643E-17,1.99526231E-17,1.58489319E-17,1.25892541E-17,1.00000000E-17,7.94328235E-18,6.30957344E-18,5.01187234E-18,3.98107171E-18,3.16227766E-18,2.51188643E-18,1.99526231E-18,1.58489319E-18,1.25892541E-18,1.00000000E-18,7.94328235E-19,6.30957344E-19,5.01187234E-19,3.98107171E-19,3.16227766E-19,2.51188643E-19,1.99526231E-19,1.58489319E-19,1.25892541E-19,1.00000000E-19,7.94328235E-20,6.30957344E-20,5.01187234E-20,3.98107171E-20,3.16227766E-20,2.51188643E-20,1.99526231E-20,1.58489319E-20,1.25892541E-20,1.00000000E-20,7.94328235E-21,6.30957344E-21,5.01187234E-21,3.98107171E-21,3.16227766E-21,2.51188643E-21,1.99526231E-21,1.58489319E-21,1.25892541E-21,1.00000000E-21,7.94328235E-22,6.30957344E-22,5.01187234E-22,3.98107171E-22,3.16227766E-22,2.51188643E-22,1.99526231E-22,1.58489319E-22,1.25892541E-22,1.00000000E-22,7.94328235E-23,6.30957344E-23,5.01187234E-23,3.98107171E-23,3.16227766E-23,2.51188643E-23,1.99526231E-23,1.58489319E-23,1.25892541E-23,1.00000000E-23,7.94328235E-24,6.30957344E-24,5.01187234E-24,3.98107171E-24,3.16227766E-24,2.51188643E-24,1.99526231E-24,1.58489319E-24,1.25892541E-24,1.00000000E-24,7.94328235E-25,6.30957344E-25,5.01187234E-25,3.98107171E-25,3.16227766E-25,2.51188643E-25,1.99526231E-25,1.58489319E-25,1.25892541E-25,1.00000000E-25,7.94328235E-26,6.30957344E-26,5.01187234E-26,3.98107171E-26};

double PhredErrorAThird[255] = {3.33333333E-01,2.64776078E-01,2.10319115E-01,1.67062411E-01,1.32702390E-01,1.05409255E-01,8.37295477E-02,6.65087438E-02,5.28297731E-02,4.19641804E-02,3.33333333E-02,2.64776078E-02,2.10319115E-02,1.67062411E-02,1.32702390E-02,1.05409255E-02,8.37295477E-03,6.65087438E-03,5.28297731E-03,4.19641804E-03,3.33333333E-03,2.64776078E-03,2.10319115E-03,1.67062411E-03,1.32702390E-03,1.05409255E-03,8.37295477E-04,6.65087438E-04,5.28297731E-04,4.19641804E-04,3.33333333E-04,2.64776078E-04,2.10319115E-04,1.67062411E-04,1.32702390E-04,1.05409255E-04,8.37295477E-05,6.65087438E-05,5.28297731E-05,4.19641804E-05,3.33333333E-05,2.64776078E-05,2.10319115E-05,1.67062411E-05,1.32702390E-05,1.05409255E-05,8.37295477E-06,6.65087438E-06,5.28297731E-06,4.19641804E-06,3.33333333E-06,2.64776078E-06,2.10319115E-06,1.67062411E-06,1.32702390E-06,1.05409255E-06,8.37295477E-07,6.65087438E-07,5.28297731E-07,4.19641804E-07,3.33333333E-07,2.64776078E-07,2.10319115E-07,1.67062411E-07,1.32702390E-07,1.05409255E-07,8.37295477E-08,6.65087438E-08,5.28297731E-08,4.19641804E-08,3.33333333E-08,2.64776078E-08,2.10319115E-08,1.67062411E-08,1.32702390E-08,1.05409255E-08,8.37295477E-09,6.65087438E-09,5.28297731E-09,4.19641804E-09,3.33333333E-09,2.64776078E-09,2.10319115E-09,1.67062411E-09,1.32702390E-09,1.05409255E-09,8.37295477E-10,6.65087438E-10,5.28297731E-10,4.19641804E-10,3.33333333E-10,2.64776078E-10,2.10319115E-10,1.67062411E-10,1.32702390E-10,1.05409255E-10,8.37295477E-11,6.65087438E-11,5.28297731E-11,4.19641804E-11,3.33333333E-11,2.64776078E-11,2.10319115E-11,1.67062411E-11,1.32702390E-11,1.05409255E-11,8.37295477E-12,6.65087438E-12,5.28297731E-12,4.19641804E-12,3.33333333E-12,2.64776078E-12,2.10319115E-12,1.67062411E-12,1.32702390E-12,1.05409255E-12,8.37295477E-13,6.65087438E-13,5.28297731E-13,4.19641804E-13,3.33333333E-13,2.64776078E-13,2.10319115E-13,1.67062411E-13,1.32702390E-13,1.05409255E-13,8.37295477E-14,6.65087438E-14,5.28297731E-14,4.19641804E-14,3.33333333E-14,2.64776078E-14,2.10319115E-14,1.67062411E-14,1.32702390E-14,1.05409255E-14,8.37295477E-15,6.65087438E-15,5.28297731E-15,4.19641804E-15,3.33333333E-15,2.64776078E-15,2.10319115E-15,1.67062411E-15,1.32702390E-15,1.05409255E-15,8.37295477E-16,6.65087438E-16,5.28297731E-16,4.19641804E-16,3.33333333E-16,2.64776078E-16,2.10319115E-16,1.67062411E-16,1.32702390E-16,1.05409255E-16,8.37295477E-17,6.65087438E-17,5.28297731E-17,4.19641804E-17,3.33333333E-17,2.64776078E-17,2.10319115E-17,1.67062411E-17,1.32702390E-17,1.05409255E-17,8.37295477E-18,6.65087438E-18,5.28297731E-18,4.19641804E-18,3.33333333E-18,2.64776078E-18,2.10319115E-18,1.67062411E-18,1.32702390E-18,1.05409255E-18,8.37295477E-19,6.65087438E-19,5.28297731E-19,4.19641804E-19,3.33333333E-19,2.64776078E-19,2.10319115E-19,1.67062411E-19,1.32702390E-19,1.05409255E-19,8.37295477E-20,6.65087438E-20,5.28297731E-20,4.19641804E-20,3.33333333E-20,2.64776078E-20,2.10319115E-20,1.67062411E-20,1.32702390E-20,1.05409255E-20,8.37295477E-21,6.65087438E-21,5.28297731E-21,4.19641804E-21,3.33333333E-21,2.64776078E-21,2.10319115E-21,1.67062411E-21,1.32702390E-21,1.05409255E-21,8.37295477E-22,6.65087438E-22,5.28297731E-22,4.19641804E-22,3.33333333E-22,2.64776078E-22,2.10319115E-22,1.67062411E-22,1.32702390E-22,1.05409255E-22,8.37295477E-23,6.65087438E-23,5.28297731E-23,4.19641804E-23,3.33333333E-23,2.64776078E-23,2.10319115E-23,1.67062411E-23,1.32702390E-23,1.05409255E-23,8.37295477E-24,6.65087438E-24,5.28297731E-24,4.19641804E-24,3.33333333E-24,2.64776078E-24,2.10319115E-24,1.67062411E-24,1.32702390E-24,1.05409255E-24,8.37295477E-25,6.65087438E-25,5.28297731E-25,4.19641804E-25,3.33333333E-25,2.64776078E-25,2.10319115E-25,1.67062411E-25,1.32702390E-25,1.05409255E-25,8.37295477E-26,6.65087438E-26,5.28297731E-26,4.19641804E-26,3.33333333E-26,2.64776078E-26,2.10319115E-26,1.67062411E-26,1.32702390E-26};

double Tol = 1.0E-8; // Tolerance

double l_check = 15;
int ncalls =0;
int ncalls_grad =0;

std::vector<nuclist> comp_nuc_llik;

tsk_struct *my_tsk_struct = NULL;

int tsk_nthreads = -1;

double **deamRateCT;
double **deamRateGA;
int number;
double *freqCT, *freqGA, *scaleCT, *scaleGA, *seqError;
double Contam_eps; //Contamination rate defined as the proportion of the number of the contaminated reads

int BinNum = -1;
double* Bin_Frag_len = new double[STRLENS];
double* Bin_Frag_freq = new double[STRLENS];

char *refName2 = NULL; // or initialize with appropriate value
char *fname2 = NULL; // or initialize with appropriate value
char *model2 = NULL; // or initialize with appropriate value
const char* chromname2, *bedname2;
int mapped_only2, se_only2, mapq2, len_limit2, len_min2;
double eps2, lambda2, delta2, delta_s2, nv2;

faidx_t *seq_ref2;
std::string s2;

// defining our main ngsBriggs function
int main(int argc, char **argv){
    argStruct *mypars = NULL;
    if(argc==1||(argc==2&&(strcasecmp(argv[1],"--help")==0||strcasecmp(argv[1],"-h")==0))){
        HelpPage(stderr);
        return 0;
    }
    
    std::string s = "./ngsBriggs ";

    for(int i=0;i<argc;i++)
    {  s += " ";
        s += argv[i];
    }
    
    //printf("The resulting string is: %s\n", s.c_str());
    mypars = pars_briggs(argc,argv);


    tsk_nthreads = mypars->nthread;
    my_tsk_struct = new tsk_struct[tsk_nthreads];
    assert(mypars);
    char* fname = mypars->hts;
    char* tabname = mypars->tab;
    const char* fastafile = mypars->ref;
    const char* lenname = mypars->len;
    const char* bedname = mypars->bed;
    const char* chromname = mypars->chr;
    char* model = mypars->model; // takes a value of b or nb.
    Contam_eps = mypars->eps; // Modern contamination rate \in [0,1)
    int len_limit = 150;
    
    faidx_t *seq_ref = NULL;
    seq_ref = fai_load(fastafile);
    fprintf(stderr,"\t-> fasta load \n");
    
    mm5p = (double**) malloc(MAXLENGTH * sizeof(double*));
    mm3p = (double**) malloc(MAXLENGTH * sizeof(double*));
    for (int i = 0; i < MAXLENGTH-1; i++){
        mm5p[i] =(double *) malloc(16 * sizeof(double));
        mm3p[i] =(double *) malloc(16 * sizeof(double));
    }
    mm5p[MAXLENGTH-1] =(double *) malloc(16 * sizeof(double));
    mm3p[MAXLENGTH-1] =(double *) malloc(16 * sizeof(double));
    for (int i=0; i<MAXLENGTH-1;i++){
        for (int j=0; j<16;j++){
            mm5p[i][j]=0;
            mm3p[i][j]=0;
        }
    }
    for (int j=0; j<16;j++){
        mm5p[MAXLENGTH-1][j]=0;
        mm3p[MAXLENGTH-1][j]=0;
    }
    
    freqCT = (double*) malloc(2*MAXLENGTH * sizeof(double));
    freqGA = (double*) malloc(2*MAXLENGTH * sizeof(double));
    scaleCT = (double*) malloc(2*MAXLENGTH * sizeof(double));
    scaleGA = (double*) malloc(2*MAXLENGTH * sizeof(double));
    seqError = (double*) malloc(2*MAXLENGTH * sizeof(double));
    int len_min;
    if ((fname != NULL && tabname != NULL) || (fname == NULL && tabname == NULL)){
        fprintf(stdout,"Please provide either a bamfile or a table!\n");
        return 0;
    }else if (fname != NULL){
        fprintf(stderr,"Loading the bamfile...\n");
        bamreader(fname,chromname,bedname,seq_ref,len_limit,len_min);
        //cout<<"Minimum length is "<<len_min<<"\n";
    }else if(tabname != NULL && lenname != NULL){
        fprintf(stderr,"Loading the table file...\n");
        tabreader(tabname);
    }else{
        fprintf(stdout,"Please provide a fragment length distribution file, if table file is provided!\n");
        return 0;
    }
    len_min = len_min > 30 ? len_min : 30;
    double max5=0;
    double max3=0;
    double maxall=0;
    for (int i=0; i<MAXLENGTH;i++){
        // Double-check this part
        // scaleCT[i] = mm5p[i][5]+mm5p[i][7];
        scaleCT[i] = mm5p[i][4]+mm5p[i][5]+mm5p[i][6]+mm5p[i][7];
        // scaleCT[2*MAXLENGTH-1-i] = mm3p[i][5]+mm3p[i][7];
        scaleCT[2*MAXLENGTH-1-i] = mm3p[i][4]+mm3p[i][5]+mm3p[i][6]+mm3p[i][7];
        // scaleGA[i] = mm3p[i][8]+mm3p[i][10];
        scaleGA[i] = mm3p[i][8]+mm3p[i][9]+mm3p[i][10]+mm3p[i][11];
        // scaleGA[2*MAXLENGTH-1-i] = mm5p[i][8]+mm5p[i][10];
        scaleGA[2*MAXLENGTH-1-i] = mm5p[i][8]+mm5p[i][9]+mm5p[i][10]+mm5p[i][11];
        
        // Overall sequencing errors are position specific, estimated by 1 - [N(AA)+N(TT)]/[N(AA)+N(AC)+N(AG)+N(AT)+N(TA)+N(TC)+N(TG)+N(TT)]
        seqError[i] = 1 - (mm5p[i][0]+mm5p[i][15])/(mm5p[i][0]+mm5p[i][1]+mm5p[i][2]+mm5p[i][3]+mm5p[i][12]+mm5p[i][13]+mm5p[i][14]+mm5p[i][15]);
        seqError[2*MAXLENGTH-1-i] = 1 - (mm3p[i][0]+mm3p[i][15])/(mm3p[i][0]+mm3p[i][1]+mm3p[i][2]+mm3p[i][3]+mm3p[i][12]+mm3p[i][13]+mm3p[i][14]+mm3p[i][15]);
        //cout << "Error "<<seqError[i]<<" "<<seqError[2*MAXLENGTH-1-i]<<"\n";
        freqCT[i] = mm5p[i][7]/scaleCT[i];
        freqCT[2*MAXLENGTH-1-i] = mm3p[i][7]/scaleCT[2*MAXLENGTH-1-i];
        freqGA[i] = mm3p[i][8]/scaleGA[i];
        freqGA[2*MAXLENGTH-1-i] = mm5p[i][8]/scaleGA[2*MAXLENGTH-1-i];
        max5 = max(max5,max(scaleCT[i],scaleCT[2*MAXLENGTH-1-i]));
        max3 = max(max3,max(scaleGA[i],scaleGA[2*MAXLENGTH-1-i]));
    }
    maxall = max(max5,max3);
    //    cout<<"max5 is "<<max5<<", max3 is "<<max3<<".\n";
    for (int i=0; i<MAXLENGTH;i++){
        scaleCT[i] = scaleCT[i]/maxall;
        scaleGA[i] = scaleGA[i]/maxall;
        scaleCT[2*MAXLENGTH-1-i] = scaleCT[2*MAXLENGTH-1-i]/maxall;
        scaleGA[2*MAXLENGTH-1-i] = scaleGA[2*MAXLENGTH-1-i]/maxall;
        //cout<<scaleCT[i]<<" "<<scaleGA[i]<<"\n";
    }
    
    
    fprintf(stderr,"\nThe misincorporting matrix is as follows:\n");
    fprintf(stderr,"Dir.\tPos.\tFreqCT\tFreqGA\n");
    for (int i=0; i<5;i++){
        fprintf(stderr,"5'\t%d\t",i+1);
        fprintf(stderr,"%f\t%f\n",freqCT[i],freqGA[2*MAXLENGTH-1-i]);
    }
    for (int i=0; i<5;i++){
        fprintf(stderr,"3'\t%d\t",i+1);
        fprintf(stderr,"%f\t%f\n",freqCT[2*MAXLENGTH-1-i],freqGA[i]);
    }
    
    //Inference
    clock_t t=clock();
    time_t t2=time(NULL);
    fprintf(stderr,"\nThe inference is under progress...\n");
    //double *x = new double[3];
    double lbd[3] = {0,0,0};
    double ubd[3] = {1,1,1};
    int nbd[3] = {2,2,2};
    double invec1[3] = {0.2,0.1,0.1};
    ncalls=0;
    ncalls_grad = 0;
    double withgrad = findmax_bfgs(3,invec1,NULL,b_loglike,b_loglike_grad,lbd,ubd,nbd,-1);
    
    double lbd1[4] = {1e-8,1e-8,1e-8,1e-8};
    double ubd1[4] = {1-1e-8,1-1e-8,1-1e-8,1-1e-8};
    int nbd1[4] = {2,2,2,2};
    double **z2;
    z2 = (double**) malloc(4 * sizeof(double*));
    for (int i = 0; i < 4; i++){
        z2[i] =(double *) malloc(4 * sizeof(double));
    }
    
    if (lenname != NULL){
        fprintf(stderr,"The provided fragment length distribution file is used!\n");
        FragArrayReader(len_limit, number, Frag_len, Frag_freq, lenname);
    }else{
        fprintf(stderr,"The fragment length distribution is calculated from the bam file!\n");
    }
    //BinNum = -1;
    FragArrayBin(number, BinNum, Frag_len, Frag_freq, Bin_Frag_len, Bin_Frag_freq);
    double invec2[4] = {invec1[0],invec1[1],invec1[2],0.01};
    ncalls=0;
    ncalls_grad = 0;
    double withgrad2 = 0;
    if (model != NULL){
        if (!strcasecmp("b",model)){
            //fprintf(stderr,"%s\n","The chosen model is a biotin model.");
            withgrad2 = findmax_bfgs(4,invec2,NULL,b_loglike_complex3_full,b_loglike_complex3_grad_full,lbd1,ubd1,nbd1,-1);
            loglike_complex3_hessian_full_b(invec2, z2, freqCT, freqGA, scaleCT, scaleGA, seqError, BinNum, Bin_Frag_len, Bin_Frag_freq, Contam_eps);
        }else if(!strcasecmp("nb",model)){
            //fprintf(stderr,"%s\n","The chosen model is non-biotin model.");
            withgrad2 = findmax_bfgs(4,invec2,NULL,nb_loglike_complex3_full,nb_loglike_complex3_grad_full,lbd1,ubd1,nbd1,-1);
            loglike_complex3_hessian_full_nb(invec2, z2, freqCT, freqGA, scaleCT, scaleGA, seqError, BinNum, Bin_Frag_len, Bin_Frag_freq, Contam_eps);
        }else{
            fprintf(stderr,"Please provide a meaningful deamination model for further calculations.\n");
            return -1;
        }
    }else{
        fprintf(stderr,"Please provide a deamination model for further calculations.\n");
        return -1;
    }
    double stdvec2[4];
    for (int i=0;i<4;i++){
        stdvec2[i] = sqrt(z2[i][i]/maxall);
    }
    
    fprintf(stderr,
            "\n"
            "\t[ALL done] cpu-time used =  %.4f sec\n"
            "\t[ALL done] walltime used =  %.4f sec\n"
            ,(float)(clock() - t) / CLOCKS_PER_SEC, (float)(time(NULL) - t2));
    if (model != NULL){
        if (!strcasecmp("b",model)){
            fprintf(stderr,"%s","The chosen model is a biotin model, ");
        }else if(!strcasecmp("nb",model)){
            fprintf(stderr,"%s","The chosen model is non-biotin model, ");
        }else{
            fprintf(stderr,"Please provide a meaningful deamination model for further calculations.\n");
            return -1;
        }
    }else{
        fprintf(stderr,"Please provide a deamination model for further calculations.\n");
        return -1;
    }
    fprintf(stderr,"the inferred parameters are given as follows:\n");
    //fprintf(stderr,"with grad val: %f=(%f,%f,%f,%f) nfunctioncalls: %d ngradientcalls: %d\n",withgrad2,invec2[0],invec2[1],invec2[2],invec2[3],ncalls,ncalls_grad);
    fprintf(stderr,"lambda: %f (%f,%f), delta: %f (%f,%f), delta_s: %f (%f,%f), nu: %f (%f,%f), nfunctioncalls: %d ngradientcalls: %d, llh: %f.\n", invec2[0], invec2[0]-1.96*stdvec2[0], invec2[0]+1.96*stdvec2[0] , invec2[1], invec2[1]-1.96*stdvec2[1], invec2[1]+1.96*stdvec2[1], invec2[2], invec2[2]-1.96*stdvec2[2], invec2[2]+1.96*stdvec2[2], invec2[3], invec2[3]-1.96*stdvec2[3], invec2[3]+1.96*stdvec2[3], ncalls,ncalls_grad, withgrad2);
    
    // Output count table
    char *otabname = mypars->otab;
    if (otabname!=NULL){
        BGZF *fp = NULL;
        fp = bgzf_open(otabname,"wb");
        kstring_t *kstr = new kstring_t;
        kstr->s = NULL;
        kstr->l = kstr->m = 0;
        ksprintf(kstr,"%s\n","The misincorporting matrix is as follows:");
        ksprintf(kstr,"%s\n","Dir.\tPos.\tAA\tAC\tAG\tAT\tCA\tCC\tCG\tCT\tGA\tGC\tGG\tGT\tTA\tTC\tTG\tTT");
        for (int i=0; i<MAXLENGTH;i++){
            ksprintf(kstr,"%s\t%d\t","5'",i+1);
            for (int j=0; j<16;j++){
                ksprintf(kstr,"%f\t",mm5p[i][j]);
            }
            ksprintf(kstr,"%f\t%f\n",freqCT[i],freqGA[2*MAXLENGTH-1-i]);
        }
        for (int i=0; i<MAXLENGTH;i++){
            ksprintf(kstr,"%s\t%d\t","3'",i+1);
            for (int j=0; j<16;j++){
                ksprintf(kstr,"%f\t",mm3p[i][j]);
            }
            ksprintf(kstr,"%f\t%f\n",freqCT[2*MAXLENGTH-1-i],freqGA[i]);
        }
        //my_bgzf_write(fp,kstr->s,kstr->l);
        assert(bgzf_write(fp,kstr->s,kstr->l)==kstr->l);

        kstr->l = 0;
        bgzf_close(fp);
    }
    
    //Output length distribution file
    char *olenname = mypars->olen;
    if (olenname!=NULL){
        BGZF *fp = NULL;
        fp = bgzf_open(olenname,"wb");
        kstring_t *kstr = new kstring_t;
        kstr->s = NULL;
        kstr->l = kstr->m = 0;
        ksprintf(kstr,"%s\n","The raw length distribution is as follows:");
        ksprintf(kstr,"%s%d\n","The number of length values is ",number);
        for (int k=0;k<number;k++){
            ksprintf(kstr,"%d\t%f\n",Frag_len[k],Frag_freq[k]);
        }
        ksprintf(kstr,"%s\n","The binned length distribution is as follows:");
        ksprintf(kstr,"%s%d\n","The number of binned length values is ",BinNum);
        for (int k=0;k<BinNum;k++){
            ksprintf(kstr,"%f\t%f\n",Bin_Frag_len[k],Bin_Frag_freq[k]);
        }
        //my_bgzf_write(fp,kstr->s,kstr->l);
        assert(bgzf_write(fp,kstr->s,kstr->l)==kstr->l);
        kstr->l = 0;
        bgzf_close(fp);
    }
    
    // Output inferred parameters
    char* oinfname = mypars->oinf;
    if (oinfname!=NULL){
        BGZF *fp = bgzf_open(oinfname,"wb");
        kstring_t *kstr = new kstring_t;
        kstr->s = NULL;
        kstr->l = kstr->m = 0;
        if (model != NULL){
            if (!strcasecmp("b",model)){
                ksprintf(kstr,"%s","The chosen model is a biotin model, ");
            }else if(!strcasecmp("nb",model)){
                ksprintf(kstr,"%s","The chosen model is non-biotin model, ");
            }else{
                fprintf(stderr,"Please provide a meaningful deamination model for further calculations.\n");
                return -1;
            }
        }else{
            fprintf(stderr,"Please provide a deamination model for further calculations.\n");
            return -1;
        }
        ksprintf(kstr,"the inferred parameters are given as follows:\n");
        ksprintf(kstr,"lambda: %f (%f,%f), delta: %f (%f,%f), delta_s: %f (%f,%f), nu: %f (%f,%f), nfunctioncalls: %d ngradientcalls: %d, llh: %f.\n", invec2[0], invec2[0]-1.96*stdvec2[0], invec2[0]+1.96*stdvec2[0] , invec2[1], invec2[1]-1.96*stdvec2[1], invec2[1]+1.96*stdvec2[1], invec2[2], invec2[2]-1.96*stdvec2[2], invec2[2]+1.96*stdvec2[2], invec2[3], invec2[3]-1.96*stdvec2[3], invec2[3]+1.96*stdvec2[3], ncalls,ncalls_grad, withgrad2);
        //my_bgzf_write(fp,kstr->s,kstr->l);
        assert(bgzf_write(fp,kstr->s,kstr->l)==kstr->l);
        kstr->l = 0;
        bgzf_close(fp);
    }
    

    char* refName;
    refName = NULL;
    int mapq =-1;
    int mapped_only = 0;
    int se_only = 1;
    len_limit = 150;
    char *fname1 = mypars->ihts; //the file needs to be recalibrated
    char *olik1 = mypars->olik; //the output nucleotide likelihood file
    const char* chromname1 = mypars->ichr;
    const char* bedname1 = mypars->ibed;
    char * ofname1 = mypars->ohts;
    int isrecal = mypars->isrecal;
    
    deamRateCT = (double**) malloc((len_limit-30) * sizeof(double*));
    deamRateGA = (double**) malloc((len_limit-30) * sizeof(double*));
    for (int i = 0; i < (len_limit-30); i++){
        deamRateCT[i] =(double *) malloc(30 * sizeof(double));
        deamRateGA[i] =(double *) malloc(30 * sizeof(double));
    }
    if (isrecal>=0){
        if (!strcasecmp("b",model)){
            CaldeamRate_b(invec2[0],invec2[1],invec2[2],invec2[3],len_limit);
        }else if (!strcasecmp("nb",model)){
            CaldeamRate_nb(invec2[0],invec2[1],invec2[2],invec2[3],len_limit);
        }

    }

    double distparam[4] = {60,45,100,45};
    //    fprintf(stderr,"fname1: %s fname2: %s\n",fname1,fname2);
    
    sam_hdr_t *hdr;
    if (fname1!=NULL && isrecal == 1){
        //tsk start
        std::vector<bam1_t*> tsk_reads;
        hdr = read_all_reads(fname1,bedname1,refName,tsk_reads);
        //tsk stop
        fprintf(stderr,"BRANCH!\n");fflush(stderr);
        refName2 = refName;
        chromname2 = chromname1;
        bedname2 = bedname1;
        mapped_only2 = mapped_only;
        se_only2 = se_only;
        mapq2 = mapq;
        seq_ref2 = seq_ref;
        len_limit2 = len_limit;
        len_min2 = len_min;
        model2 = model;
        eps2 = Contam_eps;
        s2 = s;
        lambda2 = invec2[0];
        delta2 = invec2[1];
        delta_s2 = invec2[2];
        nv2 = invec2[3];
        ncalls = 0;
        ncalls_grad = 0;
        //double lbd2[4] = {30,1e-8,30,1e-8};
        double lbd2[4] = {30,15,30,1e-8};
        double ubd2[4] = {(double)len_limit-1,100,(double)len_limit-1,100};
        int nbd2[4] = {2,2,2,2};
        int mynsites = tsk_reads.size()/tsk_nthreads;
        fprintf(stderr,"mynSites: %d\n",mynsites);
        for(int ii=0;ii<tsk_nthreads;ii++){
            my_tsk_struct[ii].from = my_tsk_struct[ii].to = -1;
            my_tsk_struct[ii].reads = &tsk_reads;
            my_tsk_struct[ii].hdr = hdr;
            my_tsk_struct[ii].seq_ref = seq_ref;
            my_tsk_struct[ii].len_limit = len_limit;
            my_tsk_struct[ii].len_min = len_min;
            my_tsk_struct[ii].model = model;
            my_tsk_struct[ii].eps = Contam_eps;
            my_tsk_struct[ii].lambda = invec2[0];
            my_tsk_struct[ii].delta = invec2[1];
            my_tsk_struct[ii].delta_s = invec2[2];
            my_tsk_struct[ii].nv = invec2[3];
            my_tsk_struct[ii].x = new double [4];
            my_tsk_struct[ii].llh_result_grad = new double [4];
            my_tsk_struct[ii].llh_result_hess = new double* [4];
            for(int i=0;i<4;i++){
                my_tsk_struct[ii].x[i] = distparam[i];
                my_tsk_struct[ii].llh_result_hess[i] = new double [4];
            }
            my_tsk_struct[ii].threadid = ii;
            my_tsk_struct[ii].from = ii==0?0:my_tsk_struct[ii-1].to;
            my_tsk_struct[ii].to =   my_tsk_struct[ii].from+mynsites;
        }
        my_tsk_struct[tsk_nthreads-1].to = tsk_reads.size();
        double withgrad3;
        if(tsk_nthreads==1)
            withgrad3 = findmax_bfgs(4,distparam,&my_tsk_struct[0],tsk_All_loglike_recalibration,tsk_All_loglike_recalibration_grad,lbd2,ubd2,nbd2,-1);
        else{
            //withgrad3 = findmax_bfgs(4,distparam,NULL,like_master,tsk_All_loglike_recalibration_grad,lbd2,ubd2,nbd2,-1);
            withgrad3 = findmax_bfgs(4,distparam,NULL,like_master,like_grad_master,lbd2,ubd2,nbd2,-1);
        }
        double **covpar = new double* [4];
        for (int i=0;i<4;i++){
            covpar[i] = new double[4];
        }
        like_hess_master(distparam,covpar);
        double stdpar[4];
        for (int i=0;i<4;i++){
            stdpar[i] = sqrt(covpar[i][i]);
            //cout<<z4[i][i]<<" std "<<stdvec4[i]<<"\n";
        }
        fprintf(stderr,"mu_anc: %f (%f,%f), sigma_anc: %f (%f,%f), mu_mod: %f (%f,%f), sigma_mod: %f (%f,%f), nfunctioncalls: %d ngradientcalls: %d, llh: %f.\n", distparam[0], distparam[0]-1.96*stdpar[0], distparam[0]+1.96*stdpar[0], distparam[1], distparam[1]-1.96*stdpar[1], distparam[1]+1.96*stdpar[1], distparam[2], distparam[2]-1.96*stdpar[2], distparam[2]+1.96*stdpar[2], distparam[3], distparam[3]-1.96*stdpar[3], distparam[3]+1.96*stdpar[3], ncalls,ncalls_grad, withgrad3);
        sam_hdr_destroy(hdr);
    }
    //fprintf(stderr,"BRANCH 2000!\n");fflush(stderr);
    if ((fname1==NULL) && (olik1!=NULL)){
        fprintf(stdout,"Please provide bam file to calculate nucleotide likelihoods!\n");
        return 0;
    }
    
    BGZF *fp3; // Compressed nucliklihood file;
    kstring_t *kstr3 = new kstring_t;
    kstr3->s = NULL;
    kstr3->l = kstr3->m = 0;
 
    if (fname1!=NULL){
        bam_hdr_t* hdr = CalPostPMDProb(refName,fname1,chromname1,bedname1,ofname1,olik1,mapped_only,se_only,mapq,seq_ref,len_limit,len_min,model,Contam_eps,invec2[0],invec2[1],invec2[2],invec2[3],distparam[0],distparam[1],distparam[2], distparam[3], isrecal,s);
        sort(comp_nuc_llik.begin(),comp_nuc_llik.end(),cmp);
        merge(comp_nuc_llik);
        if (olik1!=NULL){
            const char* comp = "comp_"; // RASMUS  
            char * concat = (char*) malloc(strlen(comp) + strlen(olik1) + 1);
            if (concat){
                strcpy(concat, comp);
                strcat(concat, olik1);
            }
            fp3 = bgzf_open(concat,"wb");
            //           fprintf(stderr,"%s\t%s\t%s\t%s\t%s\t%s\n","ChrName","NucPos","A","C","G","T");
            ksprintf(kstr3,"%s\t%s\t%s\t%s\t%s\t%s\n","ChrName","NucPos","A","C","G","T");
            //           fprintf(stderr,"%s\t%s\t%s\t%s\t%s\t%s\n","ChrName","NucPos","A","C","G","T");
            //my_bgzf_write(fp3,kstr3->s,kstr3->l);
            assert(bgzf_write(fp3,kstr3->s,kstr3->l)==kstr3->l);
            //           fprintf(stderr,"%s\t%s\t%s\t%s\t%s\t%s\n","ChrName","NucPos","A","C","G","T");
            kstr3->l = 0;
            //           fprintf(stderr,"%s\t%s\t%s\t%s\t%s\t%s\n","ChrName","NucPos","A","C","G","T");
            for (size_t s = 0; s < comp_nuc_llik.size(); s++){
                //                fprintf(stderr,"%s\n",comp_nuc_llik[s].chr);
                //fprintf(stderr,"%lld\n",comp_nuc_llik[s].pos);
                //fprintf(stderr,"%f\t%f\t%f\t%f\n",comp_nuc_llik[s].nuclik[0],comp_nuc_llik[s].nuclik[1],comp_nuc_llik[s].nuclik[2],comp_nuc_llik[s].nuclik[3]);
                //                fprintf(stderr,"%s\t%lld\t%f\t%f\t%f\t%f\n",hdr1->target_name[comp_nuc_llik[s].chrid],comp_nuc_llik[s].pos,comp_nuc_llik[s].nuclik[0],comp_nuc_llik[s].nuclik[1],comp_nuc_llik[s].nuclik[2],comp_nuc_llik[s].nuclik[3]);
                ksprintf(kstr3,"%s\t%lld\t%f\t%f\t%f\t%f\n",hdr->target_name[comp_nuc_llik[s].chrid],comp_nuc_llik[s].pos,comp_nuc_llik[s].nuclik[0],comp_nuc_llik[s].nuclik[1],comp_nuc_llik[s].nuclik[2],comp_nuc_llik[s].nuclik[3]);
                //my_bgzf_write(fp3,kstr3->s,kstr3->l);
                assert(bgzf_write(fp3,kstr3->s,kstr3->l)==kstr3->l);
                kstr3->l = 0;
            }
            bgzf_close(fp3);
        }
        sam_hdr_destroy(hdr);
    }
    
    free(Frag_len);
    free(Frag_freq);
    free(Bin_Frag_len);
    free(Bin_Frag_freq);
    
    for (int i = 0; i < 4; i++){
        free(z2[i]);
    }
    free(z2);
    for (int i = 0; i < MAXLENGTH; i++){
        free(mm5p[i]);
        free(mm3p[i]);
    }
    for (int i = 0; i < (len_limit-30); i++){
        free(deamRateCT[i]);
        free(deamRateGA[i]);
    }
    free(deamRateCT);
    free(deamRateGA);
    free(seqError);
    free(freqCT);
    free(freqGA);
    free(scaleCT);
    free(scaleGA);
    free(mm5p);
    free(mm3p);
    argStruct_destroy(mypars);
    comp_nuc_llik.clear();

    // check deallocation of my_tsk_struct values
    delete my_tsk_struct;

//    if(hdr){
//        sam_hdr_destroy(hdr);
//    }
    return 0;
}


