#include "Likelihood.h"
#include <cstring>
#include <cstdlib>
#include <cstdio>
#include <cassert>
#include <htslib/hts.h>
#include <htslib/sam.h>
#include <htslib/kstring.h>
#include <htslib/faidx.h>
#include <zlib.h>
#include <cmath>
#include <iomanip>
#include <Eigen/Core>
#include <Eigen/Eigenvalues> //sort and merge


#include <ctime>
#include <getopt.h>
#include <iostream>
#include "profile.h"
#include "bfgs.h"
#include "htslib/bgzf.h"
#include "briggs_writer.h"
#include "read_all_reads.h"

#include "misc.h"
#include "Recalibration.h"
#include "ngsBriggs.h"
#include "PosteriorProb.h"

double PhredError[255];
double PhredErrorAThird[255];

double MAX0 = 1-1e-8;
double MIN0 = 1e-8;

extern tsk_struct *my_tsk_struct;

// Naive likelihood without nick frequencies
double loglike(const double *x, double * freqCT, double * freqGA, double * scaleCT, double * scaleGA){
    double lambda = x[0];
    double delta = x[1];
    double delta_s = x[2];
    double ll = 0;
    for(int i=0; i<MAXLENGTH; i++){
        ll += scaleCT[i]*(freqCT[i]*log(delta+pow(1-lambda,i+1)/2*(delta_s-delta))+(1-freqCT[i])*log(1-delta-pow(1-lambda,i+1)/2*(delta_s-delta)))+scaleGA[i]*(freqGA[i]*log(pow(1-lambda,i+1)/2*delta_s)+(1-freqGA[i])*log(1-pow(1-lambda,i+1)/2*delta_s));
    }
    return -ll;
}

double b_loglike(const double *x, const void *ptr){
    ncalls++;
    const wrapOne *wo =(const wrapOne *) ptr;
    return loglike(x, wo->freqCT, wo->freqGA, wo->scaleCT, wo->scaleGA);
}

void loglike_grad(const double *x,double *y, double * freqCT, double * freqGA, double * scaleCT, double * scaleGA){
    double lambda = x[0];
    double delta = x[1];
    double delta_s = x[2];
    y[0] = 0;
    y[1] = 0;
    y[2] = 0;
    for(int i=0; i<MAXLENGTH; i++){
        y[0] -= scaleCT[i]*(freqCT[i]*(-(i+1)*pow(1-lambda,i)/2*(delta_s-delta))/(delta+pow(1-lambda,i+1)/2*(delta_s-delta))+(1-freqCT[i])*((i+1)*pow(1-lambda,i)/2*(delta_s-delta))/(1-delta-pow(1-lambda,i+1)/2*(delta_s-delta)))+scaleGA[i]*(freqGA[i]*(-(i+1))/(1-lambda)+(1-freqGA[i])*((i+1)*pow(1-lambda,i)/2*delta_s)/(1-pow(1-lambda,i+1)/2*delta_s));
        y[1] -= scaleCT[i]*(freqCT[i]*(1-pow(1-lambda,i+1)/2)/(delta+pow(1-lambda,i+1)/2*(delta_s-delta))+(1-freqCT[i])*(-1+pow(1-lambda,i+1)/2)/(1-delta-pow(1-lambda,i+1)/2*(delta_s-delta)));
        y[2] -= scaleCT[i]*(freqCT[i]*(pow(1-lambda,i+1)/2)/(delta+pow(1-lambda,i+1)/2*(delta_s-delta))+(1-freqCT[i])*(-pow(1-lambda,i+1)/2)/(1-delta-pow(1-lambda,i+1)/2*(delta_s-delta)))+scaleGA[i]*(freqGA[i]*1/(delta_s)+(1-freqGA[i])*(-pow(1-lambda,i+1)/2)/(1-pow(1-lambda,i+1)/2*delta_s));
    }
}

void b_loglike_grad(const double *x,double *y,const void*ptr){
    ncalls_grad++;
    const wrapOne *wo =(const wrapOne *) ptr;
    loglike_grad(x,y, wo->freqCT, wo->freqGA, wo->scaleCT, wo->scaleGA);
}


//Consider contamination rate eps + sequencing error
double loglike_complex3_full_b(const double *x, double * freqCT, double * freqGA, double * scaleCT, double * scaleGA, double * seqError, int BinNum, double* LEN, double* freqLEN, double eps){
    double lambda = x[0];
    double delta = x[1];
    double delta_s = x[2];
    double nu = x[3];
    // Change the length from a fixed value to a distribution
    double ll = 0;
    for(int n=0; n<MAXLENGTH; n++){
        double freqCT1 = 0;
        double freqGA1 = 0;
        double freqCT2 = 0;
        double freqGA2 = 0;
        for(int i=0; i<BinNum; i++){
            double L = LEN[i];
            double f = freqLEN[i];
            double p1_l = pow(1+lambda,2)*n*nu/(1+(L-2)*nu);
            for(int l=1;l<=min(MAXORDER,n);l++){
                p1_l += lambda*(1+lambda)*pow(1-lambda,l)*(n-l)*nu/(1+(L-2-l)*nu);
                for (int r=1;r<=min((double)MAXORDER,(double)(L-n-1));r++){
                    if (l+r<=MAXORDER){
                        p1_l += pow(lambda,2)*pow(1-lambda,l+r)*(n-l)*nu/(1+(L-2-l-r)*nu);
                    }
                }
            }
            for(int r=1;r<=min((double)MAXORDER,(double)(L-n-1));r++){
                p1_l += lambda*(1+lambda)*pow(1-lambda,r)*n*nu/(1+(L-2-r)*nu);
            }
            p1_l = p1_l/4;
            double p2_l = 0;
            double p3_l = 1 - pow(1-lambda,n+1)/2-p1_l;
            double p4_l = pow(1-lambda,n+1)/2;
            
            double p1_r = pow(1+lambda,2)*(L-n-1)*nu/(1+(L-2)*nu);
            for(int r=1;r<=min(MAXORDER,n);r++){
                p1_r += lambda*(1+lambda)*pow(1-lambda,r)*(L-n-1)*nu/(1+(L-2-r)*nu);
                for (int l=1;l<=min((double)MAXORDER,(double)(L-n-1));l++){
                    if (l+r<=MAXORDER){
                        p1_r += pow(lambda,2)*pow(1-lambda,l+r)*(L-n-1-l)*nu/(1+(L-2-l-r)*nu);
                    }
                }
            }
            for(int l=1;l<=min((double)MAXORDER,(double)(L-n-1));l++){
                p1_r += lambda*(1+lambda)*pow(1-lambda,l)*(L-n-1-l)*nu/(1+(L-2-l)*nu);
            }
            p1_r = p1_r/4;
            double p2_r = pow(1-lambda,n+1)/2;
            double p3_r = 1 - pow(1-lambda,n+1)/2-p1_r;
            double p4_r = 0;
            freqCT1 += (p3_l*delta+p4_l*delta_s)*f;
            freqGA1 += (p1_r*delta+p2_r*delta_s)*f;
            freqCT2 += (p3_r*delta)*f;
            freqGA2 += (p1_l*delta)*f;
        }
        freqCT1 = (1-eps)*freqCT1;
        freqGA1 = (1-eps)*freqGA1;
        freqCT2 = (1-eps)*freqCT2;
        freqGA2 = (1-eps)*freqGA2;
        double freqCT3 = freqCT1*(1-seqError[n]) + (1-freqCT1)*seqError[n]/3;
        double freqCC3 = (1-freqCT1)*(1-seqError[n]) + freqCT1*seqError[n]/3;
        double freqCT4 = freqCT2*(1-seqError[2*MAXLENGTH-1-n]) + (1-freqCT2)*seqError[2*MAXLENGTH-1-n]/3;
        double freqCC4 = (1-freqCT2)*(1-seqError[2*MAXLENGTH-1-n]) + freqCT2*seqError[2*MAXLENGTH-1-n]/3;
        double freqGA3 = freqGA1*(1-seqError[2*MAXLENGTH-1-n]) + (1-freqGA1)*seqError[2*MAXLENGTH-1-n]/3;
        double freqGG3 = (1-freqGA1)*(1-seqError[2*MAXLENGTH-1-n]) + freqGA1*seqError[2*MAXLENGTH-1-n]/3;
        double freqGA4 = freqGA2*(1-seqError[n]) + (1-freqGA2)*seqError[n]/3;
        double freqGG4 = (1-freqGA2)*(1-seqError[n]) + freqGA2*seqError[n]/3;
        double freq[8], count[8];
        freq[0] = freqCT3; freq[1] = freqCC3; freq[2] = freqCT4; freq[3] = freqCC4;
        freq[4] = freqGA3; freq[5] = freqGG3; freq[6] = freqGA4; freq[7] = freqGG4;
        count[0] = scaleCT[n]*freqCT[n]; count[1] = scaleCT[n]*(1-freqCT[n]);
        count[2] = scaleCT[2*MAXLENGTH-1-n]*freqCT[2*MAXLENGTH-1-n]; count[3] = scaleCT[2*MAXLENGTH-1-n]*(1-freqCT[2*MAXLENGTH-1-n]);
        count[4] = scaleGA[n]*freqGA[n]; count[5] = scaleGA[n]*(1-freqGA[n]);
        count[6] = scaleGA[2*MAXLENGTH-1-n]*freqGA[2*MAXLENGTH-1-n]; count[7] = scaleGA[2*MAXLENGTH-1-n]*(1-freqGA[2*MAXLENGTH-1-n]);
        for (int j=0;j<8;j++){
            if (freq[j]>0){
                ll += count[j]*log(freq[j]);
            }
            // Consider freq[j] > 1
        }
        //        ll += scaleCT[n]*(freqCT[n]*log(freqCT3)+(1-freqCT[n])*log(freqCC3))+scaleCT[2*MAXLENGTH-1-n]*(freqCT[2*MAXLENGTH-1-n]*log(freqCT4)+(1-freqCT[2*MAXLENGTH-1-n])*log(freqCC4))+scaleGA[n]*(freqGA[n]*log(freqGA3)+(1-freqGA[n])*log(freqGG3))+scaleGA[2*MAXLENGTH-1-n]*(freqGA[2*MAXLENGTH-1-n]*log(freqGA4)+(1-freqGA[2*MAXLENGTH-1-n])*log(freqGG4));
    }
    //    cout << "Likelihood "<<ll<<"\n";
    return -ll;
}


// Considering Seq Error + Non-biotin model: symmetric pattern
double loglike_complex3_full_nb(const double *x, double * freqCT, double * freqGA, double * scaleCT, double * scaleGA, double * seqError, int BinNum, double* LEN, double* freqLEN, double eps){
    double lambda = x[0];
    double delta = x[1];
    double delta_s = x[2];
    double nu = x[3];
    // Change the length from a fixed value to a distribution
    double ll = 0;
    for(int n=0; n<MAXLENGTH; n++){
        double freqCT1 = 0;
        //double freqGA1 = 0;
        double freqCT2 = 0;
        // double freqGA2 = 0;
        for(int i=0; i<BinNum; i++){
            double L = LEN[i];
            double f = freqLEN[i];
            double p1_l = pow(1+lambda,2)*n*nu/(1+(L-2)*nu);
            for(int l=1;l<=min(MAXORDER,n);l++){
                p1_l += lambda*(1+lambda)*pow(1-lambda,l)*(n-l)*nu/(1+(L-2-l)*nu);
                for (int r=1;r<=min((double)MAXORDER,(double)(L-n-1));r++){
                    if (l+r<=MAXORDER){
                        p1_l += pow(lambda,2)*pow(1-lambda,l+r)*(n-l)*nu/(1+(L-2-l-r)*nu);
                    }
                }
            }
            for(int r=1;r<=min((double)MAXORDER,(double)(L-n-1));r++){
                p1_l += lambda*(1+lambda)*pow(1-lambda,r)*n*nu/(1+(L-2-r)*nu);
            }
            p1_l = p1_l/4;
            double p2_l = 0;
            double p3_l = 1 - pow(1-lambda,n+1)/2-p1_l;
            double p4_l = pow(1-lambda,n+1)/2;
            
            double p1_r = pow(1+lambda,2)*(L-n-1)*nu/(1+(L-2)*nu);
            for(int r=1;r<=min(MAXORDER,n);r++){
                p1_r += lambda*(1+lambda)*pow(1-lambda,r)*(L-n-1)*nu/(1+(L-2-r)*nu);
                for (int l=1;l<=min((double)MAXORDER,(double)(L-n-1));l++){
                    if (l+r<=MAXORDER){
                        p1_r += pow(lambda,2)*pow(1-lambda,l+r)*(L-n-1-l)*nu/(1+(L-2-l-r)*nu);
                    }
                }
            }
            for(int l=1;l<=min((double)MAXORDER,(double)(L-n-1));l++){
                p1_r += lambda*(1+lambda)*pow(1-lambda,l)*(L-n-1-l)*nu/(1+(L-2-l)*nu);
            }
            p1_r = p1_r/4;
            double p2_r = pow(1-lambda,n+1)/2;
            double p3_r = 1 - pow(1-lambda,n+1)/2-p1_r;
            double p4_r = 0;
            double dd1 = ((p1_r+p3_l)/2*delta+(p2_r+p4_l)/2*delta_s)*f; //Lateral increments 1 for C to T (5') and G to A (3')
            double dd2 = (p1_l+p3_r)/2*delta*f; //Lateral increments 2 for C to T (3') and G to A (5')
            //freqCT1 += (p3_l*delta+p4_l*delta_s)*f;
            //freqGA1 += (p1_r*delta+p2_r*delta_s)*f;
            // freqCT2 += (p3_r*delta)*f;
            // freqGA2 += (p1_l*delta)*f;
            freqCT1 += dd1;
            freqCT2 += dd2;
        }
        freqCT1 = (1-eps)*freqCT1;
        freqCT2 = (1-eps)*freqCT2;
        double freqGA1 = freqCT1;
        double freqGA2 = freqCT2;
        double freqCT3 = freqCT1*(1-seqError[n]) + (1-freqCT1)*seqError[n]/3;
        double freqCC3 = (1-freqCT1)*(1-seqError[n]) + freqCT1*seqError[n]/3;
        double freqCT4 = freqCT2*(1-seqError[2*MAXLENGTH-1-n]) + (1-freqCT2)*seqError[2*MAXLENGTH-1-n]/3;
        double freqCC4 = (1-freqCT2)*(1-seqError[2*MAXLENGTH-1-n]) + freqCT2*seqError[2*MAXLENGTH-1-n]/3;
        double freqGA3 = freqGA1*(1-seqError[2*MAXLENGTH-1-n]) + (1-freqGA1)*seqError[2*MAXLENGTH-1-n]/3;
        double freqGG3 = (1-freqGA1)*(1-seqError[2*MAXLENGTH-1-n]) + freqGA1*seqError[2*MAXLENGTH-1-n]/3;
        double freqGA4 = freqGA2*(1-seqError[n]) + (1-freqGA2)*seqError[n]/3;
        double freqGG4 = (1-freqGA2)*(1-seqError[n]) + freqGA2*seqError[n]/3;
        double freq[8], count[8];
        freq[0] = freqCT3; freq[1] = freqCC3; freq[2] = freqCT4; freq[3] = freqCC4;
        freq[4] = freqGA3; freq[5] = freqGG3; freq[6] = freqGA4; freq[7] = freqGG4;
        count[0] = scaleCT[n]*freqCT[n]; count[1] = scaleCT[n]*(1-freqCT[n]);
        count[2] = scaleCT[2*MAXLENGTH-1-n]*freqCT[2*MAXLENGTH-1-n]; count[3] = scaleCT[2*MAXLENGTH-1-n]*(1-freqCT[2*MAXLENGTH-1-n]);
        count[4] = scaleGA[n]*freqGA[n]; count[5] = scaleGA[n]*(1-freqGA[n]);
        count[6] = scaleGA[2*MAXLENGTH-1-n]*freqGA[2*MAXLENGTH-1-n]; count[7] = scaleGA[2*MAXLENGTH-1-n]*(1-freqGA[2*MAXLENGTH-1-n]);
        for (int j=0;j<8;j++){
            if (freq[j]>0){
                ll += count[j]*log(freq[j]);
            }
            // Consider freq[j] > 1
        }
        //        ll += scaleCT[n]*(freqCT[n]*log(freqCT3)+(1-freqCT[n])*log(freqCC3))+scaleCT[2*MAXLENGTH-1-n]*(freqCT[2*MAXLENGTH-1-n]*log(freqCT4)+(1-freqCT[2*MAXLENGTH-1-n])*log(freqCC4))+scaleGA[n]*(freqGA[n]*log(freqGA3)+(1-freqGA[n])*log(freqGG3))+scaleGA[2*MAXLENGTH-1-n]*(freqGA[2*MAXLENGTH-1-n]*log(freqGA4)+(1-freqGA[2*MAXLENGTH-1-n])*log(freqGG4));
    }
    return -ll;
}

double b_loglike_complex3_full(const double *x, const void *ptr){
    ncalls++;
    const wrapOne *wo =(const wrapOne *) ptr;
    return loglike_complex3_full_b(x, wo->freqCT, wo->freqGA, wo->scaleCT, wo->scaleGA, wo->seqError, wo->BinNum, wo->Bin_Frag_len, wo->Bin_Frag_freq, wo->Contam_eps);
}

double nb_loglike_complex3_full(const double *x, const void *ptr){
    ncalls++;
    const wrapOne *wo =(const wrapOne *) ptr;
    return loglike_complex3_full_nb(x, wo->freqCT, wo->freqGA, wo->scaleCT, wo->scaleGA, wo->seqError, wo->BinNum, wo->Bin_Frag_len, wo->Bin_Frag_freq, wo->Contam_eps);
}

void loglike_complex3_grad_full_b(const double *x,double *y, double * freqCT, double * freqGA, double * scaleCT, double * scaleGA, double * seqError, int BinNum, double* LEN, double* freqLEN, double eps){
    double lambda = x[0];
    double delta = x[1];
    double delta_s = x[2];
    double nu = x[3];
    y[0] = 0;
    y[1] = 0;
    y[2] = 0;
    y[3] = 0;
    for(int n=0; n<MAXLENGTH; n++){
        double freqCT1= 0;
        double freqCT1_lambda= 0;
        double freqCT1_nv= 0;
        double freqGA1 = 0;
        double freqGA1_lambda = 0;
        double freqGA1_nv = 0;
        double freqCT2= 0;
        double freqCT2_lambda = 0;
        double freqCT2_nv = 0;
        double freqGA2= 0;
        double freqGA2_lambda = 0;
        double freqGA2_nv = 0;
        double pf3_l = 0;
        double pf4_l = 0;
        double pf1_r = 0;
        double pf2_r = 0;
        double pf3_r = 0;
        double pf1_l = 0;
        for(int i=0; i<BinNum; i++){
            double L = LEN[i];
            double f = freqLEN[i];
            double p1_l = pow(1+lambda,2)*n*nu/(1+(L-2)*nu); //double p1_l = pow(1+lambda,2)*(n+1)*nu/(1+(L-2)*nu);
            double p1_l_lambda = 2*(1+lambda)*n*nu/(1+(L-2)*nu);
            double p1_l_nv = pow(1+lambda,2)*n/pow(1+(L-2)*nu,2);
            for(int l=1;l<=min(MAXORDER,n);l++){
                p1_l += lambda*(1+lambda)*pow(1-lambda,l)*(n-l)*nu/(1+(L-2-l)*nu); //p1_l += lambda*(1+lambda)*pow(1-lambda,l)*(n-l)*nu/(1+(L-2-l)*nu);
                p1_l_lambda += ((1+2*lambda)*pow(1-lambda,l)-l*lambda*(1+lambda)*pow(1-lambda,l-1))*(n-l)*nu/(1+(L-2-l)*nu);
                p1_l_nv += lambda*(1+lambda)*pow(1-lambda,l)*(n-l)/pow(1+(L-2-l)*nu,2);
                for (int r=1;r<=min((double)MAXORDER,(double)(L-n-1));r++){
                    if (l+r<=MAXORDER){
                        p1_l += pow(lambda,2)*pow(1-lambda,l+r)*(n-l)*nu/(1+(L-2-l-r)*nu); //p1_l += pow(lambda,2)*pow(1-lambda,l+r)*(n-l)*nu/(1+(L-2-l-r)*nu);
                        p1_l_lambda += (2*lambda*pow(1-lambda,l+r)-(l+r)*pow(lambda,2)*pow(1-lambda,l+r-1))*(n-l)*nu/(1+(L-2-l-r)*nu);
                        p1_l_nv += pow(lambda,2)*pow(1-lambda,l+r)*(n-l)/pow(1+(L-2-l-r)*nu,2);
                    }
                }
            }
            for(int r=1;r<=min((double)MAXORDER,(double)(L-n-1));r++){
                p1_l += lambda*(1+lambda)*pow(1-lambda,r)*n*nu/(1+(L-2-r)*nu); // p1_l += lambda*(1+lambda)*pow(1-lambda,r)*n*nu/(1+(L-2-r)*nu);
                p1_l_lambda += ((1+2*lambda)*pow(1-lambda,r)-r*lambda*(1+lambda)*pow(1-lambda,r-1))*n*nu/(1+(L-2-r)*nu);
                p1_l_nv += lambda*(1+lambda)*pow(1-lambda,r)*n/pow(1+(L-2-r)*nu,2);
            }
            p1_l = p1_l/4;
            p1_l_lambda = p1_l_lambda/4;
            p1_l_nv = p1_l_nv/4;
            double p2_l = 0;
            double p2_l_lambda = 0;
            double p2_l_nv = 0;
            double p3_l = 1 - pow(1-lambda,n+1)/2-p1_l;
            double p3_l_lambda = (n+1)*pow(1-lambda,n)/2-p1_l_lambda;
            double p3_l_nv = -p1_l_nv;
            double p4_l = pow(1-lambda,n+1)/2;
            double p4_l_lambda = -(n+1)*pow(1-lambda,n)/2;
            double p4_l_nv = 0;
            
            double p1_r = pow(1+lambda,2)*(L-n-1)*nu/(1+(L-2)*nu); // p1_r = pow(1+lambda,2)*(L-n)*nu/(1+(L-2)*nu);
            double p1_r_lambda = 2*(1+lambda)*(L-n-1)*nu/(1+(L-2)*nu);
            double p1_r_nv = pow(1+lambda,2)*(L-n-1)/pow(1+(L-2)*nu,2);
            for(int r=1;r<=min(MAXORDER,n);r++){
                p1_r += lambda*(1+lambda)*pow(1-lambda,r)*(L-n-1)*nu/(1+(L-2-r)*nu); // p1_r += lambda*(1+lambda)*pow(1-lambda,r)*(L-n-1-r)*nu/(1+(L-2-r)*nu);
                p1_r_lambda += ((1+2*lambda)*pow(1-lambda,r)-r*lambda*(1+lambda)*pow(1-lambda,r-1))*(L-n-1)*nu/(1+(L-2-r)*nu);
                p1_r_nv += lambda*(1+lambda)*pow(1-lambda,r)*(L-n-1)/pow(1+(L-2-r)*nu,2);
                for (int l=1;l<=min((double)MAXORDER,(double)(L-n-1));l++){
                    if (l+r<=MAXORDER){
                        p1_r += pow(lambda,2)*pow(1-lambda,l+r)*(L-n-1-l)*nu/(1+(L-2-l-r)*nu); // p1_r += pow(lambda,2)*pow(1-lambda,l+r)*(L-n-1-r)*nu/(1+(L-2-l-r)*nu);
                        p1_r_lambda += (2*lambda*pow(1-lambda,l+r)-(l+r)*pow(lambda,2)*pow(1-lambda,l+r-1))*(L-n-1-l)*nu/(1+(L-2-l-r)*nu);
                        p1_r_nv += pow(lambda,2)*pow(1-lambda,l+r)*(L-n-1-l)/pow(1+(L-2-l-r)*nu,2);
                    }
                }
            }
            for(int l=1;l<=min((double)MAXORDER,(double)(L-n-1));l++){
                p1_r += lambda*(1+lambda)*pow(1-lambda,l)*(L-n-1-l)*nu/(1+(L-2-l)*nu); // p1_r += lambda*(1+lambda)*pow(1-lambda,l)*(L-n-1)*nu/(1+(L-2-l)*nu);
                p1_r_lambda += ((1+2*lambda)*pow(1-lambda,l)-l*lambda*(1+lambda)*pow(1-lambda,l-1))*(L-n-1-l)*nu/(1+(L-2-l)*nu);
                p1_r_nv += lambda*(1+lambda)*pow(1-lambda,l)*(L-n-1-l)/pow(1+(L-2-l)*nu,2);
            }
            p1_r = p1_r/4;
            p1_r_lambda = p1_r_lambda/4;
            p1_r_nv = p1_r_nv/4;
            double p2_r = pow(1-lambda,n+1)/2;
            double p2_r_lambda = -(n+1)*pow(1-lambda,n)/2;
            double p2_r_nv = 0;
            double p3_r = 1 - pow(1-lambda,n+1)/2-p1_r;
            double p3_r_lambda = (n+1)*pow(1-lambda,n)/2-p1_r_lambda;
            double p3_r_nv = -p1_r_nv;
            double p4_r = 0;
            double p4_r_lambda = 0;
            double p4_r_nv = 0; // cout<<p1_l<<" "<<p2_l<<" "<<p3_l<<" "<<p4_l<<"\n";   cout<<p1_r<<" "<<p2_r<<" "<<p3_r<<" "<<p4_r<<"\n";
            pf3_l += p3_l*f;
            pf4_l += p4_l*f;
            pf1_r += p1_r*f;
            pf2_r += p2_r*f;
            pf3_r += p3_r*f;
            pf1_l += p1_l*f;
            freqCT1 += (p3_l*delta+p4_l*delta_s)*f;
            freqCT1_lambda += (p3_l_lambda*delta+p4_l_lambda*delta_s)*f;
            freqCT1_nv += (p3_l_nv*delta+p4_l_nv*delta_s)*f;
            freqGA1 += (p1_r*delta+p2_r*delta_s)*f;
            freqGA1_lambda += (p1_r_lambda*delta+p2_r_lambda*delta_s)*f;
            freqGA1_nv += (p1_r_nv*delta+p2_r_nv*delta_s)*f;
            freqCT2 += (p3_r*delta)*f;
            freqCT2_lambda += (p3_r_lambda*delta)*f;
            freqCT2_nv += (p3_r_nv*delta)*f;
            freqGA2 += (p1_l*delta)*f;
            freqGA2_lambda += (p1_l_lambda*delta)*f;
            freqGA2_nv += (p1_l_nv*delta)*f;
        }
        freqCT1 = (1-eps)*freqCT1;
        freqCT1_lambda = (1-eps)*freqCT1_lambda;
        freqCT1_nv = (1-eps)*freqCT1_nv;
        freqCT2 = (1-eps)*freqCT2;
        freqCT2_lambda = (1-eps)*freqCT2_lambda;
        freqCT2_nv = (1-eps)*freqCT2_nv;
        freqGA1 = (1-eps)*freqGA1;
        freqGA1_lambda = (1-eps)*freqGA1_lambda;
        freqGA1_nv = (1-eps)*freqGA1_nv;
        freqGA2 = (1-eps)*freqGA2;
        freqGA2_lambda = (1-eps)*freqGA2_lambda;
        freqGA2_nv = (1-eps)*freqGA2_nv;
        pf3_l = (1-eps)*pf3_l;
        pf4_l = (1-eps)*pf4_l;
        pf1_r = (1-eps)*pf1_r;
        pf2_r = (1-eps)*pf2_r;
        pf3_r = (1-eps)*pf3_r;
        pf1_l = (1-eps)*pf1_l;
        
        double freqCT3 = freqCT1*(1-seqError[n]) + (1-freqCT1)*seqError[n]/3;
        double freqCC3 = (1-freqCT1)*(1-seqError[n]) + freqCT1*seqError[n]/3;
        double freqCT4 = freqCT2*(1-seqError[2*MAXLENGTH-1-n]) + (1-freqCT2)*seqError[2*MAXLENGTH-1-n]/3;
        double freqCC4 = (1-freqCT2)*(1-seqError[2*MAXLENGTH-1-n]) + freqCT2*seqError[2*MAXLENGTH-1-n]/3;
        double freqGA3 = freqGA1*(1-seqError[2*MAXLENGTH-1-n]) + (1-freqGA1)*seqError[2*MAXLENGTH-1-n]/3;
        double freqGG3 = (1-freqGA1)*(1-seqError[2*MAXLENGTH-1-n]) + freqGA1*seqError[2*MAXLENGTH-1-n]/3;
        double freqGA4 = freqGA2*(1-seqError[n]) + (1-freqGA2)*seqError[n]/3;
        double freqGG4 = (1-freqGA2)*(1-seqError[n]) + freqGA2*seqError[n]/3;
        double freqCT3_lambda = freqCT1_lambda*(1-4.0/3.0*seqError[n]); //double freqCC3_lambda =-freqCT1_lambda*(1-4.0/3.0*seqError[n]);
        double freqCT4_lambda = freqCT2_lambda*(1-4.0/3.0*seqError[2*MAXLENGTH-1-n]); //double freqCC4_lambda =-freqCT2_lambda*(1-4.0/3.0*seqError[2*MAXLENGTH-1-n]);
        double freqGA3_lambda = freqGA1_lambda*(1-4.0/3.0*seqError[2*MAXLENGTH-1-n]); //double freqGG3_lambda =-freqGA1_lambda*(1-4.0/3.0*seqError[2*MAXLENGTH-1-n]);
        double freqGA4_lambda = freqGA2_lambda*(1-4.0/3.0*seqError[n]); //double freqGG4_lambda =-freqGA2_lambda*(1-4.0/3.0*seqError[n]);
        double freqCT3_nv = freqCT1_nv*(1-4.0/3.0*seqError[n]); //double freqCC3_nv =-freqCT1_nv*(1-4.0/3.0*seqError[n]);
        double freqCT4_nv = freqCT2_nv*(1-4.0/3.0*seqError[2*MAXLENGTH-1-n]); //double freqCC4_nv =-freqCT2_nv*(1-4.0/3.0*seqError[2*MAXLENGTH-1-n]);
        double freqGA3_nv = freqGA1_nv*(1-4.0/3.0*seqError[2*MAXLENGTH-1-n]); //double freqGG3_nv =-freqGA1_nv*(1-4.0/3.0*seqError[2*MAXLENGTH-1-n]);
        double freqGA4_nv = freqGA2_nv*(1-4.0/3.0*seqError[n]); //double freqGG4_nv =-freqGA2_nv*(1-4.0/3.0*seqError[n]);
        
        
        double pf3_l1 = pf3_l*(1-4.0/3.0*seqError[n]); // d freqCT3/d delta
        double pf3_r1 = pf3_r*(1-4.0/3.0*seqError[2*MAXLENGTH-1-n]); // d freqCT4/d delta
        double pf1_r1 = pf1_r*(1-4.0/3.0*seqError[2*MAXLENGTH-1-n]); // d freqGA3/d delta
        double pf1_l1 = pf1_l*(1-4.0/3.0*seqError[n]); // d freqGA4/d delta
        double pf4_l1 = pf4_l*(1-4.0/3.0*seqError[n]); // d freqCT3/d delta_s
        double pf2_r1 = pf2_r*(1-4.0/3.0*seqError[2*MAXLENGTH-1-n]); // d freqGA3/ d delta_s
        
        double freq[8], freq_lambda[8], freq_delta[8], freq_deltas[8], freq_nv[8], count[8];
        freq[0] = freqCT3; freq[1] = freqCC3; freq[2] = freqCT4; freq[3] = freqCC4;
        freq[4] = freqGA3; freq[5] = freqGG3; freq[6] = freqGA4; freq[7] = freqGG4;
        freq_lambda[0] = freqCT3_lambda; freq_lambda[1] = -freqCT3_lambda; freq_lambda[2] = freqCT4_lambda; freq_lambda[3] = -freqCT4_lambda;
        freq_lambda[4] = freqGA3_lambda; freq_lambda[5] = -freqGA3_lambda; freq_lambda[6] = freqGA4_lambda; freq_lambda[7] = -freqGA4_lambda;
        freq_delta[0] = pf3_l1; freq_delta[1] = -pf3_l1; freq_delta[2] = pf3_r1; freq_delta[3] = -pf3_r1;
        freq_delta[4] = pf1_r1; freq_delta[5] = -pf1_r1; freq_delta[6] = pf1_l1; freq_delta[7] = -pf1_l1;
        freq_deltas[0] = pf4_l1; freq_deltas[1] = -pf4_l1; freq_deltas[2] = 0; freq_deltas[3] = 0;
        freq_deltas[4] = pf2_r1; freq_deltas[5] = -pf2_r1; freq_deltas[6] = 0; freq_deltas[7] = 0;
        freq_nv[0] = freqCT3_nv; freq_nv[1] = -freqCT3_nv; freq_nv[2] = freqCT4_nv; freq_nv[3] = -freqCT4_nv;
        freq_nv[4] = freqGA3_nv; freq_nv[5] = -freqGA3_nv; freq_nv[6] = freqGA4_nv; freq_nv[7] = -freqGA4_nv;
        count[0] = scaleCT[n]*freqCT[n]; count[1] = scaleCT[n]*(1-freqCT[n]);
        count[2] = scaleCT[2*MAXLENGTH-1-n]*freqCT[2*MAXLENGTH-1-n]; count[3] = scaleCT[2*MAXLENGTH-1-n]*(1-freqCT[2*MAXLENGTH-1-n]);
        count[4] = scaleGA[n]*freqGA[n]; count[5] = scaleGA[n]*(1-freqGA[n]);
        count[6] = scaleGA[2*MAXLENGTH-1-n]*freqGA[2*MAXLENGTH-1-n]; count[7] = scaleGA[2*MAXLENGTH-1-n]*(1-freqGA[2*MAXLENGTH-1-n]);
        for (int j=0;j<8;j++){
            if (freq[j]>0){
                y[0] -= count[j]/freq[j]*freq_lambda[j]; // Derivative of lambda
                y[1] -= count[j]/freq[j]*freq_delta[j]; // Derivative of delta
                y[2] -= count[j]/freq[j]*freq_deltas[j]; // Derivative of delta_s
                y[3] -= count[j]/freq[j]*freq_nv[j]; // Derivative of nv
            }
        }
    }
}

void loglike_complex3_grad_full_nb(const double *x,double *y, double * freqCT, double * freqGA, double * scaleCT, double * scaleGA, double * seqError, int BinNum, double* LEN, double* freqLEN, double eps){
    double lambda = x[0];
    double delta = x[1];
    double delta_s = x[2];
    double nu = x[3];
    y[0] = 0;
    y[1] = 0;
    y[2] = 0;
    y[3] = 0;
    for(int n=0; n<MAXLENGTH; n++){
        double freqCT1= 0;
        double freqCT1_lambda= 0;
        double freqCT1_nv= 0;
        double freqCT2= 0;
        double freqCT2_lambda = 0;
        double freqCT2_nv = 0;
        // double freqGA1 = 0;
        // double freqGA1_lambda = 0;
        // double freqGA1_nv = 0;
        // double freqGA2= 0;
        // double freqGA2_lambda = 0;
        // double freqGA2_nv = 0;
        double pf3_l = 0;
        double pf4_l = 0;
        double pf1_r = 0;
        double pf2_r = 0;
        double pf3_r = 0;
        double pf1_l = 0;
        for(int i=0; i<BinNum; i++){
            double L = LEN[i];
            double f = freqLEN[i];
            double p1_l = pow(1+lambda,2)*n*nu/(1+(L-2)*nu); //double p1_l = pow(1+lambda,2)*(n+1)*nu/(1+(L-2)*nu);
            double p1_l_lambda = 2*(1+lambda)*n*nu/(1+(L-2)*nu);
            double p1_l_nv = pow(1+lambda,2)*n/pow(1+(L-2)*nu,2);
            for(int l=1;l<=min(MAXORDER,n);l++){
                p1_l += lambda*(1+lambda)*pow(1-lambda,l)*(n-l)*nu/(1+(L-2-l)*nu); //p1_l += lambda*(1+lambda)*pow(1-lambda,l)*(n-l)*nu/(1+(L-2-l)*nu);
                p1_l_lambda += ((1+2*lambda)*pow(1-lambda,l)-l*lambda*(1+lambda)*pow(1-lambda,l-1))*(n-l)*nu/(1+(L-2-l)*nu);
                p1_l_nv += lambda*(1+lambda)*pow(1-lambda,l)*(n-l)/pow(1+(L-2-l)*nu,2);
                for (int r=1;r<=min((double)MAXORDER,(double)(L-n-1));r++){
                    if (l+r<=MAXORDER){
                        p1_l += pow(lambda,2)*pow(1-lambda,l+r)*(n-l)*nu/(1+(L-2-l-r)*nu); //p1_l += pow(lambda,2)*pow(1-lambda,l+r)*(n-l)*nu/(1+(L-2-l-r)*nu);
                        p1_l_lambda += (2*lambda*pow(1-lambda,l+r)-(l+r)*pow(lambda,2)*pow(1-lambda,l+r-1))*(n-l)*nu/(1+(L-2-l-r)*nu);
                        p1_l_nv += pow(lambda,2)*pow(1-lambda,l+r)*(n-l)/pow(1+(L-2-l-r)*nu,2);
                    }
                }
            }
            for(int r=1;r<=min((double)MAXORDER,(double)(L-n-1));r++){
                p1_l += lambda*(1+lambda)*pow(1-lambda,r)*n*nu/(1+(L-2-r)*nu); // p1_l += lambda*(1+lambda)*pow(1-lambda,r)*n*nu/(1+(L-2-r)*nu);
                p1_l_lambda += ((1+2*lambda)*pow(1-lambda,r)-r*lambda*(1+lambda)*pow(1-lambda,r-1))*n*nu/(1+(L-2-r)*nu);
                p1_l_nv += lambda*(1+lambda)*pow(1-lambda,r)*n/pow(1+(L-2-r)*nu,2);
            }
            p1_l = p1_l/4;
            p1_l_lambda = p1_l_lambda/4;
            p1_l_nv = p1_l_nv/4;
            double p2_l = 0;
            double p2_l_lambda = 0;
            double p2_l_nv = 0;
            double p3_l = 1 - pow(1-lambda,n+1)/2-p1_l;
            double p3_l_lambda = (n+1)*pow(1-lambda,n)/2-p1_l_lambda;
            double p3_l_nv = -p1_l_nv;
            double p4_l = pow(1-lambda,n+1)/2;
            double p4_l_lambda = -(n+1)*pow(1-lambda,n)/2;
            double p4_l_nv = 0;
            
            double p1_r = pow(1+lambda,2)*(L-n-1)*nu/(1+(L-2)*nu); // p1_r = pow(1+lambda,2)*(L-n)*nu/(1+(L-2)*nu);
            double p1_r_lambda = 2*(1+lambda)*(L-n-1)*nu/(1+(L-2)*nu);
            double p1_r_nv = pow(1+lambda,2)*(L-n-1)/pow(1+(L-2)*nu,2);
            for(int r=1;r<=min(MAXORDER,n);r++){
                p1_r += lambda*(1+lambda)*pow(1-lambda,r)*(L-n-1)*nu/(1+(L-2-r)*nu); // p1_r += lambda*(1+lambda)*pow(1-lambda,r)*(L-n-1-r)*nu/(1+(L-2-r)*nu);
                p1_r_lambda += ((1+2*lambda)*pow(1-lambda,r)-r*lambda*(1+lambda)*pow(1-lambda,r-1))*(L-n-1)*nu/(1+(L-2-r)*nu);
                p1_r_nv += lambda*(1+lambda)*pow(1-lambda,r)*(L-n-1)/pow(1+(L-2-r)*nu,2);
                for (int l=1;l<=min((double)MAXORDER,(double)(L-n-1));l++){
                    if (l+r<=MAXORDER){
                        p1_r += pow(lambda,2)*pow(1-lambda,l+r)*(L-n-1-l)*nu/(1+(L-2-l-r)*nu); // p1_r += pow(lambda,2)*pow(1-lambda,l+r)*(L-n-1-r)*nu/(1+(L-2-l-r)*nu);
                        p1_r_lambda += (2*lambda*pow(1-lambda,l+r)-(l+r)*pow(lambda,2)*pow(1-lambda,l+r-1))*(L-n-1-l)*nu/(1+(L-2-l-r)*nu);
                        p1_r_nv += pow(lambda,2)*pow(1-lambda,l+r)*(L-n-1-l)/pow(1+(L-2-l-r)*nu,2);
                    }
                }
            }
            for(int l=1;l<=min((double)MAXORDER,(double)(L-n-1));l++){
                p1_r += lambda*(1+lambda)*pow(1-lambda,l)*(L-n-1-l)*nu/(1+(L-2-l)*nu); // p1_r += lambda*(1+lambda)*pow(1-lambda,l)*(L-n-1)*nu/(1+(L-2-l)*nu);
                p1_r_lambda += ((1+2*lambda)*pow(1-lambda,l)-l*lambda*(1+lambda)*pow(1-lambda,l-1))*(L-n-1-l)*nu/(1+(L-2-l)*nu);
                p1_r_nv += lambda*(1+lambda)*pow(1-lambda,l)*(L-n-1-l)/pow(1+(L-2-l)*nu,2);
            }
            p1_r = p1_r/4;
            p1_r_lambda = p1_r_lambda/4;
            p1_r_nv = p1_r_nv/4;
            double p2_r = pow(1-lambda,n+1)/2;
            double p2_r_lambda = -(n+1)*pow(1-lambda,n)/2;
            double p2_r_nv = 0;
            double p3_r = 1 - pow(1-lambda,n+1)/2-p1_r;
            double p3_r_lambda = (n+1)*pow(1-lambda,n)/2-p1_r_lambda;
            double p3_r_nv = -p1_r_nv;
            double p4_r = 0;
            double p4_r_lambda = 0;
            double p4_r_nv = 0; // cout<<p1_l<<" "<<p2_l<<" "<<p3_l<<" "<<p4_l<<"\n";   cout<<p1_r<<" "<<p2_r<<" "<<p3_r<<" "<<p4_r<<"\n";
            pf3_l += p3_l*f;
            pf4_l += p4_l*f;
            pf1_r += p1_r*f;
            pf2_r += p2_r*f;
            pf3_r += p3_r*f;
            pf1_l += p1_l*f;
            double dd1 = ((p1_r+p3_l)/2*delta+(p2_r+p4_l)/2*delta_s)*f;
            double dd1_lambda = ((p1_r_lambda+p3_l_lambda)/2*delta+(p2_r_lambda+p4_l_lambda)/2*delta_s)*f;
            double dd1_nv = ((p1_r_nv+p3_l_nv)/2*delta+(p2_r_nv+p4_l_nv)/2*delta_s)*f;
            double dd2 = (p1_l+p3_r)/2*delta*f;
            double dd2_lambda = (p1_l_lambda+p3_r_lambda)/2*delta*f;
            double dd2_nv = (p1_l_nv+p3_r_nv)/2*delta*f;
            freqCT1 += dd1;
            freqCT1_lambda += dd1_lambda;
            freqCT1_nv += dd1_nv;
            freqCT2 += dd2;
            freqCT2_lambda += dd2_lambda;
            freqCT2_nv += dd2_nv;
          
        }
        freqCT1 = (1-eps)*freqCT1;
        freqCT1_lambda = (1-eps)*freqCT1_lambda;
        freqCT1_nv = (1-eps)*freqCT1_nv;
        freqCT2 = (1-eps)*freqCT2;
        freqCT2_lambda = (1-eps)*freqCT2_lambda;
        freqCT2_nv = (1-eps)*freqCT2_nv;
        pf3_l = (1-eps)*pf3_l;
        pf4_l = (1-eps)*pf4_l;
        pf1_r = (1-eps)*pf1_r;
        pf2_r = (1-eps)*pf2_r;
        pf3_r = (1-eps)*pf3_r;
        pf1_l = (1-eps)*pf1_l;
        
        double freqGA1 = freqCT1;
        double freqGA1_lambda = freqCT1_lambda;
        double freqGA1_nv = freqCT1_nv;
        double freqGA2 = freqCT2;
        double freqGA2_lambda = freqCT2_lambda;
        double freqGA2_nv = freqCT2_nv;
        // double freqGA2 = 0;
        // double freqGA2_lambda = 0;
        // double freqGA2_nv = 0;
        double freqCT3 = freqCT1*(1-seqError[n]) + (1-freqCT1)*seqError[n]/3;
        double freqCC3 = (1-freqCT1)*(1-seqError[n]) + freqCT1*seqError[n]/3;
        double freqCT4 = freqCT2*(1-seqError[2*MAXLENGTH-1-n]) + (1-freqCT2)*seqError[2*MAXLENGTH-1-n]/3;
        double freqCC4 = (1-freqCT2)*(1-seqError[2*MAXLENGTH-1-n]) + freqCT2*seqError[2*MAXLENGTH-1-n]/3;
        double freqGA3 = freqGA1*(1-seqError[2*MAXLENGTH-1-n]) + (1-freqGA1)*seqError[2*MAXLENGTH-1-n]/3;
        double freqGG3 = (1-freqGA1)*(1-seqError[2*MAXLENGTH-1-n]) + freqGA1*seqError[2*MAXLENGTH-1-n]/3;
        double freqGA4 = freqGA2*(1-seqError[n]) + (1-freqGA2)*seqError[n]/3;
        double freqGG4 = (1-freqGA2)*(1-seqError[n]) + freqGA2*seqError[n]/3;
        double freqCT3_lambda = freqCT1_lambda*(1-4.0/3.0*seqError[n]); //double freqCC3_lambda =-freqCT1_lambda*(1-4.0/3.0*seqError[n]);
        double freqCT4_lambda = freqCT2_lambda*(1-4.0/3.0*seqError[2*MAXLENGTH-1-n]); //double freqCC4_lambda =-freqCT2_lambda*(1-4.0/3.0*seqError[2*MAXLENGTH-1-n]);
        double freqGA3_lambda = freqGA1_lambda*(1-4.0/3.0*seqError[2*MAXLENGTH-1-n]); //double freqGG3_lambda =-freqGA1_lambda*(1-4.0/3.0*seqError[2*MAXLENGTH-1-n]);
        double freqGA4_lambda = freqGA2_lambda*(1-4.0/3.0*seqError[n]); //double freqGG4_lambda =-freqGA2_lambda*(1-4.0/3.0*seqError[n]);
        double freqCT3_nv = freqCT1_nv*(1-4.0/3.0*seqError[n]); //double freqCC3_nv =-freqCT1_nv*(1-4.0/3.0*seqError[n]);
        double freqCT4_nv = freqCT2_nv*(1-4.0/3.0*seqError[2*MAXLENGTH-1-n]); //double freqCC4_nv =-freqCT2_nv*(1-4.0/3.0*seqError[2*MAXLENGTH-1-n]);
        double freqGA3_nv = freqGA1_nv*(1-4.0/3.0*seqError[2*MAXLENGTH-1-n]); //double freqGG3_nv =-freqGA1_nv*(1-4.0/3.0*seqError[2*MAXLENGTH-1-n]);
        double freqGA4_nv = freqGA2_nv*(1-4.0/3.0*seqError[n]); //double freqGG4_nv =-freqGA2_nv*(1-4.0/3.0*seqError[n]);
        
        double pf3_l1 = (pf1_r+pf3_l)/2*(1-4.0/3.0*seqError[n]); // d freqCT3/d delta
        double pf3_r1 = (pf1_l+pf3_r)/2*(1-4.0/3.0*seqError[2*MAXLENGTH-1-n]); // d freqCT4/d delta
        double pf1_r1 = (pf1_r+pf3_l)/2*(1-4.0/3.0*seqError[2*MAXLENGTH-1-n]); // d freqGA3/d delta
        double pf1_l1 = (pf1_l+pf3_r)/2*(1-4.0/3.0*seqError[n]); // d freqGA4/d delta
        double pf4_l1 = (pf2_r+pf4_l)/2*(1-4.0/3.0*seqError[n]); // d freqCT3/d delta_s
        double pf2_r1 = (pf2_r+pf4_l)/2*(1-4.0/3.0*seqError[2*MAXLENGTH-1-n]); // d freqGA3/ d delta_s
        
        double freq[8], freq_lambda[8], freq_delta[8], freq_deltas[8], freq_nv[8], count[8];
        freq[0] = freqCT3; freq[1] = freqCC3; freq[2] = freqCT4; freq[3] = freqCC4;
        freq[4] = freqGA3; freq[5] = freqGG3; freq[6] = freqGA4; freq[7] = freqGG4;
        freq_lambda[0] = freqCT3_lambda; freq_lambda[1] = -freqCT3_lambda; freq_lambda[2] = freqCT4_lambda; freq_lambda[3] = -freqCT4_lambda;
        freq_lambda[4] = freqGA3_lambda; freq_lambda[5] = -freqGA3_lambda; freq_lambda[6] = freqGA4_lambda; freq_lambda[7] = -freqGA4_lambda;
        freq_delta[0] = pf3_l1; freq_delta[1] = -pf3_l1; freq_delta[2] = pf3_r1; freq_delta[3] = -pf3_r1;
        freq_delta[4] = pf1_r1; freq_delta[5] = -pf1_r1; freq_delta[6] = pf1_l1; freq_delta[7] = -pf1_l1;
        freq_deltas[0] = pf4_l1; freq_deltas[1] = -pf4_l1; freq_deltas[2] = 0; freq_deltas[3] = 0;
        freq_deltas[4] = pf2_r1; freq_deltas[5] = -pf2_r1; freq_deltas[6] = 0; freq_deltas[7] = 0;
        freq_nv[0] = freqCT3_nv; freq_nv[1] = -freqCT3_nv; freq_nv[2] = freqCT4_nv; freq_nv[3] = -freqCT4_nv;
        freq_nv[4] = freqGA3_nv; freq_nv[5] = -freqGA3_nv; freq_nv[6] = freqGA4_nv; freq_nv[7] = -freqGA4_nv;
        count[0] = scaleCT[n]*freqCT[n]; count[1] = scaleCT[n]*(1-freqCT[n]);
        count[2] = scaleCT[2*MAXLENGTH-1-n]*freqCT[2*MAXLENGTH-1-n]; count[3] = scaleCT[2*MAXLENGTH-1-n]*(1-freqCT[2*MAXLENGTH-1-n]);
        count[4] = scaleGA[n]*freqGA[n]; count[5] = scaleGA[n]*(1-freqGA[n]);
        count[6] = scaleGA[2*MAXLENGTH-1-n]*freqGA[2*MAXLENGTH-1-n]; count[7] = scaleGA[2*MAXLENGTH-1-n]*(1-freqGA[2*MAXLENGTH-1-n]);
        for (int j=0;j<8;j++){
            if (freq[j]>0){
                y[0] -= count[j]/freq[j]*freq_lambda[j]; // Derivative of lambda
                y[1] -= count[j]/freq[j]*freq_delta[j]; // Derivative of delta
                y[2] -= count[j]/freq[j]*freq_deltas[j]; // Derivative of delta_s
                y[3] -= count[j]/freq[j]*freq_nv[j]; // Derivative of nv
            }
        }
       
    }
}


void b_loglike_complex3_grad_full(const double *x,double *y,const void *ptr){
    ncalls_grad++;
    const wrapOne *wo =(const wrapOne *) ptr;
    loglike_complex3_grad_full_b(x, y, wo->freqCT, wo->freqGA, wo->scaleCT, wo->scaleGA, wo->seqError, wo->BinNum, wo->Bin_Frag_len, wo->Bin_Frag_freq, wo->Contam_eps);
}

void nb_loglike_complex3_grad_full(const double *x,double *y,const void *ptr){
    ncalls_grad++;
    const wrapOne *wo =(const wrapOne *) ptr;
    loglike_complex3_grad_full_nb(x, y, wo->freqCT, wo->freqGA, wo->scaleCT, wo->scaleGA, wo->seqError, wo->BinNum, wo->Bin_Frag_len, wo->Bin_Frag_freq, wo->Contam_eps);
}

void loglike_complex3_hessian_full_b(const double *x, double ** z, double * freqCT, double * freqGA, double * scaleCT, double * scaleGA, double * seqError, int BinNum, double * LEN, double * freqLEN, double eps){
    Eigen::Matrix4d A, B;
    A = Eigen::Matrix4d::Zero();
    double lambda = x[0];
    double delta = x[1];
    double delta_s = x[2];
    double nu = x[3];
    
    for(int n=0; n<MAXLENGTH; n++){
        double freqCT1 = 0;
        double freqCT1_lambda = 0;
        double freqCT1_nv = 0;
        double freqCT1_lambda2 = 0;
        double freqCT1_lambda_nv= 0;
        double freqCT1_nv2= 0;
        double freqGA1 = 0;
        double freqGA1_lambda = 0;
        double freqGA1_nv = 0;
        double freqGA1_lambda2 = 0;
        double freqGA1_lambda_nv = 0;
        double freqGA1_nv2 = 0;
        double freqCT2 = 0;
        double freqCT2_lambda = 0;
        double freqCT2_nv = 0;
        double freqCT2_lambda2 = 0;
        double freqCT2_lambda_nv= 0;
        double freqCT2_nv2= 0;
        double freqGA2 = 0;
        double freqGA2_lambda = 0;
        double freqGA2_nv = 0;
        double freqGA2_lambda2 = 0;
        double freqGA2_lambda_nv = 0;
        double freqGA2_nv2 = 0;
        double pf3_l = 0;
        double pf3_l_lambda = 0;
        double pf3_l_nv = 0;
        double pf4_l = 0;
        double pf4_l_lambda = 0;
        double pf4_l_nv = 0;
        double pf1_r = 0;
        double pf1_r_lambda = 0;
        double pf1_r_nv = 0;
        double pf2_r = 0;
        double pf2_r_lambda = 0;
        double pf2_r_nv = 0;
        double pf3_r = 0;
        double pf3_r_lambda = 0;
        double pf3_r_nv = 0;
        double pf1_l = 0;
        double pf1_l_lambda = 0;
        double pf1_l_nv = 0;
        for(int i=0; i<BinNum; i++){
            double L = LEN[i];
            double f = freqLEN[i];
            double p1_l = pow(1+lambda,2)*n*nu/(1+(L-2)*nu); //double p1_l = pow(1+lambda,2)*(n+1)*nu/(1+(L-2)*nu);
            double p1_l_lambda = 2*(1+lambda)*n*nu/(1+(L-2)*nu);
            double p1_l_nv = pow(1+lambda,2)*n/pow(1+(L-2)*nu,2);
            double p1_l_lambda2 = 2*n*nu/(1+(L-2)*nu);
            double p1_l_lambda_nv = 2*(1+lambda)*n/pow(1+(L-2)*nu,2);
            double p1_l_nv2 = -2*pow(1+lambda,2)*n*(L-2)/pow(1+(L-2)*nu,3);
            for(int l=1;l<=min(MAXORDER,n);l++){
                p1_l += lambda*(1+lambda)*pow(1-lambda,l)*(n-l)*nu/(1+(L-2-l)*nu); //p1_l += lambda*(1+lambda)*pow(1-lambda,l)*(n-l)*nu/(1+(L-2-l)*nu);
                p1_l_lambda += ((1+2*lambda)*pow(1-lambda,l)-l*lambda*(1+lambda)*pow(1-lambda,l-1))*(n-l)*nu/(1+(L-2-l)*nu);
                p1_l_nv += lambda*(1+lambda)*pow(1-lambda,l)*(n-l)/pow(1+(L-2-l)*nu,2);
                p1_l_lambda2 += (2*pow(1-lambda,l)-2*l*(1+2*lambda)*pow(1-lambda,l-1)+l*(l-1)*lambda*(1+lambda)*pow(1-lambda,l-2))*(n-l)*nu/(1+(L-2-l)*nu);
                p1_l_lambda_nv += ((1+2*lambda)*pow(1-lambda,l)-l*lambda*(1+lambda)*pow(1-lambda,l-1))*(n-l)/pow(1+(L-2-l)*nu,2);
                p1_l_nv2 += -2*lambda*(1+lambda)*pow(1-lambda,l)*(n-l)*(L-2-l)/pow(1+(L-2-l)*nu,3);
                for (int r=1;r<=min((double)MAXORDER,(double)(L-n-1));r++){
                    if (l+r<=MAXORDER){
                        p1_l += pow(lambda,2)*pow(1-lambda,l+r)*(n-l)*nu/(1+(L-2-l-r)*nu); //p1_l += pow(lambda,2)*pow(1-lambda,l+r)*(n-l)*nu/(1+(L-2-l-r)*nu);
                        p1_l_lambda += (2*lambda*pow(1-lambda,l+r)-(l+r)*pow(lambda,2)*pow(1-lambda,l+r-1))*(n-l)*nu/(1+(L-2-l-r)*nu);
                        p1_l_nv += pow(lambda,2)*pow(1-lambda,l+r)*(n-l)/pow(1+(L-2-l-r)*nu,2);
                        p1_l_lambda2 += (2*pow(1-lambda,l+r)-4*(l+r)*lambda*pow(1-lambda,l+r-1)+(l+r)*(l+r-1)*pow(lambda,2)*pow(1-lambda,l+r-2))*(n-l)*nu/(1+(L-2-l-r)*nu);
                        p1_l_lambda_nv += (2*lambda*pow(1-lambda,l+r)-(l+r)*pow(lambda,2)*pow(1-lambda,l+r-1))*(n-l)/pow(1+(L-2-l-r)*nu,2);
                        p1_l_nv2 += -2*pow(lambda,2)*pow(1-lambda,l+r)*(n-l)*(L-2-l-r)/pow(1+(L-2-l-r)*nu,3);
                    }
                }
            }
            for(int r=1;r<=min((double)MAXORDER,(double)(L-n-1));r++){
                p1_l += lambda*(1+lambda)*pow(1-lambda,r)*n*nu/(1+(L-2-r)*nu); // p1_l += lambda*(1+lambda)*pow(1-lambda,r)*n*nu/(1+(L-2-r)*nu);
                p1_l_lambda += ((1+2*lambda)*pow(1-lambda,r)-r*lambda*(1+lambda)*pow(1-lambda,r-1))*n*nu/(1+(L-2-r)*nu);
                p1_l_nv += lambda*(1+lambda)*pow(1-lambda,r)*n/pow(1+(L-2-r)*nu,2);
                p1_l_lambda2 += (2*pow(1-lambda,r)-2*r*(1+2*lambda)*pow(1-lambda,r-1)+r*(r-1)*lambda*(1+lambda)*pow(1-lambda,r-2))*n*nu/(1+(L-2-r)*nu);
                p1_l_lambda_nv += ((1+2*lambda)*pow(1-lambda,r)-r*lambda*(1+lambda)*pow(1-lambda,r-1))*n/pow(1+(L-2-r)*nu,2);
                p1_l_nv2 += -2*lambda*(1+lambda)*pow(1-lambda,r)*n*(L-2-r)/pow(1+(L-2-r)*nu,3);
            }
            p1_l = p1_l/4;
            p1_l_lambda = p1_l_lambda/4;
            p1_l_nv = p1_l_nv/4;
            p1_l_lambda2 = p1_l_lambda2/4;
            p1_l_lambda_nv = p1_l_lambda_nv/4;
            p1_l_nv2 = p1_l_nv2/4;
            double p2_l = 0;
            double p2_l_lambda = 0;
            double p2_l_nv = 0;
            double p2_l_lambda2 = 0;
            double p2_l_lambda_nv = 0;
            double p2_l_nv2 = 0;
            double p3_l = 1 - pow(1-lambda,n+1)/2-p1_l;
            double p3_l_lambda = (n+1)*pow(1-lambda,n)/2-p1_l_lambda;
            double p3_l_nv = -p1_l_nv;
            double p3_l_lambda2 = -(n+1)*n*pow(1-lambda,n-1)/2-p1_l_lambda2;
            double p3_l_lambda_nv = -p1_l_lambda_nv;
            double p3_l_nv2 = -p1_l_nv2;
            double p4_l = pow(1-lambda,n+1)/2;
            double p4_l_lambda = -(n+1)*pow(1-lambda,n)/2;
            double p4_l_nv = 0;
            double p4_l_lambda2 = (n+1)*n*pow(1-lambda,n-1)/2;
            double p4_l_lambda_nv = 0;
            double p4_l_nv2 = 0;
            
            double p1_r = pow(1+lambda,2)*(L-n-1)*nu/(1+(L-2)*nu); // p1_r = pow(1+lambda,2)*(L-n)*nu/(1+(L-2)*nu);
            double p1_r_lambda = 2*(1+lambda)*(L-n-1)*nu/(1+(L-2)*nu);
            double p1_r_nv = pow(1+lambda,2)*(L-n-1)/pow(1+(L-2)*nu,2);
            double p1_r_lambda2 = 2*(L-n-1)*nu/(1+(L-2)*nu);
            double p1_r_lambda_nv = 2*(1+lambda)*(L-n-1)/pow(1+(L-2)*nu,2);
            double p1_r_nv2 = -2*pow(1+lambda,2)*(L-n-1)*(L-2)/pow(1+(L-2)*nu,3);
            for(int r=1;r<=min(MAXORDER,n);r++){
                p1_r += lambda*(1+lambda)*pow(1-lambda,r)*(L-n-1)*nu/(1+(L-2-r)*nu); // p1_r += lambda*(1+lambda)*pow(1-lambda,r)*(L-n-1-r)*nu/(1+(L-2-r)*nu);
                p1_r_lambda += ((1+2*lambda)*pow(1-lambda,r)-r*lambda*(1+lambda)*pow(1-lambda,r-1))*(L-n-1)*nu/(1+(L-2-r)*nu);
                p1_r_nv += lambda*(1+lambda)*pow(1-lambda,r)*(L-n-1)/pow(1+(L-2-r)*nu,2);
                p1_r_lambda2 += (2*pow(1-lambda,r)-2*r*(1+2*lambda)*pow(1-lambda,r-1)+r*(r-1)*lambda*(1+lambda)*pow(1-lambda,r-2))*(L-n-1)*nu/(1+(L-2-r)*nu);
                p1_r_lambda_nv += ((1+2*lambda)*pow(1-lambda,r)-r*lambda*(1+lambda)*pow(1-lambda,r-1))*(L-n-1)/pow(1+(L-2-r)*nu,2);
                p1_r_nv2 += -2*lambda*(1+lambda)*pow(1-lambda,r)*(L-n-1)*(L-2-r)/pow(1+(L-2-r)*nu,3);
                for (int l=1;l<=min((double)MAXORDER,(double)(L-n-1));l++){
                    if (l+r<=MAXORDER){
                        p1_r += pow(lambda,2)*pow(1-lambda,l+r)*(L-n-1-l)*nu/(1+(L-2-l-r)*nu); // p1_r += pow(lambda,2)*pow(1-lambda,l+r)*(L-n-1-r)*nu/(1+(L-2-l-r)*nu);
                        p1_r_lambda += (2*lambda*pow(1-lambda,l+r)-(l+r)*pow(lambda,2)*pow(1-lambda,l+r-1))*(L-n-1-l)*nu/(1+(L-2-l-r)*nu);
                        p1_r_nv += pow(lambda,2)*pow(1-lambda,l+r)*(L-n-1-l)/pow(1+(L-2-l-r)*nu,2);
                        p1_r_lambda2 += (2*pow(1-lambda,l+r)-4*(l+r)*lambda*pow(1-lambda,l+r-1)+(l+r)*(l+r-1)*pow(lambda,2)*pow(1-lambda,l+r-2))*(L-n-1-l)*nu/(1+(L-2-l-r)*nu);
                        p1_r_lambda_nv += (2*lambda*pow(1-lambda,l+r)-(l+r)*pow(lambda,2)*pow(1-lambda,l+r-1))*(L-n-1-l)/pow(1+(L-2-l-r)*nu,2);
                        p1_r_nv2 += -2*pow(lambda,2)*pow(1-lambda,l+r)*(L-n-1-l)*(L-2-l-r)/pow(1+(L-2-l-r)*nu,3);
                    }
                }
            }
            for(int l=1;l<=min((double)MAXORDER,(double)(L-n-1));l++){
                p1_r += lambda*(1+lambda)*pow(1-lambda,l)*(L-n-1-l)*nu/(1+(L-2-l)*nu); // p1_r += lambda*(1+lambda)*pow(1-lambda,l)*(L-n-1)*nu/(1+(L-2-l)*nu);
                p1_r_lambda += ((1+2*lambda)*pow(1-lambda,l)-l*lambda*(1+lambda)*pow(1-lambda,l-1))*(L-n-1-l)*nu/(1+(L-2-l)*nu);
                p1_r_nv += lambda*(1+lambda)*pow(1-lambda,l)*(L-n-1-l)/pow(1+(L-2-l)*nu,2);
                p1_r_lambda2 += (2*pow(1-lambda,l)-2*l*(1+2*lambda)*pow(1-lambda,l-1)+l*(l-1)*lambda*(1+lambda)*pow(1-lambda,l-2))*(L-n-1-l)*nu/(1+(L-2-l)*nu);
                p1_r_lambda_nv += ((1+2*lambda)*pow(1-lambda,l)-l*lambda*(1+lambda)*pow(1-lambda,l-1))*(L-n-1-l)/pow(1+(L-2-l)*nu,2);
                p1_r_nv2 += -2*lambda*(1+lambda)*pow(1-lambda,l)*(L-n-1-l)*(L-2-l)/pow(1+(L-2-l)*nu,3);
            }
            p1_r = p1_r/4;
            p1_r_lambda = p1_r_lambda/4;
            p1_r_nv = p1_r_nv/4;
            p1_r_lambda2 = p1_r_lambda2/4;
            p1_r_lambda_nv = p1_r_lambda_nv/4;
            p1_r_nv2 = p1_r_nv2/4;
            double p2_r = pow(1-lambda,n+1)/2;
            double p2_r_lambda = -(n+1)*pow(1-lambda,n)/2;
            double p2_r_nv = 0;
            double p2_r_lambda2 = (n+1)*n*pow(1-lambda,n-1)/2;
            double p2_r_lambda_nv = 0;
            double p2_r_nv2 = 0;
            double p3_r = 1 - pow(1-lambda,n+1)/2-p1_r;
            double p3_r_lambda = (n+1)*pow(1-lambda,n)/2-p1_r_lambda;
            double p3_r_nv = -p1_r_nv;
            double p3_r_lambda2 = -(n+1)*n*pow(1-lambda,n-1)/2-p1_r_lambda2;
            double p3_r_lambda_nv = -p1_r_lambda_nv;
            double p3_r_nv2 = -p1_r_nv2;
            double p4_r = 0;
            double p4_r_lambda = 0;
            double p4_r_nv = 0;
            double p4_r_lambda2 = 0;
            double p4_r_lambda_nv = 0;
            double p4_r_nv2 = 0;
            //        cout<<p1_l<<" "<<p2_l<<" "<<p3_l<<" "<<p4_l<<"\n";
            //        cout<<p1_r<<" "<<p2_r<<" "<<p3_r<<" "<<p4_r<<"\n";
            freqCT1 += (p3_l*delta+p4_l*delta_s)*f;
            freqCT1_lambda += (p3_l_lambda*delta+p4_l_lambda*delta_s)*f;
            freqCT1_nv += (p3_l_nv*delta+p4_l_nv*delta_s)*f;
            freqCT1_lambda2 += (p3_l_lambda2*delta+p4_l_lambda2*delta_s)*f;
            freqCT1_lambda_nv += (p3_l_lambda_nv*delta+p4_l_lambda_nv*delta_s)*f;
            freqCT1_nv2 += (p3_l_nv2*delta+p4_l_nv2*delta_s)*f;
            freqGA1 += (p1_r*delta+p2_r*delta_s)*f;
            freqGA1_lambda += (p1_r_lambda*delta+p2_r_lambda*delta_s)*f;
            freqGA1_nv += (p1_r_nv*delta+p2_r_nv*delta_s)*f;
            freqGA1_lambda2 += (p1_r_lambda2*delta+p2_r_lambda2*delta_s)*f;
            freqGA1_lambda_nv += (p1_r_lambda_nv*delta+p2_r_lambda_nv*delta_s)*f;
            freqGA1_nv2 += (p1_r_nv2*delta+p2_r_nv2*delta_s)*f;
            freqCT2 += (p3_r*delta)*f;
            freqCT2_lambda += (p3_r_lambda*delta)*f;
            freqCT2_nv += (p3_r_nv*delta)*f;
            freqCT2_lambda2 += (p3_r_lambda2*delta)*f;
            freqCT2_lambda_nv += (p3_r_lambda_nv*delta)*f;
            freqCT2_nv2 += (p3_r_nv2*delta)*f;
            freqGA2 += (p1_l*delta)*f;
            freqGA2_lambda += (p1_l_lambda*delta)*f;
            freqGA2_nv += (p1_l_nv*delta)*f;
            freqGA2_lambda2 += (p1_l_lambda2*delta)*f;
            freqGA2_lambda_nv += (p1_l_lambda_nv*delta)*f;
            freqGA2_nv2 += (p1_l_nv2*delta)*f;
            
            pf3_l += p3_l*f;
            pf3_l_lambda += p3_l_lambda*f;
            pf3_l_nv += p3_l_nv*f;
            pf4_l += p4_l*f;
            pf4_l_lambda += p4_l_lambda*f;
            pf4_l_nv += pf4_l_nv*f;
            pf1_r += p1_r*f;
            pf1_r_lambda += p1_r_lambda*f;
            pf1_r_nv += p1_r_nv*f;
            pf2_r += p2_r*f;
            pf2_r_lambda += p2_r_lambda*f;
            pf2_r_nv += p2_r_nv*f;
            pf3_r += p3_r*f;
            pf3_r_lambda += p3_r_lambda*f;
            pf3_r_nv += p3_r_nv*f;
            pf1_l += p1_l*f;
            pf1_l_lambda += p1_l_lambda*f;
            pf1_l_nv += p1_l_nv*f;
        }
        freqCT1 = (1-eps)*freqCT1;
        freqCT1_lambda = (1-eps)*freqCT1_lambda;
        freqCT1_nv = (1-eps)*freqCT1_nv;
        freqCT2 = (1-eps)*freqCT2;
        freqCT2_lambda = (1-eps)*freqCT2_lambda;
        freqCT2_nv = (1-eps)*freqCT2_nv;
        freqGA1 = (1-eps)*freqGA1;
        freqGA1_lambda = (1-eps)*freqGA1_lambda;
        freqGA1_nv = (1-eps)*freqGA1_nv;
        freqGA2 = (1-eps)*freqGA2;
        freqGA2_lambda = (1-eps)*freqGA2_lambda;
        freqGA2_nv = (1-eps)*freqGA2_nv;
        pf3_l = (1-eps)*pf3_l;
        pf4_l = (1-eps)*pf4_l;
        pf1_r = (1-eps)*pf1_r;
        pf2_r = (1-eps)*pf2_r;
        pf3_r = (1-eps)*pf3_r;
        pf1_l = (1-eps)*pf1_l;
        freqCT1_lambda2 = (1-eps)*freqCT1_lambda2;
        freqCT1_lambda_nv = (1-eps)*freqCT1_lambda_nv;
        freqCT1_nv2 = (1-eps)*freqCT1_nv2;
        freqCT2_lambda2 = (1-eps)*freqCT2_lambda2;
        freqCT2_lambda_nv = (1-eps)*freqCT2_lambda_nv;
        freqCT2_nv2 = (1-eps)*freqCT2_nv2;
        freqGA1_lambda2 = (1-eps)*freqGA1_lambda2;
        freqGA1_lambda_nv = (1-eps)*freqGA1_lambda_nv;
        freqGA1_nv2 = (1-eps)*freqGA1_nv2;
        freqGA2_lambda2 = (1-eps)*freqGA2_lambda2;
        freqGA2_lambda_nv = (1-eps)*freqGA2_lambda_nv;
        freqGA2_nv2 = (1-eps)*freqGA2_nv2;
        pf3_l_lambda = (1-eps)*pf3_l_lambda;
        pf3_l_nv = (1-eps)*pf3_l_nv;
        pf4_l_lambda = (1-eps)*pf4_l_lambda;
        pf4_l_nv = (1-eps)*pf4_l_nv;
        pf1_r_lambda = (1-eps)*pf1_r_lambda;
        pf1_r_nv = (1-eps)*pf1_r_nv;
        pf2_r_lambda = (1-eps)*pf2_r_lambda;
        pf2_r_nv = (1-eps)*pf2_r_nv;
        pf3_r_lambda = (1-eps)*pf3_r_lambda;
        pf3_r_nv = (1-eps)*pf3_r_nv;
        pf1_l_lambda = (1-eps)*pf1_l_lambda;
        pf1_l_nv = (1-eps)*pf1_l_nv;
        
        double freqCT3 = freqCT1*(1-seqError[n]) + (1-freqCT1)*seqError[n]/3;
        double freqCC3 = (1-freqCT1)*(1-seqError[n]) + freqCT1*seqError[n]/3;
        double freqCT4 = freqCT2*(1-seqError[2*MAXLENGTH-1-n]) + (1-freqCT2)*seqError[2*MAXLENGTH-1-n]/3;
        double freqCC4 = (1-freqCT2)*(1-seqError[2*MAXLENGTH-1-n]) + freqCT2*seqError[2*MAXLENGTH-1-n]/3;
        double freqGA3 = freqGA1*(1-seqError[2*MAXLENGTH-1-n]) + (1-freqGA1)*seqError[2*MAXLENGTH-1-n]/3;
        double freqGG3 = (1-freqGA1)*(1-seqError[2*MAXLENGTH-1-n]) + freqGA1*seqError[2*MAXLENGTH-1-n]/3;
        double freqGA4 = freqGA2*(1-seqError[n]) + (1-freqGA2)*seqError[n]/3;
        double freqGG4 = (1-freqGA2)*(1-seqError[n]) + freqGA2*seqError[n]/3;
        double freqCT3_lambda = freqCT1_lambda*(1-4.0/3.0*seqError[n]); //double freqCC3_lambda =-freqCT1_lambda*(1-4.0/3.0*seqError[n]);
        double freqCT4_lambda = freqCT2_lambda*(1-4.0/3.0*seqError[2*MAXLENGTH-1-n]); //double freqCC4_lambda =-freqCT2_lambda*(1-4.0/3.0*seqError[2*MAXLENGTH-1-n]);
        double freqGA3_lambda = freqGA1_lambda*(1-4.0/3.0*seqError[2*MAXLENGTH-1-n]); //double freqGG3_lambda =-freqGA1_lambda*(1-4.0/3.0*seqError[2*MAXLENGTH-1-n]);
        double freqGA4_lambda = freqGA2_lambda*(1-4.0/3.0*seqError[n]); //double freqGG4_lambda =-freqGA2_lambda*(1-4.0/3.0*seqError[n]);
        double freqCT3_nv = freqCT1_nv*(1-4.0/3.0*seqError[n]); //double freqCC3_nv =-freqCT1_nv*(1-4.0/3.0*seqError[n]);
        double freqCT4_nv = freqCT2_nv*(1-4.0/3.0*seqError[2*MAXLENGTH-1-n]); //double freqCC4_nv =-freqCT2_nv*(1-4.0/3.0*seqError[2*MAXLENGTH-1-n]);
        double freqGA3_nv = freqGA1_nv*(1-4.0/3.0*seqError[2*MAXLENGTH-1-n]); //double freqGG3_nv =-freqGA1_nv*(1-4.0/3.0*seqError[2*MAXLENGTH-1-n]);
        double freqGA4_nv = freqGA2_nv*(1-4.0/3.0*seqError[n]); //double freqGG4_nv =-freqGA2_nv*(1-4.0/3.0*seqError[n]);
        
        double pf3_l1 = pf3_l*(1-4.0/3.0*seqError[n]);
        double pf3_r1 = pf3_r*(1-4.0/3.0*seqError[2*MAXLENGTH-1-n]);
        double pf1_r1 = pf1_r*(1-4.0/3.0*seqError[2*MAXLENGTH-1-n]);
        double pf1_l1 = pf1_l*(1-4.0/3.0*seqError[n]);
        double pf4_l1 = pf4_l*(1-4.0/3.0*seqError[n]);
        double pf2_r1 = pf2_r*(1-4.0/3.0*seqError[2*MAXLENGTH-1-n]);
        
        double freqCT3_lambda2 = freqCT1_lambda2*(1-4.0/3.0*seqError[n]);
        double freqCT4_lambda2 = freqCT2_lambda2*(1-4.0/3.0*seqError[2*MAXLENGTH-1-n]);
        double freqGA3_lambda2 = freqGA1_lambda2*(1-4.0/3.0*seqError[2*MAXLENGTH-1-n]); //double freqGG3_lambda =-freqGA1_lambda*(1-4.0/3.0*seqError[2*MAXLENGTH-1-n]);
        double freqGA4_lambda2 = freqGA2_lambda2*(1-4.0/3.0*seqError[n]); //double freqGG4_lambda =-freqGA2_lambda*(1-4.0/3.0*seqError[n]);
        double pf3_l1_lambda = pf3_l_lambda*(1-4.0/3.0*seqError[n]);
        double pf3_r1_lambda = pf3_r_lambda*(1-4.0/3.0*seqError[2*MAXLENGTH-1-n]);
        double pf1_r1_lambda = pf1_r_lambda*(1-4.0/3.0*seqError[2*MAXLENGTH-1-n]);
        double pf1_l1_lambda = pf1_l_lambda*(1-4.0/3.0*seqError[n]);
        double pf4_l1_lambda = pf4_l_lambda*(1-4.0/3.0*seqError[n]);
        double pf2_r1_lambda = pf2_r_lambda*(1-4.0/3.0*seqError[2*MAXLENGTH-1-n]);
        double pf3_l1_nv = pf3_l_nv*(1-4.0/3.0*seqError[n]);
        double pf3_r1_nv = pf3_r_nv*(1-4.0/3.0*seqError[2*MAXLENGTH-1-n]);
        double pf1_r1_nv = pf1_r_nv*(1-4.0/3.0*seqError[2*MAXLENGTH-1-n]);
        double pf1_l1_nv = pf1_l_nv*(1-4.0/3.0*seqError[n]);
        double pf4_l1_nv = pf4_l_nv*(1-4.0/3.0*seqError[n]);
        double pf2_r1_nv = pf2_r_nv*(1-4.0/3.0*seqError[2*MAXLENGTH-1-n]);
        double freqCT3_lambda_nv = freqCT1_lambda_nv*(1-4.0/3.0*seqError[n]);
        double freqCT4_lambda_nv = freqCT2_lambda_nv*(1-4.0/3.0*seqError[2*MAXLENGTH-1-n]);
        double freqGA3_lambda_nv = freqGA1_lambda_nv*(1-4.0/3.0*seqError[2*MAXLENGTH-1-n]); //double freqGG3_lambda =-freqGA1_lambda*(1-4.0/3.0*seqError[2*MAXLENGTH-1-n]);
        double freqGA4_lambda_nv = freqGA2_lambda_nv*(1-4.0/3.0*seqError[n]); //double freqGG4_lambda =-freqGA2_lambda*(1-4.0/3.0*seqError[n]);
        double freqCT3_nv2 = freqCT1_nv2*(1-4.0/3.0*seqError[n]);
        double freqCT4_nv2 = freqCT2_nv2*(1-4.0/3.0*seqError[2*MAXLENGTH-1-n]);
        double freqGA3_nv2 = freqGA1_nv2*(1-4.0/3.0*seqError[2*MAXLENGTH-1-n]); //double freqGG3_lambda =-freqGA1_lambda*(1-4.0/3.0*seqError[2*MAXLENGTH-1-n]);
        double freqGA4_nv2 = freqGA2_nv2*(1-4.0/3.0*seqError[n]); //double freqGG4_lambda =-freqGA2_lambda*(1-4.0/3.0*seqError[n]);
        double freq[8], freq_lambda[8], freq_delta[8], freq_deltas[8], freq_nv[8], count[8], freq_lambda2[8], freq_nv2[8], freq_lambda_delta[8], freq_lambda_deltas[8], freq_delta_nv[8], freq_deltas_nv[8], freq_lambda_nv[8];
        freq[0] = freqCT3; freq[1] = freqCC3; freq[2] = freqCT4; freq[3] = freqCC4;
        freq[4] = freqGA3; freq[5] = freqGG3; freq[6] = freqGA4; freq[7] = freqGG4;
        freq_lambda[0] = freqCT3_lambda; freq_lambda[1] = -freqCT3_lambda; freq_lambda[2] = freqCT4_lambda; freq_lambda[3] = -freqCT4_lambda;
        freq_lambda[4] = freqGA3_lambda; freq_lambda[5] = -freqGA3_lambda; freq_lambda[6] = freqGA4_lambda; freq_lambda[7] = -freqGA4_lambda;
        freq_delta[0] = pf3_l1; freq_delta[1] = -pf3_l1; freq_delta[2] = pf3_r1; freq_delta[3] = -pf3_r1;
        freq_delta[4] = pf1_r1; freq_delta[5] = -pf1_r1; freq_delta[6] = pf1_l1; freq_delta[7] = -pf1_l1;
        freq_deltas[0] = pf4_l1; freq_deltas[1] = -pf4_l1; freq_deltas[2] = 0; freq_deltas[3] = 0;
        freq_deltas[4] = pf2_r1; freq_deltas[5] = -pf2_r1; freq_deltas[6] = 0; freq_deltas[7] = 0;
        freq_nv[0] = freqCT3_nv; freq_nv[1] = -freqCT3_nv; freq_nv[2] = freqCT4_nv; freq_nv[3] = -freqCT4_nv;
        freq_nv[4] = freqGA3_nv; freq_nv[5] = -freqGA3_nv; freq_nv[6] = freqGA4_nv; freq_nv[7] = -freqGA4_nv;
        count[0] = scaleCT[n]*freqCT[n]; count[1] = scaleCT[n]*(1-freqCT[n]);
        count[2] = scaleCT[2*MAXLENGTH-1-n]*freqCT[2*MAXLENGTH-1-n]; count[3] = scaleCT[2*MAXLENGTH-1-n]*(1-freqCT[2*MAXLENGTH-1-n]);
        count[4] = scaleGA[n]*freqGA[n]; count[5] = scaleGA[n]*(1-freqGA[n]);
        count[6] = scaleGA[2*MAXLENGTH-1-n]*freqGA[2*MAXLENGTH-1-n]; count[7] = scaleGA[2*MAXLENGTH-1-n]*(1-freqGA[2*MAXLENGTH-1-n]);
        freq_lambda2[0] = freqCT3_lambda2; freq_lambda2[1] = -freqCT3_lambda2; freq_lambda2[2] = freqCT4_lambda2; freq_lambda2[3] = -freqCT4_lambda2;
        freq_lambda2[4] = freqGA3_lambda2; freq_lambda2[5] = -freqGA3_lambda2; freq_lambda2[6] = freqGA4_lambda2; freq_lambda2[7] = -freqGA4_lambda2;
        freq_nv2[0] = freqCT3_nv2; freq_nv2[1] = -freqCT3_nv2; freq_nv2[2] = freqCT4_nv2; freq_nv2[3] = -freqCT4_nv2;
        freq_nv2[4] = freqGA3_nv2; freq_nv2[5] = -freqGA3_nv2; freq_nv2[6] = freqGA4_nv2; freq_nv2[7] = -freqGA4_nv2;
        freq_lambda_delta[0] = pf3_l1_lambda; freq_lambda_delta[1] = -pf3_l1_lambda; freq_lambda_delta[2] = pf3_r1_lambda; freq_lambda_delta[3] = -pf3_r1_lambda;
        freq_lambda_delta[4] = pf1_r1_lambda; freq_lambda_delta[5] = -pf1_r1_lambda; freq_lambda_delta[6] = pf1_l1_lambda; freq_lambda_delta[7] = -pf1_l1_lambda;
        freq_lambda_deltas[0] = pf4_l1_lambda; freq_lambda_deltas[1] = -pf4_l1_lambda; freq_lambda_deltas[2] = 0; freq_lambda_deltas[3] = 0;
        freq_lambda_deltas[4] = pf2_r1_lambda; freq_lambda_deltas[5] = -pf2_r1_lambda; freq_lambda_deltas[6] = 0; freq_lambda_deltas[7] = 0;
        freq_delta_nv[0] = pf3_l1_nv; freq_delta_nv[1] = -pf3_l1_nv; freq_delta_nv[2] = pf3_r1_nv; freq_delta_nv[3] = -pf3_r1_nv;
        freq_delta_nv[4] = pf1_r1_nv; freq_delta_nv[5] = -pf1_r1_nv; freq_delta_nv[6] = pf1_l1_nv; freq_delta_nv[7] = -pf1_l1_nv;
        freq_deltas_nv[0] = pf4_l1_nv; freq_deltas_nv[1] = -pf4_l1_nv; freq_deltas_nv[2] = 0; freq_deltas_nv[3] = 0;
        freq_deltas_nv[4] = pf2_r1_nv; freq_deltas_nv[5] = -pf2_r1_nv; freq_deltas_nv[6] = 0; freq_deltas_nv[7] = 0;
        freq_lambda_nv[0] = freqCT3_lambda_nv; freq_lambda_nv[1] = -freqCT3_lambda_nv; freq_lambda_nv[2] = freqCT4_lambda_nv; freq_lambda_nv[3] = -freqCT4_lambda_nv;
        freq_lambda_nv[4] = freqGA3_lambda_nv; freq_lambda_nv[5] = -freqGA3_lambda_nv; freq_lambda_nv[6] = freqGA4_lambda_nv; freq_lambda_nv[7] = -freqGA4_lambda_nv;
        for (int j=0;j<8;j++){
            if (freq[j]>0){
                A(0,0) += -count[j]/pow(freq[j],2)*pow(freq_lambda[j],2)+count[j]/freq[j]*freq_lambda2[j];
                A(0,1) += -count[j]/pow(freq[j],2)*freq_lambda[j]*freq_delta[j]+count[j]/freq[j]*freq_lambda_delta[j];
                A(0,2) += -count[j]/pow(freq[j],2)*freq_lambda[j]*freq_deltas[j]+count[j]/freq[j]*freq_lambda_deltas[j];
                A(0,3) += -count[j]/pow(freq[j],2)*freq_lambda[j]*freq_nv[j]+count[j]/freq[j]*freq_lambda_nv[j];
                A(1,1) += -count[j]/pow(freq[j],2)*pow(freq_delta[j],2);
                A(1,2) += -count[j]/pow(freq[j],2)*freq_delta[j]*freq_deltas[j];
                A(1,3) += -count[j]/pow(freq[j],2)*freq_delta[j]*freq_nv[j]+count[j]/freq[j]*freq_delta_nv[j];
                A(2,2) += -count[j]/pow(freq[j],2)*pow(freq_deltas[j],2);
                A(2,3) += -count[j]/pow(freq[j],2)*freq_deltas[j]*freq_nv[j]+count[j]/freq[j]*freq_deltas_nv[j];
                A(3,3) += -count[j]/pow(freq[j],2)*pow(freq_nv[j],2)+count[j]/freq[j]*freq_nv2[j];
            }
        }
        //        A(0,0) += scaleCT[n]*(-freqCT[n]/pow(freqCT3,2)*pow(freqCT3_lambda,2)+freqCT[n]/freqCT3*freqCT3_lambda2-(1-freqCT[n])/pow(freqCC3,2)*pow(freqCT3_lambda,2)-(1-freqCT[n])/freqCC3*freqCT3_lambda2)+scaleCT[2*MAXLENGTH-1-n]*(-freqCT[2*MAXLENGTH-1-n]/pow(freqCT4,2)*pow(freqCT4_lambda,2)+freqCT[2*MAXLENGTH-1-n]/freqCT4*freqCT4_lambda2-(1-freqCT[2*MAXLENGTH-1-n])/pow(freqCC4,2)*pow(freqCT4_lambda,2)-(1-freqCT[2*MAXLENGTH-1-n])/freqCC4*freqCT4_lambda2)+scaleGA[n]*(-freqGA[n]/pow(freqGA3,2)*pow(freqGA3_lambda,2)+freqGA[n]/freqGA3*freqGA3_lambda2-(1-freqGA[n])/pow(freqGG3,2)*pow(freqGA3_lambda,2)-(1-freqGA[n])/freqGG3*freqGA3_lambda2)+scaleGA[2*MAXLENGTH-1-n]*(-freqGA[2*MAXLENGTH-1-n]/pow(freqGA4,2)*pow(freqGA4_lambda,2)+freqGA[2*MAXLENGTH-1-n]/freqGA4*freqGA4_lambda2-(1-freqGA[2*MAXLENGTH-1-n])/pow(freqGG4,2)*pow(freqGA4_lambda,2)-(1-freqGA[2*MAXLENGTH-1-n])/freqGG4*freqGA4_lambda2);
        //        A(0,1) += scaleCT[n]*(freqCT[n]/freqCT3*pf3_l1_lambda-freqCT[n]/pow(freqCT3,2)*freqCT3_lambda*pf3_l1-(1-freqCT[n])/freqCC3*pf3_l1_lambda-(1-freqCT[n])/pow(freqCC3,2)*freqCT3_lambda*pf3_l1)+scaleCT[2*MAXLENGTH-1-n]*(freqCT[2*MAXLENGTH-1-n]/freqCT4*pf3_r1_lambda-freqCT[2*MAXLENGTH-1-n]/pow(freqCT4,2)*freqCT4_lambda*pf3_r1-(1-freqCT[2*MAXLENGTH-1-n])/freqCC4*pf3_r1_lambda-(1-freqCT[2*MAXLENGTH-1-n])/pow(freqCC4,2)*freqCT4_lambda*pf3_r1)+scaleGA[n]*(freqGA[n]/freqGA3*pf1_r1_lambda-freqGA[n]/pow(freqGA3,2)*freqGA3_lambda*pf1_r1-(1-freqGA[n])/freqGG3*pf1_r1_lambda-(1-freqGA[n])/pow(freqGG3,2)*freqGA3_lambda*pf1_r1)+scaleGA[2*MAXLENGTH-1-n]*(freqGA[2*MAXLENGTH-1-n]/freqGA4*pf1_l1_lambda-freqGA[2*MAXLENGTH-1-n]/pow(freqGA4,2)*freqGA4_lambda*pf1_l1-(1-freqGA[2*MAXLENGTH-1-n])/freqGG4*pf1_l1_lambda-(1-freqGA[2*MAXLENGTH-1-n])/pow(freqGG4,2)*freqGA4_lambda*pf1_l1);
        //        A(0,2) += scaleCT[n]*(freqCT[n]/freqCT3*pf4_l1_lambda-freqCT[n]/pow(freqCT3,2)*freqCT3_lambda*pf4_l1-(1-freqCT[n])/freqCC3*pf4_l1_lambda-(1-freqCT[n])/pow(freqCC3,2)*freqCT3_lambda*pf4_l1)+scaleGA[n]*(freqGA[n]/freqGA3*pf2_r1_lambda-freqGA[n]/pow(freqGA3,2)*freqGA3_lambda*pf2_r1-(1-freqGA[n])/freqGG3*pf2_r1_lambda-(1-freqGA[n])/pow(freqGG3,2)*freqGA3_lambda*pf2_r1);
        //        A(0,3) +=scaleCT[n]*(freqCT[n]/freqCT3*freqCT3_lambda_nv-freqCT[n]/pow(freqCT3,2)*freqCT3_lambda*freqCT3_nv-(1-freqCT[n])/freqCC3*freqCT3_lambda_nv-(1-freqCT[n])/pow(freqCC3,2)*freqCT3_lambda*freqCT3_nv)+scaleCT[2*MAXLENGTH-1-n]*(freqCT[2*MAXLENGTH-1-n]/freqCT4*freqCT4_lambda_nv-freqCT[2*MAXLENGTH-1-n]/pow(freqCT4,2)*freqCT4_lambda*freqCT4_nv-(1-freqCT[2*MAXLENGTH-1-n])/freqCC4*freqCT4_lambda_nv-(1-freqCT[2*MAXLENGTH-1-n])/pow(freqCC4,2)*freqCT4_lambda*freqCT4_nv)+scaleGA[n]*(freqGA[n]/freqGA3*freqGA3_lambda_nv-freqGA[n]/pow(freqGA3,2)*freqGA3_lambda*freqGA3_nv-(1-freqGA[n])/freqGG3*freqGA3_lambda_nv-(1-freqGA[n])/pow(freqGG3,2)*freqGA3_lambda*freqGA3_nv)+scaleGA[2*MAXLENGTH-1-n]*(freqGA[2*MAXLENGTH-1-n]/freqGA4*freqGA4_lambda_nv-freqGA[2*MAXLENGTH-1-n]/pow(freqGA4,2)*freqGA4_lambda*freqGA4_nv-(1-freqGA[2*MAXLENGTH-1-n])/freqGG4*freqGA4_lambda_nv-(1-freqGA[2*MAXLENGTH-1-n])/pow(freqGG4,2)*freqGA4_lambda*freqGA4_nv);
        //        A(1,1) += scaleCT[n]*(-freqCT[n]/pow(freqCT3,2)*pow(pf3_l1,2)-(1-freqCT[n])/pow(freqCC3,2)*pow(pf3_l1,2))+scaleCT[2*MAXLENGTH-1-n]*(-freqCT[2*MAXLENGTH-1-n]/pow(freqCT4,2)*pow(pf3_r1,2)-(1-freqCT[2*MAXLENGTH-1-n])/pow(freqCC4,2)*pow(pf3_r1,2))+scaleGA[n]*(-freqGA[n]/pow(freqGA3,2)*pow(pf1_r1,2)-(1-freqGA[n])/pow(freqGG3,2)*pow(pf1_r1,2))+scaleGA[2*MAXLENGTH-1-n]*(-freqGA[2*MAXLENGTH-1-n]/pow(freqGA4,2)*pow(pf1_l1,2)-(1-freqGA[2*MAXLENGTH-1-n])/pow(freqGG4,2)*pow(pf1_l1,2));
        //        A(1,2) += scaleCT[n]*(-freqCT[n]/pow(freqCT3,2)*pf3_l1*pf4_l1-(1-freqCT[n])/pow(freqCC3,2)*pf3_l1*pf4_l1)+scaleGA[n]*(-freqGA[n]/pow(freqGA3,2)*pf1_r1*pf2_r1-(1-freqGA[n])/pow(freqGG3,2)*pf1_r1*pf2_r1);
        //        A(1,3) += scaleCT[n]*(freqCT[n]/freqCT3*pf3_l1_nv-freqCT[n]/pow(freqCT3,2)*pf3_l1*freqCT3_nv-(1-freqCT[n])/freqCC3*pf3_l1_nv-(1-freqCT[n])/pow(freqCC3,2)*pf3_l1*freqCT3_nv)+scaleCT[2*MAXLENGTH-1-n]*(freqCT[2*MAXLENGTH-1-n]/freqCT4*pf3_r1_nv-freqCT[2*MAXLENGTH-1-n]/pow(freqCT4,2)*pf3_r1*freqCT4_nv-(1-freqCT[2*MAXLENGTH-1-n])/freqCC4*pf3_r1_nv-(1-freqCT[2*MAXLENGTH-1-n])/pow(freqCC4,2)*pf3_r1*freqCT4_nv)+scaleGA[n]*(freqGA[n]/freqGA3*pf1_r1_nv-freqGA[n]/pow(freqGA3,2)*pf1_r1*freqGA3_nv-(1-freqGA[n])/freqGG3*pf1_r1_nv-(1-freqGA[n])/pow(freqGG3,2)*pf1_r1*freqGA3_nv)+scaleGA[2*MAXLENGTH-1-n]*(freqGA[2*MAXLENGTH-1-n]/freqGA4*pf1_l1_nv-freqGA[2*MAXLENGTH-1-n]/pow(freqGA4,2)*pf1_l1*freqGA4_nv-(1-freqGA[2*MAXLENGTH-1-n])/freqGG4*pf1_l1_nv-(1-freqGA[2*MAXLENGTH-1-n])/pow(freqGG4,2)*pf1_l1*freqGA4_nv);
        //        A(2,2) += scaleCT[n]*(-freqCT[n]/pow(freqCT3,2)*pow(pf4_l1,2)-(1-freqCT[n])/pow(freqCC3,2)*pow(pf4_l1,2))+scaleGA[n]*(-freqGA[n]/pow(freqGA3,2)*pow(pf2_r1,2)-(1-freqGA[n])/pow(freqGG3,2)*pow(pf2_r1,2));
        //        A(2,3) += scaleCT[n]*(freqCT[n]/freqCT3*pf4_l1_nv-freqCT[n]/pow(freqCT3,2)*pf4_l1*freqCT3_nv-(1-freqCT[n])/freqCC3*pf4_l1_nv-(1-freqCT[n])/pow(freqCC3,2)*pf4_l1*freqCT3_nv)+scaleGA[n]*(freqGA[n]/freqGA3*pf2_r1_nv-freqGA[n]/pow(freqGA3,2)*pf2_r1*freqGA3_nv-(1-freqGA[n])/(freqGG3)*pf2_r1_nv-(1-freqGA[n])/pow(freqGG3,2)*pf2_r1*freqGA3_nv);
        //        A(3,3) += scaleCT[n]*(freqCT[n]/freqCT3*freqCT3_nv2-freqCT[n]/pow(freqCT3,2)*pow(freqCT3_nv,2)-(1-freqCT[n])/freqCC3*freqCT3_nv2-(1-freqCT[n])/pow(freqCC3,2)*pow(freqCT3_nv,2))+scaleCT[2*MAXLENGTH-1-n]*(freqCT[2*MAXLENGTH-1-n]/freqCT4*freqCT4_nv2-freqCT[2*MAXLENGTH-1-n]/pow(freqCT4,2)*pow(freqCT4_nv,2)-(1-freqCT[2*MAXLENGTH-1-n])/freqCC4*freqCT4_nv2-(1-freqCT[2*MAXLENGTH-1-n])/pow(freqCC4,2)*pow(freqCT4_nv,2))+scaleGA[n]*(freqGA[n]/freqGA3*freqGA3_nv2-freqGA[n]/pow(freqGA3,2)*pow(freqGA3_nv,2)-(1-freqGA[n])/freqGG3*freqGA3_nv2-(1-freqGA[n])/pow(freqGG3,2)*pow(freqGA3_nv,2))+scaleGA[2*MAXLENGTH-1-n]*(freqGA[2*MAXLENGTH-1-n]/freqGA4*freqGA4_nv2-freqGA[2*MAXLENGTH-1-n]/pow(freqGA4,2)*pow(freqGA4_nv,2)-(1-freqGA[2*MAXLENGTH-1-n])/freqGG4*freqGA4_nv2-(1-freqGA[2*MAXLENGTH-1-n])/pow(freqGG4,2)*pow(freqGA4_nv,2));
    }
    A(1,0) = A(0,1); A(2,0) = A(0,2); A(3,0) = A(0,3); A(2,1) = A(1,2); A(3,1) = A(1,3); A(3,2) = A(2,3);
    //cout<<A<<"\n\n";
    B = -A.inverse();
    //cout<<B<<"\n";
    for (int i=0;i<4;i++){
        for (int j=0;j<4;j++){
            z[i][j] = B(i,j);
        }
    }
}

void loglike_complex3_hessian_full_nb(const double *x, double ** z, double * freqCT, double * freqGA, double * scaleCT, double * scaleGA, double * seqError, int BinNum, double * LEN, double * freqLEN, double eps){
    Eigen::Matrix4d A, B;
    A = Eigen::Matrix4d::Zero();
    double lambda = x[0];
    double delta = x[1];
    double delta_s = x[2];
    double nu = x[3];
    
    for(int n=0; n<MAXLENGTH; n++){
        double freqCT1 = 0;
        double freqCT1_lambda = 0;
        double freqCT1_nv = 0;
        double freqCT1_lambda2 = 0;
        double freqCT1_lambda_nv= 0;
        double freqCT1_nv2= 0;
        double freqCT2 = 0;
        double freqCT2_lambda = 0;
        double freqCT2_nv = 0;
        double freqCT2_lambda2 = 0;
        double freqCT2_lambda_nv= 0;
        double freqCT2_nv2= 0;

        double pf3_l = 0;
        double pf3_l_lambda = 0;
        double pf3_l_nv = 0;
        double pf4_l = 0;
        double pf4_l_lambda = 0;
        double pf4_l_nv = 0;
        double pf1_r = 0;
        double pf1_r_lambda = 0;
        double pf1_r_nv = 0;
        double pf2_r = 0;
        double pf2_r_lambda = 0;
        double pf2_r_nv = 0;
        double pf3_r = 0;
        double pf3_r_lambda = 0;
        double pf3_r_nv = 0;
        double pf1_l = 0;
        double pf1_l_lambda = 0;
        double pf1_l_nv = 0;
        for(int i=0; i<BinNum; i++){
            double L = LEN[i];
            double f = freqLEN[i];
            double p1_l = pow(1+lambda,2)*n*nu/(1+(L-2)*nu); //double p1_l = pow(1+lambda,2)*(n+1)*nu/(1+(L-2)*nu);
            double p1_l_lambda = 2*(1+lambda)*n*nu/(1+(L-2)*nu);
            double p1_l_nv = pow(1+lambda,2)*n/pow(1+(L-2)*nu,2);
            double p1_l_lambda2 = 2*n*nu/(1+(L-2)*nu);
            double p1_l_lambda_nv = 2*(1+lambda)*n/pow(1+(L-2)*nu,2);
            double p1_l_nv2 = -2*pow(1+lambda,2)*n*(L-2)/pow(1+(L-2)*nu,3);
            for(int l=1;l<=min(MAXORDER,n);l++){
                p1_l += lambda*(1+lambda)*pow(1-lambda,l)*(n-l)*nu/(1+(L-2-l)*nu); //p1_l += lambda*(1+lambda)*pow(1-lambda,l)*(n-l)*nu/(1+(L-2-l)*nu);
                p1_l_lambda += ((1+2*lambda)*pow(1-lambda,l)-l*lambda*(1+lambda)*pow(1-lambda,l-1))*(n-l)*nu/(1+(L-2-l)*nu);
                p1_l_nv += lambda*(1+lambda)*pow(1-lambda,l)*(n-l)/pow(1+(L-2-l)*nu,2);
                p1_l_lambda2 += (2*pow(1-lambda,l)-2*l*(1+2*lambda)*pow(1-lambda,l-1)+l*(l-1)*lambda*(1+lambda)*pow(1-lambda,l-2))*(n-l)*nu/(1+(L-2-l)*nu);
                p1_l_lambda_nv += ((1+2*lambda)*pow(1-lambda,l)-l*lambda*(1+lambda)*pow(1-lambda,l-1))*(n-l)/pow(1+(L-2-l)*nu,2);
                p1_l_nv2 += -2*lambda*(1+lambda)*pow(1-lambda,l)*(n-l)*(L-2-l)/pow(1+(L-2-l)*nu,3);
                for (int r=1;r<=min((double)MAXORDER,(double)(L-n-1));r++){
                    if (l+r<=MAXORDER){
                        p1_l += pow(lambda,2)*pow(1-lambda,l+r)*(n-l)*nu/(1+(L-2-l-r)*nu); //p1_l += pow(lambda,2)*pow(1-lambda,l+r)*(n-l)*nu/(1+(L-2-l-r)*nu);
                        p1_l_lambda += (2*lambda*pow(1-lambda,l+r)-(l+r)*pow(lambda,2)*pow(1-lambda,l+r-1))*(n-l)*nu/(1+(L-2-l-r)*nu);
                        p1_l_nv += pow(lambda,2)*pow(1-lambda,l+r)*(n-l)/pow(1+(L-2-l-r)*nu,2);
                        p1_l_lambda2 += (2*pow(1-lambda,l+r)-4*(l+r)*lambda*pow(1-lambda,l+r-1)+(l+r)*(l+r-1)*pow(lambda,2)*pow(1-lambda,l+r-2))*(n-l)*nu/(1+(L-2-l-r)*nu);
                        p1_l_lambda_nv += (2*lambda*pow(1-lambda,l+r)-(l+r)*pow(lambda,2)*pow(1-lambda,l+r-1))*(n-l)/pow(1+(L-2-l-r)*nu,2);
                        p1_l_nv2 += -2*pow(lambda,2)*pow(1-lambda,l+r)*(n-l)*(L-2-l-r)/pow(1+(L-2-l-r)*nu,3);
                    }
                }
            }
            for(int r=1;r<=min((double)MAXORDER,(double)(L-n-1));r++){
                p1_l += lambda*(1+lambda)*pow(1-lambda,r)*n*nu/(1+(L-2-r)*nu); // p1_l += lambda*(1+lambda)*pow(1-lambda,r)*n*nu/(1+(L-2-r)*nu);
                p1_l_lambda += ((1+2*lambda)*pow(1-lambda,r)-r*lambda*(1+lambda)*pow(1-lambda,r-1))*n*nu/(1+(L-2-r)*nu);
                p1_l_nv += lambda*(1+lambda)*pow(1-lambda,r)*n/pow(1+(L-2-r)*nu,2);
                p1_l_lambda2 += (2*pow(1-lambda,r)-2*r*(1+2*lambda)*pow(1-lambda,r-1)+r*(r-1)*lambda*(1+lambda)*pow(1-lambda,r-2))*n*nu/(1+(L-2-r)*nu);
                p1_l_lambda_nv += ((1+2*lambda)*pow(1-lambda,r)-r*lambda*(1+lambda)*pow(1-lambda,r-1))*n/pow(1+(L-2-r)*nu,2);
                p1_l_nv2 += -2*lambda*(1+lambda)*pow(1-lambda,r)*n*(L-2-r)/pow(1+(L-2-r)*nu,3);
            }
            p1_l = p1_l/4;
            p1_l_lambda = p1_l_lambda/4;
            p1_l_nv = p1_l_nv/4;
            p1_l_lambda2 = p1_l_lambda2/4;
            p1_l_lambda_nv = p1_l_lambda_nv/4;
            p1_l_nv2 = p1_l_nv2/4;
            double p2_l = 0;
            double p2_l_lambda = 0;
            double p2_l_nv = 0;
            double p2_l_lambda2 = 0;
            double p2_l_lambda_nv = 0;
            double p2_l_nv2 = 0;
            double p3_l = 1 - pow(1-lambda,n+1)/2-p1_l;
            double p3_l_lambda = (n+1)*pow(1-lambda,n)/2-p1_l_lambda;
            double p3_l_nv = -p1_l_nv;
            double p3_l_lambda2 = -(n+1)*n*pow(1-lambda,n-1)/2-p1_l_lambda2;
            double p3_l_lambda_nv = -p1_l_lambda_nv;
            double p3_l_nv2 = -p1_l_nv2;
            double p4_l = pow(1-lambda,n+1)/2;
            double p4_l_lambda = -(n+1)*pow(1-lambda,n)/2;
            double p4_l_nv = 0;
            double p4_l_lambda2 = (n+1)*n*pow(1-lambda,n-1)/2;
            double p4_l_lambda_nv = 0;
            double p4_l_nv2 = 0;
            
            double p1_r = pow(1+lambda,2)*(L-n-1)*nu/(1+(L-2)*nu); // p1_r = pow(1+lambda,2)*(L-n)*nu/(1+(L-2)*nu);
            double p1_r_lambda = 2*(1+lambda)*(L-n-1)*nu/(1+(L-2)*nu);
            double p1_r_nv = pow(1+lambda,2)*(L-n-1)/pow(1+(L-2)*nu,2);
            double p1_r_lambda2 = 2*(L-n-1)*nu/(1+(L-2)*nu);
            double p1_r_lambda_nv = 2*(1+lambda)*(L-n-1)/pow(1+(L-2)*nu,2);
            double p1_r_nv2 = -2*pow(1+lambda,2)*(L-n-1)*(L-2)/pow(1+(L-2)*nu,3);
            for(int r=1;r<=min(MAXORDER,n);r++){
                p1_r += lambda*(1+lambda)*pow(1-lambda,r)*(L-n-1)*nu/(1+(L-2-r)*nu); // p1_r += lambda*(1+lambda)*pow(1-lambda,r)*(L-n-1-r)*nu/(1+(L-2-r)*nu);
                p1_r_lambda += ((1+2*lambda)*pow(1-lambda,r)-r*lambda*(1+lambda)*pow(1-lambda,r-1))*(L-n-1)*nu/(1+(L-2-r)*nu);
                p1_r_nv += lambda*(1+lambda)*pow(1-lambda,r)*(L-n-1)/pow(1+(L-2-r)*nu,2);
                p1_r_lambda2 += (2*pow(1-lambda,r)-2*r*(1+2*lambda)*pow(1-lambda,r-1)+r*(r-1)*lambda*(1+lambda)*pow(1-lambda,r-2))*(L-n-1)*nu/(1+(L-2-r)*nu);
                p1_r_lambda_nv += ((1+2*lambda)*pow(1-lambda,r)-r*lambda*(1+lambda)*pow(1-lambda,r-1))*(L-n-1)/pow(1+(L-2-r)*nu,2);
                p1_r_nv2 += -2*lambda*(1+lambda)*pow(1-lambda,r)*(L-n-1)*(L-2-r)/pow(1+(L-2-r)*nu,3);
                for (int l=1;l<=min((double)MAXORDER,(double)(L-n-1));l++){
                    if (l+r<=MAXORDER){
                        p1_r += pow(lambda,2)*pow(1-lambda,l+r)*(L-n-1-l)*nu/(1+(L-2-l-r)*nu); // p1_r += pow(lambda,2)*pow(1-lambda,l+r)*(L-n-1-r)*nu/(1+(L-2-l-r)*nu);
                        p1_r_lambda += (2*lambda*pow(1-lambda,l+r)-(l+r)*pow(lambda,2)*pow(1-lambda,l+r-1))*(L-n-1-l)*nu/(1+(L-2-l-r)*nu);
                        p1_r_nv += pow(lambda,2)*pow(1-lambda,l+r)*(L-n-1-l)/pow(1+(L-2-l-r)*nu,2);
                        p1_r_lambda2 += (2*pow(1-lambda,l+r)-4*(l+r)*lambda*pow(1-lambda,l+r-1)+(l+r)*(l+r-1)*pow(lambda,2)*pow(1-lambda,l+r-2))*(L-n-1-l)*nu/(1+(L-2-l-r)*nu);
                        p1_r_lambda_nv += (2*lambda*pow(1-lambda,l+r)-(l+r)*pow(lambda,2)*pow(1-lambda,l+r-1))*(L-n-1-l)/pow(1+(L-2-l-r)*nu,2);
                        p1_r_nv2 += -2*pow(lambda,2)*pow(1-lambda,l+r)*(L-n-1-l)*(L-2-l-r)/pow(1+(L-2-l-r)*nu,3);
                    }
                }
            }
            for(int l=1;l<=min((double)MAXORDER,(double)(L-n-1));l++){
                p1_r += lambda*(1+lambda)*pow(1-lambda,l)*(L-n-1-l)*nu/(1+(L-2-l)*nu); // p1_r += lambda*(1+lambda)*pow(1-lambda,l)*(L-n-1)*nu/(1+(L-2-l)*nu);
                p1_r_lambda += ((1+2*lambda)*pow(1-lambda,l)-l*lambda*(1+lambda)*pow(1-lambda,l-1))*(L-n-1-l)*nu/(1+(L-2-l)*nu);
                p1_r_nv += lambda*(1+lambda)*pow(1-lambda,l)*(L-n-1-l)/pow(1+(L-2-l)*nu,2);
                p1_r_lambda2 += (2*pow(1-lambda,l)-2*l*(1+2*lambda)*pow(1-lambda,l-1)+l*(l-1)*lambda*(1+lambda)*pow(1-lambda,l-2))*(L-n-1-l)*nu/(1+(L-2-l)*nu);
                p1_r_lambda_nv += ((1+2*lambda)*pow(1-lambda,l)-l*lambda*(1+lambda)*pow(1-lambda,l-1))*(L-n-1-l)/pow(1+(L-2-l)*nu,2);
                p1_r_nv2 += -2*lambda*(1+lambda)*pow(1-lambda,l)*(L-n-1-l)*(L-2-l)/pow(1+(L-2-l)*nu,3);
            }
            p1_r = p1_r/4;
            p1_r_lambda = p1_r_lambda/4;
            p1_r_nv = p1_r_nv/4;
            p1_r_lambda2 = p1_r_lambda2/4;
            p1_r_lambda_nv = p1_r_lambda_nv/4;
            p1_r_nv2 = p1_r_nv2/4;
            double p2_r = pow(1-lambda,n+1)/2;
            double p2_r_lambda = -(n+1)*pow(1-lambda,n)/2;
            double p2_r_nv = 0;
            double p2_r_lambda2 = (n+1)*n*pow(1-lambda,n-1)/2;
            double p2_r_lambda_nv = 0;
            double p2_r_nv2 = 0;
            double p3_r = 1 - pow(1-lambda,n+1)/2-p1_r;
            double p3_r_lambda = (n+1)*pow(1-lambda,n)/2-p1_r_lambda;
            double p3_r_nv = -p1_r_nv;
            double p3_r_lambda2 = -(n+1)*n*pow(1-lambda,n-1)/2-p1_r_lambda2;
            double p3_r_lambda_nv = -p1_r_lambda_nv;
            double p3_r_nv2 = -p1_r_nv2;
            double p4_r = 0;
            double p4_r_lambda = 0;
            double p4_r_nv = 0;
            double p4_r_lambda2 = 0;
            double p4_r_lambda_nv = 0;
            double p4_r_nv2 = 0;
            //        cout<<p1_l<<" "<<p2_l<<" "<<p3_l<<" "<<p4_l<<"\n";
            //        cout<<p1_r<<" "<<p2_r<<" "<<p3_r<<" "<<p4_r<<"\n";
            double dd1 = ((p1_r+p3_l)/2*delta+(p2_r+p4_l)/2*delta_s)*f;
            double dd1_lambda = ((p1_r_lambda+p3_l_lambda)/2*delta+(p2_r_lambda+p4_l_lambda)/2*delta_s)*f;
            double dd1_nv = ((p1_r_nv+p3_l_nv)/2*delta+(p2_r_nv+p4_l_nv)/2*delta_s)*f;
            double dd1_lambda2 = ((p1_r_lambda2+p3_l_lambda2)/2*delta+(p2_r_lambda2+p4_l_lambda2)/2*delta_s)*f;
            double dd1_lambda_nv = ((p1_r_lambda_nv+p3_l_lambda_nv)/2*delta+(p2_r_lambda_nv+p4_l_lambda_nv)/2*delta_s)*f;
            double dd1_nv2 = ((p1_r_nv2+p3_l_nv2)/2*delta+(p2_r_nv2+p4_l_nv2)/2*delta_s)*f;
            double dd2 = (p1_l+p3_r)/2*delta*f;
            double dd2_lambda = (p1_l_lambda+p3_r_lambda)/2*delta*f;
            double dd2_nv = (p1_l_nv+p3_r_nv)/2*delta*f;
            double dd2_lambda2 = (p1_l_lambda2+p3_r_lambda2)/2*delta*f;
            double dd2_lambda_nv = (p1_l_lambda_nv+p3_r_lambda_nv)/2*delta*f;
            double dd2_nv2 = (p1_l_nv2+p3_r_nv2)/2*delta*f;
            
            freqCT1 += dd1;
            freqCT1_lambda += dd1_lambda;
            freqCT1_nv += dd1_nv;
            freqCT1_lambda2 += dd1_lambda2;
            freqCT1_lambda_nv += dd1_lambda_nv;
            freqCT1_nv2 += dd1_nv2;
            freqCT2 += dd2;
            freqCT2_lambda += dd2_lambda;
            freqCT2_nv += dd2_nv;
            freqCT2_lambda2 += dd2_lambda2;
            freqCT2_lambda_nv += dd2_lambda_nv;
            freqCT2_nv2 += dd2_nv2;

            pf3_l += p3_l*f;
            pf3_l_lambda += p3_l_lambda*f;
            pf3_l_nv += p3_l_nv*f;
            pf4_l += p4_l*f;
            pf4_l_lambda += p4_l_lambda*f;
            pf4_l_nv += pf4_l_nv*f;
            pf1_r += p1_r*f;
            pf1_r_lambda += p1_r_lambda*f;
            pf1_r_nv += p1_r_nv*f;
            pf2_r += p2_r*f;
            pf2_r_lambda += p2_r_lambda*f;
            pf2_r_nv += p2_r_nv*f;
            pf3_r += p3_r*f;
            pf3_r_lambda += p3_r_lambda*f;
            pf3_r_nv += p3_r_nv*f;
            pf1_l += p1_l*f;
            pf1_l_lambda += p1_l_lambda*f;
            pf1_l_nv += p1_l_nv*f;
        }
        freqCT1 = (1-eps)*freqCT1;
        freqCT1_lambda = (1-eps)*freqCT1_lambda;
        freqCT1_nv = (1-eps)*freqCT1_nv;
        freqCT2 = (1-eps)*freqCT2;
        freqCT2_lambda = (1-eps)*freqCT2_lambda;
        freqCT2_nv = (1-eps)*freqCT2_nv;
        pf3_l = (1-eps)*pf3_l;
        pf4_l = (1-eps)*pf4_l;
        pf1_r = (1-eps)*pf1_r;
        pf2_r = (1-eps)*pf2_r;
        pf3_r = (1-eps)*pf3_r;
        pf1_l = (1-eps)*pf1_l;
        freqCT1_lambda2 = (1-eps)*freqCT1_lambda2;
        freqCT1_lambda_nv = (1-eps)*freqCT1_lambda_nv;
        freqCT1_nv2 = (1-eps)*freqCT1_nv2;
        freqCT2_lambda2 = (1-eps)*freqCT2_lambda2;
        freqCT2_lambda_nv = (1-eps)*freqCT2_lambda_nv;
        freqCT2_nv2 = (1-eps)*freqCT2_nv2;
        pf3_l_lambda = (1-eps)*pf3_l_lambda;
        pf3_l_nv = (1-eps)*pf3_l_nv;
        pf4_l_lambda = (1-eps)*pf4_l_lambda;
        pf4_l_nv = (1-eps)*pf4_l_nv;
        pf1_r_lambda = (1-eps)*pf1_r_lambda;
        pf1_r_nv = (1-eps)*pf1_r_nv;
        pf2_r_lambda = (1-eps)*pf2_r_lambda;
        pf2_r_nv = (1-eps)*pf2_r_nv;
        pf3_r_lambda = (1-eps)*pf3_r_lambda;
        pf3_r_nv = (1-eps)*pf3_r_nv;
        pf1_l_lambda = (1-eps)*pf1_l_lambda;
        pf1_l_nv = (1-eps)*pf1_l_nv;
        
        double freqGA1 = freqCT1;
        double freqGA1_lambda = freqCT1_lambda;
        double freqGA1_nv = freqCT1_nv;
        double freqGA1_lambda2 = freqCT1_lambda2;
        double freqGA1_lambda_nv = freqCT1_lambda_nv;
        double freqGA1_nv2 = freqCT1_nv2;
        double freqGA2 = freqCT2;
        double freqGA2_lambda = freqCT2_lambda;
        double freqGA2_nv = freqCT2_nv;
        double freqGA2_lambda2 = freqCT2_lambda2;
        double freqGA2_lambda_nv = freqCT2_lambda_nv;
        double freqGA2_nv2 = freqCT2_nv2;
        double freqCT3 = freqCT1*(1-seqError[n]) + (1-freqCT1)*seqError[n]/3;
        double freqCC3 = (1-freqCT1)*(1-seqError[n]) + freqCT1*seqError[n]/3;
        double freqCT4 = freqCT2*(1-seqError[2*MAXLENGTH-1-n]) + (1-freqCT2)*seqError[2*MAXLENGTH-1-n]/3;
        double freqCC4 = (1-freqCT2)*(1-seqError[2*MAXLENGTH-1-n]) + freqCT2*seqError[2*MAXLENGTH-1-n]/3;
        double freqGA3 = freqGA1*(1-seqError[2*MAXLENGTH-1-n]) + (1-freqGA1)*seqError[2*MAXLENGTH-1-n]/3;
        double freqGG3 = (1-freqGA1)*(1-seqError[2*MAXLENGTH-1-n]) + freqGA1*seqError[2*MAXLENGTH-1-n]/3;
        double freqGA4 = freqGA2*(1-seqError[n]) + (1-freqGA2)*seqError[n]/3;
        double freqGG4 = (1-freqGA2)*(1-seqError[n]) + freqGA2*seqError[n]/3;
        double freqCT3_lambda = freqCT1_lambda*(1-4.0/3.0*seqError[n]); //double freqCC3_lambda =-freqCT1_lambda*(1-4.0/3.0*seqError[n]);
        double freqCT4_lambda = freqCT2_lambda*(1-4.0/3.0*seqError[2*MAXLENGTH-1-n]); //double freqCC4_lambda =-freqCT2_lambda*(1-4.0/3.0*seqError[2*MAXLENGTH-1-n]);
        double freqGA3_lambda = freqGA1_lambda*(1-4.0/3.0*seqError[2*MAXLENGTH-1-n]); //double freqGG3_lambda =-freqGA1_lambda*(1-4.0/3.0*seqError[2*MAXLENGTH-1-n]);
        double freqGA4_lambda = freqGA2_lambda*(1-4.0/3.0*seqError[n]); //double freqGG4_lambda =-freqGA2_lambda*(1-4.0/3.0*seqError[n]);
        double freqCT3_nv = freqCT1_nv*(1-4.0/3.0*seqError[n]); //double freqCC3_nv =-freqCT1_nv*(1-4.0/3.0*seqError[n]);
        double freqCT4_nv = freqCT2_nv*(1-4.0/3.0*seqError[2*MAXLENGTH-1-n]); //double freqCC4_nv =-freqCT2_nv*(1-4.0/3.0*seqError[2*MAXLENGTH-1-n]);
        double freqGA3_nv = freqGA1_nv*(1-4.0/3.0*seqError[2*MAXLENGTH-1-n]); //double freqGG3_nv =-freqGA1_nv*(1-4.0/3.0*seqError[2*MAXLENGTH-1-n]);
        double freqGA4_nv = freqGA2_nv*(1-4.0/3.0*seqError[n]); //double freqGG4_nv =-freqGA2_nv*(1-4.0/3.0*seqError[n]);
        //        double pf3_l1 = pf3_l*(1-4.0/3.0*seqError[n]);
        //        double pf3_r1 = pf3_r*(1-4.0/3.0*seqError[2*MAXLENGTH-1-n]);
        //        double pf1_r1 = pf1_r*(1-4.0/3.0*seqError[2*MAXLENGTH-1-n]);
        //        double pf1_l1 = pf1_l*(1-4.0/3.0*seqError[n]);
        //        double pf4_l1 = pf4_l*(1-4.0/3.0*seqError[n]);
        //        double pf2_r1 = pf2_r*(1-4.0/3.0*seqError[2*MAXLENGTH-1-n]);
        double pf3_l1 = (pf1_r+pf3_l)/2*(1-4.0/3.0*seqError[n]); // d freqCT3/d delta
        double pf3_r1 = (pf1_l+pf3_r)/2*(1-4.0/3.0*seqError[2*MAXLENGTH-1-n]); // d freqCT4/d delta
        double pf1_r1 = (pf1_r+pf3_l)/2*(1-4.0/3.0*seqError[2*MAXLENGTH-1-n]); // d freqGA3/d delta
        double pf1_l1 = (pf1_l+pf3_r)/2*(1-4.0/3.0*seqError[n]); // d freqGA4/d delta
        double pf4_l1 = (pf2_r+pf4_l)/2*(1-4.0/3.0*seqError[n]); // d freqCT3/d delta_s
        double pf2_r1 = (pf2_r+pf4_l)/2*(1-4.0/3.0*seqError[2*MAXLENGTH-1-n]); // d freqGA3/ d delta_s
        double freqCT3_lambda2 = freqCT1_lambda2*(1-4.0/3.0*seqError[n]);
        double freqCT4_lambda2 = freqCT2_lambda2*(1-4.0/3.0*seqError[2*MAXLENGTH-1-n]);
        double freqGA3_lambda2 = freqGA1_lambda2*(1-4.0/3.0*seqError[2*MAXLENGTH-1-n]); //double freqGG3_lambda =-freqGA1_lambda*(1-4.0/3.0*seqError[2*MAXLENGTH-1-n]);
        double freqGA4_lambda2 = freqGA2_lambda2*(1-4.0/3.0*seqError[n]); //double freqGG4_lambda =-freqGA2_lambda*(1-4.0/3.0*seqError[n]);
        double pf3_l1_lambda = (pf1_r_lambda+pf3_l_lambda)/2*(1-4.0/3.0*seqError[n]); // d^2 freqCT3/d delta d lambda
        double pf3_r1_lambda = (pf1_l_lambda+pf3_r_lambda)/2*(1-4.0/3.0*seqError[2*MAXLENGTH-1-n]); // d^2 freqCT4/d delta d lambda
        double pf1_r1_lambda = (pf1_r_lambda+pf3_l_lambda)/2*(1-4.0/3.0*seqError[2*MAXLENGTH-1-n]); // d^2 freqGA3/d delta d lambda
        double pf1_l1_lambda = (pf1_l_lambda+pf3_r_lambda)/2*(1-4.0/3.0*seqError[n]); // d^2 freqGA4/d delta d lambda
        double pf4_l1_lambda = (pf2_r_lambda+pf4_l_lambda)/2*(1-4.0/3.0*seqError[n]); // d^2 freqCT3/d delta_s d lambda
        double pf2_r1_lambda = (pf2_r_lambda+pf4_l_lambda)/2*(1-4.0/3.0*seqError[2*MAXLENGTH-1-n]); // d^2 freqGA3/d delta_s d lambda
        double pf3_l1_nv = (pf1_r_nv+pf3_l_nv)/2*(1-4.0/3.0*seqError[n]); // d^2 freqCT3/d delta d nv
        double pf3_r1_nv = (pf1_l_nv+pf3_r_nv)/2*(1-4.0/3.0*seqError[2*MAXLENGTH-1-n]); // d^2 freqCT4/d delta d nv
        double pf1_r1_nv = (pf1_r_nv+pf3_l_nv)/2*(1-4.0/3.0*seqError[2*MAXLENGTH-1-n]); // d^2 freqGA3/d delta d nv
        double pf1_l1_nv = (pf1_l_nv+pf3_r_nv)/2*(1-4.0/3.0*seqError[n]); // d^2 freqGA4/d delta d nv
        double pf4_l1_nv = (pf2_r_nv+pf4_l_nv)/2*(1-4.0/3.0*seqError[n]); // d^2 freqCT3/d delta_s d nv
        double pf2_r1_nv = (pf2_r_nv+pf4_l_nv)/2*(1-4.0/3.0*seqError[2*MAXLENGTH-1-n]); // d^2 freqGA3/d delta_s d nv
        double freqCT3_lambda_nv = freqCT1_lambda_nv*(1-4.0/3.0*seqError[n]);
        double freqCT4_lambda_nv = freqCT2_lambda_nv*(1-4.0/3.0*seqError[2*MAXLENGTH-1-n]);
        double freqGA3_lambda_nv = freqGA1_lambda_nv*(1-4.0/3.0*seqError[2*MAXLENGTH-1-n]); //double freqGG3_lambda =-freqGA1_lambda*(1-4.0/3.0*seqError[2*MAXLENGTH-1-n]);
        double freqGA4_lambda_nv = freqGA2_lambda_nv*(1-4.0/3.0*seqError[n]); //double freqGG4_lambda =-freqGA2_lambda*(1-4.0/3.0*seqError[n]);
        double freqCT3_nv2 = freqCT1_nv2*(1-4.0/3.0*seqError[n]);
        double freqCT4_nv2 = freqCT2_nv2*(1-4.0/3.0*seqError[2*MAXLENGTH-1-n]);
        double freqGA3_nv2 = freqGA1_nv2*(1-4.0/3.0*seqError[2*MAXLENGTH-1-n]); //double freqGG3_lambda =-freqGA1_lambda*(1-4.0/3.0*seqError[2*MAXLENGTH-1-n]);
        double freqGA4_nv2 = freqGA2_nv2*(1-4.0/3.0*seqError[n]); //double freqGG4_lambda =-freqGA2_lambda*(1-4.0/3.0*seqError[n]);
        double freq[8], freq_lambda[8], freq_delta[8], freq_deltas[8], freq_nv[8], count[8], freq_lambda2[8], freq_nv2[8], freq_lambda_delta[8], freq_lambda_deltas[8], freq_delta_nv[8], freq_deltas_nv[8], freq_lambda_nv[8];
        freq[0] = freqCT3; freq[1] = freqCC3; freq[2] = freqCT4; freq[3] = freqCC4;
        freq[4] = freqGA3; freq[5] = freqGG3; freq[6] = freqGA4; freq[7] = freqGG4;
        freq_lambda[0] = freqCT3_lambda; freq_lambda[1] = -freqCT3_lambda; freq_lambda[2] = freqCT4_lambda; freq_lambda[3] = -freqCT4_lambda;
        freq_lambda[4] = freqGA3_lambda; freq_lambda[5] = -freqGA3_lambda; freq_lambda[6] = freqGA4_lambda; freq_lambda[7] = -freqGA4_lambda;
        freq_delta[0] = pf3_l1; freq_delta[1] = -pf3_l1; freq_delta[2] = pf3_r1; freq_delta[3] = -pf3_r1;
        freq_delta[4] = pf1_r1; freq_delta[5] = -pf1_r1; freq_delta[6] = pf1_l1; freq_delta[7] = -pf1_l1;
        freq_deltas[0] = pf4_l1; freq_deltas[1] = -pf4_l1; freq_deltas[2] = 0; freq_deltas[3] = 0;
        freq_deltas[4] = pf2_r1; freq_deltas[5] = -pf2_r1; freq_deltas[6] = 0; freq_deltas[7] = 0;
        freq_nv[0] = freqCT3_nv; freq_nv[1] = -freqCT3_nv; freq_nv[2] = freqCT4_nv; freq_nv[3] = -freqCT4_nv;
        freq_nv[4] = freqGA3_nv; freq_nv[5] = -freqGA3_nv; freq_nv[6] = freqGA4_nv; freq_nv[7] = -freqGA4_nv;
        count[0] = scaleCT[n]*freqCT[n]; count[1] = scaleCT[n]*(1-freqCT[n]);
        count[2] = scaleCT[2*MAXLENGTH-1-n]*freqCT[2*MAXLENGTH-1-n]; count[3] = scaleCT[2*MAXLENGTH-1-n]*(1-freqCT[2*MAXLENGTH-1-n]);
        count[4] = scaleGA[n]*freqGA[n]; count[5] = scaleGA[n]*(1-freqGA[n]);
        count[6] = scaleGA[2*MAXLENGTH-1-n]*freqGA[2*MAXLENGTH-1-n]; count[7] = scaleGA[2*MAXLENGTH-1-n]*(1-freqGA[2*MAXLENGTH-1-n]);
        freq_lambda2[0] = freqCT3_lambda2; freq_lambda2[1] = -freqCT3_lambda2; freq_lambda2[2] = freqCT4_lambda2; freq_lambda2[3] = -freqCT4_lambda2;
        freq_lambda2[4] = freqGA3_lambda2; freq_lambda2[5] = -freqGA3_lambda2; freq_lambda2[6] = freqGA4_lambda2; freq_lambda2[7] = -freqGA4_lambda2;
        freq_nv2[0] = freqCT3_nv2; freq_nv2[1] = -freqCT3_nv2; freq_nv2[2] = freqCT4_nv2; freq_nv2[3] = -freqCT4_nv2;
        freq_nv2[4] = freqGA3_nv2; freq_nv2[5] = -freqGA3_nv2; freq_nv2[6] = freqGA4_nv2; freq_nv2[7] = -freqGA4_nv2;
        freq_lambda_delta[0] = pf3_l1_lambda; freq_lambda_delta[1] = -pf3_l1_lambda; freq_lambda_delta[2] = pf3_r1_lambda; freq_lambda_delta[3] = -pf3_r1_lambda;
        freq_lambda_delta[4] = pf1_r1_lambda; freq_lambda_delta[5] = -pf1_r1_lambda; freq_lambda_delta[6] = pf1_l1_lambda; freq_lambda_delta[7] = -pf1_l1_lambda;
        freq_lambda_deltas[0] = pf4_l1_lambda; freq_lambda_deltas[1] = -pf4_l1_lambda; freq_lambda_deltas[2] = 0; freq_lambda_deltas[3] = 0;
        freq_lambda_deltas[4] = pf2_r1_lambda; freq_lambda_deltas[5] = -pf2_r1_lambda; freq_lambda_deltas[6] = 0; freq_lambda_deltas[7] = 0;
        freq_delta_nv[0] = pf3_l1_nv; freq_delta_nv[1] = -pf3_l1_nv; freq_delta_nv[2] = pf3_r1_nv; freq_delta_nv[3] = -pf3_r1_nv;
        freq_delta_nv[4] = pf1_r1_nv; freq_delta_nv[5] = -pf1_r1_nv; freq_delta_nv[6] = pf1_l1_nv; freq_delta_nv[7] = -pf1_l1_nv;
        freq_deltas_nv[0] = pf4_l1_nv; freq_deltas_nv[1] = -pf4_l1_nv; freq_deltas_nv[2] = 0; freq_deltas_nv[3] = 0;
        freq_deltas_nv[4] = pf2_r1_nv; freq_deltas_nv[5] = -pf2_r1_nv; freq_deltas_nv[6] = 0; freq_deltas_nv[7] = 0;
        freq_lambda_nv[0] = freqCT3_lambda_nv; freq_lambda_nv[1] = -freqCT3_lambda_nv; freq_lambda_nv[2] = freqCT4_lambda_nv; freq_lambda_nv[3] = -freqCT4_lambda_nv;
        freq_lambda_nv[4] = freqGA3_lambda_nv; freq_lambda_nv[5] = -freqGA3_lambda_nv; freq_lambda_nv[6] = freqGA4_lambda_nv; freq_lambda_nv[7] = -freqGA4_lambda_nv;
        for (int j=0;j<8;j++){
            if (freq[j]>0){
                A(0,0) += -count[j]/pow(freq[j],2)*pow(freq_lambda[j],2)+count[j]/freq[j]*freq_lambda2[j];
                A(0,1) += -count[j]/pow(freq[j],2)*freq_lambda[j]*freq_delta[j]+count[j]/freq[j]*freq_lambda_delta[j];
                A(0,2) += -count[j]/pow(freq[j],2)*freq_lambda[j]*freq_deltas[j]+count[j]/freq[j]*freq_lambda_deltas[j];
                A(0,3) += -count[j]/pow(freq[j],2)*freq_lambda[j]*freq_nv[j]+count[j]/freq[j]*freq_lambda_nv[j];
                A(1,1) += -count[j]/pow(freq[j],2)*pow(freq_delta[j],2);
                A(1,2) += -count[j]/pow(freq[j],2)*freq_delta[j]*freq_deltas[j];
                A(1,3) += -count[j]/pow(freq[j],2)*freq_delta[j]*freq_nv[j]+count[j]/freq[j]*freq_delta_nv[j];
                A(2,2) += -count[j]/pow(freq[j],2)*pow(freq_deltas[j],2);
                A(2,3) += -count[j]/pow(freq[j],2)*freq_deltas[j]*freq_nv[j]+count[j]/freq[j]*freq_deltas_nv[j];
                A(3,3) += -count[j]/pow(freq[j],2)*pow(freq_nv[j],2)+count[j]/freq[j]*freq_nv2[j];
            }
        }
    }
    A(1,0) = A(0,1); A(2,0) = A(0,2); A(3,0) = A(0,3); A(2,1) = A(1,2); A(3,1) = A(1,3); A(3,2) = A(2,3);
    //cout<<A<<"\n\n";
    B = -A.inverse();
    //cout<<B<<"\n";
    for (int i=0;i<4;i++){
        for (int j=0;j<4;j++){
            z[i][j] = B(i,j);
        }
    }
}


double like_master(const double *xs,const void *){
    //  fprintf(stderr,"like_master\n");
    int nThreads = tsk_nthreads;
    for(int i=0;i<nThreads;i++)
        for(int ii=0;ii<4;ii++)
            my_tsk_struct[i].x[ii] = xs[ii];
    cout<<"mu_anc\tsig_anc\tmu_mod\tsig_mod\n";
    cout<<xs[0]<<"\t"<<xs[1]<<"\t"<<xs[2]<<"\t"<<xs[3]<<"\n";
    pthread_t thd[nThreads];
    for(size_t i=0;i<nThreads;i++){
        int rc = pthread_create(&thd[i],NULL,tsk_All_loglike_recalibration_slave,(void*) i);
        if(rc)
            fprintf(stderr,"Error creating thread\n");
        
    }
    for(int i=0;i<nThreads;i++)
        pthread_join(thd[i], NULL);
    
    double res=0;
    for(int i=0;i<nThreads;i++){
        fprintf(stderr,"lik[%d,%d]=%f\n",my_tsk_struct[i].from,my_tsk_struct[i].to,my_tsk_struct[i].llh_result);
        res += my_tsk_struct[i].llh_result;
    }
    fprintf(stderr,"total lik[%d,%d]: %f\n",my_tsk_struct[0].from,my_tsk_struct[nThreads-1].to,res);
    
    return res;
}

void like_grad_master(const double *xs,double *y,const void *){
    //  fprintf(stderr,"like_master\n");
    int nThreads = tsk_nthreads;
    for(int i=0;i<nThreads;i++)
        for(int ii=0;ii<4;ii++)
            my_tsk_struct[i].x[ii] = xs[ii];
    pthread_t thd[nThreads];
    for(size_t i=0;i<nThreads;i++){
        //int rc = pthread_create(&thd[i],NULL,tsk_All_loglike_recalibration_slave,(void*) i);
        int rc = pthread_create(&thd[i],NULL,tsk_All_loglike_recalibration_grad_slave,(void*) i);
        if(rc)
            fprintf(stderr,"Error creating thread\n");
        
    }
    for(int i=0;i<nThreads;i++)
        pthread_join(thd[i], NULL);
    
    for(int j=0;j<4;j++){
        y[j] = 0;
    }
    for(int i=0;i<nThreads;i++){
        fprintf(stderr,"lik_grad[%d,%d]\n",my_tsk_struct[i].from,my_tsk_struct[i].to);
        for (int j=0;j<4;j++){
            y[j] += my_tsk_struct[i].llh_result_grad[j];
        }
    }
    fprintf(stderr,"total lik_grad[%d,%d]\n",my_tsk_struct[0].from,my_tsk_struct[nThreads-1].to);
    //cout<<y[0]<<"\t"<<y[1]<<"\t"<<y[2]<<"\t"<<y[3]<<"\n";
}
// Check until here.

// Check here???
void like_hess_master(const double *xs,double **y){
    Eigen::Matrix4d A, B;
    //  fprintf(stderr,"like_master\n");???
    int nThreads = tsk_nthreads;
    for(int i=0;i<nThreads;i++)
        for(int ii=0;ii<4;ii++)
            my_tsk_struct[i].x[ii] = xs[ii];
    pthread_t thd[nThreads];
    for(size_t i=0;i<nThreads;i++){
        //int rc = pthread_create(&thd[i],NULL,tsk_All_loglike_recalibration_slave,(void*) i);
        int rc = pthread_create(&thd[i],NULL,tsk_All_loglike_recalibration_hess_slave,(void*) i);
        if(rc)
            fprintf(stderr,"Error creating thread\n");
        
    }
    for(int i=0;i<nThreads;i++)
        pthread_join(thd[i], NULL);
    
    for(int j=0;j<4;j++){
        for (int k=0;k<4;k++){
            A(j,k) = 0;
        }
    }
    for(int i=0;i<nThreads;i++){
        fprintf(stderr,"lik_hess[%d,%d]\n",my_tsk_struct[i].from,my_tsk_struct[i].to);
        for (int j=0;j<4;j++){
            for (int k=0;k<4;k++){
                A(j,k) += my_tsk_struct[i].llh_result_hess[j][k];
            }
            //y[j] += my_tsk_struct[i].llh_result_grad[j];
        }
    }
    //cout<<A<<"\n";
    fprintf(stderr,"total lik_hess[%d,%d]\n",my_tsk_struct[0].from,my_tsk_struct[nThreads-1].to);
    B = -A.inverse();
    for (int i=0;i<4;i++){
        for (int j=0;j<4;j++){
            y[i][j] = B(i,j);
        }
    }
}

double ErrorLik(char reffrag[], char frag[], int L, int seqError[]){
    double l1 = 0;

    for (int i=0; i<l_check;i++){
        //cout<<nuc[(int)reffrag[L-i-1]];
        if (reffrag[i]<4 && frag[i]<4){
            double error1 = PhredError[seqError[i]];
            double error1AThird = PhredErrorAThird[seqError[i]];
            if (reffrag[i]==1 && frag[i]==3){ // C to T change, left part
                l1 += log(error1AThird);
            }else if (reffrag[i]==1 && frag[i]==1){ // C to C, left part
                l1 += log(1-error1);
            }else if(reffrag[i]==2 && frag[i]==0){ // G to A change, left part
                l1 += log(error1AThird);
            }else if (reffrag[i]==2 && frag[i]==2){ // G to G, left part
                l1 += log(1-error1);
            }else if (reffrag[i]==frag[i]){ // All other possible no changes, left part
                l1 += log(1-error1);
            }else if (reffrag[i]!=frag[i]){ // All other possible changes, left part
                l1 += log(error1AThird);
            }
        }
        //cout<<(char)(seqError[i]+33);
        if (reffrag[L-i-1]<4 && frag[L-i-1]<4){
            double error2 = PhredError[seqError[L-i-1]];
            double error2AThird = PhredErrorAThird[seqError[L-i-1]];
            if (reffrag[L-i-1]==1 && frag[L-i-1]==3){ // C to T change, right part
                l1 += log(error2AThird);
            }else if (reffrag[L-i-1]==1 && frag[L-i-1]==1){ // C to C, right part
                l1 += log(1-error2);
            }else if (reffrag[L-i-1]==2 && frag[L-i-1]==0){ // G to A change, right part
                l1 += log(error2AThird);
            }else if (reffrag[L-i-1]==2 && frag[L-i-1]==2){ // G to G, right part
                l1 += log(1-error2);
            }else if (reffrag[L-i-1]!=frag[L-i-1]){ // All other possible no changes, right part
                l1 += log(error2AThird);
            }else if (reffrag[L-i-1]==frag[L-i-1]){ // All other possible changes, right part
                l1 += log(1-error2);
            }
        }
        //cout<<"Position "<<"L-i-1 "<<nuc[(int)reffrag[L-i-1]]<<nuc[(int)frag[L-i-1]]<<" "<<seqError[L-i-1]<<" "<<l1;
    }
    return exp(l1);
}

// Calculate the observation likelihood based on the Ancient model/ biotin model
// additive mode -> multiplicity model + speed up
double PMDLik_b(char reffrag[], char frag[], int L, double lambda, double delta, double delta_s, double nv, int seqError[]){
    double l_pmd=0; //Likelihood
    double p = 0;
    // To change the additive effects to multiplicative effects
    double exterm_s = 0; // to replace delta_s
    double exterm = 0; // to replace delta
    //Investigate each possible left and right 5' overhang pair with (l,r) <= l_check (15)
    for (int l=0; l<=l_check; l++){
        double pl = 0.5*lambda*pow(1-lambda,l)/(1-pow(1-lambda,L-1));
        if (l==0){
            pl += 0.5;
        }
        for (int r=0; r<=l_check; r++){
            double pr = 0.5*lambda*pow(1-lambda,r)/(1-pow(1-lambda,L-1));
            if (r == 0){
                pr += 0.5;
            }
            // Joint probability of (l,r)
            double dp = 4*pl*pr/(3+(1-pow(1-lambda,L-1)-(L-1)*lambda*pow(1-lambda,L-1))/(pow(1-pow(1-lambda,L-1),2)));
            if (dp > Tol){
                // Probability
                double pnick = 0;
                for (int npos=l; npos<=l_check-1; npos++){ //Nick is within left l_check region
                    double l2 = 0;
                    for (int i=0; i<l_check; i++){
                        double error1 = PhredError[seqError[i]];
                        double error1AThird = PhredErrorAThird[seqError[i]];
                        exterm_s = delta_s-4*delta_s*error1AThird;
                        exterm = delta-4*delta*error1AThird;
                        if (reffrag[i] < 4 && frag[i] < 4){
                            if (reffrag[i]==1 && frag[i]==3 && i<=l-1){// C to T change, left single strand part
                                l2 += log(min(error1AThird+exterm_s,MAX0));
                            }else if (reffrag[i]==2 && frag[i]==0 && i<=l-1){ // G to A change, left single strand part
                                l2 += log(min(error1AThird,MAX0));
                            }else if (reffrag[i]==1 && frag[i]==1 && i<=l-1){ // C to C, left single strand part
                                l2 += log(max(1-error1-exterm_s,MIN0));
                            }else if (reffrag[i]==2 && frag[i]==2 && i<=l-1){ // G to G, left single strand part
                                l2 += log(max(1-error1,MIN0));
                            }else if (reffrag[i]==1 && frag[i]==3 && i<=npos){ // C to T change, left double strand part before nick
                                l2 += log(min(error1AThird+exterm,MAX0));
                            }else if (reffrag[i]==2 && frag[i]==0 && i<=npos){ // G to A change, left double strand part before nick
                                l2 += log(min(error1AThird,MAX0));
                            }else if (reffrag[i]==1 && frag[i]==1 && i<=npos){ // C to C, left double strand part before nick
                                l2 += log(max(1-error1-exterm,MIN0));
                            }else if (reffrag[i]==2 && frag[i]==2 && i<=npos){ // G to G, left double strand part before nick
                                l2 += log(max(1-error1,MIN0));
                            }else if (reffrag[i]==1 && frag[i]==3 && i>npos){ // C to T change, left double strand part after nick
                                l2 += log(min(error1AThird,MAX0));
                            }else if (reffrag[i]==2 && frag[i]==0 && i>npos){ // G to A change, left double strand part after nick
                                l2 += log(min(error1AThird+exterm,MAX0));
                            }else if (reffrag[i]==1 && frag[i]==1 && i>npos){ // C to C, left double strand part after nick
                                l2 += log(max(1-error1,MIN0));
                            }else if (reffrag[i]==2 && frag[i]==2 && i>npos){ // G to G, left double strand part after nick
                                l2 += log(max(1-error1-exterm,MIN0));
                            }else if (reffrag[i]==frag[i]){ // All other possible no changes
                                l2 += log(max(1-error1,MIN0));
                            }else if (reffrag[i]!=frag[i]){ // All other possible changes
                                l2 += log(min(error1AThird,MAX0));
                            }
                        }
                        //cout << "l2 " <<l2<<"\n";
                        double error2 = PhredError[seqError[L-i-1]];
                        double error2AThird = PhredErrorAThird[seqError[L-i-1]];
                        exterm_s = delta_s-4*delta_s*error2AThird;
                        exterm = delta-4*delta*error2AThird;
                        if (reffrag[L-i-1] < 4 && frag[L-i-1] < 4){
                            if (reffrag[L-i-1]==2 && frag[L-i-1]==0 && i <= r-1){ // G to A change, right single strand part
                                l2 += log(min(error2AThird+exterm_s,MAX0));
				//cout<<log(min(error2AThird+exterm_s,MAX0))<<"\n";
                            }else if(reffrag[L-i-1]==1 && frag[L-i-1]==3 && i <= r-1){ // C to T change, right single strand part
                                l2 += log(min(error2AThird,MAX0));
                                //cout<<log(min(error2AThird,MAX0))<<"\n";
                            }else if(reffrag[L-i-1]==2 && frag[L-i-1]==2 && i <= r-1){ // G to G, right single strand part
                                l2 += log(max(1-error2-exterm_s,MIN0));
                            }else if(reffrag[L-i-1]==1 && frag[L-i-1]==1 && i <= r-1){ // C to C, right single strand part
                                l2 += log(max(1-error2,MIN0));
                            }else if (reffrag[L-i-1]==2 && frag[L-i-1]==0){ // G to A change, right double strand part
                                l2 += log(min(error2AThird+exterm,MAX0));
                            }else if (reffrag[L-i-1]==1 && frag[L-i-1]==3){ // C to T change, right double strand part
                                l2 += log(min(error2AThird,MAX0));
                            }else if (reffrag[L-i-1]==2 && frag[L-i-1]==2){ // G to G, right double strand part
                                l2 += log(max(1-error2-exterm,MIN0));
                            }else if (reffrag[L-i-1]==1 && frag[L-i-1]==1){ // C to C, right double strand part
                                l2 += log(max(1-error2,MIN0));
                            }else if (reffrag[L-i-1]==frag[L-i-1]){ // All other possible no changes
                                l2 += log(max(1-error2,MIN0));
                            }else if (reffrag[L-i-1]!=frag[L-i-1]){ // All other changes
                                l2 += log(min(error2AThird,MAX0));
                            }
                        }
		//cout<<i<<" "<<nuc[reffrag[i]]<<" "<<nuc[frag[i]]<<" "<<nuc[reffrag[L-i-1]]<<" "<<nuc[frag[L-i-1]]<<" Hello l2 "<<l2<<"\n";
                    }
                    double dpnick = nv/(1+(L-l-r-2)*nv);
		    //cout<<i<<" l2 "<<MAX0<<" "<<MIN0<<" "<<l2<<"\n";
                    l_pmd += exp(l2)*dp*dpnick;
                    pnick += dpnick;
                }
                
                for (int npos = L-l_check-1; npos <= L-r-1; npos++){ //Nick is within right l_check region
                    double l2 = 0;
                    for (int i=0; i<l_check; i++){
                        double error2 = PhredError[seqError[L-i-1]];
                        double error2AThird = PhredErrorAThird[seqError[L-i-1]];
                        exterm_s = delta_s-4*delta_s*error2AThird;
                        exterm = delta-4*delta*error2AThird;
                        if (reffrag[L-i-1] < 4 && frag[L-i-1] < 4){
                            if (reffrag[L-i-1]==2 && frag[L-i-1]==0 && i <= r-1){ // G to A change, right single strand part
                                l2 += log(min(error2AThird+exterm_s,MAX0));
                            }else if(reffrag[L-i-1]==1 && frag[L-i-1]==3 && i <= r-1){ // C to T change, right single strand part
                                l2 += log(min(error2AThird,MAX0));
                            }else if(reffrag[L-i-1]==2 && frag[L-i-1]==2 && i <= r-1){ // G to G, right single strand part
                                l2 += log(max(1-error2-exterm_s,MIN0));
                            }else if(reffrag[L-i-1]==1 && frag[L-i-1]==1 && i <= r-1){ // C to C, right single strand part
                                l2 += log(max(1-error2,MIN0));
                            }else if (reffrag[L-i-1]==2 && frag[L-i-1]==0 && i <= L-npos-1){ // G to A change, right double strand part, right of the nick
                                l2 += log(min(error2AThird+exterm,MAX0));
                            }else if (reffrag[L-i-1]==1 && frag[L-i-1]==3 && i <= L-npos-1){ // C to T change, right double strand part, right of the nick
                                l2 += log(min(error2AThird,MAX0));
                            }else if (reffrag[L-i-1]==2 && frag[L-i-1]==2 && i <= L-npos-1){ // G to G, right double strand part, right of the nick
                                l2 += log(max(1-error2-exterm,MIN0));
                            }else if (reffrag[L-i-1]==1 && frag[L-i-1]==1 && i <= L-npos-1){ // C to C, right double strand part, right of the nick
                                l2 += log(max(1-error2,MIN0));
                            }else if (reffrag[L-i-1]==2 && frag[L-i-1]==0 && i > L-npos-1){ // G to A, right double strand part, left of the nick
                                l2 += log(min(error2AThird,MAX0));
                            }else if (reffrag[L-i-1]==1 && frag[L-i-1]==3 && i > L-npos-1){ // C to T, right double strand part, left of the nick
                                l2 += log(min(error2AThird+exterm,MAX0));
                            }else if (reffrag[L-i-1]==2 && frag[L-i-1]==2 && i > L-npos-1){ // G to G, right double strand part, left of the nick
                                l2 += log(max(1-error2,MIN0));
                            }else if (reffrag[L-i-1]==1 && frag[L-i-1]==1 && i > L-npos-1){ // C to C, right double strand part, left of the nick
                                l2 += log(max(1-error2-exterm,MIN0));
                            }else if (reffrag[L-i-1]==frag[L-i-1]){ // All other possible no changes
                                l2 += log(max(1-error2,MIN0));
                            }else if (reffrag[L-i-1]!=frag[L-i-1]){ // All other possible changes
                                l2 += log(min(error2AThird,MAX0));
                            }
                        }
                        
                        double error1 = PhredError[seqError[i]];
                        double error1AThird = PhredErrorAThird[seqError[i]];
                        exterm_s = delta_s-4*delta_s*error1AThird;
                        exterm = delta-4*delta*error1AThird;
                        if (reffrag[i] < 4 && frag[i] < 4){
                            if (reffrag[i]==1 && frag[i]==3 && i<=l-1){ // C to T change, left single strand part
                                l2 += log(min(error1AThird+exterm_s,MAX0));
                            }else if (reffrag[i]==2 && frag[i]==0 && i<=l-1){ // G to A change, left single strand part
                                l2 += log(min(error1AThird,MAX0));
                            }else if (reffrag[i]==1 && frag[i]==1 && i<=l-1){ // C to C, left single strand part
                                l2 += log(max(1-error1-exterm_s,MIN0));
                            }else if (reffrag[i]==2 && frag[i]==2 && i<=l-1){ // G to G, left single strand part
                                l2 += log(max(1-error1,MIN0));
                            }else if (reffrag[i]==1 && frag[i]==3){ // C to T change, left double strand part
                                l2 += log(min(error1AThird+exterm,MAX0));
                            }else if (reffrag[i]==2 && frag[i]==0){ // G to A change, left double strand part
                                l2 += log(min(error1AThird,MAX0));
                            }else if (reffrag[i]==1 && frag[i]==1){ // C to C, left double strand part
                                l2 += log(max(1-error1-exterm,MIN0));
                            }else if (reffrag[i]==2 && frag[i]==2){ // G to G, left double strand part
                                l2 += log(max(1-error1,MIN0));
                            }else if (reffrag[i]==frag[i]){ // All other possible no changes
                                l2 += log(max(1-error1,MIN0));
                            }else if (reffrag[i]!=frag[i]){ // All other possible changes
                                l2 += log(min(error1AThird,MAX0));
                            }
                        }
                    }
                    double dpnick = npos < L-r-1 ? nv/(1+(L-l-r-2)*nv) : (1-nv)/(1+(L-l-r-2)*nv);
                    l_pmd += exp(l2)*dp*dpnick;
                    pnick += dpnick;
                }
                
                //nick occurring outside both l_check region
                double l2 = 0;
                for (int i=0; i<l_check; i++){
                    double error1 = PhredError[seqError[i]];
                    double error1AThird = PhredErrorAThird[seqError[i]];
                    exterm_s = delta_s-4*delta_s*error1AThird;
                    exterm = delta-4*delta*error1AThird;
                    if (reffrag[i] < 4 && frag[i] < 4){
                        if (reffrag[i]==1 && frag[i]==3 && i<=l-1){ // C to T change, left single strand part
                            l2 += log(min(error1AThird+exterm_s,MAX0));
                        }else if (reffrag[i]==2 && frag[i]==0 && i<=l-1){ // G to A change, left single strand part
                            l2 += log(min(error1AThird,MAX0));
                        }else if (reffrag[i]==1 && frag[i]==1 && i<=l-1){ // C to C, left single strand part
                            l2 += log(max(1-error1-exterm_s,MIN0));
                        }else if (reffrag[i]==2 && frag[i]==2 && i<=l-1){ // G to G, left single strand part
                            l2 += log(max(1-error1,MIN0));
                        }else if (reffrag[i]==1 && frag[i]==3 && i>l-1){ // C to T change, left double strand part
                            l2 += log(min(error1AThird+exterm,MAX0));
                        }else if (reffrag[i]==2 && frag[i]==0 && i>l-1){ // G to A change, left double strand part
                            l2 += log(min(error1AThird,MAX0));
                        }else if (reffrag[i]==1 && frag[i]==1 && i>l-1){ // C to C, left double strand part
                            l2 += log(max(1-error1-exterm,MIN0));
                        }else if (reffrag[i]==2 && frag[i]==2 && i>l-1){ // G to G, left double strand part
                            l2 += log(max(1-error1,MIN0));
                        }else if (reffrag[i]==frag[i]){ // All other possible no changes
                            l2 += log(max(1-error1,MIN0));
                        }else if (reffrag[i]!=frag[i]){ // All other possible changes
                            l2 += log(min(error1AThird,MAX0));
                        }
                    }
                    
                    double error2 = PhredError[seqError[L-i-1]];
                    double error2AThird = PhredErrorAThird[seqError[L-i-1]];
                    exterm_s = delta_s-4*delta_s*error2AThird;
                    exterm = delta-4*delta*error2AThird;
                    if (reffrag[L-i-1] < 4 && frag[L-i-1] < 4){
                        if (reffrag[L-i-1]==2 && frag[L-i-1]==0 && i <= r-1){ // G to A change, right single strand part
                            l2 += log(min(error2AThird+exterm_s,MAX0));
                        }else if(reffrag[L-i-1]==1 && frag[L-i-1]==3 && i <= r-1){ // C to T change, right single strand part
                            l2 += log(min(error2AThird,MAX0));
                        }else if(reffrag[L-i-1]==2 && frag[L-i-1]==2 && i <= r-1){ // G to G, right single strand part
                            l2 += log(max(1-error2-exterm_s,MIN0));
                        }else if(reffrag[L-i-1]==1 && frag[L-i-1]==1 && i <= r-1){ // C to C, right single strand part
                            l2 += log(max(1-error2,MIN0));
                        }else if (reffrag[L-i-1]==2 && frag[L-i-1]==0 && i > r-1){ // G to A change, right double strand part
                            l2 += log(min(error2AThird+exterm,MAX0));
                        }else if (reffrag[L-i-1]==1 && frag[L-i-1]==3 && i > r-1){ // C to T change, right double strand part
                            l2 += log(min(error2AThird,MAX0));
                        }else if (reffrag[L-i-1]==2 && frag[L-i-1]==2 && i > r-1){ // G to G, right double strand part
                            l2 += log(max(1-error2-exterm,MIN0));
                        }else if (reffrag[L-i-1]==1 && frag[L-i-1]==1 && i > r-1){ // C to C, right double strand part
                            l2 += log(max(1-error2,MIN0));
                        }else if (reffrag[L-i-1]==frag[L-i-1]){ // All other possible no changes
                            l2 += log(max(1-error2,MIN0));
                        }else if (reffrag[L-i-1]!=frag[L-i-1]){ // All other possible changes
                            l2 += log(min(error2AThird,MAX0));
                        }
                    }
                }
                l_pmd += exp(l2)*dp*(1-pnick);
            }
            p += dp;
        }
        
        for (int r=l_check+1; r<=L-2-l; r++){ // To guarantee l+r < = L-2
            double pr = 0.5*lambda*pow(1-lambda,r)/(1-pow(1-lambda,L-1));
            if (r == 0){
                pr += 0.5;
            }
            double dp = 4*pl*pr/(3+(1-pow(1-lambda,L-1)-(L-1)*lambda*pow(1-lambda,L-1))/(pow(1-pow(1-lambda,L-1),2)));
            if (dp >= Tol){
                double pnick = 0;
                for (int npos=l; npos<=(l_check-1 < L-r-1 ? l_check-1 : L-r-1); npos++){ //Nick is within left l_check region
                    double l2 = 0;
                    for (int i=0; i<l_check; i++){
                        double error1 = PhredError[seqError[i]];
                        double error1AThird = PhredErrorAThird[seqError[i]];
                        exterm_s = delta_s-4*delta_s*error1AThird;
                        exterm = delta-4*delta*error1AThird;
                        if (reffrag[i] < 4 && frag[i] < 4){
                            if (reffrag[i]==1 && frag[i]==3 && i<=l-1){// C to T change, left single strand part
                                l2 += log(min(error1AThird+exterm_s,MAX0));
                            }else if (reffrag[i]==2 && frag[i]==0 && i<=l-1){ // G to A change, left single strand part
                                l2 += log(min(error1AThird,MAX0));
                            }else if (reffrag[i]==1 && frag[i]==1 && i<=l-1){ // C to C, left single strand part
                                l2 += log(max(1-error1-exterm_s,MIN0));
                            }else if (reffrag[i]==2 && frag[i]==2 && i<=l-1){ // G to G, left single strand part
                                l2 += log(max(1-error1,MIN0));
                            }else if (reffrag[i]==1 && frag[i]==3 && i<=npos){ // C to T change, left double strand part before nick
                                l2 += log(min(error1AThird+exterm,MAX0));
                            }else if (reffrag[i]==2 && frag[i]==0 && i<=npos){ // G to A change, left double strand part before nick
                                l2 += log(min(error1AThird,MAX0));
                            }else if (reffrag[i]==1 && frag[i]==1 && i<=npos){ // C to C, left double strand part before nick
                                l2 += log(max(1-error1-exterm,MIN0));
                            }else if (reffrag[i]==2 && frag[i]==2 && i<=npos){ // G to G, left double strand part before nick
                                l2 += log(max(1-error1,MIN0));
                            }else if (reffrag[i]==1 && frag[i]==3 && i>npos){ // C to T change, left double strand part after nick
                                l2 += log(min(error1AThird,MAX0));
                            }else if (reffrag[i]==2 && frag[i]==0 && i>npos){ // G to A change, left double strand part after nick
                                l2 += log(min(error1AThird+exterm,MAX0));
                            }else if (reffrag[i]==1 && frag[i]==1 && i>npos){ // C to C, left double strand part after nick
                                l2 += log(max(1-error1,MIN0));
                            }else if (reffrag[i]==2 && frag[i]==2 && i>npos){ // G to G, left double strand part after nick
                                l2 += log(max(1-error1-exterm,MIN0));
                            }else if (reffrag[i]==frag[i]){ // All other possible no changes
                                l2 += log(max(1-error1,MIN0));
                            }else if (reffrag[i]!=frag[i]){ // All other possible changes
                                l2 += log(min(error1AThird,MAX0));
                            }
                        }
                        
                        double error2 = PhredError[seqError[L-i-1]];
                        double error2AThird = PhredErrorAThird[seqError[L-i-1]];
                        exterm_s = delta_s-4*delta_s*error2AThird;
                        if (reffrag[L-i-1] < 4 && frag[L-i-1] < 4){
                            if (reffrag[L-i-1]==2 && frag[L-i-1]==0){ // G to A change, right single strand part
                                l2 += log(min(error2AThird+exterm_s,MAX0));
                            }else if(reffrag[L-i-1]==1 && frag[L-i-1]==3){ // C to T change, right single strand part
                                l2 += log(min(error2AThird,MAX0));
                            }else if(reffrag[L-i-1]==2 && frag[L-i-1]==2){ // G to G, right single strand part
                                l2 += log(max(1-error2-exterm_s,MIN0));
                            }else if(reffrag[L-i-1]==1 && frag[L-i-1]==1){ // C to C, right single strand part
                                l2 += log(max(1-error2,MIN0));
                            }else if (reffrag[L-i-1]==frag[L-i-1]){ // All other possible no changes
                                l2 += log(max(1-error2,MIN0));
                            }else if (reffrag[L-i-1]!=frag[L-i-1]){ // All other changes
                                l2 += log(min(error2AThird,MAX0));
                            }
                        }
                    }
                    double dpnick = npos < L-r-1 ? nv/(1+(L-l-r-2)*nv) : (1-nv)/(1+(L-l-r-2)*nv);
                    l_pmd += exp(l2)*dp*dpnick;
                    pnick += dpnick;
                }
                //nick occurring outside both l_check region
                double l2 = 0;
                for (int i=0; i<l_check; i++){
                    double error1 = PhredError[seqError[i]];
                    double error1AThird = PhredErrorAThird[seqError[i]];
                    exterm_s = delta_s-4*delta_s*error1AThird;
                    exterm = delta-4*delta*error1AThird;
                    if (reffrag[i] < 4 && frag[i] < 4){
                        if (reffrag[i]==1 && frag[i]==3 && i<=l-1){ // C to T change, left single strand part
                            l2 += log(min(error1AThird+exterm_s,MAX0));
                        }else if (reffrag[i]==2 && frag[i]==0 && i<=l-1){ // G to A change, left single strand part
                            l2 += log(min(error1AThird,MAX0));
                        }else if (reffrag[i]==1 && frag[i]==1 && i<=l-1){ // C to C, left single strand part
                            l2 += log(max(1-error1-exterm_s,MIN0));
                        }else if (reffrag[i]==2 && frag[i]==2 && i<=l-1){ // G to G, left single strand part
                            l2 += log(max(1-error1,MIN0));
                        }else if (reffrag[i]==1 && frag[i]==3 && i>l-1){ // C to T change, left double strand part
                            l2 += log(min(error1AThird+exterm,MAX0));
                        }else if (reffrag[i]==2 && frag[i]==0 && i>l-1){ // G to A change, left double strand part
                            l2 += log(min(error1AThird,MAX0));
                        }else if (reffrag[i]==1 && frag[i]==1 && i>l-1){ // C to C, left double strand part
                            l2 += log(max(1-error1-exterm,MIN0));
                        }else if (reffrag[i]==2 && frag[i]==2 && i>l-1){ // G to G, left double strand part
                            l2 += log(max(1-error1,MIN0));
                        }else if (reffrag[i]==frag[i]){ // All other possible no changes
                            l2 += log(max(1-error1,MIN0));
                        }else if (reffrag[i]!=frag[i]){ // All other possible changes
                            l2 += log(min(error1AThird,MAX0));
                        }
                    }
                    
                    double error2 = PhredError[seqError[L-i-1]];
                    double error2AThird = PhredErrorAThird[seqError[L-i-1]];
                    exterm_s = delta_s-4*delta_s*error2AThird;
                    if (reffrag[L-i-1] < 4 && frag[L-i-1] < 4){
                        if (reffrag[L-i-1]==2 && frag[L-i-1]==0 && i <= r-1){ // G to A change, right single strand part
                            l2 += log(min(error2AThird+exterm_s,MAX0));
                        }else if(reffrag[L-i-1]==1 && frag[L-i-1]==3 && i <= r-1){ // C to T change, right single strand part
                            l2 += log(min(error2AThird,MAX0));
                        }else if(reffrag[L-i-1]==2 && frag[L-i-1]==2 && i <= r-1){ // G to G, right single strand part
                            l2 += log(max(1-error2-exterm_s,MIN0));
                        }else if(reffrag[L-i-1]==1 && frag[L-i-1]==1 && i <= r-1){ // C to C, right single strand part
                            l2 += log(max(1-error2,MIN0));
                        }else if (reffrag[L-i-1]==frag[L-i-1]){ // All other possible no changes
                            l2 += log(max(1-error2,MIN0));
                        }else if (reffrag[L-i-1]!=frag[L-i-1]){ // All other possible changes
                            l2 += log(min(error2AThird,MAX0));
                        }
                    }
                }
                l_pmd += exp(l2)*dp*(1-pnick);
            }
            p += dp;
        }
    }
    
    for (int r=0; r<=l_check; r++){
        double pr = 0.5*lambda*pow(1-lambda,r)/(1-pow(1-lambda,L-1));
        if (r == 0){
            pr += 0.5;
        }
        for (int l=l_check+1;l<=L-2-r; l++){
            double pl = 0.5*lambda*pow(1-lambda,l)/(1-pow(1-lambda,L-1));
            if (l==0){
                pl += 0.5;
            }
            // Joint probability of (l,r)
            double dp = 4*pl*pr/(3+(1-pow(1-lambda,L-1)-(L-1)*lambda*pow(1-lambda,L-1))/(pow(1-pow(1-lambda,L-1),2)));
            if (dp > Tol){
                double pnick = 0;
                for (int npos = (L-l_check-1 > l ? L-l_check-1 : l); npos <= L-r-1; npos++){ //Nick is within right l_check region
                    double l2 = 0;
                    for (int i=0; i<l_check; i++){
                        double error2 = PhredError[seqError[L-i-1]];
                        double error2AThird = PhredErrorAThird[seqError[L-i-1]];
                        exterm_s = delta_s-4*delta_s*error2AThird;
                        exterm = delta-4*delta*error2AThird;
                        if (reffrag[L-i-1] < 4 && frag[L-i-1] < 4){
                            if (reffrag[L-i-1]==2 && frag[L-i-1]==0 && i <= r-1){ // G to A change, right single strand part
                                l2 += log(min(error2AThird+exterm_s,MAX0));
                            }else if(reffrag[L-i-1]==1 && frag[L-i-1]==3 && i <= r-1){ // C to T change, right single strand part
                                l2 += log(min(error2AThird,MAX0));
                            }else if(reffrag[L-i-1]==2 && frag[L-i-1]==2 && i <= r-1){ // G to G, right single strand part
                                l2 += log(max(1-error2-exterm_s,MIN0));
                            }else if(reffrag[L-i-1]==1 && frag[L-i-1]==1 && i <= r-1){ // C to C, right single strand part
                                l2 += log(max(1-error2,MIN0));
                            }else if (reffrag[L-i-1]==2 && frag[L-i-1]==0 && i <= L-npos-1){ // G to A change, right double strand part, right of the nick
                                l2 += log(min(error2AThird+exterm,MAX0));
                            }else if (reffrag[L-i-1]==1 && frag[L-i-1]==3 && i <= L-npos-1){ // C to T change, right double strand part, right of the nick
                                l2 += log(min(error2AThird,MAX0));
                            }else if (reffrag[L-i-1]==2 && frag[L-i-1]==2 && i <= L-npos-1){ // G to G, right double strand part, right of the nick
                                l2 += log(max(1-error2-exterm,MIN0));
                            }else if (reffrag[L-i-1]==1 && frag[L-i-1]==1 && i <= L-npos-1){ // C to C, right double strand part, right of the nick
                                l2 += log(max(1-error2,MIN0));
                            }else if (reffrag[L-i-1]==2 && frag[L-i-1]==0 && i > L-npos-1){ // G to A, right double strand part, left of the nick
                                l2 += log(min(error2AThird,MAX0));
                            }else if (reffrag[L-i-1]==1 && frag[L-i-1]==3 && i > L-npos-1){ // C to T, right double strand part, left of the nick
                                l2 += log(min(error2AThird+exterm,MAX0));
                            }else if (reffrag[L-i-1]==2 && frag[L-i-1]==2 && i > L-npos-1){ // G to G, right double strand part, left of the nick
                                l2 += log(max(1-error2,MIN0));
                            }else if (reffrag[L-i-1]==1 && frag[L-i-1]==1 && i > L-npos-1){ // C to C, right double strand part, left of the nick
                                l2 += log(max(1-error2-exterm,MIN0));
                            }else if (reffrag[L-i-1]==frag[L-i-1]){ // All other possible no changes
                                l2 += log(max(1-error2,MIN0));
                            }else if (reffrag[L-i-1]!=frag[L-i-1]){ // All other possible changes
                                l2 += log(min(error2AThird,MAX0));
                            }
                        }
                        
                        double error1 = PhredError[seqError[i]];
                        double error1AThird = PhredErrorAThird[seqError[i]];
                        exterm_s = delta_s-4*delta_s*error1AThird;
                        if (reffrag[i] < 4 && frag[i] < 4){
                            if (reffrag[i]==1 && frag[i]==3){ // C to T change, left single strand part
                                l2 += log(min(error1AThird+exterm_s,MAX0));
                            }else if (reffrag[i]==2 && frag[i]==0){ // G to A change, left single strand part
                                l2 += log(min(error1AThird,MAX0));
                            }else if (reffrag[i]==1 && frag[i]==1){ // C to C, left single strand part
                                l2 += log(max(1-error1-exterm_s,MIN0));
                            }else if (reffrag[i]==2 && frag[i]==2){ // G to G, left single strand part
                                l2 += log(max(1-error1,MIN0));
                            }else if (reffrag[i]==frag[i]){ // All other possible no changes
                                l2 += log(max(1-error1,MIN0));
                            }else if (reffrag[i]!=frag[i]){ // All other possible changes
                                l2 += log(min(error1AThird,MAX0));
                            }
                        }
                    }
                    double dpnick = npos < L-r-1 ? nv/(1+(L-l-r-2)*nv) : (1-nv)/(1+(L-l-r-2)*nv);
                    l_pmd += exp(l2)*dp*dpnick;
                    pnick += dpnick;
                }
                // nick occurring outside both l_check region
                double l2 = 0;
                for (int i=0; i<l_check; i++){
                    double error1 = PhredError[seqError[i]];
                    double error1AThird = PhredErrorAThird[seqError[i]];
                    exterm_s = delta_s-4*delta_s*error1AThird;
                    if (reffrag[i] < 4 && frag[i] < 4){
                        if (reffrag[i]==1 && frag[i]==3){ // C to T change, left single strand part
                            l2 += log(min(error1AThird+exterm_s,MAX0));
                        }else if (reffrag[i]==2 && frag[i]==0){ // G to A change, left single strand part
                            l2 += log(min(error1AThird,MAX0));
                        }else if (reffrag[i]==1 && frag[i]==1){ // C to C, left single strand part
                            l2 += log(max(1-error1-exterm_s,MIN0));
                        }else if (reffrag[i]==2 && frag[i]==2){ // G to G, left single strand part
                            l2 += log(max(1-error1,MIN0));
                        }else if (reffrag[i]==frag[i]){ // All other possible no changes
                            l2 += log(max(1-error1,MIN0));
                        }else if (reffrag[i]!=frag[i]){ // All other possible changes
                            l2 += log(min(error1AThird,MAX0));
                        }
                    }
                    
                    double error2 = PhredError[seqError[L-i-1]];
                    double error2AThird = PhredErrorAThird[seqError[L-i-1]];
                    exterm_s = delta_s-4*delta_s*error2AThird;
                    exterm = delta-4*delta*error2AThird;
                    if (reffrag[L-i-1] < 4 && frag[L-i-1] < 4){
                        if (reffrag[L-i-1]==2 && frag[L-i-1]==0 && i <= r-1){ // G to A change, right single strand part
                            l2 += log(min(error2AThird+exterm_s,MAX0));
                        }else if(reffrag[L-i-1]==1 && frag[L-i-1]==3 && i <= r-1){ // C to T change, right single strand part
                            l2 += log(min(error2AThird,MAX0));
                        }else if(reffrag[L-i-1]==2 && frag[L-i-1]==2 && i <= r-1){ // G to G, right single strand part
                            l2 += log(max(1-error2-exterm_s,MIN0));
                        }else if(reffrag[L-i-1]==1 && frag[L-i-1]==1 && i <= r-1){ // C to C, right single strand part
                            l2 += log(max(1-error2,MIN0));
                        }else if (reffrag[L-i-1]==2 && frag[L-i-1]==0 && i > r-1){ // G to A change, right double strand part
                            l2 += log(min(error2AThird+exterm,MAX0));
                        }else if (reffrag[L-i-1]==1 && frag[L-i-1]==3 && i > r-1){ // C to T change, right double strand part
                            l2 += log(min(error2AThird,MAX0));
                        }else if (reffrag[L-i-1]==2 && frag[L-i-1]==2 && i > r-1){ // G to G, right double strand part
                            l2 += log(max(1-error2-exterm,MIN0));
                        }else if (reffrag[L-i-1]==1 && frag[L-i-1]==1 && i > r-1){ // C to C, right double strand part
                            l2 += log(max(1-error2,MIN0));
                        }else if (reffrag[L-i-1]==frag[L-i-1]){ // All other possible no changes
                            l2 += log(max(1-error2,MIN0));
                        }else if (reffrag[L-i-1]!=frag[L-i-1]){ // All other possible changes
                            l2 += log(min(error2AThird,MAX0));
                        }
                    }
                }
                l_pmd += exp(l2)*dp*(1-pnick);
            }
            p += dp;
        }
    }
    
    //Investigate each possible left and right 5' overhang pair with (l,r) > l_check (15)
    if (1-p > Tol){
        double l2 = 0;
        for (int i=0; i<l_check; i++){
            double error1 = PhredError[seqError[i]];
            double error1AThird = PhredErrorAThird[seqError[i]];
            exterm_s = delta_s-4*delta_s*error1AThird;
            if (reffrag[i] < 4 && frag[i] < 4){
                if (reffrag[i]==1 && frag[i]==3){ // C to T change, left single strand part
                    l2 += log(min(error1AThird+delta_s,MAX0));
                }else if (reffrag[i]==2 && frag[i]==0){ // G to A change, left single strand part
                    l2 += log(min(error1AThird,MAX0));
                }else if (reffrag[i]==1 && frag[i]==1){ // C to C, left single strand part
                    l2 += log(max(1-error1-exterm_s,MIN0));
                }else if (reffrag[i]==2 && frag[i]==2){ // G to G, left single strand part
                    l2 += log(max(1-error1,MIN0));
                }else if (reffrag[i]==frag[i]){ // All other possible no changes
                    l2 += log(max(1-error1,MIN0));
                }else if (reffrag[i]!=frag[i]){  // All other possible changes
                    l2 += log(min(error1AThird,MAX0));
                }
            }
            
            double error2 = PhredError[seqError[L-i-1]];
            double error2AThird = PhredErrorAThird[seqError[L-i-1]];
            exterm_s = delta_s-4*delta_s*error2AThird;
            if (reffrag[L-i-1] < 4 && frag[L-i-1] < 4){
                if (reffrag[L-i-1]==2 && frag[L-i-1]==0){ // G to A change, right single strand part
                    l2 += log(min(error2AThird+exterm_s,MAX0));
                }else if(reffrag[L-i-1]==1 && frag[L-i-1]==3){ // C to A change, right single strand part
                    l2 += log(min(error2AThird,MAX0));
                }else if(reffrag[L-i-1]==2 && frag[L-i-1]==2){ // G to G, right single strand part
                    l2 += log(max(1-error2-exterm_s,MIN0));
                }else if(reffrag[L-i-1]==1 && frag[L-i-1]==1){ // C to C, right single strand part
                    l2 += log(max(1-error2,MIN0));
                }else if (reffrag[L-i-1]==frag[L-i-1]){ // All other possible no changes
                    l2 += log(max(1-error2,MIN0));
                }else if (reffrag[L-i-1]!=frag[L-i-1]){ // All other possible changes
                    l2 += log(min(error2AThird,MAX0));
                }
            }
        }
        l_pmd += exp(l2)*(1-p);
    }
    return l_pmd;
}


// The function below is for calculating likelihood of the reverse-complementary strand of the "biotin model" strand, the name is just for simplicity.
// additive -> mulplicative
double PMDLik_nb(char reffrag[], char frag[], int L, double lambda, double delta, double delta_s, double nv, int seqError[]){
    double l_pmd=0; // Likelihood
    double p = 0; // Accumulated prob of (l,r)
    double exterm_s = 0;
    double exterm = 0;
    //cout<<"nuc_llik lambda "<<reffrag[0]+0<<" "<<frag[0]+0<<"\n";
    //Investigate each possible left and right 5' overhang pair with (l,r) <= l_check (15)
    //--------- - - - -                        - - - ------------// --- single strand part
    for (int l=0; l<=l_check; l++){
        double pl = 0.5*lambda*pow(1-lambda,l)/(1-pow(1-lambda,L-1));
        if (l==0){
            pl += 0.5;
        }
        for (int r=0; r<=l_check; r++){
            double pr = 0.5*lambda*pow(1-lambda,r)/(1-pow(1-lambda,L-1));
            if (r == 0){
                pr += 0.5;
            }
            // Joint probability of (l,r)
            double dp = 4*pl*pr/(3+(1-pow(1-lambda,L-1)-(L-1)*lambda*pow(1-lambda,L-1))/(pow(1-pow(1-lambda,L-1),2)));
            if (dp > Tol){
                double pnick = 0; // Probability
                //--------- - - - -                        - - - ------------//
                //            /|\                                            //
                //             |___ Nick                                     //
                for (int npos=l; npos<=l_check-1; npos++){ //Nick is within left l_check region
                    double l2 = 0;
                    for (int i=0; i<l_check; i++){
                        double error1 = PhredError[seqError[L-i-1]]; //rev-comp but error rate should be the same
                        double error1AThird = PhredErrorAThird[seqError[L-i-1]]; //rev-comp
                        exterm_s = delta_s-4*delta_s*error1AThird;
                        exterm = delta-4*delta*error1AThird;
                        if (reffrag[L-i-1] < 4 && frag[L-i-1] < 4){
                            if (reffrag[L-i-1]==2 && frag[L-i-1]==0 && i<=l-1){// C to T change, left single strand part, rev-comp
                                l2 += log(min(error1AThird+exterm_s,MAX0));
                            }else if (reffrag[L-i-1]==1 && frag[L-i-1]==3 && i<=l-1){ // G to A change, left single strand part, rev-comp
                                l2 += log(min(error1AThird,MAX0));
                            }else if (reffrag[L-i-1]==2 && frag[L-i-1]==2 && i<=l-1){ // C to C, left single strand part, rev-comp
                                l2 += log(max(1-error1-exterm_s,MIN0));
                            }else if (reffrag[L-i-1]==1 && frag[L-i-1]==1 && i<=l-1){ // G to G, left single strand part, rev-comp
                                l2 += log(max(1-error1,MIN0));
                            }else if (reffrag[L-i-1]==2 && frag[L-i-1]==0 && i<=npos){ // C to T change, left double strand part before nick, rev-comp
                                l2 += log(min(error1AThird+exterm,MAX0));
                            }else if (reffrag[L-i-1]==1 && frag[L-i-1]==3 && i<=npos){ // G to A change, left double strand part before nick, rev-comp
                                l2 += log(min(error1AThird,MAX0));
                            }else if (reffrag[L-i-1]==2 && frag[L-i-1]==2 && i<=npos){ // C to C, left double strand part before nick, rev-comp
                                l2 += log(max(1-error1-exterm,MIN0));
                            }else if (reffrag[L-i-1]==1 && frag[L-i-1]==1 && i<=npos){ // G to G, left double strand part before nick, rev-comp
                                l2 += log(max(1-error1,MIN0));
                            }else if (reffrag[L-i-1]==2 && frag[L-i-1]==0 && i>npos){ // C to T change, left double strand part after nick, rev-comp
                                l2 += log(min(error1AThird,MAX0));
                            }else if (reffrag[L-i-1]==1 && frag[L-i-1]==3 && i>npos){ // G to A change, left double strand part after nick, rev-comp
                                l2 += log(min(error1AThird+exterm,MAX0));
                            }else if (reffrag[L-i-1]==2 && frag[L-i-1]==2 && i>npos){ // C to C, left double strand part after nick, rev-comp
                                l2 += log(max(1-error1,MIN0));
                            }else if (reffrag[L-i-1]==1 && frag[L-i-1]==1 && i>npos){ // G to G, left double strand part after nick, rev-comp
                                l2 += log(max(1-error1-exterm,MIN0));
                            }else if (reffrag[L-i-1]==frag[L-i-1]){ // All other possible no changes, rev-comp
                                l2 += log(max(1-error1,MIN0));
                            }else if (reffrag[L-i-1]!=frag[L-i-1]){ // All other possible changes, rev-comp
                                l2 += log(min(error1AThird,MAX0));
                            }
                        }
                        
                        double error2 = PhredError[seqError[i]];
                        double error2AThird = PhredErrorAThird[seqError[i]];
                        exterm_s = delta_s-4*delta_s*error2AThird;
                        exterm = delta-4*delta*error2AThird;
                        if (reffrag[i] < 4 && frag[i] < 4){
                            if (reffrag[i]==1 && frag[i]==3 && i <= r-1){ // G to A change, right single strand part, rev-comp
                                l2 += log(min(error2AThird+exterm_s,MAX0));
                            }else if(reffrag[i]==2 && frag[i]==0 && i <= r-1){ // C to T change, right single strand part, rev-comp
                                l2 += log(min(error2AThird,MAX0));
                            }else if(reffrag[i]==1 && frag[i]==1 && i <= r-1){ // G to G, right single strand part, rev-comp
                                l2 += log(max(1-error2-exterm_s,MIN0));
                            }else if(reffrag[i]==2 && frag[i]==2 && i <= r-1){ // C to C, right single strand part, rev-comp
                                l2 += log(max(1-error2,MIN0));
                            }else if (reffrag[i]==1 && frag[i]==3){ // G to A change, right double strand part, rev-comp
                                l2 += log(min(error2AThird+exterm,MAX0));
                            }else if (reffrag[i]==2 && frag[i]==0){ // C to T change, right double strand part, rev-comp
                                l2 += log(min(error2AThird,MAX0));
                            }else if (reffrag[i]==1 && frag[i]==1){ // G to G, right double strand part, rev-comp
                                l2 += log(max(1-error2-exterm,MIN0));
                            }else if (reffrag[i]==2 && frag[i]==2){ // C to C, right double strand part, rev-comp
                                l2 += log(max(1-error2,MIN0));
                            }else if (reffrag[i]==frag[i]){ // All other possible no changes, rev-comp
                                l2 += log(max(1-error2,MIN0));
                            }else if (reffrag[i]!=frag[i]){ // All other changes, rev-comp
                                l2 += log(min(error2AThird,MAX0));
                            }
                        }
                    }
                    double dpnick = nv/(1+(L-l-r-2)*nv);
                    l_pmd += exp(l2)*dp*dpnick;
                    pnick += dpnick;
                }
                
                //--------- - - - -                        - - - ------------//
                //                                          /|\              //
                //                                           |___ Nick       //
                for (int npos = L-l_check-1; npos <= L-r-1; npos++){ //Nick is within right l_check region
                    double l2 = 0;
                    for (int i=0; i<l_check; i++){
                        double error2 = PhredError[seqError[i]];
                        double error2AThird = PhredErrorAThird[seqError[i]];
                        exterm_s = delta_s-4*delta_s*error2AThird;
                        exterm = delta-4*delta*error2AThird;
                        if (reffrag[i] < 4 && frag[i] < 4){
                            if (reffrag[i]==1 && frag[i]==3 && i <= r-1){ // G to A change, right single strand part, rev-comp
                                l2 += log(min(error2AThird+exterm_s,MAX0));
                            }else if(reffrag[i]==2 && frag[i]==0 && i <= r-1){ // C to T change, right single strand part, rev-comp
                                l2 += log(min(error2AThird,MAX0));
                            }else if(reffrag[i]==1 && frag[i]==1 && i <= r-1){ // G to G, right single strand part, rev-comp
                                l2 += log(max(1-error2-exterm_s,MIN0));
                            }else if(reffrag[i]==2 && frag[i]==2 && i <= r-1){ // C to C, right single strand part, rev-comp
                                l2 += log(max(1-error2,MIN0));
                            }else if (reffrag[i]==1 && frag[i]==3 && i <= L-npos-1){ // G to A change, right double strand part, right of the nick, rev-comp
                                l2 += log(min(error2AThird+exterm,MAX0));
                            }else if (reffrag[i]==2 && frag[i]==0 && i <= L-npos-1){ // C to T change, right double strand part, right of the nick, rev-comp
                                l2 += log(min(error2AThird,MAX0));
                            }else if (reffrag[i]==1 && frag[i]==1 && i <= L-npos-1){ // G to G, right double strand part, right of the nick, rev-comp
                                l2 += log(max(1-error2-exterm,MIN0));
                            }else if (reffrag[i]==2 && frag[i]==2 && i <= L-npos-1){ // C to C, right double strand part, right of the nick, rev-comp
                                l2 += log(max(1-error2,MIN0));
                            }else if (reffrag[i]==1 && frag[i]==3 && i > L-npos-1){ // G to A, right double strand part, left of the nick, rev-comp
                                l2 += log(min(error2AThird,MAX0));
                            }else if (reffrag[i]==2 && frag[i]==0 && i > L-npos-1){ // C to T, right double strand part, left of the nick, rev-comp
                                l2 += log(min(error2AThird+exterm,MAX0));
                            }else if (reffrag[i]==1 && frag[i]==1 && i > L-npos-1){ // G to G, right double strand part, left of the nick, rev-comp
                                l2 += log(max(1-error2,MIN0));
                            }else if (reffrag[i]==2 && frag[i]==2 && i > L-npos-1){ // C to C, right double strand part, left of the nick, rev-comp
                                l2 += log(max(1-error2-exterm,MIN0));
                            }else if (reffrag[i]==frag[i]){ // All other possible no changes, rev-comp
                                l2 += log(max(1-error2,MIN0));
                            }else if (reffrag[i]!=frag[i]){ // All other possible changes, rev-comp
                                l2 += log(min(error2AThird,MAX0));
                            }
                        }
                        
                        double error1 = PhredError[seqError[L-i-1]];
                        double error1AThird = PhredErrorAThird[seqError[L-i-1]];
                        exterm_s = delta_s-4*delta_s*error1AThird;
                        exterm = delta-4*delta*error1AThird;
                        if (reffrag[L-i-1] < 4 && frag[L-i-1] < 4){
                            if (reffrag[L-i-1]==2 && frag[L-i-1]==0 && i<=l-1){ // C to T change, left single strand part, rev-comp
                                l2 += log(min(error1AThird+exterm_s,MAX0));
                            }else if (reffrag[L-i-1]==1 && frag[L-i-1]==3 && i<=l-1){ // G to A change, left single strand part, rev-comp
                                l2 += log(min(error1AThird,MAX0));
                            }else if (reffrag[L-i-1]==2 && frag[L-i-1]==2 && i<=l-1){ // C to C, left single strand part, rev-comp
                                l2 += log(max(1-error1-exterm_s,MIN0));
                            }else if (reffrag[L-i-1]==1 && frag[L-i-1]==1 && i<=l-1){ // G to G, left single strand part, rev-comp
                                l2 += log(max(1-error1,MIN0));
                            }else if (reffrag[L-i-1]==2 && frag[L-i-1]==0){ // C to T change, left double strand part, rev-comp
                                l2 += log(min(error1AThird+exterm,MAX0));
                            }else if (reffrag[L-i-1]==1 && frag[L-i-1]==3){ // G to A change, left double strand part, rev-comp
                                l2 += log(min(error1AThird,MAX0));
                            }else if (reffrag[L-i-1]==2 && frag[L-i-1]==2){ // C to C, left double strand part, rev-comp
                                l2 += log(max(1-error1-exterm,MIN0));
                            }else if (reffrag[L-i-1]==1 && frag[L-i-1]==1){ // G to G, left double strand part, rev-comp
                                l2 += log(max(1-error1,MIN0));
                            }else if (reffrag[L-i-1]==frag[L-i-1]){ // All other possible no changes, rev-comp
                                l2 += log(max(1-error1,MIN0));
                            }else if (reffrag[L-i-1]!=frag[L-i-1]){ // All other possible changes, rev-comp
                                l2 += log(min(error1AThird,MIN0));
                            }
                        }
                    }
                    double dpnick = npos < L-r-1 ? nv/(1+(L-l-r-2)*nv) : (1-nv)/(1+(L-l-r-2)*nv);
                    l_pmd += exp(l2)*dp*dpnick;
                    pnick += dpnick;
                }
                //nick occurring outside both l_check region
                //--------- - - - -                        - - - ------------//
                //                         /|\                               //
                //                          |___ Nick                        //
                double l2 = 0;
                for (int i=0; i<l_check; i++){
                    double error1 = PhredError[seqError[L-i-1]];
                    double error1AThird = PhredErrorAThird[seqError[L-i-1]];
                    exterm_s = delta_s-4*delta_s*error1AThird;
                    exterm = delta-4*delta*error1AThird;
                    if (reffrag[L-i-1] < 4 && frag[L-i-1] < 4){
                        if (reffrag[L-i-1]==2 && frag[L-i-1]==0 && i<=l-1){ // C to T change, left single strand part, rev-comp
                            l2 += log(min(error1AThird+exterm_s,MAX0));
                        }else if (reffrag[L-i-1]==1 && frag[L-i-1]==3 && i<=l-1){ // G to A change, left single strand part, rev-comp
                            l2 += log(min(error1AThird,MAX0));
                        }else if (reffrag[L-i-1]==2 && frag[L-i-1]==2 && i<=l-1){ // C to C, left single strand part, rev-comp
                            l2 += log(max(1-error1-exterm_s,MIN0));
                        }else if (reffrag[L-i-1]==1 && frag[L-i-1]==1 && i<=l-1){ // G to G, left single strand part, rev-comp
                            l2 += log(max(1-error1,MIN0));
                        }else if (reffrag[L-i-1]==2 && frag[L-i-1]==0 && i>l-1){ // C to T change, left double strand part, rev-comp
                            l2 += log(min(error1AThird+exterm,MAX0));
                        }else if (reffrag[L-i-1]==1 && frag[L-i-1]==3 && i>l-1){ // G to A change, left double strand part, rev-comp
                            l2 += log(min(error1AThird,MAX0));
                        }else if (reffrag[L-i-1]==2 && frag[L-i-1]==2 && i>l-1){ // C to C, left double strand part, rev-comp
                            l2 += log(max(1-error1-exterm,MIN0));
                        }else if (reffrag[L-i-1]==1 && frag[L-i-1]==1 && i>l-1){ // G to G, left double strand part,rev-comp
                            l2 += log(max(1-error1,MIN0));
                        }else if (reffrag[L-i-1]==frag[L-i-1]){ // All other possible no changes, rev-comp
                            l2 += log(max(1-error1,MIN0));
                        }else if (reffrag[L-i-1]!=frag[L-i-1]){ // All other possible changes, rev-comp
                            l2 += log(min(error1AThird,MAX0));
                        }
                    }
                    
                    double error2 = PhredError[seqError[i]];
                    double error2AThird = PhredErrorAThird[seqError[i]];
                    exterm_s = delta_s-4*delta_s*error2AThird;
                    exterm = delta-4*delta*error2AThird;
                    if (reffrag[i] < 4 && frag[i] < 4){
                        if (reffrag[i]==1 && frag[i]==3 && i <= r-1){ // G to A change, right single strand part, rev-comp
                            l2 += log(min(error2AThird+exterm_s,MAX0));
                        }else if(reffrag[i]==2 && frag[i]==0 && i <= r-1){ // C to T change, right single strand part, rev-comp
                            l2 += log(min(error2AThird,MAX0));
                        }else if(reffrag[i]==1 && frag[i]==1 && i <= r-1){ // G to G, right single strand part, rev-comp
                            l2 += log(max(1-error2-exterm_s,MIN0));
                        }else if(reffrag[i]==2 && frag[i]==2 && i <= r-1){ // C to C, right single strand part, rev-comp
                            l2 += log(max(1-error2,MIN0));
                        }else if (reffrag[i]==1 && frag[i]==3 && i > r-1){ // G to A change, right double strand part, rev-comp
                            l2 += log(min(error2AThird+exterm,MAX0));
                        }else if (reffrag[i]==2 && frag[i]==0 && i > r-1){ // C to T change, right double strand part, rev-comp
                            l2 += log(min(error2AThird,MAX0));
                        }else if (reffrag[i]==1 && frag[i]==1 && i > r-1){ // G to G, right double strand part, rev-comp
                            l2 += log(max(1-error2-exterm,MIN0));
                        }else if (reffrag[i]==2 && frag[i]==2 && i > r-1){ // C to C, right double strand part, rev-comp
                            l2 += log(max(1-error2,MIN0));
                        }else if (reffrag[i]==frag[i]){ // All other possible no changes, rev-comp, rev-comp
                            l2 += log(max(1-error2,MIN0));
                        }else if (reffrag[i]!=frag[i]){ // All other possible changes, rev-comp
                            l2 += log(min(error2AThird,MAX0));
                        }
                    }
                }
                //cout<<"l_check "<<l_check<<", l "<<l<<", r "<<r<<", "<< pnick<<"\n";
                l_pmd += exp(l2)*dp*(1-pnick);
            }
            p += dp;
        }
        
        for (int r=l_check+1; r<=L-2-l; r++){ // To guarantee l+r < = L-2
            double pr = 0.5*lambda*pow(1-lambda,r)/(1-pow(1-lambda,L-1));
            if (r == 0){
                pr += 0.5;
            }
            double dp = 4*pl*pr/(3+(1-pow(1-lambda,L-1)-(L-1)*lambda*pow(1-lambda,L-1))/(pow(1-pow(1-lambda,L-1),2)));
            if (dp >= Tol){
                double pnick = 0;
                //--------- - - - -                        ------------------//
                //            /|\                                            //
                //             |___ Nick                                     //
                for (int npos=l; npos<=(l_check-1 < L-r-1 ? l_check-1 : L-r-1); npos++){ //Nick is within left l_check region ???
                    double l2 = 0;
                    for (int i=0; i<l_check; i++){
                        double error1 = PhredError[seqError[L-i-1]];
                        double error1AThird = PhredErrorAThird[seqError[L-i-1]];
                        exterm_s = delta_s-4*delta_s*error1AThird;
                        exterm = delta-4*delta*error1AThird;
                        if (reffrag[L-i-1] < 4 && frag[L-i-1] < 4){
                            if (reffrag[L-i-1]==2 && frag[L-i-1]==0 && i<=l-1){ // C to T change, left single strand part, rev-comp
                                l2 += log(min(error1AThird+exterm_s,MAX0));
                            }else if (reffrag[L-i-1]==1 && frag[L-i-1]==3 && i<=l-1){ // G to A change, left single strand part, rev-comp
                                l2 += log(min(error1AThird,MAX0));
                            }else if (reffrag[L-i-1]==2 && frag[L-i-1]==2 && i<=l-1){ // C to C, left single strand part, rev-comp
                                l2 += log(max(1-error1-exterm_s,MIN0));
                            }else if (reffrag[L-i-1]==1 && frag[L-i-1]==1 && i<=l-1){ // G to G, left single strand part, rev-comp
                                l2 += log(max(1-error1,MIN0));
                            }else if (reffrag[L-i-1]==2 && frag[L-i-1]==0 && i<=npos){ // C to T change, left double strand part before nick, rev-comp
                                l2 += log(min(error1AThird+exterm,MAX0));
                            }else if (reffrag[L-i-1]==1 && frag[L-i-1]==3 && i<=npos){ // G to A change, left double strand part before nick, rev-comp
                                l2 += log(min(error1AThird,MAX0));
                            }else if (reffrag[L-i-1]==2 && frag[L-i-1]==2 && i<=npos){ // C to C, left double strand part before nick, rev-comp
                                l2 += log(max(1-error1-exterm,MIN0));
                            }else if (reffrag[L-i-1]==1 && frag[L-i-1]==1 && i<=npos){ // G to G, left double strand part before nick, rev-comp
                                l2 += log(max(1-error1,MIN0));
                            }else if (reffrag[L-i-1]==2 && frag[L-i-1]==0 && i>npos){ // C to T change, left double strand part after nick, rev-comp
                                l2 += log(min(error1AThird,MAX0));
                            }else if (reffrag[L-i-1]==1 && frag[L-i-1]==3 && i>npos){ // G to A change, left double strand part after nick, rev-comp
                                l2 += log(min(error1AThird+exterm,MAX0));
                            }else if (reffrag[L-i-1]==2 && frag[L-i-1]==2 && i>npos){ // C to C, left double strand part after nick, rev-comp
                                l2 += log(max(1-error1,MIN0));
                            }else if (reffrag[L-i-1]==1 && frag[L-i-1]==1 && i>npos){ // G to G, left double strand part after nick, rev-comp
                                l2 += log(max(1-error1-exterm,MIN0));
                            }else if (reffrag[L-i-1]==frag[L-i-1]){ // All other possible no changes, rev-comp
                                l2 += log(max(1-error1,MIN0));
                            }else if (reffrag[L-i-1]!=frag[L-i-1]){ // All other possible changes, rev-comp
                                l2 += log(min(error1AThird,MAX0));
                            }
                        }
                        
                        double error2 = PhredError[seqError[i]];
                        double error2AThird = PhredErrorAThird[seqError[i]];
                        exterm_s = delta_s-4*delta_s*error2AThird;
                        if (reffrag[i] < 4 && frag[i] < 4){
                            if (reffrag[i]==1 && frag[i]==3){ // G to A change, right single strand part, rev-comp
                                l2 += log(min(error2AThird+exterm_s,MAX0));
                            }else if(reffrag[i]==2 && frag[i]==0){ // C to T change, right single strand part, rev-comp
                                l2 += log(min(error2AThird,MAX0));
                            }else if(reffrag[i]==1 && frag[i]==1){ // G to G, right single strand part, rev-comp
                                l2 += log(max(1-error2-exterm_s,MIN0));
                            }else if(reffrag[i]==2 && frag[i]==2){ // C to C, right single strand part, rev-comp
                                l2 += log(max(1-error2,MIN0));
                            }else if (reffrag[i]==frag[i]){ // All other possible no changes, rev-comp
                                l2 += log(max(1-error2,MIN0));
                            }else if (reffrag[i]!=frag[i]){ // All other changes, rev-comp
                                l2 += log(min(error2AThird,MAX0));
                            }
                        }
                    }
                    double dpnick = npos < L-r-1 ? nv/(1+(L-l-r-2)*nv) : (1-nv)/(1+(L-l-r-2)*nv);
                    l_pmd += exp(l2)*dp*dpnick;
                    pnick += dpnick;
                }
                //nick occurring outside both l_check region
                //--------- - - - -                        ------------------//
                //                           /|\                             //
                //                            |___ Nick                      //
                double l2 = 0;
                for (int i=0; i<l_check; i++){
                    double error1 = PhredError[seqError[L-i-1]];
                    double error1AThird = PhredErrorAThird[seqError[L-i-1]];
                    exterm_s = delta_s-4*delta_s*error1AThird;
                    exterm = delta-4*delta*error1AThird;
                    if (reffrag[L-i-1] < 4 && frag[L-i-1] < 4){
                        if (reffrag[L-i-1]==2 && frag[L-i-1]==0 && i<=l-1){ // C to T change, left single strand part, rev-comp
                            l2 += log(min(error1AThird+exterm_s,MAX0));
                        }else if (reffrag[L-i-1]==1 && frag[L-i-1]==3 && i<=l-1){ // G to A change, left single strand part, rev-comp
                            l2 += log(min(error1AThird,MAX0));
                        }else if (reffrag[L-i-1]==2 && frag[L-i-1]==2 && i<=l-1){ // C to C, left single strand part, rev-comp
                            l2 += log(max(1-error1-exterm_s,MIN0));
                        }else if (reffrag[L-i-1]==1 && frag[L-i-1]==1 && i<=l-1){ // G to G, left single strand part, rev-comp
                            l2 += log(max(1-error1,MIN0));
                        }else if (reffrag[L-i-1]==2 && frag[L-i-1]==0 && i>l-1){ // C to T change, left double strand part, rev-comp
                            l2 += log(min(error1AThird+exterm,MAX0));
                        }else if (reffrag[L-i-1]==1 && frag[L-i-1]==3 && i>l-1){ // G to A change, left double strand part, rev-comp
                            l2 += log(min(error1AThird,MAX0));
                        }else if (reffrag[L-i-1]==2 && frag[L-i-1]==2 && i>l-1){ // C to C, left double strand part, rev-comp
                            l2 += log(max(1-error1-exterm,MIN0));
                        }else if (reffrag[L-i-1]==1 && frag[L-i-1]==1 && i>l-1){ // G to G, left double strand part, rev-comp
                            l2 += log(max(1-error1,MIN0));
                        }else if (reffrag[L-i-1]==frag[L-i-1]){ // All other possible no changes, rev-comp
                            l2 += log(max(1-error1,MIN0));
                        }else if (reffrag[L-i-1]!=frag[L-i-1]){ // All other possible changes, rev-comp
                            l2 += log(min(error1AThird,MAX0));
                        }
                    }
                    
                    double error2 = PhredError[seqError[i]];
                    double error2AThird = PhredErrorAThird[seqError[i]];
                    exterm_s = delta_s-4*delta_s*error2AThird;
                    if (reffrag[i] < 4 && frag[i] < 4){
                        if (reffrag[i]==1 && frag[i]==3){ // G to A change, right single strand part, rev-comp
                            l2 += log(min(error2AThird+exterm_s,MAX0));
                        }else if(reffrag[i]==2 && frag[i]==0){ // C to T change, right single strand part, rev-comp
                            l2 += log(min(error2AThird,MAX0));
                        }else if(reffrag[i]==1 && frag[i]==1){ // G to G, right single strand part, rev-comp
                            l2 += log(max(1-error2-exterm_s,MIN0));
                        }else if(reffrag[i]==2 && frag[i]==2){ // C to C, right single strand part, rev-comp
                            l2 += log(max(1-error2,MIN0));
                        }else if (reffrag[i]==frag[i]){ // All other possible no changes, rev-comp
                            l2 += log(max(1-error2,MIN0));
                        }else if (reffrag[i]!=frag[i]){ // All other possible changes, rev-comp
                            l2 += log(min(error2AThird,MAX0));
                        }
                    }
                }
                l_pmd += exp(l2)*dp*(1-pnick);
            }
            p += dp;
        }
    }
    //-----------------                        - - - - ----------//
    //                                            /|\            //
    //                                             |___ Nick     //
    for (int r=0; r<=l_check; r++){
        double pr = 0.5*lambda*pow(1-lambda,r)/(1-pow(1-lambda,L-1));
        if (r == 0){
            pr += 0.5;
        }
        for (int l=l_check+1;l<=L-2-r; l++){
            double pl = 0.5*lambda*pow(1-lambda,l)/(1-pow(1-lambda,L-1));
            if (l==0){
                pl += 0.5;
            }
            // Joint probability of (l,r)
            double dp = 4*pl*pr/(3+(1-pow(1-lambda,L-1)-(L-1)*lambda*pow(1-lambda,L-1))/(pow(1-pow(1-lambda,L-1),2)));
            if (dp > Tol){
                double pnick = 0;
                for (int npos = (L-l_check-1 > l ? L-l_check-1 : l); npos <= L-r-1; npos++){ //Nick is within right l_check region // nuc_llik added Please be noted that npos has the potential to go to left overhang region
                    double l2 = 0;
                    for (int i=0; i<l_check; i++){
                        double error2 = PhredError[seqError[i]];
                        double error2AThird = PhredErrorAThird[seqError[i]];
                        exterm_s = delta_s-4*delta_s*error2AThird;
                        exterm = delta-4*delta*error2AThird;
                        if (reffrag[i] < 4 && frag[i] < 4){
                            if (reffrag[i]==1 && frag[i]==3 && i <= r-1){ // G to A change, right single strand part, rev-comp
                                l2 += log(min(error2AThird+exterm_s,MAX0));
                            }else if(reffrag[i]==2 && frag[i]==0 && i <= r-1){ // C to T change, right single strand part, rev-comp
                                l2 += log(min(error2AThird,MAX0));
                            }else if(reffrag[i]==1 && frag[i]==1 && i <= r-1){ // G to G, right single strand part, rev-comp
                                l2 += log(max(1-error2-exterm_s,MIN0));
                            }else if(reffrag[i]==2 && frag[i]==2 && i <= r-1){ // C to C, right single strand part, rev-comp
                                l2 += log(max(1-error2,MIN0));
                            }else if (reffrag[i]==1 && frag[i]==3 && i <= L-npos-1){ // G to A change, right double strand part, right of the nick, rev-comp
                                l2 += log(min(error2AThird+exterm,MAX0));
                            }else if (reffrag[i]==2 && frag[i]==0 && i <= L-npos-1){ // C to T change, right double strand part, right of the nick, rev-comp
                                l2 += log(min(error2AThird,MAX0));
                            }else if (reffrag[i]==1 && frag[i]==1 && i <= L-npos-1){ // G to G, right double strand part, right of the nick, rev-comp
                                l2 += log(max(1-error2-exterm,MIN0));
                            }else if (reffrag[i]==2 && frag[i]==2 && i <= L-npos-1){ // C to C, right double strand part, right of the nick, rev-comp
                                l2 += log(max(1-error2,MIN0));
                            }else if (reffrag[i]==1 && frag[i]==3 && i > L-npos-1){ // G to A, right double strand part, left of the nick, rev-comp
                                l2 += log(min(error2AThird,MAX0));
                            }else if (reffrag[i]==2 && frag[i]==0 && i > L-npos-1){ // C to T, right double strand part, left of the nick, rev-comp
                                l2 += log(min(error2AThird+exterm,MAX0));
                            }else if (reffrag[i]==1 && frag[i]==1 && i > L-npos-1){ // G to G, right double strand part, left of the nick, rev-comp
                                l2 += log(max(1-error2,MIN0));
                            }else if (reffrag[i]==2 && frag[i]==2 && i > L-npos-1){ // C to C, right double strand part, left of the nick, rev-comp
                                l2 += log(max(1-error2-exterm,MIN0));
                            }else if (reffrag[i]==frag[i]){ // All other possible no changes, rev-comp
                                l2 += log(max(1-error2,MIN0));
                            }else if (reffrag[i]!=frag[i]){ // All other possible changes, rev-comp
                                l2 += log(min(error2AThird,MAX0));
                            }
                        }
                        
                        double error1 = PhredError[seqError[L-i-1]];
                        double error1AThird = PhredErrorAThird[seqError[L-i-1]];
                        exterm_s = delta_s-4*delta_s*error1AThird;
                        if (reffrag[L-i-1] < 4 && frag[L-i-1] < 4){
                            if (reffrag[L-i-1]==2 && frag[L-i-1]==0){ // C to T change, left single strand part, rev-comp
                                l2 += log(min(error1AThird+exterm_s,MAX0));
                            }else if (reffrag[L-i-1]==1 && frag[L-i-1]==3){ // G to A change, left single strand part, rev-comp
                                l2 += log(min(error1AThird,MAX0));
                            }else if (reffrag[L-i-1]==2 && frag[L-i-1]==2){ // C to C, left single strand part, rev-comp
                                l2 += log(max(1-error1-exterm_s,MIN0));
                            }else if (reffrag[L-i-1]==1 && frag[L-i-1]==1){ // G to G, left single strand part, rev-comp
                                l2 += log(max(1-error1,MIN0));
                            }else if (reffrag[L-i-1]==frag[L-i-1]){ // All other possible no changes, rev-comp
                                l2 += log(max(1-error1,MIN0));
                            }else if (reffrag[L-i-1]!=frag[L-i-1]){ // All other possible changes, rev-comp
                                l2 += log(min(error1AThird,MAX0));
                            }
                        }
                    }
                    double dpnick = npos < L-r-1 ? nv/(1+(L-l-r-2)*nv) : (1-nv)/(1+(L-l-r-2)*nv);
                    l_pmd += exp(l2)*dp*dpnick;
                    //cout << "Length " << L <<", l "<<l<<", r "<<r<<", nick position " << npos << ", prob " << dpnick << "\n";
                    pnick += dpnick;
                }
                //cout<<
                //-----------------                        - - - - ----------//
                //                             /|\                           //
                //                              |___ Nick                    //
                // nick occurring outside both l_check region
                double l2 = 0;
                for (int i=0; i<l_check; i++){
                    double error1 = PhredError[seqError[L-i-1]];
                    double error1AThird = PhredErrorAThird[seqError[L-i-1]];
                    exterm_s = delta_s-4*delta_s*error1AThird;
                    if (reffrag[L-i-1] < 4 && frag[L-i-1] < 4){
                        if (reffrag[L-i-1]==2 && frag[L-i-1]==0){ // C to T change, left single strand part, rev-comp
                            l2 += log(min(error1AThird+exterm_s,MAX0));
                        }else if (reffrag[L-i-1]==1 && frag[L-i-1]==3){ // G to A change, left single strand part, rev-comp
                            l2 += log(min(error1AThird,MAX0));
                        }else if (reffrag[L-i-1]==2 && frag[L-i-1]==2){ // C to C, left single strand part, rev-comp
                            l2 += log(max(1-error1-exterm_s,MIN0));
                        }else if (reffrag[L-i-1]==1 && frag[L-i-1]==1){ // G to G, left single strand part, rev-comp
                            l2 += log(max(1-error1,MIN0));
                        }else if (reffrag[L-i-1]==frag[L-i-1]){ // All other possible no changes, rev-comp
                            l2 += log(max(1-error1,MIN0));
                        }else if (reffrag[L-i-1]!=frag[L-i-1]){ // All other possible changes, rev-comp
                            l2 += log(min(error1AThird,MAX0));
                        }
                    }
                    
                    double error2 = PhredError[seqError[i]];
                    double error2AThird = PhredErrorAThird[seqError[i]];
                    exterm_s = delta_s-4*delta_s*error2AThird;
                    exterm = delta_s-4*delta_s*error2AThird;
                    if (reffrag[i] < 4 && frag[i] < 4){
                        if (reffrag[i]==1 && frag[i]==3 && i <= r-1){ // G to A change, right single strand part, rev-comp
                            l2 += log(min(error2AThird+exterm_s,MAX0));
                        }else if(reffrag[i]==2 && frag[i]==0 && i <= r-1){ // C to T change, right single strand part, rev-comp
                            l2 += log(min(error2AThird,MAX0));
                        }else if(reffrag[i]==1 && frag[i]==1 && i <= r-1){ // G to G, right single strand part, rev-comp
                            l2 += log(max(1-error2-exterm_s,MIN0));
                        }else if(reffrag[i]==2 && frag[i]==2 && i <= r-1){ // C to C, right single strand part, rev-comp
                            l2 += log(max(1-error2,MIN0));
                        }else if (reffrag[i]==1 && frag[i]==3 && i > r-1){ // G to A change, right double strand part, rev-comp
                            l2 += log(min(error2AThird+exterm,MAX0));
                        }else if (reffrag[i]==2 && frag[i]==0 && i > r-1){ // C to T change, right double strand part, rev-comp
                            l2 += log(min(error2AThird,MAX0));
                        }else if (reffrag[i]==1 && frag[i]==1 && i > r-1){ // G to G, right double strand part, rev-comp
                            l2 += log(max(1-error2-exterm,MIN0));
                        }else if (reffrag[i]==2 && frag[i]==2 && i > r-1){ // C to C, right double strand part, rev-comp
                            l2 += log(max(1-error2,MIN0));
                        }else if (reffrag[i]==frag[i]){ // All other possible no changes, rev-comp
                            l2 += log(max(1-error2,MIN0));
                        }else if (reffrag[i]!=frag[i]){ // All other possible changes, rev-comp
                            l2 += log(min(error2AThird,MAX0));
                        }
                    }
                }
                //cout<<"l_check "<<l_check<<", l "<<l<<", r "<<r<<", "<< pnick<<"\n";
                //cout<<l2<<" l2\n";
		l_pmd += exp(l2)*dp*(1-pnick);
            }
            p += dp;
        }
    }
    //cout << "nuc_llik_nuc_llik_" << l_pmd << "\n";
    //-----------------                        ------------------//
    //                             /|\                           //
    //                              |___ Nick                    //
    //Investigate each possible left and right 5' overhang pair with (l,r) > l_check (15)
    if (1-p > Tol){
        double l2 = 0;
        for (int i=0; i<l_check; i++){
            double error1 = PhredError[seqError[L-i-1]];
            double error1AThird = PhredErrorAThird[seqError[L-i-1]];
            exterm_s = delta_s-4*delta_s*error1AThird;
            if (reffrag[L-i-1] < 4 && frag[L-i-1] < 4){
                if (reffrag[L-i-1]==2 && frag[L-i-1]==0){ // C to T change, left single strand part, rev-comp
                    l2 += log(min(error1AThird+exterm_s,MAX0));
                }else if (reffrag[L-i-1]==1 && frag[L-i-1]==3){ // G to A change, left single strand part, rev-comp
                    l2 += log(min(error1AThird,MAX0));
                }else if (reffrag[L-i-1]==2 && frag[L-i-1]==2){ // C to C, left single strand part, rev-comp
                    l2 += log(max(1-error1-exterm_s,MIN0));
                }else if (reffrag[L-i-1]==1 && frag[L-i-1]==1){ // G to G, left single strand part, rev-comp
                    l2 += log(max(1-error1,MIN0));
                }else if (reffrag[L-i-1]==frag[L-i-1]){ // All other possible no changes, rev-comp
                    l2 += log(max(1-error1,MIN0));
                }else if (reffrag[L-i-1]!=frag[L-i-1]){ // All other possible changes, rev-comp
                    l2 += log(min(error1AThird,MAX0));
                }
            }
            
            double error2 = PhredError[seqError[i]];
            double error2AThird = PhredErrorAThird[seqError[i]];
            exterm_s = delta_s-4*delta_s*error2AThird;
            if (reffrag[i] < 4 && frag[i] < 4){
                if (reffrag[i]==1 && frag[i]==3){ // G to A change, right single strand part, rev-comp
                    l2 += log(min(error2AThird+exterm_s,MAX0));
                }else if(reffrag[i]==2 && frag[i]==0){ // C to T change, right single strand part, rev-comp
                    l2 += log(min(error2AThird,MAX0));
                }else if(reffrag[i]==1 && frag[i]==1){ // G to G, right single strand part, rev-comp
                    l2 += log(max(1-error2-exterm_s,MIN0));
                }else if(reffrag[i]==2 && frag[i]==2){ // C to C change, right single strand part, rev-comp
                    l2 += log(max(1-error2,MIN0));
                }else if (reffrag[i]==frag[i]){ // All other possible no changes, rev-comp
                    l2 += log(max(1-error2,MIN0));
                }else if (reffrag[i]!=frag[i]){ // All other possible changes, rev-comp
                    l2 += log(min(error2AThird,MAX0));
                }
            }
        }
        l_pmd += exp(l2)*(1-p);
    }
    //cout<<l_pmd<<"nuc_lliktest\n";
    return l_pmd;
}
