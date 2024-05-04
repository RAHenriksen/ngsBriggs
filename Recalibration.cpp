#include <cstring>
#include <cstdlib>
#include <cstdio>
#include <cassert>
#include <sstream>
#include <string>

#include <htslib/hts.h>
#include <htslib/sam.h>
#include <htslib/kstring.h>
#include <htslib/faidx.h>
#include <zlib.h>
#include <cmath>
#include <iomanip>
#include <ctime>
#include <getopt.h>
#include <iostream>
#include <array>

#include "profile.h"
#include "bfgs.h"
#include "htslib/bgzf.h"

#include "misc.h"
#include "Likelihood.h"
#include "ngsBriggs.h"

extern tsk_struct *my_tsk_struct;

//The log-likelihood for recalibration the ancient probs
double tsk_loglike_recalibration(const double *x, std::vector<bam1_t *> *reads,int from,int to,sam_hdr_t *hdr, faidx_t *seq_ref,int len_limit, int len_min, char * model, double eps, double lambda, double delta, double delta_s, double nv,int threadid,double Tol){
    //  fprintf(stderr,"(%f,%f,%f,%f) seq_ref: %p\n",x[0],x[1],x[2],x[3],seq_ref);
    double anc_mu = x[0];
    double anc_si = x[1];
    double mod_mu = x[2];
    double mod_si = x[3];
    double ll = 0;
    char reconstructedRef[512];
    char myread[512];
    char myrefe[512];
    char yourread[512];
    char yourrefe[512];
    uint8_t  yourqual[512];
    std::pair< kstring_t*, std::vector<int> >  mypair;
    kstring_t *kstr =new kstring_t;
    kstr->l=kstr->m=0;
    kstr->s=NULL;
    mypair.first = kstr;
    
    char refeBase, readBase;
    
    
    unsigned char * indref = NULL;
    double x_max1 = ((double)len_limit-1+0.5-anc_mu)/anc_si;
    double x_min1 = ((double)len_min-0.5-anc_mu)/anc_si;
    double x_max2 = ((double)len_limit-1+0.5-mod_mu)/mod_si;
    double x_min2 = ((double)len_min-0.5-mod_mu)/mod_si;
    if(from==-1||to==-1){
        from =0;
        to = reads->size();
    }
    //    fprintf(stderr,"from: %d to:%d\n",from,to);
    for(int i=from;i<to;i++) {
        //   fprintf(stderr,"ll\t[%d]\t[%d]\t%f\tin\n",threadid,i,ll);
        bam1_t *b = (*reads)[i];
        //then we simply write it to the output
        memset(reconstructedRef,0,512);
        memset(myread,'N',512);
        memset(myrefe,'N',512);
        memset(yourread,'N',512);
        memset(yourrefe,'N',512);
        
        if (seq_ref != NULL){
            wrapperwithref(b,hdr,myread,myrefe,seq_ref);
        }else{
            reconstructRefWithPosHTS(b,mypair,reconstructedRef);
            wrapper(b,mypair.first->s,mypair.second,0,0,NULL,NULL,1,myread,myrefe);
        }
        
        if (len_limit <=0){
            len_limit = 512;
        };
        //cout << b->core.l_qseq << "Test nuc_llik\n";
        //	fprintf(stderr,"ll\t[%d]\t[%d]\t%f\tin2\t %d %d\n",threadid,i,ll,b->core.l_qseq,len_limit);
        if (b->core.l_qseq>=30 && b->core.l_qseq<len_limit){
            for (int cycle=0;cycle<b->core.l_qseq;cycle++){
                yourqual[cycle] = 254;
            }
            
            for (int cycle=0;cycle<b->core.l_qseq;cycle++){
                refeBase = refToChar[myrefe[cycle]];
                readBase = refToChar[myread[cycle]];
                size_t pos = b->core.pos+cycle;
                
                //if(refeBase!=4 && readBase!=4) {
                int dist5p=cycle;
                int dist3p=b->core.l_qseq-1-cycle;
                if( bam_is_rev(b) ){
                    refeBase=com[refeBase];
                    readBase=com[readBase];
                    dist5p=int(b->core.l_qseq)-1-cycle;
                    dist3p=cycle;
                }
                yourread[dist5p] = readBase;
                yourrefe[dist5p] = refeBase;
                yourqual[dist5p] = bam_get_qual(b)[cycle];
                //}
            }
            
            double l_anc, l_err;
            int L = b->core.l_qseq;
            l_err = ErrorLik(yourrefe, yourread, L, yourqual);
            if (!strcasecmp("b",model)){
	      l_anc = PMDLik_b(yourrefe, yourread, L, lambda, delta, delta_s, nv, yourqual,Tol); // Ancient Likelihood based on biotin model
            }else if(!strcasecmp("nb",model)){
	      l_anc = 0.5*PMDLik_b(yourrefe, yourread, L, lambda, delta, delta_s, nv, yourqual,Tol)+0.5*PMDLik_nb(yourrefe, yourread, L, lambda, delta, delta_s, nv, yourqual,Tol); // Ancient Likelihood based on non-biotin model
	    }else{
                fprintf(stderr,"Please specify a deamination model for further calculations.\n");
                return -1;
            }
            double y_max1 = (std::min((double)len_limit-1,std::max((double)len_min,(double)L))+0.5-anc_mu)/anc_si;
            double y_min1 = (std::max((double)len_min,std::min((double)len_limit-1,(double)L))-0.5-anc_mu)/anc_si;
            double y_max2 = (std::min((double)len_limit-1,std::max((double)len_min,(double)L))+0.5-mod_mu)/mod_si;
            double y_min2 = (std::max((double)len_min,std::min((double)len_limit-1,(double)L))-0.5-mod_mu)/mod_si;
	    ll += log(l_anc*NormalINC(y_max1, y_min1, x_max1, x_min1)*(1-eps)+l_err*NormalINC(y_max2, y_min2, x_max2, x_min2)*eps);
	}
        
    }
    return -ll;
}

double tsk_all_loglike_recalibration(const double *x, const void *dats){
    //  fprintf(stderr,"[%s]\n",__FUNCTION__);
    ncalls++;
    tsk_struct *ts = (tsk_struct *) dats;
    double check = tsk_loglike_recalibration(x, ts->reads,ts->from,ts->to,ts->hdr, ts->seq_ref, ts->len_limit, ts->len_min, ts->model, ts->eps, ts->lambda, ts->delta, ts->delta_s, ts->nv,ts->threadid,ts->Tol);
    //  fprintf(stderr,"[%s]: total like: %f\n",__FUNCTION__,check);
    return check;
}

void *tsk_all_loglike_recalibration_slave(void *dats){
    //fprintf(stderr,"[%s]\n",__FUNCTION__);
    ncalls++;
    tsk_struct *ts = (tsk_struct *) &(my_tsk_struct[(size_t) dats]);
    ts->llh_result = tsk_loglike_recalibration(ts->x, ts->reads,ts->from,ts->to,ts->hdr, ts->seq_ref, ts->len_limit, ts->len_min, ts->model, ts->eps, ts->lambda, ts->delta, ts->delta_s, ts->nv,ts->threadid,ts->Tol);
    
    pthread_exit(NULL);
}

void tsk_loglike_recalibration_grad(const double *x, double *y, std::vector<bam1_t*> *reads,int from, int to,sam_hdr_t *hdr, faidx_t *seq_ref,int len_limit, int len_min, char * model, double eps, double lambda, double delta, double delta_s, double nv,double Tol){
    //fprintf(stderr,"[%s]\n",__FUNCTION__);
    htsFormat *dingding7 =(htsFormat*) calloc(1,sizeof(htsFormat));
    double anc_mu = x[0];
    double anc_si = x[1];
    double mod_mu = x[2];
    double mod_si = x[3];
    y[0] = 0;
    y[1] = 0;
    y[2] = 0;
    y[3] = 0;
    double ll = 0;
    char reconstructedRef[512];
    char myread[512];
    char myrefe[512];
    char yourread[512];
    char yourrefe[512];
    uint8_t  yourqual[512];
    std::pair< kstring_t*, std::vector<int> >  mypair;
    kstring_t *kstr =new kstring_t;
    kstr->l=kstr->m=0;
    kstr->s=NULL;
    mypair.first = kstr;
    char refeBase, readBase;
    unsigned char * indref = NULL;
    double x_max1 = ((double)len_limit-1+0.5-anc_mu)/anc_si;
    double x_min1 = ((double)len_min-0.5-anc_mu)/anc_si;
    double x_max2 = ((double)len_limit-1+0.5-mod_mu)/mod_si;
    double x_min2 = ((double)len_min-0.5-mod_mu)/mod_si;
    if(1||from==-1||to==-1){//fix this uncomment out.
        from = 0;
        to = reads->size();
    }
    //    fprintf(stderr,"from: %d to:%d\n",from,to):qw
    for(int i=from;i<to;i++){
        bam1_t *b = (*reads)[i];
         
        memset(reconstructedRef,0,512);
        memset(myread,'N',512);
        memset(myrefe,'N',512);
        
        if (seq_ref != NULL){
            wrapperwithref(b,hdr,myread,myrefe,seq_ref);
        }else{
            reconstructRefWithPosHTS(b,mypair,reconstructedRef);
            wrapper(b,mypair.first->s,mypair.second,0,0,NULL,NULL,1,myread,myrefe);
        }
        
        if (len_limit <=0){
            len_limit = 512;
        };
        
        if (b->core.l_qseq>=30 && b->core.l_qseq<len_limit) {
            
            for (int cycle=0;cycle<b->core.l_qseq;cycle++)
                yourqual[cycle] = 254;
            
            for (int cycle=0;cycle<b->core.l_qseq;cycle++){
                refeBase = refToChar[myrefe[cycle]];
                readBase = refToChar[myread[cycle]];
                
                //                if(refeBase!=4 && readBase!=4){
                int dist5p=cycle;
                int dist3p=b->core.l_qseq-1-cycle;
                if( bam_is_rev(b) ){
                    refeBase=com[refeBase];
                    readBase=com[readBase];
                    
                    dist5p=int(b->core.l_qseq)-1-cycle;
                    dist3p=cycle;
                }
                yourread[dist5p] = readBase;
                yourrefe[dist5p] = refeBase;
                yourqual[dist5p] = bam_get_qual(b)[cycle];

            }
            
            double l_anc, l_err;
            int L = b->core.l_qseq;
            l_err = ErrorLik(yourrefe, yourread, L, yourqual);
            if (!strcasecmp("b",model)){
	      l_anc = PMDLik_b(yourrefe, yourread, L, lambda, delta, delta_s, nv, yourqual,Tol); // Ancient Likelihood based on biotin model
            }else if(!strcasecmp("nb",model)){
	      l_anc = 0.5*PMDLik_b(yourrefe, yourread, L, lambda, delta, delta_s, nv, yourqual,Tol)+0.5*PMDLik_nb(yourrefe, yourread, L, lambda, delta, delta_s, nv, yourqual,Tol); // Ancient Likelihood based on non-biotin model
            }else{
                fprintf(stderr,"Please specify a deamination model for further calculations.\n");
                exit(0);
            }
            double y_max1 = (std::min((double)len_limit-1,std::max((double)len_min,(double)L))+0.5-anc_mu)/anc_si;
            double y_min1 = (std::max((double)len_min,std::min((double)len_limit-1,(double)L))-0.5-anc_mu)/anc_si;
            double y_max2 = (std::min((double)len_limit-1,std::max((double)len_min,(double)L))+0.5-mod_mu)/mod_si;
            double y_min2 = (std::max((double)len_min,std::min((double)len_limit-1,(double)L))-0.5-mod_mu)/mod_si;
            
            //cout<<"nuc_lliktest y_max1 "<<y_max1<<" y_min1 "<<y_min1<<" "<<NormalINC_grad_mu(y_max1, y_min1, x_max1, x_min1, anc_mu, anc_si)<<"\n";
            //ll += log(l_anc*normalINC(y_max1, y_min1, x_max1, x_min1)*(1-eps)+l_err*normalINC(y_max2, y_min2, x_max2, x_min2)*eps)
	    y[0] -= l_anc*(1-eps)*NormalINC_grad_mu(y_max1, y_min1, x_max1, x_min1, anc_mu, anc_si)/(l_anc*NormalINC(y_max1, y_min1, x_max1, x_min1)*(1-eps)+l_err*NormalINC(y_max2, y_min2, x_max2, x_min2)*eps);
            y[1] -= l_anc*(1-eps)*NormalINC_grad_si(y_max1, y_min1, x_max1, x_min1, anc_mu, anc_si)/(l_anc*NormalINC(y_max1, y_min1, x_max1, x_min1)*(1-eps)+l_err*NormalINC(y_max2, y_min2, x_max2, x_min2)*eps);
            y[2] -= l_err*eps*NormalINC_grad_mu(y_max2, y_min2, x_max2, x_min2, mod_mu, mod_si)/(l_anc*NormalINC(y_max1, y_min1, x_max1, x_min1)*(1-eps)+l_err*NormalINC(y_max2, y_min2, x_max2, x_min2)*eps);
            y[3] -= l_err*eps*NormalINC_grad_si(y_max2, y_min2, x_max2, x_min2, mod_mu, mod_si)/(l_anc*NormalINC(y_max1, y_min1, x_max1, x_min1)*(1-eps)+l_err*NormalINC(y_max2, y_min2, x_max2, x_min2)*eps);
            //double PostAncProb = AncProb(yourrefe, yourread, b->core.l_qseq, lambda, delta, delta_s, nv, yourqual, model ,eps);
            //double PostPMDProb = PMDProb(yourrefe, yourread, b->core.l_qseq, lambda, delta, delta_s, nv, yourqual, model);
            //cout<<"\t"<<"AN: "<<PostAncProb<<" \t"<<"PD: "<<PostPMDProb<<"\n";
            
            
        }
        //    fprintf(stderr,"\nmyread:\n%.*s\nmyReference:\n%.*s\n",b->core.l_qseq,myread,b->core.l_qseq,myrefe);
        //    fprintf(stderr,"---read[%d]----\n",nproc1);
    }
}

void tsk_all_loglike_recalibration_grad(const double *x, double *y,const void *dats){
    ncalls_grad++;
    //tsk_struct *ts = &my_tsk_struct[0];//assuming single thread
    tsk_struct *ts = (tsk_struct *) dats;//assuming multiple threads
    tsk_loglike_recalibration_grad(x,y,ts->reads,ts->from,ts->to,ts->hdr,ts->seq_ref,ts->len_limit, ts->len_min, ts->model,ts->eps,ts->lambda,ts->delta,ts->delta_s,ts->nv,ts->Tol);
    //fprintf(stderr,"[%s]\n",__FUNCTION__);
}

// Check here
void *tsk_all_loglike_recalibration_grad_slave(void *dats){
    //fprintf(stderr,"[%s]\n",__FUNCTION__);
    ncalls_grad++;
    tsk_struct *ts = (tsk_struct *) &(my_tsk_struct[(size_t) dats]);
    //ts->llh_grad_result;
    tsk_loglike_recalibration_grad(ts->x,ts->llh_result_grad,ts->reads,ts->from,ts->to,ts->hdr,ts->seq_ref,ts->len_limit, ts->len_min,ts->model,ts->eps,ts->lambda,ts->delta,ts->delta_s,ts->nv,ts->Tol);
    //ts->llh_result = tsk_loglike_recalibration(ts->x, ts->reads,ts->from,ts->to,ts->hdr, ts->seq_ref, ts->len_limit, ts->model, ts->eps, ts->lambda, ts->delta, ts->delta_s, ts->nv,ts->threadid);
    pthread_exit(NULL);
}

void tsk_loglike_recalibration_hess(const double *x, double **y, std::vector<bam1_t*> *reads,int from, int to,sam_hdr_t *hdr, faidx_t *seq_ref,int len_limit, int len_min, char * model, double eps, double lambda, double delta, double delta_s, double nv,double Tol){
    //fprintf(stderr,"[%s]\n",__FUNCTION__);
    htsFormat *dingding7 =(htsFormat*) calloc(1,sizeof(htsFormat));
    double anc_mu = x[0];
    double anc_si = x[1];
    double mod_mu = x[2];
    double mod_si = x[3];
    for(int i=0; i<3; i++){
        for(int j=0; j<3; j++){
            y[i][j] = 0;
        }
    }
    
    char reconstructedRef[512];
    char myread[512];
    char myrefe[512];
    char yourread[512];
    char yourrefe[512];
    uint8_t  yourqual[512];
    std::pair< kstring_t*, std::vector<int> >  mypair;
    kstring_t *kstr =new kstring_t;
    kstr->l=kstr->m=0;
    kstr->s=NULL;
    mypair.first = kstr;
    char refeBase, readBase;
    unsigned char * indref = NULL;
    double x_max1 = ((double)len_limit-1+0.5-anc_mu)/anc_si;
    double x_min1 = ((double)len_min-0.5-anc_mu)/anc_si;
    double x_max2 = ((double)len_limit-1+0.5-mod_mu)/mod_si;
    double x_min2 = ((double)len_min-0.5-mod_mu)/mod_si;
    if(1||from==-1||to==-1){//fix this uncomment out.
        from = 0;
        to = reads->size();
    }
    //    fprintf(stderr,"from: %d to:%d\n",from,to)
    for(int i=from;i<to;i++){
        bam1_t *b = (*reads)[i];
        
        memset(reconstructedRef,0,512);
        memset(myread,'N',512);
        memset(myrefe,'N',512);
        
        if (seq_ref != NULL){
            wrapperwithref(b,hdr,myread,myrefe,seq_ref);
        }else{
            reconstructRefWithPosHTS(b,mypair,reconstructedRef);
            wrapper(b,mypair.first->s,mypair.second,0,0,NULL,NULL,1,myread,myrefe);
        }
        
        if (len_limit <=0){
            len_limit = 512;
        };
        
        if (b->core.l_qseq>=30 && b->core.l_qseq<len_limit) {
            
            for (int cycle=0;cycle<b->core.l_qseq;cycle++)
                yourqual[cycle] = 254;
            
            for (int cycle=0;cycle<b->core.l_qseq;cycle++){
                refeBase = refToChar[myrefe[cycle]];
                readBase = refToChar[myread[cycle]];
                
                //if(refeBase!=4 && readBase!=4){
                int dist5p=cycle;
                int dist3p=b->core.l_qseq-1-cycle;
                if( bam_is_rev(b) ){
                    refeBase=com[refeBase];
                    readBase=com[readBase];
                    
                    dist5p=int(b->core.l_qseq)-1-cycle;
                    dist3p=cycle;
                }
                yourread[dist5p] = readBase;
                yourrefe[dist5p] = refeBase;
                yourqual[dist5p] = bam_get_qual(b)[cycle];
                //}
            }
            
            double l_anc, l_err;
            int L = b->core.l_qseq;
            l_err = ErrorLik(yourrefe, yourread, L, yourqual);
            if (!strcasecmp("b",model)){
	      l_anc = PMDLik_b(yourrefe, yourread, L, lambda, delta, delta_s, nv, yourqual,Tol); // Ancient Likelihood based on biotin model
            }else if(!strcasecmp("nb",model)){
	      l_anc = 0.5*PMDLik_b(yourrefe, yourread, L, lambda, delta, delta_s, nv, yourqual,Tol)+0.5*PMDLik_nb(yourrefe, yourread, L, lambda, delta, delta_s, nv, yourqual,Tol); // Ancient Likelihood based on non-biotin model
            }else{
                fprintf(stderr,"Please specify a deamination model for further calculations.\n");
                exit(0);
            }
            double y_max1 = (std::min((double)len_limit-1,std::max((double)len_min,(double)L))+0.5-anc_mu)/anc_si;
            double y_min1 = (std::max((double)len_min,std::min((double)len_limit-1,(double)L))-0.5-anc_mu)/anc_si;
            double y_max2 = (std::min((double)len_limit-1,std::max((double)len_min,(double)L))+0.5-mod_mu)/mod_si;
            double y_min2 = (std::max((double)len_min,std::min((double)len_limit-1,(double)L))-0.5-mod_mu)/mod_si;
            
            double denom = l_anc*NormalINC(y_max1, y_min1, x_max1, x_min1)*(1-eps)+l_err*NormalINC(y_max2, y_min2, x_max2, x_min2)*eps;
            y[0][0] += l_anc*(1-eps)*NormalINC_hess_mu2(y_max1, y_min1, x_max1, x_min1, anc_mu, anc_si)/denom - pow(l_anc*(1-eps)*NormalINC_grad_mu(y_max1, y_min1, x_max1, x_min1, anc_mu, anc_si)/denom,2);
            y[2][2] += l_err*eps*NormalINC_hess_mu2(y_max2, y_min2, x_max2, x_min2, mod_mu, mod_si)/denom - pow(l_err*eps*NormalINC_grad_mu(y_max2, y_min2, x_max2, x_min2, mod_mu, mod_si)/denom,2);
            y[1][1] += l_anc*(1-eps)*NormalINC_hess_si2(y_max1, y_min1, x_max1, x_min1, anc_mu, anc_si)/denom - pow(l_anc*(1-eps)*NormalINC_grad_si(y_max1, y_min1, x_max1, x_min1, anc_mu, anc_si)/denom,2);
            y[3][3] += l_err*eps*NormalINC_hess_si2(y_max2, y_min2, x_max2, x_min2, mod_mu, mod_si)/denom - pow(l_err*eps*NormalINC_grad_si(y_max2, y_min2, x_max2, x_min2, mod_mu, mod_si)/denom,2);
            y[0][1] += l_anc*(1-eps)*NormalINC_hess_mu_si(y_max1, y_min1, x_max1, x_min1, anc_mu, anc_si)/denom - l_anc*(1-eps)*NormalINC_grad_mu(y_max1, y_min1, x_max1, x_min1, anc_mu, anc_si)*l_anc*(1-eps)*NormalINC_grad_si(y_max1, y_min1, x_max1, x_min1, anc_mu, anc_si)/pow(denom,2);
            y[2][3] += l_err*eps*NormalINC_hess_mu_si(y_max2, y_min2, x_max2, x_min2, mod_mu, mod_si)/denom - l_err*eps*NormalINC_grad_mu(y_max2, y_min2, x_max2, x_min2, mod_mu, mod_si)*l_err*eps*NormalINC_grad_si(y_max2, y_min2, x_max2, x_min2, mod_mu, mod_si)/pow(denom,2);
            y[0][2] += - l_anc*(1-eps)*NormalINC_grad_mu(y_max1, y_min1, x_max1, x_min1, anc_mu, anc_si)*l_err*eps*NormalINC_grad_mu(y_max2, y_min2, x_max2, x_min2, mod_mu, mod_si)/pow(denom,2);
            y[0][3] += - l_anc*(1-eps)*NormalINC_grad_mu(y_max1, y_min1, x_max1, x_min1, anc_mu, anc_si)*l_err*eps*NormalINC_grad_si(y_max2, y_min2, x_max2, x_min2, mod_mu, mod_si)/pow(denom,2);
            y[1][2] += - l_anc*(1-eps)*NormalINC_grad_si(y_max1, y_min1, x_max1, x_min1, anc_mu, anc_si)*l_err*eps*NormalINC_grad_mu(y_max2, y_min2, x_max2, x_min2, mod_mu, mod_si)/pow(denom,2);
            y[1][3] += - l_anc*(1-eps)*NormalINC_grad_si(y_max1, y_min1, x_max1, x_min1, anc_mu, anc_si)*l_err*eps*NormalINC_grad_si(y_max2, y_min2, x_max2, x_min2, mod_mu, mod_si)/pow(denom,2);
        }
        //    fprintf(stderr,"\nmyread:\n%.*s\nmyReference:\n%.*s\n",b->core.l_qseq,myread,b->core.l_qseq,myrefe);
        //    fprintf(stderr,"---read[%d]----\n",nproc1);
    }
    y[1][0] = y[0][1];
    y[3][2] = y[2][3];
    y[2][0] = y[0][2];
    y[3][0] = y[0][3];
    y[2][1] = y[1][2];
    y[3][1] = y[1][3];
}

void *tsk_all_loglike_recalibration_hess_slave(void *dats){
    tsk_struct *ts = (tsk_struct *) &(my_tsk_struct[(size_t) dats]);
    tsk_loglike_recalibration_hess(ts->x,ts->llh_result_hess,ts->reads,ts->from,ts->to,ts->hdr,ts->seq_ref,ts->len_limit, ts->len_min,ts->model,ts->eps,ts->lambda,ts->delta,ts->delta_s,ts->nv,ts->Tol);
    pthread_exit(NULL);
}
