#include <cstring>
#include <cstdlib>
#include <cstdio>
#include <cassert>

#include "profile.h"
#include "bfgs.h"
#include "htslib/bgzf.h"

#include "misc.h"
#include "Likelihood.h"
#include "ngsBriggs.h"
#include "read_all_reads.h"

extern tsk_struct *my_tsk_struct;

//The log-likelihood for recalibration the ancient probs
double tsk_loglike_recalibration(const double *x, double **mat,int from,int to,int len_limit, int len_min, double eps){
  //fprintf(stderr,"(%f,%f,%f,%f) seq_ref: %p from: %d to:%d\n\n",x[0],x[1],x[2],x[3],seq_ref,from,to);
    double anc_mu = x[0];
    double anc_si = x[1];
    double mod_mu = x[2];
    double mod_si = x[3];
    double ll = 0;
   
    double x_max1 = ((double)len_limit-1+0.5-anc_mu)/anc_si;
    double x_min1 = ((double)len_min-0.5-anc_mu)/anc_si;
    double x_max2 = ((double)len_limit-1+0.5-mod_mu)/mod_si;
    double x_min2 = ((double)len_min-0.5-mod_mu)/mod_si;
    assert(from!=-1&&to!=-1);

   
    for(int i=from;i<to;i++) {
       

	double l_err = mat[0][i];
	double l_anc = mat[1][i];
	double L = mat[2][i];
	
	//	fprintf(stderr,"[%d] l_anc: %f L_err: %f L: %f\n",i,l_anc,l_err,L);	
	double y_max1 = (std::min((double)len_limit-1,std::max((double)len_min,(double)L))+0.5-anc_mu)/anc_si;
	double y_min1 = (std::max((double)len_min,std::min((double)len_limit-1,(double)L))-0.5-anc_mu)/anc_si;
	double y_max2 = (std::min((double)len_limit-1,std::max((double)len_min,(double)L))+0.5-mod_mu)/mod_si;
	double y_min2 = (std::max((double)len_min,std::min((double)len_limit-1,(double)L))-0.5-mod_mu)/mod_si;
	ll += log(l_anc*NormalINC(y_max1, y_min1, x_max1, x_min1)*(1-eps)+l_err*NormalINC(y_max2, y_min2, x_max2, x_min2)*eps);
	
    }
    //   exit(0);
    return -ll;
}

double tsk_all_loglike_recalibration(const double *x, const void *dats){
    //  fprintf(stderr,"[%s]\n",__FUNCTION__);
    ncalls++;
    tsk_struct *ts = (tsk_struct *) dats;
    double check = tsk_loglike_recalibration(x, ts->mat,ts->from,ts->to, ts->len_limit, ts->len_min, ts->eps);
    //  fprintf(stderr,"[%s]: total like: %f\n",__FUNCTION__,check);
    return check;
}

void *tsk_all_loglike_recalibration_slave(void *dats){
    //fprintf(stderr,"[%s]\n",__FUNCTION__);
    ncalls++;
    tsk_struct *ts = (tsk_struct *) &(my_tsk_struct[(size_t) dats]);
    ts->llh_result = tsk_loglike_recalibration(ts->x, ts->mat,ts->from,ts->to, ts->len_limit, ts->len_min, ts->eps);
    
    pthread_exit(NULL);
}


void tsk_loglike_recalibration_grad(const double *x, double *y, double **mat,int from, int to,int len_limit, int len_min, double eps){
  
    double anc_mu = x[0];
    double anc_si = x[1];
    double mod_mu = x[2];
    double mod_si = x[3];
    y[0] = 0;
    y[1] = 0;
    y[2] = 0;
    y[3] = 0;
    double ll = 0;
  
    double x_max1 = ((double)len_limit-1+0.5-anc_mu)/anc_si;
    double x_min1 = ((double)len_min-0.5-anc_mu)/anc_si;
    double x_max2 = ((double)len_limit-1+0.5-mod_mu)/mod_si;
    double x_min2 = ((double)len_min-0.5-mod_mu)/mod_si;
    assert(from!=-1&&to!=-1);

    for(int i=from;i<to;i++){
      
      double l_err = mat[0][i];
      double l_anc = mat[1][i];
      double L = mat[2][i];      
      
      double y_max1 = (std::min((double)len_limit-1,std::max((double)len_min,(double)L))+0.5-anc_mu)/anc_si;
      double y_min1 = (std::max((double)len_min,std::min((double)len_limit-1,(double)L))-0.5-anc_mu)/anc_si;
      double y_max2 = (std::min((double)len_limit-1,std::max((double)len_min,(double)L))+0.5-mod_mu)/mod_si;
      double y_min2 = (std::max((double)len_min,std::min((double)len_limit-1,(double)L))-0.5-mod_mu)/mod_si;
      
      
      y[0] -= l_anc*(1-eps)*NormalINC_grad_mu(y_max1, y_min1, x_max1, x_min1, anc_mu, anc_si)/(l_anc*NormalINC(y_max1, y_min1, x_max1, x_min1)*(1-eps)+l_err*NormalINC(y_max2, y_min2, x_max2, x_min2)*eps);
      y[1] -= l_anc*(1-eps)*NormalINC_grad_si(y_max1, y_min1, x_max1, x_min1, anc_mu, anc_si)/(l_anc*NormalINC(y_max1, y_min1, x_max1, x_min1)*(1-eps)+l_err*NormalINC(y_max2, y_min2, x_max2, x_min2)*eps);
      y[2] -= l_err*eps*NormalINC_grad_mu(y_max2, y_min2, x_max2, x_min2, mod_mu, mod_si)/(l_anc*NormalINC(y_max1, y_min1, x_max1, x_min1)*(1-eps)+l_err*NormalINC(y_max2, y_min2, x_max2, x_min2)*eps);
      y[3] -= l_err*eps*NormalINC_grad_si(y_max2, y_min2, x_max2, x_min2, mod_mu, mod_si)/(l_anc*NormalINC(y_max1, y_min1, x_max1, x_min1)*(1-eps)+l_err*NormalINC(y_max2, y_min2, x_max2, x_min2)*eps);
    }
}

void tsk_all_loglike_recalibration_grad(const double *x, double *y,const void *dats){
    ncalls_grad++;
    //tsk_struct *ts = &my_tsk_struct[0];//assuming single thread
    tsk_struct *ts = (tsk_struct *) dats;//assuming multiple threads
    tsk_loglike_recalibration_grad(x,y,ts->mat,ts->from,ts->to,ts->len_limit, ts->len_min,ts->eps);
    //fprintf(stderr,"[%s]\n",__FUNCTION__);
}

// Check here
void *tsk_all_loglike_recalibration_grad_slave(void *dats){
    //fprintf(stderr,"[%s]\n",__FUNCTION__);
    ncalls_grad++;
    tsk_struct *ts = (tsk_struct *) &(my_tsk_struct[(size_t) dats]);
    //ts->llh_grad_result;
    tsk_loglike_recalibration_grad(ts->x,ts->llh_result_grad,ts->mat,ts->from,ts->to,ts->len_limit, ts->len_min,ts->eps);
    //ts->llh_result = tsk_loglike_recalibration(ts->x, ts->reads,ts->from,ts->to,ts->hdr, ts->seq_ref, ts->len_limit, ts->model, ts->eps, ts->lambda, ts->delta, ts->delta_s, ts->nv,ts->threadid);
    pthread_exit(NULL);
}

void tsk_loglike_recalibration_hess(const double *x, double **y, double **mat,int from, int to,int len_limit, int len_min, double eps){
    //fprintf(stderr,"[%s]\n",__FUNCTION__);
    
    double anc_mu = x[0];
    double anc_si = x[1];
    double mod_mu = x[2];
    double mod_si = x[3];
    for(int i=0; i<3; i++){
        for(int j=0; j<3; j++){
            y[i][j] = 0;
        }
    }
    
    double x_max1 = ((double)len_limit-1+0.5-anc_mu)/anc_si;
    double x_min1 = ((double)len_min-0.5-anc_mu)/anc_si;
    double x_max2 = ((double)len_limit-1+0.5-mod_mu)/mod_si;
    double x_min2 = ((double)len_min-0.5-mod_mu)/mod_si;
    assert(from!=-1&&to!=-1);

    //    fprintf(stderr,"from: %d to:%d\n",from,to)
    for(int i=from;i<to;i++){

      double l_err = mat[0][i];
      double l_anc = mat[1][i];
      double L = mat[2][i];
      
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
    tsk_loglike_recalibration_hess(ts->x,ts->llh_result_hess,ts->mat,ts->from,ts->to,ts->len_limit, ts->len_min,ts->eps);
    pthread_exit(NULL);
}
