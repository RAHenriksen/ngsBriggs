#include <cstring>
#include <cstdlib>
#include <cstdio>
#include <cassert>
#include <Eigen/Core>
#include <Eigen/Eigenvalues> //sort and merge

#include "profile.h"
#include "bfgs.h"
#include "htslib/bgzf.h"

#include "misc.h"
#include "likelihood.h"
#include "read_all_reads.h"

#include "recalibration.h"

//extern tsk_struct *my_tsk_struct;

//The log-likelihood for recalibration the ancient probs
double tsk_loglike_recalibration(const double *x, double **mat,int from,int to,int len_limit, int len_min, double eps,int &counter){
  counter++;
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
   
    tsk_struct *ts = (tsk_struct *) dats;
    double check = tsk_loglike_recalibration(x, ts->mat,ts->from,ts->to, ts->len_limit, ts->len_min, ts->eps,ts->counter[0]);
    //  fprintf(stderr,"[%s]: total like: %f\n",__FUNCTION__,check);
    return check;
}

void *tsk_all_loglike_recalibration_slave(void *dats){
    //fprintf(stderr,"[%s]\n",__FUNCTION__);
    tsk_struct *ts = (tsk_struct *) dats;
    ts->llh_result = tsk_loglike_recalibration(ts->x, ts->mat,ts->from,ts->to, ts->len_limit, ts->len_min, ts->eps,ts->counter[0]);
    
    pthread_exit(NULL);
}


void tsk_loglike_recalibration_grad(const double *x, double *y, double **mat,int from, int to,int len_limit, int len_min, double eps,int &counter){
  counter++;
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
  //tsk_struct *ts = &my_tsk_struct[0];//assuming single thread
    tsk_struct *ts = (tsk_struct *) dats;//assuming multiple threads
    tsk_loglike_recalibration_grad(x,y,ts->mat,ts->from,ts->to,ts->len_limit, ts->len_min,ts->eps,ts->counter[1]);
    //fprintf(stderr,"[%s]\n",__FUNCTION__);
}

// Check here
void *tsk_all_loglike_recalibration_grad_slave(void *dats){
    //fprintf(stderr,"[%s]\n",__FUNCTION__);

  tsk_struct *ts = (tsk_struct *) dats;
  //  fprintf(stderr,"nthreads: %d\n",ts->nthreads);exit(0);
    //ts->llh_grad_result;
    tsk_loglike_recalibration_grad(ts->x,ts->llh_result_grad,ts->mat,ts->from,ts->to,ts->len_limit, ts->len_min,ts->eps,ts->counter[1]);
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
    }
    y[1][0] = y[0][1];
    y[3][2] = y[2][3];
    y[2][0] = y[0][2];
    y[3][0] = y[0][3];
    y[2][1] = y[1][2];
    y[3][1] = y[1][3];
}

void *tsk_all_loglike_recalibration_hess_slave(void *dats){
    tsk_struct *ts = (tsk_struct *) dats;
    //  fprintf(stderr,"ts->nthreads: %d\n",ts->nthreads);exit(0);
    tsk_loglike_recalibration_hess(ts->x,ts->llh_result_hess,ts->mat,ts->from,ts->to,ts->len_limit, ts->len_min,ts->eps);
    pthread_exit(NULL);
}



double like_master(const double *xs,const void *ptr){
  fprintf(stderr,"[%s] like_master: (%f,%f,%f,%f) \n",__FUNCTION__,xs[0],xs[1],xs[2],xs[3]);
  //  fprintf(stderr,"[%s] like_master\n",__FUNCTION__);
  tsk_struct *ts = (tsk_struct*) ptr;
  for(int i=0;i<ts->nthreads;i++)
    for(int ii=0;ii<4;ii++)
      ts[i].x[ii] = xs[ii];
    //fprintf(stderr,"mu_anc\tsig_anc\tmu_mod\tsig_mod\n");
    //fprintf(stderr,"%f\t%f\t%f\t%f\n",xs[0],xs[1],xs[2],xs[3]);
    pthread_t thd[ts->nthreads];
    for(size_t i=0;i<ts->nthreads;i++){
        int rc = pthread_create(&thd[i],NULL,tsk_all_loglike_recalibration_slave,(void*) &ts[i]);
        if(rc)
            fprintf(stderr,"Error creating thread\n");
        
    }
    for(int i=0;i<ts->nthreads;i++)
        pthread_join(thd[i], NULL);
    
    double res=0;
    for(int i=0;i<ts->nthreads;i++){
      fprintf(stderr,"[%s] lik[%d,%d]=%f\n",__FUNCTION__,ts[i].from,ts[i].to,ts[i].llh_result);
        res += ts[i].llh_result;
    }
    fprintf(stderr,"[%s] total lik[%d,%d]: %f\n",__FUNCTION__,ts[0].from,ts[ts->nthreads-1].to,res);
    // exit(0);
    return res;
}

void like_grad_master(const double *xs,double *y,const void *ptr){
  //  fprintf(stderr,"[%s] like_master: (%f,%f,%f,%f) \n",__FUNCTION__,xs[0],xs[1],xs[2],xs[3]);
    //  fprintf(stderr,"like_master\n");
  tsk_struct *ts = (tsk_struct *) ptr;
    for(int i=0;i<ts->nthreads;i++)
        for(int ii=0;ii<4;ii++)
            ts[i].x[ii] = xs[ii];
    pthread_t thd[ts->nthreads];
    for(size_t i=0;i<ts->nthreads;i++){
        //int rc = pthread_create(&thd[i],NULL,tsk_All_loglike_recalibration_slave,(void*) i);
        int rc = pthread_create(&thd[i],NULL,tsk_all_loglike_recalibration_grad_slave,(void*) &ts[i]);
        if(rc)
            fprintf(stderr,"Error creating thread\n");
        
    }
    for(int i=0;i<ts->nthreads;i++)
        pthread_join(thd[i], NULL);
    
    for(int j=0;j<4;j++){
        y[j] = 0;
    }
    for(int i=0;i<ts->nthreads;i++){
      //    fprintf(stderr,"[%s] lik_grad[%d,%d]\n",__FUNCTION__,my_tsk_struct[i].from,my_tsk_struct[i].to);
        for (int j=0;j<4;j++){
            y[j] += ts[i].llh_result_grad[j];
        }
    }
    //  fprintf(stderr,"[%s] total lik_grad[%d,%d]\n",__FUNCTION__,my_tsk_struct[0].from,my_tsk_struct[nThreads-1].to);
    //cout<<y[0]<<"\t"<<y[1]<<"\t"<<y[2]<<"\t"<<y[3]<<"\n";
}

void like_hess_master(const double *xs,double **y,tsk_struct my_tsk_struct[]){
    Eigen::Matrix4d A, B;
    //  fprintf(stderr,"like_master\n");???
    int nThreads = my_tsk_struct[0].nthreads;
    for(int i=0;i<nThreads;i++)
        for(int ii=0;ii<4;ii++)
            my_tsk_struct[i].x[ii] = xs[ii];
    pthread_t thd[nThreads];
    for(size_t i=0;i<nThreads;i++){
        //int rc = pthread_create(&thd[i],NULL,tsk_All_loglike_recalibration_slave,(void*) i);
        int rc = pthread_create(&thd[i],NULL,tsk_all_loglike_recalibration_hess_slave,(void*) &my_tsk_struct[i]);
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
      //   fprintf(stderr,"lik_hess[%d,%d]\n",my_tsk_struct[i].from,my_tsk_struct[i].to);
        for (int j=0;j<4;j++){
            for (int k=0;k<4;k++){
                A(j,k) += my_tsk_struct[i].llh_result_hess[j][k];
            }
            //y[j] += my_tsk_struct[i].llh_result_grad[j];
        }
    }
    //cout<<A<<"\n";
    //   fprintf(stderr,"total lik_hess[%d,%d]\n",my_tsk_struct[0].from,my_tsk_struct[nThreads-1].to);
    B = -A.inverse();
    for (int i=0;i<4;i++){
        for (int j=0;j<4;j++){
            y[i][j] = B(i,j);
        }
    }
}
