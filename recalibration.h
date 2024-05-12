#ifndef RECALIBRATION_H
#define RECALIBRATION_H
extern int tsk_nthreads;
typedef struct{//this struct is used in recalibration
  int len_limit;
  int len_min;
  double eps;
  int from;
  int to;
  double llh_result;
  double *x;
  double *llh_result_grad; //Add the gradient feature
  double **llh_result_hess; //Add the hessian matrix feature
  int threadid;
  double **mat;
  int counter[2];//<- ncall,ncall_gradient
  int nthreads;
}tsk_struct;


void *tsk_all_loglike_recalibration_slave(void *dats);
void *tsk_all_loglike_recalibration_grad_slave(void *dats);
void *tsk_all_loglike_recalibration_hess_slave(void *dats);
double tsk_all_loglike_recalibration(const double *x, const void *dats);
void tsk_all_loglike_recalibration_grad(const double *x, double *y,const void *dats);

double like_master(const double *xs,const void *);
void like_grad_master(const double *xs,double *y,const void *);
void like_hess_master(const double *xs,double **y);
#endif
