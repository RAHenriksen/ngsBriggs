#ifndef RECALIBRATION_H
#define RECALIBRATION_H

void *tsk_all_loglike_recalibration_slave(void *dats);
void *tsk_all_loglike_recalibration_grad_slave(void *dats);
void *tsk_all_loglike_recalibration_hess_slave(void *dats);
double tsk_all_loglike_recalibration(const double *x, const void *dats);
void tsk_all_loglike_recalibration_grad(const double *x, double *y,const void *dats);

double like_master(const double *xs,const void *);
void like_grad_master(const double *xs,double *y,const void *);
void like_hess_master(const double *xs,double **y);
#endif
