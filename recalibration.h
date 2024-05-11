#ifndef RECALIBRATION_H
#define RECALIBRATION_H

void *tsk_all_loglike_recalibration_slave(void *dats);
void *tsk_all_loglike_recalibration_grad_slave(void *dats);
void *tsk_all_loglike_recalibration_hess_slave(void *dats);
double tsk_all_loglike_recalibration(const double *x, const void *dats);
void tsk_all_loglike_recalibration_grad(const double *x, double *y,const void *dats);

#endif
