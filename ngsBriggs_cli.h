#ifndef NGSBRIGGSCLI_H
#define NGSBRIGGSCLI_H
#include <cstdlib>
#include <ctime>
#include <cstring>
#include "HelpPage.h"

void FragArrayReader(int len_limit, int& number, int*& Length, double *& Freq, const char* filename);

typedef struct{
  //filenames
  char *hts; // bam/cram/sam input to count the mismatch table for inference
  char *tab; // mismatch table input 
  char *len; // length distribution file, mass probability distribution, not the cdf
  char *ref;  // input reference file
  char *bed; // input bed file
  char *chr; // chromosome name
  char *ihts; // bam/cram/sam for inference
  char *iref; // reference file for inference
  char *ibed; // bed file for inference
  char *ichr; // chromosome name for inference
  char *ohts; // bam/cram/sam output
  char *otab; // mismatch table output
  char *oinf; // model parameter output
  char *olen; // fragment len dist output
  char *model; // model specific, either b or nb
  char *olik; // output nucleotide likelihood file
  double eps; // contamination rate provided by external software
  int isrecal; // should we recalibrate the PD with the information of length as well as PMD model
  // I also want to add output the recalibrated genotype likelihood option.
  int nthread;
}argStruct;
argStruct *pars_briggs(int argc,char ** argv);
void argStruct_destroy(argStruct *mypars);
#endif
