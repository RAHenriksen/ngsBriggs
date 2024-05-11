#ifndef NGSBRIGGSCLI_H
#define NGSBRIGGSCLI_H
int helppage(FILE *fp);
typedef struct{
  //filenames
  char *hts; // bam/cram/sam input to count the mismatch table for inference
  char *tab; // mismatch table input 
  char *len; // length distribution file, mass probability distribution, not the cdf
  char *ref;  // input reference file
  char *ihts; // bam/cram/sam for inference
  char *iref; // reference file for inference
  char *ohts; // bam/cram/sam output
  char *otab; // mismatch table output
  char *oinf; // model parameter output
  char *olen; // fragment len dist output
  int model; // model specific, either zero or one. zero=b one=nb
  char *olik; // output nucleotide likelihood file
  double eps; // contamination rate provided by external software
  int dorecal; // should we recalibrate the PD with the information of length as well as PMD model
  // I also want to add output the recalibrated genotype likelihood option.
  int nthread;
  char *bdamage; // length distribution file, mass probability distribution, not the cdf
  char *rlens; // length distribution file, mass probability distribution, not the cdf
}argStruct;
argStruct *pars_briggs(int argc,char ** argv);
void argStruct_destroy(argStruct *mypars);
#endif
