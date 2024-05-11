#include <cstdio>
#include <cstring>
#include <cstdlib>

#include "ngsBriggs_cli.h"

int helppage(FILE *fp){
    fprintf(fp,"\t-> ./ngsBriggs -bam -tab -ref -len -ibam -iref -obam -otab -oinf -olen -model -eps -isrecal -olik -nthreads\n");
    fprintf(fp,"\t-> -bam: The bam file for inference;\n");
    fprintf(fp,"\t-> -tab: The table file for inference;\n");
    fprintf(fp,"\t-> -ref: The reference file for inference;\n");
    fprintf(fp,"\t-> -len: The length mass probability distribution file for inference;\n");
    fprintf(fp,"\t-> -ibam: The bam file for ancient strand fishing;\n");
    fprintf(fp,"\t-> -iref: The reference file for ancient strand fishing;\n");
    fprintf(fp,"\t-> -obam: The output bam file name;\n");
    fprintf(fp,"\t-> -otab: The output table file name;\n");
    fprintf(fp,"\t-> -oinf: The output inferred parameters file name;\n");
    fprintf(fp,"\t-> -olen: The output length mass probability distribution file name;\n");
    fprintf(fp,"\t-> -model: Specifying the model, either b (the biotin model) or nb (the non-biotin model);\n");
    fprintf(fp,"\t-> -eps: The overall modern contamination rate, the value should be within the interval [0,1);\n");
    fprintf(fp,"\t-> -isrecal: Choose 1 if recalibration based on length is needed, otherwise 0 (default);\n");
    fprintf(fp,"\t-> -olik: The output nucleotide likelihood file;\n");
    fprintf(fp,"\t-> -nthreads: Choose the number of threads to speed up the recalibration process.\n");
    fprintf(fp,"\t-> -bdamage: The mismatch matrix in bdamage format for metagenomic framework.\n");
    fprintf(fp,"\t-> -rlens: The read length distributions for metagenomic framework.\n");
    fprintf(fp,"\t-> Examples.\n");
    fprintf(fp,"\t-> ./ngsbriggs -model nb -bdamage Chr22_024_36_68_0097.bdamage.gz -rlens Chr22_024_36_68_0097.rlens.gz\n");
    fprintf(fp,"\t-> ./ngsbriggs -tab mismatch2.txt -len len.txt -model nb\n");
    fprintf(fp,"\t-> ./ngsbriggs -bam Chr22_024_36_68_0097.sorted.MD.bam -ref chr22.fa -model nb\n");
    exit(1);
  return 0;
}


argStruct *pars_briggs(int argc,char ** argv){
    argStruct *pars = new argStruct;

    pars->hts = NULL;
    pars->tab = NULL;
    pars->ref = NULL;
    pars->len = NULL;
    pars->ihts = NULL;
    pars->iref = NULL;
    pars->ohts = NULL;
    pars->otab = NULL;
    pars->oinf = NULL;
    pars->olen = NULL;
    pars->model = 1;
    pars->olik = NULL;
    pars->eps = 0;
    pars->dorecal = 0;
    pars->nthread = 1;
    pars->bdamage = NULL;
    pars->rlens = NULL;


    ++argv;
    while(*argv){
        if(strcasecmp("-bam",*argv)==0){
            pars->hts=strdup(*(++argv));
        }
        else if(strcasecmp("-tab",*argv)==0){
            pars->tab=strdup(*(++argv));
        }
        else if(strcasecmp("-ref",*argv)==0){
            pars->ref=strdup(*(++argv));
        }
     
        else if(strcasecmp("-len",*argv)==0){
            pars->len=strdup(*(++argv));
        }
	else if(strcasecmp("-ibam",*argv)==0){
            pars->ihts=strdup(*(++argv));
        }
        else if(strcasecmp("-iref",*argv)==0){
            pars->iref=strdup(*(++argv));
        }
      
        else if(strcasecmp("-obam",*argv)==0){
            pars->ohts=strdup(*(++argv));
        }
        else if(strcasecmp("-otab",*argv)==0){
            pars->otab=strdup(*(++argv));
        }
        else if(strcasecmp("-oinf",*argv)==0){
            pars->oinf=strdup(*(++argv));
        }
        else if(strcasecmp("-olen",*argv)==0){
            pars->olen=strdup(*(++argv));
        }
        else if(strcasecmp("-model",*argv)==0){
	  ++argv;
	  if(strcasecmp(*argv,"b")==0)
	    pars->model=0;
	  else if(strcasecmp(*argv,"nb")==0)
	    pars->model=1;
	  else{
	    fprintf(stderr,"\t-> Unknown option: \'%s\'\n",*argv);
	    exit(0);
	  }
        }
        else if(strcasecmp("-olik",*argv)==0){
            pars->olik=strdup(*(++argv));
        }
        else if(strcasecmp("-eps",*argv)==0){
            pars->eps=atof(*(++argv));
        }
        else if(strcasecmp("-isrecal",*argv)==0){
            pars->dorecal=atoi(*(++argv));
        }
        else if(strcasecmp("-nthread",*argv)==0){
            pars->nthread=atoi(*(++argv));
        }
        else if(strcasecmp("-bdamage",*argv)==0){
            pars->bdamage=strdup(*(++argv));
        }
        else if(strcasecmp("-rlens",*argv)==0){
            pars->rlens=strdup(*(++argv));
        }
        else{
            fprintf(stderr,"Unrecognized input option %s, see ngsBriggs help page\n\n",*(argv));
            return NULL;
        }
        ++argv;
    }

    return pars;
}

void argStruct_destroy(argStruct *pars){
    free(pars->hts);
    free(pars->tab);
    free(pars->ref);
    free(pars->len);
    free(pars->ihts);
    free(pars->iref);
    free(pars->ohts);
    free(pars->otab);
    free(pars->oinf);
    free(pars->olen);
    free(pars->olik);
    free(pars->bdamage);
    free(pars->rlens);
    delete(pars);
}
