
#include <cstring>
#include <cstdlib>
#include <cstdio>
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
#include <array>  //Rasmus

#include "misc.h"
#include "Recalibration.h"
#include "Likelihood.h"
#include "PosteriorProb.h"
#include "ngsBriggs_cli.h"
#include "ngsBriggs.h"

#define PI acos(-1)
using namespace std;

argStruct *pars_briggs(int argc,char ** argv){
    argStruct *Briggspars = new argStruct;

    Briggspars->hts = NULL;
    Briggspars->tab = NULL;
    Briggspars->ref = NULL;
    Briggspars->len = NULL;
    Briggspars->ihts = NULL;
    Briggspars->iref = NULL;
    Briggspars->ohts = NULL;
    Briggspars->otab = NULL;
    Briggspars->oinf = NULL;
    Briggspars->olen = NULL;
    Briggspars->model = NULL;
    Briggspars->olik = NULL;
    Briggspars->eps = 0;
    Briggspars->isrecal = 0;
    Briggspars->nthread = 1;
    Briggspars->bdamage = NULL;
    Briggspars->rlens = NULL;


    ++argv;
    while(*argv){
        if(strcasecmp("-bam",*argv)==0){
            Briggspars->hts=strdup(*(++argv));
        }
        else if(strcasecmp("-tab",*argv)==0){
            Briggspars->tab=strdup(*(++argv));
        }
        else if(strcasecmp("-ref",*argv)==0){
            Briggspars->ref=strdup(*(++argv));
        }
     
        else if(strcasecmp("-len",*argv)==0){
            Briggspars->len=strdup(*(++argv));
        }
	else if(strcasecmp("-ibam",*argv)==0){
            Briggspars->ihts=strdup(*(++argv));
        }
        else if(strcasecmp("-iref",*argv)==0){
            Briggspars->iref=strdup(*(++argv));
        }
      
        else if(strcasecmp("-obam",*argv)==0){
            Briggspars->ohts=strdup(*(++argv));
        }
        else if(strcasecmp("-otab",*argv)==0){
            Briggspars->otab=strdup(*(++argv));
        }
        else if(strcasecmp("-oinf",*argv)==0){
            Briggspars->oinf=strdup(*(++argv));
        }
        else if(strcasecmp("-olen",*argv)==0){
            Briggspars->olen=strdup(*(++argv));
        }
        else if(strcasecmp("-model",*argv)==0){
            Briggspars->model=strdup(*(++argv));
        }
        else if(strcasecmp("-olik",*argv)==0){
            Briggspars->olik=strdup(*(++argv));
        }
        else if(strcasecmp("-eps",*argv)==0){
            Briggspars->eps=atof(*(++argv));
        }
        else if(strcasecmp("-isrecal",*argv)==0){
            Briggspars->isrecal=atoi(*(++argv));
        }
        else if(strcasecmp("-nthread",*argv)==0){
            Briggspars->nthread=atoi(*(++argv));
        }
        else if(strcasecmp("-bdamage",*argv)==0){
            Briggspars->bdamage=strdup(*(++argv));
        }
        else if(strcasecmp("-rlens",*argv)==0){
            Briggspars->rlens=strdup(*(++argv));
        }
        else{
            fprintf(stderr,"Unrecognized input option %s, see ngsBriggs help page\n\n",*(argv));
            return NULL;
        }
        ++argv;
    }
    return Briggspars;
}

void argStruct_destroy(argStruct *Briggspars){
    free(Briggspars->hts);
    free(Briggspars->tab);
    free(Briggspars->ref);
    free(Briggspars->len);
    free(Briggspars->ihts);
    free(Briggspars->iref);
    free(Briggspars->ohts);
    free(Briggspars->otab);
    free(Briggspars->oinf);
    free(Briggspars->olen);
    free(Briggspars->model);
    free(Briggspars->olik);
    free(Briggspars->bdamage);
    free(Briggspars->rlens);
    delete(Briggspars);
}
