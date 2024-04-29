#include <cstdio>
#include <cstring>
#include <cstdlib>

#include "ngsBriggs_cli.h"

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
    pars->model = NULL;
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
            pars->model=strdup(*(++argv));
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
    free(pars->model);
    free(pars->olik);
    free(pars->bdamage);
    free(pars->rlens);
    delete(pars);
}
