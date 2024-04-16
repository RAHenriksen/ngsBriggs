#include <cstring>
#include <cstdlib>
#include <cstdio>
#include <cassert>
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
#include <array>

#include "profile.h"
#include "bfgs.h"
#include "htslib/bgzf.h"
#include "briggs_writer.h"
#include "read_all_reads.h"

#include "misc.h"
#include "Likelihood.h"
#include "ngsBriggs.h"

extern tsk_struct *my_tsk_struct;


//The log-likelihood for recalibration the ancient prob
double loglike_recalibration(const double *x, char *refName,char *fname, const char* chromname, const char* bedname,int mapped_only,int se_only, int mapq, faidx_t *seq_ref,int len_limit, int len_min, char * model, double eps, double lambda, double delta, double delta_s, double nv,std::string s,double Tol){
    fprintf(stderr,"mapped_only: %d [%s]\n",mapped_only,__FUNCTION__);
    htsFormat *dingding7 =(htsFormat*) calloc(1,sizeof(htsFormat));
    double anc_mu = x[0];
    double anc_si = x[1];
    double mod_mu = x[2];
    double mod_si = x[3];
    double ll = 0;
    char reconstructedRef[512];
    char myread[512];
    char myrefe[512];
    char yourread[512];
    char yourrefe[512];
    int  yourqual[512];
    std::pair< kstring_t*, std::vector<int> >  mypair;
    kstring_t *kstr =new kstring_t;
    kstr->l=kstr->m=0;
    kstr->s=NULL;
    mypair.first = kstr;
    samFile *in=NULL;

    int nthreads = 1;

    char refeBase, readBase;
    //code below only relevant if using cram files
    if(refName!=NULL){
        char *ref =(char*) malloc(10 + strlen(refName) + 1);
        sprintf(ref, "reference=%s", refName);
        hts_opt_add((hts_opt **)&dingding7->specific,ref);
        free(ref);
    }
    if(strstr(fname,".cram")!=NULL && refName==NULL){
        fprintf(stderr,"\t-> cram file requires reference with -T FILE.fa \n");
        exit(0);
    }
    if((in=sam_open_format(fname,"r",dingding7))==NULL ){
        fprintf(stderr,"[%s] nonexistant file: %s\n",__FUNCTION__,fname);
        exit(0);
    }
    bam_hdr_t  *hdr = sam_hdr_read(in);
    bam1_t *b = bam_init1();
    
    int chrom_num = hdr->n_targets;
    size_t * chrom_line = (size_t*)malloc(chrom_num*(sizeof(size_t)));
    std::vector<std::array<size_t, 2> > bedsites;
    
    // loading the bedfile
    if (bedname!=NULL){
        //fprintf(stderr,"Loading the bedfile %s ...\n",bedname);
        BGZF *fp = NULL;
        fp = bgzf_open(bedname,"rb");
        
        kstring_t *kstr1 = new kstring_t;
        kstr1->s = NULL;
        kstr1->l = kstr1->m = 0;
        int line=0;
        std::string word0, word;
        bgzf_getline(fp,'\n',kstr1);
        for (int j=0; j<chrom_num; j++){
            chrom_line[j] = line;
            do{
                //while(bgzf_getline(fp,'\n',kstr1)>0){
                std::istringstream iss(kstr1->s);
                //string word0, word;
                getline(iss,word0,'\t');
                size_t num;
                if  (strcmp(word0.c_str(),hdr->target_name[j])==0){
                    line++;
                }
                //size_t bedsite[2];
                std::array<size_t, 2> bedsite;
                int idx = 1;
                while (getline(iss,word,'\t')){
                    //sscanf(word.c_str(), "%zu", &num);
                    if (strcmp(word0.c_str(),hdr->target_name[j])==0){
                        sscanf(word.c_str(), "%zu", &num);
                        idx = 1-idx;
                        bedsite[idx] = num;
                        if (idx == 1){
                            bedsites.push_back(bedsite);
                        }
                    }
                }
                //cout << "\n";
            }while (strcmp(word0.c_str(),hdr->target_name[j])==0 && bgzf_getline(fp,'\n',kstr1)>0);
        }
    }
    
    
    int ret;
    int refId=-1;
    double num = 0.0;
    size_t max_site;
    nproc1 = 0;
    int chromId=-1;
    if (chromname!=NULL){
        for (int j=0; j<chrom_num; j++){
            if (!strcmp(chromname,hdr->target_name[j])){
                chromId = j;
                //fprintf(stderr,"We will focus on Chromosome %s!\n", hdr->target_name[j]);
            }
        }
    }
    if (chromId==-1){
        //fprintf(stderr,"No meaningful chromosome name is provided, therefore we will focus on all provided chromosomes!\n");
    }
    
    uchar * indref = NULL;
    double x_max1 = ((double)len_limit-1+0.5-anc_mu)/anc_si;
    double x_min1 = ((double)len_min-0.5-anc_mu)/anc_si;
    double x_max2 = ((double)len_limit-1+0.5-mod_mu)/mod_si;
    double x_min2 = ((double)len_min-0.5-mod_mu)/mod_si;
    while(((ret=sam_read1(in,hdr,b)))>0) {
      fprintf(stderr,"looping through bam\n");
      //nproc1++;
        //if -m was used, discard unmapped reads
      if(mapped_only!=0){
	if(b->core.flag&4)
	  continue;
      }
      
      //default -a 1
      //only use single end reads
      //which is either single end or collapsed reads
      //#if 0
      if(se_only==1){
	if(b->core.flag&1)//if paired continue
	  continue;
      }
      //#endif
      if(b->core.flag>256)
	continue; //Remove possible duplicates
      
      // if mapq threshold is set
      if(mapq!=-1 && b->core.qual<mapq)
	continue;
      if (chromId==-1||chromId==b->core.tid){
	nproc1++;
      }else if (chromId!=b->core.tid){
	continue;
      }
      if(refId==-1||refId!=b->core.tid){
	refId=b->core.tid;
	//fprintf(stderr,"\t-> Now at Chromosome: %s\n",hdr->target_name[refId]);
	if (bedname!=NULL){
	  if (indref!=NULL){
	    free(indref);
	  }
	  //checkchrome=b->core.tid;
	  max_site = 0;
	  size_t max_line =  refId < chrom_num-1 ? chrom_line[refId+1] : bedsites.size();
	  for (size_t l=chrom_line[refId]; l < max_line; l++){
	    if (max_site < bedsites[l][1]){
	      max_site = bedsites[l][1]; //Assume the provided bed file is 1-based (but will treat it as 0-based internally)
	    }
	  }
	  size_t chr_len = max_site; // The maximum position considered in the bed file along this chromosome
          
	  //faidx_seq_len(seq_ref,hdr->target_name[refId]);
	  indref = (uchar*) malloc(chr_len*sizeof(uchar));
	  for (size_t l=0; l<chr_len; l++){
	    indref[l] = (uchar)0;
	  }
          
	  for (size_t l=chrom_line[refId]; l < max_line; l++){
	    for (size_t i = bedsites[l][0]; i <= bedsites[l][1]; i++){
	      indref[i-1] = (uchar)1; // Shift by 1 so that it is treated 0 based internally.
	    }
	  }
	}
      }
      
      //then we simply write it to the output
      memset(reconstructedRef,0,512);
      memset(myread,'N',512);
      memset(myrefe,'N',512);
      
      if (seq_ref != NULL){
	wrapperwithref(b,hdr,myread,myrefe,seq_ref);
      }else{
	reconstructRefWithPosHTS(b,mypair,reconstructedRef);
	wrapper(b,mypair.first->s,mypair.second,0,0,NULL,NULL,1,myread,myrefe);
      }
      
      if (len_limit <=0){
	len_limit = 512;
      };
      //cout << b->core.l_qseq << "Test nuc_llik\n";
      if (b->core.l_qseq>=30 && b->core.l_qseq<len_limit){
	for (int cycle=0;cycle<b->core.l_qseq;cycle++){
	  yourqual[cycle] = -1;
	}
	//cout << b->core.l_qseq << "\n";
	int isbedfrag = 0;
	for (int cycle=0;cycle<b->core.l_qseq;cycle++){
	  refeBase = refToChar[myrefe[cycle]];
	  readBase = refToChar[myread[cycle]];
	  size_t pos = b->core.pos+cycle;
	  //                cout << (int)refeBase << " " << (int)readBase<<" "<<pos<<" "<<max_site<<"\n";
	  if(refeBase!=4 && readBase!=4 && bedname!=NULL &&  pos < max_site && indref[pos]==1){ //check whether the fragment has intersection with the bed file
	    isbedfrag += 1;
	  }else if(refeBase!=4 && readBase!=4 && bedname==NULL){
	    isbedfrag += 1;
	  }
	}
	if (isbedfrag > 0){
	  
	  for (int cycle=0;cycle<b->core.l_qseq;cycle++){
	    int dist5p=cycle;
	    int dist3p=b->core.l_qseq-1-cycle;
	    //cout<<"flag "<<b->core.flag<<"\n";
	    if( bam_is_rev(b) ){
	      // cout<<"rev "<<"\t";
	      refeBase=com[refeBase];
	      readBase=com[readBase];
	      //dist5p=int(al.QueryBases.size())-1-i;
	      dist5p=int(b->core.l_qseq)-1-cycle;
	      dist3p=cycle;
	    }
	    yourread[dist5p] = readBase;
	    yourrefe[dist5p] = refeBase;
	    yourqual[dist5p] = bam_get_qual(b)[cycle];
	  }
          
	  double l_anc, l_err;
	  int L = b->core.l_qseq;
	  l_err = ErrorLik(yourrefe, yourread, L, yourqual);
	  if (!strcasecmp("b",model)){
	    l_anc = PMDLik_b(yourrefe, yourread, L, lambda, delta, delta_s, nv, yourqual,Tol); // Ancient Likelihood based on biotin model
	  }else if(!strcasecmp("nb",model)){
	    l_anc = 0.5*PMDLik_b(yourrefe, yourread, L, lambda, delta, delta_s, nv, yourqual,Tol)+0.5*PMDLik_nb(yourrefe, yourread, L, lambda, delta, delta_s, nv, yourqual,Tol); // Ancient Likelihood based on non-biotin model
	  }else{
	    fprintf(stderr,"Please specify a deamination model for further calculations.\n");
	    return -1;
	  }
	  double y_max1 = (std::min((double)len_limit-1,std::max((double)len_min,(double)L))+0.5-anc_mu)/anc_si;
	  double y_min1 = (std::max((double)len_min,std::min((double)len_limit-1,(double)L))-0.5-anc_mu)/anc_si;
	  double y_max2 = (std::min((double)len_limit-1,std::max((double)len_min,(double)L))+0.5-mod_mu)/mod_si;
	  double y_min2 = (std::max((double)len_min,std::min((double)len_limit-1,(double)L))-0.5-mod_mu)/mod_si;
	  ll += log(l_anc*NormalINC(y_max1, y_min1, x_max1, x_min1)*(1-eps)+l_err*NormalINC(y_max2, y_min2, x_max2, x_min2)*eps);
	  //double PostAncProb = AncProb(yourrefe, yourread, b->core.l_qseq, lambda, delta, delta_s, nv, yourqual, model ,eps);
	  //double PostPMDProb = PMDProb(yourrefe, yourread, b->core.l_qseq, lambda, delta, delta_s, nv, yourqual, model);
	  //cout<<"\t"<<"AN: "<<PostAncProb<<" \t"<<"PD: "<<PostPMDProb<<"\n";
	  num += 1.0; //number of intersected reads.
	}
      }
    }

    if (indref!=NULL){
        free(indref);
    }

    bedsites.clear();
    bedsites.shrink_to_fit();
    bam_destroy1(b);
    
    hts_opt_free((hts_opt *)dingding7->specific);
    free(dingding7);
    sam_hdr_destroy(hdr);
    assert(sam_close(in)==0);

    return -ll;
}

//The log-likelihood for recalibration the ancient probs
double tsk_loglike_recalibration(const double *x, std::vector<bam1_t *> *reads,int from,int to,sam_hdr_t *hdr, faidx_t *seq_ref,int len_limit, int len_min, char * model, double eps, double lambda, double delta, double delta_s, double nv,int threadid,double Tol){
    //  fprintf(stderr,"(%f,%f,%f,%f) seq_ref: %p\n",x[0],x[1],x[2],x[3],seq_ref);
    double anc_mu = x[0];
    double anc_si = x[1];
    double mod_mu = x[2];
    double mod_si = x[3];
    double ll = 0;
    char reconstructedRef[512];
    char myread[512];
    char myrefe[512];
    char yourread[512];
    char yourrefe[512];
    int  yourqual[512];
    std::pair< kstring_t*, std::vector<int> >  mypair;
    kstring_t *kstr =new kstring_t;
    kstr->l=kstr->m=0;
    kstr->s=NULL;
    mypair.first = kstr;
    
    char refeBase, readBase;
    
    
    uchar * indref = NULL;
    double x_max1 = ((double)len_limit-1+0.5-anc_mu)/anc_si;
    double x_min1 = ((double)len_min-0.5-anc_mu)/anc_si;
    double x_max2 = ((double)len_limit-1+0.5-mod_mu)/mod_si;
    double x_min2 = ((double)len_min-0.5-mod_mu)/mod_si;
    if(from==-1||to==-1){
        from =0;
        to = reads->size();
    }
    //    fprintf(stderr,"from: %d to:%d\n",from,to);
    for(int i=from;i<to;i++){
        //   fprintf(stderr,"ll\t[%d]\t[%d]\t%f\tin\n",threadid,i,ll);
        bam1_t *b = (*reads)[i];
        //then we simply write it to the output
        memset(reconstructedRef,0,512);
        memset(myread,'N',512);
        memset(myrefe,'N',512);
        memset(yourread,'N',512);
        memset(yourrefe,'N',512);
        
        if (seq_ref != NULL){
            wrapperwithref(b,hdr,myread,myrefe,seq_ref);
        }else{
            reconstructRefWithPosHTS(b,mypair,reconstructedRef);
            wrapper(b,mypair.first->s,mypair.second,0,0,NULL,NULL,1,myread,myrefe);
        }
        
        if (len_limit <=0){
            len_limit = 512;
        };
        //cout << b->core.l_qseq << "Test nuc_llik\n";
        //	fprintf(stderr,"ll\t[%d]\t[%d]\t%f\tin2\t %d %d\n",threadid,i,ll,b->core.l_qseq,len_limit);
        if (b->core.l_qseq>=30 && b->core.l_qseq<len_limit){
            for (int cycle=0;cycle<b->core.l_qseq;cycle++){
                yourqual[cycle] = -1;
            }
            
            for (int cycle=0;cycle<b->core.l_qseq;cycle++){
                refeBase = refToChar[myrefe[cycle]];
                readBase = refToChar[myread[cycle]];
                size_t pos = b->core.pos+cycle;
                
                //if(refeBase!=4 && readBase!=4) {
                int dist5p=cycle;
                int dist3p=b->core.l_qseq-1-cycle;
                if( bam_is_rev(b) ){
                    refeBase=com[refeBase];
                    readBase=com[readBase];
                    dist5p=int(b->core.l_qseq)-1-cycle;
                    dist3p=cycle;
                }
                yourread[dist5p] = readBase;
                yourrefe[dist5p] = refeBase;
                yourqual[dist5p] = bam_get_qual(b)[cycle];
                //}
            }
            
            double l_anc, l_err;
            int L = b->core.l_qseq;
            l_err = ErrorLik(yourrefe, yourread, L, yourqual);
            if (!strcasecmp("b",model)){
	      l_anc = PMDLik_b(yourrefe, yourread, L, lambda, delta, delta_s, nv, yourqual,Tol); // Ancient Likelihood based on biotin model
            }else if(!strcasecmp("nb",model)){
                //cout<<"liklik "<<bam_get_qname(b)<<"\n";
	      l_anc = 0.5*PMDLik_b(yourrefe, yourread, L, lambda, delta, delta_s, nv, yourqual,Tol)+0.5*PMDLik_nb(yourrefe, yourread, L, lambda, delta, delta_s, nv, yourqual,Tol); // Ancient Likelihood based on non-biotin model
                //cout<<"l_anc "<<l_anc<<"\n";
                //cout<<"Check nuc_llik "<<PMDLik_b(yourrefe, yourread, L, lambda, delta, delta_s, nv, yourqual)<<" "<<PMDLik_nb(yourrefe, yourread, L, lambda, delta, delta_s, nv, yourqual)<<"\n";
            }else{
                fprintf(stderr,"Please specify a deamination model for further calculations.\n");
                return -1;
            }
            double y_max1 = (std::min((double)len_limit-1,std::max((double)len_min,(double)L))+0.5-anc_mu)/anc_si;
            double y_min1 = (std::max((double)len_min,std::min((double)len_limit-1,(double)L))-0.5-anc_mu)/anc_si;
            double y_max2 = (std::min((double)len_limit-1,std::max((double)len_min,(double)L))+0.5-mod_mu)/mod_si;
            double y_min2 = (std::max((double)len_min,std::min((double)len_limit-1,(double)L))-0.5-mod_mu)/mod_si;
     //cout<<"liklik "<<bam_get_qname(b)<<" "<<l_anc<<" "<<l_err<<" "<<eps<<" "<<NormalINC(y_max1, y_min1, x_max1, x_min1)<<" "<<NormalINC(y_max2, y_min2, x_max2, x_min2)<<" "<<l_anc*NormalINC(y_max1, y_min1, x_max1, x_min1)*(1-eps)+l_err*NormalINC(y_max2, y_min2, x_max2, x_min2)*eps<<" "<<log(l_anc*NormalINC(y_max1, y_min1, x_max1, x_min1)*(1-eps)+l_err*NormalINC(y_max2, y_min2, x_max2, x_min2)*eps)<<"\n";
	     ll += log(l_anc*NormalINC(y_max1, y_min1, x_max1, x_min1)*(1-eps)+l_err*NormalINC(y_max2, y_min2, x_max2, x_min2)*eps);
            //    fprintf(stderr,"ll\t[%d]\t[%d]\t%f\tyourread:%s\tyourref:%s\n",threadid,i,ll,yourread,yourrefe);
            //	    break;
            //cout<<"nuc_lliknuc_llik_check_name "<<bam_get_qname(b)<<" "<<ll<<"\n";
        }
        
    }
    return -ll;
}

double tsk_All_loglike_recalibration(const double *x, const void *dats){
    //  fprintf(stderr,"[%s]\n",__FUNCTION__);
    ncalls++;
    tsk_struct *ts = (tsk_struct *) dats;
    double check = tsk_loglike_recalibration(x, ts->reads,ts->from,ts->to,ts->hdr, ts->seq_ref, ts->len_limit, ts->len_min, ts->model, ts->eps, ts->lambda, ts->delta, ts->delta_s, ts->nv,ts->threadid,ts->Tol);
    //  fprintf(stderr,"[%s]: total like: %f\n",__FUNCTION__,check);
    return check;
}

void *tsk_All_loglike_recalibration_slave(void *dats){
    //fprintf(stderr,"[%s]\n",__FUNCTION__);
    ncalls++;
    tsk_struct *ts = (tsk_struct *) &(my_tsk_struct[(size_t) dats]);
    ts->llh_result = tsk_loglike_recalibration(ts->x, ts->reads,ts->from,ts->to,ts->hdr, ts->seq_ref, ts->len_limit, ts->len_min, ts->model, ts->eps, ts->lambda, ts->delta, ts->delta_s, ts->nv,ts->threadid,ts->Tol);
    
    pthread_exit(NULL);
}

void loglike_recalibration_grad(const double *x, double *y, char *refName,char *fname, const char* chromname, const char* bedname,int mapped_only,int se_only, int mapq, faidx_t *seq_ref,int len_limit, int len_min, char * model, double eps, double lambda, double delta, double delta_s, double nv, std::string s,double Tol){
    //fprintf(stderr,"mapped_only: %d\n",mapped_only);
    htsFormat *dingding7 =(htsFormat*) calloc(1,sizeof(htsFormat));
    double anc_mu = x[0];
    double anc_si = x[1];
    double mod_mu = x[2];
    double mod_si = x[3];
    y[0] = 0;
    y[1] = 0;
    y[2] = 0;
    y[3] = 0;
    double ll = 0;
    char reconstructedRef[512];
    char myread[512];
    char myrefe[512];
    char yourread[512];
    char yourrefe[512];
    int  yourqual[512];
    std::pair< kstring_t*, std::vector<int> >  mypair;
    kstring_t *kstr =new kstring_t;
    kstr->l=kstr->m=0;
    kstr->s=NULL;
    mypair.first = kstr;
    samFile *in=NULL;

    int nthreads = 1;
    char refeBase, readBase;
    //code below only relevant if using cram files
    if(refName!=NULL){
        char *ref =(char*) malloc(10 + strlen(refName) + 1);
        sprintf(ref, "reference=%s", refName);
        hts_opt_add((hts_opt **)&dingding7->specific,ref);
        free(ref);
    }
    if(strstr(fname,".cram")!=NULL && refName==NULL){
        fprintf(stderr,"\t-> cram file requires reference with -T FILE.fa \n");
        exit(0);
    }
    if((in=sam_open_format(fname,"r",dingding7))==NULL ){
        fprintf(stderr,"[%s] nonexistant file: %s\n",__FUNCTION__,fname);
        exit(0);
    }
    bam_hdr_t  *hdr = sam_hdr_read(in);
    bam1_t *b = bam_init1();
    
    int chrom_num = hdr->n_targets;
    size_t * chrom_line = (size_t*)malloc(chrom_num*(sizeof(size_t)));
    std::vector<std::array<size_t, 2> > bedsites;
    
    // loading the bedfile
    if (bedname!=NULL){
        //fprintf(stderr,"Loading the bedfile %s ...\n",bedname);
        BGZF *fp = NULL;
        fp = bgzf_open(bedname,"rb");
        
        kstring_t *kstr1 = new kstring_t;
        kstr1->s = NULL;
        kstr1->l = kstr1->m = 0;
        int line=0;
	std::string word0, word;
        bgzf_getline(fp,'\n',kstr1);
        for (int j=0; j<chrom_num; j++){
            chrom_line[j] = line;
            do{
                //while(bgzf_getline(fp,'\n',kstr1)>0){
                std::istringstream iss(kstr1->s);
                //string word0, word;
                getline(iss,word0,'\t');
                size_t num;
                if  (strcmp(word0.c_str(),hdr->target_name[j])==0){
                    line++;
                }
                //size_t bedsite[2];
                std::array<size_t, 2> bedsite;
                int idx = 1;
                while (getline(iss,word,'\t')){
                    //sscanf(word.c_str(), "%zu", &num);
                    if (strcmp(word0.c_str(),hdr->target_name[j])==0){
                        sscanf(word.c_str(), "%zu", &num);
                        idx = 1-idx;
                        bedsite[idx] = num;
                        if (idx == 1){
                            bedsites.push_back(bedsite);
                        }
                    }
                }
                //cout << "\n";
            }while (strcmp(word0.c_str(),hdr->target_name[j])==0 && bgzf_getline(fp,'\n',kstr1)>0);
        }
    }
    
    
    int ret;
    int refId=-1;
    double num = 0.0;
    size_t max_site;
    nproc1 = 0;
    int chromId=-1;
    if (chromname!=NULL){
        for (int j=0; j<chrom_num; j++){
            if (!strcmp(chromname,hdr->target_name[j])){
                chromId = j;
                //fprintf(stderr,"We will focus on Chromosome %s!\n", hdr->target_name[j]);
            }
        }
    }
    if (chromId==-1){
        //fprintf(stderr,"No meaningful chromosome name is provided, therefore we will focus on all provided chromosomes!\n");
    }
    
    uchar * indref = NULL;
    double x_max1 = ((double)len_limit-1+0.5-anc_mu)/anc_si;
    double x_min1 = ((double)len_min-0.5-anc_mu)/anc_si;
    double x_max2 = ((double)len_limit-1+0.5-mod_mu)/mod_si;
    double x_min2 = ((double)len_min-0.5-mod_mu)/mod_si;
    while(((ret=sam_read1(in,hdr,b)))>0) {
        //nproc1++;
        //if -m was used, discard unmapped reads
        if(mapped_only!=0){
            if(b->core.flag&4)
                continue;
        }
        
        //default -a 1
        //only use single end reads
        //which is either single end or collapsed reads
        //#if 0
        if(se_only==1){
            if(b->core.flag&1)//if paired continue
                continue;
        }
        //#endif
        if(b->core.flag>256)
            continue; //Remove possible duplicates
        
        // if mapq threshold is set
        if(mapq!=-1 && b->core.qual<mapq)
            continue;
        if (chromId==-1||chromId==b->core.tid){
            nproc1++;
        }else if (chromId!=b->core.tid){
            continue;
        }
        if(refId==-1||refId!=b->core.tid){
            refId=b->core.tid;
            //fprintf(stderr,"\t-> Now at Chromosome: %s\n",hdr->target_name[refId]);
            if (bedname!=NULL){
                if (indref!=NULL){
                    free(indref);
                }
                //checkchrome=b->core.tid;
                max_site = 0;
                size_t max_line =  refId < chrom_num-1 ? chrom_line[refId+1] : bedsites.size();
                for (size_t l=chrom_line[refId]; l < max_line; l++){
                    if (max_site < bedsites[l][1]){
                        max_site = bedsites[l][1]; //Assume the provided bed file is 1-based (but will treat it as 0-based internally)
                    }
                }
                size_t chr_len = max_site; // The maximum position considered in the bed file along this chromosome
                
                //faidx_seq_len(seq_ref,hdr->target_name[refId]);
                indref = (uchar*) malloc(chr_len*sizeof(uchar));
                for (size_t l=0; l<chr_len; l++){
                    indref[l] = (uchar)0;
                }
                
                for (size_t l=chrom_line[refId]; l < max_line; l++){
                    for (size_t i = bedsites[l][0]; i <= bedsites[l][1]; i++){
                        indref[i-1] = (uchar)1; // Shift by 1 so that it is treated 0 based internally.
                    }
                }
            }
        }
        
        //then we simply write it to the output
        memset(reconstructedRef,0,512);
        memset(myread,'N',512);
        memset(myrefe,'N',512);
        
        if (seq_ref != NULL){
            wrapperwithref(b,hdr,myread,myrefe,seq_ref);
        }else{
            reconstructRefWithPosHTS(b,mypair,reconstructedRef);
            wrapper(b,mypair.first->s,mypair.second,0,0,NULL,NULL,1,myread,myrefe);
        }
        
        if (len_limit <=0){
            len_limit = 512;
        };
        //cout << b->core.l_qseq << "Test nuc_llik\n";
        if (b->core.l_qseq>=30 && b->core.l_qseq<len_limit){
            for (int cycle=0;cycle<b->core.l_qseq;cycle++){
                yourqual[cycle] = -1;
            }
            //cout << b->core.l_qseq << "\n";
            int isbedfrag = 0;
            for (int cycle=0;cycle<b->core.l_qseq;cycle++){
                refeBase = refToChar[myrefe[cycle]];
                readBase = refToChar[myread[cycle]];
                size_t pos = b->core.pos+cycle;
                //                cout << (int)refeBase << " " << (int)readBase<<" "<<pos<<" "<<max_site<<"\n";
                if(refeBase!=4 && readBase!=4 && bedname!=NULL &&  pos < max_site && indref[pos]==1){ //check whether the fragment has intersection with the bed file
                    isbedfrag += 1;
                }else if(refeBase!=4 && readBase!=4 && bedname==NULL){
                    isbedfrag += 1;
                }
            }
            if (isbedfrag > 0){
                for (int cycle=0;cycle<b->core.l_qseq;cycle++){
                    int dist5p=cycle;
                    int dist3p=b->core.l_qseq-1-cycle;
                    //cout<<"flag "<<b->core.flag<<"\n";
                    if( bam_is_rev(b) ){
                        // cout<<"rev "<<"\t";
                        refeBase=com[refeBase];
                        readBase=com[readBase];
                        //dist5p=int(al.QueryBases.size())-1-i;
                        dist5p=int(b->core.l_qseq)-1-cycle;
                        dist3p=cycle;
                    }
                    yourread[dist5p] = readBase;
                    yourrefe[dist5p] = refeBase;
                    yourqual[dist5p] = bam_get_qual(b)[cycle];
                }
                double l_anc, l_err;
                int L = b->core.l_qseq;
                l_err = ErrorLik(yourrefe, yourread, L, yourqual);
                if (!strcasecmp("b",model)){
		  l_anc = PMDLik_b(yourrefe, yourread, L, lambda, delta, delta_s, nv, yourqual,Tol); // Ancient Likelihood based on biotin model
                }else if(!strcasecmp("nb",model)){
		  l_anc = 0.5*PMDLik_b(yourrefe, yourread, L, lambda, delta, delta_s, nv, yourqual,Tol)+0.5*PMDLik_nb(yourrefe, yourread, L, lambda, delta, delta_s, nv, yourqual,Tol); // Ancient Likelihood based on non-biotin model
                }else{
                    fprintf(stderr,"Please specify a deamination model for further calculations.\n");
                    exit;
                }
                double y_max1 = (std::min((double)len_limit-1,std::max((double)len_min,(double)L))+0.5-anc_mu)/anc_si;
                double y_min1 = (std::max((double)len_min,std::min((double)len_limit-1,(double)L))-0.5-anc_mu)/anc_si;
                double y_max2 = (std::min((double)len_limit-1,std::max((double)len_min,(double)L))+0.5-mod_mu)/mod_si;
                double y_min2 = (std::max((double)len_min,std::min((double)len_limit-1,(double)L))-0.5-mod_mu)/mod_si;
                //ll += log(l_anc*normalINC(y_max1, y_min1, x_max1, x_min1)*(1-eps)+l_err*normalINC(y_max2, y_min2, x_max2, x_min2)*eps);
                y[0] -= l_anc*(1-eps)*NormalINC_grad_mu(y_max1, y_min1, x_max1, x_min1, anc_mu, anc_si)/(l_anc*NormalINC(y_max1, y_min1, x_max1, x_min1)*(1-eps)+l_err*NormalINC(y_max2, y_min2, x_max2, x_min2)*eps);
                y[1] -= l_anc*(1-eps)*NormalINC_grad_si(y_max1, y_min1, x_max1, x_min1, anc_mu, anc_si)/(l_anc*NormalINC(y_max1, y_min1, x_max1, x_min1)*(1-eps)+l_err*NormalINC(y_max2, y_min2, x_max2, x_min2)*eps);
                y[2] -= l_err*eps*NormalINC_grad_mu(y_max2, y_min2, x_max2, x_min2, mod_mu, mod_si)/(l_anc*NormalINC(y_max1, y_min1, x_max1, x_min1)*(1-eps)+l_err*NormalINC(y_max2, y_min2, x_max2, x_min2)*eps);
                y[3] -= l_err*eps*NormalINC_grad_si(y_max2, y_min2, x_max2, x_min2, mod_mu, mod_si)/(l_anc*NormalINC(y_max1, y_min1, x_max1, x_min1)*(1-eps)+l_err*NormalINC(y_max2, y_min2, x_max2, x_min2)*eps);
                //double PostAncProb = AncProb(yourrefe, yourread, b->core.l_qseq, lambda, delta, delta_s, nv, yourqual, model ,eps);
                //double PostPMDProb = PMDProb(yourrefe, yourread, b->core.l_qseq, lambda, delta, delta_s, nv, yourqual, model);
                //cout<<"\t"<<"AN: "<<PostAncProb<<" \t"<<"PD: "<<PostPMDProb<<"\n";
                num += 1.0; //number of intersected reads.
            }
        }
        //    fprintf(stderr,"\nmyread:\n%.*s\nmyReference:\n%.*s\n",b->core.l_qseq,myread,b->core.l_qseq,myrefe);
        //    fprintf(stderr,"---read[%d]----\n",nproc1);
    }
    //free(ofname);
    //assert(sam_close(out)==0);
    if (indref!=NULL){
        free(indref);
    }
    //   for(int i=0;i<queue->m;i++)
    //       bam_destroy1(queue->d[i]);
    //   free(queue->d);
    //   free(queue);
    bedsites.clear();
    bedsites.shrink_to_fit();
    bam_destroy1(b);
    
    hts_opt_free((hts_opt *)dingding7->specific);
    free(dingding7);
    sam_hdr_destroy(hdr);
    assert(sam_close(in)==0);

}

void tsk_loglike_recalibration_grad(const double *x, double *y, std::vector<bam1_t*> *reads,int from, int to,sam_hdr_t *hdr, faidx_t *seq_ref,int len_limit, int len_min, char * model, double eps, double lambda, double delta, double delta_s, double nv,double Tol){
    //fprintf(stderr,"[%s]\n",__FUNCTION__);
    htsFormat *dingding7 =(htsFormat*) calloc(1,sizeof(htsFormat));
    double anc_mu = x[0];
    double anc_si = x[1];
    double mod_mu = x[2];
    double mod_si = x[3];
    y[0] = 0;
    y[1] = 0;
    y[2] = 0;
    y[3] = 0;
    double ll = 0;
    char reconstructedRef[512];
    char myread[512];
    char myrefe[512];
    char yourread[512];
    char yourrefe[512];
    int  yourqual[512];
    std::pair< kstring_t*, std::vector<int> >  mypair;
    kstring_t *kstr =new kstring_t;
    kstr->l=kstr->m=0;
    kstr->s=NULL;
    mypair.first = kstr;
    char refeBase, readBase;
    uchar * indref = NULL;
    double x_max1 = ((double)len_limit-1+0.5-anc_mu)/anc_si;
    double x_min1 = ((double)len_min-0.5-anc_mu)/anc_si;
    double x_max2 = ((double)len_limit-1+0.5-mod_mu)/mod_si;
    double x_min2 = ((double)len_min-0.5-mod_mu)/mod_si;
    if(1||from==-1||to==-1){//fix this uncomment out.
        from = 0;
        to = reads->size();
    }
    //    fprintf(stderr,"from: %d to:%d\n",from,to):qw
    for(int i=from;i<to;i++){
        bam1_t *b = (*reads)[i];
        
        memset(reconstructedRef,0,512);
        memset(myread,'N',512);
        memset(myrefe,'N',512);
        
        if (seq_ref != NULL){
            wrapperwithref(b,hdr,myread,myrefe,seq_ref);
        }else{
            reconstructRefWithPosHTS(b,mypair,reconstructedRef);
            wrapper(b,mypair.first->s,mypair.second,0,0,NULL,NULL,1,myread,myrefe);
        }
        
        if (len_limit <=0){
            len_limit = 512;
        };
        
        if (b->core.l_qseq>=30 && b->core.l_qseq<len_limit) {
            
            for (int cycle=0;cycle<b->core.l_qseq;cycle++)
                yourqual[cycle] = -1;
            
            for (int cycle=0;cycle<b->core.l_qseq;cycle++){
                refeBase = refToChar[myrefe[cycle]];
                readBase = refToChar[myread[cycle]];
                
                //                if(refeBase!=4 && readBase!=4){
                int dist5p=cycle;
                int dist3p=b->core.l_qseq-1-cycle;
                if( bam_is_rev(b) ){
                    refeBase=com[refeBase];
                    readBase=com[readBase];
                    
                    dist5p=int(b->core.l_qseq)-1-cycle;
                    dist3p=cycle;
                }
                yourread[dist5p] = readBase;
                yourrefe[dist5p] = refeBase;
                yourqual[dist5p] = bam_get_qual(b)[cycle];

            }
            
            double l_anc, l_err;
            int L = b->core.l_qseq;
            l_err = ErrorLik(yourrefe, yourread, L, yourqual);
            if (!strcasecmp("b",model)){
	      l_anc = PMDLik_b(yourrefe, yourread, L, lambda, delta, delta_s, nv, yourqual,Tol); // Ancient Likelihood based on biotin model
            }else if(!strcasecmp("nb",model)){
	      l_anc = 0.5*PMDLik_b(yourrefe, yourread, L, lambda, delta, delta_s, nv, yourqual,Tol)+0.5*PMDLik_nb(yourrefe, yourread, L, lambda, delta, delta_s, nv, yourqual,Tol); // Ancient Likelihood based on non-biotin model
            }else{
                fprintf(stderr,"Please specify a deamination model for further calculations.\n");
                exit;
            }
            double y_max1 = (std::min((double)len_limit-1,std::max((double)len_min,(double)L))+0.5-anc_mu)/anc_si;
            double y_min1 = (std::max((double)len_min,std::min((double)len_limit-1,(double)L))-0.5-anc_mu)/anc_si;
            double y_max2 = (std::min((double)len_limit-1,std::max((double)len_min,(double)L))+0.5-mod_mu)/mod_si;
            double y_min2 = (std::max((double)len_min,std::min((double)len_limit-1,(double)L))-0.5-mod_mu)/mod_si;
            
            //cout<<"nuc_lliktest y_max1 "<<y_max1<<" y_min1 "<<y_min1<<" "<<NormalINC_grad_mu(y_max1, y_min1, x_max1, x_min1, anc_mu, anc_si)<<"\n";
            //ll += log(l_anc*normalINC(y_max1, y_min1, x_max1, x_min1)*(1-eps)+l_err*normalINC(y_max2, y_min2, x_max2, x_min2)*eps)
	    y[0] -= l_anc*(1-eps)*NormalINC_grad_mu(y_max1, y_min1, x_max1, x_min1, anc_mu, anc_si)/(l_anc*NormalINC(y_max1, y_min1, x_max1, x_min1)*(1-eps)+l_err*NormalINC(y_max2, y_min2, x_max2, x_min2)*eps);
            y[1] -= l_anc*(1-eps)*NormalINC_grad_si(y_max1, y_min1, x_max1, x_min1, anc_mu, anc_si)/(l_anc*NormalINC(y_max1, y_min1, x_max1, x_min1)*(1-eps)+l_err*NormalINC(y_max2, y_min2, x_max2, x_min2)*eps);
            y[2] -= l_err*eps*NormalINC_grad_mu(y_max2, y_min2, x_max2, x_min2, mod_mu, mod_si)/(l_anc*NormalINC(y_max1, y_min1, x_max1, x_min1)*(1-eps)+l_err*NormalINC(y_max2, y_min2, x_max2, x_min2)*eps);
            y[3] -= l_err*eps*NormalINC_grad_si(y_max2, y_min2, x_max2, x_min2, mod_mu, mod_si)/(l_anc*NormalINC(y_max1, y_min1, x_max1, x_min1)*(1-eps)+l_err*NormalINC(y_max2, y_min2, x_max2, x_min2)*eps);
            //double PostAncProb = AncProb(yourrefe, yourread, b->core.l_qseq, lambda, delta, delta_s, nv, yourqual, model ,eps);
            //double PostPMDProb = PMDProb(yourrefe, yourread, b->core.l_qseq, lambda, delta, delta_s, nv, yourqual, model);
            //cout<<"\t"<<"AN: "<<PostAncProb<<" \t"<<"PD: "<<PostPMDProb<<"\n";
            
            
        }
        //    fprintf(stderr,"\nmyread:\n%.*s\nmyReference:\n%.*s\n",b->core.l_qseq,myread,b->core.l_qseq,myrefe);
        //    fprintf(stderr,"---read[%d]----\n",nproc1);
    }
}

void tsk_All_loglike_recalibration_grad(const double *x, double *y,const void *dats){
    ncalls_grad++;
    //tsk_struct *ts = &my_tsk_struct[0];//assuming single thread
    tsk_struct *ts = (tsk_struct *) dats;//assuming multiple threads
    tsk_loglike_recalibration_grad(x,y,ts->reads,ts->from,ts->to,ts->hdr,ts->seq_ref,ts->len_limit, ts->len_min, ts->model,ts->eps,ts->lambda,ts->delta,ts->delta_s,ts->nv,ts->Tol);
    //fprintf(stderr,"[%s]\n",__FUNCTION__);
}

// Check here
void *tsk_All_loglike_recalibration_grad_slave(void *dats){
    //fprintf(stderr,"[%s]\n",__FUNCTION__);
    ncalls_grad++;
    tsk_struct *ts = (tsk_struct *) &(my_tsk_struct[(size_t) dats]);
    //ts->llh_grad_result;
    tsk_loglike_recalibration_grad(ts->x,ts->llh_result_grad,ts->reads,ts->from,ts->to,ts->hdr,ts->seq_ref,ts->len_limit, ts->len_min,ts->model,ts->eps,ts->lambda,ts->delta,ts->delta_s,ts->nv,ts->Tol);
    //ts->llh_result = tsk_loglike_recalibration(ts->x, ts->reads,ts->from,ts->to,ts->hdr, ts->seq_ref, ts->len_limit, ts->model, ts->eps, ts->lambda, ts->delta, ts->delta_s, ts->nv,ts->threadid);
    pthread_exit(NULL);
}

void tsk_loglike_recalibration_hess(const double *x, double **y, std::vector<bam1_t*> *reads,int from, int to,sam_hdr_t *hdr, faidx_t *seq_ref,int len_limit, int len_min, char * model, double eps, double lambda, double delta, double delta_s, double nv,double Tol){
    //fprintf(stderr,"[%s]\n",__FUNCTION__);
    htsFormat *dingding7 =(htsFormat*) calloc(1,sizeof(htsFormat));
    double anc_mu = x[0];
    double anc_si = x[1];
    double mod_mu = x[2];
    double mod_si = x[3];
    for(int i=0; i<3; i++){
        for(int j=0; j<3; j++){
            y[i][j] = 0;
        }
    }
    
    char reconstructedRef[512];
    char myread[512];
    char myrefe[512];
    char yourread[512];
    char yourrefe[512];
    int  yourqual[512];
    std::pair< kstring_t*, std::vector<int> >  mypair;
    kstring_t *kstr =new kstring_t;
    kstr->l=kstr->m=0;
    kstr->s=NULL;
    mypair.first = kstr;
    char refeBase, readBase;
    uchar * indref = NULL;
    double x_max1 = ((double)len_limit-1+0.5-anc_mu)/anc_si;
    double x_min1 = ((double)len_min-0.5-anc_mu)/anc_si;
    double x_max2 = ((double)len_limit-1+0.5-mod_mu)/mod_si;
    double x_min2 = ((double)len_min-0.5-mod_mu)/mod_si;
    if(1||from==-1||to==-1){//fix this uncomment out.
        from = 0;
        to = reads->size();
    }
    //    fprintf(stderr,"from: %d to:%d\n",from,to)
    for(int i=from;i<to;i++){
        bam1_t *b = (*reads)[i];
        
        memset(reconstructedRef,0,512);
        memset(myread,'N',512);
        memset(myrefe,'N',512);
        
        if (seq_ref != NULL){
            wrapperwithref(b,hdr,myread,myrefe,seq_ref);
        }else{
            reconstructRefWithPosHTS(b,mypair,reconstructedRef);
            wrapper(b,mypair.first->s,mypair.second,0,0,NULL,NULL,1,myread,myrefe);
        }
        
        if (len_limit <=0){
            len_limit = 512;
        };
        
        if (b->core.l_qseq>=30 && b->core.l_qseq<len_limit) {
            
            for (int cycle=0;cycle<b->core.l_qseq;cycle++)
                yourqual[cycle] = -1;
            
            for (int cycle=0;cycle<b->core.l_qseq;cycle++){
                refeBase = refToChar[myrefe[cycle]];
                readBase = refToChar[myread[cycle]];
                
                //if(refeBase!=4 && readBase!=4){
                int dist5p=cycle;
                int dist3p=b->core.l_qseq-1-cycle;
                if( bam_is_rev(b) ){
                    refeBase=com[refeBase];
                    readBase=com[readBase];
                    
                    dist5p=int(b->core.l_qseq)-1-cycle;
                    dist3p=cycle;
                }
                yourread[dist5p] = readBase;
                yourrefe[dist5p] = refeBase;
                yourqual[dist5p] = bam_get_qual(b)[cycle];
                //}
            }
            
            double l_anc, l_err;
            int L = b->core.l_qseq;
            l_err = ErrorLik(yourrefe, yourread, L, yourqual);
            if (!strcasecmp("b",model)){
	      l_anc = PMDLik_b(yourrefe, yourread, L, lambda, delta, delta_s, nv, yourqual,Tol); // Ancient Likelihood based on biotin model
            }else if(!strcasecmp("nb",model)){
	      l_anc = 0.5*PMDLik_b(yourrefe, yourread, L, lambda, delta, delta_s, nv, yourqual,Tol)+0.5*PMDLik_nb(yourrefe, yourread, L, lambda, delta, delta_s, nv, yourqual,Tol); // Ancient Likelihood based on non-biotin model
            }else{
                fprintf(stderr,"Please specify a deamination model for further calculations.\n");
                exit;
            }
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
        }
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

void *tsk_All_loglike_recalibration_hess_slave(void *dats){
    tsk_struct *ts = (tsk_struct *) &(my_tsk_struct[(size_t) dats]);
    tsk_loglike_recalibration_hess(ts->x,ts->llh_result_hess,ts->reads,ts->from,ts->to,ts->hdr,ts->seq_ref,ts->len_limit, ts->len_min,ts->model,ts->eps,ts->lambda,ts->delta,ts->delta_s,ts->nv,ts->Tol);
    pthread_exit(NULL);
}
