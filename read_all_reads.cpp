#include <htslib/bgzf.h>

#include "read_all_reads.h"
#include "profile.h"
#include "misc.h"
#include "ngsBriggs.h"

char *talloc(char *ori,int len){
  char *des = (char*) malloc(len);
  des =(char*) memcpy(des,ori,len);
  return des;
}

sam_hdr_t *read_all_reads(const char *htsname,const char *refName,std::vector<bam1_t*> &ret){
  htsFormat myHtsFormat;
    samFile *in=NULL;
    if(refName!=NULL){
        char *ref =(char*) malloc(10 + strlen(refName) + 1);
        snprintf(ref,10 + strlen(refName) + 1, "reference=%s", refName);
        hts_opt_add((hts_opt **)&myHtsFormat.specific,ref);
        free(ref);
    }
    if(strstr(htsname,".cram")!=NULL && refName==NULL){
        fprintf(stderr,"\t-> cram file requires reference with -T FILE.fa \n");
        exit(0);
    }
    if((in=sam_open(htsname,"r"))==NULL ){
        fprintf(stderr,"[%s] nonexistant file: %s\n",__FUNCTION__,htsname);
        exit(0);
    }
    bam_hdr_t  *hdr = sam_hdr_read(in);
    
    bam1_t *b = bam_init1();
    while(sam_read1(in,hdr,b)>=0){
      ret.push_back(b);
      b= bam_init1();
    }
    bam_destroy1(b);

    fprintf(stderr,"[%s] ret.size(): %lu\n",__FUNCTION__,ret.size());
    return hdr;
}

sam_hdr_t *read_all_reads(const char *htsname,const char *refName,faidx_t *seq_ref,std::vector<asite> &ret,int len_limit){
  //  fprintf(stderr,"len_limit: %d\n",len_limit);exit(0);
  htsFormat myHtsFormat;
  samFile *in=NULL;
  if(refName!=NULL){
    char *ref =(char*) malloc(10 + strlen(refName) + 1);
    snprintf(ref,10 + strlen(refName) + 1, "reference=%s", refName);
    hts_opt_add((hts_opt **)&myHtsFormat.specific,ref);
    free(ref);
  }
  if(strstr(htsname,".cram")!=NULL && refName==NULL){
    fprintf(stderr,"\t-> cram file requires reference with -T FILE.fa \n");
    exit(0);
  }
  if((in=sam_open(htsname,"r"))==NULL ){
    fprintf(stderr,"[%s] nonexistant file: %s\n",__FUNCTION__,htsname);
    exit(0);
  }
    bam_hdr_t  *hdr = sam_hdr_read(in);
    char reconstructedRef[512];
    char myread[512];
    char myrefe[512];
    char yourread[512];
    char yourrefe[512];
    uint8_t  yourqual[512];
    std::pair< kstring_t*, std::vector<int> >  mypair;
    kstring_t *kstr =new kstring_t;
    kstr->l=kstr->m=0;
    kstr->s=NULL;
    mypair.first = kstr;
    char refeBase, readBase;
    if (len_limit <=0){
      len_limit = 512;
    };
    bam1_t *b = bam_init1();
    while(sam_read1(in,hdr,b)>=0) {
      memset(reconstructedRef,0,512);
      memset(myread,'N',512);
      memset(myrefe,'N',512);
      memset(yourread,'N',512);
      memset(yourrefe,'N',512);
      for (int cycle=0;cycle<b->core.l_qseq;cycle++)
	yourqual[cycle] = 254;
      
      if (seq_ref != NULL){
	wrapperwithref(b,hdr,myread,myrefe,seq_ref);
      }else{
	reconstructRefWithPosHTS(b,mypair,reconstructedRef);
	wrapper(b,mypair.first->s,mypair.second,0,0,NULL,NULL,1,myread,myrefe);
      }
      
      if (b->core.l_qseq<30 || b->core.l_qseq>=len_limit){
	//	fprintf(stderr,"SKIPPING\n");
	continue;
      }
      
      for (int cycle=0;cycle<b->core.l_qseq;cycle++){
	refeBase = refToChar[myrefe[cycle]];
	readBase = refToChar[myread[cycle]];
	if(0&&cycle<8){
	  fprintf(stderr,"raw %c %c raw2 %c %c int %d %d\n",myrefe[cycle],myread[cycle],refeBase,readBase ,(int)refeBase,(int)readBase );
	}
	size_t pos = b->core.pos+cycle;
        
	//if(refeBase!=4 && readBase!=4) {
	int dist5p=cycle;
	int dist3p=b->core.l_qseq-1-cycle;
	if( bam_is_rev(b) ){
	  //	    fprintf(stderr,"comp:\n");
	  refeBase=com[refeBase];
	  readBase=com[readBase];
	  dist5p=int(b->core.l_qseq)-1-cycle;
	  dist3p=cycle;
	}
	yourread[dist5p] = readBase;
	yourrefe[dist5p] = refeBase;
	yourqual[dist5p] =(unsigned int) bam_get_qual(b)[cycle];
	//	  fprintf(stderr,"yourqual: %d \n",(int)yourqual[dist5p]);
	if(0&&cycle<8)
	  fprintf(stderr,"rrr: %d %d \n",yourread[dist5p] ,yourrefe[dist5p]);
	//}
      }
      asite as={NULL,NULL,NULL,0};
      as.read =talloc(yourread,b->core.l_qseq);
      as.ref = talloc(yourrefe,b->core.l_qseq);
      as.qual = new uint8_t[b->core.l_qseq];
      as.qual =(uint8_t *) memcpy(as.qual,yourqual,b->core.l_qseq);
      as.len = b->core.l_qseq;
      ret.push_back(as);
    }
    bam_destroy1(b);
    fprintf(stderr,"[%s] ret.size(): %lu\n",__FUNCTION__,ret.size());
    return hdr;
}


