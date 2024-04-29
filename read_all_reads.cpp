#include <htslib/bgzf.h>

#include "read_all_reads.h"

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
