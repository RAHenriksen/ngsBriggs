#include <htslib/bgzf.h>

#include "read_all_reads.h"
//#include <iostream>

sam_hdr_t *read_all_reads(const char *htsname, const char *bedfile,const char *refName,std::vector<bam1_t*> &ret){
    fprintf(stderr,"\t-> [%s] htsname: %s bedfile: %s\n",__FUNCTION__,htsname,bedfile);
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
    fprintf(stderr,"\t-> Parsing possible bedfile: %s\n",bedfile);
    hts_idx_t *idx = NULL;
    hts_itr_t *itr = NULL;
    if(bedfile!=NULL){
        idx = sam_index_load(in,htsname);
        if (idx == NULL) { // index is unavailable
            fprintf(stderr, "\t[%s ] random alignment retrieval only works for indexed BAM or CRAM files. (file: \'%s\' )\n",__FUNCTION__,htsname);
            exit(0);
        }
    }
    
    
    if(bedfile==NULL){
        bam1_t *b = bam_init1();
        while(sam_read1(in,hdr,b)){
            ret.push_back(b);
            b= bam_init1();
        }
        bam_destroy1(b);
    }else{
        kstring_t kstr;kstr.s=NULL;kstr.l=kstr.m = 0;
        BGZF *bg = bgzf_open(bedfile,"rb");
        while(bgzf_getline(bg,'\n',&kstr)>0){
            char *tok = strtok(kstr.s,"\t\n");
            char *start = strtok(NULL,"\t\n");
            char *stop = strtok(NULL,"\t\n");
            char tmp[1024];
            snprintf(tmp,1024,"%s:%s-%s",tok,start,stop);
            //fprintf(stderr,"tmp: %s\n",tmp);
            if(itr)
                hts_itr_destroy(itr);
            itr = sam_itr_querys(idx, hdr, tmp);
            if (itr == NULL) { // reference name is not found
                fprintf(stderr, "[main_samview] region \"%s\" specifies an unknown reference name. Continue anyway.\n",tmp);
                exit(0);
            }
            
            bam1_t *b = bam_init1();
            //fprintf(stderr,"LeiLeiCheck %lld\n",b->id);
            while(sam_itr_next(in,itr,b)>0){
                int j = 0;
                for (int i = 0; i<ret.size(); i++){
                    if (!strcmp(bam_get_qname(ret[i]),bam_get_qname(b))){
                        j=1;
                        break;
                    }
                }
                //if (find(ret.begin(), ret.end(), b) != ret.end())
                //if (b->id == 0){
                if (j==0)
                {   ret.push_back(b);
                    //std::cout<<b<<"\t"<<bam_get_qname(b)<<"\n";}
                    //b->id = 1;
                }
                //fprintf(stderr,"LeiLeiCheck %lld\n",b->l_data);
                b= bam_init1();
                //fprintf(stderr,"LeiLeiCheck %lld\n",b->l_data);
            }
            bam_destroy1(b);
            
        }
        
    }
    
    
    fprintf(stderr,"ret.size(): %lu\n",ret.size());
    return hdr;
}
