//see header for description

#include <cassert>

#include "likelihood.h"
#include "read_all_reads.h"
#include "profile.h"
//#include "misc.h"


double **read_all_reads(samFile *in,sam_hdr_t *hdr,faidx_t *seq_ref,int len_limit,double lambda,double delta,double delta_s,double nv,double Tol,int &ndim, int model,int l_check){
 
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
    std::vector<double> vec;
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
      double l_err = ErrorLik(yourrefe, yourread, b->core.l_qseq, yourqual,l_check);
      double l_anc = PMDLik_b(yourrefe, yourread, b->core.l_qseq, lambda, delta, delta_s, nv, yourqual,Tol,l_check);
      //  fprintf(stderr,"pmdlik_b: %f nb: %f model: %d\n",PMDLik_b(yourrefe, yourread, b->core.l_qseq, lambda, delta, delta_s, nv, yourqual,Tol),PMDLik_nb(yourrefe, yourread, b->core.l_qseq, lambda, delta, delta_s, nv, yourqual,Tol),model);

      if(model==1)
	l_anc = 0.5*l_anc + 0.5*PMDLik_nb(yourrefe, yourread, b->core.l_qseq, lambda, delta, delta_s, nv, yourqual,Tol,l_check);
      //   fprintf(stderr,"lanc: %f\n",l_anc);exit(0);
      vec.push_back(l_err);
      vec.push_back(l_anc);
      vec.push_back(b->core.l_qseq);
    }
    bam_destroy1(b);
    fprintf(stderr,"[%s] vec.size(): %lu\n",__FUNCTION__,vec.size());
    assert((vec.size() % 3)==0);//should be divisible by 3
    double **ret = new double*[3];
    ret[0] = new double[vec.size()/3];
    ret[1] = new double[vec.size()/3];
    ret[2] = new double[vec.size()/3];
    int at =0;
    for(size_t i=0;i<vec.size();i+=3){
      ret[0][at] = vec[i];//l_err
      ret[1][at] = vec[i+1];//l_anc
      ret[2][at++] = vec[i+2];//l_anc
    }
    ndim = vec.size()/3;
    return ret;
}


