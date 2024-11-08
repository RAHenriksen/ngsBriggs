#pragma once
#include <vector>
#include <map>

#include <htslib/sam.h>
#include <htslib/faidx.h>
#include <htslib/kstring.h>


#define bam_is_sec(b)         (((b)->core.flag&BAM_FSECONDARY)        != 0)
#define bam_is_supp(b)        (((b)->core.flag&BAM_FSUPPLEMENTARY)    != 0)
#define bam_is_paired(b)      (((b)->core.flag&BAM_FPAIRED)     != 0)
#define bam_is_rmdup(b)       (((b)->core.flag&BAM_FDUP)        != 0)
#define bam_is_qcfailed(b)    (((b)->core.flag&BAM_FQCFAIL)     != 0)
#define bam_is_unmapped(b)    (((b)->core.flag&BAM_FUNMAP)      != 0)
#define bam_is_read1(b)       (((b)->core.flag&BAM_FREAD1)      != 0)
#define bam_is_failed(b)      ( bam_is_qcfailed(b) || bam_is_rmdup(b) || bam_is_supp(b) )

typedef struct{
  size_t nreads;
  float **mm5pF;
  float **mm3pF;
}triple;


extern char toIndex[4][4];


class damage{
  char *reconstructedTemp;
public:
  int nthreads;
  int MAXLENGTH;
  int minQualBase;//currently not set; should be set in init_damage
  int nclass;
  float **mm5pF;//will point to first. Just a hack
  float **mm3pF;//will point to first. Just a hack
  std::pair< kstring_t*, std::vector<int> >  reconstructedReference;
  std::map<int,triple > assoc;
  void write(char *prefix,bam_hdr_t *hdr);
  void bwrite(char *prefix,bam_hdr_t *hdr);
  int damage_analysis( bam1_t *b,int whichclass,float incval);
  void printit(FILE *fp,int l);
  int temp_len;
  damage(int maxlen,int nthd,int minqb){
    temp_len = 512;
    MAXLENGTH = maxlen;
    minQualBase = minqb;
    nthreads = nthd;
    reconstructedTemp=(char*)calloc(temp_len,1);
    kstring_t *kstr =new kstring_t;
    kstr->l=kstr->m=0;
    kstr->s=NULL;
    reconstructedReference.first = kstr;
    mm5pF=mm3pF=NULL;
  }
  ~damage(){free(reconstructedTemp);}
};

void destroy_damage(damage *dmg);

std::map<int,double *> load_bdamage3(const char *fname,int howmany);

typedef struct{
  int nreads;//this is nalignements
  double *fwD;
  double *bwD;
}mydataD;

typedef struct{
  int nreads;//this is nalignements
  double *data;
}mydata2;

void wrapper(const bam1_t *b,const char * reconstructedReference,const std::vector<int> &  reconstructedReferencePos,const int & minQualBase, int MAXLENGTH,float **mm5p,float **mm3p,float incval,char myread[512],char myref[512]);
std::map<int, mydataD> load_bdamage_full(const char* fname,int &printlength);
std::map<int, mydata2> load_lcastat(const char* fname);
void  reconstructRefWithPosHTS(const bam1_t   * b,std::pair< kstring_t *, std::vector<int> > &pp,char *reconstructedTemp);
extern char refToChar[256];
extern char com[5];
void wrapperwithref(const bam1_t   * b,const bam_hdr_t  *hdr, char myread[512], char myref[512],faidx_t *seq_ref);
