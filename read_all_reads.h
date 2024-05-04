#include <vector>
#include <htslib/sam.h>
#include <htslib/faidx.h>

sam_hdr_t * read_all_reads(const char *htsname,const char *refName,std::vector<bam1_t*> &ret);

typedef struct{
  char *read;
  char *ref;
  uint8_t *qual;
  int len;
}asite;

sam_hdr_t * read_all_reads(const char *htsname,const char *refName, faidx_t *seq_ref,std::vector<asite> &ret,int len_limit);
