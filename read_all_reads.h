#include <vector>
#include <htslib/sam.h>

sam_hdr_t * read_all_reads(const char *htsname, const char *bedfile,const char *refName,std::vector<bam1_t*> &ret);
