#include <htslib/hts.h>
#include <htslib/sam.h>
#include <htslib/thread_pool.h>
#include <htslib/kstring.h>
#include <cassert>
#include <cstdlib>
#include <cstring>
#include <cstdio>
#include <map>
#include <vector>
#include <cmath>
#include <getopt.h>
#include <ctime>

typedef struct{
  bam1_t **d;//data
  unsigned l;//at pos
  unsigned m;//maxpos;
}queue_t;

queue_t *init_queue_t(int l);
extern int nreads_per_pos;//assuming this is the number reads per pos. If larger the program will reallocate
extern htsFormat *dingding6;
void realloc_queue(queue_t *q);
void do_magic(queue_t *q,bam_hdr_t *hdr,samFile *fp);
int usage1(FILE *fp, int is_long_help);

