#ifndef BDAMAGEREADER_H
#define BDAMAGEREADER_H

#include <cstring>
#include <cstdlib>
#include <cstdio>
#include <cassert>
#include <htslib/hts.h>
#include <htslib/sam.h>
#include <htslib/kstring.h>
#include <htslib/faidx.h>
#include <htslib/faidx.h>
#include "htslib/bgzf.h"
#include <zlib.h>
#include <map>
#include <cassert>
#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <zlib.h>
#include <iostream>
#include <vector>

void removeCounts(int*& rlen, double*& count, int& rlen_length, int& count_length,int len_limit,int &no_val);

void FragArrayReaderRlen(int len_limit, int &no_val,int*& rlen, double*& count,const char* filename);

void parse_bdamage_Data(mydataD &md,double **dat,int howmany);

#endif
