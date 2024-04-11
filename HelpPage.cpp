#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <ctime>
#include <cassert>
#include <cstdint>
#include <iostream>
#include <math.h>
#include "HelpPage.h"
#include "version.h"

int HelpPage(FILE *fp){
    fprintf(fp,"\t-> ./ngsBriggs -bam -tab -ref -bed -len -chr -ibam -iref -ibed -ichr -obam -otab -oinf -olen -model -eps -isrecal -olik -nthreads\n");
    fprintf(fp,"\t-> -bam: The bam file for inference;\n");
    fprintf(fp,"\t-> -tab: The table file for inference;\n");
    fprintf(fp,"\t-> -ref: The reference file for inference;\n");
    fprintf(fp,"\t-> -bed: The bed file (1-based) for inference;\n");
    fprintf(fp,"\t-> -len: The length mass probability distribution file for inference;\n");
    fprintf(fp,"\t-> -chr: The focal chromosome for inference;\n");
    fprintf(fp,"\t-> -ibam: The bam file for ancient strand fishing;\n");
    fprintf(fp,"\t-> -iref: The reference file for ancient strand fishing;\n");
    fprintf(fp,"\t-> -ibed: The bed file (1-based) for ancient strand fishing;\n");
    fprintf(fp,"\t-> -ichr: The focal chromosome for ancient strand fishing;\n");
    fprintf(fp,"\t-> -obam: The output bam file name;\n");
    fprintf(fp,"\t-> -otab: The output table file name;\n");
    fprintf(fp,"\t-> -oinf: The output inferred parameters file name;\n");
    fprintf(fp,"\t-> -olen: The output length mass probability distribution file name;\n");
    fprintf(fp,"\t-> -model: Specifying the model, either b (the biotin model) or nb (the non-biotin model);\n");
    fprintf(fp,"\t-> -eps: The overall modern contamination rate, the value should be within the interval [0,1);\n");
    fprintf(fp,"\t-> -isrecal: Choose 1 if recalibration based on length is needed, otherwise 0 (default);\n");
    fprintf(fp,"\t-> -olik: The output nucleotide likelihood file;\n");
    fprintf(fp,"\t-> -nthreads: Choose the number of threads to speed up the recalibration process.\n");
    fprintf(fp,"\t-> -bdamage: The mismatch matrix in bdamage format for metagenomic framework.\n");
    fprintf(fp,"\t-> -rlens: The read length distributions for metagenomic framework.\n");
    fprintf(fp,"\t-> Examples.\n");
    fprintf(fp,"\t-> ./ngsbriggs -model nb -bdamage Chr22_024_36_68_0097.bdamage.gz -rlens Chr22_024_36_68_0097.rlens.gz\n");
    fprintf(fp,"\t-> ./ngsbriggs -tab mismatch2.txt -len len.txt -model nb\n");
    fprintf(fp,"\t-> ./ngsbriggs -bam Chr22_024_36_68_0097.sorted.MD.bam -ref chr22.fa -model nb\n");
    exit(1);
  return 0;
}
