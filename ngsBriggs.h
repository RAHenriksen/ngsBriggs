#ifndef NGSBRIGGS_H
#define NGSBRIGGS_H

//MAXLENGTH = 5;
typedef unsigned char uchar;
extern int nreads_per_pos;
typedef struct{
    char bp;
    int offset;
} mdField;

extern int nproc1;//number of reads processed not used
extern char refToChar[256];
extern char toIndex[4][4];
extern char com[5];
extern unsigned char toDIndex[16][16];
extern int MAXLENGTH;
extern double** mm5p, **mm3p, **dmm5p, **dmm3p;

extern double Tol; // Tolerance

static char DUMMYCHAR='#';

extern double l_check;

extern int ncalls;
extern int ncalls_grad;

// I finally decide to use the accurate but the most complex model to check my inference
extern int MAXORDER; //This can be adjusted for different tolerance of errors.
//Check the code again
// OTHER DEFINES VALUES

extern char *refName2, *fname2, *model2;
extern const char* chromname2, * bedname2;
extern int mapped_only2, se_only2, mapq2, len_limit2, len_min2;
extern double eps2, lambda2, delta2, delta_s2, nv2;
extern std::string s2;
extern faidx_t *seq_ref2;

extern std::vector<nuclist> comp_nuc_llik;

int main(int argc, char **argv);

#endif
