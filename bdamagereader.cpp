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

#include "profile.h"

#define MAX_LENGTH 1000 // Define the maximum length as required
#define NUM_COLUMNS 16  // Define the number of columns to parse
#define STR_LENS 1024    // Define the maximum length of a line in the file

void parse_bdamage_Data(mydataD &md,double **dat,int howmany) {
    
    for (int i = 0; i < 15; i++) {
        // Populate dat with counts for forward direction
        for (int j = 0; j < 16; j++) {
            dat[i][j] = md.fwD[i *  16 + j]; // Assuming dat contains counts for each nucleotide
        }

        // Populate dat with counts for backward direction
        for (int j = 0; j < 16; j++) {
            dat[i+15][j] = md.bwD[i * 16 + j]; // Assuming dat contains counts for each nucleotide
        }
    }
}



void removeCounts(int*& rlen, double*& count, int& rlen_length, int& count_length,int len_limit,int &no_val) {
    // Create vectors to store non-zero values temporarily
    std::vector<int> temp_rlen;
    std::vector<double> temp_count;

    // Iterate through the arrays and copy non-zero values to temporary vectors
    for (int i = 0; i < count_length; ++i) {
        if (rlen[i] > 0 && rlen[i] < len_limit) {
            //std::cout << " i val " << i << std::endl;
            temp_rlen.push_back(rlen[i]);
            temp_count.push_back(count[i]);
        }
    }

    // Update rlen and count arrays and their lengths
    rlen_length = temp_rlen.size();
    count_length = temp_count.size();

    // Deallocate the old memory
    delete[] rlen;
    delete[] count;

    // Allocate new memory
    rlen = new int[rlen_length];
    count = new double[count_length];

    // Copy values from vectors to arrays
    for (int i = 0; i < rlen_length; ++i) {
        rlen[i] = temp_rlen[i];
    }
    for (int i = 0; i < count_length; ++i) {
        count[i] = temp_count[i];
    }
    
    no_val = count_length;
}

void FragArrayReaderRlen(int len_limit,int &no_val,int*& rlen, double*& count,const char* filename) {
    int STRLENS = 4096;


    char line[STRLENS];
    char *token;
    int line_number = 0; // Line number counter
    int rlen_length = 0;
    int count_length = 0;
    int n = 0;
    int m = 0;
    double count_sum = 0;
    gzFile file = gzopen(filename, "rb"); // Replace "your_file_name.txt.gz" with your actual gzipped file name
    
    // Read each line from the file
    while (gzgets(file, line, STRLENS) != NULL) {
        //fprintf(stderr,"Line %d:\n", line_number);

        // Tokenize the line using tab as delimiter
        token = strtok(line, "\t");

        // Skip the first token (id column)
        token = strtok(NULL, "\t");

        // Print each subsequent token (column) in the line
        while (token != NULL) {
            if (line_number == 0){
                // Check if the line length is greater than 4 characters
                if (strlen(token) > 4) {
                    // Copy the line content starting from the 5th character
                    memmove(token, token + 4, strlen(token) - 4 + 1);
                    rlen[n] = atoi(token);
                    //fprintf(stderr,"test %d \n",atoi(token));
                    n++;
                }
            }
            else if (line_number == 1) {
                count[m] = atof(token);
                count_sum += count[m];
                m++;
            }

            //printf("%s\n", token);
            token = strtok(NULL, "\t");
        }
        line_number++; // Increment line number

    }

    no_val = n;

    removeCounts(rlen, count, n, m,len_limit,no_val);

    for(int i=0;i<no_val;i++){
        count[i] = count[i]/count_sum;
    }

    gzclose(file);
}

// Can deal with both counts and freqs.


#ifdef __WITH_MAIN__
//g++ bdamagereader.cpp ../htslib/libhts.a -std=c++11 -lz -lm -lbz2 -llzma -lpthread -lcurl -lcrypto -ggdb -D__WITH_MAIN__
int main(){
    int* rlen_length = new int[4096];
    double* rlen_count = new double[4096];
    int no_values = 0;
    int len_limit = 150;
    FragArrayReaderRlen(rlen_length,rlen_count,no_values,"Chr22_024_36_68_0097.rlens.gz",len_limit);

    /*
    for (int i = 0; i < no_values; i++) {
        fprintf(stderr,"Length %d and count %f \n",rlen_length[i],rlen_count[i]);
    }*/
    
    exit(1);
    int howmany = 16;
    int MAXLENGTH = 15;
    double** mm5p, **mm3p;
    double *freqCT, *freqGA, *scaleCT, *scaleGA, *seqError;
    int len_min;

    // Load data using load_bdamage_full
    std::map<int, mydataDtmp> retmap = load_bdamage("Chr22_024_36_68_0097.bdamage.gz",howmany);
    
    fprintf(stderr, "\t-> %lu mismatch matrices read for %d base pairs\n", retmap.size(), howmany);

    for (std::map<int, mydataDtmp>::iterator it = retmap.begin(); it != retmap.end(); it++) {
        int taxid = it->first;
        mydataDtmp md = it->second;
        if (it->second.nreads == 0)
        continue;
        
        // allocating memory for my mm5p and mm3p such that i can incorporate the tables 
        mm5p = (double**) malloc(MAXLENGTH * sizeof(double*));
        mm3p = (double**) malloc(MAXLENGTH * sizeof(double*));
        for (int i = 0; i < MAXLENGTH-1; i++){
            mm5p[i] =(double *) malloc(16 * sizeof(double));
            mm3p[i] =(double *) malloc(16 * sizeof(double));
        }
        mm5p[MAXLENGTH-1] =(double *) malloc(16 * sizeof(double));
        mm3p[MAXLENGTH-1] =(double *) malloc(16 * sizeof(double));
        for (int i=0; i<MAXLENGTH-1;i++){
            for (int j=0; j<16;j++){
                mm5p[i][j]=0;
                mm3p[i][j]=0;
            }
        }
        for (int j=0; j<16;j++){
            mm5p[MAXLENGTH-1][j]=0;
            mm3p[MAXLENGTH-1][j]=0;
        }
        
        freqCT = (double*) malloc(2*MAXLENGTH * sizeof(double));
        freqGA = (double*) malloc(2*MAXLENGTH * sizeof(double));
        scaleCT = (double*) malloc(2*MAXLENGTH * sizeof(double));
        scaleGA = (double*) malloc(2*MAXLENGTH * sizeof(double));
        seqError = (double*) malloc(2*MAXLENGTH * sizeof(double));
        
	int numpos = MAXLENGTH*2+1;
        int numcolumn = 16;
        double ** Table = (double **) malloc(numpos*(sizeof(double *))); /*I allocate memory here.  If this function is called many times it may be better to move the memmory allocation out of this function*/
        for (int i=0; i<numpos; i++){
            Table[i]=(double *) malloc(numcolumn*(sizeof(double)));
        }  

        parse_bdamage_Data(md, Table, howmany);
        
        for (int i=0; i<MAXLENGTH; i++){
            for (int j=0; j<MAXLENGTH;j++){
                mm5p[i][j] = Table[i][j];
                mm3p[i][j] = Table[i+howmany][j];
            }
        }
        
        
       
        //freqCT, freqGA, scaleCT, scaleGA
        double max5=0;
        double max3=0;
        double maxall=0;
        for (int i=0; i<MAXLENGTH;i++){
            // i -> 0 to 14
            // 2*MAXLENGTH-1-i -> 29 to 15
            scaleCT[i] = mm5p[i][4]+mm5p[i][5]+mm5p[i][6]+mm5p[i][7];
            scaleCT[2*MAXLENGTH-1-i] = mm3p[i][4]+mm3p[i][5]+mm3p[i][6]+mm3p[i][7];
            freqCT[i] = mm5p[i][7]/scaleCT[i];
            freqCT[2*MAXLENGTH-1-i] = mm3p[i][7]/scaleCT[2*MAXLENGTH-1-i];
            //std::cout << freqCT[i] << std::endl;

            scaleGA[i] = mm3p[i][8]+mm3p[i][9]+mm3p[i][10]+mm3p[i][11];
            scaleGA[2*MAXLENGTH-1-i] = mm5p[i][8]+mm5p[i][9]+mm5p[i][10]+mm5p[i][11];
            freqGA[i] = mm3p[i][8]/scaleGA[i];
            freqGA[2*MAXLENGTH-1-i] = mm5p[i][8]/scaleGA[2*MAXLENGTH-1-i];
            
            max5 = std::max(max5,std::max(scaleCT[i],scaleCT[2*MAXLENGTH-1-i]));
            max3 = std::max(max3,std::max(scaleGA[i],scaleGA[2*MAXLENGTH-1-i]));

            // Overall sequencing errors are position specific, estimated by 1 - [N(AA)+N(TT)]/[N(AA)+N(AC)+N(AG)+N(AT)+N(TA)+N(TC)+N(TG)+N(TT)]
            seqError[i] = 1 - (mm5p[i][0]+mm5p[i][15])/(mm5p[i][0]+mm5p[i][1]+mm5p[i][2]+mm5p[i][3]+mm5p[i][12]+mm5p[i][13]+mm5p[i][14]+mm5p[i][15]);
            seqError[2*MAXLENGTH-1-i] = 1 - (mm3p[i][0]+mm3p[i][15])/(mm3p[i][0]+mm3p[i][1]+mm3p[i][2]+mm3p[i][3]+mm3p[i][12]+mm3p[i][13]+mm3p[i][14]+mm3p[i][15]);
        }
        maxall = std::max(max5,max3);
        for (int i=0; i<MAXLENGTH;i++){
            scaleCT[i] = scaleCT[i]/maxall;
            scaleGA[i] = scaleGA[i]/maxall;
            scaleCT[2*MAXLENGTH-1-i] = scaleCT[2*MAXLENGTH-1-i]/maxall;
            scaleGA[2*MAXLENGTH-1-i] = scaleGA[2*MAXLENGTH-1-i]/maxall;
        }

        fprintf(stderr,"\nThe misincorporting matrix is as follows:\n");
        fprintf(stderr,"Dir.\tPos.\tFreqCT\tFreqGA\n");
        for (int i=0; i<5;i++){
            fprintf(stderr,"5'\t%d\t",i+1);
            fprintf(stderr,"%f\t%f\n",freqCT[i],freqGA[2*MAXLENGTH-1-i]);
        }
        for (int i=0; i<5;i++){
            fprintf(stderr,"3'\t%d\t",i+1);
            fprintf(stderr,"%f\t%f\n",freqCT[2*MAXLENGTH-1-i],freqGA[i]);
        }
    }

    return 0;
}
#endif
