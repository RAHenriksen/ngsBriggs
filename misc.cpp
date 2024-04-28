#include <cstring>
#include <cstdlib>
#include <cstdio>
#include <cassert>
#include <cmath>
#include <ctime>
#include <iostream>
#include <zlib.h>
#include <array>
#include <fstream>
#include <sstream>
#include <algorithm>

#include <htslib/hts.h>
#include <htslib/sam.h>
#include <htslib/kstring.h>
#include <htslib/faidx.h>
#include <htslib/faidx.h>
#include <htslib/bgzf.h>

#include "profile.h"
#include "bfgs.h"
#include "misc.h"
#include "ngsBriggs.h"

#define PI M_PI


double NormalCDF(double x) // Phi(-∞, x) aka N(x)
{
    return erfc(-x / sqrt(2))/2;
}

double NormalPDF(double x)
{
    return 1/sqrt(2*PI)*exp(-pow(x,2)/2);
}

double NormalPDF_grad(double x){
    return 1/sqrt(2*PI)*exp(-pow(x,2)/2)*(-x);
}

double NormalINC(double y, double x, double x_max, double x_min)
{
    if (y>=x){
        return (NormalCDF(y)-NormalCDF(x))/(NormalCDF(x_max)-NormalCDF(x_min));
    }else{
        return -1;
    }
}

double NormalINC_grad_mu(double y, double x, double x_max, double x_min, double mu, double sigma){
    if (y>=x){
        return ((NormalPDF(y)*(-1/sigma)-NormalPDF(x)*(-1/sigma))*(NormalCDF(x_max)-NormalCDF(x_min))-(NormalPDF(x_max)*(-1/sigma)-NormalPDF(x_min)*(-1/sigma))*(NormalCDF(y)-NormalCDF(x)))/pow(NormalCDF(x_max)-NormalCDF(x_min),2);
    }
    else{
        return -1;
    }
}


double NormalINC_grad_si(double y, double x, double x_max, double x_min,double mu, double sigma){
    if (y>=x){
        return ((NormalPDF(y)*(-y/sigma)-NormalPDF(x)*(-x/sigma))*(NormalCDF(x_max)-NormalCDF(x_min))-(NormalPDF(x_max)*(-x_max/sigma)-NormalPDF(x_min)*(-x_min/sigma))*(NormalCDF(y)-NormalCDF(x)))/pow(NormalCDF(x_max)-NormalCDF(x_min),2);
    }else{
        return -1;
    }
}

double NormalINC_hess_mu2(double y, double x, double x_max, double x_min,double mu, double sigma){
    if (y>=x){
        double D3 = 1/pow(sigma,2)*(NormalPDF_grad(y)-NormalPDF_grad(x))/(NormalCDF(x_max)-NormalCDF(x_min));
        double D2 = (-2/pow(sigma,2)*(NormalPDF(y)-NormalPDF(x))*(NormalPDF(x_max)-NormalPDF(x_min))-1/pow(sigma,2)*(NormalCDF(y)-NormalCDF(x))*(NormalPDF_grad(x_max)-NormalPDF_grad(x_min)))/pow(NormalCDF(x_max)-NormalCDF(x_min),2);
        double D1 = 2/pow(sigma,2)*(NormalCDF(y)-NormalCDF(x))*pow(NormalPDF(x_max)-NormalPDF(x_min),2)/pow(NormalCDF(x_max)-NormalCDF(x_min),3);
        return D3+D2+D1;
    }else{
        return -1;
    }
}

double NormalINC_hess_si2(double y, double x, double x_max, double x_min,double mu, double sigma){
    if (y>=x){
        double D3 = (NormalPDF_grad(y)*pow(y/sigma,2)-NormalPDF_grad(x)*pow(x/sigma,2)+NormalPDF(y)*2*y/pow(sigma,2)-NormalPDF(x)*2*x/pow(sigma,2))/(NormalCDF(x_max)-NormalCDF(x_min));
        double D2 = (-(NormalPDF_grad(x_max)*pow(x_max/sigma,2)-NormalPDF_grad(x_min)*pow(x_min/sigma,2)+NormalPDF(x_max)*2*x_max/pow(sigma,2)-NormalPDF(x_min)*2*x_min/pow(sigma,2))*(NormalCDF(y)-NormalCDF(x))-2*(NormalPDF(y)*y/sigma-NormalPDF(x)*x/sigma)*(NormalPDF(x_max)*x_max/sigma-NormalPDF(x_min)*x_min/sigma))/pow(NormalCDF(x_max)-NormalCDF(x_min),2);
        double D1 = 2*pow(NormalPDF(x_max)*x_max/sigma-NormalPDF(x_min)*x_min/sigma,2)*(NormalCDF(y)-NormalCDF(x))/pow(NormalCDF(x_max)-NormalCDF(x_min),3);
        return D3+D2+D1;
    }else{
        return -1;
    }
}

double NormalINC_hess_mu_si(double y, double x, double x_max, double x_min,double mu, double sigma){
    if (y>=x){
        double D3 = (1/pow(sigma,2)*(NormalPDF(y)-NormalPDF(x))+(NormalPDF_grad(y)*y/pow(sigma,2)-NormalPDF_grad(x)*x/pow(sigma,2)))/(NormalCDF(x_max)-NormalCDF(x_min));
        double D2 = (-1/pow(sigma,2)*(NormalCDF(y)-NormalCDF(x))*(NormalPDF(x_max)-NormalPDF(x_min))-(NormalPDF(y)*y/pow(sigma,2)-NormalPDF(x)*x/pow(sigma,2))*(NormalPDF(x_max)-NormalPDF(x_min))-(NormalCDF(y)-NormalCDF(x))*(NormalPDF_grad(x_max)*x_max/pow(sigma,2)-NormalPDF_grad(x_min)*x_min/pow(sigma,2))-(NormalPDF(y)-NormalPDF(x))*(NormalPDF(x_max)*x_max/pow(sigma,2)-NormalPDF(x_min)*x_min/pow(sigma,2)))/pow(NormalCDF(x_max)-NormalCDF(x_min),2);
        double D1 = 2*(NormalCDF(y)-NormalCDF(x))*(NormalPDF(x_max)-NormalPDF(x_min))*(NormalPDF(x_max)*x_max/pow(sigma,2)-NormalPDF(x_min)*x_min/pow(sigma,2))/pow(NormalCDF(x_max)-NormalCDF(x_min),3);
        return D3+D2+D1;
    }else{
        return -1;
    }
}

//Bin the lengths to speed up the inference
void FragArrayBin(int number, int &BinNum, int*& Length, double *& Freq, double*& BinLength, double*& BinFreq){
    if((BinNum >= number) || (BinNum)<=0){
        // Can be further optimized, imagin the distribution is not increasing by 1 each time.
        BinNum = number;
        for(int j=0;j<BinNum;j++){
            BinLength[j] = (double)Length[j];
            BinFreq[j] = Freq[j];
        }
    }else{
        double d = (double)number/(double)BinNum;
        double i, f;
        BinLength[0] = Length[0]-0.5;
        BinFreq[0] = Freq[0];
        int k = 0;
        for(int j=1; j<BinNum; j++){
            f = modf(j*d, &i);
            BinLength[j] = Length[(int)i-1]+0.5+f;
            for (int l=k+1;l<(int)i;l++){
                BinFreq[j-1] += Freq[l];
            }
            double f1=BinLength[j]-round(BinLength[j])+0.5;
            BinFreq[j-1] += f1*Freq[(int)i];
            BinFreq[j] += (1-f1)*Freq[(int)i];
            k = (int)i;
            BinLength[j-1] = (BinLength[j-1]+BinLength[j])/2;
        }
        for (int l=k+1;l<number;l++){
            BinFreq[BinNum-1] += Freq[l];
        }
        BinLength[BinNum-1] = (BinLength[BinNum-1]+Length[number-1]+0.5)/2;
    }
}

// Can deal with both counts and freqs.
void FragArrayReader(int len_limit, int& number, int*& Length, double *& Freq, const char* filename,int STRLENS){
    int n = 0;
    int m = 0;
    double Sumf = 0;
    gzFile gz = Z_NULL;
    char buf[STRLENS];
    
    gz = gzopen(filename, "r");
    assert(gz!=Z_NULL);
    if (len_limit>0){
        while(gzgets(gz,buf,STRLENS)){
            Length[n] = atoi(strtok(buf,"\n\t "));
            Freq[n] = atof(strtok(NULL,"\n\t "));
            //std::std::cout << "length "<< Length[n] << " FASERQ " << Freq[n] << std::endl;
            if(Length[n] < len_limit){
                m = n;
                Sumf += Freq[n];
            }
            n++;
        }
        m++;
        gzclose(gz);
        number = m;
    }else{
        while(gzgets(gz,buf,STRLENS)){
            Length[n] = atoi(strtok(buf,"\n\t "));
            Freq[n] = atof(strtok(NULL,"\n\t "));
            Sumf += Freq[n];
            n++;
        }
        gzclose(gz);
        number = n;
    }
    for(int i=0;i<number;i++){
        Freq[i] = Freq[i]/Sumf;
    }
}

void parse_tabledata(const char* filename,double** Table,int STRLENS){
    //int numpos = MAXLENGTH*2;
    //int numcolumn = 16 + 4;
    int i = 0;
    //double Sumf = 0;
    gzFile gz = Z_NULL;
    char buf[STRLENS];
    
    gz = gzopen(filename, "r");
    assert(gz!=Z_NULL);
    while(gzgets(gz,buf,STRLENS)){
        //Length[n] = atoi(strtok(buf,"\n\t "));
        strtok(buf,"\t\n ");
        strtok(NULL,"\t\n ");
        for (int j=0;j<16;j++){
            Table[i][j+2]=atof(strtok(NULL,"\t\n "));
        }
        i++;
        //        Sumf += Freq[n];
        //        n++;
    }

    gzclose(gz);
}

int tabreader(char *tabname,int STRLENS,double** mm5p, double **mm3p){
    int numpos = MAXLENGTH*2+1;
    int numcolumn = 16 + 4;
    double ** Table = (double **) malloc(numpos*(sizeof(double *))); /*I allocate memory here.  If this function is called many times it may be better to move the memmory allocation out of this function*/
    for (int i=0; i<numpos; i++){
        Table[i]=(double *) malloc(numcolumn*(sizeof(double)));
    }
    
    clock_t t=clock();
    time_t t2=time(NULL);
    
    if(tabname){
        //parse_tabledata(tabname, Table);
      parse_tabledata(tabname,Table,STRLENS);
        fprintf(stderr,"The tab file %s has been inputted successfully \t %d\n",tabname,MAXLENGTH);
    }
    
    for (int i=0; i<MAXLENGTH; i++){
        for (int j=0; j<16;j++){
            mm5p[i][j] = Table[i+1][j+2];
            mm3p[i][j] = Table[i+MAXLENGTH+1][j+2];
        }
    }
    for (int i = 0; i < 2*MAXLENGTH+1; i++){
        free(Table[i]);
    }
    free(Table);
    return 0;
}


int bamreader(char *fname, const char* chromname,faidx_t * seq_ref, int len_limit, int &len_min,int *Frag_len, double *Frag_freq,int &number,double **mm5p,double **mm3p){
    char *refName = NULL;
    clock_t t=clock();
    time_t t2=time(NULL);
    
    int mapq =-1;
    int mapped_only = 0;
    int se_only = 1;
    
    if(fname){
        fprintf(stderr,"The bam file %s has been inputted successfully\n",fname);
        if (seq_ref!=NULL){
            fprintf(stderr, "Reference is provided!\n");
        }else{
            fprintf(stderr, "Reference is not provided, and it will be reconstructed according to the MD tags!\n");
        }
        parse_sequencingdata1(refName,fname,chromname,mapped_only,se_only,mapq,seq_ref,len_limit,len_min,Frag_len,Frag_freq,number,mm5p,mm3p);
    }
    
    fprintf(stderr,
            "\n"
            "\t[ALL done] cpu-time used =  %.2f sec\n"
            "\t[ALL done] walltime used =  %.2f sec\n"
            ,(float)(clock() - t) / CLOCKS_PER_SEC, (float)(time(NULL) - t2));
    return 0;
}

void wrapperwithref(const bam1_t   * b,const bam_hdr_t  *hdr, char myread[512], char myref[512],faidx_t *seq_ref){
    const char *alphabetHTSLIB = "NACNGNNNTNNNNNNN";
    char reconstructedTemp[512];
    memset(reconstructedTemp,0,512);
    
    //skip unmapped
    if( ((b)->core.flag&BAM_FUNMAP) != 0 ){
        fprintf(stderr,"The function reconstructRefWithPosOnReadHTS() cannot be called for unmapped reads\n");
        exit(1);
    }
    
    int32_t   n_cigar_op = b->core.n_cigar;
    uint32_t *cigar      = bam_get_cigar(b);
    
    //std::cout<<"Hello "<<n_cigar_op<<"\n";
    int at =0;
    for(int32_t i = 0; i < n_cigar_op; i++){
        char opchr = bam_cigar_opchr(cigar[i]);
        int32_t oplen = bam_cigar_oplen(cigar[i]);
        //std::cout<<opchr<<" "<< oplen <<"\n";
        memset(reconstructedTemp+at,opchr,oplen);
        at += oplen;
    }
    //std::cout<<"\n";
    int initialPositionControl=b->core.pos;
    
    for (unsigned int k=0;k<(int)b->core.l_qseq;k++){
        myread[k] = toupper(alphabetHTSLIB[ bam_seqi(bam_get_seq(b),k) ]);
    }
    
    //combine the CIGAR and REF to reconstruct both Read and REF
    unsigned int i = 0;
    unsigned int j = 0;
    int start = initialPositionControl;
    int end = start-1;
    while (i<strlen(reconstructedTemp)){
        while (reconstructedTemp[i] == 'M' ){
            i++;
            end++;
        }
        int chr_len = faidx_seq_len(seq_ref,hdr->target_name[b->core.tid]);
        char* data = faidx_fetch_seq(seq_ref, hdr->target_name[b->core.tid], start, end, &chr_len);
        for (unsigned int k=0;k<end-start+1;k++){
            myref[j] = toupper(data[k]);
            j++;
        }
        free(data);
        //        i++;
        start = end+1;
        end = start-1;
        while (reconstructedTemp[i] == 'D' ){
            i++;
            end++;
        }
        //        i++;
        start = end+1;
        end = start-1;
        while(reconstructedTemp[i] == 'S' || reconstructedTemp[i] == 'I'){
            i=i+1;
            //            Ns have already been added.
            //            myref[j] = 'N';
            myread[j] = 'N';
            j++;
        }
    }
}


void CaldeamRate_b(double lambda, double delta, double delta_s, double nu, int len_limit,double **deamRateCT,double **deamRateGA){
    for(int n=0; n<MAXLENGTH; n++){
        //    int n = 0;
        for(int L=30; L<len_limit; L++){
            //double f = freqLEN[i];
            double p1_l = pow(1+lambda,2)*n*nu/(1+(L-2)*nu);
            for(int l=1;l<=std::min(MAXLENGTH,n);l++){
                p1_l += lambda*(1+lambda)*pow(1-lambda,l)*(n-l)*nu/(1+(L-2-l)*nu);
                for (int r=1;r<=std::min((double)MAXLENGTH,(double)(L-n-1));r++){
                    if (l+r<=MAXLENGTH){
                        p1_l += pow(lambda,2)*pow(1-lambda,l+r)*(n-l)*nu/(1+(L-2-l-r)*nu);
                    }
                }
            }
            for(int r=1;r<=std::min((double)MAXLENGTH,(double)(L-n-1));r++){
                p1_l += lambda*(1+lambda)*pow(1-lambda,r)*n*nu/(1+(L-2-r)*nu);
            }
            p1_l = p1_l/4;
            double p2_l = 0;
            double p3_l = 1 - pow(1-lambda,n+1)/2-p1_l;
            double p4_l = pow(1-lambda,n+1)/2;
            
            double p1_r = pow(1+lambda,2)*(L-n-1)*nu/(1+(L-2)*nu);
            for(int r=1;r<=std::min(MAXLENGTH,n);r++){
                p1_r += lambda*(1+lambda)*pow(1-lambda,r)*(L-n-1)*nu/(1+(L-2-r)*nu);
                for (int l=1;l<=std::min((double)MAXLENGTH,(double)(L-n-1));l++){
                    if (l+r<=MAXLENGTH){
                        p1_r += pow(lambda,2)*pow(1-lambda,l+r)*(L-n-1-l)*nu/(1+(L-2-l-r)*nu);
                    }
                }
            }
            for(int l=1;l<=std::min((double)MAXLENGTH,(double)(L-n-1));l++){
                p1_r += lambda*(1+lambda)*pow(1-lambda,l)*(L-n-1-l)*nu/(1+(L-2-l)*nu);
            }
            p1_r = p1_r/4;
            double p2_r = pow(1-lambda,n+1)/2;
            double p3_r = 1 - pow(1-lambda,n+1)/2-p1_r;
            double p4_r = 0;
            deamRateCT[L-30][n] = (p3_l*delta+p4_l*delta_s);
            deamRateGA[L-30][n] = (p1_r*delta+p2_r*delta_s);
            deamRateCT[L-30][30-1-n] = (p3_r*delta);
            deamRateGA[L-30][30-1-n] = (p1_l*delta);
        }
    }
}

void CaldeamRate_nb(double lambda, double delta, double delta_s, double nu, int len_limit,double **deamRateCT,double **deamRateGA){
    for(int n=0; n<MAXLENGTH; n++){
        //double freqCT1 = 0;
        //double freqGA1 = 0;
        //double freqCT2 = 0;
        //double freqGA2 = 0;
        for(int L=30; L<len_limit; L++){
            double p1_l = pow(1+lambda,2)*n*nu/(1+(L-2)*nu);
            for(int l=1;l<=std::min(MAXLENGTH,n);l++){
                p1_l += lambda*(1+lambda)*pow(1-lambda,l)*(n-l)*nu/(1+(L-2-l)*nu);
                for (int r=1;r<=std::min((double)MAXLENGTH,(double)(L-n-1));r++){
                    if (l+r<=MAXLENGTH){
                        p1_l += pow(lambda,2)*pow(1-lambda,l+r)*(n-l)*nu/(1+(L-2-l-r)*nu);
                    }
                }
            }
            for(int r=1;r<=std::min((double)MAXLENGTH,(double)(L-n-1));r++){
                p1_l += lambda*(1+lambda)*pow(1-lambda,r)*n*nu/(1+(L-2-r)*nu);
            }
            p1_l = p1_l/4;
            double p2_l = 0;
            double p3_l = 1 - pow(1-lambda,n+1)/2-p1_l;
            double p4_l = pow(1-lambda,n+1)/2;
            
            double p1_r = pow(1+lambda,2)*(L-n-1)*nu/(1+(L-2)*nu);
            for(int r=1;r<=std::min(MAXLENGTH,n);r++){
                p1_r += lambda*(1+lambda)*pow(1-lambda,r)*(L-n-1)*nu/(1+(L-2-r)*nu);
                for (int l=1;l<=std::min((double)MAXLENGTH,(double)(L-n-1));l++){
                    if (l+r<=MAXLENGTH){
                        p1_r += pow(lambda,2)*pow(1-lambda,l+r)*(L-n-1-l)*nu/(1+(L-2-l-r)*nu);
                    }
                }
            }
            for(int l=1;l<=std::min((double)MAXLENGTH,(double)(L-n-1));l++){
                p1_r += lambda*(1+lambda)*pow(1-lambda,l)*(L-n-1-l)*nu/(1+(L-2-l)*nu);
            }
            p1_r = p1_r/4;
            double p2_r = pow(1-lambda,n+1)/2;
            double p3_r = 1 - pow(1-lambda,n+1)/2-p1_r;
            double p4_r = 0;
            deamRateCT[L-30][n] = ((p1_r+p3_l)/2*delta+(p2_r+p4_l)/2*delta_s); //Lateral increments 1 for C to T (5') and G to A (3')
            deamRateCT[L-30][30-n-1] = (p1_l+p3_r)/2*delta; //Lateral increments 2 for C to T (3') and G to A (5')
            deamRateGA[L-30][n] = deamRateCT[L-30][n];
            deamRateGA[L-30][30-n-1] = deamRateCT[L-30][30-n-1];
    
        }

    }
}


//QUALITY CONTROL and LENGTH DISTRIBUTION
void parse_sequencingdata1(char *refName,char *fname,const char* chromname, int mapped_only,int se_only,int mapq, faidx_t *seq_ref,int len_limit, int & len_min,int *Frag_len, double *Frag_freq,int &number,double** mm5p, double **mm3p){
    fprintf(stderr,"mapped_only: %d\n",mapped_only);
    htsFormat *dingding3 =(htsFormat*) calloc(1,sizeof(htsFormat));

    char reconstructedRef[512];
    char myread[512];
    char myrefe[512];
    std::pair< kstring_t*, std::vector<int> >  mypair;
    kstring_t *kstr =new kstring_t;
    kstr->l=kstr->m=0;
    kstr->s=NULL;
    mypair.first = kstr;
    samFile *in=NULL;
    char refeBase, temprefeBase1, refeBase1, refeBase2;
    char readBase, tempreadBase1, readBase1, readBase2;
    
    //code below only relevant if using cram files
    if(refName){
        char *ref =(char*) malloc(10 + strlen(refName) + 1);
        snprintf(ref,10 + strlen(refName) + 1, "reference=%s", refName);
        hts_opt_add((hts_opt **)&dingding3->specific,ref);
        free(ref);
    }
    if(strstr(fname,".cram")!=NULL && refName==NULL){
        fprintf(stderr,"\t-> cram file requires reference with -T FILE.fa \n");
        exit(0);
    }
    
    if((in=sam_open_format(fname,"r",dingding3))==NULL ){
        fprintf(stderr,"[%s] nonexistant file: %s\n",__FUNCTION__,fname);
        exit(0);
    }
    
    
    bam_hdr_t  *hdr = sam_hdr_read(in);
    
    bam1_t *b = bam_init1();
    
    int chrom_num = hdr->n_targets;
    size_t * chrom_line = (size_t*)malloc(chrom_num*(sizeof(size_t)));
    std::vector<std::array<size_t, 2> > bedsites;
      
    int ret;
    int refId=-1;
    int chromId = -1;
    if (chromname!=NULL){
        for (int j=0; j<chrom_num; j++){
            if (!strcmp(chromname,hdr->target_name[j])){
                chromId = j;
                fprintf(stderr,"We will focus on Chromosome %s!\n", hdr->target_name[j]);
            }
        }
    }
    if (chromId==-1){
        fprintf(stderr,"No meaningful chromosome name is provided, therefore we will focus on all provided chromosomes!\n");
    }
    unsigned char * indref = NULL;
    double Frag_len_count[512];
    for (int j=0;j<512;j++){
        Frag_len_count[j] = 0.0;
    }
    double num = 0.0;
    size_t max_site;
    len_min = 512;
    while(((ret=sam_read1(in,hdr,b)))>0) {
        //if -m was used, discard unmapped reads
        if(mapped_only!=0){
            if(b->core.flag&4)
                continue;
        }
        
        //default -a 1
        //only use single end reads
        //which is either single end or collapsed reads
        //#if 0
        if(se_only==1){
            if(b->core.flag&1)//if paired continue
                continue;
        }
        //#endif
        if(b->core.flag>256)
            continue; //Remove possible duplicates
        
        // if mapq threshold is set
        if(mapq!=-1 && b->core.qual<mapq)
            continue;
        
        //Only count the sites belongs to a specific chromosome
        if (chromId==-1||chromId==b->core.tid){
            nproc1++;
        }else if (chromId!=b->core.tid){
            continue;
        }
        
        if(refId==-1||refId!=b->core.tid){
            refId=b->core.tid;
            fprintf(stderr,"\t-> Now at Chromosome: %s\n",hdr->target_name[refId]);
 
        }
        //fprintf(stderr,"readid:%s len: %d\nREAD:\t\n",bam_get_qname(b),b->core.l_qseq);
        // for(int i=0;i<b->core.l_qseq;i++)
        // fprintf(stderr,"%c",seq_nt16_str[bam_seqi(bam_get_seq(b),i)]);
        //fprintf(stderr,"\nQSCOre: \n");
        //for(int i=0;i<b->core.l_qseq;i++)
        //  fprintf(stderr,"%c",33+bam_get_qual(b)[i]);
        
        //then we simply write it to the output
        memset(reconstructedRef,0,512);
        memset(myread,'N',512);
        memset(myrefe,'N',512);
        
        if (seq_ref != NULL){
            //If reference is provided, calculate the reference based on the reference, the positions, the read and the CIGAR
            wrapperwithref(b,hdr,myread,myrefe,seq_ref);
        }else{
            //If reference is not provided, calculate the reference based on the MD tag and the read
            reconstructRefWithPosHTS(b,mypair,reconstructedRef);
            wrapper(b,mypair.first->s,mypair.second,0,0,NULL,NULL,1,myread,myrefe);
        }
        //      for(int i=0;i<b->core.l_qseq;i++){
        //          std::cout<<myread[i];
        //      }
        //      std::cout<<"\n";
        int temp = 0;
        if (len_limit <=0){
            len_limit = 512;
        };
        // If the fragment length is within the interval [30,len_limit)
        if (b->core.l_qseq>=30 && b->core.l_qseq<len_limit){
            len_min = (b->core.l_qseq > len_min ? len_min : b->core.l_qseq);
            //Using the global length distribution to mimic the local distribution (if bed file if provided), and this might lead to some biases and can be fixed by counting the length distribution cycle position wisely, but I am lazy for now.
            Frag_len_count[b->core.l_qseq] = Frag_len_count[b->core.l_qseq] + 1.00;
            num = num + 1.00;
            
            
            for (int cycle=0;cycle<b->core.l_qseq;cycle++){
                refeBase = refToChar[myrefe[cycle]];
                readBase = refToChar[myread[cycle]];
                size_t pos = b->core.pos+cycle;
                //std::cout << max_site << "\n";
                //if (pos<max_site){
                //   std::cout<<pos<<" "<<max_site<<"\n";
                //}
                
		if(refeBase!=4 && readBase!=4 ){
                    int dist5p=cycle;
                    int dist3p=b->core.l_qseq-1-cycle;
                    if( bam_is_rev(b) ){
                        // std::cout<<"rev "<<"\t";
                        refeBase=com[refeBase];
                        readBase=com[readBase];
                        //dist5p=int(al.QueryBases.size())-1-i;
                        dist5p=int(b->core.l_qseq)-1-cycle;
                        dist3p=cycle;
                    }
                    if(dist5p<MAXLENGTH)
                        mm5p[dist5p][toIndex[refeBase][readBase]] += 1.00;
                    if(dist3p<MAXLENGTH)
                        mm3p[dist3p][toIndex[refeBase][readBase]] += 1.00;
                    temp = temp + 1;
                }
            }
        }
    }
    int l0 = 0;
    while (Frag_len_count[l0] == 0){
        l0++;
    }
    int l1 = 512-1;
    while (Frag_len_count[l1] == 0){
        l1--;
    }
    
    int k = 0;
    for (int l2 = l0; l2<=l1; l2++){
        Frag_len[k] = l2;
        Frag_freq[k] = Frag_len_count[l2]/num;
        k = k+1;
    }
    number = k;
    
    if (indref!=NULL){
        free(indref);
    }
    bedsites.clear();
    bedsites.shrink_to_fit();
    bam_destroy1(b);
    hts_opt_free((hts_opt *)dingding3->specific);
    free(dingding3);
    sam_hdr_destroy(hdr);
    sam_close(in);
    //free(fname);
    
}

//STOP HERE
