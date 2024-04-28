#include <cstdlib>
#include <cstdio>
#include <zlib.h>
#include <cmath>
#include <ctime>


#include <htslib/hts.h>
#include <htslib/sam.h>
#include <htslib/kstring.h>
#include <htslib/faidx.h>

#include "profile.h"
#include "bfgs.h"
#include "htslib/bgzf.h"

#include "read_all_reads.h"

#include "misc.h"
#include "Recalibration.h"
#include "Likelihood.h"
#include "PosteriorProb.h"
#include "ngsBriggs_cli.h"
#include "ngsBriggs.h"
#include "bdamagereader.h"

// definining all global variables used across multiple scripts

int nproc1 = 0;//number of reads processed maybe not used
double l_check = 15;
int ncalls =0;
int ncalls_grad =0;
tsk_struct *my_tsk_struct = NULL;

int tsk_nthreads = -1;

// defining our main ngsBriggs function
int main(int argc, char **argv){

  double** mm5p, **mm3p;// **dmm5p, **dmm3p;
  double Tol = 1.0E-8; // Tolerance
  double **deamRateCT;
  double **deamRateGA;
  int number;
  double *freqCT, *freqGA, *scaleCT, *scaleGA, *seqError;
  double Contam_eps; //Contamination rate defined as the proportion of the number of the contaminated reads

  for(int i=0;i<255;i++){
    double di = i;
    PhredError[i] = pow(10,-di/10.0);
    PhredErrorAThird[i] = PhredError[i]/3.0;
  }

  int STRLENS = 4096;
  int* Frag_len = new int[STRLENS];
  double* Frag_freq = new double[STRLENS];
  int BinNum = -1;
  double* Bin_Frag_len = new double[STRLENS];
  double* Bin_Frag_freq = new double[STRLENS];
  
  argStruct *mypars = NULL;
   if(argc==1||(argc==2&&(strcasecmp(argv[1],"--help")==0||strcasecmp(argv[1],"-h")==0))){
     HelpPage(stderr);
     return 0;
   }

   kstring_t str_cli; str_cli.s=NULL;str_cli.l=str_cli.m=0;
   ksprintf(&str_cli,"./ngsBriggs");;


    for(int i=0;i<argc;i++)
      ksprintf(&str_cli," %s",argv[0]);

    mypars = pars_briggs(argc,argv);


    tsk_nthreads = mypars->nthread;
    my_tsk_struct = new tsk_struct[tsk_nthreads];
    assert(mypars);
    char* fname = mypars->hts;
    char* tabname = mypars->tab;
    const char* fastafile = mypars->ref;
    const char* lenname = mypars->len;
    const char* bedname = mypars->bed;
    const char* chromname = mypars->chr;
    char* model = mypars->model; // takes a value of b or nb.
    Contam_eps = mypars->eps; // Modern contamination rate \in [0,1)
    int len_limit = 150;
    
    char* bdamage = mypars->bdamage;
    char* rlens = mypars->rlens;


    faidx_t *seq_ref = NULL;
    seq_ref = fai_load(fastafile);
    fprintf(stderr,"\t-> fasta load \n");
    
    int inference_type = 0;
    if ((fname != NULL && tabname != NULL && bdamage != NULL) || (fname == NULL && tabname == NULL && bdamage == NULL)){
        fprintf(stdout,"Please provide either a bamfile or a table!\n");
        return 0;
    }
    else if (fname != NULL){
        inference_type = 1;
    }
    else if(tabname != NULL && lenname != NULL){
        inference_type = 2;
    }
    else if(bdamage != NULL && rlens != NULL){
        inference_type = 3;
    }

    double max5=0;
    double max3=0;
    double maxall=0;
    int len_min;
    if(inference_type < 3){
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
        
        

        if (fname != NULL){
            fprintf(stderr,"Loading the bamfile\n");
            bamreader(fname,chromname,seq_ref,len_limit,len_min,Frag_len,Frag_freq,number,mm5p,mm3p);
            //cout<<"Minimum length is "<<len_min<<"\n";
        }else if(tabname != NULL && lenname != NULL){
	  tabreader(tabname,STRLENS,mm5p,mm3p);
            fprintf(stderr,"Loading the table file with MAXLENGTH %d\n",MAXLENGTH);
        }
        freqCT = (double*) malloc(2*MAXLENGTH * sizeof(double));
        freqGA = (double*) malloc(2*MAXLENGTH * sizeof(double));
        scaleCT = (double*) malloc(2*MAXLENGTH * sizeof(double));
        scaleGA = (double*) malloc(2*MAXLENGTH * sizeof(double));
        seqError = (double*) malloc(2*MAXLENGTH * sizeof(double));
        
        len_min = len_min > 30 ? len_min : 30;

        
        for (int i=0; i<MAXLENGTH;i++){
        // Double-check this part
        // scaleCT[i] = mm5p[i][5]+mm5p[i][7];
        scaleCT[i] = mm5p[i][4]+mm5p[i][5]+mm5p[i][6]+mm5p[i][7];
        // scaleCT[2*MAXLENGTH-1-i] = mm3p[i][5]+mm3p[i][7];
        scaleCT[2*MAXLENGTH-1-i] = mm3p[i][4]+mm3p[i][5]+mm3p[i][6]+mm3p[i][7];
        // scaleGA[i] = mm3p[i][8]+mm3p[i][10];
        scaleGA[i] = mm3p[i][8]+mm3p[i][9]+mm3p[i][10]+mm3p[i][11];
        // scaleGA[2*MAXLENGTH-1-i] = mm5p[i][8]+mm5p[i][10];
        scaleGA[2*MAXLENGTH-1-i] = mm5p[i][8]+mm5p[i][9]+mm5p[i][10]+mm5p[i][11];
        
        // Overall sequencing errors are position specific, estimated by 1 - [N(AA)+N(TT)]/[N(AA)+N(AC)+N(AG)+N(AT)+N(TA)+N(TC)+N(TG)+N(TT)]
        seqError[i] = 1 - (mm5p[i][0]+mm5p[i][15])/(mm5p[i][0]+mm5p[i][1]+mm5p[i][2]+mm5p[i][3]+mm5p[i][12]+mm5p[i][13]+mm5p[i][14]+mm5p[i][15]);
        seqError[2*MAXLENGTH-1-i] = 1 - (mm3p[i][0]+mm3p[i][15])/(mm3p[i][0]+mm3p[i][1]+mm3p[i][2]+mm3p[i][3]+mm3p[i][12]+mm3p[i][13]+mm3p[i][14]+mm3p[i][15]);
        //cout << "Error "<<seqError[i]<<" "<<seqError[2*MAXLENGTH-1-i]<<"\n";
        freqCT[i] = mm5p[i][7]/scaleCT[i];
        freqCT[2*MAXLENGTH-1-i] = mm3p[i][7]/scaleCT[2*MAXLENGTH-1-i];
        freqGA[i] = mm3p[i][8]/scaleGA[i];
        freqGA[2*MAXLENGTH-1-i] = mm5p[i][8]/scaleGA[2*MAXLENGTH-1-i];
        max5 = std::max(max5,std::max(scaleCT[i],scaleCT[2*MAXLENGTH-1-i]));
        max3 = std::max(max3,std::max(scaleGA[i],scaleGA[2*MAXLENGTH-1-i]));
        }
        maxall = std::max(max5,max3);
        //    cout<<"max5 is "<<max5<<", max3 is "<<max3<<".\n";
        for (int i=0; i<MAXLENGTH;i++){
            scaleCT[i] = scaleCT[i]/maxall;
            scaleGA[i] = scaleGA[i]/maxall;
            scaleCT[2*MAXLENGTH-1-i] = scaleCT[2*MAXLENGTH-1-i]/maxall;
            scaleGA[2*MAXLENGTH-1-i] = scaleGA[2*MAXLENGTH-1-i]/maxall;
            //cout<<scaleCT[i]<<" "<<scaleGA[i]<<"\n";
        }
        
        
        fprintf(stderr,"\nThe misincorporting matrix is as follows:\n");
        fprintf(stderr,"Dir.\tPos.\tFreqCT\tFreqGA\n");
        for (int i=0; i<5;i++){	// 
            fprintf(stderr,"5'\t%d\t",i+1);
            fprintf(stderr,"%f\t%f\n",freqCT[i],freqGA[2*MAXLENGTH-1-i]);
        }
        for (int i=0; i<5;i++){
            fprintf(stderr,"3'\t%d\t",i+1);
            fprintf(stderr,"%f\t%f\n",freqCT[2*MAXLENGTH-1-i],freqGA[i]);
        }
    }
    else if(inference_type == 3){
        std::cout << " bdamage file "<< bdamage << " rlens " << rlens << std::endl;
        std::map<int, mydataD> retmap = load_bdamage_full(bdamage,MAXLENGTH);
        fprintf(stderr, "\t-> loading bdamage file in %lu mismatch matrices read for %d base pairs\n", retmap.size(), MAXLENGTH);
        for (std::map<int, mydataD>::iterator it = retmap.begin(); it != retmap.end(); it++) {
            int taxid = it->first;
            mydataD md = it->second;
            if (it->second.nreads == 0){
                continue;
            }
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
            
            //tabreader(tabname);

            int numpos = MAXLENGTH*2+1;
            int numcolumn = 16;
            double ** Table = (double **) malloc(numpos*(sizeof(double *))); /*I allocate memory here.  If this function is called many times it may be better to move the memmory allocation out of this function*/
            for (int i=0; i<numpos; i++){
                Table[i]=(double *) malloc(numcolumn*(sizeof(double)));
            }  

            parse_bdamage_Data(md, Table, MAXLENGTH);
            
            for (int i=0; i<MAXLENGTH; i++){
                for (int j=0; j<MAXLENGTH;j++){
                    mm5p[i][j] = Table[i][j];
                    mm3p[i][j] = Table[i+MAXLENGTH][j];
                }
            }

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
    }
    
    //Inference
    clock_t t=clock();
    time_t t2=time(NULL);
    fprintf(stderr,"\nThe inference is under progress...\n");
    //double *x = new double[3];
    double lbd[3] = {0,0,0};
    double ubd[3] = {1,1,1};
    int nbd[3] = {2,2,2};
    double invec1[3] = {0.2,0.1,0.1};
    ncalls=0;
    ncalls_grad = 0;
    wrapOne wo;
    wo.freqCT = freqCT;
    wo.freqGA = freqGA;
    wo.scaleCT = scaleCT;
    wo.scaleGA = scaleGA;
    wo.seqError = seqError;
    double withgrad = findmax_bfgs(3,invec1,(void *)&wo,b_loglike,b_loglike_grad,lbd,ubd,nbd,-1);
    
    double lbd1[4] = {1e-8,1e-8,1e-8,1e-8};
    double ubd1[4] = {1-1e-8,1-1e-8,1-1e-8,1-1e-8};
    int nbd1[4] = {2,2,2,2};
    double **z2;
    z2 = (double**) malloc(4 * sizeof(double*));
    for (int i = 0; i < 4; i++){
        z2[i] =(double *) malloc(4 * sizeof(double));
    }
    
    if(inference_type < 3){
        if (lenname != NULL){
            fprintf(stderr,"The provided fragment length distribution file is used!\n");
	    //len_limit
	    //number of different bins
	    //frag_len is really the fragment length
	    //Frag_freq is the relative frequency of frag_len
	    //const char * containing filename (-len input)
            FragArrayReader(len_limit, number, Frag_len, Frag_freq, lenname,STRLENS);
        }
        else{
            fprintf(stderr,"The fragment length distribution is calculated from the bam file!\n");
        }
    }
    else if(inference_type == 3){//this input is directly from metadamage, the rlens.gz
        FragArrayReaderRlen(len_limit,number,Frag_len,Frag_freq,rlens);
   
    }
    
    //BinNum = -1;
    FragArrayBin(number, BinNum, Frag_len, Frag_freq, Bin_Frag_len, Bin_Frag_freq);
    //    invec2 = {lambda,deltad,deltas,nu};
    double invec2[4] = {invec1[0],invec1[1],invec1[2],0.01};
    ncalls=0;
    ncalls_grad = 0;
    double withgrad2 = 0;
    wo.BinNum = BinNum;
    wo.Bin_Frag_len = Bin_Frag_len;
    wo.Bin_Frag_freq = Bin_Frag_freq;
    wo.Contam_eps = Contam_eps;
    
    if (model != NULL){
        if (!strcasecmp("b",model)){
	  fprintf(stderr,"%s\n","The chosen model is a biotin model.");
	  withgrad2 = findmax_bfgs(4,invec2,(void *) &wo,b_loglike_complex3_full,b_loglike_complex3_grad_full,lbd1,ubd1,nbd1,-1);
            loglike_complex3_hessian_full_b(invec2, z2, freqCT, freqGA, scaleCT, scaleGA, seqError, BinNum, Bin_Frag_len, Bin_Frag_freq, Contam_eps);
        }else if(!strcasecmp("nb",model)){
	  fprintf(stderr,"%s\n","The chosen model is non-biotin model.");
            withgrad2 = findmax_bfgs(4,invec2,(void *) &wo,nb_loglike_complex3_full,nb_loglike_complex3_grad_full,lbd1,ubd1,nbd1,-1);
            loglike_complex3_hessian_full_nb(invec2, z2, freqCT, freqGA, scaleCT, scaleGA, seqError, BinNum, Bin_Frag_len, Bin_Frag_freq, Contam_eps);
        }else{
            fprintf(stderr,"Please provide a meaningful deamination model for further calculations.\n");
            return -1;
        }
    }else{
        fprintf(stderr,"Please provide a deamination model for further calculations.\n");
        return -1;
    }
    double stdvec2[4];
    for (int i=0;i<4;i++){
        stdvec2[i] = sqrt(z2[i][i]/maxall);
    }
    
    fprintf(stderr,
            "\n"
            "\t[ALL done] cpu-time used =  %.4f sec\n"
            "\t[ALL done] walltime used =  %.4f sec\n"
            ,(float)(clock() - t) / CLOCKS_PER_SEC, (float)(time(NULL) - t2));
    if (model != NULL){
        if (!strcasecmp("b",model)){
            fprintf(stderr,"%s","The chosen model is a biotin model, ");
        }else if(!strcasecmp("nb",model)){
            fprintf(stderr,"%s","The chosen model is non-biotin model, ");
        }else{
            fprintf(stderr,"Please provide a meaningful deamination model for further calculations.\n");
            return -1;
        }
    }else{
        fprintf(stderr,"Please provide a deamination model for further calculations.\n");
        return -1;
    }
    fprintf(stderr,"the inferred parameters are given as follows:\n");
    fprintf(stderr,"lambda: %f (%f,%f), delta: %f (%f,%f), delta_s: %f (%f,%f), nu: %f (%f,%f), nfunctioncalls: %d ngradientcalls: %d, llh: %f.\n", invec2[0], invec2[0]-1.96*stdvec2[0], invec2[0]+1.96*stdvec2[0] , invec2[1], invec2[1]-1.96*stdvec2[1], invec2[1]+1.96*stdvec2[1], invec2[2], invec2[2]-1.96*stdvec2[2], invec2[2]+1.96*stdvec2[2], invec2[3], invec2[3]-1.96*stdvec2[3], invec2[3]+1.96*stdvec2[3], ncalls,ncalls_grad, withgrad2);
    
    // Output count table
    char *otabname = mypars->otab;
    if (otabname!=NULL){
        BGZF *fp = NULL;
        fp = bgzf_open(otabname,"wb");
        kstring_t *kstr = new kstring_t;
        kstr->s = NULL;
        kstr->l = kstr->m = 0;
        ksprintf(kstr,"%s\n","The misincorporting matrix is as follows:");
        ksprintf(kstr,"%s\n","Dir.\tPos.\tAA\tAC\tAG\tAT\tCA\tCC\tCG\tCT\tGA\tGC\tGG\tGT\tTA\tTC\tTG\tTT");
        for (int i=0; i<MAXLENGTH;i++){
            ksprintf(kstr,"%s\t%d\t","5'",i+1);
            for (int j=0; j<16;j++){
                ksprintf(kstr,"%f\t",mm5p[i][j]);
            }
            ksprintf(kstr,"%f\t%f\n",freqCT[i],freqGA[2*MAXLENGTH-1-i]);
        }
        for (int i=0; i<MAXLENGTH;i++){
            ksprintf(kstr,"%s\t%d\t","3'",i+1);
            for (int j=0; j<16;j++){
                ksprintf(kstr,"%f\t",mm3p[i][j]);
            }
            ksprintf(kstr,"%f\t%f\n",freqCT[2*MAXLENGTH-1-i],freqGA[i]);
        }
        //my_bgzf_write(fp,kstr->s,kstr->l);
        assert(bgzf_write(fp,kstr->s,kstr->l)==kstr->l);

        kstr->l = 0;
        bgzf_close(fp);
    }
    
    //Output length distribution file
    char *olenname = mypars->olen;
    if (olenname!=NULL){
        BGZF *fp = NULL;
        fp = bgzf_open(olenname,"wb");
        kstring_t *kstr = new kstring_t;
        kstr->s = NULL;
        kstr->l = kstr->m = 0;
        ksprintf(kstr,"%s\n","The raw length distribution is as follows:");
        ksprintf(kstr,"%s%d\n","The number of length values is ",number);
        for (int k=0;k<number;k++){
            ksprintf(kstr,"%d\t%f\n",Frag_len[k],Frag_freq[k]);
        }
        ksprintf(kstr,"%s\n","The binned length distribution is as follows:");
        ksprintf(kstr,"%s%d\n","The number of binned length values is ",BinNum);
        for (int k=0;k<BinNum;k++){
            ksprintf(kstr,"%f\t%f\n",Bin_Frag_len[k],Bin_Frag_freq[k]);
        }
        //my_bgzf_write(fp,kstr->s,kstr->l);
        assert(bgzf_write(fp,kstr->s,kstr->l)==kstr->l);
        kstr->l = 0;
        bgzf_close(fp);
    }
    
    // Output inferred parameters
    char* oinfname = mypars->oinf;
    if (oinfname!=NULL){
        BGZF *fp = bgzf_open(oinfname,"wb");
        kstring_t *kstr = new kstring_t;
        kstr->s = NULL;
        kstr->l = kstr->m = 0;
        if (model != NULL){
            if (!strcasecmp("b",model)){
                ksprintf(kstr,"%s","The chosen model is a biotin model, ");
            }else if(!strcasecmp("nb",model)){
                ksprintf(kstr,"%s","The chosen model is non-biotin model, ");
            }else{
                fprintf(stderr,"Please provide a meaningful deamination model for further calculations.\n");
                return -1;
            }
        }else{
            fprintf(stderr,"Please provide a deamination model for further calculations.\n");
            return -1;
        }
        ksprintf(kstr,"the inferred parameters are given as follows:\n");
        ksprintf(kstr,"lambda: %f (%f,%f), delta: %f (%f,%f), delta_s: %f (%f,%f), nu: %f (%f,%f), nfunctioncalls: %d ngradientcalls: %d, llh: %f.\n", invec2[0], invec2[0]-1.96*stdvec2[0], invec2[0]+1.96*stdvec2[0] , invec2[1], invec2[1]-1.96*stdvec2[1], invec2[1]+1.96*stdvec2[1], invec2[2], invec2[2]-1.96*stdvec2[2], invec2[2]+1.96*stdvec2[2], invec2[3], invec2[3]-1.96*stdvec2[3], invec2[3]+1.96*stdvec2[3], ncalls,ncalls_grad, withgrad2);
	//my_bgzf_write(fp,kstr->s,kstr->l);
        assert(bgzf_write(fp,kstr->s,kstr->l)==kstr->l);
        kstr->l = 0;
        bgzf_close(fp);
    }
    
    char* refName;
    refName = NULL;
    int mapq =-1;
    int mapped_only = 0;
    int se_only = 1;
    len_limit = 150;
    char *fname1 = mypars->ihts; //the file needs to be recalibrated
    char *olik1 = mypars->olik; //the output nucleotide likelihood file
    const char* chromname1 = mypars->ichr;
    const char* bedname1 = mypars->ibed;

    int isrecal = mypars->isrecal;
    
    deamRateCT = (double**) malloc((len_limit-30) * sizeof(double*));
    deamRateGA = (double**) malloc((len_limit-30) * sizeof(double*));
    for (int i = 0; i < (len_limit-30); i++){
        deamRateCT[i] =(double *) malloc(30 * sizeof(double));
        deamRateGA[i] =(double *) malloc(30 * sizeof(double));
    }
    if (isrecal>=0){
        if (!strcasecmp("b",model)){
	  CaldeamRate_b(invec2[0],invec2[1],invec2[2],invec2[3],len_limit,deamRateCT,deamRateGA);
        }else if (!strcasecmp("nb",model)){
	  CaldeamRate_nb(invec2[0],invec2[1],invec2[2],invec2[3],len_limit,deamRateCT,deamRateGA);
        }

    }

    double distparam[4] = {60,45,100,45};
    //    fprintf(stderr,"fname1: %s fname2: %s\n",fname1,fname2);
    
    sam_hdr_t *hdr;
    if (fname1!=NULL && isrecal == 1){
        //tsk start
        std::vector<bam1_t*> tsk_reads;
        hdr = read_all_reads(fname1,bedname1,refName,tsk_reads);
        //tsk stop
        ncalls = 0;
        ncalls_grad = 0;
        //double lbd2[4] = {30,1e-8,30,1e-8};
        double lbd2[4] = {30,15,30,1};
        double ubd2[4] = {(double)len_limit-1,100,(double)len_limit-1,100};
        int nbd2[4] = {2,2,2,2};
        int mynsites = tsk_reads.size()/tsk_nthreads;
        fprintf(stderr,"mynSites: %d\n",mynsites);
        for(int ii=0;ii<tsk_nthreads;ii++){
            my_tsk_struct[ii].from = my_tsk_struct[ii].to = -1;
            my_tsk_struct[ii].reads = &tsk_reads;
            my_tsk_struct[ii].hdr = hdr;
            my_tsk_struct[ii].seq_ref = seq_ref;
            my_tsk_struct[ii].len_limit = len_limit;
            my_tsk_struct[ii].len_min = len_min;
            my_tsk_struct[ii].model = model;
            my_tsk_struct[ii].eps = Contam_eps;
            my_tsk_struct[ii].lambda = invec2[0];
            my_tsk_struct[ii].delta = invec2[1];
            my_tsk_struct[ii].delta_s = invec2[2];
            my_tsk_struct[ii].nv = invec2[3];
            my_tsk_struct[ii].x = new double [4];
            my_tsk_struct[ii].llh_result_grad = new double [4];
            my_tsk_struct[ii].llh_result_hess = new double* [4];
	    my_tsk_struct[ii].Tol = Tol;
            for(int i=0;i<4;i++){
                my_tsk_struct[ii].x[i] = distparam[i];
                my_tsk_struct[ii].llh_result_hess[i] = new double [4];
            }
            my_tsk_struct[ii].threadid = ii;
            my_tsk_struct[ii].from = ii==0?0:my_tsk_struct[ii-1].to;
            my_tsk_struct[ii].to =   my_tsk_struct[ii].from+mynsites;
        }
        my_tsk_struct[tsk_nthreads-1].to = tsk_reads.size();
        double withgrad3;
        if(tsk_nthreads==1)
            withgrad3 = findmax_bfgs(4,distparam,&my_tsk_struct[0],tsk_All_loglike_recalibration,tsk_All_loglike_recalibration_grad,lbd2,ubd2,nbd2,-1);
        else{
            //withgrad3 = findmax_bfgs(4,distparam,NULL,like_master,tsk_All_loglike_recalibration_grad,lbd2,ubd2,nbd2,-1);
            withgrad3 = findmax_bfgs(4,distparam,NULL,like_master,like_grad_master,lbd2,ubd2,nbd2,-1);
        }
        double **covpar = new double* [4];
        for (int i=0;i<4;i++){
            covpar[i] = new double[4];
        }
        like_hess_master(distparam,covpar);
        double stdpar[4];
        for (int i=0;i<4;i++){
            stdpar[i] = sqrt(covpar[i][i]);
            //cout<<z4[i][i]<<" std "<<stdvec4[i]<<"\n";
        }
        fprintf(stderr,"mu_anc: %f (%f,%f), sigma_anc: %f (%f,%f), mu_mod: %f (%f,%f), sigma_mod: %f (%f,%f), nfunctioncalls: %d ngradientcalls: %d, llh: %f.\n", distparam[0], distparam[0]-1.96*stdpar[0], distparam[0]+1.96*stdpar[0], distparam[1], distparam[1]-1.96*stdpar[1], distparam[1]+1.96*stdpar[1], distparam[2], distparam[2]-1.96*stdpar[2], distparam[2]+1.96*stdpar[2], distparam[3], distparam[3]-1.96*stdpar[3], distparam[3]+1.96*stdpar[3], ncalls,ncalls_grad, withgrad3);
        sam_hdr_destroy(hdr);
    }
    //fprintf(stderr,"BRANCH 2000!\n");fflush(stderr);
    if ((fname1==NULL) && (olik1!=NULL)){
        fprintf(stdout,"Please provide bam file to calculate nucleotide likelihoods!\n");
        return 0;
    }
    
    kstring_t *kstr3 = new kstring_t;
    kstr3->s = NULL;
    kstr3->l = kstr3->m = 0;
 
    if (fname1!=NULL){
      bam_hdr_t* hdr = calc_pp_pmd_prob(refName,fname1,mypars->ohts,olik1,mapped_only,se_only,mapq,seq_ref,len_limit,len_min,model,Contam_eps,invec2[0],invec2[1],invec2[2],invec2[3],distparam[0],distparam[1],distparam[2], distparam[3], isrecal,&str_cli,deamRateCT,deamRateGA,Tol);
      sam_hdr_destroy(hdr);
    }
    
    free(Frag_len);
    free(Frag_freq);
    free(Bin_Frag_len);
    free(Bin_Frag_freq);
    
    for (int i = 0; i < 4; i++){
        free(z2[i]);
    }
    free(z2);
    for (int i = 0; i < MAXLENGTH; i++){
        free(mm5p[i]);
        free(mm3p[i]);
    }
    for (int i = 0; i < (len_limit-30); i++){
        free(deamRateCT[i]);
        free(deamRateGA[i]);
    }
    free(deamRateCT);
    free(deamRateGA);
    free(seqError);
    free(freqCT);
    free(freqGA);
    free(scaleCT);
    free(scaleGA);
    free(mm5p);
    free(mm3p);
    argStruct_destroy(mypars);

    // check deallocation of my_tsk_struct values
    delete my_tsk_struct;

//    if(hdr){
//        sam_hdr_destroy(hdr);
//    }
    return 0;
}


