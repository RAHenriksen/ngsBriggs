

#include <cstring>
#include <cstdlib>
#include <cstdio>
#include <cassert>
#include <fstream>
#include <sstream>
#include <iostream>

#include <zlib.h>
#include <cmath>

#include <ctime>



#include <htslib/hts.h>
#include <htslib/sam.h>
#include <htslib/kstring.h>
#include <htslib/faidx.h>
#include <htslib/thread_pool.h>

#include "profile.h"
#include "bfgs.h"
#include "htslib/bgzf.h"
#include <array>  //Rasmus

#include "misc.h"
#include "Recalibration.h"
#include "ngsBriggs.h"
#include "Likelihood.h"
#include "PosteriorProb.h"

typedef unsigned char uchar;


double NoPMDGivenAnc_b(char reffrag[], char frag[], int L, double lambda, double delta, double delta_s, double nv,double Tol){
    double np_pmd_anc=0; // probability of no pmd given ancient
    double p = 0;
    //Investigate each possible left and right 5' overhang pair with (l,r) <= l_check (15)
    for (int l=0; l<=l_check; l++){
        double pl = 0.5*lambda*pow(1-lambda,l)/(1-pow(1-lambda,L-1));
        if (l==0){
            pl += 0.5;
        }
        for (int r=0; r<=l_check; r++){
            double pr = 0.5*lambda*pow(1-lambda,r)/(1-pow(1-lambda,L-1));
            if (r == 0){
                pr += 0.5;
            }
            // Joint probability of (l,r)
            double dp = 4*pl*pr/(3+(1-pow(1-lambda,L-1)-(L-1)*lambda*pow(1-lambda,L-1))/(pow(1-pow(1-lambda,L-1),2)));
            if (dp > Tol){
                // Probability
                double pnick = 0;
                for (int npos=l; npos<=l_check-1; npos++){ //Nick is within left l_check region
                    double l2 = 0;
                    for (int i=0; i<l_check; i++){
                        if (reffrag[i]<4){
                            if (reffrag[i]==1 && i<=l-1){// C to C, left single strand part
                                l2 += log(1-delta_s);
                            }else if (reffrag[i]==2 && i<=l-1){ // G to G, left single strand part
                                l2 += 0;
                            }else if (reffrag[i]==1 && i<=npos){ // C to C, left double strand part before nick
                                l2 += log(1-delta);
                            }else if (reffrag[i]==2 && i<=npos){ // G to G, left double strand part before nick
                                l2 += 0;
                            }else if (reffrag[i]==1 && i>npos){ // C to C, left double strand part after nick
                                l2 += 0;
                            }else if (reffrag[i]==2 && i>npos){ // G to G, left double strand part after nick
                                l2 += log(1-delta);
                            }else if (reffrag[i]==frag[i]){ // All other possible no changes
                                l2 += 0;
                            }else if (reffrag[i]!=frag[i]){ // All other possible changes
                                l2 += 0;
                            }
                        }
                        if (reffrag[L-i-1]<4){
                            if (reffrag[L-i-1]==2 && i <= r-1){ // G to G, right single strand part
                                l2 += log(1-delta_s);
                            }else if(reffrag[L-i-1]==1 && i <= r-1){ // C to C, right single strand part
                                l2 += 0;
                            }else if (reffrag[L-i-1]==2){ // G to G, right double strand part
                                l2 += log(1-delta);
                            }else if (reffrag[L-i-1]==1){ // C to C, right double strand part
                                l2 += 0;
                            }else if (reffrag[L-i-1]==frag[L-i-1]){ // All other possible no changes
                                l2 += 0;
                            }else if (reffrag[L-i-1]!=frag[L-i-1]){ // All other changes
                                l2 += 0;
                            }
                        }
                    }
                    double dpnick = nv/(1+(L-l-r-2)*nv);
                    np_pmd_anc += exp(l2)*dp*dpnick;
                    pnick += dpnick;
                }
                
                for (int npos = L-l_check; npos <= L-r-1; npos++){ //Nick is within right l_check region
                    double l2 = 0;
                    for (int i=0; i<l_check; i++){
                        //double error2 = PhredError[seqError[L-i-1]];
                        //double error2AThird = PhredErrorAThird[seqError[L-i-1]];
                        if (reffrag[L-i-1]<4){
                            if (reffrag[L-i-1]==2 && i <= r-1){ // G to G, right single strand part
                                l2 += log(1-delta_s);
                            }else if(reffrag[L-i-1]==1 && i <= r-1){ // C to C, right single strand part
                                l2 += 0;
                            }else if (reffrag[L-i-1]==2 && i <= L-npos-1){ // G to G, right double strand part, right of the nick
                                l2 += log(1-delta);
                            }else if (reffrag[L-i-1]==1 && i <= L-npos-1){ // C to C, right double strand part, right of the nick
                                l2 += 0;
                            }else if (reffrag[L-i-1]==2 && i > L-npos-1){ // G to G, right double strand part, left of the nick
                                l2 += 0;
                            }else if (reffrag[L-i-1]==1 && i > L-npos-1){ // C to C, right double strand part, left of the nick
                                l2 += log(1-delta);
                            }else if (reffrag[L-i-1]==frag[L-i-1]){ // All other possible no changes
                                l2 += 0;
                            }else if (reffrag[L-i-1]!=frag[L-i-1]){ // All other possible changes
                                l2 += 0;
                            }
                        }
                        
                        // double error1 = PhredError[seqError[i]];
                        // double error1AThird = PhredErrorAThird[seqError[i]];
                        if (reffrag[i]<4){
                            if (reffrag[i]==1 && i<=l-1){ // C to C, left single strand part
                                l2 += log(1-delta_s);
                            }else if (reffrag[i]==2 && i<=l-1){ // G to G, left single strand part
                                l2 += 0;
                            }else if (reffrag[i]==1){ // C to C, left double strand part
                                l2 += log(1-delta);
                            }else if (reffrag[i]==2){ // G to G, left double strand part
                                l2 += 0;
                            }else if (reffrag[i]==frag[i]){ // All other possible no changes
                                l2 += 0;
                            }else if (reffrag[i]!=frag[i]){ // All other possible changes
                                l2 += 0;
                            }
                        }
                    }
                    double dpnick = npos < L-r-1 ? nv/(1+(L-l-r-2)*nv) : (1-nv)/(1+(L-l-r-2)*nv);
                    np_pmd_anc += exp(l2)*dp*dpnick;
                    pnick += dpnick;
                }
                
                //nick occurring outside both l_check region
                double l2 = 0;
                for (int i=0; i<l_check; i++){
                    // double error1 = PhredError[seqError[i]];
                    // double error1AThird = PhredErrorAThird[seqError[i]];
                    if (reffrag[i]<4){
                        if (reffrag[i]==1 && i<=l-1){ // C to C, left single strand part
                            l2 += log(1-delta_s);
                        }else if (reffrag[i]==2 && i<=l-1){ // G to G, left single strand part
                            l2 += 0;
                        }else if (reffrag[i]==1 && i>l-1){ // C to C, left double strand part
                            l2 += log(1-delta);
                        }else if (reffrag[i]==2 && i>l-1){ // G to G, left double strand part
                            l2 += 0;
                        }else if (reffrag[i]==frag[i]){ // All other possible no changes
                            l2 += 0;
                        }else if (reffrag[i]!=frag[i]){ // All other possible changes
                            l2 += 0;
                        }
                    }
                    
                    // double error2 = PhredError[seqError[L-i-1]];
                    // double error2AThird = PhredErrorAThird[seqError[L-i-1]];
                    if (reffrag[L-i-1]<4){
                        if (reffrag[L-i-1]==2 && i <= r-1){ // G to G, right single strand part
                            l2 += log(1-delta_s);
                        }else if(reffrag[L-i-1]==1 && i <= r-1){ // C to C, right single strand part
                            l2 += 0;
                        }else if (reffrag[L-i-1]==2 && i > r-1){ // G to G, right double strand part
                            l2 += log(1-delta);
                        }else if (reffrag[L-i-1]==1 && i > r-1){ // C to T change, right double strand part
                            l2 += 0;
                        }else if (reffrag[L-i-1]==frag[L-i-1]){ // All other possible no changes
                            l2 += 0;
                        }else if (reffrag[L-i-1]!=frag[L-i-1]){ // All other possible changes
                            l2 += 0;
                        }
                    }
                }
                np_pmd_anc += exp(l2)*dp*(1-pnick);
            }
            p += dp;
        }
        
        for (int r=l_check+1; r<=L-2-l; r++){ // To guarantee l+r < = L-2
            double pr = 0.5*lambda*pow(1-lambda,r)/(1-pow(1-lambda,L-1));
            if (r == 0){
                pr += 0.5;
            }
            double dp = 4*pl*pr/(3+(1-pow(1-lambda,L-1)-(L-1)*lambda*pow(1-lambda,L-1))/(pow(1-pow(1-lambda,L-1),2)));
            if (dp >= Tol){
                double pnick = 0;
                for (int npos=l; npos<=(l_check-1 < L-r-1 ? l_check-1 : L-r-1); npos++){ //Nick is within left l_check region
                    double l2 = 0;
                    for (int i=0; i<l_check; i++){
                        // double error1 = PhredError[seqError[i]];
                        // double error1AThird = PhredErrorAThird[seqError[i]];
                        if (reffrag[i]<4){
                            if (reffrag[i]==1 && i<=l-1){// C to C, left single strand part
                                l2 += log(1-delta_s);
                            }else if (reffrag[i]==2 && i<=l-1){ // G to G, left single strand part
                                l2 += 0;
                            }else if (reffrag[i]==1 && i<=npos){ // C to C, left double strand part before nick
                                l2 += log(1-delta);
                            }else if (reffrag[i]==2 && i<=npos){ // G to G, left double strand part before nick
                                l2 += 0;
                            }else if (reffrag[i]==1 && i>npos){ // C to C, left double strand part after nick
                                l2 += 0;
                            }else if (reffrag[i]==2 && i>npos){ // G to A change, left double strand part after nick
                                l2 += log(1-delta);
                            }else if (reffrag[i]==frag[i]){ // All other possible no changes
                                l2 += 0;
                            }else if (reffrag[i]!=frag[i]){ // All other possible changes
                                l2 += 0;
                            }
                        }
                        
                        // double error2 = PhredError[seqError[L-i-1]];
                        // double error2AThird = PhredErrorAThird[seqError[L-i-1]];
                        if (reffrag[L-i-1]<4){
                            if (reffrag[L-i-1]==2){ // G to G, right single strand part
                                l2 += log(1-delta_s);
                            }else if(reffrag[L-i-1]==1){ // C to C, right single strand part
                                l2 += 0;
                            }else if (reffrag[L-i-1]==frag[L-i-1]){ // All other possible no changes
                                l2 += 0;
                            }else if (reffrag[L-i-1]!=frag[L-i-1]){ // All other changes
                                l2 += 0;
                            }
                        }
                    }
                    double dpnick = nv/(1+(L-l-r-2)*nv);
                    np_pmd_anc += exp(l2)*dp*dpnick;
                    pnick += dpnick;
                }
                //nick occurring outside both l_check region
                double l2 = 0;
                for (int i=0; i<l_check; i++){
                    // double error1 = PhredError[seqError[i]];
                    // double error1AThird = PhredErrorAThird[seqError[i]];
                    if (reffrag[i]<4){
                        if (reffrag[i]==1 && i<=l-1){ // C to C, left single strand part
                            l2 += log(1-delta_s);
                        }else if (reffrag[i]==2 && i<=l-1){ // G to G, left single strand part
                            l2 += 0;
                        }else if (reffrag[i]==1 && i>l-1){ // C to C, left double strand part
                            l2 += log(1-delta);
                        }else if (reffrag[i]==2 && i>l-1){ // G to G, left double strand part
                            l2 += 0;
                        }else if (reffrag[i]==frag[i]){ // All other possible no changes
                            l2 += 0;
                        }else if (reffrag[i]!=frag[i]){ // All other possible changes
                            l2 += 0;
                        }
                    }
                    
                    // double error2 = PhredError[seqError[L-i-1]];
                    // double error2AThird = PhredErrorAThird[seqError[L-i-1]];
                    if (reffrag[L-i-1]<4){
                        if (reffrag[L-i-1]==2 && i <= r-1){ // G to G, right single strand part
                            l2 += log(1-delta_s);
                        }else if(reffrag[L-i-1]==1 && i <= r-1){ // C to C, right single strand part
                            l2 += 0;
                        }else if (reffrag[L-i-1]==frag[L-i-1]){ // All other possible no changes
                            l2 += 0;
                        }else if (reffrag[L-i-1]!=frag[L-i-1]){ // All other possible changes
                            l2 += 0;
                        }
                    }
                }
                np_pmd_anc += exp(l2)*dp*(1-pnick);
            }
            p += dp;
        }
    }
    
    for (int r=0; r<=l_check; r++){
        double pr = 0.5*lambda*pow(1-lambda,r)/(1-pow(1-lambda,L-1));
        if (r == 0){
            pr += 0.5;
        }
        for (int l=l_check+1;l<=L-2-r; l++){
            double pl = 0.5*lambda*pow(1-lambda,l)/(1-pow(1-lambda,L-1));
            if (l==0){
                pl += 0.5;
            }
            // Joint probability of (l,r)
            double dp = 4*pl*pr/(3+(1-pow(1-lambda,L-1)-(L-1)*lambda*pow(1-lambda,L-1))/(pow(1-pow(1-lambda,L-1),2)));
            if (dp > Tol){
                double pnick = 0;
                
                for (int npos = (L-l_check-1 > l ? L-l_check-1 : l); npos <= L-r-1; npos++){ //Nick is within right l_check region
                    double l2 = 0;
                    for (int i=0; i<l_check; i++){
                        //double error2 = PhredError[seqError[L-i-1]];
                        //double error2AThird = PhredErrorAThird[seqError[L-i-1]];
                        if (reffrag[L-i-1]<4){
                            if (reffrag[L-i-1]==2 && i <= r-1){ // G to G, right single strand part
                                l2 += log(1-delta_s);
                            }else if(reffrag[L-i-1]==1 && i <= r-1){ // C to C, right single strand part
                                l2 += 0;
                            }else if (reffrag[L-i-1]==2 && i <= L-npos-1){ // G to G, right double strand part, right of the nick
                                l2 += log(1-delta);
                            }else if (reffrag[L-i-1]==1 && i <= L-npos-1){ // C to C, right double strand part, right of the nick
                                l2 += 0;
                            }else if (reffrag[L-i-1]==2 && i > L-npos-1){ // G to G, right double strand part, left of the nick
                                l2 += 0;
                            }else if (reffrag[L-i-1]==1 && i > L-npos-1){ // C to C, right double strand part, left of the nick
                                l2 += log(1-delta);
                            }else if (reffrag[L-i-1]==frag[L-i-1]){ // All other possible no changes
                                l2 += 0;
                            }else if (reffrag[L-i-1]!=frag[L-i-1]){ // All other possible changes
                                l2 += 0;
                            }
                        }
                        
                        // double error1 = PhredError[seqError[i]];
                        // double error1AThird = PhredErrorAThird[seqError[i]];
                        if (reffrag[i]<4){
                            if (reffrag[i]==1){ // C to C, left single strand part
                                l2 += log(1-delta_s);
                            }else if (reffrag[i]==2){ // G to G, left single strand part
                                l2 += 0;
                            }else if (reffrag[i]==frag[i]){ // All other possible no changes
                                l2 += 0;
                            }else if (reffrag[i]!=frag[i]){ // All other possible changes
                                l2 += 0;
                            }
                        }
                    }
                    double dpnick = npos < L-r-1 ? nv/(1+(L-l-r-2)*nv) : (1-nv)/(1+(L-l-r-2)*nv);
                    np_pmd_anc += exp(l2)*dp*dpnick;
                    pnick += dpnick;
                }
                // nick occurring outside both l_check region
                double l2 = 0;
                for (int i=0; i<l_check; i++){
                    // double error1 = PhredError[seqError[i]];
                    // double error1AThird = PhredErrorAThird[seqError[i]];
                    if (reffrag[i]<4){
                        if (reffrag[i]==1){ // C to C, left single strand part
                            l2 += log(1-delta_s);
                        }else if (reffrag[i]==2){ // G to G, left single strand part
                            l2 += 0;
                        }else if (reffrag[i]==frag[i]){ // All other possible no changes
                            l2 += 0;
                        }else if (reffrag[i]!=frag[i]){ // All other possible changes
                            l2 += 0;
                        }
                    }
                    
                    // double error2 = PhredError[seqError[L-i-1]];
                    // double error2AThird = PhredErrorAThird[seqError[L-i-1]];
                    if (reffrag[L-i-1]<4){
                        if (reffrag[L-i-1]==2 && i <= r-1){ // G to G, right single strand part
                            l2 += log(1-delta_s);
                        }else if(reffrag[L-i-1]==1 && i <= r-1){ // C to C, right single strand part
                            l2 += 0;
                        }else if (reffrag[L-i-1]==2 && i > r-1){ // G to G, right double strand part
                            l2 += log(1-delta);
                        }else if (reffrag[L-i-1]==1 && i > r-1){ // C to T change, right double strand part
                            l2 += 0;
                        }else if (reffrag[L-i-1]==frag[L-i-1]){ // All other possible no changes
                            l2 += 0;
                        }else if (reffrag[L-i-1]!=frag[L-i-1]){ // All other possible changes
                            l2 += 0;
                        }
                    }
                }
                np_pmd_anc += exp(l2)*dp*(1-pnick);
            }
            p += dp;
        }
    }
    
    //Investigate each possible left and right 5' overhang pair with (l,r) > l_check (15)
    if (1-p > Tol){
        double l2 = 0;
        for (int i=0; i<l_check; i++){
            // double error1 = PhredError[seqError[i]];
            // double error1AThird = PhredErrorAThird[seqError[i]];
            if (reffrag[i]<4){
                if (reffrag[i]==1){ // C to C, left single strand part
                    l2 += log(1-delta_s);
                }else if (reffrag[i]==2){ // G to G, left single strand part
                    l2 += 0;
                }else if (reffrag[i]==frag[i]){ // All other possible no changes
                    l2 += 0;
                }else if (reffrag[i]!=frag[i]){  // All other possible changes
                    l2 += 0;
                }
            }
            
            // double error2 = PhredError[seqError[L-i-1]];
            // double error2AThird = PhredErrorAThird[seqError[L-i-1]];
            if (reffrag[L-i-1]<4){
                if (reffrag[L-i-1]==2){ // G to G, right single strand part
                    l2 += log(1-delta_s);
                }else if(reffrag[L-i-1]==1){ // C to C, right single strand part
                    l2 += 0;
                }else if (reffrag[L-i-1]==frag[L-i-1]){ // All other possible no changes
                    l2 += 0;
                }else if (reffrag[L-i-1]!=frag[L-i-1]){ // All other possible changes
                    l2 += 0;
                }
            }
        }
        np_pmd_anc += exp(l2)*(1-p);
    }
    return np_pmd_anc;
}


double NoPMDGivenAnc_nb(char reffrag[], char frag[], int L, double lambda, double delta, double delta_s, double nv,double Tol){
    double np_pmd_anc=0; // Likelihood
    double p = 0; // Accumulated prob of (l,r)
    //Investigate each possible left and right 5' overhang pair with (l,r) <= l_check (15)
    //--------- - - - -                        - - - ------------// --- single strand part
    for (int l=0; l<=l_check; l++){
        double pl = 0.5*lambda*pow(1-lambda,l)/(1-pow(1-lambda,L-1));
        if (l==0){
            pl += 0.5;
        }
        for (int r=0; r<=l_check; r++){
            double pr = 0.5*lambda*pow(1-lambda,r)/(1-pow(1-lambda,L-1));
            if (r == 0){
                pr += 0.5;
            }
            // Joint probability of (l,r)
            double dp = 4*pl*pr/(3+(1-pow(1-lambda,L-1)-(L-1)*lambda*pow(1-lambda,L-1))/(pow(1-pow(1-lambda,L-1),2)));
            if (dp > Tol){
                double pnick = 0; // Probability
                //--------- - - - -                        - - - ------------//
                //            /|\                                            //
                //             |___ Nick                                     //
                for (int npos=l; npos<=l_check-1; npos++){ //Nick is within left l_check region
                    double l2 = 0;
                    for (int i=0; i<l_check; i++){
                        //double error1 = PhredError[seqError[L-i-1]]; //rev-comp but error rate should be the same
                        //double error1AThird = PhredErrorAThird[seqError[L-i-1]]; //rev-comp
                        if (reffrag[L-i-1]<4){
                            if (reffrag[L-i-1]==2 && i<=l-1){// C to C, left single strand part, rev-comp
                                l2 += log(1-delta_s);
                            }else if (reffrag[L-i-1]==1 && i<=l-1){ // G to G, left single strand part, rev-comp
                                l2 += 0;
                            }else if (reffrag[L-i-1]==2 && i<=npos){ // C to C, left double strand part before nick, rev-comp
                                l2 += log(1-delta);
                            }else if (reffrag[L-i-1]==1 && i<=npos){ // G to G, left double strand part before nick, rev-comp
                                l2 += 0;
                            }else if (reffrag[L-i-1]==2 && i>npos){ // C to C, left double strand part after nick, rev-comp
                                l2 += 0;
                            }else if (reffrag[L-i-1]==1 && i>npos){ // G to G, left double strand part after nick, rev-comp
                                l2 += log(1-delta);
                            }else if (reffrag[L-i-1]==frag[L-i-1]){ // All other possible no changes, rev-comp
                                l2 += 0;
                            }else if (reffrag[L-i-1]!=frag[L-i-1]){ // All other possible changes, rev-comp
                                l2 += 0;
                            }
                        }
                        
                        //double error2 = PhredError[seqError[i]];
                        //double error2AThird = PhredErrorAThird[seqError[i]];
                        if (reffrag[i]<4){
                            if (reffrag[i]==1 && i <= r-1){ // G to G, right single strand part, rev-comp
                                l2 += log(1-delta_s);
                            }else if(reffrag[i]==2 && i <= r-1){ // C to C, right single strand part, rev-comp
                                l2 += 0;
                            }else if (reffrag[i]==1){ // G to G, right double strand part, rev-comp
                                l2 += log(1-delta);
                            }else if (reffrag[i]==2){ // C to C, right double strand part, rev-comp
                                l2 += 0;
                            }else if (reffrag[i]==frag[i]){ // All other possible no changes, rev-comp
                                l2 += 0;
                            }else if (reffrag[i]!=frag[i]){ // All other changes, rev-comp
                                l2 += 0;
                            }
                        }
                    }
                    double dpnick = nv/(1+(L-l-r-2)*nv);
                    np_pmd_anc += exp(l2)*dp*dpnick;
                    pnick += dpnick;
                }
                
                //--------- - - - -                        - - - ------------//
                //                                          /|\              //
                //                                           |___ Nick       //
                for (int npos = L-l_check; npos <= L-r-1; npos++){ //Nick is within right l_check region
                    double l2 = 0;
                    for (int i=0; i<l_check; i++){
                        //double error2 = PhredError[seqError[i]];
                        //double error2AThird = PhredErrorAThird[seqError[i]];
                        if (reffrag[i]<4){
                            if (reffrag[i]==1 && i <= r-1){ // G to G, right single strand part, rev-comp
                                l2 += log(1-delta_s);
                            }else if(reffrag[i]==2 && i <= r-1){ // C to C, right single strand part, rev-comp
                                l2 += 0;
                            }else if (reffrag[i]==1 && i <= L-npos-1){ // G to G, right double strand part, right of the nick, rev-comp
                                l2 += log(1-delta);
                            }else if (reffrag[i]==2 && i <= L-npos-1){ // C to C, right double strand part, right of the nick, rev-comp
                                l2 += 0;
                            }else if (reffrag[i]==1 && i > L-npos-1){ // G to G, right double strand part, left of the nick, rev-comp
                                l2 += 0;
                            }else if (reffrag[i]==2 && i > L-npos-1){ // C to C, right double strand part, left of the nick, rev-comp
                                l2 += log(1-delta);
                            }else if (reffrag[i]==frag[i]){ // All other possible no changes, rev-comp
                                l2 += 0;
                            }else if (reffrag[i]!=frag[i]){ // All other possible changes, rev-comp
                                l2 += 0;
                            }
                        }
                        
                        // double error1 = PhredError[seqError[L-i-1]];
                        // double error1AThird = PhredErrorAThird[seqError[L-i-1]];
                        if (reffrag[L-i-1]<4){
                            if (reffrag[L-i-1]==2 && i<=l-1){ // C to C, left single strand part, rev-comp
                                l2 += log(1-delta_s);
                            }else if (reffrag[L-i-1]==1 && i<=l-1){ // G to G, left single strand part, rev-comp
                                l2 += 0;
                            }else if (reffrag[L-i-1]==2){ // C to C, left double strand part, rev-comp
                                l2 += log(1-delta);
                            }else if (reffrag[L-i-1]==1){ // G to G, left double strand part, rev-comp
                                l2 += 0;
                            }else if (reffrag[L-i-1]==frag[L-i-1]){ // All other possible no changes, rev-comp
                                l2 += 0;
                            }else if (reffrag[L-i-1]!=frag[L-i-1]){ // All other possible changes, rev-comp
                                l2 += 0;
                            }
                        }
                    }
                    double dpnick = npos < L-r-1 ? nv/(1+(L-l-r-2)*nv) : (1-nv)/(1+(L-l-r-2)*nv);
                    np_pmd_anc += exp(l2)*dp*dpnick;
                    pnick += dpnick;
                }
                
                //nick occurring outside both l_check region
                //--------- - - - -                        - - - ------------//
                //                         /|\                               //
                //                          |___ Nick                        //
                double l2 = 0;
                for (int i=0; i<l_check; i++){
                    // double error1 = PhredError[seqError[L-i-1]];
                    // double error1AThird = PhredErrorAThird[seqError[L-i-1]];
                    if (reffrag[L-i-1]<4){
                        if (reffrag[L-i-1]==2 && i<=l-1){ // C to C, left single strand part, rev-comp
                            l2 += log(1-delta_s);
                        }else if (reffrag[L-i-1]==1 && i<=l-1){ // G to G, left single strand part, rev-comp
                            l2 += 0;
                        }else if (reffrag[L-i-1]==2 && i>l-1){ // C to C, left double strand part, rev-comp
                            l2 += log(1-delta);
                        }else if (reffrag[L-i-1]==1 && i>l-1){ // G to G, left double strand part, rev-comp
                            l2 += 0;
                        }else if (reffrag[L-i-1]==frag[L-i-1]){ // All other possible no changes, rev-comp
                            l2 += 0;
                        }else if (reffrag[L-i-1]!=frag[L-i-1]){ // All other possible changes, rev-comp
                            l2 += 0;
                        }
                    }
                    
                    // double error2 = PhredError[seqError[i]];
                    // double error2AThird = PhredErrorAThird[seqError[i]];
                    if (reffrag[i]<4){
                        if (reffrag[i]==1 && i <= r-1){ // G to G, right single strand part, rev-comp
                            l2 += log(1-delta_s);
                        }else if(reffrag[i]==2 && i <= r-1){ // C to C, right single strand part, rev-comp
                            l2 += 0;
                        }else if (reffrag[i]==1 && i > r-1){ // G to G, right double strand part, rev-comp
                            l2 += log(1-delta);
                        }else if (reffrag[i]==2 && i > r-1){ // C to C, right double strand part, rev-comp
                            l2 += 0;
                        }else if (reffrag[i]==frag[i]){ // All other possible no changes, rev-comp, rev-comp
                            l2 += 0;
                        }else if (reffrag[i]!=frag[i]){ // All other possible changes, rev-comp
                            l2 += 0;
                        }
                    }
                }
                np_pmd_anc += exp(l2)*dp*(1-pnick);
            }
            p += dp;
        }
        
        for (int r=l_check+1; r<=L-2-l; r++){ // To guarantee l+r < = L-2
            double pr = 0.5*lambda*pow(1-lambda,r)/(1-pow(1-lambda,L-1));
            if (r == 0){
                pr += 0.5;
            }
            double dp = 4*pl*pr/(3+(1-pow(1-lambda,L-1)-(L-1)*lambda*pow(1-lambda,L-1))/(pow(1-pow(1-lambda,L-1),2)));
            if (dp >= Tol){
                double pnick = 0;
                //--------- - - - -                        ------------------//
                //            /|\                                            //
                //             |___ Nick                                     //
                for (int npos=l; npos<=(l_check-1 < L-r-1 ? l_check-1 : L-r-1); npos++){ //Nick is within left l_check region
                    double l2 = 0;
                    for (int i=0; i<l_check; i++){
                        // double error1 = PhredError[seqError[L-i-1]];
                        // double error1AThird = PhredErrorAThird[seqError[L-i-1]];
                        if (reffrag[L-i-1]<4){
                            if (reffrag[L-i-1]==2 && i<=l-1){ // C to C, left single strand part, rev-comp
                                l2 += log(1-delta_s);
                            }else if (reffrag[L-i-1]==1 && i<=l-1){ // G to G, left single strand part, rev-comp
                                l2 += 0;
                            }else if (reffrag[L-i-1]==2 && i<=npos){ // C to C, left double strand part before nick, rev-comp
                                l2 += log(1-delta);
                            }else if (reffrag[L-i-1]==1 && i<=npos){ // G to G, left double strand part before nick, rev-comp
                                l2 += 0;
                            }else if (reffrag[L-i-1]==2 && i>npos){ // C to C, left double strand part after nick, rev-comp
                                l2 += 0;
                            }else if (reffrag[L-i-1]==1 && i>npos){ // G to G, left double strand part after nick, rev-comp
                                l2 += log(1-delta);
                            }else if (reffrag[L-i-1]==frag[L-i-1]){ // All other possible no changes, rev-comp
                                l2 += 0;
                            }else if (reffrag[L-i-1]!=frag[L-i-1]){ // All other possible changes, rev-comp
                                l2 += 0;
                            }
                        }
                        
                        // double error2 = PhredError[seqError[i]];
                        // double error2AThird = PhredErrorAThird[seqError[i]];
                        if (reffrag[i]<4){
                            if (reffrag[i]==1){ // G to G, right single strand part, rev-comp
                                l2 += log(1-delta_s);
                            }else if(reffrag[i]==2){ // C to C, right single strand part, rev-comp
                                l2 += 0;
                            }else if (reffrag[i]==frag[i]){ // All other possible no changes, rev-comp
                                l2 += 0;
                            }else if (reffrag[i]!=frag[i]){ // All other changes, rev-comp
                                l2 += 0;
                            }
                        }
                    }
                    double dpnick = nv/(1+(L-l-r-2)*nv);
                    np_pmd_anc += exp(l2)*dp*dpnick;
                    pnick += dpnick;
                }
                //nick occurring outside both l_check region
                //--------- - - - -                        ------------------//
                //                           /|\                             //
                //                            |___ Nick                      //
                double l2 = 0;
                for (int i=0; i<l_check; i++){
                    // double error1 = PhredError[seqError[L-i-1]];
                    // double error1AThird = PhredErrorAThird[seqError[L-i-1]];
                    if (reffrag[L-i-1]<4){
                        if (reffrag[L-i-1]==2 && i<=l-1){ // C to C, left single strand part, rev-comp
                            l2 += log(1-delta_s);
                        }else if (reffrag[L-i-1]==1 && i<=l-1){ // G to G, left single strand part, rev-comp
                            l2 += 0;
                        }else if (reffrag[L-i-1]==2 && i>l-1){ // C to C, left double strand part, rev-comp
                            l2 += log(1-delta);
                        }else if (reffrag[L-i-1]==1 && i>l-1){ // G to G, left double strand part, rev-comp
                            l2 += 0;
                        }else if (reffrag[L-i-1]==frag[L-i-1]){ // All other possible no changes, rev-comp
                            l2 += 0;
                        }else if (reffrag[L-i-1]!=frag[L-i-1]){ // All other possible changes, rev-comp
                            l2 += 0;
                        }
                    }
                    
                    // double error2 = PhredError[seqError[i]];
                    // double error2AThird = PhredErrorAThird[seqError[i]];
                    if (reffrag[i]<4){
                        if (reffrag[i]==1){ // G to G, right single strand part, rev-comp
                            l2 += log(1-delta_s);
                        }else if(reffrag[i]==2){ // C to C, right single strand part, rev-comp
                            l2 += 0;
                        }else if (reffrag[i]==frag[i]){ // All other possible no changes, rev-comp
                            l2 += 0;
                        }else if (reffrag[i]!=frag[i]){ // All other possible changes, rev-comp
                            l2 += 0;
                        }
                    }
                }
                np_pmd_anc += exp(l2)*dp*(1-pnick);
            }
            p += dp;
        }
    }
    
    //-----------------                        - - - - ----------//
    //                                            /|\            //
    //                                             |___ Nick     //
    for (int r=0; r<=l_check; r++){
        double pr = 0.5*lambda*pow(1-lambda,r)/(1-pow(1-lambda,L-1));
        if (r == 0){
            pr += 0.5;
        }
        for (int l=l_check+1;l<=L-2-r; l++){
            double pl = 0.5*lambda*pow(1-lambda,l)/(1-pow(1-lambda,L-1));
            if (l==0){
                pl += 0.5;
            }
            // Joint probability of (l,r)
            double dp = 4*pl*pr/(3+(1-pow(1-lambda,L-1)-(L-1)*lambda*pow(1-lambda,L-1))/(pow(1-pow(1-lambda,L-1),2)));
            if (dp > Tol){
                double pnick = 0;
                for (int npos = (L-l_check-1 > l ? L-l_check-1 : l); npos <= L-r-1; npos++){ //Nick is within right l_check region
                    double l2 = 0;
                    for (int i=0; i<l_check; i++){
                        // double error2 = PhredError[seqError[i]];
                        // double error2AThird = PhredErrorAThird[seqError[i]];
                        if (reffrag[i]<4){
                            if (reffrag[i]==1 && i <= r-1){ // G to G, right single strand part, rev-comp
                                l2 += log(1-delta_s);
                            }else if(reffrag[i]==2 && i <= r-1){ // C to C, right single strand part, rev-comp
                                l2 += 0;
                            }else if (reffrag[i]==1 && i <= L-npos-1){ // G to G, right double strand part, right of the nick, rev-comp
                                l2 += log(1-delta);
                            }else if (reffrag[i]==2 && i <= L-npos-1){ // C to C, right double strand part, right of the nick, rev-comp
                                l2 += 0;
                            }else if (reffrag[i]==1 && i > L-npos-1){ // G to A, right double strand part, left of the nick, rev-comp
                                l2 += 0;
                            }else if (reffrag[i]==2 && i > L-npos-1){ // C to C, right double strand part, left of the nick, rev-comp
                                l2 += log(1-delta);
                            }else if (reffrag[i]==frag[i]){ // All other possible no changes, rev-comp
                                l2 += 0;
                            }else if (reffrag[i]!=frag[i]){ // All other possible changes, rev-comp
                                l2 += 0;
                            }
                        }
                        
                        // double error1 = PhredError[seqError[L-i-1]];
                        // double error1AThird = PhredErrorAThird[seqError[L-i-1]];
                        if (reffrag[L-i-1]<4){
                            if (reffrag[L-i-1]==2){ // C to C, left single strand part, rev-comp
                                l2 += log(1-delta_s);
                            }else if (reffrag[L-i-1]==1){ // G to G, left single strand part, rev-comp
                                l2 += 0;
                            }else if (reffrag[L-i-1]==frag[L-i-1]){ // All other possible no changes, rev-comp
                                l2 += 0;
                            }else if (reffrag[L-i-1]!=frag[L-i-1]){ // All other possible changes, rev-comp
                                l2 += 0;
                            }
                        }
                    }
                    double dpnick = npos < L-r-1 ? nv/(1+(L-l-r-2)*nv) : (1-nv)/(1+(L-l-r-2)*nv);
                    np_pmd_anc += exp(l2)*dp*dpnick;
                    pnick += dpnick;
                }
                
                //-----------------                        - - - - ----------//
                //                             /|\                           //
                //                              |___ Nick                    //
                // nick occurring outside both l_check region
                double l2 = 0;
                for (int i=0; i<l_check; i++){
                    // double error1 = PhredError[seqError[L-i-1]];
                    // double error1AThird = PhredErrorAThird[seqError[L-i-1]];
                    if (reffrag[L-i-1]<4){
                        if (reffrag[L-i-1]==2){ // C to C, left single strand part, rev-comp
                            l2 += log(1-delta_s);
                        }else if (reffrag[L-i-1]==1){ // G to G, left single strand part, rev-comp
                            l2 += 0;
                        }else if (reffrag[L-i-1]==frag[L-i-1]){ // All other possible no changes, rev-comp
                            l2 += 0;
                        }else if (reffrag[L-i-1]!=frag[L-i-1]){ // All other possible changes, rev-comp
                            l2 += 0;
                        }
                    }
                    
                    // double error2 = PhredError[seqError[i]];
                    // double error2AThird = PhredErrorAThird[seqError[i]];
                    if (reffrag[i]<4){
                        if (reffrag[i]==1 && i <= r-1){ // G to G, right single strand part, rev-comp
                            l2 += log(1-delta_s);
                        }else if(reffrag[i]==2 && frag[i]==0 && i <= r-1){ // C to T change, right single strand part, rev-comp
                            l2 += 0;
                        }else if (reffrag[i]==1 && i > r-1){ // G to A change, right double strand part, rev-comp
                            l2 += log(1-delta);
                        }else if (reffrag[i]==2 && i > r-1){ // C to T change, right double strand part, rev-comp
                            l2 += 0;
                        }else if (reffrag[i]==frag[i]){ // All other possible no changes, rev-comp
                            l2 += 0;
                        }else if (reffrag[i]!=frag[i]){ // All other possible changes, rev-comp
                            l2 += 0;
                        }
                    }
                }
                np_pmd_anc += exp(l2)*dp*(1-pnick);
            }
            p += dp;
        }
    }
    
    //-----------------                        ------------------//
    //                             /|\                           //
    //                              |___ Nick                    //
    //Investigate each possible left and right 5' overhang pair with (l,r) > l_check (15)
    if (1-p > Tol){
        double l2 = 0;
        for (int i=0; i<l_check; i++){
            // double error1 = PhredError[seqError[L-i-1]];
            // double error1AThird = PhredErrorAThird[seqError[L-i-1]];
            if (reffrag[L-i-1]<4){
                if (reffrag[L-i-1]==2){ // C to C, left single strand part, rev-comp
                    l2 += log(1-delta_s);
                }else if (reffrag[L-i-1]==1){ // G to G, left single strand part, rev-comp
                    l2 += 0;
                }else if (reffrag[L-i-1]==frag[L-i-1]){ // All other possible no changes, rev-comp
                    l2 += 0;
                }else if (reffrag[L-i-1]!=frag[L-i-1]){ // All other possible changes, rev-comp
                    l2 += 0;
                }
            }
            
            // double error2 = PhredError[seqError[i]];
            // double error2AThird = PhredErrorAThird[seqError[i]];
            if (reffrag[i]<4){
                if (reffrag[i]==1){ // G to G, right single strand part, rev-comp
                    l2 += log(1-delta_s);
                }else if(reffrag[i]==2){ // C to C, right single strand part, rev-comp
                    l2 += 0;
                }else if (reffrag[i]==frag[i]){ // All other possible no changes, rev-comp
                    l2 += 0;
                }else if (reffrag[i]!=frag[i]){ // All other possible changes, rev-comp
                    l2 += 0;
                }
            }
        }
        np_pmd_anc += exp(l2)*(1-p);
    }
    return np_pmd_anc;
}

// The function below is for calculating the posterior prob of being ancient
double AncProb(char reffrag[], char frag[], int L, double lambda, double delta, double delta_s, double nv, int seqError[], char* model, double eps, double anc_mu, double anc_si, double mod_mu, double mod_si, int isrecal, int len_limit, int len_min,double Tol){
    double l_err = ErrorLik(reffrag, frag,  L,  seqError); // Modern Likelihood/Only Sequencing-error Likelihood
    double prior_anc = 1-eps;
    double post_anc;
    double l_anc_l = 1;
    double l_err_l = 1;
    if (isrecal == 1){
        double x_max1 = ((double)len_limit-1+0.5-anc_mu)/anc_si;
        double x_min1 = ((double)len_min-0.5-anc_mu)/anc_si;
        double x_max2 = ((double)len_limit-1+0.5-mod_mu)/mod_si;
        double x_min2 = ((double)len_min-0.5-mod_mu)/mod_si;
        double y_max1 = (std::min((double)len_limit-1,std::max((double)len_min,(double)L))+0.5-anc_mu)/anc_si;
        double y_min1 = (std::max((double)len_min,std::min((double)len_limit-1,(double)L))-0.5-anc_mu)/anc_si;
        double y_max2 = (std::min((double)len_limit-1,std::max((double)len_min,(double)L))+0.5-mod_mu)/mod_si;
        double y_min2 = (std::max((double)len_min,std::min((double)len_limit-1,(double)L))-0.5-mod_mu)/mod_si;
        l_anc_l = NormalINC(y_max1, y_min1, x_max1, x_min1);
        l_err_l = NormalINC(y_max2, y_min2, x_max2, x_min2);
    }
    if (!strcasecmp("b",model)){
      double l_anc_b = PMDLik_b(reffrag, frag, L, lambda, delta, delta_s, nv, seqError,Tol); // Ancient Likelihood based on biotin model
        post_anc = prior_anc * l_anc_b * l_anc_l/(prior_anc * l_anc_b * l_anc_l+ (1-prior_anc) * l_err * l_err_l);
	//double test = exp(-10);
        //cout<<prior_anc<<" "<<l_anc_b<<" "<<l_anc_l<<" "<<prior_anc * l_anc_b * l_anc_l<<" "<<test<<"Lei Lei Lei\n";
    }else if(!strcasecmp("nb",model)){
      double l_anc_nb = 0.5*PMDLik_b(reffrag, frag, L, lambda, delta, delta_s, nv, seqError,Tol)+0.5*PMDLik_nb(reffrag, frag, L, lambda, delta, delta_s, nv, seqError,Tol); // Ancient Likelihood based on non-biotin model
        //        cout << "l_anc_nb is "<<l_anc_nb<<"\n";
        //        cout << "PMDLik_b is "<<PMDLik_b(reffrag, frag, L, lambda, delta, delta_s, nv, seqError)<<"\n";
        //        cout << "PMDLik_nb is "<<PMDLik_nb(reffrag, frag, L, lambda, delta, delta_s, nv, seqError)<<"\n";
        post_anc = prior_anc * l_anc_nb * l_anc_l/(prior_anc * l_anc_nb * l_anc_l + (1-prior_anc) * l_err * l_err_l);
    }else{
        fprintf(stderr,"Please specify a deamination model for further calculations.\n");
        return -1;
    }
    if (post_anc < 0){
        post_anc = 0;
    }else if(post_anc > 1){
        post_anc = 1;
    }
    return post_anc;
}


// The function below is for calculating the posterior prob of being damaged given ancient
double PMDProb(char reffrag[], char frag[], int L, double lambda, double delta, double delta_s, double nv, int seqError[], char* model,double Tol){
    double l_err = ErrorLik(reffrag, frag,  L,  seqError); // Modern Likelihood/Only Sequencing-error Likelihood
    //double prior_pmd = 0.25+0.5*lambda/(1-pow(1-lambda,L-1))*(1-pow((1-lambda)*(1-delta_s/4)/(1-delta/4),L-1))/(1-(1-lambda)*(1-delta_s/4)/(1-delta/4))+0.25*pow(lambda,2)/pow(1-pow(1-lambda,L-1),2)*(1-pow((1-lambda)*(1-delta_s/4)/(1-delta/4),L-1))/pow(1-(1-lambda)*(1-delta_s/4)/(1-delta/4),2)-0.25*pow(lambda,2)/pow(1-pow(1-lambda,L-1),2)*(L-1)*pow((1-lambda)*(1-delta_s/4)/(1-delta/4),L-1)/(1-(1-lambda)*(1-delta_s/4)/(1-delta/4));
    //prior_pmd = 1 - prior_pmd * pow(1-delta,L)/(0.75+0.25*(1-pow(1-lambda,L-1)-(L-1)*lambda*pow(1-lambda,L-1))/pow(1-pow(1-lambda,L-1),2));
    double post_pmd;
    if (!strcasecmp("b",model)){
      double l_anc_b = PMDLik_b(reffrag, frag, L, lambda, delta, delta_s, nv, seqError,Tol);
      double np_pmd_anc_b = NoPMDGivenAnc_b(reffrag, frag, L, lambda, delta, delta_s, nv,Tol);
        post_pmd = 1-np_pmd_anc_b*l_err/l_anc_b;
    }else if(!strcasecmp("nb",model)){
      double l_anc_nb = 0.5*PMDLik_b(reffrag, frag, L, lambda, delta, delta_s, nv, seqError,Tol) + 0.5*PMDLik_nb(reffrag, frag, L, lambda, delta, delta_s, nv, seqError,Tol);
      double np_pmd_anc_nb = 0.5*NoPMDGivenAnc_b(reffrag, frag, L, lambda, delta, delta_s, nv,Tol)+0.5*NoPMDGivenAnc_nb(reffrag, frag, L, lambda, delta, delta_s, nv,Tol);
        post_pmd = 1-np_pmd_anc_nb*l_err/l_anc_nb;
    }else{
        fprintf(stderr,"Please specify a deamination model for further calculations.\n");
        return -1;
    }
    if (post_pmd < 0){
        post_pmd = 0;
    }else if(post_pmd > 1){
        post_pmd = 1;
    }
    return post_pmd;
}

bam_hdr_t* calc_pp_pmd_prob(char *refName,char *fname, char* ofname, char* olik, int mapped_only,int se_only, int mapq, faidx_t *seq_ref, int len_limit, int len_min, char * model, double eps, double lambda, double delta, double delta_s, double nv, double anc_mu, double anc_si, double mod_mu, double mod_si, int isrecal, kstring_t *str_cli,double **deamRateCT,double **deamRateGA,double Tol)
{
  char nuc[6] = "ACGTN";
    fprintf(stderr,"print msg mapped_only: %d\n",mapped_only);
    htsFormat *dingding5 =(htsFormat*) calloc(1,sizeof(htsFormat));
    htsFormat *dingding6 =(htsFormat*) calloc(1,sizeof(htsFormat));

    char reconstructedRef[512];
    char myread[512];
    char myrefe[512];
    char yourread[512];
    char yourrefe[512];
    int  yourqual[512];
    std::pair< kstring_t*, std::vector<int> >  mypair;
    kstring_t *kstr =new kstring_t;
    kstr->l=kstr->m=0;
    kstr->s=NULL;
    mypair.first = kstr;
    samFile *in=NULL;
    samFile *out=NULL;
    
    if (ofname!=NULL){
        out = sam_open_format(ofname, strdup("wb"), dingding6);
        int nthreads = 1;
        htsThreadPool p = {NULL, 0};
        if(nthreads>1){
            if (!(p.pool = hts_tpool_init(nthreads))) {
                fprintf(stderr, "Error creating thread pool\n");
                //return 0;
            }
            //hts_set_opt(in,  HTS_OPT_THREAD_POOL, &p);
            if (out) hts_set_opt(out, HTS_OPT_THREAD_POOL, &p);
        }
        //queue_t *queue1 = init_queue_t(nreads_per_pos);
    }
    //queue_t *queue1 = init_queue_t(nreads_per_pos);
    //queue_t *queue = init_queue_t(nreads_per_pos);
    char refeBase, readBase;
    //code below only relevant if using cram files
    if(refName!=NULL){
        char *ref =(char*) malloc(10 + strlen(refName) + 1);
        snprintf(ref,10 + strlen(refName) + 1, "reference=%s", refName);
        hts_opt_add((hts_opt **)&dingding5->specific,ref);
        free(ref);
    }
    if(strstr(fname,".cram")!=NULL && refName==NULL){
        fprintf(stderr,"\t-> cram file requires reference with -T FILE.fa \n");
        exit(0);
    }
    if((in=sam_open_format(fname,"r",dingding5))==NULL ){
        fprintf(stderr,"[%s] nonexistant file: %s\n",__FUNCTION__,fname);
        exit(0);
    }
    
    bam_hdr_t  *hdr = sam_hdr_read(in);
    sam_hdr_add_pg(hdr,"metadamage_briggs","VN","1.0","CL",str_cli->s, NULL);
    assert(sam_hdr_write(out, hdr) >= 0);
    
    bam1_t *b = bam_init1();
        
    int ret;
    int refId=-1;
    double num = 0.0;
    size_t max_site;
    nproc1 = 0;
    
    uchar * indref = NULL;
    while(((ret=sam_read1(in,hdr,b)))>0) {
      if(mapped_only!=0){
	if(b->core.flag&4)
	  continue;
      }
        
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

      nproc1++;

      if(refId==-1||refId!=b->core.tid){
	refId=b->core.tid;
	fprintf(stderr,"\t-> Now at Chromosome: %s\n",hdr->target_name[refId]);
      }

      //then we simply write it to the output
      memset(reconstructedRef,0,512);
      memset(myread,'N',512);
      memset(myrefe,'N',512);
      
      if (seq_ref != NULL){
	wrapperwithref(b,hdr,myread,myrefe,seq_ref);
      }else{
	reconstructRefWithPosHTS(b,mypair,reconstructedRef);
	wrapper(b,mypair.first->s,mypair.second,0,0,NULL,NULL,1,myread,myrefe);
      }
      
      if (len_limit <=0){
	len_limit = 512;
      };

      if (b->core.l_qseq>=30 && b->core.l_qseq<len_limit){
	for (int cycle=0;cycle<b->core.l_qseq;cycle++){
	  yourqual[cycle] = -1;
	}

	for (int cycle=0;cycle<b->core.l_qseq;cycle++){
	  refeBase = refToChar[myrefe[cycle]];
	  readBase = refToChar[myread[cycle]];
	  int dist5p=cycle;
	  int dist3p=b->core.l_qseq-1-cycle;
	  //cout<<"flag "<<b->core.flag<<"\n";
	  if( bam_is_rev(b) ){
	    // cout<<"rev "<<"\t";
	    refeBase=com[refeBase];
	    readBase=com[readBase];
	    //dist5p=int(al.QueryBases.size())-1-i;
	    dist5p=int(b->core.l_qseq)-1-cycle;
	    dist3p=cycle;
	  }
	  yourread[dist5p] = readBase;
	  yourrefe[dist5p] = refeBase;
	  yourqual[dist5p] = bam_get_qual(b)[cycle];
	}
        
	
	double PostAncProb1 = AncProb(yourrefe, yourread, b->core.l_qseq, lambda, delta, delta_s, nv, yourqual, model ,eps, anc_mu, anc_si, mod_mu, mod_si, 0, len_limit, len_min,Tol);
	double PostAncProb2 = 0;
	
	if (isrecal==1) {
	  PostAncProb2 = AncProb(yourrefe, yourread, b->core.l_qseq, lambda, delta, delta_s, nv, yourqual, model ,eps, anc_mu, anc_si, mod_mu, mod_si, 1, len_limit, len_min,Tol);
	}
	
	double PostPMDProb = PMDProb(yourrefe, yourread, b->core.l_qseq, lambda, delta, delta_s, nv, yourqual, model,Tol);
	//cout << "length 0" << " " << b->core.l_qseq << "\n";

#if 0
	// Read Name
	std::cout << bam_get_qname(b) << "\t";
	cout << bam_aux_get(b,"FL") << "\t";
	for (int cycle=0;cycle<b->core.l_qseq;cycle++)
	  std::cout<<nuc[(int)refToChar[myread[cycle]]];
	std::cout<<"\t";
	cout << "length 2" << b->core.l_qseq << "\n";
	for (int cycle=0;cycle<b->core.l_qseq;cycle++)
	  std::cout<<nuc[(int)refToChar[myrefe[cycle]]];
	std::cout<<"\t";
	
	for (int cycle=0;cycle<b->core.l_qseq;cycle++)
	  std::cout<<(char)(bam_get_qual(b)[cycle]+33);
	if (isrecal==1)
	  std::cout<<"\t"<<"AO:f:"<<PostAncProb1<<"\t"<<"AN:f:"<<PostAncProb2<<"\t"<<"PD:f:"<<PostPMDProb<<"\n";
	else
	  std::cout<<"\t"<<"AN:f:"<<PostAncProb1<<"\t"<<"PD:f:"<<PostPMDProb<<"\n";
#endif
	//Bam output is the bam file orientation
	if (isrecal==1){
	  bam_aux_update_float(b,"AO",PostAncProb1);
	  bam_aux_update_float(b,"AN",PostAncProb2);
	}else{
	  bam_aux_update_float(b,"AN",PostAncProb1);
	}
	bam_aux_update_float(b,"PD",PostPMDProb);
	
	assert(sam_write1(out, hdr, b)>=0);
      }
      num += 1.0; //number of intersected reads.
      //cout<<"Number is "<<num<<"\n";

    }
    
    if (ofname != NULL)
      assert(sam_close(out)==0);
    hts_opt_free((hts_opt *)dingding6->specific);
    free(dingding6);
    
    bam_destroy1(b);
    hts_opt_free((hts_opt *)dingding5->specific);
    free(dingding5);
    
    assert(sam_close(in)==0);
 
    return hdr;
}
