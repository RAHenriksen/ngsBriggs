karensvaneborg@imac ngsBriggs % time ./ngsbriggs -model nb -bdamage Chr22_024_36_68_0097.bdamage.gz -rlens Chr22_024_36_68_0097.rlens.gz 2>&1|grep  lambda|md5sum
ff00321026e0479f4e5121f01ce95232  -
./ngsbriggs -model nb -bdamage Chr22_024_36_68_0097.bdamage.gz -rlens  2>&1  0.54s user 0.00s system 99% cpu 0.545 total
grep lambda  0.00s user 0.00s system 0% cpu 0.544 total
md5sum  0.00s user 0.00s system 0% cpu 0.544 total

karensvaneborg@imac ngsBriggs % time ./ngsbriggs -bam Chr22_024_36_68_0097_eps10.sorted.MD.bam -ref chr22.fa -eps 0.1 -model nb 2>&1|grep lambda|md5sum          
7c70930ec2f1349a52f42973f5b0bcb3  -
./ngsbriggs -bam Chr22_024_36_68_0097_eps10.sorted.MD.bam -ref chr22.fa -eps   0.54s user 0.01s system 99% cpu 0.554 total
grep lambda  0.00s user 0.00s system 0% cpu 0.553 total
md5sum  0.00s user 0.00s system 0% cpu 0.552 total
karensvaneborg@imac ngsBriggs %


karensvaneborg@imac ngsBriggs % time ./ngsbriggs -bam Chr22_024_36_68_0097_eps10.sorted.MD.coord.bam -model nb -eps 0.1 -ibam Chr22_024_36_68_0097_eps10.sorted.MD.coord.bam -obam Chr22_024_36_68_0097_eps10.sorted.MD.scores.r0.bam -isrecal 0 -ibed chr22.bed -chr chr22 -nthread 8
	-> fasta load 
Loading the bamfile
The bam file Chr22_024_36_68_0097_eps10.sorted.MD.coord.bam has been inputted successfully
Reference is not provided, and it will be reconstructed according to the MD tags!
mapped_only: 0
We will focus on Chromosome chr22!
	-> Now at Chromosome: chr22

	[ALL done] cpu-time used =  0.00 sec
	[ALL done] walltime used =  0.00 sec

The misincorporting matrix is as follows:
Dir.	Pos.	FreqCT	FreqGA
5'	1	0.232472	0.000000case compile in
esac

5'	2	0.145299	0.000000
5'	3	0.088353	0.000000
5'	4	0.071111	0.008621
5'	5	0.049808	0.000000
3'	1	0.000000	0.188976
3'	2	0.004274	0.149194
3'	3	0.004274	0.111111
3'	4	0.004274	0.052174
3'	5	0.000000	0.039216

The inference is under progress...
The fragment length distribution is calculated from the bam file!
The chosen model is non-biotin model.

	[ALL done] cpu-time used =  0.5229 sec
	[ALL done] walltime used =  0.0000 sec
The chosen model is non-biotin model, the inferred parameters are given as follows:
lambda: 0.348019 (0.300381,0.395657), delta: 0.008409 (0.004223,0.012595), delta_s: 0.730853 (0.595089,0.866617), nu: 0.004318 (-0.009951,0.018587), nfunctioncalls: 33 ngradientcalls: 33, llh: -4.633483.
print msg mapped_only: 0
Loading the bedfile chr22.bed ...
No meaningful chromosome name is provided, therefore we will focus on all provided chromosomes!
	-> Now at Chromosome: chr22
./ngsbriggs -bam Chr22_024_36_68_0097_eps10.sorted.MD.coord.bam -model nb -ep  8.35s user 0.02s system 99% cpu 8.379 total
karensvaneborg@imac ngsBriggs %

time ./ngsbriggs -bam Chr22_024_36_68_0097_eps10.sorted.MD.coord.bam -model nb -eps 0.1 -ibam Chr22_024_36_68_0097_eps10.sorted.MD.coord.bam -obam Chr22_024_36_68_0097_eps10.sorted.MD.scores.r1.bam -isrecal 1 -ibed chr22.bed -chr chr22 -nthread 8

karensvaneborg@imac ngsBriggs % samtools view Chr22_024_36_68_0097_eps10.sorted.MD.scores.r0.bam|md5sum
a6097a9e953c012b99124efb4371743a  -
karensvaneborg@imac ngsBriggs % samtools view Chr22_024_36_68_0097_eps10.sorted.MD.scores.r1.bam|md5sum
dffc345f6a0d0f6786950270b786b3b1  -
karensvaneborg@imac ngsBriggs % 


karensvaneborg@imac test % for i in *.bam;do samtools view $i|md5sum;done
a6097a9e953c012b99124efb4371743a  -
dffc345f6a0d0f6786950270b786b3b1  -
karensvaneborg@imac test % cd ..
karensvaneborg@imac ngsBriggs %                                               
karensvaneborg@imac ngsBriggs % samtools view Chr22_024_36_68_0097_eps10.sorted.MD.scores.r0.bam|md5sum
a6097a9e953c012b99124efb4371743a  -
karensvaneborg@imac ngsBriggs % samtools view Chr22_024_36_68_0097_eps10.sorted.MD.scores.r1.bam|md5sum
dffc345f6a0d0f6786950270b786b3b1  -
karensvaneborg@imac ngsBriggs % 



##new version
time ./ngsbriggs -bam Chr22_024_36_68_0097_eps10.sorted.MD.coord.bam -model nb -eps 0.1 -ibam Chr22_024_36_68_0097_eps10.sorted.MD.coord.bam -obam Chr22_024_36_68_0097_eps10.sorted.MD.scores.r1.bam -isrecal 1 -nthread 8
#./ngsbriggs -bam Chr22_024_36_68_0097_eps10.sorted.MD.coord.bam -model nb -ep  1652.89s user 4.60s system 754% cpu 3:39.59 total
#samtools view Chr22_024_36_68_0097_eps10.sorted.MD.scores.r1.bam|md5sum
#dffc345f6a0d0f6786950270b786b3b1  -


time ./ngsbriggs -bam Chr22_024_36_68_0097_eps10.sorted.MD.coord.bam -model nb -eps 0.1 -ibam Chr22_024_36_68_0097_eps10.sorted.MD.coord.bam -obam Chr22_024_36_68_0097_eps10.sorted.MD.scores.r0.bam -isrecal 0 -nthread 8
#./ngsbriggs -bam Chr22_024_36_68_0097_eps10.sorted.MD.coord.bam -model nb -ep  8.57s user 0.01s system 99% cpu 8.600 total
#samtools view Chr22_024_36_68_0097_eps10.sorted.MD.scores.r0.bam|md5sum
#a6097a9e953c012b99124efb4371743a  -

time ./ngsbriggs -model nb -bdamage Chr22_024_36_68_0097.bdamage.gz -rlens Chr22_024_36_68_0097.rlens.gz 2>&1|grep  lambda|md5sum
#ff00321026e0479f4e5121f01ce95232  -
#./ngsbriggs -model nb -bdamage Chr22_024_36_68_0097.bdamage.gz -rlens  2>&1  0.57s user 0.01s system 99% cpu 0.575 total
#grep lambda  0.00s user 0.00s system 0% cpu 0.575 total
#md5sum  0.00s user 0.00s system 0% cpu 0.574 total
#karensvaneborg@imac ngsBriggs % 

time ./ngsbriggs -bam Chr22_024_36_68_0097_eps10.sorted.MD.bam -ref chr22.fa -eps 0.1 -model nb 2>&1|grep lambda|md5sum         
#7c70930ec2f1349a52f42973f5b0bcb3  -
#./ngsbriggs -bam Chr22_024_36_68_0097_eps10.sorted.MD.bam -ref chr22.fa -eps   0.56s user 0.02s system 79% cpu 0.734 total


####10may
time ./ngsbriggs -bam Chr22_024_36_68_0097_eps10.sorted.MD.coord.bam -model nb -eps 0.1 -ibam Chr22_024_36_68_0097_eps10.sorted.MD.coord.bam -obam Chr22_024_36_68_0097_eps10.sorted.MD.scores.r0.bam -isrecal 0 -nthread 8


fvr124@SUN1024817 ngsBriggs % md5sum *.bam test/*.bam|grep scores
2f84a18c83fb5b49d70d84ae4eea31bd  Chr22_024_36_68_0097_eps10.sorted.MD.scores.r0.bam
a1e9a1bfd7aea09f4cafb906e778399c  Chr22_024_36_68_0097_eps10.sorted.MD.scores.r1.bam
2f84a18c83fb5b49d70d84ae4eea31bd  test/Chr22_024_36_68_0097_eps10.sorted.MD.scores.r0.bam
a1e9a1bfd7aea09f4cafb906e778399c  test/Chr22_024_36_68_0097_eps10.sorted.MD.scores.r1.bam
fvr124@SUN1024817 ngsBriggs %


fvr124@SUN1024817 ngsBriggs % paste <(samtools view Chr22_024_36_68_0097_eps10.sorted.MD.scores.r1.bam|cut -f15|cut -f3 -d:) <(samtools view test/old/Chr22_024_36_68_0097_eps10.sorted.MD.scores.r1.bam|cut -f15|cut -f3 -d:)|awk '{print $1-$2}'|datamash range 1
0
fvr124@SUN1024817 ngsBriggs % paste <(samtools view Chr22_024_36_68_0097_eps10.sorted.MD.scores.r1.bam|cut -f16|cut -f3 -d:) <(samtools view test/old/Chr22_024_36_68_0097_eps10.sorted.MD.scores.r1.bam|cut -f16|cut -f3 -d:)|awk '{print $1-$2}'|datamash range 1
3.2e-05
fvr124@SUN1024817 ngsBriggs % paste <(samtools view Chr22_024_36_68_0097_eps10.sorted.MD.scores.r1.bam|cut -f17|cut -f3 -d:) <(samtools view test/old/Chr22_024_36_68_0097_eps10.sorted.MD.scores.r1.bam|cut -f17|cut -f3 -d:)|awk '{print $1-$2}'|datamash range 1
1e-11
fvr124@SUN1024817 ngsBriggs % paste <(samtools view Chr22_024_36_68_0097_eps10.sorted.MD.scores.r0.bam|cut -f17|cut -f3 -d:) <(samtools view test/old/Chr22_024_36_68_0097_eps10.sorted.MD.scores.r0.bam|cut -f17|cut -f3 -d:)|awk '{print $1-$2}'|datamash range 1
0
fvr124@SUN1024817 ngsBriggs % paste <(samtools view Chr22_024_36_68_0097_eps10.sorted.MD.scores.r0.bam|cut -f16|cut -f3 -d:) <(samtools view test/old/Chr22_024_36_68_0097_eps10.sorted.MD.scores.r0.bam|cut -f16|cut -f3 -d:)|awk '{print $1-$2}'|datamash range 1
1e-11
fvr124@SUN1024817 ngsBriggs % paste <(samtools view Chr22_024_36_68_0097_eps10.sorted.MD.scores.r0.bam|cut -f15|cut -f3 -d:) <(samtools view test/old/Chr22_024_36_68_0097_eps10.sorted.MD.scores.r0.bam|cut -f15|cut -f3 -d:)|awk '{print $1-$2}'|datamash range 1
0
fvr124@SUN1024817 ngsBriggs % 
