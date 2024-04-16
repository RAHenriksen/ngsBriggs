# ngsBriggs
a fast regression method to estimate postmortem damage


## INSTALLATION & REQUIREMENTS
* (optional eigenlib)

```
git clone https://gitlab.com/libeigen/eigen
```

Then remember to add
```
FLAGS="-I ../eigen/"
```
to make commands below

* Use local installation of htslib

git clone https://github.com/RAHenriksen/ngsBriggs.git

git clone https://github.com/samtools/htslib.git

cd htslib; make; cd ../ngsBriggs; make HTSSRC=../htslib

* Use systemwide installation of htslib

git clone https://github.com/RAHenriksen/ngsBriggs.git

cd ngsBriggs; make

**NOTE:** Newer version of htslib which includes bam_set1 is required


## GENERAL
ngsBriggs - inferring damage patterns and decontaminating ancient samples - 0.0.1 

~~~~bash
Usage
./ngsBriggs -bam -tab -ref -bed -len -chr -ibam -iref -ibed -ichr -obam -otab -oinf -olen -model -eps -isrecal -olik -nthreads
	-> -bam: The bam file for inference;
	-> -tab: The table file for inference;
	-> -ref: The reference file for inference;
	-> -bed: The bed file (1-based) for inference;
	-> -len: The length mass probability distribution file for inference;
	-> -chr: The focal chromosome for inference;
	-> -ibam: The bam file for ancient strand fishing;
	-> -iref: The reference file for ancient strand fishing;
	-> -ibed: The bed file (1-based) for ancient strand fishing;
	-> -ichr: The focal chromosome for ancient strand fishing;
	-> -obam: The output bam file name;
	-> -otab: The output table file name;
	-> -oinf: The output inferred parameters file name;
	-> -olen: The output length mass probability distribution file name;
	-> -model: Specifying the model, either b (the biotin model) or nb (the non-biotin model);
	-> -eps: The overall modern contamination rate, the value should be within the interval [0,1);
	-> -isrecal: Choose 1 if recalibration based on length is needed, otherwise 0 (default);
	-> -olik: The output nucleotide likelihood file;
	-> -nthreads: Choose the number of threads to speed up the recalibration process.

Example

Inference of Briggs parameters with three different ways
./ngsbriggs -bam Chr22_024_36_68_0097.sorted.MD.bam -ref chr22.fa -model nb
./ngsbriggs -model nb -bdamage Chr22_024_36_68_0097.bdamage.gz -rlens Chr22_024_36_68_0097.rlens.gz
./ngsbriggs -tab mismatch2.txt -len len.txt -model nb

Inference of Briggs parameters with epsilon
./ngsbriggs -bam Chr22_024_36_68_0097_eps10.sorted.MD.bam -ref chr22.fa -eps 0.1 -model nb

Decontamination with epsilon
./ngsbriggs -bam Chr22_024_36_68_0097_eps10.sorted.MD.coord.bam -model nb -eps 0.1 -ibam Chr22_024_36_68_0097_eps10.sorted.MD.coord.bam -obam Chr22_024_36_68_0097_eps10.sorted.MD.scores.bam -isrecal 1 -ibed chr22.bed -chr chr22 -nthread 1

~~~~
