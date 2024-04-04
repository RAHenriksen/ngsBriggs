# ngsBriggs
a fast regression method to estimate postmortem damage


## INSTALLATION & REQUIREMENTS
* Use local installation of htslib

git clone https://github.com/RAHenriksen/ngsBriggs.git

git clone https://github.com/samtools/htslib.git

cd htslib; make; cd ../NGSNGS; make HTSSRC=../htslib

* Use systemwide installation of htslib

git clone https://github.com/RAHenriksen/ngsBriggs.git

cd NGSNGS; make

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
./ngsbriggs -bam /projects/korneliussen/people/wql443/NGSNGS/HgSimTest1.bam -ref /projects/korneliussen/people/wql443/NGSNGS/Test_Examples/Mycobacterium_leprae.fa.gz -model b

~~~~
