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

