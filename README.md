# ANGSD-asso-scripts

These are scripts used for the simulations for the article "Efficient approaches for large-scale GWAS with genotype uncertainty"

They produce data for:
Figure 2
Figure 3
Supplementary Figure 1
Supplementary Figure 2
Supplementary Figure 3 and 4
Supplementary Figure 5 and 6
Supplementary Figure 7 and 8
Supplementary Figure 12
Table 2 and 3


It is assumed that the following R-packages are installed (BE AWARE THAT THE CURRENT CODE USES MULTI-THREADING VIA "parallel::mclapply"):
parallel
SQUAREM

Furthermore that the binaries of "angsd" and "snptest_v2.5.4-beta3" are in the directory, to make these scripts work.\
( angsd can be obtained from its github: https://github.com/ANGSD/angsd )\
( snptest can be be obtained here: https://www.well.ox.ac.uk/~gav/snptest/#download )\

Example of how to run:

Rscript simulations.R 1000
Rscript simulTable2and3.R 1000 1000 4 1 1.2 1 1000
Rscript plots.R

The "ANGSD-asso" directory holds R code for running ANGSD-asso's latent model and dosage model in R.



## NOTE ON HOW TO COMPILE ANGSD 

#compile with local version of htslib
make HTSSRC=../htslib/

#Then recently I have been getting an error (/usr/bin/ld: /home/emil/git/htslib/libhts.a(hfile_s3.o): undefined reference to symbol 'EVP_sha256@@OPENSSL_1_1_0')
Can be helped by adding -lcrypto as a LIB
(Also see this link: https://github.com/ANGSD/angsd/pull/151 GO TO "Commits")
