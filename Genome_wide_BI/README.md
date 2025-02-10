## Real Data Application

This folder contains the codes required to perform the real data application in Chapter 5. 

A cvs file containing all the baseline characteristics, specifically including fields "eid", "31-0.0", "34-0.0", "52-0.0", "53-0.0", "53-2.0", "22001-0.0", "22006-0.0", "22009-0.1", "22009-0.2", "22009-0.3", "22009-0.4", "22009-0.5", is needed (in our case, named as ukb50595.csv), together with all the plink data files (.bed, .bim, and .fam) of the genomic data of Chromosome 1-22, and the phenotype data of Field 25753. These can be obtained from UK Biobank via their standard data access procedure. After obtaining the data,  follow the workflow:

1. Run "Initialized_geno_blkwise.R" and "initialize_pheno.R" 

2. For each SNP, compute the ball impurity reduction as demonstrated in "Genome-wide-BI.R". 

3. Run "summarize_final.R" to summarize and plot the manhattan plots. 


For the Ball Tree of Chromosome 19, first filter the .bed, .bim, .fam files through plink command: 
plink --bfile your_input_file --indep-pairwise 500 50 0.2 --out pruned_data

Then compute the ball impurity reduction of the 3740 SNPS individually as  demonstrated in "Genome-wide-BI.R". Iteratively repeat the scanning of each SNP to obtain the tree.  