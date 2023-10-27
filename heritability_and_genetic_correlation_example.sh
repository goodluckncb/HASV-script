# In this study, we have a total of 24 hippocampal and subfield volumes (HASV) traits, see paper for detail information. 
#Therefore, in the following code, we will use "Left_hippocampus" (one of HASV traits) and "Left_fimbria" (one of HASV traits) as examples for illustration


#To calculate heritability and genetic correlations, we utilized the ldsc software (https://www.nature.com/articles/ng.3211) (https://github.com/bulik/ldsc)
#prepare summary file [ES_IS_summary.txt]
#must have these column 
# .snpid rs number
# .chr chromosome
# .bp physical (base pair) position
# .al  first allele , used as the effect allele
# .a2  second allele in bim file, used as the reference allele
# .beta Effect size (OR or BETA)
# .pval p value

cd path/to/ldsc

#prepare file
path/to/ldsc/munge_sumstats.py --sumstats path/to/gwas_summary_dir/Left_hippocampus_GWAS.txt  --N 41197 --out Left_hippocampus_summary
#This step will generate a file named Left_hippocampus_summary.sumstats.gz

path/to/ldsc/munge_sumstats.py --sumstats path/to/gwas_summary_dir/Left_fimbria_GWAS.txt  --N 41057 --out Left_fimbria_summary
#This step will generate a file named Left_fimbria_summary.sumstats.gz

#calculate heritability
path/to/ldsc/ldsc.py  \
--ref-ld-chr path/to/ldsc/eur_w_ld_chr/ \
--out Left_hippocampus.h2  \
--h2 Left_fimbria_summary.sumstats.gz  \
--w-ld-chr path/to/ldsc/eur_w_ld_chr/ 


#calculate genetic correlations
path/to/ldsc/ldsc.py  \
--ref-ld-chr path/to/ldsc/eur_w_ld_chr/  \
--out Left_hippocampus_Left_fimbria.corr  \
--rg Left_hippocampus_summary.sumstats.gz,Left_fimbria_summary.sumstats.gz \
--w-ld-chr path/to/ldsc/eur_w_ld_chr/ 
