# In this study, we have a total of 24 hippocampal and subfield volumes (HASV) traits, see paper for detail information. 
#Therefore, in the following code, we will use "Left_hippocampus" (one of HASV traits) as examples for illustration


#quality control
plink \
  --bfile path/to/genotype/ukb_genotype \
  --geno 0.05 \
  --hwe 0.000001 \
  --keep Left_hippocampus_ID.txt \
  --maf 0.01 \
  --make-bed \
  --mind 0.1 \
  --out Left_hippocampus \

#modelSnps
plink \
  --bfile Left_hippocampus \
  --indep-pairwise 500 50 0.9 \
  --maf 0.005 \
  --out Left_hippocampus_maf_0.005

#bolt-lmm ( "Perform the GWAS for Left_hippocampus)
#(https://www.nature.com/articles/s41588-018-0144-6)
#The specific meanings of BOLT-related parameters can be found at https://storage.googleapis.com/broad-alkesgroup-public/BOLT-LMM/BOLT-LMM_manual.html
bolt \
    --maxModelSnps=8463542 \
    --bfile=Left_hippocampus \
    --geneticMapFile=path/to/bolt-lmm/BOLT-LMM_v2.3.6/tables/genetic_map_hg19_withX.txt.gz \
    --modelSnps=Left_hippocampus_maf_0.005.prune.in\
    --phenoFile=Left_hippocampus_phenofile.txt \
    --phenoCol=Left_hippocampus_INT \
    --covarFile=Left_hippocampus_phenofile.txt \
    --covarCol=assessment \
    --covarCol=sex \
    --qCovarCol=age \
    --qCovarCol=bmi \
    --qCovarCol=pca{1:10} \
    --lmm \
    --lmmForceNonInf \
    --LDscoresFile=path/to/bolt-lmm/BOLT-LMM_v2.3.6/tables/LDSCORE.1000G_EUR.tab.gz \
    --numThreads=8 \
    --statsFile=path/to/gwas_summary_dir/Left_hippocampus_lmm.stats 