# In this study, we have a total of 24 hippocampal and subfield volumes (HASV) traits, see paper for detail information. 
# At the same time, we have a total of 10 brian diorders, see paper for detail information. 
#Therefore, in the following code, we will use "Left_hippocampus" (one of HASV traits) and "PD" (Parkinson's disease) as examples for illustration

#install_github("dongjunchung/GPA")  #https://github.com/dongjunchung/GPA
library(GPA) 
library(tidyverse)

#Prepare the required columns ：SNP and PVAL
#x---HASV
Left_hippocampus <- fread("path/to/gwas_summary_dir/Left_hippocampus_GWAS.txt")
Left_hippocampus <- Left_hippocampus[,c("SNP","P_BOLT_LMM_INF")]
names(Left_hippocampu) <- c("SNP","P.x")

#y----brain disorders
PD <- fread("path/to/gwas_summary_dir/PD_GWAS.txt")
PD <- PD[,c("SNP","P")]
names(PD) <- c("SNP","P.y")

data <- inner_join(Left_hippocampus,PD,by='SNP')#inner_join by SNP column

pmat = data[,c("P.x","P.y")]
fit.GPA.noAnn <- GPA( pmat, NULL )
fit.GPA.pleiotropy.H0 <- GPA( pmat, NULL, pleiotropyH0=T)
res = pTest( fit.GPA.noAnn, fit.GPA.pleiotropy.H0 )
#save(res,file="Left_hippocampus_PD_GPA.rds")


#example
# res
# $pi
#         00         10         01         11 
# 0.76257196 0.04855922 0.16153690 0.02733192 

# $piSE
#                   pi_10       pi_01       pi_11 
# 0.017039279 0.008967324 0.016848570 0.008654316 

# $statistics
# iteration_2000 
#       3.078879 

# $pvalue
# iteration_2000 
#     0.07931519 

PM11 <- res$pi[4]
PAR <-  res$pi[4]/(res$pi[2]+res$pi[3]+res$pi[4])
p_gpa <- res$pvalue

