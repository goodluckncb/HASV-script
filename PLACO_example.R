
# In this study, we have a total of 24 hippocampal and subfield volumes (HASV) traits, see paper for detail information. 
# At the same time, we have a total of 10 brian diorders, see paper for detail information. 
#Therefore, in the following code, we will use "Left_hippocampus" (one of HASV traits) and "PD" (Parkinson's disease) as examples for illustration

#https://github.com/RayDebashree/PLACO
require(devtools)
source_url("https://github.com/RayDebashree/PLACO/blob/master/PLACO_v0.1.1.R?raw=TRUE")
library(tidyverse)

#Prepare the required columns ：SNP,PVAL,Z
#x---HASV
Left_hippocampus <- fread("path/to/gwas_summary_dir/Left_hippocampus_GWAS.txt")
Left_hippocampus <- Left_hippocampus[,c("SNP","P_BOLT_LMM_INF","BETA","SE")]
names(Left_hippocampus) <- c("SNP","PVAL","BETA","SE")
Left_hippocampus$Z <- Left_hippocampus$BETA/Left_hippocampus$SE
Left_hippocampus <- Left_hippocampus[,c("SNP","PVAL","Z")]
names(Left_hippocampus) <- c("SNP","P.x","Z,x")

#y----brain disorders
PD <- fread("path/to/gwas_summary_dir/PD_GWAS.txt")
PD <- PD[,c("SNP","PVAL","Z")]
names(PD) <- c("SNP","P.y","Z,y")

data <- inner_join(Left_hippocampus,PD,by='SNP') #inner_join by SNP column

Z.matrix <- data[,c("Z.x","Z.y")]
Z.matrix <- as.matrix.data.frame(Z.matrix)
P.matrix <- data[,c("PVAL.x","PVAL.y")]
P.matrix <- as.matrix.data.frame(P.matrix)


#calculate VarZ
k <- 2
colnames(Z.matrix) <- paste("Z",1:k,sep="")
colnames(P.matrix) <- paste("P",1:k,sep="")
p <- nrow(Z.matrix)
VarZ <- var.placo(Z.matrix, P.matrix, p.threshold=1e-4)


out <- map_dfr(1:p,function(i) placo(Z=Z.matrix[i,], VarZ=VarZ))
result <- cbind(data,out)
sig <- filter(result,p.placo < 5e-8)

#example
# SNP	PVAL.x	Z.x	PVAL.y	Z.y	T.placo	p.placo
# rs55775495	7.5e-05	3.96053644559871	8.376e-07	4.92	19.4858393123457	2.89971467619034e-08
# rs2564389	1.8e-05	4.28667479836026	1.018e-05	4.424	18.9642493079458	4.50786702996999e-08
# rs12621129	1.2e-05	4.37112069491133	4.595e-06	-4.602	-20.1158974379819	1.70267551191391e-08


fwrite(sig,"Left_hippocampus_PD_placo_sig_snp.txt",sep = "\t",col.names = T)
fwrite(result,"Left_hippocampus_PD_placo_whole_snp.txt",sep = "\t",col.names = T)