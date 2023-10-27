
# In this study, we have a total of 24 hippocampal and subfield volumes (HASV) traits, see paper for detail information. 
# At the same time, we have a total of 10 brian diorders, see paper for detail information. 
#Therefore, in the following code, we will use "Left_hippocampus" (one of HASV traits) and "PD" (Parkinson's disease) as examples for illustration

library("mr.raps")
library(TwoSampleMR)
library(MRPRESSO)
library(MendelianRandomization)
library(penalized)
library(data.table)
library(tidyverse)
library(openxlsx)

#exposure (Left_hippocampus) --Read and process instrument variables
instrument <- fread("Left_hippocampus_IV.txt",header = F) %>% rename(SNP=V1)
exp_dat <- fread("path/to/gwas_summary_dir/Left_hippocampus_GWAS ")
exp_dat <- left_join(instrument,exp_dat,by=c("SNP"))

E <- format_data(
  exp_dat,
  type='exposure',
  snp_col = "SNP",
  beta_col = "BETA",
  se_col = "SE",
  effect_allele_col ="ALLELE1",
  other_allele_col = "ALLELE0",
  pval_col = "P_BOLT_LMM_INF",
  eaf_col = "A1FREQ",
)

#outcome----PD
disease <- fread("path/to/gwas_summary_dir/PD_GWAS.txt")
O <- left_join(instrument ,disease)

SNP_OUTCOME <- format_data(
  O,
  type='outcome',
  snp_col = "SNP",
  beta_col = "BETA",
  se_col = "SE",
  effect_allele_col ="A1",
  other_allele_col = "A2",
  pval_col = "PVAL",
  #eaf_col = "A1",
  samplesize_col = "N"
)



#
dat <- harmonise_data(
  exposure_dat = SNP_EXPOSE, 
  outcome_dat = SNP_OUTCOME)

dat$exposure <- "Left_hippocampus"

dat$outcome <- "PD"

dat <- subset(dat,mr_keep == "TRUE")


#Remove pleiotropic SNPS
presso_loop <- function(dat=dat,NbDistribution=1000){
  if(nrow(dat)<4){
    res_presso=list(1)
    message(paste0("dat < 4 SNP"))
    return(list(dat,res_presso))
  }else{
    set.seed(888) 
    res_presso <- TwoSampleMR::run_mr_presso(dat = dat, NbDistribution = NbDistribution)
    p1 <- res_presso[[1]][["MR-PRESSO results"]][["Global Test"]][["Pvalue"]]
    p2 <- gsub(p1,pattern = "<",replacement = "") %>% as.numeric();p2
  }
  if(p2 >=0.05){
    return(list(dat,res_presso))
  }
  if(nrow(dat)==4){
    message(paste0("dat only have 4 SNPs"))
    return(list(dat,res_presso))
  }
  if(p2 < 0.05){
    message(paste0("第一次的MRPRESSO的P<0.05"))
    if(is.null(res_presso[[1]][["MR-PRESSO results"]][["Outlier Test"]])){
      return(list(dat,res_presso))
    }
    outliers <- data.frame(dat$SNP,res_presso[[1]][["MR-PRESSO results"]][["Outlier Test"]])
    names(outliers)[1] <- "SNP"
    outliers <- dplyr::arrange(outliers,Pvalue)
    noutliers <- nrow(outliers)
    for (ii in c(1:noutliers)) {
      message(paste0("_",ii,"_"))
      dat_new <- subset(dat,!(dat$SNP %in% outliers$SNP[1]))
      set.seed(888)
      res_presso <- TwoSampleMR::run_mr_presso(dat = dat_new, NbDistribution = NbDistribution)
      p1 <- res_presso[[1]][["MR-PRESSO results"]][["Global Test"]][["Pvalue"]]
      p <- res_presso[[1]][[2]][[1]][["Pvalue"]];p
      p2 <- gsub(p,pattern = "<",replacement = "") %>% as.numeric();p2
      if(p2 >= 0.05){
        message(paste0("_",ii,"_","P>=0.05"))
        return(list(dat_new,res_presso))
      }
      if(nrow(dat_new)<=4){
        res_presso=list(1)
        message(paste0("_",ii,"_","dat_new is less than 4 SNPS"))
        return(list(dat,res_presso))
      }
      dat <- dat_new
      if(is.null(res_presso[[1]][["MR-PRESSO results"]][["Outlier Test"]])){
        return(list(dat,res_presso))
      }
        outliers <- data.frame(dat$SNP,res_presso[[1]][["MR-PRESSO results"]][["Outlier Test"]])
        names(outliers)[1] <- "SNP"
        outliers <- dplyr::arrange(outliers,Pvalue)
        noutliers <- nrow(outliers)
    }
  }
}

res <- presso_loop(dat,2000)  


dat <- res[[1]]

#pleiotropy test 
mr_pleiotropy_test <- mr_pleiotropy_test(dat)

x_y <- subset(dat,mr_keep == "TRUE")
###Calculate PVE and F values
dataz <- x_y[,c("SNP","beta.exposure","se.exposure","samplesize.exposure")]
names(dataz) <- c("rsid","beta","se","N")
Ntotal <- max(dataz$N) 
k <- nrow(dataz) 
with(dataz,{
  dataz$PVE <<- 2*(dataz$beta)^2/(2*(dataz$beta)^2+2*dataz$N*(dataz$se)^2)
  dataz$F <<- dataz$PVE*(dataz$N-1-1)/(1-dataz$PVE)
})
Ntotal <- max(dataz$N) 
k <- nrow(dataz) 
R2 <- sum(dataz$PVE) 
overallF <- R2*(Ntotal-k-1)/((1-R2)*k) 
print(c(R2,overallF,length(x_y$rsid)))

mr_frame <- x_y[,c("SNP","beta.exposure","se.exposure","beta.outcome","se.outcome")]
names(mr_frame) <- c("rsid","beta.x","se.x","beta.y","se.y")

mr_object <- mr_input(bx = mr_frame$beta.x, bxse = mr_frame$se.x, 
                      by = mr_frame$beta.y, byse = mr_frame$se.y,
                      snps = mr_frame$rsid)

#######perform MR analysis####

### MR-ivw
mr_ivw <- mr_ivw(mr_object)

### MR-Egger
MR_Egger <- mr_egger(mr_object)

### MR-Robust
MR_Robust <- mr_ivw(mr_object,"random", robust = TRUE)

### MR-median
mr_median <- mr_median(mr_object, weighting = "weighted")

### MR-mbe
mr_mbe <- mr_mbe(mr_object, weighting = "weighted")

### MR-RAPS
MR_RAPS <- mr.raps.overdispersed.robust(mr_frame$beta.x, mr_frame$beta.y, mr_frame$se.x, mr_frame$se.y,
                                        loss.function = "huber", k = 1.345,initialization = c("l2"), 
                                        suppress.warning = FALSE, diagnosis = FALSE, niter = 20,
                                        tol = .Machine$double.eps^0.5)

#perform MR-PRESSO

mr_presso <- mr_presso(BetaOutcome = "beta.y", BetaExposure = "beta.x",
                       SdOutcome = "se.y", SdExposure = "se.x", 
                       OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = mr_frame, 
                       NbDistribution = 3000, SignifThreshold = 0.05)

#Result sorting and saving----------
result <- data.frame(matrix(NA,11,9))
colnames(result) <- c("method","beta","se","CILower","CIUpper","OR","ORLower","ORUpper","pval")
n=4 #The result retains four decimals

#IVW
result[1,1] <- "Inverse-variance weighted"
result[1,2] <- round(mr_ivw@Estimate,n)
result[1,3] <- round(mr_ivw@StdError,n)
result[1,4] <- round(mr_ivw@CILower,n)
result[1,5] <- round(mr_ivw@CIUpper,n)
result[1,6] <- round(exp(mr_ivw@Estimate),n)
result[1,7] <- round(exp(mr_ivw@CILower),n)
result[1,8] <- round(exp(mr_ivw@CIUpper),n)
result[1,9] <- mr_ivw@Pvalue

#MR-Egger
result[2,1] <- "MR-Egger"
result[2,2] <- round(MR_Egger@Estimate,n)
result[2,3] <- round(MR_Egger@StdError.Est,n)
result[2,4] <- round(MR_Egger@CILower.Est,n)
result[2,5] <- round(MR_Egger@CIUpper.Est,n)
result[2,6] <- round(exp(MR_Egger@Estimate),n)
result[2,7] <- round(exp(MR_Egger@CILower.Est),n)
result[2,8] <- round(exp(MR_Egger@CIUpper.Est),n)
result[2,9] <- MR_Egger@Pvalue.Est

result[3,1] <- "intercept"
result[3,2] <- round(MR_Egger@Intercept,n)
result[3,3] <- round(MR_Egger@StdError.Int,n)
result[3,4] <- round(MR_Egger@CILower.Int,n)
result[3,5] <- round(MR_Egger@CIUpper.Int,n)
result[3,6] <- round(exp(MR_Egger@Intercept),n)
result[3,7] <- round(exp(MR_Egger@CILower.Int),n)
result[3,8] <- round(exp(MR_Egger@CIUpper.Int),n)
result[3,9] <- MR_Egger@Pvalue.Int


#MR-Robust
result[4,1] <- "MR-Robust"
result[4,2] <- round(MR_Robust@Estimate,n)
result[4,3] <- round(MR_Robust@StdError,n)
result[4,4] <- round(MR_Robust@CILower,n)
result[4,5] <- round(MR_Robust@CIUpper,n)
result[4,6] <- round(exp(MR_Robust@Estimate),n)
result[4,7] <- round(exp(MR_Robust@CILower),n)
result[4,8] <- round(exp(MR_Robust@CIUpper),n)
result[4,9] <- MR_Robust@Pvalue

#Weighted median
result[5,1] <- "Weighted median"
result[5,2] <- round(mr_median@Estimate,n)
result[5,3] <- round(mr_median@StdError,n)
result[5,4] <- round(mr_median@CILower,n)
result[5,5] <- round(mr_median@CIUpper,n)
result[5,6] <- round(exp(mr_median@Estimate),n)
result[5,7] <- round(exp(mr_median@CILower),n)
result[5,8] <- round(exp(mr_median@CIUpper),n)
result[5,9] <- mr_median@Pvalue

#Weighted mode-based
result[6,1] <- "Weighted mode-based"
result[6,2] <- round(mr_mbe@Estimate,n)
result[6,3] <- round(mr_mbe@StdError,n)
result[6,4] <- round(mr_mbe@CILower,n)
result[6,5] <- round(mr_mbe@CIUpper,n)
result[6,6] <- round(exp(mr_mbe@Estimate),n)
result[6,7] <- round(exp(mr_mbe@CILower),n)
result[6,8] <- round(exp(mr_mbe@CIUpper),n)
result[6,9] <- mr_mbe@Pvalue

#MR-RAPS
result[7,1] <- "MR-RAPS"
result[7,2] <- round(MR_RAPS$beta.hat,n)
result[7,3] <- round(MR_RAPS$beta.se,n)
result[7,4] <- round(MR_RAPS$beta.hat-1.96*MR_RAPS$beta.se,n)
result[7,5] <- round(MR_RAPS$beta.hat+1.96*MR_RAPS$beta.se,n)
result[7,6] <- round(exp(MR_RAPS$beta.hat),n)
result[7,7] <- round(exp(MR_RAPS$beta.hat-1.96*MR_RAPS$beta.se),n)
result[7,8] <- round(exp(MR_RAPS$beta.hat+1.96*MR_RAPS$beta.se),n)
result[7,9] <- MR_RAPS$beta.p.value


#MRPRESSO
result[8,1] <- "MRPRESSO"
result[8,2] <- round(mr_presso$`Main MR results`$`Causal Estimate`[1],n)
result[8,3] <- round(mr_presso$`Main MR results`$Sd[1],n)
result[8,4] <- round(mr_presso$`Main MR results`$`Causal Estimate`[1]-1.96*mr_presso$`Main MR results`$Sd[1],n)
result[8,5] <- round(mr_presso$`Main MR results`$`Causal Estimate`[1]+1.96*mr_presso$`Main MR results`$Sd[1],n)
result[8,6] <- round(exp(mr_presso$`Main MR results`$`Causal Estimate`[1]),n)
result[8,7] <- round(exp(mr_presso$`Main MR results`$`Causal Estimate`[1]-1.96*mr_presso$`Main MR results`$Sd[1]),n)
result[8,8] <- round(exp(mr_presso$`Main MR results`$`Causal Estimate`[1]+1.96*mr_presso$`Main MR results`$Sd[1]),n)
result[8,9] <- mr_presso$`Main MR results`$`P-value`[1]

#MRPRESSO
result[9,1] <- "MRPRESSO_outlier-corrected"
result[9,2] <- round(mr_presso$`Main MR results`$`Causal Estimate`[2],n)
result[9,3] <- round(mr_presso$`Main MR results`$Sd[2],n)
result[9,4] <- round(mr_presso$`Main MR results`$`Causal Estimate`[2]-1.96*mr_presso$`Main MR results`$Sd[2],n)
result[9,5] <- round(mr_presso$`Main MR results`$`Causal Estimate`[2]+1.96*mr_presso$`Main MR results`$Sd[2],n)
result[9,6] <- round(exp(mr_presso$`Main MR results`$`Causal Estimate`[2]),n)
result[9,7] <- round(exp(mr_presso$`Main MR results`$`Causal Estimate`[2]-1.96*mr_presso$`Main MR results`$Sd[2]),n)
result[9,8] <- round(exp(mr_presso$`Main MR results`$`Causal Estimate`[2]+1.96*mr_presso$`Main MR results`$Sd[2]),n)
result[9,9] <- mr_presso$`Main MR results`$`P-value`[2]

##pleiotropy test 

result[10,1] <- "mr_pleiotropy_test"
result[10,2] <- round(mr_pleiotropy_test[1,5],n)
result[10,3] <- round(mr_pleiotropy_test[1,6],n)
result[10,4] <- round(mr_pleiotropy_test[1,5]-1.96*mr_pleiotropy_test[1,6],n)
result[10,5] <- round(mr_pleiotropy_test[1,5]+1.96*mr_pleiotropy_test[1,6],n)
result[10,6] <- round(exp(mr_pleiotropy_test[1,5]),n)
result[10,7] <- round(exp(mr_pleiotropy_test[1,5]-1.96*mr_pleiotropy_test[1,6]),n)
result[10,8] <- round(exp(mr_pleiotropy_test[1,5]+1.96*mr_pleiotropy_test[1,6]),n)
result[10,9] <- mr_pleiotropy_test[1,7]

result[11,1] <- "globle_test"
result[11,9] <- mr_presso[["MR-PRESSO results"]][["Global Test"]][["Pvalue"]]



result$expoure <- as.character(as.character(ivs_name[i,2])) 
result$outcome <- as.character(as.character(ivs_name[i,3]))
result$R2 = R2
result$fstat = overallF
result$num = nrow(x_y)

write.xlsx(result,"Left_hippocampus_PD_MR_result.xlsx"))
