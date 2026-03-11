# -*- coding: utf-8 -*-
# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .r
#       format_name: nomarker
#       format_version: '1.0'
#       jupytext_version: 1.15.2
#   kernelspec:
#     display_name: R
#     language: R
#     name: ir
# ---

library(XMR)
library(MRAPSS)
library(TwoSampleMR)

setwd("./real_data_CSA_AFR")

# ## SHBG to T2D in EAS

# ## XMR

# #### step1: prepare GWAS data (finished，with prepared dara in ./formatted_data; details in format_data.ipynb); copy BBJ_250.2 in previous experiments and rename it as EAS_250.2

# #### step2: estimate bgparas (background parameters)

library(XMAP)
library(data.table)
library(Matrix)

set.seed(1)

# read LD scores
ldscore <- fread("../nc/LDscore_eas_brit_afr.txt")

# pre-process: remove ambiguous SNPs
idx_amb <- which(ldscore$allele1 == comple(ldscore$allele2))
ldscore <- ldscore[-idx_amb,]

ldscore$EUR_AFR = ldscore$AFR_EUR
ldscore$EAS_AFR = ldscore$AFR_EAS
ldscore$EAS_EUR = ldscore$EUR_EAS
ldscore$EUR_EUR = ldscore$EUR
ldscore$EAS_EAS = ldscore$EAS
ldscore$AFR_AFR = ldscore$AFR
ldscore = as.data.frame(ldscore)

ldscore0 = ldscore

bgparas <- estimate_bgparas(
                ldscore0      = ldscore0,
                exposure_code = '30830', 
                outcome_code  = '250.2',      
                pop_ref       = "EUR",         
                pop_target    = "EAS",  
                snp_col       = "rsid",
                gwas_dir      = "./formatted_data" 
            )

bgparas

write.table(bgparas, file = paste0("./EAS_bgparas_", '30830', "_", '250.2'), quote = F, col.names = F, row.names = F)

# #### step3: LD clumping; run XMR

library(XMAP)
library(data.table)
library(Matrix)

ldscore <- fread("../nc/LDscore_eas_brit_afr.txt")

# pre-process: remove ambiguous SNPs
idx_amb <- which(ldscore$allele1 == comple(ldscore$allele2))
ldscore <- ldscore[-idx_amb,]

ldscore$EAS_EUR = ldscore$EUR_EAS
ldscore$EUR_EUR = ldscore$EUR
ldscore$EAS_EAS = ldscore$EAS
ldscore = as.data.frame(ldscore)

ldscore = ldscore[, c("rsid", "EUR", "EAS", "EUR_EAS")]
colnames(ldscore) = c("SNP", "L2.pop1", "L2.pop2", "L12")

ldscore0 = ldscore

base_dir = getwd()

threshold <- 5e-5

t = 1
s = 1

ts1 = '30830'
ts2 = '250.2'

colnames(bgparas) = c('trait1','trait2','pop1','pop2','Omega11','Omega12','Omega22','c1','c12','c2')

bg <- qc_bgparas(
          bgparas       = bgparas,
          exposure_code = ts1[t],
          outcome_code  = ts2[s],
          pop_ref       = "EUR",
          pop_target    = "EAS",
          truncation_corr = 0.95
          )

Omega0 <- bg$Omega0
C      <- bg$C


exp.pop1 = fread(paste0(base_dir, "/formatted_data/EUR_", ts1[t]))
exp.pop2 = fread(paste0(base_dir, "/formatted_data/EAS_", ts1[t]))
out.pop2 = fread(paste0(base_dir, "/formatted_data/EAS_", ts2[s]))

N1 = median(exp.pop1$N)
N2 = median(exp.pop2$N)

paras1 = align_data(dat1 = exp.pop2,
                    dat2 = out.pop2,
                    trait1.name = paste0("EAS_", ts1[t]),
                    trait2.name = paste0("EAS_", ts2[s]),
                    h2.fix.intercept = F,
                    LDSC = F)

paras2 = align_data(dat1 = exp.pop2,
                    dat2 = exp.pop1,
                    trait1.name = paste0("EAS_", ts1[t]),
                    trait2.name = paste0("EUR_", ts1[t]),
                    h2.fix.intercept = F,
                    LDSC = F)

dat1 = paras1$dat
colnames(dat1) = c("SNP", "A1", "A2", "b.exp.pop2", "b.out.pop2", "se.exp.pop2", "se.out.pop2", "pval.exp.pop2", "pval.out.pop2")

dat2 = paras2$dat
colnames(dat2) = c("SNP", "A1", "A2", "b.exp.pop2", "b.exp.pop1", "se.exp.pop2", "se.exp.pop1", "pval.exp.pop2", "pval.exp.pop1")

dat = merge(dat1, dat2[, c("SNP", "A1", "A2", "b.exp.pop1", "se.exp.pop1", "pval.exp.pop1")])
dat = merge(dat, ldscore0, by = "SNP")

clumped.xpMR = MRAPSS::clump(subset(dat, pval.exp.pop1 < threshold),
                              IV.Threshold = 1,
                              SNP_col = "SNP",
                              pval_col = "pval.exp.pop2",
                              clump_kb = 1000,
                              clump_r2 = 0.001,
                              bfile = paste0(dirname(base_dir),"/nc/1kg_pops/EAS_phase3_hg38_autosomes"),
                              plink_bin = paste0(dirname(base_dir),"/nc/plink_Linux"))
clumped.xpMR$Threshold = threshold

############# run XPMR ################

t0 = abs(qnorm(threshold / 2))
dt = 0.13 / (sqrt(N2 / N1))
faket = 2 * (1 - pnorm(abs(t0 + dt)))

if(dim(clumped.xpMR)[1] < 100){
  XPMR_res = fit_XMR(data = clumped.xpMR, C = C, Omega0 = Omega0, Threshold = faket, tol1 = 1e-08,
                      tol2 = 1e-08, min_thres = 1e-2)
} else {
  XPMR_res = fit_XMR(data = clumped.xpMR, C = C, Omega0 = Omega0, Threshold = faket, tol1 = 1e-07,
                      tol2 = 1e-07, min_thres = 1e-2)
}

res_xmr = data.frame(exposure = ts1[t],
                      outcome = ts2[s],
                      beta = XPMR_res$beta,
                      beta.se = XPMR_res$beta.se,
                      beta.pvalue = XPMR_res$beta.pvalue,
                      tau.sq = XPMR_res$tau.sq,
                      sigmax11 = XPMR_res$SigmaX[1,1],
                      sigmax12 = XPMR_res$SigmaX[1,2],
                      sigmax22 = XPMR_res$SigmaX[2,2],
                      nIV = XPMR_res$nIV,
                      IVsignal.sum = XPMR_res$post$IVsignal.sum,
                      nvalid = XPMR_res$nvalid,
                      pi0 = XPMR_res$pi0,
                      Threshold = threshold,
                      exp_N1 = N1,
                      exp_N2 = N2,
                      modified_threshold = faket)

# ## MRAPSS, Egger, Weighted-mode

ts1 = '30830'
ts2 = '250.2'

s = 1
t = 1 

print(paste('outcome:', ts2[s]))
print(paste('exposure:', ts1[t]))


exp.pop2 = fread(paste0(base_dir, "/formatted_data/EAS_", ts1[t]))
out.pop2 = fread(paste0(base_dir, "/formatted_data/EAS_", ts2[s]))

###
threshold = 5e-5

paras1 = align_data(dat1 = exp.pop2,
                   dat2 = out.pop2,
                   trait1.name = paste0("EAS_", ts1[t]),
                   trait2.name = paste0("EAS_", ts2[s]),
                   h2.fix.intercept = F,
                   LDSC = F
                   )


dat = paras1$dat
colnames(dat) = c("SNP", "A1", "A2", "b.exp.pop2", "b.out.pop2", "se.exp.pop2", "se.out.pop2", "pval.exp.pop2", "pval.out.pop2" )

dat = merge(dat, ldscore0, by="SNP")
dat$Threshold = threshold


clumped.MR = MRAPSS::clump(dat,
                           IV.Threshold = threshold,
                           SNP_col = "SNP",
                           pval_col = "pval.exp.pop2",
                           clump_kb = 1000,
                           clump_r2 = 0.001,
                           bfile = paste0(dirname(base_dir),"/nc/1kg_pops/EAS_phase3_hg38_autosomes"),
                           plink_bin = paste0(dirname(base_dir),"/nc/plink_Linux"))
clumped.MR$Threshold = threshold


print('MR-APSS')

############# run MR-APSS #################
MRdat = clumped.MR[, c( "SNP", "A1", "A2", "b.exp.pop2", "se.exp.pop2", "pval.exp.pop2",
                                           "b.out.pop2", "se.out.pop2", "pval.out.pop2", "L2.pop2", "Threshold")]
colnames(MRdat) = c( "SNP", "A1", "A2", "b.exp", "se.exp", "pval.exp", "b.out", "se.out", "pval.out", "L2", "Threshold")

fit.apss = MRAPSS::MRAPSS(MRdat, C=C[2:3,2:3], Omega = Omega0[2:3,2:3], Cor.SelectionBias = T)

res_mrapss = data.frame(exposure=ts1[t],
                       outcome = ts2[s],
                       beta = fit.apss$beta,
                       beta.se = fit.apss$beta.se,
                       beta.pvalue = fit.apss$pvalue,
                       tau.sq = fit.apss$tau.sq,
                       nIV = nrow(MRdat),
                       pi0 = fit.apss$pi0,
                       IVsignal.sum = fit.apss$post$IVsignal.sum,
                       Threshold = fit.apss$Threshold)


betaXG = clumped.MR[, 'b.exp.pop2']
betaYG = clumped.MR[, 'b.out.pop2']
sebetaXG = clumped.MR[, 'se.exp.pop2']
sebetaYG = clumped.MR[, 'se.out.pop2']

print('Egger')
egger_result = TwoSampleMR::mr_egger_regression(b_exp = betaXG, b_out = betaYG, 
                                   se_exp = sebetaXG, se_out =sebetaYG)

print('Weighted-mode')
weightmod_result = TwoSampleMR::mr_weighted_mode(b_exp = betaXG, b_out = betaYG, 
                                            se_exp = sebetaXG, se_out =sebetaYG)


res_egger = data.frame(exposure = ts1[t],
                       outcome = ts2[s],
                       beta = egger_result$b,
                       beta.se = egger_result$se,
                       beta.pvalue = egger_result$pval,
                       nIV = dim(clumped.MR)[1],
                       Threshold = threshold)

res_weightmod = data.frame(exposure = ts1[t],
                           outcome = ts2[s],
                           beta = weightmod_result$b,
                           beta.se = weightmod_result$se,
                           beta.pvalue = weightmod_result$pval,
                           nIV = dim(clumped.MR)[1],
                           Threshold = threshold)

res_egger

res_weightmod

# ### CAUSE

library(cause)

threshold = 5e-5

ts1 = '30830'
ts2 = '250.2'

s = 1
t = 1 

print(paste('outcome:', ts2[s]))
print(paste('exposure:', ts1[t]))

X1 = fread(paste0(base_dir, "/formatted_data/EAS_", ts1[t]))
X2 = fread(paste0(base_dir, "/formatted_data/EAS_", ts2[s]))

X1$b = X1$Z/sqrt(X1$N)
X2$b = X2$Z/sqrt(X2$N)
X1$se = 1/sqrt(X1$N)
X2$se = 1/sqrt(X2$N)

X <- try(gwas_merge(X1, X2, 
                    snp_name_cols = c("SNP", "SNP"),
                    beta_hat_cols = c("b", "b"),
                    se_cols = c("se", "se"),
                    A1_cols = c("A1", "A1"),
                    A2_cols = c("A2", "A2")))


d0 = X1[, c("SNP", "P")]
colnames(d0) = c("snp", "pval.exp")
X0 = merge(X, d0, by="snp")

#############clumping#######################    
# pre-process: remove ambiguous SNPs
exp_dat1 <- X1[X1$A1 != comple(X1$A2),]
out_dat1 <- X2[X2$A1 != comple(X2$A2),]

paras1 = align_data(dat1 = exp_dat1,
                   dat2 = out_dat1,
                   trait1.name = paste0("EAS_", ts1[t]),
                   trait2.name = paste0("EAS_", ts2[s]),
                   h2.fix.intercept = F,
                   LDSC = F
               )

dat = paras1$dat
colnames(dat) = c("SNP", "A1", "A2", "b.exp", "b.out", "se.exp", "se.out", "pval.exp", "pval.out")

dat = merge(dat, ldscore0, by="SNP")

clumped = MRAPSS::clump(dat,
                           IV.Threshold = threshold,
                           SNP_col = "SNP",
                           pval_col = "pval.exp",
                           clump_kb = 1000,
                           clump_r2 = 0.001,
                           bfile = paste0(dirname(base_dir),"/nc/1kg_pops/EAS_phase3_hg38_autosomes"),
                           plink_bin = paste0(dirname(base_dir),"/nc/plink_Linux"))
clumped$Threshold = threshold


SNPs = clumped$SNP

clumped0 = subset(X0, snp %in% SNPs)

varlist <- with(X, sample(snp, size=min(nrow(X), 1000000), replace=FALSE))
params <- try(est_cause_params(X, varlist))


clumped = subset(clumped0, pval.exp <= threshold)

# cause
if(!is.null(clumped)){
    top_ldl_pruned_vars =intersect(as.character(X$snp), as.character(clumped$snp))

    cause_res <- try(cause(X=X, variants = top_ldl_pruned_vars , param_ests = params, force=TRUE))


    res_elpd <- data.frame(ts1[t], ts2[s],
                           threshold,
                           length(top_ldl_pruned_vars),
                           cause_res$elpd)

    res.cause.est = summary(cause_res, ci_size=0.95)

    res_est = data.frame(ts1[t], ts2[s],
                         threshold,length(top_ldl_pruned_vars),
                         matrix(c(res.cause.est$quants[[2]][,1],
                                  res.cause.est$quants[[2]][,2],
                                  res.cause.est$quants[[2]][,3]), nrow=1))

    colnames(res_elpd) = c("exposure","outcome","Threshold","nsnp","model1","model2","delta_elpd", "se_delta_elpd", "Z")
    colnames(res_est) = c("exposure","outcome","Threshold","nsnp", "beta.hat","b_l","b_u","eta","eta_l","eta_u","q","q_l","q_u")
    res_elpd = unique(subset(res_elpd, model1=="sharing"&model2=="causal"))
    res_elpd$pval = pnorm(res_elpd$Z)
    res_est = unique(res_est[, c("exposure","outcome","Threshold","nsnp", "beta.hat","b_l","b_u")])
    res_est$se = (res_est$b_u - res_est$b_l)/2/1.96
    cause_res = unique(merge(unique(res_elpd[, c("exposure","outcome","pval")]),
                             res_est[, c("exposure","outcome","Threshold","nsnp", "beta.hat","se")],
                             by=c("exposure","outcome")))
    cause_res$Method = "CAUSE"

}

print(cause_res)



case_study = NULL
case_study = rbind(case_study, 
                   res_xmr[,c('exposure','outcome','beta','beta.se','beta.pvalue','nIV','nvalid','Threshold')])
case_study = cbind(case_study, data.frame(method='XMR',population='EAS'))

case_study = rbind(case_study, data.frame(exposure = res_mrapss$exposure, 
                                          outcome = res_mrapss$outcome, 
                                          beta = res_mrapss$beta,
                                          beta.se = res_mrapss$beta.se,
                                          beta.pvalue = res_mrapss$beta.pvalue,
                                          nIV =  res_mrapss$nIV,
                                          nvalid = res_mrapss$nIV*res_mrapss$pi0, 
                                          Threshold = res_mrapss$Threshold,
                                          method = "MRAPSS", 
                                          population = "EAS"))

case_study = rbind(case_study, data.frame(exposure = res_egger$exposure, 
                                          outcome = res_egger$outcome, 
                                          beta = res_egger$beta,
                                          beta.se = res_egger$beta.se,
                                          beta.pvalue = res_egger$beta.pvalue,
                                          nIV =  res_egger$nIV,
                                          nvalid = res_egger$nIV, 
                                          Threshold = res_egger$Threshold,
                                          method = "Egger", 
                                          population = "EAS"))

case_study = rbind(case_study, data.frame(exposure = res_weightmod$exposure, 
                                          outcome = res_weightmod$outcome, 
                                          beta = res_weightmod$beta,
                                          beta.se = res_weightmod$beta.se,
                                          beta.pvalue = res_weightmod$beta.pvalue,
                                          nIV =  res_weightmod$nIV,
                                          nvalid = res_weightmod$nIV, 
                                          Threshold = res_weightmod$Threshold,
                                          method = "Weighted-mode", 
                                          population = "EAS"))

case_study = rbind(case_study, data.frame(exposure = cause_res$exposure, 
                                          outcome = cause_res$outcome, 
                                          beta = cause_res$beta.hat,
                                          beta.se = cause_res$se,
                                          beta.pvalue = cause_res$pval,
                                          nIV = cause_res$nsnp,
                                          nvalid = cause_res$nsnp, 
                                          Threshold = cause_res$Threshold,
                                          method = "CAUSE", 
                                          population = "EAS"))

case_study

write.csv(case_study, 'case_study_EAS.csv')


