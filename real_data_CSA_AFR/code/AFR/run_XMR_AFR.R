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

setwd("./real_data_CSA_AFR")

library(XMR)

# ### step1: prepare GWAS data (finished，with prepared dara in ./formatted_data; details in format_data.ipynb)

# ### step2: estimate bgparas (background parameters)

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

info_afr = read.csv('AFR_data_info.csv', row.names = 1)

data_summary_exp = info_afr[is.na(info_afr$AFR.controls),]
data_summary_out = info_afr[!is.na(info_afr$AFR.controls),]

dir.create("./AFR_bgparas/", recursive = TRUE)

ts1 = data_summary_exp$EUR.phecode
ts2 = data_summary_out$EUR.phecode

pairs <- rbind(
    expand.grid(exposure = ts1, outcome = ts2, stringsAsFactors = FALSE),
    expand.grid(exposure = ts2, outcome = ts1, stringsAsFactors = FALSE),
    subset(expand.grid(exposure = ts2, outcome = ts2, stringsAsFactors = FALSE), exposure != outcome)
)


for (i in 1:nrow(pairs)) {
    exp_code <- pairs$exposure[i]
    out_code <- pairs$outcome[i]
    
    print(paste0(i, "/", nrow(pairs), " exposure: ", exp_code, " -> outcome: ", out_code))

    bgparas <- estimate_bgparas(
                    ldscore0      = ldscore0,
                    exposure_code = exp_code, 
                    outcome_code  = out_code,      
                    pop_ref       = "EUR",         
                    pop_target    = "AFR",  
                    snp_col       = "rsid",
                    gwas_dir      = "./formatted_data" 
                )

    write.table(bgparas, file = paste0("./AFR_bgparas/bgparas_", exp_code, "_", out_code), quote = F, col.names = F, row.names = F)
                           
}



# ### step3: LD clumping; run XMR

library(XMAP)
library(data.table)
library(Matrix)

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

ldscore = ldscore[, c("rsid", "EUR", "AFR", "EUR_AFR")]
colnames(ldscore) = c("SNP", "L2.pop1", "L2.pop2", "L12")

ldscore0 = ldscore

library(doParallel)
library(data.table)
library(reshape2)
library(MRAPSS)
library(XMAP) 
library(dplyr)
library(magrittr)
library(Matrix)

info_afr = read.csv('AFR_data_info.csv', row.names = 1)
data_summary_exp = info_afr[is.na(info_afr$AFR.controls),]
data_summary_out = info_afr[!is.na(info_afr$AFR.controls),]

ts1 = data_summary_exp$EUR.phecode
ts2 = data_summary_out$EUR.phecode

pairs <- rbind(
    expand.grid(exposure = ts1, outcome = ts2, stringsAsFactors = FALSE),
    expand.grid(exposure = ts2, outcome = ts1, stringsAsFactors = FALSE),
    subset(expand.grid(exposure = ts2, outcome = ts2, stringsAsFactors = FALSE), exposure != outcome)
)


temp_dir <- "./AFR_results/XMR"
dir.create(temp_dir, recursive = TRUE, showWarnings = FALSE)

# all parameter combinations
param_grid <- expand.grid(
  threshold = c(5e-5),
  l = 1:dim(pairs)[1],
  stringsAsFactors = FALSE
)


n_cores <- detectCores() - 1
cl <- makeCluster(n_cores)
registerDoParallel(cl)
print(cl)

base_dir <- getwd()

res_all = foreach(i = 1:nrow(param_grid),
        .packages = c("data.table", "reshape2", "MRAPSS", "XMAP", "dplyr", "magrittr", "Matrix", "XMR")
        #.export = c("ts1", "ts2", "base_dir", "param_grid", "temp_dir") #"ldscore0",
) %dopar% {

  setwd(base_dir)

  threshold <- param_grid$threshold[i]
  l <- param_grid$l[i]

  s = 1
  t = 1
    
  ts1 = pairs$exposure[l]
  ts2 = pairs$outcome[l]

  task_label <- paste0("threshold=", threshold, " | outcome=", ts2[s], " | exposure=", ts1[t])

  out_file <- file.path(temp_dir, paste0("/res_", threshold, "_", ts2[s], "_", ts1[t], ".MRres"))

  if(file.exists(out_file)) return(NULL)

  result <- tryCatch({
      
    paras = read.table(paste0(base_dir, "/AFR_bgparas/bgparas_", ts1[t], "_", ts2[s]))
    colnames(paras) = c('trait1','trait2','pop1','pop2','Omega11','Omega12','Omega22','c1','c12','c2')

    bg <- qc_bgparas(
              bgparas       = paras,
              exposure_code = ts1[t],
              outcome_code  = ts2[s],
              pop_ref       = "EUR",
              pop_target    = "AFR",
              truncation_corr = 0.95
              )

    Omega0 <- bg$Omega0
    C      <- bg$C

      
    exp.pop1 = fread(paste0(base_dir, "/formatted_data/EUR_", ts1[t]))
    exp.pop2 = fread(paste0(base_dir, "/formatted_data/AFR_", ts1[t]))
    out.pop2 = fread(paste0(base_dir, "/formatted_data/AFR_", ts2[s]))

    N1 = median(exp.pop1$N)
    N2 = median(exp.pop2$N)

    paras1 = align_data(dat1 = exp.pop2,
                        dat2 = out.pop2,
                        trait1.name = paste0("AFR_", ts1[t]),
                        trait2.name = paste0("AFR_", ts2[s]),
                        h2.fix.intercept = F,
                        LDSC = F)

    paras2 = align_data(dat1 = exp.pop2,
                        dat2 = exp.pop1,
                        trait1.name = paste0("AFR_", ts1[t]),
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
                                  bfile = paste0(dirname(base_dir),"/nc/1kg_pops/AFR_phase3_hg38_autosomes"),
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

    write.table(res_xmr, file = out_file, quote=F, col.names = F, append = T, row.names = F)
      
    return(list(status = "success", task = task_label, message = "OK"))

  }, error = function(e){
    return(list(status = "error", task = task_label, message = conditionMessage(e)))
  })

  return(result)
}


stopCluster(cl)

# ========== merge all results into one file and delete those files ==========

all_files <- list.files(temp_dir, pattern = "^res_.*\\.MRres$", full.names = TRUE)

if(length(all_files) > 0){
  all_results <- do.call(rbind, lapply(all_files, fread))
  colnames(all_results) = c("exposure", "outcome", "beta", "beta.se", "beta.pvalue", "tau.sq", 
                            "sigmax11", "sigmax12", "sigmax22", "nIV", "IVsignal.sum", "nvalid", 
                            "pi0", "Threshold", "exp_N1", "exp_N2", "modified_threshold")

  if(nrow(all_results) > 0){
    write.table(all_results,
                file = paste0("./AFR_results/XMR.MRres"),
                quote = FALSE, col.names = FALSE, row.names = FALSE)
  }

  cat("number of experiments: ", nrow(all_results), "\n")
} else {
  cat("no file found\n")
}

unlink(temp_dir, recursive = TRUE)




