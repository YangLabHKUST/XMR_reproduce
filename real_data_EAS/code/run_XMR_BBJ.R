# -*- coding: utf-8 -*-
# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .R
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.15.2
#   kernelspec:
#     display_name: R
#     language: R
#     name: ir
# ---

# +
setwd("./real_data_EAS")

library(XMR)
# -

# ### step1: prepare GWAS data (finished，with prepared dara in ./formatted_data; details in format_data.ipynb)

# ### step2: estimate bgparas (background parameters)

# +
library(XMAP)
library(data.table)
library(Matrix)

set.seed(1)

# read LD scores
ldscore <- fread("./LDscore_eas_brit_afr.txt")

# pre-process: remove ambiguous SNPs
idx_amb <- which(ldscore$allele1 == comple(ldscore$allele2))
ldscore <- ldscore[-idx_amb,]

ldscore$BBJ = ldscore$EAS
ldscore$BBJ_EUR = ldscore$EUR_EAS
ldscore$EUR_BBJ = ldscore$EUR_EAS
ldscore$EUR_EUR = ldscore$EUR
ldscore$BBJ_BBJ = ldscore$EAS

ldscore = as.data.frame(ldscore)
ldscore0 = ldscore

# +
dir.create("./BBJ_bgparas/", recursive = TRUE)

ts1 = data_summary_exp$UKB.code
ts2 = data_summary_out$phecode_unique

for(s in 1:length(ts2)){
    
    print(paste0("outcome: ",ts2[s]))
    
    for(t in 1:length(ts1)){
    
        print(paste0("trait: ",ts1[t]))

        result <- try({
            
            pairs = data.frame(trait1 = c(paste0("EUR_", ts1[t]), paste0("EUR_", ts1[t]), paste0("EUR_", ts1[t]), paste0("BBJ_", ts1[t]), paste0("BBJ_", ts1[t]), paste0("BBJ_", ts2[s])),
                               trait2 = c(paste0("EUR_", ts1[t]), paste0("BBJ_", ts1[t]), paste0("BBJ_", ts2[s]), paste0("BBJ_", ts1[t]), paste0("BBJ_", ts2[s]), paste0("BBJ_", ts2[s])),
                               pop1 = c("EUR", "EUR", "EUR", "BBJ", "BBJ", "BBJ"), 
                               pop2 = c("EUR",  "BBJ", "BBJ", "BBJ", "BBJ", "BBJ") )



            for(i in 1: nrow(pairs)){

               trait1 = as.character(pairs[i, 1])
               trait2 = as.character(pairs[i, 2])
               pop1 = as.character(pairs[i, 3])
               pop2 = as.character(pairs[i, 4])

                # read GWAS summary statistics
                dat1 <- fread(paste0("./formatted_data/", trait1))

                # read GWAS summary statistics
                dat2<- fread(paste0("./formatted_data/", trait2))

                # pre-process: overlap SNPs
                snps <- Reduce(intersect, list(ldscore0$rsid, dat1$SNP, dat2$SNP))
                dat1_ldsc <- dat1[match(snps, dat1$SNP),]
                dat2_ldsc <- dat2[match(snps, dat2$SNP),]
                ldscore <- ldscore0[match(snps, ldscore0$rsid),]

                # pre-process: flip alleles; z-score (beta/se)
                z_pop1 <- dat1_ldsc$Z
                z_pop2 <- dat2_ldsc$Z

                idx_flip <- which(dat1_ldsc$A1 != ldscore$allele1 & dat1_ldsc$A1 != comple(ldscore$allele1))
                z_pop1[idx_flip] <- -z_pop1[idx_flip]

                idx_flip <- which(dat2_ldsc$A1 != ldscore$allele1 & dat2_ldsc$A1 != comple(ldscore$allele1))
                z_pop2[idx_flip] <- -z_pop2[idx_flip]

                # select SNPs with smaller Z values; filter out extreme values
                idx1 <- which(z_pop1^2 < 30 & z_pop2^2 < 30)
                # weights of the two traits in their respective populations are calculated based on the ldscore, ensuring that SNPs with larger LD scores have a smaller impact on the regression.
                ld_pop1_w <- 1 / sapply(ldscore[, pop1], function(x) max(x, 1))
                ld_pop2_w <- 1 / sapply(ldscore[, pop2], function(x) max(x, 1))

                # bi-variate LDSC:
                # stage 1 of LDSC: estimate intercepts
                fit_step1 <- XMAP::estimate_gc(data.frame(Z = z_pop1[idx1], N = dat1_ldsc$N[idx1]), 
                                                 data.frame(Z = z_pop2[idx1], N = dat2_ldsc$N[idx1]),
                                                 ldscore[idx1, pop1], ldscore[idx1, pop2], ldscore[idx1, paste0(pop1, "_", pop2)],
                                                 reg_w1 = ld_pop1_w[idx1], reg_w2 = ld_pop2_w[idx1], reg_wx = sqrt(ld_pop1_w[idx1] * ld_pop2_w[idx1]),
                                                 constrain_intercept = F)

                # stage 2 of LDSC: fix intercepts and estimate slopes
                fit_step2 <- XMAP::estimate_gc(data.frame(Z = z_pop1, N = dat1_ldsc$N), 
                                                 data.frame(Z = z_pop2, N = dat2_ldsc$N),
                                                 ldscore[, pop1], ldscore[, pop2], ldscore[, paste0(pop1, "_", pop2)],
                                                 reg_w1 = ld_pop1_w, reg_w2 = ld_pop2_w, reg_wx = sqrt(ld_pop1_w * ld_pop2_w),
                                                 constrain_intercept = T, fit_step1$tau1$coefs[1], fit_step1$tau2$coefs[1], fit_step1$theta$coefs[1])

                write.table(data.frame(trait1, trait2, pop1, pop2, omega11=fit_step2$tau1$coefs[2], omega12=fit_step2$theta$coefs[2], omega22=fit_step2$tau2$coefs[2], 
                                       c1= fit_step2$tau1$coefs[1], c12=fit_step2$theta$coefs[1] , c2=fit_step2$tau2$coefs[1]),
                            file = paste0("./BBJ_bgparas/bgparas_", ts1[t], "_", ts2[s]), quote = F, col.names = F, row.names = F, append = T)

            }

        }, silent = TRUE)

        if (inherits(result, "try-error")) {
            print(paste("Error in round", t, "- skipping to next round."))
            next
        }
    }                           
}
# -

# ### step3: LD clumping; run XMR

# +
library(XMAP)
library(data.table)
library(Matrix)

# read LD scores
ldscore <- fread("./LDscore_eas_brit_afr.txt")

# pre-process: remove ambiguous SNPs
idx_amb <- which(ldscore$allele1 == comple(ldscore$allele2))
ldscore <- ldscore[-idx_amb,]

ldscore$BBJ_EUR = ldscore$EUR_EAS
ldscore$EUR_BBJ = ldscore$EUR_EAS
ldscore$EUR_EUR = ldscore$EUR
ldscore$BBJ_BBJ = ldscore$EAS

ldscore = as.data.frame(ldscore)

ldscore = ldscore[, c("rsid", "BBJ_BBJ", "EUR_EUR", "BBJ_EUR")]
colnames(ldscore) = c("SNP", "L2.pop2", "L2.pop1", "L12")

ldscore0 = ldscore

# +
library(doParallel)
library(data.table)
library(reshape2)
library(MRAPSS)
library(XMAP) 
library(dplyr)
library(magrittr)
library(Matrix)

file_path <- './raw_data/BBJ/EUR_BBJ_exposure_description.csv'
data_summary_exp <- read.csv(file_path)
data_summary_exp$UKB.code = as.character(data_summary_exp$UKB.code)

data_summary_out = read.csv('./raw_data/BBJ/EUR_BBJ_outcome_description.csv')
data_summary_out$UKB.phecode = as.character(data_summary_out$UKB.phecode)
data_summary_out$phecode_unique = as.character(data_summary_out$phecode_unique)

ts1 = data_summary_exp$UKB.code
ts2 = data_summary_out$phecode_unique


temp_dir <- "./BBJ_results/XMR"
dir.create(temp_dir, recursive = TRUE, showWarnings = FALSE)

# all parameter combinations
param_grid <- expand.grid(
  threshold = c(5e-5),
  s = 1:length(ts2),
  t = 1:length(ts1),
  stringsAsFactors = FALSE
)

n_cores <- detectCores() - 1
cl <- makeCluster(n_cores)
registerDoParallel(cl)
print(cl)

base_dir <- getwd()

res_all = foreach(i = 1:nrow(param_grid),
        .packages = c("data.table", "reshape2", "MRAPSS", "XMAP", "dplyr", "magrittr", "Matrix")
        #.export = c("ts1", "ts2", "base_dir", "param_grid", "temp_dir") #"ldscore0",
) %dopar% {

  setwd(base_dir)

  threshold <- param_grid$threshold[i]
  s <- param_grid$s[i]
  t <- param_grid$t[i]

  task_label <- paste0("threshold=", threshold, " | outcome=", ts2[s], " | exposure=", ts1[t])

  out_file <- file.path(temp_dir, paste0("/res_", threshold, "_", ts2[s], "_", ts1[t], ".MRres"))

  if(file.exists(out_file)) return(NULL)

  result <- tryCatch({

    paras = read.table(paste0(base_dir, "/BBJ_bgparas/bgparas_", ts1[t], "_", ts2[s]))
    colnames(paras) = c("trait1", "trait2", "pop2", "pop1", "Omega11", "Omega12","Omega22", "c1","c12","c12")
    paras = subset(paras, trait1 %in% c(paste0("EUR_", ts1[t]), paste0("BBJ_", ts1[t]), paste0("BBJ_", ts2[s])) & 
                          trait2 %in% c(paste0("EUR_", ts1[t]), paste0("BBJ_", ts1[t]), paste0("BBJ_", ts2[s])))[,  c("trait1", "trait2", "Omega12","c12")]
      
    paras_opt = paras[, c("trait2", "trait1", "Omega12","c12")]
    colnames(paras_opt) = c("trait1", "trait2", "Omega12","c12")
    paras_opt = subset(paras_opt, trait1 != trait2)
    paras_all = rbind(paras, paras_opt)

    Omega0 = reshape2::acast(unique(paras_all[, c("trait1", "trait2", "Omega12")]), trait1 ~ trait2, value.var = "Omega12")
    C = reshape2::acast(unique(paras_all[, c("trait1", "trait2", "c12")]), trait1 ~ trait2, value.var = "c12")

    Omega0 = Omega0[c(paste0("EUR_", ts1[t]), paste0("BBJ_", ts1[t]), paste0("BBJ_", ts2[s])),
                     c(paste0("EUR_", ts1[t]), paste0("BBJ_", ts1[t]), paste0("BBJ_", ts2[s]))]
    C = C[c(paste0("EUR_", ts1[t]), paste0("BBJ_", ts1[t]), paste0("BBJ_", ts2[s])),
          c(paste0("EUR_", ts1[t]), paste0("BBJ_", ts1[t]), paste0("BBJ_", ts2[s]))]

      
    ### quality control of background parameters
    if((Omega0[1,1] < 0) | (Omega0[2,2] < 0) | (Omega0[3,3] < 0)){
      print('Omega variance negative')
    }

    flag2 <- TRUE
    tryCatch(
      {chol(C)},
      error = function(e){
        flag2 <<- FALSE
      })

    if(!flag2){
      print('C not PSD')
    }

    truncation_corr = 0.95

    if(abs(Omega0[1,3] / sqrt(Omega0[1,1] * Omega0[3,3])) > 1){
      Omega0[1,3] = sqrt(Omega0[1,1] * Omega0[3,3]) * truncation_corr * sign(Omega0[1,3])
      Omega0[3,1] = Omega0[1,3]
    }

    if(abs(Omega0[2,3] / sqrt(Omega0[2,2] * Omega0[3,3])) > 1){
      Omega0[2,3] = sqrt(Omega0[2,2] * Omega0[3,3]) * truncation_corr * sign(Omega0[2,3])
      Omega0[3,2] = Omega0[2,3]
    }

    if(abs(Omega0[1,2] / sqrt(Omega0[1,1] * Omega0[2,2])) > 1){
      Omega0[1,2] = sqrt(Omega0[1,1] * Omega0[2,2]) * truncation_corr * sign(Omega0[1,2])
      Omega0[2,1] = Omega0[1,2]
    }

    exp.pop1 = fread(paste0(base_dir, "/formatted_data/EUR_", ts1[t]))
    exp.pop2 = fread(paste0(base_dir, "/formatted_data/BBJ_", ts1[t]))
    out.pop2 = fread(paste0(base_dir, "/formatted_data/BBJ_", ts2[s]))

    N1 = median(exp.pop1$N)
    N2 = median(exp.pop2$N)

    paras1 = align_data(dat1 = exp.pop2,
                        dat2 = out.pop2,
                        trait1.name = paste0("BBJ_", ts1[t]),
                        trait2.name = paste0("BBJ_", ts2[s]),
                        h2.fix.intercept = F,
                        LDSC = F)

    paras2 = align_data(dat1 = exp.pop2,
                        dat2 = exp.pop1,
                        trait1.name = paste0("BBJ_", ts1[t]),
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

    write.table(res_xmr, file = out_file, quote=F, col.names = F, append = T, row.names = F)
      
    return(list(status = "success", task = task_label, message = "OK"))

  }, error = function(e){
    return(list(status = "error", task = task_label, message = conditionMessage(e)))
  })

  return(result)
}


stopCluster(cl)

# +
# ========== merge all results into one file and delete those files ==========

all_files <- list.files(temp_dir, pattern = "^res_.*\\.MRres$", full.names = TRUE)

if(length(all_files) > 0){
  all_results <- do.call(rbind, lapply(all_files, fread))
  colnames(all_results) = c("exposure", "outcome", "beta", "beta.se", "beta.pvalue", "tau.sq", 
                            "sigmax11", "sigmax12", "sigmax22", "nIV", "IVsignal.sum", "nvalid", 
                            "pi0", "Threshold", "exp_N1", "exp_N2", "modified_threshold")

  if(nrow(all_results) > 0){
    write.table(all_results,
                file = paste0("./BBJ_results/XMR.MRres"),
                quote = FALSE, col.names = FALSE, row.names = FALSE)
  }

  cat("number of experiments: ", nrow(all_results), "\n")
} else {
  cat("no file found\n")
}

unlink(temp_dir, recursive = TRUE)
# -



# ### IV number when selected solely by EAS

# +
base_dir <- getwd()

temp_dir <- "./BBJ_results"
dir.create(temp_dir, recursive = TRUE, showWarnings = FALSE)

library(doParallel)
library(data.table)
library(reshape2)
library(MRAPSS)
library(XMAP) 
library(dplyr)
library(magrittr)
library(Matrix)

file_path <- './raw_data/BBJ/EUR_BBJ_exposure_description.csv'
data_summary_exp <- read.csv(file_path)
data_summary_exp$UKB.code = as.character(data_summary_exp$UKB.code)

data_summary_out = read.csv('./raw_data/BBJ/EUR_BBJ_outcome_description.csv')
data_summary_out$UKB.phecode = as.character(data_summary_out$UKB.phecode)
data_summary_out$phecode_unique = as.character(data_summary_out$phecode_unique)

ts1 = data_summary_exp$UKB.code
ts2 = data_summary_out$phecode_unique


temp_dir <- "./BBJ_results/IV_num"
dir.create(temp_dir, recursive = TRUE, showWarnings = FALSE)

# all parameter combinations
param_grid <- expand.grid(
  threshold = c(5e-5, 5e-6, 5e-7, 5e-8),
  s = 1:length(ts2),
  t = 1:length(ts1),
  stringsAsFactors = FALSE
)

n_cores <- detectCores() - 1
cl <- makeCluster(n_cores)
registerDoParallel(cl)
print(cl)

base_dir <- getwd()

res_all = foreach(i = 1:nrow(param_grid),
        .packages = c("data.table", "reshape2", "MRAPSS", "XMAP", "dplyr", "magrittr", "Matrix")
        #.export = c("ts1", "ts2", "base_dir", "param_grid", "temp_dir") #"ldscore0",
) %dopar% {

  setwd(base_dir)

  threshold <- param_grid$threshold[i]
  s <- param_grid$s[i]
  t <- param_grid$t[i]

  out_file <- file.path(temp_dir, paste0("/res_", threshold, "_", ts2[s], "_", ts1[t], ".MRres"))

    exp.pop2 = fread(paste0(base_dir, "/formatted_data/BBJ_", ts1[t]))
    out.pop2 = fread(paste0(base_dir, "/formatted_data/BBJ_", ts2[s]))


    paras1 = align_data(dat1 = exp.pop2,
                       dat2 = out.pop2,
                       trait1.name = paste0("BBJ_", ts1[t]),
                       trait2.name = paste0("BBJ_", ts2[s]),
                       h2.fix.intercept = F,
                       LDSC = F
                       )

    dat = paras1$dat
    colnames(dat) = c("SNP", "A1", "A2", "b.exp.pop2", "b.out.pop2", "se.exp.pop2", "se.out.pop2", "pval.exp.pop2", "pval.out.pop2" )

    dat = merge(dat, ldscore0, by="SNP")
    dat$Threshold = threshold


    result <- tryCatch({
        clumped.MR = MRAPSS::clump(dat,
                                   IV.Threshold = threshold,
                                   SNP_col = "SNP",
                                   pval_col = "pval.exp.pop2",
                                   clump_kb = 1000,
                                   clump_r2 = 0.001,
                                   bfile = paste0(dirname(base_dir),"/nc/1kg_pops/EAS_phase3_hg38_autosomes"),
                                   plink_bin = paste0(dirname(base_dir),"/nc/plink_Linux"))
        clumped.MR$Threshold = threshold

        write.table(data.frame(exposure=ts1[t],
                               outcome=ts2[s],
                               Threshold=threshold,
                               nIV=dim(clumped.MR)[1]),
                    file=out_file, append=T, col.names = F, row.names = F, quote = F)

        TRUE
    }, error = function(e) {
        # Print the actual error message before skipping
        print(paste("Error for exposure", ts1[t], "outcome", ts2[s],
                    "threshold", threshold, ":", e$message))
        FALSE
    })

    if(!result) next

}

# +
# ========== merge all results into one file and delete those files ==========

all_files <- list.files(temp_dir, pattern = "^res_.*\\.MRres$", full.names = TRUE)

if(length(all_files) > 0){
  all_results <- do.call(rbind, lapply(all_files, fread))
  colnames(all_results) = c("exposure", "outcome", "Threshold", "nIV")

  if(nrow(all_results) > 0){
    write.table(all_results,
                file = paste0("./BBJ_results/nIV_by_EAS.MRres"),
                quote = FALSE, col.names = FALSE, row.names = FALSE)
  }

  cat("number of experiments: ", nrow(all_results), "\n")
} else {
  cat("no file found\n")
}

unlink(temp_dir, recursive = TRUE)
# -


