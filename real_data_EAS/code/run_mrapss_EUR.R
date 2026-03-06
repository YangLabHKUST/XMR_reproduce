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
library(XMR)

setwd("./real_data_EAS")

# +
file_path <- './raw_data/BBJ/EUR_BBJ_exposure_description.csv'
data_summary_exp <- read.csv(file_path)
data_summary_exp$UKB.code = as.character(data_summary_exp$UKB.code)

data_summary_out = read.csv('./raw_data/BBJ/EUR_BBJ_outcome_description.csv')
data_summary_out$UKB.phecode = as.character(data_summary_out$UKB.phecode)
data_summary_out$phecode_unique = as.character(data_summary_out$phecode_unique)

ts1 = data_summary_exp$UKB.code
ts2 = data_summary_out$phecode_unique
# -

library(data.table)
library(reshape2)
library(MRAPSS)
library(XMAP) 
library(dplyr)
library(magrittr)
library(Matrix)
library(MASS)
library(TEMR)

# ### MRAPSS

# #### estimate background parameters

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
dir.create("./EUR_bgparas/", recursive = TRUE)

ts1 = data_summary_exp$UKB.code
ts2 = data_summary_out$phecode_unique

for(s in 1:length(ts2)){
    
    print(paste0("outcome: ",ts2[s]))
    
    for(t in 1:length(ts1)){
    
        print(paste0("trait: ",ts1[t]))

        result <- try({
            
            pairs = data.frame(trait1 = c(paste0("BBJ_", ts1[t]), paste0("BBJ_", ts1[t]), paste0("BBJ_", ts1[t]), paste0("EUR_", ts1[t]), paste0("EUR_", ts1[t]), paste0("EUR_", ts2[s])),
                           trait2 = c(paste0("BBJ_", ts1[t]), paste0("EUR_", ts1[t]), paste0("EUR_", ts2[s]), paste0("EUR_", ts1[t]), paste0("EUR_", ts2[s]), paste0("EUR_", ts2[s])),
                           pop1 = c("BBJ", "BBJ", "BBJ", "EUR", "EUR", "EUR"), 
                           pop2 = c("BBJ",  "EUR", "EUR", "EUR", "EUR", "EUR") )


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
                            file = paste0("./EUR_bgparas/bgparas_", ts1[t], "_", ts2[s]), quote = F, col.names = F, row.names = F, append = T)

            }

        }, silent = TRUE)

        if (inherits(result, "try-error")) {
            print(paste("Error in round", t, "- skipping to next round."))
            next
        }
    }                           
}
# -

# #### run MRAPSS

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
#colnames(ldscore) = c("SNP", "L2.pop2", "L2.pop1", "L12")
colnames(ldscore) = c("SNP", "L2.pop1", "L2.pop2", "L12")

ldscore0 = ldscore

# +
base_dir <- getwd()

temp_dir <- "./EUR_results"
dir.create(temp_dir, recursive = TRUE, showWarnings = FALSE)


for(s in 1:length(ts2)){

    print(paste('outcome:', ts2[s]))

        for(t in 1:length(ts1)){


        print(paste('exposure:', ts1[t]))

        paras = read.table(paste0(base_dir, "/EUR_bgparas/bgparas_", ts1[t], "_", ts2[s]))
        colnames(paras) = c("trait1", "trait2", "pop2", "pop1", "Omega11", "Omega12","Omega22", "c1","c12","c12")
        paras = paras = subset(paras, trait1 %in% c(paste0("BBJ_", ts1[t]), paste0("EUR_", ts1[t]), paste0("EUR_", ts2[s])) &
                                      trait2 %in% c(paste0("BBJ_", ts1[t]), paste0("EUR_", ts1[t]), paste0("EUR_", ts2[s])))[,  c("trait1", "trait2", "Omega12","c12")]

        paras_opt = paras[, c("trait2", "trait1", "Omega12","c12")]
        colnames(paras_opt) = c("trait1", "trait2", "Omega12","c12")
        paras_opt = subset(paras_opt, trait1 != trait2)
        paras_all = rbind(paras, paras_opt)

        Omega0 = reshape2::acast(unique(paras_all[, c("trait1", "trait2", "Omega12")]), trait1 ~ trait2, value.var = "Omega12")
        C = reshape2::acast(unique(paras_all[, c("trait1", "trait2", "c12")]), trait1 ~ trait2, value.var = "c12")

        Omega0 = Omega0[c(paste0("BBJ_", ts1[t]), paste0("EUR_", ts1[t]), paste0("EUR_", ts2[s])),
                         c(paste0("BBJ_", ts1[t]), paste0("EUR_", ts1[t]), paste0("EUR_", ts2[s]))]
        C = C[c(paste0("BBJ_", ts1[t]), paste0("EUR_", ts1[t]), paste0("EUR_", ts2[s])),
              c(paste0("BBJ_", ts1[t]), paste0("EUR_", ts1[t]), paste0("EUR_", ts2[s]))]

        ### quality control of background parameters
        if((Omega0[1,1] < 0) | (Omega0[2,2] < 0) | (Omega0[3,3] < 0)){
          print('Omega variance negative')
          return(NULL)
        }

        flag2 <- TRUE
        tryCatch(
          {chol(C)},
          error = function(e){
            flag2 <<- FALSE
          })

        if(!flag2){
          print('C not PSD')
          return(NULL)
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


        exp.pop2 = fread(paste0(base_dir, "/formatted_data/EUR_", ts1[t]))
        out.pop2 = fread(paste0(base_dir, "/formatted_data/EUR_", ts2[s]))


        ###
        threshold = 5e-5
        
        
        paras1 = align_data(dat1 = exp.pop2,
                           dat2 = out.pop2,
                           trait1.name = paste0("EUR_", ts1[t]),
                           trait2.name = paste0("EUR_", ts2[s]),
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
                                       bfile = paste0(dirname(base_dir),"/nc/1kg_pops/EUR_phase3_hg38_autosomes"),
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

            write.table(res_mrapss, file="./EUR_results/MRAPSS.MRres", append=T, col.names = F, row.names = F, quote = F)

            TRUE
        }, error = function(e) {
            # Print the actual error message before skipping
            print(paste("Error for exposure", ts1[t], "outcome", ts2[s],
                        "threshold", threshold, ":", e$message))
            FALSE
        })

        if(!result) next

    }
}
# -








