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
library(TEMR)
library(cause)
library(MRAPSS)
library(TwoSampleMR)
library(mr.divw)
library(MRcML)
library(MRMix)
library(MRPRESSO)
library(robustbase)

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

# ### TEMR

# +
base_dir <- getwd()

temp_dir <- "./BBJ_results"
dir.create(temp_dir, recursive = TRUE, showWarnings = FALSE)


for(threshold in c(5e-8)){
    
    for(s in 1:length(ts2)){

        print(paste('outcome:', ts2[s]))

        for(t in 1:length(ts1)){

            print(paste('exposure:', ts1[t]))

            exp_dat1 <- fread(paste0(base_dir, "/formatted_data/BBJ_", ts1[t]))
            exp_dat2 <- fread(paste0(base_dir, "/formatted_data/EUR_", ts1[t]))
            out_dat1 <- fread(paste0(base_dir, "/formatted_data/BBJ_", ts2[s]))
            out_dat2 <- fread(paste0(base_dir, "/formatted_data/EUR_", ts2[s]))

            # pre-process: remove ambiguous SNPs
            exp_dat1 <- exp_dat1[exp_dat1$A1 != comple(exp_dat1$A2),]
            exp_dat2 <- exp_dat2[exp_dat2$A1 != comple(exp_dat2$A2),]
            out_dat1 <- out_dat1[out_dat1$A1 != comple(out_dat1$A2),]
            out_dat2 <- out_dat2[out_dat2$A1 != comple(out_dat2$A2),]

            # pre-process: overlap SNPs
            snps <- Reduce(intersect, list(exp_dat1$SNP, exp_dat2$SNP, out_dat1$SNP, out_dat2$SNP))
            exp_dat1 <- exp_dat1[match(snps, exp_dat1$SNP),]
            exp_dat2 <- exp_dat2[match(snps, exp_dat2$SNP),]
            out_dat1 <- out_dat1[match(snps, out_dat1$SNP),]
            out_dat2 <- out_dat2[match(snps, out_dat2$SNP),]

            paras1 = align_data(dat1 = exp_dat1,
                               dat2 = exp_dat2,
                               trait1.name = paste0("BBJ_", ts1[t]),
                               trait2.name = paste0("EUR_", ts1[t]),
                               h2.fix.intercept = F,
                               LDSC = F
                               )


            paras2 = align_data(dat1 = out_dat1,
                               dat2 = out_dat2,
                               trait1.name = paste0("BBJ_", ts2[s]),
                               trait2.name = paste0("EUR_", ts2[s]),
                               h2.fix.intercept = F,
                               LDSC = F
                               )


            exp_dat = paras1$dat
            colnames(exp_dat) = c("SNP", "A1", "A2", "b.exp.pop1", "b.exp.pop2", "se.exp.pop1", "se.exp.pop2", "pval.exp.pop1", "pval.exp.pop2")

            out_dat = paras2$dat
            colnames(out_dat) = c("SNP", "A1", "A2", "b.out.pop1", "b.out.pop2", "se.out.pop1", "se.out.pop2", "pval.out.pop1", "pval.out.pop2")

            dat = merge(exp_dat, out_dat, by=c("SNP"))

            ### pre-process: flip alleles; z-score (beta/se)
            flip.index = which((dat$A1.x == dat$A2.y & dat$A1.y == dat$A2.x) |
                               (dat$A1.x ==comple(dat$A2.y) & dat$A1.y == comple(dat$A2.x)))
            dat[,"A1.y"] = dat[,"A1.x"]
            dat[,"A2.y"] = dat[,"A2.x"]
            dat$A1 = dat[,"A1.x"]
            dat$A2 = dat[,"A2.x"]
            dat[flip.index ,"b.out.pop1"] = -dat[flip.index ,"b.out.pop1"]
            dat[flip.index ,"b.out.pop2"] = -dat[flip.index ,"b.out.pop2"]

            print(paste0('dim of dat: ', dim(dat)[1]))
            print(paste0('length of flip.index: ', length(flip.index)))


            # Choose SNPs that are significant in at least 1 population
            dat <- dat[(dat$pval.exp.pop1 < threshold) | (dat$pval.exp.pop2 < threshold),]

            temp = qnorm(1-dat$pval.exp.pop1/2) + qnorm(1-dat$pval.exp.pop2/2)
            dat$P.combined = 2*(1 - pnorm(temp))

            dat = merge(dat, ldscore0, by="SNP")

            result <- tryCatch({
                clumped = MRAPSS::clump(dat,
                                         IV.Threshold = 1,
                                         SNP_col = "SNP",
                                         pval_col = "P.combined",
                                         clump_kb = 1000,
                                         clump_r2 = 0.001,
                                         bfile = paste0(dirname(base_dir),"/nc/1kg_pops/EAS_phase3_hg38_autosomes"),
                                         plink_bin = paste0(dirname(base_dir),"/nc/plink_Linux"))
                clumped$Threshold = threshold
                                        
                
                print('Run TEMR')
                ############# run TEMR ################
                betaXG = clumped[,c('b.exp.pop1', 'b.exp.pop2')]
                betaYG = clumped[,c('b.out.pop1', 'b.out.pop2')]
                sebetaXG = clumped[,c('se.exp.pop1', 'se.exp.pop2')]
                sebetaYG = clumped[,c('se.out.pop1', 'se.out.pop2')] 

                race_name <- c('target1','auxiliary')
                colnames(betaXG) <- race_name
                colnames(betaYG) <- race_name
                colnames(sebetaXG) <- race_name
                colnames(sebetaYG) <- race_name

                rho <- NULL
                meth <- 'CG'
                RESULT <- TEMR(betaXG,betaYG,sebetaXG,sebetaYG,rho,meth)

                RESULT$resultall

                res_temr = data.frame(exposure = ts1[t],
                                      outcome = ts2[s],
                                      beta = RESULT$resultall$target1$beta[['bbe1']],
                                      beta.pvalue = RESULT$resultall$target1$pvalue,
                                      stat = RESULT$resultall$target1$stat,
                                      nIV = dim(clumped)[1],
                                      Threshold = threshold)

                write.table(res_temr, file="./BBJ_results/TEMR.MRres", append=T, col.names = F, row.names = F, quote = F)

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
}
# -

# ### MRAPSS

# +
base_dir <- getwd()

temp_dir <- "./BBJ_results"
dir.create(temp_dir, recursive = TRUE, showWarnings = FALSE)


for(s in 1:length(ts2)){

    print(paste('outcome:', ts2[s]))

    for(t in 1:length(ts1)){

        print(paste('exposure:', ts1[t]))

        paras = read.table(paste0(base_dir, "/BBJ_bgparas/bgparas_", ts1[t], "_", ts2[s]))
        colnames(paras) = c("trait1", "trait2", "pop2", "pop1", "Omega11", "Omega12","Omega22", "c1","c12","c12")
        paras = subset(paras, trait1 %in% c(paste0("EUR_", ts1[t]), paste0("BBJ_", ts1[t]), paste0("BBJ_", ts2[s])) &
                             trait2 %in% c(paste0("EUR_", ts1[t]), paste0("BBJ_", ts1[t]), paste0("BBJ_", ts2[s])))[, c("trait1", "trait2", "Omega12","c12")]

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


        exp.pop2 = fread(paste0(base_dir, "/formatted_data/BBJ_", ts1[t]))
        out.pop2 = fread(paste0(base_dir, "/formatted_data/BBJ_", ts2[s]))


        ###
        threshold = 5e-5
        
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

            write.table(res_mrapss, file="./BBJ_results/MRAPSS.MRres", append=T, col.names = F, row.names = F, quote = F)

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

# ### TwosampleMR (IVW, Egger, Weighted-mode, Weighted-median)

# +
base_dir <- getwd()

temp_dir <- "./BBJ_results"
dir.create(temp_dir, recursive = TRUE, showWarnings = FALSE)

    
for(s in 1:length(ts2)){

    print(paste('outcome:', ts2[s]))

    for(t in 1:length(ts1)){

        print(paste('exposure:', ts1[t]))


        exp.pop2 = fread(paste0(base_dir, "/formatted_data/BBJ_", ts1[t]))
        out.pop2 = fread(paste0(base_dir, "/formatted_data/BBJ_", ts2[s]))

        ###
        N2 = median(exp.pop2$N)
        if(N2 > 20000){
            threshold = 5e-8
        }else{
            threshold = 5e-5
        }

        paras1 = align_data(dat1 = exp.pop2,
                           dat2 = out.pop2,
                           trait1.name = paste0("BBJ_", ts1[t]),
                           trait2.name = paste0("BBJ_", ts2[s]),
                           h2.fix.intercept = F,
                           LDSC = F
                           )


        dat = paras1$dat
        colnames(dat) = c("SNP", "A1", "A2", "b.exp", "b.out", "se.exp", "se.out", "pval.exp", "pval.out" )

        dat = merge(dat, ldscore0, by="SNP")
        dat$Threshold = threshold


        result <- tryCatch({
            clumped.MR = MRAPSS::clump(dat,
                                       IV.Threshold = threshold,
                                       SNP_col = "SNP",
                                       pval_col = "pval.exp",
                                       clump_kb = 1000,
                                       clump_r2 = 0.001,
                                       bfile = paste0(dirname(base_dir),"/nc/1kg_pops/EAS_phase3_hg38_autosomes"),
                                       plink_bin = paste0(dirname(base_dir),"/nc/plink_Linux"))
            clumped.MR$Threshold = threshold


            betaXG = clumped.MR[, 'b.exp']
            betaYG = clumped.MR[, 'b.out']
            sebetaXG = clumped.MR[, 'se.exp']
            sebetaYG = clumped.MR[, 'se.out']


            IVW_result = TwoSampleMR::mr_ivw(b_exp = betaXG, b_out = betaYG, 
                                se_exp = sebetaXG, se_out =sebetaYG)

            egger_result = TwoSampleMR::mr_egger_regression(b_exp = betaXG, b_out = betaYG, 
                                               se_exp = sebetaXG, se_out =sebetaYG)

            weightmed_result = TwoSampleMR::mr_weighted_median(b_exp = betaXG, b_out = betaYG, 
                                                  se_exp = sebetaXG, se_out =sebetaYG)

            weightmod_result = TwoSampleMR::mr_weighted_mode(b_exp = betaXG, b_out = betaYG, 
                                                se_exp = sebetaXG, se_out =sebetaYG)


            write.table(data.frame(exposure = ts1[t],
                                   outcome = ts2[s],
                                   beta = IVW_result$b,
                                   beta.se = IVW_result$se,
                                   beta.pvalue = IVW_result$pval,
                                   nIV = dim(clumped.MR)[1],
                                   Threshold = threshold), 
                        file="./BBJ_results/IVW.MRres", append=T, col.names = F, row.names = F, quote = F)

            write.table(data.frame(exposure = ts1[t],
                                   outcome = ts2[s],
                                   beta = egger_result$b,
                                   beta.se = egger_result$se,
                                   beta.pvalue = egger_result$pval,
                                   nIV = dim(clumped.MR)[1],
                                   Threshold = threshold), 
                        file="./BBJ_results/Egger.MRres", append=T, col.names = F, row.names = F, quote = F)

            write.table(data.frame(exposure = ts1[t],
                                   outcome = ts2[s],
                                   beta = weightmed_result$b,
                                   beta.se = weightmed_result$se,
                                   beta.pvalue = weightmed_result$pval,
                                   nIV = dim(clumped.MR)[1],
                                   Threshold = threshold), 
                        file="./BBJ_results/Weightmed.MRres", append=T, col.names = F, row.names = F, quote = F)

            write.table(data.frame(exposure = ts1[t],
                                   outcome = ts2[s],
                                   beta = weightmod_result$b,
                                   beta.se = weightmod_result$se,
                                   beta.pvalue = weightmod_result$pval,
                                   nIV = dim(clumped.MR)[1],
                                   Threshold = threshold), 
                        file="./BBJ_results/Weightmode.MRres", append=T, col.names = F, row.names = F, quote = F)


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

# ### IVW_fe, MRMix, MR-Robust, MR-Lasso, cML-MA, MR-PRESSO, ConMix, dIVW

# +
base_dir <- getwd()

temp_dir <- "./BBJ_results"
dir.create(temp_dir, recursive = TRUE, showWarnings = FALSE)


methods = c("IVW_fe", "MRMix", "MR-Robust", "MR-Lasso", "cML-MA", "MR-PRESSO")

    
for(s in 1:length(ts2)){

    print(paste('outcome:', ts2[s]))

    for(t in 1:length(ts1)){

        print(paste('exposure:', ts1[t]))


        exp.pop2 = fread(paste0(base_dir, "/formatted_data/BBJ_", ts1[t]))
        out.pop2 = fread(paste0(base_dir, "/formatted_data/BBJ_", ts2[s]))


        ###
        N2 = median(exp.pop2$N)
        if(N2 > 20000){
            threshold = 5e-8
        }else{
            threshold = 5e-5
        }

        paras1 = align_data(dat1 = exp.pop2,
                           dat2 = out.pop2,
                           trait1.name = paste0("BBJ_", ts1[t]),
                           trait2.name = paste0("BBJ_", ts2[s]),
                           h2.fix.intercept = F,
                           LDSC = F
                           )


        dat = paras1$dat
        colnames(dat) = c("SNP", "A1", "A2", "b.exp", "b.out", "se.exp", "se.out", "pval.exp", "pval.out" )

        dat = merge(dat, ldscore0, by="SNP")
        dat$Threshold = threshold


        result <- tryCatch({
            clumped.MR = MRAPSS::clump(dat,
                                       IV.Threshold = threshold,
                                       SNP_col = "SNP",
                                       pval_col = "pval.exp",
                                       clump_kb = 1000,
                                       clump_r2 = 0.001,
                                       bfile = paste0(dirname(base_dir),"/nc/1kg_pops/EAS_phase3_hg38_autosomes"),
                                       plink_bin = paste0(dirname(base_dir),"/nc/plink_Linux"))
            clumped.MR$Threshold = threshold


            ldmat = NULL

            for(i in 1:length(methods)){

              ### only run when nIV > 3
              res = runMRmethods(dat=clumped.MR, ldmat=ldmat, Threshold=threshold, 
                               exposure=ts1[t], outcome=ts2[s], methods.list = methods[i])

              write.table(res, paste0("./BBJ_results/",methods[i],".MRres"), quote=F, col.names = F, append = T,row.names = F)

            }


            res_conmix = runMRmethods(dat=clumped.MR, ldmat=ldmat, Threshold=threshold, 
                               exposure=ts1[t], outcome=ts2[s], ConMixonly = T)

            write.table(res_conmix, paste0("./BBJ_results/","ConMix",".MRres"), quote=F, col.names = F, append = T,row.names = F)


            res_divw = mr.divw::mr.divw(beta.exposure=clumped.MR$b.exp, beta.outcome=clumped.MR$b.out, 
                                        se.exposure=clumped.MR$se.exp, se.outcome=clumped.MR$se.out, 
                                        diagnostics=F, over.dispersion=T)

            divw.res = data.frame(exposure = ts1[t], outcome=ts2[s],
                                    method = "dIVW", Threshold = threshold,
                                    nsnp = res_divw$n.IV , beta = res_divw$beta.hat, se = res_divw$beta.se, 
                                    pval = pchisq(res_divw$beta.hat^2/res_divw$beta.se^2, 1, lower.tail = F))

            write.table(divw.res, paste0("./BBJ_results/","dIVW",".MRres"), quote=F, col.names = F, append = T,row.names = F)


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

# ### CAUSE

# +
library(readr)
library(cause)

base_dir <- getwd()

temp_dir <- "./BBJ_results"
dir.create(temp_dir, recursive = TRUE, showWarnings = FALSE)

for(threshold in c(5e-5)){
    
    for(s in 1:length(ts2)){

        print(paste('outcome:', ts2[s]))

        for(t in 1:1){

            print(paste('exposure:', ts1[t]))

            X1 = fread(paste0(base_dir, "/formatted_data/BBJ_", ts1[t]))
            X2 = fread(paste0(base_dir, "/formatted_data/BBJ_", ts2[s]))

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

            if(inherits(X , 'try-error')) next

            d0 = X1[, c("SNP", "P")]
            colnames(d0) = c("snp", "pval.exp")
            X0 = merge(X, d0, by="snp")

            #############clumping#######################    
            # pre-process: remove ambiguous SNPs
            exp_dat1 <- X1[X1$A1 != comple(X1$A2),]
            out_dat1 <- X2[X2$A1 != comple(X2$A2),]

            paras1 = align_data(dat1 = exp_dat1,
                               dat2 = out_dat1,
                               trait1.name = paste0("BBJ_", ts1[t]),
                               trait2.name = paste0("BBJ_", ts2[s]),
                               h2.fix.intercept = F,
                               LDSC = F
                           )

            dat = paras1$dat
            colnames(dat) = c("SNP", "A1", "A2", "b.exp", "b.out", "se.exp", "se.out", "pval.exp", "pval.out")

            dat = merge(dat, ldscore0, by="SNP")

            result <- tryCatch({
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
                if(inherits(params , 'try-error')) next
                if(inherits(clumped0 , 'try-error')) next


                  clumped = subset(clumped0, pval.exp <= threshold)
                
                  # cause
                  if(!is.null(clumped)){
                    top_ldl_pruned_vars =intersect(as.character(X$snp), as.character(clumped$snp))

                    cause_res <- try(cause(X=X, variants = top_ldl_pruned_vars , param_ests = params, force=TRUE))

                    if(inherits( cause_res , 'try-error')) next

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
                    write.table(cause_res, file=paste0("./BBJ_results/CAUSE.MRres"), append=T, col.names = F, row.names = F, quote = F)

                    rm(top_ldl_pruned_vars)
                    rm(res_est)
                    rm(res_elpd)
                    rm(res.cause.est)
                    rm(cause_res)
                      
                }

                
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
}
# -














