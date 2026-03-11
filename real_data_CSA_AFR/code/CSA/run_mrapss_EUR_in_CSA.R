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

setwd("./real_data_CSA_AFR")

info_csa = read.csv('CSA_data_info.csv', row.names = 1)
data_summary_exp = info_csa[is.na(info_csa$n_controls_CSA),]
data_summary_out = info_csa[!is.na(info_csa$n_controls_CSA),]

ts1 = data_summary_exp$code
ts2 = data_summary_out$code

library(data.table)
library(reshape2)
library(MRAPSS)
library(XMAP) 
library(dplyr)
library(magrittr)
library(Matrix)
library(MASS)
library(TEMR)

pairs <- rbind(
    expand.grid(exposure = ts1, outcome = ts2, stringsAsFactors = FALSE),
    expand.grid(exposure = ts2, outcome = ts1, stringsAsFactors = FALSE),
    subset(expand.grid(exposure = ts2, outcome = ts2, stringsAsFactors = FALSE), exposure != outcome)
)

# ### MRAPSS

# #### estimate background parameters

library(XMAP)
library(data.table)
library(Matrix)

set.seed(1)

# read LD scores
ldscore <- read.csv("../nc/LDscore_SAS_EUR.csv", sep='\t')

# pre-process: remove ambiguous SNPs
#idx_amb <- which(ldscore$allele1 == comple(ldscore$allele2))
#ldscore <- ldscore[-idx_amb,]

ldscore$EUR_CSA = ldscore$CSA_EUR
ldscore$EUR_EUR = ldscore$EUR
ldscore$CSA_CSA = ldscore$CSA
ldscore = as.data.frame(ldscore)

ldscore0 = ldscore

dir.create("./EUR_bgparas_in_CSA/", recursive = TRUE)

for(l in 1:dim(pairs)[1]){

    ts1 = pairs$exposure[l]
    ts2 = pairs$outcome[l]

    s = 1
    t = 1
    
    print(paste0("outcome: ",ts2[s]))
    
    print(paste0("trait: ",ts1[t]))

    bgparas <- estimate_bgparas(
                    ldscore0      = ldscore0,
                    exposure_code = ts1[t], 
                    outcome_code  = ts2[s],      
                    pop_ref       = "CSA",         
                    pop_target    = "EUR",  
                    snp_col       = "SNP",
                    gwas_dir      = "./formatted_data" 
                )

    write.table(bgparas, file = paste0("./EUR_bgparas_in_CSA/bgparas_", ts1[t], "_", ts2[s]), quote = F, col.names = F, row.names = F)
                      
}

# #### run MRAPSS

library(XMAP)
library(data.table)
library(Matrix)

# read LD scores
ldscore <- read.csv("../nc/LDscore_SAS_EUR.csv", sep='\t')

# pre-process: remove ambiguous SNPs
#idx_amb <- which(ldscore$allele1 == comple(ldscore$allele2))
#ldscore <- ldscore[-idx_amb,]

ldscore$EUR_CSA = ldscore$CSA_EUR
ldscore$EUR_EUR = ldscore$EUR
ldscore$CSA_CSA = ldscore$CSA
ldscore = as.data.frame(ldscore)

ldscore = ldscore[, c("SNP", "EUR", "CSA", "EUR_CSA")]
#colnames(ldscore) = c("SNP", "L2.pop1", "L2.pop2", "L12")
colnames(ldscore) = c("SNP", "L2.pop2", "L2.pop1", "L12")

ldscore0 = ldscore

base_dir <- getwd()

temp_dir <- "./EUR_results_in_CSA"
dir.create(temp_dir, recursive = TRUE, showWarnings = FALSE)


for(l in 1:dim(pairs)[1]){

    ts1 = pairs$exposure[l]
    ts2 = pairs$outcome[l]

    s = 1
    t = 1

    print(paste('outcome:', ts2[s]))

    print(paste('exposure:', ts1[t]))

    paras = read.table(paste0(base_dir, "/EUR_bgparas_in_CSA/bgparas_", ts1[t], "_", ts2[s]))
    colnames(paras) = c('trait1','trait2','pop1','pop2','Omega11','Omega12','Omega22','c1','c12','c2')

    bg <- qc_bgparas(
              bgparas       = paras,
              exposure_code = ts1[t],
              outcome_code  = ts2[s],
              pop_ref       = "CSA",
              pop_target    = "EUR",
              truncation_corr = 0.95
              )

    Omega0 <- bg$Omega0
    C      <- bg$C


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

        write.table(res_mrapss, file="./EUR_results_in_CSA/MRAPSS.MRres", append=T, col.names = F, row.names = F, quote = F)

        TRUE
    }, error = function(e) {
        # Print the actual error message before skipping
        print(paste("Error for exposure", ts1[t], "outcome", ts2[s],
                    "threshold", threshold, ":", e$message))
        FALSE
    })

    if(!result) next

}


