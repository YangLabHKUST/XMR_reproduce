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

setwd("./nc")

library(MRAPSS)
library(data.table)
basdir = getwd()
files = dir("raw_data")
variants = fread("full_variant_qc_metrics.txt")
variants = subset(variants, rsid %in% MRAPSS::w_hm3.snplist$SNP)[, c("chrom", "pos", "ref", "alt","rsid","af_AFR", "af_EUR", "af_CSA")]
pheno.info = read.csv("info.csv")

# ### one specific example

# +
for(i in 10:10){
  
  dat_raw <- try(data.table::fread(paste0("./raw_data/", files[i])))
  dat_raw = data.frame(dat_raw)
  dat_raw <- dat_raw[, !(names(dat_raw) %in% c("af_AFR", "af_CSA", "af_EUR"))]
  
  dat_raw = merge(dat_raw, variants,  by.x=c("chr", "pos", "ref", "alt"), by.y=c("chrom", "pos", "ref", "alt"))
  
  filename = sub("-both_sexes-irnt\\.tsv\\.bgz_uncompressed","",files[i])
  filename = sub("-both_sexes-1677\\.tsv\\.bgz_uncompressed","",filename)
  filename = sub("-both_sexes\\.tsv\\.bgz_uncomperssed","",filename)
  filename = sub("biomarkers-","",filename) 
  filename = sub("continuous-","",filename) 
  filename = sub("categorical-","",filename) 
  filename = sub(".tsv.bgz_uncompressed","",filename) 
  filename = ifelse(filename=="1747", "1747_5", filename)
  
  if(inherits(dat_raw,"try-error")) next
  
  dat_EUR = format_data(dat_raw,
                        snp_col = "rsid",
                        b_col = "beta_EUR",
                        se_col = "se_EUR",
                        freq_col = "af_EUR",
                        A1_col = "alt",
                        A2_col = "ref",
                        p_col = "neglog10_pval_EUR",
                        log_pval = T,
                        n =  subset(pheno.info, phenocode==filename)$n_EUR,
                        min_freq = 0.01)
  
  dat_AFR = format_data(dat_raw,
                        snp_col = "rsid",
                        b_col = "beta_AFR",
                        se_col = "se_AFR",
                        freq_col = "af_AFR",
                        A1_col = "alt",
                        A2_col = "ref",
                        p_col = "neglog10_pval_AFR",
                        log_pval = T,
                        n =  subset(pheno.info, phenocode==filename)$n_AFR,
                        min_freq = 0.01)
  
  dat_CSA = format_data(dat_raw,
                        snp_col = "rsid",
                        b_col = "beta_CSA",
                        se_col = "se_CSA",
                        freq_col = "af_CSA",
                        A1_col = "alt",
                        A2_col = "ref",
                        p_col = "neglog10_pval_CSA",
                        log_pval = T,
                        n =  subset(pheno.info, phenocode==filename)$n_CSA,
                        min_freq = 0.01)
  
  
}
# -

# ### code to formate all traits

# +
dir.create("./formatted_data", recursive = TRUE)

for(i in 1:37){
  
  dat_raw <- try(data.table::fread(paste0("./raw_data/", files[i])))
  dat_raw = data.frame(dat_raw)
  dat_raw <- dat_raw[, !(names(dat_raw) %in% c("af_AFR", "af_CSA", "af_EUR"))]
  
  dat_raw = merge(dat_raw, variants,  by.x=c("chr", "pos", "ref", "alt"), by.y=c("chrom", "pos", "ref", "alt"))
  
  filename = sub("-both_sexes-irnt\\.tsv\\.bgz_uncompressed","",files[i])
  filename = sub("-both_sexes-1677\\.tsv\\.bgz_uncompressed","",filename)
  filename = sub("-both_sexes\\.tsv\\.bgz_uncomperssed","",filename)
  filename = sub("biomarkers-","",filename) 
  filename = sub("continuous-","",filename) 
  filename = sub("categorical-","",filename) 
  filename = sub(".tsv.bgz_uncompressed","",filename) 
  filename = ifelse(filename=="1747", "1747_5", filename)
  
  if(inherits(dat_raw,"try-error")) next
  
  dat_EUR = format_data(dat_raw,
                        snp_col = "rsid",
                        b_col = "beta_EUR",
                        se_col = "se_EUR",
                        freq_col = "af_EUR",
                        A1_col = "alt",
                        A2_col = "ref",
                        p_col = "neglog10_pval_EUR",
                        log_pval = T,
                        n =  subset(pheno.info, phenocode==filename)$n_EUR,
                        min_freq = 0.01)
  
  dat_AFR = format_data(dat_raw,
                        snp_col = "rsid",
                        b_col = "beta_AFR",
                        se_col = "se_AFR",
                        freq_col = "af_AFR",
                        A1_col = "alt",
                        A2_col = "ref",
                        p_col = "neglog10_pval_AFR",
                        log_pval = T,
                        n =  subset(pheno.info, phenocode==filename)$n_AFR,
                        min_freq = 0.01)
  
  dat_CSA = format_data(dat_raw,
                        snp_col = "rsid",
                        b_col = "beta_CSA",
                        se_col = "se_CSA",
                        freq_col = "af_CSA",
                        A1_col = "alt",
                        A2_col = "ref",
                        p_col = "neglog10_pval_CSA",
                        log_pval = T,
                        n =  subset(pheno.info, phenocode==filename)$n_CSA,
                        min_freq = 0.01)
    
    write.table(dat_AFR, 
              paste0("./formatted_data/AFR_", filename),
              sep="\t", quote = F, row.names = F, col.names = T)
    write.table(dat_CSA, 
              paste0("./formatted_data/CSA_", filename),
              sep="\t", quote = F, row.names = F, col.names = T)
    write.table(dat_EUR, 
              paste0("./formatted_data/EUR_", filename),
              sep="\t", quote = F, row.names = F, col.names = T)

}
# -


