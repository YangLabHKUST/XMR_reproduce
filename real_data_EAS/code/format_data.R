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

setwd("./real_data_EAS")

# +
library(MRAPSS)
library(data.table)
library(MASS)
library(ggpubr)

basdir = getwd()
# -

# ### format TPMI data

pheno.info = read.csv("./raw_data/TPMI/summary.csv")
pheno.info$EUR.phecode = as.character(pheno.info$EUR.phecode)

# +
# Read in dbSNP153 hg38 data
dbSNP38 <- readr::read_delim("./raw_data/TPMI/dbSnp153Common.bed", col_names = FALSE)[, c(1:4)]
colnames(dbSNP38) <- c("CHR", "START", "END", "rsID")

setDT(dbSNP38)

setkey(dbSNP38, CHR, END)
# -

head(pheno.info)

# +
for(i in 1:dim(pheno.info)[1]){
    
    phecode =  pheno.info[['EAS.phecode']][i]
    
    print(paste0("Trait: ", phecode, "..."))
    
    tw_address = paste0('./raw_data/TPMI/', pheno.info[['EAS.address']][i])
    data_tw = fread(tw_address)
    
    ### assign rsID
    # some with the suffix 'chr'; some not
    if(nchar(data_tw$CHR[1]) < 3){
        data_tw$CHR_modified = paste0("chr", data_tw$CHR)
    }else{
        data_tw$CHR_modified = data_tw$CHR
    }
    
    # match with keys; assign 'rsID' to data_tw in column 'SNPID'
    data_tw$SNPID = data_tw$MarkerID
    data_tw[dbSNP38, on = .(CHR_modified = CHR, POS = END), SNPID := i.rsID]

    data_tw[is.na(data_tw$SNPID)]$SNPID = data_tw[is.na(data_tw$SNPID)]$MarkerID
    data_tw[(data_tw$SNPID=='NA')]$SNPID = data_tw[(data_tw$SNPID=='NA')]$MarkerID
    
    
    result <- try({
        
        dat_tw = MRAPSS::format_data(data_tw,
                                      snp_col = "SNPID",
                                      b_col = "BETA",
                                      se_col = "SE",
                                      A1_col = "Allele2", ### effect allele
                                      A2_col = "Allele1", ### non-effect allele
                                      p_col = "p.value",
                                      freq_col = "AF_Allele2",
                                      n = pheno.info[['EAS.samples']][i],
                                      min_freq = 0.01)

        ########################### name rule: 'TW_' + its EUR.phecode ############################
        write.table(dat_tw, file=paste0("./formatted_data/TW_", pheno.info[['EUR.phecode']][i]), quote = F, col.names = T, row.names = F)

    }, silent = TRUE)

    if (inherits(result, "try-error")) {
        print(paste("Error in round", i, "- skipping to next round."))
        next
    }
    
}
# -



# ### format BBJ data

# #### exposure traits

file_path <- './raw_data/BBJ/EUR_BBJ_exposure_description.csv'
data_summary_exp <- read.csv(file_path)

data_summary_exp

# +
for(i in 1:dim(data_summary_exp)[1]){
    
    phecode =  data_summary_exp[['BBJ.code']][i]
    
    print(paste0("Trait: ", phecode, "..."))
    
    BBJ_address = paste0("./raw_data/BBJ/hum0197.v3.BBJ.", phecode, ".v1/GWASsummary_", phecode, "_Japanese_SakaueKanai2020.auto.txt.gz")
    
    data_BBJ = fread(BBJ_address)
    
    result <- try({
        
        dat_BBJ = format_data(data_BBJ,
                              snp_col = "SNP",
                              b_col = "BETA",
                              se_col = "SE",
                              A1_col = "ALLELE1", ### effect allele
                              A2_col = "ALLELE0", ### non-effect allele
                              p_col = "P_BOLT_LMM_INF",
                              freq_col = "A1FREQ",
                              n = data_summary_exp[['BBJ.samples']][i],
                              min_freq = 0.01)

        ########################### name rule: 'BBJ_' + its EUR.phecode ############################
        write.table(dat_BBJ, file=paste0("./formatted_data/BBJ_", data_summary_exp[['UKB.code']][i]), quote = F, col.names = T, row.names = F)

    }, silent = TRUE)

    if (inherits(result, "try-error")) {
        print(paste("Error in round", i, "- skipping to next round."))
        next
    }
    
}
# -



# #### outcome traits

data_summary_out = read.csv('./raw_data/BBJ/EUR_BBJ_outcome_description.csv')
data_summary_out$UKB.phecode = as.character(data_summary_out$UKB.phecode)
data_summary_out$phecode_unique = as.character(data_summary_out$phecode_unique)

head(data_summary_out)

# +
for(i in 1:dim(data_summary_out)[1]){
    phecode =  data_summary_out[['phecode_unique']][i]
    
    print(paste0("Trait: ", phecode, "..."))
    
    BBJ_address = paste0("./raw_data/BBJ/", data_summary_out[['BBJ.address']][i])

    data_BBJ = data.table::fread(BBJ_address)
    
    result <- try({
    
        dat_BBJ = format_data(data_BBJ,
                              snp_col = "SNPID",
                              b_col = "BETA",
                              se_col = "SE",
                              A1_col = "Allele2", ### effect allele
                              A2_col = "Allele1", ### non-effect allele
                              p_col = "p.value",
                              freq_col = "AF_Allele2",
                              n = data_summary_out[['BBJ.samples']][i],
                              min_freq = 0.01)
        
        ########################### name rule: 'BBJ_' + its EUR.phecode ############################
        write.table(dat_BBJ, file=paste0("./formatted_data/BBJ_", data_summary_out[['UKB.phecode']][i]), quote = F, col.names = T, row.names = F)

        
    }, silent = TRUE)

    if (inherits(result, "try-error")) {
        print(paste("Error in round", i, "- skipping to next round."))
        next
    }
    
}
# -



# ### format EUR data

# #### exposure traits

variants = fread('./raw_data/EUR/variants.tsv')

file_path <- './raw_data/BBJ/EUR_BBJ_exposure_description.csv'
data_summary_exp <- read.csv(file_path)
data_summary_exp

# +
for(i in 1:dim(data_summary_exp)[1]){
    
    phecode = data_summary_exp[['UKB.code']][i]
    
    print(paste0("Trait: ", phecode, "..."))
    
    UKB_address = paste0("./raw_data/EUR/", phecode, ".tsv")
    
    data_UKB = fread(UKB_address)
    
    check1 = sum(data_UKB$variant != variants$variant)
    
    if(check1 != 0){
        print('Need more careful alignment with variants.tsv... skip')
        next
    }
    
    data_UKB$SNP = variants$rsid
    
    data_UKB$A1 = variants$alt # effect allele
    data_UKB$A2 = variants$ref # non-effect allele
    
    result <- try({
    
        dat_UKB = format_data(data_UKB,
                      snp_col = "SNP",
                      b_col = "beta",
                      se_col = "se",
                      A1_col = "A1",
                      A2_col = "A2",
                      p_col = "pval",
                      freq_col = "minor_AF",
                      #n = data_summary_exp[['UKB.samples']][i],
                      n_col = 'n_complete_samples',
                      min_freq = 0.01)
        
        write.table(dat_UKB, file=paste0("./formatted_data/EUR_",phecode), quote = F, col.names = T, row.names = F)
        
    }, silent = TRUE)

    if (inherits(result, "try-error")) {
        print(paste("Error in round", i, "- skipping to next round."))
        next
    }
    
}
# -



# #### outcome traits

data_summary_out = read.csv('./raw_data/BBJ/EUR_BBJ_outcome_description.csv')
data_summary_out$UKB.phecode = as.character(data_summary_out$UKB.phecode)
data_summary_out$phecode_unique = as.character(data_summary_out$phecode_unique)

data_summary_out

# +
# Read in dbSNP153 data
dbSNP <- readr::read_delim("./raw_data/EUR/dbSnp153Common.hg19.bed", col_names = FALSE)[, c(1:4)]
colnames(dbSNP) <- c("CHR", "START", "END", "rsID")

setDT(dbSNP)

setkey(dbSNP, CHR, END)

# +
for(i in 1:dim(data_summary_out)[1]){
    phecode =  data_summary_out[['phecode_unique']][i]
    
    print(paste0("Trait: ", phecode, "..."))
    
    data_UKB = data.table::fread(paste0('./raw_data/EUR/', data_summary_out[['UKB.address']][i]))
    
    result <- try({
        
        data_UKB$CHR_modified = paste0("chr", data_UKB$CHR)
        data_UKB$SNPID = data_UKB$v

        # match with keys; assign 'rsID' to data_tw in column 'SNPID'
        data_UKB[dbSNP, on = .(CHR_modified = CHR, POS = END), SNPID := i.rsID]

        data_UKB[is.na(data_UKB$SNPID)]$SNPID = data_UKB[is.na(data_UKB$SNPID)]$v
        data_UKB[(data_UKB$SNPID=='NA')]$SNPID = data_UKB[(data_UKB$SNPID=='NA')]$v

        data_UKB$AF_Allele2 = data_UKB$AF_Allele2_UKB
        data_UKB[is.na(data_UKB$AF_Allele2)]$AF_Allele2 = data_UKB[is.na(data_UKB$AF_Allele2)]$AF_Allele2_FG
        
        data_UKB = format_data(data_UKB,
                              snp_col = "SNPID",
                              b_col = "BETA",
                              se_col = "SE",
                              A1_col = "Allele2", # effect allele
                              A2_col = "Allele1", # effect allele
                              p_col = "p.value",
                              freq_col = "AF_Allele2",
                              n = data_summary_out[['UKB.samples']][i],
                              min_freq = 0.01)

        write.table(data_UKB, file=paste0("./formatted_data/EUR_",phecode), quote = F, col.names = T, row.names = F)

    }, silent = TRUE)

    if (inherits(result, "try-error")) {
        print(paste("Error in round", i, "- skipping to next round."))
        next
    }
    
}
# -


