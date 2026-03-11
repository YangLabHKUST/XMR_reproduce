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

library(MRAPSS)
library(data.table)
library(MASS)
library(ggpubr)

basdir = getwd()

# ### format CSA/EUR data

info_csa = read.csv('CSA_data_info.csv', row.names = 1)

head(info_csa)

variants = data.table::fread('./full_variant_qc_metrics.txt')

for(i in 1:dim(info_csa)[1]){
    
    phecode =  info_csa[['code']][i]
    
    print(paste0("Trait: ", phecode, "..."))
    
    data = data.table::fread(paste0('./raw_data/', info_csa[['address']][i]))
    
    if(sum(!(data$chr == variants$chrom & data$pos == variants$pos)) != 0){
        print('Data and variants not fully matched. Skip...')
        next
    }
    
    data$rsid = variants$rsid
    
    if(is.na(info_csa[['n_controls_CSA']][i])){
        
        data = data[, c('chr','pos','rsid','ref','alt','af_CSA','af_EUR','beta_CSA','beta_EUR','se_CSA','se_EUR','neglog10_pval_CSA','neglog10_pval_EUR','low_confidence_CSA','low_confidence_EUR')]
        data_sub = data[!(is.na(data$low_confidence_CSA) | is.na(data$low_confidence_EUR)),]

        data_sub = data_sub[(data_sub$low_confidence_CSA == FALSE) & (data_sub$low_confidence_EUR == FALSE),]

        data_sub$pval_CSA = 10^(-data_sub$neglog10_pval_CSA)
        data_sub$pval_EUR = 10^(-data_sub$neglog10_pval_EUR)

        data_sub = as.data.frame(data_sub)

        dat_EUR = format_data(data_sub,
                              snp_col = "rsid",
                              b_col = "beta_EUR",
                              se_col = "se_EUR",
                              A1_col = "alt",
                              A2_col = "ref",
                              p_col = "pval_EUR",
                              freq_col = "af_EUR",
                              n = info_csa[['n_EUR']][i],
                              min_freq = 0.01)

        dat_CSA = format_data(data_sub,
                              snp_col = "rsid",
                              b_col = "beta_CSA",
                              se_col = "se_CSA",
                              A1_col = "alt",
                              A2_col = "ref",
                              p_col = "pval_CSA",
                              freq_col = "af_CSA",
                              n = info_csa[['n_CSA']][i],
                              min_freq = 0.01)

    }else{
        
        data = data[, c('chr','pos','rsid','ref','alt','af_cases_CSA','af_cases_EUR','beta_CSA','beta_EUR','se_CSA','se_EUR','neglog10_pval_CSA','neglog10_pval_EUR','low_confidence_CSA','low_confidence_EUR')]
        data_sub = data[!(is.na(data$low_confidence_CSA) | is.na(data$low_confidence_EUR)),]

        data_sub = data_sub[(data_sub$low_confidence_CSA == FALSE) & (data_sub$low_confidence_EUR == FALSE),]

        data_sub$pval_CSA = 10^(-data_sub$neglog10_pval_CSA)
        data_sub$pval_EUR = 10^(-data_sub$neglog10_pval_EUR)

        data_sub = as.data.frame(data_sub)

        dat_EUR = format_data(data_sub,
                              snp_col = "rsid",
                              b_col = "beta_EUR",
                              se_col = "se_EUR",
                              A1_col = "alt",
                              A2_col = "ref",
                              p_col = "pval_EUR",
                              freq_col = "af_cases_EUR",
                              n = info_csa[['n_EUR']][i],
                              min_freq = 0.01)

        dat_CSA = format_data(data_sub,
                              snp_col = "rsid",
                              b_col = "beta_CSA",
                              se_col = "se_CSA",
                              A1_col = "alt",
                              A2_col = "ref",
                              p_col = "pval_CSA",
                              freq_col = "af_cases_CSA",
                              n = info_csa[['n_CSA']][i],
                              min_freq = 0.01)

    }
                  
    write.table(dat_EUR, file=paste0("./formatted_data/EUR_",phecode), quote = F, col.names = T, row.names = F)
    write.table(dat_CSA, file=paste0("./formatted_data/CSA_",phecode), quote = F, col.names = T, row.names = F)

}



# ### format AFR data

info_afr = read.csv('AFR_data_info.csv', row.names = 1)

info_afr$address =sub("\\.bgz$", "", basename(info_afr$AFR.data.link))
info_afr[dim(info_afr)[1], 'address'] = 'AFR_alb_new.txt'

info_afr

for(i in 1:(dim(info_afr)[1]-1)){
    
    phecode =  info_afr[['EUR.phecode']][i]
    
    print(paste0("Trait: ", phecode, "..."))
    
    data = data.table::fread(paste0('./raw_data/', info_afr[['address']][i]))
    
    if(sum(!(data$chr == variants$chrom & data$pos == variants$pos)) != 0){
        print('Data and variants not fully matched. Skip...')
        next
    }
    
    data$rsid = variants$rsid
    
    if(is.na(info_afr[['AFR.controls']][i])){
        
        data = data[, c('chr','pos','rsid','ref','alt','af_AFR','af_EUR','beta_AFR','beta_EUR','se_AFR','se_EUR','neglog10_pval_AFR','neglog10_pval_EUR','low_confidence_AFR','low_confidence_EUR')]
        data_sub = data[!(is.na(data$low_confidence_AFR) | is.na(data$low_confidence_EUR)),]

        data_sub = data_sub[(data_sub$low_confidence_AFR == FALSE) & (data_sub$low_confidence_EUR == FALSE),]

        data_sub$pval_AFR = 10^(-data_sub$neglog10_pval_AFR)
        data_sub$pval_EUR = 10^(-data_sub$neglog10_pval_EUR)

        data_sub = as.data.frame(data_sub)

        dat_AFR = format_data(data_sub,
                              snp_col = "rsid",
                              b_col = "beta_AFR",
                              se_col = "se_AFR",
                              A1_col = "alt",
                              A2_col = "ref",
                              p_col = "pval_AFR",
                              freq_col = "af_AFR",
                              n = info_afr[['AFR.samples']][i],
                              min_freq = 0.01)

    }else{
        
        data = data[, c('chr','pos','rsid','ref','alt','af_cases_AFR','af_cases_EUR','beta_AFR','beta_EUR','se_AFR','se_EUR','neglog10_pval_AFR','neglog10_pval_EUR','low_confidence_AFR','low_confidence_EUR')]
        data_sub = data[!(is.na(data$low_confidence_AFR) | is.na(data$low_confidence_EUR)),]

        data_sub = data_sub[(data_sub$low_confidence_AFR == FALSE) & (data_sub$low_confidence_EUR == FALSE),]

        data_sub$pval_AFR = 10^(-data_sub$neglog10_pval_AFR)

        data_sub = as.data.frame(data_sub)
        
        dat_AFR = format_data(data_sub,
                              snp_col = "rsid",
                              b_col = "beta_AFR",
                              se_col = "se_AFR",
                              A1_col = "alt",
                              A2_col = "ref",
                              p_col = "pval_AFR",
                              freq_col = "af_cases_AFR",
                              n = info_afr[['AFR.samples']][i],
                              min_freq = 0.01)
              
    }
    
    
    write.table(dat_AFR, file=paste0("./formatted_data/AFR_",phecode), quote = F, col.names = T, row.names = F)

    
}

# #### albumin (AFR)

# Read in dbSNP153 hg19 data
dbSNP <- readr::read_delim("../real_data_EAS/raw_data/EUR/dbSnp153Common.hg19.bed", col_names = FALSE)[, c(1:4)]
colnames(dbSNP) <- c("CHR", "START", "END", "rsID")

setDT(dbSNP)

setkey(dbSNP, CHR, END)

for(i in dim(info_afr)[1]){
    
    phecode =  info_afr[['EUR.phecode']][i]
    
    print(paste0("Trait: ", phecode, "..."))

    #data = data.table::fread(paste0('./raw_data/', 'AFR_alb_new.txt'))
    data = data.table::fread(paste0('./raw_data/', info_afr[['address']][i]))
    
    all_cols = colnames(data)

    data = data[!is.na(data$beta_re),]

    data[, c("CHR", "POS", "allele1", "allele2") := tstrsplit(snpid, ":", fixed = TRUE)]

    data$POS = as.numeric(data$POS)

    data$af = NA

    n = 0
    if('no_uganda' %in% all_cols){
        n1 = unique(data$no_uganda)
        if(length(n1)==1){
            n = n + n1
        }else{
            print('sample size not unique in uganda...')
            next
        }

        data$af = data$af_uganda
    }

    if('no_AADM' %in% all_cols){
        n2 = unique(data$no_AADM)
        if(length(n2)==1){
            n = n + n2
        }else{
            print('sample size not unique i AADM...')
            next
        }

        data[is.na(data$af),]$af = data[is.na(data$af),]$af_AADM
    }

    if('no_DCC' %in% all_cols){
        n3 = unique(data$no_DCC)
        if(length(n3)==1){
            n = n + n3
        }else{
            print('sample size not unique in DCC...')
            next
        }

        data[is.na(data$af),]$af = data[is.na(data$af),]$af_DCC
    }

    if('no_DDS' %in% all_cols){
        n4 = unique(data$no_DDS)
        if(length(n4)==1){
            n = n + n4
        }else{
            print('sample size not unique in DDS...')
            next
        }

        data[is.na(data$af),]$af = data[is.na(data$af),]$af_DDS

    }


    data$CHR_modified = paste0("chr", data$CHR)
    data$SNPID = data$snpid

    setDT(data)
    data[dbSNP, on = .(CHR_modified = CHR, POS = END), SNPID := i.rsID]

    data[is.na(data$SNPID)]$SNPID = data[is.na(data$SNPID)]$snpid
    data[(data$SNPID=='NA')]$SNPID = data[(data$SNPID=='NA')]$snpid


    data = format_data(data,
                      snp_col = "SNPID",
                      b_col = "beta_re",
                      se_col = "se_re",
                      A1_col = "allele1",
                      A2_col = "allele2",
                      p_col = "pval_re",
                      freq_col = "af",
                      n = n,
                      min_freq = 0.01)

    write.table(data, file=paste0("./formatted_data/AFR_",phecode), quote = F, col.names = T, row.names = F)
    
}

# #### albumin (EUR)

variants_eas = fread('../real_data_EAS/raw_data/EUR/variants.tsv')

head(variants_eas)

head(data_UKB)

phecode = '30600'

print(paste0("Trait: ", phecode, "..."))

UKB_address = "./raw_data/30600.tsv"
data_UKB = fread(UKB_address)

check1 = sum(data_UKB$variant != variants_eas$variant)

if(check1 != 0){
    print('Need more careful alignment with variants.tsv... skip')
}

data_UKB$SNP = variants_eas$rsid

data_UKB$A1 = variants_eas$alt
data_UKB$A2 = variants_eas$ref

dat_UKB = format_data(data_UKB,
              snp_col = "SNP",
              b_col = "beta",
              se_col = "se",
              A1_col = "A1",
              A2_col = "A2",
              p_col = "pval",
              freq_col = "minor_AF",
              n = 315268,
              min_freq = 0.01)

write.table(dat_UKB, file="./formatted_data/EUR_30600", quote = F, col.names = T, row.names = F)





# ### Case study: SHBG in EAS

i = which(info_csa$code == '30830')

phecode =  info_csa[['code']][i]

print(paste0("Trait: ", phecode, "..."))

data = data.table::fread(paste0('./raw_data/', info_csa[['address']][i]))

if(sum(!(data$chr == variants$chrom & data$pos == variants$pos)) != 0){
    print('Data and variants not fully matched. Skip...')
    next
}

data$rsid = variants$rsid
data = data[, c('chr','pos','rsid','ref','alt','af_EAS','af_EUR','beta_EAS','beta_EUR','se_EAS','se_EUR','neglog10_pval_EAS','neglog10_pval_EUR','low_confidence_EAS','low_confidence_EUR')]
data_sub = data[!(is.na(data$low_confidence_EAS) | is.na(data$low_confidence_EUR)),]

data_sub = data_sub[(data_sub$low_confidence_EAS == FALSE) & (data_sub$low_confidence_EUR == FALSE),]

data_sub$pval_EAS = 10^(-data_sub$neglog10_pval_EAS)

data_sub = as.data.frame(data_sub)


dat_EAS = format_data(data_sub,
                      snp_col = "rsid",
                      b_col = "beta_EAS",
                      se_col = "se_EAS",
                      A1_col = "alt",
                      A2_col = "ref",
                      p_col = "pval_EAS",
                      freq_col = "af_EAS",
                      n = 2327,
                      min_freq = 0.01)

write.table(dat_EAS, file=paste0("./formatted_data/EAS_",phecode), quote = F, col.names = T, row.names = F)


