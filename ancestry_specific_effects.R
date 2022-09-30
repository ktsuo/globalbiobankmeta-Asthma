#######################################################
########### Ancestry-specific heterogeneity ###########
#######################################################

library(readxl)
library(data.table)
library(tidyverse)
library(tidyr)
library(ggplot2)
library(writexl)

# set working directory
output_dir=""
setwd(output_dir)

######################
###### Load data #####
######################

# ancestry-specific meta-analyses
meta_files <- list.files("", pattern=glob2rx("*ANC.txt"), full.names=TRUE)
meta_list <- lapply(meta_files,function(x) fread(x, select=c('SNP_ID', 'CHR', 'POS', 'anc_inv_var_meta_beta', 'anc_inv_var_meta_sebeta')))

# all-asthma meta-analysis
all <- read.table("", sep="\t", header=T)

######################
###### Functions #####
######################

# STEP0: compute inverse-variance weighted meta-analysis beta estimates, SE, and P-value
meta.1 <- function(data){
  data <- data %>% mutate(num = anc_inv_var_meta_beta/(anc_inv_var_meta_sebeta^2),
                          denom = 1/(anc_inv_var_meta_sebeta^2))
}

meta.2 <- function(data){
  num = as.numeric(data[seq(2,10,2)])
  denom = as.numeric(data[seq(3,12,2)])
  meta = sum(num)/sum(denom)
  se = sqrt(1/(sum(denom)))
  return(list(b.meta = meta))
}

# STEP1: compute (1/ANC_SE^2*(ANC_BETA-ALL_BETA)^2) per SNP per ancestry meta-analysis
computeQ_step1 <- function(data, b.meta, beta, sebeta){
  data <- cbind(data, b.meta)
  names(data)[8] <- 'meta_beta'
  data <- data %>% mutate(het = (1/(sebeta^2))*((beta - meta_beta)^2))
}

anc_list <- list('AFR', 'AMR', 'EAS', 'EUR', 'CSA')
rename_cols <- function(data, anc_names){
  data <- data %>%
    rename_with(., function(x){paste(x,anc_names,sep="_")})
  return(data)
}

# STEP2: compute cochran's Q per SNP across ancestry meta-analyses
computeQ_step2 <- function(data){
  data <- data %>% mutate(Q = rowSums(select(.,starts_with('het'))), pval_het = pchisq(Q, df=5-1, lower=F)) 
  return(data)
}

########################
###### Cochran's Q #####
########################

meta_list <- lapply(meta_list, meta.1)
meta_list <- meta_list %>% lapply(., function(x) x %>% select(SNP_ID, num, denom)) %>% reduce(full_join, by='SNP_ID') %>% na.omit()
meta_list.betas <- apply(meta_list, 1, function(x) meta.2(x))

snps <- meta_list %>% reduce(full_join, by='SNP_ID') %>% na.omit() %>% select(SNP_ID)
meta_list <- meta_list %>% lapply(., function(x) x %>% filter(SNP_ID %in% SNPs$SNP_ID))
meta_list_q1 <- lapply(meta_list, function (x) computeQ_step1(data=x, b.meta=unlist(meta_list.betas), 
                                                              beta=anc_inv_var_meta_beta,
                                                              sebeta=anc_inv_var_meta_sebeta))

meta_list_q1 <- Map(rename_cols, meta_list_q1, anc_list)

# merge ancestry meta-analyses together
meta_all_q1 <- meta_list_q1 %>% reduce(full_join, by='SNP_ID')

meta_all_q2 <- computeQ_step2(meta_all_q1)


######################
####### Plot #########
######################

# get number of top hits that have p-values lower than Bonferroni threshold (significantly different across ancestries)
top_snps <- meta_all_q2 %>% filter(SNP_ID == '16:27344041:G:A' | SNP_ID == '10:9010779:G:A')

meta_betas <- data.frame(unlist(meta_list.betas))
meta.2.se <- function(data){
  num = as.numeric(data[seq(2,10,2)])
  denom = as.numeric(data[seq(3,12,2)])
  meta = sum(num)/sum(denom)
  se = sqrt(1/(sum(denom)))
  return(list(b.semeta = se))
}
meta_sebetas <- apply(meta_list, 1, function(x) meta.2.se(x))
meta_sebetas <- data.frame(unlist(meta_sebetas))
meta <- cbind(meta_all_q2$SNP_ID, meta_betas,meta_sebetas)
names(meta) <- c('SNP_ID', 'all_inv_var_meta_beta', 'all_inv_var_meta_sebeta')

# format top hits
format_tophits_for_forestplots <- function(tophits, all_meta){
  tophits <- tophits %>% left_join(., all_meta, by='SNP_ID')
  # format betas
  pivot1 <- tophits %>% dplyr:::select(c('SNP_ID', starts_with(c('anc_inv_var_meta_beta', 'all_inv_var_meta_beta')))) %>%
    pivot_longer(., cols=-c('SNP_ID'), names_to = "ancestry", values_to = "beta") %>%
    mutate(ancestry = ifelse(grepl("anc_inv_var_meta_beta_", ancestry),gsub("anc_inv_var_meta_beta_", "", ancestry),"ALL"))
  # format SEs
  pivot2 <- tophits %>% dplyr:::select(c('SNP_ID', starts_with(c('anc_inv_var_meta_sebeta', 'all_inv_var_meta_sebeta')))) %>%
    pivot_longer(., cols=-c('SNP_ID'), names_to = "ancestry", values_to = "SE") %>%
    mutate(ancestry = ifelse(grepl("anc_inv_var_meta_sebeta_", ancestry),gsub("anc_inv_var_meta_sebeta_", "", ancestry),"ALL"))
  
  tophits_long <- full_join(pivot1, pivot2, by=c("SNP_ID", "ancestry"))
  
  # order according to sample size of ancestry
  tophits_long <- tophits_long %>% group_by(SNP_ID) %>% mutate(ancestry=factor(ancestry,levels=c('AMR', 'CSA', 'AFR', 'EAS', 'EUR', "ALL"))) 
  
  # add column to indicate all_ancestries (for separate color/shape in forest plots)
  tophits_long$group <- ifelse(tophits_long$ancestry == "ALL", "2", "1") 
  
  # add new SNP name column
  tophits_long$SNP_ID.name <- paste('chr', tophits_long$SNP_ID, sep='')
  
  return(tophits_long)
}

top_snps <- format_tophits_for_forestplots(top_snps, meta)

forestplots <- top_snps %>% group_by(SNP_ID.name) %>% 
  ggplot(., aes(y = ancestry, x = beta, xmin=beta-1.96*SE, xmax=beta+1.96*SE)) +
  geom_point(aes(shape = group, color = group), size=5) +
  geom_errorbarh(height=.1, alpha=0.7, size=1, aes(color=group)) + 
  geom_vline(xintercept=0, color="black", linetype="dashed", alpha=.5) +
  scale_shape_manual(values = c(20,18)) +
  scale_color_manual(values = c("steelblue3", "orangered3")) +
  facet_wrap(~SNP_ID.name) +
  xlab('Effect size') + ylab('Ancestry') +
  theme_classic() +
  theme(legend.position = "none", text=element_text(size=20), panel.spacing=unit(2, "lines")) 

forestplots_output_name='.jpeg'
ggsave(forestplots, height = 6, width = 10, dpi = 300, filename=forestplots_output_name)
