########################
##### COA vs. AOA ######
########################

library(data.table)
library(ggplot2)
library(dplyr)
library(broom)
library(writexl)
library(readr)
library(devtools)
library(deming)

# set working directory
dir=""
setwd(dir)

######################
###### Load data #####
######################

# index variants in AOA and COA meta-analysis
AOA <- read_delim("AsthmaAOA_Bothsex_inv_var_meta.asthma_all_tophits.txt", 
                  delim = "\t", escape_double = FALSE, 
                  trim_ws = TRUE)
COA <- read_delim("AsthmaCOA_Bothsex_inv_var_meta.asthma_all_tophits.txt", 
                  delim = "\t", escape_double = FALSE, 
                  trim_ws = TRUE)

AOA <- AOA %>% rename(AOA_beta = meta_BETA, AOA_sebeta = meta_SE, AOA_pval = meta_p.value)
COA <- COA %>% rename(COA_beta = meta_BETA, COA_sebeta = meta_SE, COA_pval = meta_p.value)

########################
###### Format data #####
########################

# risk allele beta
add_risk_allele <- function(x, beta_col, ALT_col, REF_col){
  x %>% mutate(risk_allele = ifelse(beta_col > 0, ALT_col, REF_col),
               risk_allele_beta = ifelse(beta_col > 0, beta_col, -(beta_col)))
}

AOA <- add_risk_allele(AOA, AOA_beta, Allele2, Allele1)

# align risk alelles
align_alleles <- function(AOA, COA){
  merge <- right_join(AOA, COA, by='SNP_ID')
  merge <- merge %>% mutate(matched_allele = ifelse(Allele2.y == risk_allele, Allele2.y, Allele1.y),
                            matched_allele_beta = ifelse(Allele2.y == risk_allele, COA_beta, -(COA_beta)))
}

AOA_COA_aligned <- align_alleles(AOA, COA)

###################################
###### Plot Deming Regression #####
###################################

# deming regression
demingfit <- deming(matched_allele_beta ~ risk_allele_beta + 0, data=AOA_COA_aligned, xstd=AOA_sebeta, ystd=COA_sebeta)

slope1 = unlist(demingfit$coefficients[2])
intercept1 = 0

deming <- ggplot(AOA_COA_aligned, aes(x = risk_allele_beta, y = matched_allele_beta)) +
  geom_point(size=3) +
  geom_errorbarh(aes(xmin = risk_allele_beta - 1.96*AOA_sebeta, xmax = risk_allele_beta + 1.96*AOA_sebeta), size=0.3, alpha=0.5) +
  geom_errorbar(aes(ymin = matched_allele_beta- 1.96*COA_sebeta, ymax = matched_allele_beta + 1.96*COA_sebeta), size=0.3, alpha=0.5) +
  labs(title = paste(
    " Slope =",signif(slope1,2))) +
  geom_abline(slope=slope1, intercept=intercept1, color="orangered3", size=0.8) + 
  geom_abline(slope = 1, intercept = 0, linetype='dashed', alpha=0.5) +
  theme_classic() +
  xlim(c(0, 0.4)) + ylim(c(-0.1,0.4)) +
  xlab('AOA risk allele effect size') + ylab('COA matched allele effect size') +
  theme(axis.text = element_text(size=15), axis.title = element_text(size=15), title=element_text(size=15)) 

deming_output_name = '.jpeg'
ggsave(deming, height = 7, width = 9, dpi = 300, filename=deming_output_name)