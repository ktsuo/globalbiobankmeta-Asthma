##############################################################
##### Biobanks VS. Leave-that-biobank-out meta-analyses ######
##############################################################

library(data.table)
library(ggplot2)
library(dplyr)
library(broom)
library(readxl)
library(ggrepel)
library(writexl)
library(tidyverse)
library(purrr)

# set working directory
output_dir=""
setwd(output_dir)

######################
###### Load data #####
######################

results_dir = ""

# load biobank results
biobank_files <- list.files(results_dir, 
                            pattern=glob2rx("*BBONLY.txt"), full.names=TRUE)
list.biobanks <- lapply(biobank_files,function(x) fread(x, select=c('SNP_ID', 'CHR', 'POS', 'Allele1', 'Allele2', 'biobank_BETA', 'biobank_SE', 'biobank_p.value')))
list.biobanks.names <- list('BBJ', 'CKB', 'DECODE', 'ESTBB', 'FinnGen', 'GNH', 'GS', 'HUNT', 'Lifelines', 'MGB', 'QSKIN', 'TWB')
names(list.biobanks) <- list.biobanks.names

biobank_meta_files <- list.files(results_dir, pattern=glob2rx("*BBMETA.txt"), full.names=TRUE)
list.bb_meta <- lapply(biobank_meta_files, function(x) fread(x, select=c(1:8)))
list.bb_meta.names <- list('BioMe', 'BioVU', 'CCPM', 'MGI', 'UCLA', 'UKBB')
names(list.bb_meta) <- list.bb_meta.names

# load LOO meta-analyses
loo_files <- list.files(results_dir, pattern=glob2rx("*LOO.txt"), full.names=TRUE)
list.loo <- lapply(loo_files,function(x) fread(x, select=c('SNP_ID','CHR', 'POS', 'REF', 'ALT', 'inv_var_meta_beta', 'inv_var_meta_sebeta',
                                                           'inv_var_meta_p')))
list.loo.names <- list('BBJ', 'BIOME', 'BIOVU', 'CCPM', 'CKB', 'DECODE', 'ESTBB', 'FG', 'GNH', 'GS', 'HUNT', 'LIFELINES', 'MGB', 'MGI', 'QSKIN', 'TWB', 'UCLA', 'UKBB')
names(list.loo) <- list.loo.names

# load top hits
top_hits <- read_excel("")
top_hits <- top_hits %>% 
  select(SNP, CHR, POS, REF, ALT, all_inv_var_meta_p)

# load sample sizes 
sample_sizes <- read_excel("biobanks_sampling_prevalence.xlsx")

########################
###### Format data #####
########################

# combine all biobank-specific meta-analyses into one list
list.all.biobanks <- c(list.bb_meta, list.biobanks)
list.all.biobanks <- lapply(list.all.biobanks, setNames, nm = c('SNP_ID','CHR', 'POS', 'REF', 'ALT', 'biobank_beta', 'biobank_sebeta', 'biobank_p'))

list.all.biobanks <- list.all.biobanks[order(names(list.all.biobanks))]
list.loo <- list.loo[order(names(list.loo))]

# align to LOO risk allele
add_risk_allele <- function(x, beta_col, ALT_col, REF_col){
  x %>% mutate(risk_allele = ifelse(beta_col > 0, ALT_col, REF_col),
               risk_allele_beta = ifelse(beta_col > 0, beta_col, -(beta_col)))
}

list.loo <- lapply(list.loo, add_risk_allele, beta_col=inv_var_meta_beta, ALT_col=ALT, REF_col=REF)
list.loo <- lapply(list.loo, setNames, nm = c('SNP_ID', 'CHR', 'POS', 'REF', 'ALT', 'inv_var_meta_beta', 'inv_var_meta_sebeta', 'inv_var_meta_p ', 'risk_allele_LOO', 'risk_allele_beta_LOO'))

align_alleles <- function(loo, biobank, beta_col, ALT_col, REF_col){
  merge <- right_join(loo,biobank, by='SNP_ID')
  merge <- merge %>% mutate(matched_allele_BIOBANK = ifelse(ALT_col == risk_allele_LOO, ALT_col, REF_col),
                            matched_allele_beta_BIOBANK = ifelse(ALT_col == risk_allele_LOO, beta_col, -(beta_col)))
}

list.aligned <- Map(align_alleles, list.loo, list.all.biobanks, beta_col=biobank_beta, ALT_col=ALT.y, REF_col=REF.y)

# add column with biobank names
biobank_names <- c(names(list.all.biobanks)[order(names(list.all.biobanks))])
list.aligned <- Map(function(x,y){x <- x %>% mutate(biobank = y)}, list.aligned, biobank_names)

df <- Reduce('rbind',list.aligned)

###################################################
##### Plot ratios of biobank vs. LOO effects ######
###################################################

# compute beta ratios
df <- df %>% mutate(beta_ratio = matched_allele_beta_BIOBANK / risk_allele_beta_LOO)
averages <- df %>% group_by(biobank) %>% summarise(average = mean(beta_ratio), n=n())

# plot average beta ratios
output_name='.jpeg'
text_size=25

barplot_ratios <- ggplot(averages, aes(x=biobank, y=average)) + 
  geom_bar(stat="identity", fill="steelblue3", width=0.5) +
  geom_hline(yintercept = 1.0, linetype='dashed', size=1) +
  ylab(paste0("Average SNP effect in biobank over", "\n", 
              "leave-that-biobank-out meta-analysis")) +
  theme_classic() + 
  theme(text=element_text(size=text_size),axis.text.x=element_text(angle=65, hjust=1), axis.title.y=element_text(size=20))

ggsave(barplot_ratios, height = 9, width = 12, dpi = 300, filename=output_name)

###################################
##### Plot Deming Regressions #####
###################################

df1 <- lapply(list.aligned, as.data.frame)

ggplotDemingRegression <- function(fit, df, param) { 
  require(ggplot2)

  deming_slope = unlist(fit$coefficients[2])
  p=ggplot(fit$model, aes_string(x = names(fit$model)[2], y = names(fit$model)[1])) + 
    geom_errorbar(aes(ymin=df$lowV, ymax=df$highV), width=.005, alpha=0.8) + 
    geom_errorbarh(aes( xmin=df$lowH,xmax=df$highH), height=.005, alpha=0.8) + 
    geom_abline(slope=deming_slope, intercept=0, color = "orangered3") + 
    theme_bw() + 
    theme(plot.title = element_text(hjust = 0.5,size=18), axis.text=element_text(size=10),
          axis.title=element_text(size=10,face="bold"),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          legend.position = "none",
          plot.margin=grid::unit(c(5,5,5,5), "mm")) +
    geom_point(size=0.5) +
    geom_abline(intercept = 0, slope=1, col="grey3", show.legend = T, linetype = "dashed") +
    xlab(param$xlabtext) + ylab(param$ylabtext)+ 

    scale_x_continuous(expand=c(0,0), limits=c(param$MinVal-0.1,param$MaxVal+0.1)) +
    scale_y_continuous(expand=c(0,0), limits=c(param$MinVal-0.1,param$MaxVal+0.1)) +
    coord_cartesian(xlim=c(param$MinVal,param$MaxVal), ylim=c(param$MinVal,param$MaxVal)) +
    labs(title = paste(
      " Slope =",signif(deming_slope, 2)
    )
    )
  ggsave(p,filename=param$pdffile, dpi=300,width = 5, height = 5)
}

# deming regression models for each biobank
fit_deming <- function(dat){
  require(deming)
  require(devtools)
  demingfit <- deming(risk_allele_beta_LOO ~ matched_allele_beta_BIOBANK + 0, data=dat, xstd=biobank_sebeta, ystd=inv_var_meta_sebeta)
}

# plot parameters
ggplotDemingRegression_input_step1 <- function(data){
  data <- data %>% mutate(lowH = matched_allele_beta_BIOBANK - 1.96*biobank_sebeta,
                          highH = matched_allele_beta_BIOBANK + 1.96*biobank_sebeta,
                          lowV = risk_allele_beta_LOO - 1.96*inv_var_meta_sebeta,
                          highV = risk_allele_beta_LOO + 1.96*inv_var_meta_sebeta)
}

ggplotDemingRegression_input_step2 <- function(data, biobank_name){
  xlabtext <- paste("Effect sizes reported by", biobank_name, sep=" ")
  ylabtext <- paste("Effect sizes reported by LOO meta-analysis excluding", biobank_name, sep=" ")
  MinVal <- min(0, min(data$lowH)-0.01, min(data$lowV)-0.01)
  MaxVal <- max(max(data$highH)+0.01, max(data$highV)+0.01)
  pdffile <- paste0('tophits_bb_loo_effectsize_comparison_DemingRegression_allSNPs', biobank_name, '.jpeg')
  input.df <- data.frame(xlabtext, ylabtext, as.numeric(MaxVal), as.numeric(MinVal), pdffile)
  names(input.df) <- c('xlabtext', 'ylabtext', 'MaxVal', 'MinVal', 'pdffile')
  return(input.df)
}

# plot

list.deming <- lapply(df1, FUN=fit_deming)
df.deming <- lapply(df1, ggplotDemingRegression_input_step1)
plot.param <- Map(ggplotDemingRegression_input_step2, df.deming, biobank_names)

Map(ggplotDemingRegression, fit=list.deming, df=df.deming, param=plot.param)
