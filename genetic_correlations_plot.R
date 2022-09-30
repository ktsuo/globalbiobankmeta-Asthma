#################################################
##### Genetic correlations in UKBB and BBJ ######
#################################################

library(readxl)
library(tidyverse)
library(dplyr)
library(corrplot)
library(qvalue)

# set working directory
dir=""
setwd(dir)

######################
###### Load data #####
######################

# import BBJ rg
BBJ_rg <- read.delim('BBJ_clinical.LDSC.txt')

# import BBJ heritability estimates 
BBJ_h2 <- read.delim('BBJ_heritabilities.txt')

# merge
BBJ <- left_join(BBJ_rg, BBJ_h2, by='trait')

# significantly heritable BBJ phenotypes 
bonf_thresh = 0.05/48
BBJ_sig <- BBJ %>% filter(p < bonf_thresh) 

# import UKB phenotypes that are significantly heritable in BBJ
UKB_matched_to_BBJ <- read.delim('GBMI-UKB_rg_asthma_phenos_in_BBJ.tsv')
# import UKB phenotypes that are significantly heritable in UKB
UKB <- read.delim('GBMI-UKB_rg_asthma.tsv')

# import UKB vs. UKB rgs
UKB_ref <- read.delim('UKB-UKB_rg_asthma.tsv')

########################
###### Format data #####
########################

# phenotypes heritable in BBJ
heritable_BBJ_phenos <- BBJ %>% 
  select(pheweb_disease_name, phecode1, phecode2) %>% 
  filter(!is.na(pheweb_disease_name))

# phenotypes heritable in UKB that are not heritable in BBJ
heritable_UKB_phenos <- UKB %>% filter(phenocode %in% BBJ$phecode1 | phenocode %in% BBJ$phecode2) %>% 
  filter(!(phenocode %in% heritable_BBJ_phenos$phecode1 | phenocode %in% heritable_BBJ_phenos$phecode2))
heritable_UKB_phenos <- BBJ %>% 
  select(pheweb_disease_name, phecode1) %>% right_join(., heritable_UKB_phenos, by=c('phecode1' = 'phenocode')) %>%
  select(pheweb_disease_name, phecode1) %>% mutate(phecode2 = NA)

# get union of heritable phenotypes in BBJ and UKB
heritable_phenos <- rbind(heritable_BBJ_phenos, heritable_UKB_phenos)

# GBMI-UKB rgs
UKB_rg1 <- UKB %>% filter(phenocode %in% heritable_phenos$phecode1 | phenocode %in% heritable_phenos$phecode2) %>%
  Select(phenocode, rg, rg_p, h2_z) %>% left_join(., heritable_phenos, by=c('phenocode' = 'phecode1'))
UKB_rg2 <- UKB_matched_to_BBJ %>% mutate(phenocode=as.character(phenocode)) %>% filter(phenocode %in% heritable_phenos$phecode1 | phenocode %in% heritable_phenos$phecode2) %>%
  select(phenocode, rg, rg_p, h2_z) %>% left_join(., heritable_phenos, by=c('phenocode' = 'phecode1')) %>% filter(!is.na(pheweb_disease_name))
UKB_rg <- rbind(UKB_rg1, UKB_rg2)
UKB_rg <- UKB_rg[!duplicated(UKB_rg$pheweb_disease_name),]
UKB_rg[UKB_rg == "NaN"] <- NA
UKB_rg <- UKB_rg %>% select(pheweb_disease_name, rg, rg_p, h2_z)
names(UKB_rg) <- c('pheweb_disease_name', 'GBMI_UKB_rg', 'GBMI_UKB_p', 'UKB_heritability_z')

# UKB-UKB rgs
UKB_vs_UKB_rg <- UKB_ref %>% mutate(phenocode = as.character(phenocode)) %>% 
  filter(phenocode %in% heritable_phenos$phecode1 | phenocode %in% heritable_phenos$phecode2) %>%
  select(phenocode, rg, rg_p, h2_z) %>% left_join(., heritable_phenos, by=c('phenocode' = 'phecode1')) %>% 
  filter(!is.na(pheweb_disease_name))
UKB_vs_UKB_rg[UKB_vs_UKB_rg == "NaN"] <- NA
UKB_vs_UKB_rg <- UKB_vs_UKB_rg %>% select(pheweb_disease_name, rg, rg_p, h2_z)
names(UKB_vs_UKB_rg) <- c('pheweb_disease_name', 'UKB_UKB_rg', 'UKB_UKB_p', 'UKB_UKB_heritability_z')

# BBJ
BBJ_rg <- BBJ %>% filter(pheweb_disease_name %in% heritable_phenos$pheweb_disease_name) %>% 
  select(pheweb_disease_name, GBMI_EAS_rg, GBMI_EAS_p, BBJ_EAS_rg, BBJ_EAS_p, p) %>% rename(.,'BBJ_heritability_p' = 'p')

# all rg
all_rg <- full_join(UKB_rg, BBJ_rg, by='pheweb_disease_name') %>% full_join(., UKB_vs_UKB_rg, by='pheweb_disease_name')
names(all_rg) <- c('pheweb_disease_name', 'GBMI_UKB_rg', 'GBMI_UKB_p', 'UKB_heritability_z', 'GBMI_BBJ_rg', 'GBMI_BBJ_p', 'BBJ_BBJ_rg', 'BBJ_BBJ_p', 'BBJ_heritability_p', 'UKB_UKB_rg', 'UKB_UKB_p', 'UKB_UKB_heritability_z')

########################
###### Corrplot ########
########################

all_rg_long <- all_rg %>% gather('GBMI_BBJ_rg', 'GBMI_UKB_rg', 'BBJ_BBJ_rg', 'GBMI_UKB_rg', 'UKB_UKB_rg',key = p1, value = rg) %>% 
  select(pheweb_disease_name, p1, rg)
all_rg_long2 <- all_rg %>% gather('GBMI_BBJ_p', 'GBMI_UKB_p', 'BBJ_BBJ_p', 'GBMI_UKB_p', 'UKB_UKB_p', key = p1, value = p) %>% 
  select(p)
all_rg_long_final <- cbind(all_rg_long, all_rg_long2) %>% mutate(p1 = gsub('_rg', '', p1))
names(all_rg_long_final) <- c('p1', 'p2', 'rg', 'p')

q <- qvalue(all_rg_long_final$p, fdr.level=0.05, pi0=1) # compute q values

all_rg_long_final$q <- q$qvalues # add column with q values to input2

rg = all_rg_long_final
tmp = rg
tmp$p1 = rg$p2
tmp$p2 = rg$p1
rg = rbind(rg, tmp)


x2 = reshape2::dcast(rg, p1 ~ p2, value.var = "rg")
mat2 = as.matrix(x2[, 2:ncol(x2)])
rownames(mat2) = x2$p1

mat2[mat2 > 1] = 1
mat2[mat2 < -1] = -1
mat2[is.na(mat2)] = NA

x2 = dcast(rg, p1 ~ p2, value.var = "q")
qmat2 = as.matrix(x2[, 2:ncol(x2)])
rownames(qmat2) = x2$p1
qmat2[is.na(qmat2)] = NA

if (nrow(mat2) == ncol(mat2)) {
  diag(mat2) = 1
  diag(qmat2) = -1
}

mat2_final <- mat2[-c(1,3,10,11,23),c('BBJ_BBJ', 'GBMI_BBJ', 'GBMI_UKB', 'UKB_UKB'), drop=FALSE]
qmat2_final <- qmat2[-c(1,3,10,11,23),c('BBJ_BBJ', 'GBMI_BBJ', 'GBMI_UKB', 'UKB_UKB'), drop=FALSE]

phenos_order <- c('Chronic obstructive pulmonary disease', 'Glaucoma', 'Myocardial infarction', 'Prostate cancer', 'Rheumatoid arthritis', 'Stable angina pectoris',
                  'Type 2 diabetes', 'Urolithiasis', 'Colorectal cancer', "Grave's disease", 'Ventricular arrhythmia', 'Chronic heart failure', 'Ischemic stroke', 'Cerebral aneurysm', 'Peripheral arterial disease',
                  'Pollinosis', 'Atopic dermatitis (including contact dermatitis)', 'Breast cancer', 'Cataract', 'Osteoporosis')

mat2_tmp <- as.data.frame(mat2_final) %>% mutate(phenos = factor(row.names(mat2_final), levels = phenos_order)) %>% arrange(phenos) %>% select(-phenos)
qmat2_tmp <- as.data.frame(qmat2_final) %>% mutate(phenos = factor(row.names(mat2_final), levels = phenos_order)) %>% arrange(phenos) %>% select(-phenos)
mat2_final2 <- as.matrix(mat2_tmp)
qmat2_final2 <- as.matrix(qmat2_tmp)


png('.png', width = 7, height = 8, units = 'in', res = 300, family = "Helvetica")
corrplot(corr=mat2_final2, method = "square",insig = "label_sig", order = 'original',
         pch = "*", pch.cex = 1.5, p.mat = qmat2_final2, sig.level = 0.05/20,
         tl.col="cyan3", na.label = "square", na.label.col = "grey80", cl.ratio = 0.4)


phenos_notHeritable <- c('Colorectal cancer', "Grave's disease", 'Ventricular arrhythmia', 'Chronic heart failure', 'Ischemic stroke', 'Cerebral aneurysm', 'Peripheral arterial disease',
                         'Pollinosis', 'Atopic dermatitis (including contact dermatitis)')
phenos_notHeritable_index <- which(rownames(mat2_final2) %in% phenos_notHeritable)
rownames(mat2_final2)[phenos_notHeritable_index] <- ""
rownames(qmat2_final2)[phenos_notHeritable_index] <- ""
corrplot(corr=mat2_final2, method = "square",insig = "label_sig", order = 'original',
         pch = "*", pch.cex = 1.5, p.mat = qmat2_final2, sig.level = 0.05/20,
         tl.col="palegreen3", cl.pos='n',na.label = "square", na.label.col = "grey80", add=TRUE)

phenos_HeritableUKBonly <- c('Breast cancer', 'Cataract', 'Osteoporosis')
phenos_HeritableUKBonly_index <- which(rownames(mat2_final2) %in% phenos_HeritableUKBonly)
rownames(mat2_final2)[phenos_HeritableUKBonly_index] <- ""
rownames(qmat2_final2)[phenos_HeritableUKBonly_index] <- ""
corrplot(corr=mat2_final2, method = "square",insig = "label_sig", order = 'original',
         pch = "*", pch.cex = 1.5, p.mat = qmat2_final2, sig.level = 0.05/20,
         tl.col="black", cl.pos='n',na.label = "square", na.label.col = "grey80", add=TRUE)

corrplot_legend <- c('Heritable in all', 'Heritable in BBJ', 'Heritable in UKB')
corrplot_legendcol <- c('black', 'cyan3', 'palegreen3')
legend(x='bottomleft',legend=corrplot_legend,lty=1, col=corrplot_legendcol, lwd=4)


dev.off()