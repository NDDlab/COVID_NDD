library(data.table)
library(tidyverse)

##################################################
#                    Finn                        
##################################################
# ASD --------------------------------------------
asd <- fread("finngen_R9_KRA_PSY_AUTISM_EXMORE",fill=T)
head(asd)
asd <- asd %>% 
  dplyr::mutate(n="NA", freq="NA") %>% 
  dplyr::select(SNP=5,A1=3,A2=4,freq,b=beta,se=sebeta,p=pval,n)
head(asd)
# SNP	A1	A2	freq	b	se	p	n
# rs2691328 G A NA -0.626368 0.741497 0.398259 NA
# rs878915777 C T NA 0.249283 0.514392 0.627948 NA
fwrite(asd,"smr_input_finn_asd_gwas.txt",sep = "\t")

# BIP --------------------------------------------
bp <- fread("finngen_R9_F5_BIPO",fill=T)
head(bp)
bp <- bp %>% dplyr::mutate(n="NA", freq="NA") %>% 
  dplyr::select(SNP=5,A1=3,A2=4,freq,b=beta,se=sebeta,p=pval,n)
head(bp)
fwrite(bp,"smr_input_finn_bip_gwas.txt",sep = "\t")

# DEP --------------------------------------------
dep <- fread("finngen_R9_F5_DEPRESSIO",fill=T)
head(dep)
dep <- dep %>% dplyr::mutate(n="NA", freq="NA") %>% 
  dplyr::select(SNP=5,A1=3,A2=4,freq,b=beta,se=sebeta,p=pval,n)
head(dep)
fwrite(dep,"smr_input_finn_dep_gwas.txt",sep = "\t")

# SCZ --------------------------------------------
scz <- fread("finngen_R9_F5_SCHZPHR",fill=T)
head(scz)
scz <- scz %>% dplyr::mutate(n="NA", freq="NA") %>% 
  dplyr::select(SNP=5,A1=3,A2=4,freq,b=beta,se=sebeta,p=pval,n)
fwrite(scz,"smr_input_finn_scz_gwas.txt",sep = "\t")

##################################################
#                    PGC                         
##################################################
# ASD --------------------------------------------
asd <- fread("iPSYCH-PGC_ASD_Nov2017",fill=T)
head(asd)
asd <- asd %>% dplyr::mutate(b=log(OR),n="NA", freq="NA") %>% 
  dplyr::select(SNP,A1,A2,freq,b,se=SE,p=P,n)
head(asd)
fwrite(asd,"smr_input_pgc_asd_gwas.txt",sep = "\t")

# BIP --------------------------------------------
bd <- fread("pgc-bip2021-all.vcf.tsv",fill=T)
head(bd)
bd <- bd %>% dplyr::mutate(n=NCAS+NCON, freq="NA") %>% 
  dplyr::select(SNP=ID,A1,A2,freq,b=BETA,se=SE,p=PVAL,n)
head(bd)
fwrite(bd,"smr_input_pgc_bip_gwas.txt",sep = "\t")

# MDD --------------------------------------------
mdd <- fread("PGC_UKB_depression_genome-wide.txt",fill=T)
head(mdd)
mdd <- mdd %>% dplyr::mutate(n="NA") %>% 
  dplyr::select(SNP=1,A1,A2,freq=Freq,b=LogOR,se=StdErrLogOR,p=P,n)
head(mdd)
fwrite(mdd,"smr_input_pgc_mdd_gwas.txt",sep = "\t")

# SCZ --------------------------------------------
scz <- fread("PGC3_SCZ_wave3.core.autosome.public.v3.vcf.tsv",fill=T)
dim(scz)
scz <- scz %>% dplyr::mutate(n=NCAS+NCON, freq="NA") %>% 
  dplyr::select(SNP=ID,A1,A2,freq,b=BETA,se=SE,p=PVAL,n)
fwrite(scz,"smr_input_pgc_scz_gwas.txt",sep = "\t")