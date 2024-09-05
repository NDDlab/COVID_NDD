library(tidyverse)
library(data.table)

######################################################
#                     1. RNA-Seq brain cell types deconvolution
# ref: # https://voineagulab.shinyapps.io/BrainDeconvShiny/
######################################################
library(CIBERSORT)

get_file<-function(){
  VL <<- "BrainDeconv/VL.txt"
}  
get_file()

# 1.1 TCx ----
lnames=load("2.output/4.ndd_wgcna/data/02_Consensus_dataInput_TC.RData")

for (i in 1:exprSize_TC$nSets) {
  
  data <- as.data.frame(t(multiExpr_TC[[i]]$data))
  data <- 2^data-1
  
  filename <- paste0("2.output/4.ndd_wgcna/data/02_dataInput_Expmat/TC_",setLabels_TC[i],".txt")
  fwrite(data,filename,sep="\t",row.names = T)
  
  cibersort_VL = cibersort(VL, filename, perm = 1000, QN = T)
  filenameVL <- paste0("2.output/4.ndd_wgcna/data/02_dataInput_Expmat/VL_TC_",setLabels_TC[i],".txt")
  fwrite(cibersort_VL %>% as.data.frame(),filenameVL,sep="\t",row.names = T)
}

# 1.2 PCx ----
lnames=load("2.output/4.ndd_wgcna/data/02_Consensus_dataInput_PC.RData")

for (i in 1:exprSize_PC$nSets) {
  
  data <- as.data.frame(t(multiExpr_PC[[i]]$data))
  data <- 2^data-1
  
  filename <- paste0("2.output/4.ndd_wgcna/data/02_dataInput_Expmat/PC_",setLabels_PC[i],".txt")
  fwrite(data,filename,sep="\t",row.names = T)
  
  cibersort_VL = cibersort(VL, filename, perm = 1000, QN = T)
  filenameVL <- paste0("2.output/4.ndd_wgcna/data/02_dataInput_Expmat/VL_PC_",setLabels_PC[i],".txt")
  fwrite(cibersort_VL %>% as.data.frame(),filenameVL,sep="\t",row.names = T)
}

# 1.3 CB ----
lnames=load("2.output/4.ndd_wgcna/data/02_Consensus_dataInput_CB-re.RData")

for (i in 1:exprSize_CB$nSets) {
  
  data <- as.data.frame(t(multiExpr_CB[[i]]$data))
  data <- 2^data-1
  
  filename <- paste0("2.output/4.ndd_wgcna/data/02_dataInput_Expmat/CB_",setLabels_CB[i],".txt")
  fwrite(data,filename,sep="\t",row.names = T)
  
  cibersort_VL = cibersort(VL, filename, perm = 1000, QN = T)
  filenameVL <- paste0("2.output/4.ndd_wgcna/data/02_dataInput_Expmat/VL_CB_",setLabels_CB[i],".txt")
  fwrite(cibersort_VL %>% as.data.frame(),filenameVL,sep="\t",row.names = T)
}

# 1.4 HIP ----
lnames=load("2.output/4.ndd_wgcna/data/02_Consensus_dataInput_HIP.RData")

for (i in 1:exprSize_HIP$nSets) {
  
  data <- as.data.frame(t(multiExpr_HIP[[i]]$data))
  data <- 2^data-1
  
  filename <- paste0("2.output/4.ndd_wgcna/data/02_dataInput_Expmat/HIP_",setLabels_HIP[i],".txt")
  fwrite(data,filename,sep="\t",row.names = T)
  
  cibersort_VL = cibersort(VL, filename, perm = 1000, QN = T)
  filenameVL <- paste0("2.output/4.ndd_wgcna/data/02_dataInput_Expmat/VL_HIP_",setLabels_HIP[i],".txt")
  fwrite(cibersort_VL %>% as.data.frame(),filenameVL,sep="\t",row.names = T)
}

# 1.5 OCx ----
lnames=load("2.output/4.ndd_wgcna/data/02_NonConsensus_dataInput_OC_asd.RData")

data <- as.data.frame(t(asd_OC))
data <- 2^data-1

filename <- paste0("2.output/4.ndd_wgcna/data/02_dataInput_Expmat/OC_ASD.txt")
fwrite(data,filename,sep="\t",row.names = T)

cibersort_VL = cibersort(VL, filename, perm = 1000, QN = T)
filenameVL <- "2.output/4.ndd_wgcna/data/02_dataInput_Expmat/VL_OC_ASD.txt"
fwrite(cibersort_VL %>% as.data.frame(),filenameVL,sep="\t",row.names = T)

######################################################
#                     2. Module Cibersort Correlation
# Calculate module enginegene (MEs) -> Correlation
######################################################

## 2.1 CB ----
## 2.1.1 import expression and cibersort result ----
load(file = "2.output/4.ndd_wgcna/data/02_Consensus_dataInput_CB-re.RData")
load(file = "2.output/4.ndd_wgcna/data/03_Consensus-NetworkConstruction-CB.RData")

vl_cb_bd <-fread("2.output/4.ndd_wgcna/data/02_dataInput_Expmat/VL_CB_BD.txt")
vl_cb_dep <-fread("2.output/4.ndd_wgcna/data/02_dataInput_Expmat/VL_CB_DEP.txt")
vl_cb_scz <-fread("2.output/4.ndd_wgcna/data/02_dataInput_Expmat/VL_CB_SCZ.txt")

vl_cb <- list()
vl_cb[[1]] <- vl_cb_bd %>% column_to_rownames("V1") %>% dplyr::select(1:8)
vl_cb[[2]] <- vl_cb_dep %>% column_to_rownames("V1") %>% dplyr::select(1:8)
vl_cb[[3]] <- vl_cb_scz %>% column_to_rownames("V1") %>% dplyr::select(1:8)

## 2.1.2 Calculate the correlations -----
# Set up variables to contain the module-trait correlations
moduleTraitCor_CB<-list()
moduleTraitPvalue_CB<-list()
moduleTraitPadj_CB<-list()
# Calculate the correlations
for (set in 1:nSets_CB)
{
  moduleTraitCor_CB[[set]] = cor(consMEs_CB[[set]]$data, vl_cb[[set]], use = "p")
  moduleTraitPvalue_CB[[set]] = corAndPvalue(consMEs_CB[[set]]$data, vl_cb[[set]])$p
  moduleTraitPadj_CB[[set]] = matrix(p.adjust(moduleTraitPvalue_CB[[set]], method="BH"),ncol = 8)
  
  colnames(moduleTraitCor_CB[[set]])<-colnames(moduleTraitPvalue_CB[[set]])<-colnames(moduleTraitPadj_CB[[set]])<-colnames(vl_cb[[1]])
}

save(moduleTraitCor_CB,moduleTraitPvalue_CB,moduleTraitPadj_CB,
        file = "2.output/4.ndd_wgcna/data/08_ModuleCibersortCor_CB.Rdata")

## 2.2 HIP ----
## 2.2.1 import expression and cibersort result ----
load(file = "2.output/4.ndd_wgcna/data/02_Consensus_dataInput_HIP.RData")
load(file = "2.output/4.ndd_wgcna/data/03_Consensus-NetworkConstruction-HIP.RData")

vl_hip_bd <-fread("2.output/4.ndd_wgcna/data/02_dataInput_Expmat/VL_HIP_BD.txt")
vl_hip_dem <-fread("2.output/4.ndd_wgcna/data/02_dataInput_Expmat/VL_HIP_DEM.txt")
vl_hip_mdd <-fread("2.output/4.ndd_wgcna/data/02_dataInput_Expmat/VL_HIP_MDD.txt")
vl_hip_scz <-fread("2.output/4.ndd_wgcna/data/02_dataInput_Expmat/VL_HIP_SCZ.txt")

vl_hip <- list()
vl_hip[[1]] <- vl_hip_bd %>% column_to_rownames("V1") %>% dplyr::select(1:8)
vl_hip[[2]] <- vl_hip_dem %>% column_to_rownames("V1") %>% dplyr::select(1:8)
vl_hip[[3]] <- vl_hip_mdd %>% column_to_rownames("V1") %>% dplyr::select(1:8)
vl_hip[[4]] <- vl_hip_scz %>% column_to_rownames("V1") %>% dplyr::select(1:8)

## 4.2.2 Calculate the correlations -----
# Set up variables to contain the module-trait correlations
moduleTraitCor_HIP<-list()
moduleTraitPvalue_HIP<-list()
moduleTraitPadj_HIP<-list()
# Calculate the correlations
for (set in 1:nSets_HIP)
{
  moduleTraitCor_HIP[[set]] = cor(consMEs_HIP[[set]]$data, vl_hip[[set]], use = "p")
  moduleTraitPvalue_HIP[[set]] = corAndPvalue(consMEs_HIP[[set]]$data, vl_hip[[set]])$p
  moduleTraitPadj_HIP[[set]] = matrix(p.adjust(moduleTraitPvalue_HIP[[set]], method="BH"),ncol = 8)
  
  colnames(moduleTraitCor_HIP[[set]])<-colnames(moduleTraitPvalue_HIP[[set]])<-colnames(moduleTraitPadj_HIP[[set]])<-colnames(vl_hip[[1]])
}

save(moduleTraitCor_HIP,moduleTraitPvalue_HIP,moduleTraitPadj_HIP,
        file = "2.output/4.ndd_wgcna/data/08_ModuleCibersortCor_HIP.Rdata")

## 2.3 OC ----
## 2.3.1 import expression and cibersort result ----
load(file = "2.output/4.ndd_wgcna/data/02_NonConsensus_dataInput_OC_asd.RData")
load(file = "2.output/4.ndd_wgcna/data/03_NonConsensus-NetworkConstruction-OCx.RData")

vl_oc_asd <-fread("2.output/4.ndd_wgcna/data/02_dataInput_Expmat/VL_OC_ASD.txt")

## 4.3.2 Calculate the correlations -----
# Calculate the correlations
moduleTraitCor_OC = cor(MEs_OC, vl_oc_asd%>% column_to_rownames("V1") %>% dplyr::select(1:8), use = "p")
moduleTraitPvalue_OC = corAndPvalue(MEs_OC, vl_oc_asd%>% column_to_rownames("V1") %>% dplyr::select(1:8))$p
moduleTraitPadj_OC = matrix(p.adjust(moduleTraitPvalue_OC, method="BH"),ncol = 8)

save(moduleTraitCor_OC,moduleTraitPvalue_OC,moduleTraitPadj_OC,
        file = "2.output/4.ndd_wgcna/data/08_ModuleCibersortCor_OC.Rdata")

## 2.4 PC ----
## 2.4.1 import expression and cibersort result ----
load(file = "2.output/4.ndd_wgcna/data/02_Consensus_dataInput_PC.RData")
load(file = "2.output/4.ndd_wgcna/data/03_Consensus-NetworkConstruction-PCx.RData")

vl_pc_bd <-fread("2.output/4.ndd_wgcna/data/02_dataInput_Expmat/VL_PC_BD.txt")
vl_pc_dem <-fread("2.output/4.ndd_wgcna/data/02_dataInput_Expmat/VL_PC_DEM.txt")
vl_pc_dep <-fread("2.output/4.ndd_wgcna/data/02_dataInput_Expmat/VL_PC_DEP.txt")
vl_pc_scz <-fread("2.output/4.ndd_wgcna/data/02_dataInput_Expmat/VL_PC_SCZ.txt")

vl_pc <- list()
vl_pc[[1]] <- vl_pc_bd %>% column_to_rownames("V1") %>% dplyr::select(1:8)
vl_pc[[2]] <- vl_pc_dep %>% column_to_rownames("V1") %>% dplyr::select(1:8)
vl_pc[[3]] <- vl_pc_dem %>% column_to_rownames("V1") %>% dplyr::select(1:8)
vl_pc[[4]] <- vl_pc_scz %>% column_to_rownames("V1") %>% dplyr::select(1:8)

## 4.4.2 Calculate the correlations -----
# Set up variables to contain the module-trait correlations
moduleTraitCor_PC<-list()
moduleTraitPvalue_PC<-list()
moduleTraitPadj_PC<-list()
# Calculate the correlations
for (set in 1:nSets_PC)
{
  moduleTraitCor_PC[[set]] = cor(consMEs_PC[[set]]$data, vl_pc[[set]], use = "p")
  moduleTraitPvalue_PC[[set]] = corAndPvalue(consMEs_PC[[set]]$data, vl_pc[[set]])$p
  moduleTraitPadj_PC[[set]] = matrix(p.adjust(moduleTraitPvalue_PC[[set]], method="BH"),ncol = 8)
  
  colnames(moduleTraitCor_PC[[set]])<-colnames(moduleTraitPvalue_PC[[set]])<-colnames(moduleTraitPadj_PC[[set]])<-colnames(vl_pc[[1]])
}

save(moduleTraitCor_PC,moduleTraitPvalue_PC,moduleTraitPadj_PC,
        file = "2.output/4.ndd_wgcna/data/08_ModuleCibersortCor_PC.Rdata")

## 2.5 TC ----
## 2.5.1 import expression and cibersort result ----
load(file = "2.output/4.ndd_wgcna/data/02_Consensus_dataInput_TC.RData")
load(file = "2.output/4.ndd_wgcna/data/03_Consensus-NetworkConstruction-TCx.RData")

vl_tc_asd <-fread("2.output/4.ndd_wgcna/data/02_dataInput_Expmat/VL_TC_ASD.txt")
vl_tc_dem <-fread("2.output/4.ndd_wgcna/data/02_dataInput_Expmat/VL_TC_DEM.txt")
vl_tc_scz <-fread("2.output/4.ndd_wgcna/data/02_dataInput_Expmat/VL_TC_SCZ.txt")

vl_tc <- list()
vl_tc[[1]] <- vl_tc_asd %>% column_to_rownames("V1") %>% dplyr::select(1:8)
vl_tc[[2]] <- vl_tc_dem %>% column_to_rownames("V1") %>% dplyr::select(1:8)
vl_tc[[3]] <- vl_tc_scz %>% column_to_rownames("V1") %>% dplyr::select(1:8)

## 2.5.2 Calculate the correlations -----
# Set up variables to contain the module-trait correlations
moduleTraitCor_TC<-list()
moduleTraitPvalue_TC<-list()
moduleTraitPadj_TC<-list()
# Calculate the correlations
for (set in 1:nSets_TC)
{
  moduleTraitCor_TC[[set]] = cor(consMEs_TC[[set]]$data, vl_tc[[set]], use = "p")
  moduleTraitPvalue_TC[[set]] = corAndPvalue(consMEs_TC[[set]]$data, vl_tc[[set]])$p
  moduleTraitPadj_TC[[set]] = matrix(p.adjust(moduleTraitPvalue_TC[[set]], method="BH"),ncol = 8)
  
  colnames(moduleTraitCor_TC[[set]])<-colnames(moduleTraitPvalue_TC[[set]])<-colnames(moduleTraitPadj_TC[[set]])<-colnames(vl_tc[[1]])
}

save(moduleTraitCor_TC,moduleTraitPvalue_TC,moduleTraitPadj_TC,
        file = "2.output/4.ndd_wgcna/data/08_ModuleCibersortCor_TC.Rdata")

# 2.6 collect ----
## construct data ----
load("2.output/4.ndd_wgcna/data/08_ModuleCibersortCor_CB.Rdata")
load("2.output/4.ndd_wgcna/data/08_ModuleCibersortCor_HIP.Rdata")
load("2.output/4.ndd_wgcna/data/08_ModuleCibersortCor_OC.Rdata")
load("2.output/4.ndd_wgcna/data/08_ModuleCibersortCor_PC.Rdata")
load("2.output/4.ndd_wgcna/data/08_ModuleCibersortCor_TC.Rdata")
### CB ----
consensusCor_CB = matrix(NA, nrow(moduleTraitCor_CB[[1]]), ncol(moduleTraitCor_CB[[1]]))
consensusPvalue_CB = matrix(NA, nrow(moduleTraitCor_CB[[1]]), ncol(moduleTraitCor_CB[[1]]))
consensusPadj_CB = matrix(NA, nrow(moduleTraitCor_CB[[1]]), ncol(moduleTraitCor_CB[[1]]))

colnames(consensusCor_CB)<-colnames(consensusPvalue_CB)<-colnames(consensusPadj_CB)<-colnames(moduleTraitCor_CB[[1]])
rownames(consensusCor_CB)<-rownames(consensusPvalue_CB)<-rownames(consensusPadj_CB)<-rownames(moduleTraitCor_CB[[1]])

negative = apply(as.data.frame(sapply(moduleTraitCor_CB, function(x) x<0 )),1,function(y){all(y)})


if(length(which(negative==TRUE)) > 1)
{
  consensusCor_CB[negative] = pmax(moduleTraitCor_CB[[1]][negative],moduleTraitCor_CB[[2]][negative],moduleTraitCor_CB[[3]][negative])
  consensusPvalue_CB[negative] = pmax(moduleTraitPvalue_CB[[1]][negative],moduleTraitPvalue_CB[[2]][negative],moduleTraitPvalue_CB[[3]][negative])
  consensusPadj_CB[negative] = pmax(moduleTraitPadj_CB[[1]][negative],moduleTraitPadj_CB[[2]][negative],moduleTraitPadj_CB[[3]][negative])
}

# Find consensus positive correlations
positive = apply(as.data.frame(sapply(moduleTraitCor_CB, function(x) x>0 )),1,function(y){all(y)})

if(length(which(positive==TRUE)) > 1)
{
  consensusCor_CB[positive] = pmin(moduleTraitCor_CB[[1]][positive],moduleTraitCor_CB[[2]][positive],moduleTraitCor_CB[[3]][positive])
  consensusPvalue_CB[positive] = pmax(moduleTraitPvalue_CB[[1]][positive],moduleTraitPvalue_CB[[2]][positive],moduleTraitPvalue_CB[[3]][positive])
  consensusPadj_CB[positive] = pmax(moduleTraitPadj_CB[[1]][positive],moduleTraitPadj_CB[[2]][positive],moduleTraitPadj_CB[[3]][positive])
}

### HIP ----
consensusCor_HIP = matrix(NA, nrow(moduleTraitCor_HIP[[1]]), ncol(moduleTraitCor_HIP[[1]]))
consensusPvalue_HIP = matrix(NA, nrow(moduleTraitCor_HIP[[1]]), ncol(moduleTraitCor_HIP[[1]]))
consensusPadj_HIP = matrix(NA, nrow(moduleTraitCor_HIP[[1]]), ncol(moduleTraitCor_HIP[[1]]))
colnames(consensusCor_HIP)<-colnames(consensusPvalue_HIP)<-colnames(consensusPadj_HIP)<-colnames(moduleTraitCor_HIP[[1]])
rownames(consensusCor_HIP)<-rownames(consensusPvalue_HIP)<-rownames(consensusPadj_HIP)<-rownames(moduleTraitCor_HIP[[1]])

negative = apply(as.data.frame(sapply(moduleTraitCor_HIP, function(x) x<0 )),1,function(y){all(y)})
negative[is.na(negative)]<-FALSE

if(length(which(negative==TRUE)) > 1)
{
  consensusCor_HIP[negative] = pmax(moduleTraitCor_HIP[[1]][negative],moduleTraitCor_HIP[[2]][negative],moduleTraitCor_HIP[[3]][negative],moduleTraitCor_HIP[[4]][negative])
  consensusPvalue_HIP[negative] = pmax(moduleTraitPvalue_HIP[[1]][negative],moduleTraitPvalue_HIP[[2]][negative],moduleTraitPvalue_HIP[[3]][negative],moduleTraitPvalue_HIP[[4]][negative])
  consensusPadj_HIP[negative] = pmax(moduleTraitPadj_HIP[[1]][negative],moduleTraitPadj_HIP[[2]][negative],moduleTraitPadj_HIP[[3]][negative],moduleTraitPadj_HIP[[4]][negative])
}

# Find consensus positive correlations
positive = apply(as.data.frame(sapply(moduleTraitCor_HIP, function(x) x>0 )),1,function(y){all(y)})
positive[is.na(positive)]<-FALSE

if(length(which(positive==TRUE)) > 1)
{
  consensusCor_HIP[positive] = pmin(moduleTraitCor_HIP[[1]][positive],moduleTraitCor_HIP[[2]][positive],moduleTraitCor_HIP[[3]][positive],moduleTraitCor_HIP[[4]][positive])
  consensusPvalue_HIP[positive] = pmax(moduleTraitPvalue_HIP[[1]][positive],moduleTraitPvalue_HIP[[2]][positive],moduleTraitPvalue_HIP[[3]][positive],moduleTraitPvalue_HIP[[4]][positive])
  consensusPadj_HIP[positive] = pmax(moduleTraitPadj_HIP[[1]][positive],moduleTraitPadj_HIP[[2]][positive],moduleTraitPadj_HIP[[3]][positive],moduleTraitPadj_HIP[[4]][positive])
}

### PC ----
consensusCor_PC = matrix(NA, nrow(moduleTraitCor_PC[[1]]), ncol(moduleTraitCor_PC[[1]]))
consensusPvalue_PC = matrix(NA, nrow(moduleTraitCor_PC[[1]]), ncol(moduleTraitCor_PC[[1]]))
consensusPadj_PC = matrix(NA, nrow(moduleTraitCor_PC[[1]]), ncol(moduleTraitCor_PC[[1]]))
colnames(consensusCor_PC)<-colnames(consensusPvalue_PC)<-colnames(consensusPadj_PC)<-colnames(moduleTraitCor_PC[[1]])
rownames(consensusCor_PC)<-rownames(consensusPvalue_PC)<-rownames(consensusPadj_PC)<-rownames(moduleTraitCor_PC[[1]])

negative = apply(as.data.frame(sapply(moduleTraitCor_PC, function(x) x<0 )),1,function(y){all(y)})
negative[is.na(negative)]<-FALSE

if(length(which(negative==TRUE)) > 1)
{
  consensusCor_PC[negative] = pmax(moduleTraitCor_PC[[1]][negative],moduleTraitCor_PC[[2]][negative],moduleTraitCor_PC[[3]][negative],moduleTraitCor_PC[[4]][negative])
  consensusPvalue_PC[negative] = pmax(moduleTraitPvalue_PC[[1]][negative],moduleTraitPvalue_PC[[2]][negative],moduleTraitPvalue_PC[[3]][negative],moduleTraitPvalue_PC[[4]][negative])
  consensusPadj_PC[negative] = pmax(moduleTraitPadj_PC[[1]][negative],moduleTraitPadj_PC[[2]][negative],moduleTraitPadj_PC[[3]][negative],moduleTraitPadj_PC[[4]][negative])
}

# Find consensus positive correlations
positive = apply(as.data.frame(sapply(moduleTraitCor_PC, function(x) x>0 )),1,function(y){all(y)})
positive[is.na(positive)]<-FALSE

if(length(which(positive==TRUE)) > 1)
{
  consensusCor_PC[positive] = pmin(moduleTraitCor_PC[[1]][positive],moduleTraitCor_PC[[2]][positive],moduleTraitCor_PC[[3]][positive],moduleTraitCor_PC[[4]][positive])
  consensusPvalue_PC[positive] = pmax(moduleTraitPvalue_PC[[1]][positive],moduleTraitPvalue_PC[[2]][positive],moduleTraitPvalue_PC[[3]][positive],moduleTraitPvalue_PC[[4]][positive])
  consensusPadj_PC[positive] = pmax(moduleTraitPadj_PC[[1]][positive],moduleTraitPadj_PC[[2]][positive],moduleTraitPadj_PC[[3]][positive],moduleTraitPadj_PC[[4]][positive])
}

### TC ----
consensusCor_TC = matrix(NA, nrow(moduleTraitCor_TC[[1]]), ncol(moduleTraitCor_TC[[1]]))
consensusPvalue_TC = matrix(NA, nrow(moduleTraitCor_TC[[1]]), ncol(moduleTraitCor_TC[[1]]))
consensusPadj_TC = matrix(NA, nrow(moduleTraitCor_TC[[1]]), ncol(moduleTraitCor_TC[[1]]))
colnames(consensusCor_TC)<-colnames(consensusPvalue_TC)<-colnames(consensusPadj_TC)<-colnames(moduleTraitCor_TC[[1]])
rownames(consensusCor_TC)<-rownames(consensusPvalue_TC)<-rownames(consensusPadj_TC)<-rownames(moduleTraitCor_TC[[1]])

negative = apply(as.data.frame(sapply(moduleTraitCor_TC, function(x) x<0 )),1,function(y){all(y)})
negative[is.na(negative)]<-FALSE

if(length(which(negative==TRUE)) > 1)
{
  consensusCor_TC[negative] = pmax(moduleTraitCor_TC[[1]][negative],moduleTraitCor_TC[[2]][negative],moduleTraitCor_TC[[3]][negative])
  consensusPvalue_TC[negative] = pmax(moduleTraitPvalue_TC[[1]][negative],moduleTraitPvalue_TC[[2]][negative],moduleTraitPvalue_TC[[3]][negative])
  consensusPadj_TC[negative] = pmax(moduleTraitPadj_TC[[1]][negative],moduleTraitPadj_TC[[2]][negative],moduleTraitPadj_TC[[3]][negative])
}

# Find consensus positive correlations
positive = apply(as.data.frame(sapply(moduleTraitCor_TC, function(x) x>0 )),1,function(y){all(y)})
positive[is.na(positive)]<-FALSE

if(length(which(positive==TRUE)) > 1)
{
  consensusCor_TC[positive] = pmin(moduleTraitCor_TC[[1]][positive],moduleTraitCor_TC[[2]][positive],moduleTraitCor_TC[[3]][positive])
  consensusPvalue_TC[positive] = pmax(moduleTraitPvalue_TC[[1]][positive],moduleTraitPvalue_TC[[2]][positive],moduleTraitPvalue_TC[[3]][positive])
  consensusPadj_TC[positive] = pmax(moduleTraitPadj_TC[[1]][positive],moduleTraitPadj_TC[[2]][positive],moduleTraitPadj_TC[[3]][positive])
}

## retidy data ----

rbind(cbind(region="CB",module=rownames(consensusCor_CB),consensusCor_CB),
      cbind(region="HIP",module=rownames(consensusCor_HIP),consensusCor_HIP),
      cbind(region="PC",module=rownames(consensusCor_PC),consensusCor_PC),
      cbind(region="TC",module=rownames(consensusCor_TC),consensusCor_TC),
      cbind(region="OC",module=rownames(moduleTraitCor_OC),moduleTraitCor_OC)) %>% as.data.frame()-> moduleTraitCor_total

rbind(cbind(region="CB",module=rownames(consensusPvalue_CB),consensusPvalue_CB),
      cbind(region="HIP",module=rownames(consensusPvalue_HIP),consensusPvalue_HIP),
      cbind(region="PC",module=rownames(consensusPvalue_PC),consensusPvalue_PC),
      cbind(region="TC",module=rownames(consensusPvalue_TC),consensusPvalue_TC),
      cbind(region="OC",module=rownames(moduleTraitPvalue_OC),moduleTraitPvalue_OC)) %>% as.data.frame()-> moduleTraitPvalue_total

rbind(cbind(region="CB",module=rownames(consensusPadj_CB),consensusPadj_CB),
      cbind(region="HIP",module=rownames(consensusPadj_HIP),consensusPadj_HIP),
      cbind(region="PC",module=rownames(consensusPadj_PC),consensusPadj_PC),
      cbind(region="TC",module=rownames(consensusPadj_TC),consensusPadj_TC),
      cbind(region="OC",module=rownames(moduleTraitPvalue_OC),moduleTraitPadj_OC)) %>% as.data.frame()-> moduleTraitPadj_total

save(moduleTraitCor_total,moduleTraitPvalue_total,moduleTraitPadj_total,
        file = "2.output/4.ndd_wgcna/data/08_ModuleCibersortCor_Total.Rdata")
