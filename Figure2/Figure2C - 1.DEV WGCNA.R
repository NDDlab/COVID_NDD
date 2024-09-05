library(WGCNA)
library(tidyverse)
library(data.table)

get_expmat <- function(myfiles,type){
  
  myfiles <- unique(gsub("_.*","",myfiles))
  
  all_files <- list.files(paste0("1.rawdata/expmat_QC/",type),pattern = "txt",full.names = T)
  
  sapply(
    myfiles,
    function(x){
      idx = str_detect(
        toupper(all_files),
        toupper(x)
      )
      idx = all_files[idx]#which(idx==TRUE)
      return(idx)
    }
  ) -> purpose_files
  
  expmat <- map(purpose_files,fread,sep="\t") %>% reduce(full_join,"gene")
  
  expmat[is.na(expmat)] <- 0
  
  idx_80<-apply(expmat[,-1], 1, function(x){ length(x[as.numeric(x)==0]) / length(x) < 0.2}) # 在80%的样本里面都>0
  
  expmat_filt<-expmat[idx_80,]
  
  return(expmat_filt)
}

read_sample.ndd<-function(){
  sample.disease<-fread("./1.rawdata/Disease_sample_info.txt")
  return(sample.disease)
}

read_sample.dev<-function(){
  sample.develop<-fread("./1.rawdata/Develop_sample_info_20220816.txt")
  return(sample.develop)
}

######################################################
#                     1. import sample and Development(DEV) expression data
######################################################
sample <- read_sample.dev()
sample %<>% 
  dplyr::mutate(GSE = case_when(GSE=="Brainspan" ~ paste0(GSE,"-",GPL), TRUE ~ GSE))

samplelist_brain <- sample %>% split(.$brain_region2)

names(samplelist_brain) <- gsub("\\(.*","",gsub(".*/","",names(samplelist_brain)))

gselist_brain <- map(samplelist_brain,function(x){unique(x$GSE)})

dtlist_brain <- lapply(gselist_brain,get_expmat,type="dev")
names(dtlist_brain) <- gsub(".*/","",names(dtlist_brain))

for (i in 1:length(dtlist_brain)) {
  
  name_i = names(dtlist_brain)[i]
  sample_i = sample$GSM[str_detect(sample$brain_region2,name_i)]
  cat( paste0(length(sample_i), "/", ncol(dtlist_brain[[i]])-1) )
  dtlist_brain[[i]] <- dtlist_brain[[i]] %>% dplyr::select(gene,all_of(sample_i))
  
  cat(" -> ", ncol(dtlist_brain[[i]])-1,"\n")
}

save(dtlist_brain,samplelist_brain,file="./2.output/9.dev_wgcna/data/1_dev_ExpmatAndSample_of_each_region.Rdata")

######################################################
#                     2. WGCNA data preparation
# batch correction -> filter genes with mad above  0 -> detect sample outlier
# FCx -> TCx -> OCx -> PCx -> CB -> HIP -> STR/VZ -> AMY
######################################################
get_gooddt_idx <- function(dtlist_FC){
  
  idx_filt <- c()
  
  ### filter the sample
  for (dt in 1:length(dtlist_FC)){
    
    window = dtlist_FC[[dt]] %>% group_by(age_stage2) %>% count()
    cat(names(dtlist_FC)[dt])
    print(window)
    
    # condition1 : filter the datasets with all windows
    if(length(window$age_stage2)==6){
      # condition2 : filter the datasets with above 15 samples
      if(sum(window$n) >= 15){
        # condition3 : filter the datasets with each window above 5 samples
        window_filt = window$age_stage2[window$n>=3]
        if(length(window_filt)==6){
          idx_filt <- append(idx_filt,dt)
        }
      }  
    }
  }
  return(idx_filt)
}

batch_correction <- function(trait,data){
  
  batch = trait$GSE
  batch2 = trait$GPL
  # considering window
  mod = model.matrix(~trait$age_stage2)
  
  # if exists batch correction
  if(length(unique(batch))>1 & length(unique(batch2))>1)
  {
    dt1 = as.data.frame(ComBat(data,batch))
    as.data.frame(ComBat(dt1,batch2,mod))
  } else if(length(unique(batch))>1 & length(unique(batch2))==1){
    as.data.frame(ComBat(data,batch,mod))
  } else if(length(unique(batch))==1 & length(unique(batch2))>1){
    as.data.frame(ComBat(data,batch2,mod))
  } else {
    data
  } 
  
}

load(file="./2.output/9.dev_wgcna/data/1_dev_ExpmatAndSample_of_each_region.Rdata")

batch_correction <- function(trait,data){
  
  batch = trait$GSE
  batch2 = trait$GPL
  # considering window
  mod = model.matrix(~trait$age_stage2)
  
  # if exists batch correction
  if(length(unique(batch))>1 & length(unique(batch2))>1)
  {
    dt1 = as.data.frame(ComBat(data,batch))
    as.data.frame(ComBat(dt1,batch2,mod))
  } else if(length(unique(batch))>1 & length(unique(batch2))==1){
    as.data.frame(ComBat(data,batch,mod))
  } else if(length(unique(batch))==1 & length(unique(batch2))>1){
    as.data.frame(ComBat(data,batch2,mod))
  } else {
    data
  } 
}

## 2.1 FCx -----
### 2.1.1 preprocess data ----
### get dt
dt_FC = dtlist_brain$FCx %>% column_to_rownames("gene")
sample_FC = samplelist_brain$FCx

### split by gse
dtlist_FC <- sample_FC %>% split(.$GSE)

### filter the effective
idx_filt <- get_gooddt_idx(dtlist_FC)
dtlist_FC = dtlist_FC[idx_filt]
dtinfo_FC <- names(dtlist_FC)

brainspanArray_FC <- dt_FC  %>% dplyr::select(all_of(dtlist_FC$`Brainspan-array`$GSM))
brainspanSeq_FC <- dt_FC  %>% dplyr::select(all_of(dtlist_FC$`Brainspan-rnaseq`$GSM))
GSE25219_FC <- dt_FC  %>% dplyr::select(all_of(dtlist_FC$GSE25219$GSM))
PsychENCODE_FC <- dt_FC  %>% dplyr::select(all_of(dtlist_FC$PsychENCODE$GSM))

idx_del_FC <- unique(c(
  which(apply(brainspanArray_FC, 1, mad)==0),
  which(apply(brainspanSeq_FC, 1, mad)==0),
  which(apply(GSE25219_FC, 1, mad)==0),
  which(apply(PsychENCODE_FC, 1, mad)==0)
))

if(length(idx_del_FC)>0){
  brainspanArray_FC <- brainspanArray_FC[-idx_del_FC, ] 
  brainspanSeq_FC <- brainspanSeq_FC[-idx_del_FC, ] 
  GSE25219_FC <- GSE25219_FC[-idx_del_FC, ] 
  PsychENCODE_FC <- PsychENCODE_FC[-idx_del_FC, ] 
}

### 2.1.2 Form data ----
nSets_FC = 4
setLabels_FC = c("BrainspanArray","BrainspanRNAseq","GSE25219","PsychENCODE")

multiExpr_FC = vector(mode = "list", length = nSets_FC)
multiExpr_FC[[1]] = list(data = as.data.frame(t(brainspanArray_FC)))
multiExpr_FC[[2]] = list(data = as.data.frame(t(brainspanSeq_FC)))
multiExpr_FC[[3]] = list(data = as.data.frame(t(GSE25219_FC)))
multiExpr_FC[[4]] = list(data = as.data.frame(t(PsychENCODE_FC)))

# Check that the data has the correct format for many functions operating on multiple sets:
exprSize_FC = checkSets(multiExpr_FC)
exprSize_FC

### 2.1.3 Rudimentary data cleaning and outlier removal ----
# Check that all genes and samples have sufficiently low numbers of missing values.
gsg <- goodSamplesGenes(dt_FC,verbose = 3)
gsg$allOK

### 2.1.4 sample clustering ----
sampleTrees_FC = list()

for (set in 1:nSets_FC){
  sampleTrees_FC[[set]] = hclust(dist(multiExpr_FC[[set]]$data), method = "average")
}

# output
for (set in 1:nSets_FC) {
  pdf(file = paste0("2.output/4.dev_wgcna/figure/FCx/dtsplit/02_SampleClustering_",setLabels_FC[set],".pdf"), width = 14, height = 7)
  plot(sampleTrees_FC[[set]], xlab="", sub="", cex = 0.7,
       main = paste("Sample clustering on all genes in", setLabels_FC[[set]])) 
  dev.off()
}

# set cutheight
setLabels_FC
cutHeights = c(0,90,40,55,90)
# plot cluster tree
pdf(file = "2.output/2.ndd_wgcna/figure/FCx/02_SampleClustering_cutHeights.pdf", width = 12, height = 12)
par(mfrow=c(3,2))
par(mar = c(0, 4, 2, 0))
for (set in 1:nSets_FC){
  plot(sampleTrees_FC[[set]], xlab="", sub="", cex = 0.7,
       main = paste("Sample clustering on all genes in", setLabels_FC[set])) 
  abline(h=cutHeights[set], col = "red")
}
dev.off()

# Now comes the actual outlier removal:
# Find clusters cut by the line
# Keep the largest one (labeled by the number 1)
setLabels_FC
# BD
labels = cutreeStatic(sampleTrees_FC[[2]],cutHeight = cutHeights[2],minSize = 10)
table(labels)
keep = (labels==1)
multiExpr_FC[[2]]$data = multiExpr_FC[[2]]$data[keep, ]

collectGarbage()
# Check the size of the leftover data
exprSize_FC = checkSets(multiExpr_FC)
exprSize_FC

### 2.1.5. import clinical information -----
sample_FC %<>% 
  dplyr::mutate(HW1 = ifelse(age_stage2=="HW1",1,0),
                HW2 = ifelse(age_stage2=="HW2",1,0),
                HW3 = ifelse(age_stage2=="HW3",1,0),
                HW4 = ifelse(age_stage2=="HW4",1,0),
                HW5 = ifelse(age_stage2=="HW5",1,0),
                HW6 = ifelse(age_stage2=="HW6",1,0)) %>% 
  dplyr::select(GSM,starts_with("HW")) %>% 
  unique() %>% 
  column_to_rownames("GSM")

Traits_FC = vector(mode="list", length = nSets_FC)
for (set in 1:nSets_FC)
{
  setSamples = rownames(multiExpr_FC[[set]]$data)
  traitRows = match(setSamples, rownames(sample_FC))
  Traits_FC[[set]] = list(data = sample_FC[traitRows, ]) 
}

# Define data set dimensions
nGenes_FC = exprSize_FC$nGenes
nSamples_FC = exprSize_FC$nSamples

### 2.1.6. reclustering samples with trait ----
sampleTrees2_FC = list()

for (set in 1:nSets_FC){
  sampleTrees2_FC[[set]] = hclust(dist(multiExpr_FC[[set]]$data), method = "average")
}

#pdf(file = "2.output/2.ndd_wgcna/figure/FCx/02_SampleClustering_cutHeights.pdf", width = 12, height = 12)
par(mar = c(0, 4, 2, 0))
for (set in 1:nSets_FC) {
  pdf(file = paste0("2.output/4.dev_wgcna/figure/FCx/dtsplit/03_SampleClustering_Trait_",setLabels_FC[set],".pdf"), width = 10, height = 5)
  traitColors = numbers2colors(Traits_FC[[set]]$data, signed = FALSE)
  plotDendroAndColors(sampleTrees2_FC[[set]], 
                      traitColors,
                      groupLabels = colnames(Traits_FC[[1]]$data),
                      autoColorHeight = T,
                      colorHeight = 0.02,
                      colorHeightBase = 0.02,
                      main = paste0("Sample dendrogram and trait heatmap-",setLabels_FC[set]))
  dev.off()
}

pdf(file = paste0("2.output/4.dev_wgcna/figure/FCx/dtsplit/03_SampleClustering_Trait_Total.pdf"), width = 10, height = 5)
for (set in 1:nSets_FC) {
  traitColors = numbers2colors(Traits_FC[[set]]$data, signed = FALSE)
  plotDendroAndColors(sampleTrees2_FC[[set]], traitColors,
                      groupLabels = "Status", 
                      main = paste0("Sample dendrogram and trait heatmap-",setLabels_FC[set]))
}
dev.off()

### 2.1.7. export exp, trait, and other info ----
save(multiExpr_FC, Traits_FC, nGenes_FC, nSamples_FC, setLabels_FC, exprSize_FC,
     file = "./2.output/4.dev_wgcna/data/02_Consensus_dataInput_FC.RData")

## 2.2 TCx----
### 2.2.1 preprocess data ----
### get dt
dt_TC = dtlist_brain$TCx %>% column_to_rownames("gene")
sample_TC = samplelist_brain$TCx

dtlist_TC <- sample_TC %>% split(.$GSE)

### filter the effective
idx_filt <- get_gooddt_idx(dtlist_TC)
dtlist_TC = dtlist_TC[idx_filt]
dtinfo_TC <- names(dtlist_TC)

brainspanArray_TC <- dt_TC  %>% dplyr::select(all_of(dtlist_TC$`Brainspan-array`$GSM))
brainspanSeq_TC <- dt_TC  %>% dplyr::select(all_of(dtlist_TC$`Brainspan-rnaseq`$GSM))
GSE25219_TC <- dt_TC  %>% dplyr::select(all_of(dtlist_TC$GSE25219$GSM))
PsychENCODE_TC <- dt_TC  %>% dplyr::select(all_of(dtlist_TC$PsychENCODE$GSM))

idx_del_TC <- unique(c(
  which(apply(brainspanArray_TC, 1, mad)==0),
  which(apply(brainspanSeq_TC, 1, mad)==0),
  which(apply(GSE25219_TC, 1, mad)==0),
  which(apply(PsychENCODE_TC, 1, mad)==0)
))

if(length(idx_del_TC)>0){
  brainspanArray_TC <- brainspanArray_TC[-idx_del_TC, ] 
  brainspanSeq_TC <- brainspanSeq_TC[-idx_del_TC, ] 
  GSE25219_TC <- GSE25219_TC[-idx_del_TC, ] 
  PsychENCODE_TC <- PsychENCODE_TC[-idx_del_TC, ] 
}

### 2.2.2 Form data ----
nSets_TC = 4
setLabels_TC = c("BrainspanArray","BrainspanRNAseq","GSE25219","PsychENCODE")

multiExpr_TC = vector(mode = "list", length = nSets_TC)
multiExpr_TC[[1]] = list(data = as.data.frame(t(brainspanArray_TC)))
multiExpr_TC[[2]] = list(data = as.data.frame(t(brainspanSeq_TC)))
multiExpr_TC[[3]] = list(data = as.data.frame(t(GSE25219_TC)))
multiExpr_TC[[4]] = list(data = as.data.frame(t(PsychENCODE_TC)))

# Check that the data has the correct format for many functions operating on multiple sets:
exprSize_TC = checkSets(multiExpr_TC)
exprSize_TC

### 2.2.3 Rudimentary data cleaning and outlier removal ----
# Check that all genes and samples have sufficiently low numbers of missing values.
gsg <- goodSamplesGenes(dt_TC,verbose = 3)
gsg$allOK

### 2.2.4 sample clustering ----
sampleTrees_TC = list()

for (set in 1:nSets_TC){
  sampleTrees_TC[[set]] = hclust(dist(multiExpr_TC[[set]]$data), method = "average")
}

# output
for (set in 1:nSets_TC) {
  pdf(file = paste0("2.output/9.dev_wgcna/figure/TCx/dtsplit/02_SampleClustering_",setLabels_TC[set],".pdf"), width = 14, height = 7)
  plot(sampleTrees_TC[[set]], xlab="", sub="", cex = 0.7,
       main = paste("Sample clustering on all genes in", setLabels_TC[[set]])) 
  dev.off()
}

# set cutheight
setLabels_TC
cutHeights = c(100,90,0,110)
# plot cluster tree
pdf(file = "2.output/9.dev_wgcna/figure/TCx/dtsplit/02_SampleClustering_cutHeights.pdf", width = 12, height = 12)
par(mfrow=c(3,2))
par(mar = c(0, 4, 2, 0))
for (set in 1:nSets_TC){
  plot(sampleTrees_TC[[set]], xlab="", sub="", cex = 0.7,
       main = paste("Sample clustering on all genes in", setLabels_TC[set])) 
  abline(h=cutHeights[set], col = "red")
}
dev.off()

# Now comes the actual outlier removal:
# Find clusters cut by the line
# Keep the largest one (labeled by the number 1)
setLabels_TC
# BrainspanArray
labels = cutreeStatic(sampleTrees_TC[[1]],cutHeight = cutHeights[1],minSize = 10)
table(labels)
keep = (labels==1)
multiExpr_TC[[1]]$data = multiExpr_TC[[1]]$data[keep, ]
# BrainspanRNAseq
labels = cutreeStatic(sampleTrees_TC[[2]],cutHeight = cutHeights[2],minSize = 10)
table(labels)
keep = (labels==1)
multiExpr_TC[[2]]$data = multiExpr_TC[[2]]$data[keep, ]
# PsychENCODE
labels = cutreeStatic(sampleTrees_TC[[4]],cutHeight = cutHeights[4],minSize = 10)
table(labels)
keep = (labels==1)
multiExpr_TC[[4]]$data = multiExpr_TC[[4]]$data[keep, ]

# Check the size of the leftover data
exprSize_TC = checkSets(multiExpr_TC)
exprSize_TC

### 2.2.5. import clinical information -----
unique(sample_TC$age_stage2)
sample_TC %<>% 
  dplyr::mutate(HW1 = ifelse(age_stage2=="HW1",1,0),
                HW2 = ifelse(age_stage2=="HW2",1,0),
                HW3 = ifelse(age_stage2=="HW3",1,0),
                HW4 = ifelse(age_stage2=="HW4",1,0),
                HW5 = ifelse(age_stage2=="HW5",1,0),
                HW6 = ifelse(age_stage2=="HW6",1,0)) %>% 
  dplyr::select(GSM,starts_with("HW")) %>% 
  unique() %>% 
  column_to_rownames("GSM")

Traits_TC = vector(mode="list", length = nSets_TC)
for (set in 1:nSets_TC)
{
  setSamples = rownames(multiExpr_TC[[set]]$data)
  traitRows = match(setSamples, rownames(sample_TC))
  Traits_TC[[set]] = list(data = sample_TC[traitRows, ]) 
}

# Define data set dimensions
nGenes_TC = exprSize_TC$nGenes
nSamples_TC = exprSize_TC$nSamples

### 2.2.6. reclustering samples with trait ----
sampleTrees2_TC = list()

for (set in 1:nSets_TC){
  sampleTrees2_TC[[set]] = hclust(dist(multiExpr_TC[[set]]$data), method = "average")
}

par(mar = c(0, 4, 2, 0))
for (set in 1:nSets_TC) {
  pdf(file = paste0("2.output/9.dev_wgcna/figure/TCx/dtsplit/03_SampleClustering_Trait_",setLabels_TC[set],".pdf"), width = 10, height = 5)
  traitColors = numbers2colors(Traits_TC[[set]]$data, signed = FALSE)
  plotDendroAndColors(sampleTrees2_TC[[set]], 
                      traitColors,
                      groupLabels = colnames(Traits_TC[[1]]$data),
                      autoColorHeight = T,
                      colorHeight = 0.02,
                      colorHeightBase = 0.02,
                      main = paste0("Sample dendrogram and trait heatmap-",setLabels_TC[set]))
  dev.off()
}

pdf(file = paste0("2.output/9.dev_wgcna/figure/TCx/dtsplit/03_SampleClustering_Trait_Total.pdf"), width = 10, height = 5)
for (set in 1:nSets_TC) {
  traitColors = numbers2colors(Traits_TC[[set]]$data, signed = FALSE)
  plotDendroAndColors(sampleTrees2_TC[[set]], traitColors,
                      groupLabels =  colnames(Traits_TC[[1]]$data),
                      main = paste0("Sample dendrogram and trait heatmap-",setLabels_TC[set]))
}
dev.off()

### 2.2.7. export exp, trait, and other info ----
save(multiExpr_TC, Traits_TC, nGenes_TC, nSamples_TC, setLabels_TC, exprSize_TC,
     file = "./2.output/9.dev_wgcna/data/02_Consensus_dataInput_TC.RData")

## 2.3 OCx ----
### 2.3.1 preprocess data ----
### get dt
dt_OC = dtlist_brain$OCx %>% column_to_rownames("gene")
sample_OC = samplelist_brain$OCx

### split by gse
dtlist_OC <- sample_OC %>% split(.$GSE)

### filter the effective
idx_filt <- get_gooddt_idx(dtlist_OC)
dtlist_OC = dtlist_OC[idx_filt]
dtinfo_OC <- names(dtlist_OC)

GSE25219_OC <- dt_OC  %>% dplyr::select(all_of(dtlist_OC$GSE25219$GSM))

idx_del_OC <- which(apply(GSE25219_OC, 1, mad)==0)

if(length(idx_del_OC)>0){
  GSE25219_OC <- GSE25219_OC[-idx_del_OC, ] 
}

### 2.3.2 Form data ----
GSE25219_OC <- as.data.frame(t(GSE25219_OC))

### 2.3.3 Rudimentary data cleaning and outlier removal ----
# Check that all genes and samples have sufficiently low numbers of missing values.
gsg <- goodSamplesGenes(GSE25219_OC,verbose = 3)
gsg$allOK

### 2.3.4 sample clustering ----
sampleTrees_OC = hclust(dist(GSE25219_OC), method = "average")

pdf(file = "2.output/9.dev_wgcna/figure/OCx/dtsplit/02_SampleClustering_OCx_cutHeights.pdf", width = 14, height = 7)
plot(sampleTrees_OC, xlab="", sub="", cex = 0.7,
     main = paste("Sample clustering on all genes in OCx")) 
#abline(h=235, col = "red")
dev.off()

### 2.3.5. import clinical information -----
sample_OC %<>% 
  dplyr::mutate(HW1 = ifelse(age_stage2=="HW1",1,0),
                HW2 = ifelse(age_stage2=="HW2",1,0),
                HW3 = ifelse(age_stage2=="HW3",1,0),
                HW4 = ifelse(age_stage2=="HW4",1,0),
                HW5 = ifelse(age_stage2=="HW5",1,0),
                HW6 = ifelse(age_stage2=="HW6",1,0)) %>% 
  dplyr::select(GSM,starts_with("HW")) %>% 
  unique() %>% 
  column_to_rownames("GSM")

setSamples = rownames(GSE25219_OC)
traitRows = match(setSamples, rownames(sample_OC))
Traits_OC = sample_OC[traitRows, ]

# Define data set dimensions
nGenes_OC = ncol(GSE25219_OC)
nSamples_OC = nrow(GSE25219_OC)

### 2.3.6. reclustering samples with trait ----
sampleTrees2_OC = hclust(dist(GSE25219_OC), method = "average")

traitColors = numbers2colors(Traits_OC, signed = FALSE)

pdf(file = paste0("2.output/9.dev_wgcna/figure/OCx/dtsplit/03_SampleClustering_OCx_Trait.pdf"), width = 10, height = 5)
plotDendroAndColors(sampleTrees2_OC, traitColors,
                    groupLabels = colnames(Traits_OC), 
                    autoColorHeight = T,
                    colorHeight = 0.02,
                    colorHeightBase = 0.02,
                    main = "Sample dendrogram and trait heatmap-OCx")
dev.off()

### 2.3.7. export exp, trait, and other info ----
save(GSE25219_OC, Traits_OC, nGenes_OC, nSamples_OC, 
     file = "./2.output/9.dev_wgcna/data/02_NonConsensus_dataInput_OC_GSE25219.RData")

## 2.4 PCx -----
### 2.4.1 preprocess data ----
### get dt
dt_PC = dtlist_brain$PCx %>% column_to_rownames("gene")
sample_PC = samplelist_brain$PCx

### split by gse
dtlist_PC <- sample_PC %>% split(.$GSE)

### filter the effective
idx_filt <- get_gooddt_idx(dtlist_PC)
dtlist_PC = dtlist_PC[idx_filt]
dtinfo_PC <- names(dtlist_PC)

GSE25219_PC <- dt_PC  %>% dplyr::select(all_of(dtlist_PC$GSE25219$GSM))

idx_del_PC <- which(apply(GSE25219_PC, 1, mad)==0)

if(length(idx_del_PC)>0){
  GSE25219_PC <- GSE25219_PC[-idx_del_PC, ] 
}

### 2.4.2 Form data ----
GSE25219_PC <- as.data.frame(t(GSE25219_PC))

### 2.4.3 Rudimentary data cleaning and outlier removal ----
# Check that all genes and samples have sufficiently low numbers of missing values.
gsg <- goodSamplesGenes(GSE25219_PC,verbose = 3)
gsg$allOK

### 2.4.4 sample clustering ----
sampleTrees_PC = hclust(dist(GSE25219_PC), method = "average")

pdf(file = "2.output/9.dev_wgcna/figure/PCx/dtsplit/02_SampleClustering_PCx_cutHeights.pdf", width = 14, height = 7)
plot(sampleTrees_PC, xlab="", sub="", cex = 0.7,
     main = paste("Sample clustering on all genes in PCx")) 
#abline(h=235, col = "red")
dev.off()

### 2.4.5. import clinical information -----
sample_PC %<>% 
  dplyr::mutate(HW1 = ifelse(age_stage2=="HW1",1,0),
                HW2 = ifelse(age_stage2=="HW2",1,0),
                HW3 = ifelse(age_stage2=="HW3",1,0),
                HW4 = ifelse(age_stage2=="HW4",1,0),
                HW5 = ifelse(age_stage2=="HW5",1,0),
                HW6 = ifelse(age_stage2=="HW6",1,0)) %>% 
  dplyr::select(GSM,starts_with("HW")) %>% 
  unique() %>% 
  column_to_rownames("GSM")

setSamples = rownames(GSE25219_PC)
traitRows = match(setSamples, rownames(sample_PC))
Traits_PC = sample_PC[traitRows, ]

# Define data set dimensions
nGenes_PC = ncol(GSE25219_PC)
nSamples_PC = nrow(GSE25219_PC)

### 2.4.6. reclustering samples with trait ----
sampleTrees2_PC = hclust(dist(GSE25219_PC), method = "average")

traitColors = numbers2colors(Traits_PC, signed = FALSE)

pdf(file = paste0("2.output/9.dev_wgcna/figure/PCx/dtsplit/03_SampleClustering_PCx_Trait.pdf"), width = 10, height = 5)
plotDendroAndColors(sampleTrees2_PC, traitColors,
                    groupLabels = colnames(Traits_PC), 
                    autoColorHeight = T,
                    colorHeight = 0.02,
                    colorHeightBase = 0.02,
                    main = "Sample dendrogram and trait heatmap-PCx")
dev.off()

### 2.4.7. export exp, trait, and other info ----
save(GSE25219_PC, Traits_PC, nGenes_PC, nSamples_PC, 
     file = "./2.output/9.dev_wgcna/data/02_NonConsensus_dataInput_PC_GSE25219.RData")

## 2.5 CB -----
### 2.5.1 preprocess data ----
### get dt
dt_CB = dtlist_brain$URL %>% column_to_rownames("gene")
sample_CB = samplelist_brain$URL

### split by gse
dtlist_CB <- sample_CB %>% split(.$GSE)

### filter the effective
idx_filt <- get_gooddt_idx(dtlist_CB)
dtlist_CB = dtlist_CB[idx_filt]
dtinfo_CB <- names(dtlist_CB)

GSE25219_CB <- dt_CB  %>% dplyr::select(all_of(dtlist_CB$GSE25219$GSM))

idx_del_CB <- which(apply(GSE25219_CB, 1, mad)==0)

if(length(idx_del_CB)>0){
  GSE25219_CB <- GSE25219_CB[-idx_del_CB, ] 
}

### 2.5.2 Form data ----
GSE25219_CB <- as.data.frame(t(GSE25219_CB))

### 2.5.3 Rudimentary data cleaning and outlier removal ----
# Check that all genes and samples have sufficiently low numbers of missing values.
gsg <- goodSamplesGenes(GSE25219_CB,verbose = 3)
gsg$allOK

### 2.5.4 sample clustering ----
sampleTrees_CB = hclust(dist(GSE25219_CB), method = "average")

pdf(file = "2.output/9.dev_wgcna/figure/CB/dtsplit/02_SampleClustering_CB_cutHeights.pdf", width = 14, height = 7)
plot(sampleTrees_CB, xlab="", sub="", cex = 0.7,
     main = paste("Sample clustering on all genes in CB")) 
#abline(h=235, col = "red")
dev.off()

### 2.4.5. import clinical information -----
sample_CB %<>% 
  dplyr::mutate(HW1 = ifelse(age_stage2=="HW1",1,0),
                HW2 = ifelse(age_stage2=="HW2",1,0),
                HW3 = ifelse(age_stage2=="HW3",1,0),
                HW4 = ifelse(age_stage2=="HW4",1,0),
                HW5 = ifelse(age_stage2=="HW5",1,0),
                HW6 = ifelse(age_stage2=="HW6",1,0)) %>% 
  dplyr::select(GSM,starts_with("HW")) %>% 
  unique() %>% 
  column_to_rownames("GSM")

setSamples = rownames(GSE25219_CB)
traitRows = match(setSamples, rownames(sample_CB))
Traits_CB = sample_CB[traitRows, ]

# Define data set dimensions
nGenes_CB = ncol(GSE25219_CB)
nSamples_CB = nrow(GSE25219_CB)

### 2.4.6. reclustering samples with trait ----
sampleTrees2_CB = hclust(dist(GSE25219_CB), method = "average")

traitColors = numbers2colors(Traits_CB, signed = FALSE)

pdf(file = paste0("2.output/9.dev_wgcna/figure/CB/dtsplit/03_SampleClustering_CB_Trait.pdf"), width = 10, height = 5)
plotDendroAndColors(sampleTrees2_CB, traitColors,
                    groupLabels = colnames(Traits_CB), 
                    autoColorHeight = T,
                    colorHeight = 0.02,
                    colorHeightBase = 0.02,
                    main = "Sample dendrogram and trait heatmap-CB")
dev.off()

### 2.4.7. export exp, trait, and other info ----
save(GSE25219_CB, Traits_CB, nGenes_CB, nSamples_CB, 
     file = "./2.output/9.dev_wgcna/data/02_NonConsensus_dataInput_CB_GSE25219.RData")

## 2.6 HIP -----
### 2.6.1 preprocess data ----
### get dt
dt_HIP = dtlist_brain$HIP %>% column_to_rownames("gene")
sample_HIP = samplelist_brain$HIP

### split by gse
dtlist_HIP <- sample_HIP %>% split(.$GSE)

### filter the effective
idx_filt <- get_gooddt_idx(dtlist_HIP)
dtlist_HIP = dtlist_HIP[idx_filt]
dtinfo_HIP <- names(dtlist_HIP)

GSE25219_HIP <- dt_HIP  %>% dplyr::select(all_of(dtlist_HIP$GSE25219$GSM))

idx_del_HIP <- which(apply(GSE25219_HIP, 1, mad)==0)

if(length(idx_del_HIP)>0){
  GSE25219_HIP <- GSE25219_HIP[-idx_del_HIP, ] 
}

### 2.6.2 Form data ----
GSE25219_HIP <- as.data.frame(t(GSE25219_HIP))

### 2.6.3 Rudimentary data cleaning and outlier removal ----
# Check that all genes and samples have sufficiently low numbers of missing values.
gsg <- goodSamplesGenes(GSE25219_HIP,verbose = 3)
gsg$allOK

### 2.6.4 sample clustering ----
sampleTrees_HIP = hclust(dist(GSE25219_HIP), method = "average")

pdf(file = "2.output/9.dev_wgcna/figure/HIP/dtsplit/02_SampleClustering_HIP_cutHeights.pdf", width = 14, height = 7)
plot(sampleTrees_HIP, xlab="", sub="", cex = 0.7,
     main = paste("Sample clustering on all genes in HIP")) 
#abline(h=235, col = "red")
dev.off()

### 2.6.5. import clinical information -----
sample_HIP %<>% 
  dplyr::mutate(HW1 = ifelse(age_stage2=="HW1",1,0),
                HW2 = ifelse(age_stage2=="HW2",1,0),
                HW3 = ifelse(age_stage2=="HW3",1,0),
                HW4 = ifelse(age_stage2=="HW4",1,0),
                HW5 = ifelse(age_stage2=="HW5",1,0),
                HW6 = ifelse(age_stage2=="HW6",1,0)) %>% 
  dplyr::select(GSM,starts_with("HW")) %>% 
  unique() %>% 
  column_to_rownames("GSM")

setSamples = rownames(GSE25219_HIP)
traitRows = match(setSamples, rownames(sample_HIP))
Traits_HIP = sample_HIP[traitRows, ]

# Define data set dimensions
nGenes_HIP = ncol(GSE25219_HIP)
nSamples_HIP = nrow(GSE25219_HIP)

### 2.6.6. reclustering samples with trait ----
sampleTrees2_HIP = hclust(dist(GSE25219_HIP), method = "average")

traitColors = numbers2colors(Traits_HIP, signed = FALSE)

pdf(file = paste0("2.output/9.dev_wgcna/figure/HIP/dtsplit/03_SampleClustering_HIP_Trait.pdf"), width = 10, height = 5)
plotDendroAndColors(sampleTrees2_HIP, traitColors,
                    groupLabels = colnames(Traits_HIP), 
                    autoColorHeight = T,
                    colorHeight = 0.02,
                    colorHeightBase = 0.02,
                    main = "Sample dendrogram and trait heatmap-HIP")
dev.off()

### 2.6.7. export exp, trait, and other info ----
save(GSE25219_HIP, Traits_HIP, nGenes_HIP, nSamples_HIP, 
     file = "./2.output/9.dev_wgcna/data/02_NonConsensus_dataInput_HIP_GSE25219.RData")

## 2.7 STR -----
### 2.7.1 preprocess data ----
### get dt
dt_STR = dtlist_brain$VZ %>% column_to_rownames("gene")
sample_STR = samplelist_brain$VZ

### split by gse
dtlist_STR <- sample_STR %>% split(.$GSE)

### filter the effective
idx_filt <- get_gooddt_idx(dtlist_STR)
dtlist_STR = dtlist_STR[idx_filt]
dtinfo_STR <- names(dtlist_STR)

GSE25219_STR <- dt_STR  %>% dplyr::select(all_of(dtlist_STR$GSE25219$GSM))

idx_del_STR <- which(apply(GSE25219_STR, 1, mad)==0)

if(length(idx_del_STR)>0){
  GSE25219_STR <- GSE25219_STR[-idx_del_STR, ] 
}

### 2.7.2 Form data ----
GSE25219_STR <- as.data.frame(t(GSE25219_STR))

### 2.7.3 Rudimentary data cleaning and outlier removal ----
# Check that all genes and samples have sufficiently low numbers of missing values.
gsg <- goodSamplesGenes(GSE25219_STR,verbose = 3)
gsg$allOK

### 2.7.4 sample clustering ----
sampleTrees_STR = hclust(dist(GSE25219_STR), method = "average")

pdf(file = "2.output/9.dev_wgcna/figure/STR/dtsplit/02_SampleClustering_STR_cutHeights.pdf", width = 14, height = 7)
plot(sampleTrees_STR, xlab="", sub="", cex = 0.7,
     main = paste("Sample clustering on all genes in STR")) 
#abline(h=235, col = "red")
dev.off()

### 2.7.5. import clinical information -----
sample_STR %<>% 
  dplyr::mutate(HW1 = ifelse(age_stage2=="HW1",1,0),
                HW2 = ifelse(age_stage2=="HW2",1,0),
                HW3 = ifelse(age_stage2=="HW3",1,0),
                HW4 = ifelse(age_stage2=="HW4",1,0),
                HW5 = ifelse(age_stage2=="HW5",1,0),
                HW6 = ifelse(age_stage2=="HW6",1,0)) %>% 
  dplyr::select(GSM,starts_with("HW")) %>% 
  unique() %>% 
  column_to_rownames("GSM")

setSamples = rownames(GSE25219_STR)
traitRows = match(setSamples, rownames(sample_STR))
Traits_STR = sample_STR[traitRows, ]

# Define data set dimensions
nGenes_STR = ncol(GSE25219_STR)
nSamples_STR = nrow(GSE25219_STR)

### 2.7.6. reclustering samples with trait ----
sampleTrees2_STR = hclust(dist(GSE25219_STR), method = "average")

traitColors = numbers2colors(Traits_STR, signed = FALSE)

pdf(file = paste0("2.output/9.dev_wgcna/figure/STR/dtsplit/03_SampleClustering_STR_Trait.pdf"), width = 10, height = 5)
plotDendroAndColors(sampleTrees2_STR, traitColors,
                    groupLabels = colnames(Traits_STR), 
                    autoColorHeight = T,
                    colorHeight = 0.02,
                    colorHeightBase = 0.02,
                    main = "Sample dendrogram and trait heatmap-STR")
dev.off()

### 2.7.7. export exp, trait, and other info ----
save(GSE25219_STR, Traits_STR, nGenes_STR, nSamples_STR, 
     file = "./2.output/9.dev_wgcna/data/02_NonConsensus_dataInput_STR_GSE25219.RData")

## 2.8 AMY -----
### 2.8.1 preprocess data ----
### get dt
dt_AMY = dtlist_brain$AMY %>% column_to_rownames("gene")
sample_AMY = samplelist_brain$AMY

### split by gse
dtlist_AMY <- sample_AMY %>% split(.$GSE)

### filter the effective
idx_filt <- get_gooddt_idx(dtlist_AMY)
dtlist_AMY = dtlist_AMY[idx_filt]
dtinfo_AMY <- names(dtlist_AMY)

GSE25219_AMY <- dt_AMY  %>% dplyr::select(all_of(dtlist_AMY$GSE25219$GSM))

idx_del_AMY <- which(apply(GSE25219_AMY, 1, mad)==0)

if(length(idx_del_AMY)>0){
  GSE25219_AMY <- GSE25219_AMY[-idx_del_AMY, ] 
}

### 2.8.2 Form data ----
GSE25219_AMY <- as.data.frame(t(GSE25219_AMY))

### 2.8.3 Rudimentary data cleaning and outlier removal ----
# Check that all genes and samples have sufficiently low numbers of missing values.
gsg <- goodSamplesGenes(GSE25219_AMY,verbose = 3)
gsg$allOK

### 2.8.4 sample clustering ----
sampleTrees_AMY = hclust(dist(GSE25219_AMY), method = "average")

pdf(file = "2.output/9.dev_wgcna/figure/AMY/dtsplit/02_SampleClustering_AMY_cutHeights.pdf", width = 14, height = 7)
plot(sampleTrees_AMY, xlab="", sub="", cex = 0.7,
     main = paste("Sample clustering on all genes in AMY")) 
#abline(h=235, col = "red")
dev.off()

### 2.8.5. import clinical information -----
sample_AMY %<>% 
  dplyr::mutate(HW1 = ifelse(age_stage2=="HW1",1,0),
                HW2 = ifelse(age_stage2=="HW2",1,0),
                HW3 = ifelse(age_stage2=="HW3",1,0),
                HW4 = ifelse(age_stage2=="HW4",1,0),
                HW5 = ifelse(age_stage2=="HW5",1,0),
                HW6 = ifelse(age_stage2=="HW6",1,0)) %>% 
  dplyr::select(GSM,starts_with("HW")) %>% 
  unique() %>% 
  column_to_rownames("GSM")

setSamples = rownames(GSE25219_AMY)
traitRows = match(setSamples, rownames(sample_AMY))
Traits_AMY = sample_AMY[traitRows, ]

# Define data set dimensions
nGenes_AMY = ncol(GSE25219_AMY)
nSamples_AMY = nrow(GSE25219_AMY)

### 2.8.6. reclustering samples with trait ----
sampleTrees2_AMY = hclust(dist(GSE25219_AMY), method = "average")

traitColors = numbers2colors(Traits_AMY, signed = FALSE)

pdf(file = paste0("2.output/9.dev_wgcna/figure/AMY/dtsplit/03_SampleClustering_AMY_Trait.pdf"), width = 10, height = 5)
plotDendroAndColors(sampleTrees2_AMY, traitColors,
                    groupLabels = colnames(Traits_AMY), 
                    autoColorHeight = T,
                    colorHeight = 0.02,
                    colorHeightBase = 0.02,
                    main = "Sample dendrogram and trait heatmap-AMY")
dev.off()

### 2.8.7. export exp, trait, and other info ----
save(GSE25219_AMY, Traits_AMY, nGenes_AMY, nSamples_AMY, 
     file = "./2.output/9.dev_wgcna/data/02_NonConsensus_dataInput_AMY_GSE25219.RData")

######################################################
#                     3. Network construction and module detection
# Coexpression networks were constructed in each brain region
# select power -> network construction
######################################################

## 3.1 FCx ----
powers = c(seq(1,10,by=1), seq(12,30, by=2)) 
lnames=load("./2.output/9.dev_wgcna/data/02_Consensus_dataInput_FC.RData")

### 3.1.1 Choosing the soft-thresholding power： analysis of network topology ----
# delete asd casue of high connectivity (temp)
nSets_FC = checkSets(multiExpr_FC)$nSets

# Initialize a list to hold the results of scale-free analysis
powerTables_FC = vector(mode = "list", length = nSets_FC)
# Call the network topology analysis function for each set in turn
for (set in 1:nSets_FC)
{
  powerTables_FC[[set]] = list(data = pickSoftThreshold(multiExpr_FC[[set]]$data, 
                                                        powerVector=powers,
                                                        corFnc ="bicor",
                                                        networkType = "signed",
                                                        verbose = 2)[[2]])  
}

# Plot the results:
colors = RColorBrewer::brewer.pal(nSets_FC,"Set2")

# Will plot these columns of the returned scale free analysis tables
plotCols = c(2,5,6,7)

colNames = c("Scale Free Topology Model Fit", "Mean connectivity", "Median connectivity","Max connectivity")

# Get the minima and maxima of the plotted points
ylim = matrix(NA, nrow = nSets_FC, ncol = 4)
for (set in 1:nSets_FC)
{
  for (col in 1:length(plotCols))
  {
    ylim[set, col] = min(ylim[set, col], powerTables_FC[[set]]$data[, plotCols[col]], na.rm = TRUE)
  }
}

# Plot the quantities in the chosen columns vs. the soft thresholding power
sizeGrWindow(8, 6)

pdf("./2.output/9.dev_wgcna/figure/FCx/dtsplit/04_Summary network indices (y-axes) as functions of the soft thresholding power (x-axes).pdf")

layout(matrix(c(1,2,3,4,5,5),ncol = 2,byrow = T),heights = c(3,3,1))
cex1 = 1
for (col in 1:length(plotCols))
{
  # set is power, col is column of matrix:ylim
  for (set in 1:nSets_FC)
  {
    if (set==1)
    {
      if(col==1)
      {
        plot(powerTables_FC[[set]]$data[,1], -sign(powerTables_FC[[set]]$data[,3])*powerTables_FC[[set]]$data[,2],
             xlab="Soft Threshold (power)",ylab=colNames[col],type="n", ylim = c(0,1),
             main = colNames[col])
        abline(h=0.8, col = "red")        
      }else{
        plot(powerTables_FC[[set]]$data[,1], -sign(powerTables_FC[[set]]$data[,3])*powerTables_FC[[set]]$data[,2],
             xlab="Soft Threshold (power)",ylab=colNames[col],type="n", ylim = range(powerTables_FC[[set]]$data[,plotCols[col]]),
             main = colNames[col])
      }
      addGrid()  
    }
    
    if (col==1)
    {
      text(powerTables_FC[[set]]$data[,1], -sign(powerTables_FC[[set]]$data[,3])*powerTables_FC[[set]]$data[,2],
           labels=powers,cex=cex1,col=colors[set]);
    } else{
      text(powerTables_FC[[set]]$data[,1], powerTables_FC[[set]]$data[,plotCols[col]],
           labels=powers,cex=cex1,col=colors[set],ylim = range(ylim[,col]))
    }    
  } 
}
par(mar=c(0,0,0,0))
plot.new()
legend("center", legend = setLabels_FC, col = colors, pch = 20,cex = 1.5,ncol = 2)

dev.off()

### 3.1.2 Network construction and consensus module detection ----
softPower_FC = c(12,18,8,20)
net_FC = blockwiseConsensusModules(multiExpr_FC,
                                    power = softPower_FC,
                                    minModuleSize = 50,
                                    maxBlockSize =50000,
                                    corType = "bicor",
                                    networkType = "signed",
                                    TOMType = "signed",
                                    networkCalibration = "full quantile",
                                    calibrationQuantile = 0.95,
                                    #deepSplit = 2,
                                    mergeCutHeight = 0.25,
                                    minKMEtoStay = 0,
                                    saveTOMs = F,
                                    saveIndividualTOMs =F,
                                    saveConsensusTOMs =F,
                                    verbose = 5)
consMEs_FC = net_FC$multiMEs
moduleColors_FC = net_FC$colors
consTree_FC = net_FC$dendrograms[[1]]


pdf("./2.output/9.dev_wgcna/figure/FCx/dtsplit/05_Dendrogram of consensus module eigengenes obtained by clustering the eigengenes on their consensus correlation.pdf",width=7, height=4)
plotDendroAndColors(consTree_FC, moduleColors_FC,
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Consensus gene dendrogram and module colors")
dev.off()

save(consMEs_FC, moduleColors_FC, consTree_FC, file = "./2.output/9.dev_wgcna/data/03_Consensus-NetworkConstruction-FCx.RData")

## 3.2 TCx ----
powers = c(seq(1,10,by=1), seq(12,30, by=2)) 
lnames=load("./2.output/9.dev_wgcna/data/02_Consensus_dataInput_TC.RData")

### 3.2.1 Choosing the soft-thresholding power： analysis of network topology ----
nSets_TC = checkSets(multiExpr_TC)$nSets

# Choose a set of soft-thresholding powers

# Initialize a list to hold the results of scale-free analysis
powerTables_TC = vector(mode = "list", length = nSets_TC)
# Call the network topology analysis function for each set in turn
for (set in 1:nSets_TC)
{
  powerTables_TC[[set]] = list(data = pickSoftThreshold(multiExpr_TC[[set]]$data, 
                                                        powerVector=powers,
                                                        corFnc ="bicor",
                                                        networkType = "signed",
                                                        verbose = 2)[[2]])  
}

# Plot the results:
colors = RColorBrewer::brewer.pal(nSets_TC,"Set2")


# Will plot these columns of the returned scale free analysis tables
plotCols = c(2,5,6,7)

colNames = c("Scale Free Topology Model Fit", "Mean connectivity", "Median connectivity",
             "Max connectivity")

# Get the minima and maxima of the plotted points
ylim = matrix(NA, nrow = nSets_TC, ncol = 4)
for (set in 1:nSets_TC)
{
  for (col in 1:length(plotCols))
  {
    ylim[set, col] = min(ylim[set, col], powerTables_TC[[set]]$data[, plotCols[col]], na.rm = TRUE)
  }
}

# Plot the quantities in the chosen columns vs. the soft thresholding power
sizeGrWindow(8, 6)

pdf("./2.output/9.dev_wgcna/figure/TCx/dtsplit/04_Summary network indices (y-axes) as functions of the soft thresholding power (x-axes).pdf")
layout(matrix(c(1,2,3,4,5,5),ncol = 2,byrow = T),heights = c(3,3,1))
cex1 = 1

for (col in 1:length(plotCols))
{
  # set is power, col is column of matrix:ylim
  for (set in 1:nSets_TC)
  {
    if (set==1)
    {
      if(col==1)
      {
        plot(powerTables_TC[[set]]$data[,1], -sign(powerTables_TC[[set]]$data[,3])*powerTables_TC[[set]]$data[,2],
             xlab="Soft Threshold (power)",ylab=colNames[col],type="n", ylim = c(0,1),
             main = colNames[col])
        abline(h=0.8, col = "red")
      }else{
        plot(powerTables_TC[[set]]$data[,1], -sign(powerTables_TC[[set]]$data[,3])*powerTables_TC[[set]]$data[,2],
             xlab="Soft Threshold (power)",ylab=colNames[col],type="n", ylim = range(powerTables_TC[[set]]$data[,plotCols[col]]),
             main = colNames[col])
      }
      addGrid()  
    }
    if (col==1)
    {
      text(powerTables_TC[[set]]$data[,1], -sign(powerTables_TC[[set]]$data[,3])*powerTables_TC[[set]]$data[,2],
           labels=powers,cex=cex1,col=colors[set]);
    } else{
      text(powerTables_TC[[set]]$data[,1], powerTables_TC[[set]]$data[,plotCols[col]],
           labels=powers,cex=cex1,col=colors[set],ylim = range(ylim[,col]))
    }
  }
}
par(mar=c(0,0,0,0))
plot.new()
legend("center", legend = setLabels_TC, col = colors, pch = 20,cex = 1.5,ncol = 2)

dev.off()

### 3.2.2 Network construction and consensus module detection ----
setLabels_TC
softPower_TC = c(15,11,7,16)
net_TC = blockwiseConsensusModules(multiExpr_TC,
                                   power = softPower_TC,
                                   minModuleSize = 50,
                                   maxBlockSize =50000,
                                   corType = "bicor",
                                   networkType = "signed",
                                   TOMType = "signed",
                                   networkCalibration = "full quantile",
                                   calibrationQuantile = 0.95,
                                   #deepSplit = 2,
                                   mergeCutHeight = 0.25,
                                   minKMEtoStay = 0,
                                   saveTOMs = F,
                                   saveIndividualTOMs =F,
                                   saveConsensusTOMs =F,
                                   verbose = 5)
consMEs_TC = net_TC$multiMEs
moduleColors_TC = net_TC$colors
consTree_TC = net_TC$dendrograms[[1]]

pdf("./2.output/9.dev_wgcna/figure/TCx/dtsplit/05_Dendrogram of consensus module eigengenes obtained by clustering the eigengenes on their consensus correlation.pdf",width=7, height=4)
plotDendroAndColors(consTree_TC, moduleColors_TC,
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Consensus gene dendrogram and module colors")
dev.off()

save(consMEs_TC, moduleColors_TC, consTree_TC, file = "./2.output/9.dev_wgcna/data/03_Consensus-NetworkConstruction-TCx.RData")

## 3.3 OCx ----
powers = c(c(1:10), seq(from = 12, to=20, by=2))
lnames=load("./2.output/9.dev_wgcna/data/02_NonConsensus_dataInput_OC_GSE25219.RData")
### 3.3.1 Choosing the soft-thresholding power： analysis of network topology ----
sft = pickSoftThreshold(GSE25219_OC, powerVector = powers, verbose = 5)

# Plot the quantities in the chosen columns vs. the soft thresholding power
pdf("./2.output/9.dev_wgcna/figure/OCx/dtsplit/04_Summary network indices (y-axes) as functions of the soft thresholding power (x-axes).pdf",width = 6.5,height = 5.5)
par(mfrow = c(1,2))
cex1 = 0.9
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"))
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red")
# this line corresponds to using an R^2 cut-off of h
abline(h=0.80,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")

dev.off()

### 3.3.2 Network construction and consensus module detection ----
net_OC = blockwiseModules(GSE25219_OC, 
                          power = 9,
                          TOMType = "signed", 
                          corType = "bicor",
                          networkType = "signed",
                          minModuleSize = 50,
                          maxBlockSize = 50000,
                          mergeCutHeight = 0.25,
                          numericLabels = FALSE, 
                          pamRespectsDendro = FALSE,
                          saveTOMs = FALSE,
                          verbose = 5)

MEs_OC = net_OC$MEs
moduleColors_OC = net_OC$colors
Tree_OC = net_OC$dendrograms[[1]]

pdf("./2.output/9.dev_wgcna/figure/OCx/dtsplit/05_Dendrogram of genes with dissimilarity based on topological overlap together with assigned module colors.pdf",width=7, height=4)
plotDendroAndColors(Tree_OC, moduleColors_OC,
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Consensus gene dendrogram and module colors")
dev.off()
save(MEs_OC, moduleColors_OC, Tree_OC, file = "./2.output/9.dev_wgcna/data/03_NonConsensus-NetworkConstruction-OCx.RData")

## 3.4 PCx ----
powers = c(c(1:10), seq(from = 12, to=20, by=2))
lnames=load("./2.output/9.dev_wgcna/data/02_NonConsensus_dataInput_PC_GSE25219.RData")

### 3.4.1 Choosing the soft-thresholding power： analysis of network topology ----
sft = pickSoftThreshold(GSE25219_PC, powerVector = powers, verbose = 5)

# Plot the quantities in the chosen columns vs. the soft thresholding power
pdf("./2.output/9.dev_wgcna/figure/PCx/dtsplit/04_Summary network indices (y-axes) as functions of the soft thresholding power (x-axes).pdf",width = 6.5,height = 5.5)
par(mfrow = c(1,2))
cex1 = 0.9
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"))
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red")
# this line corresponds to using an R^2 cut-off of h
abline(h=0.80,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")

dev.off()

### 3.4.2 Network construction and consensus module detection ----
net_PC = blockwiseModules(GSE25219_PC, 
                          power = 6,
                          TOMType = "signed", 
                          corType = "bicor",
                          networkType = "signed",
                          minModuleSize = 50,
                          maxBlockSize = 50000,
                          mergeCutHeight = 0.25,
                          numericLabels = FALSE, 
                          pamRespectsDendro = FALSE,
                          saveTOMs = FALSE,
                          verbose = 5)

MEs_PC = net_PC$MEs
moduleColors_PC = net_PC$colors
Tree_PC = net_PC$dendrograms[[1]]

pdf("./2.output/9.dev_wgcna/figure/PCx/dtsplit/05_Dendrogram of genes with dissimilarity based on topological overlap together with assigned module colors.pdf",width=7, height=4)
plotDendroAndColors(Tree_PC, moduleColors_PC,
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Consensus gene dendrogram and module colors")
dev.off()
save(MEs_PC, moduleColors_PC, Tree_PC, file = "./2.output/9.dev_wgcna/data/03_NonConsensus-NetworkConstruction-PCx.RData")

## 3.5 CB ----
powers = c(c(1:10), seq(from = 12, to=20, by=2))
lnames=load("./2.output/9.dev_wgcna/data/02_NonConsensus_dataInput_CB_GSE25219.RData")

### 3.5.1 Choosing the soft-thresholding power： analysis of network topology ----
sft = pickSoftThreshold(GSE25219_CB, powerVector = powers, verbose = 5)

# Plot the quantities in the chosen columns vs. the soft thresholding power
pdf("./2.output/9.dev_wgcna/figure/CB/dtsplit/04_Summary network indices (y-axes) as functions of the soft thresholding power (x-axes).pdf",width = 6.5,height = 5.5)
par(mfrow = c(1,2))
cex1 = 0.9
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"))
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red")
# this line corresponds to using an R^2 cut-off of h
abline(h=0.80,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")

dev.off()

### 3.5.2 Network construction and consensus module detection ----
net_CB = blockwiseModules(GSE25219_CB, 
                          power = 12,
                          TOMType = "signed", 
                          corType = "bicor",
                          networkType = "signed",
                          minModuleSize = 50,
                          maxBlockSize = 50000,
                          mergeCutHeight = 0.25,
                          numericLabels = FALSE, 
                          pamRespectsDendro = FALSE,
                          saveTOMs = FALSE,
                          verbose = 5)

MEs_CB = net_CB$MEs
moduleColors_CB = net_CB$colors
Tree_CB = net_CB$dendrograms[[1]]

pdf("./2.output/9.dev_wgcna/figure/CB/dtsplit/05_Dendrogram of genes with dissimilarity based on topological overlap together with assigned module colors.pdf",width=7, height=4)
plotDendroAndColors(Tree_CB, moduleColors_CB,
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Consensus gene dendrogram and module colors")
dev.off()
save(MEs_CB, moduleColors_CB, Tree_CB, file = "./2.output/9.dev_wgcna/data/03_NonConsensus-NetworkConstruction-CB.RData")

## 3.6 HIP ----
powers = c(c(1:10), seq(from = 12, to=20, by=2))
lnames=load("./2.output/9.dev_wgcna/data/02_NonConsensus_dataInput_HIP_GSE25219.RData")
### 3.6.1 Choosing the soft-thresholding power： analysis of network topology ----
sft = pickSoftThreshold(GSE25219_HIP, powerVector = powers, verbose = 5)

# Plot the quantities in the chosen columns vs. the soft thresholding power
pdf("./2.output/9.dev_wgcna/figure/HIP/dtsplit/04_Summary network indices (y-axes) as functions of the soft thresholding power (x-axes).pdf",width = 6.5,height = 5.5)
par(mfrow = c(1,2))
cex1 = 0.9
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"))
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red")
# this line corresponds to using an R^2 cut-off of h
abline(h=0.80,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")

dev.off()

### 3.6.2 Network construction and consensus module detection ----
net_HIP = blockwiseModules(GSE25219_HIP, 
                          power = 7,
                          TOMType = "signed", 
                          corType = "bicor",
                          networkType = "signed",
                          minModuleSize = 50,
                          maxBlockSize = 50000,
                          mergeCutHeight = 0.25,
                          numericLabels = FALSE, 
                          pamRespectsDendro = FALSE,
                          saveTOMs = FALSE,
                          verbose = 5)

MEs_HIP = net_HIP$MEs
moduleColors_HIP = net_HIP$colors
Tree_HIP = net_HIP$dendrograms[[1]]

pdf("./2.output/9.dev_wgcna/figure/HIP/dtsplit/05_Dendrogram of genes with dissimilarity based on topological overlap together with assigned module colors.pdf",width=7, height=4)
plotDendroAndColors(Tree_HIP, moduleColors_HIP,
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Consensus gene dendrogram and module colors")
dev.off()
save(MEs_HIP, moduleColors_HIP, Tree_HIP, file = "./2.output/9.dev_wgcna/data/03_NonConsensus-NetworkConstruction-HIP.RData")

## 3.7 STR ----
powers = c(c(1:10), seq(from = 12, to=20, by=2))
lnames=load("./2.output/9.dev_wgcna/data/02_NonConsensus_dataInput_STR_GSE25219.RData")
### 3.7.1 Choosing the soft-thresholding power： analysis of network topology ----
sft = pickSoftThreshold(GSE25219_STR, powerVector = powers, verbose = 5)

# Plot the quantities in the chosen columns vs. the soft thresholding power
pdf("./2.output/9.dev_wgcna/figure/STR/dtsplit/04_Summary network indices (y-axes) as functions of the soft thresholding power (x-axes).pdf",width = 6.5,height = 5.5)
par(mfrow = c(1,2))
cex1 = 0.9
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"))
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red")
# this line corresponds to using an R^2 cut-off of h
abline(h=0.80,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")

dev.off()

### 3.7.2 Network construction and consensus module detection ----
net_STR = blockwiseModules(GSE25219_STR, 
                           power = 7,
                           TOMType = "signed", 
                           corType = "bicor",
                           networkType = "signed",
                           minModuleSize = 50,
                           maxBlockSize = 50000,
                           mergeCutHeight = 0.25,
                           numericLabels = FALSE, 
                           pamRespectsDendro = FALSE,
                           saveTOMs = FALSE,
                           verbose = 5)

MEs_STR = net_STR$MEs
moduleColors_STR = net_STR$colors
Tree_STR = net_STR$dendrograms[[1]]

pdf("./2.output/9.dev_wgcna/figure/STR/dtsplit/05_Dendrogram of genes with dissimilarity based on topological overlap together with assigned module colors.pdf",width=7, height=4)
plotDendroAndColors(Tree_STR, moduleColors_STR,
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Consensus gene dendrogram and module colors")
dev.off()
save(MEs_STR, moduleColors_STR, Tree_STR, file = "./2.output/9.dev_wgcna/data/03_NonConsensus-NetworkConstruction-STR.RData")

## 3.8 AMY ----
powers = c(c(1:10), seq(from = 12, to=20, by=2))
lnames=load("./2.output/9.dev_wgcna/data/02_NonConsensus_dataInput_AMY_GSE25219.RData")

### 3.8.1 Choosing the soft-thresholding power： analysis of network topology ----
sft = pickSoftThreshold(GSE25219_AMY, powerVector = powers, verbose = 5)

# Plot the quantities in the chosen columns vs. the soft thresholding power
pdf("./2.output/9.dev_wgcna/figure/AMY/dtsplit/04_Summary network indices (y-axes) as functions of the soft thresholding power (x-axes).pdf",width = 6.5,height = 5.5)
par(mfrow = c(1,2))
cex1 = 0.9
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"))
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red")
# this line corresponds to using an R^2 cut-off of h
abline(h=0.80,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")

dev.off()

### 3.8.2 Network construction and consensus module detection ----
net_AMY = blockwiseModules(GSE25219_AMY, 
                           power = 7,
                           TOMType = "signed", 
                           corType = "bicor",
                           networkType = "signed",
                           minModuleSize = 50,
                           maxBlockSize = 50000,
                           mergeCutHeight = 0.25,
                           numericLabels = FALSE, 
                           pamRespectsDendro = FALSE,
                           saveTOMs = FALSE,
                           verbose = 5)

MEs_AMY = net_AMY$MEs
moduleColors_AMY = net_AMY$colors
Tree_AMY = net_AMY$dendrograms[[1]]

pdf("./2.output/9.dev_wgcna/figure/AMY/dtsplit/05_Dendrogram of genes with dissimilarity based on topological overlap together with assigned module colors.pdf",width=7, height=4)
plotDendroAndColors(Tree_AMY, moduleColors_AMY,
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Consensus gene dendrogram and module colors")
dev.off()
save(MEs_AMY, moduleColors_AMY, Tree_AMY, file = "./2.output/9.dev_wgcna/data/03_NonConsensus-NetworkConstruction-AMY.RData")

######################################################
#                     4. Modules Trait Correlation
######################################################
get_consensusCorAndPvalue <- function(moduleTraitCor,moduleTraitPvalue){
  
  TraitNames = colnames(moduleTraitCor[[1]])
  ModNames = rownames(moduleTraitCor[[1]])
  # form data
  consensusCor = matrix(NA, nrow(moduleTraitCor[[1]]), ncol(moduleTraitCor[[1]]))
  consensusPvalue = matrix(NA, nrow(moduleTraitCor[[1]]), ncol(moduleTraitCor[[1]]))
  
  # Find consensus negative correlations
  nSets = length(moduleTraitCor)
  
  get_idx_neg <- function(x,y){
    # numeric comparison
    if(is.numeric(x[1,1])){
      x < 0 & y < 0
    }else{
      # change logical to numereric 
      x = apply(x,2,as.numeric)
      # numeric comparison
      x == 1 & y < 0
    }
  }
  
  get_idx_pos <- function(x,y){
    # numeric comparison
    if(is.numeric(x[1,1])){
      x > 0 & y > 0
    }else{
      # change logical to numereric 
      x = apply(x,2,as.numeric)
      # numeric comparison
      x == 1 & y > 0
    }
  }
  
  negative <-  Reduce(get_idx_neg,moduleTraitCor)
  
  consensusCor[negative] = Reduce(pmax,lapply(moduleTraitCor,function(x){x[negative]}))
  consensusPvalue[negative] = Reduce(pmax,lapply(moduleTraitPvalue,function(x){x[negative]}))
  
  positive <-  Reduce(get_idx_pos,moduleTraitCor)
  
  consensusCor[positive] = Reduce(pmin,lapply(moduleTraitCor,function(x){x[positive]}))
  consensusPvalue[positive] = Reduce(pmax,lapply(moduleTraitPvalue,function(x){x[positive]}))
  
  colnames(consensusPvalue) <- colnames(consensusCor) <- TraitNames
  rownames(consensusPvalue) <- rownames(consensusCor) <- ModNames
  
  list(cor= consensusCor, p = consensusPvalue)
}

## 4.1 FCx ----
load(file = "2.output/9.dev_wgcna/data/02_Consensus_dataInput_FC.RData")
load(file = "2.output/9.dev_wgcna/data/03_Consensus-NetworkConstruction-FCx.RData")

exprSize_FC = checkSets(multiExpr_FC)
nSets_FC = exprSize_FC$nSets
nTrait = 6

### 4.1.1 Calculate the correlations -----
# Set up variables to contain the module-trait correlations
moduleTraitCor_FC = list()
moduleTraitPvalue_FC = list()

# Calculate the correlations
for (set in 1:nSets_FC)
{
    moduleTraitCor_FC[[set]] = bicorAndPvalue(consMEs_FC[[set]]$data, Traits_FC[[set]]$data, use = "p")$bicor
    moduleTraitPvalue_FC[[set]] = bicorAndPvalue(consMEs_FC[[set]]$data, Traits_FC[[set]]$data, use = "p")$p
}

# Convert numerical lables to colors for labeling of modules in the plot
MEColorNames_FC = names(consMEs_FC[[1]]$data)

### 4.1.2 calculate the consensus correlation ----
# Initialize matrices to hold the consensus correlation and p-value
consensus_FC = get_consensusCorAndPvalue(moduleTraitCor = moduleTraitCor_FC,moduleTraitPvalue = moduleTraitPvalue_FC)

### 4.1.3 Export result of module trait correlation ----
for (set in 1:nSets_FC) {
  
  moduleTraitCorRes = cbind(cor = moduleTraitCor_FC[[set]], 
                            p = moduleTraitPvalue_FC[[set]]) %>% as.data.frame()
  
  colnames(moduleTraitCorRes) <- paste0(setLabels_FC[[set]],".",
                                        colnames(moduleTraitCorRes),
                                        rep(c(".cor",".p"),each=6)
  )
  
  fwrite(moduleTraitCorRes %>% rownames_to_column("module"),
         file = paste0("2.output/9.dev_wgcna/data/04_ModuleTraitCorRes/FCx/04_moduleTraitCorRes-FCx-",setLabels_FC[set],".txt"),
         sep = "\t")  
}

moduleTraitCorRes = cbind(consensus_FC$cor, consensus_FC$p) %>% as.data.frame()
colnames(moduleTraitCorRes) <- paste0("consensus.",colnames(moduleTraitCorRes),".",rep(c("cor","p"),each=6))
fwrite(moduleTraitCorRes %>% rownames_to_column("module"),
       file = paste0("2.output/9.dev_wgcna/data/04_ModuleTraitCorRes/FCx/04_moduleTraitCorRes-FCx-consensus.txt"),
       sep = "\t")  
  
## 4.2 TCx ----
load(file = "2.output/9.dev_wgcna/data/02_Consensus_dataInput_TC.RData")
load(file = "2.output/9.dev_wgcna/data/03_Consensus-NetworkConstruction-TCx.RData")

exprSize_TC = checkSets(multiExpr_TC)
nSets_TC = exprSize_TC$nSets
nTrait = 6

### 4.2.1 Calculate the correlations -----
# Set up variables to contain the module-trait correlations
moduleTraitCor_TC = list()
moduleTraitPvalue_TC = list()

# Calculate the correlations
for (set in 1:nSets_TC)
{
  moduleTraitCor_TC[[set]] = bicorAndPvalue(consMEs_TC[[set]]$data, Traits_TC[[set]]$data, use = "p")$bicor
  moduleTraitPvalue_TC[[set]] = bicorAndPvalue(consMEs_TC[[set]]$data, Traits_TC[[set]]$data, use = "p")$p
}

# Convert numerical lables to colors for labeling of modules in the plot
MEColorNames_TC = names(consMEs_TC[[1]]$data)

### 4.2.2 calculate the consensus correlation ----
# Initialize matrices to hold the consensus correlation and p-value
consensus_TC = get_consensusCorAndPvalue(moduleTraitCor = moduleTraitCor_TC,moduleTraitPvalue = moduleTraitPvalue_TC)

### 4.2.3 Export result of module trait correlation ----
for (set in 1:nSets_TC) {
  
  moduleTraitCorRes = cbind(cor = moduleTraitCor_TC[[set]], 
                            p = moduleTraitPvalue_TC[[set]]) %>% as.data.frame()
  
  colnames(moduleTraitCorRes) <- paste0(setLabels_TC[[set]],".",
                                        colnames(moduleTraitCorRes),
                                        rep(c(".cor",".p"),each=6)
  )
  
  fwrite(moduleTraitCorRes %>% rownames_to_column("module"),
         file = paste0("2.output/9.dev_wgcna/data/04_ModuleTraitCorRes/TCx/04_moduleTraitCorRes-TCx-",setLabels_TC[set],".txt"),
         sep = "\t")  
}

moduleTraitCorRes = cbind(consensus_TC$cor, consensus_TC$p) %>% as.data.frame()
colnames(moduleTraitCorRes) <- paste0("consensus.",colnames(moduleTraitCorRes),".",rep(c("cor","p"),each=6))
fwrite(moduleTraitCorRes %>% rownames_to_column("module"),
       file = paste0("2.output/9.dev_wgcna/data/04_ModuleTraitCorRes/TCx/04_moduleTraitCorRes-TCx-consensus.txt"),
       sep = "\t")  

## 4.3 OCx ----
rm(list=ls());gc()
load(file = "2.output/9.dev_wgcna/data/02_NonConsensus_dataInput_OC_GSE25219.RData")
load(file = "2.output/9.dev_wgcna/data/03_NonConsensus-NetworkConstruction-OCx.RData")

### 4.3.1 Calculate the correlations -----
# Set up variables to contain the module-trait correlations
moduleTraitCor_OC = bicorAndPvalue(MEs_OC, Traits_OC, use = "p")$bicor
moduleTraitPvalue_OC = bicorAndPvalue(MEs_OC, Traits_OC, use = "p")$p

# Convert numerical lables to colors for labeling of modules in the plot
MEColorNames_OC = names(MEs_OC)

### 4.3.2 Export result of module trait correlation ----
moduleTraitCorRes = cbind(cor = moduleTraitCor_OC, 
                          p = moduleTraitPvalue_OC) %>% as.data.frame()
colnames(moduleTraitCorRes) <- paste0("GSE25219",".",
                                      colnames(moduleTraitCorRes),".",
                                      rep(c("cor","p"),each=6)
                                      )
fwrite(moduleTraitCorRes %>% rownames_to_column("module"),
       file = paste0("2.output/9.dev_wgcna/data/04_ModuleTraitCorRes/OCx/04_ModuleTraitCorRes-GSE25219.txt"),
       sep = "\t")  

## 4.4 PCx ----
rm(list=ls());gc()
load(file = "2.output/9.dev_wgcna/data/02_NonConsensus_dataInput_PC_GSE25219.RData")
load(file = "2.output/9.dev_wgcna/data/03_NonConsensus-NetworkConstruction-PCx.RData")

### 4.4.1 Calculate the correlations -----
# Set up variables to contain the module-trait correlations
moduleTraitCor_PC = bicorAndPvalue(MEs_PC, Traits_PC, use = "p")$bicor
moduleTraitPvalue_PC = bicorAndPvalue(MEs_PC, Traits_PC, use = "p")$p

# Convert numerical lables to colors for labeling of modules in the plot
MEColorNames_PC = names(MEs_PC)

### 4.4.3 Export result of module trait correlation ----
moduleTraitCorRes = cbind(cor = moduleTraitCor_PC, 
                          p = moduleTraitPvalue_PC) %>% as.data.frame()
colnames(moduleTraitCorRes) <- paste0("GSE25219",".",
                                      colnames(moduleTraitCorRes),".",
                                      rep(c("cor","p"),each=6)
)
fwrite(moduleTraitCorRes %>% rownames_to_column("module"),
       file = paste0("2.output/9.dev_wgcna/data/04_ModuleTraitCorRes/PCx/04_ModuleTraitCorRes-GSE25219.txt"),
       sep = "\t")  

## 4.5 CB ----
load(file = "2.output/9.dev_wgcna/data/02_NonConsensus_dataInput_CB_GSE25219.RData")
load(file = "2.output/9.dev_wgcna/data/03_NonConsensus-NetworkConstruction-CB.RData")

### 4.5.1 Calculate the correlations -----
# Set up variables to contain the module-trait correlations
moduleTraitCor_CB = bicorAndPvalue(MEs_CB, Traits_CB, use = "p")$bicor
moduleTraitPvalue_CB = bicorAndPvalue(MEs_CB, Traits_CB, use = "p")$p

# Convert numerical lables to colors for labeling of modules in the plot
MEColorNames_CB = names(MEs_CB)

### 4.5.2 Export result of module trait correlation ----
moduleTraitCorRes = cbind(cor = moduleTraitCor_CB, 
                          p = moduleTraitPvalue_CB) %>% as.data.frame()
colnames(moduleTraitCorRes) <- paste0("GSE25219",".",
                                      colnames(moduleTraitCorRes),".",
                                      rep(c("cor","p"),each=6)
)
fwrite(moduleTraitCorRes %>% rownames_to_column("module"),
       file = paste0("2.output/9.dev_wgcna/data/04_ModuleTraitCorRes/CB/04_ModuleTraitCorRes-GSE25219.txt"),
       sep = "\t")  

## 4.6 HIP ----
load(file = "2.output/9.dev_wgcna/data/02_NonConsensus_dataInput_HIP_GSE25219.RData")
load(file = "2.output/9.dev_wgcna/data/03_NonConsensus-NetworkConstruction-HIP.RData")

### 4.6.1 Calculate the correlations -----
# Set up variables to contain the module-trait correlations
moduleTraitCor_HIP = bicorAndPvalue(MEs_HIP, Traits_HIP, use = "p")$bicor
moduleTraitPvalue_HIP = bicorAndPvalue(MEs_HIP, Traits_HIP, use = "p")$p

# Convert numerical lables to colors for labeling of modules in the plot
MEColorNames_HIP = names(MEs_HIP)

### 4.6.2 Export result of module trait correlation ----
moduleTraitCorRes = cbind(cor = moduleTraitCor_HIP, 
                          p = moduleTraitPvalue_HIP) %>% as.data.frame()
colnames(moduleTraitCorRes) <- paste0("GSE25219",".",
                                      colnames(moduleTraitCorRes),".",
                                      rep(c("cor","p"),each=6)
)
fwrite(moduleTraitCorRes %>% rownames_to_column("module"),
       file = paste0("2.output/9.dev_wgcna/data/04_ModuleTraitCorRes/HIP/04_ModuleTraitCorRes-GSE25219.txt"),
       sep = "\t")  

## 4.7 STR ----
load(file = "2.output/9.dev_wgcna/data/02_NonConsensus_dataInput_STR_GSE25219.RData")
load(file = "2.output/9.dev_wgcna/data/03_NonConsensus-NetworkConstruction-STR.RData")

### 4.7.1 Calculate the correlations -----
# Set up variables to contain the module-trait correlations
moduleTraitCor_STR = bicorAndPvalue(MEs_STR, Traits_STR, use = "p")$bicor
moduleTraitPvalue_STR = bicorAndPvalue(MEs_STR, Traits_STR, use = "p")$p

# Convert numerical lables to colors for labeling of modules in the plot
MEColorNames_STR = names(MEs_STR)

### 4.7.2 Export result of module trait correlation ----
moduleTraitCorRes = cbind(cor = moduleTraitCor_STR, 
                          p = moduleTraitPvalue_STR) %>% as.data.frame()
colnames(moduleTraitCorRes) <- paste0("GSE25219",".",
                                      colnames(moduleTraitCorRes),".",
                                      rep(c("cor","p"),each=6)
)
fwrite(moduleTraitCorRes %>% rownames_to_column("module"),
       file = paste0("2.output/9.dev_wgcna/data/04_ModuleTraitCorRes/STR/04_ModuleTraitCorRes-GSE25219.txt"),
       sep = "\t")  

## 4.8 AMY ----
load(file = "2.output/9.dev_wgcna/data/02_NonConsensus_dataInput_AMY_GSE25219.RData")
load(file = "2.output/9.dev_wgcna/data/03_NonConsensus-NetworkConstruction-AMY.RData")

### 4.8.1 Calculate the correlations -----
# Set up variables to contain the module-trait correlations
moduleTraitCor_AMY = bicorAndPvalue(MEs_AMY, Traits_AMY, use = "p")$bicor
moduleTraitPvalue_AMY = bicorAndPvalue(MEs_AMY, Traits_AMY, use = "p")$p

# Convert numerical lables to colors for labeling of modules in the plot
MEColorNames_AMY = names(MEs_AMY)

### 4.8.2 Export result of module trait correlation ----
moduleTraitCorRes = cbind(cor = moduleTraitCor_AMY, 
                          p = moduleTraitPvalue_AMY) %>% as.data.frame()
colnames(moduleTraitCorRes) <- paste0("GSE25219",".",
                                      colnames(moduleTraitCorRes),".",
                                      rep(c("cor","p"),each=6)
)
fwrite(moduleTraitCorRes %>% rownames_to_column("module"),
       file = paste0("2.output/9.dev_wgcna/data/04_ModuleTraitCorRes/AMY/04_ModuleTraitCorRes-GSE25219.txt"),
       sep = "\t")  

######################################################
#                     5.Exporting results of the network analysis
######################################################
# 5.1 Module Trait Correlation Results collection -----
tidy_cor_res <- function(data){
  data %>% 
      pivot_longer(cols = contains("cor"),names_to = "CorType",values_to = "Cor") %>% 
      pivot_longer(cols = ends_with("p"),names_to = "PType",values_to = "p") %>% 
      dplyr::filter(str_extract(CorType,"HW\\d")==str_extract(PType,"HW\\d")) %>% 
      dplyr::mutate(Stage = str_extract(CorType,"HW\\d")) %>% 
      dplyr::select(module,Stage,Cor,p) %>% 
}
## 5.1.1 fc ----
cor_fc_brainspanArray <- fread("2.output/9.dev_wgcna/data/04_ModuleTraitCorRes/FCx/04_moduleTraitCorRes-FCx-BrainspanArray.txt")
cor_fc_brainspanRnaseq <- fread("2.output/9.dev_wgcna/data/04_ModuleTraitCorRes/FCx/04_moduleTraitCorRes-FCx-BrainspanRNAseq.txt")
cor_fc_gse25219 <- fread("2.output/9.dev_wgcna/data/04_ModuleTraitCorRes/FCx/04_moduleTraitCorRes-FCx-GSE25219.txt")
cor_fc_psychencode <- fread("2.output/9.dev_wgcna/data/04_ModuleTraitCorRes/FCx/04_moduleTraitCorRes-FCx-PsychENCODE.txt")

cor_fc <- rbind(
  cor_fc_brainspanArray %>% 
    tidy_cor_res() %>% 
    dplyr::mutate(datasets = "BrainspanArray"),
  cor_fc_brainspanRnaseq %>% 
    tidy_cor_res() %>% 
    dplyr::mutate(datasets = "BrainspanRnaseq"),
  cor_fc_gse25219 %>% 
    tidy_cor_res() %>% 
    dplyr::mutate(datasets = "GSE25219"),
  cor_fc_psychencode %>% 
    tidy_cor_res() %>% 
    dplyr::mutate(datasets = "PsychENCODE")
) 

cor_fc_res <- rbind(
  cor_fc %>% 
    group_by(Stage,datasets) %>% 
    dplyr::mutate(fdr = p.adjust(p,method = "BH")) %>% 
    dplyr::select(module,Stage,Cor,p,fdr,datasets),
  
  cor_fc %>% 
    group_by(Stage,datasets) %>% 
    dplyr::mutate(fdr = p.adjust(p,method = "BH")) %>% 
    ungroup() %>% 
    group_by(module,Stage) %>% 
    dplyr::mutate(Cor = case_when(max(Cor)<=0 ~ max(Cor),
                                  min(Cor)>=0 ~ min(Cor),
                                  TRUE ~ NA)) %>% 
    dplyr::mutate(p = ifelse(is.na(Cor),NA,max(p)),
                  fdr = ifelse(is.na(Cor),NA,max(fdr))) %>% 
    dplyr::select(module,Stage,Cor,p,fdr) %>% 
    dplyr::mutate(datasets="consensus") %>% unique()
)

## 5.1.2 tc ----
cor_tc_brainspanArray <- fread("2.output/9.dev_wgcna/data/04_ModuleTraitCorRes/TCx/04_moduleTraitCorRes-TCx-BrainspanArray.txt")
cor_tc_brainspanRnaseq <- fread("2.output/9.dev_wgcna/data/04_ModuleTraitCorRes/TCx/04_moduleTraitCorRes-TCx-BrainspanRNAseq.txt")
cor_tc_gse25219 <- fread("2.output/9.dev_wgcna/data/04_ModuleTraitCorRes/TCx/04_moduleTraitCorRes-TCx-GSE25219.txt")
cor_tc_psychencode <- fread("2.output/9.dev_wgcna/data/04_ModuleTraitCorRes/TCx/04_moduleTraitCorRes-TCx-PsychENCODE.txt")

cor_tc_asd <- fread("2.output/9.dev_wgcna/data/04_ModuleTraitCorRes/TCx/04_ModuleTraitCorRes-ASD.txt")
cor_tc_dem <- fread("2.output/9.dev_wgcna/data/04_ModuleTraitCorRes/TCx/04_ModuleTraitCorRes-DEM.txt")
cor_tc_scz <- fread("2.output/9.dev_wgcna/data/04_ModuleTraitCorRes/TCx/04_ModuleTraitCorRes-SCZ.txt")

cor_tc <- rbind(
  cor_tc_brainspanArray %>% 
    tidy_cor_res() %>% 
    dplyr::mutate(datasets = "BrainspanArray"),
  cor_tc_brainspanRnaseq %>% 
    tidy_cor_res() %>% 
    dplyr::mutate(datasets = "BrainspanRnaseq"),
  cor_tc_gse25219 %>% 
    tidy_cor_res() %>% 
    dplyr::mutate(datasets = "GSE25219"),
  cor_tc_psychencode %>% 
    tidy_cor_res() %>% 
    dplyr::mutate(datasets = "PsychENCODE")
) 

cor_tc_res <- rbind(
  cor_tc %>% 
    group_by(Stage,datasets) %>% 
    dplyr::mutate(fdr = p.adjust(p,method = "BH")) %>% 
    dplyr::select(module,Stage,Cor,p,fdr,datasets),
  
  cor_tc %>% 
    group_by(Stage,datasets) %>% 
    dplyr::mutate(fdr = p.adjust(p,method = "BH")) %>% 
    ungroup() %>% 
    group_by(module,Stage) %>% 
    dplyr::mutate(Cor = case_when(max(Cor)<=0 ~ max(Cor),
                                  min(Cor)>=0 ~ min(Cor),
                                  TRUE ~ NA)) %>% 
    dplyr::mutate(p = ifelse(is.na(Cor),NA,max(p)),
                  fdr = ifelse(is.na(Cor),NA,max(fdr))) %>% 
    dplyr::select(module,Stage,Cor,p,fdr) %>% 
    dplyr::mutate(datasets="consensus") %>% unique()
)

## 5.1.3 oc ----
cor_oc_gse25219 <- fread("2.output/9.dev_wgcna/data/04_ModuleTraitCorRes/OCx/04_ModuleTraitCorRes-GSE25219.txt")

cor_oc_res <- cor_oc_gse25219 %>% 
  tidy_cor_res() %>% 
  group_by(Stage) %>% 
  dplyr::mutate(fdr = p.adjust(p,method = "BH")) %>% 
  dplyr::select(module,Stage,Cor,p,fdr) %>% 
  dplyr::mutate(datasets = "GSE25219")

## 5.1.4 pc ----
cor_pc_gse25219 <- fread("2.output/9.dev_wgcna/data/04_ModuleTraitCorRes/PCx/04_ModuleTraitCorRes-GSE25219.txt")

cor_pc_res <- cor_pc_gse25219 %>% 
  tidy_cor_res() %>% 
  group_by(Stage) %>% 
  dplyr::mutate(fdr = p.adjust(p,method = "BH")) %>% 
  dplyr::select(module,Stage,Cor,p,fdr) %>% 
  dplyr::mutate(datasets = "GSE25219")


## 5.1.5 cb ----
cor_cb_gse25219 <- fread("2.output/9.dev_wgcna/data/04_ModuleTraitCorRes/CB/04_ModuleTraitCorRes-GSE25219.txt")

cor_cb_res <- cor_cb_gse25219 %>% 
  tidy_cor_res() %>% 
  group_by(Stage) %>% 
  dplyr::mutate(fdr = p.adjust(p,method = "BH")) %>% 
  dplyr::select(module,Stage,Cor,p,fdr) %>% 
  dplyr::mutate(datasets = "GSE25219")

## 5.1.6 hip ----
cor_hip_gse25219 <- fread("2.output/9.dev_wgcna/data/04_ModuleTraitCorRes/HIP/04_ModuleTraitCorRes-GSE25219.txt")

cor_hip_res <- cor_hip_gse25219 %>% 
  tidy_cor_res() %>% 
  group_by(Stage) %>% 
  dplyr::mutate(fdr = p.adjust(p,method = "BH")) %>% 
  dplyr::select(module,Stage,Cor,p,fdr) %>% 
  dplyr::mutate(datasets = "GSE25219")

## 5.1.7 str ----
cor_str_gse25219 <- fread("2.output/9.dev_wgcna/data/04_ModuleTraitCorRes/STR/04_ModuleTraitCorRes-GSE25219.txt")

cor_str_res <- cor_str_gse25219 %>% 
  tidy_cor_res() %>% 
  group_by(Stage) %>% 
  dplyr::mutate(fdr = p.adjust(p,method = "BH")) %>% 
  dplyr::select(module,Stage,Cor,p,fdr) %>% 
  dplyr::mutate(datasets = "GSE25219")

## 5.1.8 amy ----
cor_amy_gse25219 <- fread("2.output/9.dev_wgcna/data/04_ModuleTraitCorRes/AMY/04_ModuleTraitCorRes-GSE25219.txt")

cor_amy_res <- cor_amy_gse25219 %>% 
  tidy_cor_res() %>% 
  group_by(Stage) %>% 
  dplyr::mutate(fdr = p.adjust(p,method = "BH")) %>% 
  dplyr::select(module,Stage,Cor,p,fdr) %>% 
  dplyr::mutate(datasets = "GSE25219")

## 5.1.9 save res ----
cor_dev <- mget(ls(pattern = "cor_.*_res"))
names(cor_dev) <- toupper(gsub("cor_|_res","",names(cor_dev)))

save(cor_dev,file = "./2.output/9.dev_wgcna/data/04_ModuleTraitCorRes/ModuleTraitbiCorRes_DEV.Rdata")

# 5.2 Exporting results of the network analysis ----
get_network_GS_kME = function(multiExpr,moduleColors,Traits,setLabels,nSets,outdir=NULL,region=NULL){
  # get the network GS and kME for further analysis the important genes in module
  # return:
  #   the resuts include genes, module lables, GS of genes and traits, meta Result, kME of genes and  
  
  consMEs.unord = multiSetMEs(multiExpr, universalColors = moduleColors, excludeGrey = TRUE)
  GS = list()
  kME = list()
  for (set in 1:nSets)
  {
    GS[[set]] = bicorAndPvalue(multiExpr[[set]]$data, Traits[[set]]$data)
    kME[[set]] = bicorAndPvalue(multiExpr[[set]]$data, consMEs.unord[[set]]$data)
  }
  
  # perform a very simple ”meta-analysis” by combining the Z scores of correlations from each set 
  # to form a meta-Z score and the corresponding p-value
  GS.metaZ = rowSums(do.call(cbind,map(GS,"Z")))/sqrt(nSets)
  kME.metaZ = rowSums(do.call(cbind,map(kME,"Z")))/sqrt(nSets)
  GS.metaP = 2*pnorm(abs(GS.metaZ), lower.tail = FALSE)
  kME.metaP = 2*pnorm(abs(kME.metaZ), lower.tail = FALSE)
  
  #form matrices holding the GS and kME
  nTrait = 6
  
  GSmat = as.data.frame(cbind(do.call(cbind,map(GS,"bicor")),
                              do.call(cbind,map(GS,"p")),
                              GS.metaZ,GS.metaP))
  ncol_GSmat = ncol(GSmat)
  colnames(GSmat)[-c(ncol_GSmat-1,ncol_GSmat)] <- paste0(rep(setLabels,each=nTrait),".",
                                                         colnames(GSmat)[-c(ncol_GSmat-1,ncol_GSmat)],
                                                         rep(c(".GS.bicor",".GS.p"),c(nSets*nTrait,nSets&nTrait)))
  
  # Same code for kME:
  kMEmat = as.data.frame(cbind(do.call(cbind,map(kME,"bicor")),
                               do.call(cbind,map(kME,"p")),
                               kME.metaZ,kME.metaP))
  
  ncol_kMEmat = ncol(kMEmat)
  n_ME = ncol(kME[[1]]$bicor)
  colnames(kMEmat)[-c(ncol_kMEmat-1,ncol_kMEmat)] <- paste0(rep(setLabels, each = n_ME),".",
                                                            colnames(kMEmat)[-c(ncol_kMEmat-1,ncol_kMEmat)],
                                                            rep(c(".bicor",".p"),
                                                                c(nSets*n_ME,nSets*n_ME)))
  
  info = data.frame(GeneSymbol = rownames(GSmat),
                    ModuleColor = moduleColors,
                    GSmat,
                    kMEmat)
  
  if(is.null(outdir)==FALSE){
    saveRDS(info,paste0(outdir,"/05_ConsensusAnalysis-CombinedNetworkResults-",region,".rds"))
  }
  return(info) 
}
get_1network_GS_kME = function(Expr,MEs,moduleColors,Traits,outdir=NULL,region=NULL){
  # get the network GS and kME for further analysis the important genes in module
  # return:
  #   the resuts include genes, module lables, GS of genes and traits, meta Result, kME of genes and  
  GS = bicorAndPvalue(Expr, Traits)
  kME = bicorAndPvalue(Expr, MEs)
  
  # perform a very simple ”meta-analysis” by combining the Z scores of correlations from each set 
  
  #form matrices holding the GS and kME
  GSmat = do.call(cbind,GS[1:2])
  colnames(GSmat) <- paste0(colnames(GSmat),".",rep(names(GS)[1:2],each=6))
  
  # Same code for kME:
  kMEmat = do.call(cbind,kME[1:2])
  colnames(kMEmat) <- paste0(colnames(kMEmat),".",rep(names(kME)[1:2],each=ncol(kME$bicor)))
  
  info = data.frame(GeneSymbol = rownames(GSmat),
                    ModuleColor = moduleColors,
                    GSmat,
                    kMEmat)
  
  if(is.null(outdir)==FALSE){
    fwrite(info,paste0(outdir,"/05_NonConsensusAnalysis-CombinedNetworkResults-",region,".csv"))
  }
  return(info) 
}
## 5.2.1 FCx ----
info_FC = get_network_GS_kME(multiExpr = multiExpr_FC,
                             moduleColors = moduleColors_FC,
                             Traits = Traits_FC,
                             nSets = nSets_FC,
                             setLabels = setLabels_FC,
                             outdir = "./2.output/9.dev_wgcna/data/",region = "FCx")
## 5.2.2 TCx ----
info_TC = get_network_GS_kME(multiExpr = multiExpr_TC,
                             moduleColors = moduleColors_TC,
                             Traits = Traits_TC,
                             nSets = nSets_TC,
                             setLabels = setLabels_TC,
                             outdir = "./2.output/9.dev_wgcna/data/",region = "TCx")

## 5.2.3 OCx ----
lnames=load(file = "2.output/9.dev_wgcna/data/02_NonConsensus_dataInput_OC_GSE25219.RData")
load(file = "2.output/9.dev_wgcna/data/03_NonConsensus-NetworkConstruction-OCx.RData")

info_OC = get_1network_GS_kME(Expr = GSE25219_OC,
                              MEs = MEs_OC,
                              moduleColors = moduleColors_OC,
                              Traits = Traits_OC,
                              outdir = "./2.output/9.dev_wgcna/data/",region = "OCx")
## 5.2.4 PCx ----
lnames=load(file = "2.output/9.dev_wgcna/data/02_NonConsensus_dataInput_PC_GSE25219.RData")
load(file = "2.output/9.dev_wgcna/data/03_NonConsensus-NetworkConstruction-PCx.RData")

info_PC = get_1network_GS_kME(Expr = GSE25219_PC,
                              MEs = MEs_PC,
                              moduleColors = moduleColors_PC,
                              Traits = Traits_PC,
                              outdir = "./2.output/9.dev_wgcna/data/",region = "PCx")
## 5.2.5 CB ----
lnames=load(file = "2.output/9.dev_wgcna/data/02_NonConsensus_dataInput_CB_GSE25219.RData")
load(file = "2.output/9.dev_wgcna/data/03_NonConsensus-NetworkConstruction-CB.RData")

info_CB = get_1network_GS_kME(Expr = GSE25219_CB,
                              MEs = MEs_CB,
                              moduleColors = moduleColors_CB,
                              Traits = Traits_CB,
                              outdir = "./2.output/9.dev_wgcna/data/",region = "CB")

## 5.2.6 HIP ----
lnames=load(file = "2.output/9.dev_wgcna/data/02_NonConsensus_dataInput_HIP_GSE25219.RData")
load(file = "2.output/9.dev_wgcna/data/03_NonConsensus-NetworkConstruction-HIP.RData")

info_HIP = get_1network_GS_kME(Expr = GSE25219_HIP,
                              MEs = MEs_HIP,
                              moduleColors = moduleColors_HIP,
                              Traits = Traits_HIP,
                              outdir = "./2.output/9.dev_wgcna/data/",region = "HIP")

## 5.2.7 STR ----
lnames=load(file = "2.output/9.dev_wgcna/data/02_NonConsensus_dataInput_STR_GSE25219.RData")
load(file = "2.output/9.dev_wgcna/data/03_NonConsensus-NetworkConstruction-STR.RData")

info_STR = get_1network_GS_kME(Expr = GSE25219_STR,
                               MEs = MEs_STR,
                               moduleColors = moduleColors_STR,
                               Traits = Traits_STR,
                               outdir = "./2.output/9.dev_wgcna/data/",region = "STR")
## 5.2.8 AMY ----
lnames=load(file = "2.output/9.dev_wgcna/data/02_NonConsensus_dataInput_AMY_GSE25219.RData")
load(file = "2.output/9.dev_wgcna/data/03_NonConsensus-NetworkConstruction-AMY.RData")

info_AMY = get_1network_GS_kME(Expr = GSE25219_AMY,
                               MEs = MEs_AMY,
                               moduleColors = moduleColors_AMY,
                               Traits = Traits_AMY,
                               outdir = "./2.output/9.dev_wgcna/data/",region = "AMY")