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
#                     1. import sample and NDD expression data
######################################################

# 1.1. import data ---
sample <- read_sample.ndd()

samplelist_brain <- sample %>% separate_rows(Super_Cohort,sep = ";") %>% split(.$brain_region3)
names(samplelist_brain) <- gsub(".*/","",names(samplelist_brain))

gselist_brain <- map(samplelist_brain,function(x){unique(x$GSE)})

dtlist_brain <- lapply(gselist_brain,get_expmat,type="ndd")
names(dtlist_brain) <- gsub(".*/","",names(dtlist_brain))

save(dtlist_brain,samplelist_brain,file="./2.output/4.ndd_wgcna/data/1_ndd_ExpmatAndSample_of_each_region.Rdata")


######################################################
#                     2. WGCNA data preparation
# batch correction -> filter genes with mad above  0 -> detect sample outlier
# FCx -> TCx -> OCx -> PCx -> CB -> HIP -> STR -> AMY
######################################################

batch_correction_bi <- function(disorder_list,disorder,data){
  
  # disorder sample information
  info = disorder_list[[disorder]]
  # filter the samples involved in data
  info = info %>% dplyr::filter(GSM %in% colnames(data))
  
  batch = info$GSE
  batch2 = info$GPL
  
  mod = model.matrix(~info$status)
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

load(file="./2.output/4.ndd_wgcna/data/1_ndd_ExpmatAndSample_of_each_region.Rdata")

## 2.1 FCx -----
### 2.1.1 preprocess data ----
dt_FC = dtlist_brain$FCx %>% column_to_rownames("gene")
sample_FC = samplelist_brain$FCx

disorder_list_FC <- sample_FC %>% split(.$Super_Cohort)
disorder_FC <- names(disorder_list_FC)

map(seq_along(disorder_list_FC),function(i){
  disorder_list_FC[[i]] %>% dplyr::filter(status=="Control" | status== names(disorder_list_FC)[[i]])
}) -> disorder_list_FC
names(disorder_list_FC) <- disorder_FC

bd_FC <- dt_FC  %>% dplyr::select(all_of(disorder_list_FC$BP$GSM))
dep_FC <- dt_FC  %>% dplyr::select(all_of(disorder_list_FC$DEP$GSM))
mdd_FC <- dt_FC  %>% dplyr::select(all_of(disorder_list_FC$MDD$GSM))
scz_FC <- dt_FC  %>% dplyr::select(all_of(disorder_list_FC$SCZ$GSM))

### batch correction
bd_FC <- batch_correction_bi(disorder_list_FC, "BP", bd_FC)
dep_FC <- batch_correction_bi(disorder_list_FC, "DEP", dep_FC)
mdd_FC <- batch_correction_bi(disorder_list_FC, "MDD", mdd_FC)
scz_FC <- batch_correction_bi(disorder_list_FC, "SCZ", scz_FC)

idx_del_FC <- unique(c(
  which(apply(bd_FC, 1, mad)==0),
  which(apply(dep_FC, 1, mad)==0),
  which(apply(mdd_FC, 1, mad)==0),
  which(apply(scz_FC, 1, mad)==0)
))

bd_FC <- bd_FC[-idx_del_FC,]
dep_FC <- dep_FC[-idx_del_FC,]
mdd_FC <- mdd_FC[-idx_del_FC,]
scz_FC <- scz_FC[-idx_del_FC,]

### 2.1.2 Form multi-set expression data ----
nSets_FC = 4
setLabels_FC = c("BD","DEP","MDD","SCZ")

multiExpr_FC = vector(mode = "list", length = nSets_FC)
multiExpr_FC[[1]] = list(data = as.data.frame(t(bd_FC)))
multiExpr_FC[[2]] = list(data = as.data.frame(t(dep_FC)))
multiExpr_FC[[3]] = list(data = as.data.frame(t(mdd_FC)))
multiExpr_FC[[4]] = list(data = as.data.frame(t(scz_FC)))
# Check that the data has the correct format for many functions operating on multiple sets:
exprSize_FC = checkSets(multiExpr_FC)
exprSize_FC

### 2.1.3 Rudimentary data cleaning and outlier removal ----
# Check that all genes and samples have sufficiently low numbers of missing values.
gsg <- goodSamplesGenesMS(multiExpr_FC,verbose = 3)
gsg$allOK

### 2.1.4 sample clustering ----
sampleTrees_FC = list()

for (set in 1:nSets_FC){
  sampleTrees_FC[[set]] = hclust(dist(multiExpr_FC[[set]]$data), method = "average")
}

# output
for (set in 1:nSets_FC) {
  pdf(file = paste0("2.output/2.ndd_wgcna/figure/FCx/02_SampleReClustering_",setLabels_FC[set],".pdf"), width = 14, height = 7)
  plot(sampleTrees_FC[[set]], xlab="", sub="", cex = 0.7,
       main = paste("Sample clustering on all genes in", setLabels_FC[[set]])) 
  dev.off()
}

# set cutheight
setLabels_FC
cutHeights = c(90,40,55,100)
# plot cluster tree
pdf(file = "2.output/2.ndd_wgcna/figure/FCx/02_SampleReClustering_cutHeights.pdf", width = 12, height = 8)
par(mfrow=c(2,2))
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
labels = cutreeStatic(sampleTrees_FC[[1]],cutHeight = cutHeights[1],minSize = 10)
table(labels)
keep = (labels==1)
multiExpr_FC[[1]]$data = multiExpr_FC[[1]]$data[keep, ]
# DEP
labels = cutreeStatic(sampleTrees_FC[[2]],cutHeight = cutHeights[2],minSize = 10)
table(labels)
keep = (labels==1)
multiExpr_FC[[2]]$data = multiExpr_FC[[2]]$data[keep, ]
# MDD
labels = cutreeStatic(sampleTrees_FC[[3]],cutHeight = cutHeights[3],minSize = 10)
table(labels)
keep = (labels==1)
multiExpr_FC[[3]]$data = multiExpr_FC[[3]]$data[keep, ]
# SCZ
labels = cutreeStatic(sampleTrees_FC[[4]],cutHeight = cutHeights[4],minSize = 10)
table(labels)
keep = (labels==1)
multiExpr_FC[[4]]$data = multiExpr_FC[[4]]$data[keep, ]

collectGarbage()
# Check the size of the leftover data
exprSize_FC = checkSets(multiExpr_FC)
exprSize_FC

### 2.1.5. import clinical information -----
sample_FC %<>% 
  dplyr::select(GSM,status) %>% 
  dplyr::mutate(status = ifelse(status=="Control",0,1)) %>% 
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

par(mar = c(0, 4, 2, 0))
for (set in 1:nSets_FC) {
  pdf(file = paste0("2.output/2.ndd_wgcna/figure/FCx/03_SampleReClustering_Trait_",setLabels_FC[set],".pdf"), width = 10, height = 5)
  traitColors = numbers2colors(Traits_FC[[set]]$data, signed = FALSE)
  plotDendroAndColors(sampleTrees2_FC[[set]], traitColors,
                      groupLabels = "Status", 
                      main = paste0("Sample dendrogram and trait heatmap-",setLabels_FC[set]))
  dev.off()
}

pdf(file = paste0("2.output/2.ndd_wgcna/figure/FCx/03_SampleReClustering_Trait_Total.pdf"), width = 10, height = 5)
for (set in 1:nSets_FC) {
  traitColors = numbers2colors(Traits_FC[[set]]$data, signed = FALSE)
  plotDendroAndColors(sampleTrees2_FC[[set]], traitColors,
                      groupLabels = "Status", 
                      main = paste0("Sample dendrogram and trait heatmap-",setLabels_FC[set]))
}
dev.off()

### 2.1.7. export exp, trait, and other info ----
save(multiExpr_FC, Traits_FC, nGenes_FC, nSamples_FC, setLabels_FC, exprSize_FC,
     file = "./2.output/2.ndd_wgcna/data/02_Consensus_dataInput_FC-re.RData")

## 2.2 TCx -----
### 2.2.1 preprocess data ----
### get dt
dt_TC = dtlist_brain$TCx %>% column_to_rownames("gene")
sample_TC = samplelist_brain$TCx

disorder_list_TC <- sample_TC %>% split(.$Super_Cohort)
disorder_TC <- names(disorder_list_TC)

map(seq_along(disorder_list_TC),function(i){
  disorder_list_TC[[i]] %>% dplyr::filter(status=="Control" | status== names(disorder_list_TC)[[i]])
}) -> disorder_list_TC
names(disorder_list_TC) <- disorder_TC

asd_TC <- dt_TC  %>% dplyr::select(all_of(disorder_list_TC$ASD$GSM))
dem_TC <- dt_TC  %>% dplyr::select(all_of(disorder_list_TC$DEM$GSM))
scz_TC <- dt_TC  %>% dplyr::select(all_of(disorder_list_TC$SCZ$GSM))

### batch correction
asd_TC <- batch_correction_bi(disorder_list_TC, "ASD", asd_TC)
dem_TC <- batch_correction_bi(disorder_list_TC, "DEM", dem_TC)
scz_TC <- batch_correction_bi(disorder_list_TC, "SCZ", scz_TC)

idx_del_TC <- unique(c(
  which(apply(asd_TC, 1, mad)==0),
  which(apply(dem_TC, 1, mad)==0),
  which(apply(scz_TC, 1, mad)==0)
))

asd_TC <- asd_TC[-idx_del_TC,]
dem_TC <- dem_TC[-idx_del_TC,]
scz_TC <- scz_TC[-idx_del_TC,]

### 2.2.2 Form multi-set expression data ----
nSets_TC = 3
setLabels_TC = c("ASD","DEM","SCZ")

multiExpr_TC = vector(mode = "list", length = nSets_TC)
multiExpr_TC[[1]] = list(data = as.data.frame(t(asd_TC)))
multiExpr_TC[[2]] = list(data = as.data.frame(t(dem_TC)))
multiExpr_TC[[3]] = list(data = as.data.frame(t(scz_TC)))
# Check that the data has the correct format for many functions operating on multiple sets:
exprSize_TC = checkSets(multiExpr_TC)
exprSize_TC

### 2.2.3 Rudimentary data cleaning and outlier removal ----
# Check that all genes and samples have sufficiently low numbers of missing values.
gsg <- goodSamplesGenesMS(multiExpr_TC,verbose = 3)
gsg$allOK

### 2.2.4 sample clustering ----
sampleTrees_TC = list()

for (set in 1:nSets_TC){
  sampleTrees_TC[[set]] = hclust(dist(multiExpr_TC[[set]]$data), method = "average")
}

# output
for (set in 1:nSets_TC) {
  pdf(file = paste0("2.output/2.ndd_wgcna/figure/TCx/02_SampleClustering_",setLabels_TC[set],".pdf"), width = 14, height = 7)
  plot(sampleTrees_TC[[set]], xlab="", sub="", cex = 0.7,
       main = paste("Sample clustering on all genes in", setLabels_TC[[set]])) 
  dev.off()
}

### 2.2.5. import clinical information -----
sample_TC %<>% 
  dplyr::select(GSM,status) %>% 
  dplyr::mutate(status = ifelse(status=="Control",0,1)) %>% 
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
  pdf(file = paste0("2.output/2.ndd_wgcna/figure/TCx/03_SampleClustering_Trait_",setLabels_TC[set],".pdf"), width = 10, height = 5)
  traitColors = numbers2colors(Traits_TC[[set]]$data, signed = FALSE)
  plotDendroAndColors(sampleTrees2_TC[[set]], traitColors,
                      groupLabels = "Status", 
                      main = paste0("Sample dendrogram and trait heatmap-",setLabels_TC[set]))
  dev.off()
}

pdf(file = paste0("2.output/2.ndd_wgcna/figure/TCx/03_SampleClustering_Trait_Total.pdf"), width = 10, height = 5)
for (set in 1:nSets_TC) {
  traitColors = numbers2colors(Traits_TC[[set]]$data, signed = FALSE)
  plotDendroAndColors(sampleTrees2_TC[[set]], traitColors,
                      groupLabels = "Status", 
                      main = paste0("Sample dendrogram and trait heatmap-",setLabels_TC[set]))
}
dev.off()

### 2.2.7. export exp, trait, and other info ----
save(multiExpr_TC, Traits_TC, nGenes_TC, nSamples_TC, setLabels_TC, exprSize_TC,
     file = "./2.output/2.ndd_wgcna/data/02_Consensus_dataInput_TC.RData")

## 2.3 OCx -----
### 2.3.1 preprocess data ----
### get dt
dt_OC = dtlist_brain$OCx %>% column_to_rownames("gene")
sample_OC = samplelist_brain$OCx

disorder_list_OC <- sample_OC %>% split(.$Super_Cohort)
disorder_OC <- names(disorder_list_OC)

map(seq_along(disorder_list_OC),function(i){
  disorder_list_OC[[i]] %>% dplyr::filter(status=="Control" | status== names(disorder_list_OC)[[i]])
}) -> disorder_list_OC
names(disorder_list_OC) <- disorder_OC

asd_OC <- dt_OC  %>% dplyr::select(all_of(disorder_list_OC$ASD$GSM))

### batch correction
asd_OC <- batch_correction_bi(disorder_list_OC, "ASD", asd_OC)

idx_del_OC <- which(apply(asd_OC, 1, mad)==0)

asd_OC <- asd_OC[-idx_del_OC,]

### 2.3.2 Form data ----
asd_OC <- as.data.frame(t(asd_OC))

### 2.3.3 Rudimentary data cleaning and outlier removal ----
# Check that all genes and samples have sufficiently low numbers of missing values.
gsg <- goodSamplesGenes(asd_OC,verbose = 3)
gsg$allOK

### 2.3.4 sample clustering ----
sampleTrees_OC_asd = hclust(dist(asd_OC), method = "average")

pdf(file = "2.output/2.ndd_wgcna/figure/OCx/02_SampleClustering_ASD_cutHeights.pdf", width = 14, height = 7)
plot(sampleTrees_OC_asd, xlab="", sub="", cex = 0.7,
     main = paste("Sample clustering on all genes in ASD")) 
abline(h=235, col = "red")
dev.off()

# Now comes the actual outlier removal:
# Find clusters cut by the line
# Keep the largest one (labeled by the number 1)
labels = cutreeStatic(sampleTrees_OC_asd,cutHeight = 235,minSize = 10)
keep = (labels==1)
asd_OC = asd_OC[keep, ]

### 2.3.5. import clinical information -----
sample_OC %<>% 
  dplyr::select(GSM,status) %>% 
  dplyr::mutate(status = ifelse(status=="Control",0,1)) %>% 
  unique() %>% 
  column_to_rownames("GSM")

setSamples = rownames(asd_OC)
traitRows = match(setSamples, rownames(sample_OC))
Traits_OC = sample_OC[traitRows, ]

# Define data set dimensions
nGenes_OC = ncol(asd_OC)
nSamples_OC = nrow(asd_OC)

### 2.3.6. reclustering samples with trait ----
sampleTrees2_OC_asd = hclust(dist(asd_OC), method = "average")

traitColors = numbers2colors(Traits_OC, signed = FALSE)

pdf(file = paste0("2.output/2.ndd_wgcna/figure/OCx/03_SampleClustering_ASD_Trait.pdf"), width = 10, height = 5)
plotDendroAndColors(sampleTrees2_OC_asd, traitColors,
                    groupLabels = "Status", 
                    main = "Sample dendrogram and trait heatmap-ASD")
dev.off()

### 2.3.7. export exp, trait, and other info ----
save(asd_OC, Traits_OC, nGenes_OC, nSamples_OC, 
     file = "./2.output/2.ndd_wgcna/data/02_NonConsensus_dataInput_OC_asd.RData")

## 2.4 PCx -----
### 2.4.1 preprocess data ----
### get dt
dt_PC = dtlist_brain$PCx %>% column_to_rownames("gene")
sample_PC = samplelist_brain$PCx

disorder_list_PC <- sample_PC %>% split(.$Super_Cohort)
disorder_PC <- names(disorder_list_PC)

map(seq_along(disorder_list_PC),function(i){
  disorder_list_PC[[i]] %>% dplyr::filter(status=="Control" | status== names(disorder_list_PC)[[i]])
}) -> disorder_list_PC
names(disorder_list_PC) <- disorder_PC

bd_PC <- dt_PC  %>% dplyr::select(all_of(disorder_list_PC$BP$GSM))
dem_PC <- dt_PC  %>% dplyr::select(all_of(disorder_list_PC$DEM$GSM))
dep_PC <- dt_PC  %>% dplyr::select(all_of(disorder_list_PC$DEP$GSM))
scz_PC <- dt_PC  %>% dplyr::select(all_of(disorder_list_PC$SCZ$GSM))

### batch correction
bd_PC <- batch_correction_bi(disorder_list_PC, "BP", bd_PC)
dep_PC <- batch_correction_bi(disorder_list_PC, "DEP", dep_PC)
dem_PC <- batch_correction_bi(disorder_list_PC, "DEM", dem_PC)
scz_PC <- batch_correction_bi(disorder_list_PC, "SCZ", scz_PC)

idx_del_PC <- unique(c(
  which(apply(bd_PC, 1, mad)==0),
  which(apply(dep_PC, 1, mad)==0),
  which(apply(dem_PC, 1, mad)==0),
  which(apply(scz_PC, 1, mad)==0)
))

if(length(idx_del_PC)>1){
  bd_PC <- bd_PC[-idx_del_PC,]
  dep_PC <- dep_PC[-idx_del_PC,]
  dem_PC <- dem_PC[-idx_del_PC,]
  scz_PC <- scz_PC[-idx_del_PC,]
}

### 2.4.2 Form multi-set expression data ----
nSets_PC = 4
setLabels_PC = c("BD","DEP","DEM","SCZ")

multiExpr_PC = vector(mode = "list", length = nSets_PC)
multiExpr_PC[[1]] = list(data = as.data.frame(t(bd_PC)))
multiExpr_PC[[2]] = list(data = as.data.frame(t(dep_PC)))
multiExpr_PC[[3]] = list(data = as.data.frame(t(dem_PC)))
multiExpr_PC[[4]] = list(data = as.data.frame(t(scz_PC)))
# Check that the data has the correct format for many functions operating on multiple sets:
exprSize_PC = checkSets(multiExpr_PC)
exprSize_PC

### 2.4.3 Rudimentary data cleaning and outlier removal ----
# Check that all genes and samples have sufficiently low numbers of missing values.
gsg <- goodSamplesGenesMS(multiExpr_PC,verbose = 3)
gsg$allOK

### 2.4.4 sample clustering ----
sampleTrees_PC = list()

for (set in 1:nSets_PC){
  sampleTrees_PC[[set]] = hclust(dist(multiExpr_PC[[set]]$data), method = "average")
}

# output
for (set in 1:nSets_PC) {
  pdf(file = paste0("2.output/2.ndd_wgcna/figure/PCx/02_SampleClustering_",setLabels_PC[set],".pdf"), width = 14, height = 7)
  plot(sampleTrees_PC[[set]], xlab="", sub="", cex = 0.7,
       main = paste("Sample clustering on all genes in", setLabels_PC[[set]])) 
  dev.off()
}

# set cutheight
cutHeights = c(65,65,70,65)
# plot cluster tree
pdf(file = "2.output/2.ndd_wgcna/figure/PCx/02_SampleClustering_cutHeights.pdf", width = 12, height = 8)
par(mfrow=c(2,2))
par(mar = c(0, 4, 2, 0))
for (set in 1:nSets_PC){
  plot(sampleTrees_PC[[set]], xlab="", sub="", cex = 0.7,
       main = paste("Sample clustering on all genes in", setLabels_PC[set])) 
  abline(h=cutHeights[set], col = "red")
}
dev.off()

# Now comes the actual outlier removal:
# Find clusters cut by the line
# Keep the largest one (labeled by the number 1)
setLabels_PC
# BD
labels = cutreeStatic(sampleTrees_PC[[1]],cutHeight = cutHeights[1],minSize = 10)
table(labels)
keep = (labels==1)
multiExpr_PC[[1]]$data = multiExpr_PC[[1]]$data[keep, ]
# DEP
labels = cutreeStatic(sampleTrees_PC[[2]],cutHeight = cutHeights[2],minSize = 10)
table(labels)
keep = (labels==1)
multiExpr_PC[[2]]$data = multiExpr_PC[[2]]$data[keep, ]
# DEM
labels = cutreeStatic(sampleTrees_PC[[3]],cutHeight = cutHeights[3],minSize = 10)
table(labels)
keep = (labels==1)
multiExpr_PC[[3]]$data = multiExpr_PC[[3]]$data[keep, ]
# SCZ
labels = cutreeStatic(sampleTrees_PC[[4]],cutHeight = cutHeights[4],minSize = 10)
keep = (labels==1)
multiExpr_PC[[4]]$data = multiExpr_PC[[4]]$data[keep, ]

# Check the size of the leftover data
exprSize_PC = checkSets(multiExpr_PC)
exprSize_PC

### 2.4.5. import clinical information -----
sample_PC %<>% 
  dplyr::select(GSM,status) %>% 
  dplyr::mutate(status = ifelse(status=="Control",0,1)) %>% 
  unique() %>% 
  column_to_rownames("GSM")

Traits_PC = vector(mode="list", length = nSets_PC)
for (set in 1:nSets_PC)
{
  setSamples = rownames(multiExpr_PC[[set]]$data)
  traitRows = match(setSamples, rownames(sample_PC))
  Traits_PC[[set]] = list(data = sample_PC[traitRows, ]) 
}

# Define data set dimensions
nGenes_PC = exprSize_PC$nGenes
nSamples_PC = exprSize_PC$nSamples

### 2.4.6. reclustering samples with trait ----
sampleTrees2_PC = list()

for (set in 1:nSets_PC){
  sampleTrees2_PC[[set]] = hclust(dist(multiExpr_PC[[set]]$data), method = "average")
}

par(mar = c(0, 4, 2, 0))
for (set in 1:nSets_PC) {
  pdf(file = paste0("2.output/2.ndd_wgcna/figure/PCx/03_SampleClustering_Trait_",setLabels_PC[set],".pdf"), width = 10, height = 5)
  traitColors = numbers2colors(Traits_PC[[set]]$data, signed = FALSE)
  plotDendroAndColors(sampleTrees2_PC[[set]], traitColors,
                      groupLabels = "Status", 
                      main = paste0("Sample dendrogram and trait heatmap-",setLabels_PC[set]))
  dev.off()
  
}

pdf(file = paste0("2.output/2.ndd_wgcna/figure/PCx/03_SampleClustering_Trait_Total.pdf"), width = 10, height = 5)
for (set in 1:nSets_PC) {
  traitColors = numbers2colors(Traits_PC[[set]]$data, signed = FALSE)
  plotDendroAndColors(sampleTrees2_PC[[set]], traitColors,
                      groupLabels = "Status", 
                      main = paste0("Sample dendrogram and trait heatmap-",setLabels_PC[set]))
}
dev.off()

### 2.4.7. export exp, trait, and other info ----
save(multiExpr_PC, Traits_PC, nGenes_PC, nSamples_PC, setLabels_PC, exprSize_PC,
     file = "./2.output/2.ndd_wgcna/data/02_Consensus_dataInput_PC.RData")

## 2.5 CB -----
### 2.5.1 preprocess data ----
### get dt
dt_CB = dtlist_brain$CB %>% column_to_rownames("gene")
sample_CB = samplelist_brain$CB

disorder_list_CB <- sample_CB %>% split(.$Super_Cohort)
disorder_CB <- names(disorder_list_CB)

map(seq_along(disorder_list_CB),function(i){
  disorder_list_CB[[i]] %>% dplyr::filter(status=="Control" | status== names(disorder_list_CB)[[i]])
}) -> disorder_list_CB
names(disorder_list_CB) <- disorder_CB

bd_CB <- dt_CB  %>% dplyr::select(all_of(disorder_list_CB$BP$GSM))
dep_CB <- dt_CB  %>% dplyr::select(all_of(disorder_list_CB$DEP$GSM))
scz_CB <- dt_CB  %>% dplyr::select(all_of(disorder_list_CB$SCZ$GSM))

### batch correction (cuz of high heterogeneous between datasets, divide asd into two datasets)
bd_CB <- batch_correction_bi(disorder_list_CB, "BP", bd_CB)
dep_CB <- batch_correction_bi(disorder_list_CB, "DEP", dep_CB)
scz_CB <- batch_correction_bi(disorder_list_CB, "SCZ", scz_CB)

idx_del_CB <- unique(c(
  which(apply(bd_CB, 1, mad)==0),
  which(apply(dep_CB, 1, mad)==0),
  which(apply(scz_CB, 1, mad)==0)
))

if(length(idx_del_CB)>1){
  bd_CB <- bd_CB[-idx_del_CB,]
  dep_CB <- dep_CB[-idx_del_CB,]
  scz_CB <- scz_CB[-idx_del_CB,]
}

### 2.5.2 Form multi-set expression data ----
nSets_CB = 3
setLabels_CB = c("BD","DEP","SCZ")

multiExpr_CB = vector(mode = "list", length = nSets_CB)
multiExpr_CB[[1]] = list(data = as.data.frame(t(bd_CB)))
multiExpr_CB[[2]] = list(data = as.data.frame(t(dep_CB)))
multiExpr_CB[[3]] = list(data = as.data.frame(t(scz_CB)))
# Check that the data has the correct format for many functions operating on multiple sets:
exprSize_CB = checkSets(multiExpr_CB)
exprSize_CB

### 2.5.3 Rudimentary data cleaning and outlier removal ----
# Check that all genes and samples have sufficiently low numbers of missing values.
gsg <- goodSamplesGenesMS(multiExpr_CB,verbose = 3)
gsg$allOK

### 2.5.4 sample clustering ----
sampleTrees_CB = list()

for (set in 1:nSets_CB){
  sampleTrees_CB[[set]] = hclust(dist(multiExpr_CB[[set]]$data), method = "average")
}

# output
for (set in 1:nSets_CB) {
  pdf(file = paste0("2.output/2.ndd_wgcna/figure/CB/02_ReSampleClustering_",setLabels_CB[set],".pdf"), width = 14, height = 7)
  plot(sampleTrees_CB[[set]], xlab="", sub="", cex = 0.7,
       main = paste("Sample clustering on all genes in", setLabels_CB[[set]])) 
  dev.off()
}

# set cutheight
setLabels_CB
cutHeights = c(60,50,50)
# plot cluster tree
pdf(file = "2.output/2.ndd_wgcna/figure/CB/02_ReSampleClustering_cutHeights.pdf", width = 12, height = 8)
par(mfrow=c(2,2))
par(mar = c(0, 4, 2, 0))
for (set in 1:nSets_CB){
  plot(sampleTrees_CB[[set]], xlab="", sub="", cex = 0.7,
       main = paste("Sample clustering on all genes in", setLabels_CB[set])) 
  abline(h=cutHeights[set], col = "red")
}
dev.off()

# Now comes the actual outlier removal:
# Find clusters cut by the line
# Keep the largest one (labeled by the number 1)
# BD
labels = cutreeStatic(sampleTrees_CB[[1]],cutHeight = cutHeights[1],minSize = 10)
table(labels)
keep = (labels==1)
multiExpr_CB[[1]]$data = multiExpr_CB[[1]]$data[keep, ]
# DEP
labels = cutreeStatic(sampleTrees_CB[[2]],cutHeight = cutHeights[2],minSize = 10)
table(labels)
keep = (labels==1)
multiExpr_CB[[2]]$data = multiExpr_CB[[2]]$data[keep, ]
# SCZ
labels = cutreeStatic(sampleTrees_CB[[3]],cutHeight = cutHeights[3],minSize = 10)
table(labels)
keep = (labels==1)
multiExpr_CB[[3]]$data = multiExpr_CB[[3]]$data[keep, ]

collectGarbage()
# Check the size of the leftover data
exprSize_CB = checkSets(multiExpr_CB)
exprSize_CB

### 2.5.5. import clinical information -----
sample_CB %<>% 
  dplyr::select(GSM,status) %>% 
  dplyr::mutate(status = ifelse(status=="Control",0,1)) %>% 
  unique() %>% 
  column_to_rownames("GSM")

Traits_CB = vector(mode="list", length = nSets_CB)
for (set in 1:nSets_CB)
{
  setSamples = rownames(multiExpr_CB[[set]]$data)
  traitRows = match(setSamples, rownames(sample_CB))
  Traits_CB[[set]] = list(data = sample_CB[traitRows, ]) 
}

# Define data set dimensions
nGenes_CB = exprSize_CB$nGenes
nSamples_CB = exprSize_CB$nSamples

### 2.5.6. reclustering samples with trait ----
sampleTrees2_CB = list()

for (set in 1:nSets_CB){
  sampleTrees2_CB[[set]] = hclust(dist(multiExpr_CB[[set]]$data), method = "average")
}

par(mar = c(0, 4, 2, 0))
for (set in 1:nSets_CB) {
  pdf(file = paste0("2.output/2.ndd_wgcna/figure/CB/03_ReSampleClustering_Trait_",setLabels_CB[set],".pdf"), width = 10, height = 5)
  traitColors = numbers2colors(Traits_CB[[set]]$data, signed = FALSE)
  plotDendroAndColors(sampleTrees2_CB[[set]], traitColors,
                      groupLabels = "Status", 
                      main = paste0("Sample dendrogram and trait heatmap-",setLabels_CB[set]))
  dev.off()
  
}

pdf(file = paste0("2.output/2.ndd_wgcna/figure/CB/03_ReSampleClustering_Trait_Total.pdf"), width = 10, height = 5)
for (set in 1:nSets_CB) {
  traitColors = numbers2colors(Traits_CB[[set]]$data, signed = FALSE)
  plotDendroAndColors(sampleTrees2_CB[[set]], traitColors,
                      groupLabels = "Status", 
                      main = paste0("Sample dendrogram and trait heatmap-",setLabels_CB[set]))
}
dev.off()

### 2.5.7. export exp, trait, and other info ----
save(multiExpr_CB, Traits_CB, nGenes_CB, nSamples_CB, setLabels_CB, exprSize_CB,
     file = "./2.output/2.ndd_wgcna/data/02_Consensus_dataInput_CB-re.RData")

## 2.6 HIP -----
### 2.6.1 preprocess data ----
### get dt
dt_HIP = dtlist_brain$HIP %>% column_to_rownames("gene")
sample_HIP = samplelist_brain$HIP

disorder_list_HIP <- sample_HIP %>% split(.$Super_Cohort)
disorder_HIP <- names(disorder_list_HIP)

map(seq_along(disorder_list_HIP),function(i){
  disorder_list_HIP[[i]] %>% dplyr::filter(status=="Control" | status== names(disorder_list_HIP)[[i]])
}) -> disorder_list_HIP
names(disorder_list_HIP) <- disorder_HIP

bd_HIP <- dt_HIP  %>% dplyr::select(all_of(disorder_list_HIP$BP$GSM))
dem_HIP <- dt_HIP  %>% dplyr::select(all_of(disorder_list_HIP$DEM$GSM))
mdd_HIP <- dt_HIP  %>% dplyr::select(all_of(disorder_list_HIP$MDD$GSM))
scz_HIP <- dt_HIP  %>% dplyr::select(all_of(disorder_list_HIP$SCZ$GSM))

### batch correction (cuz of high heterogeneous between datasets, divide asd into two datasets)
bd_HIP <- batch_correction_bi(disorder_list_HIP, "BP", bd_HIP)
dem_HIP <- batch_correction_bi(disorder_list_HIP, "DEM", dem_HIP)
mdd_HIP <- batch_correction_bi(disorder_list_HIP, "MDD", mdd_HIP)
scz_HIP <- batch_correction_bi(disorder_list_HIP, "SCZ", scz_HIP)

idx_del_HIP <- unique(c(
  which(apply(bd_HIP, 1, mad)==0),
  which(apply(dem_HIP, 1, mad)==0),
  which(apply(mdd_HIP, 1, mad)==0),
  which(apply(scz_HIP, 1, mad)==0)
))

if(length(idx_del_HIP)>1){
  bd_HIP <- bd_HIP[-idx_del_HIP,]
  dem_HIP <- dem_HIP[-idx_del_HIP,]
  mdd_HIP <- mdd_HIP[-idx_del_HIP,]
  scz_HIP <- scz_HIP[-idx_del_HIP,]
}

### 2.6.2 Form multi-set expression data ----
nSets_HIP = 4
setLabels_HIP = c("BD","DEM","MDD","SCZ")

multiExpr_HIP = vector(mode = "list", length = nSets_HIP)
multiExpr_HIP[[1]] = list(data = as.data.frame(t(bd_HIP)))
multiExpr_HIP[[2]] = list(data = as.data.frame(t(dem_HIP)))
multiExpr_HIP[[3]] = list(data = as.data.frame(t(mdd_HIP)))
multiExpr_HIP[[4]] = list(data = as.data.frame(t(scz_HIP)))
# Check that the data has the correct format for many functions operating on multiple sets:
exprSize_HIP = checkSets(multiExpr_HIP)
exprSize_HIP

### 2.6.3 Rudimentary data cleaning and outlier removal ----
# Check that all genes and samples have sufficiently low numbers of missing values.
gsg <- goodSamplesGenesMS(multiExpr_HIP,verbose = 3)
gsg$allOK

### 2.6.4 sample clustering ----
sampleTrees_HIP = list()

for (set in 1:nSets_HIP){
  sampleTrees_HIP[[set]] = hclust(dist(multiExpr_HIP[[set]]$data), method = "average")
}

# output
for (set in 1:nSets_HIP) {
  pdf(file = paste0("2.output/2.ndd_wgcna/figure/HIP/02_SampleClustering_",setLabels_HIP[set],".pdf"), width = 14, height = 7)
  plot(sampleTrees_HIP[[set]], xlab="", sub="", cex = 0.7,
       main = paste("Sample clustering on all genes in", setLabels_HIP[[set]])) 
  dev.off()
}

# set cutheight
setLabels_HIP
cutHeights = c(50,70,70,80)
# plot cluster tree
pdf(file = "2.output/2.ndd_wgcna/figure/HIP/02_SampleClustering_cutHeights.pdf", width = 12, height = 8)
par(mfrow=c(2,2))
par(mar = c(0, 4, 2, 0))
for (set in 1:nSets_HIP){
  plot(sampleTrees_HIP[[set]], xlab="", sub="", cex = 0.7,
       main = paste("Sample clustering on all genes in", setLabels_HIP[set])) 
  abline(h=cutHeights[set], col = "red")
}
dev.off()

# Now comes the actual outlier removal:
# Find clusters cut by the line
# Keep the largest one (labeled by the number 1)
# bd
labels = cutreeStatic(sampleTrees_HIP[[1]],cutHeight = cutHeights[1],minSize = 10)
table(labels)
keep = (labels==1)
multiExpr_HIP[[1]]$data = multiExpr_HIP[[1]]$data[keep, ]
# dem
labels = cutreeStatic(sampleTrees_HIP[[2]],cutHeight = cutHeights[2],minSize = 10)
table(labels)
keep = (labels==1)
multiExpr_HIP[[2]]$data = multiExpr_HIP[[2]]$data[keep, ]
# mdd
labels = cutreeStatic(sampleTrees_HIP[[3]],cutHeight = cutHeights[3],minSize = 10)
table(labels)
keep = (labels==1)
multiExpr_HIP[[3]]$data = multiExpr_HIP[[3]]$data[keep, ]
# SCZ
labels = cutreeStatic(sampleTrees_HIP[[4]],cutHeight = cutHeights[4],minSize = 10)
table(labels)
keep = (labels==1)
multiExpr_HIP[[4]]$data = multiExpr_HIP[[4]]$data[keep, ]

collectGarbage()
# Check the size of the leftover data
exprSize_HIP = checkSets(multiExpr_HIP)
exprSize_HIP

### 2.6.5. import clinical information -----
sample_HIP %<>% 
  dplyr::select(GSM,status) %>% 
  dplyr::mutate(status = ifelse(status=="Control",0,1)) %>% 
  unique() %>% 
  column_to_rownames("GSM")

Traits_HIP = vector(mode="list", length = nSets_HIP)
for (set in 1:nSets_HIP)
{
  setSamples = rownames(multiExpr_HIP[[set]]$data)
  traitRows = match(setSamples, rownames(sample_HIP))
  Traits_HIP[[set]] = list(data = sample_HIP[traitRows, ]) 
}

# Define data set dimensions
nGenes_HIP = exprSize_HIP$nGenes
nSamples_HIP = exprSize_HIP$nSamples

### 2.6.6. reclustering samples with trait ----
sampleTrees2_HIP = list()

for (set in 1:nSets_HIP){
  sampleTrees2_HIP[[set]] = hclust(dist(multiExpr_HIP[[set]]$data), method = "average")
}

par(mar = c(0, 4, 2, 0))
for (set in 1:nSets_HIP) {
  pdf(file = paste0("2.output/2.ndd_wgcna/figure/HIP/03_SampleClustering_Trait_",setLabels_HIP[set],".pdf"), width = 10, height = 5)
  traitColors = numbers2colors(Traits_HIP[[set]]$data, signed = FALSE)
  plotDendroAndColors(sampleTrees2_HIP[[set]], traitColors,
                      groupLabels = "Status", 
                      main = paste0("Sample dendrogram and trait heatmap-",setLabels_HIP[set]))
  dev.off()
  
}

pdf(file = paste0("2.output/2.ndd_wgcna/figure/HIP/03_SampleClustering_Trait_Total.pdf"), width = 10, height = 5)
for (set in 1:nSets_HIP) {
  traitColors = numbers2colors(Traits_HIP[[set]]$data, signed = FALSE)
  plotDendroAndColors(sampleTrees2_HIP[[set]], traitColors,
                      groupLabels = "Status", 
                      main = paste0("Sample dendrogram and trait heatmap-",setLabels_HIP[set]))
}
dev.off()

### 2.6.7. export exp, trait, and other info ----
save(multiExpr_HIP, Traits_HIP, nGenes_HIP, nSamples_HIP, setLabels_HIP, exprSize_HIP,
     file = "./2.output/2.ndd_wgcna/data/02_Consensus_dataInput_HIP.RData")

## 2.7 STR -----
### 2.7.1 preprocess data ----
### get dt
dt_STR = dtlist_brain$STR %>% column_to_rownames("gene")
sample_STR = samplelist_brain$STR

disorder_list_STR <- sample_STR %>% split(.$Super_Cohort)
disorder_STR <- names(disorder_list_STR)

map(seq_along(disorder_list_STR),function(i){
  disorder_list_STR[[i]] %>% dplyr::filter(status=="Control" | status== names(disorder_list_STR)[[i]])
}) -> disorder_list_STR
names(disorder_list_STR) <- disorder_STR

bd_STR <- dt_STR  %>% dplyr::select(all_of(disorder_list_STR$BP$GSM))
mdd_STR <- dt_STR  %>% dplyr::select(all_of(disorder_list_STR$MDD$GSM))
scz_STR <- dt_STR  %>% dplyr::select(all_of(disorder_list_STR$SCZ$GSM))

bd_STR <- batch_correction_bi(disorder_list_STR, "BP", bd_STR)
mdd_STR <- batch_correction_bi(disorder_list_STR, "MDD", mdd_STR)
scz_STR <- batch_correction_bi(disorder_list_STR, "SCZ", scz_STR)

idx_del_STR <- unique(c(
  which(apply(bd_STR, 1, mad)==0),
  which(apply(mdd_STR, 1, mad)==0),
  which(apply(scz_STR, 1, mad)==0)
))

if(length(idx_del_STR)>1){
  bd_STR <- bd_STR[-idx_del_STR,]
  mdd_STR <- mdd_STR[-idx_del_STR,]
  scz_STR <- scz_STR[-idx_del_STR,]
}

### 2.7.2 Form multi-set expression data ----
nSets_STR = 3
setLabels_STR = c("BD","MDD","SCZ")

multiExpr_STR = vector(mode = "list", length = nSets_STR)
multiExpr_STR[[1]] = list(data = as.data.frame(t(bd_STR)))
multiExpr_STR[[2]] = list(data = as.data.frame(t(mdd_STR)))
multiExpr_STR[[3]] = list(data = as.data.frame(t(scz_STR)))
# Check that the data has the correct format for many functions operating on multiple sets:
exprSize_STR = checkSets(multiExpr_STR)
exprSize_STR
### 2.7.3 Rudimentary data cleaning and outlier removal ----
# Check that all genes and samples have sufficiently low numbers of missing values.
gsg <- goodSamplesGenesMS(multiExpr_STR,verbose = 3)
gsg$allOK

### 2.7.4 sample clustering ----
sampleTrees_STR = list()

for (set in 1:nSets_STR){
  sampleTrees_STR[[set]] = hclust(dist(multiExpr_STR[[set]]$data), method = "average")
}

# output
for (set in 1:nSets_STR) {
  pdf(file = paste0("2.output/4.ndd_wgcna/figure/STR/02_ReSampleClustering_",setLabels_STR[set],".pdf"), width = 14, height = 7)
  plot(sampleTrees_STR[[set]], xlab="", sub="", cex = 0.7,
       main = paste("Sample clustering on all genes in", setLabels_STR[[set]])) 
  dev.off()
}

# set cutheight
setLabels_STR
cutHeights = c(700,60,80)
# plot cluster tree
pdf(file = "2.output/4.ndd_wgcna/figure/STR/02_ReSampleClustering_cutHeights.pdf", width = 12, height = 8)
par(mfrow=c(2,2))
par(mar = c(0, 4, 2, 0))
for (set in 1:nSets_STR){
  plot(sampleTrees_STR[[set]], xlab="", sub="", cex = 0.7,
       main = paste("Sample clustering on all genes in", setLabels_STR[set])) 
  abline(h=cutHeights[set], col = "red")
}
dev.off()

# Now comes the actual outlier removal:
# Find clusters cut by the line
# Keep the largest one (labeled by the number 1)
# BD
labels = cutreeStatic(sampleTrees_STR[[1]],cutHeight = cutHeights[1],minSize = 10)
table(labels)
keep = (labels==1)
multiExpr_STR[[1]]$data = multiExpr_STR[[1]]$data[keep, ]
# MDD
labels = cutreeStatic(sampleTrees_STR[[2]],cutHeight = cutHeights[2],minSize = 10)
table(labels)
keep = (labels==1)
multiExpr_STR[[2]]$data = multiExpr_STR[[2]]$data[keep, ]
# SCZ
labels = cutreeStatic(sampleTrees_STR[[3]],cutHeight = cutHeights[3],minSize = 10)
table(labels)
keep = (labels==1)
multiExpr_STR[[3]]$data = multiExpr_STR[[3]]$data[keep, ]

collectGarbage()
# Check the size of the leftover data
exprSize_STR = checkSets(multiExpr_STR)
exprSize_STR

### 2.7.5. import clinical information -----
sample_STR %<>% 
  dplyr::select(GSM,status) %>% 
  dplyr::mutate(status = ifelse(status=="Control",0,1)) %>% 
  unique() %>% 
  column_to_rownames("GSM")

Traits_STR = vector(mode="list", length = nSets_STR)
for (set in 1:nSets_STR)
{
  setSamples = rownames(multiExpr_STR[[set]]$data)
  traitRows = match(setSamples, rownames(sample_STR))
  Traits_STR[[set]] = list(data = sample_STR[traitRows, ]) 
}

# Define data set dimensions
nGenes_STR = exprSize_STR$nGenes
nSamples_STR = exprSize_STR$nSamples

### 2.7.6. reclustering samples with trait ----
sampleTrees2_STR = list()

for (set in 1:nSets_STR){
  sampleTrees2_STR[[set]] = hclust(dist(multiExpr_STR[[set]]$data), method = "average")
}

par(mar = c(0, 4, 2, 0))
for (set in 1:nSets_STR) {
  pdf(file = paste0("2.output/4.ndd_wgcna/figure/STR/03_ReSampleClustering_Trait_",setLabels_STR[set],".pdf"), width = 10, height = 5)
  traitColors = numbers2colors(Traits_STR[[set]]$data, signed = FALSE)
  plotDendroAndColors(sampleTrees2_STR[[set]], traitColors,
                      groupLabels = "Status", 
                      main = paste0("Sample dendrogram and trait heatmap-",setLabels_STR[set]))
  dev.off()
}
disorder_list_STR$BP %>% View()

pdf(file = paste0("2.output/4.ndd_wgcna/figure/STR/03_ReSampleClustering_Trait_Total.pdf"), width = 10, height = 5)
for (set in 1:nSets_STR) {
  traitColors = numbers2colors(Traits_STR[[set]]$data, signed = FALSE)
  plotDendroAndColors(sampleTrees2_STR[[set]], traitColors,
                      groupLabels = "Status", 
                      main = paste0("Sample dendrogram and trait heatmap-",setLabels_STR[set]))
}
dev.off()

### 2.7.7. export exp, trait, and other info ----
save(multiExpr_STR, Traits_STR, nGenes_STR, nSamples_STR, setLabels_STR, exprSize_STR,
     file = "./2.output/4.ndd_wgcna/data/02_Consensus_dataInput_STR.RData")

## 2.8 AMY -----
### 2.8.1 preprocess data ----
### get dt
dt_AMY = dtlist_brain$AMY %>% column_to_rownames("gene")
sample_AMY = samplelist_brain$AMY

disorder_list_AMY <- sample_AMY %>% split(.$Super_Cohort)
disorder_AMY <- names(disorder_list_AMY)

map(seq_along(disorder_list_AMY),function(i){
  disorder_list_AMY[[i]] %>% dplyr::filter(status=="Control" | status== names(disorder_list_AMY)[[i]])
}) -> disorder_list_AMY
names(disorder_list_AMY) <- disorder_AMY

mdd_AMY <- dt_AMY  %>% dplyr::select(all_of(disorder_list_AMY$MDD$GSM))

### batch correction
mdd_AMY <- batch_correction_bi(disorder_list_AMY, "MDD", mdd_AMY)

idx_del_AMY <- which(apply(mdd_AMY, 1, mad)==0)

#mdd_AMY <- mdd_AMY[-idx_del_AMY,]

### 2.8.2 Form data ----
mdd_AMY <- as.data.frame(t(mdd_AMY))

### 2.8.3 Rudimentary data cleaning and outlier removal ----
# Check that all genes and samples have sufficiently low numbers of missing values.
gsg <- goodSamplesGenes(mdd_AMY,verbose = 3)
gsg$allOK

### 2.8.4 sample clustering ----
sampleTrees_AMY_mdd = hclust(dist(mdd_AMY), method = "average")

pdf(file = "2.output/2.ndd_wgcna/figure/AMY/02_SampleClustering_MDD_cutHeights.pdf", width = 14, height = 7)
plot(sampleTrees_AMY_mdd, xlab="", sub="", cex = 0.7,
     main = paste("Sample clustering on all genes in MDD")) 
abline(h=42, col = "red")
dev.off()

# Now comes the actual outlier removal:
# Find clusters cut by the line
# Keep the largest one (labeled by the number 1)
labels = cutreeStatic(sampleTrees_AMY_mdd,cutHeight = 42,minSize = 10)
keep = (labels==1)
mdd_AMY = mdd_AMY[keep, ]

### 2.8.5. import clinical information -----
sample_AMY %<>% 
  dplyr::select(GSM,status) %>% 
  dplyr::mutate(status = ifelse(status=="Control",0,1)) %>% 
  unique() %>% 
  column_to_rownames("GSM")

setSamples = rownames(mdd_AMY)
traitRows = match(setSamples, rownames(sample_AMY))
Traits_AMY = sample_AMY[traitRows, ]

# Define data set dimensions
nGenes_AMY = ncol(mdd_AMY)
nSamples_AMY = nrow(mdd_AMY)

### 2.8.6. reclustering samples with trait ----
sampleTrees2_AMY_mdd = hclust(dist(mdd_AMY), method = "average")

traitColors = numbers2colors(Traits_AMY, signed = FALSE)

pdf(file = paste0("2.output/2.ndd_wgcna/figure/AMY/03_SampleClustering_MDD_Trait.pdf"), width = 12, height = 5)
plotDendroAndColors(sampleTrees2_AMY_mdd, traitColors,
                    groupLabels = "Status", 
                    main = "Sample dendrogram and trait heatmap-MDD")
dev.off()

### 2.8.7. export exp, trait, and other info ----
save(mdd_AMY, Traits_AMY, nGenes_AMY, nSamples_AMY, 
     file = "./2.output/2.ndd_wgcna/data/02_NonConsensus_dataInput_AMY_mdd.RData")

######################################################
#                     3. Network construction and module detection
# Coexpression networks were constructed in each brain region
# select power -> network construction
######################################################

powers = c(seq(4,10,by=1), seq(12,30, by=2)) 
# 3.1 FCx ----
lnames=load("./2.output/2.ndd_wgcna/data/02_Consensus_dataInput_FC-re.RData")

## 3.1.1 Choosing the soft-thresholding power： analysis of network topology ----
# delete asd casue of high connectivity (temp)
nSets_FC = checkSets(multiExpr_FC)$nSets

# Initialize a list to hold the results of scale-free analysis
powerTables_FC = vector(mode = "list", length = nSets_FC)
# Call the network topology analysis function for each set in turn
for (set in 1:nSets_FC)
{
  powerTables_FC[[set]] = list(data = pickSoftThreshold(multiExpr_FC[[set]]$data, 
                                                     powerVector=powers,corFnc ="bicor",networkType = "signed",
                                                     verbose = 2)[[2]])  
}

# Plot the results:
colors = RColorBrewer::brewer.pal(nSets_FC,"Set2")

# Will plot these columns of the returned scale free analysis tables
plotCols = c(2,5,6,7)
colNames = c("Scale Free Topology Model Fit", "Mean connectivity", "Median connectivity",
             "Max connectivity")

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

pdf("./2.output/2.ndd_wgcna/figure/FCx/04_ReSummary network indices (y-axes) as functions of the soft thresholding power (x-axes).pdf")
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
legend("center", legend = setLabels_FC, col = colors, pch = 20,cex = 1.5,ncol = nSets_FC)

dev.off()

### 3.1.2 Network construction and consensus module detection ----
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

pdf("./2.output/2.ndd_wgcna/figure/FCx/05_Dendrogram of consensus module eigengenes obtained by clustering the eigengenes on their consensus correlation.pdf",width=7, height=4)
plotDendroAndColors(consTree_FC, moduleColors_FC,
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Consensus gene dendrogram and module colors")
dev.off()
save(consMEs_FC, moduleColors_FC, consTree_FC, file = "./2.output/2.ndd_wgcna/data/03_Consensus-NetworkConstruction-FCx.RData")

## 3.2 TCx ----
powers = c(seq(4,10,by=1), seq(12,30, by=2)) 
lnames=load("./2.output/2.ndd_wgcna/data/02_Consensus_dataInput_TC.RData")
### 3.2.1 Choosing the soft-thresholding power： analysis of network topology ----
# delete asd casue of high connectivity (temp)
nSets_TC = checkSets(multiExpr_TC)$nSets

# Initialize a list to hold the results of scale-free analysis
powerTables_TC = vector(mode = "list", length = nSets_TC)
# Call the network topology analysis function for each set in turn
for (set in 1:nSets_TC)
{
  powerTables_TC[[set]] = list(data = pickSoftThreshold(multiExpr_TC[[set]]$data, 
                                                        powerVector=powers,corFnc ="bicor",networkType = "signed",
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

pdf("./2.output/2.ndd_wgcna/figure/TCx/04_Summary network indices (y-axes) as functions of the soft thresholding power (x-axes).pdf")

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
legend("center", legend = setLabels_TC, col = colors, pch = 20,cex = 1.5,ncol = nSets_TC)

dev.off()

### 3.2.2 Network construction and consensus module detection ----
setLabels_TC
softPower_TC = c(9,18,5)
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

pdf("./2.output/2.ndd_wgcna/figure/TCx/05_Dendrogram of consensus module eigengenes obtained by clustering the eigengenes on their consensus correlation.pdf",width=7, height=4)
plotDendroAndColors(consTree_TC, moduleColors_TC,
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Consensus gene dendrogram and module colors")
dev.off()
save(consMEs_TC, moduleColors_TC, consTree_TC, file = "./2.output/2.ndd_wgcna/data/03_Consensus-NetworkConstruction-TCx.RData")

## 3.3 OCx ----
powers = c(c(1:10), seq(from = 12, to=20, by=2))
lnames=load("./2.output/2.ndd_wgcna/data/02_NonConsensus_dataInput_OC_asd.RData")
### 3.3.1 Choosing the soft-thresholding power： analysis of network topology ----
sft = pickSoftThreshold(asd_OC, powerVector = powers, verbose = 5)
                       
# Plot the quantities in the chosen columns vs. the soft thresholding power
pdf("./2.output/2.ndd_wgcna/figure/OCx/04_Summary network indices (y-axes) as functions of the soft thresholding power (x-axes).pdf",width = 6.5,height = 5.5)
#par(mOCol = c(2,2))
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
net_OC = blockwiseModules(asd_OC, 
                          power = 5,
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

pdf("./2.output/2.ndd_wgcna/figure/OCx/05_Dendrogram of genes with dissimilarity based on topological overlap together with assigned module colors.pdf",width=7, height=4)
plotDendroAndColors(Tree_OC, moduleColors_OC,
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Consensus gene dendrogram and module colors")
dev.off()
save(MEs_OC, moduleColors_OC, Tree_OC, file = "./2.output/2.ndd_wgcna/data/03_NonConsensus-NetworkConstruction-OCx.RData")

## 3.4 PCx ----
powers = c(seq(4,10,by=1), seq(12,30, by=2)) 
lnames=load("./2.output/2.ndd_wgcna/data/02_Consensus_dataInput_PC.RData")
### 3.4.1 Choosing the soft-thresholding power： analysis of network topology ----
# delete asd casue of high connectivity (temp)
nSets_PC = checkSets(multiExpr_PC)$nSets

# Initialize a list to hold the results of scale-free analysis
powerTables_PC = vector(mode = "list", length = nSets_PC)
# Call the network topology analysis function for each set in turn
for (set in 1:nSets_PC)
{
  powerTables_PC[[set]] = list(data = pickSoftThreshold(multiExpr_PC[[set]]$data, 
                                                        powerVector=powers,corFnc ="bicor",networkType = "signed",
                                                        verbose = 2)[[2]])  
}

# Plot the results:
colors = RColorBrewer::brewer.pal(nSets_PC,"Set2")

# Will plot these columns of the returned scale free analysis tables
ploPCols = c(2,5,6,7)
colNames = c("Scale Free Topology Model Fit", "Mean connectivity", "Median connectivity",
             "Max connectivity")

# Get the minima and maxima of the plotted points
ylim = matrix(NA, nrow = nSets_PC, ncol = 4)
for (set in 1:nSets_PC)
{
  for (col in 1:length(ploPCols))
  {
    ylim[set, col] = min(ylim[set, col], powerTables_PC[[set]]$data[, ploPCols[col]], na.rm = TRUE)
  }
}

# Plot the quantities in the chosen columns vs. the soft thresholding power
sizeGrWindow(8, 6)

pdf("./2.output/2.ndd_wgcna/figure/PCx/04_Summary network indices (y-axes) as functions of the soft thresholding power (x-axes).pdf")
#par(mPCol = c(2,2))
layout(matrix(c(1,2,3,4,5,5),ncol = 2,byrow = T),heights = c(3,3,1))
#par(mar = c(4.2, 4.2 , 2.2, 0.5))
cex1 = 1

for (col in 1:length(ploPCols))
{
  # set is power, col is column of matrix:ylim
  for (set in 1:nSets_PC)
  {
    if (set==1)
    {
      if(col==1)
      {
        plot(powerTables_PC[[set]]$data[,1], -sign(powerTables_PC[[set]]$data[,3])*powerTables_PC[[set]]$data[,2],
             xlab="Soft Threshold (power)",ylab=colNames[col],type="n", ylim = c(0,1),
             main = colNames[col])
        abline(h=0.8, col = "red")
        
      }else{
        plot(powerTables_PC[[set]]$data[,1], -sign(powerTables_PC[[set]]$data[,3])*powerTables_PC[[set]]$data[,2],
             xlab="Soft Threshold (power)",ylab=colNames[col],type="n", ylim = range(powerTables_PC[[set]]$data[,ploPCols[col]]),
             main = colNames[col])
      }
      addGrid()  
    }
    
    if (col==1)
    {
      text(powerTables_PC[[set]]$data[,1], -sign(powerTables_PC[[set]]$data[,3])*powerTables_PC[[set]]$data[,2],
           labels=powers,cex=cex1,col=colors[set]);
    } else{
      text(powerTables_PC[[set]]$data[,1], powerTables_PC[[set]]$data[,ploPCols[col]],
           labels=powers,cex=cex1,col=colors[set],ylim = range(ylim[,col]))
    }
  }
}
par(mar=c(0,0,0,0))
plot.new()
legend("center", legend = setLabels_PC, col = colors, pch = 20,cex = 1.5,ncol = nSets_PC)

dev.off()

### 3.4.2 Network construction and consensus module detection ----
setLabels_PC
softPower_PC = c(18,18,18,16)
net_PC = blockwiseConsensusModules(multiExpr_PC,
                                   power = softPower_PC,
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
consMEs_PC = net_PC$multiMEs
moduleColors_PC = net_PC$colors
consTree_PC = net_PC$dendrograms[[1]]

pdf("./2.output/2.ndd_wgcna/figure/PCx/05_Dendrogram of consensus module eigengenes obtained by clustering the eigengenes on their consensus correlation.pdf",width=7, height=4)
plotDendroAndColors(consTree_PC, moduleColors_PC,
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Consensus gene dendrogram and module colors")
dev.off()
save(consMEs_PC, moduleColors_PC, consTree_PC, file = "./2.output/2.ndd_wgcna/data/03_Consensus-NetworkConstruction-PCx.RData")

## 3.5 CB ----
powers = c(seq(4,10,by=1), seq(12,30, by=2)) 
lnames=load("./2.output/2.ndd_wgcna/data/02_Consensus_dataInput_CB-re.RData")
### 3.5.1 Choosing the soft-thresholding power： analysis of network topology ----
# delete asd casue of high connectivity (temp)
nSets_CB = checkSets(multiExpr_CB)$nSets

# Choose a set of soft-thresholding powers

# Initialize a list to hold the results of scale-free analysis
powerTables_CB = vector(mode = "list", length = nSets_CB)
# Call the network topology analysis function for each set in turn
for (set in 1:nSets_CB)
{
  powerTables_CB[[set]] = list(data = pickSoftThreshold(multiExpr_CB[[set]]$data, 
                                                        powerVector=powers,corFnc ="bicor",networkType = "signed",
                                                        verbose = 2)[[2]])  
}

# Plot the results:
colors = RColorBrewer::brewer.pal(nSets_CB,"Set2")

# Will plot these columns of the returned scale free analysis tables
ploCBols = c(2,5,6,7)
colNames = c("Scale Free Topology Model Fit", "Mean connectivity", "Median connectivity",
             "Max connectivity")

# Get the minima and maxima of the plotted points
ylim = matrix(NA, nrow = nSets_CB, ncol = 4)
for (set in 1:nSets_CB)
{
  for (col in 1:length(ploCBols))
  {
    ylim[set, col] = min(ylim[set, col], powerTables_CB[[set]]$data[, ploCBols[col]], na.rm = TRUE)
  }
}

# Plot the quantities in the chosen columns vs. the soft thresholding power
sizeGrWindow(8, 6)

pdf("./2.output/2.ndd_wgcna/figure/CB/04_ReSummary network indices (y-axes) as functions of the soft thresholding power (x-axes).pdf")
#par(mCBol = c(2,2))
layout(matrix(c(1,2,3,4,5,5),ncol = 2,byrow = T),heights = c(3,3,1))
#par(mar = c(4.2, 4.2 , 2.2, 0.5))
cex1 = 1

for (col in 1:length(ploCBols))
{
  # set is power, col is column of matrix:ylim
  for (set in 1:nSets_CB)
  {
    if (set==1)
    {
      if(col==1)
      {
        plot(powerTables_CB[[set]]$data[,1], -sign(powerTables_CB[[set]]$data[,3])*powerTables_CB[[set]]$data[,2],
             xlab="Soft Threshold (power)",ylab=colNames[col],type="n", ylim = c(0,1),
             main = colNames[col])
        abline(h=0.8, col = "red")
        
      }else{
        plot(powerTables_CB[[set]]$data[,1], -sign(powerTables_CB[[set]]$data[,3])*powerTables_CB[[set]]$data[,2],
             xlab="Soft Threshold (power)",ylab=colNames[col],type="n", ylim = range(powerTables_CB[[set]]$data[,ploCBols[col]]),
             main = colNames[col])
      }
      addGrid()  
    }
    
    if (col==1)
    {
      text(powerTables_CB[[set]]$data[,1], -sign(powerTables_CB[[set]]$data[,3])*powerTables_CB[[set]]$data[,2],
           labels=powers,cex=cex1,col=colors[set]);
    } else{
      text(powerTables_CB[[set]]$data[,1], powerTables_CB[[set]]$data[,ploCBols[col]],
           labels=powers,cex=cex1,col=colors[set],ylim = range(ylim[,col]))
    }
  }
}
par(mar=c(0,0,0,0))
plot.new()
legend("center", legend = setLabels_CB, col = colors, pch = 20,cex = 1.5,ncol = nSets_CB)
dev.off()

### 3.5.2 Network construction and consensus module detection ----
setLabels_CB
softPower_CB = 14
net_CB = blockwiseConsensusModules(multiExpr_CB,
                                   power = softPower_CB,
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
consMEs_CB = net_CB$multiMEs
moduleColors_CB = net_CB$colors
consTree_CB = net_CB$dendrograms[[1]]

pdf("./2.output/2.ndd_wgcna/figure/CB/05_Dendrogram of consensus module eigengenes obtained by clustering the eigengenes on their consensus correlation.pdf",width=7, height=4)
plotDendroAndColors(consTree_CB, moduleColors_CB,
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Consensus gene dendrogram and module colors")
dev.off()
save(consMEs_CB, moduleColors_CB, consTree_CB, file = "./2.output/2.ndd_wgcna/data/03_Consensus-NetworkConstruction-CB.RData")

## 3.6 HIP ----
powers = c(seq(4,10,by=1), seq(12,30, by=2)) 
lnames=load("./2.output/2.ndd_wgcna/data/02_Consensus_dataInput_HIP.RData")
### 3.6.1 Choosing the soft-thresholding power： analysis of network topology ----
# delete asd casue of high connectivity (temp)
nSets_HIP = checkSets(multiExpr_HIP)$nSets

# Initialize a list to hold the results of scale-free analysis
powerTables_HIP = vector(mode = "list", length = nSets_HIP)
# Call the network topology analysis function for each set in turn
for (set in 1:nSets_HIP)
{
  powerTables_HIP[[set]] = list(data = pickSoftThreshold(multiExpr_HIP[[set]]$data, 
                                                        powerVector=powers,corFnc ="bicor",networkType = "signed",
                                                        verbose = 2)[[2]])  
}

# Plot the results:
colors = RColorBrewer::brewer.pal(nSets_HIP,"Set2")

# Will plot these columns of the returned scale free analysis tables
ploHIPols = c(2,5,6,7)

colNames = c("Scale Free Topology Model Fit", "Mean connectivity", "Median connectivity", "Max connectivity")

# Get the minima and maxima of the plotted points
ylim = matrix(NA, nrow = nSets_HIP, ncol = 4)
for (set in 1:nSets_HIP)
{
  for (col in 1:length(ploHIPols))
  {
    ylim[set, col] = min(ylim[set, col], powerTables_HIP[[set]]$data[, ploHIPols[col]], na.rm = TRUE)
  }
}

# Plot the quantities in the chosen columns vs. the soft thresholding power
sizeGrWindow(8, 6)

pdf("./2.output/2.ndd_wgcna/figure/HIP/04_Summary network indices (y-axes) as functions of the soft thresholding power (x-axes).pdf")
layout(matrix(c(1,2,3,4,5,5),ncol = 2,byrow = T),heights = c(3,3,1))
cex1 = 1

for (col in 1:length(ploHIPols))
{
  # set is power, col is column of matrix:ylim
  for (set in 1:nSets_HIP)
  {
    if (set==1)
    {
      if(col==1)
      {
        plot(powerTables_HIP[[set]]$data[,1], -sign(powerTables_HIP[[set]]$data[,3])*powerTables_HIP[[set]]$data[,2],
             xlab="Soft Threshold (power)",ylab=colNames[col],type="n", ylim = c(0,1),
             main = colNames[col])
        abline(h=0.8, col = "red")
        
      }else{
        plot(powerTables_HIP[[set]]$data[,1], -sign(powerTables_HIP[[set]]$data[,3])*powerTables_HIP[[set]]$data[,2],
             xlab="Soft Threshold (power)",ylab=colNames[col],type="n", ylim = range(powerTables_HIP[[set]]$data[,ploHIPols[col]]),
             main = colNames[col])
      }
      addGrid()  
    }
    
    if (col==1)
    {
      text(powerTables_HIP[[set]]$data[,1], -sign(powerTables_HIP[[set]]$data[,3])*powerTables_HIP[[set]]$data[,2],
           labels=powers,cex=cex1,col=colors[set]);
    } else{
      text(powerTables_HIP[[set]]$data[,1], powerTables_HIP[[set]]$data[,ploHIPols[col]],
           labels=powers,cex=cex1,col=colors[set],ylim = range(ylim[,col]))
    }
  }
}
par(mar=c(0,0,0,0))
plot.new()
legend("center", legend = setLabels_HIP, col = colors, pch = 20,cex = 1.5,ncol = nSets_HIP)
dev.off()

### 3.6.2 Network construction and consensus module detection ----
setLabels_HIP
softPower_HIP = c(4,14,10,7)
net_HIP = blockwiseConsensusModules(multiExpr_HIP,
                                   power = softPower_HIP,
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
consMEs_HIP = net_HIP$multiMEs
moduleColors_HIP = net_HIP$colors
consTree_HIP = net_HIP$dendrograms[[1]]

pdf("./2.output/2.ndd_wgcna/figure/HIP/05_Dendrogram of consensus module eigengenes obtained by clustering the eigengenes on their consensus correlation.pdf",width=7, height=4)
plotDendroAndColors(consTree_HIP, moduleColors_HIP,
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Consensus gene dendrogram and module colors")
dev.off()
save(consMEs_HIP, moduleColors_HIP, consTree_HIP, file = "./2.output/2.ndd_wgcna/data/03_Consensus-NetworkConstruction-HIP.RData")

## 3.7 STR ----
powers = c(seq(4,10,by=1), seq(12,30, by=2))
lnames=load("./2.output/4.ndd_wgcna/data/02_Consensus_dataInput_STR.RData")
### 3.7.1 Choosing the soft-thresholding power： analysis of network topology ----
# delete asd casue of high connectivity (temp)
nSets_STR = checkSets(multiExpr_STR)$nSets

# Initialize a list to hold the results of scale-free analysis
powerTables_STR = vector(mode = "list", length = nSets_STR)
# Call the network topology analysis function for each set in turn
for (set in 1:nSets_STR)
{
  powerTables_STR[[set]] = list(data = pickSoftThreshold(multiExpr_STR[[set]]$data, 
                                                         powerVector=powers,corFnc ="bicor",networkType = "signed",
                                                         verbose = 2)[[2]])  
}

# Plot the results:
colors = RColorBrewer::brewer.pal(nSets_STR,"Set2")

# Will plot these columns of the returned scale free analysis tables
plotcols = c(2,5,6,7)

colNames = c("Scale Free Topology Model Fit", "Mean connectivity", "Median connectivity","Max connectivity")

# Get the minima and maxima of the plotted points
ylim = matrix(NA, nrow = nSets_STR, ncol = 4)
for (set in 1:nSets_STR)
{
  for (col in 1:length(plotcols))
  {
    ylim[set, col] = min(ylim[set, col], powerTables_STR[[set]]$data[, plotcols[col]], na.rm = TRUE)
  }
}

# Plot the quantities in the chosen columns vs. the soft thresholding power
sizeGrWindow(8, 6)

pdf("./2.output/4.ndd_wgcna/figure/STR/04_ReSummary network indices (y-axes) as functions of the soft thresholding power (x-axes).pdf")
layout(matrix(c(1,2,3,4,5,5),ncol = 2,byrow = T),heights = c(3,3,1))
cex1 = 1

for (col in 1:length(plotcols))
{
  # set is power, col is column of matrix:ylim
  for (set in 1:nSets_STR)
  {
    if (set==1)
    {
      if(col==1)
      {
        plot(powerTables_STR[[set]]$data[,1], -sign(powerTables_STR[[set]]$data[,3])*powerTables_STR[[set]]$data[,2],
             xlab="Soft Threshold (power)",ylab=colNames[col],type="n", ylim = c(0,1),
             main = colNames[col])
        abline(h=0.8, col = "red")
        
      }else{
        plot(powerTables_STR[[set]]$data[,1], -sign(powerTables_STR[[set]]$data[,3])*powerTables_STR[[set]]$data[,2],
             xlab="Soft Threshold (power)",ylab=colNames[col],type="n", ylim = range(powerTables_STR[[set]]$data[,plotcols[col]]),
             main = colNames[col])
      }
      addGrid()  
    }
    
    if (col==1)
    {
      text(powerTables_STR[[set]]$data[,1], -sign(powerTables_STR[[set]]$data[,3])*powerTables_STR[[set]]$data[,2],
           labels=powers,cex=cex1,col=colors[set]);
    } else{
      text(powerTables_STR[[set]]$data[,1], powerTables_STR[[set]]$data[,plotcols[col]],
           labels=powers,cex=cex1,col=colors[set],ylim = range(ylim[,col]))
    }
  }
}
par(mar=c(0,0,0,0))
plot.new()
legend("center", legend = setLabels_STR, col = colors, pch = 20,cex = 1.5,ncol = nSets_STR)

dev.off()

### 3.7.2 Network construction and consensus module detection ----
setLabels_STR
softPower_STR = c(10,17,11)
net_STR = blockwiseConsensusModules(multiExpr_STR,
                                    power = softPower_STR,
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
consMEs_STR = net_STR$multiMEs
moduleColors_STR = net_STR$colors
consTree_STR = net_STR$dendrograms[[1]]

pdf("./2.output/4.ndd_wgcna/figure/STR/05_ReDendrogram of consensus module eigengenes obtained by clustering the eigengenes on their consensus correlation.pdf",width=7, height=4)
plotDendroAndColors(consTree_STR, moduleColors_STR,
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Consensus gene dendrogram and module colors")
dev.off()
save(consMEs_STR, moduleColors_STR, consTree_STR, file = "./2.output/4.ndd_wgcna/data/03_Consensus-NetworkConstruction-STR.RData")

## 3.8 AMY ----
powers = c(c(1:10), seq(from = 12, to=20, by=2))
lnames=load("./2.output/2.ndd_wgcna/data/02_NonConsensus_dataInput_AMY_mdd.RData")
### 3.8.1 Choosing the soft-thresholding power： analysis of network topology ----
sft = pickSoftThreshold(mdd_AMY, powerVector = powers, verbose = 5)

# Plot the quantities in the chosen columns vs. the soft thresholding power
pdf("./2.output/2.ndd_wgcna/figure/AMY/04_Summary network indices (y-axes) as functions of the soft thresholding power (x-axes).pdf",width = 6.5,height = 5.5)
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
net_AMY = blockwiseModules(mdd_AMY, 
                          power = 5,
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

pdf("./2.output/2.ndd_wgcna/figure/AMY/05_Dendrogram of genes with dissimilarity based on topological overlap together with assigned module colors.pdf",width=7, height=4)
plotDendroAndColors(Tree_AMY, moduleColors_AMY,
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Consensus gene dendrogram and module colors")
dev.off()
save(MEs_AMY, moduleColors_AMY, Tree_AMY, file = "./2.output/2.ndd_wgcna/data/03_NonConsensus-NetworkConstruction-AMY.RData")

######################################################
#                     4. Modules Trait Correlation
######################################################

## 4.1 FCx ----
load(file = "2.output/2.ndd_wgcna/data/02_Consensus_dataInput_FC-re.RData")
load(file = "2.output/2.ndd_wgcna/data/03_Consensus-NetworkConstruction-FCx.RData")

exprSize_FC = checkSets(multiExpr_FC)
nSets_FC = exprSize_FC$nSets

### 4.1.1 Calculate the correlations -----
# Set up variables to contain the module-trait correlations
moduleTraitCor_FC = list()
moduleTraitPvalue_FC = list()
moduleTraitPadj_FC = list()

# Calculate the correlations
tmp<-cor(consMEs_FC[[1]]$data, Traits_FC[[1]]$data, use = "p")

for (set in 1:nSets_FC)
{
  moduleTraitCor_FC[[set]] = bicorAndPvalue(consMEs_FC[[set]]$data, Traits_FC[[set]]$data, use = "p")$bicor
  moduleTraitPvalue_FC[[set]] = bicorAndPvalue(consMEs_FC[[set]]$data, Traits_FC[[set]]$data, use = "p")$p
  #bicorAndPvalue(moduleTraitCor_FC[[set]], exprSize_FC$nSamples[set])$
  moduleTraitPadj_FC[[set]] = p.adjust(moduleTraitPvalue_FC[[set]], method="BH")
}

# Convert numerical lables to colors for labeling of modules in the plot
MEColorNames_FC = names(consMEs_FC[[1]]$data)

### 4.1.2 calculate the consensus correlation ----
# Initialize matrices to hold the consensus correlation and p-value
consensusCor_FC = matrix(NA, nrow(moduleTraitCor_FC[[1]]), ncol(moduleTraitCor_FC[[1]]))
consensusPvalue_FC = matrix(NA, nrow(moduleTraitCor_FC[[1]]), ncol(moduleTraitCor_FC[[1]]))
consensusPadj_FC = matrix(NA, nrow(moduleTraitCor_FC[[1]]), ncol(moduleTraitCor_FC[[1]]))

# Find consensus negative correlations
negative = apply(as.data.frame(sapply(moduleTraitCor_FC, function(x) x<0 )),1,function(y){all(y)})

if(length(which(negative==TRUE)) > 1)
{
  consensusCor_FC[negative] = apply(do.call(cbind,moduleTraitCor_FC)[negative,],1,max)
  consensusPvalue_FC[negative] = apply(do.call(cbind,moduleTraitPvalue_FC)[negative,],1,max)
  consensusPadj_FC[negative] = apply(do.call(cbind,moduleTraitPadj_FC)[negative,],1,max)
}else if(length(which(negative==TRUE)) == 1){
  consensusCor_FC[negative] = max(do.call(cbind,moduleTraitCor_FC)[negative,])
  consensusPvalue_FC[negative] = max(do.call(cbind,moduleTraitPvalue_FC)[negative,])
  consensusPadj_FC[negative] = max(do.call(cbind,moduleTraitPadj_FC)[negative,])
}

# Find consensus positive correlations
positive = apply(as.data.frame(sapply(moduleTraitCor_FC, function(x) x>0 )),1,function(y){all(y)})

if(length(which(positive==TRUE)) > 1)
{
  consensusCor_FC[positive] = apply(do.call(cbind,moduleTraitCor_FC)[positive,],1,max)
  consensusPvalue_FC[positive] = apply(do.call(cbind,moduleTraitPvalue_FC)[positive,],1,max)
  consensusPadj_FC[positive] = apply(do.call(cbind,moduleTraitPadj_FC)[positive,],1,max)
}else if(length(which(positive==TRUE)) == 1){
  consensusCor_FC[positive] = max(do.call(cbind,moduleTraitCor_FC)[positive,])
  consensusPvalue_FC[positive] = max(do.call(cbind,moduleTraitPvalue_FC)[positive,])
  consensusPadj_FC[positive] = max(do.call(cbind,moduleTraitPadj_FC)[positive,])
}

### 4.1.3 Export result of module trait correlation ----

for (set in 1:nSets_FC) {
  
  moduleTraitCorRes = cbind(cor = moduleTraitCor_FC[[set]], 
                            p = moduleTraitPvalue_FC[[set]],
                            fdr = moduleTraitPadj_FC[[set]]) %>% as.data.frame()
  colnames(moduleTraitCorRes) <- paste0(setLabels_FC[[set]],
                                        c(".cor",".p",".fdr")
                                        )
  fwrite(moduleTraitCorRes %>% rownames_to_column("module"),
         file = paste0("2.output/2.ndd_wgcna/data/04_ModuleTraitCorRes/FCx/04_moduleTraitCorRes-FCx-",setLabels_FC[set],".txt"),
         sep = "\t")  
}

moduleTraitCorRes = cbind(consensusCor_FC, consensusPvalue_FC, consensusPadj_FC) %>% as.data.frame()
colnames(moduleTraitCorRes) <- paste0("consensus.status.",c("cor","p","fdr"))
rownames(moduleTraitCorRes) <- rownames(moduleTraitCor_FC[[1]])
fwrite(moduleTraitCorRes %>% rownames_to_column("module"),
       file = paste0("2.output/2.ndd_wgcna/data/04_ModuleTraitCorRes/FCx/04_moduleTraitCorRes-FCx-consensus.txt"),
       sep = "\t")  

## 4.2 TCx ----
load(file = "2.output/2.ndd_wgcna/data/02_Consensus_dataInput_TC.RData")
load(file = "2.output/2.ndd_wgcna/data/03_Consensus-NetworkConstruction-TCx.RData")

exprSize_TC = checkSets(multiExpr_TC)
nSets_TC = exprSize_TC$nSets

### 4.2.1 Calculate the correlations -----
# Set up variables to contain the module-trait correlations
moduleTraitCor_TC = list()
moduleTraitPvalue_TC = list()
moduleTraitPadj_TC = list()
# Calculate the correlations
for (set in 1:nSets_TC)
{
  moduleTraitCor_TC[[set]] = bicorAndPvalue(consMEs_TC[[set]]$data, Traits_TC[[set]]$data, use = "p")$bicor
  moduleTraitPvalue_TC[[set]] = bicorAndPvalue(consMEs_TC[[set]]$data, Traits_TC[[set]]$data, use = "p")$p
  #corPvalueFisher(moduleTraitCor_TC[[set]], exprSize_TC$nSamples[set])
  moduleTraitPadj_TC[[set]] = p.adjust(moduleTraitPvalue_TC[[set]], method="BH")
}

# Convert numerical lables to colors for labeling of modules in the plot
MEColorNames_TC = names(consMEs_TC[[1]]$data)

### 4.2.2 calculate the consensus correlation ----
# Initialize matrices to hold the consensus correlation and p-value
consensusCor_TC = matrix(NA, nrow(moduleTraitCor_TC[[1]]), ncol(moduleTraitCor_TC[[1]]))
consensusPvalue_TC = matrix(NA, nrow(moduleTraitCor_TC[[1]]), ncol(moduleTraitCor_TC[[1]]))
consensusPadj_TC = matrix(NA, nrow(moduleTraitCor_TC[[1]]), ncol(moduleTraitCor_TC[[1]]))

# Find consensus negative correlations
negative = apply(as.data.frame(sapply(moduleTraitCor_TC, function(x) x<0 )),1,function(y){all(y)})

consensusCor_TC[negative] = apply(do.call(cbind,moduleTraitCor_TC)[negative,],1,max)
consensusPvalue_TC[negative] = apply(do.call(cbind,moduleTraitPvalue_TC)[negative,],1,max)
consensusPadj_TC[negative] = apply(do.call(cbind,moduleTraitPadj_TC)[negative,],1,max)

# Find consensus positive correlations
positive = apply(as.data.frame(sapply(moduleTraitCor_TC, function(x) x>0 )),1,function(y){all(y)})

consensusCor_TC[positive] = apply(do.call(cbind,moduleTraitCor_TC)[positive,],1,max)
consensusPvalue_TC[positive] = apply(do.call(cbind,moduleTraitPvalue_TC)[positive,],1,max)
consensusPadj_TC[positive] = apply(do.call(cbind,moduleTraitPadj_TC)[positive,],1,max)

### 4.2.3 Export result of module trait correlation ----

for (set in 1:nSets_TC) {
  
  moduleTraitCorRes = cbind(cor = moduleTraitCor_TC[[set]], 
                            p = moduleTraitPvalue_TC[[set]],
                            fdr = moduleTraitPadj_TC[[set]]) %>% as.data.frame()
  colnames(moduleTraitCorRes) <- paste0(setLabels_TC[[set]],
                                        c(".cor",".p",".fdr")
  )
  fwrite(moduleTraitCorRes %>% rownames_to_column("module"),
         file = paste0("2.output/2.ndd_wgcna/data/04_ModuleTraitCorRes/TCx/04_ModuleTraitCorRes-",setLabels_TC[set],".txt"),
         sep = "\t")  
}

moduleTraitCorRes = cbind(consensusCor_TC, consensusPvalue_TC, consensusPadj_TC) %>% as.data.frame()
colnames(moduleTraitCorRes) <- paste0("consensus.status.",c("cor","p","fdr"))
rownames(moduleTraitCorRes) <- rownames(moduleTraitCor_TC[[1]])
fwrite(moduleTraitCorRes %>% rownames_to_column("module"),
       file = paste0("2.output/2.ndd_wgcna/data/04_ModuleTraitCorRes/TCx/04_moduleTraitCorRes-TCx-consensus.txt"),
       sep = "\t")  

## 4.3 OCx ----
load(file = "2.output/2.ndd_wgcna/data/02_NonConsensus_dataInput_OC_asd.RData")
load(file = "2.output/2.ndd_wgcna/data/03_NonConsensus-NetworkConstruction-OCx.RData")

### 4.3.1 Calculate the correlations -----
# Set up variables to contain the module-trait correlations
moduleTraitCor_OC = bicorAndPvalue(MEs_OC, Traits_OC, use = "p")$bicor
moduleTraitPvalue_OC = bicorAndPvalue(MEs_OC, Traits_OC, use = "p")$p
moduleTraitPadj_OC = p.adjust(moduleTraitPvalue_OC, method="BH")

# Convert numerical lables to colors for labeling of modules in the plot
MEColorNames_OC = names(MEs_OC)

### 4.3.2 Export result of module trait correlation ----
moduleTraitCorRes = cbind(cor = moduleTraitCor_OC, 
                          p = moduleTraitPvalue_OC,
                          fdr = moduleTraitPadj_OC) %>% as.data.frame()
colnames(moduleTraitCorRes) <- paste0("ASD",
                                      c(".cor",".p",".fdr")
)
fwrite(moduleTraitCorRes %>% rownames_to_column("module"),
       file = paste0("2.output/2.ndd_wgcna/data/04_ModuleTraitCorRes/OCx/04_ModuleTraitCorRes-ASD.txt"),
       sep = "\t")  

## 4.4 PCx ----
load(file = "2.output/2.ndd_wgcna/data/02_Consensus_dataInput_PC.RData")
load(file = "2.output/2.ndd_wgcna/data/03_Consensus-NetworkConstruction-PCx.RData")

exprSize_PC = checkSets(multiExpr_PC)
nSets_PC = exprSize_PC$nSets

### 4.4.1 Calculate the correlations -----
# Set up variables to contain the module-trait correlations
moduleTraiPCor_PC = list()
moduleTraitPvalue_PC = list()
moduleTraitPadj_PC = list()
# Calculate the correlations
for (set in 1:nSets_PC)
{
  moduleTraiPCor_PC[[set]] = bicorAndPvalue(consMEs_PC[[set]]$data, Traits_PC[[set]]$data, use = "p")$bicor
  moduleTraitPvalue_PC[[set]] = bicorAndPvalue(consMEs_PC[[set]]$data, Traits_PC[[set]]$data, use = "p")$p
  moduleTraitPadj_PC[[set]] = p.adjust(moduleTraitPvalue_PC[[set]], method="BH")
}

# Convert numerical lables to colors for labeling of modules in the plot
MEColorNames_PC = names(consMEs_PC[[1]]$data)

### 4.4.2 calculate the consensus correlation ----
# Find consensus negative correlations
consensusCor_PC = matrix(NA, nrow(moduleTraiPCor_PC[[1]]), ncol(moduleTraiPCor_PC[[1]]))
consensusPvalue_PC = matrix(NA, nrow(moduleTraiPCor_PC[[1]]), ncol(moduleTraiPCor_PC[[1]]))
consensusPadj_PC = matrix(NA, nrow(moduleTraiPCor_PC[[1]]), ncol(moduleTraiPCor_PC[[1]]))

negative = apply(as.data.frame(sapply(moduleTraiPCor_PC, function(x) x<0 )),1,function(y){all(y)})

if(length(which(negative==TRUE)) > 1)
{
  consensusCor_PC[negative] = apply(do.call(cbind,moduleTraiPCor_PC)[negative,],1,max)
  consensusPvalue_PC[negative] = apply(do.call(cbind,moduleTraitPvalue_PC)[negative,],1,max)
  consensusPadj_PC[negative] = apply(do.call(cbind,moduleTraitPadj_PC)[negative,],1,max)
}else if(length(which(negative==TRUE)) == 1){
  consensusCor_PC[negative] = max(do.call(cbind,moduleTraiPCor_PC)[negative,])
  consensusPvalue_PC[negative] = max(do.call(cbind,moduleTraitPvalue_PC)[negative,])
  consensusPadj_PC[negative] = max(do.call(cbind,moduleTraitPadj_PC)[negative,])
}

# Find consensus positive correlations
positive = apply(as.data.frame(sapply(moduleTraiPCor_PC, function(x) x>0 )),1,function(y){all(y)})

if(length(which(positive==TRUE)) > 1)
{
  consensusCor_PC[positive] = apply(do.call(cbind,moduleTraiPCor_PC)[positive,],1,max)
  consensusPvalue_PC[positive] = apply(do.call(cbind,moduleTraitPvalue_PC)[positive,],1,max)
  consensusPadj_PC[positive] = apply(do.call(cbind,moduleTraitPadj_PC)[positive,],1,max)
}else if(length(which(positive==TRUE)) == 1){
  consensusCor_PC[positive] = max(do.call(cbind,moduleTraiPCor_PC)[positive,])
  consensusPvalue_PC[positive] = max(do.call(cbind,moduleTraitPvalue_PC)[positive,])
  consensusPadj_PC[positive] = max(do.call(cbind,moduleTraitPadj_PC)[positive,])
}

### 4.4.3 Export result of module trait correlation ----

for (set in 1:nSets_PC) {
  
  moduleTraiPCorRes = cbind(cor = moduleTraiPCor_PC[[set]], 
                            p = moduleTraitPvalue_PC[[set]],
                            fdr = moduleTraitPadj_PC[[set]]) %>% as.data.frame()
  colnames(moduleTraiPCorRes) <- paste0(setLabels_PC[[set]],
                                        c(".cor",".p",".fdr")
  )
  fwrite(moduleTraiPCorRes %>% rownames_to_column("module"),
         file = paste0("2.output/2.ndd_wgcna/data/04_ModuleTraitCorRes/PCx/04_ModuleTraitCorRes-",setLabels_PC[set],".txt"),
         sep = "\t")  
}

moduleTraiPCorRes = cbind(consensusCor_PC, consensusPvalue_PC, consensusPadj_PC) %>% as.data.frame()
colnames(moduleTraiPCorRes) <- paste0("consensus.status.",c("cor","p","fdr"))
rownames(moduleTraiPCorRes) <- rownames(moduleTraiPCor_PC[[1]])
fwrite(moduleTraiPCorRes %>% rownames_to_column("module"),
       file = paste0("2.output/2.ndd_wgcna/data/04_ModuleTraitCorRes/PCx/04_moduleTraitCorRes-PCx-consensus.txt"),
       sep = "\t")  

## 4.5 CB ----
load(file = "2.output/2.ndd_wgcna/data/02_Consensus_dataInput_CB-re.RData")
load(file = "2.output/2.ndd_wgcna/data/03_Consensus-NetworkConstruction-CB.RData")

exprSize_CB = checkSets(multiExpr_CB)
nSets_CB = exprSize_CB$nSets

### 4.5.1 Calculate the correlations -----
# Set up variables to contain the module-trait correlations
moduleTraitCor_CB = list()
moduleTraitPvalue_CB = list()
moduleTraitPadj_CB = list()
# Calculate the correlations
for (set in 1:nSets_CB)
{
  moduleTraitCor_CB[[set]] = bicorAndPvalue(consMEs_CB[[set]]$data, Traits_CB[[set]]$data, use = "p")$bicor
  moduleTraitPvalue_CB[[set]] = bicorAndPvalue(consMEs_CB[[set]]$data, Traits_CB[[set]]$data, use = "p")$p
  #corPvalueFisher(moduleTraitCor_CB[[set]], exprSize_CB$nSamples[set])
  moduleTraitPadj_CB[[set]] = p.adjust(moduleTraitPvalue_CB[[set]], method="BH")
}

# Convert numerical lables to colors for labeling of modules in the plot
MEColorNames_CB = names(consMEs_CB[[1]]$data)

### 4.5.2 calculate the consensus correlation ----
# Find consensus negative correlations
consensusCor_CB = matrix(NA, nrow(moduleTraitCor_CB[[1]]), ncol(moduleTraitCor_CB[[1]]))
consensusPvalue_CB = matrix(NA, nrow(moduleTraitCor_CB[[1]]), ncol(moduleTraitCor_CB[[1]]))
consensusPadj_CB = matrix(NA, nrow(moduleTraitCor_CB[[1]]), ncol(moduleTraitCor_CB[[1]]))

negative = apply(as.data.frame(sapply(moduleTraitCor_CB, function(x) x<0 )),1,function(y){all(y)})

if(length(which(negative==TRUE)) > 1)
{
  consensusCor_CB[negative] = apply(do.call(cbind,moduleTraitCor_CB)[negative,],1,max)
  consensusPvalue_CB[negative] = apply(do.call(cbind,moduleTraitPvalue_CB)[negative,],1,max)
  consensusPadj_CB[negative] = apply(do.call(cbind,moduleTraitPadj_CB)[negative,],1,max)
}else if(length(which(negative==TRUE)) == 1){
  consensusCor_CB[negative] = max(do.call(cbind,moduleTraitCor_CB)[negative,])
  consensusPvalue_CB[negative] = max(do.call(cbind,moduleTraitPvalue_CB)[negative,])
  consensusPadj_CB[negative] = max(do.call(cbind,moduleTraitPadj_CB)[negative,])
}

# Find consensus positive correlations
positive = apply(as.data.frame(sapply(moduleTraitCor_CB, function(x) x>0 )),1,function(y){all(y)})

if(length(which(positive==TRUE)) > 1)
{
  consensusCor_CB[positive] = apply(do.call(cbind,moduleTraitCor_CB)[positive,],1,max)
  consensusPvalue_CB[positive] = apply(do.call(cbind,moduleTraitPvalue_CB)[positive,],1,max)
  consensusPadj_CB[positive] = apply(do.call(cbind,moduleTraitPadj_CB)[positive,],1,max)
}else if(length(which(positive==TRUE)) == 1){
  consensusCor_CB[positive] = max(do.call(cbind,moduleTraitCor_CB)[positive,])
  consensusPvalue_CB[positive] = max(do.call(cbind,moduleTraitPvalue_CB)[positive,])
  consensusPadj_CB[positive] = max(do.call(cbind,moduleTraitPadj_CB)[positive,])
}

### 4.5.3 Export result of module trait correlation ----
for (set in 1:nSets_CB) {
  
  moduleTraiCBorRes = cbind(cor = moduleTraitCor_CB[[set]], 
                            p = moduleTraitPvalue_CB[[set]],
                            fdr = moduleTraitPadj_CB[[set]]) %>% as.data.frame()
  colnames(moduleTraiCBorRes) <- paste0(setLabels_CB[[set]],
                                        c(".cor",".p",".fdr")
  )
  fwrite(moduleTraiCBorRes %>% rownames_to_column("module"),
         file = paste0("2.output/2.ndd_wgcna/data/04_ModuleTraitCorRes/CB/04_ModuleTraitCorRes-",setLabels_CB[set],".txt"),
         sep = "\t")  
}

moduleTraiCBorRes = cbind(consensusCor_CB, consensusPvalue_CB, consensusPadj_CB) %>% as.data.frame()
colnames(moduleTraiCBorRes) <- paste0("consensus.status.",c("cor","p","fdr"))
rownames(moduleTraiCBorRes) <- rownames(moduleTraitCor_CB[[1]])
fwrite(moduleTraiCBorRes %>% rownames_to_column("module"),
       file = paste0("2.output/2.ndd_wgcna/data/04_ModuleTraitCorRes/CB/04_moduleTraitCorRes-CB-consensus.txt"),
       sep = "\t")  

## 4.6 HIP ----
load(file = "2.output/2.ndd_wgcna/data/02_Consensus_dataInput_HIP.RData")
load(file = "2.output/2.ndd_wgcna/data/03_Consensus-NetworkConstruction-HIP.RData")

exprSize_HIP = checkSets(multiExpr_HIP)
nSets_HIP = exprSize_HIP$nSets

### 4.6.1 Calculate the correlations -----
# Set up variables to contain the module-trait correlations
moduleTraitCor_HIP = list()
moduleTraitPvalue_HIP = list()
moduleTraitPadj_HIP = list()
# Calculate the correlations
for (set in 1:nSets_HIP)
{
  moduleTraitCor_HIP[[set]] = bicorAndPvalue(consMEs_HIP[[set]]$data, Traits_HIP[[set]]$data, use = "p")$bicor
  moduleTraitPvalue_HIP[[set]] = bicorAndPvalue(consMEs_HIP[[set]]$data, Traits_HIP[[set]]$data, use = "p")$p
  #corPvalueFisher(moduleTraitCor_HIP[[set]], exprSize_HIP$nSamples[set])
  moduleTraitPadj_HIP[[set]] = p.adjust(moduleTraitPvalue_HIP[[set]], method="BH")
}

# Convert numerical lables to colors for labeling of modules in the plot
MEColorNames_HIP = names(consMEs_HIP[[1]]$data)

### 4.6.2 calculate the consensus correlation ----
# Find consensus negative correlations
consensusCor_HIP = matrix(NA, nrow(moduleTraitCor_HIP[[1]]), ncol(moduleTraitCor_HIP[[1]]))
consensusPvalue_HIP = matrix(NA, nrow(moduleTraitCor_HIP[[1]]), ncol(moduleTraitCor_HIP[[1]]))
consensusPadj_HIP = matrix(NA, nrow(moduleTraitCor_HIP[[1]]), ncol(moduleTraitCor_HIP[[1]]))

negative = apply(as.data.frame(sapply(moduleTraitCor_HIP, function(x) x<0 )),1,function(y){all(y)})

if(length(which(negative==TRUE)) > 1)
{
  consensusCor_HIP[negative] = apply(do.call(cbind,moduleTraitCor_HIP)[negative,],1,max)
  consensusPvalue_HIP[negative] = apply(do.call(cbind,moduleTraitPvalue_HIP)[negative,],1,max)
  consensusPadj_HIP[negative] = apply(do.call(cbind,moduleTraitPadj_HIP)[negative,],1,max)
}else if(length(which(negative==TRUE)) == 1){
  consensusCor_HIP[negative] = max(do.call(cbind,moduleTraitCor_HIP)[negative,])
  consensusPvalue_HIP[negative] = max(do.call(cbind,moduleTraitPvalue_HIP)[negative,])
  consensusPadj_HIP[negative] = max(do.call(cbind,moduleTraitPadj_HIP)[negative,])
}

# Find consensus positive correlations
positive = apply(as.data.frame(sapply(moduleTraitCor_HIP, function(x) x>0 )),1,function(y){all(y)})

if(length(which(positive==TRUE)) > 1)
{
  consensusCor_HIP[positive] = apply(do.call(cbind,moduleTraitCor_HIP)[positive,],1,max)
  consensusPvalue_HIP[positive] = apply(do.call(cbind,moduleTraitPvalue_HIP)[positive,],1,max)
  consensusPadj_HIP[positive] = apply(do.call(cbind,moduleTraitPadj_HIP)[positive,],1,max)
}else if(length(which(positive==TRUE)) == 1){
  consensusCor_HIP[positive] = max(do.call(cbind,moduleTraitCor_HIP)[positive,])
  consensusPvalue_HIP[positive] = max(do.call(cbind,moduleTraitPvalue_HIP)[positive,])
  consensusPadj_HIP[positive] = max(do.call(cbind,moduleTraitPadj_HIP)[positive,])
}

### 4.6.3 Export result of module trait correlation ----
for (set in 1:nSets_HIP) {
  
  moduleTraiHIPorRes = cbind(cor = moduleTraitCor_HIP[[set]], 
                            p = moduleTraitPvalue_HIP[[set]],
                            fdr = moduleTraitPadj_HIP[[set]]) %>% as.data.frame()
  colnames(moduleTraiHIPorRes) <- paste0(setLabels_HIP[[set]],
                                        c(".cor",".p",".fdr")
  )
  fwrite(moduleTraiHIPorRes %>% rownames_to_column("module"),
         file = paste0("2.output/2.ndd_wgcna/data/04_ModuleTraitCorRes/HIP/04_ModuleTraitCorRes-",setLabels_HIP[set],".txt"),
         sep = "\t")  
}

moduleTraiHIPorRes = cbind(consensusCor_HIP, consensusPvalue_HIP, consensusPadj_HIP) %>% as.data.frame()
colnames(moduleTraiHIPorRes) <- paste0("consensus.status.",c("cor","p","fdr"))
rownames(moduleTraiHIPorRes) <- rownames(moduleTraitCor_HIP[[1]])
fwrite(moduleTraiHIPorRes %>% rownames_to_column("module"),
       file = paste0("2.output/2.ndd_wgcna/data/04_ModuleTraitCorRes/HIP/04_moduleTraitCorRes-HIP-consensus.txt"),
       sep = "\t")  

## 4.7 STR ----
load(file = "2.output/4.ndd_wgcna/data/02_Consensus_dataInput_STR.RData")
load(file = "2.output/4.ndd_wgcna/data/03_Consensus-NetworkConstruction-STR.RData")

exprSize_STR = checkSets(multiExpr_STR)
nSets_STR = exprSize_STR$nSets

### 4.7.1 Calculate the correlations -----
# Set up variables to contain the module-trait correlations
moduleTraitCor_STR = list()
moduleTraitPvalue_STR = list()
moduleTraitPadj_STR = list()
# Calculate the correlations
for (set in 1:nSets_STR)
{
  moduleTraitCor_STR[[set]] = bicorAndPvalue(consMEs_STR[[set]]$data, Traits_STR[[set]]$data, use = "p")$bicor
  moduleTraitPvalue_STR[[set]] = bicorAndPvalue(consMEs_STR[[set]]$data, Traits_STR[[set]]$data, use = "p")$p
  moduleTraitPadj_STR[[set]] = p.adjust(moduleTraitPvalue_STR[[set]], method="BH")
}

# Convert numerical lables to colors for labeling of modules in the plot
MEColorNames_STR = names(consMEs_STR[[1]]$data)

### 4.7.2 calculate the consensus correlation ----
# Find consensus negative correlations
consensusCor_STR = matrix(NA, nrow(moduleTraitCor_STR[[1]]), ncol(moduleTraitCor_STR[[1]]))
consensusPvalue_STR = matrix(NA, nrow(moduleTraitCor_STR[[1]]), ncol(moduleTraitCor_STR[[1]]))
consensusPadj_STR = matrix(NA, nrow(moduleTraitCor_STR[[1]]), ncol(moduleTraitCor_STR[[1]]))

negative = apply(as.data.frame(sapply(moduleTraitCor_STR, function(x) x<0 )),1,function(y){all(y)})

if(length(which(negative==TRUE)) > 1)
{
  consensusCor_STR[negative] = apply(do.call(cbind,moduleTraitCor_STR)[negative,],1,max)
  consensusPvalue_STR[negative] = apply(do.call(cbind,moduleTraitPvalue_STR)[negative,],1,max)
  consensusPadj_STR[negative] = apply(do.call(cbind,moduleTraitPadj_STR)[negative,],1,max)
}else if(length(which(negative==TRUE)) == 1){
  consensusCor_STR[negative] = max(do.call(cbind,moduleTraitCor_STR)[negative,])
  consensusPvalue_STR[negative] = max(do.call(cbind,moduleTraitPvalue_STR)[negative,])
  consensusPadj_STR[negative] = max(do.call(cbind,moduleTraitPadj_STR)[negative,])
}

# Find consensus positive correlations
positive = apply(as.data.frame(sapply(moduleTraitCor_STR, function(x) x>0 )),1,function(y){all(y)})

if(length(which(positive==TRUE)) > 1)
{
  consensusCor_STR[positive] = apply(do.call(cbind,moduleTraitCor_STR)[positive,],1,max)
  consensusPvalue_STR[positive] = apply(do.call(cbind,moduleTraitPvalue_STR)[positive,],1,max)
  consensusPadj_STR[positive] = apply(do.call(cbind,moduleTraitPadj_STR)[positive,],1,max)
}else if(length(which(positive==TRUE)) == 1){
  consensusCor_STR[positive] = max(do.call(cbind,moduleTraitCor_STR)[positive,])
  consensusPvalue_STR[positive] = max(do.call(cbind,moduleTraitPvalue_STR)[positive,])
  consensusPadj_STR[positive] = max(do.call(cbind,moduleTraitPadj_STR)[positive,])
}

### 4.7.3 Export result of module trait correlation ----
for (set in 1:nSets_STR) {
  
  moduleTraiSTRorRes = cbind(cor = moduleTraitCor_STR[[set]], 
                             p = moduleTraitPvalue_STR[[set]],
                             fdr = moduleTraitPadj_STR[[set]]) %>% as.data.frame()
  colnames(moduleTraiSTRorRes) <- paste0(setLabels_STR[[set]],
                                         c(".cor",".p",".fdr")
  )
  fwrite(moduleTraiSTRorRes %>% rownames_to_column("module"),
         file = paste0("2.output/4.ndd_wgcna/data/04_ModuleTraitCorRes/STR/04_ModuleTraitCorRes-",setLabels_STR[set],".txt"),
         sep = "\t")  
}

moduleTraiSTRorRes = cbind(consensusCor_STR, consensusPvalue_STR, consensusPadj_STR) %>% as.data.frame()
colnames(moduleTraiSTRorRes) <- paste0("consensus.status.",c("cor","p","fdr"))
rownames(moduleTraiSTRorRes) <- rownames(moduleTraitCor_STR[[1]])
fwrite(moduleTraiSTRorRes %>% rownames_to_column("module"),
       file = paste0("2.output/4.ndd_wgcna/data/04_ModuleTraitCorRes/STR/04_moduleTraitCorRes-STR-consensus.txt"),
       sep = "\t")  

## 4.8 AMY ----
load(file = "2.output/2.ndd_wgcna/data/02_NonConsensus_dataInput_AMY_mdd.RData")
load(file = "2.output/2.ndd_wgcna/data/03_NonConsensus-NetworkConstruction-AMY.RData")

### 4.8.1 Calculate the correlations -----
# Set up variables to contain the module-trait correlations
moduleTraitCor_AMY = bicorAndPvalue(MEs_AMY, Traits_AMY, use = "p")$bicor
moduleTraitPvalue_AMY = bicorAndPvalue(MEs_AMY, Traits_AMY, use = "p")$p
moduleTraitPadj_AMY = p.adjust(moduleTraitPvalue_AMY, method="BH")

# Convert numerical lables to colors for labeling of modules in the plot
MEColorNames_AMY = names(MEs_AMY)

### 4.8.2 Export result of module trait correlation ----
moduleTraitCorRes = cbind(cor = moduleTraitCor_AMY, 
                          p = moduleTraitPvalue_AMY,
                          fdr = moduleTraitPadj_AMY) %>% as.data.frame()
colnames(moduleTraitCorRes) <- paste0("MDD",
                                      c(".cor",".p",".fdr")
)
fwrite(moduleTraitCorRes %>% rownames_to_column("module"),
       file = paste0("2.output/2.ndd_wgcna/data/04_ModuleTraitCorRes/AMY/04_ModuleTraitCorRes-MDD.txt"),
       sep = "\t")  

######################################################
#                     5.Exporting results of the network analysis
######################################################
# 5.1 Module Trait Correlation Results collection -----
## 5.1.1 fc ----
cor_fc_bd <- fread("2.output/4.ndd_wgcna/data/04_ModuleTraitCorRes/FCx/04_moduleTraitCorRes-FCx-BD.txt")
cor_fc_dep <- fread("2.output/4.ndd_wgcna/data/04_ModuleTraitCorRes/FCx/04_moduleTraitCorRes-FCx-DEP.txt")
cor_fc_mdd <- fread("2.output/4.ndd_wgcna/data/04_ModuleTraitCorRes/FCx/04_moduleTraitCorRes-FCx-MDD.txt")
cor_fc_scz <- fread("2.output/4.ndd_wgcna/data/04_ModuleTraitCorRes/FCx/04_moduleTraitCorRes-FCx-SCZ.txt")

cor_fc <- rbind(
  cor_fc_bd %>% dplyr::select(1,cor=2,p=3,fdr=4) %>% dplyr::mutate(disorder="BD"),
  cor_fc_dep %>% dplyr::select(1,cor=2,p=3,fdr=4) %>% dplyr::mutate(disorder="DEP"),
  cor_fc_mdd %>% dplyr::select(1,cor=2,p=3,fdr=4) %>% dplyr::mutate(disorder="MDD"),
  cor_fc_scz %>% dplyr::select(1,cor=2,p=3,fdr=4) %>% dplyr::mutate(disorder="SCZ")
) 

cor_fc_res <- rbind(
  cor_fc,
  cor_fc %>% 
    group_by(disorder,module) %>% 
    dplyr::mutate(fdr = p.adjust(p,method = "BH")) %>% 
    ungroup() %>% 
    group_by(module) %>% 
    dplyr::mutate(cor = case_when(max(cor)<=0 ~ max(cor),
                                  min(cor)>=0 ~ min(cor),
                                  TRUE ~ NA)) %>% 
    dplyr::mutate(p = ifelse(is.na(cor),NA,max(p)),
                  fdr = ifelse(is.na(cor),NA,max(fdr))) %>% 
    dplyr::select(module,cor,p,fdr) %>% 
    dplyr::mutate(disorder="consensus") %>% unique()
)

## 5.1.2 tc ----
cor_tc_asd <- fread("2.output/4.ndd_wgcna/data/04_ModuleTraitCorRes/TCx/04_ModuleTraitCorRes-ASD.txt")
cor_tc_dem <- fread("2.output/4.ndd_wgcna/data/04_ModuleTraitCorRes/TCx/04_ModuleTraitCorRes-DEM.txt")
cor_tc_scz <- fread("2.output/4.ndd_wgcna/data/04_ModuleTraitCorRes/TCx/04_ModuleTraitCorRes-SCZ.txt")

cor_tc <- rbind(
  cor_tc_asd %>% dplyr::select(1,cor=2,p=3,fdr=4) %>% dplyr::mutate(disorder="ASD"),
  cor_tc_dem %>% dplyr::select(1,cor=2,p=3,fdr=4) %>% dplyr::mutate(disorder="DEM"),
  cor_tc_scz %>% dplyr::select(1,cor=2,p=3,fdr=4) %>% dplyr::mutate(disorder="SCZ")
) 

cor_tc_res <- rbind(
  cor_tc,
  cor_tc %>% 
    group_by(disorder,module) %>% 
    dplyr::mutate(fdr = p.adjust(p,method = "BH")) %>% 
    ungroup() %>% 
    group_by(module) %>% 
    dplyr::mutate(cor = case_when(max(cor)<=0 ~ max(cor),
                                  min(cor)>=0 ~ min(cor),
                                  TRUE ~ NA)) %>% 
    dplyr::mutate(p = ifelse(is.na(cor),NA,max(p)),
                  fdr = ifelse(is.na(cor),NA,max(fdr))) %>% 
    dplyr::select(module,cor,p,fdr) %>% 
    dplyr::mutate(disorder="consensus") %>% unique()
)

## 5.1.3 oc ----
cor_oc_asd <- fread("2.output/4.ndd_wgcna/data/04_ModuleTraitCorRes/OCx/04_ModuleTraitCorRes-ASD.txt")
cor_oc_res <- cor_oc_asd %>% dplyr::select(1,cor=2,p=3,fdr=4) %>% dplyr::mutate(disorder="ASD")

## 5.1.4 pc ----
cor_pc_bd <- fread("2.output/4.ndd_wgcna/data/04_ModuleTraitCorRes/PCx/04_ModuleTraitCorRes-BD.txt")
cor_pc_dem <- fread("2.output/4.ndd_wgcna/data/04_ModuleTraitCorRes/PCx/04_ModuleTraitCorRes-DEM.txt")
cor_pc_dep <- fread("2.output/4.ndd_wgcna/data/04_ModuleTraitCorRes/PCx/04_ModuleTraitCorRes-DEP.txt")
cor_pc_scz <- fread("2.output/4.ndd_wgcna/data/04_ModuleTraitCorRes/PCx/04_ModuleTraitCorRes-SCZ.txt")

cor_pc <- rbind(
  cor_pc_bd %>% dplyr::select(1,cor=2,p=3,fdr=4) %>% dplyr::mutate(disorder="BD"),
  cor_pc_dem %>% dplyr::select(1,cor=2,p=3,fdr=4) %>% dplyr::mutate(disorder="DEM"),
  cor_pc_dep %>% dplyr::select(1,cor=2,p=3,fdr=4) %>% dplyr::mutate(disorder="DEP"),
  cor_pc_scz %>% dplyr::select(1,cor=2,p=3,fdr=4) %>% dplyr::mutate(disorder="SCZ")
) 

cor_pc_res <- rbind(
  cor_pc,
  cor_pc %>% 
    group_by(disorder,module) %>% 
    dplyr::mutate(fdr = p.adjust(p,method = "BH")) %>% 
    ungroup() %>% 
    group_by(module) %>% 
    dplyr::mutate(cor = case_when(max(cor)<=0 ~ max(cor),
                                  min(cor)>=0 ~ min(cor),
                                  TRUE ~ NA)) %>% 
    dplyr::mutate(p = ifelse(is.na(cor),NA,max(p)),
                  fdr = ifelse(is.na(cor),NA,max(fdr))) %>% 
    dplyr::select(module,cor,p,fdr) %>% 
    dplyr::mutate(disorder="consensus") %>% unique()
)

## 5.1.5 cb ----
cor_cb_bd <- fread("2.output/4.ndd_wgcna/data/04_ModuleTraitCorRes/CB/04_ModuleTraitCorRes-BD.txt")
cor_cb_dep <- fread("2.output/4.ndd_wgcna/data/04_ModuleTraitCorRes/CB/04_ModuleTraitCorRes-DEP.txt")
cor_cb_scz <- fread("2.output/4.ndd_wgcna/data/04_ModuleTraitCorRes/CB/04_ModuleTraitCorRes-SCZ.txt")

cor_cb <- rbind(
  cor_cb_bd %>% dplyr::select(1,cor=2,p=3,fdr=4) %>% dplyr::mutate(disorder="BD"),
  cor_cb_dep %>% dplyr::select(1,cor=2,p=3,fdr=4) %>% dplyr::mutate(disorder="DEP"),
  cor_cb_scz %>% dplyr::select(1,cor=2,p=3,fdr=4) %>% dplyr::mutate(disorder="SCZ")
) 

cor_cb_res <- rbind(
  cor_cb,
  cor_cb %>% 
    group_by(disorder,module) %>% 
    dplyr::mutate(fdr = p.adjust(p,method = "BH")) %>% 
    ungroup() %>% 
    group_by(module) %>% 
    dplyr::mutate(cor = case_when(max(cor)<=0 ~ max(cor),
                                  min(cor)>=0 ~ min(cor),
                                  TRUE ~ NA)) %>% 
    dplyr::mutate(p = ifelse(is.na(cor),NA,max(p)),
                  fdr = ifelse(is.na(cor),NA,max(fdr))) %>% 
    dplyr::select(module,cor,p,fdr) %>% 
    dplyr::mutate(disorder="consensus") %>% unique()
)

## 5.1.6 hip ----
cor_hip_bd <- fread("2.output/4.ndd_wgcna/data/04_ModuleTraitCorRes/HIP/04_ModuleTraitCorRes-BD.txt")
cor_hip_dem <- fread("2.output/4.ndd_wgcna/data/04_ModuleTraitCorRes/HIP/04_ModuleTraitCorRes-DEM.txt")
cor_hip_mdd <- fread("2.output/4.ndd_wgcna/data/04_ModuleTraitCorRes/HIP/04_ModuleTraitCorRes-MDD.txt")
cor_hip_scz <- fread("2.output/4.ndd_wgcna/data/04_ModuleTraitCorRes/HIP/04_ModuleTraitCorRes-SCZ.txt")

cor_hip <- rbind(
  cor_hip_bd %>% dplyr::select(1,cor=2,p=3,fdr=4) %>% dplyr::mutate(disorder="BD"),
  cor_hip_dem %>% dplyr::select(1,cor=2,p=3,fdr=4) %>% dplyr::mutate(disorder="DEM"),
  cor_hip_mdd %>% dplyr::select(1,cor=2,p=3,fdr=4) %>% dplyr::mutate(disorder="MDD"),
  cor_hip_scz %>% dplyr::select(1,cor=2,p=3,fdr=4) %>% dplyr::mutate(disorder="SCZ")
) 

cor_hip_res <- rbind(
  cor_hip,
  cor_hip %>% 
    group_by(disorder,module) %>% 
    dplyr::mutate(fdr = p.adjust(p,method = "BH")) %>% 
    ungroup() %>% 
    group_by(module) %>% 
    dplyr::mutate(cor = case_when(max(cor)<=0 ~ max(cor),
                                  min(cor)>=0 ~ min(cor),
                                  TRUE ~ NA)) %>% 
    dplyr::mutate(p = ifelse(is.na(cor),NA,max(p)),
                  fdr = ifelse(is.na(cor),NA,max(fdr))) %>% 
    dplyr::select(module,cor,p,fdr) %>% 
    dplyr::mutate(disorder="consensus") %>% unique()
)

## 5.1.7 str ----
cor_str_bd <- fread("2.output/4.ndd_wgcna/data/04_ModuleTraitCorRes/STR/04_ModuleTraitCorRes-BD.txt")
cor_str_mdd <- fread("2.output/4.ndd_wgcna/data/04_ModuleTraitCorRes/STR/04_ModuleTraitCorRes-MDD.txt")
cor_str_scz <- fread("2.output/4.ndd_wgcna/data/04_ModuleTraitCorRes/STR/04_ModuleTraitCorRes-SCZ.txt")

cor_str <- rbind(
  cor_str_bd %>% dplyr::select(1,cor=2,p=3,fdr=4) %>% dplyr::mutate(disorder="BD"),
  cor_str_mdd %>% dplyr::select(1,cor=2,p=3,fdr=4) %>% dplyr::mutate(disorder="MDD"),
  cor_str_scz %>% dplyr::select(1,cor=2,p=3,fdr=4) %>% dplyr::mutate(disorder="SCZ")
) 

cor_str_res <- rbind(
  cor_str,
  cor_str %>% 
    group_by(disorder,module) %>% 
    dplyr::mutate(fdr = p.adjust(p,method = "BH")) %>% 
    ungroup() %>% 
    group_by(module) %>% 
    dplyr::mutate(cor = case_when(max(cor)<=0 ~ max(cor),
                                  min(cor)>=0 ~ min(cor),
                                  TRUE ~ NA)) %>% 
    dplyr::mutate(p = ifelse(is.na(cor),NA,max(p)),
                  fdr = ifelse(is.na(cor),NA,max(fdr))) %>% 
    dplyr::select(module,cor,p,fdr) %>% 
    dplyr::mutate(disorder="consensus") %>% unique()
)
## 5.1.8 amy ----
cor_amy_mdd <- fread("2.output/4.ndd_wgcna/data/04_ModuleTraitCorRes/AMY/04_ModuleTraitCorRes-MDD.txt")
cor_amy_res <- cor_amy_mdd %>% dplyr::select(1,cor=2,p=3,fdr=4) %>% dplyr::mutate(disorder="MDD")

## 5.1.9 save results ----
cor_ndd <- mget(ls(pattern = "cor_.*_res"))
names(cor_ndd) <- toupper(gsub("cor_|_res","",names(cor_ndd)))

save(cor_ndd,file = "./2.output/4.ndd_wgcna/data/04_ModuleTraitCorRes/ModuleTraitbiCorRes_NDD.Rdata")


# 5.2 Exporting results of the network analysis -----
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
  GSmat = as.data.frame(cbind(do.call(cbind,map(GS,"bicor")),
                              do.call(cbind,map(GS,"p")),
                              GS.metaZ,GS.metaP))
  colnames(GSmat)[1:(nSets*2)] <- paste0(setLabels,rep(c(".GS.bicor",".GS.p"),c(nSets,nSets)))
  
  # Same code for kME:
  kMEmat = as.data.frame(cbind(do.call(cbind,map(kME,"bicor")),
                               do.call(cbind,map(kME,"p")),
                               kME.metaZ,kME.metaP))
  
  ncol_kMEmat = ncol(kMEmat)
  n_ME = ncol(kME[[1]]$bicor)
  colnames(kMEmat)[-c(ncol_kMEmat-1,ncol_kMEmat)] <- paste0(rep(setLabels,
                                                                each = n_ME),
                                                            rep(c(".bicor.",".p."),
                                                                c(nSets*n_ME,nSets*n_ME)),
                                                            colnames(kMEmat)[-c(ncol_kMEmat-1,ncol_kMEmat)])
  
  info = data.frame(GeneSymbol = rownames(GSmat),
                    ModuleColor = moduleColors,
                    GSmat,
                    kMEmat)
  
  if(is.null(outdir)==FALSE){
    fwrite(info,paste0(outdir,"/05_ConsensusAnalysis-CombinedNetworkResults-",region,".csv"))
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
  GSmat = do.call(cbind,GS)
  colnames(GSmat) <- names(GS)
  
  # Same code for kME:
  kMEmat = do.call(cbind,kME[1:3])
  
  
  colnames(kMEmat) <- paste0(colnames(kMEmat),
                             rep(c(".bicor",".p",".Z"),
                                 c(ncol(kME$bicor),ncol(kME$bicor),ncol(kME$bicor))
                                 )
                             )
  
  info = data.frame(GeneSymbol = rownames(GSmat),
                    ModuleColor = moduleColors,
                    GSmat,
                    kMEmat)
  
  if(is.null(outdir)==FALSE){
    fwrite(info,paste0(outdir,"/05_NonConsensusAnalysis-CombinedNetworkResults-",region,".csv"))
  }
  return(info) 
}
## 5.1 FCx ----
load(file = "2.output/2.ndd_wgcna/data/02_Consensus_dataInput_FC-re.RData")
load(file = "2.output/2.ndd_wgcna/data/03_Consensus-NetworkConstruction-FCx.RData")
nSets_FC = exprSize_FC$nSets

info_FC = get_network_GS_kME(multiExpr = multiExpr_FC,
                   moduleColors = moduleColors_FC,
                   Traits = Traits_FC,
                   nSets = nSets_FC,
                   setLabels = setLabels_FC,
                   outdir = "./2.output/2.ndd_wgcna/data/",region = "FCx")

## 5.2 TCx ----
load(file = "2.output/2.ndd_wgcna/data/02_Consensus_dataInput_TC.RData")
load(file = "2.output/2.ndd_wgcna/data/03_Consensus-NetworkConstruction-TCx.RData")
nSets_TC = exprSize_TC$nSets

info_TC = get_network_GS_kME(multiExpr = multiExpr_TC,
                   moduleColors = moduleColors_TC,
                   Traits = Traits_TC,
                   nSets = nSets_TC,
                   setLabels = setLabels_TC,
                   outdir = "./2.output/2.ndd_wgcna/data/",region = "TCx")

## 5.3 OCx ----
load(file = "2.output/2.ndd_wgcna/data/02_NonConsensus_dataInput_OC_asd.RData")
load(file = "2.output/2.ndd_wgcna/data/03_NonConsensus-NetworkConstruction-OCx.RData")

info_OC = get_1network_GS_kME(Expr = asd_OC,
                              MEs = MEs_OC,
                              moduleColors = moduleColors_OC,
                              Traits = Traits_OC,
                              outdir = "./2.output/2.ndd_wgcna/data/",region = "OCx")

## 5.4 PCx ----
load(file = "2.output/2.ndd_wgcna/data/02_Consensus_dataInput_PC.RData")
load(file = "2.output/2.ndd_wgcna/data/03_Consensus-NetworkConstruction-PCx.RData")
nSets_PC = exprSize_PC$nSets

info_PC = get_network_GS_kME(multiExpr = multiExpr_PC,
                             moduleColors = moduleColors_PC,
                             Traits = Traits_PC,
                             nSets = nSets_PC,
                             setLabels = setLabels_PC,
                             outdir = "./2.output/2.ndd_wgcna/data/",region = "PCx")

## 5.5 CB ----
load(file = "2.output/2.ndd_wgcna/data/02_Consensus_dataInput_CB-re.RData")
load(file = "2.output/2.ndd_wgcna/data/03_Consensus-NetworkConstruction-CB.RData")
nSets_CB = exprSize_CB$nSets

info_CB = get_network_GS_kME(multiExpr = multiExpr_CB,
                             moduleColors = moduleColors_CB,
                             Traits = Traits_CB,
                             nSets = nSets_CB,
                             setLabels = setLabels_CB,
                             outdir = "./2.output/2.ndd_wgcna/data/",region = "CB")

## 5.6 HIP ----
load(file = "2.output/2.ndd_wgcna/data/02_Consensus_dataInput_HIP.RData")
load(file = "2.output/2.ndd_wgcna/data/03_Consensus-NetworkConstruction-HIP.RData")
nSets_HIP = exprSize_HIP$nSets

info_HIP = get_network_GS_kME(multiExpr = multiExpr_HIP,
                             moduleColors = moduleColors_HIP,
                             Traits = Traits_HIP,
                             nSets = nSets_HIP,
                             setLabels = setLabels_HIP,
                             outdir = "./2.output/2.ndd_wgcna/data/",region = "HIP")
## 5.7 STR ----
load(file = "2.output/4.ndd_wgcna/data/02_Consensus_dataInput_STR.RData")
load(file = "2.output/4.ndd_wgcna/data/03_Consensus-NetworkConstruction-STR.RData")
nSets_STR = exprSize_STR$nSets

info_STR = get_network_GS_kME(multiExpr = multiExpr_STR,
                              moduleColors = moduleColors_STR,
                              Traits = Traits_STR,
                              nSets = nSets_STR,
                              setLabels = setLabels_STR,
                              outdir = "./2.output/4.ndd_wgcna/data/",region = "STR")
## 5.8 AMY ----
load(file = "2.output/2.ndd_wgcna/data/02_NonConsensus_dataInput_AMY_mdd.RData")
load(file = "2.output/2.ndd_wgcna/data/03_NonConsensus-NetworkConstruction-AMY.RData")

info_AMY = get_1network_GS_kME(Expr = mdd_AMY,
                              MEs = MEs_AMY,
                              moduleColors = moduleColors_AMY,
                              Traits = Traits_AMY,
                              outdir = "./2.output/2.ndd_wgcna/data/",region = "AMY")

