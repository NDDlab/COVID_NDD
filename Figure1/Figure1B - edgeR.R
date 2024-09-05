library(tidyverse)
library(data.table)

##########################################
#                                1. PCA
##########################################

library(FactoMineR)
library(factoextra)

# 1. import ----
tpm <- fread("./1.rawdata/COVID/GSE185557/2.output/GSE185557_tpm.txt")
sample <- fread("./1.rawdata/COVID/GSE185557/2.output/sample_GSE185557_tidy.txt")

# 2. selete the unlabor maternal blood, preprocess data ----
sample_filt <- sample %>% dplyr::filter(str_detect(source_name,"Maternal "))
sample_filt$deliveryroute[is.na(sample_filt$deliveryroute)] <-"Cesarean section"
tpm <- tpm %>% dplyr::select(gene,sample_filt$title) %>% column_to_rownames("gene")

tpm.pca <- PCA(t(tpm), graph = T,scale.unit = T)

p1<-fviz_pca_ind(tpm.pca,
                 mean.point=F,
                 lable="none",
                 addEllipses = TRUE,
                 ellipse.type = "convex",geom = "point",
                 habillage=factor(sample_filt$group[match(colnames(tpm),sample_filt$title)]),
                 palette = RColorBrewer::brewer.pal(8,"Set2")[c(1:2)])

p2<-fviz_pca_ind(tpm.pca,
                 mean.point=F,
                 lable="none",
                 addEllipses = TRUE,
                 ellipse.type = "convex",geom = "point",
                 habillage=factor(sample_filt$deliveryroute[match(colnames(tpm),sample_filt$title)]),
                 palette = RColorBrewer::brewer.pal(8,"Set2")[c(3:4)])

library(patchwork)
p1+theme_bw()+p2+theme_bw()+plot_layout(guides = "collect")

##########################################
#                                2. data processing
##########################################
# 1. import data ----
count <- read_csv("./1.rawdata/COVID/GSE185557/1.rawdata/GSE185557_count_matrix.csv")
sample <- fread("./1.rawdata/COVID/GSE185557/2.output/sample_GSE185557_tidy.txt")

# 2. selete the unlabor maternal blood, filter outlier sample (by PCA), preprocess data ----
sample_filt <- sample %>% dplyr::filter(str_detect(source_name,"Maternal "))
sample_filt %<>% dplyr::filter(title!="MB-C-6")

sample_filt$deliveryroute[is.na(sample_filt$deliveryroute)] <-"Cesarean section"
count_filt <- count %>% dplyr::select(ID_REF,sample_filt$title)

# 3. transfer ensg id to symbol
hg38 <- fread("./1.rawdata/COVID/GSE185557/1.rawdata/gencode.v39.annotation.gtf.gz") %>% 
  dplyr::select(V3,V9) %>% 
  dplyr::filter(V3=="gene")

ensg2symbol = cbind(ID_REF = str_extract(hg38$V9,'ENSG\\d+'),
                    gene = gsub('"',
                                '',
                                gsub('gene_name "','',gsub(";.*","",str_extract(hg38$V9,'gene_name\\s".*";\\s')),fixed = T),
                                fixed = T)
                    ) %>% as.data.frame() %>% dplyr::filter(!str_detect(gene,"ENSG")) %>% unique() %>% na.omit()

count_filt <- merge(ensg2symbol,count_filt,by='ID_REF',all.y = T) %>% 
  na.omit() %>% 
  distinct(gene,.keep_all = T) %>% 
  dplyr::select(-ID_REF) %>% 
  column_to_rownames("gene")

save(count_filt,sample_filt,file="./2.output/1.covid_maternal_diff/data/4_expMat_and_sample for DEanalysis.Rdata")

##########################################
#                                3. differential analysis
##########################################
library(edgeR)

load("./2.output/1.covid_maternal_diff/data/4_expMat_and_sample for DEanalysis.Rdata")

# 1. make obj
group_info = factor(sample_filt$group,levels = c("Control","COVID-19"))

dge.list.obj <- DGEList(counts = count_filt, group = group_info)
keep <- filterByExpr(dge.list.obj) # filter low exp level gene
dge.list.obj <- dge.list.obj[keep,,keep.lib.sizes=FALSE]

# 2. normalization
dge.list.obj <- calcNormFactors(dge.list.obj,method = "TMM")
#design.mat <- model.matrix(~group_info+sample_filt$deliveryroute+sample_filt$age+sample_filt$bmi)
design.mat <- model.matrix(~group_info)
dge.list.obj <- estimateDisp(dge.list.obj,design.mat)

# 3.likelihood ratio test
dge.list.res <- exactTest(dge.list.obj)
DEGs.res <- as.data.frame(topTags(dge.list.res,n=nrow(count),sort.by = "logFC"))

# 4. merge diff result with expmat
count_diff<-merge(count_filt %>% rownames_to_column("gene"),
                  DEGs.res %>% rownames_to_column("gene"),
                  by="gene")

Nup<-count_diff %>% dplyr::filter(logFC > 1 & PValue < 0.05) %>% nrow()
Ndown<-count_diff %>% dplyr::filter(logFC < -1 & PValue < 0.05) %>% nrow()
count_diff$type <- with(count_diff,case_when(logFC > 1 & PValue < 0.05 ~ paste0('Up ',Nup),
                                             logFC < -1 & PValue < 0.05 ~ paste0('Down ',Ndown),
                                             TRUE ~ paste0('')))

fwrite(count_diff,"./2.output/1.covid_maternal_diff/data/5_diffRes(status+age+bmi+deliveryroute).txt",sep = "\t")