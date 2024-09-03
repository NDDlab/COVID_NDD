##########################################
#
# Title : Get The Maternal COVID-19 Exposure DEGs, by edgeR
# Author : Yin-Huamin
# Date : 2023-02-25
#
######################

# 1. import data ----

count <- read_csv("./1.rawdata/COVID/GSE185557/1.rawdata/GSE185557_count_matrix.csv")

sample <- fread("./1.rawdata/COVID/GSE185557/2.output/sample_GSE185557_tidy.txt")

# 2. selete the unlabor maternal blood, preprocess data ----
sample_filt <- sample %>% dplyr::filter(str_detect(source_name,"Maternal "))
sample_filt %<>% dplyr::filter(title!="MB-C-6")

sample_filt$deliveryroute[is.na(sample_filt$deliveryroute)] <-"Cesarean section"
count_filt <- count %>% dplyr::select(ID_REF,sample_filt$title)

# transfer ensg id to symbol
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


# 3. differential analysis ----
rm(list=ls());gc()
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

fwrite(count_diff,"./2.output/1.covid_maternal_diff/data/5_diffRes.txt",sep = "\t")
#fwrite(count_diff,"./2.output/1.covid_maternal_diff/data/5_diffRes(status+age+bmi+deliveryroute).txt",sep = "\t")

# 5. Visualization( Volcano + Heatmap) ----
(list = ls());gc()
library(ggrepel)
library(ggplotify)

# 1. volcano
#diff<-fread("./2.output/1.covid_maternal_diff/data/5_diffRes.txt") %>% dplyr::select(gene,logFC,PValue,type)
diff<-fread("./2.output/1.covid_maternal_diff/data/5_diffRes(status+age+bmi+deliveryroute).txt") %>% dplyr::select(gene,logFC,PValue,type)

diff.sig<-diff %>%
  dplyr::filter(str_detect(type,"\\w")) %>% # case or control
  arrange(PValue) %>%
  group_by(type) %>% dplyr::filter(row_number()<6) # top 5 sig genes in up and down, respectively

ggplot(diff, aes(x = logFC,y = -log10(PValue), colour=type)) +
  
  geom_point(alpha=0.5, size=2) +
  
  scale_color_manual(values=c("grey","blue","red"))+
  
  geom_vline(xintercept=c(-1,1),lty=4,col="black",lwd=0.8) +
  geom_hline(yintercept = -log10(0.05),lty=4,col="black",lwd=0.8) +
  
  labs(x="log2(fold change)",
       y="-log10 (FDR)")+
  
  coord_cartesian(xlim = c(-8,8),ylim = c(0,12))+
  
  theme_bw()+
  
  theme(plot.title = element_text(hjust = 0.5),
        legend.position="right",
        legend.title = element_blank(),
        legend.text = element_text(size=12)
  ) +
  
  geom_label_repel(data = diff.sig,mapping = aes(x=logFC,v=-log10(PValue),label=gene),size = 5,force = 3,key_glyph="point")

ggsave("./2.output/1.covid_maternal_diff/figure/11_diff_volcano.pdf",width = 5)
ggsave("./2.output/1.covid_maternal_diff/figure/11_diff_volcano.png",width = 5,dpi=600)

# 2. heatmap
library(ComplexHeatmap)
diffmat<-fread("./2.output/1.covid_maternal_diff/data/5_diffRes(status+age+bmi+deliveryroute).txt")
diffmat[,2:21]<-as.data.frame(cpm(diffmat[,2:21]))

#PValue<0.05,!str_detect(type,"NoSig")
diffmat %<>% dplyr::filter(type!="")  %>% dplyr::select(gene,1:21)
diffmat %<>% column_to_rownames("gene")

for (i in 1:nrow(diffmat)) {
  diffmat[i,]<-scale(as.numeric(diffmat[i,]),center = T,scale = T)
}


column_anno<-HeatmapAnnotation(
  DataType=sample_filt$group,
  col=list(
    DataType=c("COVID-19"="#EE7E30","Control"="#5D9AD3")
  )
)

Heatmap(diffmat,cluster_rows = T,cluster_columns = F,
        show_row_names = F,show_column_names = T,name="Expression",
        top_annotation = column_anno,border = T) %>% as.ggplot()


ggsave("./2.output/1.covid_maternal_diff/figure/12_diff_heatmap.pdf",width = 5,height = 4)
ggsave("./2.output/1.covid_maternal_diff/figure/12_diff_heatmap.png",width = 5,height = 4,dpi=300)



diff<-fread("./2.output/1.covid_maternal_diff/data/5_diffRes(status+age+bmi+deliveryroute).txt") %>% dplyr::select(gene,logFC,PValue,type)

p_diff <- ggplot(diff, aes(x = logFC,y = -log10(PValue), fill=type,color=type)) +
  
  geom_point(alpha=0.6, size=2,pch=21) +
  
  scale_fill_manual(values=c("grey","#00468B","#ED0000"))+
  scale_color_manual(values=c("grey","#00468B","#ED0000"))+
  
  geom_vline(xintercept=c(-1,1),lty=4,col="black",lwd=0.8) +
  geom_hline(yintercept = -log10(0.05),lty=4,col="black",lwd=0.8) +
  
  labs(x="log2(fold change)",
       y="-log10 (FDR)")+
  
  coord_cartesian(xlim = c(-5,5),ylim = c(0,10))+
  
  theme_bw()+
  
  theme(plot.title = element_text(hjust = 0.5),
        legend.position="right",
        legend.title = element_blank(),
        legend.text = element_text(size=12)
  )

sig_igg <- diff %>% 
  dplyr::filter(type!="") %>% 
  dplyr::filter(str_detect(gene,"^IG")) %>% 
  dplyr::filter(gene!="IGF2BP3") %>% 
  dplyr::pull(gene)

# Immunoglobulin 

p_igg <- diff %>% 
  dplyr::mutate(iggtype = case_when(gene %in% sig_igg ~ "IGG",
                                    TRUE ~ "other")) %>% 
  ggplot(aes(x = logFC,y = -log10(PValue), fill=iggtype, colour=iggtype)) +
  
  geom_point(alpha=0.7, size=2) +
  
  scale_color_manual(values=c("purple","grey"))+
  scale_fill_manual(values=c("purple","grey"))+
  
  geom_vline(xintercept=c(-1,1),lty=4,col="black",lwd=0.8) +
  geom_hline(yintercept = -log10(0.05),lty=4,col="black",lwd=0.8) +
  
  labs(x="log2(fold change)",
       y="-log10 (FDR)")+
  
  coord_cartesian(xlim = c(-5,5),ylim = c(0,10))+
  
  theme_bw()+
  
  theme(plot.title = element_text(hjust = 0.5),
        legend.position="right",
        legend.title = element_blank(),
        legend.text = element_text(size=12)
  )

library(patchwork)
p_diff + p_igg+plot_layout(guides = "collect")

ggsave("./2.output/1.covid_maternal_diff/figure/11_diff_volcano_with_igg.pdf",width = 7.5,height = 3)
ggsave("./2.output/1.covid_maternal_diff/figure/11_diff_volcano_with_igg.png",width = 7.5,height = 3,dpi=300)

