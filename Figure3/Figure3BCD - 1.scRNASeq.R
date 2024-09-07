library(Seurat)
library(tidyverse)
library(tidydr)

######################################################
#                     1. import data, preprocess
######################################################

process_seurat <- function(data.dir, project_name, group_name){
  data <- Read10X(data.dir = data.dir) %>% 
    CreateSeuratObject(counts = ., project = project_name, min.cells = 5) %>% 
    NormalizeData() %>% 
    FindVariableFeatures(selection.method = "vst", nfeatures = 2000)
  data$group <- group_name
  return(data)
}

HC11 <- process_seurat(data.dir = "1.rawdata/MomBloodScRNA/GSE192693_RAW/GSM5761208_HC11_processed/",
                                     project_name = "HC11", group_name = "Gestational")

HC12 <- process_seurat(data.dir = "1.rawdata/MomBloodScRNA/GSE192693_RAW/GSM5761209_HC12_processed/",
                                     project_name = "HC12", group_name = "Postpartum")

HC21 <- process_seurat(data.dir = "1.rawdata/MomBloodScRNA/GSE192693_RAW/GSM5761210_HC21_processed/",
                                     project_name = "HC21", group_name = "Gestational")

HC22 <- process_seurat(data.dir = "1.rawdata/MomBloodScRNA/GSE192693_RAW/GSM5761211_HC22_processed/",
                                     project_name = "HC22", group_name = "Postpartum")

######################################################
#                     2. Perform integration, dimension reduction, cluster
######################################################

anchors <- FindIntegrationAnchors(object.list = list(HC11, HC12,HC21, HC22), dims = 1:20)
combined <- IntegrateData(anchorset = anchors, dims = 1:20)

# Perform an integrated analysis
DefaultAssay(combined) <- "integrated"

# Run the standard workflow for visualization and clustering
combined <- ScaleData(combined, verbose = FALSE)
combined <- RunPCA(combined, npcs = 30, verbose = FALSE)
combined <- RunUMAP(combined, reduction = "pca", dims = 1:20)
combined <- FindNeighbors(combined, reduction = "pca", dims = 1:20)
combined <- FindClusters(combined, resolution = 0.5)

######################################################
#                     3. Cell type annotation
######################################################

# find markers ----
DefaultAssay(combined)<-"RNA"
Idents(combined)<-combined@meta.data$integrated_snn_res.0.5
combined.markers <- FindAllMarkers(combined)
#saveRDS(combined.markers,file = "./2.output/19.scRNAscore/data/findmarkers.rds")

# auto ----
library(scMayoMap)

scMayoMap.obj <- scMayoMap(data = combined.markers, database = scMayoMapDatabase, tissue = "blood")
# saveRDS(scMayoMap.obj,file = "./2.output/19.scRNAscore/data/scMayoMap_obj.rds")

scMayoMap.plot(scMayoMap.object = scMayoMap.obj)+
  scale_size(range = c(0,15))+
  theme(axis.text.x = element_text(color = "black", size = 14,
                                   angle = -45,vjust = 1,hjust = 0))

scMayoMap.obj[["markers"]] %>% 
  group_by(cluster) %>% 
  dplyr::arrange(desc(score)) %>% 
  dplyr::filter(row_number()<3) %>% View()

scMayoMap.obj$markers$celltype %>% unique()

Idents(combined)<-combined@meta.data$integrated_snn_res.0.5

# define ----
combined <- RenameIdents(combined,
                         '0' = 'CD4 Naive T cell', 
                         '1' = 'CD4 Central Memory T cell', 
                         '2' = 'CD14 Monocyte', 
                         '3' = 'CD8 Effector Memory T cell', 
                         '4' = 'CD56-dim natural killer cell', 
                         '5' = 'Neutrophil', #other
                         '6' = 'CD4 Effector Memory T cell', 
                         '7' = 'Monocyte',  
                         '8' = 'Naive B cell', 
                         '9' = 'CD8 Naive T cell',
                         '10' = 'Intermediate B cell',  
                         '11' = 'Erythroid precursor cell', 
                         '12' = 'CD16 Monocyte', 
                         '13' = 'Gamma delta T cell', 
                         '14' = 'Platelet', 
                         '15' = 'Neutrophil', 
                         '16' = 'CD56-dim natural killer cell', 
                         '17' = 'Platelet', 
                         '18' = 'Dendritic cell',
                         '19' = 'CD4 Proliferating T cell',
                         '20' = 'CD4 T cell',# CD4 Effector Memory T cell
                         '21' = 'Erythroid precursor cell',
                         '22' = 'Plasmablast',
                         '23' = 'Neutrophil',
                         '24' = 'Hematopoietic stem cell')

combined@meta.data$annotation1 <- Idents(combined)

# check ----
DimPlot(combined, reduction = "umap")+
  theme_dr(xlength = 0.2,
           ylength = 0.2,
           arrow = grid::arrow(length = unit(0.1,"inches"),
                               ends = "last",
                               type = "closed"))+
  theme(panel.grid = element_blank()) + 
  ggtitle("")+
  guides(color=guide_legend(ncol = 2,override.aes = list(size=2)))

scMayoMap.obj[["markers"]] %>% 
  dplyr::filter(celltype %in% combined@meta.data$annotation1) %>% 
  group_by(cluster) %>% 
  dplyr::arrange(desc(score)) %>% 
  dplyr::filter(row_number()<2) %>% 
  ungroup() %>% 
  dplyr::select(celltype,genes) %>% 
  dplyr::filter(celltype %in% combined@meta.data$annotation2) %>% 
  separate_rows(genes,sep = ",") %>% 
  unique() %>% 
  group_by(genes) %>% 
  dplyr::mutate(n=n()) %>% 
  dplyr::mutate(celltype=ifelse(n>1,"Common",celltype)) %>% 
  split(.$celltype) %>% 
  map(.,"genes") -> marker_list

glist <- unique(unlist(marker_list))
glist_filt <- c("TRDC","IL1B","FCGR1A","CCR1","TREM1","BST1","C5AR1","ADPGK","PSTPIP1","AIM2",
                "CD4","TOP2A","ASPM","TPX2","CENPF","NUSAP1","OXNAD1","TRAC","TRDC","TNFRSF13B",
                "ITGAM","ITGAX","TLR2","CD33","CLEC4E","PILRA","KLRC1","BTG1","RGS10",
                "LYST","S100A4","SERPINB1","CD7","CCL5","GZMA","DOCK8","NCF1","PECAM1",
                "CXCR4","GZMH","LTB","IL32","IGHD","LYZ","S100A9","IGHM","BANK1","RALGPS2",
                "TCF7","PIK3IP1")

DotPlot(combined, group.by = "annotation1",
        features = glist,
        col.min = 0,
        cluster.idents = F) +
  scale_color_continuous(low="grey",high =  "red")+
  theme(axis.text.x  = element_text(color="black",size=10,angle = 45,vjust = 1, hjust=1),
        axis.text.y  = element_text(color="black",size=11),
        panel.border = element_rect(color="black"), 
        panel.spacing = unit(0, "mm"), 
        axis.line = element_blank(),
  )

combined@meta.data <- combined@meta.data %>% 
  dplyr::mutate(annotation2 = case_when(annotation1=="Platelet"~"Other",
                                        annotation1=="Erythroid precursor cell"~"Other",
                                        annotation1=="Dendritic cell"~"Other",
                                        annotation1=="Plasmablast"~"Other",
                                        annotation1=="Erythroblast"~"Other",
                                        annotation1=="Hematopoietic stem cell"~"Other",
                                        str_detect(annotation1,"killer cell")~"NK cell",
                                        str_detect(annotation1,"Monocyte")~"Monocyte",
                                        TRUE~annotation1)) 

DimPlot(combined, reduction = "umap",group.by = "annotation2")+
  scale_color_manual(values = c("CD4 Central Memory T cell"="#1A8B40",
                                "CD8 Effector Memory T cell"="#1A8B40",
                                "Neutrophil"="#D34A22",
                                "CD4 T cell"="#C16CAC",
                                "Other"="grey",
                                "Intermediate B cell"="#F37C7E",
                                "Naive B cell"="#D9A867",
                                "CD4 Effector Memory T cell"="#222A6A",
                                "CD4 Naive T cell"="#04737D",
                                "CD4 Proliferating T cell"="#D61921",
                                "Gamma delta T cell"="#8BA0D2",
                                "CD8 Naive T cell"="#E7C3DD",
                                "Monocyte"="#F47E27",
                                "NK cell"="#38BDA9"))+
  theme_dr(xlength = 0.2,
           ylength = 0.2,
           arrow = grid::arrow(length = unit(0.1,"inches"),
                               ends = "last",
                               type = "closed"))+
  theme(panel.grid = element_blank()) + ggtitle("")

Idents(combined)<-combined@meta.data$annotation2

DotPlot(combined, group.by = "annotation2",
        features = glist,
        col.min = 0,
        cluster.idents = F) +
  scale_color_continuous(low="grey",high =  "red")+
  theme(axis.text.x  = element_text(color="black",size=10,angle = 45,vjust = 1, hjust=1),
        axis.text.y  = element_text(color="black",size=11),
        panel.border = element_rect(color="black"), #面板边框
        panel.spacing = unit(0, "mm"), #面板间距
        axis.line = element_blank(),
)

######################################################
#                     4. Gene Set Score
######################################################
library(UCell)

quadruplet <- fread("2.output/15.covid_ndd_dev_gene/table/covid_br_ndd_wd_quadruplet.txt")
signature <- quadruplet$gene %>% unique()

set.seed(1)
combined <- AddModuleScore_UCell(combined,
                                                          features = list("RiskGenes"=signature),
                                                          assay = "RNA",
                                                          slot = "data",
                                                          chunk.size=1000)

meta.merge <- FetchData(combined,vars="RiskGenes")
combined <- AddMetaData(combined,meta.merge)
