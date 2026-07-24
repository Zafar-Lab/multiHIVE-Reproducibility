library(Seurat)
library('EnhancedVolcano')
library(dplyr)
library(viridis)
library(fgsea)
library(tidyverse)
library(ComplexHeatmap)

sceasy::convertFormat("./Results/breast_cancer_z1_cnmf.h5ad", from="anndata", to="seurat",
                      outFile='./Results/breast_cancer_z1_cnmf.rds')

data = readRDS("./Results/breast_cancer_z1_cnmf.rds")
Idents(data) <- data$P56
sum(data@assays$RNA@counts)
data <- NormalizeData(data, normalization.method = "LogNormalize", scale.factor = 10000)
data <- ScaleData(data, features = rownames(data)) #needed for heatmap
data <- subset(data, subset = celltype_major == "CAFs")
sum(data@assays$RNA@data)


png("./Results/p56-expression.png", width = 7, height = 7, res = 200, units = 'in')
data$P56_CAF = as.character(data$P56)
data$P56_CAF[data$P56_CAF=="Program-5"] = "CAF_P5"
data$P56_CAF[data$P56_CAF=="Program-6"] = "CAF_P6"
DimPlot(data, group.by = 'P56_CAF')
dev.off()
de <- FindMarkers(data,ident.1 = "Program-5", ident.2 = "Program-6",
                  logfc.threshold = 0, min.pct = 0)

de <- de %>%  mutate_at(vars(p_val_adj, pct.1, pct.2), as.numeric)%>% filter(p_val_adj < 0.05) %>%
  filter(pct.1 > 0.01 | pct.2 > 0.01 )

png("./Results/Volcano-p5.png", width = 7, height = 7, res = 200, units = 'in')
fig <- EnhancedVolcano(de,
                lab = rownames(de),
                x = 'avg_log2FC',
                y = 'p_val_adj',
                FCcutoff = 2,
                pCutoff = 0.05)
print(fig)
dev.off()
fig
de %>% filter(p_val_adj < 0.05)  %>% arrange(desc(avg_log2FC)) %>%
  slice_tail(n = 10) %>%
  ungroup() -> bottom10
de  %>% filter(p_val_adj < 0.05)  %>% arrange(desc(avg_log2FC)) %>%
  slice_head(n = 10) %>%
  ungroup() -> top10
top.de = rbind(top10, bottom10)




scores <- de["avg_log2FC"]
scores = scores %>% arrange(desc(scores))
ranks = as.numeric(scores$avg_log2FC)
names(ranks) = rownames(scores)
results_path = "./Results/"
#hallmark

pathways <- gmtPathways("/GSEA/databases/human/h.all.v2023.2.Hs.symbols.gmt")
fgseaRes <- fgsea(pathways = pathways,
                  stats    = ranks,
                  eps      = 0.0, # get p values more accurate
                  minSize  = 15,
                  maxSize  = 500)
fgseaRes <- fgseaRes %>% arrange(NES) %>% filter(padj<0.05) %>% arrange(NES)
write.csv(apply(fgseaRes, 2, as.character), file = paste0(results_path,"hallmark_p5.csv"))


#####################################################################################
pathways = c("GAVISH_3CA_METAPROGRAM_FIBROBLASTS_CAF_2", "GAVISH_3CA_METAPROGRAM_FIBROBLASTS_LIPID_METABOLISM", "GAVISH_3CA_METAPROGRAM_FIBROBLASTS_PERICYTE_LIKE",
             "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION", "GAVISH_3CA_METAPROGRAM_FIBROBLASTS_HYPOXIA", "GAVISH_3CA_METAPROGRAM_FIBROBLASTS_CAF_7", "GAVISH_3CA_METAPROGRAM_FIBROBLASTS_CAF_1")
pathway_prog <- data.frame("pathway" = pathways) 

fgseaRes <- read.csv2("./Results/all_p5.csv", sep = ',', row.names = 1)
fgseaRes = fgseaRes[fgseaRes$pathway %in% pathways, c("pathway","NES")] %>% rename(!!"CAF_P5":="NES")
fgseaRes[,"CAF_P5"] = as.numeric(fgseaRes[,"CAF_P5"])
pathway_prog = merge(pathway_prog,fgseaRes, by="pathway", all.x = T)

pathway_prog[,"CAF_P6"] = -as.numeric(pathway_prog[,"CAF_P5"])

pathway_prog = pathway_prog %>% remove_rownames %>% column_to_rownames("pathway") %>% as.matrix()
pathway_prog = pathway_prog[pathways,]
# row_an = c(rep("Program-5",3),rep("Program-6",4))
# row_an = rowAnnotation(Program = as.factor(row_an))
rownames(pathway_prog) = sapply(rownames(pathway_prog), function (x){str_replace(x,"GAVISH_3CA_","")})
rownames(pathway_prog) = sapply(rownames(pathway_prog), function (x){str_replace(x,"HALLMARK_","")})

png("./Results/P56-pathway-heatmap.png",  width = 5.5, height = 3.3, units = 'in', fonts = "Arial", res = 300)
set.seed(0)
Heatmap(pathway_prog, cluster_rows = F, cluster_columns = F, 
        name = "NES",
        row_names_max_width =  1.2*max_text_width(rownames(pathway_prog), gp = gpar(fontsize = 10) )
) 
dev.off()
#####################################################################################
