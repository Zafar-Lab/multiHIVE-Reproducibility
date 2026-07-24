# https://github.com/ctlab/fgsea
# https://stephenturner.github.io/deseq-to-fgsea/
library("extrafont")
library(data.table)
library(fgsea)
library(ggplot2)
library(dplyr)
library(tidyverse)
library(ComplexHeatmap)
gseaenrichment_plot <- function(fgseaRes, save, topbottom=20){
  fgseaRes_plot <-  rbind(fgseaRes %>% top_n(topbottom, NES), fgseaRes %>% top_n(-topbottom, NES))
  gseaenrichment <- ggplot(data = fgseaRes_plot) +
    geom_bar(aes(x = reorder(pathway, NES), y = NES, fill = padj), stat = "identity", show.legend = T) +
    xlab("Pathway") +
    scale_fill_gradient(low = "blue", high = "red") + 
    coord_flip()+
    theme(plot.title = element_text(hjust = 0.5))
  # ggtitle(paste0("k_",i))+
  print(gseaenrichment)
  ggsave(save, plot = gseaenrichment, width = 15, height = 10, units = 'in')
}

gene_spectrum_scores_path = "/Biological_Analysis/Breast_Cancer/cNMF/breast_cancer-Z1/bc_cNMF/bc_cNMF.gene_spectra_score.k_8.dt_0_02.txt"

results_path = "./Results/Z1/"

gene_spectrum_scores = read.csv2(gene_spectrum_scores_path, sep="\t",row.names = 1)
gene_spectrum_scores = as.data.frame(sapply(gene_spectrum_scores,as.numeric))

for (i in 1:dim(gene_spectrum_scores)[1]){
  
  scores <- as.data.frame(t(gene_spectrum_scores[i,]))
  colnames(scores)= c('scores')
  scores = scores %>% arrange(desc(scores))
  
  ranks = as.numeric(scores$scores)
  names(ranks) = rownames(scores)
  
  #gobp pathways
  pathways = gmtPathways("c5.go.bp.v2023.2.Hs.symbols.gmt")
  fgseaRes <- fgsea(pathways = pathways, 
                    stats    = ranks,
                    eps      = 0.0, # get p values more accurate
                    minSize  = 15,
                    maxSize  = 500)
  fgseaRes <- fgseaRes %>% arrange(NES) %>% filter(padj<0.05) %>% arrange(NES)
  write.csv(apply(fgseaRes, 2, as.character), file = paste0(results_path,"GOBP_k",i,".csv"))
  gseaenrichment_plot(fgseaRes, paste0(results_path,"GOBP_k", i, ".png"))
  
}


#####################################################################################
fgseaRes <- read.csv2("./Results/Z1/GOBP_k8.csv", sep = ',', row.names = 1)
res_pos <- fgseaRes[fgseaRes$NES>0, ] %>% arrange(desc(NES))
res_pos$pathway[grep("different", res_pos$pathway, ignore.case = T)]
pathways = c( "GOBP_T_CELL_ACTIVATION", "GOBP_T_CELL_DIFFERENTIATION", "GOBP_ALPHA_BETA_T_CELL_ACTIVATION",
             "GOBP_MACROPHAGE_ACTIVATION","GOBP_MYELOID_LEUKOCYTE_ACTIVATION","GOBP_MYELOID_LEUKOCYTE_DIFFERENTIATION", 
             "GOBP_MUSCLE_CELL_DEVELOPMENT", "GOBP_STRIATED_MUSCLE_CELL_DEVELOPMENT", "GOBP_MUSCLE_CELL_DIFFERENTIATION",
             "GOBP_ENDOTHELIAL_CELL_DEVELOPMENT", "GOBP_ENDOTHELIUM_DEVELOPMENT",
             "GOBP_CELL_MATRIX_ADHESION", "GOBP_COLLAGEN_FIBRIL_ORGANIZATION", "GOBP_COLLAGEN_METABOLIC_PROCESS",
             "GOBP_CELL_SUBSTRATE_ADHESION", "GOBP_NEGATIVE_REGULATION_OF_BMP_SIGNALING_PATHWAY",
             "GOBP_CELL_CYCLE_PROCESS", "GOBP_REGULATION_OF_CELL_CYCLE_PROCESS",
             "GOBP_RNA_CATABOLIC_PROCESS"
             )
row_an = c(rep("Program-1",3),rep("Program-2",3),rep("Program-3",3),rep("Program-4",2),rep("Program-5",3),rep("Program-6",2),rep("Program-7",2),rep("Program-8",1))
pathway_prog <- data.frame("pathway" = pathways) 
for (k in 1:8){
  fgseaRes <- read.csv2(paste0("./Results/Z1/GOBP_k",k,".csv"), sep = ',', row.names = 1)
  col_name = paste0("Program-",k)
  fgseaRes = fgseaRes[fgseaRes$pathway %in% pathways, c("pathway","NES")] %>% rename(!!col_name:="NES")
  fgseaRes[,col_name] = as.numeric(fgseaRes[,col_name])
  pathway_prog = merge(pathway_prog,fgseaRes, by="pathway", all.x = T)
}
pathway_prog = pathway_prog %>% remove_rownames %>% column_to_rownames("pathway") %>% as.matrix()
pathway_prog = pathway_prog[pathways,]
# pathway_prog[is.na(pathway_prog)] = NA
row_an = rowAnnotation(Program = as.factor(row_an))
rownames(pathway_prog) = sapply(rownames(pathway_prog), function (x){str_replace(x,"GOBP_","")})
png("./Results/Z1.png",  width = 8, height = 7, units = 'in', fonts = "Arial", res = 300)
set.seed(0)
Heatmap(pathway_prog, cluster_rows = F, cluster_columns = F, 
        left_annotation = row_an, name = "NES",
        row_names_max_width =  max_text_width(rownames(pathway_prog), gp = gpar(fontsize = 10) )
       ) 
dev.off()
#####################################################################################

# head(fgseaRes[order(pval), ])
# http://127.0.0.1:14289/graphics/plot_zoom_png?width=762&height=401
# plotEnrichment(examplePathways[["5991130_Programmed_Cell_Death"]],
#                exampleRanks) + labs(title="Programmed Cell Death")
# 
# topPathwaysUp <- fgseaRes[ES > 0][head(order(pval), n=10), pathway]
# topPathwaysDown <- fgseaRes[ES < 0][head(order(pval), n=10), pathway]
# topPathways <- c(topPathwaysUp, rev(topPathwaysDown))
# 
# plotGseaTable(examplePathways[topPathways], exampleRanks, fgseaRes, 
#               gseaParam=0.5)
# 
# 
# ggplot(fgseaResTidy, aes(reorder(pathway, NES), NES)) +
#   geom_col(aes(fill=padj<0.05)) +
#   coord_flip() +
#   labs(x="Pathway", y="Normalized Enrichment Score",
#        title="Hallmark pathways NES from GSEA") + 
#   theme_minimal()