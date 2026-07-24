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

gene_spectrum_scores_path = "/home/anirudhn/Krushna/breast_cancer/cNMF/breast_cancer-Z2/bc_cNMF/bc_cNMF.gene_spectra_score.k_7.dt_0_03.txt"

results_path = "./Results/Z2/"


gene_spectrum_scores = read.csv2(gene_spectrum_scores_path, sep="\t",row.names = 1)
gene_spectrum_scores = as.data.frame(sapply(gene_spectrum_scores,as.numeric))

for (i in 1:dim(gene_spectrum_scores)[1]){
  print(paste0("for i ", i))
  scores <- as.data.frame(t(gene_spectrum_scores[i,]))
  colnames(scores)= c('scores')
  scores = scores %>% arrange(desc(scores))
  
  ranks = as.numeric(scores$scores)
  names(ranks) = rownames(scores)
  
  #hallmark pathways
  pathways <- gmtPathways("/home/anirudhn/Krushna/GSEA/databases/h.all.v2023.2.Hs.symbols.gmt")
  fgseaRes <- fgsea(pathways = pathways,
                    stats    = ranks,
                    eps      = 0.0, # get p values more accurate
                    minSize  = 15,
                    maxSize  = 500)
  fgseaRes <- fgseaRes %>% arrange(NES) %>% filter(padj<0.05) %>% arrange(NES)
  write.csv(apply(fgseaRes, 2, as.character), file = paste0(results_path,"hallmark_k",i,".csv"))
  gseaenrichment_plot(fgseaRes, paste0(results_path,"hallmark_k", i, ".png"))
  
  #all pathways
  pathways = gmtPathways("/home/anirudhn/Krushna/GSEA/databases/msigdb.v2023.2.Hs.symbols.gmt")
  fgseaRes <- fgsea(pathways = pathways,
                    stats    = ranks,
                    eps      = 0.0, # get p values more accurate
                    minSize  = 15,
                    maxSize  = 500)
  fgseaRes <- fgseaRes %>% arrange(NES) %>% filter(padj<0.05) %>% arrange(NES)
  write.csv(apply(fgseaRes, 2, as.character), file = paste0(results_path,"all_k",i,".csv"))
  gseaenrichment_plot(fgseaRes, paste0(results_path,"all_k", i, ".png"))

  # gobp pathways
  pathways = gmtPathways("/home/anirudhn/Krushna/GSEA/databases/c5.go.bp.v2023.2.Hs.symbols.gmt")
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
fgseaRes <- read.csv2("./Results/Z2/GOBP_k3.csv", sep = ',', row.names = 1)
res_pos <- fgseaRes[fgseaRes$NES>0, ] %>% arrange(desc(NES))
res_pos$pathway[grep("vasc", res_pos$pathway, ignore.case = T)]
pathways = c("GOBP_T_CELL_RECEPTOR_SIGNALING_PATHWAY", "GOBP_ANTIGEN_RECEPTOR_MEDIATED_SIGNALING_PATHWAY", "GOBP_T_CELL_ACTIVATION", "GOBP_ALPHA_BETA_T_CELL_ACTIVATION",
             "GOBP_POSITIVE_REGULATION_OF_TRANSCRIPTION_BY_RNA_POLYMERASE_II","GOBP_POSITIVE_REGULATION_OF_RNA_METABOLIC_PROCESS", "GOBP_MRNA_METABOLIC_PROCESS", "GOBP_RNA_CATABOLIC_PROCESS",
              "GOBP_ARTERY_MORPHOGENESIS", "GOBP_CELLULAR_COMPONENT_ASSEMBLY_INVOLVED_IN_MORPHOGENESIS",
             "GOBP_ANTIGEN_PROCESSING_AND_PRESENTATION","GOBP_PHAGOCYTOSIS","GOBP_INFLAMMATORY_RESPONSE", "GOBP_MACROPHAGE_ACTIVATION",
             "GOBP_ENDODERMAL_CELL_DIFFERENTIATION","GOBP_ENDODERM_FORMATION", "GOBP_COLLAGEN_FIBRIL_ORGANIZATION", "GOBP_REGULATION_OF_BMP_SIGNALING_PATHWAY",
             "GOBP_EXTRACELLULAR_MATRIX_DISASSEMBLY", "GOBP_RESPONSE_TO_WOUNDING", "GOBP_BLOOD_VESSEL_MORPHOGENESIS", "GOBP_VASCULATURE_DEVELOPMENT",
             "GOBP_SPINDLE_ORGANIZATION","GOBP_CELL_CYCLE_CHECKPOINT_SIGNALING", "GOBP_CELL_DIVISION", "GOBP_CELL_CYCLE_PROCESS"
             )
row_an = c(rep("Program-1",4),rep("Program-2",4),rep("Program-3",2),rep("Program-4",4),rep("Program-5",4),rep("Program-6",4),rep("Program-7",4))
pathway_prog <- data.frame("pathway" = pathways) 
for (k in 1:7){
  fgseaRes <- read.csv2(paste0("./Results/Z2/GOBP_k",k,".csv"), sep = ',', row.names = 1)
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


library("extrafont")
font_import()

png("./Results/Z2.png", width = 9, height = 7, units = 'in', fonts = "Arial", res = 300)
set.seed(0)
Heatmap(pathway_prog, cluster_rows = F, cluster_columns = F, 
        left_annotation = row_an, name = "NES",
        row_names_max_width =  max_text_width(rownames(pathway_prog), gp = gpar(fontsize = 12) ),
) 

dev.off()