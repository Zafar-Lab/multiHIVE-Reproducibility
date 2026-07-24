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

gene_spectrum_scores_path = "../cNMF/Thymus-Z2/bc_cNMF/bc_cNMF.gene_spectra_score.k_8.dt_0_02.txt"
# gene_mapping_path = "../cNMF/Thymus-Z1/gene-mapping.csv"


results_path = "./Results/Thymus-Z2/"


gene_spectrum_scores = read.csv2(gene_spectrum_scores_path, sep="\t",row.names = 1, header = T)
gene_spectrum_scores = as.data.frame(sapply(gene_spectrum_scores,as.numeric))

# gene_mapping_t = read.csv2(gene_mapping_path, sep = ',', header=F)
# gene_mapping = gene_mapping_t$V2
# names(gene_mapping) = gene_mapping_t$V1

# colnames(gene_spectrum_scores) = sapply(colnames(gene_spectrum_scores), function(x){gene_mapping[as.character(substr(x,2,nchar(x)))]})

for (i in 1:dim(gene_spectrum_scores)[1]){
  
  scores <- as.data.frame(t(gene_spectrum_scores[i,]))
  colnames(scores)= c('scores')
  scores = scores %>% arrange(desc(scores))
  
  ranks = as.numeric(scores$scores)
  names(ranks) = rownames(scores)
  
  #hallmark
  
  pathways = gmtPathways("../../GSEA/databases/mouse/mh.all.v2023.2.Mm.symbols.gmt")
  fgseaRes <- fgsea(pathways = pathways,
                    stats    = ranks,
                    eps      = 0.0, # get p values more accurate
                    minSize  = 15,
                    maxSize  = 500)
  fgseaRes <- fgseaRes %>% arrange(NES) %>% filter(padj<0.05) %>% arrange(NES)
  write.csv(apply(fgseaRes, 2, as.character), file = paste0(results_path,"hallmark_k",i,".csv"))
  gseaenrichment_plot(fgseaRes, paste0(results_path,"hallmark_k", i, ".png"))
  
  #all pathways
  pathways <- gmtPathways("../../GSEA/databases/mouse/msigdb.v2023.2.Mm.symbols.gmt")
  
  fgseaRes <- fgsea(pathways = pathways,
                    stats    = ranks,
                    eps      = 0.0, # get p values more accurate
                    minSize  = 15,
                    maxSize  = 500)
  fgseaRes <- fgseaRes %>% arrange(NES) %>% filter(padj<0.05) %>% arrange(NES)
  write.csv(apply(fgseaRes, 2, as.character), file = paste0(results_path,"all_k",i,".csv"))
  gseaenrichment_plot(fgseaRes, paste0(results_path,"all_k", i, ".png"))
  
  #gobp pathways
  pathways = gmtPathways("../../GSEA/databases/mouse/m5.go.bp.v2023.2.Mm.symbols.gmt")
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
fgseaRes <- read.csv2("./Results/Thymus-Z2/GOBP_k7.csv", sep = ',', row.names = 1)
res_pos <- fgseaRes[fgseaRes$NES>0, ] %>% arrange(desc(NES))
res_pos[grep("dna", res_pos$pathway, ignore.case = T),c("pathway", "padj", "NES")]
res_pos$pathway[grep("natur", res_pos$pathway, ignore.case = T)]
#2 "GOBP_T_CELL_DIFFERENTIATION", 
# 3 "GOBP_RESPONSE_TO_TYPE_II_INTERFERON", "GOBP_RESPONSE_TO_TYPE_I_INTERFERON", "GOBP_INTERFERON_MEDIATED_SIGNALING_PATHWAY",
# 4  "GOBP_MITOCHONDRIAL_TRANSLATION", "GOBP_RNA_PROCESSING",  "GOBP_PROTEIN_FOLDING",  "GOBP_PROTEIN_STABILIZATION", "GOBP_RNA_METHYLATION",
# 5 "GOBP_MITOTIC_CELL_CYCLE_PROCESS",  "GOBP_POSITIVE_REGULATION_OF_TRANSCRIPTION_BY_RNA_POLYMERASE_I", "GOBP_REGULATION_OF_CELL_CYCLE", "GOBP_TRANSCRIPTION_BY_RNA_POLYMERASE_I", "GOBP_NATURAL_KILLER_CELL_ACTIVATION",
# 6 "GOBP_POSITIVE_REGULATION_OF_NATURAL_KILLER_CELL_MEDIATED_IMMUNITY", "GOBP_REGULATION_OF_NATURAL_KILLER_CELL_MEDIATED_IMMUNITY", "GOBP_POSITIVE_REGULATION_OF_CELL_KILLING",
# 7 "GOBP_DNA_RECOMBINATION",
#8 "GOBP_CHROMOSOME_SEPARATION","GOBP_CHROMOSOME_SEGREGATION", "GOBP_MITOTIC_SISTER_CHROMATID_SEGREGATION",  "GOBP_NUCLEAR_CHROMOSOME_SEGREGATION", "GOBP_SISTER_CHROMATID_SEGREGATION", 
pathways = c( "GOBP_ALDEHYDE_BIOSYNTHETIC_PROCESS" ,  "GOBP_CELLULAR_ALDEHYDE_METABOLIC_PROCESS",  "GOBP_ATP_METABOLIC_PROCESS", #"GOBP_CYTOPLASMIC_TRANSLATION", "GOBP_TRANSLATION_AT_SYNAPSE",
              "GOBP_MONONUCLEAR_CELL_DIFFERENTIATION", "GOBP_T_CELL_ACTIVATION", 
              "GOBP_RESPONSE_TO_INTERFERON_BETA",   "GOBP_TYPE_I_INTERFERON_PRODUCTION", "GOBP_TYPE_II_INTERFERON_PRODUCTION",
              "GOBP_MITOCHONDRIAL_GENE_EXPRESSION", "GOBP_RRNA_METABOLIC_PROCESS",  "GOBP_TRNA_PROCESSING",  "GOBP_PROTEIN_MATURATION",
               "GOBP_MITOTIC_CELL_CYCLE", "GOBP_CELL_DIVISION", 
                "GOBP_NATURAL_KILLER_CELL_MEDIATED_IMMUNITY",  "GOBP_CELL_KILLING", 
              "GOBP_DNA_REPLICATION", "GOBP_DNA_REPAIR",  "GOBP_PROTEIN_DNA_COMPLEX_ASSEMBLY",
                "GOBP_MITOTIC_NUCLEAR_DIVISION" , "GOBP_CHROMOSOME_LOCALIZATION"

)
row_an = c(rep("Program-1",3),rep("Program-2",2),rep("Program-3",3),rep("Program-4",4),rep("Program-5",2),rep("Program-6",2),rep("Program-7",3),rep("Program-8",2))
pathway_prog <- data.frame("pathway" = pathways) 
for (k in 1:8){
  fgseaRes <- read.csv2(paste0("./Results/Thymus-Z2/GOBP_k",k,".csv"), sep = ',', row.names = 1)
  col_name = paste0("Program-",k)
  fgseaRes = fgseaRes[fgseaRes$pathway %in% pathways, c("pathway","NES")] %>% rename(!!col_name:="NES")
  fgseaRes[,col_name] = as.numeric(fgseaRes[,col_name])
  pathway_prog = merge(pathway_prog,fgseaRes, by="pathway", all.x = T)
}
pathway_prog = pathway_prog %>% distinct(pathway, .keep_all = T) %>% remove_rownames %>% column_to_rownames("pathway") %>% as.matrix()
pathway_prog = pathway_prog[pathways,]
# pathway_prog[is.na(pathway_prog)] = NA
row_an = rowAnnotation(Program = as.factor(row_an))
rownames(pathway_prog) = sapply(rownames(pathway_prog), function (x){str_replace(x,"GOBP_","")})
png("./Results/Z2.png", width = 7, height = 5.5, fonts = "Arial", units = 'in', res = 300)
set.seed(0)
Heatmap(pathway_prog, cluster_rows = F, cluster_columns = F, 
        left_annotation = row_an, name = "NES",
        row_names_max_width =  max_text_width(rownames(pathway_prog), gp = gpar(fontsize = 12) )
) 
dev.off()
#####################################################################################

