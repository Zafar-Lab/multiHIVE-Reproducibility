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

gene_spectrum_scores_path = "../cNMF/Thymus-Z1/cNMF/bc_cNMF/bc_cNMF.gene_spectra_score.k_11.dt_0_02.txt"
gene_mapping_path = "../cNMF/Thymus-Z1/gene-mapping.csv"


results_path = "./Results/Thymus-Z1/11/"


gene_spectrum_scores = read.csv2(gene_spectrum_scores_path, sep="\t",row.names = 1, header = T)
gene_spectrum_scores = as.data.frame(sapply(gene_spectrum_scores,as.numeric))

gene_mapping_t = read.csv2(gene_mapping_path, sep = ',', header=F)
gene_mapping = gene_mapping_t$V2
names(gene_mapping) = gene_mapping_t$V1

colnames(gene_spectrum_scores) = sapply(colnames(gene_spectrum_scores), function(x){gene_mapping[as.character(substr(x,2,nchar(x)))]})

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
  
  
  #reactomes pathways
  pathways = gmtPathways("../../misc/GSEA/databases/mouse/m2.cp.v2024.1.Mm.symbols.gmt")
  fgseaRes <- fgsea(pathways = pathways, 
                    stats    = ranks,
                    eps      = 0.0, # get p values more accurate
                    minSize  = 15,
                    maxSize  = 500)
  fgseaRes <- fgseaRes %>% arrange(NES) %>% filter(padj<0.05) %>% arrange(NES)
  write.csv(apply(fgseaRes, 2, as.character), file = paste0(results_path,"CP_k",i,".csv"))
  gseaenrichment_plot(fgseaRes, paste0(results_path,"CP_k", i, ".png"))
  
}

#####################################################################################
fgseaRes <- read.csv2("./Results/Thymus-Z1/11/GOBP_k4.csv", sep = ',', row.names = 1)
res_pos <- fgseaRes[fgseaRes$NES>0, ] %>% arrange(desc(NES))
res_pos[grep("GOBP_NCRNA_PROCESSING", res_pos$pathway, ignore.case = T),c("pathway", "padj", "NES")]
res_pos$pathway[grep("myeloid", res_pos$pathway, ignore.case = T)]

# 2 "GOBP_REGULATION_OF_B_CELL_ACTIVATION","REACTOME_FCERI_MEDIATED_NF_KB_ACTIVATION", "REACTOME_CELLULAR_RESPONSE_TO_HYPOXIA", "HALLMARK_MYC_TARGETS_V1","HALLMARK_MYC_TARGETS_V2","HALLMARK_MTORC1_SIGNALING", 
#3 "GOBP_REGULATION_OF_T_CELL_CYTOKINE_PRODUCTION","GOBP_NATURAL_KILLER_CELL_CYTOKINE_PRODUCTION", "GOBP_CYTOKINE_PRODUCTION", "GOBP_CYTOKINE_PRODUCTION_INVOLVED_IN_IMMUNE_RESPONSE", "GOBP_ANTIBACTERIAL_HUMORAL_RESPONSE", "GOBP_DEFENSE_RESPONSE_TO_BACTERIUM", "GOBP_NEGATIVE_REGULATION_OF_IMMUNE_RESPONSE", "GOBP_RESPONSE_TO_VIRUS", "GOBP_POSITIVE_REGULATION_OF_DEFENSE_RESPONSE", "GOBP_REGULATION_OF_INNATE_IMMUNE_RESPONSE","HALLMARK_INTERFERON_GAMMA_RESPONSE","HALLMARK_INFLAMMATORY_RESPONSE",
# 4 "GOBP_RESPONSE_TO_CYTOKINE", "HALLMARK_MTORC1_SIGNALING", "GOBP_RRNA_METABOLIC_PROCESS", "GOBP_NCRNA_PROCESSING", "GOBP_TRNA_METABOLIC_PROCESS", "REACTOME_DECTIN_1_MEDIATED_NONCANONICAL_NF_KB_SIGNALING", 
#6 "GOBP_DNA_REPLICATION", "GOBP_POSITIVE_REGULATION_OF_CELL_CYCLE_PROCESS", "HALLMARK_E2F_TARGETS", 
#10 "GOBP_ACTIVATION_OF_IMMUNE_RESPONSE", "GOBP_ANTIGEN_PROCESSING_AND_PRESENTATION",
#11 "GOBP_NEUTROPHIL_MIGRATION",
# 1,2,3,4,5,6,7,8,9,10,11
pathways = c("REACTOME_GENERATION_OF_SECOND_MESSENGER_MOLECULES","REACTOME_TRANSLOCATION_OF_ZAP_70_TO_IMMUNOLOGICAL_SYNAPSE","REACTOME_PHOSPHORYLATION_OF_CD3_AND_TCR_ZETA_CHAINS",
             "GOBP_THYMIC_T_CELL_SELECTION","GOBP_LYMPHOCYTE_COSTIMULATION","HALLMARK_TNFA_SIGNALING_VIA_NFKB","HALLMARK_HYPOXIA",
             "GOBP_CD8_POSITIVE_ALPHA_BETA_T_CELL_ACTIVATION","GOBP_T_CELL_MEDIATED_CYTOTOXICITY","GOBP_NATURAL_KILLER_CELL_MEDIATED_IMMUNITY",  "GOBP_REGULATION_OF_ACUTE_INFLAMMATORY_RESPONSE", "HALLMARK_IL2_STAT5_SIGNALING",
             "REACTOME_ANTIGEN_PROCESSING_CROSS_PRESENTATION","GOBP_RNA_PROCESSING", "GOBP_TRNA_PROCESSING",
              "GOBP_DNA_REPAIR" , 
              "GOBP_CELL_CYCLE_DNA_REPLICATION","GOBP_DNA_DAMAGE_RESPONSE",
              "GOBP_CELL_DIVISION", "GOBP_REGULATION_OF_CELL_DIVISION",
  "GOBP_TRANSMEMBRANE_RECEPTOR_PROTEIN_TYROSINE_KINASE_SIGNALING_PATHWAY",
                 "GOBP_MYELOID_CELL_ACTIVATION_INVOLVED_IN_IMMUNE_RESPONSE", "GOBP_REGULATION_OF_LEUKOCYTE_MEDIATED_IMMUNITY",
                "GOBP_MAST_CELL_ACTIVATION_INVOLVED_IN_IMMUNE_RESPONSE",  "GOBP_BLASTOCYST_GROWTH",
              "GOBP_NEUTROPHIL_ACTIVATION",  "GOBP_NEUTROPHIL_MEDIATED_IMMUNITY"
)
row_an = c(rep("Program-1",3),rep("Program-2",4),rep("Program-3",5), rep("Program-4",3),rep("Program-5",1),rep("Program-6",2),rep("Program-7",2),rep("Program-8",1),rep("Program-9",2),rep("Program-10",2),rep("Program-11",2))
row_an = factor(row_an,levels = c("Program-1", "Program-2", "Program-3", "Program-4", "Program-5", "Program-6", "Program-7", "Program-8", "Program-9", "Program-10", "Program-11"))
pathway_prog <- data.frame("pathway" =pathways) 
for (k in c(1:11)){
  fgseaRes_all <- read.csv2(paste0("./Results/Thymus-Z1/11/all_k",k,".csv"), sep = ',', row.names = 1)
  fgseaRes1 <- read.csv2(paste0("./Results/Thymus-Z1/11/hallmark_k",k,".csv"), sep = ',', row.names = 1)
  fgseaRes2 <- read.csv2(paste0("./Results/Thymus-Z1/11/GOBP_k",k,".csv"), sep = ',', row.names = 1)
  fgseaRes3 <- read.csv2(paste0("./Results/Thymus-Z1/11/CP_k",k,".csv"), sep = ',', row.names = 1)
  fgseaRes <- rbind(fgseaRes1,fgseaRes2,fgseaRes3)
  fgseaRes_all <- fgseaRes_all %>% filter(!(pathway %in% fgseaRes$pathway))
  fgseaRes <- rbind(fgseaRes,fgseaRes_all)
  
  col_name = paste0("Program-",k)
  fgseaRes = fgseaRes[fgseaRes$pathway %in% pathways, c("pathway","NES")] %>% rename(!!col_name:="NES")
  fgseaRes[,col_name] = as.numeric(fgseaRes[,col_name])
  pathway_prog = merge(pathway_prog,fgseaRes, by="pathway", all.x = T)
}

pathway_prog = pathway_prog %>% distinct(pathway, .keep_all = T) %>% remove_rownames %>% column_to_rownames("pathway") %>% as.matrix()
pathway_prog = pathway_prog[c(pathways),]
# pathway_prog[is.na(pathway_prog)] = NA
row_an = rowAnnotation(Program = as.factor(row_an))
rownames(pathway_prog) = sapply(rownames(pathway_prog), function (x){str_replace(x,"GOBP_","")})
rownames(pathway_prog) = sapply(rownames(pathway_prog), function (x){str_replace(x,"REACTOME_","")})
rownames(pathway_prog) = sapply(rownames(pathway_prog), function (x){str_replace(x,"HALLMARK_","")})

png("./Results/Z1.png", width = 10.5, height = 7, fonts = "Arial", units = 'in', res = 300)
set.seed(0)
Heatmap(pathway_prog, cluster_rows = F, cluster_columns = F, 
        left_annotation = row_an, name = "NES",
        row_names_max_width =  max_text_width(rownames(pathway_prog), gp = gpar(fontsize = 12) )
) 
dev.off()