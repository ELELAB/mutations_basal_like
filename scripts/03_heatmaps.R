# Initilization------------------------------------------------------
library(tidyverse)
library(ComplexHeatmap)
library(tidyHeatmap)

# Load Data ---------------------------------------------------------
DEGene_All_Mutation_Annotations <- read_csv("results/DEGene_All_Mutation_Annotations.csv")
DEGene_Mutation_Summary <- read_csv("results/DEGene_Mutation_Summary.csv")
basal_spsecific_ura <- readRDS("data/basal_specific_ura.RDS")
ura <- as_tibble(basal_specific_ura, rownames = NA) %>% 
  rownames_to_column(var = "Genes")

# Wrangle data ------------------------------------------------------

# Heatmap1: Driver mutations CScape-classification heatmap:
DEGs <- DEGene_All_Mutation_Annotations %>% select(Hugo_Symbol, logFC) %>% distinct()
DEGene_Mutation_Summary_wrangled <- DEGene_Mutation_Summary %>% 
  pivot_longer(cols = c(Driver, Passenger, Unclassified), 
               names_to = "Mutation_type",
               values_to = "Count") %>% 
  replace_na(list(Count = 0)) %>% 
  filter(Total_Mutations > 0) %>% #, Driver_type == "OCG")
  left_join(DEGs)

driver_mut_heatmap <- DEGene_Mutation_Summary_wrangled %>%  group_by(Driver_type) %>% 
  heatmap(.row = Mutation_type,
          .column = Hugo_Symbol,
          .value = Count,
          cluster_rows = FALSE, cluster_columns = TRUE, show_column_dend = FALSE, 
          column_names_gp = gpar(fontsize = 3), column_names_rot = 45,  
          palette_value = c("white", "blue", "darkblue")) %>% 
  #heatmap_legend_param = list(at = c(0,1,6))) %>% 
  add_tile(Driver_type, palette = c("goldenrod2", "dodgerblue3")) %>% 
  add_bar(Total_Mutations) %>% 
  add_tile(logFC, palette = c("chartreuse4","firebrick3"))

save_pdf(driver_mut_heatmap, height = 15, width = 35, units = "cm",
         filename = "results/03_driver_mut_heatmap.pdf")


# Heatmap 2: Biological process vs driver genes ------------------------------

# The differentially expressed genes, that are annotated as TSG/OCG
DEGs <- DEGene_All_Mutation_Annotations %>% 
  select(Hugo_Symbol, Moonlight_gene_z_score, logFC, AveExpr, t, P.Value, adj.P.Val, B) %>% 
  unique() %>% 
  drop_na(Moonlight_gene_z_score)

ura_wrangled <- ura %>%  
  pivot_longer(cols = !c('Genes'), 
               names_to = 'Biological_Process',
               values_to = 'Moonlight_score') %>% 
  right_join(DEGene_Mutation_Summary, by = c("Genes" = "Hugo_Symbol")) %>% 
  right_join(DEGs, by = c("Genes" = "Hugo_Symbol")) %>%
  slice_max(Moonlight_gene_z_score, n = 60) %>% 
  replace_na(list(Driver = 0, Passenger = 0, Unclassified = 0))


bp_heatmap <- heatmap(ura_wrangled,
        .row = Biological_Process,
        .column = Genes,
        .value = Moonlight_score,
        clustering_distance_columns = "euclidean",
        clustering_method_columns = "complete",
        cluster_rows = FALSE) %>%
  add_tile(Driver_type, palette = c("goldenrod2", "dodgerblue3")) %>%
  add_tile(logFC, palette = c("chartreuse4","firebrick3")) %>%
  #add_tile(Moonlight_gene_z_score, palette = c("chartreuse4","firebrick3")) %>%
  #add_tile(Passenger) %>%
  #add_tile(Unclassified) %>%
  add_tile(Driver, palette = c("white", "dodgerblue3")) %>%
  add_bar(Total_Mutations) 

save_pdf(bp_heatmap, height = 15, width = 35, units = "cm",
         filename = "results/03_bp_driver_genes_heatmap.pdf")
