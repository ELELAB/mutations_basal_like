#' plotMoonlightResults 
#' 
#' This function creates a heatmap 
#'
#' @param DEG_Mutations_Annotations A tibble 
#' @param Oncogenic_mediators_mutation_summary A tibble 
#' @param URA A loaded RDA object
#' @param gene_type A character string either "mediators" or "drivers". If NULL defaults to "drivers".
#' @param n The top number of genes plotted. If NULL defaults to 50.
#' @param genelist A vector of strings containing Hugo Symbols of genes.s
#'
#' @import dplyr  
#' @importFrom magrittr "%>%"
#' @import ComlpexHeatmap
#' @import tidyHeatmap
#'
#' @return
#' @export
#'
#' @examples
#' plotMoonlightResults(DEG_Mutations_Annotations, Oncogenic_mediators_mutation_summary,URA, gene_type = "drivers", n = 50)
#' 
plotMoonlight <- function(DEG_Mutations_Annotations, 
                          Oncogenic_mediators_mutation_summary,
                          URA,
                          gene_type = "drivers",
                          n = 50, 
                          genelist = c()){
  
  # The differentially expressed genes, that are annotated as TSG/OCG
  DEGs <- DEG_Mutations_Annotations %>% 
    select(Hugo_Symbol, Moonlight_gene_z_score, logFC, AveExpr, 
           t, P.Value, adj.P.Val, B) %>% 
    unique() %>% 
    drop_na(Moonlight_gene_z_score)
  
  # restructure URA to tibble
  ura <- as_tibble(dataURA, rownames = NA) %>% 
    rownames_to_column(var = "Genes")
  
  ura_wrangled <- ura %>%  
    pivot_longer(cols = !c('Genes'), 
                 names_to = 'Biological_Process',
                 values_to = 'Moonlight_score') %>% 
    right_join(Oncogenic_mediators_mutation_summary, 
               by = c("Genes" = "Hugo_Symbol")) %>% 
    right_join(DEGs, by = c("Genes" = "Hugo_Symbol")) %>%
    replace_na(list(CScape_Driver = 0, 
                    CScape_Passenger = 0, 
                    CScape_Unclassified = 0)) #%>% 
  
  # Type of plot:
  if (gene_type == "mediators"){
    ura_wrangled <- ura_wrangled %>% 
      slice_max(Total_Mutations, n = n, with_ties = FALSE)
    } else if (length(genelist) > 0 ){
    ura_wrangled <- ura_wrangled %>% 
      filter(Genes %in% genelist)
    } else{
    ura_wrangled <- ura_wrangled %>% 
      slice_max(CScape_Driver, n = n, with_ties = FALSE)
    }
  
  # Variable for color scaling in legend
  max_driver <- ura_wrangled %>% arrange(desc(CScape_Driver)) %>% 
    select(CScape_Driver) %>%  head(1) %>% pull
  
  # Plot Heatmap
  bp_heatmap <- heatmap(ura_wrangled,
                        .row = Biological_Process,
                        .column = Genes,
                        .value = Moonlight_score,
                        .scale = "none",
                        clustering_distance_columns = "euclidean",
                        clustering_method_columns = "complete",
                        cluster_rows = FALSE) %>%
    add_tile(Moonlight_Oncogenic_Mediator, palette = c("goldenrod2", "dodgerblue3")) %>%
    add_tile(logFC, palette = c("chartreuse4","firebrick3")) %>%
    add_tile(CScape_Driver, palette = c("dodgerblue1", "steelblue4")) %>% 
    #Use below when tidyHeatmap is updated
    #add_tile(CScape_Driver, palette = colorRamp2(c(0,max_driver), c("white", "dodgerblue3"))) %>%
    add_bar(Total_Mutations) 
  
  save_pdf(bp_heatmap, height = 15, width = 35, units = "cm",
           filename = paste(additionalFilename,"bp_heat_all_drivers.pdf", sep =""))
}