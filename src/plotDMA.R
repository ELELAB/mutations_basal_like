#' plotDMA
#' 
#' This function creates one or more heatmaps on the output from DMA. 
#' It visualises the CScape-Somatic annotations per oncogenic mediator either
#' in a single heatmap or split into several different ones.
#' It is also possible to provide a personalised genelist to visualise.
#'
#' @param DEG_Mutations_Annotations Output file from DMA FORMAT
#' @param Oncogenic_mediators_mutation_summary Output file from DMA FORMAT
#' @param type A character string. It can take the values "split" or "complete". 
#' "split" will split the entire dataset into sections of 40 genes and create individual plots. 
#' These plots will be merged into one pdf. The genes will be sorted alphabeatically.
#' "complete" will create one plot, though it will not be possible to see the individual gene names.
#' The heatmap will be clustered hierarchically. 
#' If both type and genelist is NULL, the function will default to "split".
#'  
#' @param genelist A character vector containing hugo symbols of genes. 
#' A single heatmap will be created with only these genes.
#' The heatmap will be hierarchically clustered. If NULL 
#' Do not provide both type and genelist. 
#' 
#' @param additionalFilename A character string. Adds prefix to the filename of the pdf.
#' 
#' @import dplyr  
#' @importFrom magrittr "%>%"
#' @import ComplexHeatmap
#' @import tidyHeatmap
#' @importFrom qpdf pdf_combine
#' 
#' @return
#' @export
#'
#' @examples
#' plotDMA(DEG_Mutations_Annotions, Oncogenic_mediators_mutation_summary, type = "split", additionalFilename = "myplots_")
#' plotDMA(DEG_Mutations_Annotions, Oncogenic_mediators_mutation_summary, type = "complete", additionalFilename = "myplot_")
#' plotDMA(DEG_Mutations_Annotions, Oncogenic_mediators_mutation_summary, genelist = c("BRCA1", "BRCA2", "FOX1", "GATA3", "TP53"), additionalFilename = "myplot_")

plotDMA <- function(DEG_Mutations_Annotations, 
                    Oncogenic_mediators_mutation_summary,
                    type = "split",
                    genelist = c(),
                    additionalFilename = ""){ 
  # Modify input
  DEGs <- DEG_Mutations_Annotations %>% 
    select(Hugo_Symbol, logFC) %>% 
    distinct()
  
  Summary_wrangled <- Oncogenic_mediators_mutation_summary %>% 
    pivot_longer(cols = c(CScape_Driver, CScape_Passenger, CScape_Unclassified), 
                 names_to = "Mutation_type",
                 values_to = "Count") %>% 
    replace_na(list(Count = 0)) %>% 
    filter(Total_Mutations > 0) %>% 
    left_join(DEGs)
  
  
  if (length(genelist)>0 ){
    
    Summary_wrangled <-  Summary_wrangled %>% 
      filter(Hugo_Symbol %in% genelist) %>% 
      group_by(Moonlight_Oncogenic_Mediator)# %>% 
    
    #Check that there are both OCGs and TGS
    n <- Summary_wrangled %>% summarise() %>% count() %>% pull()
    
    if(n > 1){ 
      driver_mut_heatmap <- Summary_wrangled %>% 
        heatmap(.row = Mutation_type,
                .column = Hugo_Symbol,
                .value = Count,
                .scale = "none",
                cluster_rows = FALSE, cluster_columns = TRUE, show_column_dend = FALSE, 
                show_column_names = TRUE,
                palette_value = c("white", "blue", "darkblue"),
                column_title = paste("Heatmap Driver annotation by CScape-Somatic",
                                     "\n Hugo_Symbol")) %>% 
        add_tile(Moonlight_Oncogenic_Mediator, palette = c("goldenrod2", "dodgerblue3")) %>% 
        add_bar(Total_Mutations) %>% 
        add_tile(logFC, palette = c("chartreuse4","firebrick3"))
      
      #Save plot
      save_pdf(driver_mut_heatmap, height = 15, width = 35, units = "cm",
               filename = paste(additionalFilename, "heatmap_genelist.pdf",
                                sep = "")) 
    }else{
      print("The genelist must contain at least one OCG and one TSG")
    }   
    
  } else if(type == "complete"){

      driver_mut_heatmap <-  Summary_wrangled %>% 
      group_by(Moonlight_Oncogenic_Mediator) %>% 
      heatmap(.row = Mutation_type,
              .column = Hugo_Symbol,
              .value = Count,
              .scale = "none",
              cluster_rows = FALSE, cluster_columns = TRUE, show_column_dend = FALSE, 
              show_column_names = TRUE,
              palette_value = c("white", "blue", "darkblue"),
              column_title = paste("Heatmap Driver annotation by CScape-Somatic",
                                   "\n Hugo_Symbol")) %>% 
      add_tile(Moonlight_Oncogenic_Mediator, palette = c("goldenrod2", "dodgerblue3")) %>% 
      add_bar(Total_Mutations) %>% 
      add_tile(logFC, palette = c("chartreuse4","firebrick3"))
    
    #Save plot
    save_pdf(driver_mut_heatmap, height = 15, width = 35, units = "cm",
             filename = paste(additionalFilename, "heatmap_complete.pdf",
                              sep = ""))
    
    
  }  else if(type == "split"){ 
    # Make/Check heatmap folder!
    
    #Make vector with groups (40 genes in each plot x 3 (driver, pas, unclas) per gene) = 120
    split_vector <- rep(seq(1,ceiling(nrow(Summary_wrangled)/120), by = 1),each=120)
    split_vector <- split_vector[1:nrow(Summary_wrangled)] 
    
    grouped_data <- Summary_wrangled %>% arrange(Hugo_Symbol) %>% 
      mutate(gr = split_vector) %>% 
      group_by(gr) %>%
      tidyr::nest(data = -gr) 
    #%>% head(2)
    
    test<- grouped_data %>% ungroup %>%  
      mutate(First_gene = map(grouped_data$data, ~plotHeatmap(.))) %>% 
      mutate(heatmaps = "heatmaps/heatmap_",
             pdf = ".pdf") %>% 
      unite(heatmaps, c(heatmaps, First_gene, pdf), sep ="")
    
    # Stable pdf into one pdf
    pdfs <- test %>% pull(heatmaps)
    pdf_combine(input = pdfs, 
                output = paste(additionalFilename,"heatmaps_split.pdf", 
                               sep = ""))
    
    #remove heatmap-folder with redudant pdfs
    unlink("heatmaps", recursive = TRUE)
  } 
  
}
