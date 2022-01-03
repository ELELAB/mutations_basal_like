#' DMA
#'
#' This function carries out the driver mutation analysis
#' @param MafFile Maf file containing mutations
#' @param DEGs Differentially expressed genes RDS object
#' @param Drivers driver genes obtained from PRA
#' @param coding_file 
#' @param noncoding_file
#' @param results_folder
#' @importFrom dplyr  
#' @importFrom fuzzyjoin genome_left_join
#' @importFrom magrittr "%>%"
#' @importFrom readr read_tsv read_csv write_csv
#' @import maftools read.maf plotmafSummary rainfallPlot somaticInteractions
#' @import MAFtoCscape PRAtoTibble
#' @import LiftMAF RunCscape_somatic
#'
#'
#' @return no return files and plots are saved
#' @export
#' @examples
#' 
#' DMA(MafFile, DEGs, Drivers, coding_file, noncoding_file, results_folder)

MutationAnalysis <- function(MafFile, DEGs, Drivers, coding_file, noncoding_file, cosmic_file, results_folder){
  
  # Create Output folder
  if (dir.exists(results_folder)){
    print("Output folder already exits")
  }
  else {
    dir.create(path = results_folder, showWarnings = TRUE, recursive = TRUE)
  }
  
  # Load Data --------------------------------
  # read maf and add ID number to each mutation
  mutations <- read_tsv(MafFile) %>% #Save maf as Tibble
    mutate(ID = row_number()) %>% 
    relocate(ID)  
  
  drivers_moonlight <- PRAtoTibble(Drivers)
  DEGs <- DEGs %>% rownames_to_column(var = 'Hugo_Symbol')
  
  # Load homemade mutations effect on transcription table
  transcription_binary <- read_tsv("data/transcription_mutations.tsv") %>% 
    pivot_longer(cols = !Variant_Classification, 
                 names_to = "Variant_Type", 
                 values_to = "Potential_Effect_on_Transcription")
  
  translation_binary <- read_tsv("data/translation_mutations.tsv") %>% 
    pivot_longer(cols = !Variant_Classification, 
                 names_to = "Variant_Type", 
                 values_to = "Potential_Effect_on_Translation")
  
  protein_binary <- read_tsv("data/protein_structure_mutations.tsv") %>% 
    pivot_longer(cols = !Variant_Classification, 
                 names_to = "Variant_Type", 
                 values_to = "Potential_Effect_on_Protein")
  
  # Load ENCODE files (promoter and enhancer)
  promoters <- read_tsv("data/ENCFF140XLU.bed.gz", col_names = FALSE) %>% 
    mutate(Annotation = 'Promoter')
  enhancers <- read_tsv("data/ENCFF212UAV.bed.gz", col_names = FALSE) %>% 
    mutate(Annotation = 'Enhancer')
  
  # Wrangle Data -----------------------------
  # Keep only mutations in DEGs 
  DEGs_mut <- DEGs %>% left_join(mutations, by = 'Hugo_Symbol') 
  
  
  # Run Cscape ---------------------
  #Lifting from one build to another
  mut_only <- DEGs_mut %>% filter(!is.na(ID)) 
  DEGs_mut_hg19 <- LiftMAF(Infile = mut_only, Current_Build = 'GRCh38')
  
  cscape_in <- MAFtoCscape(DEGs_mut_hg19) 
  
  print("Cscape will now run. It takes some time")
  cscape_out <- RunCscape_somatic(input = cscape_in,
                                  coding_file = coding_file,
                                  noncoding_file = noncoding_file)
  
  # merge cscape results
  print('Cscape finished. Output file is saved in result folder.')
  cscape_out <- cscape_out %>% 
    mutate(Variant_Type = "SNP")
  
  DEGs_mut_annotated_19 <- DEGs_mut_hg19 %>% 
    separate(Chromosome, into = c(NA, "Chr"), sep = 3, remove = FALSE)%>% 
    left_join(cscape_out, 
              by = c("Start_Position" = "Position", 
                     "Variant_Type", 
                     "Chr")) %>%
    mutate(Driver_Mutation = case_when((Coding_score > 0.5 | Noncoding_score > 0.5) ~ "Driver",    #Driver
                                       (Coding_score <= 0.5 | Noncoding_score <= 0.5 ~ "Passenger"), #Passenger
                                       TRUE ~ "Unclassified")) #When no score is found
  
  #Lift back to 38
  DEGs_mut_annotated <- LiftMAF(Infile = DEGs_mut_annotated_19, 
                                Current_Build = "GRCh37")
  #(Add here for SNPs and DELETETIONS Rating, with other tool)
  
  
  # Annotate file --------------------------------
  
  # Add level of consequence ----------
  DEGs_mut_annotated <- DEGs_mut_annotated %>% 
    left_join(transcription_binary,
              by = c("Variant_Classification",
                     "Variant_Type")) %>% 
    left_join(translation_binary,
              by = c("Variant_Classification",
                     "Variant_Type")) %>% 
    left_join(protein_binary,
              by = c("Variant_Classification",
                     "Variant_Type")) 
  
  # Add again the DEGs with no mutations --
    # and merge driver status from moonlight
  no_muts <- DEGs_mut %>% filter(is.na(ID))
  DEGs_mut_annotated <- bind_rows(DEGs_mut_annotated, no_muts) %>% 
    replace_na(list(Driver_Mutation = 'No_mutations'))
  
  # Add drivers from moonlight ------------
  DEGs_mut_annotated <- DEGs_mut_annotated %>% 
    left_join(drivers_moonlight,
              by = "Hugo_Symbol")
  
  
  ## Promoters and enhancers -------------
  
  #Merge Encode files into one
  encode_tibble <- full_join(promoters, enhancers) %>% 
    dplyr::select(X1, X2, X3, Annotation) %>% 
    dplyr::rename("Chromosome_annot" = "X1",
                  "Annotation_Start" = "X2",
                  "Annotation_End" = "X3")
  
  # Join when mutation overlap with an annotations position
  mutations_encode <- genome_left_join(x = DEGs_mut_annotated,
                                       y = encode_tibble,
                                       by = c("Chromosome" = "Chromosome_annot",
                                              "Start_Position" = "Annotation_Start",
                                              "End_Position" = "Annotation_End")) %>% 
    dplyr::select(-c("Chromosome_annot"))
  
  # If a mutation overlap with several annotation, keep them in the same row
  DEGs_mut_annotated <- mutations_encode %>% group_by(ID) %>%  
    mutate(Annotation = toString(Annotation),
           Annotation_Start = toString(Annotation_Start),
           Annotation_End = toString(Annotation_End)) %>% 
    ungroup() %>% unique() 
  
  
  # Transcription factors -----------------------
  # Not yet implemented
  
  
  # Make Summary Table --------------------------
  Summary_per_gene_1 <- DEGs_mut_annotated %>% 
    group_by(Hugo_Symbol, Driver_type, Driver_Mutation) %>% 
    summarise(n = n()) %>% 
    pivot_wider(names_from = Driver_Mutation,
                values_from = n) 
  
  Summary_per_gene_2 <- DEGs_mut_annotated %>% 
    group_by(Hugo_Symbol, Driver_type) %>% 
    summarise(Transcription_mut_sum = sum(Potential_Effect_on_Transcription, na.rm = TRUE),
              Translation_mut_sum = sum(Potential_Effect_on_Translation, na.rm = TRUE),
              Protein_mut_sum = sum(Potential_Effect_on_Protein, na.rm = TRUE),
              Total_Mutations = sum(!is.na(ID))) 
  
  Summary_per_gene <- full_join(Summary_per_gene_1, Summary_per_gene_2) %>% 
    filter(!is.na(Driver_type)) %>% 
    arrange(desc(Driver, Total_Mutations)) %>% 
    dplyr::select(!No_mutations)
  
  write_csv(x = Summary_per_gene,
            path = paste(results_folder,"DEGene_Mutation_Summary.csv", sep ='/'),
            col_names = TRUE)
  
  
  # Make 'Raw' Table ---------------------------
  #This table is just a cleaned-up version of the annotated table
  #(when other tools to estimate score are implemented, pivot cscape scores and other score)
  DEGs_mut_Raw_out <- DEGs_mut_annotated %>% 
    relocate(any_of(c("Moonlight_gene_z_score", 
                      "Driver_type",
                      "Coding", "Noncoding", 
                      "Driver_Mutation", 
                      "Potential_Effect_on_Transcription",
                      "Potential_Effect_on_Translation", 
                      "Potential_Effect_on_Protein")), .after = "B") %>% 
    mutate(Driver_Mutation = replace_na(Driver_Mutation, 'Unclassified'))
  
  write_csv(x = DEGs_mut_Raw_out,
            path = paste(results_folder,"/DEGene_All_Mutation_Annotations.csv", sep = ''),
            col_names = TRUE)
  
  
  # Maftools functions --------------
  MafFile <- read.maf(MafFile)
  png(filename = paste(results_folder, "/PlotSumMAF.png", sep = ""))
  plotmafSummary(MafFile)
  dev.off()
  
  png(filename = paste(results_folder, "/PlotRainfallMAF.png", sep = ""))
  rainfallPlot(MafFile)
  dev.off()
  
  png(filename = paste(results_folder, "/PlotSomaticInteractionsMAF.png", sep = ""))
  somaticInteractions(MafFile)
  dev.off()
  
} # End of function --------------------------------------------
