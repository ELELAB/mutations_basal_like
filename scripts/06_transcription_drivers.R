# Initilization------------------------------------------------------
library(tidyverse)


# Load Data ---------------------------------------------------------
DEGene_All_Mutation_Annotations <- read_csv("results/DEGene_All_Mutation_Annotations.csv")
DEGene_Mutation_Summary <- read_csv("results/DEGene_Mutation_Summary.csv")
MafFile = read_tsv("data/mutations.tsv")



# Wrangle data ------------------------------------------------------

# Are the driver mutations in the regulatory region - transcriptional level
transcription_drivers <- DEGene_All_Mutation_Annotations %>% filter(!is.na(Driver_type), 
                                           Driver_Mutation  == "Driver",
                                           Potential_Effect_on_Transcription == 1) %>%  
  group_by(Hugo_Symbol, Driver_type, logFC,
           Variant_Classification, Variant_Type, 
           One_Consequence, Start_Position, Annotation,
           COSMIC) %>% 
  summarise(n = n()) 

write_csv(transcription_drivers, file = "results/06_transcription_drivers.csv")


# How many can/cannot be caused by mutations: 
DEGene_Mutation_Summary %>% filter(Passenger == Total_Mutations) %>% group_by(Driver_type) %>% summarise(n = n())
DEGene_Mutation_Summary %>% filter(Driver == Total_Mutations) %>% group_by(Driver_type) %>% summarise(n = n())


