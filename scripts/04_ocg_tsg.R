# Initilization------------------------------------------------------
library(tidyverse)


# Load Data ---------------------------------------------------------
DEGene_All_Mutation_Annotations <- read_csv("results/DEGene_All_Mutation_Annotations.csv")
DEGene_Mutation_Summary <- read_csv("results/DEGene_Mutation_Summary.csv")



# Wrangle data ------------------------------------------------------
# get tumor suppressors and oncogenes 
tsg <- DEGene_All_Mutation_Annotations %>% 
  filter(Driver_type == 'TSG') #%>% select(Hugo_Symbol) %>% distinct()
ocg <- DEGene_All_Mutation_Annotations %>% 
  filter(Driver_type == 'OCG') 


# barplot variant classification count
both_drivers <- DEGene_All_Mutation_Annotations %>% 
  filter(Driver_type == 'OCG' |
           Driver_type == 'TSG') %>% 
  group_by(Driver_type, Variant_Classification) %>% 
  summarise(count = n()) #%>% 
  #mutate(freq = case_when(Driver_type == "OCG" ~ n/272,
  #                      Driver_type == "TSG" ~ n/94))


barplot_tsg_ocg_class <- ggplot(both_drivers, aes(x = Variant_Classification, y = count, 
                 fill = Driver_type), alpha = 0.8)+
  geom_col(position = "dodge2") +
  coord_flip()

ggsave(barplot_tsg_ocg_class, 
       filename = "results/04_barplot_tsg_ocg_classification.png",
       dpi=300)


# Mutations in driver genes -----------------------------------------

# Driver/Passenger scatter plot for TSG/OCGs
deg_moonlight <- DEGene_All_Mutation_Annotations %>% 
  select(Hugo_Symbol, Moonlight_gene_z_score) %>% 
  drop_na(Moonlight_gene_z_score) %>% 
  distinct()

DEGene_Mutation_Summary %>% 
  left_join(deg_moonlight) %>% # Moonlight score for all..
  replace_na(list(Driver = 0, Passenger = 0, Unclassified = 0)) %>% 
  mutate(Freq = Driver/Passenger) %>% #filter(Driver != 0 | Passenger != 0) %>% 
  ggplot(aes(x = Driver, 
             y = Passenger, 
             color = Driver_type), alpha = 0.5) +
  #color = Moonlight_gene_z_score)) +
  #color = Unclassified)) +
  geom_point() +
  geom_jitter() #+
#scale_size(trans = "reverse")



freq_ocg_tsg_dtypes <- DEGene_Mutation_Summary %>% 
  replace_na(list(Driver = 0, Passenger = 0, Unclassified = 0)) %>% 
  pivot_longer(cols = c(Driver, Unclassified, Passenger),
               names_to = "Type",
               values_to = "Count") %>% 
  group_by(Driver_type, Type) %>% summarise(n = sum(Count)) %>% 
  mutate(freq = case_when(Driver_type == 'OCG' ~ n/272,
                          Driver_type == 'TSG'~ n/94)) %>% 
  ggplot(aes(x = Driver_type, y = freq, fill = Type)) +
  geom_col(position = "dodge2") +
  theme(text = element_text(size = 30))+
  labs(x = "Driver Type", y = "Frequency")
freq_ocg_tsg_dtypes

ggsave(freq_ocg_tsg_dtypes,
       filename = "04_freq_drivers_barplot.png", 
       height=6,width=12, units = "in") #, dpi=300)


# Distribtion of cscape-annotation in top 10 mutated TSG/OCGGS ------------
TSG_10 <- DEGene_Mutation_Summary %>%
  filter(Driver_type == "TSG") %>% 
  slice_max(Total_Mutations, n = 10) %>% 
  pivot_longer(cols = c(Driver, Unclassified, Passenger),
               names_to = "Type",
               values_to = "Count") %>% 
  replace_na(list(Count = 0))

tsg_10 <- ggplot(TSG_10, aes(x = Hugo_Symbol, y = Count, fill = Type))+
  geom_col() + #position = "dodge2") +
  labs(title = 'Top 10 mutated TSG') +
  scale_fill_brewer(palette = "Set1")

ggsave(tsg_10,
       filename = "results/04_TSG_mutations.png", 
       dpi = 300)

OCG_10 <- DEGene_Mutation_Summary %>%
  filter(Driver_type == "OCG") %>% 
  slice_max(Total_Mutations, n = 10) %>% 
  pivot_longer(cols = c(Driver, Unclassified, Passenger),
               names_to = "Type",
               values_to = "Count") %>% 
  replace_na(list(Count = 0))

ocg_10 <- ggplot(OCG_10, aes(x = Hugo_Symbol, y = Count, fill = Type))+
  geom_col() + #position = "dodge2") +
  labs(title = 'Top 10 mutated OCG') +
  scale_fill_brewer(palette = "Set1")

ggsave(ocg_10,
       filename = "results/04_OCG_mutations.png", 
       dpi = 300)


# Distribution of Driver mutations in variant Classification ------------
#Count:
driver_genes_class_number <- DEGene_All_Mutation_Annotations %>% 
  filter(!is.na(Driver_type), 
         Driver_Mutation == 'Driver') %>% 
  group_by(Driver_type, Variant_Classification) %>% count() %>% 
  mutate(freq = case_when(Driver_type == 'OCG' ~ n/272,
                          Driver_type == 'TSG'~ n/94)) %>% 
  ggplot(aes(x = Driver_type, y = n, fill = Variant_Classification))+
  geom_col(position = "dodge2") +
  theme(text = element_text(size = 30))
ggsave(driver_genes_class_number,
       filename = "results/04_driver_gene_class_num.png", 
       dpi=300)

#Frequency: 
driver_genes_class_freq <- DEGene_All_Mutation_Annotations %>% 
  filter(!is.na(Driver_type),
         Driver_Mutation == 'Driver') %>% 
  group_by(Driver_type, Variant_Classification) %>% count() %>% 
  mutate(freq = case_when(Driver_type == 'OCG' ~ n/272,
                          Driver_type == 'TSG'~ n/94)) %>% 
  ggplot(aes(x = Driver_type, y = freq, fill = Variant_Classification))+
  geom_col(position = "dodge2") +
  theme(text = element_text(size = 30)) +
  labs(x= "Driver Type", y = "Frequency")

ggsave(driver_genes_class_freq,
       filename = "results/04_driver_gene_class_freq.png", 
       dpi=300)


#Up and down regulation in TSG and OCGs 
regulation_tsg_ocg <- DEGene_All_Mutation_Annotations %>% 
  filter(!is.na(Driver_type)) %>% 
  group_by(Hugo_Symbol, logFC, Driver_type) %>% count() %>% 
  ggplot(aes(x=logFC, y = n, color = Driver_type),alpha = 0.5) +
  geom_point() +
  theme(text = element_text(size = 20))

ggsave(regulation_tsg_ocg,
       filename = "results/04_regulation_tsg_ocg.png", 
       height=6, width=12, units = "in") #, dpi=300)
