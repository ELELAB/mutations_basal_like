# Initilization------------------------------------------------------
library(tidyverse)

# Load Data ---------------------------------------------------------
DEGene_All_Mutation_Annotations <- read_csv("results/DEGene_All_Mutation_Annotations.csv")


# Wrangle data ------------------------------------------------------

# Distribution of the Promoter/Enhancer annotation ----
Deg_dist_only_class <- DEGene_All_Mutation_Annotations %>% 
  group_by(Variant_Classification) %>% 
  summarise(Class_count = n())


# Plot both annoations 
clean_annotation <- DEGene_All_Mutation_Annotations %>% filter(!is.na(Annotation),
                                                !is.na(Variant_Classification),
                                                !is.na(Variant_Type),
                                                Annotation != 'NA, NA') %>% 
  dplyr::select(Variant_Classification, Variant_Type, Annotation)

both_annotations <- ggplot(clean_annotation, aes(x = Annotation, fill = Variant_Classification), alpha= 0.7) +
  geom_bar(position = "dodge2") +
  theme(text =element_text(size =20))
#labs(title = "Distribution of the enhancers and promoters position in genomic regions")
ggsave(plot = both_annotations, 
       filename = "results/05_promoter_enhancer_count.png")



# Promoters --------------------------------------------------------
promoter_sum <- DEGene_All_Mutation_Annotations %>% 
  filter(Annotation == "Promoter") %>% 
  dplyr::select(Variant_Classification, Variant_Type) %>% 
  group_by(Variant_Classification) %>% 
  summarise(Pro_count = n()) %>% 
  left_join(Deg_dist_only_class) %>% 
  mutate(Freq = Pro_count/Class_count)

promoter_plot <- ggplot(promoter_sum, aes(x = Variant_Classification, 
                                          y = Freq,
                                          fill = Variant_Classification), alpha = 0.7)+
  geom_col(show.legend = FALSE)+
  coord_flip() +
  geom_label(aes(label = Pro_count), nudge_y= 0.07, show.legend = FALSE, size=6)+
  theme(text = element_text(size = 20))
#labs(title = "Promoters", 
#     caption = "Normalised acording to number of observation \n in each classifiaction")
ggsave(plot = promoter_plot, 
       filename = "results/05_promoter_freq.png")




# Enhancers -------------------------------------------------------
enhancer_sum <- DEGene_All_Mutation_Annotations %>% 
  filter(Annotation == "Enhancer") %>% 
  dplyr::select(Variant_Classification, Variant_Type) %>% 
  group_by(Variant_Classification) %>% 
  summarise(Enh_count = n()) %>%
  left_join(Deg_dist_only_class) %>% 
  mutate(Freq = Enh_count/Class_count) %>% 
  add_row(Variant_Classification = "Transcriptional_start_site", Enh_count = 0,
          Class_count = 0,Freq = 0)


enhancer_plot <- ggplot(enhancer_sum, aes(x = Variant_Classification, 
                                          y = Freq, 
                                          fill = Variant_Classification), alpha = 0.7)+
  geom_col(position = position_dodge(preserve = "single"), show.legend = FALSE)+
  coord_flip() + 
  geom_label(aes(label = Enh_count), nudge_y= 0.07, show.legend = FALSE, size=6)+
  theme(text = element_text(size = 20))
#labs(title = "Enhancers", 
#     caption = "Normalised acording to number of observation \n in each classifiaction")

ggsave(plot = enhancer_plot, 
       filename = "results/05_enhancer_freq.png")
