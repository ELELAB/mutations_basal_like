# Initilization------------------------------------------------------
#devtools::install_github("yanlinlin82/ggvenn")
library(ggvenn)
library(tidyverse)

# Functions
is_outlier <- function(x) {
  return(x < quantile(x, 0.001) - 1.5 * IQR(x) | x > quantile(x, 0.999) + 1.5 * IQR(x))
}

# Load Data ---------------------------------------------------------
MafFile <- read_tsv("data/mutations.tsv")
DEGene_All_Mutation_Annotations <- read_csv("results/DEGene_All_Mutation_Annotations.csv")
#DEGene_Mutation_Summary <- read_csv("results/DEGene_Mutation_Summary.csv")



# Wrangle data ------------------------------------------------------

# Find dimensions/numbers of mutatations and genes
MafFile %>% dplyr::select(Hugo_Symbol) %>% distinct() %>% count
dim(MafFile)
DEGene_All_Mutation_Annotations %>% select(Hugo_Symbol) %>% distinct() %>% dim()
dim(DEGene_All_Mutation_Annotations)



# Visualisations of Mutations Distribution
notDEGs <- MafFile %>% 
  anti_join(DEGene_All_Mutation_Annotations, by = 'Hugo_Symbol')# %>% select(Hugo_Symbol) %>%  distinct()# %>% dim()

t <- notDEGs %>% 
  group_by(Hugo_Symbol) %>%  summarise(n = n())%>% 
  mutate(Type = 'nonDEG') %>%  mutate(logn = log10(n))

v <- DEGene_All_Mutation_Annotations %>% 
  group_by(Hugo_Symbol) %>% count() %>%
  mutate(Type = 'DEG') %>%  mutate(logn = log10(n))

both <- full_join(t,v) %>% 
  mutate(logn = log10(n), 
         outlier = ifelse(is_outlier(n), Hugo_Symbol, as.numeric(NA))) 

# Boxplot with outliers
ggplot(data =both, aes(y = n, x = Type,  color = Type), alpha = 0.5, show.legend = FALSE) +
  geom_boxplot() +
  theme(text = element_text(size = 20))+
  geom_text(aes(label = outlier, size = 12), na.rm = TRUE, check_overlap = TRUE, nudge_x = 0.15)#, hjust = -1.2, vjust = -1.6)

# Bar plot Degs vs NonDEGS no outliers
denisty_deg_nodeg <- both %>% filter(is.na(outlier)) %>% 
  ggplot(aes(x = n,  fill = Type)) +
  geom_bar(position = "dodge2")  +
  theme(text = element_text(size = 20))
ggsave(denisty_deg_nodeg,
       filename = "results/02_density_degs_vs_nodegs.png", dpi = 300) #units = "in")

# Dotplot of outliers, which was not included in above barplot
outliers_dotplot <- both %>% filter(!is.na(outlier)) %>% mutate(y = 0) %>% 
  ggplot(aes(x = n, y = 0, fill = Type), alpha = 0.6) +
  geom_dotplot(stackgroups = TRUE) +
  theme(text = element_text(size = 20)) +
  scale_y_continuous(NULL, breaks = NULL) +
  geom_text(aes(label = Hugo_Symbol, size = 12), show.legend = FALSE, na.rm = TRUE, 
            check_overlap = TRUE, vjust = -3) #nudge_x = 0.15, nodge_y = 0.3)#, hjust = -1.2, vjust = -1.6)
ggsave(outliers_dotplot,
       filename = "results/02_outliers_dotplot.png", dpi = 300)



# Venn-diagram (not used)------------------------
uniq_nondegs <- notDEGs %>% select(Hugo_Symbol) %>% 
  unique() %>% 
  mutate(Type = "nonDEG", 
         Genes = TRUE)

uniq_degs <- DEGene_All_Mutation_Annotations %>% select(Hugo_Symbol) %>% 
  unique() %>% 
  mutate(Type = 'DEG', Genes = TRUE)

ocg <- DEGene_All_Mutation_Annotations %>% filter(Driver_type == 'OCG') %>% 
  select(Hugo_Symbol) %>% unique() %>% mutate(Type = 'OCG', Genes = TRUE)

tsg <- DEGene_All_Mutation_Annotations %>% filter(Driver_type == 'TSG') %>% 
  select(Hugo_Symbol) %>% unique() %>% mutate(Type = 'TSG', Genes = TRUE)

x <- full_join(uniq_nondegs,uniq_degs) %>% full_join(ocg) %>% full_join(tsg) %>% 
  pivot_wider(id_cols = Hugo_Symbol, 
              names_from = Type, 
              values_from = Genes) %>% 
  mutate( across(everything(), ~replace_na(.x, FALSE)))

ggvenn(data = x, columns = c("DEG", "OCG", "TSG"))


# Types and Classifications in DEGS ---------------------------------------
# Distribution of variant Types in variant classifications in the DEGs
DEGs_type_class <- ggplot(DEGene_All_Mutation_Annotations, 
       aes(x = Variant_Classification, fill = Variant_Type), alpha = 0.8) +
  geom_bar() +
  coord_flip() +
  labs(title = "The distribution of mutations in DEGs")
ggsave(DEGs_type_class,
       filename = "results/02_DEGs_type_class.png", dpi = 300)


# Frequency of distribution 
# The distribution with freq compared to the MAF
Deg_dist <- DEGene_All_Mutation_Annotations %>% 
  group_by(Variant_Classification, Variant_Type) %>% 
  summarise(Class_count = n())

maf_freq <- MafFile %>% 
  group_by(Variant_Classification) %>% 
  summarise(MAF_count = n()) %>% 
  left_join(Deg_dist) %>% 
  mutate(Freq = Class_count/MAF_count)

#Expected freq
total_genes <- MafFile %>% select(Hugo_Symbol) %>% unique() %>% 
  count() %>% as.numeric()
deg_genes <- DEGene_All_Mutation_Annotations %>%  
  select(Hugo_Symbol) %>% unique() %>% 
  count() %>% as.numeric()
exp_freq <- deg_genes/total_genes

DEGs_type_class_freq <- ggplot(maf_freq, aes(x = Variant_Classification, 
                     y = Freq,
                     fill = Variant_Type), alpha = 0.7)+
  geom_col(show.legend = FALSE)+
  coord_flip() +
  #geom_text(aes(label = Class_count), nudge_y= 0.02)+
  labs(title = "DEGs Frequency compared to total in MAF", 
       caption = paste("The expected Freqency", round(exp_freq, digits = 2)))
ggsave(DEGs_type_class_freq,
       filename = "results/02_DEGs_type_class_frequency.png", dpi = 300)


# Mutations ------------------------------------------------
#Driver mutations vs Variant class --------------------
mutations_class <- DEGene_All_Mutation_Annotations %>% 
  filter(Driver_Mutation !="No_mutations") %>% 
  group_by(Driver_Mutation, Variant_Classification) %>% count() %>% 
  ggplot(aes(x = Variant_Classification, y = n, fill = Driver_Mutation))+
  geom_col(position = "dodge2") +
  theme(text = element_text(size = 20)) +
  coord_flip()

ggsave(mutations_class,
       filename = "results/02_driver_mutations_class.png", 
        dpi=300)


#Up and down regultation against number of mutations -------
regulation_mutations <- DEGene_All_Mutation_Annotations %>% 
  group_by(Hugo_Symbol, logFC, Driver_type) %>% count() %>% 
  ggplot(aes(x=logFC, y = n, color = Driver_type),alpha = 0.5) +
  geom_point() +
  theme(text = element_text(size = 20)) 

ggsave(regulation_mutations,
       filename = "results/02_regulation_mutations.png", 
       height=6, width=12, units = "in") #, dpi=300)
