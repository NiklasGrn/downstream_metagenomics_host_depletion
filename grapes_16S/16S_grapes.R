library(tidyverse)
library(phyloseq)
library(stringr)
library(readr)
install.packages(
  "microViz",
  repos = c(davidbarnett = "https://david-barnett.r-universe.dev", getOption("repos"))
)
library(microViz)


# Importieren der ddaten aus ./data/*.tsv und der metadaten aus ./metadaten.csv
file_names <- list.files("data/", pattern = "\\.tsv$", full.names = TRUE)


all_emu_data <- file_names %>%
  map_dfr(function(x) {
    read_tsv(x) %>%
      mutate(sample_id = basename(x) %>% str_remove(".tsv")) 
  })


metadata <- read_csv2("metadata.csv") %>% 
  column_to_rownames("sample_id")


# Erstellung der OTU Matrix sowie TAXonomy Matrix für das pyloseq object 
otu_mat <- all_emu_data %>%
  select(tax_id, sample_id, `estimated counts`) %>%
  group_by(tax_id, sample_id) %>%
  summarise(estimated_counts_sum = sum(`estimated counts`, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = sample_id, values_from = estimated_counts_sum, values_fill = 0) %>%
  column_to_rownames("tax_id") %>%
  as.matrix()


tax_mat <- all_emu_data %>%
  select(tax_id, superkingdom, phylum, class, order, family, genus, species) %>%
  distinct(tax_id, .keep_all = TRUE) %>%
  column_to_rownames("tax_id") %>%
  as.matrix()

ps <- phyloseq(
  otu_table(otu_mat, taxa_are_rows = TRUE),
  tax_table(tax_mat),
  sample_data(metadata)
)


# Fixen des ps objectes und dann dastellung in korrekter Reihenfolge (BC1->BC20) auf Family taxa
# Auftrennung nach sampling methode
ps_final <- ps %>%
  tax_fix() %>% 
  tax_mutate(Species_Full = paste(genus, species, sep = " "))

ps_final %>% tax_fix()

ps_final <- ps_final %>%
  ps_mutate(sample_id_col = sample_names(.)) %>% 
  ps_arrange(sample_id_col)

p1 <- ps_final %>%
  comp_barplot(
    tax_level = "family", 
    n_taxa = 12, 
    facet_by = "sampling",
    sample_order = rev(sort(sample_names(ps_final)))
  ) +
  coord_flip()

p1


# Extrahieren der daten aus phyloseq object um dannach die Mitochondira, Chloroplasten und Other
# Gruppen zu erstellen und dann wieder als phyloseq object zu speichern
tax_df <- as.data.frame(tax_table(ps_final)@.Data)


tax_df$custom_group <- tax_df$family 

tax_df$custom_group[tax_df$custom_group %in% c("Mitochondria")] <- "Mitochondria"
tax_df$custom_group[tax_df$custom_group %in% c("Incertae_Sedis")] <- "Chloroplast"
tax_df$custom_group[!tax_df$custom_group %in% c("Mitochondria", "Chloroplast")] <- "Other"


tax_table(ps_final) <- tax_table(as.matrix(tax_df))

ps_final %>% tax_fix()


ps_custom <- ps_final
tax_table(ps_custom) <- tax_table(ps_final)[, "custom_group", drop = FALSE]

ps_custom <- ps_custom %>%
  ps_mutate(sample_id_col = sample_names(.)) %>% 
  ps_arrange(sample_id_col)


# Erstellung des Plots zur darstellung Mitochondria vs Chloroplasten vs Other
p2 <- ps_custom %>%
  comp_barplot(
    tax_level = "custom_group",  
    n_taxa = 12, 
    label = "barcode",
    facet_by = "sampling",
    sample_order = rev(sort(sample_names(ps_custom)))
  ) +
  coord_flip() +
  labs(fill = "Taxonomic Group", title = "Relative taxonomic abundance of mitochondrial and chloroplast DNA in AS vs. bulk sequencing.")

p2

p3 <- ps_custom %>%
  ps_filter(sampling == "bulk") %>%
  comp_barplot(
    tax_level = "custom_group",  
    n_taxa = 12, 
    label = "barcode",
    facet_by = "pna",
    sample_order = rev(sort(sample_names(.)))
  ) +
  coord_flip() +
  labs(fill = "Taxonomic Group", title = "Relative taxonomic abundance of mitochondrial and chloroplast DNA", subtitle = "Effect of PNA clamp concentration on host depletion")

ps_custom_mapped <- ps_custom %>%
  ps_mutate(
    temp = factor(temp, levels = c("normal", "60", "65", "70"))
  )

p4 <- ps_custom_mapped %>%
  ps_filter(sampling == "bulk") %>%
  comp_barplot(
    tax_level = "custom_group",  
    n_taxa = 12, 
    label = "barcode",
    facet_by = "temp",
    sample_order = rev(sort(sample_names(.)))
  ) +
  coord_flip() +
  labs(fill = "Taxonomic Group", title = "Relative taxonomic abundance of mitochondrial and chloroplast DNA", subtitle = "Effect of PNA clamping temperature on host depletion")




# Filtering der nicht bacterial taxa um diese noch weiter darzustellen
host_terms <- c("Mitochondria", "mapped_filtered", "mapped_unclassified", "Incertae_Sedis")

ps_microbiome <- subset_taxa(ps_final, 
                               !family %in% host_terms & 
                               !order %in% host_terms &
                               !genus %in% host_terms)



p5 <- ps_microbiome %>%
  comp_barplot(
    tax_level = "genus",
    n_taxa = 12, 
    label = "barcode",
    facet_by = "sampling",
    sample_order = rev(sort(sample_names(.)))
    ) +
  coord_flip() +
  labs(fill = "Genus", title = "Relative taxonomic abundance of microbial community", subtitle = "Changes in the relative abundance while using AS and bulk sequencing")

# Darstellung der PCA "einfluss von AS auf microbial profile"
p6 <- ps_microbiome %>%
  tax_transform("clr", rank = "genus") %>%
  ord_calc(method = "PCA") %>%
  ord_plot(
    axes = c(1, 2),
    color = "sampling",
    shape = "pna",
    size = 4, 
    alpha = 0.8
  ) +
  stat_ellipse(aes(color = sampling)) + 
  theme_bw() +
  labs(title = "Effect of Adaptive Sampling on the microbial profile")

ggsave("plots/relative_abundance_mitochondira_chloroplast_AS_bulk.png", plot = p2, width = 1200, height = 700, units="px", dpi=100)
ggsave("plots/relative_taxonomic_abundance_PNA_conc.png", plot = p3, width = 1200, height = 700, units="px", dpi=100)
ggsave("plots/relative_taxonomic_abundance_PNA_temp.png", plot = p4, width = 1200, height = 700, units="px", dpi=100)
ggsave("plots/rel_tax_change_microbes_AS.png", plot = p5, width = 1200, height = 700, units="px", dpi=100)
ggsave("plots/PCA_AS.png", plot = p6, width = 1200, height = 700, units="px", dpi=100)

