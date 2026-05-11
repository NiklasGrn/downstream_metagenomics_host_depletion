library(phyloseq)
library(ggplot2)

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("phyloseq", "microbiome", "ComplexHeatmap"), update = FALSE)

install.packages(
  "microViz",
  repos = c(davidbarnett = "https://david-barnett.r-universe.dev", getOption("repos"))
)
library(microViz)
library(tidyverse)

## Importieren und Zusammenfügen der daten sowie erstellung des phyloseq objects

ps_1 <- import_biom("data/grapes_sample1.biom")
ps_2 <- import_biom("data/grapes_sampleW3_1.biom")

sample_names(ps_1) <- paste0("S1_", c("Control", "AS", "Enrichment", "Enrichment + AS"))
sample_names(ps_2) <- paste0("W3_1_", c("Control", "AS", "Enrichment", "Enrichment + AS"))


df_1 <- data.frame(
  Behandlung = c("Control", "AS", "Enrichment", "Enrichment + AS"),
  Batch = "Sample 1",
  row.names = sample_names(ps_1)
)

df_2 <- data.frame(
  Behandlung = c("Control", "AS", "Enrichment", "Enrichment + AS"),
  Batch = "Sample W3.1",
  row.names = sample_names(ps_2)
)

sample_data(ps_1) <- sample_data(df_1)
sample_data(ps_2) <- sample_data(df_2)

tax_cols <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
colnames(tax_table(ps_1)) <- tax_cols
colnames(tax_table(ps_2)) <- tax_cols

ps_merged <- merge_phyloseq(ps_1, ps_2)

ps <- ps_merged %>%
  tax_fix(unknowns = c("", " ", "unclassified"))

sample_data(ps)$Behandlung <- factor(
  sample_data(ps)$Behandlung, 
  levels = c("Control", "AS", "Enrichment", "Enrichment + AS")
)


tax_table(ps) <- apply(tax_table(ps), 2, function(x) gsub("^[a-zA-Z]_+", "", x))
tax_table(ps)[, "Species"] <- paste(tax_table(ps)[, "Genus"], tax_table(ps)[, "Species"], sep = " ")

ps <- ps %>%
  tax_fix(
    unknowns = c(" ", "", "NA", "unclassified", "unknown", "Ambiguous_taxa"),
    sep = " ", 
    min_length = 4
  )

# Erstellung des plots zur visualisierung der relativen abundance aller proben
p1 <- ps %>%
  ps_mutate(
    Behandlung = factor(Behandlung, levels = c("Enrichment + AS", "Enrichment",  "AS", "Control"))
  ) %>%
  comp_barplot(
    tax_level = "Species", 
    n_taxa = 15, 
    label = "Behandlung", 
    other_name = "Other",
    taxon_renamer = function(x) stringr::str_remove(x, " [ae]t rel."),
    palette = distinct_palette(n = 15, add = "grey90"),
    merge_other = TRUE, 
    bar_outline_colour = "darkgrey"
  ) +
  facet_grid(Batch ~ ., scales = "free", space = "free") + 
  coord_flip() +
  labs(
    title = "Taxonomic Composition at Species Level",
    subtitle = "Comparison between Sample 1 and Sample W3.1 (relative abundance)",
    x = "Treatment", 
    y = "Relative Abundance"
  ) +
  theme(
    strip.text.y = element_text(angle = 0, face = "bold"),
    legend.text = element_text(size = 8, face = "italic"),
    panel.spacing = unit(2, "lines")
  )

p1 
###### Absolute abundance

# Manuelle eingabe der gesamtanzahl der reads nach filtering 
raw_read_counts <- c(
  "S1_Control" = 750637, 
  "S1_AS" = 436354,
  "S1_Enrichment" = 578145,
  "S1_Enrichment + AS" = 575248,
  "W3_1_Control" = 793526,
  "W3_1_AS" = 323806,
  "W3_1_Enrichment" = 635193,
  "W3_1_Enrichment + AS" = 375193
)

# extrahieren der top 15 taxa und normalisierung der reads und extraktion der otu daten aus dem ps object 
# absolute abundance nur schwer mittels phyloseq darstellbar daher manuelle extraktion der anzahl der reads aus diesem
# ps object
top15_names <- names(sort(taxa_sums(ps), decreasing = TRUE)[1:15])
otu_mat <- as(otu_table(ps), "matrix")
if (taxa_are_rows(ps)) { otu_mat <- t(otu_mat) }

plot_data_rpm <- as.data.frame(otu_mat) %>%
  rownames_to_column("sample_id") %>%
  pivot_longer(cols = -sample_id, names_to = "OTU", values_to = "Reads") %>%
  
  mutate(Raw_Total = raw_read_counts[sample_id]) %>%

  left_join(as.data.frame(as(tax_table(ps), "matrix")) %>% rownames_to_column("OTU"), by = "OTU") %>%
  mutate(Species_Plot = if_else(OTU %in% top15_names, Species, "Other")) %>%

  group_by(sample_id, Species_Plot, Raw_Total) %>%
  summarise(Reads_Sum = sum(Reads), .groups = "drop") %>%

  mutate(RPM = (Reads_Sum / Raw_Total) * 1000000) %>%
  left_join(as.data.frame(as(sample_data(ps), "matrix")) %>% rownames_to_column("sample_id"), by = "sample_id")


species_order <- plot_data_rpm %>%
  group_by(Species_Plot) %>%
  summarise(total = sum(RPM)) %>%
  arrange(desc(total)) %>%
  pull(Species_Plot)
species_order <- c(setdiff(species_order, "Other"), "Other")

plot_data_rpm <- plot_data_rpm %>%
  mutate(
    Behandlung = str_replace(Behandlung, "Enrichment_AS", "Enrichment + AS"),
    Behandlung = factor(Behandlung, levels = rev(c( "Enrichment + AS", "Enrichment", "AS", "Control"))),
    Species_Plot = factor(Species_Plot, levels = species_order)
  )


my_colors <- c(distinct_palette(n = length(species_order)-1, add = NULL), "grey90")
names(my_colors) <- species_order

p2 <- ggplot(plot_data_rpm, aes(y = Behandlung, x = RPM, fill = Species_Plot)) +
  geom_bar(stat = "identity", position = position_stack(reverse = TRUE), 
           color = "darkgrey", linewidth = 0.2, width = 1) +
  facet_grid(Batch ~ ., scales = "free_y", space = "free_y") +
  scale_fill_manual(values = my_colors) +
  scale_x_continuous(labels = scales::label_number(big.mark = " "), 
                     expand = expansion(mult = c(0, 0.05))) +
  theme_bw() +
  labs(
    title = "Taxonomic Composition at Species Level",
    subtitle = "Comparison between Sample 1 and Sample W3.1 (Normalized to total reads before host filtering)",
    x = "Reads Per Million (relative to data before host filtering)",
    y = "Treatment",
    fill = "Species"
  ) + 
  theme(strip.text.y = element_text(angle = 0, face = "bold"),
        legend.text = element_text(size = 8, face = "italic"),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"))
 

plot_data_rpm %>%
  filter(sample_id == "W3_1_Control") %>%
  summarize(total_sum = sum(RPM, na.rm = TRUE))

# Speichern der erstellten Plots 
ggsave("plots/grapes_relative_abundnance.png", plot = p1, height = 740, width = 1200,units = "px", dpi = 100)
ggsave("plots/grapes_absolute_abundance.png", plot = p2, height = 740, width = 1200,units = "px", dpi = 100)
