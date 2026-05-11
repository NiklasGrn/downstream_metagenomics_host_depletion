library(scales)
library(ggplot2)
library(tidyverse)
library(phyloseq)
library(textshape)

# Zuweisung der verwendeten Farben für die Plots
my_colors <- c("#1F78B4", "#964B00", "#6A3D9A", "#FF7F00", "#E31A1C", "#33A02C", "#999999", "#E6B800")
mein_header <- c("Target_Name", "Length", "Reads", "Unmapped")


#Importieren der dateien aus ./data/mapping/... 
bulk_files <- c("data/mapping/barcode01_bulk.tsv", "data/mapping/barcode02_bulk.tsv", "data/mapping/barcode03_bulk.tsv", "data/mapping/barcode04_bulk.tsv")
as_files   <- c("data/mapping/barcode01_AS.tsv", "data/mapping/barcode02_AS.tsv", "data/mapping/barcode03_AS.tsv", "data/mapping/barcode04_AS.tsv")


# Manuelle zuweisung der Metadaten und der gesamtanzahl der Reads pro probe
df_reads_wgs <- data.frame(
  reads = c(294.101,295.063,486.39,469.139,297.939,299.126,495.3,477.499),
  Condition = c("Bulk","Bulk","Bulk","Bulk","AS","AS","AS","AS"),
  Replicate = c("Rep_1","Rep_2","Rep_3","Rep_4","Rep_1","Rep_2","Rep_3","Rep_4")
)




# Importieren der daten 
load_all_replicates <- function(files, label) {
  lapply(seq_along(files), function(i) {
    df <- read.delim(files[i], header = FALSE, skip = 0, col.names = c("Target_Name", "Length", "Reads", "Unmapped"))
    
    df %>%
      mutate(Group = case_when(
        str_detect(Target_Name, "tig") ~ "S._cerevisiae",
        str_detect(Target_Name, "BS.pilon.polished") ~ "B._subtilis",
        str_detect(Target_Name, "Enterococcus_faecalis") ~ "E._faecalis",
        str_detect(Target_Name, "Escherichia_coli") ~ "E._coli",
        str_detect(Target_Name, "Listeria_monocytogenes") ~ "L._monocytogenes",
        str_detect(Target_Name, "Pseudomonas_aeruginosa") ~ "P._aeruginosa",
        str_detect(Target_Name, "Salmonella_enterica") ~ "S._enterica",
        str_detect(Target_Name, "Staphylococcus_aureus") ~ "S._aureus",
        TRUE ~ "Other"
      )) %>%
      
      group_by(Group) %>%
      summarise(Total_Reads = sum(as.numeric(Reads), na.rm = TRUE)) %>%
      mutate(Condition = label, Replicate = paste0("Rep_", i))
  }) %>% bind_rows()
}

all_data_alt %>% 
  filter(Group == "S._cerevisiae") %>% 
  select(Total_Reads)

all_data_alt <- bind_rows(
  load_all_replicates(bulk_files, "Bulk"),
  load_all_replicates(as_files, "AS")
)

# Berechnung der normalisierten reads für alle proben
all_data <- all_data_alt %>%
  left_join(df_reads, by = c("Condition", "Replicate")) %>%
  mutate(
    CPM = (Total_Reads / reads) * 1e3
  )

# Sätzung der gewünschten reihenfolge innerhalb der plots
reihenfolge <- c("S. cerevisiae", 
                 "B. subtilis", 
                 "E. faecalis", 
                 "E. coli", 
                 "L. monocytogenes", 
                 "P. aeruginosa", 
                 "S. enterica", 
                 "S. aureus",
                 "Other",
                 "Unmapped")

all_data$Group <- gsub("_"," ",x=all_data$Group)

all_data$Group

all_data$Group <- factor(all_data$Group, levels = reihenfolge)
all_data <- all_data %>% filter(Group != "Other") %>% mutate(Group = droplevels(Group))
all_data <- all_data %>% select(Group,Condition,Replicate,CPM)

# erstellung des plots zur darstellung wie sich die anzahl der reads verändert hat vor und nach AS
p1 <- ggplot(all_data, aes(x = Group, y = CPM, fill = Condition)) +
  stat_summary(fun = mean, geom = "bar", 
               position = position_dodge(width = 0.9), 
               color = "black", width = 0.8) +
  
  stat_summary(fun = mean, geom = "text", 
               aes(label = label_number(big.mark = " ")(after_stat(round(y, 0)))), 
               position = position_dodge(width = 0.9), 
               vjust = -0.5,
               hjust = -0.2,
               angle = 90,
               size = 4.2) + 
  
  stat_summary(fun.data = mean_sdl, fun.args = list(mult = 1), 
               geom = "errorbar", 
               position = position_dodge(width = 0.9), 
               width = 0.2) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  scale_fill_manual(values = c("Bulk" = "#999999", "AS" = "#e6e600")) +
  labs(title = "Read distribution (n=4) of normalized reads", y = "Reads per million (Mean)", x = "") +
  theme_classic() +
  theme(axis.text.x = element_text(size = 12, angle = 45, hjust = 1, vjust = 1 , face = "italic"))


# Berechnung des fold changes 
fold_change_data <- all_data %>%
  pivot_wider(names_from = Condition, values_from = CPM) %>%
  mutate(Fold_Change = AS / ifelse(Bulk == 0, 1, Bulk)) %>%
  mutate(log2_FC = log2(Fold_Change))

fc_summary <- fold_change_data %>%
  group_by(Group) %>%
  summarise(
    Mean_FC = mean(Fold_Change, na.rm = TRUE),
    SD_FC = sd(Fold_Change, na.rm = TRUE),
    n = n()
  )

# Erstellung des plots für die darstellung des Fold changes nach der Verwendung AS
p2 <- ggplot(fold_change_data, aes(x = Group, y = Fold_Change, fill = Group)) +
  stat_summary(fun = mean, geom = "bar", color = "black", width = 0.7) +
  stat_summary(fun.data = mean_sdl, fun.args = list(mult = 1), geom = "errorbar", width = 0.2) +
  # geom_jitter(width = 0.1, size = 2, alpha = 0.8) + 
  geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
  scale_fill_manual(values = my_colors)  +
  labs(title = "X-fold change AS/Bulk (n=4)", 
       y = "X-fold Change", 
       x = "") +
  theme_classic() +
  theme(axis.text.x = element_text(size = 12, angle = 45, hjust = 1, face = "italic"), legend.position = "none")

df_wide <- all_data_alt %>%
  unite("Sample_ID", Condition, Replicate, sep = "_") %>%
  pivot_wider(
    names_from = Group,         
    values_from = Total_Reads   
  )

# Einlesen der daten für den standard zur darstellung der relativen abundance im vergleich zu dem standard
standard_names <- c("S. cerevisiae", "B. subtilis", "E. faecalis", "E. coli", "L. monocytogenes", "P. aeruginosa", "S. enterica", "S. aureus")
standard_amount <- c(2, 14, 14, 14, 14, 14, 14, 14)

standard <- data.frame(Group = standard_names, CPM = standard_amount, Condition = c("Standard"), Replicate = c("Standard"))
all_plus_standard <- bind_rows(all_data, standard)


# Erstellen eines phyloseq objects zur darstellung der relativen abundance aus den normalisieren reads
otu_data <- all_plus_standard %>%
  mutate(SampleID = paste(Condition, Replicate, sep = "_")) %>%
  select(Group, SampleID, CPM) %>%
  pivot_wider(names_from = SampleID, values_from = CPM) %>%
  column_to_rownames("Group") %>%
  as.matrix()


tax_data <- data.frame(
  Species = rownames(otu_data),
  row.names = rownames(otu_data)
) %>% as.matrix()


sample_data_df <- data.frame(
  SampleID = colnames(otu_data)) %>%
  separate(SampleID, into = c("Condition", "Replicate"), sep = "_", remove = FALSE) %>%
  column_to_rownames("SampleID")

OTU = otu_table(otu_data, taxa_are_rows = TRUE)
TAX = tax_table(tax_data)
SAM = sample_data(sample_data_df)


ps <- phyloseq(OTU, TAX, SAM)


print(ps)

ps_rel <- transform_sample_counts(ps, function(x) x / sum(x) * 100)

my_colors_2 <- c("#964B00", "#FF7F00", "#6A3D9A","#E31A1C", "#33A02C","#E6B800", "#1F78B4", "#999999")
metadata$Display_Name <- c("Replicate 1","Replicate 2","Replicate 3","Replicate 4","Replicate 1","Replicate 2","Replicate 3","Replicate 4","Standard")

sample_data(ps_rel) <- sample_data(metadata)
sample_data(ps_rel)

# Erstellung des plots für die darstellung der relativen abundance
p3 <- plot_bar(ps_rel,x = "Display_Name", fill = "Species") + 
  facet_wrap(~Condition, scales = "free_x") + 
  
  scale_fill_manual(values = my_colors_2) +
  theme_bw() + 
  labs(y = "Relative Abundanz (%)", x = "") +
  
  theme(
    strip.background = element_rect(fill = "white"), 
    strip.text = element_text(face = "bold", size = 12),
    axis.text.x = element_text(angle = 45, hjust = 1), 
    legend.text = element_text(face = "italic") 
  )

# Speichern der erstellten plots
ggsave("plots_wgs/zymo_read_distribution_including_SD.png" , plot=p1, width = 740, height = 740, units = "px", dpi = 100)
ggsave("plots_wgs/X-fold_change_zymo_WGS.png" , plot=p2, width = 740, height = 740, units = "px", dpi = 100)
ggsave("plots_wgs/zymo_relative_abundance_AS_bulk.png" , plot=p3, width = 740, height = 740, units = "px", dpi = 100)
