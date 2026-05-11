library(ggplot2)
library(tidyverse)

# importieren der dateien aus ./data/mapping/16S_... .tsv
# Ändern auf den entsprechenden unterordner für die jeweiligen AS run
file_paths <- list.files(path = "data/mapping/16S_pseudo", pattern = "\\.tsv$", full.names = TRUE)

# importieren und korrigieren der namen der bacterien
df <- file_paths %>%
  set_names(basename(.)) %>%
  map_dfr(read_tsv, show_col_types = FALSE, .id = "filename") %>%
  mutate(
    extracted_barcode = str_extract(filename, "barcode\\d+"),
    extracted_condition = case_when(
      str_detect(filename, "_AS") ~ "AS",
      str_detect(filename, "_bulk") ~ "bulk",
      TRUE ~ "unknown"
    )
  ) %>%
  
  select(
    name = extracted_barcode,         
    condition = extracted_condition,  
    species = species,                   
    rel_abundance = abundance,        
    total_reads = "estimated counts"    
  ) %>%
  filter(!is.na(species))

df <- df %>%
  mutate(species = recode(species, 
                          "Salmonella enterica" = "S. enterica",
                          "Bacillus subtilis" = "B. subtilis",
                          "Enterococcus faecalis" = "E. faecalis",
                          "Staphylococcus aureus" = "S. aureus",
                          "Escherichia coli" = "E. coli",
                          "Pseudomonas aeruginosa" = "P. aeruginosa",
                          "Listeria monocytogenes" = "L. monocytogenes"))


# Normalisieren der reads aller spezies
df_final <- df %>%
  group_by(name, condition) %>%
  mutate(
    rpm = (total_reads / sum(total_reads, na.rm = TRUE)) * 1000000
  ) %>%
  ungroup()

all_data_16S <- df_final %>% select(
  Group = species,
  CPM = rpm,
  Condition = condition
)

# Erstellen des plots für die darstellung 
p1 <- ggplot(all_data_16S, aes(x = Group, y = CPM, fill = Condition)) +
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
  scale_fill_manual(values = c("bulk" = "#999999", "AS" = "#e6e600")) +
  labs(title = "Read distribution (n=4) of normalized reads", y = "Reads per million (Mean)", x = "") +
  theme_classic() +
  theme(axis.text.x = element_text(size = 12, angle = 45, hjust = 1, face = "italic"))


# Berechnen des Fold_changes 
fold_change_data <- all_data_16S %>%
  group_by(Group, Condition) %>%
  summarise(mean_CPM = mean(CPM, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = Condition, values_from = mean_CPM) %>%
  mutate(Fold_Change = AS / ifelse(bulk == 0, 1, bulk)) %>%
  mutate(log2_FC = log2(Fold_Change))

fc_summary <- fold_change_data %>%
  group_by(Group) %>%
  summarise(
    Mean_FC = mean(Fold_Change, na.rm = TRUE),
    SD_FC = sd(Fold_Change, na.rm = TRUE),
    n = n()
  )

my_colors_2 <- c("#964B00", "#FF7F00", "#6A3D9A","#E31A1C", "#33A02C","#E6B800",  "#999999")

# Erstellung des Plotes zur visualisierung des Fold Changes
p2 <- ggplot(fold_change_data, aes(x = Group, y = Fold_Change, fill = Group)) +
  stat_summary(fun = mean, geom = "bar", color = "black", width = 0.7) +
  stat_summary(fun.data = mean_sdl, fun.args = list(mult = 1), geom = "errorbar", width = 0.2) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
  scale_fill_manual(values = my_colors_2)  +
  labs(title = "X-Fold Change AS/Bulk (n=4)", 
       y = "X-Fold Change", 
       x = "") +
  theme_classic() +
  theme(axis.text.x = element_text(size = 12, angle = 45, hjust = 1, face = "italic"), legend.position = "none")

ggsave("plots_16S/16S_zymo_Read_Distribution.png" , plot=p1, width = 740, height = 740, units = "px", dpi = 100)
ggsave("plots_16S/16S_zymo_Fold_change.png" , plot=p2, width = 740, height = 740, units = "px", dpi = 100)

