setwd("working directory")

GTEX <- read.csv("csv dataset1", header = FALSE)
HPA <- read.delim("tsv dataset ")

#install.packages("tidyverse")
#install.packages("ggplot2")

library(tidyverse)
library(ggplot2)


#converting protein abundance dataset to long format, 
#removing metadata rows
header_row <- GTEX [3,]
colnames(GTEX) <- c("gene_id", as.vector(unlist(header_row[ ,-1])))

GTEX_abundance <- GTEX[-c(1:3),]

#converted the GTEX file to long format
GTEX_abun_long <- GTEX_abundance %>%
  pivot_longer(
    cols = -gene_id,
    names_to = "Tissue",
    values_to = "Abundance"
  )




#filtering protein amount data by levels and in order
HPA_filtered <- HPA %>% 
  select(!IHC.tissue.name) %>%
  filter(!Level %in% c("N/A", "Not representative")) %>%
  filter(!Reliability %in% c("Uncertain")) %>%
  mutate(Level = factor(Level, levels = c("Not detected","Descending","Ascending","Low", "Medium", "High"), ordered = TRUE)) %>%
  mutate(Reliability = factor(Reliability, levels = c("Supported","Approved", "Enhanced"), ordered = TRUE)) %>%
  mutate(
    Level.num = case_when(
      Level == "Not detected" ~ 0,
      Level == "Descending" ~ 1,
      Level == "Ascending" ~ 2,
      Level == "Low" ~ 3,
      Level == "Medium" ~ 4,
      Level == "High" ~ 5,
      TRUE ~ NA_real_
    )
  )

# Grouping HPA to show only 1 gene-tissue entry 
HPA_Grouped <- HPA_filtered %>%
  group_by(Gene, Tissue, Gene.name) %>%
  summarise(
    Abundance = max(Level.num, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    Abundance = case_when(
      Abundance == 0 ~ "Not detected",
      Abundance == 1 ~ "Descending",
      Abundance == 2 ~ "Ascending",
      Abundance == 3 ~ "Low",
      Abundance == 4 ~ "Medium",
      Abundance == 5 ~ "High",
      TRUE ~ NA_character_
    )
  ) %>%
    mutate(
      presence = case_when(
        Abundance == "Not detected" ~ "absent",
        .default = "present"
      ) 
    )





#comparing tissues between dataset
GTEXSORTED <- unique(GTEX_abun_long$Tissue)
HPASORTED <- unique(HPA_Grouped$Tissue)

GTEXSORTED       
HPASORTED

shared_tissues <- intersect(GTEXSORTED, HPASORTED)
OnlyGTEX <- setdiff(GTEXSORTED, HPASORTED)
OnlyHPA <- setdiff(HPASORTED, GTEXSORTED)

shared_tissues
OnlyGTEX
OnlyHPA




#good one, renames GTEX tissues into HPA format
GTEX_final <- GTEX_abun_long %>%
  rename(Gene = gene_id) %>%
  mutate(
    Tissue = recode_values(
      Tissue, 
      "Adrenal Gland" ~ "Adrenal gland", 
      "Brain Cerebellum" ~ "Cerebellum",
      "Brain Cortex" ~ "Cerebral cortex",
      "Colon Sigmoid" ~  "Colon",
      "Colon Transverse" ~ "Colon",
      "GE junction" ~ "Esophagus",
      "Esophagus Mucosa" ~ "Esophagus",
      "Esophagus Muscle" ~ "Esophagus",
      "Heart Atrial" ~ "Heart muscle",
      "Heart muscle" ~ "Heart muscle",
      "Heart Ventricle" ~ "Heart muscle",
      "Minor Salivary" ~ "Salivary gland",
      "Muscle Skeletal" ~ "Skeletal muscle",
      "Pituitary" ~ "Pituitary gland", 
      "Skin SunExpo" ~ "Skin",
      "Skin Unexpo" ~ "Skin",
      "Small Intestine" ~ "Small intestine", 
      "Thyroid" ~ "Thyroid gland",
      "Breast" ~ "Breast",
      "Liver" ~ "Liver",
      "Lung"~ "Lung",
      "Ovary" ~ "Ovary",
      "Pancreas" ~ "Pancreas",
      "Prostate" ~ "Prostate",
      "Spleen" ~ "Spleen",
      "Stomach" ~ "Stomach",
      "Testis" ~ "Testis",
      "Vagina" ~ "Vagina",
      default = "no match"
    )
  ) %>%
  filter(!Tissue %in% "no match") %>%
  mutate(Abundance = as.numeric(Abundance)) %>%
  group_by(Gene) %>%
  mutate(Heart = if_else(
    Tissue == "Heart muscle",
    mean(Abundance[Tissue %in% "Heart muscle"], na.rm = TRUE),
    NA_real_
    )
  ) %>%
  mutate(Esophagus = if_else(
    Tissue == "Esophagus",
    mean(Abundance[Tissue %in% "Esophagus"], na.rm = TRUE),
    NA_real_
    )
  ) %>%
  mutate(Colon = if_else(
    Tissue == "Colon",
    mean(Abundance[Tissue %in% "Colon"], na.rm = TRUE),
    NA_real_
    )
  ) %>%
  mutate(Skin = if_else(
    Tissue == "Skin",
    mean(Abundance[Tissue %in% "Skin"], na.rm = TRUE),
    NA_real_
    )
  ) %>%
  mutate(
    presence = case_when(
      is.na(Abundance) ~ "absent",
      .default = "present"
    ) 
  )

GTEX_final_corrected <- GTEX_final %>%
  mutate(Abundance = case_when(
    !is.na(Abundance) ~ as.character(Abundance),
    is.na(Abundance) ~ "NA"
    )
  )



matching_tissues <- c(
  "Adrenal gland",
  "Cerebellum",
  "Cerebral cortex",
  "Breast",
  "Colon",
  "Esophagus",      
  "Heart muscle",
  "Liver",
  "Lung",
  "Salivary gland",
  "Skeletal muscle",
  "Ovary",
  "Pancreas",
  "Pituitary gland",
  "Prostate",
  "Skin",
  "Small intestine",
  "Spleen",
  "Stomach",
  "Testis",
  "Thyroid gland",
  "Vagina"
)

HPA_final <- HPA_Grouped %>%
  mutate(
    Tissue = case_when(
      Tissue %in% matching_tissues ~ Tissue,
      TRUE ~ "no match"
    )
  ) %>%
  filter(!Tissue %in% "no match")

#final combined dataframe
combined <- bind_rows(GTEX_final_corrected, HPA_final) %>%
  mutate(HPAorGTEX = case_when(
    Gene.name %in% NA_character_ ~ "GTEX",
    Gene.name %in% Gene.name ~ "HPA"
    )
  )





#calculations and analysis 
#histogram and density plot
com_pres <- combined %>% count(Gene, presence)
com_pres_pres <- com_pres %>%
  filter(presence == "present")
com_pres_abs <-com_pres %>%
  filter(presence == "absent")

combined %>% count(Tissue, presence)
combined %>% count(Tissue, HPAorGTEX, presence)
tissue_percent <- combined %>%
  group_by(Tissue) %>%
  summarise(
    n_present = sum(presence == "present"),
    n_total = n(),
    percent_present = n_present/n_total
  )

sourced_tissue_percent <- combined %>%
  group_by(Tissue, HPAorGTEX) %>%
  summarise(
    n_present = sum(presence == "present"),
    n_total = n(),
    percent_present = n_present/n_total,
    .groups = "drop"
  )

# observations: HPA sourced tissues seem to have a lower presence detected overall 
#               while GTEX tissues have a higher presence percentage. only except 
#               is the Pituitary gland proteins with a 95% presence percentage from HPA. 
#               their specific protein structure might make it easier to detect (shorter proteins maybe? )
  
#same shit but with tissue resolution now
tissue_status <- combined %>%
  group_by(Gene, Tissue, HPAorGTEX) %>%
  summarise(
    present_any = any(presence == "present"),
    .groups = "drop"
  )

gene_summary <- tissue_status %>%
  group_by(Gene, Tissue) %>%
  summarise(
    n_platforms_present = sum(present_any), 
    platforms_present = paste(HPAorGTEX[present_any], collapse = ", "),
    .groups = "drop"
  )




ggplot(sourced_tissue_percent, aes(x = percent_present, ))
















#attempted plots 
p1_pres <- com_pres_pres
ggplot(com_pres_pres,aes(x = Gene, fill = n)) +
  geom_histogram(binwidth = 20,
                 show.legend = F,
                 alpha = .5) +
  labs(title = "Histogram",
       x = "gene",
       y = "Count")

ggplot(com_pres, aes(x = gene, y = presence, fill = n)) + 
  geom_bar(stat = "n")
  
  
  
combined_presence <- combined %>%
  select(Gene, Tissue, presence, Gene.name) %>%
  # mutate(presence = case_when(
  #   presence %in% "present" ~ 1,
  #   presence %in% "absent" ~ 0,
  #   )
  # ) %>%
  mutate(HPAorGTEX = case_when(
    Gene.name %in% NA_character_ ~ "GTEX",
    Gene.name %in% Gene.name ~ "HPA"
    )
  ) %>%
  mutate(grouped_presence = case_when(
    presence == "present" & HPAorGTEX == "GTEX" ~ "present_GTEX",
    presence == "present" & HPAorGTEX == "HPA" ~ "present_HPA",
    presence == "absent" & HPAorGTEX == "GTEX" ~ "absent_GTEX",
    presence == "absent" & HPAorGTEX == "HPA" ~ "absent_HPA",
    )
  )

#creating a heat map
ggplot(combined_presence, aes(x = Tissue, y = Gene, fill = presence)) +
  geom_tile() +
  scale_fill_manual(
    values = c(
      present_GTEX = "blue",
      present_HPA = "red",
      absent_GTEX = "grey",
      absent_HPA = "white"
      )
  )


#gene only abundance
gene_abundance <- combined %>%
  group_by(Gene)
  
ggplot(combined_presence, aes(x = presence))

boxplot
