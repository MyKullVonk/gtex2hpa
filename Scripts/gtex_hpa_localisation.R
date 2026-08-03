# =============================================================================
# GTEX (paper protein abundance) + HPA IHC  +  HPA subcellular location (lastest addition)
#
# Inputs: gtex.csv
#   normal_ihc_data.tsv      HPA IHC
#   subcellular_location.tsv HPA IHC
# =============================================================================

setwd("C:/Users/blueg/Documents/binf practice projects")

data_dir <- "HPA vs GTEX" # for now. change it

library(dplyr)
library(ggplot2)

# ---- shared tissue vocabulary ------------------------------------------------
# Single source of truth: names = GTEX sub-tissues, values = HPA tissue.
# Used both to harmonise GTEX names AND to define the shared tissue set.
# you decalre it once and reuse it later. if you need to copy something it means you don't do it right.
tissue_map <- c(
  "Adrenal Gland"    = "Adrenal gland",
  "Brain Cerebellum" = "Cerebellum",
  "Brain Cortex"     = "Cerebral cortex",
  "Colon Sigmoid"    = "Colon",
  "Colon Transverse" = "Colon",
  "GE junction"      = "Esophagus",
  "Esophagus Mucosa" = "Esophagus",
  "Esophagus Muscle" = "Esophagus",
  "Heart Atrial"     = "Heart muscle",
  "Heart muscle"     = "Heart muscle",
  "Heart Ventricle"  = "Heart muscle",
  "Minor Salivary"   = "Salivary gland",
  "Muscle Skeletal"  = "Skeletal muscle",
  "Pituitary"        = "Pituitary gland",
  "Skin SunExpo"     = "Skin",
  "Skin Unexpo"      = "Skin",
  "Small Intestine"  = "Small intestine",
  "Thyroid"          = "Thyroid gland",
  "Breast" = "Breast", "Liver" = "Liver", "Lung" = "Lung", "Ovary" = "Ovary",
  "Pancreas" = "Pancreas", "Prostate" = "Prostate", "Spleen" = "Spleen",
  "Stomach" = "Stomach", "Testis" = "Testis", "Vagina" = "Vagina"
)
shared_tissues <- unname(tissue_map)

# ---- 1. load -----------------------------------------------------------------
gtex_raw <- read.csv(file.path(data_dir, "NIHMS1624446-supplement-2 (1).csv"), header = TRUE, skip = 2)
hpa_raw  <- read.csv(file.path(data_dir, "normal_ihc_data.tsv"), sep = "\t")
cc_raw   <- read.csv(file.path(data_dir, "subcellular_location.tsv"), sep = "\t") # new stuff, tsv format so same deal with sep="\t"

# ---- 2. GTEX: long, harmonise tissues, collapse split sub-tissues, presence --
# The four hand-written Heart/Esophagus/Colon/Skin columns are replaced by one
# group_by(gene, tissue) %>% summarise(mean): it collapses any set of split
# subtissues. Abundance stays numeric throughout.
gtex <- gtex_raw %>%
  rename(gene = gene.id) %>%
  tidyr::pivot_longer(-gene, names_to = "tissue", values_to = "abundance") %>%
  mutate(
    tissue    = trimws(gsub("[.]", " ", tissue)),  # GTEX dots -> HPA spaces READ about trimws and gsub. 
    tissue    = unname(tissue_map[tissue]),        # harmonise; NA if no match
    abundance = as.numeric(abundance)
  ) %>%
  filter(!is.na(tissue)) %>%
  group_by(gene, tissue) %>%
  summarise(abundance = mean(abundance, na.rm = TRUE), .groups = "drop") %>%
  mutate(present = !is.na(abundance)) # NA means didn't detect.

# ---- 3. HPA IHC: one row per gene x tissue, present if detected by any entry --
# presence is driven off "Not detected" directly, so the exact Level vocabulary
# matters less.
hpa <- hpa_raw %>%

    filter(!Level %in% c("N/A", "Not representative"),
         Reliability != "Uncertain") %>%
  mutate(detected = !Level %in% c("Not detected", "")) %>%
  group_by(gene = Gene, gene_name = Gene.name, tissue = Tissue) %>%
  summarise(present = any(detected), .groups = "drop") %>%
  filter(tissue %in% shared_tissues)

# ---- 4. combined (comparable presence, same rule both sides) -----------------
# Tag the source at bind time instead of reverse-engineering it from NA gene_name.
combined <- bind_rows(gtex = gtex, hpa = hpa, .id = "source") %>%
  select(source, gene, gene_name, tissue, present)



# That's it, this part above reproduced your script.


# New part
# ---- 5. presence summary + barplot (replaces the broken bottom plots) --------
presence_by_tissue <- combined %>%
  group_by(tissue, source) %>%
  summarise(pct_present = mean(present), n = n(), .groups = "drop")

# I don't want to feed you the ready code but also don't want to turture you. The template is ready, just insert the relevant pieces
# coming from presence_by_tissue
# this should give you the barplot to illustrate the difference between HPA and GTEx in terms of detecting proteins.
# Build it, check pituitary gland, your finding should hold

png(filename = "tissue_detection_plot.png", 
    width = 700, height = 500, units = "px", pointsize = 12)
p_presence <- ggplot(presence_by_tissue,
                     aes(x = pct_present,
                         y = reorder(tissue, pct_present),
                         fill = source )
                     ) +
  geom_col(position = position_dodge()) +
  labs(title = "Detection rate by tissue and platform",
       x = "fraction present", y = NULL, fill = NULL) 


print(p_presence)

dev.off()
# Try to save this in proc/pic with any format png(), tiff() whatever. play with resolutions and pixels.




# =============================================================================
# SUBCELLULAR LOCALISATION
# =============================================================================

# ---- 6. locations (gene-level). Main.location is ';'-separated. --------------

# same as in HPA IHC get rid of Uncertain entries in reliability + NA in Main.location
# create new column 'primary_location that takes first entry from Main.location. hint: trimws(sub()) 

cc <- cc_raw %>% 
  filter(!Reliability %in% "Uncertain", !is.na(Main.location), !Main.location %in% "") %>%
  mutate(
    primary_location = trimws(sub(";.*", "", Main.location)),
  ) %>%
  rename(gene = Gene)

# HPA presence + primary location, per gene x tissue (localisation is HPA-native)
# inner join with hpa by gene. you'll have merged dataset with detection + localisation

hpa_loc <- inner_join(hpa, cc, by = "gene")

# ---- 7. barplot: proteins per subcellular compartment ------------------------
# Build a barplot on loc primary location to see the distribution

location_distributed <- hpa_loc %>%
  group_by(primary_location) %>%
  summarise(loc_present = mean(present), n = n(), .groups = "drop")

png(filename = "proteins per subcellular compartment.png", width = 600, height = 600, units = "px")

loc_distroplot <- ggplot(location_distributed,
                     aes(x = n,
                         y = reorder(primary_location,n),
                         fill = n )) +
  geom_col(position = position_dodge()) +
  scale_x_log10() + 
  coord_flip() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  labs(title = "Proteins per subcellular compartment",
       x = "count", y = "Primary location")

print(loc_distroplot)

dev.off()


# ---- 8. barplot: detection rate by compartment (pooled tissues, HPA) ---------
# More interesting barplot that uses hpa_loc


png(filename = "detection rate by compartment.png", width = 600, height = 600, units = "px")

protein_comp <- ggplot(location_distributed,
                       aes(x = loc_present,
                           y = reorder(primary_location, loc_present),
                           fill = loc_present )) +
  geom_col(position = position_dodge()) +
  coord_flip() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  labs(title = "Detection rate by compartments",
       x = "fraction present", y = "Primary location", fill = NULL)

print(protein_comp)

dev.off()

# ---- 9. Fisher per tissue: is detection associated with localisation? --------
# One contingency table per tissue: primary location x present/absent.
# R x 2 with ~30 location levels -> use simulate.p.value (exact Fisher would hit
# the workspace limit), then BH/FDR across tissues.

tissue_table <- hpa_loc %>%
  group_by(tissue) %>%
  reframe(
    p_value = fisher.test(table(primary_location, present), simulate.p.value = TRUE
    )$p.value
  ) %>% 
  ungroup() %>%
  mutate(p_val_adjusted = p.adjust(p_value, method = "BH"))
  
write.csv(
  x = tissue_table,
  file = file.path("fisher_results.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
# =============================================================================
# optional: which compartment drives it (2x2 enrichment per tissue x location)
# The omnibus above says whether localisation matters in a tissue, this one says
# which compartments are over/under-detected.
# =============================================================================



# optinal: summary barplot: in how many tissues is each compartment significantly
# enriched vs depleted among detected proteins


