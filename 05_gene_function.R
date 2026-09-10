setwd("D:/C/Desktop/Leccino24/Landscape_156WWE")

library(fit.models)
library(adegenet)
library(vegan)
library(geosphere)
library(vcfR)
library(data.table)
library(dplyr)
library(LEA)
library(ggplot2)
library(ggpubr)
library(ggrepel)
library(data.table)

gff <- fread(
  "D:/D/vcf_file_GEA_leccino/Olea_europaea_cv_Leccino.gff3.gz",
  sep = "\t",
  header = FALSE,
  comment.char = "#"
)

library(data.table)

# ------------------------------------------------------------
# 1. Make sure GFF has proper column names
# ------------------------------------------------------------

colnames(gff) <- c(
  "seqid",
  "source",
  "type",
  "start",
  "end",
  "score",
  "strand",
  "phase",
  "attributes"
)

# ------------------------------------------------------------
# 2. Keep gene annotations
# ------------------------------------------------------------

genes <- gff[type == "gene"]

# Extract gene ID, name and functional annotation
genes[, gene_id := sub(".*ID=([^;]+).*", "\\1", attributes)]
genes[, gene_name := sub(".*Name=([^;]+).*", "\\1", attributes)]
genes[, note := sub(".*Note=([^;]+).*", "\\1", attributes)]


# ------------------------------------------------------------
# 3. Your SNP marker list
# ------------------------------------------------------------

snp_markers <- read.csv("D:/C/Desktop/Leccino24/Landscape_156WWE/GEA_bonferroni.csv")

snp_dt <- data.table(
  marker = snp_markers$X,
  pvalue = snp_markers$pvalue
)

# Remove "Response " prefix
snp_dt[, marker := sub("^Response ", "", marker)]

# Remove ".0.value" suffix
snp_dt[, marker := sub("\\.0\\.value$", "", marker)]

# Check
head(snp_dt)
# ============================================================
# 2. Extract chromosome, scaffold and SNP position
# ============================================================

# Extract chromosome number from GWHEUUU00000001
snp_dt[, chromosome := as.integer(
  sub("GWHEUUU0+", "", sub("_1_.*", "", marker))
)]

# Extract GFF scaffold ID
snp_dt[, seqid := sub("_1_.*", ".1", marker)]

# Extract SNP physical position
snp_dt[, pos := as.integer(
  sub(".*_1_", "", marker)
)]

# Check
head(snp_dt)

# ============================================================
# 3. Define ±100 kb window around each SNP
# ============================================================

window_size <- 100000

snp_dt[, `:=`(
  window_start = pmax(1, pos - window_size),
  window_end   = pos + window_size
)]


# ============================================================
# 4. Find genes within ±100 kb of each SNP
# ============================================================

# Create a separate SNP-window table
snp_windows <- snp_dt[
  ,
  .(
    marker,
    pvalue,
    chromosome,
    seqid,
    SNP_position = pos,
    window_start,
    window_end
  )
]

# Rename gene coordinates to avoid confusion
genes2 <- copy(genes)

genes2[
  ,
  `:=`(
    gene_start = start,
    gene_end = end
  )
]

# Perform overlap join
candidate_genes <- genes2[
  snp_windows,
  on = .(
    seqid,
    gene_start <= window_end,
    gene_end >= window_start
  ),
  nomatch = 0,
  allow.cartesian = TRUE
]




# ============================================================
# ABIOTIC-STRESS CANDIDATE GENE ANNOTATION
# Wild olive GEA SNPs ±100 kb
# ============================================================

library(data.table)
library(stringr)


# ============================================================
# 1. PARAMETERS
# ============================================================

window_size <- 100000L


# ============================================================
# 2. ABIOTIC-STRESS KEYWORDS
# ============================================================

adaptation_keywords <- c(
  
  # ----------------------------------------------------------
  # 1. GENERAL ABIOTIC STRESS
  # ----------------------------------------------------------
  "abiotic stress",
  "abiotic stimulus",
  "response to abiotic stimulus",
  "response to abiotic stress",
  "environmental stress",
  "environmental stimulus",
  "stress response",
  "stress tolerance",
  "stress resistance",
  
  # ----------------------------------------------------------
  # 2. TEMPERATURE
  # ----------------------------------------------------------
  "temperature",
  "temperature stress",
  "temperature response",
  "response to temperature",
  "heat stress",
  "heat response",
  "response to heat",
  "heat acclimation",
  "heat acclimatization",
  "cold stress",
  "cold response",
  "response to cold",
  "cold acclimation",
  "cold acclimatization",
  "high temperature",
  "low temperature",
  "thermal stress",
  "thermotolerance",
  "temperature homeostasis",
  
  # ----------------------------------------------------------
  # 3. DROUGHT / WATER
  # ----------------------------------------------------------
  "drought",
  "drought stress",
  "response to drought",
  "water stress",
  "water deprivation",
  "response to water deprivation",
  "water deficit",
  "water limitation",
  "water availability",
  "dehydration",
  "response to dehydration",
  "desiccation",
  "desiccation tolerance",
  "water homeostasis",
  "water balance",
  
  # ----------------------------------------------------------
  # 4. OSMOTIC STRESS
  # ----------------------------------------------------------
  "osmotic stress",
  "osmotic response",
  "response to osmotic stress",
  "osmotic adjustment",
  "osmoprotection",
  "osmotic tolerance",
  
  # ----------------------------------------------------------
  # 5. SALINITY / IONIC STRESS
  # ----------------------------------------------------------
  "salt stress",
  "salinity stress",
  "response to salt stress",
  "response to salinity",
  "salt tolerance",
  "salinity tolerance",
  "salt stress response",
  "ionic stress",
  "ion toxicity",
  
  # ----------------------------------------------------------
  # 6. OXIDATIVE STRESS / ROS
  # ----------------------------------------------------------
  "oxidative stress",
  "response to oxidative stress",
  "reactive oxygen",
  "reactive oxygen species",
  "hydrogen peroxide",
  "superoxide",
  "superoxide radical",
  "superoxide dismutase",
  "oxidative damage",
  "ROS",
  "ROS homeostasis",
  "reactive nitrogen",
  "nitric oxide",
  
  # ----------------------------------------------------------
  # 7. REDOX
  # ----------------------------------------------------------
  "redox homeostasis",
  "cell redox homeostasis",
  "cellular redox homeostasis",
  "redox regulation",
  "redox process",
  "redox balance",
  "redox signaling",
  "peroxiredoxin",
  "thioredoxin",
  "thioredoxin reductase",
  "glutaredoxin",
  "glutathione",
  "glutathione metabolism",
  
  # ----------------------------------------------------------
  # 8. ABA / HORMONAL ABIOTIC RESPONSE
  # ----------------------------------------------------------
  "abscisic acid",
  "response to abscisic acid",
  "abscisic acid signaling",
  "abscisic acid-activated signaling",
  "abscisic acid-activated",
  "ABA signaling",
  "response to hormone",
  "hormone signaling",
  
  # ----------------------------------------------------------
  # 9. OTHER STRESS-RELATED HORMONES
  # ----------------------------------------------------------
  "auxin",
  "ethylene",
  "jasmonic acid",
  "jasmonate",
  "salicylic acid",
  "brassinosteroid",
  "cytokinin",
  "gibberellin",
  
  # ----------------------------------------------------------
  # 10. LIGHT / PHOTOPERIOD / CIRCADIAN
  # ----------------------------------------------------------
  "photoperiod",
  "photoperiodism",
  "photoperiodic",
  "response to light",
  "light response",
  "response to light stimulus",
  "light stimulus",
  "light signaling",
  "circadian rhythm",
  "circadian regulation",
  "circadian clock",
  "biological rhythm",
  "seasonal",
  "seasonality",
  
  # ----------------------------------------------------------
  # 11. PHOTOSYNTHESIS / PHOTOPROTECTION
  # ----------------------------------------------------------
  "photosynthesis",
  "photosynthetic",
  "photosystem",
  "photosynthetic electron transport",
  "photoprotection",
  "photoinhibition",
  "response to high light",
  "high light",
  "ultraviolet",
  "UV response",
  "UV stress",
  
  # ----------------------------------------------------------
  # 12. ION HOMEOSTASIS
  # ----------------------------------------------------------
  "ion homeostasis",
  "cation homeostasis",
  "anion homeostasis",
  "metal ion homeostasis",
  "potassium ion homeostasis",
  "sodium ion homeostasis",
  "calcium ion homeostasis",
  "magnesium ion homeostasis",
  "chloride ion homeostasis",
  "transmembrane transporter activity",
  "potassium ion transport",
  "sodium ion transport",
  "calcium ion transport",
  "magnesium ion transport",
  "chloride ion transport",
  
  # ----------------------------------------------------------
  # 13. MEMBRANE / WATER TRANSPORT
  # ----------------------------------------------------------
  "membrane stability",
  "membrane organization",
  "membrane lipid",
  "lipid remodeling",
  "lipid homeostasis",
  "aquaporin",
  "water transport",
  "water channel activity",
  
  # ----------------------------------------------------------
  # 14. CELL WALL
  # ----------------------------------------------------------
  "cell wall modification",
  "cell wall organization",
  "cell wall remodeling",
  "cell wall integrity",
  "cell wall biogenesis",
  
  # ----------------------------------------------------------
  # 15. PROTEIN PROTECTION
  # ----------------------------------------------------------
  "protein folding",
  "protein refolding",
  "response to unfolded protein",
  "unfolded protein response",
  "heat shock protein",
  "heat shock",
  "molecular chaperone",
  "chaperone",
  "protein stability",
  
  # ----------------------------------------------------------
  # 16. PROTEIN DEGRADATION / QUALITY CONTROL
  # ----------------------------------------------------------
  "protein ubiquitination",
  "protein deubiquitination",
  "ubiquitin",
  "proteasome",
  "protein quality control",
  "protein degradation",
  
  # ----------------------------------------------------------
  # 17. DETOXIFICATION
  # ----------------------------------------------------------
  "detoxification",
  "xenobiotic",
  "xenobiotic detoxification",
  "glutathione",
  "glutathione transferase",
  "glutathione metabolism",
  
  # ----------------------------------------------------------
  # 18. NUTRIENT LIMITATION
  # ----------------------------------------------------------
  "nitrogen starvation",
  "nitrogen limitation",
  "phosphate starvation",
  "phosphate limitation",
  "nutrient starvation",
  "nutrient limitation",
  "nutrient deficiency",
  "nutrient homeostasis",
  
  # ----------------------------------------------------------
  # 19. CARBON / ENERGY STRESS
  # ----------------------------------------------------------
  "carbon starvation",
  "carbon limitation",
  "energy deprivation",
  "energy stress",
  "ATP homeostasis",
  
  # ----------------------------------------------------------
  # 20. METABOLIC ADJUSTMENT
  # ----------------------------------------------------------
  "osmolyte",
  "proline",
  "sugar",
  "trehalose",
  "sucrose",
  "starch metabolism",
  "carbohydrate metabolism",
  "secondary metabolite",
  
  # ----------------------------------------------------------
  # 21. SIGNALING
  # ----------------------------------------------------------
  "stress signaling",
  "stress signal",
  "calcium signaling",
  "calcium-mediated signaling",
  "protein kinase",
  "protein phosphorylation",
  "MAP kinase",
  "MAPK",
  "small GTPase",
  "signal transduction",
  
  # ----------------------------------------------------------
  # 22. TRANSCRIPTIONAL REGULATION
  # ----------------------------------------------------------
  "stress-responsive transcription",
  "stress response transcription",
  "transcription factor activity",
  "transcriptional regulation",
  
  # ----------------------------------------------------------
  # 23. CELLULAR PROTECTION / REPAIR
  # ----------------------------------------------------------
  "DNA repair",
  "DNA damage",
  "DNA damage response",
  "membrane repair",
  "cellular protection",
  
  # ----------------------------------------------------------
  # 24. AUTOPHAGY / CELLULAR RECYCLING
  # ----------------------------------------------------------
  "autophagy",
  "macroautophagy",
  "cellular recycling",
  
  # ----------------------------------------------------------
  # 25. SENESCENCE / CELL DEATH UNDER STRESS
  # ----------------------------------------------------------
  "stress-induced senescence",
  "senescence",
  "programmed cell death",
  "cell death",
  "apoptosis"
)




# ------------------------------------------------------------
# 1. Combine all keywords into ONE regex
# ------------------------------------------------------------
# ------------------------------------------------------------
# 1. Function to escape regex characters
# ------------------------------------------------------------

escape_regex <- function(x) {
  stringr::str_replace_all(
    x,
    "([.\\^$|()\\[\\]{}*+?\\\\])",
    "\\\\\\1"
  )
}
keyword_pattern <- paste(
  escape_regex(tolower(adaptation_keywords)),
  collapse = "|"
)

keyword_pattern <- paste0(
  "(?<![[:alnum:]_])(",
  keyword_pattern,
  ")(?![[:alnum:]_])"
)


# ------------------------------------------------------------
# 2. Find whether each annotation contains ANY keyword
# ------------------------------------------------------------

candidate_genes[
  ,
  adaptation_match := vapply(
    note,
    function(x) {
      
      if (is.na(x) || x == "") {
        return(NA_character_)
      }
      
      x <- tolower(x)
      
      matches <- str_extract_all(
        x,
        regex(keyword_pattern)
      )[[1]]
      
      if (length(matches) == 0) {
        return(NA_character_)
      }
      
      unique(matches) |>
        paste(collapse = "; ")
    },
    character(1)
  )
]


# ------------------------------------------------------------
# 3. Keep matched genes
# ------------------------------------------------------------

abiotic_genes <- candidate_genes[
  !is.na(adaptation_match)
]


# ------------------------------------------------------------
# 4. Fix gene coordinates
# ------------------------------------------------------------

abiotic_genes[
  ,
  `:=`(
    gene_start = start,
    gene_end   = end
  )
]


# ------------------------------------------------------------
# 5. Calculate true distance from SNP to gene
# ------------------------------------------------------------

abiotic_genes[
  ,
  distance_bp := fifelse(
    SNP_position < gene_start,
    gene_start - SNP_position,
    fifelse(
      SNP_position > gene_end,
      SNP_position - gene_end,
      0L
    )
  )
]


# ------------------------------------------------------------
# 6. SNP position relative to gene
# ------------------------------------------------------------

abiotic_genes[
  ,
  SNP_location := fifelse(
    SNP_position < gene_start,
    "upstream",
    fifelse(
      SNP_position > gene_end,
      "downstream",
      "inside_gene"
    )
  )
]


# ------------------------------------------------------------
# 7. Select final columns
# ------------------------------------------------------------

abiotic_genes <- abiotic_genes[
  ,
  .(
    marker,
    pvalue,
    chromosome,
    seqid,
    SNP_position,
    gene_id,
    gene_name,
    gene_start,
    gene_end,
    strand,
    distance_bp,
    SNP_location,
    note,
    adaptation_match
  )
]


# ------------------------------------------------------------
# 8. Remove duplicated SNP-gene combinations
# ------------------------------------------------------------

abiotic_genes <- unique(
  abiotic_genes,
  by = c("marker", "gene_id")
)


# ------------------------------------------------------------
# 9. Sort
# ------------------------------------------------------------

setorder(
  abiotic_genes,
  chromosome,
  SNP_position,
  distance_bp,
  gene_id
)


# ------------------------------------------------------------
# 10. Check
# ------------------------------------------------------------

head(
  abiotic_genes,
  20
)

nrow(abiotic_genes)
uniqueN(abiotic_genes$gene_id)

write.table(
  abiotic_genes,
  "genes_GEA_window100kb.csv",
  sep = ";",
  dec = ",",
  row.names = FALSE,
  col.names = TRUE,
  quote = TRUE,
  fileEncoding = "UTF-8"
)

a <- read.csv2(
  "genes_GEA_window100kb.csv",
  header = TRUE,
  stringsAsFactors = FALSE
)


library(writexl)

write_xlsx(
  abiotic_genes,
  "genes_GEA_window100kb.xlsx"
)
