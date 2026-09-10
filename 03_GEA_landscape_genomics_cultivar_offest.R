
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

#-------------
#enter vcf file
#--------------
geno155 <- read.vcfR("D:/D/vcf_file_GEA_leccino/WC156_lec24_DP10_100_miss090_ind085_mac1_MAF005.vcf.recode.vcf")#import vcf file
GI <- vcfR2genind(geno155)#transfrom file in genind object
geno155<-as.data.frame(GI)
geno155<-geno155%>% dplyr::select(ends_with(".0"))
#imputation
for (i in 1:ncol(geno155))
{
  geno155[which(is.na(geno155[,i])),i] <- median(geno155[-which(is.na(geno155[,i])),i], na.rm=TRUE)
}
geno155_data<- write.table(geno155, "geno_155.txt")
# faster input
geno155<-fread("D:/D/vcf_file_GEA_leccino/geno_155.txt")
geno155 <- as.data.frame(geno155)  # rownames only work on data.frame, not data.table
rownames(geno155) <- geno155[[1]]  # assign first column as row names
geno155 <- geno155[ , -1]          # remove the first column

##### Wilde East
listWE<-read.table("list_WE.txt")
genoWE<- geno155[rownames(geno155)%in% listWE$V1, ]

#---------------------
#Environlent Wild Weast
#---------------------

data_wild<- read.csv("Env_155_WWE.csv", header = TRUE)
data_wild <- data_wild %>%
  mutate(LAT_classes = cut(lat,
                           breaks = c(-Inf, 35, 40, 45),
                           labels = c("low_lat", "med_lat", "high_lat"),
                           right = FALSE))

# centering 
list142WW<-read.table("list142WW.txt")
dataWild_142W <- data_wild[data_wild$id%in% list142WW$V1, ]
test_env <- dataWild_142W[, c("bio2", "bio10", "bio11", "bio15", "bio18", "bio19", "clay", "N", "pH", "sand")]
Env <- scale(test_env, center=TRUE, scale=TRUE)
# Extract the centering values
env_center <- attr(Env, "scaled:center")
# Extract the scaling values
env_scale <- attr(Env, "scaled:scale")
#transform into dataset
Env <- as.data.frame(Env)
Variables_142WW<-data.frame(geno=dataWild_142W$id,group = dataWild_142W$group, region = dataWild_142W$region, lat_classes = dataWild_142W$LAT_classes, lat = dataWild_142W$lat, long = dataWild_142W$long,  Env )

genoWW<- geno155[rownames(geno155)%in% list142WW$V1, ]

#-------------
#maf filtering
#------------

Y <- genoWW
# Function to calculate MAF for each column (SNP)
calculate_maf <- function(geno_col) {
  geno_col <- na.omit(geno_col)
  allele_freq <- sum(geno_col) / (2 * length(geno_col))  # assumes diploid, genotypes 0/1/2
  maf <- min(allele_freq, 1 - allele_freq)
  return(maf)
}

# Apply function to each SNP (column)
maf_values <- apply(genoWW, 2, calculate_maf)
# Filter threshold, e.g., keep SNPs with MAF >= 0.05
maf_threshold <- 0.05
genoWW_maf <- genoWW[, maf_values >= maf_threshold]

write.table(genoWW_maf,"genoWW_maf.txt")

#-------------
#run LFMM GEA
#-------------

## Use latent factor for covariable correction
# latent factor temperature variable
Y <- genoWW_maf
Y <- as.matrix(genoWW_maf)

sel_latent<- data.frame(Variables_142WW%>% dplyr::select(bio2, bio10, bio11, bio15, bio18, bio19,clay, N, pH, sand))
write.env(sel_latent, "latent_all_variable.env")
X = read.table("latent_all_variable.env")
X <- as.matrix(X)

mod.lfmm2 <- lfmm2(input = Y, env = X, K = 3, effect.sizes = TRUE)
#get environment effect sizes
mod.lfmm2@B
#Define GEA
pv = lfmm2.test(mod.lfmm2, input = Y, env = X, full = T)
pvals <- pv$pvalue

hist(pv$pvalues)

# Estimate FDR-adjusted p-values
fdr_values <- p.adjust(pv$pvalues, method = "BH")
# Define FDR threshold
fdr_threshold <- 0.05

# Get indices of significant tests
signif_indices <- which(fdr_values < fdr_threshold)
GEA_lfmm <- data.frame(index = signif_indices, 
                       pvalue = pv$pvalues[signif_indices], 
                       fdr = fdr_values[signif_indices])

##Bonferroni threshold
thres <- 0.05/ncol(genoWW_maf)
signif_bonf <- which(pv$pvalues < thres)
GEA_bonferroni <- data.frame(index = signif_bonf, 
                             pvalue = pv$pvalues[signif_bonf])

PvaluesGEA_lfmm<-data.frame(pv$pvalues)

#define cadidate mod.lfmm2#define cadidate loci for GO 

GEA_lfmm <- data.frame(pvalue = pv$pvalues[-log10(pv$pvalue) >5])

write.csv(GEA_bonferroni, "GEA_bonferroni.csv")#selected 50 SNPs
write.csv(GEA_lfmm, "GEA_lfmm_all_var_log5.csv")#selected 255 SNPs
write.csv(pv$pvalues, "GEA_all_var_lfmm.csv")# all SNPs


#plotting Mhanattan plot using the library qqman

library(qqman)
Manhattan_all <- read.csv(file = "GEA_all_var_lfmm.csv", header=TRUE) #import the p value result for precipitation
jpeg(file = "/lustre/rocchettil/Manh_RDA_prec.jpeg")
manhattan(Manhattan_all, col = c("darkgreen", "gray60"),genomewideline = 5)

hist(Manhattan_all$P)
dev.off()
#-----------------
# Filter GEA
#-----------------
list_255<-read.table("listGEA_255.txt",  header=TRUE)
GEA_lfmm_255<-  genoWW_maf[, colnames(genoWW_maf)%in% list_255$SNP]
write.table(GEA_lfmm_255, "GEA_lfmm_all_var_255.txt")

GEA_lfmm_all_var<-read.table("GEA_lfmm_all_var_255.txt")


#----------------------------------------
#Cultivars genetic data and GEA filtering
#------------------------------------------

#upload genotypic file whole collection
geno_cultivar<- read.vcfR("D:/vcf_file_GEA_leccino/Cultivar_319_lec24_DP10_100_miss090_ind085_mac1.vcf.recode.vcf")


GI <- vcfR2genind(geno_cultivar)#transform file in genind object
geno_cultivar<-as.data.frame(GI)
geno_cultivar <- dplyr::select(geno_cultivar, ends_with(".0"))



GEA_lfmm_list<- read.table("GEA_lfmm_all_var_log5 (2).txt", header = T)
GEA_lfmm_all_var<-  genoWW_maf[, colnames(genoWW_maf)%in% GEA_lfmm_list$SNP]
write.table(GEA_lfmm_all_var, "GEA_lfmm_all_var.txt")
GEA_lfmm_all_var<-read.table("GEA_lfmm_all_var.txt")
list_bonf<-read.table("GEA_lfmm_all_var_Bonferroni.txt")
GEA_lfmm_bonf<-  GEA_lfmm_all_var[, colnames(GEA_lfmm_all_var)%in% list_bonf$x]

GEA <-colnames(GEA_lfmm_all_var)
GEA_geno_cultivar<-dplyr::select(geno_cultivar, all_of(GEA))

#imputation
for (i in 1:ncol(GEA_geno_cultivar))
{
  GEA_geno_cultivar[which(is.na(GEA_geno_cultivar[,i])),i] <- median(GEA_geno_cultivar[-which(is.na(GEA_geno_cultivar[,i])),i], na.rm=TRUE)
}

## save GEA all varible
write.table(GEA_geno_cultivar, "GEA_allWW_lfmm_all_cultivars.txt")


GEA_cultivars<-read.table("GEA_allWW_lfmm_all_cultivars.txt")

### filtering for MAF in cultivars

# Function to calculate MAF for each column (SNP)
calculate_maf <- function(geno_col) {
  geno_col <- na.omit(geno_col)
  allele_freq <- sum(geno_col) / (2 * length(geno_col))  # assumes diploid, genotypes 0/1/2
  maf <- min(allele_freq, 1 - allele_freq)
  return(maf)
}

# Apply function to each SNP (column)
maf_values <- apply(GEA_cultivars, 2, calculate_maf)
# Filter threshold, e.g., keep SNPs with MAF >= 0.05
maf_threshold <- 0.05
GEA_cultivars_maf <- GEA_cultivars[, maf_values >= maf_threshold]

GEA_124<- GEA_lfmm_all_var[, colnames(GEA_lfmm_all_var)%in% colnames(GEA_cultivars_maf)]
write.table(GEA_124,"GEA_124_WW.txt")
GEA_124<-read.table("GEA_124_WW.txt")

#-----------------------------------------
#enriched RDA 124 GEA polymorphic in wild and cultivars
#-----------------------------------------
RDA_all_enriched<-rda(GEA_124 ~ bio2 + bio10 + bio11 + bio15	+ bio18 + bio19 + clay + N+ pH + sand , Variables_142WW)
summary(eigenvals(RDA_all_enriched, model = "constrained"))
RsquareAdj(RDA_all_enriched)
plot(RDA_all_enriched)
sqrt(vif.cca(RDA_all_enriched))


# plot Geographic regions


TAB_gen <- data.frame(geno = row.names(scores(RDA_all_enriched , display = "sites")), scores(RDA_all_enriched, display = "sites", scaling = "sites"))

Geno <- merge(TAB_gen, Variables_142WW[, 1:5] ,by="geno")
TAB_var <- as.data.frame(scores(RDA_all_enriched, choices=c(1,2), display="bp"))
loading_geno_all_enriched_lat<-ggplot() +
  geom_hline(yintercept=0, linetype="dashed", color = gray(.80), size=0.6) +
  geom_vline(xintercept=0, linetype="dashed", color = gray(.80), size=0.6) +
  geom_point(data = Geno, aes(x=RDA1, y=RDA2, fill = lat), linewidth = 2.5, shape = 21, color = "black", stroke = 0.6, size = 3) +
  #scale_fill_manual(values = c("lightblue","darkgreen", "darkorange")) +
  scale_fill_gradientn(colors = c("#d73027","#fc8d59", "#fee090", "#91bfdb", "#4575b4"), 
                       name = "Latitude")+
  geom_segment(data = TAB_var, aes(xend = RDA1, yend = RDA2, x = 0, y = 0), 
               colour = "black", size = 0.15, linetype = 1, 
               arrow = arrow(length = unit(0.20, "cm"), type = "closed")) +
  geom_label_repel(data = TAB_var, aes(x = RDA1, y = RDA2, label = row.names(TAB_var)), 
                   size = 2, family = "Times") +
  xlab("RDA 1: 68%") + ylab("RDA 2: 11%") +
  guides(color=guide_legend(title="Latitude gradient")) +
  theme_bw(base_size = 12, base_family = "Times") +
  scale_x_continuous(breaks = seq(-5, 5, by = 0.5)) +
  scale_y_continuous(breaks = seq(-5, 5, by = 0.5)) +
  theme(
    axis.text = element_text(size = 7),
    axis.title = element_text(size = 7),# tick labels
    legend.title = element_text(size = 7),
    panel.background = element_blank(), legend.background = element_blank(), panel.grid = element_blank(), plot.background = element_blank(), legend.text=element_text(size=rel(.5)), strip.text = element_text(size=10))
#labs(title = "enriched RDA")
loading_geno_all_enriched_lat

loading_geno_all_enriched_region<-ggplot() +
  geom_hline(yintercept=0, linetype="dashed", color = gray(.80), size=0.6) +
  geom_vline(xintercept=0, linetype="dashed", color = gray(.80), size=0.6) +
  geom_point(data = Geno, aes(x=RDA1, y=RDA2, fill = region), linewidth = 2.5, shape = 21, color = "black", stroke = 0.8, size = 4) +
  scale_fill_manual(values = c("lightblue","darkgreen", "darkorange")) +
  geom_segment(data = TAB_var, aes(xend=RDA1, yend=RDA2, x=0, y=0), colour="black", linewidth =0.15, linetype=1, arrow=arrow(length = unit(0.02, "npc"))) +
  geom_label_repel(data = TAB_var, aes(x=RDA1, y=RDA2, label = row.names(TAB_var)), size = 3.2, family = "Times") +
  xlab("RDA 1: 68%") + ylab("RDA 2: 11%") +
  guides(color=guide_legend(title="Latitude gradient")) +
  theme_bw(base_size = 11, base_family = "Times") +
  theme(panel.background = element_blank(), legend.background = element_blank(), panel.grid = element_blank(), plot.background = element_blank(), legend.text=element_text(size=rel(.8)), strip.text = element_text(size=11))
#labs(title = "enriched RDA")
loading_geno_all_enriched_region
library(ggpubr)
combined_loading<-ggarrange(loading_geno_all_enriched_region, loading_geno_all_enriched_lat, nrow=1, ncol=2)


#Color-blind plot
loading_geno_all_enriched_lat <- ggplot() +
  # Dashed origin lines (updated 'size' to 'linewidth')
  geom_hline(yintercept = 0, linetype = "dashed", color = gray(0.80), linewidth = 0.5) +
  geom_vline(xintercept = 0, linetype = "dashed", color = gray(0.80), linewidth = 0.5) +
  
  # Data points with black outline
  geom_point(
    data = Geno, 
    aes(x = RDA1, y = RDA2, fill = lat), 
    shape = 21, 
    color = "black", 
    stroke = 0.5, 
    size = 3
  ) +
  
  # Colorblind-friendly perceptually uniform gradient
  scale_fill_distiller(
    palette = "PuOr", # Purple-Orange (High contrast, colorblind safe)
    name = "Latitude",
    direction = 1
  )+
  
  # Environmental vector arrows
  geom_segment(
    data = TAB_var, 
    aes(xend = RDA1, yend = RDA2, x = 0, y = 0), 
    colour = "black", 
    linewidth = 0.4, 
    linetype = 1, 
    arrow = arrow(length = unit(0.18, "cm"), type = "closed")
  ) +
  
  # Non-overlapping vector labels with subtle background box
  geom_label_repel(
    data = TAB_var, 
    aes(x = RDA1, y = RDA2, label = row.names(TAB_var)), 
    size = 2.5, 
    family = "Times",
    label.padding = unit(0.15, "lines"),
    box.padding = unit(0.2, "lines"),
    point.padding = unit(0.1, "lines"),
    segment.color = "grey30",
    segment.size = 0.3
  ) +
  
  # Axis labels
  labs(x = "RDA 1: 68%", y = "RDA 2: 11%") +
  
  # Coordinate breaks
  scale_x_continuous(breaks = seq(-5, 5, by = 0.5)) +
  scale_y_continuous(breaks = seq(-5, 5, by = 0.5)) +
  
  # Publication-ready minimalist theme
  theme_bw(base_size = 12, base_family = "Times") +
  theme(
    axis.text = element_text(size = 7, color = "black"),
    axis.title = element_text(size = 8, face = "bold"),
    legend.title = element_text(size = 8, face = "bold"),
    legend.text = element_text(size = 7),
    panel.background = element_blank(), 
    legend.background = element_blank(), 
    panel.grid = element_blank(), 
    plot.background = element_blank(),
    strip.text = element_text(size = 10)
  )

loading_geno_all_enriched_lat



ggsave(
  filename = "RDA_biplot_lat_CB.tiff",
  plot = loading_geno_all_enriched_lat,      # Optional if last plot was your desired one
  width = 3.2,             # In inches (default)
  height = 2.1,             # In inches
  units = "in",           # Can be "in", "cm", or "mm"
  dpi = 300               # Resolution (important for publications)
)

#Partial RDA
pRDA_all_enriched <- rda(
  GEA_124 ~ bio2 + bio10 + bio11 + bio15 + bio18 + bio19 +
    clay + N + pH + sand +
    Condition(lat + long),
  data = Variables_142WW
)

# Summary
summary(pRDA_all_enriched)

# Constrained eigenvalues
eigenvals(pRDA_all_enriched, model = "constrained")

# Adjusted R2
RsquareAdj(pRDA_all_enriched)

# Variance inflation factors
sqrt(vif.cca(pRDA_all_enriched))
anova.cca(pRDA_all_enriched, permutations = 999)

library(vegan)

# Environmental variables
ENV <- Variables_142WW[, c(
  "bio2", "bio10", "bio11", "bio15", "bio18", "bio19",
  "clay", "N", "pH", "sand"
)]

# Geographic variables
GEO <- Variables_142WW[, c("lat", "long")]

# Variance partitioning
vp <- varpart(
  GEA_124,
  ENV,
  GEO
)

vp

anova.cca(
  rda(
    GEA_124 ~ lat + long +
      Condition(
        bio2 + bio10 + bio11 + bio15 + bio18 + bio19 +
          clay + N + pH + sand
      ),
    data = Variables_142WW
  ),
  permutations = 999
)
RsquareAdj(rda(
  GEA_124 ~ lat + long +
    Condition(
      bio2 + bio10 + bio11 + bio15 + bio18 + bio19 +
        clay + N + pH + sand
    ),
  data = Variables_142WW
))
# ============================================================
# Variance partitioning table for enriched RDA
# 124 LFMM-derived GEA SNPs
# ============================================================

library(gridExtra)
library(grid)

# ------------------------------------------------------------
# Values directly from varpart()
# ------------------------------------------------------------

# Total inertia from varpart()
total_inertia <- 49.831

# Adjusted R2 fractions from varpart()
adjR2_env    <- 0.10992   # [a] ENV | GEO
adjR2_geo    <- 0.01908   # [b] GEO | ENV
adjR2_shared <- 0.22911   # [c] ENV ∩ GEO
adjR2_resid  <- 0.64189   # [d] residual

# Total adjusted R2 explained
adjR2_full <- 0.35811


# ------------------------------------------------------------
# Significance from partial RDA permutation tests
# ------------------------------------------------------------

p_env <- "0.001***"
p_geo <- "0.001***"


# ------------------------------------------------------------
# Proportion of total variance
#
# These are the variance-partitioning fractions reported
# by varpart() and sum to 1.
# ------------------------------------------------------------

prop_total_env    <- adjR2_env
prop_total_geo    <- adjR2_geo
prop_total_shared <- adjR2_shared
prop_total_resid  <- adjR2_resid


# ------------------------------------------------------------
# Proportion of explained variance
#
# Among the variance explained by ENV + GEO:
# ------------------------------------------------------------

prop_explained_env <- adjR2_env / adjR2_full
prop_explained_geo <- adjR2_geo / adjR2_full
prop_explained_shared <- adjR2_shared / adjR2_full


# ------------------------------------------------------------
# Create table
# ------------------------------------------------------------

partition_table <- data.frame(
  
  `Variance component` = c(
    "Full model: ENV + GEO model",
    "Pure environment:ENV | GEO",
    "Pure geography:GEO | ENV",
    "Shared environment–geography",
    "Residual"
  ),
  
  `Adjusted R²` = c(
    adjR2_full,
    adjR2_env,
    adjR2_geo,
    adjR2_shared,
    NA
  ),
  
  `p(>F)` = c(
    "0.001***",
    p_env,
    p_geo,
    "—",
    "—"
  ),
  
  `Proportion of explained variance` = c(
    "1.00",
    sprintf("%.2f", prop_explained_env),
    sprintf("%.2f", prop_explained_geo),
    sprintf("%.2f", prop_explained_shared),
    "—"
  ),
  
  `Proportion of total variance` = c(
    sprintf("%.2f", adjR2_full),
    sprintf("%.2f", prop_total_env),
    sprintf("%.2f", prop_total_geo),
    sprintf("%.2f", prop_total_shared),
    sprintf("%.2f", prop_total_resid)
  ),
  
  check.names = FALSE
)


# ------------------------------------------------------------
# Format adjusted R2
# ------------------------------------------------------------

partition_table$`Adjusted R²` <- ifelse(
  is.na(partition_table$`Adjusted R²`),
  "—",
  sprintf("%.3f", partition_table$`Adjusted R²`)
)


# ------------------------------------------------------------
# Table theme
# ------------------------------------------------------------

tab_theme <- ttheme_minimal(
  base_size = 12,
  base_family = "Times",
  
  colhead = list(
    fg_params = list(
      fontface = "bold",
      col = "black"
    ),
    bg_params = list(
      fill = "grey85",
      col = NA
    )
  ),
  
  core = list(
    fg_params = list(
      col = "black"
    ),
    bg_params = list(
      fill = "white",
      col = NA
    )
  )
)


# ------------------------------------------------------------
# Create table grob
# ------------------------------------------------------------

tbl <- tableGrob(
  partition_table,
  rows = NULL,
  theme = tab_theme
)


# ------------------------------------------------------------
# Column widths
# ------------------------------------------------------------

tbl$widths <- unit(
  c(
    6.4,   # Variance component
    3.4,   # Adjusted R2
    2.4,   # p value
    6.6,   # Proportion explained
    6    # Proportion total
  ),
  "cm"
)


# # ------------------------------------------------------------
# # Footnote
# # ------------------------------------------------------------
# 
# footnote <- textGrob(
#   paste0(
#     "*** p ≤ 0.001. Pure environmental variation represents ",
#     "the environmental component remaining after conditioning ",
#     "on latitude and longitude; pure geographic variation ",
#     "represents the geographic component remaining after ",
#     "conditioning on the environmental variables. ",
#     "Permutation tests based on 999 permutations."
#   ),
#   x = 0,
#   hjust = 0,
#   gp = gpar(
#     fontsize = 8,
#     fontfamily = "Times"
#   )
# )


# ------------------------------------------------------------
# Combine table + footnote
# ------------------------------------------------------------

tbl_final <- arrangeGrob(
  tbl,
  #footnote,
  ncol = 1,
  heights = unit.c(
    unit(1, "npc") - unit(0.65, "cm"),
    unit(0.65, "cm")
  )
)


# ------------------------------------------------------------
# Export TIFF
# ------------------------------------------------------------

tiff(
  "Supplementary_Variance_Partitioning_Enriched_RDA.tiff",
  width = 11,
  height = 3.2,
  units = "in",
  res = 600,
  compression = "lzw"
)

grid.newpage()
grid.draw(tbl_final)

dev.off()


#--------------------
#Spatial analysis
#--------------------

# enter spatial pixel values

library(raster)
library("readxl")


bio2<- raster(paste("D:/D/raster files/ENM_bioclim_soil_25/bio2_ENM_def_clip.tif"))
bio10<- raster(paste("D:/D/raster files/ENM_bioclim_soil_25/bio10_ENM_def_clip.tif"))
bio11<- raster(paste("D:/D/raster files/ENM_bioclim_soil_25/bio11_ENM_def_clip.tif"))
bio15<- raster(paste("D:/D/raster files/ENM_bioclim_soil_25/bio15_ENM_def_clip.tif"))
bio18<- raster(paste("D:/D/raster files/ENM_bioclim_soil_25/bio18_ENM_def_clip.tif"))
bio19<- raster(paste("D:/D/raster files/ENM_bioclim_soil_25/bio19_ENM_def_clip.tif"))
soilN<- raster(paste("D:/D/raster files/ENM_bioclim_soil_25/soilN_ENM_def_clip.tif"))
soilpH<- raster(paste("D:/D/raster files/ENM_bioclim_soil_25/soilpH_ENM_def_clip.tif"))

soilclay<- raster(paste("D:/D/raster files/ENM_bioclim_soil_25/soilclay_ENM_def_clip.tif"))
soilsand<- raster(paste("D:/D/raster files/ENM_bioclim_soil_25/soilsand_ENM_def_clip.tif"))

names(bio2) = 'bio2'
names(bio10) = 'bio10'
names(bio11) = 'bio11'
names(bio15) = 'bio15'
names(bio18) = 'bio18'
names(bio19) = 'bio19'
names(soilN ) = 'N'
names(soilpH) = 'pH'
names(soilclay) = 'clay'
names(soilsand) = 'sand'



#stack the different raster file
ras_current_var<-stack(c(bio2,bio10, bio11, bio15, bio18, bio19, soilclay,soilN,soilpH, soilsand))
plot(ras_current_var, 
     xlim = c(-10, 12), 
     ylim = c(27, 50))



pixel <- as.data.frame(rasterToPoints(ras_current_var))
pixel <- data.frame(x=pixel$x, y=pixel$y, bio2=pixel$bio2,bio10=pixel$bio10,bio11=pixel$bio11,bio15=pixel$bio15,bio18=pixel$bio18,bio19=pixel$bio19,clay=pixel$clay/10,N=pixel$N/100,pH=pixel$pH/10,sand=pixel$sand/10)
pixel<-na.omit(pixel)
pixel<- pixel[pixel$x>-10, ]
pixel_env<- pixel%>% dplyr::select(bio2, bio10, bio11, bio15, bio18, bio19,clay, N, pH, sand)

scaled_pixel <- scale(pixel_env, center = env_center, scale = env_scale)
scaled_pixel<-as.data.frame(scaled_pixel)


#prediction of pixel in the RDA space
scaled_pixel_LC <- predict(RDA_all_enriched, newdata=scaled_pixel, type="lc")
TAB_pixel_LC<- data.frame(lat = pixel$y, long = pixel$x, scaled_pixel_LC)
TAB_var <- as.data.frame(scores(RDA_all_enriched, choices=c(1,2), display="bp"))


### mapping with palette

library(ggplot2)
library(ggrepel)

# --- 1. Extract RDA axes ---
a1 <- TAB_pixel_LC$RDA1
a2 <- TAB_pixel_LC$RDA2

# --- 2. Calculate distance from origin ---
dist_origin <- sqrt(a1^2 + a2^2)
dist_norm <- (dist_origin - min(dist_origin)) / (max(dist_origin) - min(dist_origin))

# --- 3. Base color mapping ---
# Red for negative RDA1
red <- pmax(-a1, 0)
# Blue for positive RDA1
blue <- pmax(a1, 0)
# Green for positive RDA2
green <- pmax(a2, 0)

# Normalize channels
red <- (red - min(red)) / (max(red) - min(red) + 1e-6) * 255
green <- (green - min(green)) / (max(green) - min(green) + 1e-6) * 255
blue <- (blue - min(blue)) / (max(blue) - min(blue) + 1e-6) * 255

# --- 4. Grey factor (distance-based saturation) --- 
grey_factor <- dist_norm        # 0 = near origin, 1 = far
grey_target <- 200             # lighter grey target (closer to white) used 220
baseline_grey <- 0.3         # ensures center points are not black used 0.4

# Apply a power transformation for vibrancy
#red <- (red / 255) ^ 0.4 * 255
#green <- (green / 255) ^ 0.4 * 255
#blue <- (blue / 255) ^ 0.4 * 255


# Exaggerate color intensity 1.2
boost_factor <- 1.05
red   <- pmin(((red   / 255) ^ 0.5) * 255 * boost_factor, 255)
green <- pmin(((green / 255) ^ 0.5) * 255 * boost_factor, 255)
blue  <- pmin(((blue  / 255) ^ 0.5) * 255 * boost_factor, 255)

# Blend RGB toward grey
r_adj <- (1 - grey_factor) * ((1 - baseline_grey) * red   + baseline_grey * grey_target) +
  grey_factor * grey_target
g_adj <- (1 - grey_factor) * ((1 - baseline_grey) * green + baseline_grey * grey_target) +
  grey_factor * grey_target
b_adj <- (1 - grey_factor) * ((1 - baseline_grey) * blue  + baseline_grey * grey_target) +
  grey_factor * grey_target

# Final RGB colors
colors <- rgb(r_adj, g_adj, b_adj, maxColorValue = 255)

# --- 5. Plot ---
pp <- ggplot(as.data.frame(TAB_pixel_LC)) +
  geom_point(aes(x = RDA1, y = RDA2),
             color = colors,
             size = 2, shape = 21, fill = colors, stroke = 0.2) +
  geom_hline(yintercept = 0, linetype = "dashed", color = gray(0.80), size = 0.6) +
  geom_vline(xintercept = 0, linetype = "dashed", color = gray(0.80), size = 0.6) +
  geom_segment(data = TAB_var, aes(xend = RDA1, yend = RDA2, x = 0, y = 0), 
               colour = "black", size = 0.15, linetype = 1, 
               arrow = arrow(length = unit(0.1, "cm"), type = "closed")) +
  geom_label_repel(data = TAB_var, aes(x = RDA1, y = RDA2, label = row.names(TAB_var)), 
                   size = 1.7, family = "Times") +
  xlab("RDA 1: 68%") + 
  ylab("RDA 2: 11%") +
  scale_x_continuous(breaks = seq(-5, 5, by = 0.5)) +
  scale_y_continuous(breaks = seq(-5, 5, by = 0.5)) +
  theme_bw(base_size = 12, base_family = "Times") +
  theme(axis.text = element_text(size = 6),
        axis.title = element_text(size = 6),
        panel.background = element_blank(), 
        legend.background = element_blank(), 
        panel.grid = element_blank(), 
        plot.background = element_blank(), 
        legend.text = element_text(size = rel(0.5)), 
        strip.text = element_text(size = 10)) +
  scale_color_identity()

pp


#color blind version

library(ggplot2)
library(ggrepel)
# --- 1. Extract RDA axes ---
a1 <- TAB_pixel_LC$RDA1
a2 <- TAB_pixel_LC$RDA2
# --- 2. Calculate distance from origin ---
dist_origin <- sqrt(a1^2 + a2^2)
dist_norm <- (dist_origin - min(dist_origin)) / (max(dist_origin) - min(dist_origin))

# --- 3. Base color mapping (Red / Gold / Blue, no green) ---
w_red  <- pmax(-a1, 0)
w_gold <- pmax(a2, 0)
w_blue <- pmax(a1, 0)

w_red  <- (w_red  - min(w_red))  / (max(w_red)  - min(w_red)  + 1e-6)
w_gold <- (w_gold - min(w_gold)) / (max(w_gold) - min(w_gold) + 1e-6)
w_blue <- (w_blue - min(w_blue)) / (max(w_blue) - min(w_blue) + 1e-6)

red_rgb  <- as.vector(col2rgb("darkred"))  # try "#FF0000" for pure red instead
gold_rgb <- as.vector(col2rgb("#FFD700"))

red   <- w_red * red_rgb[1] + w_gold * gold_rgb[1] + w_blue * 0
green <- w_red * red_rgb[2] + w_gold * gold_rgb[2] + w_blue * 0
blue  <- w_red * red_rgb[3] + w_gold * gold_rgb[3] + w_blue * 255

# --- 4. Grey factor (distance-based saturation) — unchanged ---
grey_factor <- dist_norm
grey_target <- 200
baseline_grey <- 0.3

boost_factor <- 1.05
red   <- pmin(((red   / 255) ^ 0.5) * 255 * boost_factor, 255)
green <- pmin(((green / 255) ^ 0.5) * 255 * boost_factor, 255)
blue  <- pmin(((blue  / 255) ^ 0.5) * 255 * boost_factor, 255)

r_adj <- (1 - grey_factor) * ((1 - baseline_grey) * red   + baseline_grey * grey_target) +
  grey_factor * grey_target
g_adj <- (1 - grey_factor) * ((1 - baseline_grey) * green + baseline_grey * grey_target) +
  grey_factor * grey_target
b_adj <- (1 - grey_factor) * ((1 - baseline_grey) * blue  + baseline_grey * grey_target) +
  grey_factor * grey_target

colors <- rgb(r_adj, g_adj, b_adj, maxColorValue = 255)

# --- 5. Plot ---
pp <- ggplot(as.data.frame(TAB_pixel_LC)) +
  geom_point(aes(x = RDA1, y = RDA2),
             color = colors,
             size = 2, shape = 21, fill = colors, stroke = 0.2) +
  geom_hline(yintercept = 0, linetype = "dashed", color = gray(0.80), size = 0.6) +
  geom_vline(xintercept = 0, linetype = "dashed", color = gray(0.80), size = 0.6) +
  geom_segment(data = TAB_var, aes(xend = RDA1, yend = RDA2, x = 0, y = 0), 
               colour = "black", size = 0.15, linetype = 1, 
               arrow = arrow(length = unit(0.1, "cm"), type = "closed")) +
  geom_label_repel(data = TAB_var, aes(x = RDA1, y = RDA2, label = row.names(TAB_var)),
                   size = 1.7, family = "Times",
                   min.segment.length = Inf,   # never draw the leader line
                   max.overlaps = Inf,         # never drop a label for overlapping
                   box.padding = 0.1,          # allow labels to sit closer to points
                   point.padding = 0,
                   force = 0.5) +               # weaker repulsion → less movement, more overlap +
  xlab("RDA 1: 68%") + 
  ylab("RDA 2: 11%") +
  scale_x_continuous(breaks = seq(-5, 5, by = 0.5)) +
  scale_y_continuous(breaks = seq(-5, 5, by = 0.5)) +
  # Publication-ready minimalist theme
  theme_bw(base_size = 12, base_family = "Times") +
  theme(
    axis.text = element_text(size = 7, color = "black"),
    axis.title = element_text(size = 8, face = "bold"),
    legend.title = element_text(size = 8, face = "bold"),
    legend.text = element_text(size = 7),
    panel.background = element_blank(), 
    legend.background = element_blank(), 
    panel.grid = element_blank(), 
    plot.background = element_blank(),
    strip.text = element_text(size = 10)
  )+
  scale_color_identity()

pp






ggsave(
  filename = "biplot_landscape_CB.png",
  plot = pp,      # Optional if last plot was your desired one
  width = 2.3,             # In inches (default)
  height = 2.1,             # In inches
  units = "in",           # Can be "in", "cm", or "mm"
  dpi = 600               # Resolution (important for publications)
)



## plot in geographic map

library(sf)
library(rnaturalearth)
library(rnaturalearthdata)

# Load geographic boundaries of France, Spain, Morocco, Portugal, and Algeria
countries <- ne_countries(scale = "medium", country = c("France", "Spain", "Morocco", "Portugal", "Algeria"), returnclass = "sf")

# Remove French Guiana and Atlantic French territories
countries <- countries[!(countries$geounit %in% c("French Guiana", "Guadeloupe", "Martinique", "Saint Pierre and Miquelon", 
                                                  "Reunion", "Mayotte", "New Caledonia", "French Polynesia", 
                                                  "Wallis and Futuna", "Saint Barthelemy", "Saint Martin")), ]

# Convert TAB_pixel_LC to an sf object
TAB_pixel_LC_sf <- st_as_sf(TAB_pixel_LC,
                            coords = c("long", "lat"),
                            crs = 4326,
                            remove = FALSE)  # Keeps original long/lat columns
# Create the map
library(scales)

library(scales)

# Embed alpha directly into your colors
TAB_pixel_LC_sf$colors_alpha <- scales::alpha(colors, 1)

map <- ggplot(data = countries) +
  geom_sf(fill = "#EBEBEB", color = "black") +
  geom_point(data = TAB_pixel_LC_sf, 
             x = TAB_pixel_LC_sf$long, 
             y = TAB_pixel_LC_sf$lat, 
             color = TAB_pixel_LC_sf$colors_alpha,
             size = 0.05, show.legend = FALSE) +
  scale_color_identity() +
  coord_sf(xlim = c(-15, 15), ylim = c(28, 52), expand = FALSE,
  ) +
  theme_minimal() +
  # Publication-ready minimalist theme
  theme_bw(base_size = 12, base_family = "Times") +
  theme(
    axis.text = element_text(size = 7, color = "black"),
    axis.title = element_text(size = 8, face = "bold"),
    legend.title = element_text(size = 8, face = "bold"),
    legend.text = element_text(size = 7),
    panel.background = element_blank(), 
    legend.background = element_blank(), 
    panel.grid = element_blank(), 
    plot.background = element_blank(),
    strip.text = element_text(size = 10)
  )

map

library(ggplot2)
library(dplyr)
library(patchwork)


combined <- (
  (loading_geno_all_enriched_lat | pp) / (map | plot_spacer())
) +
  plot_layout(widths = c(1, 1), heights = c(1, 1))


ggsave(
  filename = "adaptive_landscape_color_CB.png",
  plot=map,
  dpi=600,
  width = 3,
  height =3,
  units = 'in'
)


ggsave(
  filename = "figure3.pdf",
  plot=combined,
  dpi=300,
  width = 6.2,
  height = 4,
  units = 'in'
)

#-----------------------------
#Estimation of cultivar offset
#------------------------------

#upload genotypic file whole collection
geno_cultivar<- read.vcfR("D:/vcf_file_GEA_leccino/Cultivar_319_lec24_DP10_100_miss090_ind085_mac1.vcf.recode.vcf")


GI <- vcfR2genind(geno_cultivar)#transform file in genind object
geno_cultivar<-as.data.frame(GI)
geno_cultivar <- dplyr::select(geno_cultivar, ends_with(".0"))



GEA_lfmm_list<- read.table("GEA_lfmm_all_var_log5 (2).txt", header = T)
GEA_lfmm_all_var<-  genoWW_maf[, colnames(genoWW_maf)%in% GEA_lfmm_list$SNP]
write.table(GEA_lfmm_all_var, "GEA_lfmm_all_var.txt")
GEA_lfmm_all_var<-read.table("GEA_lfmm_all_var.txt")
list_bonf<-read.table("GEA_lfmm_all_var_Bonferroni.txt")
GEA_lfmm_bonf<-  GEA_lfmm_all_var[, colnames(GEA_lfmm_all_var)%in% list_bonf$x]

GEA <-colnames(GEA_lfmm_all_var)
GEA_geno_cultivar<-dplyr::select(geno_cultivar, all_of(GEA))

#imputation
for (i in 1:ncol(GEA_geno_cultivar))
{
  GEA_geno_cultivar[which(is.na(GEA_geno_cultivar[,i])),i] <- median(GEA_geno_cultivar[-which(is.na(GEA_geno_cultivar[,i])),i], na.rm=TRUE)
}

## save GEA all varible
write.table(GEA_geno_cultivar, "GEA_allWW_lfmm_all_cultivars.txt")


GEA_cultivars<-read.table("GEA_allWW_lfmm_all_cultivars.txt")

### filtering for MAF in cultivars

# Function to calculate MAF for each column (SNP)
calculate_maf <- function(geno_col) {
  geno_col <- na.omit(geno_col)
  allele_freq <- sum(geno_col) / (2 * length(geno_col))  # assumes diploid, genotypes 0/1/2
  maf <- min(allele_freq, 1 - allele_freq)
  return(maf)
}

# Apply function to each SNP (column)
maf_values <- apply(GEA_cultivars, 2, calculate_maf)
# Filter threshold, e.g., keep SNPs with MAF >= 0.05
maf_threshold <- 0.05
GEA_cultivars_maf <- GEA_cultivars[, maf_values >= maf_threshold]







RDAscore_cul <- predict(RDA_all_enriched, newdata=GEA_cultivars_maf, type="wa")
RDAscore_cul<-as.data.frame(RDAscore_cul)
write.table(RDAscore_cul, "D:/C/Desktop/Leccino24/landscape_cultivar_offset/RDAscore_cul.txt")

### predictetion of WildEast based on genotypes

GEA <-colnames(GEA_124)
GEA_WE<-dplyr::select(genoWE, all_of(GEA))
RDAscore_WE <- predict(RDA_all_enriched, newdata=GEA_WE, type="wa", scaling = 2)


##predicte Wild West
RDAscore_WW <- predict(RDA_all_enriched, newdata=GEA_lfmm_all_var, type="wa")

anova.cca(RDA_all_enriched, permutations = 999)
anova.cca(RDA_all_enriched, by = "axis", permutations = 999)

###### plot predicted cultivar and spatial pixels in the RDA space

Tab_cultivar<- data.frame(ID = row.names(RDAscore_cul),RDAscore_cul[,1:2] )
Tab_cultivar$group<-"cultivars"
RDA_wild<-predict(RDA_all_enriched,newdata =  scaled_pixel, type = "lc")
Tab_wild<-data.frame(ID = TAB_pixel_LC$lat,RDA_wild[,1:2])
Tab_wild$group<-"spatial_point"


TAB_var <- as.data.frame(scores(RDA_all_enriched, choices=c(1,2), display="bp"))

wild_cult_pred<-rbind(Tab_wild, Tab_cultivar)
wild_cult_pred$group <- as.factor(wild_cult_pred$group)

arrow_scale <- 0.5  

ll <- ggplot() +
  geom_hline(yintercept = 0, linetype = "dashed", color = gray(0.8), size = 0.6) +
  geom_vline(xintercept = 0, linetype = "dashed", color = gray(0.8), size = 0.6) +
  geom_point(data = wild_cult_pred, aes(x = RDA1, y = RDA2, color = group), size = 2.5) +
  #scale_shape_manual(values = c(24,21))+
  scale_color_manual(values = c('#E69F00', "lightgrey")) +
  geom_segment(data = TAB_var, 
               aes(x = 0, y = 0, xend = RDA1 * arrow_scale, yend = RDA2 * arrow_scale),
               arrow = arrow(length = unit(0.15, "cm"), type = "closed"),
               color = "black", size = 0.15) +
  geom_label_repel(data = TAB_var, 
                   aes(x = RDA1 * arrow_scale, y = RDA2 * arrow_scale, label = row.names(TAB_var)),
                   size = 2, family = "Times") +
  xlab("RDA 1: 68%") + 
  ylab("RDA 2: 11%") +
  #labs(title = "RDA Cultivar and Spatial Pixel Predictions") +
  theme_bw(base_size = 11, base_family = "Times") +
  scale_x_continuous(breaks = seq(-5, 5, by = 0.2)) +
  scale_y_continuous(breaks = seq(-5, 5, by = 0.2)) +
  theme(
    axis.text = element_text(size = 7),
    axis.title = element_text(size = 7),# tick labels
    legend.title = element_text(size = 7),
    panel.background = element_blank(), legend.background = element_blank(), panel.grid = element_blank(), plot.background = element_blank(), legend.text=element_text(size=rel(.5)), strip.text = element_text(size=10))
ll

ggsave(
  filename = "biplot_CUL-spatial.png",
  plot = ll,      # Optional if last plot was your desired one
  width = 3.5,             # In inches (default)
  height = 2,             # In inches
  units = "in",           # Can be "in", "cm", or "mm"
  dpi = 300               # Resolution (important for publications)
)


#### max GO within the niche

coords <- Tab_wild[, 2:3]
n <- nrow(coords)
max_dist <- 0
pair <- c(NA, NA)

for (i in 1:(n - 1)) {
  for (j in (i + 1):n) {
    d <- sum((coords[i, ] - coords[j, ])^2)
    if (d > max_dist) {
      max_dist <- d
      pair <- c(i, j)
    }
  }
}
print(max_dist)
print(pair)

################################# map cultivar and wild in the RDA space
Tab_cultivar<- data.frame(geno = row.names(RDAscore_cul),RDAscore_cul[,1:2] )
Tab_cultivar$group<-"cultivars"

RDAscore_WE<-data.frame(geno = row.names(RDAscore_WE),RDAscore_WE[,1:2])
RDAscore_WE$group<-"genoWE"

RDA_wild<-predict(RDA_all_enriched,newdata =  GEA_124, type = "wa")
Tab_wild<-data.frame(geno = Variables_142WW[,1],RDA_wild[,1:2])
Tab_wild$group<-"wild"

WE<-data.frame(data_wild[143:155,7:16])
scaled_WE <- scale(WE, center = env_center, scale = env_scale)
scaled_WE<-as.data.frame(scaled_WE)
RDA_WE<-predict(RDA_all_enriched,newdata =  scaled_WE, type = "lc")
Tab_WE<-data.frame(geno = row.names(WE),RDA_WE[,1:2])
Tab_WE$group<-"WE"

TAB_var <- as.data.frame(scores(RDA_all_enriched, choices=c(1,2), display="bp"))

wild_cult_pred<-rbind(Tab_wild, Tab_cultivar)
wild_cult_pred$group <- as.factor(wild_cult_pred$group)

hh <- ggplot() +
  geom_hline(yintercept=0, linetype="dashed", color = gray(.80), size=0.6) +
  geom_vline(xintercept=0, linetype="dashed", color = gray(.80), size=0.6) +
  geom_point(data = wild_cult_pred, aes(x = RDA1, y = RDA2, fill = group, shape = group),size = 2, color = "black", stroke = 0.6)+
  scale_shape_manual(values = c(24,21,21))+
  scale_fill_manual(values=c('#E69F00',"grey48","lightblue"))+
  #scale_size_manual(values=c(3,3))+
  geom_segment(data = TAB_var, aes(xend=RDA1*0.3, yend=RDA2*0.3, x=0, y=0), colour="black", size=0.15, linetype=1, arrow = arrow(length=unit(0.15,"cm"),type = "closed")) +
  geom_label_repel(data = TAB_var, aes(x = RDA1*0.3, y = RDA2*0.3, label = row.names(TAB_var)),
                   size = 1.7, family = "Times",
                   min.segment.length = Inf,   # never draw the leader line
                   max.overlaps = Inf,         # never drop a label for overlapping
                   box.padding = 0.1,          # allow labels to sit closer to points
                   point.padding = 0,
                   force = 0.5) +          
  xlab("RDA 1: 68%") + ylab("RDA 2: 11%") +
  theme_bw(base_size = 11, base_family = "Times") +
  scale_x_continuous(breaks = seq(-5, 5, by = 0.2)) +
  scale_y_continuous(breaks = seq(-5, 5, by = 0.2)) +
  # Publication-ready minimalist theme
  theme_bw(base_size = 12, base_family = "Times") +
  theme(
    axis.text = element_text(size = 7, color = "black"),
    axis.title = element_text(size = 8, face = "bold"),
    legend.title = element_text(size = 8, face = "bold"),
    legend.text = element_text(size = 7),
    panel.background = element_blank(), 
    legend.background = element_blank(), 
    panel.grid = element_blank(), 
    plot.background = element_blank(),
    strip.text = element_text(size = 10)
  )
hh
library(ggpubr)
proj<-ggarrange(ll, hh,nrow=2, ncol=1)
ggsave(
  filename = "proj_cultivar.tif",
  plot=hh,
  dpi=600,
  width = 4,
  height = 2.5,
  units = 'in'
)

#----------------------------------
# Map specific cultivar offset
#----------------------------------

F <- GEA_cultivars_maf[rownames(GEA_cultivars) == "Picholine2", ]

FRDA <- predict(RDA_all_enriched, newdata=F, type="wa")
FRDA<-as.data.frame(FRDA)



TAB_pixel_LC$offset <- (FRDA$RDA1 - TAB_pixel_LC$RDA1)^2 * 0.684 + 
  (FRDA$RDA2 - TAB_pixel_LC$RDA2)^2 * 0.1095 + 
  (FRDA$RDA3 - TAB_pixel_LC$RDA3)^2 * 0.06071+
  (FRDA$RDA4 - TAB_pixel_LC$RDA4)^2 * 0.03739 
hist(TAB_pixel_LC$offset)


#TAB_pixel_LC$Frantoio_offset<-scale(TAB_pixel_LC$Frantoio_offset, center = intersection_x, scale = sd(results_df$GO))
#hist(TAB_pixel_LC$Frantoio_offset)

level1 <- quantile(TAB_pixel_LC$offset, probs = 0.10, na.rm = TRUE)
level2 <- quantile(TAB_pixel_LC$offset, probs = 0.20, na.rm = TRUE)
level3 <- quantile(TAB_pixel_LC$offset, probs = 0.30, na.rm = TRUE)
level4 <- quantile(TAB_pixel_LC$offset, probs = 0.40, na.rm = TRUE)
level5 <- quantile(TAB_pixel_LC$offset, probs = 0.50, na.rm = TRUE)
level6 <- quantile(TAB_pixel_LC$offset, probs = 0.60, na.rm = TRUE)
level7 <- quantile(TAB_pixel_LC$offset, probs = 0.70, na.rm = TRUE)
level8 <- quantile(TAB_pixel_LC$offset, probs = 0.80, na.rm = TRUE)
level9 <- quantile(TAB_pixel_LC$offset, probs = 0.90, na.rm = TRUE)



# Compute breaks for the column
sd_breaks <- c( min(TAB_pixel_LC$offset, na.rm = TRUE), level1, level2, level3, level4, level5, level6, level7, level8, level9, max(TAB_pixel_LC$offset, na.rm = TRUE))



# Create a color palette from blue to yellow
library(RColorBrewer)
color_palette <- brewer.pal(10, "PuOr")
# color_palette <- c(
#   "#004d00",  # very dark green
#   "#228B22",  # forest green
#   "#66C200",  # yellow-green
#   "#CCCC00",  # mustard yellow
#   "#FFD700",  # golden yellow
#   "#FFA500",  # orange
#   "#FF8C00",  # dark orange
#   "#FF4500",  # orange-red
#   "#B22222",  # firebrick
#   "#8B0000"   # dark red
# )

# Assign colors based on quantiles
TAB_pixel_LC$Foffset <- cut(TAB_pixel_LC$offset, breaks = sd_breaks, labels = color_palette)

library(ggplot2)
library(sf)
library(rnaturalearth)

# Load geographic boundaries of France, Spain, Morocco, Portugal, and Algeria
countries <- ne_countries(scale = "medium", country = c("France", "Spain", "Morocco", "Portugal", "Algeria"), returnclass = "sf")

# Remove French Guiana and Atlantic French territories
countries <- countries[!(countries$geounit %in% c("French Guiana", "Guadeloupe", "Martinique", "Saint Pierre and Miquelon", 
                                                  "Reunion", "Mayotte", "New Caledonia", "French Polynesia", 
                                                  "Wallis and Futuna", "Saint Barthelemy", "Saint Martin")), ]

# Convert TAB_pixel_LC to an sf object
TAB_pixel_LC_sf <- st_as_sf(TAB_pixel_LC, coords = c("long", "lat"), crs = 4326)

# --- Convert offset to percentile rank (0-100) ---
TAB_pixel_LC$offset_percentile <- rank(TAB_pixel_LC$offset, na.last = "keep") / 
  sum(!is.na(TAB_pixel_LC$offset)) * 100

# --- Plot with continuous colorblind-safe gradient + colorbar legend ---
P_map <- ggplot(data = countries) +
  geom_sf(fill = "#EBEBEB", color = "black") +
  geom_sf(data = TAB_pixel_LC_sf, aes(color = TAB_pixel_LC$offset_percentile), 
          size = 0.05, show.legend = TRUE) +
  scale_color_gradientn(
    colours = c("#2D004B", "#542788", "#8073AC", "#B2ABD2",
                "grey85", "grey85",
                "#FDB863", "#E08214", "#B35806", "#7F3B08"),
    values = c(0, 0.15, 0.30, 0.42,
               0.48, 0.52,
               0.58, 0.70, 0.85, 1),
    limits = c(0, 100),
    breaks = c(0, 25, 50, 75, 100),
    name = "Offset\nPercentile (%)",
    guide = guide_colorbar(
      barheight = unit(3, "cm"),
      barwidth  = unit(0.4, "cm"),
      ticks.colour = "black",
      frame.colour = "black"
    )
  )+
  coord_sf(xlim = c(-15, 15), ylim = c(28, 52), expand = FALSE) +
  theme_bw(base_size = 10) +
  labs(title = "Adaptive Landscape Picholine") +
  theme(
    panel.background = element_blank(),
    legend.position = "none",
    plot.title = element_text(size = 10)
    # legend.title = element_text(size = 8),
    # legend.text = element_text(size = 7)
  )






cultivar_offset<- ggarrange(Picholine2_map, Picual2_map, Manzanilla_map, Picholine_Marocaine2_map, nrow=2, ncol=2 )
cultivar_offset<- (Picholine2_map+Picual2_map)/ (Manzanilla_Cacerena_map + Picholine_Marocaine2_map)
ggsave(
  filename = " Picholine.tif",
  plot=  P_map,
  dpi=600,
  width = 3,
  height = 3,
  units = 'in'
)

#-------------------------
# scaled cultivar based on min max of F9
#---------------------------------------

# F9 reference
F9 <- GEA_124[rownames(GEA_124) == "OES_F9_10_S62_L004", ]

FRDA <- predict(RDA_all_enriched, newdata = F9, type = "wa")
FRDA <- as.data.frame(FRDA)

TAB_pixel_LC$offsetF9 <-
  (FRDA$RDA1 - TAB_pixel_LC$RDA1)^2 * 0.684 +
  (FRDA$RDA2 - TAB_pixel_LC$RDA2)^2 * 0.1095 +
  (FRDA$RDA3 - TAB_pixel_LC$RDA3)^2 * 0.06071 +
  (FRDA$RDA4 - TAB_pixel_LC$RDA4)^2 * 0.03739

# ------------------------------------------------------------
# Define quantile breaks from F9 (these are the reference scale
# that Aggezi_Akse1 will be projected onto)
# ------------------------------------------------------------

percentile_probs <- seq(0, 1, by = 0.10)          # 0, 0.10, ..., 1
sd_breaks <- quantile(TAB_pixel_LC$offsetF9, probs = percentile_probs, na.rm = TRUE)
percentile_labels <- percentile_probs * 100        # 0, 10, ..., 100

# ------------------------------------------------------------
# cultivar offset
# ------------------------------------------------------------

F <- GEA_cultivars_maf[
  rownames(GEA_cultivars_maf) == "Beladi-577",
]

FRDA <- predict(RDA_all_enriched, newdata = F, type = "wa")
FRDA <- as.data.frame(FRDA)

TAB_pixel_LC$offset <-
  (FRDA$RDA1 - TAB_pixel_LC$RDA1)^2 * 0.684 +
  (FRDA$RDA2 - TAB_pixel_LC$RDA2)^2 * 0.1095 +
  (FRDA$RDA3 - TAB_pixel_LC$RDA3)^2 * 0.06071 +
  (FRDA$RDA4 - TAB_pixel_LC$RDA4)^2 * 0.03739



TAB_pixel_LC$offset_percentile <- approx(
  x    = sd_breaks,
  y    = percentile_labels,
  xout = TAB_pixel_LC$offset,
  rule = 2
)$y

# ------------------------------------------------------------
# Convert to sf (only needed once)
# ------------------------------------------------------------

TAB_pixel_LC_sf <- st_as_sf(
  TAB_pixel_LC,
  coords = c("long", "lat"),
  crs = 4326
)

library(ggplot2)
library(sf)
library(rnaturalearth)

# Load geographic boundaries of France, Spain, Morocco, Portugal, and Algeria
countries <- ne_countries(
  scale = "medium",
  country = c("France", "Spain", "Morocco", "Portugal", "Algeria"),
  returnclass = "sf"
)

# Remove French Guiana and Atlantic French territories
countries <- countries[!(countries$geounit %in% c(
  "French Guiana", "Guadeloupe", "Martinique", "Saint Pierre and Miquelon",
  "Reunion", "Mayotte", "New Caledonia", "French Polynesia",
  "Wallis and Futuna", "Saint Barthelemy", "Saint Martin"
)), ]

# --- Plot with continuous colorblind-safe gradient + colorbar legend ---
P_map <- ggplot(data = countries) +
  geom_sf(fill = "#EBEBEB", color = "black") +
  geom_sf(
    data = TAB_pixel_LC_sf,
    aes(color = offset_percentile),
    size = 0.05, show.legend = TRUE
  ) +
  scale_color_gradientn(
    colours = c("#2D004B", "#542788", "#8073AC", "#B2ABD2",
                "grey85", "grey85",
                "#FDB863", "#E08214", "#B35806", "#7F3B08"),
    values = c(0, 0.15, 0.30, 0.42,
               0.48, 0.52,
               0.58, 0.70, 0.85, 1),
    limits = c(0, 100),
    breaks = c(0, 25, 50, 75, 100),
    name = "Offset\nPercentile (%)\n(F9 scale)",
    guide = guide_colorbar(
      barheight = unit(3, "cm"),
      barwidth  = unit(0.4, "cm"),
      ticks.colour = "black",
      frame.colour = "black"
    )
  ) +
  coord_sf(xlim = c(-15, 15), ylim = c(28, 52), expand = FALSE) +
  theme_bw(base_size = 10) +
  labs(title = "Adaptive Beladi") +
  theme(
    panel.background = element_blank(),
    legend.position = "none",
    plot.title = element_text(size = 10)
  )

ggsave(
  filename = "Beladi_scaled.tif",
  plot = P_map,
  dpi = 600,
  width = 3,
  height = 3,
  units = 'in'
)
















#--------------------
#Permutation
#----------------------

#----------------------------------
# Real Leccino offset (as you already have)
#----------------------------------
F <- GEA_cultivars_maf[rownames(GEA_cultivars) == "Leccino", ]
FRDA <- predict(RDA_all_enriched, newdata = F, type = "wa")
FRDA <- as.data.frame(FRDA)

TAB_pixel_LC$offset <- (FRDA$RDA1 - TAB_pixel_LC$RDA1)^2 * 0.684 +
  (FRDA$RDA2 - TAB_pixel_LC$RDA2)^2 * 0.1095 +
  (FRDA$RDA3 - TAB_pixel_LC$RDA3)^2 * 0.06071 +
  (FRDA$RDA4 - TAB_pixel_LC$RDA4)^2 * 0.03739

# Define the "adapted" zone using the same quantile cut you already use for mapping
level1 <- quantile(TAB_pixel_LC$offset, probs = 0.20, na.rm = TRUE)
adapted_idx <- which(TAB_pixel_LC$offset <= level1)     # dark-green pixels
adapted_pixels <- TAB_pixel_LC[adapted_idx, ]

real_mean_adapted <- mean(TAB_pixel_LC$offset[adapted_idx], na.rm = TRUE)
# (this will trivially be low, since the zone is defined from it — that's expected,
#  it's just the fixed target region for the null comparison below)

#--------------------
# Permutation restricted to the adapted zone
#--------------------
n_perm <- 1000
geno_cul <- GEA_cultivars_maf
null_mean_adapted <- numeric(n_perm)
null_offsets_adapted <- vector("list", n_perm)   # keep full null values in-zone if you want distribution overlay

for (i in 1:n_perm) {
  geno_cul_perm <- as.data.frame(
    lapply(geno_cul, function(col) sample(col))
  )
  rownames(geno_cul_perm) <- rownames(geno_cul)
  
  F_null <- geno_cul_perm[rownames(geno_cul_perm) == "Leccino", , drop = FALSE]
  FRDA_null <- predict(RDA_all_enriched, newdata = F_null, type = "wa")
  
  # only compute offset for the fixed adapted-zone pixels, not the whole landscape
  offset_null_adapted <- (FRDA_null[,1] - adapted_pixels$RDA1)^2 * 0.684 +
    (FRDA_null[,2] - adapted_pixels$RDA2)^2 * 0.1095 +
    (FRDA_null[,3] - adapted_pixels$RDA3)^2 * 0.06071 +
    (FRDA_null[,4] - adapted_pixels$RDA4)^2 * 0.03739
  
  null_offsets_adapted[[i]] <- offset_null_adapted
  null_mean_adapted[i] <- mean(offset_null_adapted, na.rm = TRUE)
}

#--------------------
# Compare
#--------------------
hist(null_mean_adapted, main = "Null: mean offset in Leccino's adapted zone",
     xlab = "mean offset")
abline(v = real_mean_adapted, col = "red", lwd = 2)

p_empirical <- mean(null_mean_adapted <= real_mean_adapted)
p_empirical

#--------------------
#Permutation
#----------------------
#----------------------------------
# Real Aggezi_Akse1 offset
#----------------------------------
F <- GEA_cultivars_maf[rownames(GEA_cultivars) == "Aggezi_Akse1", ]
FRDA <- predict(RDA_all_enriched, newdata = F, type = "wa")
FRDA <- as.data.frame(FRDA)

TAB_pixel_LC$offset <- (FRDA$RDA1 - TAB_pixel_LC$RDA1)^2 * 0.684 +
  (FRDA$RDA2 - TAB_pixel_LC$RDA2)^2 * 0.1095 +
  (FRDA$RDA3 - TAB_pixel_LC$RDA3)^2 * 0.06071 +
  (FRDA$RDA4 - TAB_pixel_LC$RDA4)^2 * 0.03739

# Define the "adapted" zone using the same quantile cut you already use for mapping
level1 <- quantile(TAB_pixel_LC$offset, probs = 0.20, na.rm = TRUE)
adapted_idx <- which(TAB_pixel_LC$offset <= level1)     # dark-green pixels
adapted_pixels <- TAB_pixel_LC[adapted_idx, ]

real_mean_adapted <- mean(TAB_pixel_LC$offset[adapted_idx], na.rm = TRUE)

#--------------------
# Permutation restricted to the adapted zone
#--------------------
n_perm <- 1000
geno_cul <- GEA_cultivars_maf
null_mean_adapted <- numeric(n_perm)
null_offsets_adapted <- vector("list", n_perm)

for (i in 1:n_perm) {
  geno_cul_perm <- as.data.frame(
    lapply(geno_cul, function(col) sample(col))
  )
  rownames(geno_cul_perm) <- rownames(geno_cul)
  
  F_null <- geno_cul_perm[rownames(geno_cul_perm) == "Aggezi_Akse1", , drop = FALSE]
  FRDA_null <- predict(RDA_all_enriched, newdata = F_null, type = "wa")
  
  offset_null_adapted <- (FRDA_null[,1] - adapted_pixels$RDA1)^2 * 0.684 +
    (FRDA_null[,2] - adapted_pixels$RDA2)^2 * 0.1095 +
    (FRDA_null[,3] - adapted_pixels$RDA3)^2 * 0.06071 +
    (FRDA_null[,4] - adapted_pixels$RDA4)^2 * 0.03739
  
  null_offsets_adapted[[i]] <- offset_null_adapted
  null_mean_adapted[i] <- mean(offset_null_adapted, na.rm = TRUE)
}

#--------------------
# Compare
#--------------------
hist(null_mean_adapted, main = "Null: mean offset in Aggezi_Akse1's adapted zone",
     xlab = "mean offset")
abline(v = real_mean_adapted, col = "red", lwd = 2)

p_empirical <- mean(null_mean_adapted <= real_mean_adapted)
p_empirical




# ============================================================
# PERMUTATION TEST: SPATIAL OVERLAP OF LOW-CGO REGIONS
# Leccino
#
# Question:
# Is the geographic distribution of Leccino's lowest-CGO
# environments more spatially specific than expected under
# random genotype permutations?
#
# Null hypothesis:
# A randomized genotype can generate a low-CGO geographic
# region with similar overlap to the observed Leccino region.
# ============================================================

set.seed(123)

# ------------------------------------------------------------
# 1. Real Leccino genomic position in RDA space
# ------------------------------------------------------------

F <- GEA_cultivars_maf[
  rownames(GEA_cultivars_maf) == "Leccino",
  ,
  drop = FALSE
]

FRDA <- predict(
  RDA_all_enriched,
  newdata = F,
  type = "wa"
)

FRDA <- as.data.frame(FRDA)


# ------------------------------------------------------------
# 2. Calculate observed Leccino CGO across the landscape
# ------------------------------------------------------------

TAB_pixel_LC$offset <- 
  (FRDA$RDA1 - TAB_pixel_LC$RDA1)^2 * 0.684 +
  (FRDA$RDA2 - TAB_pixel_LC$RDA2)^2 * 0.1095 +
  (FRDA$RDA3 - TAB_pixel_LC$RDA3)^2 * 0.06071 +
  (FRDA$RDA4 - TAB_pixel_LC$RDA4)^2 * 0.03739


# ------------------------------------------------------------
# 3. Define observed LOW-CGO region
#    Lowest 20% of Leccino CGO
# ------------------------------------------------------------

level20 <- quantile(
  TAB_pixel_LC$offset,
  probs = 0.20,
  na.rm = TRUE
)

observed_low_idx <- which(
  TAB_pixel_LC$offset <= level20
)

observed_low_pixels <- TAB_pixel_LC[
  observed_low_idx,
  ,
  drop = FALSE
]

N_low <- length(observed_low_idx)

cat("\nObserved Leccino low-CGO region\n")
cat("--------------------------------\n")
cat("20% CGO threshold:", level20, "\n")
cat("Number of low-CGO pixels:", N_low, "\n")


# ------------------------------------------------------------
# 4. Permutation
# ------------------------------------------------------------

n_perm <- 1000

geno_cul <- as.data.frame(
  GEA_cultivars_maf
)

# Store overlap statistics
null_overlap <- numeric(n_perm)

# Optional: store geographic centroid of random low-CGO regions
null_lat_centroid <- numeric(n_perm)
null_long_centroid <- numeric(n_perm)


# ------------------------------------------------------------
# 5. Permutation loop
# ------------------------------------------------------------

for (i in 1:n_perm) {
  
  # ----------------------------------------------------------
  # Randomize each SNP independently
  # ----------------------------------------------------------
  
  geno_cul_perm <- as.data.frame(
    lapply(
      geno_cul,
      function(col) sample(col)
    )
  )
  
  rownames(geno_cul_perm) <- rownames(geno_cul)
  
  
  # ----------------------------------------------------------
  # Extract randomized Leccino genotype
  # ----------------------------------------------------------
  
  F_null <- geno_cul_perm[
    rownames(geno_cul_perm) == "Leccino",
    ,
    drop = FALSE
  ]
  
  
  # ----------------------------------------------------------
  # Project randomized genotype into FIXED wild RDA
  # ----------------------------------------------------------
  
  FRDA_null <- predict(
    RDA_all_enriched,
    newdata = F_null,
    type = "wa"
  )
  
  
  # ----------------------------------------------------------
  # Calculate CGO across ALL landscape pixels
  # ----------------------------------------------------------
  
  offset_null <- 
    (FRDA_null[,1] - TAB_pixel_LC$RDA1)^2 * 0.684 +
    (FRDA_null[,2] - TAB_pixel_LC$RDA2)^2 * 0.1095 +
    (FRDA_null[,3] - TAB_pixel_LC$RDA3)^2 * 0.06071 +
    (FRDA_null[,4] - TAB_pixel_LC$RDA4)^2 * 0.03739
  
  
  # ----------------------------------------------------------
  # Define LOW-CGO region for this random genotype
  # ----------------------------------------------------------
  
  null_level20 <- quantile(
    offset_null,
    probs = 0.20,
    na.rm = TRUE
  )
  
  null_low_idx <- which(
    offset_null <= null_level20
  )
  
  
  # ----------------------------------------------------------
  # Calculate geographic overlap
  # ----------------------------------------------------------
  
  # Number of pixels shared between:
  #
  # observed Leccino low-CGO region
  #             AND
  # random low-CGO region
  
  shared_pixels <- length(
    intersect(
      observed_low_idx,
      null_low_idx
    )
  )
  
  
  # ----------------------------------------------------------
  # Overlap coefficient
  #
  # Both regions contain ~20% of the landscape,
  # therefore this measures the proportion of the
  # observed Leccino region reproduced by the random genotype.
  # ----------------------------------------------------------
  
  null_overlap[i] <- shared_pixels / N_low
  
  
  # ----------------------------------------------------------
  # Geographic centroid of random low-CGO region
  #
  # Change 'lat' and 'long' below if your coordinate
  # columns have different names.
  # ----------------------------------------------------------
  
  null_lat_centroid[i] <- mean(
    TAB_pixel_LC$lat[null_low_idx],
    na.rm = TRUE
  )
  
  null_long_centroid[i] <- mean(
    TAB_pixel_LC$long[null_low_idx],
    na.rm = TRUE
  )
}


# ============================================================
# 6. Observed geographic centroid
# ============================================================

observed_lat_centroid <- mean(
  observed_low_pixels$lat,
  na.rm = TRUE
)

observed_long_centroid <- mean(
  observed_low_pixels$long,
  na.rm = TRUE
)


# ============================================================
# 7. Null distribution of spatial overlap
# ============================================================

hist(
  null_overlap,
  breaks = 40,
  main = "Null distribution of spatial overlap",
  xlab = "Overlap with observed Leccino low-CGO region"
)

# There is no 'observed overlap' with itself to plot here.
# The important comparison is whether random overlap
# is systematically low/high relative to the observed
# geographic pattern.


# ============================================================
# 8. Null distribution of latitude centroid
# ============================================================

hist(
  null_lat_centroid,
  breaks = 40,
  main = "Leccino",
  xlab = "Mean latitude of low offest pixels"
)

abline(
  v = observed_lat_centroid,
  col = "red",
  lwd = 3
)
tiff("Leccino_null_lat_centroid.tif", width = 6, height = 5, units = "in", res = 600)

hist(
  null_lat_centroid,
  breaks = 40,
  main = "Leccino",
  xlab = "Mean latitude of low offset pixels",
  cex.lab = 1.4,
  cex.axis = 1.2,
  cex.main = 1.4
)

abline(
  v = observed_lat_centroid,
  col = "red",
  lwd = 3
)

dev.off()

# ============================================================
# 9. Summary of null distributions
# ============================================================

overlap_mean <- mean(
  null_overlap,
  na.rm = TRUE
)

overlap_sd <- sd(
  null_overlap,
  na.rm = TRUE
)

overlap_05 <- quantile(
  null_overlap,
  0.05,
  na.rm = TRUE
)

overlap_50 <- quantile(
  null_overlap,
  0.50,
  na.rm = TRUE
)

overlap_95 <- quantile(
  null_overlap,
  0.95,
  na.rm = TRUE
)


# ============================================================
# 10. Latitude statistics
# ============================================================

lat_null_mean <- mean(
  null_lat_centroid,
  na.rm = TRUE
)

lat_null_sd <- sd(
  null_lat_centroid,
  na.rm = TRUE
)

lat_null_05 <- quantile(
  null_lat_centroid,
  0.05,
  na.rm = TRUE
)

lat_null_95 <- quantile(
  null_lat_centroid,
  0.95,
  na.rm = TRUE
)


# ============================================================
# 11. Empirical P-value for geographic concentration
#
# If Leccino's low-CGO region is unusually NORTHERN,
# test whether random low-CGO regions have a centroid
# at least as far north as Leccino.
# ============================================================

p_lat_north <- (
  sum(
    null_lat_centroid >= observed_lat_centroid,
    na.rm = TRUE
  ) + 1
) / (
  n_perm + 1
)


# ============================================================
# 12. Z-score for latitude
# ============================================================

z_lat <- (
  observed_lat_centroid -
    lat_null_mean
) / lat_null_sd


# ============================================================
# 13. Results table
# ============================================================

permutation_spatial_Leccino <- data.frame(
  
  Cultivar = "Leccino",
  
  Low_CGO_threshold = level20,
  
  N_low_CGO_pixels = N_low,
  
  Observed_lat_centroid = observed_lat_centroid,
  
  Null_lat_mean = lat_null_mean,
  
  Null_lat_SD = lat_null_sd,
  
  Null_lat_05 = lat_null_05,
  
  Null_lat_95 = lat_null_95,
  
  Empirical_P_lat_north = p_lat_north,
  
  Z_lat = z_lat,
  
  Null_overlap_mean = overlap_mean,
  
  Null_overlap_SD = overlap_sd,
  
  Null_overlap_05 = overlap_05,
  
  Null_overlap_median = overlap_50,
  
  Null_overlap_95 = overlap_95
)


permutation_spatial_Leccino


# ============================================================
# PERMUTATION TEST: SPATIAL OVERLAP OF LOW-CGO REGIONS
# Aggezi_Akse1
#
# Question:
# Is the geographic distribution of Aggezi_Akse1's lowest-CGO
# environments different from that expected under random
# genotype permutations?
#
# Two complementary statistics are evaluated:
#
# 1. Spatial overlap:
#    How often does a randomized genotype reproduce the same
#    geographic low-CGO region as the observed cultivar?
#
# 2. Geographic centroid:
#    Is the observed low-CGO region unusually SOUTHERN
#    compared with randomized genotypes?
#
# The RDA model remains FIXED and was fitted using wild olive.
# ============================================================

set.seed(123)

# ============================================================
# 1. CULTIVAR
# ============================================================

cultivar_name <- "Aggezi_Akse1"


# ============================================================
# 2. Real cultivar genomic position in RDA space
# ============================================================

F <- GEA_cultivars_maf[
  rownames(GEA_cultivars_maf) == cultivar_name,
  ,
  drop = FALSE
]

FRDA <- predict(
  RDA_all_enriched,
  newdata = F,
  type = "wa"
)

FRDA <- as.data.frame(FRDA)


# ============================================================
# 3. Calculate observed CGO across the entire landscape
# ============================================================

TAB_pixel_LC$offset <- 
  (FRDA$RDA1 - TAB_pixel_LC$RDA1)^2 * 0.684 +
  (FRDA$RDA2 - TAB_pixel_LC$RDA2)^2 * 0.1095 +
  (FRDA$RDA3 - TAB_pixel_LC$RDA3)^2 * 0.06071 +
  (FRDA$RDA4 - TAB_pixel_LC$RDA4)^2 * 0.03739


# ============================================================
# 4. Define observed LOW-CGO region
#    Lowest 20% of CGO
# ============================================================

level20 <- quantile(
  TAB_pixel_LC$offset,
  probs = 0.20,
  na.rm = TRUE
)

observed_low_idx <- which(
  TAB_pixel_LC$offset <= level20
)

observed_low_pixels <- TAB_pixel_LC[
  observed_low_idx,
  ,
  drop = FALSE
]

N_low <- length(observed_low_idx)


cat("\nObserved", cultivar_name, "low-CGO region\n")
cat("---------------------------------------------\n")
cat("20% CGO threshold:", level20, "\n")
cat("Number of low-CGO pixels:", N_low, "\n")


# ============================================================
# 5. Permutation settings
# ============================================================

n_perm <- 1000

geno_cul <- as.data.frame(
  GEA_cultivars_maf
)

# Store overlap statistics
null_overlap <- numeric(n_perm)

# Store geographic centroid of random low-CGO regions
null_lat_centroid <- numeric(n_perm)
null_long_centroid <- numeric(n_perm)


# ============================================================
# 6. PERMUTATION LOOP
# ============================================================

for (i in 1:n_perm) {
  
  # ----------------------------------------------------------
  # Randomize each SNP independently
  # ----------------------------------------------------------
  
  geno_cul_perm <- as.data.frame(
    lapply(
      geno_cul,
      function(col) sample(col)
    )
  )
  
  rownames(geno_cul_perm) <- rownames(geno_cul)
  
  
  # ----------------------------------------------------------
  # Extract randomized cultivar genotype
  # ----------------------------------------------------------
  
  F_null <- geno_cul_perm[
    rownames(geno_cul_perm) == cultivar_name,
    ,
    drop = FALSE
  ]
  
  
  # ----------------------------------------------------------
  # Project randomized genotype into the FIXED wild RDA
  # ----------------------------------------------------------
  
  FRDA_null <- predict(
    RDA_all_enriched,
    newdata = F_null,
    type = "wa"
  )
  
  
  # ----------------------------------------------------------
  # Calculate CGO across ALL landscape pixels
  # ----------------------------------------------------------
  
  offset_null <- 
    (FRDA_null[,1] - TAB_pixel_LC$RDA1)^2 * 0.684 +
    (FRDA_null[,2] - TAB_pixel_LC$RDA2)^2 * 0.1095 +
    (FRDA_null[,3] - TAB_pixel_LC$RDA3)^2 * 0.06071 +
    (FRDA_null[,4] - TAB_pixel_LC$RDA4)^2 * 0.03739
  
  
  # ----------------------------------------------------------
  # Define LOW-CGO region for randomized genotype
  # ----------------------------------------------------------
  
  null_level20 <- quantile(
    offset_null,
    probs = 0.20,
    na.rm = TRUE
  )
  
  null_low_idx <- which(
    offset_null <= null_level20
  )
  
  
  # ==========================================================
  # 6A. SPATIAL OVERLAP
  # ==========================================================
  
  shared_pixels <- length(
    intersect(
      observed_low_idx,
      null_low_idx
    )
  )
  
  # Fraction of the observed cultivar's low-CGO region
  # reproduced by the randomized genotype
  
  null_overlap[i] <- shared_pixels / N_low
  
  
  # ==========================================================
  # 6B. GEOGRAPHIC CENTROID
  # ==========================================================
  
  null_lat_centroid[i] <- mean(
    TAB_pixel_LC$lat[null_low_idx],
    na.rm = TRUE
  )
  
  null_long_centroid[i] <- mean(
    TAB_pixel_LC$long[null_low_idx],
    na.rm = TRUE
  )
}


# ============================================================
# 7. OBSERVED GEOGRAPHIC CENTROID
# ============================================================

observed_lat_centroid <- mean(
  observed_low_pixels$lat,
  na.rm = TRUE
)

observed_long_centroid <- mean(
  observed_low_pixels$long,
  na.rm = TRUE
)


cat("\nObserved geographic centroid\n")
cat("----------------------------\n")
cat("Latitude:", observed_lat_centroid, "\n")
cat("Longitude:", observed_long_centroid, "\n")


# ============================================================
# 8. NULL DISTRIBUTION OF SPATIAL OVERLAP
# ============================================================

hist(
  null_overlap,
  breaks = 40,
  main = paste(
    "Null spatial overlap:",
    cultivar_name
  ),
  xlab = "Overlap with observed low-CGO region"
)


# ============================================================
# 9. NULL DISTRIBUTION OF LATITUDE
# ============================================================

hist(
  null_lat_centroid,
  breaks = 40,
  main = "Aggezi_Akse",
  xlab = "Mean latitude of low offest pixels",
  xlim = c(36, max(null_lat_centroid, observed_lat_centroid))
  
)

abline(
  v = observed_lat_centroid,
  col = "red",
  lwd = 3
)

tiff("Aggezi_null_lat_centroid.tif", width = 6, height = 5, units = "in", res = 600)

hist(
  null_lat_centroid,
  breaks = 40,
  main = "Aggezi_Akse",
  xlab = "Mean latitude of low offset pixels",
  cex.lab = 1.4,
  cex.axis = 1.2,
  cex.main = 1.4,
  xlim = c(36, max(null_lat_centroid, observed_lat_centroid))
)

abline(
  v = observed_lat_centroid,
  col = "red",
  lwd = 3
)

dev.off()
# ============================================================
# 10. NULL OVERLAP STATISTICS
# ============================================================

overlap_mean <- mean(
  null_overlap,
  na.rm = TRUE
)

overlap_sd <- sd(
  null_overlap,
  na.rm = TRUE
)

overlap_05 <- quantile(
  null_overlap,
  0.05,
  na.rm = TRUE
)

overlap_50 <- quantile(
  null_overlap,
  0.50,
  na.rm = TRUE
)

overlap_95 <- quantile(
  null_overlap,
  0.95,
  na.rm = TRUE
)


# ============================================================
# 11. NULL LATITUDE STATISTICS
# ============================================================

lat_null_mean <- mean(
  null_lat_centroid,
  na.rm = TRUE
)

lat_null_sd <- sd(
  null_lat_centroid,
  na.rm = TRUE
)

lat_null_05 <- quantile(
  null_lat_centroid,
  0.05,
  na.rm = TRUE
)

lat_null_50 <- quantile(
  null_lat_centroid,
  0.50,
  na.rm = TRUE
)

lat_null_95 <- quantile(
  null_lat_centroid,
  0.95,
  na.rm = TRUE
)


# ============================================================
# 12. EMPIRICAL P-VALUE: SOUTHERN SHIFT
#
# Aggezi_Akse1 is expected to have a southern low-CGO region.
#
# Test:
# How often does a randomized genotype generate a low-CGO
# centroid at least as far SOUTH as the observed Aggezi region?
# ============================================================

p_lat_south <- (
  sum(
    null_lat_centroid <= observed_lat_centroid,
    na.rm = TRUE
  ) + 1
) / (
  n_perm + 1
)


# ============================================================
# 13. Z-SCORE FOR LATITUDE
#
# Negative Z = observed region is SOUTH of null expectation
# Positive Z = observed region is NORTH of null expectation
# ============================================================

z_lat <- (
  observed_lat_centroid -
    lat_null_mean
) / lat_null_sd


# ============================================================
# 14. OPTIONAL: NORTH/SOUTH INTERPRETATION
# ============================================================

lat_shift <- (
  observed_lat_centroid -
    lat_null_mean
)

cat("\nGeographic shift\n")
cat("----------------\n")
cat(
  "Observed latitude - null mean:",
  lat_shift,
  "degrees\n"
)

if (lat_shift < 0) {
  
  cat(
    "Observed low-CGO region is SOUTH of the null expectation.\n"
  )
  
} else {
  
  cat(
    "Observed low-CGO region is NORTH of the null expectation.\n"
  )
}


# ============================================================
# 15. FINAL SUMMARY TABLE
# ============================================================

permutation_spatial_Aggezi <- data.frame(
  
  Cultivar = cultivar_name,
  
  Low_CGO_threshold = level20,
  
  N_low_CGO_pixels = N_low,
  
  Observed_lat_centroid = observed_lat_centroid,
  
  Observed_long_centroid = observed_long_centroid,
  
  Null_lat_mean = lat_null_mean,
  
  Null_lat_SD = lat_null_sd,
  
  Null_lat_05 = lat_null_05,
  
  Null_lat_median = lat_null_50,
  
  Null_lat_95 = lat_null_95,
  
  Latitude_shift = lat_shift,
  
  Empirical_P_lat_south = p_lat_south,
  
  Z_lat = z_lat,
  
  Null_overlap_mean = overlap_mean,
  
  Null_overlap_SD = overlap_sd,
  
  Null_overlap_05 = overlap_05,
  
  Null_overlap_median = overlap_50,
  
  Null_overlap_95 = overlap_95
)



# ============================================================
# PERMUTATION TEST: SPATIAL OVERLAP OF LOW-CGO REGIONS
# Beladi
#
# Question:
# Is the geographic distribution of Beladi's lowest-CGO
# environments more spatially specific than expected under
# random genotype permutations?
#
# Null hypothesis:
# A randomized genotype can generate a low-CGO geographic
# region with similar overlap to the observed Beladi region.
# ============================================================

set.seed(123)

# ------------------------------------------------------------
# 1. Real Beladi genomic position in RDA space
# ------------------------------------------------------------

F <- GEA_cultivars_maf[
  rownames(GEA_cultivars_maf) == "Beladi",
  ,
  drop = FALSE
]

FRDA <- predict(
  RDA_all_enriched,
  newdata = F,
  type = "wa"
)

FRDA <- as.data.frame(FRDA)


# ------------------------------------------------------------
# 2. Calculate observed Beladi CGO across the landscape
# ------------------------------------------------------------

TAB_pixel_LC$offset <- 
  (FRDA$RDA1 - TAB_pixel_LC$RDA1)^2 * 0.684 +
  (FRDA$RDA2 - TAB_pixel_LC$RDA2)^2 * 0.1095 +
  (FRDA$RDA3 - TAB_pixel_LC$RDA3)^2 * 0.06071 +
  (FRDA$RDA4 - TAB_pixel_LC$RDA4)^2 * 0.03739


# ------------------------------------------------------------
# 3. Define observed LOW-CGO region
#    Lowest 20% of Beladi CGO
# ------------------------------------------------------------

level20 <- quantile(
  TAB_pixel_LC$offset,
  probs = 0.20,
  na.rm = TRUE
)

observed_low_idx <- which(
  TAB_pixel_LC$offset <= level20
)

observed_low_pixels <- TAB_pixel_LC[
  observed_low_idx,
  ,
  drop = FALSE
]

N_low <- length(observed_low_idx)

cat("\nObserved Beladi low-CGO region\n")
cat("--------------------------------\n")
cat("20% CGO threshold:", level20, "\n")
cat("Number of low-CGO pixels:", N_low, "\n")


# ------------------------------------------------------------
# 4. Permutation
# ------------------------------------------------------------

n_perm <- 1000

geno_cul <- as.data.frame(
  GEA_cultivars_maf
)

# Store overlap statistics
null_overlap <- numeric(n_perm)

# Optional: store geographic centroid of random low-CGO regions
null_lat_centroid <- numeric(n_perm)
null_long_centroid <- numeric(n_perm)


# ------------------------------------------------------------
# 5. Permutation loop
# ------------------------------------------------------------

for (i in 1:n_perm) {
  
  # ----------------------------------------------------------
  # Randomize each SNP independently
  # ----------------------------------------------------------
  
  geno_cul_perm <- as.data.frame(
    lapply(
      geno_cul,
      function(col) sample(col)
    )
  )
  
  rownames(geno_cul_perm) <- rownames(geno_cul)
  
  
  # ----------------------------------------------------------
  # Extract randomized Beladi genotype
  # ----------------------------------------------------------
  
  F_null <- geno_cul_perm[
    rownames(geno_cul_perm) == "Beladi",
    ,
    drop = FALSE
  ]
  
  
  # ----------------------------------------------------------
  # Project randomized genotype into FIXED wild RDA
  # ----------------------------------------------------------
  
  FRDA_null <- predict(
    RDA_all_enriched,
    newdata = F_null,
    type = "wa"
  )
  
  
  # ----------------------------------------------------------
  # Calculate CGO across ALL landscape pixels
  # ----------------------------------------------------------
  
  offset_null <- 
    (FRDA_null[,1] - TAB_pixel_LC$RDA1)^2 * 0.684 +
    (FRDA_null[,2] - TAB_pixel_LC$RDA2)^2 * 0.1095 +
    (FRDA_null[,3] - TAB_pixel_LC$RDA3)^2 * 0.06071 +
    (FRDA_null[,4] - TAB_pixel_LC$RDA4)^2 * 0.03739
  
  
  # ----------------------------------------------------------
  # Define LOW-CGO region for this random genotype
  # ----------------------------------------------------------
  
  null_level20 <- quantile(
    offset_null,
    probs = 0.20,
    na.rm = TRUE
  )
  
  null_low_idx <- which(
    offset_null <= null_level20
  )
  
  
  # ----------------------------------------------------------
  # Calculate geographic overlap
  # ----------------------------------------------------------
  
  # Number of pixels shared between:
  #
  # observed Beladi low-CGO region
  #             AND
  # random low-CGO region
  
  shared_pixels <- length(
    intersect(
      observed_low_idx,
      null_low_idx
    )
  )
  
  
  # ----------------------------------------------------------
  # Overlap coefficient
  #
  # Both regions contain ~20% of the landscape,
  # therefore this measures the proportion of the
  # observed Beladi region reproduced by the random genotype.
  # ----------------------------------------------------------
  
  null_overlap[i] <- shared_pixels / N_low
  
  
  # ----------------------------------------------------------
  # Geographic centroid of random low-CGO region
  #
  # Change 'lat' and 'long' below if your coordinate
  # columns have different names.
  # ----------------------------------------------------------
  
  null_lat_centroid[i] <- mean(
    TAB_pixel_LC$lat[null_low_idx],
    na.rm = TRUE
  )
  
  null_long_centroid[i] <- mean(
    TAB_pixel_LC$long[null_low_idx],
    na.rm = TRUE
  )
}


# ============================================================
# 6. Observed geographic centroid
# ============================================================

observed_lat_centroid <- mean(
  observed_low_pixels$lat,
  na.rm = TRUE
)

observed_long_centroid <- mean(
  observed_low_pixels$long,
  na.rm = TRUE
)


# ============================================================
# 7. Null distribution of spatial overlap
# ============================================================

hist(
  null_overlap,
  breaks = 40,
  main = "Null distribution of spatial overlap",
  xlab = "Overlap with observed Beladi low-CGO region"
)

# There is no 'observed overlap' with itself to plot here.
# The important comparison is whether random overlap
# is systematically low/high relative to the observed
# geographic pattern.


# ============================================================
# 8. Null distribution of latitude centroid
# ============================================================

hist(
  null_lat_centroid,
  breaks = 40,
  main = "Beladi",
  xlab = "Mean latitude of low offest pixels"
)

abline(
  v = observed_lat_centroid,
  col = "red",
  lwd = 3
)


tiff("Beladi_null_lat_centroid.tif", width = 6, height = 5, units = "in", res = 600)

hist(
  null_lat_centroid,
  breaks = 40,
  main = "Beladi_Akse",
  xlab = "Mean latitude of low offset pixels",
  cex.lab = 1.4,
  cex.axis = 1.2,
  cex.main = 1.4,
  #xlim = c(36, max(null_lat_centroid, observed_lat_centroid))
)

abline(
  v = observed_lat_centroid,
  col = "red",
  lwd = 3
)

dev.off()

#============================================================
# PERMUTATION TEST: SPATIAL OVERLAP OF LOW-CGO REGIONS
# Razzaio1
#
# Question:
# Is the geographic distribution of Razzaio1's lowest-CGO
# environments more spatially specific than expected under
# random genotype permutations?
#
# Null hypothesis:
# A randomized genotype can generate a low-CGO geographic
# region with similar overlap to the observed Razzaio1 region.
# ============================================================

set.seed(123)

# ------------------------------------------------------------
# 1. Real Razzaio1 genomic position in RDA space
# ------------------------------------------------------------

F <- GEA_cultivars_maf[
  rownames(GEA_cultivars_maf) == "Razzaio1",
  ,
  drop = FALSE
]

FRDA <- predict(
  RDA_all_enriched,
  newdata = F,
  type = "wa"
)

FRDA <- as.data.frame(FRDA)


# ------------------------------------------------------------
# 2. Calculate observed Razzaio1 CGO across the landscape
# ------------------------------------------------------------

TAB_pixel_LC$offset <- 
  (FRDA$RDA1 - TAB_pixel_LC$RDA1)^2 * 0.684 +
  (FRDA$RDA2 - TAB_pixel_LC$RDA2)^2 * 0.1095 +
  (FRDA$RDA3 - TAB_pixel_LC$RDA3)^2 * 0.06071 +
  (FRDA$RDA4 - TAB_pixel_LC$RDA4)^2 * 0.03739


# ------------------------------------------------------------
# 3. Define observed LOW-CGO region
#    Lowest 20% of Razzaio1 CGO
# ------------------------------------------------------------

level20 <- quantile(
  TAB_pixel_LC$offset,
  probs = 0.20,
  na.rm = TRUE
)

observed_low_idx <- which(
  TAB_pixel_LC$offset <= level20
)

observed_low_pixels <- TAB_pixel_LC[
  observed_low_idx,
  ,
  drop = FALSE
]

N_low <- length(observed_low_idx)

cat("\nObserved Razzaio1 low-CGO region\n")
cat("--------------------------------\n")
cat("20% CGO threshold:", level20, "\n")
cat("Number of low-CGO pixels:", N_low, "\n")


# ------------------------------------------------------------
# 4. Permutation
# ------------------------------------------------------------

n_perm <- 1000

geno_cul <- as.data.frame(
  GEA_cultivars_maf
)

# Store overlap statistics
null_overlap <- numeric(n_perm)

# Optional: store geographic centroid of random low-CGO regions
null_lat_centroid <- numeric(n_perm)
null_long_centroid <- numeric(n_perm)


# ------------------------------------------------------------
# 5. Permutation loop
# ------------------------------------------------------------

for (i in 1:n_perm) {
  
  # ----------------------------------------------------------
  # Randomize each SNP independently
  # ----------------------------------------------------------
  
  geno_cul_perm <- as.data.frame(
    lapply(
      geno_cul,
      function(col) sample(col)
    )
  )
  
  rownames(geno_cul_perm) <- rownames(geno_cul)
  
  
  # ----------------------------------------------------------
  # Extract randomized Razzaio1 genotype
  # ----------------------------------------------------------
  
  F_null <- geno_cul_perm[
    rownames(geno_cul_perm) == "Razzaio1",
    ,
    drop = FALSE
  ]
  
  
  # ----------------------------------------------------------
  # Project randomized genotype into FIXED wild RDA
  # ----------------------------------------------------------
  
  FRDA_null <- predict(
    RDA_all_enriched,
    newdata = F_null,
    type = "wa"
  )
  
  
  # ----------------------------------------------------------
  # Calculate CGO across ALL landscape pixels
  # ----------------------------------------------------------
  
  offset_null <- 
    (FRDA_null[,1] - TAB_pixel_LC$RDA1)^2 * 0.684 +
    (FRDA_null[,2] - TAB_pixel_LC$RDA2)^2 * 0.1095 +
    (FRDA_null[,3] - TAB_pixel_LC$RDA3)^2 * 0.06071 +
    (FRDA_null[,4] - TAB_pixel_LC$RDA4)^2 * 0.03739
  
  
  # ----------------------------------------------------------
  # Define LOW-CGO region for this random genotype
  # ----------------------------------------------------------
  
  null_level20 <- quantile(
    offset_null,
    probs = 0.20,
    na.rm = TRUE
  )
  
  null_low_idx <- which(
    offset_null <= null_level20
  )
  
  
  # ----------------------------------------------------------
  # Calculate geographic overlap
  # ----------------------------------------------------------
  
  # Number of pixels shared between:
  #
  # observed Razzaio1 low-CGO region
  #             AND
  # random low-CGO region
  
  shared_pixels <- length(
    intersect(
      observed_low_idx,
      null_low_idx
    )
  )
  
  
  # ----------------------------------------------------------
  # Overlap coefficient
  #
  # Both regions contain ~20% of the landscape,
  # therefore this measures the proportion of the
  # observed Razzaio1 region reproduced by the random genotype.
  # ----------------------------------------------------------
  
  null_overlap[i] <- shared_pixels / N_low
  
  
  # ----------------------------------------------------------
  # Geographic centroid of random low-CGO region
  #
  # Change 'lat' and 'long' below if your coordinate
  # columns have different names.
  # ----------------------------------------------------------
  
  null_lat_centroid[i] <- mean(
    TAB_pixel_LC$lat[null_low_idx],
    na.rm = TRUE
  )
  
  null_long_centroid[i] <- mean(
    TAB_pixel_LC$long[null_low_idx],
    na.rm = TRUE
  )
}


# ============================================================
# 6. Observed geographic centroid
# ============================================================

observed_lat_centroid <- mean(
  observed_low_pixels$lat,
  na.rm = TRUE
)

observed_long_centroid <- mean(
  observed_low_pixels$long,
  na.rm = TRUE
)


# ============================================================
# 7. Null distribution of spatial overlap
# ============================================================

hist(
  null_overlap,
  breaks = 40,
  main = "Null distribution of spatial overlap",
  xlab = "Overlap with observed Razzaio1 low-CGO region"
)

# There is no 'observed overlap' with itself to plot here.
# The important comparison is whether random overlap
# is systematically low/high relative to the observed
# geographic pattern.


# ============================================================
# 8. Null distribution of latitude centroid
# ============================================================

hist(
  null_lat_centroid,
  breaks = 40,
  main = "Razzaio",
  xlab = "Mean latitude of low offest pixels",
  xlim = c(38, 42)
)

abline(
  v = observed_lat_centroid,
  col = "red",
  lwd = 3
)

tiff("Razzaio_null_lat_centroid.tif", width = 6, height = 5, units = "in", res = 600)

hist(
  null_lat_centroid,
  breaks = 40,
  main = "Razzaio",
  xlab = "Mean latitude of low offset pixels",
  cex.lab = 1.4,
  cex.axis = 1.2,
  cex.main = 1.4,
  xlim = c(38, 42)
)

abline(
  v = observed_lat_centroid,
  col = "red",
  lwd = 3
)

dev.off()


p_lat_north <- (
  sum(
    null_lat_centroid >= observed_lat_centroid,
    na.rm = TRUE
  ) + 1
) / (
  n_perm + 1
)






















#---------------------------------
#Cultivar Spatial GO for Rshinyapp
#---------------------------------

# Create empty list to store each cultivar's offset vector
offset_results <- list()

# Loop over all cultivars
for (cultivar in rownames(GEA_cultivars_maf)) {
  
  # Extract SNP data for that cultivar
  F <- GEA_cultivars_maf[rownames(GEA_cultivars_maf) == cultivar, ]
  
  # Predict RDA scores for the cultivar
  FRDA <- predict(RDA_all_enriched, newdata = F, type = "wa")
  FRDA <- as.data.frame(FRDA)
  
  # Compute offset for each pixel
  offset_results[[cultivar]] <- (FRDA$RDA1 - TAB_pixel_LC$RDA1)^2 * 0.684 + 
    (FRDA$RDA2 - TAB_pixel_LC$RDA2)^2 * 0.1095 + 
    (FRDA$RDA3 - TAB_pixel_LC$RDA3)^2 * 0.06071 +
    (FRDA$RDA4 - TAB_pixel_LC$RDA4)^2 * 0.03739
}

# Combine results into one dataframe: rows = TAB_pixel_LC (pixels), columns = cultivars
offset_df <- as.data.frame(offset_results)


offset_df <- data.frame(lat = TAB_pixel_LC$lat, long = TAB_pixel_LC$long, offset_df)


# Write to CSV
write.csv(offset_df, "D:/C/Desktop/Olive_GO_paper/Spatial_cultivar_GO.csv", row.names = FALSE)

cat("✅ File written: offset_by_pixel_and_cultivar.csv\n")

#-------------------------------------------
#Phenotypic Evaluation
#---------------------------------------


# Define the target latitude and longitude
target_lat <- 31.816095


target_long <- -7.60

# Define a tolerance value (small range of acceptable difference)
tolerance <- 1e-2

# Filter the data using a tolerance for matching
Marrakech_row <- TAB_pixel_LC %>%
  filter(abs(lat - target_lat) < tolerance & abs(long - target_long) < tolerance)

Marrakech_row<-as.data.frame(Marrakech_row)
### offset at marakesh


Marrakech_row <- Marrakech_row[1, ]

# Cultivar offset in Marrakesh
RDAscore_cul$offsetM <- (Marrakech_row$RDA1 - RDAscore_cul$RDA1)^2 * 0.684 + 
  (Marrakech_row$RDA2 - RDAscore_cul$RDA2)^2 * 0.1095 + 
  (Marrakech_row$RDA3 - RDAscore_cul$RDA3)^2 * 0.06071 + 
  (Marrakech_row$RDA4 - RDAscore_cul$RDA4)^2 * 0.03739


GoMar<-write.table(RDAscore_cul, "Go_Mar.txt")



## Wild offset in Marrakesh
RDA_wild<-as.data.frame(RDA_wild)

RDA_wild$offsetM <- (Marrakech_row$RDA1 - RDA_wild$RDA1)^2 * 0.684 + 
  (Marrakech_row$RDA2 - RDA_wild$RDA2)^2 * 0.1095 + 
  (Marrakech_row$RDA3 - RDA_wild$RDA3)^2 * 0.06071 + 
  (Marrakech_row$RDA4 - RDA_wild$RDA4)^2 * 0.03739


GoMar<-read.csv("GO_marakesh.csv")
GoMar <- na.omit(GoMar)

boxplot(GO ~ class, data = GoMar) 
library(ggplot2)

# Set factor levels in the desired order
GoMar$class <- factor(GoMar$class, levels = c("early_blooming", "mid_blooming", "late_blooming"))

# Plot
a<-ggplot(GoMar, aes(x = class, y = GO, fill = class)) +
  geom_boxplot() +
  scale_fill_manual(values = c("purple", "lightgrey", "darkblue")) +
  theme_minimal() +
  labs(title = "GO by Class",
       x = "Class",
       y = "GO") +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )
a
model <- lm(logGO ~ class, data = GoMar)

summary(model)

GoMar$class <- as.factor(GoMar$class)
model <- lm(GO ~ class, data = GoMar)

library(multcomp)
summary(glht(model, linfct = mcp(class = "Tukey")))

hist(GoMar$BLUE_FFD)
abline(v = c(120, 126), col = "red", lwd = 2, lty = 2)
b<-ggplot(GoMar, aes(x = FFD)) +
  geom_histogram(binwidth = 1, fill = "grey", color = "black") +
  #geom_vline(xintercept = c(120, 126), color = "red", linetype = "dashed", size = 1) +
  theme_minimal() +
  labs(title = "FFD distribution",
       x = "FFD",
       y = "Count") +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 12)
  )

b
plot(logGO ~ FFD, data = GoMar)

model <- lm(GO ~ BLUE_FFD, data = GoMar)
#abline(model, col = "red", lwd = 2)
summary(model)

c <- ggplot(GoMar, aes(x = BLUE_FFD, y = GO)) +
  
  geom_point(
    shape = 21,
    size = 3.5,
    stroke = 0.5,
    colour = "black",
    fill = "grey70",
    alpha = 0.9
  ) +
  
  geom_smooth(
    method = "lm",
    colour = "darkred",
    fill = "grey85",
    linewidth = 1,
    se = TRUE
  ) +
  
  annotate(
    "text",
    x = -Inf,
    y = Inf,
    hjust = -0.1,
    vjust = 1.2,
    size = 5,
    label = "R² = 0.04\np = 0.002"
  ) +
  
  labs(
    x = "BLUE full flowering days (FFD)",
    y = "Genomic Offset",
    title = NULL
  ) +
  
  theme_classic(base_size = 14) +
  
  theme(
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 13, colour = "black"),
    axis.line = element_line(linewidth = 0.6),
    axis.ticks = element_line(linewidth = 0.6)
  )

c

library(ggpubr)
d<-ggarrange(b,c, nrow = 1, ncol=2)
d
ggsave(
  filename = "GO_vs_BLUEFFD.tiff",
  plot = c,
  width = 4.5,
  height = 5,
  units = "in",
  dpi = 300,
  compression = "lzw"
)
ggsave("GO_vs_BLUEFFD.tiff", plot = c, width = 5, height = 5, dpi = 300)
ggsave("pheno_ev.jpg", plot = d, width = 6, height = 3, dpi = 300)








## plot in geographic map

library(sf)
library(rnaturalearth)
library(rnaturalearthdata)

# Load geographic boundaries of France, Spain, Morocco, Portugal, and Algeria
countries <- ne_countries(scale = "medium", country = c("France", "Spain", "Morocco", "Portugal", "Algeria"), returnclass = "sf")

# Remove French Guiana and Atlantic French territories
countries <- countries[!(countries$geounit %in% c("French Guiana", "Guadeloupe", "Martinique", "Saint Pierre and Miquelon", 
                                                  "Reunion", "Mayotte", "New Caledonia", "French Polynesia", 
                                                  "Wallis and Futuna", "Saint Barthelemy", "Saint Martin")), ]


# Define the specific latitude and longitude
latitude <- 31.80403
longitude <- -7.60 

# Create a data frame with the point
highlight_point <- data.frame(
  lon = longitude,
  lat = latitude
)

# Convert to sf object
highlight_point_sf <- st_as_sf(highlight_point, coords = c("lon", "lat"), crs = 4326)


# Create the map with continuous legend
map <- ggplot(data = countries) +
  geom_sf(fill = "#EBEBEB", color = "black") +
  geom_sf(data = highlight_point_sf, color = "red", shape = 17, size = 4)+  # shape = 17 for triangle
  coord_sf(xlim = c(-15, 15), ylim = c(28, 52), expand = FALSE) +
  theme_minimal() +
  labs(title = "WOGBM collection sites") +
  theme(
    panel.background = element_blank(),
    legend.position = "right",
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 8)
  )
map


### Spatial current offset Marrakesh

TAB_pixel_LC$offset <- (Marrakech_row$RDA1 - TAB_pixel_LC$RDA1)^2 * 0.684 + 
  (Marrakech_row$RDA2 - TAB_pixel_LC$RDA2)^2 * 0.1095 + 
  (Marrakech_row$RDA3 - TAB_pixel_LC$RDA3)^2 * 0.06071+
  (Marrakech_row$RDA4 - TAB_pixel_LC$RDA4)^2 * 0.03739 
hist(TAB_pixel_LC$offset)


#TAB_pixel_LC$Frantoio_offset<-scale(TAB_pixel_LC$Frantoio_offset, center = intersection_x, scale = sd(results_df$GO))
#hist(TAB_pixel_LC$Frantoio_offset)

level1 <- quantile(TAB_pixel_LC$offset, probs = 0.10, na.rm = TRUE)
level2 <- quantile(TAB_pixel_LC$offset, probs = 0.20, na.rm = TRUE)
level3 <- quantile(TAB_pixel_LC$offset, probs = 0.30, na.rm = TRUE)
level4 <- quantile(TAB_pixel_LC$offset, probs = 0.40, na.rm = TRUE)
level5 <- quantile(TAB_pixel_LC$offset, probs = 0.50, na.rm = TRUE)
level6 <- quantile(TAB_pixel_LC$offset, probs = 0.60, na.rm = TRUE)
level7 <- quantile(TAB_pixel_LC$offset, probs = 0.70, na.rm = TRUE)
level8 <- quantile(TAB_pixel_LC$offset, probs = 0.80, na.rm = TRUE)
level9 <- quantile(TAB_pixel_LC$offset, probs = 0.90, na.rm = TRUE)



# Compute breaks for the column
sd_breaks <- c( min(TAB_pixel_LC$offset, na.rm = TRUE), level1, level2, level3, level4, level5, level6, level7, level8, level9, max(TAB_pixel_LC$offset, na.rm = TRUE))



# Create a color palette from blue to yellow

color_palette <- c(
  "#004d00",  # very dark green
  "#228B22",  # forest green
  "#66C200",  # yellow-green
  "#CCCC00",  # mustard yellow
  "#FFD700",  # golden yellow
  "#FFA500",  # orange
  "#FF8C00",  # dark orange
  "#FF4500",  # orange-red
  "#B22222",  # firebrick
  "#8B0000"   # dark red
)

# Assign colors based on quantiles
TAB_pixel_LC$Foffset <- cut(TAB_pixel_LC$offset, breaks = sd_breaks, labels = color_palette)

library(ggplot2)
library(sf)
library(rnaturalearth)

# Load geographic boundaries of France, Spain, Morocco, Portugal, and Algeria
countries <- ne_countries(scale = "medium", country = c("France", "Spain", "Morocco", "Portugal", "Algeria"), returnclass = "sf")

# Remove French Guiana and Atlantic French territories
countries <- countries[!(countries$geounit %in% c("French Guiana", "Guadeloupe", "Martinique", "Saint Pierre and Miquelon", 
                                                  "Reunion", "Mayotte", "New Caledonia", "French Polynesia", 
                                                  "Wallis and Futuna", "Saint Barthelemy", "Saint Martin")), ]

# Convert TAB_pixel_LC to an sf object
TAB_pixel_LC_sf <- st_as_sf(TAB_pixel_LC, coords = c("long", "lat"), crs = 4326)

highlight_point <- st_as_sf(
  data.frame(lon = -7.60, lat = 31.80403),
  coords = c("lon", "lat"),
  crs = st_crs(countries)  # Use the same CRS as your map
)

# Step 4: Plot the map with quantile-based color scale
map <- ggplot(data = countries) +
  geom_sf(fill = "#EBEBEB", color = "black") +
  geom_sf(data = TAB_pixel_LC_sf, aes(color = TAB_pixel_LC$Foffset), size = 0.5, show.legend = FALSE) +
  scale_color_manual(values = color_palette, name = "Offset Quantile") +
  geom_sf(data = highlight_point, color = "black", size = 3, shape = 21, fill = "lightblue", stroke = 1) +
  coord_sf(xlim = c(-15, 15), ylim = c(28, 52), expand = FALSE) +
  theme_minimal() +
  labs(title = "Spatial offset Tassaout") +
  theme(
    panel.background = element_blank(),
    legend.position = "right",
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 8)
  )

map








plot(GOMar$GOM, GOMar$BLUP_FFD, 
     xlab = "GOMar", 
     ylab = "BLUP_FFD", 
     pch = 19, col = "steelblue")

# Add regression line
abline(model, col = "red", lwd = 2)

summary(model)




#------------------------------------------------------------------------------------------------------
#Spatial genomic offset of wild sample at the niche extreames. North Corse (F9) and south Morocco (M30)
#------------------------------------------------------------------------------------------------------
F <- GEA_124[rownames(GEA_124) == "OES_F9_10_S62_L004", ] #OES_M30_16_S58_L004 OES_F9_10_S62_L004 

FRDA <- predict(RDA_all_enriched, newdata=F, type="wa")
FRDA<-as.data.frame(FRDA)



TAB_pixel_LC$offsetF9 <- (FRDA$RDA1 - TAB_pixel_LC$RDA1)^2 * 0.684 + 
  (FRDA$RDA2 - TAB_pixel_LC$RDA2)^2 * 0.1095 + 
  (FRDA$RDA3 - TAB_pixel_LC$RDA3)^2 * 0.06071+
  (FRDA$RDA4 - TAB_pixel_LC$RDA4)^2 * 0.03739 
hist(TAB_pixel_LC$offsetF9)

level1 <- quantile(TAB_pixel_LC$offset, probs = 0.10, na.rm = TRUE)
level2 <- quantile(TAB_pixel_LC$offset, probs = 0.20, na.rm = TRUE)
level3 <- quantile(TAB_pixel_LC$offset, probs = 0.30, na.rm = TRUE)
level4 <- quantile(TAB_pixel_LC$offset, probs = 0.40, na.rm = TRUE)
level5 <- quantile(TAB_pixel_LC$offset, probs = 0.50, na.rm = TRUE)
level6 <- quantile(TAB_pixel_LC$offset, probs = 0.60, na.rm = TRUE)
level7 <- quantile(TAB_pixel_LC$offset, probs = 0.70, na.rm = TRUE)
level8 <- quantile(TAB_pixel_LC$offset, probs = 0.80, na.rm = TRUE)
level9 <- quantile(TAB_pixel_LC$offset, probs = 0.90, na.rm = TRUE)



# Compute breaks for the column
sd_breaks <- c( min(TAB_pixel_LC$offsetM30, na.rm = TRUE), level1, level2, level3, level4, level5, level6, level7, level8, level9, max(TAB_pixel_LC$offset, na.rm = TRUE))



# Create a color palette from blue to yellow
library(RColorBrewer)
color_palette <- brewer.pal(10, "PuOr")
# color_palette <- c(
#   "#004d00",  # very dark green
#   "#228B22",  # forest green
#   "#66C200",  # yellow-green
#   "#CCCC00",  # mustard yellow
#   "#FFD700",  # golden yellow
#   "#FFA500",  # orange
#   "#FF8C00",  # dark orange
#   "#FF4500",  # orange-red
#   "#B22222",  # firebrick
#   "#8B0000"   # dark red
# )

# Assign colors based on quantiles
TAB_pixel_LC$Foffset <- cut(TAB_pixel_LC$offsetM30, breaks = sd_breaks, labels = color_palette)



library(ggplot2)
library(sf)
library(rnaturalearth)


# Load geographic boundaries of France, Spain, Morocco, Portugal, and Algeria
countries <- ne_countries(scale = "medium", country = c("France", "Spain", "Morocco", "Portugal", "Algeria"), returnclass = "sf")

# Remove French Guiana and Atlantic French territories
countries <- countries[!(countries$geounit %in% c("French Guiana", "Guadeloupe", "Martinique", "Saint Pierre and Miquelon", 
                                                  "Reunion", "Mayotte", "New Caledonia", "French Polynesia", 
                                                  "Wallis and Futuna", "Saint Barthelemy", "Saint Martin")), ]

# Convert TAB_pixel_LC to an sf object
TAB_pixel_LC_sf <- st_as_sf(TAB_pixel_LC, coords = c("long", "lat"), crs = 4326)

# --- Convert offset to percentile rank (0-100) ---
TAB_pixel_LC$offset_percentile <- rank(TAB_pixel_LC$offsetM30, na.last = "keep") / 
  sum(!is.na(TAB_pixel_LC$offsetM30)) * 100

# --- Plot with continuous colorblind-safe gradient + colorbar legend ---
M30 <- ggplot(data = countries) +
  geom_sf(fill = "#EBEBEB", color = "black") +
  geom_sf(data = TAB_pixel_LC_sf, aes(color = TAB_pixel_LC$offset_percentile), 
          size = 0.05, show.legend = TRUE) +
  scale_color_gradientn(
    colours = c("#2D004B", "#542788", "#8073AC", "#B2ABD2",
                "grey85", "grey85",
                "#FDB863", "#E08214", "#B35806", "#7F3B08"),
    values = c(0, 0.15, 0.30, 0.42,
               0.48, 0.52,
               0.58, 0.70, 0.85, 1),
    limits = c(0, 100),
    breaks = c(0, 25, 50, 75, 100),
    name = "Offset\nPercentile (%)",
    guide = guide_colorbar(
      barheight = unit(3, "cm"),
      barwidth  = unit(0.4, "cm"),
      ticks.colour = "black",
      frame.colour = "black"
    )
  )+
  coord_sf(xlim = c(-15, 15), ylim = c(28, 52), expand = FALSE) +
  theme_bw(base_size = 10) +
  labs(title = "Adaptive Landscape M30") +
  theme(
    panel.background = element_blank(),
    legend.position = "none",
    plot.title = element_text(size = 10)
    # legend.title = element_text(size = 8),
    # legend.text = element_text(size = 7)
  )
ggsave(
  filename = "M30.jpeg",
  plot = M30,
  width = 3,       # inches
  height = 3,      # inches
  dpi = 600,       # publication-quality resolution
  units = "in",
  device = "jpeg"
)


# Load required libraries
library(ggplot2)
library(ggridges)
library(dplyr)
library(tidyr)

library(ggplot2)
library(ggridges)

ridge_data <- data.frame(
  Offset = c(TAB_pixel_LC$offsetF9, TAB_pixel_LC$offsetM30),
  Variable = rep(c("offsetF9", "offsetM30"), each = nrow(TAB_pixel_LC))
)

# Rename variables for clarity
ridge_data$Variable <- dplyr::recode(
  ridge_data$Variable,
  offsetF9  = "Spatial offset F9",
  offsetM30 = "Spatial offset M30"
)

ridge_plot <- ggplot(ridge_data, aes(x = Offset, y = Variable, fill = Variable)) +
  geom_density_ridges(scale = 1.2, alpha = 0.8, color = "black") +
  scale_fill_manual(values = c("Spatial offset F9" = "#3d5a80",
                               "Spatial offset M30" = "#e07a5f")) +
  theme_minimal(base_size = 12) +
  labs(
    title = "Distribution of Spatial Offset Values",
    x = "Offset value",
    y = NULL
  ) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 12),
    axis.text.y = element_text(face = "bold", size = 12),
    axis.text.x = element_text(size = 12),
    legend.position = "none"
  ) +
  coord_cartesian(xlim = c(0, 0.1))  # limit x-axis range
# Save as high-resolution JPEG
ggsave(
  filename = "ridge_offset_plot.jpeg",
  plot = ridge_plot,
  width = 7,       # inches
  height = 5,      # inches
  dpi = 600,       # publication-quality resolution
  units = "in",
  device = "jpeg"
)
#-----------------------------------------------
#PCA of cultivars just with GEA QTL from lfmm WW
#-----------------------------------------------
#PCA
library(FactoMineR)
library(factoextra)

GEA_wild_cultivar<-rbind(GEA_124, GEA_cultivars_maf)

res.pcacultivar<-PCA(GEA_wild_cultivar, scale.unit = FALSE, ncp = 5, graph = TRUE)
ind <- get_pca_ind(res.pcacultivar)
pca_data <- as.data.frame(ind$coord)
pca_data_wild<-cbind(pca_data[1:142,], group = Variables_142WW[,3])
pca_data_wild<-data.frame( geno = rownames(pca_data_wild), pca_data_wild)
pca_data_cul<-pca_data[143:461, ]
pca_data_cul$group = "cultivar"

write.table(pca_data_cul,'pca_data_cul.txt')
pca_data_cul<-read.csv("pca_data_cul.csv")
pca_data_cul <- pca_data_cul %>%
  filter_all(all_vars(. != ""))

pca_wild_cult<-rbind(pca_data_wild,pca_data_cul )
qq<-ggplot() +
  geom_hline(yintercept=0, linetype="dashed", color = gray(.80), linewidth=0.6) +
  geom_vline(xintercept=0, linetype="dashed", color = gray(.80), linewidth=0.6) +
  geom_point(data = pca_data, aes(x=Dim.1, y=Dim.2), size = 2.5) +
  xlab("PC1: 27%") + ylab("PC2: 10%") +
  guides(color=guide_legend(title="Group")) +
  theme_bw(base_size = 11, base_family = "Times") +
  theme(panel.background = element_blank(), legend.background = element_blank(), panel.grid = element_blank(), plot.background = element_blank(), legend.text=element_text(size=rel(.8)), strip.text = element_text(size=11))
qq

data_pca_flow <- read.delim("PCA_GEA_flow_class_cultivar.txt", 
                            header = TRUE, 
                            na.strings = c("NA", ""), 
                            stringsAsFactors = FALSE, 
                            fill = TRUE)
data_pca_flow <- na.omit(data_pca_flow)

qq <- ggplot() +
  geom_hline(yintercept=0, linetype="dashed", color = gray(.80), linewidth=0.6) +
  geom_vline(xintercept=0, linetype="dashed", color = gray(.80), linewidth=0.6) +
  geom_point(data = pca_wild_cult, 
             aes(x = Dim.1, y = Dim.2, fill = group), 
             size = 3.5, shape = 21, color = "black", stroke = 0.5) +
  scale_fill_manual(values = c("purple", "darkblue", "lightgrey", "lightblue", "darkgreen", "darkorange")) +
  xlab("PC1: 15%") + ylab("PC2: 7%") +
  guides(fill = guide_legend(title = "Group")) +
  theme_bw(base_size = 11, base_family = "Times") +
  theme(
    panel.background = element_blank(), 
    legend.background = element_blank(), 
    panel.grid = element_blank(), 
    plot.background = element_blank(), 
    legend.text = element_text(size = rel(.8)), 
    strip.text = element_text(size = 11)
  )

qq

boxplot(Dim.2 ~ group, data = pca_data_cul)
model<-lm(Dim.2 ~ group, data = pca_data_cul)
summary(model)
