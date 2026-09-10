
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
geno155<-fread("D:/D/vcf_file_GEA_leccino/geno_155.txt")
geno155 <- as.data.frame(geno155)  # rownames only work on data.frame, not data.table
rownames(geno155) <- geno155[[1]]  # assign first column as row names
geno155 <- geno155[ , -1]          # remove the first column

##### Wilde East
listWE<-read.table("list_WE.txt")
genoWE<- geno155[rownames(geno155)%in% listWE$V1, ]

geno229 <- read.vcfR("D:/D/vcf_file_GEA_leccino/WC229_Admixed_lec24_DP10_100_miss090_ind085_mac1.vcf.recode.vcf")#import vcf file
GI <- vcfR2genind(geno229)#transfrom file in genind object
geno229 <- as.data.frame(GI)

geno229 <- geno229 %>%
  dplyr::select(ends_with(".0"))

# list of wild with ancestry q>0.6
list_q06<-read.table("D:/C/Desktop/Leccino24/PopulationStructure/list_154WW.txt")
geno06 <- geno229[
  rownames(geno229) %in% list_q06$ rownames.geno154_WW.,
  ,
  drop = FALSE
]

geno06 <- geno06[, colnames(geno06) %in% colnames(geno155), drop = FALSE]
geno06 <- apply(geno06, 2, function(x) {
  x[is.na(x)] <- median(x, na.rm = TRUE)
  x
})

geno06<- rbind(geno06, geno155)


# environmental variable
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


locations = read.csv("D:/C/Desktop/Leccino24/417_wild_EW.csv")
# Extract environmental variables
Env <- data.frame(
  extract(ras_current_var, locations[, c(3, 4)])
)

# Combine locations + environmental variables
Env_all <- cbind(locations, Env)
Env_all[, c("clay", "N", "pH", "sand")] <- 
  Env_all[, c("clay", "N", "pH", "sand")] / 10


# centering 
dataWild_q06W <- Env_all[Env_all$vcf_name%in% list_q06$rownames.geno154_WW, ]
dataWild_q06W <- dataWild_q06W %>%
  mutate(LAT_classes = cut(lat,
                           breaks = c(-Inf, 35, 40, 45),
                           labels = c("low_lat", "med_lat", "high_lat"),
                           right = FALSE))
test_env <- dataWild_q06W[, c("bio2", "bio10", "bio11", "bio15", "bio18", "bio19", "clay", "N", "pH", "sand")]
Env <- scale(test_env, center=TRUE, scale=TRUE)
# Extract the centering values
env_center <- attr(Env, "scaled:center")
# Extract the scaling values
env_scale <- attr(Env, "scaled:scale")
#transform into dataset
Env <- as.data.frame(Env)
Variables_q06W<-data.frame(geno=dataWild_q06W$vcf_name, lat_classes = dataWild_q06W$LAT_classes, lat = dataWild_q06W$lat,  Env )
Variables_q06W <- na.omit(Variables_q06W)
geno_q06WW <- geno06[
  match(Variables_q06W$geno, rownames(geno06)),
  ,
  drop = FALSE
]

rownames(geno_q06WW) <- Variables_q06W$geno

#-------------
#maf filtering
#------------

Y <- geno_q06WW
# Function to calculate MAF for each column (SNP)
calculate_maf <- function(geno_col) {
  geno_col <- na.omit(geno_col)
  allele_freq <- sum(geno_col) / (2 * length(geno_col))  # assumes diploid, genotypes 0/1/2
  maf <- min(allele_freq, 1 - allele_freq)
  return(maf)
}

# Apply function to each SNP (column)
maf_values <- apply(geno_q06WW, 2, calculate_maf)
# Filter threshold, e.g., keep SNPs with MAF >= 0.05
maf_threshold <- 0.05
geno_q06WW_maf <- geno_q06WW[, maf_values >= maf_threshold]

write.table(geno_q06WW_maf,"geno_q06WW_maf.txt")

#-------------
#run LFMM GEA
#-------------

## Use latent factor for covariable correction
# latent factor temperature variable
Y <- geno_q06WW_maf
Y <- as.matrix(geno_q06WW_maf)

sel_latent<- data.frame(Variables_q06W%>% dplyr::select(bio2, bio10, bio11, bio15, bio18, bio19,clay, N, pH, sand))
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
thres <- 0.05/ncol(geno_q06WW_maf)
signif_bonf <- which(pv$pvalues < thres)
GEA_bonferroni <- data.frame(index = signif_bonf, 
                             pvalue = pv$pvalues[signif_bonf])

PvaluesGEA_lfmm<-data.frame(pv$pvalues)

#define cadidate mod.lfmm2#define cadidate loci for GO 

GEA_lfmm <- data.frame(pvalue = pv$pvalues[-log10(pv$pvalue) >5])

write.csv(GEA_bonferroni, "GEA_bonferroni_q06.csv")#selected 50 SNPs
write.csv(GEA_lfmm, "GEA_lfmm_all_var_log5_q06.csv")#selected 255 SNPs
write.csv(pv$pvalues, "GEA_all_var_lfmm_q06.csv")# all SNPs


#plotting Mhanattan plot using the library qqman

library(qqman)

Manhattan_q06 <- read.csv(file = "GEA_all_var_lfmm_q06.csv", header=TRUE) #import the p value result for precipitation
# =========================================================
# Function to prepare LFMM results for qqman
# =========================================================

prepare_manhattan <- function(pvalue_file) {
  
  # Read LFMM p-values
  dat <- read.csv(
    file = pvalue_file,
    header = TRUE,
    stringsAsFactors = FALSE
  )
  
  # The first column contains the LFMM SNP names
  # The second column contains the p-values
  colnames(dat)[1:2] <- c("SNP_raw", "P")
  
  # -------------------------------------------------------
  # Clean SNP names
  # -------------------------------------------------------
  
  dat$SNP <- dat$SNP_raw
  
  # Remove "Response "
  dat$SNP <- gsub(
    "^Response ",
    "",
    dat$SNP
  )
  
  # Remove ".0.value"
  dat$SNP <- gsub(
    "\\.0\\.value$",
    "",
    dat$SNP
  )
  
  # -------------------------------------------------------
  # Extract chromosome
  # -------------------------------------------------------
  
  # GWHEUUU00000001 -> 1
  # GWHEUUU00000002 -> 2
  # GWHEUUU00000010 -> 10
  
  dat$CHR <- as.numeric(
    sub(
      "_.*$",
      "",
      dat$SNP
    ) |>
      sub(
        "^GWHEUUU000000",
        "",
        x = _
      )
  )
  
  # -------------------------------------------------------
  # Extract physical position
  # -------------------------------------------------------
  
  # GWHEUUU00000001_1_4178813
  #                         ^^^^^^^
  
  dat$BP <- as.numeric(
    sub(
      ".*_[0-9]+_([0-9]+)$",
      "\\1",
      dat$SNP
    )
  )
  
  # -------------------------------------------------------
  # Keep qqman columns
  # -------------------------------------------------------
  
  dat <- dat[, c("SNP", "CHR", "BP", "P")]
  
  # Sort by chromosome and position
  dat <- dat[
    order(dat$CHR, dat$BP),
  ]
  
  rownames(dat) <- NULL
  
  return(dat)
}
Manhattan_q06 <- prepare_manhattan(
  "GEA_all_var_lfmm_q06.csv"
)

tiff(
  filename = "Manhattan_all_06.tiff",
  width = 12,
  height = 4,
  units = "in",
  res = 600,
  compression = "lzw"
)

manhattan(
  Manhattan_q06,
  col = c("darkgreen", "gray60"),
  genomewideline = 5,
  cex.axis = 1.0,
  cex.lab = 1.1
)

dev.off()

#-----------------
# Filter GEA
#-----------------
rownames(GEA_lfmm) <- gsub("^Response |\\.value$", "", rownames(GEA_lfmm))
list_GEA <- data.frame(SNP = rownames(GEA_lfmm))

GEA_lfmm <- geno_q06WW_maf[, colnames(geno_q06WW_maf) %in% list_GEA$SNP, drop = FALSE]
write.table(GEA_lfmm, "GEA_lfmm_all_var_q06.txt")

GEA_lfmm_all_var<-read.table("GEA_lfmm_all_var_q06.txt")


#----------------------------------------
#Cultivars genetic data and GEA filtering
#------------------------------------------

#upload genotypic file whole collection
geno_cultivar<- read.vcfR("D:/D/vcf_file_GEA_leccino/Cultivar_319_lec24_DP10_100_miss090_ind085_mac1.vcf.recode.vcf")


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
RDA_all_enriched<-rda(GEA_lfmm_all_var ~ bio2 + bio10 + bio11 + bio15	+ bio18 + bio19 + clay + N+ pH + sand , Variables_q06W)
summary(eigenvals(RDA_all_enriched, model = "constrained"))
RsquareAdj(RDA_all_enriched)
plot(RDA_all_enriched)
sqrt(vif.cca(RDA_all_enriched))

# plot Geographic regions


TAB_gen <- data.frame(geno = row.names(scores(RDA_all_enriched , display = "sites")), scores(RDA_all_enriched, display = "sites", scaling = "sites"))

Geno <- merge(TAB_gen, Variables_q06W[, 1:5] ,by="geno")
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
  xlab("RDA 1: 62%") + ylab("RDA 2: 13%") +
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

ggsave(
  filename = "RDA_biplot_lat_q06.png",
  plot = loading_geno_all_enriched_lat,      # Optional if last plot was your desired one
  width = 3,             # In inches (default)
  height = 2.1,             # In inches
  units = "in",           # Can be "in", "cm", or "mm"
  dpi = 300               # Resolution (important for publications)
)












# list of wild with ancestry q>0.8
list_q08<-read.table("D:/C/Desktop/Leccino24/PopulationStructure/list_126WW.txt")
geno08 <- geno229[
  rownames(geno229) %in% list_q08$ rownames.geno126_WW,
  ,
  drop = FALSE
]

geno08 <- geno08[, colnames(geno08) %in% colnames(geno155), drop = FALSE]
geno08 <- apply(geno08, 2, function(x) {
  x[is.na(x)] <- median(x, na.rm = TRUE)
  x
})

geno08<- rbind(geno08, geno155)


# environmental variable
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


locations = read.csv("D:/C/Desktop/Leccino24/417_wild_EW.csv")
# Extract environmental variables
Env <- data.frame(
  extract(ras_current_var, locations[, c(3, 4)])
)

# Combine locations + environmental variables
Env_all <- cbind(locations, Env)
Env_all[, c("clay", "N", "pH", "sand")] <- 
  Env_all[, c("clay", "N", "pH", "sand")] / 10


# centering 
dataWild_q08W <- Env_all[Env_all$vcf_name%in% list_q08$rownames.geno126_WW., ]
dataWild_q08W <- dataWild_q08W %>%
  mutate(LAT_classes = cut(lat,
                           breaks = c(-Inf, 35, 40, 45),
                           labels = c("low_lat", "med_lat", "high_lat"),
                           right = FALSE))
test_env <- dataWild_q08W[, c("bio2", "bio10", "bio11", "bio15", "bio18", "bio19", "clay", "N", "pH", "sand")]
Env <- scale(test_env, center=TRUE, scale=TRUE)
# Extract the centering values
env_center <- attr(Env, "scaled:center")
# Extract the scaling values
env_scale <- attr(Env, "scaled:scale")
#transform into dataset
Env <- as.data.frame(Env)
Variables_q08W<-data.frame(geno=dataWild_q08W$vcf_name, lat_classes = dataWild_q08W$LAT_classes, lat = dataWild_q08W$lat,  Env )
Variables_q08W <- na.omit(Variables_q08W)
geno_q08WW <- geno08[
  match(Variables_q08W$geno, rownames(geno08)),
  ,
  drop = FALSE
]

rownames(geno_q08WW) <- Variables_q08W$geno
#-------------
#maf filtering
#------------

Y <- geno_q08WW
# Function to calculate MAF for each column (SNP)
calculate_maf <- function(geno_col) {
  geno_col <- na.omit(geno_col)
  allele_freq <- sum(geno_col) / (2 * length(geno_col))  # assumes diploid, genotypes 0/1/2
  maf <- min(allele_freq, 1 - allele_freq)
  return(maf)
}

# Apply function to each SNP (column)
maf_values <- apply(geno_q08WW, 2, calculate_maf)
# Filter threshold, e.g., keep SNPs with MAF >= 0.05
maf_threshold <- 0.05
geno_q08WW_maf <- geno_q08WW[, maf_values >= maf_threshold]

write.table(geno_q08WW_maf,"geno_q08WW_maf.txt")

#-------------
#run LFMM GEA
#-------------

## Use latent factor for covariable correction
# latent factor temperature variable
Y <- geno_q08WW_maf
Y <- as.matrix(geno_q08WW_maf)

sel_latent<- data.frame(Variables_q08W%>% dplyr::select(bio2, bio10, bio11, bio15, bio18, bio19,clay, N, pH, sand))
write.env(sel_latent, "latent_all_variable.env")
X = read.table("latent_all_variable.env")
X <- as.matrix(X)

mod.lfmm2 <- lfmm2(input = Y, env = X, K = 2, effect.sizes = TRUE)
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
thres <- 0.05/ncol(geno_q08WW_maf)
signif_bonf <- which(pv$pvalues < thres)
GEA_bonferroni <- data.frame(index = signif_bonf, 
                             pvalue = pv$pvalues[signif_bonf])

PvaluesGEA_lfmm<-data.frame(pv$pvalues)

#define cadidate mod.lfmm2#define cadidate loci for GO 

GEA_lfmm <- data.frame(pvalue = pv$pvalues[-log10(pv$pvalue) >5])

write.csv(GEA_bonferroni, "GEA_bonferroni_q08.csv")#selected 50 SNPs
write.csv(GEA_lfmm, "GEA_lfmm_all_var_log5_q08.csv")#selected 255 SNPs
write.csv(pv$pvalues, "GEA_all_var_lfmm_q08.csv")# all SNPs


#plotting Mhanattan plot using the library qqman

library(qqman)
Manhattan_q08 <- prepare_manhattan(
  "GEA_all_var_lfmm_q08.csv"
)


tiff(
  filename = "Manhattan_all_08.tiff",
  width = 12,
  height = 4,
  units = "in",
  res = 600,
  compression = "lzw"
)

manhattan(
  Manhattan_q08,
  col = c("darkgreen", "gray60"),
  genomewideline = 5,
  cex.axis = 1.0,
  cex.lab = 1.1
)

dev.off()

#-----------------
# Filter GEA
#-----------------
rownames(GEA_lfmm) <- gsub("^Response |\\.value$", "", rownames(GEA_lfmm))
list_GEA <- data.frame(SNP = rownames(GEA_lfmm))

GEA_lfmm <- geno_q08WW_maf[, colnames(geno_q08WW_maf) %in% list_GEA$SNP, drop = FALSE]
write.table(GEA_lfmm, "GEA_lfmm_all_var_q08.txt")

GEA_lfmm_all_var<-read.table("GEA_lfmm_all_var_q08.txt")


#----------------------------------------
#Cultivars genetic data and GEA filtering
#------------------------------------------

#upload genotypic file whole collection
geno_cultivar<- read.vcfR("D:/D/vcf_file_GEA_leccino/Cultivar_319_lec24_DP10_100_miss090_ind085_mac1.vcf.recode.vcf")


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
RDA_all_enriched<-rda(GEA_lfmm_all_var ~ bio2 + bio10 + bio11 + bio15	+ bio18 + bio19 + clay + N+ pH + sand , Variables_q08W)
summary(eigenvals(RDA_all_enriched, model = "constrained"))
RsquareAdj(RDA_all_enriched)
plot(RDA_all_enriched)
sqrt(vif.cca(RDA_all_enriched))

# plot Geographic regions


TAB_gen <- data.frame(geno = row.names(scores(RDA_all_enriched , display = "sites")), scores(RDA_all_enriched, display = "sites", scaling = "sites"))

Geno <- merge(TAB_gen, Variables_q08W[, 1:5] ,by="geno")
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
  xlab("RDA 1: 61%") + ylab("RDA 2: 16%") +
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

ggsave(
  filename = "RDA_biplot_lat_q08.png",
  plot = loading_geno_all_enriched_lat,      # Optional if last plot was your desired one
  width = 3,             # In inches (default)
  height = 2.1,             # In inches
  units = "in",           # Can be "in", "cm", or "mm"
  dpi = 300               # Resolution (important for publications)
)






file <- "D:/C/Desktop/Olive_GO_paper/Submission_Horticulture_reasearch/revision_HR/revision_2/GEA_sensitivity_analysis.xlsx"

GEA_sensitivity_analysis <- read_excel(
  file,
  sheet = "Sheet1",
  col_names = TRUE
)

dim(GEA_sensitivity_analysis)
head(GEA_sensitivity_analysis, 20)
library(dplyr)

GEA_sets <- GEA_sensitivity_analysis %>%
  mutate(across(
    everything(),
    ~ gsub("^Response |\\.value$", "", .)
  ))

sets <- lapply(GEA_sets, function(x) {
  x <- x[!is.na(x) & x != ""]
  unique(x)
})

names(sets)

overlap_coefficient <- function(A, B) {
  length(intersect(A, B)) / min(length(A), length(B))
}

jaccard <- function(A, B) {
  length(intersect(A, B)) / length(union(A, B))
}

# Pairwise comparisons
comparisons <- list(
  q07_vs_q06 = c("GEA_bonferroni_q07", "GEA_bonferroni_q06"),
  q07_vs_q08 = c("GEA_bonferroni_q07", "GEA_bonferroni_q08"),
  q06_vs_q08 = c("GEA_bonferroni_q06", "GEA_bonferroni_q08")
)

sensitivity_summary <- do.call(rbind, lapply(comparisons, function(x) {
  
  A <- sets[[x[1]]]
  B <- sets[[x[2]]]
  
  data.frame(
    Set1 = x[1],
    Set2 = x[2],
    N_Set1 = length(A),
    N_Set2 = length(B),
    Shared_SNPs = length(intersect(A, B)),
    Overlap_coefficient = overlap_coefficient(A, B),
    Jaccard = jaccard(A, B)
  )
}))

sensitivity_summary
q07 <- sets$GEA_bonferroni_q07

data.frame(
  Comparison = c("q > 0.6", "q > 0.8"),
  N_SNPs = c(
    length(sets$GEA_bonferroni_q06),
    length(sets$GEA_bonferroni_q08)
  ),
  Shared_with_q07 = c(
    length(intersect(q07, sets$GEA_bonferroni_q06)),
    length(intersect(q07, sets$GEA_bonferroni_q08))
  ),
  Percent_q07_recovered = c(
    100 * length(intersect(q07, sets$GEA_bonferroni_q06)) / length(q07),
    100 * length(intersect(q07, sets$GEA_bonferroni_q08)) / length(q07)
  )
)

sensitivity_summary2 <- data.frame(
  Threshold = c("q > 0.6", "q > 0.7", "q > 0.8"),
  Individuals = c(150, 142, 126),
  Significant_SNPs = c(
    length(sets$GEA_bonferroni_q06),
    length(sets$GEA_bonferroni_q07),
    length(sets$GEA_bonferroni_q08)
  )
)

sensitivity_summary2
sensitivity_summary2$SNPs_per_individual <-
  sensitivity_summary2$Significant_SNPs /
  sensitivity_summary2$Individuals

sensitivity_summary2

# =========================================================
# Publication-quality sensitivity analysis table
# =========================================================

library(dplyr)
library(gridExtra)
library(grid)
library(gtable)

# ---------------------------------------------------------
# 1. Clean SNP identifiers
# ---------------------------------------------------------

GEA_sets <- GEA_sensitivity_analysis %>%
  mutate(
    across(
      everything(),
      ~ gsub("^Response |\\.value$", "", .)
    )
  )

sets <- lapply(GEA_sets, function(x) {
  x <- x[!is.na(x) & x != ""]
  unique(x)
})


# ---------------------------------------------------------
# 2. Functions
# ---------------------------------------------------------

overlap_coefficient <- function(A, B) {
  
  A <- unique(na.omit(A))
  B <- unique(na.omit(B))
  
  if (length(A) == 0 || length(B) == 0) {
    return(NA_real_)
  }
  
  length(intersect(A, B)) / min(length(A), length(B))
}


jaccard <- function(A, B) {
  
  A <- unique(na.omit(A))
  B <- unique(na.omit(B))
  
  if (length(A) == 0 || length(B) == 0) {
    return(NA_real_)
  }
  
  length(intersect(A, B)) / length(union(A, B))
}


# ---------------------------------------------------------
# 3. Pairwise overlap analysis
# ---------------------------------------------------------

comparisons <- list(
  "q > 0.7 vs q > 0.6" =
    c("GEA_bonferroni_q07", "GEA_bonferroni_q06"),
  
  "q > 0.7 vs q > 0.8" =
    c("GEA_bonferroni_q07", "GEA_bonferroni_q08"),
  
  "q > 0.6 vs q > 0.8" =
    c("GEA_bonferroni_q06", "GEA_bonferroni_q08")
)


overlap_table <- do.call(
  rbind,
  lapply(names(comparisons), function(comp) {
    
    x <- comparisons[[comp]]
    
    A <- sets[[x[1]]]
    B <- sets[[x[2]]]
    
    data.frame(
      Comparison = comp,
      `SNPs (Set 1)` = length(A),
      `SNPs (Set 2)` = length(B),
      `Shared SNPs` = length(intersect(A, B)),
      `Overlap coefficient` =
        overlap_coefficient(A, B),
      Jaccard = jaccard(A, B),
      check.names = FALSE
    )
  })
)

# Round statistics
overlap_table$`Overlap coefficient` <-
  sprintf("%.3f", overlap_table$`Overlap coefficient`)

overlap_table$Jaccard <-
  sprintf("%.3f", overlap_table$Jaccard)

overlap_table

# =========================================================
# 4. Create publication-quality table
# =========================================================

table_theme <- ttheme_minimal(
  
  base_size = 10,
  base_family = "serif",
  
  core = list(
    fg_params = list(
      fontsize = 10,
      fontfamily = "serif"
    ),
    bg_params = list(
      fill = "white",
      col = NA
    ),
    padding = unit(
      c(4, 8),
      "mm"
    )
  ),
  
  colhead = list(
    fg_params = list(
      fontsize = 10,
      fontface = "bold",
      fontfamily = "serif"
    ),
    bg_params = list(
      fill = "white",
      col = NA
    ),
    padding = unit(
      c(5, 8),
      "mm"
    )
  )
)


tbl <- tableGrob(
  overlap_table,
  rows = NULL,
  theme = table_theme
)


# ---------------------------------------------------------
# 5. Journal-style horizontal rules
# ---------------------------------------------------------

# Top rule
tbl <- gtable::gtable_add_grob(
  tbl,
  grobs = segmentsGrob(
    x0 = unit(0, "npc"),
    x1 = unit(1, "npc"),
    y0 = unit(1, "npc"),
    y1 = unit(1, "npc"),
    gp = gpar(lwd = 1)
  ),
  t = 1,
  l = 1,
  r = ncol(tbl)
)


# Header separator
tbl <- gtable::gtable_add_grob(
  tbl,
  grobs = segmentsGrob(
    x0 = unit(0, "npc"),
    x1 = unit(1, "npc"),
    y0 = unit(0, "npc"),
    y1 = unit(0, "npc"),
    gp = gpar(lwd = 0.7)
  ),
  t = 2,
  l = 1,
  r = ncol(tbl)
)


# Bottom rule
tbl <- gtable::gtable_add_grob(
  tbl,
  grobs = segmentsGrob(
    x0 = unit(0, "npc"),
    x1 = unit(1, "npc"),
    y0 = unit(0, "npc"),
    y1 = unit(0, "npc"),
    gp = gpar(lwd = 1)
  ),
  t = nrow(tbl),
  l = 1,
  r = ncol(tbl)
)


# =========================================================
# 6. Export as high-resolution TIFF
# =========================================================

tiff(
  filename = "GEA_sensitivity_overlap_table.tiff",
  width = 190,
  height = 65,
  units = "mm",
  res = 600,
  compression = "lzw"
)

grid.newpage()
grid.draw(tbl)

dev.off()


