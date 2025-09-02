rasterlayer<-raster(paste("C:/Users/rocchetti/Desktop/Leccino24/ENM_403/proj_CurrentEM_ent_West.med.clim.soil_ensemble.tif"))
r_vals <- getValues(rasterlayer)

threshold <- quantile(r_vals, probs = 0.80, na.rm = TRUE)
hist(r_vals)
hist(r_vals, main = "weighted mean distribution", xlab = "Suitability", col = "grey", border = "black")


  # Add a vertical line at the threshold
 abline(v = threshold, col = "red", lwd = 2, lty = 2)

 # Optionally add text label
 text(threshold, par("usr")[4]*0.9, labels = paste0("80th percentile = ", round(threshold, 3)), pos = 4, col = "red")








setwd("D:/C/Desktop/Leccino24/Landscape_156WWE")

library(fit.models)
library(adegenet)
library(vegan)
library(geosphere)
library(vcfR)
library(data.table)
#enter vcf file
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
geno155<-read.table("D:/vcf_file_GEA_leccino/geno_155.txt")

##### Wilde East
listWE<-read.table("list_WE.txt")
genoWE<- geno155[rownames(geno155)%in% listWE$V1, ]

############################## Environlent Wild Weast

data_wild<- read.csv("Env_155_WWE.csv", header = TRUE)
data_wild <- data_wild %>%
  mutate(LAT_classes = cut(lat,
                           breaks = c(-Inf, 35, 40, 45),
                           labels = c("low_lat", "med_lat", "high_lat"),
                           right = FALSE))

#################################### GEA all together
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
Variables_142WW<-data.frame(geno=dataWild_142W$id,group = dataWild_142W$group, region = dataWild_142W$region, lat_classes = dataWild_142W$LAT_classes, lat = dataWild_142W$lat,  Env )

genoWW<- geno155[rownames(geno155)%in% list142WW$V1, ]



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
genoWW_maf<-read.table("genoWW_maf.txt")


####  run LFMM GEA


## Use latent factor for covariable correction
# latent factor temperature variable
Y <- genoWW_maf
sel_latent<- data.frame(Variables_142WW%>% dplyr::select(bio2, bio10, bio11, bio15, bio18, bio19,clay, N, pH, sand))
write.env(sel_latent, "latent_all_variable.env")
X = read.table("latent_all_variable.env")

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

PvaluesGEA_lfmm<-data.frame(pv$pvalues)

#define cadidatemod.lfmm2#define cadidate loci for GO 

GEA_lfmm <- data.frame(pvalue = pv$pvalues[-log10(pv$pvalue) >5])

write.csv(GEA_lfmm, "GEA_lfmm_all_var_log5.csv")#selected 255 SNPs
write.csv(pv$pvalues, "GEA_all_var_lfmm.csv")# all SNPs


#plotting Mhanattan plot using the library qqman

library(qqman)
Manhattan_all <- read.csv(file = "GEA_all_var_lfmm.csv", header=TRUE) #import the p value result for precipitation
jpeg(file = "/lustre/rocchettil/Manh_RDA_prec.jpeg")
manhattan(Manhattan_all, col = c("darkgreen", "gray60"),genomewideline = 5)

hist(Manhattan_all$P)
dev.off()

########################################    enriched redundancy analysis (RDA only with 255 GEA QTL)  #######################
list_255<-read.table("listGEA_255.txt",  header=TRUE)
GEA_lfmm_255<-  genoWW_maf[, colnames(genoWW_maf)%in% list_255$SNP]
write.table(GEA_lfmm_255, "GEA_lfmm_all_var_255.txt")

GEA_lfmm_all_var<-read.table("GEA_lfmm_all_var_255.txt")


### tentative selecting GEA that have variations among cultivars

GEA_124<- GEA_lfmm_all_var[, colnames(GEA_lfmm_all_var)%in% colnames(GEA_cultivars_maf)]
write.table(GEA_124,"GEA_124_WW.txt")
GEA_124<-read.table("GEA_124_WW.txt")

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
  geom_point(data = Geno, aes(x=RDA1, y=RDA2, fill = lat), linewidth = 2.5, shape = 21, color = "black", stroke = 0.8, size = 4) +
  #scale_fill_manual(values = c("lightblue","darkgreen", "darkorange")) +
  scale_fill_gradientn(colors = c("#d73027","#fc8d59", "#fee090", "#91bfdb", "#4575b4"), 
                       name = "Latitude")+
  geom_segment(data = TAB_var, aes(xend=RDA1, yend=RDA2, x=0, y=0), colour="black", linewidth =0.15, linetype=1, arrow=arrow(length = unit(0.02, "npc"))) +
  geom_label_repel(data = TAB_var, aes(x=RDA1, y=RDA2, label = row.names(TAB_var)), size = 3.2, family = "Times") +
  xlab("RDA 1: 68%") + ylab("RDA 2: 11%") +
  guides(color=guide_legend(title="Latitude gradient")) +
  theme_bw(base_size = 11, base_family = "Times") +
  theme(panel.background = element_blank(), legend.background = element_blank(), panel.grid = element_blank(), plot.background = element_blank(), legend.text=element_text(size=rel(.8)), strip.text = element_text(size=11))
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


ggsave(
  filename = "RDA_biplot_lat.png",
  plot = loading_geno_all_enriched_lat,      # Optional if last plot was your desired one
  width = 5,             # In inches (default)
  height = 4,             # In inches
  units = "in",           # Can be "in", "cm", or "mm"
  dpi = 300               # Resolution (important for publications)
)

#### enter spatial pixel values

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

###### mapping adaptive landscape### map of GEA

# Extract RDA values
a1 <- TAB_pixel_LC$RDA1
a2 <- TAB_pixel_LC$RDA2

# Compute the distance from the origin
distance <- sqrt(a1^2 + a2^2)

# Assign colors based on quadrants and the 5th sector (circle radius < 0.5)
TAB_pixel_LC$color <- ifelse(distance < 0.25, "#717171",  # 5th sector - Purple
                             ifelse(a1 > 0 & a2 > 0, "#377EB8",  # Quadrant 1 - Red
                                    ifelse(a1 < 0 & a2 > 0, "red",  # Quadrant 2 - Blue
                                           ifelse(a1 < 0 & a2 < 0, "#FF7F00",  # Quadrant 3 - Green
                                                  "limegreen"))))  # Quadrant 4 - Orange

# Update ggplot with quadrant-based colors and 5th sector
pp <- ggplot() +
  geom_hline(yintercept = 0, linetype = "dashed", color = gray(0.80), size = 0.6) +
  geom_vline(xintercept = 0, linetype = "dashed", color = gray(0.80), size = 0.6) +
  geom_point(data = TAB_pixel_LC, aes(x = RDA1, y = RDA2, color = color), size = 2) +  # Use sector colors
  geom_segment(data = TAB_var, aes(xend = RDA1, yend = RDA2, x = 0, y = 0), 
               colour = "black", size = 0.15, linetype = 1, 
               arrow = arrow(length = unit(0.20, "cm"), type = "closed")) +
  geom_label_repel(data = TAB_var, aes(x = RDA1, y = RDA2, label = row.names(TAB_var)), 
                   size = 4, family = "Times") +
  xlab("RDA 1: 68%") + 
  ylab("RDA 2: 11%") +
  theme_bw(base_size = 9, base_family = "Times") +
  theme(panel.background = element_blank(), 
        legend.background = element_blank(), 
        panel.grid = element_blank(), 
        plot.background = element_blank(), 
        legend.text = element_text(size = rel(0.8)), 
        strip.text = element_text(size = 9)) +
  scale_color_identity()  # Use predefined colors directly

pp

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
  coord_sf(xlim = c(-15, 15), ylim = c(28, 52), expand = FALSE) +
  theme_minimal() +
  #labs(title = "Adaptive Landscape") +
  theme(panel.background = element_blank())
map

ggsave(
  filename = "adaptive_landscape_color.jpg",
  plot=map,
  dpi=300,
  width = 5,
  height = 7,
  units = 'in'
)

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
             size = 4, shape = 21, fill = colors, stroke = 0.2) +
  geom_hline(yintercept = 0, linetype = "dashed", color = gray(0.80), size = 0.6) +
  geom_vline(xintercept = 0, linetype = "dashed", color = gray(0.80), size = 0.6) +
  geom_segment(data = TAB_var, aes(xend = RDA1, yend = RDA2, x = 0, y = 0), 
               colour = "black", size = 0.15, linetype = 1, 
               arrow = arrow(length = unit(0.20, "cm"), type = "closed")) +
  geom_label_repel(data = TAB_var, aes(x = RDA1, y = RDA2, label = row.names(TAB_var)), 
                   size = 4, family = "Times") +
  xlab("RDA 1: 68%") + 
  ylab("RDA 2: 11%") +
  theme_bw(base_size = 9, base_family = "Times") +
  theme(panel.background = element_blank(), 
        legend.background = element_blank(), 
        panel.grid = element_blank(), 
        plot.background = element_blank(), 
        legend.text = element_text(size = rel(0.8)), 
        strip.text = element_text(size = 9)) +
  scale_color_identity()

pp

ggsave(
  filename = "biplot_landscape.png",
  plot = pp,      # Optional if last plot was your desired one
  width = 5,             # In inches (default)
  height = 5,             # In inches
  units = "in",           # Can be "in", "cm", or "mm"
  dpi = 300               # Resolution (important for publications)
)


library(patchwork)

layout <- "
ABB
CBB
"


final_plot <- pp+map+plot_spacer()+
  plot_layout(design = layout)+
   plot_annotation(
    title = "Adaptive Landscape",
    tag_levels = 'a'
  )


  final_plot
  
  ((pp / plot_spacer() + plot_layout(guides = 'keep')) | map) + plot_layout(guides = 'collect')


library(scales)  # for rescale()
library(ggplot2)
library(ggrepel)

# Rescale RDA1 and RDA2 to [0, 1]
TAB_pixel_LC$R <- rescale(TAB_pixel_LC$RDA1, to = c(0, 1))
TAB_pixel_LC$G <- rescale(TAB_pixel_LC$RDA2, to = c(0, 1))
TAB_pixel_LC$B <- rescale(sqrt(TAB_pixel_LC$RDA1^2 + TAB_pixel_LC$RDA2^2), to = c(0, 1))  # Optional: use magnitude as blue channel

# Convert to RGB hex code
TAB_pixel_LC$color <- rgb(TAB_pixel_LC$R, TAB_pixel_LC$G, TAB_pixel_LC$B)

# Plot with continuous RGB
pp <- ggplot() +
  geom_hline(yintercept = 0, linetype = "dashed", color = gray(0.80), size = 0.6) +
  geom_vline(xintercept = 0, linetype = "dashed", color = gray(0.80), size = 0.6) +
  geom_point(data = TAB_pixel_LC, aes(x = RDA1, y = RDA2, color = color), size = 2) +
  geom_segment(data = TAB_var, aes(xend = RDA1, yend = RDA2, x = 0, y = 0), 
               colour = "black", size = 0.15, linetype = 1, 
               arrow = arrow(length = unit(0.20, "cm"), type = "closed")) +
  geom_label_repel(data = TAB_var, aes(x = RDA1, y = RDA2, label = row.names(TAB_var)), 
                   size = 4, family = "Times") +
  xlab("RDA 1: 68%") + 
  ylab("RDA 2: 11%") +
  theme_bw(base_size = 9, base_family = "Times") +
  theme(panel.background = element_blank(), 
        legend.background = element_blank(), 
        panel.grid = element_blank(), 
        plot.background = element_blank(), 
        legend.text = element_text(size = rel(0.8)), 
        strip.text = element_text(size = 9)) +
  scale_color_identity()

pp

    ######################################################### Estimation of cultivar offset #######################

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
               arrow = arrow(length = unit(0.20, "cm"), type = "closed"),
               color = "black", size = 0.15) +
  geom_label_repel(data = TAB_var, 
                   aes(x = RDA1 * arrow_scale, y = RDA2 * arrow_scale, label = row.names(TAB_var)),
                   size = 3, family = "Times") +
  xlab("RDA 1: 68%") + 
  ylab("RDA 2: 11%") +
  labs(title = "RDA Cultivar and Spatial Pixel Predictions") +
  theme_bw(base_size = 11, base_family = "Times") +
  theme(
    panel.background  = element_blank(),
    plot.background   = element_blank(),
    legend.background = element_blank(),
    panel.grid        = element_blank(),
    legend.text       = element_text(size = rel(0.8)),
    strip.text        = element_text(size = 11)
  )
ll

ggsave(
  filename = "biplot_CUL-spaatial.png",
  plot = ll,      # Optional if last plot was your desired one
  width = 7,             # In inches (default)
  height = 6,             # In inches
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
  geom_point(data = wild_cult_pred, aes(x = RDA1, y = RDA2, fill = group, shape = group),size = 3, color = "black", stroke = 0.8)+
  scale_shape_manual(values = c(24,21,21))+
  scale_fill_manual(values=c('#E69F00',"grey48","lightblue"))+
  #scale_size_manual(values=c(3,3))+
  geom_segment(data = TAB_var, aes(xend=RDA1*0.3, yend=RDA2*0.3, x=0, y=0), colour="black", size=0.15, linetype=1, arrow = arrow(length=unit(0.20,"cm"),type = "closed")) +
  geom_label_repel(data = TAB_var, aes(x=RDA1*0.3, y=RDA2*0.3, label = row.names(TAB_var)), size = 3, family = "Times") +
  xlab("RDA 1: 68%") + ylab("RDA 2: 11%") +
  #guides(legend(title="Group")) +
  theme_bw(base_size = 11, base_family = "Times") +
  theme(panel.background = element_blank(), legend.background = element_blank(), panel.grid = element_blank(), plot.background = element_blank(), legend.text=element_text(size=rel(.8)), strip.text = element_text(size=11))+
  labs(title = "predicted Wild and Cultivar with GEAs")
hh
library(ggpubr)
proj<-ggarrange(ll, hh,nrow=2, ncol=1)
ggsave(
  filename = "proj_cultivar.jpg",
  plot=proj,
  dpi=300,
  width = 5,
  height = 7,
  units = 'in'
)

### Map specific cultivar offset

###### Frantoio

F <- GEA_cultivars[rownames(GEA_cultivars) == "Arbequina", ]

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
  
  # Step 4: Plot the map with quantile-based color scale
  Picholine_Marocaine_map <- ggplot(data = countries) +
    geom_sf(fill = "#EBEBEB", color = "black") +
    geom_sf(data = TAB_pixel_LC_sf, aes(color = TAB_pixel_LC$Foffset), size = 0.5, show.legend = FALSE) +
    scale_color_manual(values = color_palette, name = "Offset Quantile") +
    coord_sf(xlim = c(-15, 15), ylim = c(28, 52), expand = FALSE) +
    theme_minimal() +
    labs(title = "Adaptive Landscape  Picholine_Marocaine") +
    theme(
      panel.background = element_blank(),
      legend.position = "right",
      legend.title = element_text(size = 10),
      legend.text = element_text(size = 8)
    )
  


cultivar_offset<- ggarrange(Picholine2_map, Picual2_map, Manzanilla_map, Picholine_Marocaine2_map, nrow=2, ncol=2 )
cultivar_offset<- (Picholine2_map+Picual2_map)/ (Manzanilla_Cacerena_map + Picholine_Marocaine2_map)
ggsave(
  filename = " Picholine_Marocaine_map.jpg",
  plot= Picholine_Marocaine_map,
  dpi=300,
  width = 8,
  height = 8,
  units = 'in'
)

  #### PCA of cultivars just with GEA QTL from lfmm WW

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


#### Comparison with data laila marakesh


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

hist(GoMar$FFD)
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

model <- lm(GO ~ logFFD, data = GoMar)
#abline(model, col = "red", lwd = 2)
summary(model)

c<-ggplot(GoMar, aes(x = logFFD, y = GO)) +
  geom_point(shape = 21, size = 5, color = "black", fill = "grey", alpha = 0.8) +
  geom_abline(intercept = coef(model)[1], slope = coef(model)[2],
              color = "red", linetype = "solid", linewidth = 1.2) +
  #scale_fill_manual(values = c("purple", "lightgrey", "darkblue")) +
  theme_minimal() +
  labs(title = "GO vs FFD by Class",
       x = "log(FFD)",
       y = "GO",
       fill = "Class") +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 12)
  )
c

library(ggpubr)
d<-ggarrange(b,c, nrow = 1, ncol=2)
d

ggsave("GO_vs_FFD_by_Class.jpg", plot = d, width = 7, height = 5, dpi = 300)
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





##### extreams wild

  F <- GEA_124[rownames(GEA_124) == "OES_F9_08_S16_L004", ]

FRDA <- predict(RDA_all_enriched, newdata=F, type="wa")
FRDA<-as.data.frame(FRDA)



TAB_pixel_LC$offset <- (FRDA$RDA1 - TAB_pixel_LC$RDA1)^2 * 0.684 + 
  (FRDA$RDA2 - TAB_pixel_LC$RDA2)^2 * 0.1095 + 
  (FRDA$RDA3 - TAB_pixel_LC$RDA3)^2 * 0.06071+
  (FRDA$RDA4 - TAB_pixel_LC$RDA4)^2 * 0.03739 
hist(TAB_pixel_LC$offset)

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

# Step 4: Plot the map with quantile-based color scale
m30 <- ggplot(data = countries) +
  geom_sf(fill = "#EBEBEB", color = "black") +
  geom_sf(data = TAB_pixel_LC_sf, aes(color = TAB_pixel_LC$Foffset), size = 0.5, show.legend = FALSE) +
  scale_color_manual(values = color_palette, name = "Offset Quantile") +
  coord_sf(xlim = c(-15, 15), ylim = c(28, 52), expand = FALSE) +
  theme_minimal() +
  labs(title = "Adaptive Landscape F9") +
  theme(
    panel.background = element_blank(),
    legend.position = "right",
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 8)
  )

m30

##### 2100 scenarios

#### enter spatial pixel values

library(raster)
library("readxl")


bio2_2100<- raster(paste("D:/raster files/MPI_esm/ssp585/bio2_2100_masked.tif"))
bio10_2100<- raster(paste("D:/raster files/MPI_esm/ssp585/bio10_2100_masked.tif"))
bio11_2100<- raster(paste("D:/raster files/MPI_esm/ssp585/bio11_2100_masked.tif"))
bio15_2100<- raster(paste("D:/raster files/MPI_esm/ssp585/bio15_2100_masked.tif"))
bio18_2100<- raster(paste("D:/raster files/MPI_esm/ssp585/bio18_2100_masked.tif"))
bio19_2100<- raster(paste("D:/raster files/MPI_esm/ssp585/bio19_2100_masked.tif"))
soilN<- raster(paste("D:/raster files/Current_ENM_clipped_biova/resampled_soilN.tif"))
soilpH<- raster(paste("D:/raster files/Current_ENM_clipped_biova/resampled_soilpH.tif"))

soilclay<- raster(paste("D:/raster files/Current_ENM_clipped_biova/resampled_soilclay.tif"))
soilsand<- raster(paste("D:/raster files/Current_ENM_clipped_biova/resampled_soilsand_1.tif"))

names(bio2_2100) = 'bio2_2100'
names(bio10_2100) = 'bio10_2100'
names(bio11_2100) = 'bio11_2100'
names(bio15_2100) = 'bio15_2100'
names(bio18_2100) = 'bio18_2100'
names(bio19_2100) = 'bio19_2100'
names(soilN ) = 'N'
names(soilpH) = 'pH'
names(soilclay) = 'clay'
names(soilsand) = 'sand'



#stack the different raster file
ras_future_var<-stack(c(bio2_2100,bio10_2100, bio11_2100, bio15_2100, bio18_2100, bio19_2100, soilclay,soilN,soilpH, soilsand))

pixel_2100 <- as.data.frame(rasterToPoints(ras_future_var))
pixel_2100 <- data.frame(x=pixel_2100$x, y=pixel_2100$y, bio2=pixel_2100$bio2_2100,bio10=pixel_2100$bio10_2100,bio11=pixel_2100$bio11_2100,bio15=pixel_2100$bio15_2100,bio18=pixel_2100$bio18_2100,bio19=pixel_2100$bio19_2100,clay=pixel_2100$clay/10,N=pixel_2100$N/100,pH=pixel_2100$pH/10,sand=pixel_2100$sand/10)
pixel_2100<-na.omit(pixel_2100)
pixel_2100<- pixel_2100[pixel_2100$x>-10, ]
pixel_2100env<- pixel_2100%>% dplyr::select(bio2, bio10, bio11, bio15, bio18, bio19,clay, N, pH, sand)

scaled_pixel_2100 <- scale(pixel_2100env, center = env_center, scale = env_scale)
scaled_pixel_2100<-as.data.frame(scaled_pixel_2100)




#prediction of pixel in the RDA space
scaled_pixel_2100_LC <- predict(RDA_all_enriched, newdata=scaled_pixel_2100, type="lc")
TAB_pixel_2100_LC<- data.frame(lat = pixel$y, long = pixel$x, scaled_pixel_2100_LC)
TAB_var <- as.data.frame(scores(RDA_all_enriched, choices=c(1,2), display="bp"))


# Extract RDA values
a1 <- TAB_pixel_2100_LC$RDA1
a2 <- TAB_pixel_2100_LC$RDA2

# Compute the distance from the origin
distance <- sqrt(a1^2 + a2^2)

# Assign colors based on quadrants and the 5th sector (circle radius < 0.5)
TAB_pixel_2100_LC$color <- ifelse(distance < 0.05, "#717171",  # 5th sector - Purple
                             ifelse(a1 > 0 & a2 > 0, "#377EB8",  # Quadrant 1 - Red
                                    ifelse(a1 < 0 & a2 > 0, "red",  # Quadrant 2 - Blue
                                           ifelse(a1 < 0 & a2 < 0, "#FF7F00",  # Quadrant 3 - Green
                                                  "limegreen"))))  # Quadrant 4 - Orange

# Update ggplot with quadrant-based colors and 5th sector
pp <- ggplot() +
  geom_hline(yintercept = 0, linetype = "dashed", color = gray(0.80), size = 0.6) +
  geom_vline(xintercept = 0, linetype = "dashed", color = gray(0.80), size = 0.6) +
  geom_point(data = TAB_pixel_2100_LC, aes(x = RDA1, y = RDA2, color = color), size = 2) +  # Use sector colors
  geom_segment(data = TAB_var, aes(xend = RDA1*2, yend = RDA2*2, x = 0, y = 0), 
               colour = "black", size = 0.15, linetype = 1, 
               arrow = arrow(length = unit(0.20, "cm"), type = "closed")) +
  geom_label_repel(data = TAB_var, aes(x = RDA1*2, y = RDA2*2, label = row.names(TAB_var)), 
                   size = 4, family = "Times") +
  xlab("RDA 1: 64%") + 
  ylab("RDA 2: 11%") +
  theme_bw(base_size = 9, base_family = "Times") +
  theme(panel.background = element_blank(), 
        legend.background = element_blank(), 
        panel.grid = element_blank(), 
        plot.background = element_blank(), 
        legend.text = element_text(size = rel(0.8)), 
        strip.text = element_text(size = 9)) +
  scale_color_identity()  # Use predefined colors directly

pp

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
TAB_pixel_LC_sf <- st_as_sf(TAB_pixel_2100_LC, coords = c("long", "lat"), crs = 4326)  # Assumes 'longitude' and 'latitude' columns exist
# Create the map
map <- ggplot(data = countries) +
  geom_sf(fill = "#EBEBEB", color = "black") +  # Countries' borders
  geom_sf(data = TAB_pixel_LC_sf, aes(color = color), size = 0.05, show.legend = FALSE) +  # Points with custom colors
  scale_color_identity() +  # Use exact colors from the 'color' column
  coord_sf(xlim = c(-15, 15), ylim = c(28, 52), expand = FALSE) +  # Set geographic limits
  theme_minimal() +
  labs(title = "Adaptive Landscape") +
  theme(panel.background = element_blank())

#jpeg(file = "C:/Users/rocchetti/Desktop/running RDA GO/adaptive_landscape_rgb.jpeg",width = 18, height = 14, units = "cm", res = 800)
map

##### Marrakesh 2100

latitude <- 31.80403
longitude <- -7.60 

# Define a tolerance value (small range of acceptable difference)
tolerance <- 1e-2

# Filter the data using a tolerance for matching
Marrakech_2100_row <- TAB_pixel_2100_LC %>%
  filter(abs(lat - target_lat) < tolerance & abs(long - target_long) < tolerance)

Marrakech_2100_row<-as.data.frame(Marrakech_2100_row)
### offset at marakesh


Marrakech_2100_row <- Marrakech_2100_row[1, ]

# Perform the calculation
RDAscore_cul$offsetM2100 <- (Marrakech_2100_row$RDA1 - RDAscore_cul$RDA1)^2 * 0.684 + 
  (Marrakech_2100_row$RDA2 - RDAscore_cul$RDA2)^2 * 0.1095 + 
  (Marrakech_2100_row$RDA3 - RDAscore_cul$RDA3)^2 * 0.06071 + 
  (Marrakech_2100_row$RDA4 - RDAscore_cul$RDA4)^2 * 0.03739


GoMar<-write.table(RDAscore_cul, "Go_Mar_2100.txt")

GoMar<-read.csv("GO_marakesh.csv")
GoMar <- na.omit(GoMar)

boxplot(GO_2100 ~ class, data = GoMar)

model <- lm(offset_WW ~ class, data = GoMar)

summary(model)

# Perform Tukey's Honest Significant Difference test
tukey_result <- TukeyHSD(aov(model))

# View results
print(tukey_result)

hist(GoMar$FFD)
abline(v = c(120.62, 125.78), col = "red", lwd = 2, lty = 2)


plot(offset_WW ~ FFD, data = GoMar)

model <- lm(offset_WW ~ FFD, data = GoMar)
abline(model, col = "red", lwd = 2)
summary(model)


Gopresent_future<-read.csv("GO_present_future_morocco.csv")

Gopresent_future$class <- as.factor(Gopresent_future$class)
Gopresent_future$time <- as.factor(Gopresent_future$time)
ggplot(Gopresent_future, aes(x = class, y = offset, fill = time)) +
  geom_boxplot() +
  labs(title = "Offset by Class and Time", x = "Class", y = "Offset") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text.x = element_text(angle = 0, hjust = 0.5)
  )


# Fit the linear model with the interaction term
model <- lm(offset ~ class + time + class*time, data = Gopresent_future)
emmip(model, ~ class + time + class*time, CIs = TRUE)  # or time ~ class
tab_model(model, show.ci = FALSE, show.se = TRUE)
# Show the summary of the model



model<-lm(offset ~ FFD + time + FFD*time, data = Gopresent_future)

summary(model)

# Add the fitted values to the original dataset
Gopresent_future$predicted_offset <- predict(model, newdata = Gopresent_future)

# Create the ggplot
ggplot(Gopresent_future, aes(x = FFD, y = predicted_offset, color = time)) +
  geom_line() +  # Adds the lines for predicted values
  geom_point(aes(x = FFD, y = offset), alpha = 0.6, size = 3) +  # Increased point size
  labs(title = "Predicted Offset vs FFD by Time",
       x = "FFD",
       y = "Predicted Offset") +
  theme_minimal() +
  theme(legend.position = "top")

# cultivar offset 2100

F <- GEA_cultivars[rownames(GEA_cultivars) == "Aggezi_Akse1", ]

FRDA <- predict(RDA_all_enriched, newdata=F, type="wa")
FRDA<-as.data.frame(FRDA)



TAB_pixel_2100_LC$offset <- (FRDA$RDA1 - TAB_pixel_2100_LC$RDA1)^2 * 0.684 + 
  (FRDA$RDA2 - TAB_pixel_2100_LC$RDA2)^2 * 0.1095 + 
  (FRDA$RDA3 - TAB_pixel_2100_LC$RDA3)^2 * 0.06071+
  (FRDA$RDA4 - TAB_pixel_2100_LC$RDA4)^2 * 0.03739 
hist(TAB_pixel_2100_LC$offset)


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
sd_breaks_2100 <- c( min(TAB_pixel_LC$offset, na.rm = TRUE), level1, level2, level3, level4, level5, level6, level7, level8, level9, max(TAB_pixel_LC$offset, na.rm = TRUE))



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
TAB_pixel_2100_LC$Foffset <- cut(TAB_pixel_2100_LC$offset, breaks = sd_breaks, labels = color_palette)

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
TAB_pixel_2100_LC_sf <- st_as_sf(TAB_pixel_2100_LC, coords = c("long", "lat"), crs = 4326)

# Step 4: Plot the map with quantile-based color scale
map <- ggplot(data = countries) +
  geom_sf(fill = "#EBEBEB", color = "black") +
  geom_sf(data = TAB_pixel_LC_sf, aes(color = TAB_pixel_2100_LC$Foffset), size = 0.5, show.legend = FALSE) +
  scale_color_manual(values = color_palette, name = "Offset Quantile") +
  coord_sf(xlim = c(-15, 15), ylim = c(28, 52), expand = FALSE) +
  theme_minimal() +
  labs(title = "Adaptive Landscape Aggezi_Akse1_2100") +
  theme(
    panel.background = element_blank(),
    legend.position = "right",
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 8)
  )

map


############################################################### spatial offset marrakesh

TAB_pixel_2100_LC$offset <- (Marrakech_row$RDA1 - TAB_pixel_2100_LC$RDA1)^2 * 0.684 + 
  (Marrakech_row$RDA2 - TAB_pixel_2100_LC$RDA2)^2 * 0.1095 + 
  (Marrakech_row$RDA3 - TAB_pixel_2100_LC$RDA3)^2 * 0.06071+
  (Marrakech_row$RDA4 - TAB_pixel_2100_LC$RDA4)^2 * 0.03739 
hist(TAB_pixel_2100_LC$offset)


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
TAB_pixel_2100_LC$Foffset <- cut(TAB_pixel_2100_LC$offset, breaks = sd_breaks, labels = color_palette)

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
TAB_pixel_2100_LC_sf <- st_as_sf(TAB_pixel_2100_LC, coords = c("long", "lat"), crs = 4326)

highlight_point <- st_as_sf(
  data.frame(lon = -7.60, lat = 31.80403),
  coords = c("lon", "lat"),
  crs = st_crs(countries)  # Use the same CRS as your map
)

# Step 4: Plot the map with quantile-based color scale
map <- ggplot(data = countries) +
  geom_sf(fill = "#EBEBEB", color = "black") +
  geom_sf(data = TAB_pixel_2100_LC_sf, aes(color = TAB_pixel_2100_LC$Foffset), size = 0.5, show.legend = FALSE) +
  scale_color_manual(values = color_palette, name = "Offset Quantile") +
  geom_sf(data = highlight_point, color = "black", size = 3, shape = 21, fill = "lightblue", stroke = 1) +
  coord_sf(xlim = c(-15, 15), ylim = c(28, 52), expand = FALSE) +
  theme_minimal() +
  labs(title = "Spatial offset Tassaout_2100") +
  theme(
    panel.background = element_blank(),
    legend.position = "right",
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 8)
  )

map



##### future spatial offset wild extream

F <- GEA_124[rownames(GEA_124) == "OES_F9_08_S16_L004", ]

FRDA <- predict(RDA_all_enriched, newdata=F, type="wa")
FRDA<-as.data.frame(FRDA)


TAB_pixel_2100_LC$offset <- (FRDA$RDA1 - TAB_pixel_2100_LC$RDA1)^2 * 0.684 + 
  (FRDA$RDA2 - TAB_pixel_2100_LC$RDA2)^2 * 0.1095 + 
  (FRDA$RDA3 - TAB_pixel_2100_LC$RDA3)^2 * 0.06071+
  (FRDA$RDA4 - TAB_pixel_2100_LC$RDA4)^2 * 0.03739 
hist(TAB_pixel_2100_LC$offset)


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
sd_breaks_2100 <- c( min(TAB_pixel_LC$offset, na.rm = TRUE), level1, level2, level3, level4, level5, level6, level7, level8, level9, max(TAB_pixel_LC$offset, na.rm = TRUE))



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
TAB_pixel_2100_LC$Foffset <- cut(TAB_pixel_2100_LC$offset, breaks = sd_breaks, labels = color_palette)

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
TAB_pixel_2100_LC_sf <- st_as_sf(TAB_pixel_2100_LC, coords = c("long", "lat"), crs = 4326)

# Step 4: Plot the map with quantile-based color scale
map <- ggplot(data = countries) +
  geom_sf(fill = "#EBEBEB", color = "black") +
  geom_sf(data = TAB_pixel_LC_sf, aes(color = TAB_pixel_2100_LC$Foffset), size = 0.5, show.legend = FALSE) +
  scale_color_manual(values = color_palette, name = "Offset Quantile") +
  coord_sf(xlim = c(-15, 15), ylim = c(28, 52), expand = FALSE) +
  theme_minimal() +
  labs(title = "Adaptive Landscape OES_F9_08_S16_L004_2100") +
  theme(
    panel.background = element_blank(),
    legend.position = "right",
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 8)
  )

map























##### future vs present RDA plot
scaled_pixel_2100_LC <- predict(RDA_all_enriched, newdata=scaled_pixel_2100, type="lc", scaling = "sites")
TAB_pixel_2100_LC<- data.frame(lat = pixel$y, long = pixel$x, scaled_pixel_2100_LC)
TAB_pixel_2100_LC$time <-"future_2100"
scaled_pixel_LC <- predict(RDA_all_enriched, newdata=scaled_pixel, type="lc", scaling = "sites")
TAB_pixel_LC<- data.frame(lat = pixel$y, long = pixel$x, scaled_pixel_LC)
TAB_pixel_LC$time <-"current"
Tab<-rbind(TAB_pixel_LC, TAB_pixel_2100_LC)

TAB_var <- as.data.frame(scores(RDA_all_enriched, choices=c(1,2), display="bp"))

aa <- ggplot() +
  geom_hline(yintercept = 0, linetype = "dashed", color = gray(0.80), size = 0.6) +
  geom_vline(xintercept = 0, linetype = "dashed", color = gray(0.80), size = 0.6) +
  geom_point(data = Tab, aes(x = RDA1, y = RDA2, color = time), size = 2) +
  geom_segment(data = TAB_var, aes(xend = RDA1, yend = RDA2, x = 0, y = 0), 
               colour = "black", size = 0.15, linetype = 1, 
               arrow = arrow(length = unit(0.20, "cm"), type = "closed")) +
  geom_label_repel(data = TAB_var, aes(x = RDA1, y = RDA2, label = row.names(TAB_var)), 
                   size = 4, family = "Times") +
  xlab("RDA 1: 68%") + 
  ylab("RDA 2: 11%") +
  theme_bw(base_size = 9, base_family = "Times") +
  theme(panel.background = element_blank(), 
        legend.background = element_blank(), 
        panel.grid = element_blank(), 
        plot.background = element_blank(), 
        legend.text = element_text(size = rel(0.8)), 
        strip.text = element_text(size = 9)) +
  scale_color_manual(values = c("current" = "lightblue", "future_2100" = "orange"))  # blue & red

aa

boxplot(RDA2 ~ time, data = Tab,
        col = c("lightblue", "orange"),
        xlab = "Time",
        ylab = "RDA2")



########################## check adaptive value of admixed

genoadm <- read.vcfR("D:/vcf_file_GEA_leccino/WC229_Admixed_lec24_DP10_100_miss090_ind085_mac1.vcf.recode.vcf")#import vcf file
GI <- vcfR2genind(genoadm)#transfrom file in genind object
genoadm<-as.data.frame(GI)
genoadm<-genoadm%>% dplyr::select(ends_with(".0"))

GEA_admixed<- genoadm[, colnames(genoadm)%in% colnames(GEA_cultivars_maf)]

#imputation
for (i in 1:ncol(GEA_admixed))
{
  GEA_admixed[which(is.na(GEA_admixed[,i])),i] <- median(GEA_admixed[-which(is.na(GEA_admixed[,i])),i], na.rm=TRUE)
}


write.table(GEA_admixed, "GEA_admixed.txt")
GEA_admixed<-read.table("GEA_admixed.txt", header = T)


RDAscore_adm <- predict(RDA_all_enriched, newdata=GEA_admixed, type="wa")
RDAscore_adm<-as.data.frame(RDAscore_adm)


Tab_adm<- data.frame(geno = row.names(RDAscore_adm),RDAscore_adm[,1:2] )
Tab_adm$group<-"geno_Admixed"

### enter environmental data for lc prediction

adm_env<-read.csv("Admixed_env.csv")
adm_env <- data.frame(bio2=adm_env$bio2,bio10=adm_env$bio10,bio11=adm_env$bio11,bio15=adm_env$bio15,bio18=adm_env$bio18,bio19=adm_env$bio19,clay=adm_env$clay/10,N=adm_env$N/100,pH=adm_env$pH/10,sand=adm_env$sand/10)
scaled_adm_env <- scale(adm_env, center = env_center, scale = env_scale)
scaled_adm_env<-as.data.frame(scaled_adm_env)



RDA_env_adm<- predict(RDA_all_enriched, newdata=scaled_adm_env, type="lc")


Tab_env_adm<-data.frame(geno = row.names(GEA_admixed),RDA_env_adm[,1:2])
Tab_env_adm$group<-"env_Admixed"

TAB_var <- as.data.frame(scores(RDA_all_enriched, choices=c(1,2), display="bp"))

adm_geno_env_pred<-rbind(Tab_adm, Tab_env_adm)
adm_geno_env_pred$group <- as.factor(adm_geno_env_pred$group)

hh <- ggplot() +
  geom_hline(yintercept=0, linetype="dashed", color = gray(.80), size=0.6) +
  geom_vline(xintercept=0, linetype="dashed", color = gray(.80), size=0.6) +
  geom_point(data = adm_geno_env_pred, aes(x = RDA1, y = RDA2, fill = group, shape = group),size = 2.5, color = "black", stroke = 0.8)+
  scale_shape_manual(values = c(21,24))+
  scale_fill_manual(values=c("grey48", '#E69F00'))+
  #scale_size_manual(values=c(3,3))+
  geom_segment(data = TAB_var, aes(xend=RDA1*0.3, yend=RDA2*0.3, x=0, y=0), colour="black", size=0.15, linetype=1, arrow = arrow(length=unit(0.20,"cm"),type = "closed")) +
  geom_label_repel(data = TAB_var, aes(x=RDA1*0.3, y=RDA2*0.3, label = row.names(TAB_var)), size = 3, family = "Times") +
  xlab("RDA 1: 68%") + ylab("RDA 2: 11%") +
  #guides(legend(title="Group")) +
  theme_bw(base_size = 11, base_family = "Times") +
  theme(panel.background = element_blank(), legend.background = element_blank(), panel.grid = element_blank(), plot.background = element_blank(), legend.text=element_text(size=rel(.8)), strip.text = element_text(size=11))+
  labs(title = "predicted Wild and Cultivar with GEAs")
hh



















write.csv(adm_geno_env_pred, "prediction_adm_geno_env.csv")
adm_geno_env_pred<-read.csv("prediction_adm_geno_env.csv")
# Split the data by group
circle_data <- subset(adm_geno_env_pred, group == "geno_Admixed")
triangle_data <- subset(adm_geno_env_pred, group == "env_Admixed")

# Shared color vector
color_values <- c("darkorange", "lightblue",  "darkgreen")
names(color_values) <- unique(adm_geno_env_pred$country)  # make sure the levels match

library(ggplot2)
library(ggrepel)

hh <- ggplot() +
  geom_hline(yintercept = 0, linetype = "dashed", color = gray(.80), size = 0.6) +
  geom_vline(xintercept = 0, linetype = "dashed", color = gray(.80), size = 0.6) +
  
  # ▶️ Circles (filled)
  geom_point(data = circle_data,
             aes(x = RDA1, y = RDA2, fill = country),
             shape = 21, size = 2.5, stroke = 0.8, color = "black") +
  
  # 🔺 Triangles (outline only, colored by country)
  geom_point(data = triangle_data,
             aes(x = RDA1, y = RDA2, color = country),
             shape = 24, size = 2.5, stroke = 1.2, fill = NA) +
  
  scale_fill_manual(values = color_values) +
  scale_color_manual(values = color_values) +
  
  geom_segment(data = TAB_var,
               aes(xend = RDA1 * 0.3, yend = RDA2 * 0.3, x = 0, y = 0),
               colour = "black", size = 0.15, linetype = 1,
               arrow = arrow(length = unit(0.20, "cm"), type = "closed")) +
  
  geom_label_repel(data = TAB_var,
                   aes(x = RDA1 * 0.3, y = RDA2 * 0.3, label = row.names(TAB_var)),
                   size = 3, family = "Times") +
  
  xlab("RDA 1: 68%") + ylab("RDA 2: 11%") +
  theme_bw(base_size = 11, base_family = "Times") +
  theme(
    panel.background = element_blank(),
    legend.background = element_blank(),
    panel.grid = element_blank(),
    plot.background = element_blank(),
    legend.text = element_text(size = rel(.8)),
    strip.text = element_text(size = 11)
  ) +
  labs(title = "Predicted Wild and Cultivar with GEAs")

hh

library(dplyr)
library(tidyr)
summary_df <- adm_geno_env_pred %>%
  dplyr::select(geno, country, RDA1, RDA2, group) %>%
  pivot_wider(names_from = group, values_from = c(RDA1, RDA2)) %>%
  mutate(
    RDA1_diff = RDA1_geno_Admixed - RDA1_env_Admixed,
    RDA2_diff = RDA2_geno_Admixed - RDA2_env_Admixed,
    euclidean_squared = RDA1_diff^2*0.684 + RDA2_diff^2*0.1095
  )

summary_df <- summary_df %>%
  filter(geno != "OES_M23_02_S5_L004")


write.table(summary_df, "adaptation_adm_lat.txt")

adaptation_adm<-read.csv("Adaptation_ADM_lat.csv")



library(ggplot2)
library(sf)
library(rnaturalearth)
library(dplyr)

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
adaptation_adm$Foffset <- cut(adaptation_adm$euclidean_squared, breaks = sd_breaks, labels = color_palette)

# 1. Load geographic boundaries
countries <- ne_countries(
  scale = "medium",
  country = c("France", "Spain", "Morocco", "Portugal", "Algeria"),
  returnclass = "sf"
)

# 2. Remove overseas French territories
countries <- countries[!(countries$geounit %in% c(
  "French Guiana", "Guadeloupe", "Martinique", "Saint Pierre and Miquelon", 
  "Reunion", "Mayotte", "New Caledonia", "French Polynesia", 
  "Wallis and Futuna", "Saint Barthelemy", "Saint Martin"
)), ]

# 3. Convert your adaptation dataset to sf object
adaptation_adm_sf <- st_as_sf(adaptation_adm, coords = c("long", "lat"), crs = 4326)


# Create the map with the custom color palette
map <- ggplot(data = countries) +
  geom_sf(fill = "#EBEBEB", color = "black") +
  geom_sf(data = adaptation_adm_sf, aes(color = adaptation_adm$Foffset), size = 4, show.legend = FALSE) +
  scale_color_manual(values = color_palette, name = "Offset Quantile") +
    coord_sf(xlim = c(-15, 15), ylim = c(28, 52), expand = FALSE) +
  theme_minimal() +
  labs(title = "Offset of Admixed populations") +
  theme(
    panel.background = element_blank(),
    legend.position = "right",
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 8)
  )

map

boxplot(
  TAB_pixel_LC$offset, # calcolated from extream wild pop F9
  adaptation_adm$euclidean_squared,
  names = c("wild Offset", "Admixed offset"),
  main = "Comparison of wild and admixed Offset",
  ylab = "Values",
  col = c("lightblue", "darkorange")
)


library(ggplot2)

ggplot(adaptation_adm, aes(x = "Local offset admixed", y = euclidean_squared)) +
  geom_boxplot(fill = "darkorange") +
  geom_hline(yintercept = 0.01, linetype = "dashed", color = "blue") +
  geom_hline(yintercept = 0.019, linetype = "dashed", color = "red") +
  annotate("text", x = 1, y = 0.011, label = "GO early blooming variety", color = "blue", hjust = 0, size = 4) +
  annotate("text", x = 1, y = 0.02, label = "GO late blooming variety", color = "red", hjust = 0, size = 4) +
  labs(
    x = "",
    y = "Local GO",
    title = "Comparison of Wild and Admixed Offset"
  ) +
  ylim(0, 0.022) +
  theme_minimal() +
  theme(
    axis.title.y = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    axis.text.x = element_text(size = 12),
    plot.title = element_text(size = 12, face = "bold", hjust = 0.5)
  )
  

ggarrange(map, m30, nrow=1, ncol=2)

boxplot(
  TAB_pixel_LC$offset,
  col = "purple",
  ylab = "Spatial GO Value",
  main = "Distribution of Spatial GO F9",
  ylim = c(0, 0.08),
  axes = FALSE  # suppress default axes
)

# Add custom axes
axis(2, at = seq(0, 0.08, by = 0.01))  # Y-axis with steps of 0.01

box()




























################ Plan B###################################################################################################################################

# Wild Environment datafile

#standardize bioclim variable
data_wild<- read.csv("Env_155_WWE.csv", header = TRUE)
test_env <- data_wild%>% dplyr::select(bio2, bio10, bio11, bio15, bio18, bio19,clay, N, pH, sand)
Env <- scale(test_env, center=TRUE, scale=TRUE) #each comumn will be subtracted by its mean and divided by its sd
# Extract the centering values
env_center <- attr(Env, "scaled:center") #mean of each variable
# Extract the scaling values
env_scale <- attr(Env, "scaled:scale") #standard deviation of each variable
#transform into dataset
Env <- as.data.frame(Env)


#combining geographic, Popstructure, environmental (scaled) variables
Variables <- data.frame(data_wild$id, data_wild$group,data_wild$region, data_wild$lat, data_wild$long,  Env)
names(Variables)[1]<- paste("geno")
names(Variables)[2]<- paste("group")
names(Variables)[3]<- paste("region")
names(Variables)[4]<- paste("lat")
names(Variables)[5]<- paste("long")

prcomp

RDAgeo_env <- rda(geno155 ~ bio2+bio10+bio11+	bio15	+ bio18 + bio19 + clay+ N+ pH+ sand , Variables)
sqrt(vif.cca(RDAgeo_env))
plot(RDAgeo_env)

RDAT <- rda(geno155 ~ bio2+bio10+bio11,Variables)
sqrt(vif.cca(RDAT))

#PCA to visualize how many latent factor I should use

library(FactoMineR)
library(factoextra)

res.pca155<-PCA(geno155, scale.unit = TRUE, ncp = 6, graph = TRUE)
eig_values <-res.pca155$eig
# Convert to a data frame
eig_df <- data.frame(PC = seq_along(eig_values[, 1]), 
                     Variance = eig_values[, 2])

eig_values10 <- eig_df %>% filter(PC < 10)

# Create the scree plot using ggplot2
ggplot(eig_values10, aes(x = PC, y = Variance)) +
  geom_bar(stat = "identity") +
  geom_line(aes(y = Variance), color = "red", size = 1) +
  geom_point(color = "red", size = 3) +
  labs(title = "Scree Plot", x = "Principal Components", y = "Percentage of Variance Explained") +
  theme_minimal()



## RDA for Genotype Environment Associations (GEA)

#Redundancy analysis can be used to identify GEA based on the Mhallanoise distance of SNPs in the RDA-biplot. Within the RDA model we can effectively correct for population structure  and geography (latitude and longitude) using them as covariates in the RDA model. As population structure correction we used latent factor derived from the LEA package.

#As first attempt I decided to run the anlysis seperate for temperature and precipitation variables.

#Temperature


## Use latent factor for covariable correction
# latent factor temperature variable
Y <- geno155
sel_temp<- data.frame(Env%>% dplyr::select(bio2, bio10, bio11))
write.env(sel_temp, "Temp_variable.env")
X = read.table("Temp_variable.env")

mod.lfmm2 <- lfmm2(input = Y, env = X, K = 3)
str(mod.lfmm2)
mod.lfmm2@U
#Merge latent factor to Variable
latent_temp<-data.frame(rownames(geno155), mod.lfmm2@U)
Temp_Var<-cbind(Variables,latent_temp)

#GEA Temperature
RDA_temp <- rda(geno155 ~ bio2+bio10+bio11 +  Condition(X1 + X2 + X3), Temp_Var)
summary(eigenvals(RDA_temp, model = "constrained"))
library(robust)
rdadapt<-function(rda,K)
{
  zscores<-rda$CCA$v[,1:as.numeric(K)]
  resscale <- apply(zscores, 2, scale)
  resmaha <- covRob(resscale, distance = TRUE, na.action= na.omit, estim="pairwiseGK")$dist
  lambda <- median(resmaha)/qchisq(0.5,df=K)
  reschi2test <- pchisq(resmaha/lambda,K,lower.tail=FALSE)
  qval <- qvalue(reschi2test)
  q.values_rdadapt<-qval$qvalues
  return(data.frame(p.values=reschi2test, q.values=q.values_rdadapt))
}

rdadapt_temp<- rdadapt(RDA_temp, 2)
## P-values threshold after Bonferroni correction
thres_env <- 0.05/length(rdadapt_temp$p.values)
## Identifying the loci that are below the p-value threshold
top_outliers <- data.frame(Loci = colnames(geno155)[which(rdadapt_temp$p.values<thres_env)], p.value = rdadapt_temp$p.values[which(rdadapt_temp$p.values<thres_env)], contig = unlist(lapply(strsplit(colnames(geno155)[which(rdadapt_temp$p.values<thres_env)], split = "_"), function(x) x[1])))
qvalue <- data.frame(Loci = colnames(geno155), p.value = rdadapt_temp$p.values, q.value = rdadapt_temp$q.value)
outliers <- data.frame(Loci = colnames(geno155)[which(rdadapt_temp$q.values<0.05)], p.value = rdadapt_temp$p.values[which(rdadapt_temp$q.values<0.05)])


#plot GEA temp

locus_scores <- scores(RDA_temp, choices=c(1:2), display="species", scaling="sites")
TAB_loci <- data.frame(names = row.names(locus_scores), locus_scores)
TAB_loci$type <- "Not associated"
TAB_loci$type[TAB_loci$names%in%outliers$Loci] <- "FDR"
TAB_loci$type[TAB_loci$names%in%top_outliers$Loci] <- "Bonferroni"
TAB_loci$type <- factor(TAB_loci$type, levels = c("Not associated", "FDR", "Bonferroni"))
TAB_var <- as.data.frame(scores(RDA_temp, choices=c(1,2), display="bp"))
loading_temp<-ggplot() +
  geom_hline(yintercept=0, linetype="dashed", color = gray(.80), linewidth=0.6) +
  geom_vline(xintercept=0, linetype="dashed", color = gray(.80), linewidth=0.6) +
  geom_point(data = TAB_loci, aes(x=RDA1, y=RDA2, colour = type), size = 2.5) +
  scale_color_manual(values = c("gray90", "#F9A242FF", "#6B4596FF")) +
  geom_segment(data = TAB_var, aes(xend=1.1*RDA1, yend=1.1*RDA2, x=0, y=0), colour="black", size=0.15, linetype=1, arrow=arrow(length = unit(0.02, "npc"))) +
  geom_label_repel(data = TAB_var, aes(x=1.1*RDA1, y=1.1*RDA2, label = row.names(TAB_var)), size = 3.8, family = "Times") +
  xlab("RDA 1: 57.8%") + ylab("RDA 2: 21.6%") +
  guides(color=guide_legend(title="Locus type")) +
  theme_bw(base_size = 11, base_family = "Times") +
  theme(panel.background = element_blank(), legend.background = element_blank(), panel.grid = element_blank(), plot.background = element_blank(), legend.text=element_text(size=rel(.8)), strip.text = element_text(size=11))
loading_temp
jpeg(file = "/lustre/rocchettil/RDA_temp_biplot.jpeg")
plot(loading_temp)
dev.off()

write.table(qvalue, "Temp_GEA_WWE_3.csv", append = FALSE, quote = TRUE, sep = " ",
            eol = "\n", na = "NA", dec = ".", row.names = FALSE,
            col.names = TRUE, qmethod = c("escape", "double"),
            fileEncoding = "")

Manhattan_temp <- read.csv(file = "Temp_GEA_WWE_3.csv", header=TRUE) #import the p value result for temperature
manhattan(Manhattan_temp, col = c("darkred", "gray60"),suggestiveline = -log10(0.001158020), genomewideline = -log10(1.924104e-07))
jpeg(file = "/lustre/rocchettil/Manh_RDA_temp.jpeg")
manhattan(Manhattan_temp, col = c("darkred", "gray60"),suggestiveline = -log10(0.001158020), genomewideline = -log10(1.924104e-07))
dev.off()

#P distribution
jpeg(file = "/lustre/rocchettil/Phist_Manh_RDA_temp")
hist(Manhattan_temp$P)
dev.off()

hist(qvalue$p.value)


  #Precipitation
  ```
  #latent factor precipitation variable
  
  Y <- geno155
  sel_prec<- data.frame(Env%>% dplyr::select(bio15, bio18, bio19))
  write.env(sel_prec, "prec_variable.env")
  X = read.table("prec_variable.env")
  mod.lfmm2 <- lfmm2(input = Y, env = X, K = 3)
  str(mod.lfmm2)
  mod.lfmm2@U
  #Merge latent factor to Variable
  latent_prec<-data.frame(rownames(geno155), mod.lfmm2@U)
  Prec_Var<-cbind(Variables,latent_prec)
  
  
  
  ## GEA Precipitation
  RDA_prec <- rda(geno155 ~ 	bio15	+ bio18 + bio19 +  Condition(X1 + X2 + X3 ), Prec_Var)
   
  summary(eigenvals(RDA_prec, model = "constrained"))
  
  rdadapt_prec<- rdadapt(RDA_prec, 2)
  ## P-values threshold after Bonferroni correction
  thres_env <- 0.05/length(rdadapt_prec$p.values)
  ## Identifying the loci that are below the p-value threshold
  top_outliers <- data.frame(Loci = colnames(geno155)[which(rdadapt_prec$p.values<thres_env)], p.value = rdadapt_prec$p.values[which(rdadapt_prec$p.values<thres_env)], contig = unlist(lapply(strsplit(colnames(geno155)[which(rdadapt_prec$p.values<thres_env)], split = "_"), function(x) x[1])))
  
  qvalue <- data.frame(Loci = colnames(geno155), p.value = rdadapt_prec$p.values, q.value = rdadapt_prec$q.value)
  outliers <- data.frame(Loci = colnames(geno155)[which(rdadapt_prec$q.values<0.05)], p.value = rdadapt_prec$p.values[which(rdadapt_prec$q.values<0.05)])
  
  #plot GEA precipitation
  
  locus_scores <- scores(RDA_prec, choices=c(1:2), display="species", scaling="sites")
  TAB_loci <- data.frame(names = row.names(locus_scores), locus_scores)
  TAB_loci$type <- "Not associated"
  TAB_loci$type[TAB_loci$names%in%outliers$Loci] <- "FDR"
  TAB_loci$type[TAB_loci$names%in%top_outliers$Loci] <- "Bonferroni"
  TAB_loci$type <- factor(TAB_loci$type, levels = c("Not associated", "FDR", "Bonferroni"))
  TAB_var <- as.data.frame(scores(RDA_prec, choices=c(1,2), display="bp"))
  loading_prec<-ggplot() +
    geom_hline(yintercept=0, linetype="dashed", color = gray(.80), size=0.6) +
    geom_vline(xintercept=0, linetype="dashed", color = gray(.80), size=0.6) +
    geom_point(data = TAB_loci, aes(x=RDA1, y=RDA2, colour = type), size = 2.5) +
    scale_color_manual(values = c("gray90", "#F9A242FF", "#6B4596FF")) +
    geom_segment(data = TAB_var, aes(xend=RDA1, yend=RDA2, x=0, y=0), colour="black", size=0.15, linetype=1, arrow=arrow(length = unit(0.02, "npc"))) +
    geom_label_repel(data = TAB_var, aes(x=1.1*RDA1, y=1.1*RDA2, label = row.names(TAB_var)), size = 2.8, family = "Times") +
    xlab("RDA 1: 48.5%") + ylab("RDA 2: 27.0%") +
    guides(color=guide_legend(title="Locus type")) +
    theme_bw(base_size = 11, base_family = "Times") +
    theme(panel.background = element_blank(), legend.background = element_blank(), panel.grid = element_blank(), plot.background = element_blank(), legend.text=element_text(size=rel(.8)), strip.text = element_text(size=11))
  loading_prec
  jpeg(file = "/lustre/rocchettil/RDA_prec_biplot.jpeg")
  plot(loading_prec)
  dev.off()
  
  
  write.table(qvalue, "Prec_GEA_WWE_3.csv", append = FALSE, quote = TRUE, sep = " ",
              eol = "\n", na = "NA", dec = ".", row.names = FALSE,
              col.names = TRUE, qmethod = c("escape", "double"),
              fileEncoding = "")
  
  #plotting Mhanattan plot using the library qqman
  
  library(qqman)
  Manhattan_prec <- read.csv(file = "Prec_GEA_WWE_3.csv", header=TRUE) #import the p value result for precipitation
  jpeg(file = "/lustre/rocchettil/Manh_RDA_prec.jpeg")
  manhattan(Manhattan_prec, col = c("blue", "gray60"),suggestiveline = -log10(0.0010233235), genomewideline = -log10(1.914289e-07))
  dev.off()
  
  #P distribution
  jpeg(file = "/lustre/rocchettil/Phist_Manh_RDA_prec.jpeg")
  hist(Manhattan_prec$P)
  dev.off()


# Soil variable
  

  #latent factor precipitation variable
  
  Y <- geno155
  sel_soil<- data.frame(Env%>% dplyr::select(clay, N, pH, sand))
  write.env(sel_soil, "soil_variable.env")
  X = read.table("soil_variable.env")
  mod.lfmm2 <- lfmm2(input = Y, env = X, K = 3)
  str(mod.lfmm2)
  mod.lfmm2@U
  #Merge latent factor to Variable
  latent_soil<-data.frame(rownames(geno155), mod.lfmm2@U)
  soil_Var<-cbind(Variables,latent_soil)
  
  
  
  ## GEA soil
  RDA_soil <- rda(geno155 ~ 	clay+ N+ pH+ sand+ Condition(X1 + X2 + X3 ), soil_Var)
  summary(eigenvals(RDA_soil, model = "constrained"))
  
  rdadapt_soil<- rdadapt(RDA_soil, 2)
  ## P-values threshold after Bonferroni correction
  thres_env <- 0.05/length(rdadapt_soil$p.values)
  ## Identifying the loci that are below the p-value threshold
  top_outliers <- data.frame(Loci = colnames(geno155)[which(rdadapt_soil$p.values<thres_env)], p.value = rdadapt_soil$p.values[which(rdadapt_soil$p.values<thres_env)], contig = unlist(lapply(strsplit(colnames(geno155)[which(rdadapt_soil$p.values<thres_env)], split = "_"), function(x) x[1])))
  
  qvalue <- data.frame(Loci = colnames(geno155), p.value = rdadapt_soil$p.values, q.value = rdadapt_soil$q.value)
  outliers <- data.frame(Loci = colnames(geno155)[which(rdadapt_soil$q.values<0.05)], p.value = rdadapt_soil$p.values[which(rdadapt_soil$q.values<0.05)])
  
  #plot GEA precipitation
  
  locus_scores <- scores(RDA_soil, choices=c(1:2), display="species", scaling="sites")
  TAB_loci <- data.frame(names = row.names(locus_scores), locus_scores)
  TAB_loci$type <- "Not associated"
  TAB_loci$type[TAB_loci$names%in%outliers$Loci] <- "FDR"
  TAB_loci$type[TAB_loci$names%in%top_outliers$Loci] <- "Bonferroni"
  TAB_loci$type <- factor(TAB_loci$type, levels = c("Not associated", "FDR", "Bonferroni"))
  TAB_var <- as.data.frame(scores(RDA_soil, choices=c(1,2), display="bp"))
  loading_soil<-ggplot() +
    geom_hline(yintercept=0, linetype="dashed", color = gray(.80), size=0.6) +
    geom_vline(xintercept=0, linetype="dashed", color = gray(.80), size=0.6) +
    geom_point(data = TAB_loci, aes(x=RDA1, y=RDA2, colour = type), size = 2.5) +
    scale_color_manual(values = c("gray90", "#F9A242FF", "#6B4596FF")) +
    geom_segment(data = TAB_var, aes(xend=RDA1, yend=RDA2, x=0, y=0), colour="black", size=0.15, linetype=1, arrow=arrow(length = unit(0.02, "npc"))) +
    geom_label_repel(data = TAB_var, aes(x=1.1*RDA1, y=1.1*RDA2, label = row.names(TAB_var)), size = 3.8, family = "Times") +
    xlab("RDA 1: 39.1%") + ylab("RDA 2: 24.2%") +
    guides(color=guide_legend(title="Locus type")) +
    theme_bw(base_size = 11, base_family = "Times") +
    theme(panel.background = element_blank(), legend.background = element_blank(), panel.grid = element_blank(), plot.background = element_blank(), legend.text=element_text(size=rel(.8)), strip.text = element_text(size=11))
  loading_soil
  jpeg(file = "/lustre/rocchettil/RDA_soil_biplot.jpeg")
  plot(loading_soil)
  dev.off()
  
  
  write.table(qvalue, "Soil_GEA_WWE_3.csv", append = FALSE, quote = TRUE, sep = " ",
              eol = "\n", na = "NA", dec = ".", row.names = FALSE,
              col.names = TRUE, qmethod = c("escape", "double"),
              fileEncoding = "")
  
  #plotting Mhanattan plot using the library qqman
  
  library(qqman)
  Manhattan_soil <- read.csv(file = "Soil_GEA_WWE_3.csv", header=TRUE) #import the p value result for precipitation
  jpeg(file = "/lustre/rocchettil/Manh_RDA_prec.jpeg")
  manhattan(Manhattan_soil, col = c("#ab7e4c", "gray60"),suggestiveline = -log10(0.001236238), genomewideline = -log10(1.961360e-07))
  dev.off()
  
  #P distribution
  jpeg(file = "/lustre/rocchettil/Phist_Manh_RDA_prec.jpeg")
  hist(Manhattan_soil$P)
  dev.off()


  
  
  
#################################### GEA all together RDA
  Y <- geno155
  sel_latent<- data.frame(Variables%>% dplyr::select(bio2, bio10, bio11, bio15, bio18, bio19,clay, N, pH, sand))
  write.env(sel_latent, "latent_all_variable.env")
  X = read.table("latent_all_variable.env")
  
  mod.lfmm2 <- lfmm2(input = Y, env = X, K = 3)
  str(mod.lfmm2)
  mod.lfmm2@U
  #Merge latent factor to Variable
  latent_all_var<-data.frame(rownames(geno155), mod.lfmm2@U)
  latent_all_var<-cbind(Variables,latent_all_var)  
  
  RDA_all <- rda(geno155 ~ 	bio2+ bio10+ bio11+ bio15+ bio18+ bio19+clay+ N+ pH+ sand+ Condition(X1 + X2 + X3 ), latent_all_var)
  summary(eigenvals(RDA_all, model = "constrained"))
  plot(RDA_all)
  
  # plot Geographic regions
  
  TAB_gen <- data.frame(geno = row.names(scores(RDA_all , display = "sites")), scores(RDA_all, display = "sites"))
  Geno <- merge(TAB_gen, Variables[, 1:5] ,by="geno")
  TAB_var <- as.data.frame(scores(RDA_all, choices=c(1,2), display="bp"))
  loading_geno_all_enriched_region<-ggplot() +
    geom_hline(yintercept=0, linetype="dashed", color = gray(.80), size=0.6) +
    geom_vline(xintercept=0, linetype="dashed", color = gray(.80), size=0.6) +
    geom_point(data = Geno, aes(x=RDA1, y=RDA2, fill = region), size = 2.5, shape = 21, color = "black", stroke = 0.8) +
    scale_fill_manual(values = c("lightblue","darkgreen", "darkorange", "darkblue")) +
    geom_segment(data = TAB_var, aes(xend=RDA1*5, yend=RDA2*5, x=0, y=0), colour="black", linewidth =0.15, linetype=1, arrow=arrow(length = unit(0.02, "npc"))) +
    geom_label_repel(data = TAB_var, aes(x=5*RDA1, y=5*RDA2, label = row.names(TAB_var)), size = 3.2, family = "Times") +
    xlab("RDA 1: 23%") + ylab("RDA 2: 12%") +
    guides(color=guide_legend(title="Latitude gradient")) +
    theme_bw(base_size = 11, base_family = "Times") +
    theme(panel.background = element_blank(), legend.background = element_blank(), panel.grid = element_blank(), plot.background = element_blank(), legend.text=element_text(size=rel(.8)), strip.text = element_text(size=11))
  loading_geno_all_enriched_region
  
  rdadapt_all<- rdadapt(RDA_all, 2) 
  ## P-values threshold after Bonferroni correction
  thres_env <- 0.05/length(rdadapt_all$p.values)
  ## Identifying the loci that are below the p-value threshold
  top_outliers <- data.frame(Loci = colnames(geno155)[which(rdadapt_all$p.values<thres_env)], p.value = rdadapt_all$p.values[which(rdadapt_all$p.values<thres_env)], contig = unlist(lapply(strsplit(colnames(geno155)[which(rdadapt_all$p.values<thres_env)], split = "_"), function(x) x[1])))
  
  qvalue <- data.frame(Loci = colnames(geno155), p.value = rdadapt_all$p.values, q.value = rdadapt_all$q.value)
  outliers <- data.frame(Loci = colnames(geno155)[which(rdadapt_all$q.values<0.05)], p.value = rdadapt_all$p.values[which(rdadapt_all$q.values<0.05)])
  hist(outliers$p.value)
  
  #plot GEA all
  
  locus_scores <- scores(RDA_all, choices=c(1:2), display="species", scaling="sites")
  TAB_loci <- data.frame(names = row.names(locus_scores), locus_scores)
  TAB_loci$type <- "Not associated"
  TAB_loci$type[TAB_loci$names%in%outliers$Loci] <- "FDR"
  TAB_loci$type[TAB_loci$names%in%top_outliers$Loci] <- "Bonferroni"
  TAB_loci$type <- factor(TAB_loci$type, levels = c("Not associated", "FDR", "Bonferroni"))
  TAB_var <- as.data.frame(scores(RDA_all, choices=c(1,2), display="bp"))
  loading_all<-ggplot() +
    geom_hline(yintercept=0, linetype="dashed", color = gray(.80), size=0.6) +
    geom_vline(xintercept=0, linetype="dashed", color = gray(.80), size=0.6) +
    geom_point(data = TAB_loci, aes(x=RDA1, y=RDA2, colour = type), size = 2.5) +
    scale_color_manual(values = c("gray90", "#F9A242FF", "#6B4596FF")) +
    geom_segment(data = TAB_var, aes(xend=RDA1, yend=RDA2, x=0, y=0), colour="black", size=0.15, linetype=1, arrow=arrow(length = unit(0.02, "npc"))) +
    geom_label_repel(data = TAB_var, aes(x=1.1*RDA1, y=1.1*RDA2, label = row.names(TAB_var)), size = 3.8, family = "Times") +
    xlab("RDA 1: 23%") + ylab("RDA 2: 12%") +
    guides(color=guide_legend(title="Locus type")) +
    theme_bw(base_size = 11, base_family = "Times") +
    theme(panel.background = element_blank(), legend.background = element_blank(), panel.grid = element_blank(), plot.background = element_blank(), legend.text=element_text(size=rel(.8)), strip.text = element_text(size=11))
  loading_all
  
  
  write.table(qvalue, "All_GEA_WWE_RDA.csv", append = FALSE, quote = TRUE, sep = " ",
              eol = "\n", na = "NA", dec = ".", row.names = FALSE,
              col.names = TRUE, qmethod = c("escape", "double"),
              fileEncoding = "")
  
  
  
  library(qqman)
  Manhattan_all <- read.csv(file = "All_GEA_WWE_RDA.csv", header=TRUE) #import the p value result for precipitation
 
  manhattan(Manhattan_all, col = c("darkgreen", "gray60"),suggestiveline = -log10(0.001227563), genomewideline = -log10(1.961360e-07))
 
  
  #P distribution
  jpeg(file = "/lustre/rocchettil/Phist_Manh_RDA_prec.jpeg")
  hist(Manhattan_all$P)
  dev.off()
  
  ############## 155 WWE GEA all together LFMM
  
  Y <- geno155
  sel_latent<- data.frame(Variables%>% dplyr::select(bio2, bio10, bio11, bio15, bio18, bio19,clay, N, pH, sand))
  write.env(sel_latent, "latent_all_variable.env")
  X = read.table("latent_all_variable.env")
  
  mod_lfmm <- LEA::lfmm2(input = Y, env = X, K = 3, effect.sizes = TRUE)
  #get environment effect sizes
  mod_lfmm@B
  #Define GEA
  pv = lfmm2.test(mod_lfmm, input = geno155, env = X, full = T)
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
  thres <- 0.05/ncol(geno155)
  
  PvaluesGEA_lfmm<-data.frame(pv$pvalues)
  
  #define cadidatemod.lfmm2#define cadidate loci for GO 
  
  GEA_lfmm <- data.frame(pvalue = pv$pvalues[pv$pvalue < thres])
  write.csv(GEA_lfmm, "GEA_lfmm_WWE_bonferroni.csv")
  write.csv(pv$pvalues, "GEA_lfmm_WWE.csv")
  
  #plotting Mhanattan plot using the library qqman
  
  library(qqman)
  Manhattan_all <- read.csv(file = "GEA_lfmm_WWE.csv", header=TRUE) #import the p value result for precipitation
  jpeg(file = "/lustre/rocchettil/Manh_RDA_prec.jpeg")
  manhattan(Manhattan_all, col = c("darkgreen", "gray60"),genomewideline = 6.705138)
  
  hist(Manhattan_all$P)
  dev.off()
  
  ### Save lfmm GEA results 
  GEA_lfmm_bonferroni<-read.csv(file='GEA_lfmm_WWE_bonferroni.csv', header = T)
  qtl_gea <- as.character(GEA_lfmm_bonferroni$QTL)
  GEA_WestEast_bonf<-dplyr::select(geno155, all_of(qtl_gea))
  GEA_lfmm_WWE_bonf<-write.table(GEA_WestEast_bonf, "GEA_lfmm_WWE_bonf.txt")
  GEA_lfmm_WWE_bonf<-read.table("GEA_lfmm_WWE_bonf.txt")
  ########################################     partial redundancy analysis (RDA only with GEA QTL)  #######################
  ## A total of 12143 GEA QTL (FDR threshold) have been identified and used for RDA
  ## A total of 438 GEA QTL (Bonferroni threshold) have been identified and used for RDA
  ## A total of 11161 GEA QTL (FDR threshold) identified with pRDA using all env variables (prdadapt used with 4 dimentions)
  ## A total of 336 GEA QTL (Bonferroni threshold) identified with pRDA using all env variables (prdadapt used with 4 dimentions)
  
  
  geno_Wild_GEA<-geno155[which((rdadapt_temp$p.values<thres_env)|(rdadapt_prec$p.values<thres_env)|(rdadapt_soil$p.values<thres_env))]
  
  geno_Wild_allGEA_RDA<-genoWW[which((rdadapt_all$p.values<thres_env))]
  geno_Wild_allGEA<-geno155[which((rdadapt_all$q.values<0.05))]
  
  write.table(geno_Wild_GEA, "geno_Wild_GEA_WWE.txt") #save the new GEA genotype data
  geno_Wild_GEA<-read.table("geno_Wild_GEA_WWE.txt")
  write.table(geno_Wild_GEA, "geno_Wild_GEA_Bonferroni_WWE.txt") #save the new GEA genotype data
  geno_Wild_GEA<-read.table("geno_Wild_GEA_Bonferroni_WWE.txt")
  write.table(geno_Wild_allGEA, "geno_Wild_GEA_AllVariable_Bonferroni_WWE.txt")
  geno_Wild_allGEA<-read.table("geno_Wild_GEA_AllVariable_Bonferroni_WWE.txt")
  write.table(geno_Wild_allGEA, "geno_Wild_GEA_AllVariable_FDR_WWE.txt")

  

  
  ### tentative selecting GEA that have variations among cultivars
  
  GEA_191<- geno_Wild_allGEA[, colnames(geno_Wild_allGEA)%in% colnames(GEA_cultivars_maf)]

  
  
  
  
  RDA_all_enriched<-rda(GEA_191 ~ bio2 + bio10 + bio11 + bio15	+ bio18 + bio19 + clay + N+ pH + sand , Variables)
  summary(eigenvals(RDA_all_enriched, model = "constrained"))
  RsquareAdj(RDA_all_enriched)
plot(RDA_all_enriched)
sqrt(vif.cca(RDA_all_enriched))

# plot Geographic regions


TAB_gen <- data.frame(geno = row.names(scores(RDA_all_enriched , display = "sites")), scores(RDA_all_enriched, display = "sites"))

Geno <- merge(TAB_gen, Variables[, 1:5] ,by="geno")
TAB_var <- as.data.frame(scores(RDA_all_enriched, choices=c(1,2), display="bp"))
loading_geno_all_enriched_region<-ggplot() +
  geom_hline(yintercept=0, linetype="dashed", color = gray(.80), size=0.6) +
  geom_vline(xintercept=0, linetype="dashed", color = gray(.80), size=0.6) +
  geom_point(data = Geno, aes(x=RDA1, y=RDA2, fill = region), size = 2.5, shape = 21, color = "black", stroke = 0.8) +
  scale_fill_manual(values = c("lightblue","darkgreen", "darkorange", "darkblue")) +
  geom_segment(data = TAB_var, aes(xend=RDA1, yend=RDA2, x=0, y=0), colour="black", linewidth =0.15, linetype=1, arrow=arrow(length = unit(0.02, "npc"))) +
  geom_label_repel(data = TAB_var, aes(x=RDA1, y=RDA2, label = row.names(TAB_var)), size = 3.2, family = "Times") +
  xlab("RDA 1: 56%") + ylab("RDA 2: 22%") +
  guides(color=guide_legend(title="Latitude gradient")) +
  theme_bw(base_size = 11, base_family = "Times") +
  theme(panel.background = element_blank(), legend.background = element_blank(), panel.grid = element_blank(), plot.background = element_blank(), legend.text=element_text(size=rel(.8)), strip.text = element_text(size=11))
loading_geno_all_enriched_region
jpeg(file = "/lustre/rocchettil/RDA_all_geno_biplot_region.jpeg")
plot(loading_geno_all_enriched_region)
dev.off()

### RDA plot based on GEA considering only West Mediterrenean
#filter west med individuals

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
Variables_142WW<-data.frame(geno=dataWild_142W$id,group = dataWild_142W$group, region = dataWild_142W$region, Env )



geno_wild_GEA_142WW<- geno_Wild_allGEA[rownames(geno_Wild_allGEA)%in% list142WW$V1, ]

RDA_142WW_enriched<-rda(geno_wild_GEA_142WW ~ bio2 + bio10 + bio11 + bio15	+ bio18 + bio19 +clay + N + pH +  sand, Variables_142WW)
plot(RDA_142WW_enriched)
summary(eigenvals(RDA_142WW_enriched, model = "constrained"))
plot(eigenvals(RDA_142WW_enriched, model = "constrained"))

TAB_gen <- data.frame(geno = row.names(scores(RDA_142WW_enriched , display = "sites")), scores(RDA_142WW_enriched, choices=c(1,2,3), display = "sites"))

Geno <- merge(TAB_gen, Variables_142WW[, 1:3] ,by="geno")
TAB_var <- as.data.frame(scores(RDA_142WW_enriched, choices=c(1,2), display="bp"))
loading_geno_142W_enriched_region<-ggplot() +
  geom_hline(yintercept=0, linetype="dashed", color = gray(.80), size=0.6) +
  geom_vline(xintercept=0, linetype="dashed", color = gray(.80), size=0.6) +
  geom_point(data = Geno, aes(x=RDA1, y=RDA2, fill = region), size = 2.5, shape = 21, color = "black", stroke = 0.8) +
  scale_fill_manual(values = c("lightblue","darkgreen", "darkorange")) +
  geom_segment(data = TAB_var, aes(xend=RDA1*5, yend=RDA2*5, x=0, y=0), colour="black", linewidth =0.15, linetype=1, arrow=arrow(length = unit(0.02, "npc"))) +
  geom_label_repel(data = TAB_var, aes(x=5*RDA1, y=5*RDA2, label = row.names(TAB_var)), size = 3.2, family = "Times") +
  xlab("RDA 1: 69.7%") + ylab("RDA 2: 12.2%") +
  guides(color=guide_legend(title="Latitude gradient")) +
  theme_bw(base_size = 11, base_family = "Times") +
  theme(panel.background = element_blank(), legend.background = element_blank(), panel.grid = element_blank(), plot.background = element_blank(), legend.text=element_text(size=rel(.8)), strip.text = element_text(size=11))
loading_geno_142W_enriched_region

ggarrange(loading_geno_all_enriched_region, loading_geno_142W_enriched_region, nrow=1,ncol=2)



### Model of the adaptive Landscape
# estimation of adaptive value of each environmental pixel

library(raster)
library("readxl")


bio2<- raster(paste("D:/raster files/Current_ENM_clipped_biova/bio2_current_masked.tif"))
bio10<- raster(paste("D:/raster files/Current_ENM_clipped_biova/bio10_current_masked.tif"))
bio11<- raster(paste("D:/raster files/Current_ENM_clipped_biova/bio11_current_masked.tif"))
bio15<- raster(paste("D:/raster files/Current_ENM_clipped_biova/bio15_current_masked.tif"))
bio18<- raster(paste("D:/raster files/Current_ENM_clipped_biova/bio18_current_masked.tif"))
bio19<- raster(paste("D:/raster files/Current_ENM_clipped_biova/bio19_current_masked.tif"))
soilN<- raster(paste("D:/raster files/Current_ENM_clipped_biova/resampled_soilN.tif"))
soilpH<- raster(paste("D:/raster files/Current_ENM_clipped_biova/resampled_soilpH.tif"))

soilclay<- raster(paste("D:/raster files/Current_ENM_clipped_biova/resampled_soilclay.tif"))
soilsand<- raster(paste("D:/raster files/Current_ENM_clipped_biova/resampled_soilsand.tif"))

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


#alignment of soil rasters with bioclimatic variables

soilN <- resample(soilN, bio2, method="bilinear")
writeRaster(soilN, "resampled_soilN.tif", format="GTiff", overwrite=TRUE)


soilpH<- resample(soilpH, bio2, method="bilinear")
writeRaster(soilpH, "D:/raster files/Current_ENM_clipped_biova/resampled_soilpH.tif", format="GTiff", overwrite=TRUE)

soilSOC<- resample(soilSOC, bio2, method="bilinear")
writeRaster(soilSOC, "D:/raster files/Current_ENM_clipped_biova/resampled_soilSOC.tif", format="GTiff", overwrite=TRUE)

soilclay <- resample(soilclay, bio2, method="bilinear")
writeRaster(soilclay, "D:/raster files/Current_ENM_clipped_biova/resampled_soilclay.tif", format="GTiff", overwrite=TRUE)

soilsand <- resample(soilsand, bio2, method="bilinear")
writeRaster(soilsand, "D:/raster files/Current_ENM_clipped_biova/resampled_soilsand_1.tif", format="GTiff", overwrite=TRUE)


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

#plot
#create color palette
a1 <- TAB_pixel_LC$RDA1
a2 <- TAB_pixel_LC$RDA2
a3 <- TAB_pixel_LC$RDA3


# Colorblind-friendly RGB calculation based on the Color Universal Design (CUD) palette
# CUD palette uses blue, orange, green, purple, and red tones
r <- (a1 + 1) * 30  # Linear scaling of RDA1 with an offset for better contrast
g <- abs(a2) * 60  # Scaling absolute RDA2 with a factor for visible green shades
b <- (a3^2 + a1^2)^(1/2) * 40  # Euclidean distance of RDA1 and RDA3 for blue intensity

# Normalize RGB values to range 0-255
r <- (r - min(r)) / (max(r) - min(r)) * 255
g <- (g - min(g)) / (max(g) - min(g)) * 255
b <- (b - min(b)) / (max(b) - min(b)) * 255



# Combine into hexadecimal color
TAB_pixel_LC$color <- rgb(r, g, b, maxColorValue = 255)

# Update ggplot to use custom RGB colors
pp <- ggplot() +
  geom_hline(yintercept = 0, linetype = "dashed", color = gray(0.80), size = 0.6) +
  geom_vline(xintercept = 0, linetype = "dashed", color = gray(0.80), size = 0.6) +
  geom_point(data = TAB_pixel_LC, aes(x = RDA1, y = RDA2), color = TAB_pixel_LC$color) +  # Use custom colors
  geom_segment(data = TAB_var, aes(xend = RDA1, yend = RDA2, x = 0, y = 0), 
               colour = "black", size = 0.15, linetype = 1, 
               arrow = arrow(length = unit(0.20, "cm"), type = "closed")) +
  geom_label_repel(data = TAB_var, aes(x = RDA1, y = RDA2, label = row.names(TAB_var)), 
                   size = 4, family = "Times") +
  xlab("RDA 1: 69.7%") + 
  ylab("RDA 2: 12.2%") +
  theme_bw(base_size = 9, base_family = "Times") +
  theme(panel.background = element_blank(), 
        legend.background = element_blank(), 
        panel.grid = element_blank(), 
        plot.background = element_blank(), 
        legend.text = element_text(size = rel(0.8)), 
        strip.text = element_text(size = 9))
pp




# Extract RDA values
a1 <- TAB_pixel_LC$RDA1
a2 <- TAB_pixel_LC$RDA2

# Compute the distance from the origin
distance <- sqrt(a1^2 + a2^2)

# Assign colors based on quadrants and the 5th sector (circle radius < 0.5)
TAB_pixel_LC$color <- ifelse(distance < 0.25, "#717171",  # 5th sector - Purple
                             ifelse(a1 > 0 & a2 > 0, "#FF7F00",  # Quadrant 1 - Red
                                    ifelse(a1 < 0 & a2 > 0, "limegreen",  # Quadrant 2 - Blue
                                           ifelse(a1 < 0 & a2 < 0, "#377EB8",  # Quadrant 3 - Green
                                                  "red"))))  # Quadrant 4 - Orange

# Update ggplot with quadrant-based colors and 5th sector
pp <- ggplot() +
  geom_hline(yintercept = 0, linetype = "dashed", color = gray(0.80), size = 0.6) +
  geom_vline(xintercept = 0, linetype = "dashed", color = gray(0.80), size = 0.6) +
  geom_point(data = TAB_pixel_LC, aes(x = RDA1, y = RDA2, color = color), size = 2) +  # Use sector colors
  geom_segment(data = TAB_var, aes(xend = RDA1, yend = RDA2, x = 0, y = 0), 
               colour = "black", size = 0.15, linetype = 1, 
               arrow = arrow(length = unit(0.20, "cm"), type = "closed")) +
  geom_label_repel(data = TAB_var, aes(x = RDA1, y = RDA2, label = row.names(TAB_var)), 
                   size = 4, family = "Times") +
  xlab("RDA 1: 61.6%") + 
  ylab("RDA 2: 19.1%") +
  theme_bw(base_size = 9, base_family = "Times") +
  theme(panel.background = element_blank(), 
        legend.background = element_blank(), 
        panel.grid = element_blank(), 
        plot.background = element_blank(), 
        legend.text = element_text(size = rel(0.8)), 
        strip.text = element_text(size = 9)) +
  scale_color_identity()  # Use predefined colors directly

pp

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
TAB_pixel_LC_sf <- st_as_sf(TAB_pixel_LC, coords = c("long", "lat"), crs = 4326)  # Assumes 'longitude' and 'latitude' columns exist
# Create the map
map <- ggplot(data = countries) +
  geom_sf(fill = "#EBEBEB", color = "black") +  # Countries' borders
  geom_sf(data = TAB_pixel_LC_sf, aes(color = color), size = 0.05, show.legend = FALSE) +  # Points with custom colors
  scale_color_identity() +  # Use exact colors from the 'color' column
  coord_sf(xlim = c(-15, 15), ylim = c(28, 52), expand = FALSE) +  # Set geographic limits
  theme_minimal() +
  labs(title = "Adaptive Landscape") +
  theme(panel.background = element_blank())

#jpeg(file = "C:/Users/rocchetti/Desktop/running RDA GO/adaptive_landscape_rgb.jpeg",width = 18, height = 14, units = "cm", res = 800)
map
dev.off()


######################################################### Estimation of cultivar offset #######################

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

GEA <-colnames(GEA_lfmm_bonf)
GEA_geno_cultivar<-dplyr::select(geno_cultivar, all_of(GEA))

#imputation
for (i in 1:ncol(GEA_geno_cultivar))
{
  GEA_geno_cultivar[which(is.na(GEA_geno_cultivar[,i])),i] <- median(GEA_geno_cultivar[-which(is.na(GEA_geno_cultivar[,i])),i], na.rm=TRUE)
}

 ## save GEA all varible
write.table(GEA_geno_cultivar, "GEA_lfmm_Bonferroni_all_cultivars.txt")


GEA_cultivars<-read.table("GEA_all_Variable_Bonferroni_all_cultivars.txt")

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

##### cultivar offset

###### Frantoio

F <- GEA_cultivars[rownames(GEA_cultivars) == "Berri_Meslal-3971", ]

FRDA <- predict(RDA_all_enriched, newdata=F, type="wa")
FRDA<-as.data.frame(FRDA)



TAB_pixel_LC$offset <- (FRDA$RDA1 - TAB_pixel_LC$RDA1)^2 * 0.5687 + 
  (FRDA$RDA2 - TAB_pixel_LC$RDA2)^2 * 0.2218 + 
  (FRDA$RDA3 - TAB_pixel_LC$RDA3)^2 * 0.07756 
  
hist(TAB_pixel_LC$offset)


#TAB_pixel_LC$Frantoio_offset<-scale(TAB_pixel_LC$Frantoio_offset, center = intersection_x, scale = sd(results_df$GO))
#hist(TAB_pixel_LC$Frantoio_offset)

TAB_pixel_LC$offset_percentile <- ecdf(TAB_pixel_LC$offset)(TAB_pixel_LC$offset) * 100

# Re-convert to sf object with updated column
TAB_pixel_LC_sf <- st_as_sf(TAB_pixel_LC, coords = c("long", "lat"), crs = 4326)

# Define continuous palette function
continuous_palette <- colorRampPalette(c("#004d00",  # very dark green
                                         "#228B22",  # forest green
                                         "#66C200",  # yellow-green
                                         "#CCCC00",  # mustard yellow
                                         "#FFD700",  # golden yellow
                                         "#FFA500",  # orange
                                         "#FF8C00",  # dark orange
                                         "#FF4500",  # orange-red
                                         "#B22222",  # firebrick
                                         "#8B0000"))

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
TAB_pixel_LC_sf <- st_as_sf(TAB_pixel_LC, coords = c("long", "lat"), crs = 4326)
# Create the map

# Create the map with continuous legend
map <- ggplot(data = countries) +
  geom_sf(fill = "#EBEBEB", color = "black") +
  geom_sf(data = TAB_pixel_LC_sf, aes(color = offset_percentile), size = 0.05) +
  scale_color_gradientn(
    colors = continuous_palette(100),
    limits = c(0, 100),
    name = "GO Percentile"
  ) +
  coord_sf(xlim = c(-15, 15), ylim = c(28, 52), expand = FALSE) +
  theme_minimal() +
  labs(title = "Adaptive Landscape Berri_Meslal-3971") +
  theme(
    panel.background = element_blank(),
    legend.position = "right",
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 8)
  )
map


######## Offset at marrakesh

#### Comparison with data laila marakesh


# Define the target latitude and longitude
target_lat <- 31.80403
target_long <- -7.254303

# Define a tolerance value (small range of acceptable difference)
tolerance <- 1e-5

# Filter the data using a tolerance for matching
Marrakech_row <- TAB_pixel_LC %>%
  filter(abs(lat - target_lat) < tolerance & abs(long - target_long) < tolerance)

Marrakech_row<-as.data.frame(Marrakech_row)
### offset at marakesh


Marrakech_row <- Marrakech_row[1, ]

# Perform the calculation
RDAscore_cul$offsetM <- (Marrakech_row$RDA1 - RDAscore_cul$RDA1)^2 * 0.5687 + 
  (Marrakech_row$RDA2 - RDAscore_cul$RDA2)^2 * 0.2218 + 
  (Marrakech_row$RDA3 - RDAscore_cul$RDA3)^2 * 0.07756 

write.table(RDAscore_cul, "offset_marakkesh_wildE_RDAGEA.txt")

GoMar<-read.csv("GO_marakesh.csv")
GoMar <- na.omit(GoMar)

boxplot(GOM ~ class, data = GoMar)

model <- lm(GOM ~ class, data = GoMar)

summary(model)

# Perform Tukey's Honest Significant Difference test
tukey_result <- TukeyHSD(aov(model))

# View results
print(tukey_result)

pairwise.t.test(GoMar$GOM, GoMar$class, p.adjust.method = "bonferroni")

plot(GO_We ~ CHL, data = GoMar)

model <- lm(GOM ~ CHL, data = GoMar)
abline(model, col = "red", lwd = 2)
summary(model)


###### plot predicted cultivar and wild genotypes in the RDA space



Tab_cultivar<- data.frame(geno = row.names(RDAscore_cul),RDAscore_cul[,1:2] )
Tab_cultivar$group<-"cultivars"

RDA_genowild<-predict(RDA_all_enriched,newdata =  GEA_191, type = "wa")
Tab_genowild<-data.frame(geno = Variables[,1],RDA_genowild[,1:2])
Tab_genowild$group<-"wild_geno"
RDA_wild<-predict(RDA_all_enriched,newdata =  Variables[, 6:15], type = "lc")
Tab_wild<-data.frame(geno = Variables[,1],RDA_wild[,1:2])
Tab_wild$group<-"wild_env"


TAB_var <- as.data.frame(scores(RDA_all_enriched, choices=c(1,2), display="bp"))

wild_cult_pred<-rbind(Tab_cultivar, Tab_genowild)
wild_cult_pred$group <- as.factor(wild_cult_pred$group)

hh <- ggplot() +
  geom_hline(yintercept=0, linetype="dashed", color = gray(.80), size=0.6) +
  geom_vline(xintercept=0, linetype="dashed", color = gray(.80), size=0.6) +
  geom_point(data = wild_cult_pred, aes(x = RDA1, y = RDA2, fill = group, shape = group),size = 2.5, color = "black", stroke = 0.8)+
  scale_shape_manual(values = c(24,21))+
  scale_fill_manual(values=c('#E69F00',"grey48"))+
  #scale_size_manual(values=c(3,3))+
  #geom_segment(data = TAB_var, aes(xend=RDA1*0.5, yend=RDA2*0.5, x=0, y=0), colour="black", size=0.15, linetype=1, arrow = arrow(length=unit(0.20,"cm"),type = "closed")) +
  #geom_label_repel(data = TAB_var, aes(x=RDA1*0.5, y=RDA2*0.5, label = row.names(TAB_var)), size = 3, family = "Times") +
  xlab("RDA 1: 56%") + ylab("RDA 2: 22%") +
  #guides(legend(title="Group")) +
  theme_bw(base_size = 11, base_family = "Times") +
  theme(panel.background = element_blank(), legend.background = element_blank(), panel.grid = element_blank(), plot.background = element_blank(), legend.text=element_text(size=rel(.8)), strip.text = element_text(size=11))+
labs(title = "RDA Wild and Cultivar with GEAs")
hh


#### PCA GEA cultivar and Wild
### genetic variation of GEA QTLs beween wild and cultivated


wildwest<-data.frame(geno = row.names(geno_Wild_allGEA[1:142, ]), group = "Wildwest", geno_Wild_allGEA[1:142,1:232])
cultivars<-data.frame(geno = row.names(GEA_cultivars), group = "Cultivar", GEA_cultivars[,1:232])
wildeast<- data.frame(geno = row.names(geno_Wild_allGEA[143:155,]), group = "Wildeast", geno_Wild_allGEA[143:155,1:232])
df<-rbind(wildwest, wildeast, cultivars)
library(FactoMineR)
library(factoextra)
res.pca_df<-PCA(df[,3:234], scale.unit = TRUE, ncp = 5, graph = TRUE)
inddf <- get_pca_ind(res.pca_df)
pca_data_df <- as.data.frame(inddf$coord)
pca_data_df<-cbind(df[,1:2], pca_data_df)
qq<-ggplot() +
  geom_hline(yintercept=0, linetype="dashed", color = gray(.80), linewidth=0.6) +
  geom_vline(xintercept=0, linetype="dashed", color = gray(.80), linewidth=0.6) +
  geom_point(data = pca_data_df, aes(x=Dim.1, y=Dim.2, fill = group, shape = group),size = 2.5, color = "black", stroke = 0.8) +
  scale_shape_manual(values = c(24,21,21))+
  scale_fill_manual(values=c('#E69F00',"darkblue","grey48"))+
  xlab("PC1: 23.9%") + ylab("PC2: 12.0%") +
  guides(color=guide_legend(title="Group")) +
  theme_bw(base_size = 11, base_family = "Times") +
  theme(panel.background = element_blank(), legend.background = element_blank(), panel.grid = element_blank(), plot.background = element_blank(), legend.text=element_text(size=rel(.8)), strip.text = element_text(size=11))+
  labs(title = "PCA Wild and Cultivar with GEAs")
qq


library(ggpubr)
ggarrange(hh,qq,nrow=1, ncol=2)



C2_GEA <- Tab_cultivar[rownames(Tab_cultivar) %in% c("Olivastra_di_Populonia1", "Sivigliana_da_Olio1", "Lazzero_di_prata1","Lastrino1", "Ottobratica1", "Tabelout1","Berri_Meslal-3971", "Vaddarica1"), ]
wild_cult_C2<-rbind(Tab_wild, C2_GEA)
aa <- ggplot() +
  geom_hline(yintercept=0, linetype="dashed", color = gray(.80), size=0.6) +
  geom_vline(xintercept=0, linetype="dashed", color = gray(.80), size=0.6) +
  geom_point(data = wild_cult_C2, aes(x = RDA1, y = RDA2, fill = type, shape = type),size = 2.5, color = "black", stroke = 0.8)+
  scale_shape_manual(values = c(24,21))+
  scale_fill_manual(values=c('#E69F00','grey48'))+
  #scale_size_manual(values=c(3,3))+
  geom_segment(data = TAB_var, aes(xend=RDA1, yend=RDA2, x=0, y=0), colour="black", size=0.15, linetype=1, arrow = arrow(length=unit(0.20,"cm"),type = "closed")) +
  geom_label_repel(data = TAB_var, aes(x=RDA1, y=RDA2, label = row.names(TAB_var)), size = 3, family = "Times") +
  geom_label_repel(data = C2_GEA, aes(x=RDA1, y=RDA2, label = row.names(C2_GEA)), size = 3, family = "Times") +
  xlab("RDA 1: 69.7%") + ylab("RDA 2: 12.2%") +
  #guides(legend(title="Group")) +
  theme_bw(base_size = 11, base_family = "Times") +
  theme(panel.background = element_blank(), legend.background = element_blank(), panel.grid = element_blank(), plot.background = element_blank(), legend.text=element_text(size=rel(.8)), strip.text = element_text(size=11))+
  labs(title = "RDA Wild and Cultivar with GEAs")
aa


#PCA Cultivars with GEA

library(FactoMineR)
library(factoextra)

res.pca<-PCA(GEA_geno_3_cul, scale.unit = TRUE, ncp = 5, graph = TRUE)
fviz_eig(res.pca, addlabels = TRUE, ylim = c(0, 50))#  scree plot
library(tibble)
ind <- get_pca_ind(res.pca)
pca_data <- as.data.frame(ind$coord)
eig.val <- get_eigenvalue(res.pca) #selected 5 PCs
# Create a data frame for PCA results
ind135 <- get_pca_ind(res.pca)
ind135
pca_data135 <- as.data.frame(ind135$coord)
qq<-ggplot() +
  geom_hline(yintercept=0, linetype="dashed", color = gray(.80), size=0.6) +
  geom_vline(xintercept=0, linetype="dashed", color = gray(.80), size=0.6) +
  geom_point(data = pca_data135, aes(x=Dim.1, y=Dim.2), size = 2.5) +
  geom_label_repel(data = pca_data135, aes(x=Dim.1, y=Dim.2, label = row.names(pca_data135)), size = 3.2, family = "Times") +
    xlab("PC1: 28.1%") + ylab("PC2: 24.8%") +
  guides(color=guide_legend(title="Group")) +
  theme_bw(base_size = 11, base_family = "Times") +
  theme(panel.background = element_blank(), legend.background = element_blank(), panel.grid = element_blank(), plot.background = element_blank(), legend.text=element_text(size=rel(.8)), strip.text = element_text(size=11))
qq



####### Loop cycle for GO estimation all cultivars


# empty dataframe
results_df <- data.frame(lat = numeric(0), long = numeric(0), 
                         RDA1 = numeric(0), RDA2 = numeric(0), 
                         RDA3 = numeric(0), GO = numeric(0), row_name = character(0))

for (i in row.names(RDAscore_cul)) {
  GO <- data.frame(
    lat = TAB_pixel_LC$lat, 
    long = TAB_pixel_LC$long, 
    RDA1 = TAB_pixel_LC$RDA1, 
    RDA2 = TAB_pixel_LC$RDA2, 
    RDA3 = TAB_pixel_LC$RDA3, 
    GO = (RDAscore_cul[i,1] - TAB_pixel_LC$RDA1)^2 + 
      (RDAscore_cul[i,2] - TAB_pixel_LC$RDA2)^2 + 
      (RDAscore_cul[i,3] - TAB_pixel_LC$RDA3)^2
  )
  
  GO$row_name <- i  
  results_df <- dplyr::bind_rows(results_df, GO) 
}

write.table(results_df, "D:/Leccino24_genome/cultivar_pop_geom_RDAGO.txt")
results_df<-read.table("D:/Leccino24_genome/cultivar_pop_geom_RDAGO.txt")
hist(results_df$GO, main = "geometric GO distribution cultivar population", las = 1, yaxt = "n")



TAB_pixel_LC$arbequina_offset<- (arbequina_RDAscore$RDA1 - TAB_pixel_LC$RDA1)^2 + (arbequina_RDAscore$RDA2 - TAB_pixel_LC$RDA2)^2 + (arbequina_RDAscore$RDA3 - TAB_pixel_LC$RDA3)^2
hist(TAB_pixel_LC$arbequina_offset)
install.packages("fitdistrplus")
library(fitdistrplus)
# Estimate an initial guess for the degrees of freedom (mean of the data)
initial_df <- mean(results_df$GO)

# Fit a chi-square distribution with the manually specified starting value for df
fit <- fitdist(results_df$GO, "chisq", start = list(df = initial_df))

hist(TAB_pixel_LC$arbequina_offset, main = "Histogram of Arbequina Offset", las = 1, 
     yaxt = "n", col = "lightgray", border = "darkgray", probability = TRUE)


# Generate values for plotting the fitted chi-square PDF
x_values <- seq(min(TAB_pixel_LC$arbequina_offset), max(TAB_pixel_LC$arbequina_offset), length.out = 1000)
pdf_values <- dchisq(x_values, df = fit$estimate["df"])

# Overlay the fitted PDF on the histogram
lines(x_values, pdf_values, col = "blue", lwd = 2)

# Estimate and plot the CDF of the fitted chi-square distribution
cdf_values <- pchisq(x_values, df = fit$estimate["df"])
lines(x_values, cdf_values, col = "red", lwd = 2)


# Find the x-value where PDF and CDF intersect
intersection_index <- which.min(abs(pdf_values - cdf_values))  # Find closest match
intersection_x <- x_values[intersection_index]

cat("The intersection point is approximately at x =", intersection_x, "\n")

# Add a vertical line at the intersection point for visualization
abline(v = intersection_x, col = "purple", lty = 2, lwd = 2)
legend("topright", legend = c("PDF", "CDF", "Intersection"), 
       col = c("blue", "red", "purple"), lty = c(1, 1, 2), lwd = c(2, 2, 2))





















culitivar_pop <- scale(results_df$GO, center=TRUE, scale=TRUE)
# Extract the centering values
culitivar_pop_center <- attr(culitivar_pop, "scaled:center")
# Extract the scaling values
culitivar_pop_scale <- attr(culitivar_pop, "scaled:scale")

#save parameters
write.csv(data.frame(center = culitivar_pop_center, scale = culitivar_pop_scale), "culitivar_scaling.csv", row.names = FALSE)
scaling_params <- read.csv("culitivar_scaling.csv")
culitivar_pop_center <- scaling_params$center
culitivar_pop_scale <- scaling_params$scale







####### Loop cycle for GO estimation 142 reference wild west


TAB_gen <- data.frame(geno = row.names(scores(RDA_142WW_enriched , display = "sites")), scores(RDA_142WW_enriched, choices=c(1,2,3), display = "sites", scaling = "sites"))

# empty dataframe
results_df <- data.frame(lat = numeric(0), long = numeric(0),   
                         RDA1 = numeric(0), RDA2 = numeric(0), 
                         RDA3 = numeric(0), GO = numeric(0))

for (i in 1:142) {
  GO <- data.frame(
    lat = TAB_pixel_LC$lat, 
    long = TAB_pixel_LC$long, 
    RDA1 = TAB_pixel_LC$RDA1, 
    RDA2 = TAB_pixel_LC$RDA2, 
    RDA3 = TAB_pixel_LC$RDA3, 
    GO = (TAB_gen[i,2] - TAB_pixel_LC$RDA1)^2 + 
      (TAB_gen[i,3] - TAB_pixel_LC$RDA2)^2 + 
      (TAB_gen[i,4] - TAB_pixel_LC$RDA3)^2
  )
  
  results_df <- dplyr::bind_rows(results_df, GO) 
}

hist(results_df$GO, main = "geometric GO distribution wild population", las = 1, yaxt = "n")



#save parameters
write.csv(data.frame(intersection = intersection_x, sdev = sd(results_df$GO)), "wild_pop_geomGOthreshold_param.csv", row.names = FALSE)
params <- read.csv("wild_pop_geomGOthreshold_param.csv")

#####Fit Chisquare PDF

library(fitdistrplus)

initial_df <- mean(results_df$GO)

# Fit a chi-square distribution with the manually specified starting value for df
fit <- fitdist(results_df$GO, "chisq", start = list(df = initial_df))

hist(results_df$GO, main = "wild geometric GO", las = 1, 
     yaxt = "n", col = "lightgray", border = "darkgray", probability = TRUE)
axis(2, at = seq(0, 1, by = 0.1))

# Generate x values for plotting the fitted chi-square PDF
x_values <- seq(min(results_df$GO), max(results_df$GO), length.out = 100)
pdf_values <- dchisq(x_values, df = fit$estimate["df"])

# plot the fitted PDF on the histogram
lines(x_values, pdf_values, col = "blue", lwd = 2)

# Estimate and plot the CDF of the fitted chi-square distribution
cdf_values <- pchisq(x_values, df = fit$estimate["df"])
lines(x_values, cdf_values, col = "red", lwd = 2)


# Find the x-value where PDF and CDF intersect
intersection_index <- which.min(abs(pdf_values - cdf_values))  # Find closest match
intersection_x <- x_values[intersection_index]
intersection = intersection_x
cat("The intersection point is approximately at x =", intersection_x, "\n")

# Add a vertical line at the intersection point for visualization
abline(v = intersection_x, col = "purple", lty = 2, lwd = 2)
legend("topright", legend = c("PDF", "CDF", "Intersection"), 
       col = c("blue", "red", "purple"), lty = c(1, 1, 2), lwd = c(2, 2, 2))















wild142_pop <- scale(results_df$GO, center=TRUE, scale=TRUE)
# Extract the centering values
wild142_pop_center <- attr(wild142_pop, "scaled:center")
# Extract the scaling values
wild142_pop_scale <- attr(wild142_pop, "scaled:scale")

#save parameters
write.csv(data.frame(center = culitivar_pop_center, scale = culitivar_pop_scale), "culitivar_scaling.csv", row.names = FALSE)
scaling_params <- read.csv("culitivar_scaling.csv")
culitivar_pop_center <- scaling_params$center
culitivar_pop_scale <- scaling_params$scale







library(raster)
library("readxl")


bio2<- raster(paste("D:/raster files/Env_var_west_med/bio2_west_med.tif"))
bio10<- raster(paste("D:/raster files/Current_ENM_clipped_biova/bio10_current_masked.tif"))
bio11<- raster(paste("D:/raster files/Current_ENM_clipped_biova/bio11_current_masked.tif"))
bio15<- raster(paste("D:/raster files/Current_ENM_clipped_biova/bio15_current_masked.tif"))
bio18<- raster(paste("D:/raster files/Current_ENM_clipped_biova/bio18_current_masked.tif"))
bio19<- raster(paste("D:/raster files/Current_ENM_clipped_biova/bio19_current_masked.tif"))
soilN<- raster(paste("D:/raster files/Current_ENM_clipped_biova/resampled_soilN.tif"))
soilpH<- raster(paste("D:/raster files/Current_ENM_clipped_biova/resampled_soilpH.tif"))

soilclay<- raster(paste("D:/raster files/soil_clay_5_15_west_med.tif"))
soilsand<- raster(paste("D:/raster files/Current_ENM_clipped_biova/resampled_soilsand.tif"))

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


#alignment of soil rasters with bioclimatic variables

soilN <- resample(soilN, bio2, method="bilinear")
writeRaster(soilN, "resampled_soilN.tif", format="GTiff", overwrite=TRUE)


soilpH<- resample(soilpH, bio2, method="bilinear")
writeRaster(soilpH, "D:/raster files/Current_ENM_clipped_biova/resampled_soilpH.tif", format="GTiff", overwrite=TRUE)

soilSOC<- resample(soilSOC, bio2, method="bilinear")
writeRaster(soilSOC, "D:/raster files/Current_ENM_clipped_biova/resampled_soilSOC.tif", format="GTiff", overwrite=TRUE)

soilclay <- resample(soilclay, bio2, method="bilinear")
writeRaster(soilclay, "D:/raster files/Env_var_west_med/resampled_west_med_soilclay.tif", format="GTiff", overwrite=TRUE)

soilsand <- resample(soilsand, bio2, method="bilinear")
writeRaster(soilclay, "D:/raster files/Current_ENM_clipped_biova/resampled_soilsand.tif", format="GTiff", overwrite=TRUE)


#stack the different raster file
ras_current_var<-stack(c(bio2,bio10, bio11, bio15, bio18, bio19, soilclay,soilN,soilpH, soilsand))







####  run LFMM GEA


## Use latent factor for covariable correction
# latent factor temperature variable
Y <- genoWW_maf
sel_latent<- data.frame(Variables_142WW%>% dplyr::select(bio2, bio10, bio11, bio15, bio18, bio19,clay, N, pH, sand))
write.env(sel_latent, "latent_all_variable.env")
X = read.table("latent_all_variable.env")

mod.lfmm2 <- lfmm2(input = Y, env = X, K = 3, effect.sizes = TRUE)
#get environment effect sizes
mod.lfmm2@B
#Define GEA
pv = lfmm2.test(mod.lfmm2, input = Y, env = X, full = T)
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

PvaluesGEA_lfmm<-data.frame(pv$pvalues)

#define cadidatemod.lfmm2#define cadidate loci for GO 

GEA_lfmm <- data.frame(pvalue = pv$pvalues[pv$pvalue < thres])
write.csv(GEA_lfmm, "GEA_lfmm_all_var_log5.csv")
write.table(row.names(GEA_lfmm), "GEA_lfmm_all_var_Bonferroni.txt")
hist(pv$pvalues)
top_outliers<- data.frame(candidates)
write.csv(pv$pvalues, "GEA_all_var_lfmm.csv")

#plotting Mhanattan plot using the library qqman

library(qqman)
Manhattan_all <- read.csv(file = "GEA_all_var_lfmm.csv", header=TRUE) #import the p value result for precipitation
jpeg(file = "/lustre/rocchettil/Manh_RDA_prec.jpeg")
manhattan(Manhattan_all, col = c("darkgreen", "gray60"),genomewideline = 5)

hist(Manhattan_all$P)
dev.off()



######### represent GEA lfmm QTL in mhanattan plot derived from pRDA method
Manhattan_all <- read.csv(file = "All_GEA_WW_3.csv", header=TRUE)

LFMM_snps <- read.csv(file = "GEA_lfmm_all_var_log5.csv", header=TRUE) 


manhattan(Manhattan_all, 
          col = c("gray", "black"),  
          highlight = LFMM_snps$SNP,  # Highlight specific SNPs
          suggestiveline = FALSE,  
          genomewideline = -log10(2.189362e-07),  
          main = "RDA_all GEA highlighting LFMM GEA SNPs", 
          ylab = "-log10(P-value)",
          xlab = "Chromosome")


################################################################################# Enriched analysis with LFMM
GEA_lfmm_list<- read.table("GEA_lfmm_all_var_log5 (2).txt", header = T)
GEA_lfmm_all_var<-  genoWW_maf[, colnames(genoWW_maf)%in% GEA_lfmm_list$SNP]
write.table(GEA_lfmm_all_var, "GEA_lfmm_all_var.txt")
GEA_lfmm_all_var<-read.table("GEA_lfmm_all_var.txt")


#GEA cultivar


GEA <- colnames(GEA_lfmm_all_var)
GEA_geno_cultivar<-dplyr::select(geno_cultivar, all_of(GEA))
GEA_geno_cultivar_clean <- GEA_geno_cultivar[complete.cases(GEA_geno_cultivar), ]

write.table(GEA_geno_cultivar_clean, "GEA_allWW_lfmm_all_cultivars.txt")


GEA_cultivars<-read.table("GEA_allWW_lfmm_all_cultivars.txt")

#### Linear model (no latent factor) + PCA
Y <- GEA_124
sel_env<- data.frame(Variables_142WW[, 4:13])
write.env(sel_env, "env_variable.env")
X = read.table("env_variable.env")

# Ensure X is a numeric matrix
X <- as.matrix(sel_env)
storage.mode(X) <- "numeric"

# Initialize matrix to store effect sizes (slopes)
slopes <- matrix(NA, ncol = ncol(Y), nrow = ncol(X))
rownames(slopes) <- colnames(X)
colnames(slopes) <- colnames(Y)

# Loop through each SNP and fit a linear model
for (j in 1:ncol(Y)) {
  geno <- Y[, j]
  if (var(geno, na.rm = TRUE) == 0) next  # skip invariant SNPs
  model <- lm(geno ~ X)
  slopes[, j] <- summary(model)$coefficients[-1, 1]  # extract slopes (exclude intercept)
}

# Transpose: rows = SNPs, columns = environmental variables
slopes_t <- t(slopes)


Ypred_linear <- X %*% t(slopes_t)
Ypred_linear <- tcrossprod(X, slopes_t)

pc_linear <- prcomp(Ypred_linear)
# Calculate the proportion of variance explained
eigenvalues <- pc_linear$sdev^2
variance_explained <- eigenvalues / sum(eigenvalues)


pc_wild<-as.data.frame(pc_linear$x)
plot(pc_linear$x, col = "darkgrey", pch = 19, main = "PCA of Predicted Genotypes")


# Project original environmental variables into PCA space
# This gives loadings of ENVIRONMENTAL variables
env_arrows <- cor(X, pc_linear$x)

# Plot PCA of individuals
plot(pc_linear$x, col = "gray", pch = 19,
     xlab = "PC1", ylab = "PC2", main = "linear model + PCA")

# Add environmental vectors
arrow_scale <- 5  # tweak for visibility
arrows(0, 0, env_arrows[,1]*arrow_scale, env_arrows[,2]*arrow_scale, col = "red", length = 0.1)
text(env_arrows[,1]*arrow_scale, env_arrows[,2]*arrow_scale,
     labels = colnames(X), col = "red", pos = 4)

##### project estimated cultivar in the PCA
Y.center <- scale(Y, scale = FALSE)  # Centering matrix for wilds
Y.cult.centered <- scale(GEA_cultivars_maf, center = attr(Y.center, "scaled:center"), scale = FALSE)

Y.cult.projected <- Y.cult.centered %*% pc_linear$rotation
Y.cult.projected<-as.data.frame(Y.cult.projected)
points(Y.cult.projected, col = "orange", pch = 19)


#######

### LFMM + PCA
Y <- GEA_124
sel_env<- data.frame(Variables_142WW[, 4:13])
write.env(sel_env, "env_variable.env")
X = read.table("env_variable.env")

mod.lfmm2 <- lfmm2(input = Y, env = X, K = 3, effect.sizes = TRUE)
#get environment effect sizes
mod.lfmm2@B

## compute cultivar prediction
X <- as.matrix(X)
storage.mode(X) <- "numeric"  # Ensures it's really numeric
Y.pred <- tcrossprod(X, mod.lfmm2@B) #or
Y.pred <- X %*% t(mod.lfmm2@B)


pc <- prcomp(Y.pred)
# Calculate the proportion of variance explained
eigenvalues <- pc$sdev^2
variance_explained <- eigenvalues / sum(eigenvalues)


pc_wild<-as.data.frame(pc$x)
plot(pc$x, col = "darkgrey", pch = 19, main = "PCA of Predicted Genotypes")


# Project original environmental variables into PCA space
# This gives loadings of ENVIRONMENTAL variables
env_arrows <- cor(X, pc$x)

# Plot PCA of individuals
plot(pc$x, col = "gray", pch = 19,
     xlab = "PC1", ylab = "PC2", main = "LFMM + PCA")

# Add environmental vectors
arrow_scale <- 3  # tweak for visibility
arrows(0, 0, env_arrows[,1]*arrow_scale, env_arrows[,2]*arrow_scale, col = "red", length = 0.1)
text(env_arrows[,1]*arrow_scale, env_arrows[,2]*arrow_scale,
     labels = colnames(X), col = "red", pos = 4)

##### project estimated cultivar in the PCA
Y.center <- scale(Y, scale = FALSE)  # Centering matrix for wilds
Y.cult.centered <- scale(GEA_cultivars_maf, center = attr(Y.center, "scaled:center"), scale = FALSE)

Y.cult.projected <- Y.cult.centered %*% pc$rotation
Y.cult.projected<-as.data.frame(Y.cult.projected)
points(Y.cult.projected, col = "orange", pch = 19)


  #project spatial pixel point in the same PC

X_pixel <- as.matrix(scaled_pixel)
storage.mode(X_pixel) <- "numeric"


Yspatial_pred <- X_pixel %*% t(slopes_t)

pc_spatial <- prcomp(Yspatial_pred)
pc_spatial<-as.data.frame(pc_spatial$x)


plot(pc_spatial$x, col = "darkgrey", main = "PCA of spatial_pixel")
points(Y.cult.projected, col = "orange", pch = 19)

env_arrows <- cor(X_pixel, pc_spatial$x)
arrow_scale <- 5
arrows(0, 0, env_arrows[,1]*arrow_scale, env_arrows[,2]*arrow_scale, col = "red", length = 0.1)
text(env_arrows[,1]*arrow_scale, env_arrows[,2]*arrow_scale,
     labels = colnames(X), col = "red", pos = 4)

points(pc$x, col = "black", pch = 19)


###### analysis just in west Med

genoWW<- geno155[rownames(geno155)%in% list142WW$V1, ]
#imputation







#### Toffahi olive cultivar from egypt


BerriPC <- Y.cult.projected[rownames(Y.cult.projected) == "Berri_Meslal-3971", ]

spatial_LC<- data.frame(lat = pixel$y, long = pixel$x)

spatial_LC$Berri_offsetpc<- (BerriPC$PC1 - pc_spatial$PC1)^2 + (BerriPC$PC2 - pc_spatial$PC2)^2 + (BerriPC$PC3 - pc_spatial$PC3)^2 

hist(spatial_LC$Berri_offsetpc)

level1 <- quantile(spatial_LC$Berri_offsetpc, probs = 0.30, na.rm = TRUE)
level2 <- quantile(spatial_LC$Berri_offsetpc, probs = 0.5, na.rm = TRUE)
level3 <- quantile(spatial_LC$Berri_offsetpc, probs = 0.75, na.rm = TRUE)

# Compute breaks for the column
sd_breaks <- c( min(spatial_LC$Berri_offsetpc, na.rm = TRUE), level1,level2, level3, max(spatial_LC$Berri_offsetpc, na.rm = TRUE))



# Create a color palette from blue to yellow

color_palette <- c("#61b8ef", "#f4d03f", "#f38b26", "#8e44ad")  # 4 quantiles
# Assign colors based on quantiles
spatial_LC$spatial_LC$Berri_offsetpc_color <- cut(spatial_LC$Berri_offsetpc, breaks = sd_breaks, labels = color_palette)

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
spatial_sf <- st_as_sf(spatial_LC, coords = c("long", "lat"), crs = 4326)  # Assumes 'longitude' and 'latitude' columns exist
# Create the map
map <- ggplot(data = countries) +
  geom_sf(fill = "#EBEBEB", color = "black") +  # Countries' borders
  geom_sf(data = spatial_sf, aes(color = spatial_LC$Berri_offsetpc_color), size = 0.05, show.legend = FALSE) +  # Points with custom colors
  scale_color_identity() +  # Use exact colors from the 'color' column
  coord_sf(xlim = c(-15, 15), ylim = c(28, 52), expand = FALSE) +  # Set geographic limits
  theme_minimal() +
  labs(title = "Adaptive Landscape Berri Meslal") +
  theme(panel.background = element_blank())

#jpeg(file = "C:/Users/rocchetti/Desktop/running RDA GO/adaptive_landscape_rgb.jpeg",width = 18, height = 14, units = "cm", res = 800)
map
















