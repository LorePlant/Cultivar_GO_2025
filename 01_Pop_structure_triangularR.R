
setwd("D:/C/Desktop/Leccino24/PopulationStructure")


library(tidyverse)
library(ggpubr)
library(vcfR)
library(LEA)

geno708 <- read.vcfR("D:/D/vcf_file_GEA_leccino/WC708_lec24_DP10_100_miss090_ind085_mac1_Thinned.recode.vcf")#import vcf file
GI <- vcfR2genind(geno708)#transfrom file in genind object
geno708<-as.data.frame(GI)
geno708<-geno708%>% select(ends_with(".0"))
list708<-data.frame(row.names(geno708))
write.table(list708, "list708.txt")#save individual order

write.geno(geno708, "Pop_stru_708.geno")

pop_stru = snmf("Pop_stru_708.geno", K = 1:10, entropy = TRUE, repetitions = 10, project = "new")

# plot cross-entropy criterion for all runs in the snmf project
jpeg(file = "/storage/replicated/cirad/projects/CLIMOLIVEMED/results/GenomicOffsets/Lorenzo/Leccino_new_genome24/genotypes.snmf/cross_entropy_decay.JPEG")
plot(pop_stru, col = "blue", pch = 19, cex = 1.2)
dev.off()

## Print Q matrixes for K runs from K2 to K4


best = which.min(cross.entropy(pop_stru, K = 2))
qmatrix_K2 = Q(pop_stru, K = 2, run = best)

best = which.min(cross.entropy(pop_stru, K = 3))
qmatrix_K3 = Q(pop_stru, K = 3, run = best)


best = which.min(cross.entropy(pop_stru, K = 4))
qmatrix_K4 = Q(pop_stru, K = 4, run = best)


pop_info_708<-read.table("708_pop_info.txt", header = T)
## K2 partition

K2_Qmatrix<-cbind(pop_info_708, qmatrix_K2)


K2_Qmatrix <- K2_Qmatrix%>%
  pivot_longer(cols = starts_with("V"), 
               names_to = "Cluster", 
               values_to = "Ancestry")
  


K2<-ggplot(K2_Qmatrix, aes(x =id, y = Ancestry, fill = Cluster)) +
  geom_bar(stat = "identity", width = 1) +
  theme_minimal() +
  labs(x = "Individuals", y = "Ancestry Proportion") +
  scale_fill_manual(values = c("deepskyblue4", "darkorange")) + # Customize colors
  theme(
    axis.text.x = element_blank(),  # Hide individual labels if too many
    axis.ticks.x = element_blank(),
    panel.grid = element_blank(),
    panel.spacing = unit(0.1, "lines")) +
  facet_grid(~POP, scales = "free_x", space = "free") # Separate by population

K2

K3_Qmatrix<-cbind(pop_info_708, qmatrix_K3)


K3_Qmatrix <- K3_Qmatrix%>%
  pivot_longer(cols = starts_with("V"), 
               names_to = "Cluster", 
               values_to = "Ancestry")

K3<-ggplot(K3_Qmatrix, aes(x =id, y = Ancestry, fill = Cluster)) +
  geom_bar(stat = "identity", width = 1) +
  theme_minimal() +
  labs(x = "Individuals", y = "Ancestry Proportion") +
  scale_fill_manual(values = c("darkorange","deepskyblue4", "darkgray")) + # Customize colors
  theme(
    axis.text.x = element_blank(),  # Hide individual labels if too many
    axis.ticks.x = element_blank(),
    panel.grid = element_blank(),
    panel.spacing = unit(0.1, "lines")) +
  facet_grid(~POP, scales = "free_x", space = "free_x") # Separate by population

K3

K4_Qmatrix<-cbind(pop_info_708, qmatrix_K4)


K4_Qmatrix <- K4_Qmatrix%>%
  pivot_longer(cols = starts_with("V"), 
               names_to = "Cluster", 
               values_to = "Ancestry")

K4<-ggplot(K4_Qmatrix, aes(x =id, y = Ancestry, fill = Cluster)) +
  geom_bar(stat = "identity", width = 1) +
  theme_minimal() +
  labs(x = "Individuals", y = "Ancestry Proportion") +
  scale_fill_manual(values = c("darkgreen", "darkorange", "gray", "deepskyblue4")) + # Customize colors
  theme(
    axis.text.x = element_blank(),  # Hide individual labels if too many
    axis.ticks.x = element_blank(),
    panel.grid = element_blank(),
    panel.spacing = unit(0.1, "lines")) +
  facet_grid(~POP, scales = "free_x", space = "free_x") # Separate by population

K4

ggarrange(K2,K3,K4,nrow=3,ncol=1)

# pure wild based on K=4

K4_Qmatrix<-cbind(pop_info_708, qmatrix_K4)
pure_wild_west <- subset(K4_Qmatrix, V1  > 0.7)
write.table(pure_wild_west, "pure_wildW_070.txt")

############################################## WILD pop structure
# filtering q>0.7
Pop_stru_708 <- load.snmfProject("Pop_stru_708.snmfProject")

best = which.min(cross.entropy(Pop_stru_708, K = 4))
qmatrix_K4 = Q(Pop_stru_708, K = 4, run = best)
K4_Qmatrix<-cbind(pop_info_708, qmatrix_K4)
pure_wild_west <- subset(K4_Qmatrix, V1  > 0.7)
pure_wildW <- pure_wild_west %>% select(id)

geno142_WW <-  geno708[rownames(geno708)%in% pure_wildW$id, ]
list142_wildW<- data.frame(rownames(geno142_WW))
write.table(list142_wildW, "list_142WW.txt")
setwd("C:/Users/rocchetti/Desktop/Leccino24/PopulationStructure/Pop_structure_142_wild")

write.geno(geno142_WW, "Pop_stru_142_WW.geno")

# filtering q>0.6
Pop_stru_708 <- load.snmfProject("Pop_stru_708.snmfProject")

best = which.min(cross.entropy(Pop_stru_708, K = 4))
qmatrix_K4 = Q(Pop_stru_708, K = 4, run = best)
K4_Qmatrix<-cbind(pop_info_708, qmatrix_K4)
pure_wild_west <- subset(K4_Qmatrix, V1  > 0.6)
pure_wildW <- pure_wild_west %>% select(id)

geno154_WW <-  geno708[rownames(geno708)%in% pure_wildW$id, ]
list154_wildW<- data.frame(rownames(geno154_WW))
write.table(list154_wildW, "list_154WW.txt")

write.geno(geno154_WW, "Pop_stru_154_WW.geno")


# filtering q>0.8
Pop_stru_708 <- load.snmfProject("Pop_stru_708.snmfProject")

best = which.min(cross.entropy(Pop_stru_708, K = 4))
qmatrix_K4 = Q(Pop_stru_708, K = 4, run = best)
K4_Qmatrix<-cbind(pop_info_708, qmatrix_K4)
pure_wild_west <- subset(K4_Qmatrix, V1  > 0.8)
pure_wildW <- pure_wild_west %>% select(id)

geno126_WW <-  geno708[rownames(geno708)%in% pure_wildW$id, ]
list126_wildW<- data.frame(rownames(geno126_WW))
write.table(list126_wildW, "list_126WW.txt")

write.geno(geno126_WW, "Pop_stru_126_WW.geno")




pop_stru_142WW = snmf("Pop_stru_142_WW.geno", K = 1:10, entropy = TRUE, repetitions = 10, project = "new")



pop_stru_142WW <- load.snmfProject("Pop_stru_142_WW.snmfProject")


# plot cross-entropy criterion for all runs in the snmf project
jpeg(file = "/storage/replicated/cirad/projects/CLIMOLIVEMED/results/GenomicOffsets/Lorenzo/Leccino_new_genome24/genotypes.snmf/cross_entropy_decay.JPEG")
plot(pop_stru_142WW, col = "blue", pch = 19, cex = 1.2)
dev.off()

best = which.min(cross.entropy(pop_stru_142WW, K = 2))
K2Q = Q(pop_stru_142WW, K = 2, run = best)
write.table(K2Q, "Qmatrix_K2_142WW.txt")

best = which.min(cross.entropy(pop_stru_142WW, K = 3))
K3Q = Q(pop_stru_142WW, K = 3, run = best)
write.table(K3Q, "Qmatrix_K3_142WW.txt")

pop_info_142WW<- pop_info_708[pop_info_708$id %in% pure_wildW$id, ]
write.csv(pop_info_142WW, "pop_info_142WW.csv", row.names = FALSE)
pop_info_142WW<-read.csv("pop_info_142WW.csv", header = T)


K2Q<-cbind(pop_info_142WW, K2Q)

K2Q <- K2Q%>%
  pivot_longer(cols = starts_with("V"), 
               names_to = "Cluster", 
               values_to = "Ancestry")

K2Q$id <- factor(K2Q$id, levels = unique(K2Q$id))

K2Q <- K2Q %>%
  arrange(POP, sort) %>%                 # Sort first by POP, then alphabetically within POP
  mutate(id = factor(id, levels = unique(id)))  # Set factor levels based on this order

# Step 3: Generate the plot
K2W <- ggplot(K2Q, aes(x = id, y = Ancestry, fill = Cluster)) +
  geom_bar(stat = "identity", width = 1) +
  theme_minimal() +
  labs(x = "Individuals", y = "Ancestry Proportion") +
  scale_fill_manual(values = c("#E69F00", "#56B4E9")) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.grid = element_blank(),
    panel.spacing = unit(0.1, "lines")
  ) +
  facet_grid(~sort, scales = "free_x", space = "free")  # Now using POP

K2W

# Step 4: Display the plot
K2W

K3Q<-cbind(pop_info_142WW, K3Q)

K3Q <- K3Q%>%
  pivot_longer(cols = starts_with("V"), 
               names_to = "Cluster", 
               values_to = "Ancestry")
K3W<-ggplot(K3Q, aes(x =id, y = Ancestry, fill = Cluster)) +
  geom_bar(stat = "identity", width = 1) +
  theme_minimal() +
  labs(x = "Individuals", y = "Ancestry Proportion") +
  scale_fill_manual(values = c("#56B4E9", "#E69F00", "darkgrey")) + # Customize colors
  theme(
    axis.text.x = element_blank(),  # Hide individual labels if too many
    axis.ticks.x = element_blank(),
    panel.grid = element_blank(),
    panel.spacing = unit(0.1, "lines")) +
  facet_grid(~POP, scales = "free_x", space = "free") # Separate by population

K3W




library(ggpubr)
ggarrange(K2W,K3W,nrow=2, ncol=1)


project <- load.snmfProject("genotypes.snmfProject")




#### hybrid

library(triangulaR)

setwd("D:/C/Desktop/Leccino24/hybridization_Wild_cult")
popmap<-read.table("pop_map_wild_adm.txt", header = T)
geno710 <- read.vcfR("D:/D/vcf_file_GEA_leccino/WC708_lec24_DP10_100_miss090_ind085_mac1_Thinned.recode.vcf")#import vcf file
vcfR.diff <- alleleFreqDiff(vcfR = geno710, pm = popmap, p1 = "WW", p2 = "WE", difference = 0.8)
hi.het <- hybridIndex(vcfR = vcfR.diff, pm = popmap, p1 = "WW", p2 = "WE")
cols <- c("#777777", "#7B3294", "#D95F02",  "#009E00")
triangle.plot(hi.het, colors = cols)

write.table(hi.het, "hybrids_classes.txt")

triangle_data <- data.frame(
  x = c(0, 1, 0.5),   # Example coordinates for a triangle
  y = c(0, 0, 1)  # Height of the triangle
)

ggplot(hi.het, aes(x = hybrid.index, y = heterozygosity, fill = factor(pop))) +
  # Plot the triangle as the background
  geom_polygon(data = triangle_data, aes(x = x, y = y), fill = "lightgrey", color = "black", alpha = 0.3) +
  # Plot the points inside the triangle
  geom_point(aes(color = factor(pop)), size = 3) +
  # Customize the fill colors for different populations
  scale_fill_manual(values = c("#777777", "#7B3294", "#D95F02",  "#009E00")) +  
  scale_color_manual(values = c("#777777", "#7B3294", "#D95F02",  "#009E00")) +  # Same color mapping for points
  theme_minimal() +
  theme(legend.position = "right")  # Ensure the legend appears on the right side



triangle_data <- data.frame(
  x = c(0, 1, 0.5),
  y = c(0, 0, 1)
)

t<-ggplot(hi.het, aes(
  x = hybrid.index,
  y = heterozygosity,
  fill = factor(pop)
)) +
  
  # Triangle background
  geom_polygon(
    data = triangle_data,
    aes(x = x, y = y),
    fill = "lightgrey",
    color = "black",
    alpha = 0.3,
    inherit.aes = FALSE
  ) +
  
  # Points
  geom_point(
    aes(color = factor(pop)),
    size = 3
  ) +
  
  # Colours
  scale_fill_manual(
    values = c(
      "#777777",
      "#7B3294",
      "#D95F02",
      "#009E00"
    ),
    name = "Groups"
  ) +
  
  scale_color_manual(
    values = c(
      "#777777",
      "#7B3294",
      "#D95F02",
      "#009E00"
    ),
    name = "Groups"
  ) +
  
  # Axis labels
  labs(
    x = "Hybrid index",
    y = "Heterozygosity"
  ) +
  
  theme_classic() +
  
  theme(
    # Axis titles
    axis.title.x = element_text(
      size = 16,
      color = "black"
    ),
    axis.title.y = element_text(
      size = 16,
      color = "black"
    ),
    
    # Axis tick labels
    axis.text.x = element_text(
      size = 13,
      color = "black"
    ),
    axis.text.y = element_text(
      size = 13,
      color = "black"
    ),
    
    # Black axis lines
    axis.line = element_line(
      color = "black",
      linewidth = 0.8
    ),
    
    # Tick marks
    axis.ticks = element_line(
      color = "black",
      linewidth = 0.6
    ),
    
    # Legend
    legend.title = element_text(
      size = 14,
      color = "black"
    ),
    legend.text = element_text(
      size = 12,
      color = "black"
    ),
    legend.position = "right"
  )

ggsave(
  filename = "triangular.tiff",
  plot = t,      # Optional if last plot was your desired one
  width = 6.8,             # In inches (default)
  height = 5.1,             # In inches
  units = "in",           # Can be "in", "cm", or "mm"
  dpi = 300               # Resolution (important for publications)
)






















gl.genoLAND <- vcfR2genind(genoLAND.VCF)#transfrom file in genind object
geno710<-as.data.frame(gl.genoLAND)
geno710<-geno710%>% select(ends_with(".0"))


library(FactoMineR)
library(factoextra)

res.pca708<-PCA(geno708, scale.unit = FALSE, ncp = 5, graph = TRUE)
ind708 <- get_pca_ind(res.pca708)
pca_data708 <- as.data.frame(ind708$coord)





pca_data708<-cbind(popmap, pca_data708)


qq<-ggplot() +
  geom_hline(yintercept=0, linetype="dashed", color = gray(.80), linewidth=0.6) +
  geom_vline(xintercept=0, linetype="dashed", color = gray(.80), linewidth=0.6) +
  geom_point(data = pca_data708, 
             aes(x = Dim.1, y = Dim.2, fill = pop), 
             size = 3.5, shape = 21, color = "black", stroke = 0.5) +
  scale_fill_manual(values = c("darkgrey", "purple", "darkorange", "darkgreen")) +
  xlab("PC1: 10%") + ylab("PC2: 5.9%") +
  guides(color=guide_legend(title="Group")) +
  theme_bw(base_size = 11, base_family = "Times") +
  theme(panel.background = element_blank(), legend.background = element_blank(), panel.grid = element_blank(), plot.background = element_blank(), legend.text=element_text(size=rel(.8)), strip.text = element_text(size=11))
qq




library(FactoMineR)
library(factoextra)


PCA_lea<-vcf2lfmm("WC710_lec24_DP10_100_miss090_ind085_mac1_Thinned.recode.vcf.recode.vcf", "genotypes.lfmm")
# Read the genotype data
genotypes <- read.lfmm("WC710_lec24_DP10_100_miss090_ind085_mac1_Thinned.recode.vcf.recode.lfmm")


# Compute PCA
pca_result <- pca("WC710_lec24_DP10_100_miss090_ind085_mac1_Thinned.recode.vcf.recode.lfmm", scale = TRUE)

# Plot the first two principal components
plot(pca_result$projections[,1], pca_result$projections[,2],
     xlab="PC1", ylab="PC2", main="PCA of Genotypes", pch=20, col="blue")
  
