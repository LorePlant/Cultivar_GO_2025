# Landscape Genomic Tool to Assess the Adaptive Value of Cultivars: A Case Study in Western Mediterranean Olive
This page is created to track progresses on my postdoctoral research in landscape genomics in a wester Mediterrenean Olive population.
## Abstract
Crop wild relatives and landraces represent locally adapted germplasm, enabling landscape genomic studies to identify genomic signatures of environmental adaptation (GEA). The derived information from such studies can be used to enhance and evaluate cultivar’s adaptive value. This approach can represent an advantageous opportunity particularly for perennial crops where multi-location field trials to test specific cultivar or new breeding clone adaptations are costly and required a significant amount of time to get results.

We applied this approach to olive trees (Olea europaea subsp. europaea var. sylvestris), leveraging truly local wild germplasm from 30 populations sampled along a 13-degree north-south gradient, covering regions in Southern France, Corsica, Spain, and Morocco.
Using Redundancy analysis (RDA), we identified candidate genotype-environment associations (GEAs), allowing us to map the adaptive genetic variation of olive trees across the estimated species range in the western Mediterranean basin. These findings can inform sapling recommendations aimed at capturing wild adaptive diversity.

Building on these GEAs, we applied the same RDA framework to predict the adaptive landscape of new genotypes, which in this case was represented by widely used cultivars. For a set of cultivars the prediction allowed to identify regions where adaptation could enhance local performance and reduce reliance on agricultural inputs. The model confirms a latitude-based correspondence between specific cultivars' growing environments and adaptive regions in the western Mediterranean. Additionally, it highlights the broad adaptive-genetic variation present in western wild genotypes, which is currently underrepresented in cultivated germplasm, emphasizing the need for conservation strategies and breeding utilization.

Given the co-occurrence of wild and cultivars in the same area and the one-way gene flow from cultivars that are vegetatively propagated to wild populations, we could pinpoint areas where the cultivar offset is minimal. In these regions, the coexistence of wild and cultivated genotypes is less likely to disrupt the GEA pattern thereby preserving the integrity of wild adaptive diversity.

Collectively, these findings offer valuable insights for developing in-situ conservation strategies for wild relatives and informing breeding programs to enhance environmental adaptation.

## vcf file preparation
I started by filtering sites quality

```
vcftools --gzvcf /storage/replicated/cirad/projects/CLIMOLIVEMED/results/GenomicOffsets/Lorenzo/Leccino_new_genome24/All_wild_cultivated_olive_2_run.vcf.gz --remove-indels --minDP 10 --max-meanDP 100 --minQ 200  --max-alleles 2 --min-alleles 2 --max-missing 0.90  --recode --recode-INFO-all --out /storage/replicated/cirad/projects/CLIMOLIVEMED/results/GenomicOffsets/Lorenzo/Leccino_new_genome24/site_filtered_wild_cult_lec24_DP10_100_miss090.vcf

remouve individual withi missingness >0.85

vcftools --vcf /storage/replicated/cirad/projects/CLIMOLIVEMED/results/GenomicOffsets/Lorenzo/Leccino_new_genome24/site_filtered_wild_cult_lec24_DP10_100_miss090.vcf.recode.vcf --missing-indv --out /storage/replicated/cirad/projects/CLIMOLIVEMED/results/GenomicOffsets/Lorenzo/Leccino_new_genome24/F_missing_individuals_DP10_100

vcftools --vcf /storage/replicated/cirad/projects/CLIMOLIVEMED/results/GenomicOffsets/Lorenzo/Leccino_new_genome24/site_filtered_wild_cult_lec24_DP10_100_miss090.vcf.recode.vcf --keep /storage/replicated/cirad/projects/CLIMOLIVEMED/results/GenomicOffsets/Lorenzo/Leccino_new_genome24/list710WD_cul.txt  --recode --recode-INFO-all --out /storage/replicated/cirad/projects/CLIMOLIVEMED/results/GenomicOffsets/Lorenzo/Leccino_new_genome24/WC710_lec24_DP10_100_miss090_ind085.vcf

further remouve OES_M29_07_S9_L004  OES_M29_02_S56_L004

vcftools --vcf /storage/replicated/cirad/projects/CLIMOLIVEMED/results/GenomicOffsets/Lorenzo/Leccino_new_genome24/WC710_lec24_DP10_100_miss090_ind085.vcf --keep /storage/replicated/cirad/projects/CLIMOLIVEMED/results/GenomicOffsets/Lorenzo/Leccino_new_genome24/list710WD_cul.txt  --recode --recode-INFO-all --out /storage/replicated/cirad/projects/CLIMOLIVEMED/results/GenomicOffsets/Lorenzo/Leccino_new_genome24/WC708_lec24_DP10_100_miss090_ind085.vcf

Found monomorphic SNPs. Filter for mac minor allele count 1

vcftools --vcf /storage/replicated/cirad/projects/CLIMOLIVEMED/results/GenomicOffsets/Lorenzo/Leccino_new_genome24/WC708_lec24_DP10_100_miss090_ind085.vcf --mac 1  --recode --recode-INFO-all --out /storage/replicated/cirad/projects/CLIMOLIVEMED/results/GenomicOffsets/Lorenzo/Leccino_new_genome24/WC708_lec24_DP10_100_miss090_ind085_mac1.vcf

```

In total we obtained a dataset of 708 individuals

## Analysis of Population Structure

Let's remouve SNPs that are in linkage using the thinning option in vcf tool. This command will select a SNPs every 1Kbp. Considering the low LD present in Olive (~250bp) we can consider this window appropriate.

```
thinning vcf
vcftools --vcf /storage/replicated/cirad/projects/CLIMOLIVEMED/results/GenomicOffsets/Lorenzo/Leccino_new_genome24/WC710_lec24_DP10_100_miss090_ind085.vcf.recode.vcf --thin 1000 --recode --out /storage/replicated/cirad/projects/CLIMOLIVEMED/results/GenomicOffsets/Lorenzo/Leccino_new_genome24/WC710_lec24_DP10_100_miss090_ind085_Thinned

```
Analysis of Population Structure using LEA package

```
geno708 <- read.vcfR("WC708_lec24_DP10_100_miss090_ind085_mac1_Thinned.recode.vcf")#import vcf file
GI <- vcfR2genind(geno708)#transfrom file in genind object
geno708<-as.data.frame(GI)
geno708<-geno708%>% select(ends_with(".0"))
list708<-data.frame(row.names(geno708))
write.table(list708, "list708.txt")#save individual order

write.geno(geno708, "Pop_stru_708.geno")

pop_stru = snmf("Pop_stru_708.geno", K = 1:10, entropy = TRUE, repetitions = 10, project = "new")
```
![image](https://github.com/user-attachments/assets/2c7f1026-c05f-4253-8d50-b910331a0c9b)


Print Q matrixes for K runs from K2 to K4 

```

best = which.min(cross.entropy(pop_stru, K = 2))
qmatrix_K2 = Q(pop_stru, K = 2, run = best)

best = which.min(cross.entropy(pop_stru, K = 3))
qmatrix_K3 = Q(pop_stru, K = 3, run = best)


best = which.min(cross.entropy(pop_stru, K = 4))
qmatrix_K4 = Q(pop_stru, K = 4, run = best)


pop_info_708<-read.table("708_pop_info.txt", header = T)

```
Plot bar plot for the K runs using ggplot

```
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

```
![image](https://github.com/user-attachments/assets/22f6cad6-556c-41ea-a115-40603761bd78)


At K4 we can clearly distinguish a specific group present in the Wil West and absent in cultivars and Wild East. Using the ancestry coefficient at K=4 we selected individual q>0.7 for the Wild Weast group. 
In total 142 individuals have been selected

```
K4_Qmatrix<-cbind(pop_info_708, qmatrix_K4)
pure_wild_west <- subset(K4_Qmatrix, V1  > 0.7)
write.table(pure_wild_west, "pure_wildW_070.txt")
```

We are going to repat the population structure analysis for the 142 individuals with the aim to define population differentiation among pure wild western olive trees.

```
# pure wild based on K=4


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
pop_stru_142WW = snmf("Pop_stru_142_WW.geno", K = 1:10, entropy = TRUE, repetitions = 10, project = "new")

# plot cross-entropy criterion for all runs in the snmf project
jpeg(file = "/storage/replicated/cirad/projects/CLIMOLIVEMED/results/GenomicOffsets/Lorenzo/Leccino_new_genome24/genotypes.snmf/cross_entropy_decay.JPEG")
plot(pop_stru_142WW, col = "blue", pch = 19, cex = 1.2)
dev.off()

```
![image](https://github.com/user-attachments/assets/44d90bd8-540d-4a1d-959b-97dc463c5223)




The cross-entropy coefficient reached a minumum at K=3. We are going to use this partition to construnct the ancestry barplot using ggplot and map the spatial interpolation of ancestry coefficients in the species niche using QGIS.

```
best = which.min(cross.entropy(pop_stru_142WW, K = 2))
K2Q = Q(pop_stru_142WW, K = 2, run = best)
write.table(K2Q, "Qmatrix_K2_142WW.txt")

best = which.min(cross.entropy(pop_stru_142WW, K = 3))
K3Q = Q(pop_stru_142WW, K = 3, run = best)
write.table(K3Q, "Qmatrix_K3_142WW.txt")

pop_info_142WW<- pop_info_708[pop_info_708$id %in% pure_wildW$id, ]

K2Q<-cbind(pop_info_142WW, K2Q)

K2Q <- K2Q%>%
  pivot_longer(cols = starts_with("V"), 
               names_to = "Cluster", 
               values_to = "Ancestry")
K2W<-ggplot(K2Q, aes(x =id, y = Ancestry, fill = Cluster)) +
  geom_bar(stat = "identity", width = 1) +
  theme_minimal() +
  labs(x = "Individuals", y = "Ancestry Proportion") +
  scale_fill_manual(values = c("#E69F00", "#56B4E9")) + # Customize colors
  theme(
    axis.text.x = element_blank(),  # Hide individual labels if too many
    axis.ticks.x = element_blank(),
    panel.grid = element_blank(),
    panel.spacing = unit(0.1, "lines")) +
  facet_grid(~POP, scales = "free_x", space = "free") # Separate by population

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
  scale_fill_manual(values = c("#E69F00", "#56B4E9", "darkgrey")) + # Customize colors
  theme(
    axis.text.x = element_blank(),  # Hide individual labels if too many
    axis.ticks.x = element_blank(),
    panel.grid = element_blank(),
    panel.spacing = unit(0.1, "lines")) +
  facet_grid(~POP, scales = "free_x", space = "free") # Separate by population

K3W

```
![image](https://github.com/user-attachments/assets/16fb6a37-bc28-4731-9270-a2f64100a53c)

![image](https://github.com/user-attachments/assets/255ddb7e-5b89-4247-ad17-0003e53ecb8a)


## hybrid index

142 truly wild genotypes were selected with ancestry q>0.70. To distinguish between historical vs recent introgression, we analyzed the hybrid index using ancestry-informative SNPs. We run the intercalls heterozygosity analysis classyfing as parental group the two wild gene pool Wild East and Wild West.
The significant admixture present in our collection between wild and cultivated material can be derived from recent crossing forming F1 hybrids or from past generations of crossing where natural selection had the possibilty to act. The GEA identification presume that the associated QTL derived from processes of local adaptation where environmental selection had the generational time to act. In our case the cultivated genome from cultivars vegetatevly progated mainly evolved in the eastern part of the mediterrenan basin, so it is paramount to identify the presence of recent F1 hybrids where the cultivated genome have been recently introduced and where selections did not have the generational time to act.

To investigate the presence of F1 hybrids I identified a recent devoped Rpackage that allow to identified ancestry-informative markers and estimate their hybrids index with the relative presence of F1, BC1, BC2 or past introgression. https://omys-omics.github.io/triangulaR/index.html

In this analysis I used the vcf file that was not filtered for MAF 5%.

```
library(triangulaR)
# make a pop map
popmap<-read.table("pop_map_wild_adm.txt", header = T)

  genoLAND.VCF <- read.vcfR("/storage/replicated/cirad/projects/CLIMOLIVEMED/results/GenomicOffsets/Lorenzo/Leccino_new_genome24/WC710_lec24_DP10_100_miss090_ind085_Thinned.recode.vcf")#import vcf file

library(triangulaR)

setwd("C:/Users/rocchetti/Desktop/Leccino24/hybridization_Wild_cult")
popmap<-read.table("pop_map_wild_adm.txt", header = T)
geno708 <- read.vcfR("WC708_lec24_DP10_100_miss090_ind085_mac1_Thinned.recode.vcf")#import vcf file
vcfR.diff <- alleleFreqDiff(vcfR = geno708, pm = popmap, p1 = "WW", p2 = "WE", difference = 0.8)
# 1405 sites passed allele freauency differences of 0.8

hi.het <- hybridIndex(vcfR = vcfR.diff, pm = popmap, p1 = "WW", p2 = "WE")
cols <- c("darkgrey", "purple", "darkorange", "darkgreen")
triangle.plot(hi.het, colors = cols)
```
![image](https://github.com/user-attachments/assets/bcb97c23-7c9f-463f-8bae-67de89a23833)

The results highlight the large presence of rencet hybrids like F1 and BC1, making the admixed population not suited for landscape genomics and GEA discovery. The potential adaptation of these individuals can be given from phenotypic plasticity and/or hybrid vigor.
The large part of cultivated material is closer to the WildEast group. Among them only a few hybrids, mainly from Italt Spain and France (Picholine) are considered F1 with the WildWest group. These result confirm the unexplored diversity of WildWest for the Cultivated germplasm.
Similar results are confirmed as well from PCA analysis

```
#PCA
library(FactoMineR)
library(factoextra)

res.pca708<-PCA(geno708, scale.unit = TRUE, ncp = 5, graph = TRUE)
ind708 <- get_pca_ind(res.pca708)
pca_data708 <- as.data.frame(ind708$coord)
pca_data708<-cbind(popmap, pca_data708)
qq<-ggplot() +
  geom_hline(yintercept=0, linetype="dashed", color = gray(.80), linewidth=0.6) +
  geom_vline(xintercept=0, linetype="dashed", color = gray(.80), linewidth=0.6) +
  geom_point(data = pca_data708, aes(x=Dim.1, y=Dim.2, color = pop), size = 2.5) +
  scale_color_manual(values = c("darkgrey", "purple", "darkorange", "darkgreen")) +
  xlab("PC1: 10%") + ylab("PC2: 5.9%") +
  guides(color=guide_legend(title="Group")) +
  theme_bw(base_size = 11, base_family = "Times") +
  theme(panel.background = element_blank(), legend.background = element_blank(), panel.grid = element_blank(), plot.background = element_blank(), legend.text=element_text(size=rel(.8)), strip.text = element_text(size=11))
qq
```
![image](https://github.com/user-attachments/assets/581a716a-2c28-46c4-a84c-35a2cab53d3f)

## RDA

To construct a landscape model with genetic and envirnmental varaible we used the Redundancy Analysis (RDA)approach. As first step we are going to enter the genotypic vcf.file (filtered MAF<0.05) and Environmental variable composed by bioclimatic variables extracted from Chelsa (https://chelsa-climate.org/downloads/) and soil variable extracted from SoilGrid 250 (https://soilgrids.org/). To ensure comparability, the environmental variable distributions will be scaled to achieve a mean of 0 and a standard deviation of 1.

```
#enter vcf file
geno155 <- read.vcfR("D:/vcf_file_GEA_leccino/WC156_lec24_DP10_100_miss090_ind085_mac1_MAF005.vcf.recode.vcf")#import vcf file
GI <- vcfR2genind(geno155)#transfrom file in genind object
geno155<-as.data.frame(GI)
geno155<-geno155%>% select(ends_with(".0"))
#imputation
for (i in 1:ncol(geno155))
{
  geno155[which(is.na(geno155[,i])),i] <- median(geno155[-which(is.na(geno155[,i])),i], na.rm=TRUE)
}
geno155_data<- write.table(geno155, "geno_155.txt")
# Wild Environment datafile

#standardize bioclim variable
data_wild<- read.csv("Env_155_WWE.csv", header = TRUE)
test_env <- data_wild%>% select(long, lat, bio2, bio10, bio11, bio15, bio18, bio19, clay, N, pH, sand)
Env <- scale(test_env, center=TRUE, scale=TRUE)
# Extract the centering values
env_center <- attr(Env, "scaled:center") #mean of each variable
# Extract the scaling values
env_scale <- attr(Env, "scaled:scale") #standard deviation of each variable
#transform into dataset
Env <- as.data.frame(Env)


#combining geographic, Popstructure, environmental (scaled) variables
Variables <- data.frame(data_wild$id, data_wild$group,data_wild$region, Env)
names(Variables)[1]<- paste("geno")
names(Variables)[2]<- paste("group")
names(Variables)[3]<- paste("region")
```
Check Variance Inflation Factor of selected environmental variable including bioclim and soil variable;

```
RDAgeo_env <- rda(geno155 ~ bio2+bio10+bio11+	bio15	+ bio18 + bio19 + clay+ N+ pH+ sand , Variables)
sqrt(vif.cca(RDAgeo_env))
```

|  bio2    |   bio10   |   bio11  |  bio15  |   bio18  |   bio19  |clay    |   N  |    pH  |     sand |
|---------|----------|---------|----------|-----------|----------|-------|-------|---------|--------|
1.648260 |2.132844| 2.717164| 3.642094| 2.693061 |2.360481|3.067049 |2.399268 |1.567922| 2.776405|
