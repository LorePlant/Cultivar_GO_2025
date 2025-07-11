# Landscape Genomic Tool to Assess the Adaptive Value of Cultivars: A Case Study in Western Mediterranean Olive
This page is created to track progresses on my postdoctoral research in landscape genomics in a wester Mediterrenean Olive population.
## Abstract
Crop wild relatives and landraces represent locally adapted germplasm, enabling landscape genomic studies to identify genomic signatures of environmental adaptation. The derived information can be used to evaluate the adaptive value of new genotypes such as breeding material or commercial varieties. This approach can represent an advantageous opportunity particularly for perennial crops where multi-location field trials to test specific cultivar or new breeding clone adaptations are costly and time-consuming.

We applied this approach to olive trees (Olea europaea subsp. europaea var. sylvestris), leveraging truly local wild germplasm from 27 populations sampled along a 13-degree north-south gradient, covering regions in Southern France, Corsica, Spain, and Morocco.
Using capture sequencing targeting over 50,000 predicted genes, we performed genotype-environment association (GEA) analyses to identify candidate polymorphisms associated with climate and soil variables. This enabled us to map adaptive genetic variation across the western Mediterranean distribution of wild olive populations and inform sampling strategies to capture adaptive diversity.

Building on these GEAs, we used redundancy analysis (RDA) to estimate spatial genomic offset for widely used domesticated cultivars (Olea europaea subsp. europaea var. europaea). This predictive model identified regions where adaptation could enhance cultivar performance and reduce dependence on external agricultural inputs. The model confirms a latitude-based correspondence between specific cultivars' growing environments and adaptive regions in the western Mediterranean. 
Moreover, our predictions aligned with phenological data from the Worldwide Olive Germplasm Bank of Marrakech (WOGBM), where early-flowering cultivars—typically requiring less winter chilling to initiate flowering—exhibited lower offset values than late-flowering ones, suggesting better local adaptation under current climatic conditions.

Collectively these findings demonstrate the utility of wild germplasm for detecting adaptive genomic variation and modeling the adaptive potential of cultivated or admixed populations.
Building on the adaptive model, we highlight the broad adaptive-genetic variation present in western wild genotypes, which remain largely underrepresented in cultivated germplasm. This underscores the importance of developing in-situ conservation strategies for wild relatives and integrating their diversity into breeding programs to enhance environmental adaptation.


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

```
![image](https://github.com/user-attachments/assets/2a983895-2ca7-4ab1-93cf-c069bec1e067)


## Landscape Dataset preparation

As first step we are going to imputate the genotypic vcf.file and apply filter of maf 0.05. To ensure comparability, the environmental variable distributions will be scaled to achieve a mean of 0 and a standard deviation of 1.

```
#I enter the vcf file I have that ionclude wild east and wild west
geno155 <- read.vcfR("D:/vcf_file_GEA_leccino/WC156_lec24_DP10_100_miss090_ind085_mac1_MAF005.vcf.recode.vcf")#import vcf file
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

# filter the 142 wild West and apply maf 0.05
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


# Wild Environment datafile

data_wild<- read.csv("Env_155_WWE.csv", header = TRUE)

list142WW<-read.table("list142WW.txt") #list with 142 wild west
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

```
Check Variance Inflation Factor (VIF) of selected environmental variable from RDA model.
VIF measure how much the variance of a regression coefficient is inflated due to collinearity with all other predictors. We keep variable with VIF<5

```
RDAgeo_env <- rda(genoWW_maf ~ bio2+bio10+bio11+	bio15	+ bio18 + bio19 + clay+ N+ pH+ sand , Variables)
sqrt(vif.cca(RDAgeo_env))
```

|  bio2    |   bio10   |   bio11  |  bio15  |   bio18  |   bio19  |clay    |   N  |    pH  |     sand |
|---------|----------|---------|----------|-----------|----------|-------|-------|---------|--------|
3.733407 |1.968997| 3.874802| 4.199532| 2.594169 |2.328000|3.054609 |2.587595|1.480445| 2.885255|


## Genotype Environment Association (GEA analysis)

The GEA analysis enabled the identification of specific SNP markers associated with multivariate environmental variation, aiming to detect genomic signatures of local adaptation. To achieve this, we used the LFMM approach (https://doi.org/10.1093/molbev/msz008). This model estimates the effects of latent factors related to demographic history and incorporates them into a linear model to identify associations with individual environmental variables. Significant values for each environmental variable were then combined using a squared-max transformation, resulting in a single P-value for each marker

run GEA analysis

```
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
```
We defined a trheshold of log10Pvalue >5

```
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
```
![image](https://github.com/user-attachments/assets/c4d42434-ab6d-40da-999c-2213da576e2a)

In total we identified 255 QTLs

# RDA Genomic offset model

We used the RDA framework to construct a genomic offset model to predict the adaptation of new spatial locations and novel genotypes (cultivars). To achive this initially we selected among the 255 QTLs those one that are polymorphic across the cultivar population (maf>0.05), obtaining a total of 124 SNPs. Subsequently we used these 124 SNPs to run an enriched RDA.

Define the 124 GEA QTL
```
list_255<-read.table("listGEA_255.txt",  header=TRUE)
GEA_lfmm_255<-  genoWW_maf[, colnames(genoWW_maf)%in% list_255$SNP]
write.table(GEA_lfmm_255, "GEA_lfmm_all_var_255.txt")
GEA_lfmm_all_var<-read.table("GEA_lfmm_all_var_255.txt")

GEA_124<- GEA_lfmm_all_var[, colnames(GEA_lfmm_all_var)%in% colnames(GEA_cultivars_maf)]
write.table(GEA_124,"GEA_124_WW.txt")
GEA_124<-read.table("GEA_124_WW.txt")
```
Selection of polymorphic GEA across cultivars
```
geno_cultivar<- read.vcfR("D:/vcf_file_GEA_leccino/Cultivar_319_lec24_DP10_100_miss090_ind085_mac1.vcf.recode.vcf")
GI <- vcfR2genind(geno_cultivar)#transform file in genind object
geno_cultivar<-as.data.frame(GI)
geno_cultivar <- dplyr::select(geno_cultivar, ends_with(".0"))

GEA_lfmm_all_var<-read.table("GEA_lfmm_all_var.txt")
GEA <-colnames(GEA_lfmm_all_var)
GEA_geno_cultivar<-dplyr::select(geno_cultivar, all_of(GEA))

#imputation
for (i in 1:ncol(GEA_geno_cultivar))
{
  GEA_geno_cultivar[which(is.na(GEA_geno_cultivar[,i])),i] <- median(GEA_geno_cultivar[-which(is.na(GEA_geno_cultivar[,i])),i], na.rm=TRUE)
}
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
```
Run RDA using the selected GEA QTL that are polymorphic across the cultivar populations
```
GEA_124<- GEA_lfmm_all_var[, colnames(GEA_lfmm_all_var)%in% colnames(GEA_cultivars_maf)]
write.table(GEA_124,"GEA_124_WW.txt")
GEA_124<-read.table("GEA_124_WW.txt")

RDA_all_enriched<-rda(GEA_124 ~ bio2 + bio10 + bio11 + bio15	+ bio18 + bio19 + clay + N+ pH + sand , Variables_142WW)
summary(eigenvals(RDA_all_enriched, model = "constrained"))
RsquareAdj(RDA_all_enriched)
plot(RDA_all_enriched)
sqrt(vif.cca(RDA_all_enriched))
```
Plot the RDA biplot. In this graphical representation I used _scaling = 2_ to graphically represent the location points.

```
# plot Geographic regions

TAB_gen <- data.frame(geno = row.names(scores(RDA_all_enriched , display = "sites")), scores(RDA_all_enriched, display = "sites", scaling = "sites"))
Geno <- merge(TAB_gen, Variables_142WW[, 1:5] ,by="geno")
TAB_var <- as.data.frame(scores(RDA_all_enriched, choices=c(1,2), display="bp"))
loading_geno_all_enriched_region<-ggplot() +
  geom_hline(yintercept=0, linetype="dashed", color = gray(.80), size=0.6) +
  geom_vline(xintercept=0, linetype="dashed", color = gray(.80), size=0.6) +
  geom_point(data = Geno, aes(x=RDA1, y=RDA2, fill = region), size = 2.5, shape = 21, color = "black", stroke = 0.8) +
  scale_fill_manual(values = c("lightblue","darkgreen", "darkorange")) +
  geom_segment(data = TAB_var, aes(xend=RDA1*2, yend=RDA2*2, x=0, y=0), colour="black", linewidth =0.15, linetype=1, arrow=arrow(length = unit(0.02, "npc"))) +
  geom_label_repel(data = TAB_var, aes(x=RDA1*2, y=RDA2*2, label = row.names(TAB_var)), size = 3.2, family = "Times") +
  xlab("RDA 1: 68%") + ylab("RDA 2: 11%") +
  guides(color=guide_legend(title="Latitude gradient")) +
  theme_bw(base_size = 11, base_family = "Times") +
  theme(panel.background = element_blank(), legend.background = element_blank(), panel.grid = element_blank(), plot.background = element_blank(), legend.text=element_text(size=rel(.8)), strip.text = element_text(size=11))+
  labs(title = "enriched RDA")
loading_geno_all_enriched_region
```
![image](https://github.com/user-attachments/assets/d13ec740-e818-4c87-a2cd-c2d0f8108033)

> The biplot shows that environmental differentiation among the 142 truly wild locations follows a latitudinal gradient, with French sites associated with higher summer precipitation, lower temperatures, and greater soil fertility, while Moroccan sites are characterized by higher temperatures, lower precipitation, and reduced fertility.

## Adaptive landscape projection
I leveraged the enriched RDA model constructed using only the 142 wild western samples to predict the adaptive value of each spatial pixel within the olive niche.

As fist step lets upload the raster files previously clipped using niche model raster in QGIS

```
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
soilsand<- raster(paste("D:/raster files/Current_ENM_clipped_biova/resampled_soilsand_1.tif"))

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

soilclay <- resample(soilclay, bio2, method="bilinear")
writeRaster(soilclay, "D:/raster files/Current_ENM_clipped_biova/resampled_soilclay.tif", format="GTiff", overwrite=TRUE)

soilsand <- resample(soilsand, bio2, method="bilinear")
writeRaster(soilclay, "D:/raster files/Current_ENM_clipped_biova/resampled_soilsand.tif", format="GTiff", overwrite=TRUE)

#stack the different raster file
ras_current_var<-stack(c(bio2,bio10, bio11, bio15, bio18, bio19, soilclay,soilN,soilpH, soilsand))
plot(ras_current_var, 
     xlim = c(-10, 12), 
     ylim = c(27, 50))
```
>in this specific code chunk the function _resample_ was only used once to generate a new alligned raster. In this case the function allowed to generate new soil raster files that were named _resample_soil_. Once generetaed this new files were opend with the function _raster_ and than stacked all together with the function _stack_


Transform the stacked raster in table and use the previous environmental scaling factor to scale the pixel table.
NB: clay, pH; N and sand were adjusted according to the previous RDA model. In particular clay, pH, sand are divided by 10 a Nitrigen is dividied by 100.
Subsequently the environlental matrix is scaled using the mean (center) and standard deviation (scale) of the environmental variables of the 142 reference population

```
pixel <- as.data.frame(rasterToPoints(ras_current_var))
pixel <- data.frame(x=pixel$x, y=pixel$y, bio2=pixel$bio2,bio10=pixel$bio10,bio11=pixel$bio11,bio15=pixel$bio15,bio18=pixel$bio18,bio19=pixel$bio19,clay=pixel$clay/10,N=pixel$N/100,pH=pixel$pH/10,sand=pixel$sand/10)
pixel<-na.omit(pixel)
pixel<- pixel[pixel$x>-10, ]
pixel_env<- pixel%>% dplyr::select(bio2, bio10, bio11, bio15, bio18, bio19,clay, N, pH, sand)

scaled_pixel <- scale(pixel_env, center = env_center, scale = env_scale)
scaled_pixel<-as.data.frame(scaled_pixel)
```
Used the RDA model (_RDA_142WW_enriched_) to predict pixel adaptive value (position in the RDA space). 
I used the _predict_ function of _vegan_ package _type="lc"_. 

This function allows to compute the site (pixel) scores  as a **linear combination of environmental variables**:

$$
LC_i = \sum_{j} (X_{ij} \cdot b_j)
$$

Where:
- `LC_i` is the linear constrained score for site `i`,
- `X_{ij}` is the value of environmental variable `j` for site `i`,
- `b_j` is the regression coefficient for environmental variable `j`.

```
#prediction of pixel in the RDA space
scaled_pixel_LC <- predict(RDA_all_enriched, newdata=scaled_pixel, type="lc", scaling  = "sites")
TAB_pixel_LC<- data.frame(lat = pixel$y, long = pixel$x, scaled_pixel_LC)
TAB_var <- as.data.frame(scores(RDA_all_enriched, choices=c(1,2), display="bp"))
```
>NB _scaling = "sites_ is only used for graphical representation. For the prediction of Genomic Offset run the _predict_ function without the _scaling_ option

Graphical representation using RDA biplot and geographic map

```
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
```

![image](https://github.com/user-attachments/assets/587c2a14-cf53-4d79-b0ef-e68801acaac7)

## Cultivar Genomic offset

Prepare genotypic datafile for cultivated accessions
Filter cultivar accessions frol vcf file
```
vcftools --vcf /storage/replicated/cirad/projects/CLIMOLIVEMED/results/GenomicOffsets/Lorenzo/Leccino_new_genome24/WC708_lec24_DP10_100_miss090_ind085_mac1.vcf.recode.vcf --keep /storage/replicated/cirad/projects/CLIMOLIVEMED/results/GenomicOffsets/Lorenzo/Leccino_new_genome24/list_cultivars.txt  --recode --recode-INFO-all --out /storage/replicated/cirad/projects/CLIMOLIVEMED/results/GenomicOffsets/Lorenzo/Leccino_new_genome24/Cultivar_319_lec24_DP10_100_miss090_ind085_mac1.vcf
```
Upload vcf file in R and filter for GEA QTLs. This procedure was run in the MesoLab cluster

```
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
```
Used the enriched RDA model to predict the adaptive value of new genotypes using the function _predict_. The function with _type = "wa"_ estimate the weighted mean value of the genotypes based on its SNP value. 

The site scores are computed as a **weighted average of species (allele) scores**:

$$
WA_i = \sum_{s} (Y_{is} \cdot u_s)
$$

Where:
- `WA_i` is the weighted average site score for site `i`,
- `Y_{is}` is the abundance (or presence-absence) of species (allele) `s` at site `i`,
- `u_s` is the species (allele) score for species (allele) `s`.

```
RDAscore_cul <- predict(RDA_all_enriched, newdata=GEA_cultivars_maf, type="wa")
RDAscore_cul<-as.data.frame(RDAscore_cul)
```
Plot prediction of Wild West and Cultivar in the same RDA biplot

```
################################# map cultivar and wild in the RDA space
Tab_cultivar<- data.frame(geno = row.names(RDAscore_cul),RDAscore_cul[,1:2] )
Tab_cultivar$group<-"cultivars"

RDA_wild<-predict(RDA_all_enriched,newdata =  GEA_124, type = "wa")
Tab_wild<-data.frame(geno = row.names(RDA_wild),RDA_wild[,1:2])
Tab_wild$group<-"wild"

TAB_var <- as.data.frame(scores(RDA_all_enriched, choices=c(1,2), display="bp"))

wild_cult_pred<-rbind(Tab_wild, Tab_cultivar)
wild_cult_pred$group <- as.factor(wild_cult_pred$group)

hh <- ggplot() +
  geom_hline(yintercept=0, linetype="dashed", color = gray(.80), size=0.6) +
  geom_vline(xintercept=0, linetype="dashed", color = gray(.80), size=0.6) +
  geom_point(data = wild_cult_pred, aes(x = RDA1, y = RDA2, fill = group, shape = group),size = 2.5, color = "black", stroke = 0.8)+
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
```
![image](https://github.com/user-attachments/assets/aa59112c-3619-4f45-a9c1-67dd569ed539)

As we plotted wildwest and cultivar in the RDA space we can plot them in the PCA space

```
#PCA
library(FactoMineR)
library(factoextra)

GEA_wild_cultivar<-rbind(GEA_124, GEA_cultivars_maf)

res.pcacultivar<-PCA(GEA_wild_cultivar, scale.unit = TRUE, ncp = 5, graph = TRUE)
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
```
![image](https://github.com/user-attachments/assets/ff9beb55-5451-460f-9973-7a77fb4b7f49)
> The plot shows that among the cultivated group, the genetic diversity explained by GEA QTLs differentiate early genotype from late genotype. Genotypes flowering classes derived from the work of Omar https://www.mdpi.com/2073-4395/12/12/2975

We can use the predicted position of each singol cultivar to predict the current spatial offset in the olive niche. The current spatial genomic offset is calcolated as quadratic distance between the two predicted point in the RDA space

```

F <- GEA_cultivars[rownames(GEA_cultivars) == "Istarska_crnica1", ]

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
  map <- ggplot(data = countries) +
    geom_sf(fill = "#EBEBEB", color = "black") +
    geom_sf(data = TAB_pixel_LC_sf, aes(color = TAB_pixel_LC$Foffset), size = 0.5, show.legend = FALSE) +
    scale_color_manual(values = color_palette, name = "Offset Quantile") +
    coord_sf(xlim = c(-15, 15), ylim = c(28, 52), expand = FALSE) +
    theme_minimal() +
    labs(title = "Adaptive Landscape Aggezi_Akse1") +
    theme(
      panel.background = element_blank(),
      legend.position = "right",
      legend.title = element_text(size = 10),
      legend.text = element_text(size = 8)
    )
  
map
```
![image](https://github.com/user-attachments/assets/c8b4c02a-5ff1-474b-a52e-451fcfb2d72f)

> Here is an example of spatial genomic offset estimated for four cultivars. The two cultivars at the top, Frantoio and Razzaio, originate from the northern shore of the Mediterranean (Tuscany, Italy), while the two at the bottom, Karme and Berri Meslal, originate from southern latitudes in Egypt and Morocco, respectively. The results show that the model effectively captured the adaptive responses along the latitudinal gradient, highlighting that cultivars originating from the northern and southern Mediterranean shores are more adaptive to the northern and southern regions of our study area, respectively


## Offset validation with flowering data from WOGBM
To validate further validate the spatial GO results of cultivars we used the flowering data derived from the WOGBM published in by Abou-Saaid et al 2022 https://doi.org/10.3390/agronomy12122975. In that study, the authors confirmed that early-flowering cultivars have lower chilling requirements, making them potentially better adapted to warmer climates.

Given the combined influence of chilling and forcing temperatures on floral initiation, our aim was to investigate whether a correlation exists between flowering time and the estimated genomic offset of cultivars at the WOGBM location.

```
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
```
Adjust the table with the flowering data from the Abou-Saaid et al 2022 publication

```
GoMar<-read.csv("GO_marakesh.csv")
GoMar <- na.omit(GoMar)

boxplot(GO ~ class, data = GoMar)
library(ggplot2)

ggplot(GoMar, aes(x = class, y = GO, fill = class)) +
  geom_boxplot() +
  scale_fill_manual(values = c("purple", "darkblue", "lightgrey")) +
  theme_minimal() +
  labs(title = "GO by Class",
       x = "Class",
       y = "GO") +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text.x = element_text(angle = 0, hjust = 0.5)
  )
model <- lm(GO ~ class, data = GoMar)

summary(model)

GoMar$class <- as.factor(GoMar$class)
 model <- lm(GO ~ class, data = GoMar)

library(multcomp)
 summary(glht(model, linfct = mcp(class = "Tukey")))

hist(GoMar$FFD)
abline(v = c(120, 126), col = "red", lwd = 2, lty = 2)


plot(GO ~ FFD, data = GoMar)

model <- lm(GO ~ FFD, data = GoMar)
abline(model, col = "red", lwd = 2)
summary(model)

ggplot(GoMar, aes(x = FFD, y = GO, fill = class)) +
  geom_point(shape = 21, size = 2.5, color = "black", alpha = 0.8) +
  geom_abline(intercept = coef(model)[1], slope = coef(model)[2],
              color = "red", linetype = "solid", linewidth = 1.2) +
  scale_fill_manual(values = c("purple", "darkblue", "lightgrey")) +
  theme_minimal() +
  labs(title = "GO vs FFD by Class",
       x = "FFD",
       y = "GO",
       fill = "Class") +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold")
  )
```
![image](https://github.com/user-attachments/assets/8458fd86-255a-4357-8a56-c5d51df12d91)

Results show a significant correlation between flowering time and Genomic offset. While the overall variance explained by the model was modest (R² ≈ 5%), we observed a strong and statistically significant difference between early- and late-flowering cultivar classes.
Although the overall variance explained was modest, the marked differentiation between early- and late-flowering cultivars supports the biological relevance of the predictions.

Importantly, from the overall GEAs identified in the wild, we filtered those (GEAs) that were also polymorphic within cultivated germplasm. The correlation between these polymorphisms and flowering time suggests that the adaptive variation captured in cultivars may result from multiple evolutionary processes. These include parallel selection under similar environmental pressures, introgression from wild relatives during post-domestication period, and/or balancing selection maintaining functional diversity across both wild and cultivated gene pools.

