
# Landscape Genomics of Wild Olive Reveals Opportunities to Inform Cultivar Adaptation

## Abstract

- Perennial crops' wild relatives provide locally adapted germplasm for identifying genomic signatures of environmental adaptation. Such insights, assuming shared adaptive variation between wild and cultivated compartments, can improve predictions of cultivar adaptive potential and guide cost-effective multi-location field trials for long-lived species such as olives.
- We analyzed 27 wild olive (*Olea europaea* subsp. *europaea* var. *sylvestris*) populations across a 13° latitudinal gradient in the western Mediterranean (southern France, Corsica, Spain, Morocco). Genetic structure analyses identified genuinely wild genotypes likely to harbor locally adapted variation, later used for genotype-environment association (GEA) studies. Environmental Niche Modelling (ENM) and Redundancy Analysis (RDA) were combined to construct a landscape genomic model anchored in the species' current ecological niche. This model was then used to estimate the spatial Genomic Offset of cultivated olives, a framework we term **Cultivar Genomic Offset**.
- Predictions revealed a latitude-dependent correspondence between cultivar environments and candidate adaptive genomic regions. In a Moroccan common garden, early-flowering cultivars requiring less winter chilling showed lower local Cultivar Genomic Offset values than late-flowering ones, indicating stronger adaptation to local conditions.
- These results underscore the value of wild germplasm for uncovering candidate adaptive genomic variation that can ultimately inform the resilience of cultivated populations.

## Repository contents

Scripts are numbered in the order they are meant to be run. Each stage takes the output of the previous one as input.

| File | Description |
|---|---|
| `00_filtering_vcf.txt` | `vcftools` command log for raw SNP filtering (depth, quality, missingness, biallelic sites), individual-missingness filtering, sample removal, and minor-allele-count filtering of the wild + cultivated olive VCF. |
| `01_Pop_structure_triangularR.R` | Population structure analysis: genotype import (`vcfR`), conversion to `.geno` format, `sNMF` ancestry estimation (K = 1–10) with cross-entropy model selection, ancestry barplots by K, identification of genuinely wild individuals from admixture proportions (q > 0.6 / 0.7 / 0.8 thresholds), hybrid index / heterozygosity triangle plots (`triangulaR`) to distinguish wild, admixed, and cultivated genotypes, and PCA of genetic variation. |
| `02_Environmental_niche_model.R` | Species Distribution / Environmental Niche Modelling with `biomod2`: bioclimatic and soil raster assembly, pseudo-absence generation, single-algorithm modelling (GLM, RF, MAXENT, MAXNET), model evaluation and variable importance, ensemble modelling, and projection of the current suitable niche across the western Mediterranean. Includes a sensitivity analysis of the suitability threshold (70th/80th/90th percentile) used to define the modelled niche extent. |
| `03_GEA_landscape_genomics_cultivar_offest.R` | Core landscape genomics pipeline: LFMM2 genotype-environment association on the wild dataset (latent-factor correction, FDR/Bonferroni candidate SNP detection, Manhattan plots), enriched RDA model (loci polymorphic in both wild and cultivated compartments) constrained by climate and soil variables, variance partitioning (environment vs. geography), prediction of wild and cultivar genotypes into RDA space, calculation of spatial genomic offset ("Cultivar Genomic Offset") for individual cultivars against the wild-derived model, geographic mapping of offset surfaces, and permutation tests assessing whether a cultivar's low-offset (best-adapted) region is more spatially/latitudinally specific than expected from randomized genotypes. |
| `04_Sensitivity_analysis.R` | Robustness checks re-running the wild-genome GEA/RDA/offset workflow under alternative wild-genotype ancestry thresholds (e.g., q > 0.6 vs. q > 0.8) to confirm that candidate loci and offset predictions are not sensitive to how "genuinely wild" individuals are defined. |
| `05_gene_function.R` | Functional annotation of candidate adaptive loci: matching GEA SNP positions to genes within a ±100 kb window using the reference GFF, and screening resulting gene annotations against an abiotic-stress-related keyword set to flag candidate genes potentially involved in climate adaptation. |

## Workflow overview

```
00_filtering_vcf.txt
        │  (site + individual filtering, MAC filtering)
        ▼
01_Pop_structure_triangularR.R
        │  (sNMF ancestry, wild/hybrid/cultivar classification)
        ▼
02_Environmental_niche_model.R
        │  (ENM ensemble, current niche extent + threshold sensitivity)
        ▼
03_GEA_landscape_genomics_cultivar_offest.R
        │  (LFMM2 GEA → enriched RDA → Cultivar Genomic Offset → permutation tests)
        ├──► 04_Sensitivity_analysis.R   (robustness to wild-sample ancestry threshold)
        └──► 05_gene_function.R          (candidate gene annotation)
```

## Data

Raw genomic data consist of whole-genome resequencing of 27 wild olive populations and a panel of cultivated olive accessions, jointly genotyped against the olive reference genome. Environmental layers are bioclimatic (WorldClim-derived: bio2, bio10, bio11, bio15, bio18, bio19) and soil (nitrogen, pH, clay, sand content) rasters covering the western Mediterranean (southern France, Corsica, Spain, Morocco, Portugal, Algeria).

Raw sequence data, VCF files, and raster layers are not included in this repository due to size; paths in the scripts point to their storage location on the CIRAD CLIMOLIVEMED project server and will need to be updated to reproduce the analysis on a different system.

## Requirements

Analyses were run in R and rely on the following packages:

```r
install.packages(c(
  "tidyverse", "dplyr", "tidyr", "stringr", "data.table",
  "vcfR", "adegenet", "LEA", "triangulaR", "vegan",
  "FactoMineR", "factoextra",
  "biomod2", "raster", "rasterVis", "geodata",
  "sf", "rnaturalearth", "rnaturalearthdata", "geosphere",
  "ggplot2", "ggpubr", "ggrepel", "ggridges", "gridExtra", "gtable", "grid",
  "patchwork", "scales", "RColorBrewer",
  "qqman", "multcomp", "fit.models", "writexl"
))
```

`LEA` is a Bioconductor package and should be installed with:

```r
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("LEA")
```

SNP filtering (`00_filtering_vcf.txt`) requires [`vcftools`](https://vcftools.github.io/) installed on the system.

## Citation

If you use this code, please cite the associated manuscript:

> Landscape Genomics of Wild Olive Reveals Opportunities to Inform Cultivar Adaptation.
