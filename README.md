
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
| `00_filtering_vcf.txt` | SNP and sample quality filtering of the wild + cultivated olive VCF using `vcftools`. |
| `01_Pop_structure_triangularR.R` | Population structure and hybridization analysis, used to identify genuinely wild individuals for the landscape genomics model. |
| `02_Environmental_niche_model.R` | Environmental Niche Modelling of the species' current distribution, including a sensitivity analysis of the niche extent threshold. |
| `03_GEA_landscape_genomics_cultivar_offest.R` | Core landscape genomics pipeline: genotype-environment association, RDA model construction, and calculation and mapping of Cultivar Genomic Offset, with permutation tests of its spatial patterns. |
| `04_Sensitivity_analysis.R` | Robustness checks of the GEA/RDA/offset workflow under alternative definitions of the wild reference sample. |
| `05_gene_function.R` | Functional annotation of candidate adaptive loci and screening for genes linked to abiotic-stress response. |

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
