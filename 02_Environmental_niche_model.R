setwd("D:/C/Desktop/Leccino24/ENM_403")
install.packages("biomod2", dependencies = TRUE)
install.packages("rasterVis")
install.packages("geodata")

library(biomod2)
library(ggplot2)
library(gridExtra)
library(raster)
library(rasterVis)
library(geodata)
library(tidyr)

setwd("/storage/replicated/cirad/projects/CLIMOLIVEMED/results/GenomicOffsets/Lorenzo/Leccino_new_genome24/ENM_medit")
install.packages("biomod2", dependencies = TRUE)
install.packages("rasterVis")
install.packages("geodata")

library(biomod2)
library(ggplot2)
library(gridExtra)
library(raster)
library(rasterVis)
library(geodata)
library(tidyr)




#open the raster from the directory without downloading it everytime


bio2<- raster(paste("/storage/replicated/cirad/projects/CLIMOLIVEMED/results/GenomicOffsets/Lorenzo/Leccino_new_genome24/ENM_medit/rasters_med/bio2_west_med.tif"))
bio10<- raster(paste("/storage/replicated/cirad/projects/CLIMOLIVEMED/results/GenomicOffsets/Lorenzo/Leccino_new_genome24/ENM_medit/rasters_med/bio10_west_med.tif"))
bio11<- raster(paste("/storage/replicated/cirad/projects/CLIMOLIVEMED/results/GenomicOffsets/Lorenzo/Leccino_new_genome24/ENM_medit/rasters_med/bio11_west_med.tif"))
bio15<- raster(paste("/storage/replicated/cirad/projects/CLIMOLIVEMED/results/GenomicOffsets/Lorenzo/Leccino_new_genome24/ENM_medit/rasters_med/bio15_west_med.tif"))
bio18<- raster(paste("/storage/replicated/cirad/projects/CLIMOLIVEMED/results/GenomicOffsets/Lorenzo/Leccino_new_genome24/ENM_medit/rasters_med/bio18_west_med.tif"))
bio19<- raster(paste("/storage/replicated/cirad/projects/CLIMOLIVEMED/results/GenomicOffsets/Lorenzo/Leccino_new_genome24/ENM_medit/rasters_med/bio19_west_med.tif"))
soilN<- raster(paste("/storage/replicated/cirad/projects/CLIMOLIVEMED/results/GenomicOffsets/Lorenzo/Leccino_new_genome24/ENM_medit/rasters_med/resampled_soilN_west_med.tif"))
soilpH<- raster(paste("/storage/replicated/cirad/projects/CLIMOLIVEMED/results/GenomicOffsets/Lorenzo/Leccino_new_genome24/ENM_medit/rasters_med/resampled_soilpH_west_med.tif"))
soilclay<- raster(paste("/storage/replicated/cirad/projects/CLIMOLIVEMED/results/GenomicOffsets/Lorenzo/Leccino_new_genome24/ENM_medit/rasters_med/resampled_soilclay_west_med.tif"))
soilsand<- raster(paste("/storage/replicated/cirad/projects/CLIMOLIVEMED/results/GenomicOffsets/Lorenzo/Leccino_new_genome24/ENM_medit/rasters_med/resampled_soilsand_west_med.tif"))
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


#stack bio clim variables

ras_current_var<-stack(c(bio2,bio10, bio11, bio15, bio18, bio19, soilclay,soilN,soilpH, soilsand))




#spatial grid

coord_r<-rasterToPoints(bioclim_current, spatial = TRUE)
map_pts<-data.frame(x = coordinates(coord_r)[,1], y=coordinates(coord_r)[,2], coord_r@data)

#format data with pseudo absence
locations = read.csv("403_wild_geoloc.csv", header=TRUE)
Olive_data = BIOMOD_FormatingData(resp.var=locations$P.A, 
                                  expl.var=ras_current_var,
                                  resp.xy= locations[,c("longitude","latitude")],
                                  resp.name = "West_med_clim_soil",
                                  PA.nb.rep = 3,
                                  PA.nb.absence = 500,
                                  PA.strategy = "random")
plot(Olive_data)
jpeg(file="PA_test.jpg")
plot(Olive_data)
dev.off()

###########run individual model

Olive_niches<- BIOMOD_Modeling(bm.format = Olive_data,
                               modeling.id = 'Single_Olive_niche_maxent',
                               models = c('GLM', 'RF', 'MAXENT', 'MAXNET'),
                               CV.strategy = 'random',
                               CV.nb.rep = 2,
                               CV.perc = 0.80,
                               OPT.strategy = "default",
                               metric.eval = c('TSS','ROC'),
                               var.import = 2,
                               seed.val = 42,
                               weights = NULL)


Olive_niches

## reopen the model
Olive_niches <- get(load("/storage/replicated/cirad/projects/CLIMOLIVEMED/results/GenomicOffsets/Lorenzo/Leccino_new_genome24/ENM_medit/West.med.clim.soil/West.med.clim.soil.Single_Olive_niche_maxent.models.out"))


# Get evaluation scores & variables importance
get_evaluations(Olive_niches)
get_variables_importance(Olive_niches)
var_imp<-data.frame(get_variables_importance(Olive_niches))
write.table(var_imp, 'variqnce.importance')


# Represent evaluation scores & variables importance
bm_PlotEvalMean(bm.out = Olive_niches)
bm_PlotEvalBoxplot(bm.out = Olive_niches, group.by = c('algo', 'algo'))
bm_PlotVarImpBoxplot(bm.out = Olive_niches, group.by = c('expl.var', 'algo', 'algo'))
#bm_PlotResponseCurves(bm.out = Olive_niches,models.chosen = get_built_models(Olive_niches)[c(121:130)], fixed.var = 'mean')


########### Project single models for all raster point using the same env variable of PA dataset

get_built_models(Olive_niches, full.name = NULL, PA = NULL, run = NULL, algo = NULL)#selecting models

myBiomodProj <- BIOMOD_Projection(bm.mod = Olive_niches,
                                  proj.name = 'Current',
                                  new.env = bioclim_current,
                                  models.chosen = 'all',
                                  metric.binary = c('TSS', 'ROC'),
                                  metric.filter = c('TSS', 'ROC'),
                                  build.clamping.mask = TRUE)
myBiomodProj
plot(myBiomodProj) #not enough RAMt

prova<- raster(paste("C:/Users/rocchetti/OneDrive - Cirad/Desktop/Olive_ENM/P.A/proj_current/proj_Current_P.A_ensemble.tif"))
plot(prova)




#################### Model ensemble models

get_built_models(Olive_niches , full.name = NULL, PA = NULL, run = NULL, algo = c('RF', 'GLM', 'MAXNET'))

Olive_EM <- BIOMOD_EnsembleModeling(bm.mod = Olive_niches,
                                    models.chosen = 'all',
                                    em.by = 'all',
                                    em.algo = c('EMwmean'),
                                    metric.select = c('TSS', 'ROC'),
                                    metric.select.thresh = c(0.85, 0.85),
                                    metric.eval = c('TSS', 'ROC'),
                                    var.import = 3,
                                    EMci.alpha = 0.05,
                                    EMwmean.decay = 'proportional')

bm_PlotVarImpBoxplot(bm.out = Olive_EM, group.by = c('expl.var', 'algo', 'merged.by.run'))

bm_PlotResponseCurves(bm.out = Olive_EM, 
                      models.chosen = get_built_models(Olive_EM)[c(1, 6, 7)],
                      fixed.var = 'median')

## reopen the model
Olive_EM <- get(load("/storage/replicated/cirad/projects/CLIMOLIVEMED/results/GenomicOffsets/Lorenzo/Leccino_new_genome24/ENM_medit/West.med.clim.soil/West.med.clim.soil.Single_Olive_niche_maxent.ensemble.models.out"))


########################### Project ensamble models
########## Current environmental variable

myBiomodEMProj   <- BIOMOD_EnsembleForecasting(bm.em = Olive_EM,
                                               proj.name = 'CurrentEM_ent',
                                               new.env = ras_current_var,
                                               models.chosen = 'all',
                                               metric.binary = 'all',
                                               metric.filter = 'all')
myBiomodEMProj
plot(myBiomodEMProj)

prova<- raster(paste("C:/Users/rocchetti/OneDrive - Cirad/Desktop/Olive_ENM/P.A/proj_current/proj_Current_P.A_ensemble.tif"))
plot(prova)


# sensitivity analysis


ENM<- raster(paste("D:/C/Desktop/Leccino24/ENM_403/proj_CurrentEM_ent_West.med.clim.soil_ensemble.tif"))
summary(values(ENM))

vals <- values(ENM)
vals <- vals[!is.na(vals)]
q70 <- quantile(vals, 0.70)
q80 <- quantile(vals, 0.80)  # your original threshold
q90 <- quantile(vals, 0.90)

mask70 <- ENM >= q70
mask80 <- ENM >= q80
mask90 <- ENM >= q90
# quick sanity check on how much area each retains
cellStats(mask70, sum) / cellStats(!is.na(ENM), sum)

cellStats(mask80, sum) / cellStats(!is.na(ENM), sum)

cellStats(mask90, sum) / cellStats(!is.na(ENM), sum)

extent(ENM)  # check this first — confirm it matches your intended region

tiff("ENM_threshold_sensitivity.tif", width = 12, height = 4.5, units = "in", res = 600, compression = "lzw")

par(mfrow = c(1,3))
plot(mask70, main = "70th percentile", legend = FALSE, col = c("grey90","forestgreen"),
     xlim = c(-15,10), ylim = c(27,55))
plot(mask80, main = "80th percentile (original)", legend = FALSE, col = c("grey90","forestgreen"),
     xlim = c(-15,10), ylim = c(27,55))
plot(mask90, main = "90th percentile", legend = FALSE, col = c("grey90","forestgreen"),
     xlim = c(-15,10), ylim = c(27,55))
par(mfrow = c(1,1))

dev.off()  # this line actually writes the file to disk



bio2<- raster(paste("/storage/replicated/cirad/projects/CLIMOLIVEMED/results/GenomicOffsets/Lorenzo/Leccino_new_genome24/ENM_medit/rasters_med/bio2_west_med.tif"))
bio10<- raster(paste("/storage/replicated/cirad/projects/CLIMOLIVEMED/results/GenomicOffsets/Lorenzo/Leccino_new_genome24/ENM_medit/rasters_med/bio10_west_med.tif"))
bio11<- raster(paste("/storage/replicated/cirad/projects/CLIMOLIVEMED/results/GenomicOffsets/Lorenzo/Leccino_new_genome24/ENM_medit/rasters_med/bio11_west_med.tif"))
bio15<- raster(paste("/storage/replicated/cirad/projects/CLIMOLIVEMED/results/GenomicOffsets/Lorenzo/Leccino_new_genome24/ENM_medit/rasters_med/bio15_west_med.tif"))
bio18<- raster(paste("/storage/replicated/cirad/projects/CLIMOLIVEMED/results/GenomicOffsets/Lorenzo/Leccino_new_genome24/ENM_medit/rasters_med/bio18_west_med.tif"))
bio19<- raster(paste("/storage/replicated/cirad/projects/CLIMOLIVEMED/results/GenomicOffsets/Lorenzo/Leccino_new_genome24/ENM_medit/rasters_med/bio19_west_med.tif"))
soilN<- raster(paste("/storage/replicated/cirad/projects/CLIMOLIVEMED/results/GenomicOffsets/Lorenzo/Leccino_new_genome24/ENM_medit/rasters_med/resampled_soilN_west_med.tif"))
soilpH<- raster(paste("/storage/replicated/cirad/projects/CLIMOLIVEMED/results/GenomicOffsets/Lorenzo/Leccino_new_genome24/ENM_medit/rasters_med/resampled_soilpH_west_med.tif"))
soilclay<- raster(paste("/storage/replicated/cirad/projects/CLIMOLIVEMED/results/GenomicOffsets/Lorenzo/Leccino_new_genome24/ENM_medit/rasters_med/resampled_soilclay_west_med.tif"))
soilsand<- raster(paste("/storage/replicated/cirad/projects/CLIMOLIVEMED/results/GenomicOffsets/Lorenzo/Leccino_new_genome24/ENM_medit/rasters_med/resampled_soilsand_west_med.tif"))
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


#stack bio clim variables

ras_current_var<-stack(c(bio2,bio10, bio11, bio15, bio18, bio19, soilclay,soilN,soilpH, soilsand))




#spatial grid

coord_r<-rasterToPoints(bioclim_current, spatial = TRUE)
map_pts<-data.frame(x = coordinates(coord_r)[,1], y=coordinates(coord_r)[,2], coord_r@data)

#format data with pseudo absence
locations = read.csv("403_wild_geoloc.csv", header=TRUE)
Olive_data = BIOMOD_FormatingData(resp.var=locations$P.A, 
                                  expl.var=ras_current_var,
                                  resp.xy= locations[,c("longitude","latitude")],
                                  resp.name = "West_med_clim_soil",
                                  PA.nb.rep = 3,
                                  PA.nb.absence = 500,
                                  PA.strategy = "random")
plot(Olive_data)
jpeg(file="PA_test.jpg")
plot(Olive_data)
dev.off()

###########run individual model

Olive_niches<- BIOMOD_Modeling(bm.format = Olive_data,
                               modeling.id = 'Single_Olive_niche_maxent',
                               models = c('GLM', 'RF', 'MAXENT', 'MAXNET'),
                               CV.strategy = 'random',
                               CV.nb.rep = 2,
                               CV.perc = 0.80,
                               OPT.strategy = "default",
                               metric.eval = c('TSS','ROC'),
                               var.import = 2,
                               seed.val = 42,
                               weights = NULL)


Olive_niches

## reopen the model
Olive_niches <- get(load("/storage/replicated/cirad/projects/CLIMOLIVEMED/results/GenomicOffsets/Lorenzo/Leccino_new_genome24/ENM_medit/West.med.clim.soil/West.med.clim.soil.Single_Olive_niche_maxent.models.out"))


# Get evaluation scores & variables importance
get_evaluations(Olive_niches)
get_variables_importance(Olive_niches)
var_imp<-data.frame(get_variables_importance(Olive_niches))
write.table(var_imp, 'variqnce.importance')


# Represent evaluation scores & variables importance
bm_PlotEvalMean(bm.out = Olive_niches)
bm_PlotEvalBoxplot(bm.out = Olive_niches, group.by = c('algo', 'algo'))
bm_PlotVarImpBoxplot(bm.out = Olive_niches, group.by = c('expl.var', 'algo', 'algo'))
#bm_PlotResponseCurves(bm.out = Olive_niches, 
models.chosen = get_built_models(Olive_niches)[c(121:130)],
fixed.var = 'mean')


########### Project single models for all raster point using the same env variable of PA dataset

get_built_models(Olive_niches, full.name = NULL, PA = NULL, run = NULL, algo = NULL)#selecting models

myBiomodProj <- BIOMOD_Projection(bm.mod = Olive_niches,
                                  proj.name = 'Current',
                                  new.env = bioclim_current,
                                  models.chosen = 'all',
                                  metric.binary = c('TSS', 'ROC'),
                                  metric.filter = c('TSS', 'ROC'),
                                  build.clamping.mask = TRUE)
myBiomodProj
plot(myBiomodProj) #not enough RAMt

prova<- raster(paste("C:/Users/rocchetti/OneDrive - Cirad/Desktop/Olive_ENM/P.A/proj_current/proj_Current_P.A_ensemble.tif"))
plot(prova)




#################### Model ensemble models

get_built_models(Olive_niches , full.name = NULL, PA = NULL, run = NULL, algo = c('RF', 'GLM', 'MAXNET'))

Olive_EM <- BIOMOD_EnsembleModeling(bm.mod = Olive_niches,
                                    models.chosen = 'all',
                                    em.by = 'all',
                                    em.algo = c('EMwmean'),
                                    metric.select = c('TSS', 'ROC'),
                                    metric.select.thresh = c(0.85, 0.85),
                                    metric.eval = c('TSS', 'ROC'),
                                    var.import = 3,
                                    EMci.alpha = 0.05,
                                    EMwmean.decay = 'proportional')

bm_PlotVarImpBoxplot(bm.out = Olive_EM, group.by = c('expl.var', 'algo', 'merged.by.run'))

bm_PlotResponseCurves(bm.out = Olive_EM, 
                      models.chosen = get_built_models(Olive_EM)[c(1, 6, 7)],
                      fixed.var = 'median')

## reopen the model
Olive_EM <- get(load("/storage/replicated/cirad/projects/CLIMOLIVEMED/results/GenomicOffsets/Lorenzo/Leccino_new_genome24/ENM_medit/West.med.clim.soil/West.med.clim.soil.Single_Olive_niche_maxent.ensemble.models.out"))


########################### Project ensamble models
########## Current environmental variable

myBiomodEMProj   <- BIOMOD_EnsembleForecasting(bm.em = Olive_EM,
                                               proj.name = 'CurrentEM_ent',
                                               new.env = ras_current_var,
                                               models.chosen = 'all',
                                               metric.binary = 'all',
                                               metric.filter = 'all')
myBiomodEMProj
plot(myBiomodEMProj)

prova<- raster(paste("C:/Users/rocchetti/OneDrive - Cirad/Desktop/Olive_ENM/P.A/proj_current/proj_Current_P.A_ensemble.tif"))
plot(prova)

#------------------------------------------------------------
# Sensitivity analysis threshold 70th,80th, 90th percentile
#--------------------------------------------------------------

ENM<- raster(paste("D:/C/Desktop/Leccino24/ENM_403/proj_CurrentEM_ent_West.med.clim.soil_ensemble.tif"))
summary(values(ENM))

vals <- values(ENM)
vals <- vals[!is.na(vals)]
q70 <- quantile(vals, 0.70)
q80 <- quantile(vals, 0.80)  # your original threshold
q90 <- quantile(vals, 0.90)

mask70 <- ENM >= q70
mask80 <- ENM >= q80
mask90 <- ENM >= q90
# quick sanity check on how much area each retains
cellStats(mask70, sum) / cellStats(!is.na(ENM), sum)

cellStats(mask80, sum) / cellStats(!is.na(ENM), sum)

cellStats(mask90, sum) / cellStats(!is.na(ENM), sum)

extent(ENM)  # check this first — confirm it matches your intended region

tiff("ENM_threshold_sensitivity.tif", width = 12, height = 4.5, units = "in", res = 600, compression = "lzw")

par(mfrow = c(1,3))
plot(mask70, main = "70th percentile", legend = FALSE, col = c("grey90","forestgreen"),
     xlim = c(-15,10), ylim = c(27,55))
plot(mask80, main = "80th percentile (original)", legend = FALSE, col = c("grey90","forestgreen"),
     xlim = c(-15,10), ylim = c(27,55))
plot(mask90, main = "90th percentile", legend = FALSE, col = c("grey90","forestgreen"),
     xlim = c(-15,10), ylim = c(27,55))
par(mfrow = c(1,1))

dev.off()  # this line actually writes the file to disk








