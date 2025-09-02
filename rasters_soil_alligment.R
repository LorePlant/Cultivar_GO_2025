library(raster)
library("readxl")

bio2<- raster(paste("D:/raster files/ENM_bioclim_soil_25/bio2_ENM_def_clip.tif"))
bio10<- raster(paste("D:/raster files/ENM_bioclim_soil_25/bio10_ENM_def_clip.tif"))
bio11<- raster(paste("D:/raster files/ENM_bioclim_soil_25/bio11_ENM_def_clip.tif"))
bio15<- raster(paste("D:/raster files/ENM_bioclim_soil_25/bio15_ENM_def_clip.tif"))
bio18<- raster(paste("D:/raster files/ENM_bioclim_soil_25/bio18_ENM_def_clip.tif"))
bio19<- raster(paste("D:/raster files/ENM_bioclim_soil_25/bio19_ENM_def_clip.tif"))
soilN<- raster(paste("D:/raster files/ENM_bioclim_soil_25/soilN_ENM_def_clip.tif"))
soilpH<- raster(paste("D:/raster files/ENM_bioclim_soil_25/soilpH_ENM_def_clip.tif"))

soilclay<- raster(paste("D:/raster files/ENM_bioclim_soil_25/soilclay_ENM_def_clip.tif"))
soilsand<- raster(paste("D:/raster files/ENM_bioclim_soil_25/soilsand_ENM_def_clip.tif"))

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
ras_current_var<-stack(c(bio2,bio10, bio11, bio15, bio18, bio19, soilN, soilpH, soilclay, soilsand))
plot(ras_current_var)




bio2<- raster(paste("D:/raster files/Env_var_west_med/bio2_west_med.tif"))
bio10<- raster(paste("D:/raster files/Env_var_west_med/bio10_west_med.tif"))
bio11<- raster(paste("D:/raster files/Env_var_west_med/bio11_west_med.tif"))
bio15<- raster(paste("D:/raster files/Env_var_west_med/bio15_west_med.tif"))
bio18<- raster(paste("D:/raster files/Env_var_west_med/bio18_west_med.tif"))
bio19<- raster(paste("D:/raster files/Env_var_west_med/bio19_west_med.tif"))
soilN<- raster(paste("D:/raster files/resampled_soilN_west_med.tif"))
soilpH<- raster(paste("D:/raster files/resampled_soilpH_west_med.tif"))

soilclay<- raster(paste("D:/raster files/resampled_soilclay_west_med.tif"))
soilsand<- raster(paste("D:/raster files/resampled_soilsand_west_med.tif"))

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
writeRaster(soilN, "D:/raster files/resampled_soilN_west_med.tif", format="GTiff", overwrite=TRUE)


soilpH <- resample(soilpH, bio2, method="bilinear")
writeRaster(soilpH, "D:/raster files/resampled_soilpH_west_med.tif", format="GTiff", overwrite=TRUE)

soilclay <- resample(soilclay, bio2, method="bilinear")
writeRaster(soilclay, "D:/raster files/resampled_soilclay_west_med.tif", format="GTiff", overwrite=TRUE)

soilsand <- resample(soilsand, bio2, method="bilinear")
writeRaster(soilsand, "D:/raster files/resampled_soilsand_west_med.tif", format="GTiff", overwrite=TRUE)


#stack the different raster file
ras_current_var<-stack(c(bio2,bio10, bio11, bio15, bio18, bio19, soilN, soilpH, soilclay, soilsand))
plot(ras_current_var)
