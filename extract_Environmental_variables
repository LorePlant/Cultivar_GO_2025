# exctact_environmental_variable
The following code show how to extract environmental variable from raster file for specific georeferenced point.
In this example I'm going to exctract 10 different bioclimatic variable from CHELSA raster

```
library(raster)
library("readxl")


bio2<- raster(paste("/storage/replicated/cirad/projects/CLIMOLIVEMED/results/GenomicOffsets/Lorenzo/Current_biovariable/West_Med_raster/bio2_studyarea_ext.tif"))
bio10<- raster(paste("/storage/replicated/cirad/projects/CLIMOLIVEMED/results/GenomicOffsets/Lorenzo/Current_biovariable/West_Med_raster/bio10_studyarea_ext.tif"))
bio11<- raster(paste("/storage/replicated/cirad/projects/CLIMOLIVEMED/results/GenomicOffsets/Lorenzo/Current_biovariable/West_Med_raster/bio11_studyarea_ext.tif"))
bio15<- raster(paste("/storage/replicated/cirad/projects/CLIMOLIVEMED/results/GenomicOffsets/Lorenzo/Current_biovariable/West_Med_raster/bio15_studyarea_ext.tif"))
bio18<- raster(paste("/storage/replicated/cirad/projects/CLIMOLIVEMED/results/GenomicOffsets/Lorenzo/Current_biovariable/West_Med_raster/bio18_studyarea_ext.tif"))
bio19<- raster(paste("/storage/replicated/cirad/projects/CLIMOLIVEMED/results/GenomicOffsets/Lorenzo/Current_biovariable/West_Med_raster/bio19_studyarea_ext.tif"))
names(bio1) = 'bio1'
names(bio2) = 'bio2'
names(bio4) = 'bio4'
names(bio5) = 'bio5'
names(bio6) = 'bio6'
names(bio8) = 'bio8'
names(bio9) = 'bio9'
names(bio10) = 'bio10'
names(bio11) = 'bio11'
names(bio12) = 'bio12'
names(bio14) = 'bio14'
names(bio15) = 'bio15'
names(bio18) = 'bio18'
names(bio19) = 'bio19'
#stack the different raster file
ras_current<-stack(c(bio1, bio2, bio4, bio5, bio6, bio8, bio9, bio10, bio11, bio12, bio14, bio15, bio18, bio19))
```

We can create a spatial grid x y with cell centered coordinates. 
The spatial grid is like a data frame with x long and y lat longitude derived from the center of each pixel and the value derived from all the raster layers
```
#spatial grid

coord_r<-rasterToPoints(ras_current, spatial = TRUE)
map_pts<-data.frame(x = coordinates(coord_r)[,1], y=coordinates(coord_r)[,2], coord_r@data)

```

upload the excel file with longitude (x) and latitude (y) data. In this case present in column 9 and 10 respectively

```
locations = read_excel("359_geo_info.xlsx")

## Extracting environmental values for each source population
Env <- data.frame(extract(ras_current, locations[,9:10]))
write.table(Env, "Current_clim_wise")

```

