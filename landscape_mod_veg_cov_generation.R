#  Prepare annually-varying human modifcation, vegetation cover and human night light intensity 
#    covariate layers 

library(dplyr)
library(sf)
library(terra)
library(ggplot2)
library(tidyterra)

setwd("D:/Puma/Olympic_Peninsula/Current_HR_analysis/")



############################################################################
load("hr_extent.Rdata")
load("hr.Rdata")
load("hr_core.Rdata")

hr_UTM <- st_transform(hr, 32610)
hr_total <- st_union(hr_UTM) %>%
	st_buffer(50000)

hr_v <- vect(hr_UTM)
hr_total_v <- vect(hr_total)



############################################################################
#  Process covariates to correct extent and projection

# Human modification land cover
lulc_human_mod <- rast("lulc/lulc-human-modification-terrestrial-systems_geographic.tif")
lulc_human_mod <- crop(lulc_human_mod, hr_extent)
names(lulc_human_mod) <- "lulc_human_mod"
# The Global Human Modification of Terrestrial Systems data set provides a cumulative 
#  measure of the human modification of terrestrial lands across the globe at a 1-km resolution. 
#  It is a continuous 0-1 metric that reflects the proportion of a landscape modified, based on 
#  modeling the physical extents of 13 anthropogenic stressors and their estimated impacts 
#  using spatially-explicit global data sets with a median year of 2016.


# Human modification - ghls built up
ghsl_human_mod <- rast("ghsl-population-built-up-estimates-degree-urban-smod_smod-9ss-2015.tif")
ghsl_human_mod <- crop(ghsl_human_mod, hr_extent)
names(ghsl_human_mod) <- "ghsl_human_mod"
# GHS-BUILT describes the percent built-up area for each 9 arc-second grid cell 
#  (based on Landsat imagery  -- 2015


#  Vegetation cover
#   Perc tree canopy
load("ptc.Rdata")
ptc <- raster_ts
ptc <- rast(ptc)
ptc[ptc > 100] <- NA
names(ptc) <- rep("ptc", 10)
ptc <- resample(ptc, hnl)

#   Perc non-tree vegetation
load("pntv.Rdata")
pntv <- raster_ts
pntv <- rast(pntv)
pntv[pntv > 100] <- NA
names(pntv) <- rep("pntv", 10)
pntv <- resample(pntv, hnl)

#   Perc non-vegetation (bare)
load("pnv.Rdata")
pnv <- raster_ts
pnv <- rast(pnv)
pnv[pnv > 100] <- NA
names(pnv) <- rep("pnv", 10)
pnv <- resample(pnv, hnl)

# Human night light
tmp <- rast("tc.tif")
hnl_ls <- list.files("D:/Puma/Olympic_Peninsula/Current_HR_analysis/hnl")
setwd("D:/Puma/Olympic_Peninsula/Current_HR_analysis/hnl")

for(i in 1:length(hnl_ls)){

	tmp1 <- rast(hnl_ls[[i]])
	tmp2 <- resample(tmp1, tmp)
	nam <- paste("hnl", i, sep = "_")
	assign(nam, tmp2)
	}

human_night_light <- c(hnl_1, hnl_2, hnl_3, hnl_4, hnl_5, hnl_6, hnl_7, hnl_8, hnl_9, hnl_10)
names(human_night_light) <- rep("human_night_light", 10)


#  Save all new covariates
setwd("D:/Puma/Olympic_Peninsula/Current_HR_analysis/")


# writeRaster(lc, file = "lulc_human_mod.tif.tif")
# writeRaster(tree_canopy, file = "tree_canopy.tif")
# writeRaster(ghsl_human_mod, file = "ghsl_human_mod.tif")
# writeRaster(human_night_light, file = "human_night_light.tif")
# writeRaster(ptc, file = "ptc.tif")
# writeRaster(pntv, file = "pntv.tif")
# writeRaster(pnv, file = "pnv.tif")



############################################################################
#  Get summarized covariate value within each home range

yr_df <- data.frame(year = seq(from = 2011, to = 2020, by = 1),
	year_ind = (1:10))

cov_df <- data.frame()
for(i in 1:nrow(hr)){

    tmp1 <- hr[i,] %>% 
        as.data.frame() %>%
        dplyr::select(-geometry)
	yr <- hr$year[i]
	ind <- yr_df %>% 
		dplyr::filter(year == yr)
    tmp2a <- terra::extract(human_night_light[[ind$year_ind]], hr[i,], fun = "mean", na.rm = TRUE) %>%
		dplyr::select(-ID) 
	tmp2b <- terra::extract(ptc[[ind$year_ind]], hr[i,], fun = "mean", na.rm = TRUE) %>%
		dplyr::select(-ID) 
	tmp2c <- terra::extract(pntv[[ind$year_ind]], hr[i,], fun = "mean", na.rm = TRUE) %>%
		dplyr::select(-ID) 
	tmp2d <- terra::extract(pnv[[ind$year_ind]], hr[i,], fun = "mean", na.rm = TRUE) %>%
		dplyr::select(-ID) 
	tmp2e <- terra::extract(lulc_human_mod, hr[i,], fun = "mean", na.rm = TRUE) %>%
		dplyr::select(-ID) 
	tmp2f <- terra::extract(ghsl_human_mod, hr[i,], fun = "mean", na.rm = TRUE) %>%
		dplyr::select(-ID) 

	tmp3 <- tmp1 %>%
        cbind(tmp2a) %>%
        cbind(tmp2b) %>% 
        cbind(tmp2c) %>%
        cbind(tmp2d) %>%
        cbind(tmp2e) %>% 
        cbind(tmp2f) 
			
    cov_df <- cov_df %>%
        rbind(tmp3)
   }

human_mod_veg_cov_dat <- cov_df %>%
	dplyr::select(-CL, -boundary, -name_year_sq_km)
#save(human_mod_veg_cov_dat, file = "human_mod_veg_cov_dat.Rdata")


##################################################################################################
##################################################################################################