#  Prepare forest/non-forest boundary ("edge") covariate layer and summarize covariate value per home range

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
#  Processing forest/non-forest boundary ("edge")

#  Initial layer processing
# Tree canopy cover
# tree_canopy <- rast("nlcd_forest/nlcd_2016_treecanopy_2019_08_31.img")
# tree_canopy <- crop(tree_canopy, vect(st_transform(hr_extent, crs(tree_canopy))))
# tree_canopy <- project(tree_canopy, crs(lulc_human_mod))
# tree_canopy[tree_canopy > 100] <- NA
# names(tree_canopy) <- "tree_canopy"


#  Original tree canopy cover layer
forest <- rast("tree_canopy.tif")
forest <- crop(forest, project(hr_total_v, "epsg: 4326"))
forest_UTM <- project(forest, "epsg:32610")

#  Reclassify so that cells with value 30 or greater are "forested" and others are "non-forested"
forest_thresh30 <- forest_UTM
forest_thresh30[forest_thresh30 >= 30] <- 100
forest_thresh30[forest_thresh30 < 30] <- 1
#  Layer values: 100 = forest; 1 = non-forest

#  Reclassify to create a layer that shows non-forest only (for masking later on)
nonforest_thresh30 <- forest_thresh30
nonforest_thresh30[nonforest_thresh30 == 100] <- NA
#  Layer values: 1 = non-forest; all others are NA

#  Get distance from each forested cell to the closest non-forested cell
dist_forest_to_nonforest <- distance(forest_thresh30, target = 100)
#  Layer values: continuous 0 - 2823.31m

#  Mask out non-forested areas
dist_forest_to_nonforest_m <- mask(dist_forest_to_nonforest, nonforest_thresh30, inverse = TRUE)
#  Layer values: continuous 23 - 2823.31m

#  Reclassify to select cells that that have a distance from forest to closest non-forested cell is 
#   less than or equal to 60 m (ie, within 60 m means it is forest/non-forest edge)
dist_forest_to_nonforest_m[dist_forest_to_nonforest_m > 60] <- NA
dist_forest_to_nonforest_m[dist_forest_to_nonforest_m <= 60] <- 1
#  Layer values: 1 = boundary between forest/non-forest ("edge"); all others are NA


#  Write raster - Jamie will manually deal with werid edges (around lakes and beaches)
#writeRaster(dist_forest_to_nonforest_m, file = "edge_between_forested_and_nonforested.tif", overwrite = TRUE)


# Read back in layer from Jamie 
edge_r <- rast("edge_between_forested_and_nonforested_60m_servalfix.tif")

# Adjusted forested layer to values 1 and 0
forest_thresh30[forest_thresh30 == 1] <- 0
forest_thresh30[forest_thresh30 == 100] <- 1
#  Layer values: 1 = forest; 0 = non-forest



############################################################################
#  Get summarized length of edge within each home range
cov_df <- data.frame()
for(i in 1:nrow(hr_UTM)){
	tmp1 <- hr[i,] %>% 
		as.data.frame() %>%
		dplyr::select(-geometry)
	tmp2a <- crop(edge_r, vect(hr_UTM[i,]))
	tmp2b <- mask(tmp2a, vect(hr_UTM[i,]))
	tmp3 <- freq(tmp2b)$count
	tmp4 <- res(tmp2b)[1]
	tmp5 <- tmp3 * tmp4
	tmp6 <- tmp5/1000
	
	tmp7a <- crop(forest_thresh30, vect(hr_UTM[i,]))
	tmp7b <- mask(tmp7a, vect(hr_UTM[i,]))
	tmp8 <- values(tmp7b)
	tmp9 <-  res(tmp7b)[1] * res(tmp7b)[2]
	tmp10 <- sum(tmp8, na.rm = TRUE) 
	tmp11 <- tmp10 * tmp9
	tmp12 <- tmp11/1000000
	
	tmp13 <- tmp1 %>%
		cbind(edge_km = tmp6) %>%
		cbind(forest_sq_km = tmp12) %>%
		mutate(edge_forest_ratio = edge_km/forest_sq_km) # larger number means more edge
	cov_df <- cov_df %>%
		rbind(tmp13)
	}

edge_cov_dat <- cov_df %>%
	dplyr::select(-CL, -boundary, -name_year_sq_km)
#save(edge_cov_dat, file = "edge_cov_dat.Rdata")



############################################################################
# TESTING ALTERNATIVE METHOD USING 'landscapemetrics'
#   OUTPUT IS THE SAME
library(landscapemetrics)

tmp <- rast("E:/New_OP_HR/puma_waterDist_op.tif")
tmp2 <- tmp
tmp2[is.na(tmp2)] <- 0
tmp3 <- project(tmp2, forest_thresh30, method = "near")
tmp3[tmp3 > 5] <- 100
tmp3[tmp3 <= 5] <- 500
la <- mask(forest_thresh30, tmp3, maskvalue = 500, updatevalue = 2)
#0 = non-forest, 1 = forest, 2 = water

hr_la <- crop(la, vect(hr_total))
hr_la <- mask(hr_la, vect(hr_total))

patch_sd <- lsm_c_area_sd(hr_la)
patch_mn <- lsm_c_area_mn(hr_la)
patch_cv <- lsm_c_area_cv(hr_la)

#  Get summarized length of edge within each home range
cov_df <- data.frame()
for(i in 1:nrow(hr_UTM)){
	tmp1 <- hr[i,] %>% 
		as.data.frame() %>%
		dplyr::select(-geometry)
	#tmp2a <- crop(forest_thresh30, vect(hr_UTM[i,]))
	tmp2a <- crop(la, vect(hr_UTM[i,]))
	tmp2b <- mask(tmp2a, vect(hr_UTM[i,]))
	
	#tmp3 <- lsm_c_ed(tmp2b)
	#tmp4 <- lsm_l_ed(tmp2b)
	#tmp5 <- lsm_c_te(tmp2b)
	#tmp6 <- lsm_l_te(tmp2b)
	tmp7 <- lsm_c_ai(tmp2b)
	tmp8 <- lsm_c_area_mn(tmp2b)
	tmp9 <- lsm_c_ca(tmp2b)
	tmp10 <- lsm_l_ed(tmp2b)
	tmp11 <- lsm_c_pland(tmp2b)
	tmp12 <- lsm_c_np(tmp2b)


	
#  Equals ED = 0 if only one patch is present (and the landscape boundary is not included) and 
#   increases, without limit, as the landscapes becomes more patchy
# Equals TE = 0 if all cells are edge cells. Increases, without limit, as landscape becomes 
#   more fragmented
# AI is an 'Aggregation metric'. It equals the number of like adjacencies divided by the theoretical 
#  maximum possible number of like adjacencies for that class.
#  Equals 0 for maximally disaggregated and 100 for maximally aggregated classes.

	tmp13 <- tmp1 %>%
		cbind(lm_perc_nonforest = tmp11$value[1]) %>%
		cbind(lm_n_patch_nonforest = tmp12$value[1]) %>%
		cbind(lm_perc_forest = tmp11$value[2]) %>%
		cbind(lm_n_patch_forest = tmp12$value[2]) %>%
		cbind(lm_edge_density = tmp10$value[1]) %>% 
		cbind(lm_agg_index = tmp7$value[2]) %>% 
		cbind(lm_mean_nonforest_patch_area = tmp8$value[1]) %>%
		cbind(lm_total_nonforest_area = tmp9$value[1]) %>%
		mutate(patch_structure_nonforest = lm_mean_nonforest_patch_area/lm_total_nonforest_area) %>% # many small patches vs. few larges patches
		cbind(lm_mean_forest_patch_area = tmp8$value[2]) %>%
		cbind(lm_total_forest_area = tmp9$value[2]) %>%
		mutate(patch_structure_forest = lm_mean_forest_patch_area/lm_total_forest_area) # many small patches vs. few larges patches
	cov_df <- cov_df %>%
		rbind(tmp13)
	}


lm_edge_cov_dat <- cov_df %>%
	dplyr::select(-CL, -boundary, -name_year_sq_km)
save(lm_edge_cov_dat, file = "lm_edge_cov_dat.Rdata")


window_circular <- matrix(c(NA, 1, NA, 1, 1, 1, NA, 1, NA), nrow = 9, ncol = 9)
tst <- window_lsm(forest, window = window_circular, what = c("lsm_l_te"))

##################################################################################################
##################################################################################################