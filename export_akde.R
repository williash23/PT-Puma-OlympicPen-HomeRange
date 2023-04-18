#  Export AKDE for each indivdiual-year as raster and polygons; visual exploration

library(dplyr)
library(ctmm)
library(sf)
library(terra)

setwd("D:/Puma/Olympic_Peninsula/")



############################################################################
#  Load ctmm components

load("UD_df.Rdata")
load("UD_ls.Rdata")
load("top_mod_ls.Rdata")

#  Load individual identifying data
load("dat.Rdata")



############################################################################
#  Process ctmm components
#   Remove individuals that do not have enough locations or sufficient degrees of freedom to 
#   accept models

UD_df <- UD_df %>%
	mutate(dur_locations = as.numeric(dur_locations)) %>%
	mutate(name_year = paste(name, year, sep = "_"))
nrow(UD_df)

rem <- UD_df %>%
	dplyr::filter(model_DOF <= 5 | num_daily_locations <= 15) 
nrow(rem)

red_UD_df <- UD_df %>%
	dplyr::filter(num_daily_locations > 15) %>%
	dplyr::filter(model_DOF > 5) 
nrow(red_UD_df)
	
	
nam_ls <- list()
for(i in 1:length(UD_ls)){
	tmp1 <- UD_ls[[i]]
	tmp2 <- tmp1@info[[1]]
	nam_ls[[i]] <- tmp2
	}
	
list_order <- unlist(nam_ls)
df_order <- red_UD_df$name_year

kp <- which(list_order %in% red_UD_df$name_year == TRUE)

red_UD_ls <- UD_ls[kp]

red_top_mod_ls <- top_mod_ls[kp]


red_nam_ls <- list()
for(i in 1:length(red_UD_ls)){
	tmp1 <- red_UD_ls[[i]]
	tmp2 <- tmp1@info[[1]]
	red_nam_ls[[i]] <- tmp2
	}
	
red_list_order <- unlist(red_nam_ls)
table(df_order == red_list_order)




############################################################################
# Export AKDE rasters

# Template for raster
r <- rast("D:/Puma/Olympic_Peninsula/op/puma_elev_op.tif")

for(i in 1:length(red_UD_ls)){

	print(i)
	print(red_list_order[i])
	
	# 95% boundary (core)
	tmp1a1 <- as.sf(red_UD_ls[[i]], level.UD = 0.95) %>%
		tidyr::separate(name, c("name_year", "CL", "boundary"), sep = " ") 
	rownames(tmp1a1) <- NULL
	tmp1a2 <- tmp1a1 %>%
		slice(2)
	UD_sf <- tmp1a1 %>%
		st_transform(4326)
	save(UD_sf, file = paste(paste("akde_export/b95/sf/sf95", red_list_order[i], sep = "-"), ".Rdata", sep = ""))
	
	tmp1b <- raster(red_UD_ls[[i]], DF = "PDF", level.UD = 0.95)
	tmp1c <- raster::mask(tmp1b, tmp1a1)
	tmp1d <- rast(climateStability::rescale0to1(tmp1c))
	tmp1<- resample(tmp1d, r)
	writeRaster(tmp1, file = paste(paste("akde_export/b95/rasters/r95", red_list_order[i], sep = "-"), ".tif", sep = ""))
	
	rm(UD_sf)
	
	# 50% boundary (core)
	tmp2a1 <- as.sf(red_UD_ls[[i]], level.UD = 0.50) %>%
		tidyr::separate(name, c("name_year", "CL", "boundary"), sep = " ") 
	rownames(tmp2a1) <- NULL
	tmp2a2 <- tmp2a1 %>%
		slice(2)
	UD_sf <- tmp2a1 %>%
		st_transform(4326)
	save(UD_sf, file = paste(paste("akde_export/b50/sf/sf50", red_list_order[i], sep = "-"), ".Rdata", sep = ""))

	tmp2b <- raster(red_UD_ls[[i]], DF = "PDF", level.UD = 0.50)
	tmp2c <- raster::mask(tmp2b, tmp2a2)
	tmp2d <- rast(climateStability::rescale0to1(tmp2c))
	tmp2 <- resample(tmp2d, r)
	writeRaster(tmp2, file = paste(paste("akde_export/b50/rasters/r50", red_list_order[i], sep = "-"), ".tif", sep = ""))

	rm(UD_sf)
	}
	



############################################################################
# Extract speed estimates
	
spd_df <- data.frame()
for(i in 1:length(red_top_mod_ls)){

	txt1 <- strsplit(rownames(spd_ls[[i]]$CI), split = " ")[[1]][2]
	
	tmp <- data.frame(speed_DOF = spd_ls[[i]]$DOF[1],
		speed_est = spd_ls[[i]]$CI[2],
		speed_LCI = spd_ls[[i]]$CI[1],
		speed_UCI = spd_ls[[i]]$CI[3],
		units = txt1)
	
	rownames(tmp) <- NULL
	spd_df <- spd_df %>% rbind(tmp)
	}

spd_df[spd_df == "Inf"] <- NA

	
red_UD_spd_df <- red_UD_df %>%
	cbind(spd_df)
	

sx_dat <- read.csv("puma_names_sex.csv")

UD_spd_sx <- red_UD_spd_df %>%
	left_join(sx_dat, by = c("name" = "Names"))
#save(UD_spd_sx, file = "UD_spd_sex.Rdata")
	
	
	
	
############################################################################	
#  Population-level averages

kp1 <- which(red_UD_df$year %in% c(2011, 2012, 2013))
kp2 <- which(red_UD_df$year %in% c(2014, 2015, 2016, 2017))
kp3 <- which(red_UD_df$year %in% c(2018, 2019, 2020))


m <- meta(red_UD_ls, sort = TRUE)

early_UD_ls <- red_UD_ls[kp1]
mid_UD_ls <- red_UD_ls[kp2]
late_UD_ls <- red_UD_ls[kp3]

m_early <- meta(early_UD_ls, sort = TRUE)
m_mid <- meta(mid_UD_ls, sort = TRUE)
m_late <- meta(late_UD_ls, sort = TRUE)



############################################################################	
#  Exploratory plots

p_spd <- ggplot() +
	geom_jitter(data = m_dat_spd, width = 0.1, size = 2,
		aes(x = Sex, y = speed_est, colour = Sex, fill = Sex)) +
	geom_boxplot(data = m_dat_spd, alpha = 0.4,
		aes(x = Sex, y = speed_est, colour = Sex, fill = Sex)) +
	theme_bw() +
	theme(legend.position = "none") +
	ylab("Estimated daily speed (km)") +
	xlab("Sex") +
	scale_colour_manual(name = "Sex", values = c("#00A08A", "#F2AD00")) +
	scale_fill_manual(name = "Sex", values = c("#00A08A", "#F2AD00")) 
p_spd	
	

p_sze <- ggplot() +
	geom_jitter(data = m_dat_sze, width = 0.1, size = 2,
		aes(x = Sex, y = aKDE_sq_km_est, colour = Sex, fill = Sex)) +
	geom_boxplot(data = m_dat_sze, alpha = 0.4,
		aes(x = Sex, y = aKDE_sq_km_est, colour = Sex, fill = Sex)) +
	theme_bw() +
	theme(legend.position = "none") +
	ylab("Home range size (sq km)") +
	xlab("Sex") +
	scale_colour_manual(name = "Sex",values = c("#00A08A", "#F2AD00"))  +
	scale_fill_manual(name = "Sex", values = c("#00A08A", "#F2AD00")) 
p_sze	


p_sze + p_spd


# Outlier for speed is is butch 2013
plot(SVF_ls[[22]])
plot(SVF_ls[[22]], fraction = 0.1)

# Outlier for size is is crash 2015
plot(SVF_ls[[35]])
plot(SVF_ls[[35]], fraction = 0.1)


##################################################################################################
##################################################################################################