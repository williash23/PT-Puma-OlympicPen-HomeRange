# Fit AKDE models for each individual-year  - Olympic Penninsula Mountain Lion home ranges

library(dplyr)
library(ctmm)
library(sf)
library(terra)

setwd("D:/Puma/Olympic_Peninsula/")



############################################################################
#  Load individual locations nad identifying data

load("dat.Rdata")


unq_inds <- unique(dat$name_year)

tst <- dat %>%
	group_by(name_year) %>%
	summarise(n_annual_locs = n(),
		annual_sampling_duration = max(new_date) - min(new_date)) %>%
	tidyr::separate(name_year, c("name", "year"), sep = "_") %>%
	as.data.frame() %>%
	dplyr::select(name, year, n_annual_locs,annual_sampling_duration)
	
	
	
############################################################################
#  Landmass boundary; used so UDs don't spill into water
bnd <- st_read("wa-state-4326/wa-state-4326.shp") %>%
	dplyr::select(name)
wa <- as(bnd, "Spatial")
centroids <- coordinates(wa)
x <- centroids[,1]
y <- centroids[,2]
wa <- SpatialPolygonsDataFrame(wa, data = data.frame(x = x, y = y))



############################################################################
#  Run CTMMs and AKDEs for each indivdiual-year and save all components to lists and 
#   data frame(each list element or data frame row is the output for an individula-year)

SVF_ls <- list()
all_mods_summ_ls <- list()
top_mod_ls <- list()
UD_ls <- list()
UD_df <- data.frame()


for(i in 1:length(unq_inds)){

	print(i)
	print(unq_inds[i])
	
	ind_dat <- dat %>%
		filter(name_year == unq_inds[i]) %>%
		distinct(new_date, .keep_all = TRUE) %>%
		mutate(name = name_year) %>%
		dplyr::select(name, new_date, date_time, latitude, longitude, elevation)
	nam <- ind_dat %>%
		tidyr::separate(name, c("name", "year"), sep = "_")
	
	ind_tel <- as.telemetry(ind_dat)
	SVF <- variogram(ind_tel)
	GUESS <- variogram.fit(SVF)
	all_mods <- ctmm.select(ind_tel, CTMM = GUESS, verbose = TRUE)
	print(summary(all_mods))

	top_mod <- all_mods[[1]]
	print(summary(top_mod))

	UD_top_mod <- akde(ind_tel, top_mod, SP = wa, SP.in = TRUE)
	print(summary(UD_top_mod))
	
	tmp <- data.frame(name = nam$name[1],
		year = nam$year[1],
		num_daily_locations = nrow(ind_dat),
		dur_locations = max(ind_dat$new_date) - min(ind_dat$new_date),
		model_DOF = summary(top_mod)$DOF[1],
		aKDE_sq_km_est = summary(UD_top_mod)$CI[2],
		aKDE_sq_km_LCI = summary(UD_top_mod)$CI[1],
		aKDE_sq_km_UCI = summary(UD_top_mod)$CI[3])
	
	rownames(tmp) <- NULL
	
	UD_df <- UD_df %>% rbind(tmp)
	
	SVF_ls[[i]] <- SVF
	all_mods_summ_ls[[i]] <- summary(all_mods)
	top_mod_ls[[i]] <- top_mod
	UD_ls[[i]] <- UD_top_mod
	
	}
	
#### ISSUE WITH individual 87 (forrestm12)
#  remove 87 from UD_ls
UD_ls <- UD_ls[[-87]]

#save(UD_df, file = "UD_df.Rdata")	
#save(SVF_ls, file = "SVF_ls.Rdata")
#save(all_mods_summ_ls, file = "all_mods_summ_ls.Rdata")
#save(top_mod_ls, file = "top_mod_ls.Rdata")
#save(UD_ls, file = "UD_ls.Rdata")


############################################################################
#  Run speed estimation
spd_df <- data.frame()
for(i in 1:length(top_mod_ls)){

	ind <- top_mod_ls[[i]]
	mod <- speed(ind)
	nam <-ind@info[[1]]
	tmp <- data.frame(name = nam,
			spd_DOF = mod$DOF,
			unit = rownames(mod$CI),
			spd_est = mod$CI[1,2],
			spd_LCI = mod$CI[1,1],
			spd_UCI = mod$CI[1,3]) %>%
		tidyr::separate(name, c("name", "year"), sep = "_")
	
	rownames(tmp) <- NULL
	
	spd_df <- spd_df %>% rbind(tmp)
	
	}
	
spd_df <- spd_df %>% 
	dplyr::filter(spd_est != "Inf")
#save(spd_df , file = "spd_df .Rdata")	


##################################################################################################
##################################################################################################