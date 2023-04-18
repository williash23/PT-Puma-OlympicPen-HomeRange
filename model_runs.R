#  Run models to evaluate effect of covariates on home range size and movement speed

library(dplyr)
library(MCMCglmm)
library(ggplot2)
library(patchwork)
library(tidyr)
library(sf)
library(broom.mixed)
library(dotwhisker)

setwd("D:/Puma/Olympic_Peninsula/Current_HR_analysis/")



############################################################################
load("hr.Rdata")
load("spd_df_w_boundary.Rdata")	
sex_dat <- read.csv("puma_names_sex.csv")
sex_dat$Sex[sex_dat$Sex == "Unknown"] <- "F"

load("human_mod_veg_cov_dat.Rdata")
load("edge_cov_dat.Rdata")
load("lm_edge_cov_dat.Rdata")



############################################################################
dat_tmp1 <- hr %>%
	dplyr::select(-CL, -boundary) %>%
	left_join(sex_dat, by = c("name" = "Names")) %>%
	left_join(spd_df) 
hist(dat_tmp1 $name_year_sq_km, breaks = 25)
hist(dat_tmp1 $spd_est)
# Crash 2015 HR is huge (alsmot 2400 sq km)
# Butch 2013 speed is 18 km/day

dat_tmp2 <- dat_tmp1 %>%
	dplyr::filter(name_year_sq_km < 1250) %>%
	rename(hr_est = name_year_sq_km) %>%
	dplyr::select(-unit)


# Join covariates and make scaled versions
dat_sz <- dat_tmp2 %>%
	left_join(human_mod_veg_cov_dat) %>%
	left_join(edge_cov_dat) %>%
	left_join(lm_edge_cov_dat) %>%
	as.data.frame() %>%
	dplyr::select(-spd_DOF, -spd_est, -spd_LCI, -spd_UCI, -geometry) 

dat_sz$Sex <- as.factor(dat_sz$Sex)
dat_sz$lm_perc_nonforest_sc <- as.numeric(scale(dat_sz$lm_perc_nonforest))
dat_sz$lm_perc_forest_sc <- as.numeric(scale(dat_sz$lm_perc_forest))
dat_sz$lm_edge_density_sc <- as.numeric(scale(dat_sz$lm_edge_density))
dat_sz$lm_agg_index_sc <- as.numeric(scale(dat_sz$lm_agg_index))
dat_sz$ratio_sc <- as.numeric(scale(dat_sz$edge_forest_ratio))
dat_sz$lulc_human_mod_sc <- as.numeric(scale(dat_sz$lulc_human_mod))
dat_sz$ghsl_human_mod_sc <- as.numeric(scale(dat_sz$ghsl_human_mod))
dat_sz$human_night_light_sc <- as.numeric(scale(dat_sz$human_night_light))
dat_sz$hr_est_lg <- log(dat_sz$hr_est)


dat_spd <- dat_tmp2 %>%
	dplyr::filter(spd_est < 16) %>%
	left_join(human_mod_veg_cov_dat) %>%
	left_join(edge_cov_dat) %>%
	left_join(lm_edge_cov_dat) %>%
	as.data.frame() %>%
	dplyr::select(-spd_DOF, -hr_est, -spd_LCI, -spd_UCI, -geometry) 

dat_spd$Sex <- as.factor(dat_spd$Sex)
dat_spd$lm_perc_nonforest_sc <- as.numeric(scale(dat_spd$lm_perc_nonforest))
dat_spd$lm_perc_forest_sc <- as.numeric(scale(dat_spd$lm_perc_forest))
dat_spd$lm_edge_density_sc <- as.numeric(scale(dat_spd$lm_edge_density))
dat_spd$lm_agg_index_sc <- as.numeric(scale(dat_spd$lm_agg_index))
dat_spd$ratio_sc <- as.numeric(scale(dat_spd$edge_forest_ratio))
dat_spd$lulc_human_mod_sc <- as.numeric(scale(dat_spd$lulc_human_mod))
dat_spd$ghsl_human_mod_sc <- as.numeric(scale(dat_spd$ghsl_human_mod))
dat_spd$human_night_light_sc <- as.numeric(scale(dat_spd$human_night_light))



############################################################################
#  Summaries
dat_sz %>% 	
	group_by(Sex) %>%
	mutate(mean_hr = mean(hr_est)) %>%
	slice(1) %>%
	as.data.frame() %>%
	dplyr::select(Sex, mean_hr)
dat_spd %>% 	
	group_by(Sex) %>%
	mutate(mean_spd = mean(spd_est, na.rm = TRUE)) %>%
	slice(1) %>%
	as.data.frame() %>%
	dplyr::select(Sex, mean_spd)

dat_f_hr <- dat_sz %>% dplyr::filter(Sex == "F") %>% dplyr::filter(!is.na(hr_est))
dat_f_spd <- dat_spd %>% dplyr::filter(Sex == "F") 

dat_m_hr <- dat_sz %>% dplyr::filter(Sex == "M") %>% dplyr::filter(!is.na(hr_est))
dat_m_spd <- dat_spd %>% dplyr::filter(Sex == "M") 



############################################################################
#  Models -- AREA

sz_cor_dat <- dat_sz %>%
	dplyr::select(lm_perc_forest, lm_perc_nonforest,
		edge_forest_ratio, lm_edge_density, lm_agg_index, 
		lulc_human_mod, ghsl_human_mod, human_night_light)
colnames(sz_cor_dat) <- c("Percent forested", "Percent non-forested",
	"Edge density (lm)", "Edge density (manual)", "Aggregation index",
	"Landscape modification", "Urbanization", "Nighttime light intensity")

corrplot::corrplot(cor(sz_cor_dat), method = "number")

#### edge_forest_ratio = edge_km/forest_sq_km 
#### ratio_sc: larger number means more edge relative to forested area

m0 <- MCMCglmm(hr_est_lg ~ Sex, random = ~name, data = dat_sz, nitt = 2500, burnin = 500, thin = 3)
m1 <- MCMCglmm(hr_est_lg ~ ratio_sc * Sex, random = ~name, data = dat_sz, nitt = 2500, burnin = 500, thin = 3)
m2 <- MCMCglmm(hr_est_lg ~ lm_edge_density_sc * Sex, random = ~name, data = dat_sz, nitt = 2500, burnin = 500, thin = 3)
m3 <- MCMCglmm(hr_est_lg ~ lm_perc_nonforest_sc * Sex, random = ~name, data = dat_sz, nitt = 2500, burnin = 500, thin = 3)
m4 <- MCMCglmm(hr_est_lg ~ lm_perc_forest_sc * Sex, random = ~name, data = dat_sz, nitt = 2500, burnin = 500, thin = 3)
m5 <- MCMCglmm(hr_est_lg ~ lm_agg_index_sc * Sex, random = ~name, data = dat_sz, nitt = 2500, burnin = 500, thin = 3)
m6 <- MCMCglmm(hr_est_lg ~ ghsl_human_mod_sc *Sex, random = ~name, data = dat_sz, nitt = 2500, burnin = 500, thin = 3)
m7 <- MCMCglmm(hr_est_lg ~ lulc_human_mod_sc * Sex, random = ~name, data = dat_sz, nitt = 2500, burnin = 500, thin = 3)
m8 <- MCMCglmm(hr_est_lg ~ human_night_light_sc * Sex, random = ~name, data = dat_sz, nitt = 2500, burnin = 500, thin = 3)

# Global	
g1 <- MCMCglmm(hr_est_lg ~ ratio_sc + 
	ghsl_human_mod_sc * Sex, 
	random = ~name, data = dat_sz, 
	nitt = 2500, burnin = 500, thin = 3)
summary(g1)

tt <- tidy(g1)  %>%
    filter(effect == "fixed")
	
dwplot(tt, dot_args = list(size = 3), line_args = list(lwd = 2)) +
	geom_vline(xintercept = 0, lty = 2) +
	xlab("Coefficient") +
	theme_bw() + 
	theme(text = element_text(size = 14),
		legend.position = "none") +
	scale_y_discrete("Parameter", 
		limits = c("ghsl_human_mod_sc:SexM", "ghsl_human_mod_sc", "ratio_sc", "SexM"),
		labels = c("Urbanization * Sex (male)", "Urbanization", "Edge density", "Sex (male)"))



############################################################################
#  Models -- SPEED

spd_cor_dat <- dat_spd %>%
	dplyr::select(lm_perc_forest, lm_perc_nonforest,
		edge_forest_ratio, lm_edge_density, lm_agg_index, 
		lulc_human_mod, ghsl_human_mod, human_night_light)
colnames(spd_cor_dat) <- c("Percent forested", "Percent non-forested",
	"Edge density (lm)", "Edge density (manual)", "Aggregation index",
	"Landscape modification", "Urbanization", "Nighttime light intensity")

corrplot::corrplot(cor(spd_cor_dat), method = "number")

#### edge_forest_ratio = edge_km/forest_sq_km 
#### ratio_sc: larger number means more edge relative to forested area

m0 <- MCMCglmm(spd_est ~ Sex, random = ~name, data = dat_spd, nitt = 2500, burnin = 500, thin = 3)
m1 <- MCMCglmm(spd_est ~ ratio_sc * Sex, random = ~name, data = dat_spd, nitt = 2500, burnin = 500, thin = 3)
m2 <- MCMCglmm(spd_est ~ lm_edge_density_sc * Sex, random = ~name, data = dat_spd, nitt = 2500, burnin = 500, thin = 3)
m3 <- MCMCglmm(spd_est ~ lm_perc_nonforest_sc * Sex, random = ~name, data = dat_spd, nitt = 2500, burnin = 500, thin = 3)
m4 <- MCMCglmm(spd_est ~ lm_perc_forest_sc * Sex, random = ~name, data = dat_spd, nitt = 2500, burnin = 500, thin = 3)
m5 <- MCMCglmm(spd_est ~ lm_agg_index_sc * Sex, random = ~name, data = dat_spd, nitt = 2500, burnin = 500, thin = 3)
m6 <- MCMCglmm(spd_est ~ ghsl_human_mod_sc * Sex, random = ~name, data = dat_spd, nitt = 2500, burnin = 500, thin = 3)
m7 <- MCMCglmm(spd_est ~ lulc_human_mod_sc * Sex, random = ~name, data = dat_spd, nitt = 2500, burnin = 500, thin = 3)
m8 <- MCMCglmm(spd_est ~ human_night_light_sc * Sex, random = ~name, data = dat_spd, nitt = 2500, burnin = 500, thin = 3)
#  No interaction terms are significant

m0 <- MCMCglmm(spd_est ~ Sex, random = ~name, data = dat_spd, nitt = 2500, burnin = 500, thin = 3)
m1 <- MCMCglmm(spd_est ~ ratio_sc + Sex, random = ~name, data = dat_spd, nitt = 2500, burnin = 500, thin = 3)
m2 <- MCMCglmm(spd_est ~ lm_edge_density_sc + Sex, random = ~name, data = dat_spd, nitt = 2500, burnin = 500, thin = 3)
m3 <- MCMCglmm(spd_est ~ lm_perc_nonforest_sc + Sex, random = ~name, data = dat_spd, nitt = 2500, burnin = 500, thin = 3)
m4 <- MCMCglmm(spd_est ~ lm_perc_forest_sc + Sex, random = ~name, data = dat_spd, nitt = 2500, burnin = 500, thin = 3)
m5 <- MCMCglmm(spd_est ~ lm_agg_index_sc + Sex, random = ~name, data = dat_spd, nitt = 2500, burnin = 500, thin = 3)
m6 <- MCMCglmm(spd_est ~ ghsl_human_mod_sc +Sex, random = ~name, data = dat_spd, nitt = 2500, burnin = 500, thin = 3)
m7 <- MCMCglmm(spd_est ~ lulc_human_mod_sc + Sex, random = ~name, data = dat_spd, nitt = 2500, burnin = 500, thin = 3)
m8 <- MCMCglmm(spd_est ~ human_night_light_sc + Sex, random = ~name, data = dat_spd, nitt = 2500, burnin = 500, thin = 3)
	
	
tt <- tidy(m7)  %>%
    filter(effect == "fixed")
	
dwplot(tt, dot_args = list(size = 3), line_args = list(lwd = 2)) +
	geom_vline(xintercept = 0, lty = 2) +
	xlab("Coefficient") +
	theme_bw() + 
	xlim(-2, 5) +
	theme(text = element_text(size = 14),
		legend.position = "none") +
	scale_y_discrete("Parameter", 
		limits = c("lulc_human_mod_sc", "SexM"),
		labels = c("Landscape modification", "Sex (male)"))



############################################################################
# Covariate plots

ggplot() +
	geom_smooth(data = dat_sz, aes(x = ghsl_human_mod_sc, y =  hr_est_lg, color = Sex, fill = Sex), method = "lm") +
	theme_bw() +
	scale_colour_manual(name = "Sex", values = c("#0B775E", "#35274A")) +
	scale_fill_manual(name = "Sex", values = c("#0B775E", "#35274A")) +
	xlab("Urbanization (Z-transformed)") +
	ylab("Log-transformed home range (sq km)") +
	theme(text = element_text(size = 14))
	
ggplot() +
	geom_smooth(data = dat_sz, aes(x = ratio_sc, y =  hr_est_lg), method = "lm", color = "black") +
	theme_bw() +
	xlab("Edge density (manual, Z-transformed)") +
	ylab("Log-transformed home range (sq km)") +
	theme(text = element_text(size = 14))
	
ggplot() +
	geom_smooth(data = dat_spd, aes(x = lulc_human_mod_sc, y =  spd_est), method = "lm", color = "black") +
	theme_bw() +
	xlab("Landscape modificiation (Z-transformed)") +
	ylab("Speed (km/day)") +
	theme(text = element_text(size = 14))



p1 <- ggplot() +
	geom_point(data = dat_sz, aes(x = lm_perc_nonforest, y = hr_est_lg)) +
	geom_smooth(data = dat_sz, aes(x = lm_perc_nonforest, y =  hr_est_lg), method = "lm") +
	theme_bw() +
	xlab("Perc.non-forested") +
	ylab("Log-transformed home range (sq km)")


p2 <- ggplot() +
	geom_point(data = dat_sz, aes(x = lm_perc_forest, y = hr_est_lg)) +
	geom_smooth(data = dat_sz, aes(x = lm_perc_forest, y =  hr_est_lg), method = "lm") +
	theme_bw() +
	xlab("Perc. forested") +
	ylab("Log-transformed home range (sq km)")

p3 <- ggplot() +
	geom_point(data = dat_sz, aes(x = lm_edge_density, y = hr_est_lg)) +
	geom_smooth(data = dat_sz, aes(x = lm_edge_density, y =  hr_est_lg), method = "lm") +
	theme_bw() +
	xlab("Edge density") +
	ylab("Log-transformed home range (sq km)")


p4 <- ggplot() +
	geom_point(data = dat_sz, aes(x = lm_perc_nonforest, y = hr_est_lg, color = Sex)) +
	geom_smooth(data = dat_sz, aes(x = lm_perc_nonforest, y =  hr_est_lg, color = Sex, fill = Sex), method = "lm") +
	theme_bw() +
	scale_colour_manual(name = "Sex", values = c("#0B775E", "#35274A")) +
	scale_fill_manual(name = "Sex", values = c("#0B775E", "#35274A")) +
	xlab("Perc.non-forested") +
	ylab("Log-transformed home range (sq km)")

p5 <- ggplot() +
	geom_point(data = dat_sz, aes(x = lm_perc_forest, y = hr_est_lg, color = Sex)) +
	geom_smooth(data = dat_sz, aes(x = lm_perc_forest, y =  hr_est_lg, color = Sex, fill = Sex), method = "lm") +
	theme_bw() +
	scale_colour_manual(name = "Sex", values = c("#0B775E", "#35274A")) +
	scale_fill_manual(name = "Sex", values = c("#0B775E", "#35274A")) +
	xlab("Perc. forested") +
	ylab("Log-transformed home range (sq km)")

p6 <- ggplot() +
	geom_point(data = dat_sz, aes(x = lm_edge_density, y = hr_est_lg, color = Sex)) +
	geom_smooth(data = dat_sz, aes(x = lm_edge_density, y =  hr_est_lg, color = Sex, fill = Sex), method = "lm") +
	theme_bw() +
	scale_colour_manual(name = "Sex", values = c("#0B775E", "#35274A")) +
	scale_fill_manual(name = "Sex", values = c("#0B775E", "#35274A")) +
	xlab("Edge density") +
	ylab("Log-transformed home range (sq km)")

p7 <- ggplot() +
	geom_point(data = dat_sz, aes(y = lm_edge_density, x = lm_perc_nonforest)) +
	geom_smooth(data = dat_sz, aes(y = lm_edge_density, x = lm_perc_nonforest), method = "lm") +
	theme_bw() +
	xlab("Perc.non-forested") +
	ylab("Edge density") +
	theme(text = element_text(size = 14))

	
p8 <- ggplot() +
	geom_point(data = dat_sz, aes(y = lm_edge_density, x = lm_perc_nonforest, color = Sex)) +
	geom_smooth(data = dat_sz, aes(y = lm_edge_density, x = lm_perc_nonforest, color = Sex, fill = Sex), method = "lm") +
	theme_bw() +
	scale_colour_manual(name = "Sex", values = c("#0B775E", "#35274A")) +
	scale_fill_manual(name = "Sex", values = c("#0B775E", "#35274A")) +
	xlab("Perc. non-forested") +
	ylab("Edge density") +
	theme(text = element_text(size = 14))



ggplot() +
	geom_point(data = dat_sz, aes(y = lm_mean_nonforest_patch_area, x = lm_perc_nonforest)) +
	#geom_smooth(data = dat_sz, aes(y = patch_structure_nonforest, x = lm_perc_nonforest), method = "lm") +
	theme_bw() +
	xlab("Perc. non-forested") +
	ylab("Patch size")
	

ggplot() +
	geom_point(data = dat_sz, aes(y = lm_n_patch_nonforest, x = lm_perc_nonforest)) +
	#geom_smooth(data = dat_sz, aes(y =  lm_n_patch_nonforest,, x = lm_perc_nonforest), 
		#formula = y ~ poly(x,2), method = "lm") +
		#method = "lm") +
	theme_bw() +
	ylim(0, 4750) +
	xlab("Perc. non-forested") +
	ylab("Number of non-forested patches") +
	theme(text = element_text(size = 14))

ggplot() +
	geom_point(data = dat_sz, aes(y = lm_n_patch_nonforest, x = lm_edge_density)) +
	#geom_smooth(data = dat_sz, aes(y =  lm_n_patch_nonforest,, x = lm_edge_density), 
		#formula = y ~ poly(x,2), method = "lm") +
		#method = "lm") +
	theme_bw() +
	ylim(0, 4750) +
	xlab("Edge density") +
	ylab("Number of non-forested patches") +
	theme(text = element_text(size = 14))


ggplot() +
	geom_point(data = dat_sz, aes(y = lm_n_patch_nonforest, x = lm_perc_nonforest, color = Sex)) +
	geom_smooth(data = dat_sz, aes(y =  lm_n_patch_nonforest,, x = lm_perc_nonforest, color = Sex, fill = Sex), 
		#formula = y ~ poly(x,2), method = "lm") +
		method = "lm") +
	theme_bw() +
	scale_colour_manual(name = "Sex", values = c("#0B775E", "#35274A")) +
	scale_fill_manual(name = "Sex", values = c("#0B775E", "#35274A")) +
	xlab("Perc. non-forested") +
	ylab("Number of non-forested patches")


ggplot() +
	geom_point(data = dat_sz, aes(y = lm_n_patch_nonforest, x = lm_total_nonforest_area)) +
	geom_smooth(data = dat_sz, aes(y =  lm_n_patch_nonforest, x = lm_total_nonforest_area), 
		#formula = y ~ poly(x,2), method = "lm") +
		method = "lm") +
	theme_bw() +
	xlab("Total area non-forested") +
	ylab("Number of patches")
	
p9 <- ggplot() +
	geom_point(data = dat_sz, aes(y = patch_structure_nonforest, x = lm_perc_nonforest)) +
	geom_smooth(data = dat_sz, aes(y = patch_structure_nonforest, x = lm_perc_nonforest), method = "lm") +
	theme_bw() +
	xlab("Perc.non-forested") +
	ylab("Patch structure")

p10 <- ggplot() +
	geom_point(data = dat_sz, aes(y = patch_structure_nonforest, x = lm_perc_nonforest, color = Sex)) +
	geom_smooth(data = dat_sz, aes(y = patch_structure_nonforest, x = lm_perc_nonforest, color = Sex, fill = Sex), method = "lm") +
	theme_bw() +
	scale_colour_manual(name = "Sex", values = c("#0B775E", "#35274A")) +
	scale_fill_manual(name = "Sex", values = c("#0B775E", "#35274A")) +
	xlab("Perc.non-forested") +
	ylab("Patch structure")
	
(p2  + p5) / (p1 + p4)

(p7 + p8)  / (p3 + p6)
 

ggplot() +
	geom_histogram(data = dat_sz, aes(x = lm_perc_forest), color = "black", fill = "grey50", bins = 25)+
	theme_bw() +
	#scale_colour_manual(name = "Sex", values = c("#0B775E", "#35274A")) +
	#scale_fill_manual(name = "Sex", values = c("#0B775E", "#35274A")) +
	xlab("Perc. forested")


p1 <- ggplot() +
	geom_point(data = dat_spd, aes(x = lm_perc_nonforest, y = spd_est)) +
	geom_smooth(data = dat_spd, aes(x = lm_perc_nonforest, y =  spd_est), method = "lm") +
	theme_bw() +
	xlab("Perc.non-forested") +
	ylab("Speed (km/day)")

p2 <- ggplot() +
	geom_point(data = dat_spd, aes(x = lm_perc_forest, y = spd_est)) +
	geom_smooth(data = dat_spd, aes(x = lm_perc_forest, y =  spd_est), method = "lm") +
	theme_bw() +
	xlab("Perc. forested") +
	ylab("Speed (km/day)")

p3 <- ggplot() +
	geom_point(data = dat_spd, aes(x = lm_edge_density, y = spd_est)) +
	geom_smooth(data = dat_spd, aes(x = lm_edge_density, y =  spd_est), method = "lm") +
	theme_bw() +
	xlab("Edge density") +
	ylab("Speed (km/day)")

p4 <- ggplot() +
	geom_point(data = dat_spd, aes(x = lm_perc_nonforest, y = spd_est, color = Sex)) +
	geom_smooth(data = dat_spd, aes(x = lm_perc_nonforest, y =  spd_est, color = Sex, fill = Sex), method = "lm") +
	theme_bw() +
	scale_colour_manual(name = "Sex", values = c("#0B775E", "#35274A")) +
	scale_fill_manual(name = "Sex", values = c("#0B775E", "#35274A")) +
	xlab("Perc.non-forested") +
	ylab("Speed (km/day)")

p5 <- ggplot() +
	geom_point(data = dat_spd, aes(x = lm_perc_forest, y = spd_est, color = Sex)) +
	geom_smooth(data = dat_spd, aes(x = lm_perc_forest, y =  spd_est, color = Sex, fill = Sex), method = "lm") +
	theme_bw() +
	scale_colour_manual(name = "Sex", values = c("#0B775E", "#35274A")) +
	scale_fill_manual(name = "Sex", values = c("#0B775E", "#35274A")) +
	xlab("Perc. forested") +
	ylab("Speed (km/day)")

p6 <- ggplot() +
	geom_point(data = dat_spd, aes(x = lm_edge_density, y = spd_est, color = Sex)) +
	geom_smooth(data = dat_spd, aes(x = lm_edge_density, y =  spd_est, color = Sex, fill = Sex), method = "lm") +
	theme_bw() +
	scale_colour_manual(name = "Sex", values = c("#0B775E", "#35274A")) +
	scale_fill_manual(name = "Sex", values = c("#0B775E", "#35274A")) +
	xlab("Edge density") +
	ylab("Speed (km/day)")
	
(p2  + p5) / (p1 + p4) / (p3 + p6) 

##################################################################################################
##################################################################################################