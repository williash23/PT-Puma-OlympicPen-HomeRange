#  Initial data processing for AKDE analysis - Olympic Penninsula Mountain Lion home ranges

library(dplyr)
library(lubridate)
library(tidyr)
library(stringr)

setwd("D:/Puma/Olympic_Peninsula/")



############################################################################
#  Clean location (GPS/VHF) data and format for AKDE

dat1 <- read.csv("ohr_gps.csv") 


#  Clean up date formatting
dat2 <- dat1 

#  Is an entry with month/day/year or is an entry with day.month.year
tst_date <- dat2$date
date_is_DMY <- str_detect(tst_date, "\\.")
date_is_DMY[date_is_DMY == TRUE] <- 1
date_is_DMY[date_is_DMY == FALSE] <- 0
#date_is_DMY[104070:104082]	

#   Create correctly formatted date columns for each format	
date1a <- strftime(mdy(dat2$date), format="%Y-%m-%d")
date1b <- strftime(mdy_hm(dat2$date) , format="%Y-%m-%d")
date1 <- date1a
date1[is.na(date1a)] <- date1b[is.na(date1a)] 
#date1[140:150]
date1c <- strftime(ymd_hms(dat2$date) , format="%Y-%m-%d")
date1[is.na(date1)] <- date1c[is.na(date1)] 

date2a <- strftime(dmy(dat2$date), format="%Y-%m-%d")
date2b <- strftime(dmy_hm(dat2$date) , format="%Y-%m-%d")
date2 <- date2a
date2[is.na(date2a)] <- date2b[is.na(date2a)] 

#  Combine correced datae formats and select proper one
dat3 <- dat1 %>%
	cbind(date_is_DMY, date1, date2) %>%
	mutate(new_date = ifelse(date_is_DMY == "0", date1, date2))
dat3$new_date <- ymd(dat3$new_date)

# Checks
#dat3[140:150,c(4, 10:13)]
#dat3[104072:104102,c(4, 10:13)]



#  Clean up time formatting
#  Some rows have NA for time and time is embedded in date
na_times <- which(is.na(dat3$time))
pull_times <- data.frame(tmp1 = dat1$date[na_times]) %>%
	separate(tmp1, c("tmp1", "time"), sep = " ")
	
dat4 <- dat3
for(i in 1:length(na_times)){
	x <- na_times[i]
	dat4$time[x] <- pull_times$time[i]
	}


time1 <- parse_date_time(dat4$time, 
	c("mdy HM", "mdy HMS", "HM", "HMS", "ymd HMS", "ymd HM"))
time2 <- parse_date_time(dat4$time, 
	c("%I:%M:%S %p", "%I:%M %p"))
time3 <- time2
time3[is.na(time3)] <- time1[is.na(time3)]


dat5 <- dat4 %>%	
	select(fid, gid, name, new_date, time, latitude, longitude, elevation, filename) %>%
	cbind(time1, time2, time3) %>%
	separate(time3, c("tmp1", "new_time"), sep = " ") %>%
	mutate(date_time = paste(new_date, new_time, sep = " "))

#dat5[1:20, c(5, 10:13)] 
#dat5[109627:109633, c(5, 10:13)] 
#which(is.na(dat5$new_time))
#which(is.na(dat5$new_date))
#which(is.na(dat5$date_time))



#  Clean up white spaces
dat <- dat5 %>%
	select(fid, gid, name, new_date, new_time, date_time, latitude, longitude, elevation, filename) %>%
	mutate(name = tolower(str_replace_all(name, fixed(" "), ""))) %>%
	mutate(name_year = paste(name, year(date_time), sep = "_"))

dat$date_time_frmt <- ymd_hms(dat$date_time)

dat <- dat %>%
	arrange(name) %>%
	group_by(name) %>%
	arrange(date_time_frmt) %>%
	as.data.frame()


################################################
#  Save	
#save(dat, file = "dat.Rdata")


# Summary per individual
ind_dat <- dat %>%
	group_by(name) %>%
	summarize(
		n_locs = n(),
		earliest_date = min(date_time_frmt),
		last_date = max(date_time_frmt)) %>%
	as.data.frame()


##################################################################################################
##################################################################################################