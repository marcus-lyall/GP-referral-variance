
#load library and set plot theme perameters.

library(lubridate)
library(ggplot2)
library(dplyr)
library(tidyverse)
library(zoo)
library(stringr)
library(tidyr)
library(data.table)
library(knitr)
library(broom)
library(ggpubr)
library(MASS)
library(finalfit)



mytheme = theme(axis.line = element_line(colour = "black"),
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                panel.background = element_blank(),
                axis.text.x = element_text(size = 6,face = "bold", angle = 90),
                axis.text.y = element_text(size = 8,face = "bold"),
                axis.title = element_text(size = 8, face = "bold"),
                legend.title = element_text(size = 8, colour = "black"),
                legend.text=element_text(size=8), 
                legend.key.size = unit(0.6, "cm"))



#Load in data - this assumes the data is in wie format with one patient referral per row and associated variables such as referring practice, patient age, gender, Charlson score, AKIN score, SIMD quintile in columns. 

ALL = read.csv(file = "yourdata.csv", 
               header = TRUE)

#convert character vectors of dates to POSIXCT format and relevant variables to factors. 

ALL$ARRIVAL_DATE = as.POSIXct(ALL$ARRIVAL_DATE,format = "%Y-%m-%d %H:%M:%S", tmz = "GMT" )
ALL$IP_WARD_ADMITTED_TO = as.factor(ALL$IP_WARD_ADMITTED_TO)
ALL$GP_PRACTICE_CODE = as.factor(ALL$GP_PRACTICE_CODE)
ALL$SOURCE_OF_REFERRAL = as.factor(ALL$SOURCE_OF_REFERRAL)

# select out timeframes and factor levels to ensure the correct population is examined as follows. 

ALL = subset(ALL, format(ARRIVAL_DATE,"%H:%M") > "08:30")
ALL = subset(ALL, format(ARRIVAL_DATE,"%H:%M") < "22:00")

ALL = ALL %>%
  filter(SOURCE_OF_REFERRAL %in% c("General Practitioner", "Flow Centre")) %>%
  filter(LOCATION_CODE %in% c("WGH", "RIE"))


ALL = ALL %>% 
  filter(LOCATION_CODE %in% c("WGH")) %>% 
  filter(SOURCE_OF_REFERRAL %in% c("General Practitioner", "Flow Centre")) %>% 
  filter(!FIRST_TRIAGE_CATEGORY %in% c("5 - Non Urgent")) %>%
  filter(ADMISSION_WARD %in% c("WGH Red Zone")) %>%
  filter(!DISCHARGE_WARD_CODE %in% c("WGHSAU", "WGHMINOR"))


# Collapse data to practice level


GPrefs = ALL %>% 
  group_by(GP_PRACTICE_CODE, admitted) %>%
  summarise(n = n()) %>% 
  mutate(prop = prop.table(n)) %>%
  filter(admitted == TRUE) 

# add total numbers again
GPrefs$Referrals = GPrefs$n/GPrefs$prop


# add practice level variables to data frame in the following way (you'll nee to to this for all practice level variables) such as %patient over 65, % in bottom 2 SIMD quintiles etc. 
GPrefs$practice_list_size = GP_data$Practice.List.Size[match(GPrefs$GP_PRACTICE_CODE, GP_data$Practice.Code)]
GPrefs$main_hospital = GPzones$Hospital[match(GPrefs$GP_PRACTICE_CODE, GPzones$Practice.Code)]


# calculate referrals per 100 patient years
GPrefs$patient_years = GPrefs$practice_list_size*(33/12)
GPrefs$ref_rate = GPrefs$Referrals/(GPrefs$patient_years/100)

#  add distances from hospital

library(geosphere)

# postcodes.csv is a data table of uk postcodes(zipcodes) with associated longditude nad latitude cooordinates - there are many web resources such as this https://www.doogal.co.uk/PostcodeDownloads.php

ukpost = read.csv("postcodes.csv", header = TRUE)

# get the coordinates for each major hospital
RIEcoord = ukpost[ukpost$Postcode %in% c("EH16 4SA"),c(4,5)]
WGHcoord = ukpost[ukpost$Postcode %in% c("EH4 2XU"),c(4,5)]

# make columns to populate Long and lat. 
GPrefs$Hospital_Lat = NA
GPrefs$Hospital_long = NA

GPrefs$Hospital_Lat = ifelse(GPrefs$main_hospital %in% c("RIE"), RIEcoord$Latitude, GPrefs$Hospital_Lat)
GPrefs$Hospital_Lat = ifelse(GPrefs$main_hospital %in% c("WGH"), WGHcoord$Latitude, GPrefs$Hospital_Lat)

GPrefs$Hospital_long = ifelse(GPrefs$main_hospital %in% c("RIE"), RIEcoord$Longitude, GPrefs$Hospital_long)
GPrefs$Hospital_long = ifelse(GPrefs$main_hospital %in% c("WGH"), WGHcoord$Longitude, GPrefs$Hospital_long)


# add longditude and lattitude for each GP practice based on postcode. 
GPrefs$GPlat = ukpost$Latitude[match(GPrefs$practice_postcode, ukpost$Postcode)]
GPrefs$GPlong = ukpost$Longitude[match(GPrefs$practice_postcode, ukpost$Postcode)]

# make data frames for disthaversine function.  
Hospital_coord = as.matrix(data.frame(longditude = GPrefs$Hospital_long , latitude = GPrefs$Hospital_Lat))
GP_coord = as.matrix(data.frame(longditude = GPrefs$GPlong, latitude = GPrefs$GPlat))

#get the distances to RIE and add to the GP referrals matri. 
GPrefs$distance_2_hospital = (distHaversine(Hospital_coord, GP_coord, r=6378137))/1000



# construct model and write results


model.1q <- glm(referrals ~  pop_in_bottom_two_quartiles + perc_over65b + percentageover65incare + distance_2_hospital + practice_list_size, offset = log(patient_years), family = quasipoisson(link = "log") , data = GPrefs)


modelgpref = tidy(model.1q)
write.csv(modelgpref, file = "modelgpref.csv")



# using the model get the predicted referral rate for each practice 
GPrefs$predicted_referrals <- exp(predict(model.1q))

# then calculate predicted referral number and adjusted referral rate. 
GPrefs$expected_referral_rate = GPrefs$predicted_referrals*(100/GPrefs$patient_years)

GPrefs$observed_over_expected_rate = GPrefs$Referrals/GPrefs$predicted_referrals

total_mean_rate = mean(GPrefs$ref_rate)

GPrefs$adjusted_referral_rate = GPrefs$observed_over_expected_rate*total_mean_rate

# use ggpubr package to general expercted versus observed plots and export. 

GPrefs$pred.ref.resid = GPrefs$referrals - GPrefs$predicted_referral

p = ggscatter(GPrefs,  x = "pred.refs", y = "referrals", 
              add = "reg.line",  # Add regressin line
              add.params = list(color = "red", linetype = "dashed"))

p = p +   stat_cor(
  aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")),
  label.x = 1
)
p = p + mytheme
p = p + xlab("predicted referrals")
p = p + ylab("observed referrals")

r = ggplot(GPrefs)
r = r + geom_point(aes(y = pred.ref.resid, x = referrals))
#you can general your own aesthetics http://ggplot2.tidyverse.org/reference/theme.html
r = r + theme(axis.text = element_text(angle = 90, hjust = 1))
r = r + theme(axis.line = element_line(colour = "black"),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.background = element_blank(),
              axis.text.x = element_text(size = 6,face = "bold", angle = 90),
              axis.text.y = element_text(size = 8,face = "bold"),
              axis.title = element_text(size = 8, face = "bold"),
              legend.title = element_text(size = 8, colour = "black"),
              legend.text=element_text(size=8), 
              legend.key.size = unit(0.6, "cm"))
r = r + ylab("observed - predicted referrals")
r = r + xlab("observed referrals")
r = r + geom_hline(yintercept = 0, color = "red", linetype = "dashed")


jpeg("gp_refs_model-scatter.tiff", width = 4, height = 4, units = 'in', res = 300)

p

dev.off()

tiff("gprefs-resid.tiff", width = 4, height = 4, units = 'in', res = 300)

r

dev.off()

grid.arrange(p, r, ncol = 2)


#  Use scripts such as the following to generate quartiles and deciles and fold changes 

GPrefs$quartile <- ntile(GPrefs$ref_rate, 4)
GPrefs$decile = ntile(GPrefs$referral_rate, 10)


deciles = GPrefs %>% 
  group_by(decile) %>% 
  summarise(referral_rate = mean(referal_rate))



# using calculted fold changes create scatter plots with horizontal lines and export

GPrefs_var = GPrefs

GPrefs_var_long = as.data.frame(GPrefs_var)
GPrefs_var_long = GPrefs_var_long[,c("GP_PRACTICE_CODE", "ref_rate", "adjusted_referral_rate")]

colnames(GPrefs_var_long)=c("GP_PRACTICE_CODE", "crude_referral_rate", "adjusted_referral_rate")

GPrefs_var_long = gather(GPrefs_var_long, key, value, -GP_PRACTICE_CODE)

GPrefs_var_long$key = factor(GPrefs_var_long$key, levels = c("crude_referral_rate", "adjusted_referral_rate"))


sumtab = GPrefs_var_long %>%
  group_by(key) %>%
  summarise(median = median(value, na.rm = TRUE))
sumtab$midfirst = c(0.8,0.9)
sumtab$midtenth = c(4.3,3.6)
sumtab$midfirstquart = c(1.29,1.43)
sumtab$midfourthquart = c(3.27, 3.08)

supp.labs <- c("crude referral rate", "adjusted referral rate")
names(supp.labs) <- c("crude_referral_rate", "adjusted_referral_rate")

b = ggplot(GPrefs_var_long) 

b = b + geom_jitter(aes(x = key, y = value), width = 0.1, size = 2.0, pch = 21)

b = b + mytheme

b = b + ylab("referrals per 100 patient years")

b = b + geom_hline(data = sumtab, aes(x = key, yintercept = median), width = 0.2)

b = b + geom_hline(data = sumtab, aes(x = key, yintercept = midfirst), linetype = "dashed", width = 0.2)
b = b + geom_hline(data = sumtab, aes(x = key, yintercept = midtenth), linetype = "dashed", width = 0.2)

b = b + geom_hline(data = sumtab, aes(x = key, yintercept = midfirstquart), linetype = "dotted", width = 0.2)
b = b + geom_hline(data = sumtab, aes(x = key, yintercept = midfourthquart), linetype = "dotted", width = 0.2)



b = b + facet_wrap(~key, scales = "free_x", labeller = labeller(key = supp.labs))
b = b + xlab("")
b = b + theme(axis.title.x=element_blank(),
              axis.text.x=element_blank(),
              axis.ticks.x=element_blank())




b

tiff("variance.tiff", width = 3.5, height = 3, units = 'in', res = 300)
b # Make plot
dev.off()





#  For patient level analysis add practice referral quartile to the individualised data set. 


ALL = referred_cohort

referred_cohort$referral_quartile = GPrefs$quartile[match(referred_cohort$GP_PRACTICE_CODE, GPrefs$GP_PRACTICE_CODE)]

# name your referral quartiles in way that suites. Do this for AKIN score any other variables as needed. 

referred_cohort$practice_referral_quartile = ifelse(referred_cohort$referral_quartile %in% "3", "second referral quartile", referred_cohort$practice_referral_quartile )

referred_cohort$practice_referral_quartile = ifelse(referred_cohort$referral_quartile %in% "2", "third referral quartile", referred_cohort$practice_referral_quartile)

referred_cohort$practice_referral_quartile = ifelse(referred_cohort$referral_quartile %in% "1", "lowest referral quartile", referred_cohort$practice_referral_quartile)

referred_cohort$practice_referral_quartile = ifelse(referred_cohort$referral_quartile %in% "4", "highest referral quartile", referred_cohort$practice_referral_quartile)

# order your factors. 

referred_cohort$practice_referral_quartile_group = factor(referred_cohort$practice_referral_quartile, levels = c( "highest referral quartile", "third referral quartile", "second referral quartile", "lowest referral quartile"))


# Perform binary logistic regression modelling using finalfit package and export output. repeat with AMIN score included for sensitivity analysis.  

library(finalfit)
explanatory = c("SEX", "AGE_ON_ARRIVAL", "SIMD","distance_to_hospital", "care_home",  "charleson_factor", "practice_referral_quartile")
random_effect = "LOCATION_CODE"
dependent = "admitted"



referred_cohort %>%
  summary_factorlist(dependent, explanatory, fit_id=TRUE) -> t1

referred_cohort %>%
  glmmixed(dependent, explanatory, random_effect = random_effect)%>%
  fit2df(estimate_suffix=" (multilevel)") -> t2


t1 %>% 
  ff_merge(t2) %>% 
  dplyr::select(-c(fit_id, index)) %>% 
  dependent_label(referred_cohort, dependent)-> t3


write.csv(t3, file = "model_output.csv" )

# generate and output OR plot at required image resolution. 


x = glmmixed( referred_cohort, explanatory = explanatory, dependent = dependent, random_effect = random_effect)


tiff("OR_plot.tiff", width = 8, height = 5, units = 'in', res = 300)

referred_cohort%>% 
  or_plot(dependent, explanatory, random_effect = random_effect,  glmfit = x, table_text_size = 3,
          title_text_size = 14, column_space = c(-0.65, 0, 0.4)) 


dev.off()





