library(tidyverse)
library(lubridate)
library(scales)
library(UpSetR)
library(naniar)
library(labelled)
#### Step 1 Import data & define study timeframe ####

smr02_data <- read_csv("1920-0125_smr02.csv")

#first need to order observations sequentially using weeks and cohort ID variables
smr02_data <- smr02_data %>%
  mutate(time.weeks = ifelse(COHORT_IDENTIFIER == "UNEXPOSED",
                             (WEEKS_BEFORE_AFTER_BABYBOX - 106) * -1,
                             WEEKS_BEFORE_AFTER_BABYBOX + 105))

#plot of observations over time
smr02_data %>% 
  group_by(time.weeks) %>%
  mutate(initial.cases = n()) %>%
  summarise(initial.cases = mean(initial.cases)) %>%
  ggplot(aes(x=time.weeks, y=initial.cases)) + geom_line()

#low observations when time.weeks = 1 & 209
#because these timepoints do not provide a full week of data
#spike due to week of SBBS introduction being coded as time.weeks = 104, as weeks before after = 0. 
#drop weeks 1, 209, and week of introduction

smr02_data <- smr02_data %>%
  filter(WEEKS_BEFORE_AFTER_BABYBOX > 0, 
         time.weeks >1, 
         time.weeks <209)

#need to re-write time.weeks to match loss of weeks 
smr02_data <- smr02_data %>%
  mutate(time.weeks = time.weeks-1)

#plot of observations again
smr02_data %>% 
  group_by(time.weeks) %>%
  mutate(initial.cases = n()) %>%
  summarise(initial.cases = mean(initial.cases)) %>%
  ggplot(aes(x=time.weeks, y=initial.cases)) + geom_line() 


#### STEP2 Initial transformation ####

#rename variables & apply exclusion criteria
smr02_data %>%
  filter(num_of_births_this_pregnancy > 1) %>%
  nrow()


smr02_data %>%
  filter(outcome_of_pregnancy_baby_1 > 1) %>%
  nrow()


smr02_data <- smr02_data %>% 
  select(c(MOTHERID, WEEKS_BEFORE_AFTER_BABYBOX, COHORT_IDENTIFIER, 
           num_of_births_this_pregnancy, 
           outcome_of_pregnancy_baby_1,
           sex_baby_1,
           condition_on_discharge,
           APGAR_5_MINUTES_BABY_1, 
           BIRTHWEIGHT_BABY_1, 
           FEED_ON_DISCHARGE_BABY1, 
           MARITAL_STATUS, 
           TOTAL_PREVIOUS_PREGNANCIES, 
           BOOKING_SMOKING_HISTORY, 
           SMOKER_DURING_PREGNANCY, 
           MODE_OF_DELIVERY_BABY_1,  
           ETHNIC_GROUP,  
           MOTHER_AGE_GRPS, 
           simd2016_sc_quintile,
           time.weeks)) %>%
  filter(num_of_births_this_pregnancy == 1, outcome_of_pregnancy_baby_1 == 1, condition_on_discharge == 3) %>%
  rename(mother.id = MOTHERID, 
         weeks.ba.bb = WEEKS_BEFORE_AFTER_BABYBOX, 
         cohort.id = COHORT_IDENTIFIER, 
         number.births = num_of_births_this_pregnancy,
         pregnancy.outcome = outcome_of_pregnancy_baby_1,
         sex = sex_baby_1,
         discharge.condition = condition_on_discharge,
         apgar.5min = APGAR_5_MINUTES_BABY_1, 
         birthweight = BIRTHWEIGHT_BABY_1, 
         feed.on.discharge = FEED_ON_DISCHARGE_BABY1, 
         marital.status = MARITAL_STATUS, 
         total.previous.pregnancies = TOTAL_PREVIOUS_PREGNANCIES, 
         book.smoking.history = BOOKING_SMOKING_HISTORY, 
         smoking.during.pregnancy = SMOKER_DURING_PREGNANCY, 
         delivery.mode = MODE_OF_DELIVERY_BABY_1, 
         ethnic.grp = ETHNIC_GROUP, 
         mother.age = MOTHER_AGE_GRPS, 
         simd.quintile.2016 = simd2016_sc_quintile)


# check for distinct observations across all variables
nrow(distinct(smr02_data))
# drop duplicates
smr02_data <- distinct(smr02_data, .keep_all = FALSE)

# check for distinct mother.id, weeks.ba.bb, & cohort.id
nrow(distinct(smr02_data, mother.id, weeks.ba.bb, cohort.id))
# it is not possible to know if one is authentic over the other
# thus ALL matching cases of this kind are dopped
smr02_data <- smr02_data %>%
  group_by(mother.id, weeks.ba.bb, cohort.id) %>%
  filter(n()==1)


# check
nrow(distinct(smr02_data, mother.id, weeks.ba.bb, cohort.id))
nrow(distinct(smr02_data))

#All observations are unique live singleton births


#---- STEP3 SMR02 check variables & missingness ####
# Here I a) summarise each variable,
# b) convert values indicating missingness to NAs,
# c) alter variables in preparation for analysis.

#mother.id
smr02_data$mother.id = str_sub(smr02_data$mother.id, -6, -1)#removes unwanted character string
sum(is.na(smr02_data$mother.id))

#week.ba.bb
sum(is.na(smr02_data$weeks.ba.bb))
nrow(filter(smr02_data, weeks.ba.bb ==0))

#cohort.id
sum(is.na(smr02_data$cohort.id))

#number.births
sum(is.na(smr02_data$number.births))
hist(smr02_data$number.births) #all singleton births (i.e. =1)

#pregnancy.outcome
sum(is.na(smr02_data$pregnancy.outcome))
hist(smr02_data$pregnancy.outcome)#all live briths (i.e. =1)

#feed.on.discharge
smr02_data$feed.on.discharge = na_if(smr02_data$feed.on.discharge, 9) #recoding 'not known' to NAs
smr02_data$feed.on.discharge = na_if(smr02_data$feed.on.discharge, 3)
table(smr02_data$feed.on.discharge, useNA = "always")

#total.previous.pregnancies
smr02_data$total.previous.pregnancies = na_if(smr02_data$total.previous.pregnancies, 99)
summary(smr02_data$total.previous.pregnancies)
smr02_data <- smr02_data %>%
  replace_with_na_at(.vars = "total.previous.pregnancies", condition = ~.x > 15)

#book.smoking.history
smr02_data$book.smoking.history = na_if(smr02_data$book.smoking.history, 9)
table(smr02_data$book.smoking.history, useNA = "always")

#smoking.during.pregnancy
smr02_data$smoking.during.pregnancy = na_if(smr02_data$smoking.during.pregnancy, 9)
table(smr02_data$smoking.during.pregnancy, useNA = "always")

#mother.age
table(smr02_data$mother.age, useNA = "always")

#simd.quintile.2016
table(smr02_data$simd.quintile.2016, useNA = "always")

#sex
smr02_data$sex = na_if(smr02_data$sex, 9) #recoding 'not specified' to NAs
smr02_data$sex = na_if(smr02_data$sex, 0) #recoding 'not known' to NAs
table(smr02_data$sex, useNA = "always")

#can now check missingness across time

#missingeness between variables
smr02_data %>%
  select(c(time.weeks, feed.on.discharge, total.previous.pregnancies,
           smoking.during.pregnancy, book.smoking.history, simd.quintile.2016)) %>%
  gg_miss_fct(., time.weeks)

smr02_data %>%
  select(c(time.weeks, feed.on.discharge, total.previous.pregnancies,
           smoking.during.pregnancy, book.smoking.history, simd.quintile.2016)) %>%
  gg_miss_upset(nsets = 10, nintersects = 50)


#---- STEP4 SMR02 complete case > write to file ####
#Creat complete case dataset
smr02_data_cc <- smr02_data %>%
  select(c(mother.id, 
           weeks.ba.bb,
           cohort.id,
           feed.on.discharge,
           total.previous.pregnancies,
           book.smoking.history,
           smoking.during.pregnancy,
           mother.age,
           simd.quintile.2016,
           time.weeks, 
           sex)) 

smr02_data_cc <- smr02_data_cc %>%
  drop_na() 

#create level variable 
smr02_data_cc <- smr02_data_cc %>% 
  mutate(level = 
           case_when(between(time.weeks, 1, 104) ~ 0,
                     between(time.weeks, 105, 207) ~ 1))

#create trend variables
smr02_data_cc <- smr02_data_cc %>%
  mutate(trend.weeks =
           case_when(between(time.weeks, 1, 104) ~ 0,
                     between(time.weeks, 105, 105) ~ 1,
                     between(time.weeks, 106, 106) ~ 2,
                     between(time.weeks, 107, 107) ~ 3,
                     between(time.weeks, 108, 108) ~ 4,
                     between(time.weeks, 109, 109) ~ 5,
                     between(time.weeks, 110, 110) ~ 6, 
                     between(time.weeks, 111, 111) ~ 7,
                     between(time.weeks, 112, 112) ~ 8,
                     between(time.weeks, 113, 113) ~ 9, 
                     between(time.weeks, 114, 114) ~ 10,
                     between(time.weeks, 115, 115) ~ 11, 
                     between(time.weeks, 116, 116) ~ 12,
                     between(time.weeks, 117, 117) ~ 13,
                     between(time.weeks, 118, 118) ~ 14,
                     between(time.weeks, 119, 119) ~ 15,
                     between(time.weeks, 120, 120) ~ 16, 
                     between(time.weeks, 121, 121) ~ 17,
                     between(time.weeks, 122, 122) ~ 18, 
                     between(time.weeks, 123, 123) ~ 19,
                     between(time.weeks, 124, 124) ~ 20,
                     between(time.weeks, 125, 125) ~ 21,
                     between(time.weeks, 126, 126) ~ 22, 
                     between(time.weeks, 127, 127) ~ 23,
                     between(time.weeks, 128, 128) ~ 24, 
                     between(time.weeks, 129, 129) ~ 25, 
                     between(time.weeks, 130, 130) ~ 26, 
                     between(time.weeks, 131, 131) ~ 27,
                     between(time.weeks, 132, 132) ~ 28, 
                     between(time.weeks, 133, 133) ~ 29, 
                     between(time.weeks, 134, 134) ~ 30,
                     between(time.weeks, 135, 135) ~ 31,
                     between(time.weeks, 136, 136) ~ 32, 
                     between(time.weeks, 137, 137) ~ 33,
                     between(time.weeks, 138, 138) ~ 34, 
                     between(time.weeks, 139, 139) ~ 35, 
                     between(time.weeks, 140, 140) ~ 36,
                     between(time.weeks, 141, 141) ~ 37,
                     between(time.weeks, 142, 142) ~ 38, 
                     between(time.weeks, 143, 143) ~ 39, 
                     between(time.weeks, 144, 144) ~ 40, 
                     between(time.weeks, 145, 145) ~ 41, 
                     between(time.weeks, 146, 146) ~ 42, 
                     between(time.weeks, 147, 147) ~ 43, 
                     between(time.weeks, 148, 148) ~ 44, 
                     between(time.weeks, 149, 149) ~ 45,
                     between(time.weeks, 150, 150) ~ 46, 
                     between(time.weeks, 151, 151) ~ 47, 
                     between(time.weeks, 152, 152) ~ 48,
                     between(time.weeks, 153, 153) ~ 49,
                     between(time.weeks, 154, 154) ~ 50, 
                     between(time.weeks, 155, 155) ~ 51,
                     between(time.weeks, 156, 156) ~ 52,
                     between(time.weeks, 157, 157) ~ 53,
                     between(time.weeks, 158, 158) ~ 54,
                     between(time.weeks, 159, 159) ~ 55,
                     between(time.weeks, 160, 160) ~ 56,
                     between(time.weeks, 161, 161) ~ 57, 
                     between(time.weeks, 162, 162) ~ 58,
                     between(time.weeks, 163, 163) ~ 59,
                     between(time.weeks, 164, 164) ~ 60, 
                     between(time.weeks, 165, 165) ~ 61,
                     between(time.weeks, 166, 166) ~ 62, 
                     between(time.weeks, 167, 167) ~ 63,
                     between(time.weeks, 168, 168) ~ 64,
                     between(time.weeks, 169, 169) ~ 65,
                     between(time.weeks, 170, 170) ~ 66,
                     between(time.weeks, 171, 171) ~ 67, 
                     between(time.weeks, 172, 172) ~ 68,
                     between(time.weeks, 173, 173) ~ 69, 
                     between(time.weeks, 174, 174) ~ 70,
                     between(time.weeks, 175, 175) ~ 71,
                     between(time.weeks, 176, 176) ~ 72,
                     between(time.weeks, 177, 177) ~ 73, 
                     between(time.weeks, 178, 178) ~ 74,
                     between(time.weeks, 179, 179) ~ 75, 
                     between(time.weeks, 180, 180) ~ 76, 
                     between(time.weeks, 181, 181) ~ 77, 
                     between(time.weeks, 182, 182) ~ 78,
                     between(time.weeks, 183, 183) ~ 79, 
                     between(time.weeks, 184, 184) ~ 80, 
                     between(time.weeks, 185, 185) ~ 81,
                     between(time.weeks, 186, 186) ~ 82,
                     between(time.weeks, 187, 187) ~ 83, 
                     between(time.weeks, 188, 188) ~ 84,
                     between(time.weeks, 189, 189) ~ 85, 
                     between(time.weeks, 190, 190) ~ 86, 
                     between(time.weeks, 191, 191) ~ 87,
                     between(time.weeks, 192, 192) ~ 88,
                     between(time.weeks, 193, 193) ~ 89, 
                     between(time.weeks, 194, 194) ~ 90, 
                     between(time.weeks, 195, 195) ~ 91, 
                     between(time.weeks, 196, 196) ~ 92, 
                     between(time.weeks, 197, 197) ~ 93, 
                     between(time.weeks, 198, 198) ~ 94, 
                     between(time.weeks, 199, 199) ~ 95, 
                     between(time.weeks, 200, 200) ~ 96,
                     between(time.weeks, 201, 201) ~ 97, 
                     between(time.weeks, 202, 202) ~ 98, 
                     between(time.weeks, 203, 203) ~ 99,
                     between(time.weeks, 204, 204) ~ 100,
                     between(time.weeks, 205, 205) ~ 101, 
                     between(time.weeks, 206, 206) ~ 102,
                     between(time.weeks, 207, 207) ~ 103))



#Write the cc dataframe to file to make it easier to work with
write_excel_csv(smr02_data_cc, 
                col_names = TRUE,
                delim = ",",
                path = "smr02_data_cc_pub.csv"
)

