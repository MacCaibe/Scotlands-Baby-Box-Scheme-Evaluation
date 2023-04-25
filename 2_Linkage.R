library(tidyverse)
library(lubridate)
library(scales)
library(UpSetR)
library(naniar)
library(labelled)

smr02 <- read_csv("smr02_data_cc_pub.csv")

# SMR01c ####

smr01c_data <- read_csv("1920-0125_smr01_children.csv")

smr01c_data <- smr01c_data %>%
  select(c(MOTHERID, WEEKS_BEFORE_AFTER_BABYBOX,
           COHORT_IDENTIFIER, EPISODES_WITHIN_26_WEEKS,
           EPISODES_WITHIN_52_WEEKS)) %>%
  rename(mother.id = MOTHERID,
         weeks.ba.bb = WEEKS_BEFORE_AFTER_BABYBOX, 
         cohort.id = COHORT_IDENTIFIER,
         episodes.26w.child = EPISODES_WITHIN_26_WEEKS,
         episodes.52w.child = EPISODES_WITHIN_52_WEEKS) 
smr01c_data$mother.id = str_sub(smr01c_data$mother.id, -6, -1)
sum(is.na(smr01c_data$weeks.ba.bb))

#keep all distinct (across all)
smr01c_data <- distinct(smr01c_data, .keep_all = FALSE)

#removing duplicates across ID variables
smr01c_data <- smr01c_data %>%
  group_by(mother.id, weeks.ba.bb, cohort.id) %>%
  filter(n()==1)

#join with smr02 
smr02 <- smr02 %>% 
  left_join(smr01c_data, by = NULL)

#replace NA with 0
smr02$episodes.26w.child[is.na(smr02$episodes.26w.child)] <- 0
smr02$episodes.52w.child[is.na(smr02$episodes.52w.child)] <- 0

# SMR01m ----

smr01m_data <- read_csv("1920-0125_smr01_mother.csv")

smr01m_data <- smr01m_data %>%
  select(c(MOTHERID, WEEKS_BEFORE_AFTER_BABYBOX,
           COHORT_IDENTIFIER, EPISODES_WITHIN_26_WEEKS,
           EPISODES_WITHIN_52_WEEKS)) %>%
  rename(mother.id = MOTHERID,
         weeks.ba.bb = WEEKS_BEFORE_AFTER_BABYBOX, 
         cohort.id = COHORT_IDENTIFIER,
         episodes.26w.mother = EPISODES_WITHIN_26_WEEKS,
         episodes.52w.mother = EPISODES_WITHIN_52_WEEKS) 
smr01m_data$mother.id = str_sub(smr01m_data$mother.id, -6, -1)
sum(is.na(smr01c_data$weeks.ba.bb))


#keep all distinct (across all)
smr01m_data <- distinct(smr01m_data, .keep_all = FALSE)

# Removing all duplicates across ID variables
smr01m_data <- smr01m_data %>%
  group_by(mother.id, weeks.ba.bb, cohort.id) %>%
  filter(n()==1)

#join 
smr02 <- smr02 %>% 
  left_join(smr01m_data, by = NULL)

#convert NA to 0 
smr02$episodes.26w.mother[is.na(smr02$episodes.26w.mother)] <- 0
smr02$episodes.52w.mother[is.na(smr02$episodes.52w.mother)] <- 0

# CHSP-PS fv ####
CHSP_fv_data <- read.csv("1920-0125_child_health_fv.csv")

#rename variables 
CHSP_fv_data <- CHSP_fv_data %>% 
  select(c(MOTHERID, WEEKS_BEFORE_AFTER_BABYBOX, COHORT_IDENTIFIER, 
           smoker_pc, second_hand_smoke,
           feed_2016)) %>%
  rename(mother.id = MOTHERID, 
         weeks.ba.bb = WEEKS_BEFORE_AFTER_BABYBOX, 
         cohort.id = COHORT_IDENTIFIER,
         smoke.pc.fv = smoker_pc, 
         smoke.sh.fv = second_hand_smoke, 
         feed.fv = feed_2016)
CHSP_fv_data$mother.id = str_sub(CHSP_fv_data$mother.id, -6, -1)

#keep all distinct (across all)
CHSP_fv_data <- distinct(CHSP_fv_data, .keep_all = FALSE)

# Removing all duplicates across ID variables
CHSP_fv_data <- CHSP_fv_data %>%
  group_by(mother.id, weeks.ba.bb, cohort.id) %>%
  filter(n()==1)

###input NAs

#smoke.pc.fv
table(CHSP_fv_data$smoke.pc.fv, useNA = "always")
CHSP_fv_data$smoke.pc.fv = na_if(CHSP_fv_data$smoke.pc.fv, "Unknown/Invalid")

#smoke.sh.fv
table(CHSP_fv_data$smoke.sh.fv, useNA = "always")
CHSP_fv_data$smoke.sh.fv = na_if(CHSP_fv_data$smoke.sh.fv, "Unknown/Invalid")

#feed.fv
table(CHSP_fv_data$feed.fv, useNA = "always")
CHSP_fv_data$feed.fv = na_if(CHSP_fv_data$feed.fv, "Unknown/Invalid")


#join 
smr02 <- smr02 %>% 
  left_join(CHSP_fv_data, by = NULL)



# CHSP-PS 6-8 ####
CHSP_6_8_data <- read_csv("1920-0125_child_health_6_8.csv")

#rename variables
CHSP_6_8_data <- CHSP_6_8_data %>% 
  select(c(MOTHERID, WEEKS_BEFORE_AFTER_BABYBOX, COHORT_IDENTIFIER, 
           smoker_pc, second_hand_smoke, feed_2016, sleep_pr,
           sleep_su, sleep_si)) %>%
  rename(mother.id = MOTHERID, 
         weeks.ba.bb = WEEKS_BEFORE_AFTER_BABYBOX, 
         cohort.id = COHORT_IDENTIFIER, 
         smoke.pc.68 = smoker_pc, 
         smoke.sh.68 = second_hand_smoke, 
         feed.68 = feed_2016,
         sleep.pr.68 = sleep_pr,
         sleep.su.68 = sleep_su,
         sleep.si.68 = sleep_si)

CHSP_6_8_data$mother.id = str_sub(CHSP_6_8_data$mother.id, -6, -1)

# keep all distinct
CHSP_6_8_data <- distinct(CHSP_6_8_data, .keep_all = FALSE)
# Removing all duplicates across ID variables
CHSP_6_8_data <- CHSP_6_8_data %>%
  group_by(mother.id, weeks.ba.bb, cohort.id) %>%
  filter(n()==1)

#smoke.pc
table(CHSP_6_8_data$smoke.pc.68, useNA = "always")
CHSP_6_8_data$smoke.pc.68 = na_if(CHSP_6_8_data$smoke.pc.68, "Unknown/Invalid")

#smoke.sh
table(CHSP_6_8_data$smoke.sh.68, useNA = "always")
CHSP_6_8_data$smoke.sh.68 = na_if(CHSP_6_8_data$smoke.sh.68, "Unknown/Invalid")

#feed
table(CHSP_6_8_data$feed.68, useNA = "always")
CHSP_6_8_data$feed.68 = na_if(CHSP_6_8_data$feed.68, "Unknown/Invalid")

#sleep
table(CHSP_6_8_data$sleep.pr.68, useNA = "always")
table(CHSP_6_8_data$sleep.su.68, useNA = "always")
table(CHSP_6_8_data$sleep.si.68, useNA = "always")


CHSP_6_8_data$sleep.pr.68[is.na(CHSP_6_8_data$sleep.pr.68)] <- 0
CHSP_6_8_data$sleep.su.68[is.na(CHSP_6_8_data$sleep.su.68)] <- 0
CHSP_6_8_data$sleep.si.68[is.na(CHSP_6_8_data$sleep.si.68)] <- 0

CHSP_6_8_data$sleep.pr.68[CHSP_6_8_data$sleep.pr.68 == "N"] <- 0
CHSP_6_8_data$sleep.su.68[CHSP_6_8_data$sleep.su.68 == "N"] <- 0
CHSP_6_8_data$sleep.si.68[CHSP_6_8_data$sleep.si.68 == "N"] <- 0

CHSP_6_8_data$sleep.pr.68[CHSP_6_8_data$sleep.pr.68 == "Y"] <- 1
CHSP_6_8_data$sleep.su.68[CHSP_6_8_data$sleep.su.68 == "Y"] <- 1
CHSP_6_8_data$sleep.si.68[CHSP_6_8_data$sleep.si.68 == "Y"] <- 1

CHSP_6_8_data$sleep.pr.68 <- as.numeric(CHSP_6_8_data$sleep.pr.68)
CHSP_6_8_data$sleep.su.68 <- as.numeric(CHSP_6_8_data$sleep.su.68)
CHSP_6_8_data$sleep.si.68 <- as.numeric(CHSP_6_8_data$sleep.si.68)

CHSP_6_8_data <- CHSP_6_8_data %>% 
  mutate(sleep.count = sleep.pr.68 + sleep.su.68 + sleep.si.68)

#create binary
CHSP_6_8_data <- CHSP_6_8_data %>%
  mutate(sleep.bin.su = if_else(sleep.count == 1 & sleep.su.68 == 1, 1, 0))

#code missing values (i.e., missing for all three sleep variables)
CHSP_6_8_data <- CHSP_6_8_data %>%
  mutate(sleep.bin.su = replace(sleep.bin.su, sleep.count==0, NA))

table(CHSP_6_8_data$sleep.bin.su, useNA = "always")

#join 
CHSP_6_8_data <- CHSP_6_8_data %>%
  select(c(mother.id, weeks.ba.bb, cohort.id, smoke.sh.68,
           smoke.pc.68, feed.68, sleep.bin.su))

smr02 <- smr02 %>% 
  left_join(CHSP_6_8_data, by = NULL)

# write linked SMR02 to file ----
write_excel_csv(smr02, 
                col_names = TRUE,
                delim = ",",
                path = "smr02_linked_pub.csv")
