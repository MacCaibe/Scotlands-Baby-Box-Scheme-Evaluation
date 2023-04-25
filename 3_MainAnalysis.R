library(foreign)
library(splines)
library(lmtest)
library(Epi)
library(tsModel)
library(vcd)
library(tidyverse)
library(lubridate)
library(scales)
library(UpSetR)
library(naniar)
library(nlme)
library(gmodels)
library(sjPlot)
library(ciTools)


data_cc <- read_csv("smr02_linked_pub.csv")
trend_pred <- read_csv("trend_pred.csv")

data_cc <- data_cc %>%
  mutate(tpp_bin = if_else(total.previous.pregnancies > 0, 1,0))

#population characteristics

table(data_cc$mother.age, useNA = 'always')
table(data_cc$simd.quintile.2016)
table(data_cc$tpp_bin)

# SMR-01 - child 26w ----
smr01c_26 <- data_cc %>%
  select(time.weeks, level, trend.weeks, episodes.26w.child)

#create denominator (person-weeks)
smr01c_26 <- smr01c_26 %>%
  group_by(time.weeks) %>%
  mutate(cases.per.week.smr01 = n()) 
smr01c_26$cases.per.week.smr01 <- smr01c_26$cases.per.week.smr01*26

#aggregate
smr01c_26 <- smr01c_26 %>%
  group_by(time.weeks) %>%
  summarise(trend.weeks = mean(trend.weeks), 
            level = mean(level), 
            cases.per.week.smr01 = mean(cases.per.week.smr01), 
            episodes.26w.child = sum(episodes.26w.child))

plot(smr01c_26$cases.per.week.smr01, type ="n", ylim = c(10000,30000), xlab = "Week of birth", 
     ylab = "Observations", xaxt="n")
points(smr01c_26$cases.per.week.smr01, cex=0.7)
axis(1, at=0:52*4, labels = TRUE)
title("Infant admissions within 26 weeks")


# Incidence rate
# This gives the average number of episodes per 1000 persons per week at risk
smr01c_26$rate <- with(smr01c_26, episodes.26w.child/cases.per.week.smr01)*1000

#plot pre-intervention trend
plot(smr01c_26$rate, type ="n", ylim = c(0,20), xlab = "Week of birth", 
     ylab = "Rate x 1000", xaxt="n")
rect(105, 0, 208, 20, col=grey(0.9), border = F)
points(smr01c_26$rate, cex=0.7)
axis(1, at=0:52*4, labels = TRUE)
title("Infant admissions within 26 weeks")

summary(smr01c_26$episodes.26w.child[smr01c_26$level==0])
summary(smr01c_26$episodes.26w.child[smr01c_26$level==1])
summary(smr01c_26$rate)
summary(smr01c_26$rate[smr01c_26$level==0])
summary(smr01c_26$rate[smr01c_26$level==1])

### model 1 - unadjusted
model1 <- glm(episodes.26w.child ~ offset(log(cases.per.week.smr01))
              + time.weeks + level + trend.weeks, family = poisson, smr01c_26)
summary(model1)  
summary(model1)$dispersion
round(ci.exp(model1), 4)

# model2 - allow for over-dispersion
model2 <- glm(episodes.26w.child ~ offset(log(cases.per.week.smr01))
              + time.weeks + level + trend.weeks, family = quasipoisson, smr01c_26)
summary(model2)  
summary(model2)$dispersion  
round(ci.exp(model2), 3)

res2 <- residuals(model2, type = "deviance")
plot(smr01c_26$time.weeks, res2, 
     ylim = c(-10,10), pch=19, cex=0.7, col=grey(0.6),
     main = "Infant admissions within 26 weeks (model 2)", ylab = "Deviance residuals",
     xlab = "Week of birth")
abline(h=0, lty=2, lwd=2)

#check acf & p-acf
acf(res2,
    main = "Infant admissions within 26 weeks (model 2)")
pacf(res2,
     main = "Infant admissions within 26 weeks (model 2)")




#adjust for seasonality
model3 <- glm(episodes.26w.child ~ offset(log(cases.per.week.smr01))
              + time.weeks + level + trend.weeks + harmonic(time.weeks,2,52), 
              family = quasipoisson, smr01c_26)
summary(model3)
summary(model3)$dispersion  
round(ci.exp(model3), 4)

res3 <- residuals(model3, type = "deviance")
plot(smr01c_26$time.weeks, res3, 
     ylim = c(-10,10), pch=19, cex=0.7, col=grey(0.6),
     main = "Infant admissions within 26 weeks (model 3)", ylab = "Deviance residuals",
     xlab = "Week of birth")
abline(h=0, lty=2, lwd=2)

#check acf & p-acf
acf(res3,
    main =  "Infant admissions within 26 weeks (model 3)")
pacf(res3,
     main = "Infant admissions within 26 weeks (model 3)")

# model 3 predictions
datanew <- data.frame(cases.per.week.smr01=mean(smr01c_26$cases.per.week.smr01),
                      level=rep(c(0,1), c(104,103)),
                      time.weeks=1:207)
datanew <- cbind(datanew,trend_pred)
smr01c_26$factual=predict(model3, type="response", datanew)/mean(smr01c_26$cases.per.week.smr01)*1000

datanew_b <- data.frame(cases.per.week.smr01=mean(smr01c_26$cases.per.week.smr01),
                        level=0,
                        trend.weeks=0,
                        time.weeks=1:207)
smr01c_26$c.factual=predict(model3, type="response", datanew_b)/mean(smr01c_26$cases.per.week.smr01)*1000

#absolute effects for model 2
#at 1 month
smr01c_26$factual[109]-smr01c_26$c.factual[109]
#at 6 months
smr01c_26$factual[129]-smr01c_26$c.factual[129]

# SMR-01 - child 52w ----
smr01c_52 <- data_cc %>%
  select(time.weeks, level, trend.weeks, episodes.52w.child)

#create denominator (person-weeks)
smr01c_52 <- smr01c_52 %>%
  group_by(time.weeks) %>%
  mutate(cases.per.week.smr01 = n()) 
smr01c_52$cases.per.week.smr01 <- smr01c_52$cases.per.week.smr01*52

#aggregate
smr01c_52 <- smr01c_52 %>%
  group_by(time.weeks) %>%
  summarise(trend.weeks = mean(trend.weeks), 
            level = mean(level), 
            cases.per.week.smr01 = mean(cases.per.week.smr01), 
            episodes.52w.child = sum(episodes.52w.child))

# Incidence rate
# This gives the average number of episodes per 1000 persons per week at risk
smr01c_52$rate <- with(smr01c_52, episodes.52w.child/cases.per.week.smr01)*1000

#create plot of pre-intervention trend
plot(smr01c_52$rate, type ="n", ylim = c(0,15), xlab = "Week of birth", 
     ylab = "Rate x 1000", xaxt="n")
rect(105, 0, 207, 15, col=grey(0.9), border = F)
points(smr01c_52$rate, cex=0.7)
axis(1, at=0:52*4, labels = TRUE)
title("Infant admissions within 52 weeks")

summary(smr01c_52$episodes.52w.child[smr01c_52$level==0])
summary(smr01c_52$episodes.52w.child[smr01c_52$level==1])
summary(smr01c_52$rate)
summary(smr01c_52$rate[smr01c_52$level==0])
summary(smr01c_52$rate[smr01c_52$level==1])

### model 1 - unadjusted
model1 <- glm(episodes.52w.child ~ offset(log(cases.per.week.smr01))
              + time.weeks + level + trend.weeks, family = poisson, smr01c_52)
summary(model1)  
summary(model1)$dispersion
round(ci.exp(model1), 4)


# model2 - allow for over-dispersion
model2 <- glm(episodes.52w.child ~ offset(log(cases.per.week.smr01))
              + time.weeks + level + trend.weeks, family = quasipoisson, smr01c_52)
summary(model2)  
summary(model2)$dispersion  
round(ci.exp(model2), 2)

res2 <- residuals(model2, type = "deviance")
plot(smr01c_52$time.weeks, res2, 
     ylim = c(-10,10), pch=19, cex=0.7, col=grey(0.6),
     main = "Infant admissions within 52 weeks (model 2)", ylab = "Deviance residuals",
     xlab = "Week of birth")
abline(h=0, lty=2, lwd=2)

#check acf & p-acf
acf(res2, 
    main = "Infant admissions within 52 weeks (model 2)")
pacf(res2,
     main = "Infant admissions within 52 weeks (model 2)")

# model 2 predictions
datanew <- data.frame(cases.per.week.smr01=mean(smr01c_52$cases.per.week.smr01),
                      level=rep(c(0,1), c(104,103)),
                      time.weeks=1:207)
datanew <- cbind(datanew,trend_pred)
smr01c_52$factual=predict(model2, type="response", datanew)/mean(smr01c_52$cases.per.week.smr01)*1000

datanew_b <- data.frame(cases.per.week.smr01=mean(smr01c_52$cases.per.week.smr01),
                        level=0,
                        trend.weeks=0,
                        time.weeks=1:207)
smr01c_52$c.factual=predict(model2, type="response", datanew_b)/mean(smr01c_52$cases.per.week.smr01)*1000

#absolute effects for model 2
#at 1 month
smr01c_52$factual[109]-smr01c_52$c.factual[109]
#at 6 months
smr01c_52$factual[129]-smr01c_52$c.factual[129]

# SMR-01 - mother 26w ----
smr01m_26 <- data_cc %>%
  select(time.weeks, level, trend.weeks, episodes.26w.mother)

#create denominator (person-weeks)
smr01m_26 <- smr01m_26 %>%
  group_by(time.weeks) %>%
  mutate(cases.per.week.smr01 = n()) 
smr01m_26$cases.per.week.smr01 <- smr01m_26$cases.per.week.smr01*26

#aggregate
smr01m_26 <- smr01m_26 %>%
  group_by(time.weeks) %>%
  summarise(trend.weeks = mean(trend.weeks), 
            level = mean(level), 
            cases.per.week.smr01 = mean(cases.per.week.smr01), 
            episodes.26w.mother = sum(episodes.26w.mother))

# Incidence rate
# This gives the average number of episodes per 1000 persons per week at risk
smr01m_26$rate <- with(smr01m_26, episodes.26w.mother/cases.per.week.smr01)*1000

#create plot of pre-intervention trend
plot(smr01m_26$rate, type ="n", ylim = c(0,10), xlab = "Week of birth", 
     ylab = "Rate x 1000", xaxt="n")
rect(105, 0, 208, 10, col=grey(0.9), border = F)
points(smr01m_26$rate, cex=0.7)
axis(1, at=0:52*4, labels = TRUE)
title("Maternal admissions within 26 weeks")

summary(smr01m_26$episodes.26w.mother[smr01m_26$level==0])
summary(smr01m_26$episodes.26w.mother[smr01m_26$level==1])
summary(smr01m_26$rate)
summary(smr01m_26$rate[smr01m_26$level==0])
summary(smr01m_26$rate[smr01m_26$level==1])

### model 1 - unadjusted
model1 <- glm(episodes.26w.mother ~ offset(log(cases.per.week.smr01))
              + time.weeks + level + trend.weeks, family = poisson, smr01m_26)
summary(model1)  
summary(model1)$dispersion
round(ci.exp(model1), 3)


# model2 - allow for over-dispersion
model2 <- glm(episodes.26w.mother ~ offset(log(cases.per.week.smr01))
              + time.weeks + level + trend.weeks, family = quasipoisson, smr01m_26)
summary(model2)  
summary(model2)$dispersion  
round(ci.exp(model2), 2)

res2 <- residuals(model2, type = "deviance")
plot(smr01m_26$time.weeks, res2, 
     ylim = c(-10,10), pch=19, cex=0.7, col=grey(0.6),
     main = "Maternal admissions within 26 weeks (model 2)", ylab = "Deviance residuals",
     xlab = "Week of birth")
abline(h=0, lty=2, lwd=2)

#check acf & p-acf
acf(res2,
    main = "Maternal admissions within 26 weeks (model 2)")
pacf(res2,
     main = "Maternal admissions within 26 weeks (model 2)")

# model 2 predictions
datanew <- data.frame(cases.per.week.smr01=mean(smr01m_26$cases.per.week.smr01),
                      level=rep(c(0,1), c(104,103)),
                      time.weeks=1:207)
datanew <- cbind(datanew,trend_pred)
smr01m_26$factual=predict(model2, type="response", datanew)/mean(smr01m_26$cases.per.week.smr01)*1000

datanew_b <- data.frame(cases.per.week.smr01=mean(smr01m_26$cases.per.week.smr01),
                        level=0,
                        trend.weeks=0,
                        time.weeks=1:207)
smr01m_26$c.factual=predict(model2, type="response", datanew_b)/mean(smr01m_26$cases.per.week.smr01)*1000

#absolute effects for model 2
#at 1 month
smr01m_26$factual[109]-smr01m_26$c.factual[109]
#at 6 months
smr01m_26$factual[129]-smr01m_26$c.factual[129]



# SMR-01 - mother 52w ----
smr01m_52 <- data_cc %>%
  select(time.weeks, level, trend.weeks, episodes.52w.mother)

#create denominator (person-weeks)
smr01m_52 <- smr01m_52 %>%
  group_by(time.weeks) %>%
  mutate(cases.per.week.smr01 = n()) 
smr01m_52$cases.per.week.smr01 <- smr01m_52$cases.per.week.smr01*52

#aggregate
smr01m_52 <- smr01m_52 %>%
  group_by(time.weeks) %>%
  summarise(trend.weeks = mean(trend.weeks), 
            level = mean(level), 
            cases.per.week.smr01 = mean(cases.per.week.smr01), 
            episodes.52w.mother = sum(episodes.52w.mother))

# Incidence rate
# This gives the average number of episodes per 1000 persons per week at risk
smr01m_52$rate <- with(smr01m_52, episodes.52w.mother/cases.per.week.smr01)*1000

#create plot of pre-intervention trend
plot(smr01m_52$rate, type ="n", ylim = c(0,10), xlab = "Week of birth", 
     ylab = "Rate x 1000", xaxt="n")
rect(105, 0, 207, 10, col=grey(0.9), border = F)
points(smr01m_52$rate, cex=0.7)
axis(1, at=0:52*4, labels = TRUE)
title("Maternal admissions within 52 weeks")

summary(smr01m_52$episodes.52w.mother[smr01m_52$level==0])
summary(smr01m_52$episodes.52w.mother[smr01m_52$level==1])
summary(smr01m_52$rate)
summary(smr01m_52$rate[smr01m_52$level==0])
summary(smr01m_52$rate[smr01m_52$level==1])

### model 1 - unadjusted
model1 <- glm(episodes.52w.mother ~ offset(log(cases.per.week.smr01))
              + time.weeks + level + trend.weeks, family = poisson, smr01m_52)
summary(model1)  
summary(model1)$dispersion
round(ci.exp(model1), 2)

# model2 - allow for over-dispersion
model2 <- glm(episodes.52w.mother ~ offset(log(cases.per.week.smr01))
              + time.weeks + level + trend.weeks, family = quasipoisson, smr01m_52)
summary(model2)  
summary(model2)$dispersion  
round(ci.exp(model2), 4)

res2 <- residuals(model2, type = "deviance")
plot(smr01m_52$time.weeks, res2, 
     ylim = c(-10,10), pch=19, cex=0.7, col=grey(0.6),
     main = "Maternal admissions within 52 weeks (model 2)", ylab = "Deviance residuals",
     xlab = "Week of birth")
abline(h=0, lty=2, lwd=2)

#check acf & p-acf
acf(res2,
    main = "Maternal admissions within 52 weeks (model 2)")
pacf(res2,
     main = "Maternal admissions within 52 weeks (model 2)")

# model 2 predictions
datanew <- data.frame(cases.per.week.smr01=mean(smr01m_52$cases.per.week.smr01),
                      level=rep(c(0,1), c(104,103)),
                      time.weeks=1:207)
datanew <- cbind(datanew,trend_pred)
smr01m_52$factual=predict(model2, type="response", datanew)/mean(smr01m_52$cases.per.week.smr01)*1000

datanew_b <- data.frame(cases.per.week.smr01=mean(smr01m_52$cases.per.week.smr01),
                        level=0,
                        trend.weeks=0,
                        time.weeks=1:207)
smr01m_52$c.factual=predict(model2, type="response", datanew_b)/mean(smr01m_52$cases.per.week.smr01)*1000

#absolute effects for model 2
#at 1 month
smr01m_52$factual[109]-smr01m_52$c.factual[109]
#at 6 months
smr01m_52$factual[129]-smr01m_52$c.factual[129]


# CHSP_fv - smoke.pc (spc) ----
chsp_fv_spc <- data_cc %>%
  select(time.weeks, level, trend.weeks, smoke.pc.fv)

#check & drop missingness
table(chsp_fv_spc$smoke.pc.fv, useNA = "always")

chsp_fv_spc %>%
  select(time.weeks, smoke.pc.fv) %>%
  gg_miss_fct(time.weeks) + 
  labs(y = "Variable", x="Week of birth") 


sum(!complete.cases(chsp_fv_spc))
182516-177746
chsp_fv_spc <- chsp_fv_spc %>%
  drop_na()

#create denominator
chsp_fv_spc <- chsp_fv_spc %>%
  group_by(time.weeks) %>%
  mutate(cases.per.week = n())

#create binary outcome
chsp_fv_spc <- chsp_fv_spc %>%
  mutate(smoke.pc.bin = ifelse(smoke.pc.fv == "Y", 1, 0))
table(chsp_fv_spc$smoke.pc.bin)

#aggregate
chsp_fv_spc <- chsp_fv_spc %>%
  group_by(time.weeks) %>%
  summarise(trend.weeks = mean(trend.weeks), 
            level = mean(level), 
            cases.per.week = mean(cases.per.week),
            smoke.pc.bin = sum(smoke.pc.bin))

plot(chsp_fv_spc$cases.per.week, type ="n", ylim = c(300,1100), xlab = "Week of birth", 
     ylab = "Observations", xaxt="n")
points(chsp_fv_spc$cases.per.week, cex=0.7)
axis(1, at=0:52*4, labels = TRUE)
title("Primary carer current smoker")

#descriptive statistics

#compute prevalence (episodes per 100 births)

chsp_fv_spc$rate <- with(chsp_fv_spc, smoke.pc.bin/cases.per.week*100)

#create plot of pre-intervention trend
plot(chsp_fv_spc$rate, type ="n", ylim = c(0,30), xlab = "Week of birth", 
     ylab = "Proportion (%)", xaxt="n")
rect(105, 0, 208, 30, col=grey(0.9), border = F)
points(chsp_fv_spc$rate[chsp_fv_spc$level==0], cex=0.7)
axis(1, at=0:52*4, labels = TRUE)
title("Primary carer current smoker")

summary(chsp_fv_spc)
summary(chsp_fv_spc$smoke.pc.bin[chsp_fv_spc$level==0])
summary(chsp_fv_spc$smoke.pc.bin[chsp_fv_spc$level==1])
summary(chsp_fv_spc$rate)
summary(chsp_fv_spc$rate[chsp_fv_spc$level==0])
summary(chsp_fv_spc$rate[chsp_fv_spc$level==1])

### model 1 - unadjusted
model1 <- glm(smoke.pc.bin ~ offset(log(cases.per.week))
              + time.weeks + level + trend.weeks, family = poisson, chsp_fv_spc)
tab_model(model1, vcov.fun = "CL", vcov.type = "HC0", show.se = TRUE, digits = 4)


### model 2 
model2 <- glm(smoke.pc.bin ~ offset(log(cases.per.week))
              + time.weeks + level + trend.weeks, family = quasipoisson, chsp_fv_spc)

tab_model(mode21, vcov.fun = "CL", vcov.type = "HC0", show.se = TRUE, digits = 4)

res2 <- residuals(model2, type = "deviance")
plot(chsp_fv_spc$time.weeks, res2, 
     ylim = c(-10,10), pch=19, cex=0.7, col=grey(0.6),
     main = "Primary carer current smoker (model2)", ylab = "Deviance residuals",
     xlab = "Week of birth")
abline(h=0, lty=2, lwd=2)

#check acf & p-acf
acf(res2,
    main = "Primary carer current smoker (model2)")
pacf(res2,
     main = "Primary carer current smoker (model2)")

#adjust for seasonality
model3 <- glm(smoke.pc.bin ~ offset(log(cases.per.week))
              + time.weeks + level + trend.weeks + harmonic(time.weeks,2,52),
              family = quasipoisson, chsp_fv_spc)

tab_model(model3, vcov.fun = "CL", vcov.type = "HC0", show.se = TRUE, digits = 4)

res3 <- residuals(model3, type = "deviance")
plot(chsp_fv_spc$time.weeks, res3, 
     ylim = c(-10,10), pch=19, cex=0.7, col=grey(0.6),
     main = "Primary carer current smoker (model 3)", ylab = "Deviance residuals",
     xlab = "Week of birth")
abline(h=0, lty=2, lwd=2)

#check acf & p-acf
acf(res3,
    main = "Primary carer current smoker (model 3)")
pacf(res3,
     main = "Primary carer current smoker (model 3)")

# model 3 predictions
datanew <- data.frame(cases.per.week=mean(chsp_fv_spc$cases.per.week),
                      level=rep(c(0,1), c(104,103)),
                      time.weeks=1:207)
datanew <- cbind(datanew,trend_pred)
factual=add_ci(datanew, model3, names = c("lcb", "ucb"), alpha = 0.05)/mean(chsp_fv_spc$cases.per.week)*100

datanew_b <- data.frame(cases.per.week=mean(chsp_fv_spc$cases.per.week),
                        level=0,
                        trend.weeks=0,
                        time.weeks=1:207)
chsp_fv_spc$c.factual=predict(model3, type="response", datanew_b)/mean(chsp_fv_spc$cases.per.week)*100


#absolute effects for model 3
#at 1 month
factual$pred[109]-chsp_fv_spc$c.factual[109]
#at 6 months
factual$pred[129]-chsp_fv_spc$c.factual[129]


# plot model
date <- read_csv("Dates.csv")
chsp_fv_spc <- cbind(chsp_fv_spc,date)
chsp_fv_spc <- chsp_fv_spc %>%
  mutate(Week.beginning = as.Date(Week.beginning, format = "%Y/%m/%d"))
  
#modify factual for segments
factual <- factual %>%
  mutate(firstseg = if_else(trend.weeks == 0, pred, 0)) %>%
  mutate(firstlcb = if_else(trend.weeks == 0, lcb, 0)) %>%
  mutate(firstucb = if_else(trend.weeks==0, ucb,0)) %>%
  mutate(seclcb = if_else(trend.weeks>0, lcb,0)) %>%
  mutate(secucb = if_else(trend.weeks>0, ucb,0)) %>%
  mutate(secseg = if_else(trend.weeks> 0, pred, 0))

factual$firstseg[factual$firstseg==0] <- NA 
factual$secseg[factual$secseg==0] <- NA
factual$firstlcb[factual$firstlcb==0] <- NA
factual$firstucb[factual$firstucb==0] <- NA
factual$secucb[factual$secucb==0] <- NA
factual$seclcb[factual$seclcb==0] <- NA


ggplot(chsp_fv_spc, aes(y=rate, x=Week.beginning, linetype = "solid")) +
  geom_ribbon(aes(y=factual$firstseg, ymin=factual$firstlcb, ymax=factual$firstucb), fill = "#CCCCCC") +
  geom_ribbon(aes(y=factual$secseg, ymin=factual$seclcb, ymax=factual$secucb), fill = "#CCCCCC") +
  geom_line() +
  geom_line(aes(y=factual$firstseg, linetype = "dashed"), colour="#990999", size = 1.3) +
  geom_line(aes(y=factual$secseg, linetype = "dashed"), colour="#990999", size = 1.3) +
  xlab(expression(bold(paste("Year")))) +
  ylab(expression(bold(paste("Prevalence (%)")))) +
  scale_linetype_manual(name = "Trend", values = c("solid", "dashed"), labels = c("Modelled", "Actual")) +
  theme(axis.text.x = element_text(size=12, face = "bold"),
        axis.text.y = element_text(size=12, face = "bold"),
        axis.title.x = element_text(size=15, face = "bold"),
        axis.title.y = element_text(size=15, face = "bold"),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        panel.border = element_rect(colour = "black", fill=NA, size=1.2),
        legend.position = "none") 



  


# CHSP_68 - smoke.pc (spc) ####
chsp_68_spc <- data_cc %>%
  select(time.weeks, level, trend.weeks, smoke.pc.68, simd.quintile.2016, mother.age, tpp_bin)

#check & drop missingness
table(chsp_68_spc$smoke.pc.68, useNA = "always")

chsp_68_spc %>%
  select(time.weeks, smoke.pc.68) %>%
  gg_miss_fct(time.weeks) + 
  labs(y = "Variable", x="Week of birth") 

sum(!complete.cases(chsp_68_spc))


#create binary indicating missingness in smoke.pc.68
chsp_68_spc <- chsp_68_spc %>%
  mutate(miss_bin = as.integer(complete.cases(chsp_68_spc$smoke.pc.68)))
table(chsp_68_spc$miss_bin)

setwd("Results")
spc_68_miss_simd <- CrossTable(chsp_68_spc$simd.quintile.2016, chsp_68_spc$miss_bin)
write.csv(spc_68_miss_simd, "spc_68_miss_simd.csv")

spc_68_miss_age <- CrossTable(chsp_68_spc$mother.age, chsp_68_spc$miss_bin)
write.csv(spc_68_miss_age, "spc_68_miss_age.csv")

spc_68_miss_npp <- CrossTable(chsp_68_spc$tpp_bin, chsp_68_spc$miss_bin)


chsp_68_spc <- chsp_68_spc %>%
  drop_na()


#create denominator
chsp_68_spc <- chsp_68_spc %>%
  group_by(time.weeks) %>%
  mutate(cases.per.week = n())

#create binary outcome
chsp_68_spc <- chsp_68_spc %>%
  mutate(smoke.pc.bin = ifelse(smoke.pc.68 == "Y", 1, 0))
table(chsp_68_spc$smoke.pc.bin)

#truncate timeframe
chsp_68_spc <- chsp_68_spc %>%
  filter(time.weeks >25)

#aggregate
chsp_68_spc <- chsp_68_spc %>%
  group_by(time.weeks) %>%
  summarise(trend.weeks = mean(trend.weeks), 
            level = mean(level), 
            cases.per.week = mean(cases.per.week),
            smoke.pc.bin = sum(smoke.pc.bin))

plot(chsp_68_spc$cases.per.week, type ="n", ylim = c(300,1100), xlab = "Week of birth", 
     ylab = "Observations", xaxt="n")
points(chsp_68_spc$cases.per.week, cex=0.7)
axis(1, at=0:52*4, labels = TRUE)
title("Primary carer current smoker")

#add truncated time variable
time_truncated <- read_csv("time_truncated.csv")
chsp_68_spc <- cbind(chsp_68_spc,time_truncated)
#descriptive statistics

#compute prevalence (episodes per 100 births)

chsp_68_spc$rate <- with(chsp_68_spc, smoke.pc.bin/cases.per.week*100)

#create plot of pre-intervention trend
plot(chsp_68_spc$rate, type ="n", ylim = c(0,30), xlab = "Week of birth", 
     ylab = "Proportion (%)", xaxt="n")
points(chsp_68_spc$rate, cex=0.7)
axis(1, at=0:52*4, labels = TRUE)
title("Primary carer current smoker")

summary(chsp_68_spc)
summary(chsp_68_spc$smoke.pc.bin[chsp_68_spc$level==0])
summary(chsp_68_spc$smoke.pc.bin[chsp_68_spc$level==1])
summary(chsp_68_spc$rate)
summary(chsp_68_spc$rate[chsp_68_spc$level==0])
summary(chsp_68_spc$rate[chsp_68_spc$level==1])

### model 1 - unadjusted
model1 <- glm(smoke.pc.bin ~ offset(log(cases.per.week))
              + time.truncated + level + trend.weeks, family = poisson, chsp_68_spc)

tab_model(model1, vcov.fun = "CL", vcov.type = "HC0", show.se = TRUE, digits = 4)


### model 2 
model2 <- glm(smoke.pc.bin ~ offset(log(cases.per.week))
              + time.truncated + level + trend.weeks, family = quasipoisson, chsp_68_spc)
summary(model2)

tab_model(model2, vcov.fun = "CL", vcov.type = "HC0", show.se = TRUE, digits = 4)

#check acf & p-acf


res2 <- residuals(model2, type = "deviance")
par(mfrow=c(1,1))

plot(chsp_68_spc$time.truncated, res2, 
     ylim = c(-7,7), pch=19, cex=0.7, col=grey(0.6),
     main = "Primary carer smoker 6-8 (model2)", ylab = "deviance residuals",
     xlab = "Week of birth") + abline(h=0, lty=2, lwd=2)
acf(res2,
    main = "Primary carer smoker 6-8 (model2)")
pacf(res2,
     main = "Primary carer smoker 6-8 (model2)")

#model 2 predictions
trend_pred <- read_csv("trend_weeks_truncated.csv")
datanew <- data.frame(cases.per.week=mean(chsp_68_spc$cases.per.week),
                      level=rep(c(0,1), c(79,103)),
                      time.truncated=1:182)
datanew <- cbind(datanew,trend_pred)
chsp_68_spc$factual=predict(model2, type="response", datanew)/mean(chsp_68_spc$cases.per.week)*100

datanew_b <- data.frame(cases.per.week=mean(chsp_68_spc$cases.per.week),
                        level=0,
                        trend.weeks=0,
                        time.truncated=1:182)
chsp_68_spc$c.factual=predict(model2, type="response", datanew_b)/mean(chsp_68_spc$cases.per.week)*100

#absolute effects for model 2
#at 1 month
chsp_68_spc$factual[84]-chsp_68_spc$c.factual[84]
#at 6 months
chsp_68_spc$factual[104]-chsp_68_spc$c.factual[104]

# CHSP_fv - smoke.sh (ssh)----
chsp_fv_ssh <- data_cc %>%
  select(time.weeks, level, trend.weeks, smoke.sh.fv)

#check & drop missingness
table(chsp_fv_ssh$smoke.sh.fv, useNA = "always")
chsp_fv_ssh %>% 
  select(smoke.sh.fv, time.weeks) %>%
  gg_miss_fct(time.weeks)+ 
  labs(y = "Variable", x="Week of birth") 

sum(!complete.cases(chsp_fv_ssh))


chsp_fv_ssh <- chsp_fv_ssh %>%
  drop_na()

#create denominator
chsp_fv_ssh <- chsp_fv_ssh %>%
  group_by(time.weeks) %>%
  mutate(cases.per.week = n())

#create binary outcome
chsp_fv_ssh <- chsp_fv_ssh %>%
  mutate(smoke.sh.bin = ifelse(smoke.sh.fv == "Y", 1, 0))
table(chsp_fv_ssh$smoke.sh.bin)

#aggregate
chsp_fv_ssh <- chsp_fv_ssh %>%
  group_by(time.weeks) %>%
  summarise(trend.weeks = mean(trend.weeks), 
            level = mean(level), 
            cases.per.week = mean(cases.per.week),
            smoke.sh.bin = sum(smoke.sh.bin))

plot(chsp_fv_ssh$cases.per.week, type ="n", ylim = c(550,1100), xlab = "Week of birth", 
     ylab = "Observations", xaxt="n")
points(chsp_fv_ssh$cases.per.week, cex=0.7)
axis(1, at=0:52*4, labels = TRUE)
title("Exposure to second-hand smoke")


### descriptive statistics
#compute rates (episodes per 1000 births)
chsp_fv_ssh$rate <- with(chsp_fv_ssh, smoke.sh.bin/cases.per.week*100)
#create plot of pre-intervention trend
plot(chsp_fv_ssh$rate, type ="n", ylim = c(0,25), xlab = "Week of birth", 
     ylab = "Proportion (%)", xaxt="n")
rect(105, 0, 208, 25, col=grey(0.9), border = F)
points(chsp_fv_ssh$rate, cex=0.7)
axis(1, at=0:52*4, labels = TRUE)
title("Exposure to second-hand smoke")

summary(chsp_fv_ssh)
summary(chsp_fv_ssh$smoke.sh.bin[chsp_fv_ssh$level==0])
summary(chsp_fv_ssh$smoke.sh.bin[chsp_fv_ssh$level==1])
summary(chsp_fv_ssh$rate[chsp_fv_ssh$level==0])
summary(chsp_fv_ssh$rate[chsp_fv_ssh$level==1])
summary(chsp_fv_ssh$rate)
### model 1 - unadjusted
model1 <- glm(smoke.sh.bin ~ offset(log(cases.per.week))
              + time.weeks + level + trend.weeks, family = poisson, chsp_fv_ssh)
summary(model1)  
tab_model(model1, vcov.fun = "CL", vcov.type = "HC0", show.se = TRUE, digits = 4)


res1 <- residuals(model1, type = "deviance")
plot(chsp_fv_ssh$time.weeks, res1, 
     ylim = c(-10,10), pch=19, cex=0.7, col=grey(0.6),
     main = "Exposure to second-hand smoke (model 1)", ylab = "Deviance residuals",
     xlab = "Week of birth")
abline(h=0, lty=2, lwd=2)

#check acf & p-acf
acf(res1,
    main = "Exposure to second-hand smoke (model 1)")
pacf(res1,
     main = "Exposure to second-hand smoke (model 1)")


#adjust for seasonality
model2 <- glm(smoke.sh.bin ~ offset(log(cases.per.week))
              + time.weeks + level + trend.weeks + harmonic(time.weeks,2,52),
              family = poisson, chsp_fv_ssh)
summary(model2)  
tab_model(model2, vcov.fun = "CL", vcov.type = "HC0", show.se = TRUE, digits = 4)

res2 <- residuals(model2, type = "deviance")
plot(chsp_fv_ssh$time.weeks, res2, 
     ylim = c(-7,7), pch=19, cex=0.7, col=grey(0.6),
     main = "Exposure to second-hand smoke (model 2)", ylab = "Deviance residuals",
     xlab = "Week of birth")
abline(h=0, lty=2, lwd=2)

#check acf & p-acf
acf(res2,
    main = "Exposure to second-hand smoke (model 2)")
pacf(res2,
     main = "Exposure to second-hand smoke (model 2)")

#model 2 predictions
datanew <- data.frame(cases.per.week=mean(chsp_fv_ssh$cases.per.week),
                      level=rep(c(0,1), c(104,103)),
                      time.weeks=1:207)
datanew <- cbind(datanew,trend_pred)
chsp_fv_ssh$factual=predict(model2, type="response", datanew)/mean(chsp_fv_ssh$cases.per.week)*100
factual=add_ci(datanew, model2, names = c("lcb", "ucb"), alpha = 0.05)/mean(chsp_fv_ssh$cases.per.week)*100


datanew_b <- data.frame(cases.per.week=mean(chsp_fv_ssh$cases.per.week),
                        level=0,
                        trend.weeks=0,
                        time.weeks=1:207)
chsp_fv_ssh$c.factual=predict(model2, type="response", datanew_b)/mean(chsp_fv_ssh$cases.per.week)*100

#absolute effects for model 2
#at 1 month
chsp_fv_ssh$factual[109]-chsp_fv_ssh$c.factual[109]
#at 6 months
chsp_fv_ssh$factual[129]-chsp_fv_ssh$c.factual[129]

#plot
date <- read_csv("Dates.csv")
chsp_fv_ssh <- cbind(chsp_fv_ssh,date)
chsp_fv_ssh <- chsp_fv_ssh %>%
  mutate(Week.beginning = as.Date(Week.beginning, format = "%Y/%m/%d"))


factual <- factual %>%
  mutate(firstseg = if_else(trend.weeks == 0, pred, 0)) %>%
  mutate(firstlcb = if_else(trend.weeks == 0, lcb, 0)) %>%
  mutate(firstucb = if_else(trend.weeks==0, ucb,0)) %>%
  mutate(seclcb = if_else(trend.weeks>0, lcb,0)) %>%
  mutate(secucb = if_else(trend.weeks>0, ucb,0)) %>%
  mutate(secseg = if_else(trend.weeks> 0, pred, 0))

factual$firstseg[factual$firstseg==0] <- NA 
factual$secseg[factual$secseg==0] <- NA
factual$firstlcb[factual$firstlcb==0] <- NA
factual$firstucb[factual$firstucb==0] <- NA
factual$secucb[factual$secucb==0] <- NA
factual$seclcb[factual$seclcb==0] <- NA


chsp_fv_ssh <- chsp_fv_ssh %>% add_column(firstseg = factual$firstseg,
                                          secseg = factual$secseg,
                                          firstlcb = factual$firstlcb,
                                          firstucb = factual$firstucb,
                                          secucb = factual$secucb,
                                          seclcb = factual$seclcb)
  

ggplot(chsp_fv_ssh, aes(y=rate, x=Week.beginning, linetype = "solid")) +
  geom_ribbon(aes(y=chsp_fv_ssh$secseg, ymin=chsp_fv_ssh$seclcb, ymax=chsp_fv_ssh$secucb), fill = "#CCCCCC") +
  geom_ribbon(aes(y=chsp_fv_ssh$firstseg, ymin=chsp_fv_ssh$firstlcb, ymax=chsp_fv_ssh$firstucb), fill = "#CCCCCC") +
  geom_line() +
  geom_line(aes(y=chsp_fv_ssh$firstseg, linetype = "dashed"), colour="#CC0000", size = 1.3) +
  geom_line(aes(y=chsp_fv_ssh$secseg, linetype = "dashed"), colour="#CC0000", size = 1.3) +
  xlab(expression(bold(paste("Year")))) +
  ylab(expression(bold(paste("Prevalence (%)")))) +
  scale_linetype_manual(name = "Trend", values = c("solid", "dashed"), labels = c("Modelled", "Actual")) +
  theme(axis.text.x = element_text(size=12, face = "bold"),
        axis.text.y = element_text(size=12, face = "bold"),
        axis.title.x = element_text(size=15, face = "bold"),
        axis.title.y = element_text(size=15, face = "bold"),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        panel.border = element_rect(colour = "black", fill=NA, size=1.2),
        legend.position = "none") 





#ggplot(chsp_fv_ssh, aes(y=rate, x=Week.beginning, linetype = "solid")) +
#  geom_ribbon(aes(y=factual_ssh$pred, ymin=factual_ssh$lcb, ymax=factual_ssh$ucb), fill = "#CCCCCC") +
#  geom_line() +
#  geom_line(aes(y=factual_ssh$pred, linetype = "dashed"), colour="#CC0000", size = 1.3) +
#  xlab(expression(bold(paste("Year")))) +
#  ylab(expression(bold(paste("Prevalence (%)")))) +
#  scale_linetype_manual(name = "Trend", values = c("solid", "dashed"), labels = c("Modelled", "Actual")) +
#  theme(axis.text.x = element_text((size=15)),
#        axis.text.y = element_text((size=15)),
#        axis.title.x = element_text((size=17)),
#        axis.title.y = element_text((size=17)),
#        panel.background = element_blank(),
#        panel.grid.major = element_blank(),
#        panel.grid.minor = element_blank(),
#        axis.line = element_line(colour = "black"),
#        panel.border = element_rect(colour = "black", fill=NA, size=1.2),
#        legend.position = "none")





# CHSP_68 - smoke.sh (ssh)####
chsp_68_ssh <- data_cc %>%
  select(time.weeks, level, trend.weeks, smoke.sh.68, simd.quintile.2016, mother.age, tpp_bin)

#check & drop missingness
table(chsp_68_ssh$smoke.sh.68, useNA = "always")
chsp_68_ssh %>% 
  select(smoke.sh.68, time.weeks) %>%
  gg_miss_fct(time.weeks)+ 
  labs(y = "Variable", x="Week of birth") 


#create binary indicating missingness in smoke.pc.68
chsp_68_ssh <- chsp_68_ssh %>%
  mutate(miss_bin = as.integer(complete.cases(chsp_68_ssh$smoke.sh.68)))
table(chsp_68_ssh$miss_bin)

setwd("Results")
ssh_68_miss_simd <- CrossTable(chsp_68_ssh$simd.quintile.2016, chsp_68_ssh$miss_bin)
write.csv(ssh_68_miss_simd, "ssh_68_miss_simd.csv")

ssh_68_miss_age <- CrossTable(chsp_68_ssh$mother.age, chsp_68_ssh$miss_bin)
write.csv(ssh_68_miss_age, "ssh_68_miss_age.csv")

ssh_68_miss_npp <- CrossTable(chsp_68_ssh$tpp_bin, chsp_68_ssh$miss_bin)


sum(!complete.cases(chsp_68_ssh))
chsp_68_ssh <- chsp_68_ssh %>%
  drop_na()

#create denominator
chsp_68_ssh <- chsp_68_ssh %>%
  group_by(time.weeks) %>%
  mutate(cases.per.week = n())

#create binary outcome
chsp_68_ssh <- chsp_68_ssh %>%
  mutate(smoke.sh.bin = ifelse(smoke.sh.68 == "Y", 1, 0))
table(chsp_68_ssh$smoke.sh.bin)

#aggregate
chsp_68_ssh <- chsp_68_ssh %>%
  group_by(time.weeks) %>%
  summarise(trend.weeks = mean(trend.weeks), 
            level = mean(level), 
            cases.per.week = mean(cases.per.week),
            smoke.sh.bin = sum(smoke.sh.bin))

plot(chsp_68_ssh$cases.per.week, type ="n", ylim = c(500,1100), xlab = "Week of birth", 
     ylab = "Observations", xaxt="n")
points(chsp_68_ssh$cases.per.week, cex=0.7)
axis(1, at=0:52*4, labels = TRUE)
title("Exposure to second-hand smoke")


### descriptive statistics
#compute rates (episodes per 1000 births)
chsp_68_ssh$rate <- with(chsp_68_ssh, smoke.sh.bin/cases.per.week*100)

#create plot of pre-intervention trend
plot(chsp_68_ssh$rate, type ="n", ylim = c(0,20), xlab = "Week of birth", 
     ylab = "Proportion (%)", xaxt="n")
rect(105, 0, 208, 20, col=grey(0.9), border = F)
points(chsp_68_ssh$rate, cex=0.7)
axis(1, at=0:52*4, labels = TRUE)
title("Exposure to second-hand smoke 6-8 weeks")

summary(chsp_68_ssh$smoke.sh.bin[chsp_68_ssh$level==0])
summary(chsp_68_ssh$smoke.sh.bin[chsp_68_ssh$level==1])
summary(chsp_68_ssh$rate[chsp_68_ssh$level==0])
summary(chsp_68_ssh$rate[chsp_68_ssh$level==1])

### model 1 - unadjusted
model1 <- glm(smoke.sh.bin ~ offset(log(cases.per.week))
              + time.weeks + level + trend.weeks, family = poisson, chsp_68_ssh)
summary(model1)
tab_model(model1, vcov.fun = "CL", vcov.type = "HC0", show.se = TRUE, digits = 4)

# model 2 
model2 <- glm(smoke.sh.bin ~ offset(log(cases.per.week))
              + time.weeks + level + trend.weeks, family = quasipoisson, chsp_68_ssh)

summary(model2)
tab_model(model2, vcov.fun = "CL", vcov.type = "HC0", show.se = TRUE, digits = 4)

#check acf & p-acf
res2 <- residuals(model2, type = "deviance")
acf(res2,
    main = "Exposure to second-hand smoke (model 2)")
pacf(res2,
     main = "Exposure to second-hand smoke (model 2)")

# model 3 
model3 <- glm(smoke.sh.bin ~ offset(log(cases.per.week))
              + time.weeks + level + trend.weeks + harmonic(time.weeks,2,52), family = quasipoisson, chsp_68_ssh)

summary(model3)
tab_model(model3, vcov.fun = "CL", vcov.type = "HC0", show.se = TRUE, digits = 4)

#check acf & p-acf
res3 <- residuals(model3, type = "deviance")
plot(chsp_68_ssh$time.weeks, res3, 
     ylim = c(-7,7), pch=19, cex=0.7, col=grey(0.6),
     main = "Exposure to second-hand smoke 6-8 (model3)", ylab = "deviance residuals",
     xlab = "Week of birth") + abline(h=0, lty=2, lwd=2)
acf(res3,
    main = "Exposure to second-hand smoke 6-8 (model 3)")
pacf(res3,
     main = "Exposure to second-hand smoke 6-8 (model 3)")

#model 3 predictions
datanew <- data.frame(cases.per.week=mean(chsp_68_ssh$cases.per.week),
                      level=rep(c(0,1), c(104,103)),
                      time.weeks=1:207)
datanew <- cbind(datanew,trend_pred)
chsp_68_ssh$factual=predict(model3, type="response", datanew)/mean(chsp_68_ssh$cases.per.week)*100

datanew_b <- data.frame(cases.per.week=mean(chsp_68_ssh$cases.per.week),
                        level=0,
                        trend.weeks=0,
                        time.weeks=1:207)
chsp_68_ssh$c.factual=predict(model3, type="response", datanew_b)/mean(chsp_68_ssh$cases.per.week)*100

#absolute effects for model 2
#at 1 month
chsp_68_ssh$factual[109]-chsp_68_ssh$c.factual[109]
#at 6 months
chsp_68_ssh$factual[129]-chsp_68_ssh$c.factual[129]


# CHSP_68 - smoke.sh (ssh) - truncated ----
chsp_68_ssh <- data_cc %>%
  select(time.weeks, level, trend.weeks, smoke.sh.68, simd.quintile.2016, mother.age)

chsp_68_ssh <- chsp_68_ssh %>%
  drop_na()

#create denominator
chsp_68_ssh <- chsp_68_ssh %>%
  group_by(time.weeks) %>%
  mutate(cases.per.week = n())

#create binary outcome
chsp_68_ssh <- chsp_68_ssh %>%
  mutate(smoke.sh.bin = ifelse(smoke.sh.68 == "Y", 1, 0))
table(chsp_68_ssh$smoke.sh.bin)

#truncate timeframe
chsp_68_ssh <- chsp_68_ssh %>%
  filter(time.weeks >25)

#aggregate
chsp_68_ssh <- chsp_68_ssh %>%
  group_by(time.weeks) %>%
  summarise(trend.weeks = mean(trend.weeks), 
            level = mean(level), 
            cases.per.week = mean(cases.per.week),
            smoke.sh.bin = sum(smoke.sh.bin))


#add truncated time variable
time_truncated <- read_csv("time_truncated.csv")
chsp_68_ssh <- cbind(chsp_68_ssh,time_truncated)
#descriptive statistics

#compute prevalence (episodes per 100 births)

chsp_68_ssh$rate <- with(chsp_68_ssh, smoke.sh.bin/cases.per.week*100)

#create plot of pre-intervention trend
plot(chsp_68_ssh$rate, type ="n", ylim = c(0,20), xlab = "Week of birth", 
     ylab = "Proportion (%)", xaxt="n")
points(chsp_68_ssh$rate, cex=0.7)
axis(1, at=0:52*4, labels = TRUE)
title("Primary carer current smoker")

summary(chsp_68_ssh$rate[chsp_68_ssh$level==0])
summary(chsp_68_ssh$rate[chsp_68_ssh$level==1])

### model 1 - unadjusted
model1 <- glm(smoke.sh.bin ~ offset(log(cases.per.week))
              + time.truncated + level + trend.weeks, family = quasipoisson, chsp_68_ssh)

summary(model1)
tab_model(model1, vcov.fun = "CL", vcov.type = "HC0", show.se = TRUE, digits = 4)
#check acf & p-acf
res1 <- residuals(model1, type = "deviance")
acf(res1,
    main = "Exposure to second-hand smoke tr (model 1)")
pacf(res1,
     main = "Exposure to second-hand smoke tr (model 1)")

#model 1 predictions
datanew <- data.frame(cases.per.week=mean(chsp_68_ssh$cases.per.week),
                      level=rep(c(0,1), c(79,103)),
                      time.truncated=1:182)
trend_pred <- read_csv("trend_weeks_truncated.csv")
datanew <- cbind(datanew,trend_pred)
chsp_68_ssh$factual=predict(model1, type="response", datanew)/mean(chsp_68_ssh$cases.per.week)*100

datanew_b <- data.frame(cases.per.week=mean(chsp_68_ssh$cases.per.week),
                        level=0,
                        trend.weeks=0,
                        time.truncated=1:182)
chsp_68_ssh$c.factual=predict(model1, type="response", datanew_b)/mean(chsp_68_ssh$cases.per.week)*100

#absolute effects for model 2
#at 1 month
chsp_68_ssh$factual[84]
chsp_68_ssh$factual[84]-chsp_68_ssh$c.factual[84]
#at 6 months
chsp_68_ssh$factual[104]-chsp_68_ssh$c.factual[104]

# SMR-02 - book.smoking.history ----
smr02_bsh <- data_cc %>%
  select(time.weeks, level, trend.weeks, book.smoking.history)
table(smr02_bsh$book.smoking.history, useNA = "always")

#create denominator
smr02_bsh <- smr02_bsh %>%
  group_by(time.weeks) %>%
  mutate(cases.per.week = n())

#create binary outcome
smr02_bsh <- smr02_bsh %>%
  mutate(bsh.bin = ifelse(book.smoking.history == 1, 1, 0))
table(smr02_bsh$bsh.bin)

#aggregate

smr02_bsh <- smr02_bsh %>%
  group_by(time.weeks) %>%
  summarise(trend.weeks = mean(trend.weeks), 
            level = mean(level), 
            cases.per.week = mean(cases.per.week),
            bsh.bin = sum(bsh.bin))


### descriptive statistics
#compute rates (episodes per 1000 births)
smr02_bsh$rate <- with(smr02_bsh, bsh.bin/cases.per.week*100)
#create plot of pre-intervention trend
plot(smr02_bsh$rate, type ="n", ylim = c(0,25), xlab = "Week of birth", 
     ylab = "Proportion (%)", xaxt="n")
rect(105, 0, 208, 25, col=grey(0.9), border = F)
points(smr02_bsh$rate, cex=0.7)
axis(1, at=0:52*4, labels = TRUE)
title("Current smoker at booking")

summary(smr02_bsh$bsh.bin[smr02_bsh$level==0])
summary(smr02_bsh$bsh.bin[smr02_bsh$level==1])
summary(smr02_bsh$rate[smr02_bsh$level==0])
summary(smr02_bsh$rate[smr02_bsh$level==1])
summary(smr02_bsh$rate)
### model 1 - unadjusted
model1 <- glm(bsh.bin ~ offset(log(cases.per.week))
              + time.weeks + level + trend.weeks, family = poisson, smr02_bsh)
summary(model1)  
tab_model(model1, vcov.fun = "CL", vcov.type = "HC0", show.se = TRUE, digits = 4)


### model 2 
model2 <- glm(bsh.bin ~ offset(log(cases.per.week))
              + time.weeks + level + trend.weeks, family = quasipoisson, smr02_bsh)
summary(model2)  
tab_model(model2, vcov.fun = "CL", vcov.type = "HC0", show.se = TRUE, digits = 4)




res2 <- residuals(model2, type = "deviance")
plot(smr02_bsh$time.weeks, res2, 
     ylim = c(-7,7), pch=19, cex=0.7, col=grey(0.6),
     main = "Current smoker at booking (model 2)", ylab = "Deviance residuals",
     xlab = "Week of birth")
abline(h=0, lty=2, lwd=2)

#check acf & p-acf
acf(res2,
    main = "Current smoker at booking (model 2)")
pacf(res2,
     main = "Current smoker at booking (model 2)")



#adjust for seasonality
model3 <- glm(bsh.bin ~ offset(log(cases.per.week))
              + time.weeks + level + trend.weeks + harmonic(time.weeks,2,52),
              family = quasipoisson, smr02_bsh)

summary(model3)
tab_model(model3, vcov.fun = "CL", vcov.type = "HC0", show.se = TRUE, digits = 4)


res3 <- residuals(model3, type = "deviance")
plot(smr02_bsh$time.weeks, res3, 
     ylim = c(-7,7), pch=19, cex=0.7, col=grey(0.6),
     main = "Current smoker at booking (model 3)", ylab = "Deviance residuals",
     xlab = "Week of birth")
abline(h=0, lty=2, lwd=2)

#check acf & p-acf
acf(res3,
    main = "Current smoker at booking (model 3)")
pacf(res3,
     main = "Current smoker at booking (model 3)")

#model 3 predictions
datanew <- data.frame(cases.per.week=mean(smr02_bsh$cases.per.week),
                      level=rep(c(0,1), c(104,103)),
                      time.weeks=1:207)
datanew <- cbind(datanew,trend_pred)
smr02_bsh$factual=predict(model3, type="response", datanew)/mean(smr02_bsh$cases.per.week)*100

datanew_b <- data.frame(cases.per.week=mean(smr02_bsh$cases.per.week),
                        level=0,
                        trend.weeks=0,
                        time.weeks=1:207)
smr02_bsh$c.factual=predict(model3, type="response", datanew_b)/mean(smr02_bsh$cases.per.week)*100

#absolute effects for model 2
#at 1 month
smr02_bsh$factual[109]-smr02_bsh$c.factual[109]
#at 6 months
smr02_bsh$factual[129]-smr02_bsh$c.factual[129]


# SMR-02 - smoking.during.pregnancy ----
smr02_sp <- data_cc %>%
  select(time.weeks, level, trend.weeks, smoking.during.pregnancy)
table(smr02_sp$smoking.during.pregnancy, useNA = "always")

#create denominator
smr02_sp <- smr02_sp %>%
  group_by(time.weeks) %>%
  mutate(cases.per.week = n())

#aggregate

smr02_sp <- smr02_sp %>%
  group_by(time.weeks) %>%
  summarise(trend.weeks = mean(trend.weeks), 
            level = mean(level), 
            cases.per.week = mean(cases.per.week),
            smoking.during.pregnancy = sum(smoking.during.pregnancy))


### descriptive statistics
#compute rates (episodes per 1000 births)
smr02_sp$rate <- with(smr02_sp, smoking.during.pregnancy/cases.per.week*100)
#create plot of pre-intervention trend
plot(smr02_sp$rate, type ="n", ylim = c(0,30), xlab = "Week of birth", 
     ylab = "Proportion (%)", xaxt="n")
rect(105, 0, 208, 30, col=grey(0.9), border = F)
points(smr02_sp$rate, cex=0.7)
axis(1, at=0:52*4, labels = TRUE)
title("Smoker during pregnancy")

summary(smr02_sp$smoking.during.pregnancy[smr02_sp$level==0])
summary(smr02_sp$smoking.during.pregnancy[smr02_sp$level==1])
summary(smr02_sp$rate[smr02_sp$level==0])
summary(smr02_sp$rate[smr02_sp$level==1])
summary(smr02_sp$rate)

### model 1 - unadjusted
model1 <- glm(smoking.during.pregnancy ~ offset(log(cases.per.week))
              + time.weeks + level + trend.weeks, family = poisson, smr02_sp)
summary(model1)  
tab_model(model1, vcov.fun = "CL", vcov.type = "HC0", show.se = TRUE, digits = 4)


### model 2 
model2 <- glm(smoking.during.pregnancy ~ offset(log(cases.per.week))
              + time.weeks + level + trend.weeks, family = quasipoisson, smr02_sp)
summary(model2)  
tab_model(model2, vcov.fun = "CL", vcov.type = "HC0", show.se = TRUE, digits = 4)


res2 <- residuals(model2, type = "deviance")
plot(smr02_sp$time.weeks, res2, 
     ylim = c(-7,7), pch=19, cex=0.7, col=grey(0.6),
     main = "Smoker during pregnancy (model 2)", ylab = "Deviance residuals",
     xlab = "Week of birth")
abline(h=0, lty=2, lwd=2)

#check acf & p-acf
acf(res2,
    main = "Smoker during pregnancy (model 2)")
pacf(res2,
     main = "Smoker during pregnancy (model 2)")

#adjust for seasonality
model3 <- glm(smoking.during.pregnancy ~ offset(log(cases.per.week))
              + time.weeks + level + trend.weeks + harmonic(time.weeks,2,52),
              family = quasipoisson, smr02_sp)
summary(model3)
tab_model(model3, vcov.fun = "CL", vcov.type = "HC0", show.se = TRUE, digits = 4)


res3 <- residuals(model3, type = "deviance")
plot(smr02_sp$time.weeks, res3, 
     ylim = c(-7,7), pch=19, cex=0.7, col=grey(0.6),
     main = "Smoker during pregnancy (model 3)", ylab = "Deviance residuals",
     xlab = "Week of birth")
abline(h=0, lty=2, lwd=2)

#check acf & p-acf
acf(res3,
    main = "Smoker during pregnancy (model 3)")
pacf(res3,
     main = "Smoker during pregnancy (model 3)")

#model 3 predictions
datanew <- data.frame(cases.per.week=mean(smr02_sp$cases.per.week),
                      level=rep(c(0,1), c(104,103)),
                      time.weeks=1:207)
datanew <- cbind(datanew,trend_pred)
smr02_sp$factual=predict(model3, type="response", datanew)/mean(smr02_sp$cases.per.week)*100

datanew_b <- data.frame(cases.per.week=mean(smr02_sp$cases.per.week),
                        level=0,
                        trend.weeks=0,
                        time.weeks=1:207)
smr02_sp$c.factual=predict(model3, type="response", datanew_b)/mean(smr02_sp$cases.per.week)*100

#absolute effects for model 2 
#at 1 month
smr02_sp$factual[109]-smr02_sp$c.factual[109]
#at 6 months
smr02_sp$factual[129]-smr02_sp$c.factual[129]

# CHSP_fv - feed ----
chsp_fv_bf <- data_cc %>%
  select(time.weeks, level, trend.weeks, feed.fv)

#check & drop missingness
table(chsp_fv_bf$feed.fv, useNA = "always")


chsp_fv_bf %>%
  select(time.weeks, feed.fv) %>%
  gg_miss_fct(time.weeks)+ 
  labs(y = "Variable", x="Week of birth") 

sum(!complete.cases(chsp_fv_bf))

chsp_fv_bf <- chsp_fv_bf %>%
  drop_na()

#create denominator
chsp_fv_bf <- chsp_fv_bf %>%
  group_by(time.weeks) %>%
  mutate(cases.per.week = n())

#create binary outcome
chsp_fv_bf <- chsp_fv_bf %>%
  mutate(feed.bin = ifelse(feed.fv == "Breast", 1, 0))
table(chsp_fv_bf$feed.bin)

#aggregate
chsp_fv_bf <- chsp_fv_bf %>%
  group_by(time.weeks) %>%
  summarise(trend.weeks = mean(trend.weeks), 
            level = mean(level), 
            cases.per.week = mean(cases.per.week),
            feed.bin = sum(feed.bin))

plot(chsp_fv_bf$cases.per.week, type ="n", ylim = c(550,1100), xlab = "Week of birth", 
     ylab = "Observations", xaxt="n")
points(chsp_fv_bf$cases.per.week, cex=0.7)
axis(1, at=0:52*4, labels = TRUE)
title("Breastfeeding at initial visit")



### descriptive statistics
#compute rates (episodes per 1000 births)
chsp_fv_bf$rate <- with(chsp_fv_bf, feed.bin/cases.per.week*100)
#create plot of pre-intervention trend
plot(chsp_fv_bf$rate, type ="n", ylim = c(20,50), xlab = "Week of birth", 
     ylab = "Proportion (%)", xaxt="n")
rect(105, 20, 208, 50, col=grey(0.9), border = F)
points(chsp_fv_bf$rate, cex=0.7)
axis(1, at=0:52*4, labels = TRUE)
title("Breastfeeding at initial visit")

summary(chsp_fv_bf)
summary(chsp_fv_bf$feed.bin[chsp_fv_bf$level==0])
summary(chsp_fv_bf$feed.bin[chsp_fv_bf$level==1])
summary(chsp_fv_bf$rate[chsp_fv_bf$level==0])
summary(chsp_fv_bf$rate[chsp_fv_bf$level==1])
summary(chsp_fv_bf$rate)


### model 1 - unadjusted
model1 <- glm(feed.bin ~ offset(log(cases.per.week))
              + time.weeks + level + trend.weeks, family = poisson, chsp_fv_bf)
summary(model1)  
tab_model(model1, vcov.fun = "CL", vcov.type = "HC0", show.se = TRUE, digits = 4)


res2 <- residuals(model1, type = "deviance")
plot(chsp_fv_bf$time.weeks, res2, 
     ylim = c(-7,7), pch=19, cex=0.7, col=grey(0.6),
     main = "Breastfeeding at initial visit (model 1)", ylab = "Deviance residuals",
     xlab = "Week of birth")
abline(h=0, lty=2, lwd=2)

#check acf & p-acf
acf(res2,
    main = "Breastfeeding at initial visit (model 1)")
pacf(res2,
     main = "Breastfeeding at initial visit (model 1)")

#adjust for seasonality
model2 <- glm(feed.bin ~ offset(log(cases.per.week))
              + time.weeks + level + trend.weeks + harmonic(time.weeks,2,52),
              family = poisson, chsp_fv_bf)

tab_model(model2, vcov.fun = "CL", vcov.type = "HC0", show.se = TRUE, digits = 4)


res2 <- residuals(model2, type = "deviance")
plot(chsp_fv_bf$time.weeks, res2, 
     ylim = c(-7,7), pch=19, cex=0.7, col=grey(0.6),
     main = "Breastfeeding at initial visit (model 2)", ylab = "Deviance residuals",
     xlab = "Week of birth")
abline(h=0, lty=2, lwd=2)

#check acf & p-acf
acf(res2,
    main = "Breastfeeding at initial visit (model 2)")
pacf(res2,
     main = "Breastfeeding at initial visit (model 2)")

#model 2 predictions
datanew <- data.frame(cases.per.week=mean(chsp_fv_bf$cases.per.week),
                      level=rep(c(0,1), c(104,103)),
                      time.weeks=1:207)
datanew <- cbind(datanew,trend_pred)
chsp_fv_bf$factual=predict(model2, type="response", datanew)/mean(chsp_fv_bf$cases.per.week)*100

datanew_b <- data.frame(cases.per.week=mean(chsp_fv_bf$cases.per.week),
                        level=0,
                        trend.weeks=0,
                        time.weeks=1:207)
chsp_fv_bf$c.factual=predict(model2, type="response", datanew_b)/mean(chsp_fv_bf$cases.per.week)*100

#absolute effects for model 2
#at 1 month
chsp_fv_bf$factual[109]-chsp_fv_bf$c.factual[109]
#at 6 months
chsp_fv_bf$factual[129]-chsp_fv_bf$c.factual[129]

# CHSP_68 - feed ----
chsp_68_bf <- data_cc %>%
  select(time.weeks, level, trend.weeks, feed.68, simd.quintile.2016, mother.age, total.previous.pregnancies, tpp_bin)

#check & drop missingness
table(chsp_68_bf$feed.68, useNA = "always")


#create binary indicating missingness in feed.68
chsp_68_bf <- chsp_68_bf %>%
  mutate(miss_bin = as.integer(complete.cases(chsp_68_bf$feed.68)))
table(chsp_68_bf$miss_bin)

chsp_68_bf %>% 
  select(feed.68, time.weeks) %>%
  gg_miss_fct(time.weeks) + labs(y = "Variable", x="Week of birth") 


setwd("Results")
feed_68_miss_simd <- CrossTable(chsp_68_bf$simd.quintile.2016, chsp_68_bf$miss_bin, missing.include = TRUE)
write.csv(feed_68_miss_simd, "feed_68_miss_simd.csv")

feed_68_miss_age <- CrossTable(chsp_68_bf$mother.age, chsp_68_bf$miss_bin, missing.include = TRUE)
write.csv(feed_68_miss_age, "feed_68_miss_age.csv")

feed_68_miss_npp <- CrossTable(chsp_68_bf$tpp_bin, chsp_68_bf$miss_bin, missing.include = TRUE)


sum(!complete.cases(chsp_68_bf))
chsp_68_bf <- chsp_68_bf %>%
  drop_na()

#create denominator
chsp_68_bf <- chsp_68_bf %>%
  group_by(time.weeks) %>%
  mutate(cases.per.week = n())

#create binary outcome
chsp_68_bf <- chsp_68_bf %>%
  mutate(feed.bin = ifelse(feed.68 == "Breast", 1, 0))
table(chsp_68_bf$feed.bin)

#aggregate
chsp_68_bf <- chsp_68_bf %>%
  group_by(time.weeks) %>%
  summarise(trend.weeks = mean(trend.weeks), 
            level = mean(level), 
            cases.per.week = mean(cases.per.week),
            feed.bin = sum(feed.bin))

plot(chsp_68_bf$cases.per.week, type ="n", ylim = c(500,1000), xlab = "Week of birth", 
     ylab = "Observations", xaxt="n")
points(chsp_68_bf$cases.per.week, cex=0.7)
axis(1, at=0:52*4, labels = TRUE)
title("Breastfeeding at 6-8 weeks")

### descriptive statistics
#compute prevalence (episodes per 1000 births)
chsp_68_bf$rate <- with(chsp_68_bf, feed.bin/cases.per.week*100)
#create plot of pre-intervention trend
plot(chsp_68_bf$rate, type ="n", ylim = c(0,60), xlab = "Week of birth", 
     ylab = "Prevalence", xaxt="n")
points(chsp_68_bf$rate, cex=0.7)
axis(1, at=0:52*4, labels = TRUE)
title("Breastfeeding at 6-8 weeks")

summary(chsp_68_bf)
summary(chsp_68_bf$feed.bin[chsp_68_bf$level==0])
summary(chsp_68_bf$feed.bin[chsp_68_bf$level==1])
summary(chsp_68_bf$rate[chsp_68_bf$level==0])
summary(chsp_68_bf$rate[chsp_68_bf$level==1])
summary(chsp_68_bf$rate)

### model 1 - unadjusted
model1 <- glm(feed.bin ~ offset(log(cases.per.week))
              + time.weeks + level + trend.weeks, family = poisson, chsp_68_bf)
summary(model1)  
summary(model1)$dispersion  

tab_model(model1, vcov.fun = "CL", vcov.type = "HC0", show.se = TRUE, digits = 4)


res1 <- residuals(model1, type = "deviance")
plot(chsp_68_bf$time.weeks, res1, 
     ylim = c(-7,7), pch=19, cex=0.7, col=grey(0.6),
     main = "Breastfeeding at 6-8 weeks (model 2)", ylab = "deviance residuals",
     xlab = "Week of birth")
abline(h=0, lty=2, lwd=2)

#check acf & p-acf
acf(res1,
    main = "Breastfeeding at 6-8 weeks (model 1)")
pacf(res1,
     main = "Breastfeeding at 6-8 weeks (model 1)")


#adjust for seasonality
model2 <- glm(feed.bin ~ offset(log(cases.per.week))
              + time.weeks + level + trend.weeks + harmonic(time.weeks,2,52),
              family = poisson, chsp_68_bf)

tab_model(model2, vcov.fun = "CL", vcov.type = "HC0", show.se = TRUE, digits = 4)




par(mfrow = c(2,2))

res2 <- residuals(model2, type = "deviance")
plot(chsp_68_bf$time.weeks, res2, 
     ylim = c(-7,7), pch=19, cex=0.7, col=grey(0.6),
     main = "Breastfeeding at 6-8 weeks (model 2)", ylab = "Deviance residuals",
     xlab = "Week of birth")
abline(h=0, lty=2, lwd=2)

#check acf & p-acf
acf(res2,
    main = "Breastfeeding at 6-8 weeks (model 2)")
pacf(res2,
     main = "Breastfeeding at 6-8 weeks (model 2)")

#model 2 predictions
datanew <- data.frame(cases.per.week=mean(chsp_68_bf$cases.per.week),
                      level=rep(c(0,1), c(104,103)),
                      time.weeks=1:207)
datanew <- cbind(datanew,trend_pred)
chsp_68_bf$factual=predict(model2, type="response", datanew)/mean(chsp_68_bf$cases.per.week)*100

datanew_b <- data.frame(cases.per.week=mean(chsp_68_bf$cases.per.week),
                        level=0,
                        trend.weeks=0,
                        time.weeks=1:207)
chsp_68_bf$c.factual=predict(model2, type="response", datanew_b)/mean(chsp_68_bf$cases.per.week)*100

#absolute effects for model 2
#at 1 month
chsp_68_bf$factual[109]-chsp_68_bf$c.factual[109]
#at 6 months
chsp_68_bf$factual[129]-chsp_68_bf$c.factual[129]


# CHSP_68 - sleep.bin.su.su ----
chsp_68_sleep <- data_cc %>%
  select(time.weeks, level, trend.weeks, sleep.bin.su)

#missing cases across all sleep.bin.su variables are dropped
table(chsp_68_sleep$sleep.bin.su, useNA = "always")

chsp_68_sleep %>%
  select(time.weeks,sleep.bin.su) %>%
  gg_miss_fct(time.weeks)+
  labs(x="Week of birth", y="Variable")

chsp_68_sleep <- chsp_68_sleep %>%
  drop_na()


#create denominator
chsp_68_sleep <- chsp_68_sleep %>%
  group_by(time.weeks) %>%
  mutate(cases.per.week = n())


#aggregate
chsp_68_sleep <- chsp_68_sleep %>%
  group_by(time.weeks) %>%
  summarise(trend.weeks = mean(trend.weeks), 
            level = mean(level), 
            cases.per.week = mean(cases.per.week),
            sleep.bin.su = sum(sleep.bin.su))

plot(chsp_68_sleep$cases.per.week, type ="n", ylim = c(600,900), xlab = "Week of birth", 
     ylab = "Observations", xaxt="n")
points(chsp_68_sleep$cases.per.week, cex=0.7)
axis(1, at=0:52*4, labels = TRUE)
title("sleep supine only at 6-8 weeks")

### descriptive statistics su
#compute prevalence
chsp_68_sleep$rate <- with(chsp_68_sleep, sleep.bin.su/cases.per.week*100)
#create plot of pre-intervention trend
plot(chsp_68_sleep$rate, type ="n", ylim = c(85,100), xlab = "Week of birth", 
     ylab = "Proportion (%)", xaxt="n")
rect(105, 85, 208, 100, col=grey(0.9), border = F)
points(chsp_68_sleep$rate, cex=0.7)
axis(1, at=0:52*4, labels = TRUE)
title("sleep supine only at 6-8 weeks")

summary(chsp_68_sleep)
summary(chsp_68_sleep$sleep.bin.su[chsp_68_sleep$level==0])
summary(chsp_68_sleep$sleep.bin.su[chsp_68_sleep$level==1])
summary(chsp_68_sleep$rate[chsp_68_sleep$level==0])
summary(chsp_68_sleep$rate[chsp_68_sleep$level==1])
summary(chsp_68_sleep$rate)


### model 1 su - unadjusted
model1 <- glm(sleep.bin.su ~ offset(log(cases.per.week))
              + time.weeks + level + trend.weeks, family = poisson, chsp_68_sleep)
summary(model1)
tab_model(model1, vcov.fun = "CL", vcov.type = "HC0", show.se = TRUE, digits = 4)


res2 <- residuals(model1, type = "deviance")
plot(chsp_68_sleep$time.weeks, res2, 
     ylim = c(-2,2), pch=19, cex=0.7, col=grey(0.6),
     main = "sleep.bin.su supine only at 6-8 weeks (model 1)", ylab = "Deviance residuals",
     xlab = "Week of birth")
abline(h=0, lty=2, lwd=2)

#check acf & p-acf
acf(res2,
    main = "sleep.bin.su supine only at 6-8 weeks (model 1)")
pacf(res2,
     main = "sleep.bin.su supine only at 6-8 weeks (model 1)")


#model 1 predictions
datanew <- data.frame(cases.per.week=mean(chsp_68_sleep$cases.per.week),
                      level=rep(c(0,1), c(104,103)),
                      time.weeks=1:207)
datanew <- cbind(datanew,trend_pred)
chsp_68_sleep$factual=predict(model1, type="response", datanew)/mean(chsp_68_sleep$cases.per.week)*100

datanew_b <- data.frame(cases.per.week=mean(chsp_68_sleep$cases.per.week),
                        level=0,
                        trend.weeks=0,
                        time.weeks=1:207)
chsp_68_sleep$c.factual=predict(model1, type="response", datanew_b)/mean(chsp_68_sleep$cases.per.week)*100

#absolute effects for model 2
#at 1 month
chsp_68_sleep$factual[109]-chsp_68_sleep$c.factual[109]
#at 6 months
chsp_68_sleep$factual[129]-chsp_68_sleep$c.factual[129]



