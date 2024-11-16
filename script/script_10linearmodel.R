library(lme4)
library(sjstats)
library(car)
library(sjPlot)
library(arm)
library(texreg)
library(dplyr)
library(RNOmni)
library(ggpubr)
library(rstatix)
library(broom)
library(arsenal)
library(epiDisplay)
library(rcompanion) #needed for Nagelkerke
library(lmtest) #needed for overall p-value of glm

demog <- read.csv("/Users/nagrodzkij/data/angry/output/table_demog_withedu_byccid.csv")
traces_v <- read.csv("/Users/nagrodzkij/data/angry/output/params_v_byparticip.csv")

#demog <- read.csv("/Users/nagrodzkij/data/angry/output/fmri/2ndLevel/main/demog_filtered.csv")
#traces_v <- read.csv("/Users/nagrodzkij/data/angry/output/fmri/2ndLevel/main/params_v_filtered.csv")
#traces_v <- traces_v[-c(which(traces_v$X==520055)),]

demog['v_angry']<-traces_v$v_angry
demog['v_neutral']<-traces_v$v_neutral
demog['t']<-traces_v$t
demog['a']<-traces_v$a

FaceScores <- read.csv("/Users/nagrodzkij/data/angry/input/demog/FaceScores.csv")
FaceScores <- FaceScores[-c(which(FaceScores$ccid==520055)),]
demog['benton']<-FaceScores$benton

reg <- glm(formula = prev_depressed ~ hads_depression + age + v_neutral + v_angry
           ,data = demog, family="binomial" )
summary(reg)

nagelkerke(reg)

anova(reg, 
      update(reg, ~1),    # update here produces null model for comparison 
      test="Chisq")

lrtest(reg)

exp(cbind(OR = coef(reg), confint(reg)))

predict <- predict(reg, demog, type = 'response')
table_mat <- table(demog$prev_depressed, predict > 0.5)
table_mat

accuracy_Test <- sum(diag(table_mat)) / sum(table_mat)
accuracy_Test

# Check if this persists - CONTROL ANALYSIS 

reg <- glm(formula = prev_depressed ~ hads_depression + age + v_neutral + v_angry + sex + benton + highest_qualification 
           + a + t,data = demog, family="binomial" )
summary(reg)


nagelkerke(reg)

anova(reg, 
      update(reg, ~1),    # update here produces null model for comparison 
      test="Chisq")

lrtest(reg)

exp(cbind(OR = coef(reg), confint(reg)))

predict <- predict(reg, demog, type = 'response')
table_mat <- table(demog$prev_depressed, predict > 0.5)
table_mat

accuracy_Test <- sum(diag(table_mat)) / sum(table_mat)
accuracy_Test

# Check v_angry and on_antidepressant + diagnosis_delay

demog_prevdepressed = demog[demog$prev_depressed==1,]

reg <- glm(formula = v_angry ~ hads_depression + age + antidepressant + diagnosisdelay
           ,data = demog, family="gaussian" )
summary(reg)
