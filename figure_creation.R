library(ggpubr)
library(ggplot2)
library(brms)
library(bayestestR)
library(bayesplot)
library(dplyr)
library(performance)
library(mice)

fgc <- read.csv('/Users/maddiethall/R Repos/replication_project/Pritchard_fGCm_baboons_R3.csv')

#Convert columns to correct class
fgc$Sex <- factor(fgc$Sex)
fgc$Group.1 <- factor(fgc$Group.1)
fgc$id <- factor(fgc$id)
fgc$Metabolites.Conc <- as.numeric(fgc$Metabolites.Conc)
fgc$MAXAvg <- as.numeric(fgc$MAXAvg)
fgc$SWDIgr <- as.numeric(fgc$SWDIgr)
fgc$MINAvg <- as.numeric(fgc$MINAvg)
fgc$RankProp <- as.numeric(fgc$RankProp)
fgc$Coping <- as.numeric(fgc$Coping)
fgc$Strg <- as.numeric(fgc$Strg)
fgc$Deg <- as.numeric(fgc$Deg)

#Scale and center numeric predictors
columns_to_scale_center <- c("MAXAvg", "MINAvg", "RankProp", "SWDIgr", "Strg", "Deg")  # Specify the column names you want to scale and center

scaled_centered_data <- fgc %>%
  mutate_at(columns_to_scale_center,
            list(~ (. - mean(.)) / (2 * sd(.))))

fg <- scaled_centered_data

#Set reference categories
fg$Sex <- relevel(fg$Sex, "M")
fg$Group.1 <- relevel(fg$Group.1, "S")

#Calculate days as a numeric, relative to date of first fecal sample collection
for(i in 1:nrow(fg)){fg$Days[i] <- as.numeric(difftime(fg$Date[i], min(fg$Date)))}

#Scale days
fg$Days.S <- (fg$Days - mean(fg$Days)) / (2 * sd(fg$Days))

#Scale coping style scores, while excluding 'NA' values, which were not subjected to experiments
fg$Coping[!is.na(fg$Coping)] <- (fg$Coping[!is.na(fg$Coping)] - mean(fg$Coping[!is.na(fg$Coping)])) / (2 * sd(fg$Coping[!is.na(fg$Coping)]))

#Obtain predictor matrix
ini <- mice(fg,maxit=0) 
pred1 <- ini$predictorMatrix 
pred1[,-c(12:13)] <- 0 

#Populate imputed dataframe
#Preliminary models use a smaller imputed dataframe to limit computation time
fgI1 <- mice(fg, pred = pred1, method="pmm", m = 1, maxit = 20, seed = 1234)





#Run more imputations before proceeding with the final models
fgI100 <- mice(fg, pred = pred1, method="pmm", m = 100, maxit = 10, seed = 5678)

#This is the reported model, which allows checking and bayesR2
M3_100_noC <- brm_multiple(data=fgI100, formula= Metabolites.Conc ~ MAXAvg + MINAvg + Sex * SWDIgr + Coping + (1|Days.S) + (1|id), combine=FALSE, chains = 4, iter = 3000, warmup = 1000, seed = 321, family = lognormal(), cores=16, thin = 2,  prior = prior(normal(0, 1), class = "b"), sample_prior = "yes")

obj <- combine_models(mlist=M3_100_noC, check_data=FALSE)





#Populate dataframes with fitted data, fixed to Min, Mean, or Average values

dfxAxR <- data.frame(Sex="F", SWDIgr = max(fg$SWDIgr), Coping = median(fg$Coping[!is.na(fg$Coping)]), MINAvg = median(fg$MINAvg), MAXAvg = median(fg$MAXAvg))
fit.FMx <- as.data.frame(fitted(obj, newdata = dfxAxR, summary = F, re_formula = NA))

dfxAaR <- data.frame(Sex="M", SWDIgr = max(fg$SWDIgr), Coping = median(fg$Coping[!is.na(fg$Coping)]), MINAvg = median(fg$MINAvg), MAXAvg = median(fg$MAXAvg))
fit.MMx <- as.data.frame(fitted(obj, newdata = dfxAaR, summary = F, re_formula = NA))

dfaAxR <- data.frame(Sex="F", SWDIgr = median(fg$SWDIgr), Coping = median(fg$Coping[!is.na(fg$Coping)]), MINAvg = median(fg$MINAvg), MAXAvg = median(fg$MAXAvg))
fit.FMn <- as.data.frame(fitted(obj, newdata = dfaAxR, summary = F, re_formula = NA))

dfaAaR <- data.frame(Sex="M", SWDIgr = median(fg$SWDIgr), Coping = median(fg$Coping[!is.na(fg$Coping)]), MINAvg = median(fg$MINAvg), MAXAvg = median(fg$MAXAvg))
fit.MMn <- as.data.frame(fitted(obj, newdata = dfaAaR, summary = F, re_formula = NA))

dfnAxR <- data.frame(Sex="F", SWDIgr = min(fg$SWDIgr), Coping = median(fg$Coping[!is.na(fg$Coping)]), MINAvg = median(fg$MINAvg), MAXAvg = median(fg$MAXAvg))
fit.FMi <- as.data.frame(fitted(obj, newdata = dfnAxR, summary = F, re_formula = NA))

dfnAxRF <- data.frame(Sex="F", SWDIgr = min(fg$SWDIgr[fg$Sex == "F"]), Coping = median(fg$Coping[!is.na(fg$Coping)]), MINAvg = median(fg$MINAvg), MAXAvg = median(fg$MAXAvg))
fit.FMiF <- as.data.frame(fitted(obj, newdata = dfnAxRF, summary = F, re_formula = NA))

dfnAaR <- data.frame(Sex="M", SWDIgr = min(fg$SWDIgr), Coping = median(fg$Coping[!is.na(fg$Coping)]), MINAvg = median(fg$MINAvg), MAXAvg = median(fg$MAXAvg))
fit.MMi <- as.data.frame(fitted(obj, newdata = dfnAaR, summary = F, re_formula = NA))

dfnAaRF <- data.frame(Sex="M", SWDIgr = min(fg$SWDIgr[fg$Sex == "F"]), Coping = median(fg$Coping[!is.na(fg$Coping)]), MINAvg = median(fg$MINAvg), MAXAvg = median(fg$MAXAvg))
fit.MMiF <- as.data.frame(fitted(obj, newdata = dfnAaRF, summary = F, re_formula = NA))

#Add factor identifiers for sex

fit.FMx$Sex <- "Female"
fit.FMn$Sex <- "Female"
fit.FMi$Sex <- "Female"
fit.FMiF$Sex <- "Female"

fit.MMx$Sex <- "Male"
fit.MMn$Sex <- "Male"
fit.MMi$Sex <- "Male"
fit.MMiF$Sex <- "Male"

#Add factor identifiers for SWDI value
fit.FMx$SWDIgr <- "High"
fit.FMn$SWDIgr <- "Median"
fit.FMi$SWDIgr <- "Low"
fit.FMiF$SWDIgr <- "Low"

fit.MMx$SWDIgr <- "High"
fit.MMn$SWDIgr <- "Median"
fit.MMi$SWDIgr <- "Low"
fit.MMiF$SWDIgr <- "Low"

#Combine to a dataframe with, either, a minimum set to female minimum
Int <- rbind(fit.FMx, fit.FMn, fit.FMiF, fit.MMx, fit.MMn, fit.MMiF)
#Or a minimum for the entire dataset
Int <- rbind(fit.FMx, fit.FMn, fit.FMi, fit.MMx, fit.MMn, fit.MMi)

#Add column identifiers
colnames(Int) <- c("fGCm Concentration", "Sex", "SWDI")

#Identify factor levels
Int$Sex <- factor(Int$Sex, levels = c("Female", "Male"))
Int$SWDI <- factor(Int$SWDI, levels = c("High", "Median", "Low"))

#Plot that matches figure in manuscript and supplementary
ggplot(Int, aes(x=`fGCm Concentration`, y=SWDI)) +
  geom_density_ridges(aes(fill=Sex, alpha=Sex)) + 
  theme_pubclean() +  
  theme(text = element_text(size = 20))  +
  coord_cartesian(xlim = c(min(Int$`fGCm Concentration`), max(Int$`fGCm Concentration`))) +
  scale_fill_manual(values = c("Female" = "purple", "Male"="goldenrod")) + 
  scale_alpha_discrete(range = c(0.75, 0.4)) +
  scale_x_continuous(breaks = seq(round(min(Int$`fGCm Concentration`),-1), round(max(Int$`fGCm Concentration`),-1), by = 25))
