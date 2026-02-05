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

#You might need to adjust the # of cores to match those in your machine
#Be wary of running all these at once unless you have time
#Initial model runs were with a reduced number of iterations to compare model fit
m0 <- brm(data=fgI1, formula= Metabolites.Conc ~ 1, chains = 2, iter = 1000, warmup = 200, family = lognormal(), cores=16)
m0re <- brm(data=fgI1, formula= Metabolites.Conc ~ 1 + (1|id), chains = 2, iter = 1000, warmup = 200, family = lognormal(), cores=16)
m0ds <- brm(data=fgI1, formula= Metabolites.Conc ~ s(Days.S) + (1|id), chains = 2, iter = 1000, warmup = 200, family = lognormal(), cores=16, prior = prior(normal(0, 1), class = "b"))
m0dl <- brm(data=fgI1, formula= Metabolites.Conc ~ Days.S + (1|id), chains = 2, iter = 1000, warmup = 200, family = lognormal(), cores=16, prior = prior(normal(0, 1), class = "b"))
m0dr <- brm(data=fgI1, formula= Metabolites.Conc ~ 1 + (1|Days.S) + (1|id), chains = 2, iter = 1000, warmup = 200, family = lognormal(), cores=16)
m0dg <- brm(data=fgI1, formula= Metabolites.Conc ~ Group.1 + Days.S + (1|id), chains = 2, iter = 1000, warmup = 200, family = lognormal(), cores=16, prior = prior(normal(0, 1), class = "b"))
m0rg <- brm(data=fgI1, formula= Metabolites.Conc ~ Group.1 + (1|Days.S) + (1|id), chains = 2, iter = 1000, warmup = 200, family = lognormal(), cores=16, prior = prior(normal(0, 1), class = "b"))
m0dgi <- brm(data=fgI1, formula= Metabolites.Conc ~ Group.1 * Days.S + (1|id), chains = 2, iter = 1000, warmup = 200, family = lognormal(), cores=16, prior = prior(normal(0, 1), class = "b"))
loo_compare(loo(m0),loo(m0re),loo(m0ds),loo(m0dl),loo(m0dr),loo(m0dg),loo(m0rg),loo(m0dgi))
#move forward with m0dr

#assess temperature
m1xn <- brm(data=fgI1, formula= Metabolites.Conc ~ MAXAvg + MINAvg + (1|Days.S) + (1|id), chains = 2, iter = 1000, warmup = 200, family = lognormal(), cores=16)
m1x <- brm(data=fgI1, formula= Metabolites.Conc ~ MAXAvg + (1|Days.S) + (1|id), chains = 2, iter = 1000, warmup = 200, family = lognormal(), cores=16)
m1n <- brm(data=fgI1, formula= Metabolites.Conc ~ MINAvg + (1|Days.S) + (1|id), chains = 2, iter = 1000, warmup = 200, family = lognormal(), cores=16)
m1xn2 <- brm(data=fgI1, formula= Metabolites.Conc ~ MAXAvg * MINAvg + (1|Days.S) + (1|id), chains = 2, iter = 1000, warmup = 200, family = lognormal(), cores=16)
loo_compare(loo(m1xn),loo(m1x),loo(m1n),loo(m1xn2))
performance::check_collinearity(m1xn)

#These are preliminary models prior to including Coping Style scores
m2s <- brm(data=fgI1, formula= Metabolites.Conc ~ MAXAvg + MINAvg + Sex + (1|Days.S) + (1|id), chains = 2, iter = 1000, warmup = 200, family = lognormal(), cores=16)
m2rs <- brm(data=fgI1, formula= Metabolites.Conc ~ MAXAvg + MINAvg + Sex + RankProp + (1|Days.S) + (1|id), chains = 2, iter = 1000, warmup = 200, family = lognormal(), cores=16)
m2rsi <- brm(data=fgI1, formula= Metabolites.Conc ~ MAXAvg + MINAvg + Sex * RankProp + (1|Days.S) + (1|id), chains = 2, iter = 1000, warmup = 200, family = lognormal(), cores=16)
loo_compare(loo(m2s),loo(m2rs),loo(m2rsi))

#Candidate models that include Coping Style scores
m3 <- brm(data=fgI1, formula= Metabolites.Conc ~ MAXAvg + MINAvg + Sex + SWDIgr + Coping + (1|Days.S) + (1|id), chains = 2, iter = 1000, warmup = 200, family = lognormal(), cores=16)
m3sc <- brm(data=fgI1, formula= Metabolites.Conc ~ MAXAvg + MINAvg + Sex + SWDIgr * Coping + (1|Days.S) + (1|id), chains = 2, iter = 1000, warmup = 200, family = lognormal(), cores=16)
m3xs <- brm(data=fgI1, formula= Metabolites.Conc ~ MAXAvg + MINAvg + Sex * SWDIgr + Coping + (1|Days.S) + (1|id), chains = 2, iter = 1000, warmup = 200, family = lognormal(), cores=16)
m3xc <- brm(data=fgI1, formula= Metabolites.Conc ~ MAXAvg + MINAvg + Sex * Coping + SWDIgr + (1|Days.S) + (1|id), chains = 2, iter = 1000, warmup = 200, family = lognormal(), cores=16)
m3xsc <- brm(data=fgI1, formula= Metabolites.Conc ~ MAXAvg + MINAvg + Sex * Coping * SWDIgr + (1|Days.S) + (1|id), chains = 2, iter = 1000, warmup = 200, family = lognormal(), cores=16)
loo_compare(loo(m3),loo(m3sc),loo(m3xs),loo(m3xc),loo(m3xsc))

#Add run date to see if it improves fit
m3xsd <- brm(data=fgI1, formula= Metabolites.Conc ~ MAXAvg + MINAvg + Sex * SWDIgr + Coping + (1|RunDate) + (1|Days.S) + (1|id), chains = 2, iter = 1000, warmup = 200, family = lognormal(), cores=16)
loo_compare(loo(m3xsd),loo(m3xs))
RL <- loo(m3xs, reloo = TRUE)
loo_compare(RL, loo(m3xs))

#Run more imputations before proceeding with the final models
fgI100 <- mice(fg, pred = pred1, method="pmm", m = 100, maxit = 10, seed = 5678)

#Originally computed this, but easier to run posterior checks without combining
#Retaining this line for posterity, but I used the M3_100_noC model
#M3_100 <- brm_multiple(data=fgI100, formula= Metabolites.Conc ~ MAXAvg + MINAvg + Sex * SWDIgr + Coping + (1|Days.S) + (1|id), chains = 4, iter = 3000, warmup = 1000, seed = 321, family = lognormal(), cores=16, thin = 2,  prior = prior(normal(0, 1), class = "b"), sample_prior = "yes")

#This is the reported model, which allows checking and bayesR2
M3_100_noC <- brm_multiple(data=fgI100, formula= Metabolites.Conc ~ MAXAvg + MINAvg + Sex * SWDIgr + Coping + (1|Days.S) + (1|id), combine=FALSE, chains = 4, iter = 3000, warmup = 1000, seed = 321, family = lognormal(), cores=16, thin = 2,  prior = prior(normal(0, 1), class = "b"), sample_prior = "yes")

# Pooled Bayes_R2 estimates obtained with code from Jeffrey Girard # https://github.com/jmgirard
# posted Sep 15, 2020
# https://github.com/paul-buerkner/brms/issues/997#issuecomment-692817166
pool_R2 <- function(mlist, probs = c(0.025, 0.975), robust = FALSE, ...) {
  r2_post <- sapply(mlist, bayes_R2, summary = FALSE, ..., simplify = TRUE)
  posterior_summary(c(r2_post), probs = probs, robust = robust)
}

pool_R2(M3_100_noC)
pool_R2(M3_100_noC, re_formula=NA)

#Plot model fit for a subset of the imputed models
set.seed(1241) 
seq <- sample(1:100,9, replace=FALSE)
a <- pp_check(M3_100_noC[[seq[1]]], ndraws=25)
b <- pp_check(M3_100_noC[[seq[2]]], ndraws=25)
c <- pp_check(M3_100_noC[[seq[3]]], ndraws=25)
d <- pp_check(M3_100_noC[[seq[4]]], ndraws=25)
e <- pp_check(M3_100_noC[[seq[5]]], ndraws=25)
f <- pp_check(M3_100_noC[[seq[6]]], ndraws=25)
g <- pp_check(M3_100_noC[[seq[7]]], ndraws=25)
h <- pp_check(M3_100_noC[[seq[8]]], ndraws=25)
i <- pp_check(M3_100_noC[[seq[9]]], ndraws=25)

ggpubr::ggarrange(a,b,c,d,e,f,g,h,i, nrow=3, ncol=3, legend="none")

seq <- sample(1:100,9, replace=FALSE)
a <- pairs(M3_100_noC[[seq[1]]])
b <- pairs(M3_100_noC[[seq[2]]])
c <- pairs(M3_100_noC[[seq[3]]])
d <- pairs(M3_100_noC[[seq[4]]])

ggpubr::ggarrange(a,b,c,d, nrow=2, ncol=2, legend="none")

test <- ggpubr::ggarrange(a,b,c,d, nrow=2, ncol=2)
ggsave("E:/fGCm_Bayes/Dryad/Imput_Pairs.jpg",test, width=8000, height=8000, units="px") 

#Now that I am happy with the basic model, I compared some small variations
#These models use the original (non-imputed) dataset, as Coping is excluded. 
Fin_mod3.2 <- brm(data=fg, formula= Metabolites.Conc ~ MAXAvg + MINAvg + SWDIgr * Sex + (1 | id) + (1|Days.S), chains = 4, iter = 6000, warmup = 1000, seed = 321, family = lognormal(), cores=16, thin = 2,  prior = prior(normal(0, 1), class = "b"), save_pars = save_pars(all = TRUE), sample_prior = "yes")

#Confirm nothing odd based on reference level
fg$Sex <- relevel(fg$Sex, "F")
Fin_mod3.22 <- brm(data=fg, formula= Metabolites.Conc ~ MAXAvg + MINAvg + SWDIgr * Sex + (1 | id) + (1|Days.S), chains = 4, iter = 6000, warmup = 1000, seed = 321, family = lognormal(), cores=16, thin = 2,  prior = prior(normal(0, 1), class = "b"), save_pars = save_pars(all = TRUE), sample_prior = "yes")

#Examine interaction's posteriors
#First, combine the listed models into a single model object
obj <- combine_models(mlist=M3_100_noC, check_data=FALSE)

#hypothesis(M3_100, c("SWDIgr = 0", "SWDIgr + SexF + SexF:SWDIgr = 0 + SexF"), alpha=0.05)
hypothesis(obj, c("SWDIgr = 0", "SWDIgr + SexF + SexF:SWDIgr = 0 + SexF"), alpha=0.05)

#hypothesis(Fin_mod3.2, c("SWDIgr = 0", "SWDIgr + SexF + SWDIgr:SexF = 0 + SexF"), alpha=0.05)
hypothesis(Fin_mod3.22, c("SWDIgr = 0", "SWDIgr + SexM + SWDIgr:SexM = 0 + SexM"), alpha=0.05)

#emmeans::emtrends(M3_100, pairwise ~ Sex, var = "SWDIgr")
emmeans::emtrends(obj, pairwise ~ Sex, var = "SWDIgr")

#What might be driving the effects of SWDI? We included both degree and strength
DegStr <- brm(data=fg, formula= Metabolites.Conc ~ MAXAvg + MINAvg + Deg + Strg + Sex + Strg:Deg + Strg:Sex + Deg:Sex + (1 | id) + (1|Days.S), chains = 4, iter = 6500, warmup = 1500, seed = 321, family = lognormal(), cores=4, thin = 2,  prior = prior(normal(0, 1), class = "b"), save_pars = save_pars(all = TRUE))
performance::check_collinearity(DegStr)
#Some multicollinearity, but not pronounced. This is to be expected with interactions.
#As Strg:Sex is not contributing, can omit

DegStr2 <- brm(data=fg, formula= Metabolites.Conc ~ MAXAvg + MINAvg + Deg + Strg + Sex + Strg:Deg + Deg:Sex + (1 | id) + (1|Days.S), chains = 4, iter = 6500, warmup = 1500, seed = 321, family = lognormal(), cores=4, thin = 2,  prior = prior(normal(0, 1), class = "b"), save_pars = save_pars(all = TRUE))
performance::check_collinearity(DegStr2)
#Resolves multicollinearity; results remain the same

emmeans::emtrends(DegStr, pairwise ~ Sex, var = "Deg")

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
