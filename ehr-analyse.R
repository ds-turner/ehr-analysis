###########################################################################
## ehr-analyse.R
## 
## Analysis of data extracted using the script 'ehr-clean.R'.
## Inital descriptive analysis then fits a multivaritae Cox model to 
## estimate the hazard ratio between treatment groups.
##
## Author: David Turner
###########################################################################

# Load libraries ----------------------------------------------------------

library(tableone)
library(sjPlot)
library(survival)
library(KMunicate)
library(broom)
library(survminer)
library(cobalt)
library(survival)
library(tidyverse)

# Plot theme --------------------------------------------------------------

# set the theme for the plots
theme_update(
  axis.title = element_text(),
  plot.caption = element_text(hjust = 0, vjust = 0),
  plot.background = element_rect(fill = "white", colour = "white"),
  panel.background = element_rect(fill = "white", colour = "white"),
  legend.background = element_rect(fill = "white", colour = "white"),
  legend.box.background = element_rect(fill = "white", colour = "white"),
  plot.title = element_text(),
  plot.subtitle = element_text()
)

cbf_pal  <- c("#E69F00", # orange
              "#009E73", # bluishgreen
              "#56B4E9", # yellow
              "#0072B2", # blue
              "#D55E00", # vermillion
              "#CC79A7") # reddishpurple

pal.538 <- c("#30a2da",
             "#fc4f30",
             "#e5ae38",
             "#6d904f",
             "#8b8b8b")


# Load analysis data ------------------------------------------------------

analysis_data <- readRDS("analysis_data.rds")

# Descriptive statistics and preliminary analysis -------------------------





## Follow-up duration ------------------------------------------------------

analysis_data %>%
  group_by(ppi) %>%
  summarise(min_followup = min(followup_dur),
            max_followup = max(followup_dur))

## Patients who have an follow-up duration of 0 are not to have been considered 
## at risk in the context of survival analysis and are therefore uninformative.  
## As there patients were prescribed PPI or H2RA on the same day as the enddate 
## we can assume that they were at risk of an amount of time less than 1 day.
## therefore I will add 0.001 (less than half a day) to the folowup_dur while 
## retaining the recorded enddate. 

## for followup_dur of 0, add 0.001 years to followup_dur.
analysis_data <- analysis_data %>%
  mutate(followup_dur = if_else(
    indexdate == enddate, followup_dur + 0.001, followup_dur))

## check
analysis_data %>%
  group_by(ppi) %>%
  summarise(min_followup = min(followup_dur),
            max_followup = max(followup_dur))


# Covariate balance -------------------------------------------------------


# checking the covariate balance between treatment groups.

bal.tab(analysis_data, treat = analysis_data$ppi)

# list of covariates
covariates <- c(
  "age",
  "gender",
  "eth5",
  "imd_person",
  "calendarperiod",
  "bmi",
  "nconsult",
  "prior_diabetes",
  "prior_gastric_cancer",
  "recent_gerd",
  "recent_peptic_ulcer"
)

# list of balance plots
balance_plots <- lapply(covariates, function(x)
  bal.plot(analysis_data, 
           var.nam = x,
           treat = analysis_data$ppi))
# name list
names(balance_plots) <- covariates

# View plots
plot_grid(balance_plots)

# Plots
p <- ggplot(analysis_data)

## density of followup duration by treatment
p +
  geom_density(aes(followup_dur,
                   colour = factor(ppi))) +
  labs(title = "Follow duration by treatment group")

## histogram of followup duration by treatment
p +
  geom_histogram(aes(followup_dur,
                     fill = factor(ppi), 
                     alpha = 0.05), 
                 position = 'identity') +
  labs(title = "Follow duration by treatment group")

## density of followup duration by outcome
p +
  geom_density(aes(followup_dur,
                   colour = factor(died))) +
  labs(title = "Follow duration by treatment outcome")

## Patient characteristics -------------------------------------------------


## Tables
varstofactor <- c(
  "ppi",
  "gender",
  "eth5",
  "imd_person",
  "calendarperiod",
  "died",
  "prior_diabetes",
  "prior_gastric_cancer",
  "recent_gerd",
  "recent_peptic_ulcer"
)
table1vars <- c(
  "ppi",
  "age",
  "gender",
  "eth5",
  "imd_person",
  "calendarperiod",
  "bmi",
  "died",
  "nconsult",
  "prior_diabetes",
  "prior_gastric_cancer",
  "recent_gerd",
  "recent_peptic_ulcer",
  "followup_dur"
)
nonnormal <- c("age",
               "bmi",
               "nconsult",
               "followup_dur")
ppi.cov_summary <- CreateTableOne(
  vars = table1vars,
  data = analysis_data,
  strata = c("ppi"),
  factorVars = varstofactor,
  includeNA = TRUE,
  test = FALSE,
  addOverall = T
)
print(ppi.cov_summary, 
      smd = T, 
      nonnormal = nonnormal,
      quote = TRUE, noSpaces = TRUE)


## Kaplan-Meier plot -------------------------------------------------------

# Survival data
## Timescale is time of exposure to outcome.
## No right censored (delayed entry) needed time origin is exposure date
## (indexdate) end of follow-up  is first of date of death, study end or
## patient leave GP (enddate), unit is years.

## right censored survival data (in years)
surv_data <- Surv(time = analysis_data$followup_dur,
                  event = analysis_data$died)


# Kaplan-Meier model
ppi.km <- survfit(surv_data ~ factor(ppi, 
                                     levels = c(0, 1), 
                                     labels = c("Histamine H2-receptor antagonists", 
                                                "Proton pump inhibitors")),
                  data = analysis_data)
print(ppi.km)

# Kaplan-Meier plot

time_scale <- seq(0, 25, by = 5)
KMunicate(fit = survfit(surv_data ~ factor(ppi, 
                                           levels = c(0, 1), 
                                           labels = c("Histamine H2-receptor antagonists", 
                                                      "Proton pump inhibitors")),
                        data = analysis_data), 
          time_scale = time_scale, 
          .theme = ggplot2::theme_minimal(),
          .title = "Plot of the Kaplan-Meier estimate comparing survival of patients prescribed proton
pump inhibitors (PPI) or histamine H2-receptor antagonists (H2RA)"
)

# Test Survival Curve Differences - Log rank
survdiff(surv_data ~ ppi,
         data = analysis_data)

# Survival analysis: Cox regression ---------------------------------------

# Using a multivariate cox model with the covariates included in the model

## Initial assessment of the proportional hazard assumption ------------------

## log-log plot to assess if the proportional hazard assumption holds for the exposure

# visual assessment using log-log plot
plot(
  survfit(surv_data ~ ppi,
          data = analysis_data),
  fun = "cloglog",
  xlab = "Time of death (log scale)",
  ylab = "logH(tjx) = log(-log S(t))",
  col = cbf_pal,
  main = "Treatment Group"
)
      
legend(
  x = "topleft",
  box.lwd = 0,
  legend = c("H2RA", "PPI"),
  col = cbf_pal,
  lty = 1,
  cex = 0.75,
  inset = 0.01,
  bty = "n"
)
      # Visually the line for each treatment group are almost parallel, 
      # providing some evidence that the ratio of the hazards between the two 
      # treatment groups is constant over time and therefore the proportional 
      # hazard assumption holds for the treatment variable. 

### Assessing correct form of continuous variable --------------------------
## Martingale residual plots
ppi.cox.lin_assesment <- coxph(
  surv_data ~ age +
    log(age) +
    I(age ^ 2) +
    I(age ^ 3) +
    nconsult +
    log(nconsult) +
    I(nconsult ^ 2) +
    I(nconsult ^ 3),
  data = analysis_data
)

ggcoxfunctional(
  ppi.cox.lin_assesment,
  data = analysis_data,
  point.col = "blue",
  point.alpha = 0.3,
  title = "Martingale Residuals"
)

      # untransformed forms of age and nconsult show signs of linearity so will 
      # be included in the model

## for BMI.  
analysis_data.bmi.cc <- drop_na(analysis_data, bmi)

ppi.cox.lin_assesment.bmi <- coxph(Surv(time = analysis_data.bmi.cc$followup_dur,
                                        event = analysis_data.bmi.cc$died) ~ bmi + 
                                     log(bmi) +
                                     I(bmi^2) +
                                     I(bmi^3),
                                   data = analysis_data.bmi.cc)
#
ggcoxfunctional(
  ppi.cox.lin_assesment.bmi,
  data = analysis_data.bmi.cc,
  point.col = "blue",
  point.alpha = 0.3,
  title = "Martingale Residuals"
)
# BMI did not show evidence of linearity untransformed or transformed so 
# grouping BMI into factors that follow the standard groups of underweight, 
# normal, overweight ans obese.

analysis_data <- analysis_data %>%
  mutate(bmi_grp = case_when(
    bmi < 18.5 ~ "< 18.5",
    bmi >= 18.5 & bmi < 25 ~ "18.5 - 25",
    bmi >= 25 & bmi < 30 ~ "25 - 30",
    bmi >= 30  ~ "> 30"),
    bmi_grp = factor(bmi_grp, levels = c("< 18.5",
                                         "18.5 - 25",
                                         "25 - 30",
                                         "> 30"),
                     ordered = TRUE))

### Potential Confounder analysis --------------------------------------------

# Potential confounders assessed by estimating the association of the variate to 
# the risk and the outcome.  No scientific judgment was made as part of this mock
# study.

# list of potential confounders
covariates <- c(
  "age",
  "gender",
  "eth5",
  "imd_person",
  "calendarperiod",
  "bmi_grp",
  "nconsult",
  "prior_diabetes",
  "prior_gastric_cancer",
  "recent_gerd",
  "recent_peptic_ulcer"
)
# List of the association of each covariate to ppi using logistic regression
ppi.association <- lapply(covariates, function(x)
  glm(
    as.formula(paste("ppi", x, sep = "~")),
    family = binomial(link = "logit"),
    data = analysis_data
  ))

# name list
names(ppi.association) <- covariates

# tidy the results
tidy.ppi.association <- lapply(covariates, function(x)
  tidy(
    ppi.association[[x]],
    conf.int = TRUE,
    exponentiate = TRUE
  ))

# name list
names(tidy.ppi.association) <- covariates

# List of the association of each covariate to death using logistic regression
died.association <- lapply(covariates, function(x)
  glm(
    as.formula(paste("died", x, sep = "~")),
    family = binomial(link = "logit"),
    data = analysis_data
  ))
# name list
names(died.association) <- covariates
# tidy the results
tidy.died.association <-lapply(covariates, function(x)
  tidy(
    died.association[[x]],
    conf.int = TRUE,
    exponentiate = TRUE
  ))
# name list
names(tidy.died.association) <- covariates

## likelihood ratio test

## likelihood ratio test to obtain a p-value testing the hypothesis that each 
## covariate is not associated with the probability of receiving a PPI prescription.

lr.ppi <- lapply(covariates, function(x)
  lrtest(
    ppi.association[[x]],
    glm(ppi ~ 1, data = drop_na(analysis_data, x), family = "binomial")
  ))

names(lr.ppi) <- covariates

lr.ppi

## likelihood ratio test to obtain a p-value testing the hypothesis that each 
## covariate is not associated with the probability of death.

lr.died <- lapply(covariates, function(x)
  lrtest(
    died.association[[x]],
    glm(died ~ 1, data = drop_na(analysis_data, x), family = "binomial")
  ))

names(lr.died) <- covariates


## bind values into a dataframe
# ppi association results
lr.ppi_Chisq <-
  sapply(covariates, function(x)
    lr.ppi[[x]]$Chisq[2])
lr.ppi_p_value <-
  sapply(covariates, function(x)
    lr.ppi[[x]]$`Pr(>Chisq)`[2])

# death association results
lr.died_Chisq <-
  sapply(covariates, function(x)
    lr.died[[x]]$Chisq[2])
lr.died_p_value <-
  sapply(covariates, function(x)
    lr.died[[x]]$`Pr(>Chisq)`[2])

# bind to data frame
lr.values <- cbind(as.data.frame(lr.ppi_Chisq),
                   as.data.frame(lr.ppi_p_value),
                   as.data.frame(lr.died_Chisq),
                   as.data.frame(lr.died_p_value)) %>%
  arrange(lr.ppi_Chisq) %>%
  round(4) 

# Recent GERD and ethnicity show lack of evidence of association of either PPI 
# prescription or of all cause mortality.  IMD shows lack of evidence of 
# association of PPI prescriptions.  All 3 variables excluded from the analysis.
## Fit One -----------------------------------------------------------------

## drop record's with missing BMI to perform complete case analysis
analysis_data.cc <- drop_na(analysis_data, bmi)

## survival time for complete cases
surv_data.2 <- Surv(time = analysis_data.cc$followup_dur,
                    event = analysis_data.cc$died)

## fit cox model
ppi.cox.multivarable.1 <- coxph(surv_data.2 ~ ppi +
                                  age +
                                  gender +
                                  calendarperiod +
                                  bmi_grp +
                                  nconsult +
                                  prior_diabetes +
                                  prior_gastric_cancer +
                                  recent_peptic_ulcer,
                                data = analysis_data.cc)

## tidy to view results
tidy(ppi.cox.multivarable.1,
     conf.int = TRUE,
     exponentiate = TRUE)

### Fit 1 Predicted survival curve ------
newdata.1 <- with(
  analysis_data,
  data.frame(
    ppi = c(0, 1),
    age = rep(mean(age, na.rm = TRUE), 2),
    bmi_grp = "18.5 - 25",
    gender = "Female",
    imd_person = "Least Deprived (1)",
    calendarperiod = "1991-1999",
    nconsult = rep(mean(nconsult, na.rm = TRUE), 2),
    prior_diabetes = 0,
    prior_gastric_cancer = 0,
    recent_peptic_ulcer = 0
  )
)

fit.1 <- survfit(ppi.cox.multivarable.1,
                 newdata = newdata.1)

## Plot predicted survival curve 
ggsurvplot(
  fit.1,
  data = analysis_data,
  censor = F,
  conf.int = T,
  #surv.median.line = "hv",
  #risk.table = "abs_pct",
  title = "Cox estimate of the survival curves of patients prescribed proton
pump inhibitors (PPI) and histamine H2-receptor antagonists (H2RA)",
  font.title = 20,
  legend.title = "",
  legend.labs = c("H2RA", "PPI"),
  legend = c(0.8, 0.8)
)

#### Fit 1 Assumptions checking #######################

#### Fit 1 Assessing linearity of covariates --------
## Martingale residual plots
par(mfrow = c(2,2), oma=c(0,0,2,0))
# Age
plot(x = analysis_data.cc$age,
     y = ppi.cox.multivarable.1$residuals,
     xlab = "Age",
     ylab = "Martingale residuals",
     main = "Age")
lines(lowess(analysis_data.cc$age, ppi.cox.multivarable.1$residuals),
      col = pal.538[2])
# Age zoomed in
plot(x = analysis_data.cc$age,
     y = ppi.cox.multivarable.1$residuals,
     xlab = "Age",
     ylab = "Martingale residuals",
     main = "Age zoomed in",
     ylim = c(-1,1))
lines(lowess(analysis_data.cc$age, ppi.cox.multivarable.1$residuals),
      col = pal.538[2])
# nconsult
plot(x = analysis_data.cc$nconsult,
     y = ppi.cox.multivarable.1$residuals,
     xlab = "Number of consultations in last year",
     ylab = "Martingale residuals",
     main = "Number of consultations in last year")

lines(lowess(analysis_data.cc$nconsult, ppi.cox.multivarable.1$residuals),
      col = pal.538[2])
# nconsult zoomed in
plot(x = analysis_data.cc$nconsult,
     y = ppi.cox.multivarable.1$residuals,
     xlab = "Number of consultations in last year",
     ylab = "Martingale residuals",
     main = "Number of consultations in last year zoomed in",
     ylim = c(-1,1))
lines(lowess(analysis_data.cc$nconsult, ppi.cox.multivarable.1$residuals),
      col = pal.538[2])
# Plot title
mtext("Plots of Martingale residuals against covariates", 
      line=0, 
      #side=3, 
      outer=TRUE, 
      cex=2,
      adj = 0.5,
      padj = 0.5)

#### Fit 1 Proportional hazard assumption checking ---------------
## Schoenfeld plot

# proportional hazards assumption test
ph.fit.1 <- cox.zph(ppi.cox.multivarable.1, 
                    transform = "km") 

# plot the results
par(mfrow = c(3,3), oma = c(0, 0, 2, 0))
plot(
  ph.fit.1,
  resid = T,
  col = pal.538,
  ylab = c(
    "Scaled Schoenfeld residuals: ppi",
    "Scaled Schoenfeld residuals: age",
    "Scaled Schoenfeld residuals: gender",
    "Scaled Schoenfeld residuals: calendarperiod",
    "Scaled Schoenfeld residuals: BMI group",
    "Scaled Schoenfeld residuals: nconsult",
    "Scaled Schoenfeld residuals: prior_diabetes",
    "Scaled Schoenfeld residuals: prior_gastric_cancer",
    "Scaled Schoenfeld residuals: recent_peptic_ulcer"
  )
)
mtext(
  "Fit 1: Scaled Schoenfeld residuals for covariates against time",
  line = 0,
  #side=3,
  outer = TRUE,
  cex = 2,
  adj = 0.5,
  padj = 0.5
)

# The calendar period date showed some visual evidence of violating the 
# proportional hazard assumption in the Schoenfeld Residuals test and plots.  

### Multivariable Fit 2 ----------

# Stratified by the calendar period to avoid violation of the 
# proportional hazard assumption.
ppi.cox.multivarable.2 <- coxph(surv_data.2 ~ ppi +
                                  age +
                                  gender +
                                  strata(calendarperiod) +
                                  bmi_grp +
                                  nconsult +
                                  prior_diabetes +
                                  prior_gastric_cancer +
                                  recent_peptic_ulcer,
                                data = analysis_data.cc
)
tidy(ppi.cox.multivarable.2,
     conf.int = TRUE,
     exponentiate = TRUE)
summary(ppi.cox.multivarable.1)
summary(ppi.cox.multivarable.2)

#### Fit 2 Proportional hazard assumption checking ---------------
## Schoenfeld plot
ph.fit.2 <- cox.zph(ppi.cox.multivarable.2, 
                    transform = "km",
                    terms = T,
                    singledf = F)
ph.fit.2
#
ggcoxzph(
  ph.fit.2,
  resid = T,
  se = F,
  point.col = cbf_pal[1],
  point.size = 0.05,
  point.shape = 19,
  point.alpha = 0.1,
)
#
par(mfrow = c(4, 2)
    ,oma = c(0, 0, 2, 0)
)
plot(
  ph.fit.2,
  resid = T,
  col = pal.538,
  ylab = c(
    "Scaled Schoenfeld residuals: ppi",
    "Scaled Schoenfeld residuals: age",
    "Scaled Schoenfeld residuals: gender",
    "Scaled Schoenfeld residuals: BMI group",
    "Scaled Schoenfeld residuals: nconsult",
    "Scaled Schoenfeld residuals: prior_diabetes",
    "Scaled Schoenfeld residuals: prior_gastric_cancer",
    "Scaled Schoenfeld residuals: recent_peptic_ulcer"
  )
)
mtext(
  "Fit 2: Scaled Schoenfeld residuals for covariates against time",
  line = 0,
  #side=3,
  outer = TRUE,
  cex = 2,
  adj = 0.5,
  padj = 0.5
)

#### Fit 2 Predicted survival curve ------

# data for a stratified plot
newdata.2 <- with(
  analysis_data,
  data.frame(
    ppi = rep(c(0, 1), 5),
    age = mean(age, na.rm = TRUE),
    gender = "Female",
    bmi_grp = "18.5 - 25",
    imd_person = "Least Deprived (1)",
    calendarperiod = c("1991-1999", "1991-1999", "2000-2004", "2000-2004", "2005-2009", "2005-2009", "2010-2014","2010-2014", "2015-2017", "2015-2017"),
    nconsult = mean(nconsult, na.rm = TRUE),
    prior_diabetes = 0,
    prior_gastric_cancer = 0,
    recent_peptic_ulcer = 0
  )
)

labs <- expand_grid(year = c("1991-1999", "2000-2004", "2005-2009", "2010-2014","2015-2017"), 
                    trt = c("H2RA", "PPI")) %>%
  mutate(labs = paste0(trt, " - ", year)) %>%
  pull(labs)

fit.2 <- survfit(ppi.cox.multivarable.2, newdata = newdata.2)

# stratified plot
ggsurvplot(fit.2,
           data = analysis_data,
           censor = F,
           conf.int = T,
           #surv.median.line = "hv",
           #risk.table = "abs_pct",
           title = "Stratified Cox estimate of the survival curves of patients prescribed proton
pump inhibitors (PPI) and histamine H2-receptor antagonists (H2RA)",
           font.title = 20,
           legend.title = "",
           legend.labs = labs,
           legend = c(0.15, 0.4)
)

# data for a 'Stratified Cox estimate of the survival curves of Female 
# patients, with a bmi of 25 - 30, prescribed proton pump inhibitors 
# (PPI) and histamine H2-receptor antagonists (H2RA) between 1991 and 
# 1999 with no underlining conditions' plot
newdata.3 <- with(
  analysis_data,
  data.frame(
    ppi = rep(c(0, 1)),
    age = mean(age, na.rm = TRUE),
    gender = "Female",
    bmi_grp = "25 - 30",
    imd_person = "Least Deprived (1)",
    calendarperiod = "1991-1999",
    nconsult = mean(nconsult, na.rm = TRUE),
    prior_diabetes = 0,
    prior_gastric_cancer = 0,
    recent_peptic_ulcer = 0
  )
)

# plot
fit.3 <- survfit(ppi.cox.multivarable.2, newdata = newdata.3)
table(analysis_data$bmi_grp)
ggsurvplot(fit.3,
           data = analysis_data,
           censor = F,
           conf.int = T,
           #surv.median.line = "hv",
           #risk.table = "abs_pct",
           title = 
             "Stratified Cox estimate of the survival curves of Female 
           patients, with a bmi of 25 - 30, prescribed proton pump inhibitors 
           (PPI) and histamine H2-receptor antagonists (H2RA) between 1991 and 
           1999 with no underlining conditions",
           font.title = 20,
           legend.title = "",
           legend.labs = c("H2RA", "PPI"),
           legend = c(0.15, 0.4)
)



