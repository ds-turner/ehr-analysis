###########################################################################
## ehr-clean.R
##
## Define and extract an appropriate analysis dataset from artificial data 
## based on CPRD GOLD format.  
##  Applies the study inclusion/exclusion criteria to include only 
##    patients that meet the criteria.
##  Derives variables to be included in the analysis.
##
## Please see readme file
###########################################################################

# Load libraries ----------------------------------------------------------
library(dplyr)
library(purrr)
library(stringr)
library(forcats)

# Read in data and codelists ----------------------------------------------

## Data
## List paths for the data files
files <-
  list.files("Data",
             full.names = TRUE,
             pattern = ".csv")

## extract the data domain name from the file name
domain_names <- files %>% 
  basename() %>% 
  tools::file_path_sans_ext() %>% 
  str_remove("Sim_") %>% 
  str_to_lower()

## read data into a list of data frames
data <- lapply(files, readr::read_csv)

## name the data frames with the relevant domain name
names(data) <- domain_names

## Codelists
## List paths of the codelists
files <-
  list.files("Codelists",
             full.names = TRUE,
             pattern = ".csv")

## extract the domain name from the file name
domain_names <- files %>% 
  basename() %>% 
  tools::file_path_sans_ext() %>% 
  str_remove("_codes") %>% 
  str_to_lower()

## read codelists into a list of data frames
codes <- lapply(files, readr::read_csv)

## name the data frames with the relevant domain name
names(codes) <- domain_names


# Factor vars -------------------------------------------------------------

## Check class of each var
sapply(data, function(x)
  sapply(x, class))

## Columns to be factored
factor_vars <- c("gender", "eth5", "constype", "imd_person")

## Factor columns in factor_vars across all data.frames in data
data <-
  map(data, ~ .x %>% mutate(across(any_of(factor_vars), as_factor)))

## Order imd_person
data$imd$imd_person <- factor(
  data$imd$imd_person,
  levels = c("Least Deprived (1)",
             "2",
             "3",
             "4",
             "Most Deprived (5)"),
  ordered = TRUE
)


# Identify cohort ---------------------------------------------------------


# filter therapy for all PPI prescriptions, add treatment indicator (ppi) and set to 1.
ppis <- data$therapy %>%
  semi_join(codes$ppi, by = "prodcode") %>%
  mutate(ppi = 1)

# Filter therapy for all H2RA prescriptions and set treatment indicator to 0.
h2ras <- data$therapy %>%
  semi_join(codes$h2ra, by = "prodcode") %>%
  mutate(ppi = 0)

# Retain number of patients extracted
pat_extracted <- nrow(data$patient)

# Build cohort of patients
cohort <-
  rbind(ppis, h2ras) %>% # combined prescriptions with indicators
  arrange(patid, eventdate) %>%
  filter(!duplicated(patid)) %>% # date of first prescription
  left_join(data$patient, by = "patid") 

# check for missing data. 
cohort %>% 
  summarise(across(everything(), ~ sum(is.na(.)))) %>% 
  select(where(~any(. != 0)))
    # tod and death date fine to have NAs present

# check values of ethnicity and gender
count(cohort, eth5) 
count(cohort, gender)

# drop unknown ethnicity
cohort <- cohort %>% 
  filter(eth5 != "Unknown")

# Retain number of patients in cohort.
pat_complete <- nrow(cohort)

# Drop patients who's first prescriptions are not after patient's registration 
# at GP plus one year
cohort <- cohort %>%
  filter(eventdate > crd + 365.25)

# retain cohort size
pat_lookback <- nrow(cohort) 

# Drop patients who's first prescription were not between 1-Jan-1991 and 17-Apr-2017
cohort <- cohort %>%
  filter(eventdate >= as.Date("1991-01-01") & eventdate <= as.Date("2017-04-17")) %>%
  select(patid, eventdate, ppi) %>%
  rename(indexdate = eventdate) # rename date of first prescription to indexdate

# retain cohort size
pat_enrol_period <- nrow(cohort)

# retain treatment group sizes
n_ppi <- sum(cohort$ppi)
n_h2ras <- nrow(cohort) - sum(cohort$ppi)



# Identify the date of the end of the study for each patient --------------

enddates <- cohort %>%
  left_join(data$patient, by = "patid") %>%
  mutate(enddate = pmin(
    #Minimum of ...
    as.Date("2017-04-17"), #...the study end date 
    # ...or ...
    pmin(tod, deathdate, na.rm = TRUE), #...the min of transfer out of practice date or death date.
    na.rm = TRUE
  ),
  died = ifelse(is.na(deathdate), 0, (deathdate == enddate))) %>% # set the date of death
  select(patid, deathdate, enddate, died)


# Identify covariates -----------------------------------------------------

## Demographics
demographics <- cohort %>%
  left_join(data$patient, by = "patid") %>%
  left_join(data$imd, by = "patid") %>%
  mutate( # derive age assuming mid year birth
    age = round(as.numeric(indexdate - as.Date(
      paste0(yob, "-06-15")
    )) / 365.25), 
    calendarperiod = case_when( # split into period indexdate occurred
      indexdate < as.Date("2000-01-01") ~ "1991-1999",
      indexdate < as.Date("2005-01-01") ~ "2000-2004",
      indexdate < as.Date("2010-01-01") ~ "2005-2009",
      indexdate < as.Date("2015-01-01") ~ "2010-2014",
      indexdate >= as.Date("2015-01-01") ~ "2015-2017"
    ), # 
    calendarperiod = factor(
      calendarperiod,
      levels = c(
        "1991-1999",
        "2000-2004",
        "2005-2009",
        "2010-2014",
        "2015-2017"
      ),
      ordered = TRUE
    )
  ) %>%
  rename(pracid = pracid.x) %>%
  select(patid, age, gender, eth5, imd_person, calendarperiod)

## Body Mass Index

# Weight data
weight_data <- data$clinical %>%
  filter(enttype == 13) %>%
  inner_join(cohort, by = "patid") %>%
  filter(between(as.numeric(indexdate - eventdate), 
                 0, 
                 5 * 365.25))   %>% # within 5 years of the index date
  inner_join(data$additional, by = "adid") %>%
  rename(weight_kg = data1,
         patid = patid.x,
         weight_date = eventdate) %>%
  filter(!is.na(weight_kg), weight_kg >= 20)  %>% # excludes unfeasible weights
  group_by(patid, weight_date) %>%
  summarize(weight_kg = mean(weight_kg))  %>%
  arrange(patid, weight_date) %>%
  filter(!duplicated(patid)) %>%
  select(patid, weight_kg, weight_date)

bmi_data <- data$clinical %>%
  filter(enttype == 13) %>%
  inner_join(cohort, by = "patid") %>%
  filter(between(as.numeric(indexdate - eventdate), 
                 0, 
                 5 * 365.25))   %>% # within 5 years of the index date
  inner_join(data$additional, by = "adid") %>%
  rename(bmi = data3,
         patid = patid.x,
         bmi_date = eventdate) %>%
  filter(!is.na(bmi), between(bmi, 5, 200)) %>%
  group_by(patid, bmi_date) %>%
  summarize(bmi = mean(bmi))  %>%
  arrange(patid, bmi_date) %>%
  filter(!duplicated(patid)) %>%
  mutate(preference = 2) %>%
  select(patid, bmi, bmi_date, preference)

# Obtain the patients' year of birth (to be used to restrict heights to post-17 years)
yob <- data$patient %>%
  select(patid, yob)

# Height measurements
height_data <- data$clinical %>%
  filter(enttype == 14)  %>% # select additional data of height
  inner_join(cohort, by = "patid") %>%
  mutate(yoh = as.numeric(format(eventdate, "%Y")))  %>%
  inner_join(yob, by = "patid") %>%
  filter(yoh - yob >= 17)   %>% # remove heights taken at approx under 17 years
  inner_join(data$additional, by = "adid") %>%
  rename(height_m = data1,
         patid = patid.x,
         height_date = eventdate) %>%
  filter(!is.na(height_m), 
         between(height_m, 1.20, 2.15)) %>% # drop missing or unfeasible values
  group_by(patid, height_date) %>%
  summarize(height_m = mean(height_m)) %>% # tame mean of heights
  arrange(patid, height_date) %>%
  filter(!duplicated(patid)) %>%
  select(patid, height_m, height_date)

# Calculate BMI from weight and height
bmi_calculate <-
  inner_join(height_data, weight_data, by = "patid")  %>%
  mutate(bmi = weight_kg / height_m ^ 2, 
         preference = 1) %>% # preference to aid selection of calculated BMIs
  rename(bmi_date = weight_date) %>%
  select(patid, bmi, bmi_date, preference)

# Take calculated BMI if available, recorded BMI if not
bmi <- rbind(bmi_calculate, bmi_data)
bmi <- bmi %>%
  arrange(patid, preference) %>%
  filter(!duplicated(patid)) %>%
  select(patid, bmi)

# Diagnosis of diabetes prior to indexdate
diabetes <- data$clinical %>%
  inner_join(codes$diabetes, by = "medcode") %>%
  arrange(patid, eventdate) %>%
  filter(!duplicated(patid)) %>% # date of first diagnosis
  inner_join(cohort, by = "patid") %>%
  filter(eventdate < indexdate) %>% # drop patients diagnosed after to indexdate
  mutate(prior_diabetes = 1) %>%
  select(patid, prior_diabetes)

## Diagnosis of gastric cancer prior to indexdate
gastric_cancer <- data$clinical %>%
  inner_join(codes$gastric_cancer, by = "medcode") %>%
  arrange(patid, eventdate) %>%
  filter(!duplicated(patid)) %>% # date of first diagnosis
  inner_join(cohort, by = "patid") %>%
  filter(eventdate < indexdate) %>% # drop patients diagnosed after to indexdate
  mutate(prior_gastric_cancer = 1) %>%
  select(patid, prior_gastric_cancer)

## Gastro-oesophageal reflux disease (GERD) in the 6 months prior to indexdate
recent_gerd <- data$clinical %>%
  inner_join(cohort, by = "patid") %>%
  arrange(patid, eventdate) %>%
  inner_join(codes$gerd, by = "medcode") %>%
  # retain patients with diagnosis in the 6 months prior to indexdate
  filter(eventdate >= indexdate - 180 & eventdate < indexdate) %>%
  filter(!duplicated(patid)) %>%
  select(patid) %>%
  mutate(recent_gerd = 1)

## Peptic ulcer in the 6 months prior to prescription
peptic_ulcer <- data$clinical %>%
  inner_join(cohort, by = "patid") %>%
  arrange(patid, eventdate) %>%
  inner_join(codes$peptic_ulcer, by = "medcode") %>%
  # retain patients with diagnosis in the 6 months prior to indexdate
  filter(eventdate >= indexdate - 180 & eventdate < indexdate) %>%
  filter(!duplicated(patid)) %>%
  select(patid) %>%
  mutate(recent_peptic_ulcer = 1)

## Number of consultations in the year prior to indexdate

# consultation type to be counted
cons <- c(
  "Surgery consultation",
  "Follow-up/routine visit",
  "Clinic",
  "Telephone call from a patient",
  "Acute visit",
  "Home Visit",
  "Emergency Consultation"
)

# count number of consultations
numbers_consult <- data$consultations %>%
  inner_join(cohort, by = "patid") %>%
  semi_join(as.data.frame(cons), # filter for consultation types of interest only
            by = c('constype' = 'cons')) %>%
  filter(eventdate >= indexdate - 365 & eventdate <= indexdate) %>%
  group_by(patid) %>%
  summarise(nconsult = n())


# Merge for final analysis dataset ----------------------------------------

analysis_data <- cohort %>%
  left_join(demographics, by = "patid") %>%
  left_join(bmi, by = "patid") %>%
  left_join(enddates, by = "patid") %>%
  left_join(numbers_consult, by = "patid") %>%
  left_join(diabetes, by = "patid") %>%
  left_join(gastric_cancer, by = "patid") %>%
  left_join(recent_gerd, by = "patid") %>%
  left_join(peptic_ulcer, by = "patid") %>%
  mutate(across(
    .cols = c(
      "prior_diabetes",
      "prior_gastric_cancer",
      "recent_gerd",
      "recent_peptic_ulcer"
    ),
    .fns = ~ ifelse(is.na(.x), 0, 1)
  ),
  followup_dur = as.numeric((enddate - indexdate) / 365.25)) %>% 
  droplevels() # drup the 'Unknown' level from eth5

## Quick check
summary(analysis_data)
glimpse(analysis_data)
analysis_data %>% filter(ppi == 1) %>% summary()
analysis_data %>% filter(ppi == 0) %>% summary()

saveRDS(analysis_data, "analysis_data.rds")
