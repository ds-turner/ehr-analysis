# Example code for a mock retrospective cohort study using EHRs.

## Introduction
The following code was written to clean and analyse artificial data that follows a format similar to CPRD GOLD, as part of a mock retrospective cohort study.  The original data is not included in this repository but the analysis dataset is.  

The mock study aimed to compare mortality among new users of proton pump inhibitors (PPIs) to new users of an alternative acid suppression drug, H2 receptor antagonists (H2RAs).  Hazard ratios were estimated using a multivariate Cox proportional hazard model.

CPRD collects fully-coded patient level EHRs from GPs in the UK. Researches can access anonymised patient records that match the criteria of CPRD approved protocols.  Anonymised patient data is provided in 2 formats, CPRD GOLD and CPRD Aurum, depending on the system used to collect the data at the GP.  More information, including data specifications and data resource profiles can be found [here](https://cprd.com/primary-care-data-public-health-research)

The following code was written to clean and analyse artificial data that follows a format similar to CPRD GOLD.  The original data is not included in this repository but the analysis dataset is.

## Purpose of scripts

'ehr-clean.R' manipulates the raw data to produce the analysis dataset.  In brief it applies the inclusion/exclusion criteria of the study to inclued only patients that meet the criteria then it derives veriables to be included in the analyis.  To derive these variables codelists are sometimes used to identify information such as a diagnosis of a disease or prscreption of a medician.  It then saves a copy of the analysis dataset, that can be found in this respository.

'ehr-analyse.R' loads the analysis dataset then fits a cox regression modle, checking the proportional hazard assumption and the correct for with which to fit the covariates

## Mock study design

### Eligibility criteria:
Patients whose records are held in the database, who have:
*	a recorded year of birth, sex and ethnicity
*	at least one prescription for either a PPI or a H2RA ++
*	acceptable quality of data for research ++

### Further, the first such prescription must be:
*	after the patient’s 18th birthday ++
*	after the date of the patient's registration at the general practice plus one year
*	between 1 January 1991 and 17 April 2017

++ Eigibility criteria imposed prior to raw data recieved

### Exposure

New prescription for a PPI compared to new prescriptions for a H2RA.

### Potential confounders:
* Age
* sex, 
* ethnicity, 
* deprivation (measured by the Index of Multiple Deprivation, IMD), 
* body mass index (BMI), 
* diagnosis of diabetes prior to prescription, 
* diagnosis of gastric cancer prior to prescription, 
* gastro-oesophageal reflux disease (GERD) in the 6 months prior to prescription, 
* peptic ulcer in the 6 months prior to prescription, 
* number of consultations in the year prior to prescription. 

### Follow-up 

For each patient, follow-up begins at the first prescription for PPI or H2RA and ends at the first of: 
* death, 
* transfer out of the general practice or 
* end of study (17 April 2017).

### Outcome 

All cause-mortality




