#---------------------------------------------------#
# Prediction models for albuminuria
# Script 1: cohort derivation
# Roemer J. Janse - 2024/01/10
#---------------------------------------------------#

# 0. Set-up ----
# Load packages
pacman::p_load("dplyr",           # Data wrangling
               "tidyr",           # Data tidying
               "magrittr",        # Efficient pipelines
               "stringr"          # Working with strings
)

# Set working directionary
setwd("C:/Users/rjjanse/OneDrive - LUMC/Research/Projects/13. dm_alb/codes/dataframes")

# 1. Creating functions ----
# Wrapper for getting number of individuals
noi <- function() n_distinct(cohort[["lopnr"]])

# Derive diagnostic codes
diag_codes <- function(icd, keep_first = TRUE){
    # Get diagnoses of interest
    dat_tmp <- diagn_roemer
        
    # Keep only diabetes type 2 diagnoses
    dat_tmp <- filter(dat_tmp, grepl(icd, diagnosis)) %>%
        # Change date to date format
        mutate(diag_dt = as.Date(bdat)) %>%
        # Keep only relevant columns
        select(lopnr, diag_dt, diagnosis)

    # Keep first if requested
    if(keep_first){
        # Keep first row
        dat_tmp <- dat_tmp %>%
            # Arrange per person on date
            arrange(lopnr, diag_dt) %>%
            # Group per person
            group_by(lopnr) %>%
            # Keep only first diagnosis
            slice(1L) %>%
            # Remove grouping
            ungroup()
    }
    
    # Return data
    return(dat_tmp)
}

# Reduce size based on diabetes codes
reduce <- function(.data) filter(.data, lopnr %in% dat_diab[["lopnr"]])

# CKD-EPI SCr
ckd_epi_scr <- function(creatinine, creat_dt, birth_dt, female, year = 2009){
    # If year not 2009 or 2021, stop
    if(!(year %in% c(2009, 2021))) stop(paste0("No formula available for CKD-EPI SCr ", year))
    
    # K coefficient
    k <- ifelse(female, 62, 80)
    
    # Alpha coefficient
    alpha <- case_when(year == 2009 & female ~ -0.329,
                       year == 2009 & !female ~ -0.411,
                       year == 2021 & female ~ -0.241,
                       year == 2021 & !female ~ -0.302)
    
    # Calculate age at date of birth
    age <- as.numeric((creat_dt - birth_dt) / 365.25)
    
    # Intercept for 2009 or 2021 formula
    intercept <- ifelse(year == 2009, 141, 142)
    
    # Power coefficient for 2009 or 2021 formula
    power_coef <- ifelse(year == 2009, -1.209, -1.200)
    
    # Age coefficient for 2009 or 2021 formula
    age_coef <- ifelse(year == 2009, 0.993, 0.9938)
    
    # Female coefficient 
    female_coef <- ifelse(year == 2009, 1.018, 1.012)
    
    # Calculate eGFR
    egfr <- intercept * (pmin(creatinine / k, 1) ^ alpha) * (pmax(creatinine / k, 1) ^ (power_coef)) * 
        (age_coef ^ age) * ifelse(female, female_coef, 1)
    
    # Return eGFR
    return(egfr)
}

# 2. Prepare data sources ----
## 2.1. Albumin measurements ----
load("data/albumin_complete.Rda")

# Prepare data (alb tests are already not inpatient)
dat_alb <- albuminuria_lab_test_clean %>%
    # Change datum to date
    mutate(alb_dt = datum) %>%
    # Keep only relevant variables
    select(lopnr, alb_dt, result, test, cat, op, ip)

## 2.2. Diabetes diagnosis ----
# Load diagnostic data
load("data/diagn.Rda")

# Derive diabetes diagnoses and bind data together
dat_diab <- diag_codes(icd = "E11")

## 2.3. Renal complications ----
dat_ren_com <- diag_codes(icd = "E112|N0[0-8]|N1[0-9]|N2[5-9]|Q6[013]|I1[23]")

## 2.4. Creatinine measurements ----
# Creatinine data
load("data/creat.Rdata")

# Keep only individuals who have diabetes to reduce size
creat <- reduce(creat)

# Arterial creatinine
load("data/artcreat.Rdata")

# Keep only individuals who have diabetes to reduce size
artcreat <- reduce(artcreat)

# Venous creatinine
load("data/vencreat.Rdata")

# Keep only individuals who have diabetes to reduce size
vencreat <- reduce(vencreat)

# Demographics data
load("data/demo.Rdata")

# Calculate eGFR from creatinine
dat_egfr <- bind_rows(mutate(creat, type = "creat"), 
                      mutate(artcreat, type = "artcreat"),
                      mutate(vencreat, type = "vencreat")) %>%
    # Keep only relevant rows
    select(lopnr, datum, resultat, ip) %>%
    # Add demographics
    left_join(select(demo, lopnr, female, dob_s3), "lopnr") %>%
    # Keep only non-in patient values
    filter(ip == 0) %>%
    # Change dates to dates
    mutate(# Date of birth
           dob = as.Date(dob_s3),
           # Creat date
           creat_dt = as.Date(datum)) %>%
    # Calculate eGFR
    mutate(egfr = ckd_epi_scr(resultat, creat_dt, dob, female, year = 2009)) %>%
    # Keep only relevant columns
    select(lopnr, creat_dt, egfr)

## 2.5. Demographic data ----
# Process demographic data
dat_demo <- demo %>%
    # Clean variables
    mutate(female = female_s3,                 # Sex
           birth_dt = as.Date(dob_s3)) %>%     # Birth date
    # Keep only relevant variables
    select(lopnr, female, birth_dt)

## 2.6. Migration data ----
# Migration data
load("data/migration.Rdata")

# Join data together
dat_migr <- migration %>%
    # Clean migration date
    mutate(migration_dt = as.Date(hdat)) %>%   
    # Keep only relevant variables
    select(lopnr, migration_dt, hkod) 

# Remove data
rm(artcreat, creat, demo, migration, diagn_roemer, albuminuria_lab_test_clean, vencreat)

## 2.7. Death data ----
# Load data
load("data/death.Rdata")

# Prepare death data
dat_death <- death %>%
    # Keep only relevant variables
    select(lopnr, dodsdat) %>%
    # Change date of death to date format
    mutate(death_dt = as.Date(dodsdat), .keep = "unused") %>%
    # Keep one row per person
    distinct(lopnr, death_dt)

# 3. Cohort derivation ----
## 3.1. Initiate cohort with normal albumin measurements ----
# Only normal albumin measurements
cohort <- dat_alb %>%
    # Arrange per person on alb date
    arrange(lopnr, alb_dt) %>%
    # Group per person
    group_by(lopnr) %>%
    # Make an indicator for albuminuria
    mutate(albuminuria = ifelse(result >= 30, 1, NA)) %>%
    # Fill the indicator downwards (to all next measurements)
    fill(albuminuria, .direction = "down") %>%
    # If albuminuria is still missing, no previous albuminuria was present
    mutate(albuminuria = ifelse(is.na(albuminuria), 0, albuminuria)) %>%
    # Remove grouping structure
    ungroup() %>%
    # Remove all albuminurias and subsequent measurements
    filter(albuminuria == 0) %>%
    # Remove all non-outpatient/primary care measurements
    filter(ip == 0) %>%
    # Remove indicators
    select(-albuminuria, -op, -ip)

# Individuals = 916,724
noi()

## 3.2. Urine albumin measurements date selection ----
cohort <- cohort %>%
    # Keep measurements between 01-01-2007 and 31-12-2021
    filter(between(alb_dt, as.Date("2007-01-01"), as.Date("2021-12-31")))

# Individuals = 887,600
noi()

## 3.3. Only measurements where the individual was 18 ----
# This also adds demographics information on baseline age and sex
cohort <- cohort %>%
    # Join age data
    left_join(dat_demo, "lopnr") %>%
    # Calculate age
    mutate(age = as.numeric(alb_dt - birth_dt) / 365.25) %>%
    # Keep only individuals above 18 years old
    filter(age >= 18)
    
# Individuals = 750,656
noi()

## 3.4. Only measurements after T2DM diagnosis ----
cohort <- cohort %>%
    # Add diabetes diagnosis
    left_join(dat_diab, "lopnr") %>%
    # Keep only rows where the diagnosis is on or before the day of the albumin measurement (this also drops NAs)
    filter(diag_dt <= alb_dt) %>%
    # Rename diag_dt
    rename(dm_dt = diag_dt)

# Individuals = 91,011
noi()

## 3.5. No renal complications prior to diabetes ----
cohort <- cohort %>%
    # Add renal complications
    left_join(select(dat_ren_com, -diagnosis), "lopnr") %>%
    # Keep only rows where there is no diagnosis or diagnosis is after albumin measurement
    filter(diag_dt > alb_dt | is.na(diag_dt)) %>%
    # Remove diagnosis date and source
    select(-diag_dt)

# Individuals = 85,332
noi()

## 3.6. Kidney function >=60 mL/min/1.73m2 ----
cohort <- cohort %>%
    # Join eGFR data
    left_join(dat_egfr, "lopnr", relationship = "many-to-many") %>%
    # Keep onlly measurements with eGFR >=60
    filter(egfr >= 60)
    
# Individuals = 83,994
noi()

## 3.7. Migrated before baseline ----
cohort <- cohort %>%
    # Join migration
    left_join(dat_migr, "lopnr", relationship = "many-to-many") %>%
    # Set any migration after albumin date to NA (also date for sorting)
    mutate(across(migration_dt:hkod, \(x) x = ifelse(migration_dt > alb_dt, NA, x))) %>%
    # Change migration_dt back to date format
    mutate(migration_dt = as.Date(migration_dt, origin = "1970-01-01")) %>%
    # Arrange per person on albumin date and migration date
    arrange(lopnr, alb_dt, migration_dt) %>%
    # Group per person and albumin date
    group_by(lopnr, alb_dt) %>%
    # If the last migration type per person is immigration (I) or missing (no migration), there is no problem
    mutate(indicator = ifelse(last(hkod, na_rm = TRUE) == "I" | is.na(last(hkod, na_rm = TRUE)), 1, 0)) %>%
    # Keep only measurements without emigration before measurement
    filter(indicator == 1) %>%
    # Remove unncessary variables
    select(-c(indicator, hkod, migration_dt)) %>%
    # Keep one row per person per albumin date
    slice(1L) %>%
    # Remove grouping structure
    ungroup()

# Individuals = 83,659
noi()

## 3.8. Died before or at baseline ----
# This is possible if lab is communicated later
cohort <- cohort %>%
    # Join death data
    left_join(dat_death, "lopnr") %>%
    # Keep only if death date missing (not died yet) or after baseline
    filter(is.na(death_dt) | death_dt > alb_dt)

# Individuals = 83,658
noi()

## 3.9. Keep first albumin measurement per person ----
# This will be their baseline, no individuals are excluded here
cohort <- cohort %>%
    # Arrange per person per albumin test
    arrange(lopnr, alb_dt) %>%
    # Group per person
    group_by(lopnr) %>%
    # Keep first row per person
    slice(1L) %>%
    # Remove grouping structure
    ungroup() %>%
    # Rename variables
    rename(index_dt = alb_dt, index_alb = result, dm_type = diagnosis) %>%
    # Create indicator for development and validation cohort
    mutate(cohort = ifelse(index_dt <= as.Date("2013-12-31"), "development", "validation")) %>%
    # Change variable order
    select(lopnr, cohort, index_dt, index_alb, egfr, female, birth_dt, age, dm_dt, dm_type)

# Save cohort
save(cohort, file = "cohort.Rdata")
