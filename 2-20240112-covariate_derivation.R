#---------------------------------------------------#
# Prediction models for albuminuria
# Script 2: covriate & predictor derivation
# Roemer J. Janse - 2024/01/12
#---------------------------------------------------#

# 0. Set-up ----
# Load packages
pacman::p_load("dplyr",           # Data wrangling
               "tidyr",           # Data tidying
               "magrittr",        # Efficient pipelines
               "stringr",         # Working with strings
               "mice"             # Multiple imputation
)

# Set working directionary
setwd("C:/Users/rjjanse.LUMCNET/OneDrive - LUMC/Research/Projects/13. dm_alb/codes/dataframes")

# Load cohort
load("cohort.Rdata")

# 1. Creating functions ----
# Derive diagnoses
diag <- function(.data, icd, name = NULL){
    # Derive diagnoses of interest
    dat_tmp <- cohort %>%
        # Join data
        left_join(.data %>%
                      # Keep only diagnosis of interest
                      filter(grepl(icd, diagnosis)), "lopnr") %>%
        # Create indicator for ICD code present
        mutate(diag = if_else(diag_dt <= index_dt, 1, 0)) %>%
        # Arrange for grouping
        arrange(lopnr, desc(diag)) %>%
        # Group per person
        group_by(lopnr) %>%
        # Keep only first row per person, which is 1 if diagnosis was present before index
        slice(1L) %>%
        # Ungroup again
        ungroup() %>%
        # Any remaining NA means no diagnosis
        replace_na(list(diag = 0)) %>%
        # Keep one row per person and select only variables of interest
        distinct(lopnr, diag)
    
    # Do not rename variable if no name is given
    if(is.null(name)) return(dat_tmp)
    
    # Rename and return if name is given
    else {
        # Rename variable
        dat_tmp[[name]] <- dat_tmp[["diag"]]
        
        # Remove variable
        dat_tmp <- select(dat_tmp, -diag)
        
        # Return data
        return(dat_tmp)
    }
}

# Derive drugs
drug <- function(.data, drg, name = NULL){
    # Derive diagnoses of interest
    dat_tmp <- cohort %>%
        # Join data
        left_join(.data %>%
                      # Keep only diagnosis of interest
                      filter(grepl(drg, atc)), "lopnr") %>%
        # Create indicator for ICD code present in the year before
        mutate(drug = if_else(drug_dt <= index_dt & drug_dt >= (index_dt - 365.25), 1, NA)) %>%
        # Arrange for grouping
        arrange(lopnr, desc(drug)) %>%
        # Group per person
        group_by(lopnr) %>%
        # Keep only first row per person, which is 1 if diagnosis was present before index
        slice(1L) %>%
        # Ungroup again
        ungroup() %>%
        # Any remaining NA means no diagnosis
        replace_na(list(drug = 0)) %>%
        # Keep one row per person and select only variables of interest
        distinct(lopnr, drug)
    
    # Do not rename variable if no name is given
    if(is.null(name)) return(dat_tmp)
    
    # Rename and return if name is given
    else {
        # Rename variable
        dat_tmp[[name]] <- dat_tmp[["drug"]]
        
        # Remove variable
        dat_tmp <- select(dat_tmp, -drug)
        
        # Return data
        return(dat_tmp)
    }
}

# Derive lab
lab <- function(.data, lab, name = NULL){
    # Derive diagnoses of interest
    dat_tmp <- cohort %>%
        # Join data
        left_join(.data %>%
                      # Keep only diagnosis of interest
                      filter(grepl(lab, test)), "lopnr") %>%
        # Create indicator for ICD code present in the year before
        mutate(resultat = if_else(lab_dt <= index_dt & lab_dt >= (index_dt - 365.25), resultat, NA)) %>%
        # Arrange for grouping
        arrange(lopnr) %>%
        # Calculate mean per person
        mutate(resultat = mean(resultat, na.rm = TRUE), .by = "lopnr") %>%
        # Keep one row per person and select only variables of interest
        distinct(lopnr, resultat)
    
    # Do not rename variable if no name is given
    if(is.null(name)) return(dat_tmp)
    
    # Rename and return if name is given
    else {
        # Rename variable
        dat_tmp[[name]] <- dat_tmp[["resultat"]]
        
        # Remove variable
        dat_tmp <- select(dat_tmp, -resultat)
        
        # Return data
        return(dat_tmp)
    }
}
    
# 2. Deriving diagnoses ----
## 2.0. Preparing data ----
# Load diagnostics data
load("data/diagn.Rda")

# Prepare diagnostics data
dat_diags <- diagn_roemer %>%
    # Keep only individuals that are in the cohort
    filter(lopnr %in% cohort[["lopnr"]]) %>%
    # Change bdat
    mutate(diag_dt = as.Date(bdat), .keep = "unused")

# Remove data
rm(diagn_roemer)

## 2.1. Diabetic retinopathy
diag_ret <- diag(dat_diags, "E113[ABC]", "retinopathy")

## 2.2. Diabetic neuropathy ----
diag_neu <- diag(dat_diags, "E114[BCD]", "neuropathy")

## 2.3. Diabetic foot ulcer ----
diag_ulc <- diag(dat_diags, "E116D", "ulcer")

## 2.4. Hypertension ----
diag_hyp <- diag(dat_diags, "I1[0-5]", "hypertension")

## 2.5. Ischemic heart disease ----
diag_ihd <- diag(dat_diags, "I2[0-5]", "ihd")

## 2.6. Cerebrovascular disease ----
diag_cvd <- diag(dat_diags, "I6[0-9]|G45[01]", "cvd")

## 2.7. Peripheral vascular disease ----
diag_pvd <- diag(dat_diags, "I7[0-2]", "pvd")

## 2.8. Congestive heart failure ----
diag_chf <- diag(dat_diags, "I50|I110|I13[02]|I255|K761|I4[23]", "chf")

## 2.9. Atrial fibrillation ----
diag_afb <- diag(dat_diags, "I48", "fibrillation")

# 3. Deriving drug use ----
## 3.0. Prepare medication data ----
# Load medication data
load("data/meds.Rda")

# Prepare medication data
dat_meds <- meds_roemer %>%
    # Keep only individuals that are in the cohort
    filter(lopnr %in% cohort[["lopnr"]]) %>%
    # Keep only relevant columns
    select(lopnr, atc, edatum) %>%
    # Rename edatum
    rename(drug_dt = edatum)

# Remove data
rm(meds_roemer)

## 3.1. Insulin ----
drug_ins <- drug(dat_meds, "A10A", "insulin")

## 3.2. Blood glucose lowering drugs (excl. insulin) ----
drug_bgl <- drug(dat_meds, "A10B|A10XA", "glucose_lowering")

## 3.3. ACEis and ARBs ----
drug_ras <- drug(dat_meds, "C09[A-D]", "rasi")

## 3.4. Statins ----
drug_stn <- drug(dat_meds, "C10AA", "statins")

## 3.5. Aspirin ----
drug_asp <- drug(dat_meds, "B01AC06|N02BA01", "aspirin")

## 3.6. Antihypertensives ----
drug_hyp <- drug(dat_meds, "C02", "antihypertensives")

## 3.7. Anticoagulants ----
drug_coa <- drug(dat_meds, "B01", "anticoagulants")

## 3.8. Beta-blockers ----
drug_bbs <- drug(dat_meds, "C07", "bblockers")

## 3.9. Calcium channel blockers ----
drug_ccb <- drug(dat_meds, "C08[CD]", "ccbs")

## 3.10. MRAs ----
drug_mra <- drug(dat_meds, "C03DA", "mras")

# 4. Deriving laboratory ----
## 4.0. Preparing laboratory data ----
# Load data
load("data/labs.Rda")

# Load spreadsheet with tests
tests <- readxl::read_xlsx("c:/users/rjjanse.lumcnet/onedrive - lumc/research/projects/13. dm_alb/spreadsheets/tests.xlsx") %>%
    # Rename variables
    rename(analys = test, test = lab)

# Prepare lab data
dat_labs <- lab_roemer %>%
    # Keep only individuals in cohort
    filter(lopnr %in% cohort[["lopnr"]]) %>%
    # Add lab
    left_join(tests, "analys") %>%
    # Some tests are not in the .xlsx file, these are mainly HbA1c, rest is LDL cholesterol ratio. Add test here
    mutate(test = case_when(is.na(test) & grepl("HbA1c", analys) ~ "hba1c",
                            is.na(test) & grepl("LDL-kol", analys) ~ "ldl",
                            .default = test)) %>%
    # Clean resultat
    mutate(resultat = str_replace(resultat, ",", "."),               # Comma's to dots for decimals
           resultat = str_replace_all(resultat, "\\<|\\>|\\*", ""),  # Additional notation
           resultat = as.numeric(resultat)) %>%                      # Change character to numeric
    # Keep only not inpatient measures
    filter(ip == 0) %>%
    # Change lab date to date format
    mutate(lab_dt = as.Date(datum), .keep = "unused") %>%
    # Convert based on units
    mutate(resultat = case_when(# HbA1c ( to mmol/mol)
                                test == "hba1c" & enhet %in% c("%", "% av totalt Hb") ~ 10.93 * resultat - 23.5,
                                test == "hba1c" & is.na(enhet) ~ 10.93 * resultat - 23.5,
                                # HDL, LDL, and total cholesterol (to mg/dL)
                                test %in% c("hdl", "ldl", "tc") & enhet %in% c("mmol/L", "mmol/l") ~ resultat * 38.67,
                                # LDL (to mg/dL)
                                test == "ldl" & is.na(enhet) ~ resultat * 38.67,
                                # Triglycerides (to mg/dL)
                                test == "triglyc" & enhet %in% c("mmol/L", "mmol/l") ~ resultat * 88.57,
                                # Urine creatinine (to mg/dL)
                                test == "ucreat" & enhet %in% c("mmol/L", "mmol/d") ~ resultat * 11.312,
                                # Remainder stays the same
                                .default = resultat)) %>%
    # Remove negative values
    filter(resultat >= 0) %>%
    # Keep only relevant columns
    select(lopnr, resultat, lab_dt, test) 

# Create data for total cholesterol:HDL cholesterol ratio
dat_tc_hdl <- dat_labs %>%
    # Keep only relevant tests
    filter(test %in% c("hdl", "tc")) %>%
    # Arrange for grouping
    arrange(lopnr, lab_dt, test) %>%
    # Group per person per date per test
    group_by(lopnr, lab_dt, test) %>%
    # Get mean per day for days with multiple tests
    mutate(resultat = mean(resultat)) %>%
    # Keep one row per test per day per person
    slice(1L) %>%
    # Remove grouping structure
    ungroup() %>%
    # Rearrange to make sure no mistakes are made in arranging
    arrange(lopnr, lab_dt, test) %>%
    # Regroup to person and date
    group_by(lopnr, lab_dt) %>%
    # Check length of different labs (should be 2 (hdl and tc) if we can calculate ratio
    mutate(calculateable = ifelse(length(unique(test)) == 2, 1, 0)) %>%
    # Keep only dates of individuals where we  can calculate tc:hdl ratio
    filter(calculateable == 1) %>%
    # Remove grouping structure
    ungroup() %>%
    # Pivot to wide format for calculation
    pivot_wider(names_from = test, values_from = resultat) %>%
    # Calculate ratio
    mutate(resultat = tc / hdl, test = "tc_hdl") %>%
    # Keep only relevant variables
    select(lopnr, lab_dt, resultat, test)

# Add data together
dat_labs <- rbind(dat_labs, dat_tc_hdl) %>%
    # Arrange per person per date per test
    arrange(lopnr, lab_dt, test)
        
# Remove data
rm(lab_roemer, dat_tc_hdl)

## 4.1. Urine creatinine ----
lab_ucrea <- lab(dat_labs, "ucrea", "ucrea")

## 4.2. HbA1c ----
lab_hba1c <- lab(dat_labs, "hba1c", "hba1c")

## 4.3. Total cholesterol ----
lab_tc <- lab(dat_labs, "tc", "total_cholesterol") 

## 4.4. HDL cholesterol ----
lab_hdl <- lab(dat_labs, "hdl", "hdl")

## 4.5. LDL cholesterol ----
lab_ldl <- lab(dat_labs, "ldl", "ldl")

## 4.6. Triglycerides ----
lab_triglyc <- lab(dat_labs, "triglyc", "triglycerides")

## 4.7. Total cholesterol:HDL ratio ----
lab_tchdl <- lab(dat_labs, "tc_hdl", "tc_hdl")

# 5. Deriving educational level ----
# Load data
load("data/education.Rdata")

# Sun2000niva.old has 7 levels, see legend of Figure 2 by Ludvigsson
# 1 = Compulsory school, less than 9 years (did not complete compulsory education [<9 years])
# 2 = Compulsory school, 9 years (completed compulsory education [9 years])
# 3 = Secondary school (upper secondary [2 years])
# 4 = Secondary school (upper secondary [3 years])
# 5 = University (college/university < 3 years)
# 6 = University (college/university >= 3 years)
# 7 = University (research education)
# Coded:
# 0 = Low category (up to 4)
# 1 = High category (from 5 onwards)
dat_educ <- education %>%
    # Create categories for education
    mutate(education = if_else(sun2000niva_old <= 4, 0, 1, missing = NA), .keep = "unused")

# Determine education
# One error on coercion leading to NAs is fine: these values where already NA
cohort <- cohort %>%
    # Join educational data
    left_join(dat_educ, "lopnr") %>%
    # Mutate variables
    mutate(# Get year of baseline
           index_year = as.numeric(format(index_dt, "%Y")),
           # If exam year is **** assume that the educational level was achieved before index
           exam_year = if_else(exam_year == "****", "0", exam_year),
           # Change exam year to numeric
           exam_year = as.numeric(exam_year),
           # Set education to NA if exam year is after baseline year
           education = if_else(exam_year > index_year, NA, education)) %>%
    # Sort for grouping and slicing (highest level up)
    arrange(lopnr, desc(education)) %>%
    # Group per person
    group_by(lopnr) %>%
    # Keep one row per person
    slice(1L) %>%
    # Remove grouping structure
    ungroup() %>%
    # Remove unnecessary variables
    select(-c(exam_year, index_year))

# 6. Joining data together ----
placeholder <- 
    # Join data frames
    plyr::join_all(lapply(ls()[str_detect(ls(), "cohort|diag_|drug_|^lab_")], get), by = "lopnr", type = "left") %>%
    # Some lab values are infinite, set to NA
    mutate(across(c(tc_hdl, hdl, total_cholesterol), \(x) x = if_else(is.infinite(x), NA, x)))

# Data to cohort
cohort <- placeholder %>%
    # Calculate time to diabetes in months
    mutate(diab_months = as.numeric(index_dt - dm_dt) / (365.25 / 12), .after = "dm_type")

# Save data
save(cohort, file = "cohort_covariates.Rdata")

# 7. Derive outcome ----
# Load laboratory data
load("data/albumin_complete.Rda")

# Load death data
load("data/death.Rdata")

# Prepare death data
dat_death <- death %>%
    # Keep only relevant variables
    select(lopnr, dodsdat) %>%
    # Change date of death to date format
    mutate(death_dt = as.Date(dodsdat), .keep = "unused") %>%
    # Keep one row per person
    distinct(lopnr, death_dt)

# Add outcomes to cohort
cohort <- cohort %>%
    # Join albumin tests
    left_join(albuminuria_lab_test_clean %>%
                  # Keep only tests of A2 and A3 which are outcomes of interest (not A1)
                  filter(cat %in% c("A2", "A3")), 
              "lopnr") %>%
    # Join death data
    left_join(dat_death, "lopnr") %>%
    # Set any test before index to NA
    mutate(lab_dt = if_else(datum <= index_dt, NA, datum)) %>%
    # Create outcome variables
    mutate(## Fine-Gray & AFT prediction models
           # Date at three years
           year3_dt = index_dt + 365.25 * 3,
           # Administrative censoring
           cens_dt = as.Date("2021-12-31"),
           # For microalbuminuria, choose earliest date as outcome
           microalb_y3_dt = pmin(year3_dt, death_dt, cens_dt, lab_dt, na.rm = TRUE),
           # If date is equal to test, event is microalbuminuria, if equal to death then death, else censoring
           microalb_y3 = case_when(microalb_y3_dt == lab_dt ~ 1,
                                   microalb_y3_dt == death_dt ~ 2,
                                   microalb_y3_dt == cens_dt ~ 3,
                                   microalb_y3_dt == year3_dt ~ 3),
           # Calculate time to event for 3 years
           tte_y3 = as.numeric(microalb_y3_dt - index_dt),
           ## Multi-state model
           # Date of microalbuminurias (if no measurement then impossibly high date)
           microalb_dt = as.Date(ifelse(cat == "A2" & !is.na(lab_dt), lab_dt, 50000), origin = "1970-01-01"),
           # Date of macroalbuminurias (if no measurement then impossibly high date)
           macroalb_dt = as.Date(ifelse(cat == "A3" & !is.na(lab_dt), lab_dt, 50000), origin = "1970-01-01")) %>%
    # Put outcomes of interest on top, then death, then censoring
    arrange(lopnr, microalb_y3_dt, microalb_y3) %>%
    # Group per person
    group_by(lopnr) %>%
    # Update variables
    mutate(# First microalbuminuria date for each person
           microalb_dt = min(microalb_dt, na.rm = TRUE),
           # First macroalbuminuria date for each person
           macroalb_dt = min(macroalb_dt, na.rm = TRUE),
           # Censoring date is smallest of 3 years and administrative censoring
           censor_dt = pmin(cens_dt, year3_dt),
           # Determine if an individual had microalbuminuria (should have occurred before censoring, death, and macroalbuminuria)
           stage_microalb = if_else(microalb_dt == pmin(microalb_dt, macroalb_dt, censor_dt, death_dt, na.rm = TRUE), 1, 0, missing = 0),
           # Determine if an individual had macroalbuminuria (should have occurred before censoring and death)
           stage_macroalb = if_else(macroalb_dt == pmin(macroalb_dt, censor_dt, death_dt, na.rm = TRUE), 1, 0, missing = 0),
           # Determine if an individual died (should have occurred before censoring)
           stage_death = if_else(death_dt == pmin(death_dt, censor_dt, na.rm = TRUE), 1, 0, missing = 0)) %>%
    # Keep one row (thus one outcome) per person
    slice(1L) %>%
    # Ungroup again
    ungroup() %>%
    # Change outcomes to factor
    mutate(microalb_y3 = factor(microalb_y3, levels = c(3, 1, 2), labels = c("censored", "microalbuminuria", "died"))) %>%
    # Keep only relevant columns
    select(-(datum:enhet), -(lab_dt:cens_dt)) %>%
    # Also change NaNs to NAs
    replace_na(list(hba1c = NA, hdl = NA, ldl = NA, total_cholesterol = NA, tc_hdl = NA, triglycerides = NA, ucrea = NA)) %>%
    # Rename lopnr to studynr
    rename(studynr = lopnr) %>%
    # Set dates for no death artificially high to bypass imputation
    mutate(death_dt = as.Date(ifelse(stage_death == 0, 50000, death_dt), origin = "1970-01-01"))

# 8. Impute data ----
# Predictor matrix check
pred_mat <- mice(cohort, m = 1, maxit = 0, seed = 1)[["predictorMatrix"]]

# Edit predictor matrix
pred_mat <- pred_mat %>%
    # Change to data frame
    as.data.frame() %>%
    # Set columns to 0
    mutate(across(c(microalb_y3_dt, microalb_dt, macroalb_dt, censor_dt, death_dt), \(x) x = 0))

# Set rows of predictor matrix to 0 too for dates
for(i in c("microalb_y3_dt", "microalb_dt", "macroalb_dt", "censor_dt", "death_dt")) pred_mat[i, ] <- 0

# Predictor matrix back to matrix
pred_mat <- as.matrix(pred_mat)

# Imputation object
cohort_imputed <- mice(cohort, m = 10, maxit = 50, seed = 1, method = "pmm", predictorMatrix = pred_mat)

# Save imputation object
save(cohort_imputed, file = "imputation.Rdata")

# Get imputed data
cohort <- complete(cohort_imputed, action = "long") %>%
    # Set dates for absent stages to missing
    mutate(# Microalbuminuria
           microalb_dt = as.Date(ifelse(stage_microalb == 0, NA, microalb_dt), origin = "1970-01-01"),
           # Macroalbuminuria
           macroalb_dt = as.Date(ifelse(stage_macroalb == 0, NA, macroalb_dt), origin = "1970-01-01"),
           # Death
           death_dt = as.Date(ifelse(stage_death == 0, NA, death_dt), origin = "1970-01-01"))

# Save imputed data
save(cohort, file = "cohort_imputed.Rdata")

# Remove everything
rm(list = ls())


