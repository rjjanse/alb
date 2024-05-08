#---------------------------------------------------#
# Prediction models for albuminuria
# Shiny app server
# Roemer J. Janse - 2024/05/01
#---------------------------------------------------#

# 0. Set-up ----
# Load packages
library(dplyr)          # Data wrangling
library(magrittr)       # Better pipelines
library(stringr)        # Working with strings
library(survival)       # Survival models
library(ggplot2)        # Data visualization
library(broom)          # Tidy data
library(splines)        # Splines
library(cowplot)        # Data viz. add-on
library(plotly)         # Interactive plots
library(mstate)         # Multistate modelling
library(tidyr)          # Tidying data
library(shiny)          # Shiny web app
library(bslib)          # Bootstrap layouts

# Load functions
source("https://raw.githubusercontent.com/rjjanse/alb/main/3x-functions-20240215.R")

# Define server logic required to draw a histogram
function(input, output, session){
    # Output prediction plot
    output[["prediction_plot"]] <- renderPlotly({
        # Check value of input model
        if(input[["model"]] == "Fine-Gray"){
            survcurv(female = as.numeric(input$female),
                     age = input$age,
                     education = as.numeric(input$education),
                     diabetes_months = input$diab_month,
                     egfr = input$egfr,
                     urine_creatinine = input$ucrea, 
                     hba1c = input$hba1c,
                     total_cholesterol = input$tc,
                     hdl = input$hdl,
                     ldl = input$ldl,
                     triglycerides = input$tgcs,
                     total_cholesterol_hdl_ratio = input$tc_hdl_ratio,
                     albumin = input$uacr,
                     atrial_fibrillation = as.numeric(input$afib),
                     congestive_heart_failure = as.numeric(input$chf),
                     cerebrovascular_disease = as.numeric(input$cvd),
                     hypertension = as.numeric(input$hyp),
                     ischemic_heart_disease = as.numeric(input$ihd),
                     neuropathy = as.numeric(input$neuro),
                     peripheral_vascular_disease = as.numeric(input$pvd),
                     retinopathy = as.numeric(input$retino),
                     ulcer = as.numeric(input$ulcer),
                     aspirin = as.numeric(input$aspirin),
                     beta_blockers = as.numeric(input$bbs),
                     glucose_lowering_drugs = as.numeric(input$gld),
                     calcium_channel_blockers = as.numeric(input$ccbs),
                     anticoagulants = as.numeric(input$coags),
                     antihypertensives = as.numeric(input$hyps),
                     insulin = as.numeric(input$insu),
                     mineralocorticoid_receptor_antagonists = as.numeric(input$mra),
                     ras_inhibitors = as.numeric(input$rasi),
                     statins = as.numeric(input$statins),
                     # Interactivity
                     interactive = TRUE
            )
        } else {
            # If area chart, plot area
            if(input[["mstate_output"]] == "Area chart"){
                # Multi-state model
                plot_mstate_prep(female = as.numeric(input$female),
                                 age = input$age,
                                 education = as.numeric(input$education),
                                 diabetes_months = input$diab_month,
                                 egfr = input$egfr,
                                 urine_creatinine = input$ucrea,
                                 hba1c = input$hba1c,
                                 total_cholesterol = input$tc,
                                 hdl = input$hdl,
                                 ldl = input$ldl,
                                 triglycerides = input$tgcs,
                                 total_cholesterol_hdl_ratio = input$tc_hdl_ratio,
                                 albumin = input$uacr,
                                 atrial_fibrillation = as.numeric(input$afib),
                                 congestive_heart_failure = as.numeric(input$chf),
                                 cerebrovascular_disease = as.numeric(input$cvd),
                                 hypertension = as.numeric(input$hyp),
                                 ischemic_heart_disease = as.numeric(input$ihd),
                                 neuropathy = as.numeric(input$neuro),
                                 peripheral_vascular_disease = as.numeric(input$pvd),
                                 retinopathy = as.numeric(input$retino),
                                 ulcer = as.numeric(input$ulcer),
                                 aspirin = as.numeric(input$aspirin),
                                 beta_blockers = as.numeric(input$bbs),
                                 glucose_lowering_drugs = as.numeric(input$gld),
                                 calcium_channel_blockers = as.numeric(input$ccbs),
                                 anticoagulants = as.numeric(input$coags),
                                 antihypertensives = as.numeric(input$hyps),
                                 insulin = as.numeric(input$insu),
                                 mineralocorticoid_receptor_antagonists = as.numeric(input$mra),
                                 ras_inhibitors = as.numeric(input$rasi),
                                 statins = as.numeric(input$statins),
                                 # Interactivity
                                 interactive = TRUE
                )
            } else {
                # Multi-state model
                plot_mstate_prep(female = as.numeric(input$female),
                                 age = input$age,
                                 education = as.numeric(input$education),
                                 diabetes_months = input$diab_month,
                                 egfr = input$egfr,
                                 urine_creatinine = input$ucrea,
                                 hba1c = input$hba1c,
                                 total_cholesterol = input$tc,
                                 hdl = input$hdl,
                                 ldl = input$ldl,
                                 triglycerides = input$tgcs,
                                 total_cholesterol_hdl_ratio = input$tc_hdl_ratio,
                                 albumin = input$uacr,
                                 atrial_fibrillation = as.numeric(input$afib),
                                 congestive_heart_failure = as.numeric(input$chf),
                                 cerebrovascular_disease = as.numeric(input$cvd),
                                 hypertension = as.numeric(input$hyp),
                                 ischemic_heart_disease = as.numeric(input$ihd),
                                 neuropathy = as.numeric(input$neuro),
                                 peripheral_vascular_disease = as.numeric(input$pvd),
                                 retinopathy = as.numeric(input$retino),
                                 ulcer = as.numeric(input$ulcer),
                                 aspirin = as.numeric(input$aspirin),
                                 beta_blockers = as.numeric(input$bbs),
                                 glucose_lowering_drugs = as.numeric(input$gld),
                                 calcium_channel_blockers = as.numeric(input$ccbs),
                                 anticoagulants = as.numeric(input$coags),
                                 antihypertensives = as.numeric(input$hyps),
                                 insulin = as.numeric(input$insu),
                                 mineralocorticoid_receptor_antagonists = as.numeric(input$mra),
                                 ras_inhibitors = as.numeric(input$rasi),
                                 statins = as.numeric(input$statins),
                                 # Area
                                 area = FALSE,
                                 # Interactivity
                                 interactive = TRUE
                )
            }
        }
    })
}
