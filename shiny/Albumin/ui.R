#---------------------------------------------------#
# Prediction models for albuminuria
# Shiny app user interface
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
library(shinythemes)    # Theming
library(bslib)          # Bootstrap layouts

# 1. Define UI ----
page_sidebar(
    # Theme
    theme = shinytheme("flatly"),
    
    # Title
    title = "Predicting the risk of albuminuria in patients with diabetes:  - Janse et al. 2024",
    
    # Sidebar with options
    sidebar = sidebar(
        # Sidebar styling
        position = "right",
        width = 400,
        
        # Model
        selectInput(
            "model",
            label = "Prediction model",
            choices = list("Fine-Gray", "Multi-state"),
            selected = "Fine-Gray"
        ),
        
        # Output for multi-state
        conditionalPanel(
            condition = "input.model == 'Multi-state'",
            selectInput("mstate_output",
                        label = "Multi-state model output",
                        choices = list("Area chart", "Line chart"),
                        selected = "Area chart"
            )
        ),
        
        ## Predictors
        # Subtitle
        h4("Demographics"),
        
        # Sex
        checkboxInput(
            "female",
            label = "Female sex"
        ),
        
        # Age
        numericInput(
            "age",
            label = "Age (years)",
            value = 63,
        ),
        
        # Education
        checkboxInput(
            "education",
            label = "Post-secondary education"
        ),
        

        # Diabetes months
        numericInput(
            "diab_month",
            label = "Time since diabetes diagnosis (months)",
            value = 36
        ),
        
        # Subtitle
        h4("Laboratory"),
        
        # Kidney function
        numericInput(
            "egfr",
            label = "Kidney function (mL/min/1.73m2)",
            value = 90
        ),
        
        # Urine creatinine
        numericInput(
            "ucrea",
            label = "Urine creatinine (mg/dL)",
            value = 109
        ),
        
        # HbA1c
        numericInput(
            "hba1c",
            label = "HbA1c (mmol/mol)",
            value = 49
        ),
        
        # Total cholesterol
        numericInput(
            "tc",
            label = "Total cholesterol (mg/dL)",
            value = 121
        ),
        
        # HDL
        numericInput(
            "hdl",
            label = "HDL (mg/dL)",
            value = 25
        ),
        
        # LDL
        numericInput(
            "ldl",
            label = "LDL (mg/dL)",
            value = 92
        ),
        
        # Triglycerides
        numericInput(
            "tgcs",
            label = "Triglycerides (mg/dL)",
            value = 165
        ),
        
        # Total cholesterol : HDL ratio
        numericInput(
            "tc_hdl_ratio",
            label = "Total cholesterol:HDL ratio",
            value = 4.2
        ),
        
        # uACR
        numericInput(
            "uacr",
            label = "urine albumin-creatinine ratio (mg/g)",
            value = 15
        ),
        
        # Subtitle
        h4("Medical history"),
        
        # Atrial fibrillation
        checkboxInput(
            "afib",
            label = "Atrial fibrillation"
        ),
        
        # Congestive heart failure
        checkboxInput(
            "chf",
            label = "Congestive heart failure"
        ),
        
        # Cerebrovascular disease
        checkboxInput(
            "cvd",
            label = "Cerebrovascular disease"
        ),
        
        # Hypertension
        checkboxInput(
            "hyp",
            label = "Hypertension"
        ),
        
        # Ischemic heart disease
        checkboxInput(
            "ihd",
            label = "Ischemic heart disease"
        ),
        
        # Neuropathy
        checkboxInput(
            "neuro",
            label = "Diabetic neuropathy"
        ),
        
        # Retinopathy
        checkboxInput(
            "retino",
            label = "Diabetic retinopathy"
        ),
        
        # Ulcer
        checkboxInput(
            "ulcer",
            label = "Diabetic foot ulcer"
        ),
        
        # Peripheral vascular disease
        checkboxInput(
            "pvd",
            label = "Peripheral vascular disease"
        ),
        
        # Subtitle
        h4("Medication use"),
        
        # Aspirin
        checkboxInput(
            "aspirin",
            label = "Aspirin"
        ),
        
        # Beta blockers
        checkboxInput(
            "bbs",
            label = "Beta-blockers"
        ),
        
        # RASi
        checkboxInput(
            "rasi",
            label = "ACEis/ARBs"
        ),
        
        # Calcium channel blockers
        checkboxInput(
            "ccbs",
            label = "Calcium channel blockers"
        ),
        
        # Antihypertensives
        checkboxInput(
            "hyps",
            label = "Antihypertensives"
        ),
        
        # MRAs
        checkboxInput(
            "mra",
            label = "Mineralocorticoid recepetor antagonists"
        ),
        
        # Anticoagulants
        checkboxInput(
            "coags",
            label = "Anticoagulants"
        ),
        
        # Statins
        checkboxInput(
            "statins",
            label = "Statins"
        ),
        
        # Glucose lowering drugs
        checkboxInput(
            "gld",
            label = "Glucose-lowering drugs (excl. insulin)"
        ),
        
        # Insulin
        checkboxInput(
            "insu",
            label = "Insulin"
        ),

    ),
    
    # Main content
    card(
        # Output plot
        plotlyOutput("prediction_plot")
        )
)
