#--------------------------------------------------------------#
# Prediction models for hospitalization in CKD
# Creating functions for development script
# Roemer J. Janse - 2023/11/29
#--------------------------------------------------------------#

# Function for quick save
quick_ggsave <- function(name, plot){ggsave(name, plot, path = paste0(path, "/figures"), width = 5, height = 5, dpi = 600)}

# Function for quiet output
quiet <- function(x){ 
    # Put data in a temporary file
    sink(tempfile()) 
    
    # Sink again on end of function
    on.exit(sink()) 
    
    # Forcefully make output invisible
    invisible(force(x)) 
} 

# Function for linearity check
lin_check <- function(df, var, event, time, annotation, model, xlab){
    # Get event
    if(is.factor(df[[event]])) ev <- as.numeric(filter(df, .imp == 1)[[event]]) - 1 
    else ev <- as.numeric(filter(df, .imp == 1)[[event]])
    
    # Fine-Gray model
    if(model == "fine-gray"){
        # Get data
        dat_plot <- do.call("cbind", quiet(rcspline.plot(x = filter(df, .imp == 1)[[var]], y = filter(df, .imp == 1)[[time]], noprint = TRUE,
                                                         event = ev, model = "cox", statloc = "none"))) %>%
            # Change to data frame
            as.data.frame() 
    }
    
    # Cox and AFT model
    if(model %in% c("cox", "aft")){
        # Get data
        dat_plot <- do.call("cbind", quiet(rcspline.plot(x = filter(df, .imp == 1)[[var]], y = filter(df, .imp == 1)[[time]], noprint = TRUE,
                                                         event = ev == 1, model = "cox", statloc = "none"))) %>%
            # Change to data frame
            as.data.frame() 
    }
    
    # Logistic model
    if(model == "logistic"){
        # Get data
        dat_plot <- do.call("cbind", quiet(rcspline.plot(x = filter(df, .imp == 1)[[var]], y = ev == 1, 
                                                         noprint = TRUE, model = "logistic", statloc = "none"))) %>%
            # Change to data frame
            as.data.frame() 
    }
    
    # Get x-axis label
    xlabel <- xlab
    
    # Get plot
    p <- ggplot(dat_plot, aes(x = x, y = V3)) +
        # Geometries
        geom_ribbon(aes(ymin = V4, ymax = V5), alpha = 0.2) +
        geom_line() +
        # Scaling
        scale_y_continuous(limits = c(min(dat_plot[["V4"]]) - 1, max(dat_plot[["V5"]]) + 1), 
                           name = ifelse(model == "logistic", "Log odds", "Log relative hazard")) +
        # Labelling
        xlab(xlabel) +
        annotate("text", x = as.numeric(annotation[2]), y = as.numeric(annotation[3]), label = annotation[1], fontface = "bold", size = 8) +
        # Aesthetics
        theme_bw() +
        theme(panel.grid = element_blank())
    
    # Return plot
    return(p)
}

# Function for proportional hazards plot
ph_plot <- function(var, ylabel, annotation){
    # Get data for assumption check
    ph_check <- 
        # Bind data together
        cbind(# Time (x-axis)
            cox.zph(fit)[["time"]],
            # Y per variable
            cox.zph(fit)[["y"]]) %>%
        # Change to data frame
        as.data.frame() %>%
        # Remove factor from column names
        set_colnames(gsub("as.factor\\(|\\)", "", colnames(.)))
    
    # Get data
    dat_plot <- ph_check %>%
        # Keep only relevant variables
        dplyr::select(V1, all_of(var)) %>%
        # Rename second variable
        rename(y = 2)
    
    # Draw plot
    p <- ggplot(dat_plot, aes(x = V1, y = y)) +
        # Geometries
        geom_point(alpha = 0.2) +
        geom_smooth(formula = y ~ x, method = "loess", colour = "black") +
        # Scaling
        scale_x_continuous(breaks = seq(0, max(dat_plot[["V1"]]), 100)) +
        # Labels
        xlab("Time (days)") +
        ylab(bquote(beta[t] *" for " * .(ylabel))) +
        annotate("text", x = as.numeric(annotation[2]), y = as.numeric(annotation[3]), label = annotation[1], fontface = "bold", size = 8) +
        # Aesthetics
        theme_bw() +
        theme(panel.grid = element_blank())
    
    # Return plot
    return(p)
}

# Function to create model per imputation
develop <- function(df, imp, formula, model, horizon, aft_dist, all_bhs){
    # Data
    dat_mod_tmp <- filter(df, .imp == imp) 
    
    # Fine-Gray model fit
    if(model == "fine-gray"){
        # Get left-hand-side of formula
        lhs <- str_extract(formula, ".*(?=~)")
        
        # Get right hand side of formula
        rhs <- str_extract(formula, "(?<=~).*")
        
        # Fit model
        fit <- coxph(as.formula(paste0("Surv(fgstart, fgstop, fgstatus) ~ ", rhs)), x = TRUE, 
                     weight = fgwt, data = finegray(as.formula(paste0(lhs, "~ .")), data = dat_mod_tmp))
    }
    
    # Cox model fit
    if(model == "cox") fit <- coxph(as.formula(formula), x = TRUE, data = dat_mod_tmp)
    
    # Accelerated failure time model fit
    if(model == "aft") fit <- survreg(as.formula(formula), data = dat_mod_tmp, dist = aft_dist)

    # Logistic model fit
    if(model == "logistic") fit <- glm(as.formula(formula), family = "binomial", data = dat_mod_tmp)
    
    # If survival model, get baseline hazard and add to coefficients
    if(model %in% c("cox", "fine-gray")){
        # Get baseline hazard for 1 year
        bh <- basehaz(fit) %>%
            # Change to data frame
            as.data.frame() %>%
            # Remove everything after prediction horizon
            filter(time <= horizon) %>%
            # Keep only last observation
            slice_tail(n = 1L) %>%
            # Keep baseline hazard
            extract2("hazard")
        
        # Get coefficients
        coeffs <- fit[["coefficients"]] %>%
            # Change to data frame
            as.data.frame() %>%
            # Transpose
            t() %>%
            # Add baseline hazard
            bind_cols(as.data.frame(bh), .)
        
        # Export all baseline hazards if needed
        if(all_bhs){
            # All times
            times <- 1:horizon
            
            # Get baseline hazards for each timepoint
            bhs <- do.call("c", c(0, lapply(times[2:length(times)], \(x) {
                # All baseline hazards
                basehaz(fit) %>% 
                    # Hazards before timepoint of interest
                    filter(time <= x) %>% 
                    # Keep hazard only
                    extract2("hazard") %>% 
                    # Keep last hazard only
                    last()})))
        }
    }
    
    # If logistic or accelerated failure time model, get intercept
    if(model %in% c("logistic", "aft")){
        # Get coefficients
        coeffs <- fit[["coefficients"]] %>%
            # Change to data frame
            as.data.frame() %>%
            # Transpose
            t() %>%
            # Rename intercept
            set_colnames(c("a", colnames(.)[2:length(colnames(.))]))
    }
    
    # Return model data and bhs if asked for
    if(!all_bhs) return(coeffs) else return(list(coeffs, bhs))
}

# Function for model development
dev <- function(df, formula, model, horizon, aft_dist = NULL, all_bhs = FALSE){
    # Fit model for each imputation and add together the final models as normal if not all baseline hazards are requested
    if(!all_bhs){
        mod <- do.call("rbind", lapply(unique(df[[".imp"]]), \(x) develop(df, x, formula, model, horizon, aft_dist, all_bhs)))
    }
    
    # Else store list and then bind models
    else {
        # Individual models
        models <- lapply(unique(df[[".imp"]]), \(x) develop(df, x, formula, model, horizon, aft_dist, all_bhs))
        
        # Bind models
        mod <- do.call("rbind", lapply(unique(df[[".imp"]]), \(x) models[[x]][[1]]))
        
        # Bind baseline hazards
        all_bh <- do.call("rbind", lapply(unique(df[[".imp"]]), \(x) models[[x]][[2]])) %>%
            # Change to data frame
            as.data.frame()
    }
    
    # Get final model fit by taking means of each model
    model_vars <- summary(mod) %>% 
        # Change to data frame
        as.data.frame() %>% 
        # Keep only mean
        filter(grepl("Mean", Freq)) %>% 
        # Change parameters to numeric
        mutate(param = as.numeric(gsub("\\D{4} *:", "", Freq))) %>%
        # Rename Var2
        rename(var = Var2) %>%
        # Keep only columns of interest
        dplyr::select(var, param) %>%
        # Transpose
        t() %>%
        # Set column names
        set_colnames(gsub(" ", "", .[1, ])) %>%
        # Change to data frame
        as.data.frame() %>%
        # Keep only second row
        filter(row_number() == 2) 
    
    # Get spline information if splines were used
    if(str_detect(formula, "ns\\(")){
        # Get all variables for which splines were used
        model_spline_info <- tibble(variables = str_extract_all(formula, "(?<=ns\\().*?(?=,)")[[1]]) 
        
        # Empty boundaries
        lower_boundaries <- c(); upper_boundaries <- c()
        
        # Get boundaries
        for(i in 1:nrow(model_spline_info)){
            # Lower boundary
            lower_boundaries[i] <- min(df[[model_spline_info[[i, 1]]]])
            
            # Upper boundary
            upper_boundaries[i] <- max(df[[model_spline_info[[i, 1]]]])
        }
        
        # Add information to data and put data in global environment
        model_spline_info <<- mutate(model_spline_info, lower_boundary = lower_boundaries, upper_boundary = upper_boundaries)
    }
    
    # Return average baseline hazards if requested
    if(all_bhs){
        # Store baseline hazards globally
        bhs <<- summary(all_bh) %>% 
            # Change to data frame
            as.data.frame() %>% 
            # Keep only mean
            filter(grepl("Mean", Freq)) %>% 
            # Change parameters to numeric
            mutate(bh = as.numeric(gsub("\\D{4} *:", "", Freq))) %>%
            # Add timepoints
            mutate(time = 1:horizon) %>%
            # Keep variables of interest
            select(time, bh)
    }
    
    # Return model vars
    return(model_vars)
}

# C-statistic function
cstat <- function(df, model){
    # Wolbers' adaptation for competing events
    if(model == "fine-gray") dat_tmp <- mutate(df, tim = ifelse(obs >= 2, Inf, tim)) else dat_tmp <- df
    
    # For logistic models
    if(model == "logistic"){
        # All predictions of individuals with event
        events <- dat_tmp %>%
            # Keep only events
            filter(obs == 1) %>%
            # Extract predictions
            extract2("prd")
        
        # All predictions of individuals without event
        no_events <- dat_tmp %>%
            # Keep only events
            filter(obs == 0) %>%
            # Extract predictions
            extract2("prd")
        
        # Compare all events with all non-events and calculate C-statistic
        cstatistic <- sum((tmp <- do.call("c", lapply(events, \(x) ifelse(x > no_events, 1, ifelse(x == no_events, 0.5, 0)))))) / length(tmp); rm(tmp)
    }
    
    # For survival models (Harrell's C-statistic)
    if(model %in% c("cox", "fine-gray")) cstatistic <- cIndex(dat_tmp[["tim"]], dat_tmp[["obs_ncr"]], dat_tmp[["prd"]])[["index"]]
    
    # For AFT models, use aft_status obs
    if(model == "aft") cstatistic <- cIndex(dat_tmp[["tim"]], dat_tmp[["aft_status"]], dat_tmp[["prd"]])[["index"]]
    
    # Return C-statistic
    return(cstatistic)
}

# Function for validation
validate <- function(.data,                                     # Data
                     observed,                                  # Observed outcome
                     predicted,                                 # Predicted outcome
                     lp,                                        # Linear predictor
                     model,                                     # Regression model used to create the prediction model,
                     time = NULL,                               # Time variable (only relevant for Cox/AFT/Fine-Gray)
                     aft_dist = "lognormal",                    # Distribution used for the AFT model
                     # Output
                     print_stratified = FALSE,                  # Print stratified C-statistics (only relevant for survival model)
                     plot = TRUE,                               # Should a calibration plot be made
                     deciles = TRUE,                            # Should deciles be added in calibration plot
                     bootstraps = 500,                          # Number of bootstraps for C-statistic
                     # Calibration plot details
                     unit = "probability",                      # Unit of prediction for axes of plot
                     annotation = c("", 0, 1),                  # Annotation to add to plot as c("annotation", x, y)
                     smoother = TRUE,                           # Should a smoother be added? This also determines whether pseudo-obs are calculated in FG models
                     smooth_colour = "darkred",                 # Colour of smoother
                     histogram_label = NULL                     # Location of event / no-event label in probabilities histogram
){
    ## Prepare data
    # Observed
    obs <- .data[[deparse(substitute(observed))]]
    
    # Set observed to numeric if it is factor
    if(is.factor(obs)) obs <- as.numeric(obs) - 1
    
    # Observed without competing risks
    if(model == "fine-gray") obs_ncr <- ifelse(obs == 1, 1, 0) else obs_ncr <- obs
    
    # Predicted
    prd <- .data[[deparse(substitute(predicted))]]
    
    # Linear predictors
    lps <- .data[[deparse(substitute(lp))]] 
    
    # Time
    if(!(deparse(substitute(time)) == "NULL")) tim <- .data[[deparse(substitute(time))]] else tim <- NA
    
    # If model is AFT, observed should be time to event and status should be stored elsewhere
    if(model == "aft") {
        # Store status elsewhere
        aft_status <- obs
        
        # Store EFT
        obs <- obs_ncr <- tim
    }
    
    # Else empty aft_status variable
    else aft_status <- NA
    
    # Create data to work with
    dat <- tibble(obs = obs,
                  obs_ncr = obs_ncr,
                  lps = lps,
                  prd = prd, 
                  tim = tim,
                  aft_status = aft_status)
    
    # Get total number of individuals
    n <- nrow(dat)
    
    ## If plotting, prepare data and make plot
    if(plot){
        # Calculate deciles
        if(deciles){
            # Create deciles data
            dat_dec <- dat %>%
                # Create decile
                mutate(dec = cut(prd, breaks = as.numeric(quantile(prd, probs = seq(0, 1, 0.1))), include.lowest = TRUE)) %>%
                # Sort for grouping
                arrange(dec) %>%
                # Group per decile
                group_by(dec) %>%
                # Get mean outcome and probability
                mutate(# Outcome
                       out_prop = mean(obs_ncr),
                       # Predicted
                       pred_prop = mean(prd),
                       # Check number of individuals
                       nmax = max(row_number())) %>%
                # Keep one row per decile
                slice(1L) %>%
                # Ungroup again
                ungroup() %>%
                # Keep only relevant columns
                select(dec, out_prop, pred_prop, nmax)
        }
        
        # Calculate pseudo-observations if smoother
        if(model == "fine-gray" && smoother){
            # Create empty id column
            dat <- mutate(dat, ids = 1:nrow(dat))
            
            # Get A-J estimate for outcome in total population
            aj_tot <- survfit(Surv(tim, obs, type = "mstate") ~ 1, data = dat) %>%
                # Change fit list to dataframe
                tidy() %>%
                # Keep only events
                filter(state == 1) %>%
                # Keep only last observation, assuming this is the time point of interest
                slice_tail(n = 1L) %>%
                # Keep baseline hazard
                extract2("estimate")
            
            # Get A-J estimates for excluding each individual with jackknife
            aj_in <- do.call("c", lapply(dat[["ids"]], \(x){
                # Create new data
                dat_tmp <- filter(dat, ids != x)
                
                # Calculate A-J estimate
                aj <- survfit(Surv(tim, obs, type = "mstate") ~ 1, data = dat_tmp) %>%
                    # Change fit list to dataframe
                    tidy() %>%
                    # Keep only events
                    filter(state == 1) %>%
                    # Keep only last observation, assuming this is the time point of interest
                    slice_tail(n = 1L) %>%
                    # Keep baseline hazard
                    extract2("estimate")
                
                # Return estimate
                return(aj)
            }))
            
            # Get jackknife A-J estimate per individual
            dat <- dat %>%
                # Calculate new variables
                mutate(# Total A-J estimate
                    aj_t = aj_tot,
                    # Jackknife A-J estimate
                    aj_i = aj_in,
                    # Pseudovalue for the outcome
                    aj_o = (n * aj_t) - ((n - 1) * aj_i))
        }
        
        # Set y variable based on model
        if(model == "fine-gray" && smoother) dat <- mutate(dat, y = aj_o) else dat <- mutate(dat, y = obs_ncr)
        
        ## Define characteristics of the plot
        # Upper limits of axes
        if(model %in% c("poisson", "linear", "aft")) xmax <- ymax <- max(obs) else xmax <- ymax <- 1
        
        # Lower limits of axes
        if(model %in% c("poisson", "linear", "aft")) xmin <- ymin <- pmin(0, min(obs)) else xmin <- ymin <- 0
        
        # Breaks of axes
        brks <- seq(xmin, xmax, (xmax - xmin) / 10)
        
        # Labels of x-axis
        if(model %in% c("poisson", "linear", "aft")) xlab <- paste0("Predicted ", unit) else xlab <- "Predicted probability"
        
        # Label of y-axis
        if(model %in% c("poisson", "linear", "aft")) ylab <- paste0("Observed ", unit) else ylab <- "Observed probability"
        
        # Determine location of event label in probabilities histogram automatically if not specified, based on least predicted probability
        if(is.null(histogram_label)) histogram_label <- as.numeric(names(sort(table(round(dat[["prd"]], 3))))[[1]])
        
        ## Create plot
        # Create base scatter plot
        plot_cal <- ggplot(dat, aes(x = prd, y = y)) +
            # Geometries
            geom_abline(colour = "black", linewidth = 2, alpha = 0.33) +
            geom_point(alpha = 0.25)

        # If smoother, add smoother
        if(smoother) plot_cal <- plot_cal + geom_smooth(colour = smooth_colour, fill = smooth_colour, method = "loess", formula = y ~ x)
        
        # If AFT model, overwrite base scatter plot to add colouring for different statuses
        if(model == "aft"){
            # Create base scatter plot
            plot_cal <- ggplot(dat, aes(x = prd, y = y)) +
                # Geometries
                geom_abline(colour = "black", linewidth = 2, alpha = 0.33) +
                geom_point(alpha = 0.25, mapping = aes(colour = factor(aft_status, levels = c(0, 1), labels = c("No-event", "Event")))) +
                geom_smooth(colour = smooth_colour, fill = smooth_colour, method = "loess", formula = y ~ x) + 
                # Scaling
                scale_colour_manual(values = c("darkorange", "darkred")) +
                # Aesthetics
                theme(legend.position = "bottom",
                      legend.title = element_blank())
        }
        
        # Add deciles
        if(deciles) plot_cal <- plot_cal + geom_point(inherit.aes = FALSE, data = dat_dec, aes(x = pred_prop, y = out_prop), shape = 0)
        
        # Finish calibration plot
        plot_cal <- plot_cal +
            # Scaling
            scale_y_continuous(breaks = brks, name = ylab) +
            scale_x_continuous(breaks = brks, name = xlab) +
            # Transformations
            coord_cartesian(xlim = c(xmin, xmax), ylim = c(ymin, ymax)) +
            # Labels
            annotate("text", x = as.numeric(annotation[[2]]), y = as.numeric(annotation[[3]]), 
                     label = annotation[[1]], fontface = "bold", size = 10) +
            # Aesthetics
            theme(panel.background = element_blank(),
                  panel.grid = element_blank(),
                  axis.line = element_blank()) +
            panel_border(colour = "black", size = 1) 
        
        # Create probability histogram for non-continuous models
        if(!(model %in% c("poisson", "linear", "aft"))){
            # Turn off x-axis for calibration plot
            plot_cal <- plot_cal + theme(axis.ticks.x = element_blank(),
                                         axis.text.x = element_blank(),
                                         axis.title.x = element_blank())
            
            # Create histogram of predicted probabilities
            plot_his <- ggplot(data = dat) +
                # Geometries
                geom_histogram(data = filter(dat, obs_ncr == 1), aes(x = prd), binwidth = 0.001, colour = "black", fill = "black") +
                geom_histogram(data = filter(dat, obs_ncr == 0), aes(x = prd, y = -after_stat(count)), binwidth = 0.001, colour = "black", 
                               fill = "black") +
                annotate("text", x = histogram_label, y = max(table(round(dat[["prd"]], 3))) / 2, label = "Event", hjust = 0) +
                annotate("text", x = histogram_label, y = -max(table(round(dat[["prd"]], 3))) / 2, label = "No event", hjust = 0) +
                # Scaling
                scale_x_continuous(breaks = brks, name = xlab, limits = c(xmin, xmax)) +
                # Transformations
                coord_cartesian(xlim = c(xmin, xmax), ylim = c(-max(table(round(dat[["prd"]], 3))), max(table(round(dat[["prd"]], 3))))) +
                # Aesthetics
                theme(panel.background = element_blank(),
                      panel.grid = element_blank(),
                      axis.line = element_blank(),
                      axis.ticks.y = element_blank(),
                      axis.text.y = element_blank(),
                      axis.title.y = element_blank()) +
                panel_border(colour = "black", size = 1)
            
            # Combine plots
            plot_cal <- suppressWarnings(plot_grid(plot_cal, plot_his, align = c("hv"), axis = c("tblr"), ncol = 1, rel_heights = c(3, 1)))
        }
    }
    
    ## Calculate performance measures
    # Get proportion of outcome for non-continuous models
    if(!(model %in% c("linear", "poisson", "aft"))) prop_out <- prop.table(table(dat[["obs"]]))[["1"]]
    
    # Get mean of outcome for continous models
    else prop_out <- mean(obs)
    
    # Calibration-in-the-large
    citl <- format(round(prop_out / mean(prd), 3), nsmall = 3)
    
    ## Calculate calibration slope
    # Generalized linear model
    if(!(model %in% c("cox", "fine-gray", "aft"))) cslope <- format(round(glm(y ~ lps, family = ifelse(model == "logistic", "binomial", model), 
                                                                            data = dat)[["coefficients"]][["lps"]], 3), nsmall = 3)
    
    # Cox model
    if(model == "cox") cslope <- format(round(coxph(Surv(tim, obs) ~ lps, data = dat)[["coefficients"]][["lps"]], 3), nsmall = 3)
    
    # Fine-Gray model
    if(model == "fine-gray") cslope <- format(round(coxph(Surv(fgstart, fgstop, fgstatus) ~ lps, weight = fgwt, 
                                                          data = finegray(Surv(tim, as.factor(obs)) ~ ., data = dat))[["coefficients"]][["lps"]], 
                                                    3), nsmall = 3)
    
    # AFT model
    if(model == "aft") cslope <- format(round(survreg(Surv(tim, aft_status) ~ lps, data = dat, dist = aft_dist)[["coefficients"]][["lps"]], 3), nsmall = 3)

    # C-statistic
    if(!(model %in% c("linear", "poisson"))){
        # Calculate C statistic
        c <- cstat(dat, model)
        
        # Calculate confidence interval around C statistic
        ci <- quantile(do.call("c", lapply(1:bootstraps, \(x){
            # Set seed
            set.seed(x)
            
            # Get random sample numbers
            nrs <- sample.int(length(prd), replace = TRUE)
            
            # Get samples
            dat_samp <- dat[nrs, ] %>%
                # Change studynr
                mutate(studynr = 1:nrow(dat))
            
            # Calculate statistic
            c_bootstrap <- cstat(dat_samp, model)
            
            # Return statistic
            return(c_bootstrap)
        })), probs = c(0.025, 0.975))
        
        # Pretty up C statistic
        c <- paste0(format(round(c, 2), nsmall = 2), " (",
                    format(round(ci[[1]], 2), nsmall = 2), "-",
                    format(round(ci[[2]], 2), nsmall = 2), ")")
    }
    
    # List with output
    output <- list(
        # Calibration statistics
        tibble("Calibration-in-the-large" = citl, "Calibration slope" = cslope)
    )
    
    # If C-statistic was calculated, add it
    if(!(model %in% c("linear", "poisson"))) output[[1]] <- mutate(output[[1]], "C statistic" = c)
    
    # If plot was desired, add it
    if(plot) output[[2]] <- plot_cal
    
    # Return final result
    return(output)
}

# Create function to get model information
model_inf <- function(mod){
    # Get model
    model_vars <- mod
    
    # Get factor info if factors were used
    if(length(names(model_vars[grepl("factor\\(", names(model_vars))])) > 0){
        # Get information on factor variables
        factor_info <- names(model_vars[grepl("factor\\(", names(model_vars))]) %>%
            # Change into data frame
            as.data.frame() %>%
            # Rename first column
            rename(input = 1) %>%
            # Get variable
            mutate(var = str_replace_all(input, "(as.)+factor\\(|\\)\\d*", "")) %>%
            # Group per variable
            group_by(var) %>%
            # Create new columns
            mutate(
                # Levels
                levels = as.numeric(str_extract(input, "(?<=\\))\\d+")),
                # Rows per group
                rows = max(row_number())) %>%
            # Ungroup again
            ungroup()
        
        # Add information to factor_info
        factor_info <- mutate(factor_info, 
                              # Coefficient per level
                              coef = as.numeric(as.vector(model_vars[grepl("factor", names(model_vars))])),
                              # Add type
                              type = "factor") %>%
            # Keep only relevant variables
            select(var, levels, coef, type) %>%
            # Create empty lower and upper level columns
            mutate(lower_level = -Inf, upper_level = Inf)
    }
    
    # Else, create empty data
    else factor_info <- tibble(var = as.character(),
                               coef = as.numeric(),
                               lower_level = as.numeric(), 
                               upper_level = as.numeric(),
                               type = as.character())
    
    # Get spline info if B-spline matrix basis splines were used
    if(length(names(model_vars[grepl("knot", names(model_vars))])) > 0){
        # Get information on splined variables
        bspline <- names(model_vars[grepl("knot", names(model_vars))]) %>%
            # Change into data frame
            as.data.frame() %>%
            # Rename first column
            rename(input = 1) %>%
            # Get variable
            mutate(var = str_extract(input, "(?<=.{1,5}\\().*(?=,.?knots=)")) %>%
            # Group per variable
            group_by(var) %>%
            # Create new columns
            mutate(
                # Knots
                knots = paste0("-Inf,", str_extract(input, "(?<=[\\(|=|,])[\\d.,]+(?=\\))"), ",Inf"),
                # Levels
                levels = as.numeric(str_extract(input, "(?<=\\))\\d+")),
                # Row number
                row = row_number()) %>%
            # Ungroup again
            ungroup()
        
        # Matrix of spline levels
        spline_levels <- str_split(bspline[["knots"]], ",", simplify = TRUE)
        
        # Get vector with lower levels
        lower_levels <- as.numeric(lapply(1:nrow(bspline), \(x) spline_levels[x, bspline[[x, "levels"]]]))
        
        # Get vector with upper levels
        upper_levels <- as.numeric(lapply(1:nrow(bspline), \(x) spline_levels[x, bspline[[x, "levels"]] + 1]))
        
        # Add information to spline_info
        transspline_info <- mutate(bspline, 
                                   # Lower level of section
                                   lower_level = lower_levels,
                                   # Upper level of section
                                   upper_level = upper_levels,
                                   # Coefficient per level
                                   coef = as.numeric(as.vector(model_vars[grepl("knot", names(model_vars))])),
                                   # Add type
                                   type = "bspline") %>%
            # Keep only relevant variables
            select(var, lower_level, upper_level, coef, type, levels) 
    }
    
    # Else, create empty data
    else transspline_info <- tibble(var = as.character(),
                                    coef = as.numeric(),
                                    lower_level = as.numeric(), 
                                    upper_level = as.numeric(),
                                    type = as.character(),
                                    levels = NA)
    
    # Create data for all variables in the model, adding onto the splined variables
    model_info <- tibble(var = names(model_vars[!grepl("knot|factor\\(", names(model_vars))][-1]),
                         coef = as.numeric(t(model_vars[!grepl("knot|factor\\(", names(model_vars))][-1])),
                         lower_level = -Inf, 
                         upper_level = Inf,
                         levels = NA,
                         type = "as_is") %>%
        # Clean up names
        mutate(var = str_replace_all(var, ".*(?<!,knot=c)\\(|\\)\\d", "")) %>%
        # Add splined and factor variables
        bind_rows(transspline_info, factor_info)
    
    # Return data
    return(model_info)
}

# Function for B-spline basis matrix transformation
bspline_trans <- function(variable, value, df = 3, mod = model_vars){
    # Get model info
    model_info <- model_inf(mod)
    
    # Get knots
    knots <- c(filter(model_info, var == variable)[["lower_level"]], 
               filter(model_info, var == variable)[["upper_level"]]) %>%
        # Keep only unique
        unique() %>%
        # Keep only not infinites
        `[`(!is.infinite(.))
    
    # Lower boundary
    lb <- filter(model_spline_info, variables == variable)[["lower_boundary"]]
    
    # Upper boundary
    ub <- filter(model_spline_info, variables == variable)[["upper_boundary"]]
    
    # Get B-spline matrix
    mat <- as.data.frame(ns(value, knots = knots, df = df, Boundary.knots = c(lb, ub)))
    
    # Return matrix
    return(mat)
}
    

# Create function to calculate sample linear predictor
lpsamp <- function(df, mod = model_vars){
    # Store data
    dat_tmp <- df
    
    # Get model information
    model_info <- model_inf(mod)
    
    # Calculate LP fractions
    lp_fracs <- bind_cols(lapply(1:nrow(model_info), \(x){
        # Get variable of interest
        var <- model_info[x, "var"][[1]]

        # Get coefficient of interest
        coef <- model_info[x, "coef"][[1]]
        
        # Get lower level of section
        lower_level <- model_info[x, "lower_level"][[1]]
        
        # Get upper level of section
        upper_level <- model_info[x, "upper_level"][[1]]
        
        # Get type of variable
        type <- model_info[x, "type"][[1]]
        
        # Get level of variable
        level <- model_info[x, "levels"][[1]]
        
        # Change var if transforming spline
        if(type == "bspline") {
            # Update var
            var_new <- paste0(var, level)
            
            # Add new var in data
            dat_tmp[[var_new]] <- bspline_trans(var, dat_tmp[[var]])[, level]
            
            # Overwrite var
            var <- var_new
        }
        
        # Determine whether variable is logical
        if(length(table(dat_tmp[[var]])) == 2 | type == "factor") logic <- TRUE else logic <- FALSE
        
        # Calculate fraction of LP
        lp_frac <- dat_tmp %>%
            # Select relevant variables
            select(.imp, all_of(var)) %>%
            # Rename varying 'var' to value
            rename(value = 2) %>%
            # Arrange for grouping
            arrange(.imp) %>%
            # Group per imputation
            group_by(.imp) %>%
            # Calculate new variables
            mutate(# Set value based on conditions
                value = case_when(type %in% c("as_is", "bspline") ~ value,                                # No changes for as_is and spline variables
                                  type == "factor" & as.character(value) == as.character(level) ~ 1,      # Factor to 1 if same level (this represents dummy variable)
                                  type == "factor" & as.character(value) != as.character(level) ~ 0,      # Factor to 0 if not the same level
                                  .default = NA),                                                         # Otherwise, set value to missing
                # Get mean of value
                mean_value = mean(value, na.rm = TRUE),
                # For factors, mean value is 0 (see details coxph)
                mean_value = ifelse(logic, 0, mean_value),
                # Multiply mean value with coefficient
                lp = mean_value * coef) %>%
            # Keep one row per group
            slice(1L) %>%
            # Ungroup again
            ungroup() %>%
            # Select only relevant variables
            select(lp) %>%
            # Rename to make name unique
            set_colnames(paste0("frac_", x))
        
        # Return data
        return(lp_frac)
    })) 
    
    # Get and return linear predictor of sample
    return(mean(rowSums(lp_fracs, na.rm = TRUE)))
}

# Create function to get predicted risks
pred <- function(df, model, observed, time = NULL, lpsamp = NULL, aft_dist = NULL, mod = model_vars){
    # Store data
    dat_tmp <- df
    
    # Get model information
    model_info <- model_inf(mod)
    
    # Calculate LP fractions
    lp_fracs <- bind_cols(lapply(1:nrow(model_info), \(x){
        # Get variable of interest
        var <- model_info[x, "var"][[1]]
        
        # Get coefficient of interest
        coef <- model_info[x, "coef"][[1]]
        
        # Get lower level of section
        lower_level <- model_info[x, "lower_level"][[1]]
        
        # Get upper level of section
        upper_level <- model_info[x, "upper_level"][[1]]
        
        # Get type of variable
        type <- model_info[x, "type"][[1]]
        
        # Get level of variable
        level <- model_info[x, "levels"][[1]]
        
        # Change var if spline
        if(type == "bspline") {
            # Update var
            var_new <- paste0(var, level)
            
            # Add new var in data
            dat_tmp[[var_new]] <- bspline_trans(var, dat_tmp[[var]])[, level]
            
            # Overwrite var
            var <- var_new
        }
        
        # Calculate fraction of LP
        lp_frac <- dat_tmp %>%
            # Select relevant variables
            select(studynr, .imp, all_of(var)) %>%
            # Rename varying 'var' to value
            rename(value = 3) %>%
            # Arrange for grouping
            arrange(studynr, .imp) %>%
            # Group per imputation
            group_by(studynr, .imp) %>%
            # Calculate fraction
            mutate(lp_frac = case_when(type %in% c("as_is", "bspline") ~ value * coef,                                # as_is & B-spline
                                       type == "factor" & as.character(value) == as.character(level) ~ coef,          # Factors
                                       .default = 0)) %>%  
            # Ungroup again
            ungroup() %>%
            # Keep only LP fraction
            select(lp_frac) %>%
            # Rename to make name unique
            set_colnames(paste0("frac_", x))
        
        # Return data
        return(lp_frac)
    })) 
    
    # Get observed event if given
    if(deparse(substitute(observed)) != "NULL") obs <- arrange(df, studynr, .imp)[[deparse(substitute(observed))]] else obs <- NA
    
    # Get time (for survival models)
    tim <- arrange(df, studynr, .imp)[[deparse(substitute(time))]]
    
    # Sum linear predictors
    lp_sums <- do.call("c", lapply(1:nrow(lp_fracs), \(x) sum(t(lp_fracs[x, ]), na.rm = TRUE)))
    
    # Get final linear predictors
    dat_tmp <- df %>%
        # Arrange on imputation
        arrange(studynr, .imp) %>%
        # Keep only imputation
        select(studynr, .imp) %>%
        # Add variables of interest
        mutate(# Linear predictor
            lp = lp_sums,
            # Observed event
            observed = obs,
            # Time (for survival models)
            time = tim)
    
    # For Cox/FIne-Gray models
    if(model %in% c("cox", "fine-gray")){
        ## Calculate individual risks
        # The predicted risks from predict.coxph() from {survival} do not correspond to our calculated probabilities, because they are
        # the survival at the observed follow-up time for that individual. Instead, predictCox() from {riskRegression} gives predictions
        # at a certain timepoint: predictCox(fit, newdata = df, times = 365, type = "survival")[["survival"]]
        # https://community.rstudio.com/t/predict-probability-of-default-in-cox-ph-model-and-get-baseline-hazard/118953 
        dat_tmp <- dat_tmp %>%
            # Calculate new variables
            mutate(# Centered linear predictors
                lp = lp - lpsamp,
                # Predicted probability
                prob = 1 - exp(-as.numeric(model_vars[["bh"]])) ^ exp(lp)) %>%
            # Ungroup again
            ungroup()
    }
    
    # For AFT models (currently only for lognormal and Weibull distributions)
    if(model == "aft"){
        # Calculate individual times
        dat_tmp <- dat_tmp %>%
            # Calculate expected time to failure
            mutate(prob = survreg.distributions[[aft_dist]][["itrans"]](as.numeric(model_vars[["a"]]) + lp))
    }
    
    # For logistic models
    if(model == "logistic"){
        # Calculate individual risks
        dat_tmp <- dat_tmp %>%
            # Calculate predicted probability
            mutate(prob = 1 / (1 + exp(-(as.numeric(model_vars[["a"]]) + lp))))
    }
    
    # Take mean of predicted value
    dat_tmp <- dat_tmp %>%
        # Sort for grouping
        arrange(studynr) %>%
        # Group on person
        group_by(studynr) %>%
        # Calculate final predicted probability
        mutate(pred = mean(prob)) %>%
        # Keep one row per person
        slice(1L) %>%
        # Ungroup again
        ungroup() %>%
        # Drop irrelevant columns
        dplyr::select(-prob, -.imp)
    
    # Return data
    return(dat_tmp)
}

# Function for fitted model
mv <- function(){
    # Create model vars based on dput()
    structure(list(bh = " 0.1313000", female = "-0.0317540", `ns(age,knots=60)1` = " 1.1301000", 
                   `ns(age,knots=60)2` = " 0.7079000", education = "-0.0584500", 
                   diab_months = " 0.0007586", `ns(egfr,knots=105)1` = " 0.8377000", 
                   `ns(egfr,knots=105)2` = " 2.2880000", `ns(ucrea,knots=180)1` = " 0.0240500", 
                   `ns(ucrea,knots=180)2` = " 0.7989000", `ns(hba1c,knots=c(40,65))1` = " 0.8660000", 
                   `ns(hba1c,knots=c(40,65))2` = " 0.4963000", `ns(hba1c,knots=c(40,65))3` = "-0.3145000", 
                   `ns(total_cholesterol,knots=100)1` = "-0.3122600", `ns(total_cholesterol,knots=100)2` = "-0.2978500", 
                   hdl = "-0.0045870", `ns(ldl,knots=80)1` = " 0.3539000", `ns(ldl,knots=80)2` = " 0.6085000", 
                   triglycerides = " 0.0004498", `ns(tc_hdl,knots=5)1` = "-0.7897600", 
                   `ns(tc_hdl,knots=5)2` = "-1.0283000", `ns(index_alb,knots=15)1` = " 1.9680000", 
                   `ns(index_alb,knots=15)2` = " 0.5054000", fibrillation = " 0.1761000", 
                   chf = " 0.0841500", cvd = " 0.0104780", hypertension = " 0.0705500", 
                   ihd = "-0.0649900", neuropathy = " 0.5330000", pvd = " 0.1403000", 
                   retinopathy = "-0.0651600", ulcer = " 0.0739400", aspirin = "-0.0441700", 
                   bblockers = " 0.1307000", glucose_lowering = " 0.1220000", 
                   ccbs = " 0.1177000", anticoagulants = " 0.1477000", antihypertensives = " 0.2271000", 
                   insulin = " 0.1510000", mras = " 0.0089220", rasi = " 0.0577200", 
                   statins = "-0.0241210"), class = "data.frame", row.names = "param")
}

# Function for baseline hazards of fitted model
bh <- function(){
    # Create model baseline hazards per day (from dput())
    structure(list(time = 1:1095, bh = c(0, 0.0005211, 0.0006053, 
                                         0.0007735, 0.0009587, 0.0009924, 0.001093, 0.001211, 0.00138, 
                                         0.001498, 0.001548, 0.00165, 0.001751, 0.001852, 0.001869, 0.001953, 
                                         0.002004, 0.002088, 0.002139, 0.002291, 0.002477, 0.002544, 0.002679, 
                                         0.00273, 0.002798, 0.002882, 0.003001, 0.003119, 0.003187, 0.003237, 
                                         0.003322, 0.003356, 0.003441, 0.003491, 0.003678, 0.003813, 0.003847, 
                                         0.003915, 0.003966, 0.003983, 0.004118, 0.004186, 0.004288, 0.004322, 
                                         0.004406, 0.004542, 0.004678, 0.004746, 0.004831, 0.004933, 0.005051, 
                                         0.005119, 0.005204, 0.005255, 0.005357, 0.005459, 0.005527, 0.005578, 
                                         0.005646, 0.005714, 0.005748, 0.005884, 0.006054, 0.006191, 0.006225, 
                                         0.00631, 0.006378, 0.006497, 0.006616, 0.006719, 0.006804, 0.006838, 
                                         0.006855, 0.006889, 0.006923, 0.006991, 0.007094, 0.007145, 0.007247, 
                                         0.007332, 0.007401, 0.007469, 0.007588, 0.007657, 0.007708, 0.007776, 
                                         0.007776, 0.007879, 0.007981, 0.008118, 0.008254, 0.00834, 0.008425, 
                                         0.008545, 0.008562, 0.008579, 0.008596, 0.008767, 0.008836, 0.008956, 
                                         0.009058, 0.00911, 0.009212, 0.009332, 0.009521, 0.009572, 0.009641, 
                                         0.009743, 0.009863, 0.009881, 0.009966, 0.010104, 0.01036, 0.01043, 
                                         0.01048, 0.01055, 0.01065, 0.01077, 0.01088, 0.01094, 0.01098, 
                                         0.01101, 0.01112, 0.01122, 0.01139, 0.01153, 0.01163, 0.01172, 
                                         0.01177, 0.01189, 0.01199, 0.01211, 0.01225, 0.01234, 0.01242, 
                                         0.01246, 0.01253, 0.01258, 0.0127, 0.01282, 0.01291, 0.01301, 
                                         0.01304, 0.01318, 0.01323, 0.01336, 0.01353, 0.01363, 0.01375, 
                                         0.01384, 0.01389, 0.01401, 0.0141, 0.01431, 0.01441, 0.01446, 
                                         0.01463, 0.01472, 0.01484, 0.015, 0.01517, 0.01527, 0.01536, 
                                         0.01541, 0.01543, 0.01553, 0.01567, 0.01585, 0.016, 0.01606, 
                                         0.01612, 0.01625, 0.01635, 0.01645, 0.01673, 0.01691, 0.01696, 
                                         0.01711, 0.01732, 0.01738, 0.01762, 0.01781, 0.01795, 0.01807, 
                                         0.01818, 0.01833, 0.01842, 0.01847, 0.01865, 0.01877, 0.01882, 
                                         0.01887, 0.01896, 0.0191, 0.01917, 0.01929, 0.01943, 0.01952, 
                                         0.01954, 0.01961, 0.01971, 0.01982, 0.01997, 0.02006, 0.0201, 
                                         0.02015, 0.02018, 0.02027, 0.02039, 0.02064, 0.0208, 0.02095, 
                                         0.02102, 0.02108, 0.02115, 0.02127, 0.02139, 0.02146, 0.02158, 
                                         0.02172, 0.0219, 0.02202, 0.02206, 0.02222, 0.0223, 0.02236, 
                                         0.02241, 0.0225, 0.02265, 0.02271, 0.02288, 0.02299, 0.02302, 
                                         0.02315, 0.02327, 0.02336, 0.02343, 0.02355, 0.02367, 0.02378, 
                                         0.02385, 0.02401, 0.02408, 0.02417, 0.02431, 0.0245, 0.02459, 
                                         0.02464, 0.02473, 0.02479, 0.02493, 0.02517, 0.02528, 0.02535, 
                                         0.02544, 0.02555, 0.02563, 0.02572, 0.02599, 0.02611, 0.02615, 
                                         0.0262, 0.02625, 0.02629, 0.02639, 0.02655, 0.02673, 0.02675, 
                                         0.02678, 0.02689, 0.02698, 0.02703, 0.02726, 0.02737, 0.02751, 
                                         0.02756, 0.02765, 0.02776, 0.02787, 0.02799, 0.02808, 0.0282, 
                                         0.02826, 0.02836, 0.02849, 0.02861, 0.02874, 0.02897, 0.02902, 
                                         0.02913, 0.02918, 0.02927, 0.02938, 0.02948, 0.02966, 0.02972, 
                                         0.02982, 0.02989, 0.03007, 0.03025, 0.03041, 0.0305, 0.03066, 
                                         0.03073, 0.03082, 0.03091, 0.03096, 0.03118, 0.03138, 0.03148, 
                                         0.03155, 0.03161, 0.03166, 0.03182, 0.03209, 0.03222, 0.03238, 
                                         0.03247, 0.03263, 0.03274, 0.03284, 0.03299, 0.03317, 0.03324, 
                                         0.03333, 0.03338, 0.03347, 0.0337, 0.03401, 0.03419, 0.03431, 
                                         0.03446, 0.03457, 0.03464, 0.03478, 0.035, 0.03527, 0.0355, 0.03559, 
                                         0.03566, 0.03577, 0.03595, 0.03626, 0.03649, 0.03662, 0.03673, 
                                         0.03693, 0.03704, 0.03729, 0.03758, 0.03778, 0.03792, 0.03798, 
                                         0.03823, 0.03848, 0.03876, 0.03914, 0.03941, 0.03963, 0.03977, 
                                         0.03985, 0.0401, 0.04048, 0.04072, 0.04103, 0.04132, 0.0415, 
                                         0.04167, 0.04187, 0.04214, 0.04252, 0.04284, 0.04307, 0.04322, 
                                         0.04346, 0.04368, 0.04399, 0.0443, 0.04448, 0.04465, 0.04474, 
                                         0.04487, 0.04516, 0.04538, 0.0458, 0.04606, 0.04619, 0.04636, 
                                         0.04656, 0.04674, 0.04695, 0.04713, 0.04735, 0.04741, 0.04755, 
                                         0.0477, 0.04787, 0.04809, 0.04835, 0.04844, 0.04853, 0.04872, 
                                         0.04879, 0.04888, 0.04896, 0.04905, 0.04916, 0.04935, 0.04951, 
                                         0.04959, 0.04979, 0.05005, 0.05023, 0.05036, 0.05042, 0.05049, 
                                         0.05068, 0.05081, 0.05092, 0.05111, 0.05124, 0.05138, 0.05148, 
                                         0.05163, 0.05168, 0.05192, 0.05205, 0.05222, 0.0523, 0.05246, 
                                         0.05267, 0.05276, 0.05289, 0.05306, 0.05313, 0.05321, 0.05334, 
                                         0.05345, 0.05364, 0.05369, 0.0539, 0.05401, 0.05407, 0.05414, 
                                         0.05427, 0.05438, 0.05457, 0.05476, 0.05494, 0.05504, 0.05509, 
                                         0.05519, 0.05528, 0.05536, 0.05553, 0.05562, 0.05571, 0.05584, 
                                         0.05597, 0.05607, 0.05626, 0.05646, 0.05656, 0.05663, 0.05669, 
                                         0.05678, 0.0568, 0.05686, 0.05703, 0.05719, 0.05725, 0.05733, 
                                         0.0574, 0.0575, 0.05759, 0.05778, 0.05789, 0.05797, 0.05808, 
                                         0.05821, 0.05836, 0.05847, 0.0587, 0.05883, 0.05891, 0.05898, 
                                         0.05913, 0.05923, 0.0594, 0.05951, 0.05957, 0.0597, 0.05976, 
                                         0.05987, 0.06002, 0.0601, 0.06025, 0.0604, 0.06047, 0.06053, 
                                         0.06055, 0.06064, 0.0608, 0.06093, 0.06104, 0.06115, 0.06127, 
                                         0.06132, 0.06144, 0.06165, 0.06182, 0.06189, 0.06199, 0.06204, 
                                         0.06214, 0.06223, 0.06242, 0.06261, 0.06277, 0.06294, 0.06299, 
                                         0.06301, 0.06309, 0.06334, 0.06347, 0.0636, 0.06364, 0.06375, 
                                         0.06385, 0.06394, 0.06406, 0.06423, 0.06438, 0.06453, 0.06467, 
                                         0.06476, 0.0648, 0.06497, 0.06516, 0.0652, 0.06526, 0.06534, 
                                         0.06545, 0.06553, 0.06574, 0.06587, 0.066, 0.06621, 0.06633, 
                                         0.0664, 0.06648, 0.0666, 0.06675, 0.06681, 0.06698, 0.06709, 
                                         0.06723, 0.06725, 0.0673, 0.0674, 0.0675, 0.06761, 0.06771, 0.06784, 
                                         0.06796, 0.06809, 0.06836, 0.06851, 0.06861, 0.06874, 0.06874, 
                                         0.06893, 0.06903, 0.06926, 0.0694, 0.06949, 0.06959, 0.06972, 
                                         0.06982, 0.0699, 0.07001, 0.07017, 0.07022, 0.07028, 0.07038, 
                                         0.07042, 0.07057, 0.07067, 0.07076, 0.07082, 0.07094, 0.07099, 
                                         0.07113, 0.07126, 0.0714, 0.07152, 0.07155, 0.07165, 0.07169, 
                                         0.07179, 0.07186, 0.07204, 0.07221, 0.07235, 0.0725, 0.0726, 
                                         0.07264, 0.07281, 0.07293, 0.07306, 0.07314, 0.07322, 0.07327, 
                                         0.07343, 0.07366, 0.07382, 0.07387, 0.07393, 0.07405, 0.07411, 
                                         0.07419, 0.07442, 0.07473, 0.07485, 0.07494, 0.07502, 0.07504, 
                                         0.0751, 0.07521, 0.07541, 0.07562, 0.07566, 0.07578, 0.07591, 
                                         0.07601, 0.07615, 0.07632, 0.07644, 0.07648, 0.07664, 0.07677, 
                                         0.07685, 0.07693, 0.07703, 0.07712, 0.07722, 0.07728, 0.0774, 
                                         0.07761, 0.07771, 0.0779, 0.07808, 0.07816, 0.07833, 0.07843, 
                                         0.07855, 0.07876, 0.0789, 0.07898, 0.0791, 0.07914, 0.07922, 
                                         0.07931, 0.07941, 0.07953, 0.07967, 0.07978, 0.07984, 0.07986, 
                                         0.07998, 0.0801, 0.08029, 0.08045, 0.08057, 0.08061, 0.08075, 
                                         0.08086, 0.08104, 0.0812, 0.08128, 0.08139, 0.08151, 0.08167, 
                                         0.08185, 0.082, 0.08226, 0.08242, 0.08254, 0.08261, 0.08273, 
                                         0.08279, 0.08297, 0.08319, 0.08328, 0.08338, 0.08358, 0.08368, 
                                         0.08376, 0.08394, 0.08406, 0.08423, 0.08431, 0.08437, 0.08445, 
                                         0.08461, 0.08481, 0.085, 0.08516, 0.08526, 0.0854, 0.0855, 0.08562, 
                                         0.08574, 0.08594, 0.08615, 0.08633, 0.08639, 0.08649, 0.08659, 
                                         0.08675, 0.08707, 0.08725, 0.08727, 0.08748, 0.0877, 0.08778, 
                                         0.08798, 0.08838, 0.08866, 0.08882, 0.08894, 0.08902, 0.08922, 
                                         0.08954, 0.0897, 0.08983, 0.08997, 0.09011, 0.09019, 0.09033, 
                                         0.09047, 0.09073, 0.09097, 0.09111, 0.09127, 0.09145, 0.09158, 
                                         0.09172, 0.09196, 0.09216, 0.09228, 0.09238, 0.09246, 0.09258, 
                                         0.09288, 0.09302, 0.09318, 0.09326, 0.09338, 0.0935, 0.09362, 
                                         0.0938, 0.09407, 0.09421, 0.09433, 0.09445, 0.09451, 0.09465, 
                                         0.09477, 0.09495, 0.09507, 0.09513, 0.09517, 0.09531, 0.09546, 
                                         0.09552, 0.09558, 0.0958, 0.09596, 0.09606, 0.09612, 0.09628, 
                                         0.09645, 0.09665, 0.09677, 0.09681, 0.09691, 0.09711, 0.09725, 
                                         0.09729, 0.09738, 0.0975, 0.09758, 0.0976, 0.09772, 0.0978, 0.098, 
                                         0.09819, 0.09831, 0.09845, 0.09847, 0.09853, 0.09861, 0.0988, 
                                         0.09892, 0.09912, 0.09924, 0.09928, 0.09938, 0.09953, 0.09973, 
                                         0.09981, 0.09987, 0.10006, 0.10014, 0.10024, 0.10038, 0.10052, 
                                         0.10075, 0.10085, 0.10095, 0.10103, 0.10111, 0.1012, 0.1013, 
                                         0.10146, 0.1016, 0.1018, 0.1019, 0.102, 0.102, 0.1022, 0.1023, 
                                         0.1023, 0.1024, 0.1026, 0.1027, 0.1028, 0.1029, 0.103, 0.1031, 
                                         0.1032, 0.1033, 0.1035, 0.1035, 0.1037, 0.1038, 0.1039, 0.104, 
                                         0.104, 0.1041, 0.1041, 0.1042, 0.1044, 0.1046, 0.1047, 0.1049, 
                                         0.1049, 0.105, 0.1052, 0.1053, 0.1054, 0.1055, 0.1056, 0.1057, 
                                         0.1058, 0.1059, 0.106, 0.1061, 0.1062, 0.1063, 0.1064, 0.1064, 
                                         0.1066, 0.1067, 0.1069, 0.1069, 0.1071, 0.1072, 0.1072, 0.1072, 
                                         0.1074, 0.1075, 0.1076, 0.1076, 0.1076, 0.1078, 0.1079, 0.1081, 
                                         0.1082, 0.1083, 0.1083, 0.1084, 0.1085, 0.1087, 0.1089, 0.1091, 
                                         0.1092, 0.1093, 0.1095, 0.1096, 0.1098, 0.11, 0.1101, 0.1102, 
                                         0.1103, 0.1103, 0.1104, 0.1105, 0.1106, 0.1107, 0.1108, 0.1109, 
                                         0.1109, 0.1111, 0.1113, 0.1114, 0.1115, 0.1115, 0.1117, 0.1118, 
                                         0.1118, 0.112, 0.1121, 0.1122, 0.1123, 0.1123, 0.1124, 0.1126, 
                                         0.1127, 0.1129, 0.1131, 0.1132, 0.1133, 0.1133, 0.1134, 0.1136, 
                                         0.1137, 0.1138, 0.114, 0.114, 0.1142, 0.1142, 0.1145, 0.1147, 
                                         0.1149, 0.1149, 0.1149, 0.115, 0.1151, 0.1152, 0.1153, 0.1154, 
                                         0.1156, 0.1157, 0.1157, 0.1158, 0.1159, 0.116, 0.1162, 0.1163, 
                                         0.1164, 0.1165, 0.1166, 0.1167, 0.1169, 0.117, 0.1171, 0.1172, 
                                         0.1173, 0.1174, 0.1176, 0.1177, 0.1178, 0.1179, 0.1181, 0.1182, 
                                         0.1182, 0.1183, 0.1185, 0.1185, 0.1187, 0.1187, 0.1188, 0.1189, 
                                         0.119, 0.119, 0.119, 0.1191, 0.1191, 0.1192, 0.1193, 0.1194, 
                                         0.1196, 0.1197, 0.1198, 0.1199, 0.1199, 0.12, 0.1201, 0.1202, 
                                         0.1204, 0.1205, 0.1206, 0.1207, 0.1208, 0.121, 0.1211, 0.1213, 
                                         0.1214, 0.1215, 0.1216, 0.1217, 0.1218, 0.1219, 0.1219, 0.1221, 
                                         0.1222, 0.1223, 0.1223, 0.1225, 0.1227, 0.1228, 0.1229, 0.1229, 
                                         0.123, 0.1231, 0.1232, 0.1234, 0.1235, 0.1236, 0.1236, 0.1238, 
                                         0.1239, 0.1241, 0.1243, 0.1243, 0.1244, 0.1245, 0.1246, 0.1247, 
                                         0.1248, 0.1249, 0.1251, 0.1252, 0.1253, 0.1253, 0.1254, 0.1255, 
                                         0.1257, 0.1258, 0.1259, 0.1261, 0.1262, 0.1263, 0.1264, 0.1266, 
                                         0.1267, 0.1268, 0.127, 0.1271, 0.1272, 0.1272, 0.1274, 0.1276, 
                                         0.1277, 0.1278, 0.1278, 0.1279, 0.128, 0.1283, 0.1284, 0.1285, 
                                         0.1286, 0.1288, 0.1288, 0.129, 0.1291, 0.1295, 0.1297, 0.1297, 
                                         0.1298, 0.1298, 0.1299, 0.1301, 0.1303, 0.1304, 0.1304, 0.1304, 
                                         0.1306, 0.1308, 0.131, 0.1311, 0.1312, 0.1313)), 
              class = "data.frame", row.names = c(NA, -1095L))
}

# Fitted model spline info
si <- function(){
    # Output from dput()
    structure(list(variables = c("age", "egfr", "ucrea", "hba1c", 
                                 "total_cholesterol", "ldl", "tc_hdl", "index_alb"), 
                   lower_boundary = c(18.031485284052, 60.0001689340495, 0, 6.011, 
                                      34.0837857142857, 5.9005, 1, 0), 
                   upper_boundary = c(99.6112251882272, 185.770461135844, 994.3248, 
                                      175, 458.2395, 314.0004, 36.7714285714286, 29.97)), 
              row.names = c(NA, -8L), class = c("tbl_df", "tbl", "data.frame"))
}

# Fit individual survival curve
survcurv <- function(# Model variables
                     female,
                     age,
                     education,
                     diabetes_months,
                     egfr,
                     urine_creatinine,
                     hba1c,
                     total_cholesterol,
                     hdl,
                     ldl,
                     triglycerides,
                     total_cholesterol_hdl_ratio,
                     albumin,
                     atrial_fibrillation,
                     congestive_heart_failure,
                     cerebrovascular_disease,
                     hypertension,
                     ischemic_heart_disease,
                     neuropathy,
                     peripheral_vascular_disease,
                     retinopathy,
                     ulcer,
                     aspirin,
                     beta_blockers,
                     glucose_lowering_drugs,
                     calcium_channel_blockers,
                     anticoagulants,
                     antihypertensives,
                     insulin,
                     mineralocorticoid_receptor_antagonists,
                     ras_inhibitors,
                     statins,
                     # Data frame
                     newdata = NULL,
                     # Model info
                     mod = mv(),
                     bhs = bh(),
                     msi = si(),
                     # Output
                     interactive = FALSE
                     ){
    # Define column names
    columns <- c("female", "age", "education", "diab_months", "egfr", "ucrea", "hba1c", "total_cholesterol", "hdl", "ldl", "triglycerides",
                 "tc_hdl", "index_alb", "fibrillation", "chf", "cvd", "hypertension", "ihd", "neuropathy", "pvd", "retinopathy", "ulcer",                         
                 "aspirin", "bblockers", "glucose_lowering", "ccbs", "anticoagulants", "antihypertensives", "insulin", "mras", "rasi", "statins")
    
    # Get model vars
    model_vars <<- mod
    
    # Get spline info
    model_spline_info <<- msi
    
    # Create data frame based on variables if newdata is not entered
    if(is.null(newdata)){
        # Create data frame based on variables
        data_new <- c(female, age, education, diabetes_months, egfr, urine_creatinine, hba1c, total_cholesterol, hdl, ldl, triglycerides,
                      total_cholesterol_hdl_ratio, albumin, atrial_fibrillation, congestive_heart_failure, cerebrovascular_disease,
                      hypertension, ischemic_heart_disease, neuropathy, peripheral_vascular_disease, retinopathy, ulcer, aspirin, beta_blockers,
                      glucose_lowering_drugs, calcium_channel_blockers, anticoagulants, antihypertensives, insulin, mineralocorticoid_receptor_antagonists,
                      ras_inhibitors, statins) %>%
            # Change to matrix
            as.matrix() %>%
            # Transpose
            t() %>% as.data.frame() %>% print()
            # To data frame
            # as.data.frame() %>%
            # # Set column names
            # set_colnames(columns)
    }
    
    # Else new data is the entered data frame
    if(!is.null(newdata)){
        # New data is data frame
        data_new <- newdata %>%
            # Select columns of interest; this returns error if not all necessary columns are present
            select(all_of(columns))
    }
    
    # Add studynr and imputation
    data_new <- mutate(data_new, studynr = 1, .imp = 1)
    
    # Calculate linear predictor
    lp <- pred(data_new, "fine-gray", NULL, lpsamp = 1.406445, mod = mod)[["lp"]]
    
    # Time points of interest
    times <- bhs[["time"]]
    
    # Data with probability of outcome
    dat <- data.frame(time = times,
                      Month = times / (365.25 / 12),
                      Probability = 1 - exp(-exp(lp) * (bhs[["bh"]])))
    
    # Plot data
    plot <- ggplot(dat, aes(x = Month, y = Probability)) +
        # Geometries
        geom_line(colour = "darkred", linewidth = 1) +
        # Scaling
        scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.1), labels = paste0(seq(0, 100, 10), "%")) +
        scale_x_continuous(limits = c(0, max(dat[["Month"]]) + 1), breaks = seq(0, max(dat[["Month"]] + 1), 3)) +
        # Labels
        xlab("Time (months)") +
        ylab("Probability of outcome") +
        # Aesthetics
        theme(panel.background = element_blank(),
              panel.grid.major = element_line(colour = "lightgrey"),
              axis.line = element_line(colour = "black"))
    
    # If interactive, return interactive plot, else just plot
    if(interactive) return(ggplotly(plot)) else return(plot)
}

# Transition matrix
tm <- function(){
    structure(c(NA, NA, NA, NA, 1L, NA, NA, NA, 2L, 4L, NA, NA, 3L, 5L, 6L, NA), dim = c(4L, 4L), 
              dimnames = list(from = c("Baseline", "Microalbuminuria", "Macroalbuminuria", "Death"), 
                              to = c("Baseline", "Microalbuminuria", "Macroalbuminuria", "Death")))
}

# Data preparation for multi state prediction
plot_mstate_prep <- function(# Model variables
                             female,
                             age,
                             education,
                             diabetes_months,
                             egfr,
                             urine_creatinine,
                             hba1c,
                             total_cholesterol,
                             hdl,
                             ldl,
                             triglycerides,
                             total_cholesterol_hdl_ratio,
                             albumin,
                             atrial_fibrillation,
                             congestive_heart_failure,
                             cerebrovascular_disease,
                             hypertension,
                             ischemic_heart_disease,
                             neuropathy,
                             peripheral_vascular_disease,
                             retinopathy,
                             ulcer,
                             aspirin,
                             beta_blockers,
                             glucose_lowering_drugs,
                             calcium_channel_blockers,
                             anticoagulants,
                             antihypertensives,
                             insulin,
                             mineralocorticoid_receptor_antagonists,
                             ras_inhibitors,
                             statins,
                             # Data frame
                             newdata = NULL,
                             # Model info
                             transition_matrix = tm(),
                             # Output
                             area = TRUE,
                             interactive = FALSE
){
    # Define column names
    columns <- c("female", "age", "education", "diab_months", "egfr", "ucrea", "hba1c", "total_cholesterol", "hdl", "ldl", "triglycerides",
                 "tc_hdl", "index_alb", "fibrillation", "chf", "cvd", "hypertension", "ihd", "neuropathy", "pvd", "retinopathy", "ulcer",                         
                 "aspirin", "bblockers", "glucose_lowering", "ccbs", "anticoagulants", "antihypertensives", "insulin", "mras", "rasi", "statins")
    
    # Create data frame based on variables if newdata is not entered
    if(is.null(newdata)){
        # Create data frame based on variables
        data_new <- c(female, age, education, diabetes_months, egfr, urine_creatinine, hba1c, total_cholesterol, hdl, ldl, triglycerides,
                      total_cholesterol_hdl_ratio, albumin, atrial_fibrillation, congestive_heart_failure, cerebrovascular_disease,
                      hypertension, ischemic_heart_disease, neuropathy, peripheral_vascular_disease, retinopathy, ulcer, aspirin, beta_blockers,
                      glucose_lowering_drugs, calcium_channel_blockers, anticoagulants, antihypertensives, insulin, mineralocorticoid_receptor_antagonists,
                      ras_inhibitors, statins) %>%
            # Change to matrix
            as.matrix() %>%
            # Transpose
            t() %>%
            # To data frame
            as.data.frame() %>%
            # Set column names
            set_colnames(columns)
    }
    
    # Else new data is the entered data frame
    if(!is.null(newdata)){
        # New data is data frame
        data_new <- newdata %>%
            # Select columns of interest; this returns error if not all necessary columns are present
            select(all_of(columns))
    }
    
    # Create new column names
    new_cols <- sort(do.call("c", lapply(1:6, \(x) paste0(columns, ".", x))))
    
    # Duplicate new data for each transition (6)
    dat_new <- rbind(data_new, data_new, data_new, data_new, data_new, data_new)
    
    # Add new (empty) columns
    for(i in new_cols) dat_new[[i]] <- NA
    
    # Set new column values to covariate value if stratum, else 0
    dat_new <- dat_new %>%
        # Set column values
        mutate(across(all_of(new_cols), \(x) x = ifelse(as.numeric(str_extract(cur_column(), "(?<=\\.)\\d")) == row_number(), 
                                                        dat_new[[str_replace(cur_column(), "\\.\\d", "")]], 0))) %>%
        # Add stratum and trans variables
        mutate(strata = 1:6,
               trans = 1:6) %>%
        # Remove unnecessary variables
        select(-all_of(columns))
 
    # Load mstate fit
    load(url("https://github.com/rjjanse/alb/raw/main/multistate/mstate_fit.Rdata"))
    
    # Load anonymous development data
    load(url("https://github.com/rjjanse/alb/raw/main/multistate/dat_dev.Rdata"))
    
    # Empty data frame
    dat_mstate_dev <- NULL
    
    # Prepare individual data
    dat_fit <- msfit(mstate_object, newdata = dat_new, trans = transition_matrix)
    
    # Get individual probabilities (predt = 365 * 3, else code does not function)
    probs <- probtrans(dat_fit, predt = 0, method = "aalen", direction = "forward", variance = FALSE)[[1]]
    
    # Create plotting data
    dat_plot <- pivot_longer(probs, cols = pstate1:pstate4) %>%
        # Change variables
        mutate(# Change state names
               state = case_match(name,
                                  "pstate1" ~ "Event-free",
                                  "pstate2" ~ "Microalbuminuria",
                                  "pstate3" ~ "Macroalbuminuria",
                                  "pstate4" ~ "Death"),
               # State as factor
               state = factor(state, levels = c("Death", "Macroalbuminuria", "Microalbuminuria", "Event-free")),
               # Time to months
               month = time / (365.25 / 12))
    
    # Create base plot
    plot <- ggplot(dat_plot, aes(x = month, y = value, fill = state, colour = state)) 
    
    # If area, add area, else add line
    if(area) plot <- plot + geom_area() else plot <- plot + geom_line()
    
    # Finish plot
    plot <- plot +
        # Scalings
        scale_colour_manual(values = c("#BC4749", "#1D3354", "#467599", "#F2E8CF")) +
        scale_fill_manual(values = c("#BC4749", "#1D3354", "#467599", "#F2E8CF")) +
        scale_x_continuous(breaks = seq(0, 36, 3), name = "Time (months)", expand = c(0, 0)) +
        scale_y_continuous(breaks = seq(0, 1, 0.1), labels = paste0(seq(0, 100, 10), "%"),
                           name = "Probability", expand = c(0, 0)) +
        # Transformations
        coord_cartesian(xlim = c(0, 36), ylim = c(0, 1)) +
        # Aesthetics
        theme(panel.border = element_rect(colour = "black", fill = "transparent"),
              legend.position = "bottom",
              legend.title = element_blank(),
              panel.background = element_rect(fill = "white"),
              panel.grid.major = element_line(colour = "darkgrey"),
              panel.grid.minor = element_blank())
    
    # If interactive, return interactive plot, else just plot
    if(interactive) return(ggplotly(plot)) else return(plot)
}
    
