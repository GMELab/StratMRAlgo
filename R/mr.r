library(data.table)
library(dplyr)
library(metafor)

#' @title Stratified Mendelian Randomization
#' @param method Either "by_x", "by_g", "resid_strat", "raw_strat", or "resid_wald"
#' @param x The exposure name
#' @param y The outcome name
#' @param c The stratifying variable name
#' @param quantile_number The number of quantiles or strata to use
#' @param data The name of the file containing the data or a data.table object. Column names need to be in a specific format; see README page for details.
#' @param input_dir The directory where the data is stored
#' @param output_dir The directory where the output should be stored, defaults to NULL meaning the results will not be saved
#' @export
mr <- function(method, exposure_name, outcome_name, stratifying_name, quantile_number, data, input_dir = getwd(), output_dir = NULL) {
  if (method != "by_x" && method != "by_g" && method != "resid_strat" && method != "raw_strat" && method != "resid_wald") {
    stop("method must be either 'by_x', 'by_g', 'resid_strat', 'raw_strat', or 'resid_wald'")
  }

  start <- Sys.time()

  quantile_number <- as.numeric(quantile_number)

  if (grepl("_EXP$", exposure_name)) {
    exposure_name_full <- exposure_name
    exposure_name <- gsub("_EXP$", "", exposure_name)
  } else {
    exposure_name_full <- paste0(exposure_name, "_EXP")
    exposure_name <- x
  }

  if (grepl("_STRAT$", stratifying_name)) {
    stratifying_name_full <- stratifying_name
    stratifying_name <- gsub("_STRAT$", "", stratifying_name)
  } else {
    stratifying_name_full <- paste0(stratifying_name, "_STRAT")
  }

  if (grepl("_CASE_STATUS$", outcome_name)) {
    outcome_name_full <- outcome_name
    outcome_name <- gsub("_CASE_STATUS$", "", outcome_name)
    outcome_type <- "binom"
  } else if (grepl("_OUTCOME$", outcome_name)) {
    outcome_name_full <- outcome_name
    outcome_name <- gsub("_OUTCOME$", "", outcome_name)
    outcome_type <- "cont"
  } else {
    stop("could not determine the outcome type.")
  }
  if (method == "resid_wald" && outcome_type == "binom") {
    stop("resid_wald method is only available for continuous outcomes.")
  }
  prs_name_full <- paste0(exposure_name, "_PRS")

  message("Starting stratified MR analysis at ", start)
  message("Exposure: ", exposure_name, ", Stratifying variable: ", stratifying_name, ", Outcome: ", outcome_name, ", Quantile number: ", quantile_number, ", Method: ", method)

  if (is.character(data)) {
    if (!file.exists(file.path(input_dir, data))) {
      stop("data file does not exist.")
    }
    df <- fread(
      file.path(input_dir, data),
      stringsAsFactors = FALSE,
      header = TRUE,
      data.table = TRUE,
    )
  } else {
    if (!is.data.table(data)) {
      stop("data must be a data.table.")
    }
    df <- data
  }

  ################################################################################################################################################
  # Step 1) Create dataframe with variables of interest
  ################################################################################################################################################

  message("STEP 1: Extracting necessary variables...")
  
  covars <- grep("_COVAR$", colnames(df), value = TRUE)
  mr_df <- df[, c(
    exposure_name_full, outcome_name_full, stratifying_name_full, prs_name_full,
    covars
  ), with = FALSE]

  pts_full <- dim(mr_df)[1]
  message("Total number of participants at start of analysis: ", dim(mr_df)[1])
  
  # Metadata ------------------------------------------------------

  # Count participants with NA in the PRS prior to removal
  prs_na_list_full <- mr_df[is.na(mr_df[[prs_name_full]]), ]
  prs_na_count_full <- dim(prs_na_list_full)[1]

  # Count participants with missing exposure variable prior to removal
  exp_na_list_full <- mr_df[is.na(mr_df[[exposure_name_full]]), ]
  exp_na_count_full <- dim(exp_na_list_full)[1]

  # Count participants with missing outcomes prior to removal
  outcome_na_list_full <- mr_df[is.na(mr_df[[outcome_name_full]]), ]
  outcome_na_count_full <- dim(outcome_na_list_full)[1]

  # Count participants with NA in the stratifying variable prior to removal
  strat_na_list_full <- mr_df[is.na(mr_df[[stratifying_name_full]]), ]
  strat_na_count_full <- dim(strat_na_list_full)[1]

  # Count number of cases prior to removal for binary outcomes
  if (outcome_type == "binom") {
    case_list_full <- mr_df[(mr_df[[outcome_name_full]]) == 1, ]
    case_count_full <- dim(case_list_full)[1]
  }

  ################################################################################################################################################
  # Step 2) Remove participants with NA in PRS, STRAT, COVAR, EXPOSURE, and OUTCOME
  ################################################################################################################################################

  message("STEP 2: Keeping only complete cases...")
  
  # Keep only complete cases
  combined <- mr_df[complete.cases(mr_df), ]

  # Count number of participants in data frame after removal
  pts_with_complete_vars <- dim(combined)[1]

  message("Total number of participants with the PRS, stratifying variable, covariates, exposure, and outcome available: ", dim(combined)[1])

  # Metadata ------------------------------------------------------
  
  # Count remaining individuals with NA in the PRS (should be 0)
  prs_na_list <- combined[is.na(combined[[prs_name_full]]), ]
  prs_na_count <- dim(prs_na_list)[1]

  # Count remaining participants with missing exposure variable (should be 0)
  exp_na_list <- combined[is.na(combined[[exposure_name_full]]), ]
  exp_na_count <- dim(exp_na_list)[1]

  # Count remaining participants with missing outcomes (should be 0)
  outcome_na_list <- combined[is.na(combined[[outcome_name_full]]), ]
  outcome_na_count <- dim(outcome_na_list)[1]

  # Count remaining participants with missing stratifying variable (should be 0)
  strat_na_list <- combined[is.na(combined[[stratifying_name_full]]), ]
  strat_na_count <- dim(strat_na_list)[1]

  # Count remaining number of cases
  if (outcome_type == "binom") {
    case_list <- combined[(combined[[outcome_name_full]]) == 1, ]
    case_count <- dim(case_list)[1]
  }

  ################################################################################################################################################
  # Step 3) Define strata
  ################################################################################################################################################

  # Code for doubly ranking steps adapted from the DRMR R package written by Tian et al: https://github.com/HDTian/DRMR/blob/main/R/Stratify.R
  # Variable naming: X = exposure, Y = outcome, X0 (or C) = stratifying variable, G = instrument
  
  message("STEP 3: Defining strata...")
  
  message("Performing first ranking step...")
  if (method == "by_x") {
    combined_1 <- combined[order(combined[[exposure_name_full]]), ] # rank by X
  } else if (method == "by_g") {
    combined_1 <- combined[order(combined[[prs_name_full]]), ] # rank by G
  } else if (method == "resid_strat") {
    temp_exp <- paste0(exposure_name_full, "_2")
    combined[[temp_exp]] <- combined[[exposure_name_full]]^2 # calculate X^2
    if (outcome_type == "binom") {
      message("Residualizing stratifying variable by X, X2, and Y (binary variable)...")
      combined$x0_resid <- resid(lm(paste0(stratifying_name_full, " ~ ", exposure_name_full, " + ", exposure_name_full, "_2 + ", outcome_name_full), data = combined)) # calculate residualized C
    } else if (outcome_type == "cont") {
      message("Residualizing stratifying variable by X, X2, Y, and Y2 (continuous variable)...")
      temp_out <- paste0(outcome_name_full, "_2")
      combined[[temp_out]] <- combined[[outcome_name_full]]^2 # calculate Y^2
      combined$x0_resid <- resid(lm(paste0(stratifying_name_full, " ~ ", exposure_name_full, " + ", exposure_name_full, "_2 + ", outcome_name_full, " + ", outcome_name_full, "_2"), data = combined)) # calculate residualized C
    }
    combined_1 <- combined[order(combined$x0_resid), ] # rank by residualized C
  } else if (method == "raw_strat" || method == "resid_wald") {
    combined_1 <- combined[order(combined[[stratifying_name_full]]), ] # rank by C
  }

  message("Performing first stratification step...")
  n <- nrow(combined_1)
  combined_2 <- combined_1
  combined_2$ID <- 1:n # create column with the rank
  if (method == "by_x" || method == "by_g") {
    combined_2$pre_stratum <- rep(1:(floor(n / quantile_number) + 1), each = quantile_number, length.out = n) # create pre-strata based on rank of either X or G
    combined_3 <- combined_2
    
    message("Performing second ranking step...") # within each pre-stratum, rank in ascending order based on stratifying variable
    combined_3 <- arrange(combined_3, combined_3[[stratifying_name_full]]) 
    combined_3 <- arrange(combined_3, pre_stratum) # nolint
    combined_3 <- as.data.table(lapply(combined_3, as.numeric)) # ensure that all columns are numeric

    message("Performing second stratification step...") # within each pre-stratum, rank in ascending order based on stratifying variable
    # Assign a vector to final_stratum that repeats numbers from 1 to q; the last may have fewer elements if not evenly divisible
    combined_3$final_stratum <- as.vector(unlist(
      sapply(
        as.numeric(table(combined_3$pre_stratum)),
        function(x) sort(rep(1:quantile_number, length.out = x))
      )
    ))
    
    combined_final <- combined_3

  } else if (method == "resid_strat" || method == "raw_strat" || method == "resid_wald") {
    combined_2$final_stratum <- cut(combined_2$ID, quantile_number, include.lowest = TRUE, labels = FALSE) # stratification based on the first ranking step by C or residualized C

    combined_final <- combined_2
  }

  if (method == "resid_wald") {
    setDT(combined_final)[, paste0(stratifying_name_full, "_2") := get(stratifying_name_full)^2] # need C^2 for later steps
  }

  ################################################################################################################################################
  # Step 4) Define key variables as vectors for stratified MR
  ################################################################################################################################################

  message("STEP 4: Define key variables as vectors...")
  
  # x: exposure
  x <- combined_final[[exposure_name_full]]

  # y: outcome
  y <- combined_final[[outcome_name_full]] # nolint

  # g: instrumental variable
  g <- combined_final[[prs_name_full]] # nolint

  covars <- grep("_COVAR$", colnames(combined_final), value = TRUE)
  
  # Remove age covariates if age is the stratifying variable
  if (grepl("^AGE_STRAT$", stratifying_name_full, ignore.case = TRUE)) {
    age_covars <- grep("^AGE", colnames(combined_final), value = TRUE, ignore.case = TRUE)
    covars <- covars[!covars %in% age_covars]
    if (length(age_covars) > 0) {
      message("Age-related variables ending with _COVAR in this list will be removed: ", paste(age_covars, collapse = ", "))
    }
  }
  
  # Remove age covariates if age is the continuous outcome (negative control analyses)
  if (grepl("^AGE_OUTCOME$", outcome_name_full, ignore.case = TRUE) && outcome_type == "cont") {
    age_covars <- grep("^AGE", colnames(combined_final), value = TRUE, ignore.case = TRUE)
    covars <- covars[!covars %in% age_covars]
    if (length(age_covars) > 0) {
      message("Age-related variables ending with _COVAR in this list will be removed: ", paste(age_covars, collapse = ", "))
    }
  }
  
  # Remove sex covariates if sex is the binary outcome (negative control analyses)
  if (grepl("^SEX_CASE_STATUS$", outcome_name_full, ignore.case = TRUE) && outcome_type == "binom") {
    sex_covars <- grep("^SEX", colnames(combined_final), value = TRUE, ignore.case = TRUE)
    covars <- covars[!covars %in% sex_covars]
    if (length(sex_covars) > 0) {
      message("Sex-related variables ending with _COVAR in this list will be removed: ", paste(sex_covars, collapse = ", "))
    }
  }
  
  # covar: covariates
  covar <- combined_final[, covars, with = FALSE] # nolint
  covar <- as.matrix(covar)

  # x0: stratifying variable (residualized for resid_strat)
  if (method == "by_x" || method == "by_g" || method == "raw_strat" || method == "resid_wald") {
    x0 <- combined_final[[stratifying_name_full]]
  } else if (method == "resid_strat") {
    x0 <- combined_final$x0_resid
  }

  # x0q: final strata ID (which stratum the participant belongs to; doubly-ranked strata for DR methods, residualized strata for resid_strat, regular strata otherwise)
  x0q <- combined_final$final_stratum

  # number of strata or quantiles
  q <- quantile_number

  ################################################################################################################################################
  # Step 5) Calculate stratum-specific estimates
  ################################################################################################################################################

  # Code for stratum-specific estimates adapted from the nlmr R package written by Burgess et al: https://github.com/jrs95/nlmr/blob/master/R/lace.R 
  
  message("STEP 5: Calculating stratum-specific estimates...")

  ycoef <- NULL
  ycoef_se <- NULL
  xcoef <- NULL
  xcoef_se <- NULL
  x0mean <- NULL
  coef <- NULL
  coef_se <- NULL
  coef_95_ci <- NULL
  
  if (method == "resid_wald") {
    y <- resid(lm(paste0(outcome_name_full, " ~ ", stratifying_name_full, " + ", stratifying_name_full, "_2"), data = combined_final))
    x <- resid(lm(paste0(exposure_name_full, " ~ ", stratifying_name_full, " + ", stratifying_name_full, "_2"), data = combined_final))
  }

  for (i in 1:q) {
    if (outcome_type == "binom") {
      ycoef_model <- glm(y[x0q == i] ~ g[x0q == i] + covar[x0q == i, , drop = FALSE], family = "binomial")
    } else if (outcome_type == "cont") {
      ycoef_model <- lm(y[x0q == i] ~ g[x0q == i] + covar[x0q == i, , drop = FALSE])
    }

    ycoef[i] <- ycoef_model$coef[2] # ycoef in stratum
    ycoef_se[i] <- summary(ycoef_model)$coef[2, 2] # SE of ycoef

    if (outcome_type == "binom") {
      xcoef_model <- lm(x[(x0q == i & y == 0)] ~ g[(x0q == i & y == 0)] + covar[(x0q == i & y == 0), , drop = FALSE]) # omits cases
    } else if (outcome_type == "cont") {
      xcoef_model <- lm(x[(x0q == i)] ~ g[(x0q == i)] + covar[(x0q == i), , drop = FALSE]) # no cases to omit
    }

    xcoef[i] <- xcoef_model$coef[2] # xcoef in stratum
    xcoef_se[i] <- summary(xcoef_model)$coef[2, 2] # SE of xcoef
    
    x0mean[i] <- mean(x0[x0q == i]) # mean x0 in stratum
    coef[i] <- ycoef[i] / xcoef[i] # stratum-specific Wald Ratio
    coef_se[i] <- ycoef_se[i] / xcoef[i] # error calculated with first order estimation using the asymptotic SE of ratio estimate
    coef_95_ci[i] <- abs(coef_se[i]) * 1.96
  }

  loc <- data.frame(ycoef = ycoef, ycoef_se = ycoef_se, xcoef = xcoef, xcoef_se = xcoef_se, x0mean = x0mean, coef = coef, coef_se = coef_se, coef_95_ci = coef_95_ci) # nolint

  ################################################################################################################################################
  # Step 6) Meta-regression of stratum estimates
  ################################################################################################################################################

  message("STEP 6: Meta-regression of stratum estimates...")

  # Regress stratum-specific estimates on x0mean per stratum using a fixed effects meta-analysis model
  rma_model <- rma.uni(coef ~ x0mean, vi = (coef_se)^2, method = "FE")

  # Output summary of model
  beta_rma <- coef(summary(rma_model))[, "estimate"][2] # binary outcome units: log-odds Y per unit X per unit C; continuous outcome units: units Y per unit X per unit C
  se_rma <- coef(summary(rma_model))[, "se"][2]
  p_rma <- coef(summary(rma_model))[, "pval"][2]

  ################################################################################################################################################
  # Step 7) Summarize all relevant results into table
  ################################################################################################################################################

  message("STEP 7: Summarizing all relevant results into table...")

  current_date <- gsub("-", "_", Sys.Date())
  if (outcome_type == "binom") {
    results <- data.frame(
      exposure_name, stratifying_name, outcome_name, quantile_number, current_date, # inputs
      pts_full, prs_na_count_full, exp_na_count_full, outcome_na_count_full, strat_na_count_full, case_count_full, # metadata on full cohort
      pts_with_complete_vars, prs_na_count, exp_na_count, outcome_na_count, strat_na_count, case_count, # metadata on participants included in stratified calculations
      beta_rma, se_rma, p_rma # meta-regression outputs
    )
  } else if (outcome_type == "cont") {
    results <- data.frame(
      exposure_name, stratifying_name, outcome_name, quantile_number, current_date, # inputs
      pts_full, prs_na_count_full, exp_na_count_full, outcome_na_count_full, strat_na_count_full, # metadata on full cohort
      pts_with_complete_vars, prs_na_count, exp_na_count, outcome_na_count, strat_na_count, # metadata on participants included in stratified calculations
      beta_rma, se_rma, p_rma # meta-regression outputs
    )
  }

  if (!is.null(output_dir)) {
    if (!dir.exists(output_dir)) {
      dir.create(output_dir)
    }
    message("Saving results...")

    fwrite(
      results,
      file = file.path(output_dir, paste0("exp_", exposure_name, "_strat_", stratifying_name, "_out_", outcome_name, "_", quantile_number, "q_", current_date, "_", method, ".csv")), quote = FALSE,
      na = "NA"
    )
  }

  message("STRATIFIED MR ANALYSIS COMPLETE")
  print(Sys.time() - start)

  return(list(p = p_rma, beta = beta_rma, se = se_rma))
}
