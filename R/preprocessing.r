library(data.table)

#' @title Preprocess data for stratified MR
#' @param input_dir The directory where the input data is stored, defaults to the current working directory.
#' @param output_dir The directory where the output data will be stored, defaults to NULL meaning that the output will not be saved.
#' @param data The name of the file containing the data, or a data.table object, defaults to "data.csv". Column names need to be in a specific format; see README page for details.
#' @param diagnostically_related The name of the file containing the diagnostically related pairs or a data.table object, defaults to NULL meaning that no diagnostically related pairs are provided. No header required; see README page for details.
#' @param binary The file name of the binary combinations, or a data frame with three columns (EXP, STRAT, CASE_STATUS). Column headers required. Defaults to NULL meaning the combinations will be automatically generated.
#' @param continuous The file name of the continuous combinations, or a data frame with three columns (EXP, STRAT, OUTCOME). Column headers required. Defaults to NULL meaning the combinations will be automatically generated.
#' @export
preprocess <- function(
    input_dir = getwd(),
    output_dir = NULL,
    data = "data.csv",
    diagnostically_related = NULL,
    binary = NULL,
    continuous = NULL) {

  start <- Sys.time()
  message("Starting preprocessing at ", start)

  ###############################################################################################
  # Step 1: Load packages and data
  ###############################################################################################

  # Load in data
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

  if (!is.data.table(df)) {
    stop("df is still not a data.table.")
  }

  message("Loaded data.")

  # Load in the diagnostically related pairs, if provided
  if (!is.null(diagnostically_related)) {
    if (is.character(diagnostically_related)) {
      if (!file.exists(file.path(input_dir, diagnostically_related))) {
        stop("diagnostically related file does not exist.")
      }
      diagnostically_related <- fread(
        file.path(input_dir, diagnostically_related),
        stringsAsFactors = FALSE,
        header = FALSE,
        data.table = TRUE,
      )
    } else {
      if (!is.data.table(diagnostically_related)) {
        stop("diagnostically_related must be a data.table.")
      }
    }
    message("Loaded diagnostically related pairs.")
  } else {
    diagnostically_related <- NULL
    message("No diagnostically related pairs provided.")
  }

  ###############################################################################################
  # Step 2: Create list of exposure - stratifying variable - outcome combinations
  ###############################################################################################

  # Separate out the variables provided
  variables <- colnames(df)
  exp <- grep("_EXP$", variables, value = TRUE)
  strat <- grep("_STRAT$", variables, value = TRUE)
  binary_outcome <- grep("_CASE_STATUS$", variables, value = TRUE)
  continuous_outcome <- grep("_OUTCOME$", variables, value = TRUE)
  covar <- grep("_COVAR$", variables, value = TRUE)

  if (!is.null(binary)) {
    if (is.character(binary)) {
      if (!file.exists(file.path(input_dir, binary))) {
        stop("binary file does not exist.")
      }
      # Load in the binary matrix file
      binary <- fread(
        file.path(input_dir, binary),
        stringsAsFactors = FALSE,
        header = TRUE,
      )
    } else {
      if (!is.data.table(binary)) {
        stop("binary must be a data.table.")
      }
    }
    message("Loaded binary matrix.")
  } else {
    # Create the binary matrix
    binary <- expand.grid(exp, strat, binary_outcome)
    if (sum(
      dim(binary)[1] == length(exp) * length(strat) * length(binary_outcome)
    )) {
      message(paste0("All ", dim(binary)[1], " combinations captured in binomial outcomes matrix."))
    } else {
      message(paste0(length(exp) * length(strat) * length(binary_outcome) - dim(binary)[1], " combinations missing in binomial outcomes matrix."))
    }
    binary <- as.data.table(binary)
    colnames(binary) <- c("EXP", "STRAT", "CASE_STATUS")
    message("Binary outcome matrix automatically generated using data provided.")
  }

  if (!is.null(continuous)) {
    if (is.character(continuous)) {
      if (!file.exists(file.path(input_dir, continuous))) {
        stop("continuous file does not exist.")
      }
      # Load in the continuous matrix file
      continuous <- fread(
        file.path(input_dir, continuous),
        stringsAsFactors = FALSE,
        header = TRUE,
      )
    } else {
      if (!is.data.table(continuous)) {
        stop("continuous must be a data.table.")
      }
    }
    message("Loaded continuous matrix.")
  } else {
    # Create the continuous matrix
    continuous <- expand.grid(exp, strat, continuous_outcome)
    if (sum(
      dim(continuous)[1] == length(exp) * length(strat) * length(continuous_outcome)
    )) {
      message(paste0("All ", dim(continuous)[1], " combinations captured in continuous outcomes matrix."))
    } else {
      message(paste0(length(exp) * length(strat) * length(continuous_outcome) - dim(continuous)[1], " combinations missing in continuous outcomes matrix."))
    }
    continuous <- as.data.table(continuous)
    colnames(continuous) <- c("EXP", "STRAT", "OUTCOME")
    message("Continuous outcome matrix automatically generated using data provided.")
  }

  ###############################################################################################
  # Step 3: Output error if there are any binary exposures or stratifying variables
  ###############################################################################################

  exp_strat_names <- c(exp, strat)
  exp_strat_df <- df[, exp_strat_names, with = FALSE]

  # Identify binary columns
  binary_cols <- which(sapply(exp_strat_df, function(x) {
    length(unique(x)) == 2
  }))

  if (any(binary_cols)) {
    binary_col_names <- names(exp_strat_df)[binary_cols]
    stop(paste0("binary exposures or stratifying variables detected: ", paste(binary_col_names, collapse = ", ")))
  }

  ###############################################################################################
  # Step 4: Filter out combinations with identical variables
  ###############################################################################################

  # Binary outcome combinations without the _EXP _STRAT and _CASE_STATUS indicators
  binary_simple <- binary
  binary_simple[[1]] <- sub("_EXP$", "", binary[[1]])
  binary_simple[[2]] <- sub("_STRAT$", "", binary[[2]])
  binary_simple[[3]] <- sub("_CASE_STATUS$", "", binary[[3]])

  # Continuous outcome combinations without the _EXP _STRAT and _OUTCOME indicators
  continuous_simple <- continuous
  continuous_simple[[1]] <- sub("_EXP$", "", continuous[[1]])
  continuous_simple[[2]] <- sub("_STRAT$", "", continuous[[2]])
  continuous_simple[[3]] <- sub("_OUTCOME$", "", continuous[[3]])

  # Identify combinations with identical variables
  binary_identical <- which(binary_simple[[1]] == binary_simple[[2]] | binary_simple[[1]] == binary_simple[[3]] | binary_simple[[2]] == binary_simple[[3]])
  continuous_identical <- which(continuous_simple[[1]] == continuous_simple[[2]] | continuous_simple[[1]] == continuous_simple[[3]] | continuous_simple[[2]] == continuous_simple[[3]])

  # Remove combinations with identical variables
  binary_simple <- binary_simple[-binary_identical, ]
  continuous_simple <- continuous_simple[-continuous_identical, ]

  if (length(binary_identical) > 0) {
    binary_simple[[1]] <- paste0(binary_simple[[1]], "_EXP")
    binary_simple[[2]] <- paste0(binary_simple[[2]], "_STRAT")
    binary_simple[[3]] <- paste0(binary_simple[[3]], "_CASE_STATUS")
    binary <- binary_simple
  }

  if (length(continuous_identical) > 0) {
    continuous_simple[[1]] <- paste0(continuous_simple[[1]], "_EXP")
    continuous_simple[[2]] <- paste0(continuous_simple[[2]], "_STRAT")
    continuous_simple[[3]] <- paste0(continuous_simple[[3]], "_OUTCOME")
    continuous <- continuous_simple
  }

  removed <- data.table(EXP = character(), STRAT = character(), CASE_STATUS = character(), OUTCOME = character(), REASON_FOR_REMOVAL = character())
  removed2 <- rbind(binary_simple[binary_identical, ], continuous_simple[continuous_identical, ], fill = TRUE)
  if (nrow(removed2) > 0) {
    removed2$REASON_FOR_REMOVAL <- "Identical variables"
  }
  removed <- rbind(removed, removed2, fill = TRUE)
  if (nrow(removed) > 0) {
    message(paste0(nrow(removed), " combinations removed due to identical variables."))
  } else {
    message("No combinations removed due to identical variables.")
  }

  ###############################################################################################
  # Step 5: Z-standardize data
  ###############################################################################################

  standardize <- function(x) {
    (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)
  }

  df_std <- df

  # Select columns to standardize
  cols_to_standardize <- grep("_CASE_STATUS$", variables, value = TRUE, invert = TRUE)

  # Standardize the data
  for (col in cols_to_standardize) {
    df_std[[col]] <- standardize(df[[col]])
  }

  df_exp_std <- df_std[, exp, with = FALSE]
  df_strat_std <- df_std[, strat, with = FALSE]
  df_continuous_outcome_std <- df_std[, continuous_outcome, with = FALSE]

  ###############################################################################################
  # Step 6: Filter out combinations with strongly correlated variables
  ###############################################################################################

  # Calculate the correlation matrix between standardized exposures and stratifying variables
  cor_exp_strat <- cor(df_exp_std, df_strat_std, use = "pairwise.complete.obs", method = "pearson")

  # Change matrices to long format
  cor_exp_strat_long <- as.data.frame(as.table(cor_exp_strat)) # correlation matrix 1
  colnames(cor_exp_strat_long) <- c("EXP", "STRAT", "R")
  
  # Calculate the correlation matrix between standardized continuous outcomes and stratifying variables
  if (nrow(continuous) > 0) {
    cor_outcome_strat <- cor(df_continuous_outcome_std, df_strat_std, use = "pairwise.complete.obs", method = "pearson")
    cor_outcome_strat_long <- as.data.frame(as.table(cor_outcome_strat)) # Correlation matrix 2
    colnames(cor_outcome_strat_long) <- c("OUTCOME", "STRAT", "R")
  }

  # Identify pairs with R > abs(sqrt(0.1))
  cor_exp_strat_high <- cor_exp_strat_long[abs(cor_exp_strat_long$R) > sqrt(0.1), ]
  if (nrow(continuous) > 0) {
    cor_outcome_strat_high <- cor_outcome_strat_long[abs(cor_outcome_strat_long$R) > sqrt(0.1), ]
  }

  # Remove any rows with those pairs in the binary and continuous dfs
  pairs_to_remove_exp_strat <- paste(cor_exp_strat_high$EXP, cor_exp_strat_high$STRAT, sep = "_")

  rows_to_remove_exp_strat_bin <- NULL
  rows_to_remove_exp_strat_cont <- NULL
  rows_to_remove_outcome_strat_cont <- NULL
  
  # nolint start
  if (nrow(binary) > 0) {
    rows_to_remove_exp_strat_bin <- binary %>%
      filter((paste(EXP, STRAT, sep = "_") %in% pairs_to_remove_exp_strat))
    if (nrow(rows_to_remove_exp_strat_bin) > 0) {
      rows_to_remove_exp_strat_bin$REASON_FOR_REMOVAL <- "Exposure and stratifying variable highly correlated (R > sqrt(0.1))"
    }
    binary <- binary %>%
      filter(!(paste(EXP, STRAT, sep = "_") %in% pairs_to_remove_exp_strat))
  }

  if (nrow(continuous) > 0) {
    rows_to_remove_exp_strat_cont <- continuous %>%
      filter((paste(EXP, STRAT, sep = "_") %in% pairs_to_remove_exp_strat))
    if (nrow(rows_to_remove_exp_strat_cont) > 0) {
      rows_to_remove_exp_strat_cont$REASON_FOR_REMOVAL <- "Exposure and stratifying variable highly correlated (R > sqrt(0.1)"
    }
    continuous <- continuous %>%
      filter(!(paste(EXP, STRAT, sep = "_") %in% pairs_to_remove_exp_strat))
    
    pairs_to_remove_outcome_strat <- paste(cor_outcome_strat_high$OUTCOME, cor_outcome_strat_high$STRAT, sep = "_")
    rows_to_remove_outcome_strat_cont <- continuous %>%
      filter((paste(OUTCOME, STRAT, sep = "_") %in% pairs_to_remove_outcome_strat))
    continuous <- continuous %>%
      filter(!(paste(OUTCOME, STRAT, sep = "_") %in% pairs_to_remove_outcome_strat))
    if (nrow(rows_to_remove_outcome_strat_cont) > 0) {
      rows_to_remove_outcome_strat_cont$REASON_FOR_REMOVAL <- "Outcome and stratifying variable highly correlated (R > sqrt(0.1)"
    }
  }
   # nolint end

  # Add filtered out combinations to the removed df only if they exist
  dfs_to_bind <- list(removed)
  if (!is.null(rows_to_remove_exp_strat_bin) && nrow(rows_to_remove_exp_strat_bin) > 0) {
    dfs_to_bind <- append(dfs_to_bind, list(rows_to_remove_exp_strat_bin))
  }
  if (!is.null(rows_to_remove_exp_strat_cont) && nrow(rows_to_remove_exp_strat_cont) > 0) {
    dfs_to_bind <- append(dfs_to_bind, list(rows_to_remove_exp_strat_cont))
  }
  if (!is.null(rows_to_remove_outcome_strat_cont) && nrow(rows_to_remove_outcome_strat_cont) > 0) {
    dfs_to_bind <- append(dfs_to_bind, list(rows_to_remove_outcome_strat_cont))
  }
  removed <- rbindlist(dfs_to_bind, fill = TRUE)

  # Calculate the number of rows removed for each variable
  rows_removed_exp_strat_bin <- if (!is.null(rows_to_remove_exp_strat_bin)) nrow(rows_to_remove_exp_strat_bin) else 0
  rows_removed_exp_strat_cont <- if (!is.null(rows_to_remove_exp_strat_cont)) nrow(rows_to_remove_exp_strat_cont) else 0
  rows_removed_outcome_strat_cont <- if (!is.null(rows_to_remove_outcome_strat_cont)) nrow(rows_to_remove_outcome_strat_cont) else 0

  # Output message
  total_removed <- rows_removed_exp_strat_bin + rows_removed_exp_strat_cont + rows_removed_outcome_strat_cont
  message(total_removed, " combinations with highly correlated outcome-stratifying variable or exposure-stratifying variable pairs removed.")

  ###############################################################################################
  # Step 7: Filter out combinations where variables are diagnostically related
  ###############################################################################################

  if (!is.null(diagnostically_related) && nrow(diagnostically_related) > 0) {
    dr_exp <- diagnostically_related
    colnames(dr_exp) <- c("EXP", "CASE_STATUS")
    dr_exp$EXP <- paste0(dr_exp$EXP, "_EXP")
    dr_exp$CASE_STATUS <- paste0(dr_exp$CASE_STATUS, "_CASE_STATUS")

    dr_strat <- diagnostically_related
    colnames(dr_strat) <- c("STRAT", "CASE_STATUS")
    dr_strat$STRAT <- paste0(dr_strat$STRAT, "_STRAT")
    dr_strat$CASE_STATUS <- paste0(dr_strat$CASE_STATUS, "_CASE_STATUS")

    pairs_to_remove_exp_case <- paste(dr_exp$EXP, dr_exp$CASE_STATUS, sep = "_")
    # nolint start
    rows_to_remove_exp_case <- binary %>%
      filter((paste(EXP, CASE_STATUS, sep = "_") %in% pairs_to_remove_exp_case))
    binary <- binary %>%
      filter(!(paste(EXP, CASE_STATUS, sep = "_") %in% pairs_to_remove_exp_case))

    pairs_to_remove_strat_case <- paste(dr_strat$STRAT, dr_strat$CASE_STATUS, sep = "_")
    rows_to_remove_strat_case <- binary %>%
      filter((paste(STRAT, CASE_STATUS, sep = "_") %in% pairs_to_remove_strat_case))
    binary <- binary %>%
      filter(!(paste(STRAT, CASE_STATUS, sep = "_") %in% pairs_to_remove_strat_case))
    # nolint end
    if (nrow(rows_to_remove_exp_case) > 0) {
      rows_to_remove_exp_case$REASON_FOR_REMOVAL <- "Exposure and binary outcome diagnostically related"
    }
    if (nrow(rows_to_remove_strat_case) > 0) {
      rows_to_remove_strat_case$REASON_FOR_REMOVAL <- "Stratifying variable and binary outcome diagnostically related"
    }

    # Add filtered out combinations to the removed df only if they exist
    removed <- rbind(removed, rows_to_remove_exp_case, rows_to_remove_strat_case, fill = TRUE)

    message(nrow(rows_to_remove_exp_case) + nrow(rows_to_remove_strat_case), " combinations removed due to diagnostically related pairs.")
  } else {
    message("No diagnostically related pairs provided.")
  }

  ###############################################################################################
  # Step 8: Calculate the overall causal effect of the exposures on binary outcomes
  ###############################################################################################

  if (nrow(binary) > 0) {
    # Create combinations of EXP and CASE_STATUS
    binary_overall <- expand.grid(exp, binary_outcome)
    binary_overall <- as.data.table(binary_overall)
    colnames(binary_overall) <- c("EXP", "CASE_STATUS")

    # Ensure that they are characters and not factors
    binary_overall$EXP <- as.character(binary_overall$EXP)
    binary_overall$CASE_STATUS <- as.character(binary_overall$CASE_STATUS)

    # Create new columns to store Wald Ratio values
    binary_overall$wald_den <- NA
    binary_overall$wald_num <- NA
    binary_overall$wald_ratio_OR <- NA
    binary_overall$wald_ratio_OR_95_LL <- NA
    binary_overall$wald_ratio_OR_95_UL <- NA
    binary_overall$X_G_R2 <- NA

    # Make sure the vector of covariates is unnamed
    covar_cols <- unname(covar)

    # Calculate overall causal effects using the Wald Ratio
    for (i in seq_len(nrow(binary_overall))) {

      # Extract EXP and OUTCOME column names from row
      exp_var <- binary_overall$EXP[i]
      case_status_var <- binary_overall$CASE_STATUS[i]
      prs_var <- paste0(sub("_EXP", "", exp_var), "_PRS")

      # Conditionally remove SEX_COVAR if outcome is SEX_CASE_STATUS (negative control analyses)
      covar_to_use <- if (case_status_var == "SEX_CASE_STATUS") covar_cols[covar_cols != "SEX_COVAR"] else covar_cols

      # Run the denominator model: EXP ~ PRS + covariates
      model_den <- lm(as.formula(paste(exp_var, "~", prs_var, "+ .")), data = df[, c(exp_var, prs_var, covar_to_use), with = FALSE]) # note that covariate columns are optional
      coef_den <- model_den$coef[2]

      # Run the numerator model: CASE_STATUS ~ PRS + covariates
      model_num <- glm(as.formula(paste(case_status_var, "~", prs_var, "+ .")), data = df[, c(case_status_var, prs_var, covar_to_use), with = FALSE],
                       family = binomial) # note that covariate columns are optional
      coef_num <- model_num$coef[2]
      se_num <- coef(summary(model_num))[2, "Std. Error"]

      # Calculate Wald Ratio
      wald_ratio <- coef_num / coef_den
      wald_ratio_se <- se_num / coef_den # uses delta approximation method

      # Calculate 95% CI for the Wald Ratio OR
      wald_ratio_or <- exp(wald_ratio)
      wald_ratio_or_95_ll <- exp(wald_ratio - 1.96 * wald_ratio_se)
      wald_ratio_or_95_ul <- exp(wald_ratio + 1.96 * wald_ratio_se)

      # Calculate instrument strength
      model_x_g <- lm(as.formula(paste(exp_var, "~", prs_var)), data = df[, c(exp_var, prs_var), with = FALSE]) # no covariates included
      x_g_r2 <- summary(model_x_g)$r.squared # variance explained as a calculation of instrument strength

      # Store the model summaries in binary_overall
      binary_overall$wald_den[i] <- coef_den
      binary_overall$wald_num[i] <- coef_num
      binary_overall$wald_ratio_OR[i] <- wald_ratio_or
      binary_overall$wald_ratio_OR_95_LL[i] <- wald_ratio_or_95_ll
      binary_overall$wald_ratio_OR_95_UL[i] <- wald_ratio_or_95_ul
      binary_overall$X_G_R2[i] <- x_g_r2

      # Print iteration
      message("Wald ratio calculation complete for exposure-binary outcome pair number ", i)
    }

    message("Wald ratio calculations complete for all exposure-binary outcome pairs.")

    # Identify rows where the Wald Ratio OR 95% CI include 1
    binary_overall_filter <- which((binary_overall$wald_ratio_OR_95_LL < 1 & binary_overall$wald_ratio_OR_95_UL > 1) | (binary_overall$wald_ratio_OR_95_UL < 1 & binary_overall$wald_ratio_OR_95_LL > 1))
    binary_overall_filter <- binary_overall[binary_overall_filter, ]

    # Filter out any of those rows out of binary with null overall causal effects
    pairs_to_remove_exp_case <- paste(binary_overall_filter$EXP, binary_overall_filter$CASE_STATUS, sep = "_")
    # nolint start
    rows_to_remove_exp_case <- binary %>%
      filter((paste(EXP, CASE_STATUS, sep = "_") %in% pairs_to_remove_exp_case))
    if (nrow(rows_to_remove_exp_case) > 0) {
      rows_to_remove_exp_case$REASON_FOR_REMOVAL <- "Overall causal effect of X on Y is null"
    }
    binary <- binary %>%
      filter(!(paste(EXP, CASE_STATUS, sep = "_") %in% pairs_to_remove_exp_case))
    # nolint end

    # Add filtered out combinations to the removed df
    removed <- rbind(removed, rows_to_remove_exp_case, fill = TRUE)

    message(nrow(rows_to_remove_exp_case), " combinations with a null causal effect removed (binary outcomes).")
  } else {
    message("No combinations with a null causal effect removed (binary outcomes), as there are no binary outcomes.")
  }
  
  ###############################################################################################
  # Step 9: Calculate the overall causal effect of the exposures on continuous outcomes
  ###############################################################################################

  if (nrow(continuous) > 0) {
    # Create combinations of EXP and OUTCOME
    continuous_overall <- expand.grid(exp, continuous_outcome)
    continuous_overall <- as.data.table(continuous_overall)
    colnames(continuous_overall) <- c("EXP", "OUTCOME")

    # Ensure that they are characters and not factors
    continuous_overall$EXP <- as.character(continuous_overall$EXP)
    continuous_overall$OUTCOME <- as.character(continuous_overall$OUTCOME)

    # Create new columns to store Wald Ratio values
    continuous_overall$wald_den <- NA
    continuous_overall$wald_num <- NA
    continuous_overall$wald_ratio <- NA
    continuous_overall$wald_ratio_se <- NA
    continuous_overall$wald_ratio_95_LL <- NA
    continuous_overall$wald_ratio_95_UL <- NA
    continuous_overall$X_G_R2 <- NA

    # Make sure the vector of covariates is unnamed
    covar_cols <- unname(covar)

    # Calculate overall causal effects using the Wald Ratio
    for (i in seq_len(nrow(continuous_overall))) {

      # Extract EXP and OUTCOME column names from row
      exp_var <- continuous_overall$EXP[i]
      outcome_var <- continuous_overall$OUTCOME[i]
      prs_var <- paste0(sub("_EXP", "", exp_var), "_PRS")

      # Conditionally remove AGE_COVAR if outcome is AGE_OUTCOME (negative control analyses)
      covar_to_use <- if (outcome_var == "AGE_OUTCOME") covar_cols[covar_cols != "AGE_COVAR"] else covar_cols

      # Run the denominator model: EXP ~ PRS + covariates
      model_den <- lm(as.formula(paste(exp_var, "~", prs_var, "+ .")), data = df[, c(exp_var, prs_var, covar_to_use), with = FALSE]) # note that covariate columns are optional
      coef_den <- model_den$coef[2]

      # Run the numerator model: OUTCOME ~ PRS + covariates
      model_num <- lm(as.formula(paste(outcome_var, "~", prs_var, "+ .")), data = df[, c(outcome_var, prs_var, covar_to_use), with = FALSE]) # note that covariate columns are optional
      coef_num <- model_num$coef[2]
      se_num <- coef(summary(model_num))[2, "Std. Error"]

      # Calculate Wald Ratio
      wald_ratio <- coef_num / coef_den
      wald_ratio_se <- se_num / coef_den # uses delta approximation method

      # Calculate 95% CI
      wald_ratio_95_ll <- wald_ratio - 1.96 * wald_ratio_se
      wald_ratio_95_ul <- wald_ratio + 1.96 * wald_ratio_se

      # Calculate instrument strength
      model_x_g <- lm(as.formula(paste(exp_var, "~", prs_var)), data = df[, c(exp_var, prs_var), with = FALSE]) # no covariates included
      x_g_r2 <- summary(model_x_g)$r.squared # variance explained as a calculation of instrument strength

      # Store the model summaries in continuous_overall
      continuous_overall$wald_den[i] <- coef_den
      continuous_overall$wald_num[i] <- coef_num
      continuous_overall$wald_ratio[i] <- wald_ratio
      continuous_overall$wald_ratio_se[i] <- wald_ratio_se
      continuous_overall$wald_ratio_95_LL[i] <- wald_ratio_95_ll
      continuous_overall$wald_ratio_95_UL[i] <- wald_ratio_95_ul
      continuous_overall$X_G_R2[i] <- x_g_r2

      message("Wald ratio calculation complete for exposure-continuous outcome pair number ", i)
    }

    message("Wald ratio calculations complete for all exposure-continuous outcome pairs.")

    # Identify rows where the Wald Ratio 95% CI include 0
    continuous_overall_filter <- which((continuous_overall$wald_ratio_95_LL < 0 & continuous_overall$wald_ratio_95_UL > 0) | (continuous_overall$wald_ratio_95_UL < 0 & continuous_overall$wald_ratio_95_LL > 0))
    continuous_overall_filter <- continuous_overall[continuous_overall_filter, ]

    # Filter out any of those rows out of continuous with null overall causal effects
    pairs_to_remove_exp_outcome <- paste(continuous_overall_filter$EXP, continuous_overall_filter$OUTCOME, sep = "_")
    # nolint start
    rows_to_remove_exp_outcome <- continuous %>%
      filter((paste(EXP, OUTCOME, sep = "_") %in% pairs_to_remove_exp_outcome))
    if (nrow(rows_to_remove_exp_outcome) > 0) {
      rows_to_remove_exp_outcome$REASON_FOR_REMOVAL <- "Overall causal effect of X on Y is null"
    }
    continuous <- continuous %>%
      filter(!(paste(EXP, OUTCOME, sep = "_") %in% pairs_to_remove_exp_outcome))
    # nolint end

    # Add filtered out combinations to the removed df
    removed <- rbind(removed, rows_to_remove_exp_outcome, fill = TRUE)

    message(nrow(rows_to_remove_exp_outcome), " combinations with a null causal effect removed (continuous outcomes).")
  } else {
    message("No combinations with a null causal effect removed (continuous outcomes), as there are no continuous outcomes.")
  }

  ###############################################################################################
  # Step 10: Calculate all relevant "biases"
  ###############################################################################################

  # Herein we use only standardized data

  # Calculate quadratic stratifying variables (C2), exposures (X2), and continuous outcomes (Y2)
  to_square <- c(exp, strat, continuous_outcome)
  for (var in to_square) {
    var2 <- paste0(var, "2") # add 2 to the end of the variable name
    df_std[[var2]] <- standardize(df_std[[var]]^2) # Z-standardize the quadratic term
  }

  bias_lm <- function(var1, var2, var1_name, var2_name, name, gstrat = FALSE) {
    # Create combinations of var1 and var2
    grid <- as.data.table(expand.grid(var1, var2))
    colnames(grid) <- c(var1_name, var2_name)

    # Ensure that they are characters and not factors
    grid[[var1_name]] <- as.character(grid[[var1_name]])
    grid[[var2_name]] <- as.character(grid[[var2_name]])

    # Create new columns to store R2 values
    beta_name <- paste0("beta_", name)
    grid[[beta_name]] <- NA
    beta_se_name <- paste0("beta_", name, "_se")
    grid[[beta_se_name]] <- NA
    beta2_name <- paste0("beta_", name, "2")
    grid[[beta2_name]] <- NA
    beta2_se_name <- paste0("beta_", name, "2_se")
    grid[[beta2_se_name]] <- NA
    r2_suffix <- if (gstrat) "_R2" else "2_R2"
    r2_name <- paste0(name, r2_suffix)
    grid[[r2_name]] <- NA

    for (i in seq_len(nrow(grid))) {
      # Extract relevant variables from row
      v1 <- grid[[var1_name]][i]
      v2 <- grid[[var2_name]][i]
      if (gstrat) {
        prs_var <- sub("_EXP$", "_PRS", v1)
        gstrat_var <- paste0(prs_var, "_", v2)

        # Create gstrat_var in df_std
        df_std[[gstrat_var]] <- df_std[[prs_var]] * df_std[[v2]]
        df_std[[gstrat_var]] <- standardize(df_std[[gstrat_var]]) # standardize
        # Run the model: var1 ~ prs_var + v2 + gstrat_var
        model <- lm(as.formula(paste(v1, "~", prs_var, "+", v2, "+", gstrat_var)), data = df_std[, c(v1, prs_var, v2, gstrat_var), with = FALSE])
      } else {
        v2_2 <- paste0(v2, "2")
        # Run the model: var1 ~ var2 + var2_2
        model <- lm(as.formula(paste(v1, "~", v2, "+", v2_2)), data = df_std[, c(v1, v2, v2_2), with = FALSE])
      }

      beta <- coef(model)[2]
      beta_se <- summary(model)$coef[2, 2]
      if (gstrat) {
        beta2 <- coef(model)[4]
        beta2_se <- summary(model)$coef[4, 2]
      } else {
        beta2 <- coef(model)[3]
        beta2_se <- summary(model)$coef[3, 2]
      }
      r2 <- beta2^2 # if all involved variables have a SD of 1, it is possible to approximate the R2 using coefficient^2

      # Store model summaries in grid
      grid[[beta_name]][i] <- beta
      grid[[beta_se_name]][i] <- beta_se
      grid[[beta2_name]][i] <- beta2
      grid[[beta2_se_name]][i] <- beta2_se
      grid[[r2_name]][i] <- r2
    }

    return(grid)
  }

  bias_glm <- function(var1, var2, var1_name, var2_name, name) {
    # Create combinations of var1 and var2
    grid <- as.data.table(expand.grid(var1, var2))
    colnames(grid) <- c(var1_name, var2_name)

    # Ensure that they are characters and not factors
    grid[[var1_name]] <- as.character(grid[[var1_name]])
    grid[[var2_name]] <- as.character(grid[[var2_name]])

    # Create new columns to store OR values
    or_name <- paste0(name, "_OR")
    grid[[or_name]] <- NA
    or_95_ll_name <- paste0(name, "_OR_95_LL")
    grid[[or_95_ll_name]] <- NA
    or_95_ul_name <- paste0(name, "_OR_95_UL")
    grid[[or_95_ul_name]] <- NA
    or2_name <- paste0(name, "2_OR")
    grid[[or2_name]] <- NA
    or2_95_ll_name <- paste0(name, "2_OR_95_LL")
    grid[[or2_95_ll_name]] <- NA
    or2_95_ul_name <- paste0(name, "2_OR_95_UL")
    grid[[or2_95_ul_name]] <- NA

    for (i in seq_len(nrow(grid))) {
      # Extract relevant variables from row
      v1 <- grid[[var1_name]][i]
      v2 <- grid[[var2_name]][i]
      v2_2 <- paste0(v2, "2")

      # Run the model: var1 ~ var2 + var2_2
      model <- glm(as.formula(paste(v1, "~", v2, "+", v2_2)), data = df_std[, c(v1, v2, v2_2), with = FALSE], family = binomial)
      beta <- coef(model)[2]
      beta_se <- summary(model)$coef[2, 2]
      beta2 <- coef(model)[3]
      beta2_se <- summary(model)$coef[3, 2]

      # Calculate OR
      or <- exp(beta)
      or_95_ll <- exp(beta - 1.96 * beta_se)
      or_95_ul <- exp(beta + 1.96 * beta_se)
      or2 <- exp(beta2)
      or2_95_ll <- exp(beta2 - 1.96 * beta2_se)
      or2_95_ul <- exp(beta2 + 1.96 * beta2_se)

      # Store model summaries in grid
      grid[[or_name]][i] <- or
      grid[[or_95_ll_name]][i] <- or_95_ll
      grid[[or_95_ul_name]][i] <- or_95_ul
      grid[[or2_name]][i] <- or2
      grid[[or2_95_ll_name]][i] <- or2_95_ll
      grid[[or2_95_ul_name]][i] <- or2_95_ul
    }

    return(grid)
  }

  x_c <- bias_lm(exp, strat, "EXP", "STRAT", "X_C")

  if (nrow(continuous) > 0) {
    c_y <- bias_lm(strat, continuous_outcome, "STRAT", "OUTCOME", "C_Y")
    c_y_c <- bias_lm(continuous_outcome, strat, "OUTCOME", "STRAT", "Y_C")
    c_y_x <- bias_lm(continuous_outcome, exp, "OUTCOME", "EXP", "Y_X")
    c_x <- bias_lm(strat, exp, "STRAT", "EXP", "C_X")
    x_gc <- bias_lm(exp, strat, "EXP", "STRAT", "X_GC", gstrat = TRUE)
  }

  if (nrow(binary) > 0) {
    b_y_c <- bias_glm(binary_outcome, strat, "CASE_STATUS", "STRAT", "Y_C")
    b_y_x <- bias_glm(binary_outcome, exp, "CASE_STATUS", "EXP", "Y_X")
    x_gc <- bias_lm(exp, strat, "EXP", "STRAT", "X_GC", gstrat = TRUE)
  }

  message("All biases calculated.")

  ###############################################################################################
  # Step 11: Append the "biases" regression results to the appropriate files and print
  ###############################################################################################

  # Binary outcomes
  if (nrow(binary) > 0) {
    # Append X_C2_R2
    binary_1 <- merge(binary, x_c[, c("EXP", "STRAT", "X_C2_R2")], by = c("EXP", "STRAT"), all.x = TRUE)

    # Append Y_C2_OR
    binary_2 <- merge(binary_1, b_y_c[, c("CASE_STATUS", "STRAT", "Y_C2_OR")], by = c("CASE_STATUS", "STRAT"), all.x = TRUE)

    # Append Y_X2_OR
    binary_3 <- merge(binary_2, b_y_x[, c("CASE_STATUS", "EXP", "Y_X2_OR")], by = c("CASE_STATUS", "EXP"), all.x = TRUE)

    # Append X_GC_R2
    x_gc <- x_gc[, c("EXP", "STRAT", "X_GC_R2")]
    binary_4 <- merge(binary_3, x_gc, by = c("EXP", "STRAT"), all.x = TRUE)

    # Append Wald ratios
    wald_ratios <- binary_overall[, c("EXP", "CASE_STATUS", "wald_ratio_OR", "wald_ratio_OR_95_LL", "wald_ratio_OR_95_UL", "X_G_R2")]
    binary_5 <- merge(binary_4, wald_ratios, by = c("EXP", "CASE_STATUS"), all.x = TRUE)

    # Reorder file
    order <- c("EXP", "STRAT", "CASE_STATUS", "X_C2_R2", "Y_C2_OR", "Y_X2_OR", "X_GC_R2", "wald_ratio_OR", "wald_ratio_OR_95_LL", "wald_ratio_OR_95_UL", "X_G_R2")
    setcolorder(binary_5, order)
    setorder(binary_5, EXP, CASE_STATUS, STRAT) # nolint
    binary_final <- binary_5
  }

  # Continuous outcomes
  if (nrow(continuous) > 0) {
    # Append C_Y2_R2
    continuous_1 <- merge(continuous, c_y[, c("STRAT", "OUTCOME", "C_Y2_R2")], by = c("STRAT", "OUTCOME"), all.x = TRUE)

    # Append X_C2_R2
    x_c <- x_c[, c("EXP", "STRAT", "X_C2_R2")]
    continuous_2 <- merge(continuous_1, x_c, by = c("EXP", "STRAT"), all.x = TRUE)

    message("Ordering columns...")
    # Append Y_C2_R2
    y_c <- c_y_c[, c("OUTCOME", "STRAT", "Y_C2_R2")]
    continuous_3 <- merge(continuous_2, y_c, by = c("OUTCOME", "STRAT"), all.x = TRUE)

    # Append Y_X2_R2
    y_x <- c_y_x[, c("OUTCOME", "EXP", "Y_X2_R2")]
    continuous_4 <- merge(continuous_3, y_x, by = c("OUTCOME", "EXP"), all.x = TRUE)

    # Append C_X2_R2
    c_x <- c_x[, c("STRAT", "EXP", "C_X2_R2")]
    continuous_5 <- merge(continuous_4, c_x, by = c("STRAT", "EXP"), all.x = TRUE)

    # Append X_GC_R2
    x_gc <- x_gc[, c("EXP", "STRAT", "X_GC_R2")]
    continuous_6 <- merge(continuous_5, x_gc, by = c("EXP", "STRAT"), all.x = TRUE)

    # Append wald ratios
    wald_ratios <- continuous_overall[, c("EXP", "OUTCOME", "wald_ratio", "wald_ratio_95_LL", "wald_ratio_95_UL", "X_G_R2")]
    continuous_7 <- merge(continuous_6, wald_ratios, by = c("EXP", "OUTCOME"), all.x = TRUE)

    # Reorder file
    order <- c("EXP", "STRAT", "OUTCOME", "C_Y2_R2", "X_GC_R2", "X_C2_R2", "Y_C2_R2", "Y_X2_R2", "C_X2_R2", "wald_ratio", "wald_ratio_95_LL", "wald_ratio_95_UL", "X_G_R2")
    setcolorder(continuous_7, order)
    setorder(continuous_7, EXP, OUTCOME, STRAT) # nolint
    continuous_final <- continuous_7
  }

  if (!is.null(output_dir)) {
    # Print files
    if (nrow(binary) > 0) {
      fwrite(binary_final, file = file.path(output_dir, "binary_comb.txt"), sep = "\t", quote = FALSE, na = "NA")
    }
    if (nrow(continuous) > 0) {
      fwrite(continuous_final, file = file.path(output_dir, "continuous_comb.txt"), sep = "\t", quote = FALSE, na = "NA")
    }

    message("Biases added to the combination list and outputted.")

    ###############################################################################################
    # Step 12: Print overall causal estimates and the strength of the instruments
    ###############################################################################################

    if (nrow(binary) > 0) {
      fwrite(binary_overall, file = file.path(output_dir, "binary_overall.txt"), sep = "\t", quote = FALSE, na = "NA")
    }

    if (nrow(continuous) > 0) {
      fwrite(continuous_overall, file = file.path(output_dir, "continuous_overall.txt"), sep = "\t", quote = FALSE, na = "NA")
    }

    if (nrow(removed) > 0) {
      fwrite(removed, file = file.path(output_dir, "removed_combinations.txt"), sep = "\t", quote = FALSE, na = "NA")
    }
  }

  message("Overall causal estimates and instrument strength outputted. Removed combinations printed as well.")

  message("PRE-PROCESSING COMPLETE")
  print(Sys.time() - start)

  if (nrow(binary) == 0) {
    binary_final <- NULL
  }
  if (nrow(continuous) == 0) {
    continuous_final <- NULL
  }

  return(list(binary = binary_final, continuous = continuous_final))
}
