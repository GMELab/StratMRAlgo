library(data.table)

#' @title Main algorithm
#' @param quantile_number The number of quantiles or strata to use.
#' @param input_dir The directory where the data is stored, defaults to the current working directory.
#' @param output_dir The directory where the output should be stored, defaults to NULL meaning no output is saved.
#' @param data The name of the data file or an input data.table, defaults to "data.csv". Column names need to be in a specific format; see README page for details.
#' @param diagnostically_related The file name of the diagnostically related variables, or a data frame with two columns. No column headers required; see README page for details. Defaults to NULL meaning no diagnostically related variables are provided. Parameter used for preprocessing step.
#' @param binary The file name of the binary combinations, or a data frame with three columns (EXP, STRAT, CASE_STATUS). Column headers required. Defaults to NULL meaning the combinations will be automatically generated. Parameter used for preprocessing step.
#' @param continuous The file name of the continuous combinations, or a data frame with three columns (EXP, STRAT, OUTCOME). Column headers required. Defaults to NULL meaning the combinations will be automatically generated. Parameter used for preprocessing step.
#' @param n The number of combinations in the dataset (used for calculating the Bonferroni threshold), defaults to NULL meaning the number of observations will be automatically calculated.
#' @param causal_relationships The file name of the causal relationships, or a data frame with two columns defining a causal relationship where C1 is known to cause C2 (C1 -> C2). Headers are required (named C1 and C2). This defines the causal and collider relationships between C and X in the actual algorithm. If NULL, or there is no entry for C1 -> C2, then the causal relationship is unknown. Defaults to NULL meaning no causal relationships are provided.
#' @param combinations Skip preprocessing and pass in data frames (combinations$binary and combinations$continuous) returned from the preprocess function. Defaults to NULL meaning it will preprocess the provided data.
#' @param p_threshold The threshold for statistical significance, defaults to 0.05.
#' @return A list containing two data tables, one for binary results and one for continuous results. See README page for details.
#' @export
algo <- function(
    quantile_number,
    input_dir = getwd(),
    output_dir = NULL,
    data = "data.csv",
    diagnostically_related = NULL,
    binary = NULL,
    continuous = NULL,
    n = NULL,
    causal_relationships = NULL,
    combinations = NULL,
    p_threshold = 0.05
  ) {
  quantile_number <- as.integer(quantile_number)
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
  } else if (!is.data.frame(data)) {
    stop("data must be a file name or a data frame.")
  } else {
    df <- data
  }

  if (is.null(combinations)) {
    combinations <- preprocess(input_dir = input_dir, data = df, diagnostically_related = diagnostically_related, binary = binary, continuous = continuous, output_dir = output_dir) # nolint
    binary <- combinations$binary
    continuous <- combinations$continuous
  } else {
    binary <- combinations$binary
    continuous <- combinations$continuous
    message("Preprocessing steps skipped, using combinations$binary and/or combinations$continuous directly from R environment.")
  }

  binary_n <- if (is.null(binary)) {
    0
  } else {
    nrow(binary)
  }
  continuous_n <- if (is.null(continuous)) {
    0
  } else {
    nrow(continuous)
  }
  n <- if (is.null(n)) {
    binary_n + continuous_n
  } else {
    n
  }
  message(paste("Total number of combinations used for Bonferroni correction =", n))
  
  # Load in the causal_relationships, if provided
  if (!is.null(causal_relationships)) {
    if (is.character(causal_relationships)) {
      if (!file.exists(file.path(input_dir, causal_relationships))) {
        stop("causal relationships file does not exist.")
      }
      causal_relationships <- fread(
        file.path(input_dir, causal_relationships),
        stringsAsFactors = FALSE,
        header = TRUE,
        data.table = TRUE,
      )
    } else {
      if (!is.data.table(causal_relationships)) {
        stop("causal_relationships must be a data.table.")
      }
    }
    message("Loaded causal relationships.")
  } else {
    causal_relationships <- NULL
    message("No causal relationships provided.")
  }

  binary_results <- data.table(
    EXP = character(),
    STRAT = character(),
    CASE_STATUS = character(),
    SIGNIFICANT = character(), # TRUE if significant, FALSE if not significant, NULL if inconclusive
    X_C2_R2 = numeric(),
    Y_C2_OR = numeric(),
    Y_X2_OR = numeric(),
    DR_BY_X_P = numeric(),
    DR_BY_X_BETA = numeric(),
    DR_BY_X_SE = numeric(),
    DR_BY_G_P = numeric(),
    DR_BY_G_BETA = numeric(),
    DR_BY_G_SE = numeric(),

    WALD_RATIO_OR = numeric(),
    WALD_RATIO_OR_95_LL = numeric(),
    WALD_RATIO_OR_95_UL = numeric(),
    X_G_R2 = numeric()
  )
  continuous_results <- data.table(
    EXP = character(),
    STRAT = character(),
    OUTCOME = character(),
    SIGNIFICANT = character(), # TRUE if significant, FALSE if not significant, NULL if inconclusive
    C_Y2_R2 = numeric(),
    X_GC_R2 = numeric(),
    X_C2_R2 = numeric(),
    Y_C2_R2 = numeric(),
    Y_X2_R2 = numeric(),
    C_X2_R2 = numeric(),
    DR_BY_X_P = numeric(),
    DR_BY_X_BETA = numeric(),
    DR_BY_X_SE = numeric(),
    DR_BY_G_P = numeric(),
    DR_BY_G_BETA = numeric(),
    DR_BY_G_SE = numeric(),
    RESID_STRAT_P = numeric(),
    RESID_STRAT_BETA = numeric(),
    RESID_STRAT_SE = numeric(),

    WALD_RATIO = numeric(),
    WALD_RATIO_95_LL = numeric(),
    WALD_RATIO_95_UL = numeric(),
    X_G_R2 = numeric()
  )

  # Binary algorithm
  if (binary_n > 0) {
    for (i in seq_len(nrow(binary))) {
      c <- binary[i, ]
      exp <- c$EXP
      strat <- c$STRAT
      outcome <- c$CASE_STATUS

      # X -> C
      c_is_collider <- if (is.null(causal_relationships)) {
        FALSE
      } else {
        any(causal_relationships$C1 == sub("_EXP$", "", exp) & 
            causal_relationships$C2 == sub("_STRAT$", "", strat)
            )
      }
      # C -> X
      c_is_causal <- if (is.null(causal_relationships)) {
        FALSE
      } else {
        any(causal_relationships$C1 == sub("_STRAT$", "", strat) & 
            causal_relationships$C2 == sub("_EXP$", "", exp)
            )
      }

      dr_by_x <- mr(
        method = "by_x",
        exposure_name = exp,
        outcome_name = outcome,
        stratifying_name = strat,
        quantile_number = quantile_number,
        data = df,
        output_dir = output_dir
      )
      if (dr_by_x$p < (p_threshold / n)) {
        if (c_is_collider) {
          binary_results <- rbind(binary_results, data.table(
            EXP = exp,
            STRAT = strat,
            CASE_STATUS = outcome,
            SIGNIFICANT = "TRUE",
            X_C2_R2 = c$X_C2_R2,
            Y_C2_OR = c$Y_C2_OR,
            Y_X2_OR = c$Y_X2_OR,
            DR_BY_X_P = dr_by_x$p,
            DR_BY_X_BETA = dr_by_x$beta,
            DR_BY_X_SE = dr_by_x$se,
            DR_BY_G_P = NA,
            DR_BY_G_BETA = NA,
            DR_BY_G_SE = NA,

            WALD_RATIO_OR = c$wald_ratio_OR,
            WALD_RATIO_OR_95_LL = c$wald_ratio_OR_95_LL,
            WALD_RATIO_OR_95_UL = c$wald_ratio_OR_95_UL,
            X_G_R2 = c$X_G_R2
          ), fill = TRUE)
        } else { # c_is_causal or unclear_causal_relationship
          m1_c2_r2 <- c$X_C2_R2
          m2_c2_or <- c$Y_C2_OR
          if (m1_c2_r2 <= 0.02 && m2_c2_or <= 1.05 && m2_c2_or >= 1 / 1.05) {
            binary_results <- rbind(binary_results, data.table(
              EXP = exp,
              STRAT = strat,
              CASE_STATUS = outcome,
              SIGNIFICANT = "TRUE",
              X_C2_R2 = c$X_C2_R2,
              Y_C2_OR = c$Y_C2_OR,
              Y_X2_OR = c$Y_X2_OR,
              DR_BY_X_P = dr_by_x$p,
              DR_BY_X_BETA = dr_by_x$beta,
              DR_BY_X_SE = dr_by_x$se,
              DR_BY_G_P = NA,
              DR_BY_G_BETA = NA,
              DR_BY_G_SE = NA,

              WALD_RATIO_OR = c$wald_ratio_OR,
              WALD_RATIO_OR_95_LL = c$wald_ratio_OR_95_LL,
              WALD_RATIO_OR_95_UL = c$wald_ratio_OR_95_UL,
              X_G_R2 = c$X_G_R2
            ), fill = TRUE)
          } else {
            m3_x2_or <- c$Y_X2_OR
            if (m3_x2_or <= 1.02 && m3_x2_or >= 1 / 1.02) {
              dr_by_g <- mr(
                method = "by_g",
                exposure_name = exp,
                outcome_name = outcome,
                stratifying_name = strat,
                quantile_number = quantile_number,
                data = df,
                output_dir = output_dir
              )
              if (dr_by_g$p < (p_threshold / n)) {
                binary_results <- rbind(binary_results, data.table(
                  EXP = exp,
                  STRAT = strat,
                  CASE_STATUS = outcome,
                  SIGNIFICANT = "TRUE",
                  X_C2_R2 = c$X_C2_R2,
                  Y_C2_OR = c$Y_C2_OR,
                  Y_X2_OR = c$Y_X2_OR,
                  DR_BY_X_P = dr_by_x$p,
                  DR_BY_X_BETA = dr_by_x$beta,
                  DR_BY_X_SE = dr_by_x$se,
                  DR_BY_G_P = dr_by_g$p,
                  DR_BY_G_BETA = dr_by_g$beta,
                  DR_BY_G_SE = dr_by_g$se,

                  WALD_RATIO_OR = c$wald_ratio_OR,
                  WALD_RATIO_OR_95_LL = c$wald_ratio_OR_95_LL,
                  WALD_RATIO_OR_95_UL = c$wald_ratio_OR_95_UL,
                  X_G_R2 = c$X_G_R2
                ), fill = TRUE)
              } else {
                binary_results <- rbind(binary_results, data.table(
                  EXP = exp,
                  STRAT = strat,
                  CASE_STATUS = outcome,
                  SIGNIFICANT = "FALSE",
                  X_C2_R2 = c$X_C2_R2,
                  Y_C2_OR = c$Y_C2_OR,
                  Y_X2_OR = c$Y_X2_OR,
                  DR_BY_X_P = dr_by_x$p,
                  DR_BY_X_BETA = dr_by_x$beta,
                  DR_BY_X_SE = dr_by_x$se,
                  DR_BY_G_P = dr_by_g$p,
                  DR_BY_G_BETA = dr_by_g$beta,
                  DR_BY_G_SE = dr_by_g$se,

                  WALD_RATIO_OR = c$wald_ratio_OR,
                  WALD_RATIO_OR_95_LL = c$wald_ratio_OR_95_LL,
                  WALD_RATIO_OR_95_UL = c$wald_ratio_OR_95_UL,
                  X_G_R2 = c$X_G_R2
                ), fill = TRUE)
              }
            } else {
              binary_results <- rbind(binary_results, data.table(
                EXP = exp,
                STRAT = strat,
                CASE_STATUS = outcome,
                SIGNIFICANT = "INCONCLUSIVE",
                X_C2_R2 = c$X_C2_R2,
                Y_C2_OR = c$Y_C2_OR,
                Y_X2_OR = c$Y_X2_OR,
                DR_BY_X_P = dr_by_x$p,
                DR_BY_X_BETA = dr_by_x$beta,
                DR_BY_X_SE = dr_by_x$se,
                DR_BY_G_P = NA,
                DR_BY_G_BETA = NA,
                DR_BY_G_SE = NA,

                WALD_RATIO_OR = c$wald_ratio_OR,
                WALD_RATIO_OR_95_LL = c$wald_ratio_OR_95_LL,
                WALD_RATIO_OR_95_UL = c$wald_ratio_OR_95_UL,
                X_G_R2 = c$X_G_R2
              ), fill = TRUE)
            }
          }
        }
      } else {
        binary_results <- rbind(binary_results, data.table(
          EXP = exp,
          STRAT = strat,
          CASE_STATUS = outcome,
          SIGNIFICANT = "FALSE",
          X_C2_R2 = c$X_C2_R2,
          Y_C2_OR = c$Y_C2_OR,
          Y_X2_OR = c$Y_X2_OR,
          DR_BY_X_P = dr_by_x$p,
          DR_BY_X_BETA = dr_by_x$beta,
          DR_BY_X_SE = dr_by_x$se,
          DR_BY_G_P = NA,
          DR_BY_G_BETA = NA,
          DR_BY_G_SE = NA,

          WALD_RATIO_OR = c$wald_ratio_OR,
          WALD_RATIO_OR_95_LL = c$wald_ratio_OR_95_LL,
          WALD_RATIO_OR_95_UL = c$wald_ratio_OR_95_UL,
          X_G_R2 = c$X_G_R2
        ), fill = TRUE)
      }
    }
  }
  if (continuous_n > 0) {
    for (i in seq_len(nrow(continuous))) {
      c <- continuous[i, ]
      exp <- c$EXP
      strat <- c$STRAT
      outcome <- c$OUTCOME

      # X -> C
      c_is_collider <- if (is.null(causal_relationships)) {
        FALSE
      } else {
        any(causal_relationships$C1 == sub("_EXP$", "", exp) & 
            causal_relationships$C2 == sub("_STRAT$", "", strat)
            )
      }
      # C -> X
      c_is_causal <- if (is.null(causal_relationships)) {
        FALSE
      } else {
        any(causal_relationships$C1 == sub("_STRAT$", "", strat) & 
        causal_relationships$C2 == sub("_EXP$", "", exp)
        )
      }

      dr_by_x <- mr(
        method = "by_x",
        exposure_name = exp,
        outcome_name = outcome,
        stratifying_name = strat,
        quantile_number = quantile_number,
        data = df,
        output_dir = output_dir
      )
      if (dr_by_x$p < (p_threshold / n)) {
        if (c_is_collider) {
          m1_y2_r2 <- c$C_Y2_R2
          if (m1_y2_r2 <= 2 * 10^-4) {
            continuous_results <- rbind(continuous_results, data.table(
              EXP = exp,
              STRAT = strat,
              OUTCOME = outcome,
              SIGNIFICANT = "TRUE",
              C_Y2_R2 = c$C_Y2_R2,
              X_GC_R2 = c$X_GC_R2,
              X_C2_R2 = c$X_C2_R2,
              Y_C2_R2 = c$Y_C2_R2,
              Y_X2_R2 = c$Y_X2_R2,
              C_X2_R2 = c$C_X2_R2,
              DR_BY_X_P = dr_by_x$p,
              DR_BY_X_BETA = dr_by_x$beta,
              DR_BY_X_SE = dr_by_x$se,
              DR_BY_G_P = NA,
              DR_BY_G_BETA = NA,
              DR_BY_G_SE = NA,
              RESID_STRAT_P = NA,
              RESID_STRAT_BETA = NA,
              RESID_STRAT_SE = NA,

              WALD_RATIO = c$wald_ratio,
              WALD_RATIO_95_LL = c$wald_ratio_95_LL,
              WALD_RATIO_95_UL = c$wald_ratio_95_UL,
              X_G_R2 = c$X_G_R2
            ), fill = TRUE)
          } else {
            m2_gc_r2 <- c$X_GC_R2
            if (m2_gc_r2 <= 5 * 10^-4) {
              resid_strat <- mr(
                method = "resid_strat",
                exposure_name = exp,
                outcome_name = outcome,
                stratifying_name = strat,
                quantile_number = quantile_number,
                data = df,
                output_dir = output_dir
              )
              if (resid_strat$p < (p_threshold / n)) {
                continuous_results <- rbind(continuous_results, data.table(
                  EXP = exp,
                  STRAT = strat,
                  OUTCOME = outcome,
                  SIGNIFICANT = "TRUE",
                  C_Y2_R2 = c$C_Y2_R2,
                  X_GC_R2 = c$X_GC_R2,
                  X_C2_R2 = c$X_C2_R2,
                  Y_C2_R2 = c$Y_C2_R2,
                  Y_X2_R2 = c$Y_X2_R2,
                  C_X2_R2 = c$C_X2_R2,
                  DR_BY_X_P = dr_by_x$p,
                  DR_BY_X_BETA = dr_by_x$beta,
                  DR_BY_X_SE = dr_by_x$se,
                  DR_BY_G_P = NA,
                  DR_BY_G_BETA = NA,
                  DR_BY_G_SE = NA,
                  RESID_STRAT_P = resid_strat$p,
                  RESID_STRAT_BETA = resid_strat$beta,
                  RESID_STRAT_SE = resid_strat$se,

                  WALD_RATIO = c$wald_ratio,
                  WALD_RATIO_95_LL = c$wald_ratio_95_LL,
                  WALD_RATIO_95_UL = c$wald_ratio_95_UL,
                  X_G_R2 = c$X_G_R2
                ), fill = TRUE)
              } else {
                continuous_results <- rbind(continuous_results, data.table(
                  EXP = exp,
                  STRAT = strat,
                  OUTCOME = outcome,
                  SIGNIFICANT = "FALSE",
                  C_Y2_R2 = c$C_Y2_R2,
                  X_GC_R2 = c$X_GC_R2,
                  X_C2_R2 = c$X_C2_R2,
                  Y_C2_R2 = c$Y_C2_R2,
                  Y_X2_R2 = c$Y_X2_R2,
                  C_X2_R2 = c$C_X2_R2,
                  DR_BY_X_P = dr_by_x$p,
                  DR_BY_X_BETA = dr_by_x$beta,
                  DR_BY_X_SE = dr_by_x$se,
                  DR_BY_G_P = NA,
                  DR_BY_G_BETA = NA,
                  DR_BY_G_SE = NA,
                  RESID_STRAT_P = resid_strat$p,
                  RESID_STRAT_BETA = resid_strat$beta,
                  RESID_STRAT_SE = resid_strat$se,

                  WALD_RATIO = c$wald_ratio,
                  WALD_RATIO_95_LL = c$wald_ratio_95_LL,
                  WALD_RATIO_95_UL = c$wald_ratio_95_UL,
                  X_G_R2 = c$X_G_R2
                ), fill = TRUE)
              }
            } else {
              continuous_results <- rbind(continuous_results, data.table(
                EXP = exp,
                STRAT = strat,
                OUTCOME = outcome,
                SIGNIFICANT = "INCONCLUSIVE",
                C_Y2_R2 = c$C_Y2_R2,
                X_GC_R2 = c$X_GC_R2,
                X_C2_R2 = c$X_C2_R2,
                Y_C2_R2 = c$Y_C2_R2,
                Y_X2_R2 = c$Y_X2_R2,
                C_X2_R2 = c$C_X2_R2,
                DR_BY_X_P = dr_by_x$p,
                DR_BY_X_BETA = dr_by_x$beta,
                DR_BY_X_SE = dr_by_x$se,
                DR_BY_G_P = NA,
                DR_BY_G_BETA = NA,
                DR_BY_G_SE = NA,
                RESID_STRAT_P = NA,
                RESID_STRAT_BETA = NA,
                RESID_STRAT_SE = NA,

                WALD_RATIO = c$wald_ratio,
                WALD_RATIO_95_LL = c$wald_ratio_95_LL,
                WALD_RATIO_95_UL = c$wald_ratio_95_UL,
                X_G_R2 = c$X_G_R2
              ), fill = TRUE)
            }
          }
        } else if (c_is_causal) {
          m2_gc_r2 <- c$X_GC_R2
          m3_c2_r2 <- c$X_C2_R2
          m4_c2_r2 <- c$Y_C2_R2
          if (m2_gc_r2 <= 2 * 10^-4 && m3_c2_r2 <= 2 * 10^-3 && m4_c2_r2 <= 2 * 10^-4) {
            continuous_results <- rbind(continuous_results, data.table(
              EXP = exp,
              STRAT = strat,
              OUTCOME = outcome,
              SIGNIFICANT = "TRUE",
              C_Y2_R2 = c$C_Y2_R2,
              X_GC_R2 = c$X_GC_R2,
              X_C2_R2 = c$X_C2_R2,
              Y_C2_R2 = c$Y_C2_R2,
              Y_X2_R2 = c$Y_X2_R2,
              C_X2_R2 = c$C_X2_R2,
              DR_BY_X_P = dr_by_x$p,
              DR_BY_X_BETA = dr_by_x$beta,
              DR_BY_X_SE = dr_by_x$se,
              DR_BY_G_P = NA,
              DR_BY_G_BETA = NA,
              DR_BY_G_SE = NA,
              RESID_STRAT_P = NA,
              RESID_STRAT_BETA = NA,
              RESID_STRAT_SE = NA,

              WALD_RATIO = c$wald_ratio,
              WALD_RATIO_95_LL = c$wald_ratio_95_LL,
              WALD_RATIO_95_UL = c$wald_ratio_95_UL,
              X_G_R2 = c$X_G_R2
            ), fill = TRUE)
          } else {
            m5_x2_r2 <- c$Y_X2_R2
            if (m5_x2_r2 <= 5 * 10^-5) {
              dr_by_g <- mr(
                method = "by_g",
                exposure_name = exp,
                outcome_name = outcome,
                stratifying_name = strat,
                quantile_number = quantile_number,
                data = df,
                output_dir = output_dir
              )
              if (dr_by_g$p < (p_threshold / n)) {
                continuous_results <- rbind(continuous_results, data.table(
                  EXP = exp,
                  STRAT = strat,
                  OUTCOME = outcome,
                  SIGNIFICANT = "TRUE",
                  C_Y2_R2 = c$C_Y2_R2,
                  X_GC_R2 = c$X_GC_R2,
                  X_C2_R2 = c$X_C2_R2,
                  Y_C2_R2 = c$Y_C2_R2,
                  Y_X2_R2 = c$Y_X2_R2,
                  C_X2_R2 = c$C_X2_R2,
                  DR_BY_X_P = dr_by_x$p,
                  DR_BY_X_BETA = dr_by_x$beta,
                  DR_BY_X_SE = dr_by_x$se,
                  DR_BY_G_P = dr_by_g$p,
                  DR_BY_G_BETA = dr_by_g$beta,
                  DR_BY_G_SE = dr_by_g$se,
                  RESID_STRAT_P = NA,
                  RESID_STRAT_BETA = NA,
                  RESID_STRAT_SE = NA,

                  WALD_RATIO = c$wald_ratio,
                  WALD_RATIO_95_LL = c$wald_ratio_95_LL,
                  WALD_RATIO_95_UL = c$wald_ratio_95_UL,
                  X_G_R2 = c$X_G_R2
                ), fill = TRUE)
              } else {
                continuous_results <- rbind(continuous_results, data.table(
                  EXP = exp,
                  STRAT = strat,
                  OUTCOME = outcome,
                  SIGNIFICANT = "FALSE",
                  C_Y2_R2 = c$C_Y2_R2,
                  X_GC_R2 = c$X_GC_R2,
                  X_C2_R2 = c$X_C2_R2,
                  Y_C2_R2 = c$Y_C2_R2,
                  Y_X2_R2 = c$Y_X2_R2,
                  C_X2_R2 = c$C_X2_R2,
                  DR_BY_X_P = dr_by_x$p,
                  DR_BY_X_BETA = dr_by_x$beta,
                  DR_BY_X_SE = dr_by_x$se,
                  DR_BY_G_P = dr_by_g$p,
                  DR_BY_G_BETA = dr_by_g$beta,
                  DR_BY_G_SE = dr_by_g$se,
                  RESID_STRAT_P = NA,
                  RESID_STRAT_BETA = NA,
                  RESID_STRAT_SE = NA,

                  WALD_RATIO = c$wald_ratio,
                  WALD_RATIO_95_LL = c$wald_ratio_95_LL,
                  WALD_RATIO_95_UL = c$wald_ratio_95_UL,
                  X_G_R2 = c$X_G_R2
                ), fill = TRUE)
              }
            } else {
              continuous_results <- rbind(continuous_results, data.table(
                EXP = exp,
                STRAT = strat,
                OUTCOME = outcome,
                SIGNIFICANT = "INCONCLUSIVE",
                C_Y2_R2 = c$C_Y2_R2,
                X_GC_R2 = c$X_GC_R2,
                X_C2_R2 = c$X_C2_R2,
                Y_C2_R2 = c$Y_C2_R2,
                Y_X2_R2 = c$Y_X2_R2,
                C_X2_R2 = c$C_X2_R2,
                DR_BY_X_P = dr_by_x$p,
                DR_BY_X_BETA = dr_by_x$beta,
                DR_BY_X_SE = dr_by_x$se,
                DR_BY_G_P = NA,
                DR_BY_G_BETA = NA,
                DR_BY_G_SE = NA,
                RESID_STRAT_P = NA,
                RESID_STRAT_BETA = NA,
                RESID_STRAT_SE = NA,

                WALD_RATIO = c$wald_ratio,
                WALD_RATIO_95_LL = c$wald_ratio_95_LL,
                WALD_RATIO_95_UL = c$wald_ratio_95_UL,
                X_G_R2 = c$X_G_R2
              ), fill = TRUE)
            }
          }
        } else { # unclear causal distinction
          m1_y2_r2 <- c$C_Y2_R2
          m2_gc_r2 <- c$X_GC_R2
          m3_c2_r2 <- c$X_C2_R2
          m4_c2_r2 <- c$Y_C2_R2
          if (m1_y2_r2 > 2 * 10^-4 && m2_gc_r2 <= 5 * 10^-5 && m3_c2_r2 <= 2 * 10^-4 && m4_c2_r2 <= 5 * 10^-3) {
            resid_strat <- mr(
              method = "resid_strat",
              exposure_name = exp,
              outcome_name = outcome,
              stratifying_name = strat,
              quantile_number = quantile_number,
              data = df,
              output_dir = output_dir
            )
            if (resid_strat$p < (p_threshold / n)) {
              continuous_results <- rbind(continuous_results, data.table(
                EXP = exp,
                STRAT = strat,
                OUTCOME = outcome,
                SIGNIFICANT = "TRUE",
                C_Y2_R2 = c$C_Y2_R2,
                X_GC_R2 = c$X_GC_R2,
                X_C2_R2 = c$X_C2_R2,
                Y_C2_R2 = c$Y_C2_R2,
                Y_X2_R2 = c$Y_X2_R2,
                C_X2_R2 = c$C_X2_R2,
                DR_BY_X_P = dr_by_x$p,
                DR_BY_X_BETA = dr_by_x$beta,
                DR_BY_X_SE = dr_by_x$se,
                DR_BY_G_P = NA,
                DR_BY_G_BETA = NA,
                DR_BY_G_SE = NA,
                RESID_STRAT_P = resid_strat$p,
                RESID_STRAT_BETA = resid_strat$beta,
                RESID_STRAT_SE = resid_strat$se,

                WALD_RATIO = c$wald_ratio,
                WALD_RATIO_95_LL = c$wald_ratio_95_LL,
                WALD_RATIO_95_UL = c$wald_ratio_95_UL,
                X_G_R2 = c$X_G_R2
              ), fill = TRUE)
            } else {
              continuous_results <- rbind(continuous_results, data.table(
                EXP = exp,
                STRAT = strat,
                OUTCOME = outcome,
                SIGNIFICANT = "FALSE",
                C_Y2_R2 = c$C_Y2_R2,
                X_GC_R2 = c$X_GC_R2,
                X_C2_R2 = c$X_C2_R2,
                Y_C2_R2 = c$Y_C2_R2,
                Y_X2_R2 = c$Y_X2_R2,
                C_X2_R2 = c$C_X2_R2,
                DR_BY_X_P = dr_by_x$p,
                DR_BY_X_BETA = dr_by_x$beta,
                DR_BY_X_SE = dr_by_x$se,
                DR_BY_G_P = NA,
                DR_BY_G_BETA = NA,
                DR_BY_G_SE = NA,
                RESID_STRAT_P = resid_strat$p,
                RESID_STRAT_BETA = resid_strat$beta,
                RESID_STRAT_SE = resid_strat$se,

                WALD_RATIO = c$wald_ratio,
                WALD_RATIO_95_LL = c$wald_ratio_95_LL,
                WALD_RATIO_95_UL = c$wald_ratio_95_UL,
                X_G_R2 = c$X_G_R2
              ), fill = TRUE)
            }
          } else if (m1_y2_r2 > 2 * 10^-4 && 
                     (m2_gc_r2 > 5 * 10^-5 || 
                      m3_c2_r2 > 2 * 10^-4 || 
                      m4_c2_r2 > 5 * 10^-3)) {
            continuous_results <- rbind(continuous_results, data.table(
              EXP = exp,
              STRAT = strat,
              OUTCOME = outcome,
              SIGNIFICANT = "INCONCLUSIVE",
              C_Y2_R2 = c$C_Y2_R2,
              X_GC_R2 = c$X_GC_R2,
              X_C2_R2 = c$X_C2_R2,
              Y_C2_R2 = c$Y_C2_R2,
              Y_X2_R2 = c$Y_X2_R2,
              C_X2_R2 = c$C_X2_R2,
              DR_BY_X_P = dr_by_x$p,
              DR_BY_X_BETA = dr_by_x$beta,
              DR_BY_X_SE = dr_by_x$se,
              DR_BY_G_P = NA,
              DR_BY_G_BETA = NA,
              DR_BY_G_SE = NA,
              RESID_STRAT_P = NA,
              RESID_STRAT_BETA = NA,
              RESID_STRAT_SE = NA,

              WALD_RATIO = c$wald_ratio,
              WALD_RATIO_95_LL = c$wald_ratio_95_LL,
              WALD_RATIO_95_UL = c$wald_ratio_95_UL,
              X_G_R2 = c$X_G_R2
            ), fill = TRUE)
          } else if (m1_y2_r2 <= 2 * 10^-4 && m2_gc_r2 <= 2 * 10^-4 && m3_c2_r2 <= 2 * 10^-3 && m4_c2_r2 <= 2 * 10^-4) {
            continuous_results <- rbind(continuous_results, data.table(
              EXP = exp,
              STRAT = strat,
              OUTCOME = outcome,
              SIGNIFICANT = "TRUE",
              C_Y2_R2 = c$C_Y2_R2,
              X_GC_R2 = c$X_GC_R2,
              X_C2_R2 = c$X_C2_R2,
              Y_C2_R2 = c$Y_C2_R2,
              Y_X2_R2 = c$Y_X2_R2,
              C_X2_R2 = c$C_X2_R2,
              DR_BY_X_P = dr_by_x$p,
              DR_BY_X_BETA = dr_by_x$beta,
              DR_BY_X_SE = dr_by_x$se,
              DR_BY_G_P = NA,
              DR_BY_G_BETA = NA,
              DR_BY_G_SE = NA,
              RESID_STRAT_P = NA,
              RESID_STRAT_BETA = NA,
              RESID_STRAT_SE = NA,

              WALD_RATIO = c$wald_ratio,
              WALD_RATIO_95_LL = c$wald_ratio_95_LL,
              WALD_RATIO_95_UL = c$wald_ratio_95_UL,
              X_G_R2 = c$X_G_R2
            ), fill = TRUE)
          } else {
            m5_x2_r2 <- c$Y_X2_R2
            m6_x2_r2 <- c$C_X2_R2
            if (m5_x2_r2 <= 2 * 10^-5 && m6_x2_r2 <= 2 * 10^-3) {
              dr_by_g <- mr(
                method = "by_g",
                exposure_name = exp,
                outcome_name = outcome,
                stratifying_name = strat,
                quantile_number = quantile_number,
                data = df,
                output_dir = output_dir
              )
              if (dr_by_g$p < (p_threshold / n)) {
                continuous_results <- rbind(continuous_results, data.table(
                  EXP = exp,
                  STRAT = strat,
                  OUTCOME = outcome,
                  SIGNIFICANT = "TRUE",
                  C_Y2_R2 = c$C_Y2_R2,
                  X_GC_R2 = c$X_GC_R2,
                  X_C2_R2 = c$X_C2_R2,
                  Y_C2_R2 = c$Y_C2_R2,
                  Y_X2_R2 = c$Y_X2_R2,
                  C_X2_R2 = c$C_X2_R2,
                  DR_BY_X_P = dr_by_x$p,
                  DR_BY_X_BETA = dr_by_x$beta,
                  DR_BY_X_SE = dr_by_x$se,
                  DR_BY_G_P = dr_by_g$p,
                  DR_BY_G_BETA = dr_by_g$beta,
                  DR_BY_G_SE = dr_by_g$se,
                  RESID_STRAT_P = NA,
                  RESID_STRAT_BETA = NA,
                  RESID_STRAT_SE = NA,

                  WALD_RATIO = c$wald_ratio,
                  WALD_RATIO_95_LL = c$wald_ratio_95_LL,
                  WALD_RATIO_95_UL = c$wald_ratio_95_UL,
                  X_G_R2 = c$X_G_R2
                ), fill = TRUE)
              } else {
                continuous_results <- rbind(continuous_results, data.table(
                  EXP = exp,
                  STRAT = strat,
                  OUTCOME = outcome,
                  SIGNIFICANT = "FALSE",
                  C_Y2_R2 = c$C_Y2_R2,
                  X_GC_R2 = c$X_GC_R2,
                  X_C2_R2 = c$X_C2_R2,
                  Y_C2_R2 = c$Y_C2_R2,
                  Y_X2_R2 = c$Y_X2_R2,
                  C_X2_R2 = c$C_X2_R2,
                  DR_BY_X_P = dr_by_x$p,
                  DR_BY_X_BETA = dr_by_x$beta,
                  DR_BY_X_SE = dr_by_x$se,
                  DR_BY_G_P = dr_by_g$p,
                  DR_BY_G_BETA = dr_by_g$beta,
                  DR_BY_G_SE = dr_by_g$se,
                  RESID_STRAT_P = NA,
                  RESID_STRAT_BETA = NA,
                  RESID_STRAT_SE = NA,

                  WALD_RATIO = c$wald_ratio,
                  WALD_RATIO_95_LL = c$wald_ratio_95_LL,
                  WALD_RATIO_95_UL = c$wald_ratio_95_UL,
                  X_G_R2 = c$X_G_R2
                ), fill = TRUE)
              }
            } else {
              continuous_results <- rbind(continuous_results, data.table(
                EXP = exp,
                STRAT = strat,
                OUTCOME = outcome,
                SIGNIFICANT = "INCONCLUSIVE",
                C_Y2_R2 = c$C_Y2_R2,
                X_GC_R2 = c$X_GC_R2,
                X_C2_R2 = c$X_C2_R2,
                Y_C2_R2 = c$Y_C2_R2,
                Y_X2_R2 = c$Y_X2_R2,
                C_X2_R2 = c$C_X2_R2,
                DR_BY_X_P = dr_by_x$p,
                DR_BY_X_BETA = dr_by_x$beta,
                DR_BY_X_SE = dr_by_x$se,
                DR_BY_G_P = NA,
                DR_BY_G_BETA = NA,
                DR_BY_G_SE = NA,
                RESID_STRAT_P = NA,
                RESID_STRAT_BETA = NA,
                RESID_STRAT_SE = NA,
                
                WALD_RATIO = c$wald_ratio,
                WALD_RATIO_95_LL = c$wald_ratio_95_LL,
                WALD_RATIO_95_UL = c$wald_ratio_95_UL,
                X_G_R2 = c$X_G_R2
              ), fill = TRUE)
            }
          }
        }
      } else {
        continuous_results <- rbind(continuous_results, data.table(
          EXP = exp,
          STRAT = strat,
          OUTCOME = outcome,
          SIGNIFICANT = "FALSE",
          C_Y2_R2 = c$C_Y2_R2,
          X_GC_R2 = c$X_GC_R2,
          X_C2_R2 = c$X_C2_R2,
          Y_C2_R2 = c$Y_C2_R2,
          Y_X2_R2 = c$Y_X2_R2,
          C_X2_R2 = c$C_X2_R2,
          DR_BY_X_P = dr_by_x$p,
          DR_BY_X_BETA = dr_by_x$beta,
          DR_BY_X_SE = dr_by_x$se,
          DR_BY_G_P = NA,
          DR_BY_G_BETA = NA,
          DR_BY_G_SE = NA,
          RESID_STRAT_P = NA,
          RESID_STRAT_BETA = NA,
          RESID_STRAT_SE = NA,

          WALD_RATIO = c$wald_ratio,
          WALD_RATIO_95_LL = c$wald_ratio_95_LL,
          WALD_RATIO_95_UL = c$wald_ratio_95_UL,
          X_G_R2 = c$X_G_R2
        ), fill = TRUE)
      }
    }
  }

  if (!is.null(output_dir)) {
    if (!dir.exists(output_dir)) {
      dir.create(output_dir)
    }
    if (nrow(binary_results) > 0) {
      fwrite(binary_results, file.path(output_dir, "binary_results.csv"))
    }
    if (nrow(continuous_results) > 0) {
      fwrite(continuous_results, file.path(output_dir, "continuous_results.csv"))
    }
  }

  return(list(binary_results = binary_results, continuous_results = continuous_results))
}

#' @title All methods of the algorithm
#' @param quantile_number The number of quantiles or strata to use.
#' @param input_dir The directory where the data is stored, defaults to the current working directory.
#' @param output_dir The directory where the output should be stored, defaults to NULL meaning no output is saved.
#' @param data The name of the data file or an input data.table, defaults to "data.csv". Column names need to be in a specific format; see README page for details.
#' @param diagnostically_related The file name of the diagnostically related variables, or a data frame with two columns. No column headers required; see README page for details. Defaults to NULL meaning no diagnostically related variables are provided. Parameter used for preprocessing step.
#' @param binary The file name of the binary combinations, or a data frame with three columns (EXP, STRAT, CASE_STATUS). Column headers required. Defaults to NULL meaning the combinations will be automatically generated. Parameter used for preprocessing step.
#' @param continuous The file name of the continuous combinations, or a data frame with three columns (EXP, STRAT, OUTCOME). Column headers required. Defaults to NULL meaning the combinations will be automatically generated. Parameter used for preprocessing step.
#' @param n The number of combinations in the dataset (used for calculating the Bonferroni threshold), defaults to 1 meaning that no Bonferroni correction will be automatically applied.
#' @param causal_relationships The file name of the causal relationships, or a data frame with two columns defining a causal relationship where C1 is known to cause C2 (C1 -> C2). Headers are required (named C1 and C2). This defines the causal and collider relationships between C and X in the actual algorithm. If NULL, or there is no entry for C1 -> C2, then the causal relationship is unknown. Defaults to NULL meaning no causal relationships are provided.
#' @param combinations Skip preprocessing and pass in data frames (combinations$binary and combinations$continuous) returned from the preprocess function. Defaults to NULL meaning it will preprocess the provided data.
#' @param p_threshold The threshold for statistical significance, defaults to 0.05.
#' @return A list containing two data tables, one for binary results and one for continuous results. See README page for details.
#' @export
all_methods <- function(
    quantile_number,
    input_dir = getwd(),
    output_dir = NULL,
    data = "data.csv",
    diagnostically_related = NULL,
    binary = NULL,
    continuous = NULL,
    n = NULL,
    causal_relationships = NULL,
    combinations = NULL,
    p_threshold = 0.05
  ) {
  biased <- "POSSIBLE_BIAS"
  unbiased <- "NO_BIAS_DETECTED"

  quantile_number <- as.integer(quantile_number)
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
  } else if (!is.data.frame(data)) {
    stop("data must be a file name or a data frame.")
  } else {
    df <- data
  }

  if (is.null(combinations)) {
    combinations <- preprocess(input_dir = input_dir, data = df, diagnostically_related = diagnostically_related, binary = binary, continuous = continuous, output_dir = output_dir) # nolint
    binary <- combinations$binary
    continuous <- combinations$continuous
  } else {
    binary <- combinations$binary
    continuous <- combinations$continuous
    message("Preprocessing steps skipped, using combinations$binary and/or combinations$continuous directly from R environment.")
  }

  binary_n <- if (is.null(binary)) {
    0
  } else {
    nrow(binary)
  }
  continuous_n <- if (is.null(continuous)) {
    0
  } else {
    nrow(continuous)
  }
  n <- if (is.null(n)) {
    # No need for Bonferroni correction
    message("No Bonferroni correction applied.")
    1
  } else {
    message(paste("Total number of combinations used for Bonferroni correction =", n))
    n
  }
  
  # Load in the causal_relationships, if provided
  if (!is.null(causal_relationships)) {
    if (is.character(causal_relationships)) {
      if (!file.exists(file.path(input_dir, causal_relationships))) {
        stop("causal relationships file does not exist.")
      }
      causal_relationships <- fread(
        file.path(input_dir, causal_relationships),
        stringsAsFactors = FALSE,
        header = TRUE,
        data.table = TRUE,
      )
    } else {
      if (!is.data.table(causal_relationships)) {
        stop("causal_relationships must be a data.table")
      }
    }
    message("Loaded causal relationships.")
  } else {
    causal_relationships <- NULL
    message("No causal relationships provided.")
  }

  binary_results <- data.table(
    EXP = character(),
    STRAT = character(),
    CASE_STATUS = character(),
    DR_BY_X_P = numeric(),
    DR_BY_X_BETA = numeric(),
    DR_BY_X_SE = numeric(),
    DR_BY_X_SIGNIFICANT = character(),
    DR_BY_G_P = numeric(),
    DR_BY_G_BETA = numeric(),
    DR_BY_G_SE = numeric(),
    DR_BY_G_SIGNIFICANT = character(),
    RESID_STRAT_P = numeric(),
    RESID_STRAT_BETA = numeric(),
    RESID_STRAT_SE = numeric(),
    RESID_STRAT_SIGNIFICANT = character(),
    RAW_STRAT_P = numeric(),
    RAW_STRAT_BETA = numeric(),
    RAW_STRAT_SE = numeric(),
    RAW_STRAT_SIGNIFICANT = character(),
    WALD_RATIO_OR = numeric(),
    WALD_RATIO_OR_95_LL = numeric(),
    WALD_RATIO_OR_95_UL = numeric(),
    X_G_R2 = numeric()
  )
  continuous_results <- data.table(
    EXP = character(),
    STRAT = character(),
    OUTCOME = character(),
    DR_BY_X_P = numeric(),
    DR_BY_X_BETA = numeric(),
    DR_BY_X_SE = numeric(),
    DR_BY_X_SIGNIFICANT = character(),
    DR_BY_G_P = numeric(),
    DR_BY_G_BETA = numeric(),
    DR_BY_G_SE = numeric(),
    DR_BY_G_SIGNIFICANT = character(),
    RESID_STRAT_P = numeric(),
    RESID_STRAT_BETA = numeric(),
    RESID_STRAT_SE = numeric(),
    RESID_STRAT_SIGNIFICANT = character(),
    RESID_WALD_P = numeric(),
    RESID_WALD_BETA = numeric(),
    RESID_WALD_SE = numeric(),
    RESID_WALD_SIGNIFICANT = character(),
    RAW_STRAT_P = numeric(),
    RAW_STRAT_BETA = numeric(),
    RAW_STRAT_SE = numeric(),
    RAW_STRAT_SIGNIFICANT = character(),
    WALD_RATIO = numeric(),
    WALD_RATIO_95_LL = numeric(),
    WALD_RATIO_95_UL = numeric(),
    X_G_R2 = numeric()
  )

  is_significant <- function(p) {
    if (p < (p_threshold / n)) {
      "TRUE"
    } else {
      "FALSE"
    }
  }

  or_biased <- function(or, test) {
    return(test > or || test < (1 / or))
  }

  if (binary_n > 0) {
    for (i in seq_len(nrow(binary))) {
      c <- binary[i, ]
      exp <- c$EXP
      strat <- c$STRAT
      outcome <- c$CASE_STATUS

      # X -> C
      c_is_collider <- if (is.null(causal_relationships)) {
        FALSE
      } else {
        any(causal_relationships$C1 == sub("_EXP$", "", exp) & 
            causal_relationships$C2 == sub("_STRAT$", "", strat)
            )
      }
      # C -> X
      c_is_causal <- if (is.null(causal_relationships)) {
        FALSE
      } else {
        any(causal_relationships$C1 == sub("_STRAT$", "", strat) & 
            causal_relationships$C2 == sub("_EXP$", "", exp)
            )
      }

      dr_by_x <- mr(
        method = "by_x",
        exposure_name = exp,
        outcome_name = outcome,
        stratifying_name = strat,
        quantile_number = quantile_number,
        data = df,
        output_dir = output_dir
      )
      dr_by_g <- mr(
        method = "by_g",
        exposure_name = exp,
        outcome_name = outcome,
        stratifying_name = strat,
        quantile_number = quantile_number,
        data = df,
        output_dir = output_dir
      )
      resid_strat <- mr(
        method = "resid_strat",
        exposure_name = exp,
        outcome_name = outcome,
        stratifying_name = strat,
        quantile_number = quantile_number,
        data = df,
        output_dir = output_dir
      )
      raw_strat <- mr(
        method = "raw_strat",
        exposure_name = exp,
        outcome_name = outcome,
        stratifying_name = strat,
        quantile_number = quantile_number,
        data = df,
        output_dir = output_dir
      )

      binary_results <- rbind(binary_results, data.table(
        EXP = exp,
        STRAT = strat,
        CASE_STATUS = outcome,
        DR_BY_X_P = dr_by_x$p,
        DR_BY_X_BETA = dr_by_x$beta,
        DR_BY_X_SE = dr_by_x$se,
        DR_BY_X_SIGNIFICANT = is_significant(dr_by_x$p),
        DR_BY_X_BIAS = if (c_is_collider) {
          unbiased
        } else { # if (c_is_causal or unknown direction)
          if (or_biased(1.05, c$Y_C2_OR) || c$X_C2_R2 > 2 * 10^-2) {
            biased
          } else {
            unbiased
          }
        },
        DR_BY_G_P = dr_by_g$p,
        DR_BY_G_BETA = dr_by_g$beta,
        DR_BY_G_SE = dr_by_g$se,
        DR_BY_G_SIGNIFICANT = is_significant(dr_by_g$p),
        DR_BY_G_BIAS = if (or_biased(1.02, c$Y_X2_OR)) {
          biased
        } else {
          unbiased
        },
        RESID_STRAT_P = resid_strat$p,
        RESID_STRAT_BETA = resid_strat$beta,
        RESID_STRAT_SE = resid_strat$se,
        RESID_STRAT_SIGNIFICANT = is_significant(resid_strat$p),
        RESID_STRAT_BIAS = if (c_is_collider) {
          unbiased
        } else {
          if (or_biased(1.02, c$Y_C2_OR) || (c$X_C2_R2 > (2 * 10^-3)) || (c$X_GC_R2 > (5 * 10^-4))) {
            biased
          } else {
            unbiased
          }
        },
        RAW_STRAT_P = raw_strat$p,
        RAW_STRAT_BETA = raw_strat$beta,
        RAW_STRAT_SE = raw_strat$se,
        RAW_STRAT_SIGNIFICANT = is_significant(raw_strat$p),
        RAW_STRAT_BIAS = if (or_biased(1.02, c$Y_X2_OR)) {
          biased
        } else {
          unbiased
        },
        WALD_RATIO_OR = c$wald_ratio_OR,
        WALD_RATIO_OR_95_LL = c$wald_ratio_OR_95_LL,
        WALD_RATIO_OR_95_UL = c$wald_ratio_OR_95_UL,
        X_G_R2 = c$X_G_R2
      ), fill = TRUE)
    }
  }

  if (continuous_n > 0) {
    for (i in seq_len(nrow(continuous))) {
      c <- continuous[i, ]
      exp <- c$EXP
      strat <- c$STRAT
      outcome <- c$OUTCOME

      # X -> C
      c_is_collider <- if (is.null(causal_relationships)) {
        FALSE
      } else {
        any(causal_relationships$C1 == sub("_EXP$", "", exp) & 
            causal_relationships$C2 == sub("_STRAT$", "", strat)
            )
      }
      # C -> X
      c_is_causal <- if (is.null(causal_relationships)) {
        FALSE
      } else {
        any(causal_relationships$C1 == sub("_STRAT$", "", strat) & 
            causal_relationships$C2 == sub("_EXP$", "", exp)
            )
      }

      dr_by_x <- mr(
        method = "by_x",
        exposure_name = exp,
        outcome_name = outcome,
        stratifying_name = strat,
        quantile_number = quantile_number,
        data = df,
        output_dir = output_dir
      )
      dr_by_g <- mr(
        method = "by_g",
        exposure_name = exp,
        outcome_name = outcome,
        stratifying_name = strat,
        quantile_number = quantile_number,
        data = df,
        output_dir = output_dir
      )
      resid_strat <- mr(
        method = "resid_strat",
        exposure_name = exp,
        outcome_name = outcome,
        stratifying_name = strat,
        quantile_number = quantile_number,
        data = df,
        output_dir = output_dir
      )
      resid_wald <- mr(
        method = "resid_wald",
        exposure_name = exp,
        outcome_name = outcome,
        stratifying_name = strat,
        quantile_number = quantile_number,
        data = df,
        output_dir = output_dir
      )
      raw_strat <- mr(
        method = "raw_strat",
        exposure_name = exp,
        outcome_name = outcome,
        stratifying_name = strat,
        quantile_number = quantile_number,
        data = df,
        output_dir = output_dir
      )

      continuous_results <- rbind(continuous_results, data.table(
        EXP = exp,
        STRAT = strat,
        OUTCOME = outcome,
        DR_BY_X_P = dr_by_x$p,
        DR_BY_X_BETA = dr_by_x$beta,
        DR_BY_X_SE = dr_by_x$se,
        DR_BY_X_SIGNIFICANT = is_significant(dr_by_x$p),
        DR_BY_X_BIAS = if (c_is_collider) {
          if (c$C_Y2_R2 > 2 * 10^-4) {
            biased
          } else {
            unbiased
          }
        } else if (c_is_causal) {
          if (c$X_GC_R2 > 2 * 10^-4 || c$X_C2_R2 > 2 * 10^-3 || c$Y_C2_R2 > 2 * 10^-4) {
            biased
          } else {
            unbiased
          }
        } else {
          if (c$C_Y2_R2 > 2 * 10^-4 || c$X_GC_R2 > 2 * 10^-4 || c$X_C2_R2 > 2 * 10^-3 || c$Y_C2_R2 > 2 * 10^-4) {
            biased
          } else {
            unbiased
          }
        },
        DR_BY_G_P = dr_by_g$p,
        DR_BY_G_BETA = dr_by_g$beta,
        DR_BY_G_SE = dr_by_g$se,
        DR_BY_G_SIGNIFICANT = is_significant(dr_by_g$p),
        DR_BY_G_BIAS = if (c_is_collider) {
          if (c$Y_X2_R2 > 2 * 10^-5 || c$C_X2_R2 > 2 * 10^-3 || c$C_Y2_R2 > 2 * 10^-4) {
            biased
          } else {
            unbiased
          }
        } else if (c_is_causal) {
          if (c$Y_X2_R2 > 5 * 10^-5) {
            biased
          } else {
            unbiased
          }
        } else {
          if (c$Y_X2_R2 > 2 * 10^-5 || c$C_X2_R2 > 2 * 10^-3 || c$C_Y2_R2 > 2 * 10^-4) {
            biased
          } else {
            unbiased
          }
        },
        RESID_STRAT_P = resid_strat$p,
        RESID_STRAT_BETA = resid_strat$beta,
        RESID_STRAT_SE = resid_strat$se,
        RESID_STRAT_SIGNIFICANT = is_significant(resid_strat$p),
        RESID_STRAT_BIAS = if (c_is_collider) {
          if (c$X_GC_R2 > 5 * 10^-4) {
            biased
          } else {
            unbiased
          }
        } else if (c_is_causal) {
          if (c$X_GC_R2 > 5 * 10^-5 || c$X_C2_R2 > 2 * 10^-4 || c$Y_C2_R2 > 5 * 10^-3) {
            biased
          } else {
            unbiased
          }
        } else {
          if (c$X_GC_R2 > 5 * 10^-5 || c$X_C2_R2 > 2 * 10^-4 || c$Y_C2_R2 > 5 * 10^-3) {
            biased
          } else {
            unbiased
          }
        },
        RESID_WALD_P = resid_wald$p,
        RESID_WALD_BETA = resid_wald$beta,
        RESID_WALD_SE = resid_wald$se,
        RESID_WALD_SIGNIFICANT = is_significant(resid_wald$p),
        RESID_WALD_BIAS = if (c_is_collider) {
          if (c$Y_X2_R2 > 5 * 10^-5 || c$C_X2_R2 > 2 * 10^-4 || c$C_Y2_R2 > 5 * 10^-4) {
            biased
          } else {
            unbiased
          }
        } else if (c_is_causal) {
          if (c$Y_X2_R2 > 1 * 10^-4) {
            biased
          } else {
            unbiased
          }
        } else {
          if (c$Y_X2_R2 > 5 * 10^-5 || c$C_X2_R2 > 2 * 10^-4 || c$C_Y2_R2 > 5 * 10^-4) {
            biased
          } else {
            unbiased
          }
        },
        RAW_STRAT_P = raw_strat$p,
        RAW_STRAT_BETA = raw_strat$beta,
        RAW_STRAT_SE = raw_strat$se,
        RAW_STRAT_SIGNIFICANT = is_significant(raw_strat$p),
        RAW_STRAT_BIAS = if (c_is_collider) {
          if (c$Y_X2_R2 > 5 * 10^-5 || c$C_X2_R2 > 2 * 10^-4 || c$C_Y2_R2 > 1 * 10^-3) {
            biased
          } else {
            unbiased
          }
        } else if (c_is_causal) {
          if (c$Y_X2_R2 > 1 * 10^-4) {
            biased
          } else {
            unbiased
          }
        } else {
          if (c$Y_X2_R2 > 5 * 10^-5 || c$C_X2_R2 > 2 * 10^-4 || c$C_Y2_R2 > 1 * 10^-3) {
            biased
          } else {
            unbiased
          }
        },
        WALD_RATIO = c$wald_ratio,
        WALD_RATIO_95_LL = c$wald_ratio_95_LL,
        WALD_RATIO_95_UL = c$wald_ratio_95_UL,
        X_G_R2 = c$X_G_R2
      ), fill = TRUE)
    }
  }

  if (!is.null(output_dir)) {
    if (!dir.exists(output_dir)) {
      dir.create(output_dir)
    }
    if (nrow(binary_results) > 0) {
      fwrite(binary_results, file.path(output_dir, "binary_results_all.csv"))
    }
    if (nrow(continuous_results) > 0) {
      fwrite(continuous_results, file.path(output_dir, "continuous_results_all.csv"))
    }
  }
  return(list(binary_results = binary_results, continuous_results = continuous_results))
}
