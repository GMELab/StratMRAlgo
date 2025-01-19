library(data.table)
library(dplyr)

#' @title Generate sample data
#' @param output_dir The directory where the sample data files will be stored
#' @return
#' Writes the following files to the specified output directory:
#' \itemize{
#' \item data.csv: A data frame with 4 exposures, 4 PRS, 4 stratifying variables, 1 continuous outcome,
#' and 1 binary outcome for 100K individuals. This setup results in 32 initial combinations.
#' \item diagnostically_related.txt: A tab-delimited text file containing diagnostically-related variables
#' for any exposure or stratifying variable associated with the outcome.
#' \item continuous_comb_interest.txt: A tab-delimited text file with columns EXP, STRAT, and OUTCOME containing
#' combinations of interest with continuous outcomes.
#' \item binary_comb_interest.txt: A tab-delimited text file with columns EXP, STRAT, and CASE_STATUS containing
#' combinations of interest with binary outcomes.
#' \item causal_relationships.txt: A tab-delimited text file with columns C1 and C2, where C1 is known to cause
#' C2.
#' }
#' @export
generate_sample_data <- function(output_dir) {
  #######################################################################################
  # Step 1: Generate data
  #######################################################################################

  set.seed(120)
  nb.ind <- 100000 # Number of participants

  # Generate genetic instrument for each exposure (n = 4)
  X1_PRS <- scale(rnorm(nb.ind)) # PRS for X1
  X2_PRS <- scale(rnorm(nb.ind)) # PRS for X2
  X3_PRS <- scale(rnorm(nb.ind)) # PRS for X3
  X4_PRS <- scale(rnorm(nb.ind)) # PRS for X4

  # Generate exposures (n = 4)
  EX1 <- scale(rnorm(nb.ind)) # the independent term for X1
  EX2 <- scale(rnorm(nb.ind)) # the independent term for X2
  EX3 <- scale(rnorm(nb.ind)) # the independent term for X3
  EX4 <- scale(rnorm(nb.ind)) # the independent term for X4

  X_G_r2 <- 0.05 # variance of X explained by G
  EX_variance <- 1 - X_G_r2 # variance of X explained by X_variance term

  X1_EXP <- scale((X_G_r2^0.5) * X1_PRS + (EX_variance^0.5) * EX1) # exposure X1
  X2_EXP <- scale((X_G_r2^0.5) * X2_PRS + (EX_variance^0.5) * EX2) # exposure X2
  X3_EXP <- scale((X_G_r2^0.5) * X3_PRS + (EX_variance^0.5) * EX3) # exposure X3
  X4_EXP <- scale((X_G_r2^0.5) * X4_PRS + (EX_variance^0.5) * EX4) # exposure X4

  # Generate stratifying variables (n = 4)
  X1_STRAT <- X1_EXP # first stratifying variable equals the exposure X1 (therefore, will be filtered out in some instances)

  C1_X1_r2 <- 0.20 # variance of C1 explained by X1 (very high, will lead to this combination being filtered out since r2 of 0.1 is the maximum)
  EC1 <- scale(rnorm(nb.ind)) # the independent term for C1
  EC1_variance <- 1 - C1_X1_r2 # variance of C1 explained by the independent term
  C1_STRAT <- scale((C1_X1_r2^0.5) * X1_EXP + (EC1_variance^0.5) * EC1) # second stratifying variable C1 is explained by a combination of X1 and an independent variable

  C2_X1_r2 <- 0.03 # variance of C2 explained by X1
  EC2 <- scale(rnorm(nb.ind)) # the independent term for C2
  EC2_variance <- 1 - C2_X1_r2 # variance of C2 explained by the independent term

  C2_STRAT <- scale((C2_X1_r2^0.5) * X1_EXP + (EC2_variance^0.5) * EC2) # third stratifying variable C2 is explained by a combination of X1 and an independent variable

  C3_STRAT <- scale(rnorm(nb.ind)) # fourth stratifying variable C3

  # Generate continuous outcomes (n = 2)
  Y1_X1_r2 <- 0.08 # variance of Y1 explained by X1
  Y1_X2_r2 <- 0.07 # variance of Y1 explained by X2
  Y1_X3_r2 <- 0 # variance of Y1 explained by X3
  Y1_X4_r2 <- 0.06 # variance of Y1 explained by X4
  Y1_X1C2_r2 <- 0.08 # variance of Y1 explained by X1*C2 (effect modification or interaction term)
  EY1 <- scale(rnorm(nb.ind)) # the independent term for Y1
  EY1_variance <- 1 - Y1_X1_r2 - Y1_X2_r2 - Y1_X3_r2 - Y1_X4_r2 - Y1_X1C2_r2 # variance of Y1 explained by the independent term

  Y1_OUTCOME <- scale((Y1_X1_r2^0.5) * X1_EXP + (Y1_X2_r2^0.5) * X2_EXP + (Y1_X3_r2^0.5) * X3_EXP + (Y1_X4_r2^0.5) * X4_EXP
    + (Y1_X1C2_r2^0.5) * scale(X1_EXP * C2_STRAT) + (EY1_variance^0.5) * EY1) # Y1 is explained by X1, X2, X4, X1*C2, and an independent variable

  # Generate binary outcome (n = 1)
  Y2_X1_OR <- 1.5 # effect of X1 on Y2 (odds ratio)
  Y2_X1_logit <- log(Y2_X1_OR) # convert to log(odds)

  Y2_X2_OR <- 1.3 # effect of X2 on Y2 (odds ratio)
  Y2_X2_logit <- log(Y2_X2_OR) # convert to log(odds)

  Y2_X3_OR <- 1 # effect of X3 on Y2 (odds ratio)
  Y2_X3_logit <- log(Y2_X3_OR) # convert to log(odds)

  Y2_X4_OR <- 1.4 # effect of X4 on Y2 (odds ratio)
  Y2_X4_logit <- log(Y2_X4_OR) # convert to log(odds)

  Y2_X1C2_OR <- 1.4 # effect of X1*C2 on Y (odds ratio)
  Y2_X1C2_logit <- log(Y2_X1C2_OR) # convert to log(odds)

  prev_logit <- -log(24) # log(odds) representing prevalence
  prevalence <- 1 / (1 + exp(-prev_logit)) # set disease prevalence to 4%

  Y2_fitted_prob <- 1 / (1 + exp(-(Y2_X1_logit * X1_EXP + Y2_X2_logit * X2_EXP + Y2_X3_logit * X3_EXP + Y2_X4_logit * X4_EXP
    + Y2_X1C2_logit * scale(X1_EXP * C2_STRAT) + prev_logit))) # Y2 is explained by X1, X2, X4, X1*C2, and an independent variable

  Y2_CASE_STATUS <- rbinom(nb.ind, 1, prob = Y2_fitted_prob) # define Y2 using the fitted probability

  # Generate covariates (n = 3)
  V1_COVAR <- scale(rnorm(nb.ind))
  V2_COVAR <- scale(rnorm(nb.ind))
  V3_COVAR <- scale(rnorm(nb.ind))
  covar_df <- data.frame(V1_COVAR = V1_COVAR, V2_COVAR = V2_COVAR, V3_COVAR = V3_COVAR)

  # Merge all of these data into one data frame as the input for the algorithm

  df <- data.frame(
    X1_PRS = X1_PRS, X2_PRS = X2_PRS, X3_PRS = X3_PRS, X4_PRS = X4_PRS,
    X1_EXP = X1_EXP, X2_EXP = X2_EXP, X3_EXP = X3_EXP, X4_EXP = X4_EXP,
    X1_STRAT = X1_STRAT, C1_STRAT = C1_STRAT, C2_STRAT = C2_STRAT, C3_STRAT = C3_STRAT,
    Y1_OUTCOME = Y1_OUTCOME, Y2_CASE_STATUS = Y2_CASE_STATUS, V1_COVAR = V1_COVAR,
    V2_COVAR = V2_COVAR, V3_COVAR = V3_COVAR
  )

  df$eid <- seq(1, nrow(df)) # create participant ID column
  df <- df %>% select(eid, everything())

  if (!dir.exists(output_dir)) {
    dir.create(output_dir)
  }

  # Write file
  fwrite(df, file = file.path(output_dir, "data.csv"), sep = ",", quote = FALSE, na = "NA")

  #######################################################################################
  # Step 2: Generate list of diagnostically-related variables
  #######################################################################################

  # Say a cut-off value for C3_STRAT determines whether you have disease Y2_CASE_STATUS
  # Applies to binary outcomes only

  # Create df
  txt <- data.frame("C3", "Y2", stringsAsFactors = FALSE) # must be in the order of continuous variable and outcome
  colnames(txt) <- NULL

  # Write file
  fwrite(txt, file = file.path(output_dir, "diagnostically_related.txt"), sep = "\t", quote = FALSE, na = "NA")

  #######################################################################################
  # Step 3: Generate list of input combinations of interest
  #######################################################################################

  # Say you are specifically interested in studying a few input combinations of interest
  # You can specify this such that only those combinations go through the preprocesing steps
  # To do this, you will need a dataframe with columns EXP, STRAT, and OUTCOME containing combinations of interest with continuous outcomes
  # You will also need a dataframe with columns EXP, STRAT, and CASE_STATUS containing combinations of interest with binary outcomes

  # Create the continuous outcomes combinations data frame
  continuous_combinations <- data.frame(
    EXP = c("X1_EXP", "X1_EXP", "X3_EXP", "X4_EXP"),
    STRAT = c("C2_STRAT", "C3_STRAT", "C3_STRAT", "X1_STRAT"),
    OUTCOME = c("Y1_OUTCOME", "Y1_OUTCOME", "Y1_OUTCOME", "Y1_OUTCOME"),
    stringsAsFactors = FALSE
  )

  # Create the binary outcomes combinations data frame
  binary_combinations <- data.frame(
    EXP = c("X1_EXP", "X3_EXP", "X4_EXP", "X4_EXP"),
    STRAT = c("C2_STRAT", "C3_STRAT", "X1_STRAT", "C1_STRAT"),
    CASE_STATUS = c("Y2_CASE_STATUS", "Y2_CASE_STATUS", "Y2_CASE_STATUS", "Y2_CASE_STATUS"),
    stringsAsFactors = FALSE
  )

  # Write files
  fwrite(continuous_combinations, file = file.path(output_dir, "continuous_comb_interest.txt"), sep = "\t", quote = FALSE, na = "NA")
  fwrite(binary_combinations, file = file.path(output_dir, "binary_comb_interest.txt"), sep = "\t", quote = FALSE, na = "NA")

  #######################################################################################
  # Step 4: Specify known directionality between X and C
  #######################################################################################

  # Say you know the directionality between some of the X and C variables
  # You can specify this for both the algo() and all_methods() functions
  # To do this, you will need a dataframe with two columns defining a causal relationship where C1 is known to cause C2 (C1 -> C2)

  # Create the causal relationships data frame
  causal_relationships <- data.frame(
    C1 = c("X1", "X1"),
    C2 = c("C1", "C2"),
    stringsAsFactors = FALSE
  )

  # Write file
  fwrite(causal_relationships, file = file.path(output_dir, "causal_relationships.txt"), sep = "\t", quote = FALSE, na = "NA")

  #######################################################################################
  # Note: Starting combinations and their expected effects
  #######################################################################################

  # X1_EXP X1_STRAT Y1_OUTCOME - Filtered out as EXP and STRAT are identical
  # X1_EXP C1_STRAT Y1_OUTCOME - Filtered out as X1 and C1 are strongly correlated
  # X1_EXP C2_STRAT Y1_OUTCOME - Effect modification explicitly modeled for this combination
  # X1_EXP C3_STRAT Y1_OUTCOME
  # X2_EXP X1_STRAT Y1_OUTCOME
  # X2_EXP C1_STRAT Y1_OUTCOME
  # X2_EXP C2_STRAT Y1_OUTCOME
  # X2_EXP C3_STRAT Y1_OUTCOME
  # X3_EXP X1_STRAT Y1_OUTCOME - Filtered out as X3 is not causal of Y1 overall
  # X3_EXP C1_STRAT Y1_OUTCOME - Filtered out as X3 is not causal of Y1 overall
  # X3_EXP C2_STRAT Y1_OUTCOME - Filtered out as X3 is not causal of Y1 overall
  # X3_EXP C3_STRAT Y1_OUTCOME - Filtered out as X3 is not causal of Y1 overall
  # X4_EXP X1_STRAT Y1_OUTCOME
  # X4_EXP C1_STRAT Y1_OUTCOME
  # X4_EXP C2_STRAT Y1_OUTCOME
  # X4_EXP C3_STRAT Y1_OUTCOME

  # X1_EXP X1_STRAT Y2_CASE_STATUS - Filtered out as EXP and STRAT are identical
  # X1_EXP C1_STRAT Y2_CASE_STATUS - Filtered out as X1 and C1 are strongly correlated
  # X1_EXP C2_STRAT Y2_CASE_STATUS - Effect modification explicitly modeled for this combination
  # X1_EXP C3_STRAT Y2_CASE_STATUS - Filtered out as STRAT and CASE STATUS are diagnostically related
  # X2_EXP X1_STRAT Y2_CASE_STATUS
  # X2_EXP C1_STRAT Y2_CASE_STATUS
  # X2_EXP C2_STRAT Y2_CASE_STATUS
  # X2_EXP C3_STRAT Y2_CASE_STATUS - Filtered out as STRAT and CASE STATUS are diagnostically related
  # X3_EXP X1_STRAT Y2_CASE_STATUS - Filtered out as X3 is not causal of Y2 overall
  # X3_EXP C1_STRAT Y2_CASE_STATUS - Filtered out as X3 is not causal of Y2 overall
  # X3_EXP C2_STRAT Y2_CASE_STATUS - Filtered out as X3 is not causal of Y2 overall
  # X3_EXP C3_STRAT Y2_CASE_STATUS - Filtered out as X3 is not causal of Y2 overall + STRAT and CASE STATUS are diagnostically related
  # X4_EXP X1_STRAT Y2_CASE_STATUS
  # X4_EXP C1_STRAT Y2_CASE_STATUS
  # X4_EXP C2_STRAT Y2_CASE_STATUS
  # X4_EXP C3_STRAT Y2_CASE_STATUS - Filtered out as STRAT and CASE STATUS are diagnostically related

  #######################################################################################
}
