# StratMRAlgo

StratMRAlgo is a stratified Mendelian randomization (MR) algorithm designed to identify effect modifiers of causal relationships, accomodating both binary and continuous outcomes. The algorithm emerged as a heuristic solution following the evaluation of five stratified MR methods across a total of 63 simulation scenarios, which considered forms of collider bias, confounding, mediation, and nonlinear relationships between variables. The derivation of the algorithm and its scope of use can be found in detail in our main manuscript.

This R package automates the algorithm and the five MR methods tested using four functions: <code>preprocess()</code>, <code>mr()</code>, <code>algo()</code>, and <code>all_methods()</code>.

### Table of Contents
- [Installation](#1)
- [Part 1: Data Pre-processing](#2)
- [Part 2: Algorithm](#3)
- [Part 3: All Stratified MR Methods](#4)
- [Part 4: Individual Stratified MR Methods](#5)
- [License](#6)
- [Contact Information](#7)
- [Citation](#8)

## Installation <a name="1"></a>

StratMRAlgo can be directly installed in R using <code>devtools::install_github</code>:

``` r
devtools::install_github("GMELab/StratMRAlgo")
```

If issues are encountered with <code>devtools::install_github</code>, the repository can be downloaded locally using the <code>git clone</code> command as follows:

``` sh
# Enter a directory to download StratMRAlgo
cd directory

# Initialize a local git repository
git init

# Execute the git clone command
git clone https://github.com/GMELab/StratMRAlgo.git
```

And then loaded in <code>R</code> using <code>devtools::load_all</code>:

``` r
setwd("/directory/R")
devtools::load_all(".")
```

This package requires the dependencies <code>data.table</code>, <code>dplyr</code>, and <code>metafor</code>.

## Workflow

### Part 1: Data Pre-processing <a name="2"></a>

The <code>preprocess()</code> function automates the pre-processing steps required for the algorithm. The function requires an input data file with columns representing all exposures, their respective polygenic risk scores (PRS), outcomes, stratifying variables, and covariates which you plan to include in your analyses. For the package to recognize the variable type properly, the variable names must be specified as follows:

| Variable Type | Required Suffix | Sample Column Name(s) | Supported Data Type(s) |
| ------ | ------ | ------ | ------ |
|  |  |  |
| Exposure | <code>_EXP</code> | <code>BMI_EXP</code> |  Continuous |
| Exposure PRS | <code>_PRS</code> | <code>BMI_PRS</code> |  Continuous |
| Stratifying Variable | <code>_STRAT</code> | <code>LDL_STRAT</code> | Continuous |
| Outcome | <code>_OUTCOME</code> if continuous or <code>_CASE_STATUS</code> if binary | <code>HBA1C_OUTCOME</code> <code>T2DM_CASE_STATUS</code> <code>AGE_OUTCOME</code> <code>SEX_CASE_STATUS</code> | Binary and continuous |
| Covariate | <code>_COVAR</code> | <code>AGE_COVAR</code> <code>SEX_COVAR</code> <code>PC1_COVAR</code> | Binary and continuous |

If you wish to consider a variable as more than one variable type, multiple copies of the column should be present with the appropriate suffices added to their names. For example, if you wish to consider BMI as both an exposure and stratifying variable, a duplicate set of columns labeled BMI_EXP and BMI_STRAT should be included in the input file. Additionally, age and sex columns (used either as covariates or negative control outcomes) should be specified in all capital letters <code>AGE</code> and <code>SEX</code> prior to their appropriate suffix.

All possible exposure, stratifying variable, and outcome combinations will be automatically generated using the main input data file before being filtered by the <code>preprocess()</code> function. Alternatively, specific combinations to process can also be added as optional input files using the <code>continuous</code> and <code>binary</code> parameters.

Another optional input file, a .txt file with diagnostically-related variables can also be specified. Each row in this file should represent a pair of variables which are diagnostically-related, such that combinations of exposures, stratifying variables, and outcomes containing these pairs can be filtered out by <code>preprocess()</code>. The first column should list the continuous variable(s) while the second column should list the binary outcome(s). For example, for all combinations containing both HbA1c and diabetes to be excluded, one row should be included in this .txt file, where HbA1c (continuous variable) is in the first column, and diabetes (binary outcome) is in the second column.

The following steps are performed by the <code>preprocess()</code> function. Sample input files can be obtained using <code>generate_sample_data()</code>. 

<div align="center">
    <img src="https://github.com/manalice/draft-algo/blob/main/preprocess.png" alt="Preprocess Steps" width="600">
</div>
<br>

Sample code for generating sample input files and pre-processing steps is shown below:

``` r
# Generate sample data
StratMRAlgo::generate_sample_data(output_dir = "/Your_Directory/Package_Testing_Sample_Data/") # modify

# Pre-process sample data
combinations <- StratMRAlgo::preprocess(
  input_dir = "/Your_Directory/Package_Testing_Sample_Data/", # modify
  output_dir = "/Your_Directory/1_preprocess/", # modify
  data = "data.csv",
  diagnostically_related = "diagnostically_related.txt",
  binary = NULL, # could specify "binary_comb_interest.txt"
  continuous = NULL) # could specify "continuous_comb_interest.txt"
```

The function outputs a list object containing two elements: <code>$binary</code> and <code>$continuous</code> corresponding to binary_comb.txt and continuous_comb.txt respectively if output directories are specified. Both data frames contain the filtered exposure, stratifying variable, and outcome combinations in each row with their biases, overall causal estimates, and the instrument strength. The data in the two data frames are summarized below:

| Column Name | Defintion | Data Frame(s) |
| ------ | ------ | ------ |
|  |  |  |
| EXP | Exposure name | <code>$binary</code> and <code>$continuous</code> |
| STRAT | Stratifying variable name | <code>$binary</code> and <code>$continuous</code> |
| CASE_STATUS | Binary outcome name | <code>$binary</code> |
| OUTCOME | Continuous outcome name | <code>$continuous</code> |
|  |  |  |
| X_C2_R2 | Estimation of the variance of X explained by C<sup>2</sup> from Binary Algorithm Model 1 | <code>$binary</code> |
| Y_C2_OR | Odds ratio estimating the effect of C<sup>2</sup> on Y from Binary Algorithm Model 2 | <code>$binary</code> |
| Y_X2_OR | Odds ratio estimating the effect of X<sup>2</sup> on Y from Binary Algorithm Model 3 | <code>$binary</code> |
| X_GC_R2 | Estimation of the variance of X explained by the GC interaction term from Continuous Algorithm Model 2 (used for the <code>all_methods()</code> function but not <code>algo()</code> for binary outcomes) | <code>$binary</code> |
|  |  |  |
| C_Y2_R2 | Estimation of the variance of C explained by Y<sup>2</sup> from Continuous Algorithm Model 1 | <code>$continuous</code> |
| X_GC_R2 | Estimation of the variance of X explained by the GC interaction term from Continuous Algorithm Model 2 | <code>$continuous</code> |
| X_C2_R2 | Estimation of the variance of X explained by C<sup>2</sup> from Continuous Algorithm Model 3 | <code>$continuous</code> |
| Y_C2_R2 | Estimation of the variance of Y explained by C<sup>2</sup> from Continuous Algorithm Model 4 | <code>$continuous</code> |
| Y_X2_R2 | Estimation of the variance of Y explained by X<sup>2</sup> from Continuous Algorithm Model 5 | <code>$continuous</code> |
| C_X2_R2 | Estimation of the variance of C explained by X<sup>2</sup> from Continuous Algorithm Model 6 | <code>$continuous</code> |
|  |  |  |
| wald_ratio_OR | Causal estimate of X on Y (OR) calculated using the allele score Wald Ratio (odds Y per unit X) | <code>$binary</code> |
| wald_ratio_OR_95_LL | Lower limit of the 95% CI of the odds ratio | <code>$binary</code> |
| wald_ratio_OR_95_UL | Upper limit of the 95% CI of the odds ratio | <code>$binary</code> |
|  |  |  |
| wald_ratio | Causal estimate of X on Y calculated using the allele score Wald Ratio (unit Y per unit X) | <code>$continuous</code> |
| wald_ratio_95_LL | Lower limit of the 95% CI of the causal estimate | <code>$continuous</code> |
| wald_ratio_95_UL | Upper limit of the 95% CI of the causal estimate | <code>$continuous</code> |
|  |  |  |
| X_G_R2 | Variance of X explained by G | <code>$binary</code> and <code>$continuous</code> |

If an output directory is defined, additional metadata with the excluded combinations (removed_combinations.txt), as well as overall allele score MR causal effects and the instrument strength for all pairs of X and Y (binary_overall.txt and continuous_overall.txt) will also be outputted.

### Part 2: Algorithm <a name="3"></a>

The <code>algo()</code> function automates the algorithm as described in our main manuscript. First, all the filtered combinations outputted by <code>preprocess()</code> are evaluated using the main DR by X method. Next, all Bonferroni significant combinations are evaluated as per the steps in either diagram below. 

*Algorithm for Binary Outcomes*
<br>
<div align="center">
    <img src="https://github.com/manalice/draft-algo/blob/main/binary.png" alt="Binary Algorithm" width="550">
</div>

<br>

*Algorithm for Continuous Outcomes*
<br>
<div align="center">
    <img src="https://github.com/manalice/draft-algo/blob/main/continuous.png" alt="Continuous Algorithm" width="850">
</div>
<br>

The <code>algo()</code> function automatically calls the <code>preprocess()</code> function. Therefore, it is not necessary to first pre-process the data before running this function. However, if you already completed the pre-processing steps and do not want to run <code>preprocess()</code> again, you can specify a non-null value for the parameter <code>combinations</code> such that it uses the final outputs from <code>preprocess()</code> directly from the R environment.

An optional input file for this function can be specified using the <code>causal_relationships</code> parameter. This is a file defining the directionality of causal relationships between the exposures and stratifying variables used. For example, age must be causal of BMI, as BMI cannot cause age. Each row in the causal relationships file should represent a pair wherein the causal direction is known. The first column should contain the causal parent, while the second column should contain the causal child. For example, <code>AGE</code> would be in the first column, while <code>BMI</code> would be in the second column. The column names should be <code>C1</code> and <code>C2</code> from left to right. There is no need to specify the suffices representing the variable type in this causal relationships file. Examples of all optional input files can also be generated using <code>generate_sample_data()</code>. 

Sample code for running the algorithm can be found below:

``` r
# Run algorithm
algo <- StratMRAlgo::algo(
  quantile_number = 50,
  input_dir = "/Your_Directory/Package_Testing_Sample_Data/", # modify
  output_dir = "/Your_Directory/2_algo/", # modify
  data = "data.csv",
  diagnostically_related = "diagnostically_related.txt",
  binary = NULL, # could specify "binary_comb_interest.txt"
  continuous = NULL, # could specify "continuous_comb_interest.txt"
  n = NULL,
  causal_relationships = "causal_relationships.txt",
  combinations = NULL, # could specify combinations = combinations if you do not wish to run the pre-processing steps again
  p_threshold = 0.05)
``` 

The function outputs a list object containing two elements: <code>$binary_results</code> and <code>$continuous_results</code> corresponding to binary_results.csv and continuous_results.csv respectively if output directories are specified. Both data frames contain the filtered exposure, stratifying variable, and outcome combinations in each row with their effect modification coefficients, p-values, and standard errors calculated using DR by X and any subsequent stratified MR methods, as well as whether they were determined to be overall Bonferroni significant by the algorithm. Data from <code>preprocess()</code> outputs pertaining to each combination remain in the two data frames. New variables which have not yet been defined in [Part 1](#2) are summarized below:

| Column Name | Defintion | Data Frame(s) |
| ------ | ------ | ------ |
|  |  |  |
| SIGNIFICANT | TRUE if determined to be overall significant by the algorithm; FALSE if determined not to be overall significant by the algorithm; INCONCLUSIVE otherwise | <code>$binary_results</code> and <code>$continuous_results</code> |
|  |  |  |
| DR_BY_X_P | P-value for the effect modification coefficient calculated using DR by X | <code>$binary_results</code> and <code>$continuous_results</code> |
| DR_BY_X_BETA | Effect modification coefficient calculated using DR by X; units described in [Part 4](#5) | <code>$binary_results</code> and <code>$continuous_results</code> |
| DR_BY_X_SE | Effect modification coefficient standard error calculated using DR by X | <code>$binary_results</code> and <code>$continuous_results</code> |
|  |  |  |
| DR_BY_G_P | P-value for the effect modification coefficient calculated using DR by G | <code>$binary_results</code> and <code>$continuous_results</code> |
| DR_BY_G_BETA | Effect modification coefficient calculated using DR by G; units described in [Part 4](#5) | <code>$binary_results</code> and <code>$continuous_results</code> |
| DR_BY_G_SE | Effect modification coefficient standard error calculated using DR by G | <code>$binary_results</code> and <code>$continuous_results</code> |
|  |  |  |
| RESID_STRAT_P | P-value for the effect modification coefficient calculated using Resid-Strat | <code>$continuous_results</code> |
| RESID_STRAT_BETA | Effect modification coefficient calculated using Resid-Strat; units described in [Part 4](#5) | <code>$continuous_results</code> |
| RESID_STRAT_SE | Effect modification coefficient standard error calculated Resid-Strat | <code>$continuous_results</code> |

If an output directory is defined, individual .csv files with additional metadata for each combination/method applied will be outputted as well.

### Part 3: All Stratified MR Methods <a name="4"></a>

All five stratified MR methods can also be run for a set of pre-processed exposure, stratifying variable, and outcome combinations. The <code>all_methods()</code> function has the same input parameters as <code>algo()</code>. However, instead of evaluating all combinations using DR by X and running Bonferroni significant combinations through the algorithm, it simply outputs the effect modification estimates, p-values, and standard errors for each combination using all methods and states whether bias is likely per method based on the simulation-derived thresholds.

Sample code for running all stratified MR methods (DR by X, DR by G, Resid-Strat, Resid-Wald, and Raw Strat) can be found below:

``` r
# Run all methods for each combination
all_methods <- StratMRAlgo::all_methods(
  quantile_number = 50,
  input_dir = "/Your_Directory/Package_Testing_Sample_Data/", # modify
  output_dir = "/Your_Directory/3_all_methods/", # modify
  data = "data.csv",
  diagnostically_related = "diagnostically_related.txt",
  binary = NULL, # could specify "binary_comb_interest.txt"
  continuous = NULL, # could specify "continuous_comb_interest.txt"
  n = NULL, # defaults to n = 1, no Bonferroni correction
  causal_relationships = "causal_relationships.txt",
  combinations = NULL, # could specify combinations = combinations if you do not wish to run the pre-processing steps again
  p_threshold = 0.05)
```

The function outputs a list object containing two elements: <code>$binary_results</code> and <code>$continuous_results</code> corresponding to binary_results_all.csv and continuous_results_all.csv respectively if output directories are specified. Both data frames contain the filtered exposure, stratifying variable, and outcome combinations in each row with their effect modification coefficients, p-values, and standard errors calculated using each method, as well as whether they were determined to be likely biased based on the thresholds shown in Table 2 of the manuscript. Data from preprocess() outputs pertaining to each combination remain in the two data frames. New variables which have not yet been defined in [Part 1](#2) are summarized below:

| Column Name | Defintion | Data Frame(s) |
| ------ | ------ | ------ |
|  |  |  |
| DR_BY_X_P | P-value for the effect modification coefficient calculated using DR by X | <code>$binary_results</code> and <code>$continuous_results</code> |
| DR_BY_X_BETA | Effect modification coefficient calculated using DR by X; units described in [Part 4](#5) | <code>$binary_results</code> and <code>$continuous_results</code> |
| DR_BY_X_SE | Effect modification coefficient standard error calculated using DR by X | <code>$binary_results</code> and <code>$continuous_results</code> |
| DR_BY_X_SIGNIFICANT | TRUE if determined to be significant (< p_threshold/n, where n defaults to 1) by DR by X; FALSE if determined not to be significant by DR by X | <code>$binary_results</code> and <code>$continuous_results</code> |
|  |  |  |
| DR_BY_G_P | P-value for the effect modification coefficient calculated using DR by G | <code>$binary_results</code> and <code>$continuous_results</code> |
| DR_BY_G_BETA | Effect modification coefficient calculated using DR by G; units described in [Part 4](#5) | <code>$binary_results</code> and <code>$continuous_results</code> |
| DR_BY_G_SE | Effect modification coefficient standard error calculated using DR by G | <code>$binary_results</code> and <code>$continuous_results</code> |
| DR_BY_G_SIGNIFICANT | TRUE if determined to be significant (< p_threshold/n, where n defaults to 1) by DR by G; FALSE if determined not to be significant by DR by G | <code>$binary_results</code> and <code>$continuous_results</code> |
|  |  |  |
| RESID_STRAT_P | P-value for the effect modification coefficient calculated using Resid-Strat | <code>$binary_results</code> and <code>$continuous_results</code> |
| RESID_STRAT_BETA | Effect modification coefficient calculated using Resid-Strat; units described in [Part 4](#5) | <code>$binary_results</code> and <code>$continuous_results</code> |
| RESID_STRAT_SE | Effect modification coefficient standard error calculated Resid-Strat | <code>$binary_results</code> and <code>$continuous_results</code> |
| RESID_STRAT_SIGNIFICANT | TRUE if determined to be significant (< p_threshold/n, where n defaults to 1) by Resid-Strat; FALSE if determined not to be significant by Resid-Strat | <code>$binary_results</code> and <code>$continuous_results</code> |
|  |  |  |
| RAW_STRAT_P | P-value for the effect modification coefficient calculated using Raw Strat | <code>$binary_results</code> and <code>$continuous_results</code> |
| RAW_STRAT_BETA | Effect modification coefficient calculated using Raw Strat; units described in [Part 4](#5) | <code>$binary_results</code> and <code>$continuous_results</code> |
| RAW_STRAT_SE | Effect modification coefficient standard error calculated Raw Strat | <code>$binary_results</code> and <code>$continuous_results</code> |
| RAW_STRAT_SIGNIFICANT | TRUE if determined to be significant (< p_threshold/n, where n defaults to 1) by Raw Strat; FALSE if determined not to be significant by Raw Strat | <code>$binary_results</code> and <code>$continuous_results</code> |
|  |  |  |
| RAW_STRAT_P | P-value for the effect modification coefficient calculated using Raw Strat | <code>$binary_results</code> and <code>$continuous_results</code> |
| RAW_STRAT_BETA | Effect modification coefficient calculated using Raw Strat; units described in [Part 4](#5) | <code>$binary_results</code> and <code>$continuous_results</code> |
| RAW_STRAT_SE | Effect modification coefficient standard error calculated Raw Strat | <code>$binary_results</code> and <code>$continuous_results</code> |
| RAW_STRAT_SIGNIFICANT | TRUE if determined to be significant (< p_threshold/n, where n defaults to 1) by Raw Strat; FALSE if determined not to be significant by Raw Strat | <code>$binary_results</code> and <code>$continuous_results</code> |
|  |  |  |
| RESID_WALD_P | P-value for the effect modification coefficient calculated using Resid-Wald | <code>$continuous_results</code> |
| RESID_WALD_BETA | Effect modification coefficient calculated using Resid-Wald; units described in [Part 4](#5) | <code>$continuous_results</code> |
| RESID_WALD_SE | Effect modification coefficient standard error calculated Resid-Wald | <code>$continuous_results</code> |
| RESID_WALD_SIGNIFICANT | TRUE if determined to be significant (< p_threshold/n, where n defaults to 1) by Resid-Wald; FALSE if determined not to be significant by Resid-Wald | <code>$continuous_results</code> |
|  |  |  |
| DR_BY_X_BIAS | POSSIBLE_BIAS indicates potential bias for DR by X; NO_BIAS_DETECTED indicates no bias detected for DR by X | <code>$binary_results</code> and <code>$continuous_results</code> |
| DR_BY_G_BIAS | POSSIBLE_BIAS indicates potential bias for DR by G; NO_BIAS_DETECTED indicates no bias detected for DR by G | <code>$binary_results</code> and <code>$continuous_results</code> |
| RESID_STRAT_BIAS | POSSIBLE_BIAS indicates potential bias for Resid-Strat; NO_BIAS_DETECTED indicates no bias detected for Resid-Strat | <code>$binary_results</code> and <code>$continuous_results</code> |
| RESID_WALD_BIAS | POSSIBLE_BIAS indicates potential bias for Resid-Wald; NO_BIAS_DETECTED indicates no bias detected for Resid-Wald | <code>$continuous_results</code> |
| RAW_STRAT_BIAS | POSSIBLE_BIAS indicates potential bias for Raw Strat; NO_BIAS_DETECTED indicates no bias detected for Raw Strat | <code>$binary_results</code> and <code>$continuous_results</code> |

Similar to <code>algo()</code>, if an output directory is defined, individual .csv files with additional metadata for each combination/method applied will be outputted as well.

### Part 4: Individual Stratified MR Methods <a name="5"></a>

Any of the five methods can also be used to evaluate any exposure, stratifying variable, and outcome combination of interest using the function <code>mr()</code>. Doubly ranked methods were adapted from [Tian et al. (2023)](https://doi.org/10.1371/journal.pgen.1010823) and residual stratification methods were adapted from [Staley and Burgess (2017)](https://doi.org/10.1002/gepi.22041) and their respective Github repositories. Sample code can be found below:

``` r
# Test one specific combination using DR by X
dr_x_sample <- StratMRAlgo::mr(method = "by_x", # can specify 'by_x', 'by_g', 'resid_strat', 'raw_strat', or 'resid_wald'
   exposure_name = "X1_EXP", 
   outcome_name = "Y2_CASE_STATUS", 
   stratifying_name = "C2_STRAT", 
   quantile_number = 50, 
   data = "data.csv", 
   input_dir = "/Your_Directory/Package_Testing_Sample_Data/", # modify
   output_dir = "/Your_Directory/4_dr_x_sample/") # modify
```

The output is a list containing three elements: <code>$p</code>, <code>$beta</code> (effect modification coefficient), and <code>$se</code> obtained from the meta-regression of stratum estimates. For binary outcomes, <code>$beta</code> is reported in log-odds Y per unit X per unit C. For continuous outcomes <code>$beta</code> is reported in units Y per unit X per unit C.

Note that when age is the outcome, sex is the outcome, or age is the stratifying variable, the respective variable will be automatically removed as a covariate for the calculation of the stratum-specific MR estimates. This applies to <code>algo()</code>, <code>all_methods()</code>, and <code>mr()</code>.

## License <a name="6"></a>
GNU General Public License v3.0

## Contact information <a name="7"></a>
Any queries pertaining to the StratMRAlgo R package or methodological framework can be addressed to either: Alice Man (mana3@mcmaster.ca) or Guillaume Paré (pareg@mcmaster.ca).

## Citation <a name="8"></a>
Citation will be made available in the future.


