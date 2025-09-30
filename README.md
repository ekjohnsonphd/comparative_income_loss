# Comparative Income Loss Analysis

A comprehensive analysis of income effects from disease onset using Danish health registry data. This project examines how different diseases affect income trajectories through  matching.

## Project Overview

This study investigates the impact of disease related health shock on individual income using longitudinal data from Danish health registries (2000-2018). The analysis employs exact matching methods and estimates the realized income loss following hospitalization and the standardized loss which accounts for differential risk of each diagnosis. We estimate losses for 10 years following exposure and report losses stratified along population characteristics.

### Key Features
- **Diseases Analyzed**: Depression, stroke, alcohol-related disorders
- **Data Sources**: Danish National Patient Registry (LPR), Danish Breast Cancer Group register (DBCG), Danish Income Statistics Register, the Civil Registration Register, the Integrated Student Register
- **Methods**: Exact matching. Variations of matching variables and lookback windows are tested in the analysis. Standardization along matching characteristics is performed to estimate income losses which control for differential disease risk. We also consider the impact of differential risk of mortality as a sensitivity analysis.
- **Time Period**: 2000-2018 with up to 10-year follow-up

## Project Structure

### Setup and Configuration
- **`environment_setup.R`** - Global environment configuration, package loading, and path definitions

### Data Preparation Pipeline (`01_data_prep/`)

#### `01_hospital_utilization.R`
**Purpose**: Calculate hospital utilization metrics for matching covariates
- Processes LPR (National Patient Registry) data across multiple years
- Aggregates inpatient and outpatient hospital contacts
- Creates 5-year lookback utilization histories for each individual
- Exports utilization counts for use in matching algorithms

#### `02_calculate_comorbidities.R`
**Purpose**: Generate comorbidity indices using standard medical classification systems
- Calculates Charlson Comorbidity Index and Elixhauser Comorbidity Index from ICD information in the registers
- Produces monthly comorbidity scores for matching

#### `03_populations.R`
**Purpose**: Identify and define study populations from health registries
- Extracts individuals with target conditions using ICD-10 regex patterns
- Combines data from multiple sources (LPR, private clinics, psychiatric hospitals)
- Handles diagnosis timing and data source prioritization
- Creates population-specific datasets for each disease category

#### `04_data_prep.R`
**Purpose**: Main data integration and preparation pipeline
- Merges health registry data with income and demographic information
- Applies inflation adjustments to monetary values using Danish consumer price indices
- Creates analysis-ready datasets with proper variable encoding
- Generates income quintiles and socioeconomic indicators

### Analysis Pipeline (`02_analysis/`)

#### `01_matching.R`
**Purpose**: Implement exact matching to create comparable treatment and control groups
- Defines multiple matching specifications with increasing covariate sets
- Performs exact matching on demographics, comorbidities, and utilization patterns
- Exports matched datasets for downstream analysis

#### `02_probability_of_exposure.R`
**Purpose**: Analyze exposure patterns and matching completeness
- Fits decision trees to predict exposure probability based on covariates
- Examines factors associated with missing exposure data
- Creates descriptive statistics and visualization of exposure patterns
- Assesses potential selection bias in the matching process

#### `03_matched_outcomes_data.R`
**Purpose**: Prepare longitudinal outcome data for causal analysis
- Links matched individuals to 10-year income follow-up data
- Tracks mortality status and creates survival-adjusted outcomes
- Estimates conditional impacts
- Handles missing data and creates analysis-ready longitudinal datasets

#### `04_mortality_stratified_outcomes.R`
**Purpose**: Conduct mortality-stratified analysis to address competing risks
- Stratifies analysis by survival status to assess competing risk bias
- Creates separate estimates for survivors vs. overall population
- Generates visualizations comparing mortality-adjusted and unadjusted effects
- Exports stratified results for policy interpretation

#### `xx_match_functions.R`
**Purpose**: Utility functions for matching and causal inference
- **`deaths_fcn()`** - Handles mortality in follow-up data with counterfactual survival assumptions
- **`deaths_cfactual()`** - Adjusts income estimates for mortality in CACE analysis
- **`calculate_ace()`** - Calculates Average Causal Effects with proper standard error propagation
- **`identify_exposed()`** - Classifies individuals as exposed, control, or excluded based on timing
- **`standardized_unstandardized()`** - Creates both real and population-standardized estimates

## Methodology

### Matching Strategy
The analysis uses exact matching on key confounders to create comparable treatment and control groups:
- **Demographics**: Age, sex, education, income history
- **Health Status**: Comorbidity indices, prior healthcare utilization
- **Geographic**: Municipality of residence
- **Temporal**: Calendar year effects

### Causal Inference Approach
1. **Exact Matching**: Creates balanced treatment and control groups
2. **CACE Estimation**: Focuses on compliers (those who would be treated when assigned to treatment)
3. **Mortality Adjustment**: Accounts for competing risk of death using counterfactual survival scenarios
4. **Sensitivity Analysis**: Multiple matching specifications and robustness checks

### Data Processing
- **Parallel Computing**: Utilizes multi-core processing for large dataset operations
- **Memory Efficiency**: Arrow/Parquet format for optimized data storage and retrieval
- **Temporal Indexing**: Standardized time indexing relative to match years for consistent analysis

## Key Outputs

- **Income Loss Estimates**: Causal effect of disease onset on subsequent income
- **Mortality-Stratified Results**: Separate estimates accounting for survival bias
- **Robustness Checks**: Multiple matching specifications and sensitivity analyses
- **Policy Metrics**: Percentage income losses and absolute monetary effects

## Dependencies

### R Packages
- **Data Processing**: `data.table`, `arrow`, `dplyr`
- **Matching**: `MatchIt` for propensity score methods
- **Modeling**: `rpart` for decision trees
- **Parallel Computing**: `parallel`, `mclapply`
- **Visualization**: `ggplot2`, `gridExtra`

### Data Requirements
- Access to Danish health registries (LPR, private clinics, psychiatric hospitals)
- Individual-level income data with long-term follow-up
- Population demographics and mortality records

## Usage Notes

1. **Sequential Execution**: Files should be run in numerical order within each directory
2. **Memory Requirements**: Large datasets require sufficient RAM for data.table operations  
3. **Parallel Processing**: Adjust core counts in scripts based on available computational resources
4. **Data Privacy**: All code assumes anonymized data meeting Danish data protection requirements

## Contact

For questions about methodology or implementation, please refer to the documented code or contact the project maintainer.
