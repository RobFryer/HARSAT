## Change log

### Version 1.0.3

This release is that used to run the OSPAR 2025 CEMP assessment.

#### ICES data extractions

`read_data` and `read_contaminants` now explicitly look for the variables `amap_arctic_lme` and `casenumber` in the contaminants data file when argument `data_format` is set to `"ICES"`. These variables were introduced in ICES extractions in mid-June 2024. A warning is printed if the variables are missing. `read_stations` also explicitly looks for `amap_arctic_lme` in the ICES station dictionary.

* `amap_arctic_lme` contains AMAP large marine ecosystem regional information 
* `casenumber` is the accession identification id in the AIMS (Accession Information Management System) system introduced in 2024; `casenumber` replaces accessionid from the now discontinued DAD system

#### AMAP assessments

There are two changes that affect data processing for AMAP assessments (when argument `purpose` in `read_data` is set to `"AMAP"`).

First, the default behaviour by which `read_data` matches stations to data from an ICES extraction is now to restrict eligible stations in those in the AMAP area. Previously it didn't matter where the stations were. The default behaviour can be changed using the `add_stations` element of the `control` argument of `read_data`. See the help file of `add_stations` for more information about station matching.

Second, the default values for the elements of control$region are now

* id: "amap_arctic_lme"
* names: "AMAP_arctic_lme"
* all: FALSE (this is because there are discrepancies between shape files for the AMAP area and the AMAP large marine ecosystems)

Users of external data should either rename the variable with regional information in the station dictionary to `amap_arctic_lme` or set `control$region` so that it picks up the existing name of this variable.

#### Summations involving censored values

An argument `sum_censored` has been added to `determinand.link.sum` which allows the user to control whether censored values are included when computing a sum.

`sum_censored = TRUE` (default) maintains the previous behaviour where the sum is computed by adding all the non-censored measurements and the censored values of the censored measurements. If all the measurements are censored and all the censored values are equal to the limit of detection, then the sum will be the sum of the limits of detection.

`sum_censored = FALSE` typically treats each censored value as zero when calculating the sum. Usually, the output is the sum of the non-censored measurements. There are two exceptions:

* when all measurements are censored, the output is the largest censored value (and is flagged as a censored measurement)
* when the sum of the non-censored measurements is less than the largest censored value, the output is the largest censored value (and is flagged as a censored measurement); this will be unusual and is most likely to occur when using weights (and small weights are applied to non-censored measurements and a large weight is applied to the censored measurements)

The same argument has also been added to customised link functions such as `determinand.link.BBKF` which call `determinand.link.sum`.

It is likely that most users will use `sum_censored = FALSE` since this is compatible with Marine Strategy advice on the treatment of censored values. If so, this might become the default option in a future release.

#### Summary tables - variable renaming

Several variables in the output from `write_summary_table` have been renamed:  

* `p_nonlinear`    --> `p_nonlinear_trend`
* `p_linear`       --> `p_linear_trend`
* `p_overall       --> `p_overall_trend`
* `p_linear_trend` --> `p_overall_change`
* `p_recent_trend` --> `p_recent change`
* `linear_trend`   --> `overall_change`
* `recent_trend`   --> `recent_change`

The first three variables report evidence of systematic trends in the data. The remaining variables all relate to changes in concentration between the start and end of either the whole time series or the 'recent' part of the time series (typically the last 20 monitoring years). The renaming was prompted by confusion about the original names `p_linear_trend` and `linear_trend` which suggest linearity that is not always the case. The harsat team struggled to come up with suitable alternatives and hopes that the new names are less confusing (if not crystal clear). 

#### Calculation of recent_change

The `overall_change` is the change in concentration over the whole monitoring period. It is only calculated if there are at least five years of data with at least one non-censored measurement.

The `recent_change` is the change in concentration in recent years, typically taken to be the last twenty monitoring years. In previous releases, `recent_change` was only calculated if `overall_change` had been calculated and there were at least five years of data in the recent period. However there was no requirement on the number of years with non-censored measurements in the recent period. This meant that the evidence base for `recent_change` could be weak (e.g. if only one or two years in the recent period had non-censored measurements). Occasionally, `recent_change` could be undefined for long time series with infrequent monitoring and many censored measurements. 

The user can now control the calculation of `recent_change` using the `control` argument of `run_assessment`. Specifically, `control$recent_change` has two components:

* `n_year_fit` - default 5L
* `n_year_positive` - default 5L

where `n_year_fit` is the required number of years of data in the recent period and `n_year_positive` is the required number of years of data with at least one non-censored measurement in the recent period.

By default, `recent_change` will now only be calculated if there are at least 5 years of data with non-censored measurements in the recent period.

The behaviour of previous releases can be replicated (almost) by seting `n_year_fit` to 5L and `n_year_positive` to 2L (the smallest value allowed to avoid pathological behaviour).

#### Imposex assessments - annual indices equal to zero 

Imposex assessments (based on individual measurements) involve the estimation of cut-points that measure the transition from one imposex stage to the next on the latent odds scale. This estimation pools the data from several (often many) timeseries to improve the precision of the cut-point estimates and is done before the call to `run_assessment`. 

Data from timeseries / year combinations where all the individual measurements are zero (equivalently, the annual index is zero) contain no information on cut-points and are now removed from the estimation procedure. This improves convergence. The estimation routine has been renamed as `ordinal_theta_est` (previously `ctsm.VDS.index.opt )

Estimation of the cut-points is followed by estimation of confidence intervals for the annual indices. This is possible because the model output also includes estimates of each index. The confidence intervals are estimated from the posterior distribution of the parameter estimates. However, this does not work for zero indices whose fitted values would be infinite on the latent scale. A good solution would be to calculate likelihood intervals, but this is a challenging numerical problem and has been left for a future release. Instead, the previous very ad-hoc approach has been replaced by a moderately ad-hoc approach. Specifically, an infinite value for a zero index is replaced by a value larger than all the fitted values for positive indices. This is done by taking all the positive indices, fitting a linear model of fitted value against square root index, and predicting what the fitted value should be when the index is zero. The square root scale is based on the typical relationship between fitted values and indices observed in the OSPAR 2025 assessment. The associated standard error is taken to be the upper 90th quantile of the standard errors associated with positive indices (with a suitable adjustment for the number of measurements used in each), which is reasonable since standard errors increase as the number of zero measurements increase but not in a very predictable manner (at least not in the OSPAR 2025 assessment). The upper confidence limit on a zero index is then estimated by simulating from the posterior distribution of the parameter estimates as before. (The lower confidence limit is zero by definition.) The estimation routine has been renamed `ordinal_theta_cl` (previously `ctsm.VDS.index.cl`).

#### Error handling 

Error handling for the `sample` variable in the contaminants data file (input to `read_data`) has been tightened up. This is only likely to affect external data.  

An error is now thrown if there are any missing values.

There are now checks for non-unique sample identifiers, for example when the same sample identifier has been used in different years, or for different stations or species. If these are found, a warning is printed exhorting the user to sort out their data. However, the code also attempts to create unique sample identifiers by pasting together `year`, `station_code`, `species` (biota only) and `sample`. 

#### Deleted functions

* `determinand.link.TEQDFP` was deprecated in release 1.0.2 and is now deleted; use `determinand.link.sum` instead.

#### Minor bug fixes

* `write_summary_table` now works for biological effects assessments with non-standard summary variables (e.g. imposex_class)
* `ctsm_symbology_OSPAR` now works when there are no non-parametric tests of status (e.g. when only imposex is assessed)
* `determinand.link.sum`: the uncertainty of the summed concentration is now independent of the order of the data; testing suggests that data order previously affected the uncertainty of a tiny proportion of samples, typically where all the measurements were censored
* `plot_assessment`: no extra page is produced in pdf plots of assessments with thresholds
* `write_summary_tabe`: a header is always produced when a new summary table is written to file; previously failed if the `append` argument was set to `TRUE` and there was no existing summary file to append to


### Version 1.0.2

This release is that used to run the OSPAR 2024 CEMP assessment.

#### Weighted sums for TEQs etc.

Weighted sums of concentrations are now calculated using `determinand.link.sum`. The weights are supplied using the `weights` argument. This function supersedes `determinand.link.TEQDFP` which is now deprecated and will be removed at the next release.

Updated and corrected World Health Organisation TEFs for dioxins, furans and planar PCBs are now available in the data object `info_TEF`: there are four versions available:  

* DFP_2005; the 2005 values  
* DFP_2022; the 2022 values  
* DFP_HOLAS3; the values used in the HOLAS3 assessment; these are the 2005 values excluding three PCBs and are included for backward compatibility  
* DFP_CEMP; the values used in CEMP assessments <= 2024; these are the 2005 values excluding three PCBs and with the TEF for CDFO ten times too small; they are included for backward compatibility  

DFP_2005 and DFP_2022 use determinand OCDF (the correct code) rather than CDFO (which is a grouped determinand code). 

#### Auxiliary variables

Key auxiliary variables can now be plotted in `plot_assessment`. The default variables are:  

* biota: concentration, LNMEA, DRYWT%, LIPIDWT%  
* sediment: non-normalised concentration, normalised concentration, AL, CORG  
* water: no plots are currently produced  

The choice of auxiliary variables can be altered using the `auxiliary` argument, although the options here are still limited.

The merging of auxiliary variables with the data is now controlled using the `control` argument of `read_data`. `control$auxiliary$by_matrix` is a list which determines which auxiliary variables are merged by sample and matrix and which are just merged by sample. The default values are `c("DRYWT%", "LIPIDWT%)` for biota and `"all"` for sediment and water. Thus, by default, dry weight and lipid weight measurements are matched with chemical concentrations in the same tissue (matrix). but all other auxiliary variables in biota are matched at the sample level. For sediment (and water) all auxiliary variables (e.g. aluminium and organic carbon measurements) are matched with chemical concentrations in the same grain size fraction.

#### Reporting

`report_assessment` now has the full functionality required for the OSPAR CEMP assessment. This includes:  

* scatterplot matrices of concentrations of related compounds using the non-exported function `plot_multidata`  
* plots of assessments of related compounds using the non-exported function `plot_multiassessment`  
* plots of contaminant ratios using the non-exported function `plot_ratio`  

There is still some work required to make `report_assessment` suitable for all purposes.

#### Minor bug fixes  

* ensures that if `uncertainty` column is present in external data then so is `unit_uncertainty` and vice versa
* `plot_assessment` now correctly plots the 90% two-sided confidence intervals on VDSI estimates from imposex assessments 
* correct treatment of censoring data in `determinand.link.sum`
* `ctsm.check.sex.biota` now works for any auxiliary variable
* `get_timeseries` now always shows the series identifier for each timeseries
* `estimate_uncertainties` now traps for the case then `DRYWT%` and `LIPIDWT%` are not specified as auxiliary variables


### Version 1.0.1

Updates (mostly) required to run the OSPAR 2024 CEMP assessment. 

#### Data import

For OSPAR and HELCOM style assessments, data from Germany are now matched to stations by name for 2023 onwards. This applies to biota, sediment and water. Note that for HELCOM, biota data from Germany are already matched by name for all years.  

#### Uncertainty processing

harsat 1.0.0 replaced implausibly large relative uncertainties ($>=$ 100%) and replaced them with imputed values. However, implausibly small relative uncertainties were not changed. The code now replaces relative uncertainties $<=$ 1% with imputed values. 

The defaults can be changed using `control$relative_uncertainty` in `read_data`. To replicate the defaults in harsat 1.0.0, set `control$relative_uncertainty = c(0, 100)`. To keep all uncertainties, regardless of how ridiculous they are, set `control$relative_uncertainty = c(0, Inf)`.

Two minor bug fixes:  

* relative uncertainties were being filtered for all distributional types, but this is only a reliable procedure for determinands with `distribution == "lognormal"`; the checks are now only applied to lognormal data   
* some biological effect data with distributions other than normal or lognormal were being incorrectly deleted; this has now been corrected  

The oddity files have been updated to show:  

* implausible_uncertainties_reported.csv - all reported uncertainties that are replaced by imputed values  
* missing_uncertainties.csv - all uncertainties (normal or lognormal data) that are not reported and can't be imputed  
* implausible_uncertaintes_calculated.csv - all uncertainties that are calculated during the data processing (e.g. during normalisation) that are implausible and are set to missing  

#### Uncertainty coefficients

The function `ctsm_uncrt_workup` and related supporting functions are used in OSPAR assessments to update the fixed and proportional standard deviations which are subsequently used to impute missing uncertainties. These functions were ignored during the initial development of harsat and are now harsat compatible.  

#### Biological effect assessments

Imposex assessments: these are now fully reproducible with seeds for random number generation provided in the calls to `ctsm.VDS.cl` and `assess_imposex` 

Assessment functions for negative binomial data have been added. Negative binomial data includes MNC - the number of micronucleated cells.

#### Reporting

`report_assessment` generates default file names. These are based on the series identifier with additional station information. It is now possible to override this behaviour for a single report by providing a different file name using the `output_file` argument.

#### Reference tables  

* new values added to method_extraction table  

#### Minor bug fixes  

* correct behaviour of argument `return_early` in `create_timeseries`  
* pass `info` component of the harsat object to `determinand.link.sum`, `determinand.link.replace`, and `determinand.link.imposex`  
* ensure early return from `ctsm_convert_basis` when there is nothing to convert (avoids issues e.g. when all the data are biological effects)  
* ensure SURVT (in pargroup B-BIO) is recognised as a biological effect in `ctsm_get_datatype` (SURVT is the only determinand in this pargroup that isn't an auxiliary variable)  
* pass `good_status` to assessment functions for data with distributions other than normal and lognormal
* trap pathological case in estimation of `prtrend`; see #436
* ensure `ctsm_OHAT_legends` uses the symbology as specified in `write_summary_table` 


### Version 1.0.0

- Initial public release

### Version 0.1.3

- Various fixes

### Version 0.1.2

- Fixed issues when packaged; see #326, #328
- Updated AMAP data and packaging; see #329

### Version 0.1.1

- Fixed issue with auxiliary variables: see #289
- Small documentation improvements
- Added build processes for package bundles

### Version 0.1.0

- Initial release
