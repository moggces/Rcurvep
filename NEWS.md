# Rcurvep 0.1.0

* initial release

# Rcurvep 0.2.0

* new: add diagnostic plot
* change: the separate the id parameters into endpoint and direction parameters in identify_basenoise_threshold()

# Rcurvep 0.2.1

* change: TLOG(curvep.r) is to -24

# Rcurvep 0.3

* new: simplify_output is added for run_curvep_job() to reduce the amount of output
* change: levels in curvep() is seperated into ECxx, Cxx, and xx
* change: combine conc_hl to act function in extract_curvep_data()
* change: input column became a list of list structure

# Rcurvep 0.3.1

* change: expose the p1 and p2 for identify_basenoise_threshold()

# Rcurvep 0.4
* new: add modifier parameter for extract_curvep_data()
* new: add reparam_curvep_job() to recalculate based on run_curvep_job() complex output

# Rcurvep 0.4.1
* new: add the abline line that is used to calculate distance in the dignostic plot
* new: add generate_diagnostic_plot() and cal_exponential_inflection()
* new: add an comment (thresDistComment) that shows if the identified baseline threshold is usable
* change: thresCurva approach is removed

# Rcurvep 0.5
* new: add a method to calculate the baseline noise threshold based on the exponential fitting (thresDist_exp column)
* new: statistical outputs from identify_basenoise_threshold() are reported ( y_exp_fit, y_lm_fit, dist2l_exp, cor_exp_fit, cor_lm_fit)
* change: identify_basenoise_threshold() output becomes a list instead of a tibble; first element contains the statistical data and the second element contains the baseline noise threshold information
* change: cal_exponential_inflection() is renamed as cal_knee_point()
* change: p1, p2 => p1_raw, p2_raw; dist2l => dist2l_raw

# Rcurvep 0.5.1
* new: summarize_curvep_output() to get confidence interval of activites by bypassing the extract_curvep_data()
* change: use SSasymp to fit the exponential decay curve

# Rcurvep 0.5.2
* change: fix comment_threshold()

# Rcurvep 0.6.0
* new: print.rcurvep_thres_stats() -> print()
* change: run_curvep_job() -> run_curvep_batch()
* change: extract_curvep_data() -> withdraw.rcurvep_out_nested() -> withdraw()
* change: reparam_curvep_job() -> reparam.rcurvep_out_nested() -> reparam()
* change: extract_curvep_data(summary = TRUE) -> summary.rcurvep_out() -> summary()
* change: identify_basenoise_threshold() -> get_baseline_threshold()
* change: generate_diagnostic_plot() -> plot.rcurvep_thres_stats() -> plot()



# Todo
* implement Curvep MASK parameter
* implement stratified bootstrap
* figure out toJSON/fromJSON conversion

