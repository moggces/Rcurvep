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


# Todo
* implement Curvep MASK parameter
* implement stratified bootstrap
* figure out toJSON/fromJSON conversion

