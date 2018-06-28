# Rcurvep 0.1.0

* initial release

# Rcurvep 0.2.0

* new: add diagnostic plot
* change: the separate the id parameters into endpoint and direction parameters in identify_basenoise_threshold function

# Rcurvep 0.2.1

* change: TLOG(curvep.r) is to -24

# Rcurvep 0.3

* new: simplify_output is added for run_curvep_job() to reduce the amount of output
* change: levels in curvep() is seperated into ECxx, Cxx, and xx
* change: combine conc_hl to act function in extract_curvep_data()
* change: input column became a list of list structure


# Todo
* implement Curvep MASK parameter
* figure out toJSON/fromJSON conversion
* allow user to define p1 and p2 in the identify_basenoise_threshold
