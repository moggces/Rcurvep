# Rcurvep 0.1.0

* initial release

# Rcurvep 0.2.0

* new: add diagnostic plot
* change: the separate the id parameters into endpoint and direction parameters in identify_basenoise_threshold function

# Rcurvep 0.2.1

* change: TLOG(curvep.r) is to -24

# Note
* Curvep MASK parameter is not yet implemented

# Todo
* fix n_sample = NULL and simulate_percent_resp for direction = -1
* add a parameter in extract_curvep_data function to get simulated response curve
* simplify the run_curvep_job output file so toJSON/fromJSON can work together
* allow user to define p1 and p2 in the identify_basenoise_threshold
