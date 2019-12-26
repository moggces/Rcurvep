


tcplHillVal <- function(logc, tp, ga, gw, bt = 0) {

  bt + (tp - bt)/(1 + 10^((ga - logc)*gw))

}

tcplHillConc <- function(val, tp, ga, gw, bt = 0) {

  ga - log10((tp - bt)/(val - bt) - 1)/gw

}

tcplHillACXX <- function(XX, tp, ga, gw, bt = 0) {

  y <- tp * XX/100
  ga - log10((tp - bt)/(y - bt) - 1)/gw

}


ObjHillnorm <-  function(p, lconc, resp) {

  mu <- p[1]/(1 + 10^((p[2] - lconc)*p[3]))
  sum(dnorm(x = resp,mean = mu,sd = p[4],log = TRUE))

}

# tcplObjCnst: Generate a constant model objective function to optimize (from tcpl package)
tcplObjCnst <- function(p, resp) {

  ### This function takes creates an objective function to be optimized using
  ### the starting constant model parameter, and response.
  ###
  ### Arguments:
  ###   p:     a numeric vector of length 1 containg the starting values for
  ###          the constant model, in order: log error term
  ###   lresp: a numeric vector containing the response values to produce the
  ###          objective function
  ###
  ### Value:
  ###   An objective function for the constant model and the given resp data

  mu <- 0
  sum(dt((resp - mu)/exp(p[1]), df = 4, log = TRUE) - p[1])

}


# tcplObjHill: Generate a hill model objective function to optimize (from tcpl package)

tcplObjHill <- function(p, lconc, resp) {

  ### This function takes creates an objective function to be optimized using
  ### the starting hill parameters, log concentration, and response.
  ###
  ### Arguments:
  ###   p:     a numeric vector of length 4 containg the starting values for
  ###          the hill model, in order: top, log AC50, hill
  ###          coefficient, and log error term
  ###   lconc: a numeric vector containing the log concentration values to
  ###          produce the objective function
  ###   lresp: a numeric vector containing the response values to produce the
  ###          objective function
  ###
  ### Value:
  ###   An objective function for the hill model and the given conc-resp data

  mu <- p[1]/(1 + 10^((p[2] - lconc)*p[3]))
  sum(dt((resp - mu)/exp(p[4]), df = 4, log = TRUE) - p[4])

}


# tcplObjGnls: Generate a gain-loss model objective function to optimize (from tcpl package)

tcplObjGnls <- function(p, lconc, resp) {

  ### This function takes creates an objective function to be optimized using
  ### the starting gain-loss parameters, log concentration, and response.
  ###
  ### Arguments:
  ###   p:     a numeric vector of length 5 containg the starting values for
  ###          the gain-loss model, in order: top, gain log AC50, gain hill
  ###          coefficient, loss log AC50, loss hill coefficient and log error
  ###          term
  ###   lconc: a numeric vector containing the log concentration values to
  ###          produce the objective function
  ###   lresp: a numeric vector containing the response values to produce the
  ###          objective function
  ###
  ### Value:
  ###   An objective function for the gain-loss model and the given conc-resp
  ###   data

  gn <- 1/(1 + 10^((p[2] - lconc)*p[3]))
  ls <- 1/(1 + 10^((lconc - p[4])*p[5]))
  mu <- p[1]*gn*ls
  sum(dt((resp - mu)/exp(p[6]), df = 4, log = TRUE) - p[6])

}
