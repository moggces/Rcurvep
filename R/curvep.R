Impute <- function(lgC, V, L, JUNK = -999, strict = TRUE)
{
  #returns estimated log conc., at which response level L is reached.
	s <- 1
	e <- length(lgC)
	aV <- abs(V)
	#block below deals with the rare case of U-shape curves with gradually (< MXDV) decreasing shoulder detected as plateau
	u <- e
	while (u > s)  if (aV[u-1] > aV[u]) u <- u - 1	else break

	#------
	if (strict)
	{
	  if (aV[s] >= L) return (lgC[s]) #07.10.17 this check has to be done first, otherwise result is voided for inverse curves
	  if (aV[u] < L) return (JUNK)
	}

	xT <- JUNK
	z <- s
	while (z < e)
	{
	  z <- z + 1
		if ( ( (aV[z] >= L) || (z == e) ) && ( aV[z] > aV[s] ) )
		{
			xT  <- lgC[s] + (lgC[z] - lgC[s])*(L - aV[s]) / abs(V[z] - V[s])
			break
		}
		s <- z
	}

	return ( xT )
}

get_monotonics	<- function (vals, vdif = 0, bads = NULL)
{
  #returns array of decreasing/constant/increasing (1/0/-1) intervals in vals.
  #Intervals are constant if values differ <= vdif
  #NB: inters[c] correspond to (vals[c-1]; vals[c]). For c = 0, it corresponds to interval from infinite dilution to the first test conc.
  #NB: for gene expression data with control included, this needs to be adjusted, because infinite dilution is explicitly given

  n <- length(vals)

  mask <- bads
  if (is.null(mask)) mask <- rep(FALSE, n)

  inters <- rep(0, n)
  prev <- 0
  for (c in 1:n)
  {
    if (mask[c]) next
    if (abs(vals[c] - prev) > vdif)
    {
      if (vals[c] < prev) inters[c] <-  1
      if (vals[c] > prev) inters[c] <- -1
    }
    prev <- vals[c]
  }

  return (inters)
}

#' The Curvep function to process one set of concentration-response data
#'
#' The relationship between concentration and response has to be 1 to 1.
#'   The function is the backbone of [run_rcurvep()] and [combi_run_rcurvep()].
#'
#' @param Conc Array of concentrations, e.g., in Molar units, can be log-transformed, in which case internal log-transformation is skipped.
#' @param Resp Array of responses at corresponding concentrations, e.g., raw measurements or normalized to controls.
#' @param Mask array of 1/0 flags indicating invalidated measurements (default = NULL).
#' @param TRSH Base(zero-)line threshold (default = 15).
#' @param RNGE Target range of responses (default = -100).
#' @param MXDV Maximum allowed deviation from monotonicity (default = 5).
#'
#' @param CARR Carryover detection threshold (default = 0, analysis skipped if set to 0)
#' @param BSFT For baseline shift issue, min.#points to detect baseline shift (default = 3, analysis skipped if set to 0).
#' @param USHP For u-shape curves, min.#points to avoid flattening (default = 4, analysis skipped if set to 0).
#' @param TrustHi For equal sets of corrections, trusts those retaining measurements at high concentrations (default = FALSE).
#' @param StrictImp It prevents extrapolating over concentration-range boundaries; used for POD, ECxx etc (default = TRUE).
#' @param DUMV A dummy value, default = -999.
#' @param TLOG A scaling factor for calculating the wAUC, default = -24.
#'
#' @return A list with corrected concentration-response measurements and several calculated curve metrics.
#' \itemize{
#'   \item resp: corrected responses
#'   \item corr: flags for corrections
#'   \item ECxx: effective concentration values at various thresholds
#'   \item Cxx: concentrations for various absolute response levels
#'   \item Emax: maximum effective concentration, slope of the mid-curve (b/w EC25 and EC75)
#'   \item wConc: response-weighted concentration
#'   \item wResp: concentration-weighed response
#'   \item POD: point-of-departure (first concentration with response >TRSH)
#'   \item AUC: area-under-curve (in units of log-concentration X response)
#'   \item wAUC: AUC weighted by concentration range and POD / TLOG (-24)
#'   \item wAUC_pre: AUC weighted by concentration range and POD
#'   \item nCorrected: number of points corrected (basically, sum of flags in corr)
#'   \item Comments: warning and notes about the dose-response curve
#'   \item Settings: input parameters for this run
#' }
#' @export
#' @references{
#'   \insertRef{PMID:20980217}{Rcurvep}\cr
#'
#'   \insertRef{PMID:27518631}{Rcurvep}
#' }
#' @seealso  [run_rcurvep()] and [combi_run_rcurvep()]
#' @examples
#'
#' curvep(Conc = c(-8, -7, -6, -5, -4) , Resp = c(0, -3, -5, -15, -30))
#'
#'


#[[Rccp::export]]
curvep <- function(Conc, Resp, Mask = NULL,
                    TRSH = 15, RNGE = -100, MXDV = 5, CARR = 0, BSFT = 3, USHP = 4,
                    TrustHi = FALSE, StrictImp = TRUE, DUMV = -999, TLOG = -24)
# Conc - array of concentrations, e.g., in Molar units, can be log-transformed, in which case internal log-transformation is skipped
# Resp - array of responses at corresponding concentrations, e.g., raw measurements or normalized to controls
# !!!NB: Conc & Resp arrays should be in order from lo to hi concentrations
# Mask - array of 1/0 flags indicating invalidated measurements (will be marked with DUMV value)
#
# TRSH - base(zero-)line threshold, RNGE - target range of responses, MXDV - maximum allowed deviation from monotonicity,
#
#-- parameters for advanced corrections of curve issues:
#   CARR - carryover detection threshold, (analysis skipped if set to 0)
#   BSFT - for baseline shift, min.#points to detect baseline shift, (analysis skipped if set to 0)
#   USHP - for u-shape curves, min.#points to avoid flattening, (analysis skipped if set to 0)
#--
#
# TrustHi - for equal sets of corrections, trusts those retaining measurements at high concentrations
# StrictImp - prevents extrapolating over concentration-range boundaries; used for POD, ECxx etc, it is used inside Impute()
#
#
#----------------
# curvep() returns a list with corrected dose-response measurements and several calculated curve metrics:
#
#     resp = corrected responses, corr = flags for corrections,
#     ECxx = effective concentration values at various thresholds,
#     Cxx  = concentrations for various absolute response levels,
#     Emax - maximum effective concentration, slope of the mid-curve (b/w EC25 and EC75),
#     wConc = response-weighted concentration, wResp = concentration-weighed response,
#     POD = point-of-departure (first concentration with response >TRSH),
#     AUC=area-under-curve (in units of log-concentration X response),
#     wAUC= AUC weighted by concentration range and POD / TLOG
#     nCorrected = number of points corrected (basically, sum of flags in corr)
#     Comments = warning and notes about the dose-response curve
#     Settings = input parameters for this run

# CurveP for R - 2017
#
#ToDo / History
#  30.04.2018 modified wAUC calculations for better scaling against diverse test ranges
#
#  10.18.2017 modified wAUC calculation to any scale of log-doses (could span negative and positive values both)
#             this now requires TLOG constant to represent the log10(dose) at infinite dilution or at the limit-of-detection
#             TLOG's default is -12 for assumed log10(M) scale) but could go as low as -24 (Avogadro number based)
#             NB: e.g., if user wishes to use milliMolar dose units as input throughout, the n TLOG should be likewise adjutsed by adding 6,
#             Additional output added (run settings)
#
#  10.16.2017 handling of non-numeric imput added (quits with message)
#
#  09.20.2017 curvecase of just above baseline: Conc <- c(-5.52 -5.00 -4.52 -4.00 -3.52), Resp <- c(0,5,0,0,5), if TRSH=MXDV=5, it will
#              result in no correction but significant response (because MXDV is not exceeded, and threshold is also satisfied)
#              Hence, TRSH should be always > MXDV. For TRSH <= MXDV, there may be non-monotonic points saved near threshold due to MXDV rescue.
#
#  07.24.2017 minor bug for slope estimation fixed
#             EC50 added to output as explicit metric
#
#  07.10.2017 minor bug fixed in Impute() which caused POD invalidation in some cases of carryover/inverse curves
#             15 additional cases of curves checked
#
#  07.05.2017 bugs fixed in get_monotonics() and in BLIP handling
#
#  07.03.2017 U-shape work check on 28099 case
#  06.16.2017 small cosmetic change in get_monotonics()
#
#  05.26.2017 testing on "tox21-er-luc-bg1-4e2-agonist-p2_luc.Rdata" curves
#
#  removed DUMV comparison for Conc, since no reading from file, instead the user must provide clean input
#
#
#
#
{
  if (!is.numeric(Conc))
  {
    print ("Non-numeric doses")
    return (NULL)
  }

  if (!is.numeric(Resp))
  {
    print ("Non-numeric response")
    return (NULL)
  }

	mxc <- max(Conc)
	mnc <- min(Conc)
	lgConc <- Conc
	if ( (mxc > 0) && (mnc > 0) )
	{
	  if ((mxc / mnc) > 10.0 )  lgConc <- log10(Conc)  #normalize
	}


  nCols <- length(lgConc)
	if (nCols != length(Resp)) return (NULL)


  #--- handle invalid points
	HTS <- Resp
	if (!is.null(Mask))
	{
	  if (length(Mask) == nCols)
	  {
	    Baddies <- (Mask == 1)
	    HTS[Baddies] <- DUMV
	  }
	}
	Baddies <- (HTS == DUMV)
	mxr <- max(0, RNGE)
	mnr <- min(0, RNGE)

	#--- pre-process curve
	for (c in 1:nCols)
	{
		if ( HTS[c] >  mxr)	HTS[c] <- mxr #may be not right for gene expression data
		if ( HTS[c] <  mnr)	HTS[c] <- mnr
		if ( abs( HTS[c] ) < TRSH )	HTS[c] <- 0
	}

	origBaddies <- Baddies

	#------- the curve should be monotonous, therefore, need to detect "violating" points
	CheckMore <- TRUE
	CheckYetMore <- TRUE

	tVals <- HTS[!Baddies]
	xRange <- 0

	Warn <- ""

	if ( is.null(tVals) || (length(tVals) == 0)  )
	{
		Warn <- "NO_VALID_POINTS"
		HTS <- rep(0, nCols)
		CheckMore <- FALSE
		CheckYetMore <- FALSE
	}
  else
  {
  	f = 1
  	while (HTS[f] == DUMV) f <- f + 1

  	c = nCols
  	while (HTS[c] == DUMV) c <- c - 1

  	xRange <- (HTS[f] - HTS[c])


  	if (abs(HTS[f]) > 0)
  	{ #treat blips
  		v <- f + 1
  		while (HTS[v] == DUMV) v <- v + 1
  		if ((v < nCols) && (HTS[v] == 0))
  		{
  		  v <- v + 1
  		  while (HTS[v] == DUMV) v <- v + 1
  		}

  		if ((v < nCols) && (HTS[v] == 0)) #erase starting point (likely blip)
  		{
  		  xRange <- xRange - HTS[f]
  		  HTS[f] <- 0
  		  Warn <- "BLIP"
  		  Baddies[f] <- TRUE
  		}
  	}

  	sdHTS <- sd(tVals)
  	if ( (abs(xRange) < TRSH) && ((sdHTS < MXDV) || (max(abs(HTS[f]), abs(HTS[c])) < TRSH)) )
  	{
  	  #redo the slope-direction analysis differently
  		xRange <- 0
  	}
  	else
  	{
  	  #check for possible partial U-curve
  	  c <- 0
  		for (v in 1:length(tVals))
  		{
  		  if ( abs(HTS[f] - tVals[v]) > (abs(xRange) + MXDV) ) c <- c + 1
  		}
  		if (c < 2) CheckMore <- FALSE     #no problems
  	}
  } #else.. if ( is.null(tVals) || (length(tVals) == 0)  )




	if ((sdHTS > MXDV) && CheckMore)
	{
		Ivals <- get_monotonics(HTS, MXDV, Baddies)

		#flat or U-shaped curves

		#pars for the optimal pivot, corrections, index, base
		errP <- nCols
		bP <- 1
		mP <- 0

    #--------------------------------------------------------
    #-------- find best pivot
    #--------------------------------------------------------
		for (c in 2:(nCols-1))
		{
			if (Ivals[c] == 0) next

			f <- sum(Ivals[1:c])
			v <- sum(Ivals[(c+1):nCols])

			cP = abs(v - f)      #reflects the width of the spike in count of non-flat intervals
			if (cP < mP) next

			vP <- c   #determines pivot
			xP <- 0   #counts corrections (skips continuously monotonic points on the slope of U-shape)

			if ( RNGE*(f - v) > 0 )
			{ #hi-lo-hi u-shapes - normally, this case should not happen! So...
			  #treat as spurious, unless wo parts of the curve move clearly in different directions
				if (f*v >= 0) next

				flat <- FALSE
				tf <- vP
				while (tf > 1)
				{
					tf <- tf - 1
					if (Baddies[tf]) next

					if (Ivals[tf] == Ivals[c])
					{
					  flat =FALSE
					  next
					}

					if ( (Ivals[tf] * Ivals[c]) < 0 )
					{
					   xP <- xP + 1
					   break
					} #counter-slop, stop

					if (flat)
					{
					  xP <- xP + 2
					  break
					} #two flat regions in a row, stop

					flat = TRUE
				} #while (tf > 1)

				xP <- xP + tf - 1

				while (tf > 1) #remove invalidated points from error count
				{
				  tf <- tf - 1
				  if (Baddies[tf]) xP <- xP - 1
				}
			}
			else
			{ #lo-hi-lo u-shapes
				#move the pivot to the last plateau point
				while (vP < nCols - 1)
				{
				  if (Ivals[vP + 1] == 0)
				  {
				    vP <- vP + 1
				  }
				  else break
				}

				tf <- vP
				while (tf < nCols)
				{
				  tf <- tf + 1
					if (Baddies[tf]) next
					if (HTS[tf] == 0) #to count baseline points as errors (lo area), not the slope
					{
					  tf <- tf - 1 #set one back for proper count of errors
					  break
					}
				}
				xP <- nCols - tf
			}


			#--- compare pivot with current best
			if (errP + cP < xP + mP) next    #compares errors / support
			mP   <- cP
			bP   <- vP
			errP <- xP
		} #for c
    #--------------------------------------------------------
    #------------ end of best pivot detection
    #--------------------------------------------------------

		#--- U-shape pivot point was found, now check if this is a correctable case
		xP <- 0
		vP <- 0
		v <- 0
		f <- 1
		for (c in 1:nCols)
		{
			if (Ivals[c] == 0)	next
			if (Ivals[c]*Ivals[f] < 0) v <- v + 1      #count spikes (changes in monotonicity)
			if (c > bP) vP <- vP + Ivals[c] else xP <- xP + Ivals[c]
			f <- c #save last non-flat interval for in-cycle spike detection
		}

		xWrk <- abs(HTS[bP])                 #may use high absolute response to rescue U-shape in some cases
		if (CARR == 0)  xWrk <- -1
		if (Ivals[bP] == 0) mP <- mP + 1    #if pivot it at a plateau, add support to base
		if (mP < 3) xWrk <- -1              #narrow base, true spike, disable potency-based rescue

		errP <- errP - 2*(mP - USHP)   #allows 2 addl corrections per extra support

		if ( (bP == 1) || ((xWrk < CARR) && (max(v, errP) > mP)) || (USHP == 0) )
		{

			if (xRange == 0) Warn <- paste(Warn, "NOISY") else Warn <- paste(Warn, "PART_U?")
			if (xWrk > CARR) Warn <- paste(Warn, "POTENT")
		}
		else
		{#treat as U-shape, its pivot represents peak/plateau
			Warn <- paste(Warn, "U_SHAPE")
			if (errP > mP) Warn <- paste(Warn, "CHECK")
			if (RNGE > 0) xRange <- -1 else xRange <- 1

			#--- invalidate the wrong part to aid outlier-detection
			if ( RNGE*(xP - vP) > 0 )
			{#hi-lo-hi
				Warn <- paste(Warn, "CARRY_OVER")
				Baddies[1:bP] <- TRUE
				HTS[bP] = 0      #important seed for corrections-to-baseline
			}
			else Baddies[(bP+1):nCols] <- TRUE      #lo-hi-lo
		}
	} #if ((sdHTS > MXDV) && CheckMore)
	else
	{
		if ((xRange == 0) && CheckMore)
		{
			avv = mean(tVals)
			if (abs(avv) >= TRSH)
			{#const curve with low variance and signif signal: can be baseline shift or carry over, do not erase
			  nonc <- abs(HTS - avv) > MXDV
			  Baddies[nonc] <- TRUE
			  HTS[nonc] <- avv
			  CheckYetMore <- FALSE
			}
		}
	}

  #Now xRange stores a range of the curve, if < 0 then it's rising
	if ((xRange == 0) && CheckYetMore)
	{#if still flat, make it so
	  nonf <- abs(HTS) > 0
	  Baddies[nonf] <- TRUE
	  HTS[nonf] <- 0
	  CheckYetMore <- FALSE
	}

	#------------------------------------------------------------------
	#---------------  Detect a minimum set of violating points --------
	#------------------------------------------------------------------
  if ((xRange != 0) && CheckYetMore)
	{#Detect a minimum set of violating points

		TrialBest <- rep(FALSE, nCols)
		tbSize <- nCols
		bdSize <- length(HTS[Baddies])

		for (v in 1:nCols)
		{ #v seeds the starting point to trust as "true"
			if (Baddies[v]) next

			Trial <- Baddies

			#detecting discrepancies from v onward
			f <- v
			c <- f
			while (c < nCols)
			{
			  c <- c + 1
				if (Baddies[c]) next
				tdff <- HTS[c] - HTS[f]
				if ((xRange*tdff > 0) && (abs(tdff) > MXDV)) Trial[c] <- TRUE else f <- c
			}

			#detecting discrepancies on the curve from v backward
			f <- v
			c <- f
			while (c > 1)
			{
			  c <- c - 1
				if (Baddies[c]) next
				tdff <- HTS[f] - HTS[c]
				if ((xRange*tdff > 0) && (abs(tdff) > MXDV)) Trial[c] <- TRUE else f <- c
			}

			f = length(HTS[Trial])
			if (tbSize < f) next

			if ( (tbSize > f) || TrustHi)
			{   #update best set of corrections
			    TrialBest <- Trial
			    tbSize <- f
			}

			if (!TrustHi)	if (tbSize == bdSize) break #the optimum is reached
		}#for v

		Baddies <- TrialBest	#output mask
	} #if ((xRange != 0) && CheckYetMore)
	#------------------------------------------------------------------
	#------------------------------------------------------------------


	#------------------------------------------------------------------
	#------- replace bad points with appropriate extrapolations -------
	#------------------------------------------------------------------
  if (CheckYetMore)
  {
		for (c in 1:nCols)
		{
			if (!Baddies[c]) next

			f <- c
			v <- c

			while (v > 0) if (Baddies[v]) v<-v-1 else break
			while (f <= nCols) if (Baddies[f]) f<-f+1 else break

			if (v == 0)
			{
				if (f > nCols)
				{
				  Warn <- paste(Warn, "NO_VALID_POINTS")
				  HTS[c] <- 0
				}	#all points are invalidated, set to baseline
			  else HTS[c] <- HTS[f]
			}
			else
			{
				if (f > nCols) HTS[c] <- HTS[v]
				else
				{
					HTS[c]	 <- (lgConc[c] - lgConc[v])/(lgConc[f] - lgConc[v])
					HTS[c]	<- HTS[c] * (HTS[f] - HTS[v]) + HTS[v]
				}
			}
		}#for c
  } #if (CheckYetMore)

  f <- 0
  v <- 0
	for (c in 1:nCols)
	{#counts non-zero signals, original and imputed
		if (HTS[c] == 0) next
		if (HTS[c] == DUMV) next
		f <- f + 1
		if (Baddies[c]) next
		v <- v + 1
	}

	if ((v == 1) && (f > 1)) Warn <- paste(Warn, "SINGLE_POINT_ACT")

	tdff <- xRange * RNGE
	if (tdff > 0) Warn <- paste(Warn, "INVERSE")



	##------- detection & handling of carry-over and related problems
	if ( (CARR > 0) && (v > 0) )
	{
		c <- 1
		f <- nCols

		xWrk = abs(HTS[c])
		cuRange = abs(HTS[c] - HTS[f])
		if (xWrk > 0)
		{
				if ( (tdff > 0) && (cuRange >= TRSH) )
				{ #decrease in signal, unconditional carryover (if inhibitor) or potency-conditional (if agonist)
					if ( (RNGE < 0) || (xWrk < CARR) )
					{
					  Warn <- paste(Warn, "CARRY_OVER")
						for (z in c:f)
						{
						    if (HTS[z] == 0) break
						    Baddies[z] <- TRUE
						    HTS[z] <- 0
						}
					} else Warn <- paste(Warn, "CHECK")
				}
				else #increasing or near-constant signal
					if (xWrk < CARR)
					{
						if ( (xRange == 0) || (cuRange < TRSH) )
						{#constant or nearly so, likely carryover or baseline shift
							Warn <- paste(Warn, "CARRY_OVER? BASE_SHIFT?")
							for (z in c:f)
							{
							  if (HTS[z] == 0) break
							  Baddies[z] <- TRUE
							  HTS[z] <- 0
							}
						}
						else
						{ #increasing, can be potent active, carry over or baseline shift
							vv <- HTS[c]
							nvv <- 1
							#baseline detection
							for (f in (c+1):nCols)
							{
								vv <- c(vv, HTS[f])
								if (sd(vv) < MXDV) nvv <- length(vv)
								if (nvv + 3 < length(vv)) break #stop scanning after several violations in a row
							}

							if (nvv < BSFT)
								Warn <- paste(Warn, "CHECK")
							else
							{#baseline shift, adjust
								Warn <- paste(Warn, "BASE_SHIFT?")
								#xWrk <- mean(vv)
								xWrk <- median(vv[1:nvv])
								for (f in 1:nCols)
								{
									nvv <- nvv - 1
									if (BSFT > 0)
									{
										HTS[f] <- HTS[f] - xWrk		#fix as base-shift
										if ( abs(HTS[f]) < TRSH )	HTS[f] <- 0
									}
									else
									{
										HTS[f] <- 0
										Baddies[f] <- TRUE      #fix as carry-over
										if (nvv == 0) break
									}
								}
							}
						}
					} #if (xWrk < CARR)
					else Warn <- paste(Warn, "TOO_POTENT")
		} #if (xWrk > 0)
	}#if ( (CARR > 0) && (v > 0) )


	##------------------------------------------------------------
	##------- output updated mask, and adjusted data -------------
	  rlevels <- c (1, 5, 10, 20, 25, 50, 75, 80, 90, 95, 99, 100)
		Corrections <- rep(0, nCols)
		Corrections[Baddies] <- 1

		Emax <- 0
		POD <- DUMV
		slope <- 0
		CAs <- rep(DUMV, length(rlevels))
		ECs <- rep(DUMV, length(rlevels))

		##-------- print C@.. and EC..
		if (RNGE > 0) Emax <- max(HTS) else Emax <- min(HTS)

		wConc <- DUMV
		wResp <- 0

		if (abs(Emax) > 0)
		{#some signal
			slope <- Emax/100
			xWrk <- abs(slope)
			for(c in 1:length(rlevels))
			{
				CAs[c] <- Impute(lgConc, HTS, rlevels[c], DUMV, StrictImp)          #NB: estimated C for a given level of response (from 1 to 99)
				ECs[c] <- Impute(lgConc, HTS,  xWrk*rlevels[c], DUMV, StrictImp)    #NB: ECs are tied to Emax
			}
			POD <- Impute(lgConc, HTS, TRSH, DUMV, StrictImp)
			slope <- slope *50/(ECs[7] - ECs[5]) #corr.07.24.17

			pdr <- sum(lgConc * HTS)
			wConc <- pdr / sum(HTS)
			wResp <- pdr / sum(lgConc)	#only reasonable if log10(M) units are used that are always the same sign
		}

		#--- area under curve calculation
		AUC <- 0
		dlgC <- lgConc[2:nCols] - lgConc[1:(nCols-1)]
		sR <- HTS[2:nCols] + HTS[1:(nCols-1)]
		AUC <- sum(sR * dlgC) / 2

		wAUC <- 0
		wAUC_old <- 0

		if (POD != DUMV)
		{
		  wAUC_old <- AUC * POD/(lgConc[1] - lgConc[nCols])
		  #old equation above assumes log10(M) units, and thus is scale-dependent (e..g, switching sign for positive log-doses)

		  #30.04.2018 --
		  wAUC <- AUC * (lgConc[1] - TLOG) / (POD - TLOG) / (lgConc[nCols] - TLOG)

		  #10.18.2017 --
		  #wAUC <- AUC / (POD - TLOG) / (lgConc[nCols] - lgConc[1])
		}


		#print AUC, wAUC, POD, Warn

		if (nchar(Warn) == 0) Warn <- "OK"
		Warn <- gsub(" ", "|", Warn)
		Warn <- sub("^\\|", "", Warn) #07.03.17

		n <- length(HTS[Baddies]) - length(HTS[origBaddies])

		#gather settings in one list:
		sett <- list(TRSH = TRSH, RNGE = RNGE, MXDV = MXDV, CARR=CARR, BSFT=BSFT, USHP=USHP, TrustHi = TrustHi, StrictImp = StrictImp, DUMV = DUMV, TLOG = TLOG)

		return
		(

		  list(resp = HTS, corr = Corrections, #levels = data.frame(xx = rlevels, ECxx =ECs, Cxx =CAs, row.names = paste(rlevels, "%", sep='')),
		       ECxx = ECs, Cxx = CAs, xx = rlevels,
		       Emax = Emax, slope = slope, wConc = wConc, wResp = wResp, EC50 = ECs[which(rlevels == 50)], C50 = CAs[which(rlevels == 50)], POD = POD, AUC=AUC, wAUC=wAUC, wAUC_prev=wAUC_old, nCorrected = n,
		       Comments = Warn, Settings = sett)
		)

} #curvep(...)















#--------------------------------------------------------------------------
#debug area
# load("C:\\DATA\\BMDExpress\\tox21-er-luc-bg1-4e2-agonist-p2_luc.Rdata")
# lrows <- 1:dim(cebs)[1]
#
# totstuff <- NULL
#
# special <- !(cebs$curvep_remark == "OK")
# length(cebs$curvep_remark[special])
#
# special1 <- (cebs$curvep_remark == "U_SHAPE")
# i <- 28099
# for (i in lrows)
# #for (i in lrows[special1])
# {
#   C <- cebs[i, 14:28]
#   R <- cebs[i, 29:43]
#   valid <- !is.na(C) #there are some NAs to handle
#   results <- curvep(C[valid], R[valid], TRSH = 25, RNGE = 1e+006, CARR = 60, TLOG = -24)
#   #totstuff <- rbind(totstuff, results$resp)
#   totstuff <- rbind(totstuff, c(results$AUC, results$wResp, results$wConc, results$EC50, results$wAUC, results$wAUC_prev) )
# }
#
# colnames(totstuff) <- c("AUC", "wResp", "wConc", "EC50", "wAUC", "wAUC_prev")
# write.table(totstuff, file="C:\\DATA\\BMDExpress\\curvep_tlog24_tox21-er-luc-bg1-4e2-agonist-p2_luc_101817.txt", row.names = FALSE, sep = "\t")
# #points(C[valid], results$resp, col = "red")
# #dbi <- (cebs$uniqueID == "N10168")
# #dbi <- (cebs$uniqueID == "N23531")
# #dbi <- (cebs$uniqueID == "N26478")
# #dbi <- (cebs$uniqueID == "N29722")
# #dbi <- (cebs$uniqueID == "N14827")
# #iter2
# #dbi <- (cebs$uniqueID == "N17409")
# #dbi <- (cebs$uniqueID == "N23332")
# #dbi <- (cebs$uniqueID == "N34085")
# #dbi <- (cebs$uniqueID == "N34757")
# #dbi <- (cebs$uniqueID == "N3906")
# #
# #dbi <- (cebs$uniqueID == "N1133")
# #dbi <- (cebs$uniqueID == "N12666")
# #dbi <- (cebs$uniqueID == "N12698")
# #dbi <- (cebs$uniqueID == "N14476")
# dbi <- (cebs$uniqueID == "N14477")
#
# cebs_dbi <- cebs[dbi,]
# C <- cebs_dbi[14:28]
# R <- cebs_dbi[29:43]
# valid <- !is.na(C)
# plot(unlist(C[valid]), unlist(R[valid]), ylim = c(0, 100))
# points(unlist(C[valid]), results$resp, col="red")
# #debug init
# Conc <- C[valid]
# Resp <- R[valid]
# Mask <- NULL
# TRSH <- 25
# RNGE <- 1000000
# MXDV <- 5
# CARR <- 60
# #CARR <- 80
# BSFT <- 3
# USHP <- 4
# TrustHi <- FALSE
# StrictImp <- TRUE
# DUMV <- -999
