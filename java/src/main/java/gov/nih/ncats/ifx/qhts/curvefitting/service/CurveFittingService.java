package gov.nih.ncats.ifx.qhts.curvefitting.service;

import java.util.Optional;

import gov.nih.ncats.ifx.qhts.curvefitting.FittingResult;

public interface CurveFittingService {

	
	public FittingResult doFit(Double[] conc, Double[] resp);
	
	public FittingResult doFit(Double[] conc, Double[] resp, Double sdPct, Double minYRange);
	
	public double curveClass2(Double[] conc, Double[] resp, Double sdPct, Double minYRange);
	
	public double curveClass2(Double[] conc, Double[] resp, Optional<Double> sdPctOptional, Optional<Double> minYRangeOptional);
	
}
