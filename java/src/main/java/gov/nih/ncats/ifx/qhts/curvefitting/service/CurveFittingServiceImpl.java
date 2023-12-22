package gov.nih.ncats.ifx.qhts.curvefitting.service;
  
import java.util.*;
import gov.nih.ncats.ifx.qhts.utils.*;
import gov.nih.ncats.ifx.qhts.curvefitting.FittingResult;
import gov.nih.ncats.ifx.qhts.curvefitting.algorithm.*;

import org.apache.commons.lang3.ArrayUtils; 
  
public class CurveFittingServiceImpl implements CurveFittingService{	
	 

	/**
	 * Use default fitting constants in HillConstants.java
	 */
	@Override
	public FittingResult doFit(Double[] conc, Double[] resp) { 
		
		return doFit(conc, resp, null, null);
	}
	
	
	/**
	 * Use customized sdPct & minYRange
	 */
	@Override
	public FittingResult doFit(Double[] conc, Double[] resp, Double classificationSD, Double minYRange) {
		 
		FittingResult fitResult = new FittingResult();
		//1.
		Optional<Double> classificationSDOptional = Optional.ofNullable(classificationSD);
		Optional<Double> minYRangeOptional = Optional.ofNullable(minYRange);
		
		//2.
		double[] concs = toLog10(conc);
    	double[] resps = ArrayUtils.toPrimitive(resp);     	 
		
    	//3.
		Object[] fitResults = _doFit(concs, resps, minYRangeOptional);
		
		if(fitResults == null || (double[])fitResults[1] == null) {
			fitResult.setCurveClass2(4d);
			return fitResult;
		}
		 	
		String flagStr = AlgotithmUtils.composeMaskFlag((boolean[])fitResults[0]);
				
		//4.
		boolean[] flagBoolean = (boolean[])fitResults[0]; 
		
		int numOfMask = 0;
		for(boolean f : flagBoolean) 
			numOfMask += !f?1:0;//count the number of the 'false'
		
		//{ x05_min, yinf_min, slope_min, r2, dev_const, dev_min, y0_min };
		double[] fitValues = (double[])fitResults[1];
		 
		
		double ac50 = Math.pow(10, fitValues[0]+6.0);
		fitResult.setAc50(ac50);
		
		fitResult.setLogAc50(fitValues[0]);
		fitResult.setInfActivity(fitValues[1]);
		fitResult.setHillCoef(fitValues[2]);
		fitResult.setR2(fitValues[3]);		
		fitResult.setZeroActivity(fitValues[6]);
		
		fitResult.setMaskFlag(flagStr);
		fitResult.setMaskNo(numOfMask);
		
		double maxResp = resp[resp.length-1];
		fitResult.setMaxResponse(maxResp);
		
		double efficacy = fitValues[1]-fitValues[6];
		fitResult.setEfficacy(efficacy);
		 
		//5.
		double y0 = fitValues[6];
		double x05 = fitValues[0];
		double yinf = fitValues[1];
		double slope = fitValues[2]; 

        CurveClass2 curveClass2 = new CurveClass2();
        
        if(classificationSDOptional.isPresent())
        	HillConstants.CLASSIFICATION_SD = classificationSDOptional.get();
        	 
        curveClass2.setSDPct(HillConstants.CLASSIFICATION_SD);
        curveClass2.setSdFactor(HillConstants.CLASSIFICATION_SD_FACTOR);
        curveClass2.setRelativeActivity(HillConstants.RELATIVE_ACTIVITY);
        

		double pValue = HillStat.calcPValue(y0, yinf, x05, slope, concs, resps, flagBoolean);
		
		try {
			double classNo2 = curveClass2.calcCurveClass(fitValues, concs, resps, flagBoolean, pValue);	
			
			fitResult.setpValue(pValue); 
			fitResult.setCurveClass2(classNo2);
			
		} catch (Exception e) { 
			e.printStackTrace();
		}
		
		return fitResult;
	}
	
	
	/**
	 * calculate the curveClass2 only
	 */
	@Override
	public double curveClass2(Double[] conc, Double[] resp, Double sdPct, Double minYRange) {
		
		Optional<Double> sdPctOptional = Optional.ofNullable(sdPct);		
		Optional<Double> minYRangeOptional = Optional.ofNullable(minYRange);
		
		return curveClass2(conc, resp, sdPctOptional, minYRangeOptional);
	}
	
	 
	/**
	 * sdPct --- used by Tox21, sdPct=5.0
	 * default --- sdPct = 10.0 
	 */
	@Override
	public double curveClass2(Double[] conc, Double[] resp, Optional<Double> sdPctOptional, Optional<Double> minYRangeOptional) {
		 
		double[] concs = toLog10(conc);
    	double[] resps = ArrayUtils.toPrimitive(resp); 
    	 
		Object[] fitResults = _doFit(concs, resps, minYRangeOptional);
		
		double classNo2 = 4d;
		if (fitResults == null)
			return classNo2;		
	 
        double[] fitValues = (double[])fitResults[1];
		if (fitValues == null)
			return classNo2;
 
		boolean[] flags = (boolean[])fitResults[0];
		
		return _curveClass2(concs, resps, fitValues, flags, sdPctOptional);
	}

	
	/**
	 * do fitting
	 * 
	 * minYRangeOptional default is 20.0
	 */
	private static Object[] _doFit(double[] concs, double[] resps, Optional<Double> minYRangeOptional) {
		
		if(minYRangeOptional.isPresent())
			HillConstants.MIN_Y_RANGE = minYRangeOptional.get();
		
		try { 
			Object[] fitResults = HillFit.doHill(concs, resps);
			return fitResults;
			
		} catch (Exception e) {
			e.printStackTrace();
		}
	
    	return null;
	}
	

	/**
	 * generate curveClass2
	 */
	private static double _curveClass2(double[] concs,double[] resps, double[] fitValues, boolean[] flags, Optional<Double> sdPctOptional) {
		
		//1.
        double p = _pValue(concs, resps, fitValues, flags);
        
        double classNo2 = 4d;
        
        //2.
        CurveClass2 curveClass2 = new CurveClass2();
        
        //config
        if(sdPctOptional.isPresent())
        	HillConstants.CLASSIFICATION_SD = sdPctOptional.get();
        	 
        curveClass2.setSDPct(HillConstants.CLASSIFICATION_SD);
        curveClass2.setSdFactor(HillConstants.CLASSIFICATION_SD_FACTOR);
        curveClass2.setRelativeActivity(HillConstants.RELATIVE_ACTIVITY);
        
        //3.
		try {
			classNo2 = curveClass2.calcCurveClass(fitValues, concs, resps, flags, p);			
			return classNo2;
			
		} catch (Exception e) { 
			e.printStackTrace();
		}
		
		return classNo2;
	}

	
	/**
	 * calculate the p value
	 */
	private static double _pValue(double[] concs,double[] resps, double[] fitValues, boolean[] flags) {
		
		double y0 = fitValues[6];
		double x05 = fitValues[0];
		double yinf = fitValues[1];
		double slope = fitValues[2]; 
		
		return HillStat.calcPValue(y0, yinf, x05, slope, concs, resps, flags);
	}

	
	/**
	 * convert to log10
	 */
	private static double[] toLog10(Double[] conc) { 
		
		//get log10 of concentrations		
		//return Arrays.asList(conc).stream().sequential().mapToDouble(c->Math.log10(c)).toArray(); 
		
		Double[] temp = new Double[conc.length];
		for(int i=0;i<conc.length;i++)
			temp[i] = Math.log10(conc[i]);
		
		return ArrayUtils.toPrimitive(temp);
	}
 	
}
