package gov.nih.ncats.ifx.qhts.curvefitting.algorithm;


import java.util.*;
import gov.nih.ncats.ifx.qhts.utils.*;

import static gov.nih.ncats.ifx.qhts.curvefitting.algorithm.AlgotithmUtils.*;
import static gov.nih.ncats.ifx.qhts.utils.HillConstants.*;

public class HillFit {

	public static double LN10 = Math.log(10.0);
	protected static boolean[] maskings_ = null;

	private static boolean fastFlag = false;
	private static boolean slopeFlag = true;
	private static boolean p4Fit = true;

	
	public static void setFastFlag(boolean flag) {
		fastFlag = flag;
	}

	
	public static void setSlopeFlag(boolean flag) {
		slopeFlag = flag;
	}
	

	public static Object[] doHill(double[] x, double[] y) throws Exception{
		return doHill(x, y, null, PARTIAL_FIT_MASK_NO, null, 1, CLASSIFICATION_SD, CLASSIFICATION_SD_FACTOR);
	}

	
	public static Object[] doHill(double[] x, double[] y, boolean flag) throws Exception{
		p4Fit = flag;
		return doHill(x, y);
	}

	
	public static Object[] doHill(double[] x, double[] y, String maskFlag, int iterNo, String fixPara) throws Exception{
		return doHill(x, y, maskFlag, iterNo, fixPara, 1, CLASSIFICATION_SD, CLASSIFICATION_SD_FACTOR);
	}

	
	public static Object[] doHill(double[] conc, double[] y, String maskFlag, int iterNo, String fixPara, int dn, double sd, double sdf) throws Exception {

		int numPoints = conc.length;
		
		if (numPoints < 3)
			throw new Exception("Concentratrion points should be equl or greater than 3");

		double bestR2 = 0.0;
		HashMap[] bestMaps = null;
		double[] bestFitValues = null;
		boolean[] bestFlags = null;
		double bestDelta = 0.0;
		int bestIterNo = 0;

		boolean[] maskFlags = parseMaskFlag(maskFlag);

		if (maskFlags != null || MASK_FLAG == false)
			iterNo = 1;

		for (int i = 0; i < iterNo; i++) {
			double[] xs = new double[numPoints - i];
			double[] ys = new double[xs.length];

			for (int j = 0; j < xs.length; j++) {
				xs[j] = conc[j];
				ys[j] = y[j];
			}

			double ys0 = ys[0];

			if (SHIFT_ACTIVITY) {
				for (int j = 0; j < xs.length; j++) 
					ys[j] -= ys0;
			}

			HashMap[] maps = new HashMap[xs.length];

			boolean noChangeFlag = false;
			boolean[] flags = null;

			if (maskFlags != null) {
				flags = maskFlags;
				noChangeFlag = true;
			} else {
				flags = Masking.maskDif(xs, ys, null, null, maps, dn, sd, sdf);
			}

			if (xs.length < 3)
				break;

			// if bell shape, no change
			if (BELL_MASK)
				noChangeFlag = true;

			double[] fitValues = null;

			if (FIT_TYPE.equals("P5"))
				fitValues = hillFit5P(PI_MAX, PS_MIN, flags, null, xs, ys, noChangeFlag, fixPara);
			else
				fitValues = hillFit(PI_MAX, PS_MIN, flags, null, xs, ys, noChangeFlag, fixPara);

			if (fitValues == null && bestFlags == null)
				bestFlags = flags;

			if (fitValues == null)
				continue;

			if (SHIFT_ACTIVITY) {
				fitValues[1] += ys0;
				fitValues[6] += ys0;
			}

			double r2 = fitValues[3];
			double delta = Math.abs(fitValues[1] - fitValues[6]);

			if (r2 * xs.length > bestR2 && delta > bestDelta) {
				bestR2 = r2 * xs.length;
				bestFitValues = fitValues;
				bestFlags = flags;
				bestMaps = maps;
				bestDelta = delta;
				bestIterNo = i;
 
				//System.out.println(i + "\t" + xs.length + "\t" + r2 + "\t" + bestR2 + "\t" + delta + "\t" + bestDelta);
				
			} else {
				//System.out.println(i + "\t" + xs.length + "\t" + r2 + "\t" + bestR2 + "\t" + delta + "\t" + bestDelta);
				break;
			}
		}

		if (bestFitValues == null) {
			Object[] values = { bestFlags, bestFitValues, bestMaps, String.valueOf(bestIterNo) };

			return values;
		}

		
		boolean[] flags = new boolean[numPoints];
		
		for (int i = 0; i < bestFlags.length; i++)
			flags[i] = bestFlags[i];

		for (int i = bestFlags.length; i < numPoints; i++)
			flags[i] = false;

		HashMap[] maps = new HashMap[numPoints];
		for (int i = 0; i < bestFlags.length; i++)
			maps[i] = bestMaps[i];

		for (int i = bestFlags.length; i < numPoints; i++)
			maps[i] = null;

		Object[] values = { flags, bestFitValues, maps, String.valueOf(bestIterNo) };

		return values;

	}
	

	public static double[] hillFit5p(double ymin, double ymax, boolean[] flags, double[] ws, double[] xs, double[] ys, String fixPara, boolean upFlag) {
		if (checkRange(ymin, ymax, flags, ys) == false)
			return null;

		double xl = -10.0;
		if (xs[0] * 1.2 < xl)
			xl = xs[0] * 1.2;

		double xr = -2.0;
		double xd = 0.5;

		double y0l = -150.0;
		double y0r = 150.0;
		double y0d = 5.0;

		double yl = -150.0;
		double yr = 150.0;
		double yd = 5.0;

		double sl = MIN_SLOPE;
		double sr = MAX_SLOPE;
		double sd = 0.2;

		double ml = 0.1;
		double mr = 10.0;
		double md = 1.0;

		double max_y = -Double.MAX_VALUE;
		double min_y = Double.MAX_VALUE;

		for (int i = 0; i < ys.length; i++)
			if (flags[i]) {
				if (ys[i] > max_y)
					max_y = ys[i];

				if (ys[i] < min_y)
					min_y = ys[i];
			}

		if (upFlag) {
			if (min_y < 0) {
				y0l = Y0_INF_COEF * min_y;
				y0r = 0.0;
			} else {
				y0l = 0.0;
				y0r = Y0_INF_COEF * min_y;
			}

			if (max_y > 0.0) {
				yl = 0.0;
				yr = Y0_INF_COEF * max_y;
			} else {
				yr = 0.0;
				yl = Y0_INF_COEF * max_y;
			}

			sl = MIN_SLOPE;
			sr = MAX_SLOPE;
		} else {
			if (max_y < 0) {
				y0l = Y0_INF_COEF * max_y;
				y0r = 0.0;
			} else {
				y0l = 0.0;
				y0r = Y0_INF_COEF * max_y;
			}

			if (min_y < 0.0) {
				yl = Y0_INF_COEF * min_y;
				yr = 0;
			} else {
				yr = Y0_INF_COEF * min_y;
				yl = 0;
			}

			sl = MIN_SLOPE;
			sr = MAX_SLOPE;
		}

		// y0d = delta_y;
		// yd = delta_y;

		double delta_y = 0.05 * (max_y - min_y);

		y0d = 0.05 * (y0r - y0l);
		y0d = (y0d < 2.0) ? 2.0 : y0d;
		y0d = (y0d < delta_y) ? delta_y : y0d;

		yd = 0.05 * (yr - yl);
		yd = (yd < 2.0) ? 2.0 : yd;
		yd = (yd < delta_y) ? delta_y : yd;

		if (delta_y < 1.0 && y0d > 1.0)
			y0d = 1.0;

		if (delta_y < 1.0 && yd > 1.0)
			yd = 1.0;
		//System.out.println(y0l + ":" + y0l + ":" + yl + ":" + yr + ":" + xl + ":" + xr);

		// System.out.println(p4Fit+":"+min_y+":"+max_y+":"+delta_y);
		double[] fitValues = null;

		fitValues = hillFit(y0l, y0r, y0d, yl, yr, yd, xl, xr, xd, sl, sr, sd, flags, ws, xs, ys, fixPara);

		if (fitValues == null)
			return null;

		xl = fitValues[0] - 1.0;
		xr = fitValues[0] + 1.0;
		xd = EC50_GRID_SIZE;

		y0l = fitValues[6] - 2 * y0d;
		y0r = fitValues[6] + 2 * y0d;

		yl = fitValues[1] - 2 * yd;
		yr = fitValues[1] + 2 * yd;

//            delta_y = 0.1*delta_y;
//            if (delta_y < 0.5)
//                delta_y = 0.5;
//
//            y0d = delta_y;
//            yd = delta_y;

		y0d = 0.2 * y0d;
		y0d = (y0d < 1.0) ? 1.0 : y0d;

		yd = 0.2 * yd;
		yd = (yd < 1.0) ? 1.0 : yd;

		if (fitValues[2] - 0.5 < 0.1)
			sl = 0.1;
		else
			sl = fitValues[2] - 0.4;

		sr = fitValues[2] + 0.4;
		sd = 0.2;

		ml = 0.4;
		mr = 2.0;
		md = 0.2;

		fitValues = hillFit5p(y0l, y0r, y0d, yl, yr, yd, xl, xr, xd, sl, sr, sd, ml, mr, md, flags, ws, xs, ys,fixPara);

		return fitValues;
	}
	


    public static double[] hillFitFast5P(double ymin, double ymax,boolean[] flags, double[] ws, double[] xs,double[] ys, String fixPara, boolean upFlag)
    {
            if (checkRange(ymin,ymax,flags,ys) == false)
                    return null;

            double xl = -10.0;
            if (xs[0]*1.2 < xl)
                xl = xs[0]*1.2;

            double xr = -2.0;
            double xd = 0.5;

            double y0l = -150.0;
            double y0r = 150.0;
            double y0d = 5.0;

            double yl = -150.0;
            double yr = 150.0;
            double yd = 5.0;

            double sl = MIN_SLOPE;
            double sr = MAX_SLOPE;
            double sd = 0.2;

            double ml = 0.1;
            double mr = 10.0;
            double md = 1.0;

            double max_y = -Double.MAX_VALUE;
            double min_y = Double.MAX_VALUE;

            for (int i=0; i<ys.length; i++)
            if (flags[i])
            {
                if (ys[i] > max_y)
                    max_y = ys[i];

                if (ys[i] < min_y)
                    min_y = ys[i];
            }

            if (upFlag)
            {
                    if (min_y<0)
                    {
                        y0l = Y0_INF_COEF*min_y;
                        y0r = 0.0;
                    }
                    else
                    {
                        y0l = 0.0;
                        y0r = Y0_INF_COEF*min_y;
                    }

                    if (max_y > 0.0)
                    {
                        yl = 0.0;
                        yr = Y0_INF_COEF*max_y;
                    }
                    else
                    {
                        yr = 0.0;
                        yl = Y0_INF_COEF*max_y;
                    }


                    sl = MIN_SLOPE;
                    sr = MAX_SLOPE;
            }
            else
            {
                    if (max_y<0)
                    {
                        y0l = Y0_INF_COEF*max_y;
                        y0r = 0.0;
                   }
                    else
                    {
                        y0l = 0.0;
                        y0r = Y0_INF_COEF*max_y;
                   }

                    if (min_y < 0.0)
                    {
                        yl = Y0_INF_COEF*min_y;
                        yr = 0;
                    }
                    else
                    {
                        yr = Y0_INF_COEF*min_y;
                        yl = 0;
                   }

                    sl = MIN_SLOPE;
                    sr = MAX_SLOPE;
            }


            //y0d = delta_y;
            //yd = delta_y;

            double delta_y = 0.05*(max_y-min_y);

            y0d = 0.05*(y0r-y0l);
            y0d = (y0d < 2.0)?2.0:y0d;
            y0d = (y0d < delta_y)?delta_y:y0d;

            yd = 0.05*(yr-yl);
            yd = (yd < 2.0)?2.0:yd;
            yd = (yd < delta_y)?delta_y:yd;

            if (delta_y < 1.0 && y0d > 1.0)
                y0d = 1.0;

            if (delta_y < 1.0 && yd > 1.0)
                yd = 1.0;

            //System.out.println(y0l+":"+y0l+":"+yl+":"+yr+":"+xl+":"+xr);

            //System.out.println(p4Fit+":"+min_y+":"+max_y+":"+delta_y);
            double[] fitValues = null;


            fitValues = hillFit(y0l,y0r,y0d,yl,yr,yd,xl,xr,xd,sl,sr,sd,flags,ws,xs,ys,fixPara);


            if (fitValues == null)
                    return null;

            xl = fitValues[0]-1.0;
            xr = fitValues[0]+1.0;
            xd = HillConstants.EC50_GRID_SIZE;

            y0l = fitValues[6]-2*y0d;
            y0r = fitValues[6]+2*y0d;

            yl = fitValues[1]-2*yd;
            yr = fitValues[1]+2*yd;

//            delta_y = 0.1*delta_y;
//            if (delta_y < 0.5)
//                delta_y = 0.5;
//
//            y0d = delta_y;
//            yd = delta_y;

            y0d = 0.2*y0d;
            y0d = (y0d < 1.0)?1.0:y0d;

            yd = 0.2*yd;
            yd = (yd < 1.0)?1.0:yd;

            if (fitValues[2]-0.5 < 0.1)
                    sl = 0.1;
            else
                    sl = fitValues[2]-0.4;

            sr = fitValues[2]+0.4;
            sd = 0.2;

            ml = 0.4;
            mr = 2.0;
            md = 0.2;

            fitValues = hillFit5p(y0l,y0r,y0d,yl,yr,yd,xl,xr,xd,sl,sr,sd,ml,mr,md,flags,ws,xs,ys,fixPara);

            return fitValues;
    }

	// 5p fit
	public static double[] hillFit5p(double y0l, double y0r, double y0d, double yl, double yr, double yd, double xl,double xr, double xd, double sl, 
					double sr, double sd, double ml, double mr, double md, boolean[] flags, double[] ws, double[] xs, double[] ys, String fixPara) {
				double y0_min = 0.0;
				double yinf_min = 0.0;
				double slope_min = 0.0;
				double x05_min = 0.0;
				double dev_min = 100000000.0;
				double sym_min = 1.0;

				double delta = sd;

				int no = 0;

				for (double y0 = y0l; y0 <= y0r; y0 += y0d) {
					for (double yinf = yl; yinf <= yr; yinf += yd) {
						for (double x05 = xl; x05 <= xr; x05 += xd) {
							for (double slope = sl; slope <= sr; slope += delta) {
								if (Math.abs(slope) < MIN_SLOPE)
									continue;

								// for curve class > 2.1 or < -2.1
								if (slope > 5.0)
									delta = 5.0;

								for (double sym = ml; sym <= mr; sym += md) {
									double dev = calcHillDeviation(x05, y0, yinf, slope, sym, flags, ws, xs, ys);

									if (dev < dev_min) {
										dev_min = dev;
										slope_min = slope;
										x05_min = x05;
										yinf_min = yinf;
										sym_min = sym;
										y0_min = y0;
									}

									delta = Math.abs(slope) * 0.1;
									if (delta < 0.1)
										delta = 0.1;

									no++;
								}
							}
						}
					}
				}

				double dev_const = calcConstantDeviation(flags, ws, xs, ys);

				// double r2 = 1.0-dev_min/(1.0+dev_const);

				double r2 = 0.0;

				// for curve class 5, compare flat fit with no fit
				if (CURVE_CLASS5) {
					double dev0 = xs.length * CLASSIFICATION_SD * CLASSIFICATION_SD;

					r2 = 1.0 - dev_const / dev0;
				} else {
					r2 = 1.0 - dev_min / (1.0 + dev_const);
				}

				if (r2 < R2)
					return null;

				double pvalue = HillStat.calcPValue(y0_min, yinf_min, x05_min, slope_min, sym_min, xs, ys, flags);
								
				double[] values = { x05_min, yinf_min, slope_min, r2, dev_const, dev_min, y0_min, sym_min, pvalue };

				return values;
			}
	
	
	public static double[] hillFit(double ymin, double ymax, boolean[] flags, double[] ws, double[] xs, double[] ys) {
		return hillFit(ymin, ymax, flags, ws, xs, ys, false, null);
	}

	public static double[] hillFit(double ymin, double ymax, boolean[] flags, double[] ws, double[] xs, double[] ys, boolean flag2, String fixPara) {
		
		double[] values = hillFitFast_1(ymin, ymax, flags, ws, xs, ys, fixPara);

		if (values == null)
			return null;

		double max = -Double.MAX_VALUE;
		double min = Double.MAX_VALUE;

		for (int i = 0; i < xs.length; i++) {
			max = Math.max(max, ys[i]);
			min = Math.min(min, ys[i]);
		}

		double range = max - min;

		// deal with near flat curves
		if (range < 100.0)
			range = 100.0;

		double y0 = values[6];
		double x05 = values[0];
		double yinf = values[1];
		double slope = values[2];

		if (flag2)
			return values;

		double[] yfit = calcHillFitCurve(x05, y0, yinf, slope, xs);

		boolean flag = false;

		// first unmask good points
		int maskNo = 0;

		for (int i = 1; i < ys.length; i++) {
			if (flags[i] == false && Math.abs(ys[i] - yfit[i]) * 100.0 / range < TPV) {
				flag = true;
				flags[i] = true;
			}

			if (flags[i] == false)
				maskNo++;
		}

		ArrayList list = new ArrayList();

		if (maskNo < xs.length / MASK) {
			for (int i = 0; i < ys.length - 1; i++) {
				if (i > 0 && flags[i] == true && Math.abs(ys[i] - yfit[i]) * 100.0 / range > TPV) {
					// flag = true;
					list.add(String.valueOf(i));
					// flags[i]= false;

				} else if (i == 0 && flags[i] == true && Math.abs(ys[i] - yfit[i]) * 100.0 / range > TP0) {
					// flag = true;
					list.add(String.valueOf(i));
					// flags[i]= false;
				}
			}
		}

		if (list.size() > 0 && MASK_FLAG) {
			flag = true;

			for (int i = maskNo; i < xs.length / MASK && (i - maskNo) < list.size(); i++) {
				int k = Integer.parseInt((String) list.get(i - maskNo));

				flags[k] = false;
			}
		}

		if (flag)
			return hillFitFast_1(ymin, ymax, flags, ws, xs, ys, fixPara); 

		return values;
	}

	public static double[] hillFit5P(double ymin, double ymax, boolean[] flags, double[] ws, double[] xs, double[] ys,boolean flag2, String fixPara) {
		double[] values = hillFitFast_1(ymin, ymax, flags, ws, xs, ys, fixPara);

		if (values == null)
			return null;

		double max = -Double.MAX_VALUE;
		double min = Double.MAX_VALUE;

		for (int i = 0; i < xs.length; i++) {
			max = Math.max(max, ys[i]);
			min = Math.min(min, ys[i]);
		}

		double range = max - min;

		// deal with near flat curves
		if (range < 100.0)
			range = 100.0;

		double y0 = values[6];
		double x05 = values[0];
		double yinf = values[1];
		double slope = values[2];
		double sym = values[7];

		if (flag2)
			return values;

		double[] yfit = calcHillFitCurve(x05, y0, yinf, slope, sym, xs);

		boolean flag = false;

		// first unmask good points
		int maskNo = 0;

		for (int i = 1; i < ys.length; i++) {
			if (flags[i] == false && Math.abs(ys[i] - yfit[i]) * 100.0 / range < TPV) {
				flag = true;
				flags[i] = true;
			}

			if (flags[i] == false)
				maskNo++;
		}

		ArrayList list = new ArrayList();

		if (maskNo < xs.length / MASK) {
			for (int i = 0; i < ys.length - 1; i++) {
				if (i > 0 && flags[i] == true && Math.abs(ys[i] - yfit[i]) * 100.0 / range > TPV) {
					// flag = true;
					list.add(String.valueOf(i));
					// flags[i]= false;
				} else if (i == 0 && flags[i] == true && Math.abs(ys[i] - yfit[i]) * 100.0 / range > TP0) {
					// flag = true;
					list.add(String.valueOf(i));
					// flags[i]= false;
				}
			}
		}

		if (list.size() > 0) {
			flag = true;

			for (int i = maskNo; i < xs.length / MASK && (i - maskNo) < list.size(); i++) {
				int k = Integer.parseInt((String) list.get(i - maskNo));

				flags[k] = false;
			}
		}

		return values;
	}
	

	// 4p fit
	public static double[] hillFit(double y0l, double y0r, double y0d, double yl, double yr, double yd, double xl,
			double xr, double xd, double sl, double sr, double sd, boolean[] flags, double[] ws, double[] xs, double[] ys, String fixPara) {
		if (fastFlag)
			return hillFitFast_0(y0l, y0r, y0d, yl, yr, yd, xl, xr, xd, sl, sr, sd, flags, ws, xs, ys, fixPara);

		double y0_min = 0.0;
		double yinf_min = 0.0;
		double slope_min = 0.0;
		double x05_min = 0.0;
		double dev_min = Double.MAX_VALUE;
		double delta = sd;

		int no = 0;

		// from console client
		if (fixPara != null) {
			if (fixPara.indexOf("y0") >= 0) {
				double[] paras = parseFixPara(fixPara, "y0");
				y0l = paras[0];
				y0r = paras[1];
			}

			if (fixPara.indexOf("yinf") >= 0) {
				double[] paras = parseFixPara(fixPara, "yinf");
				yl = paras[0];
				yr = paras[1];
			}

			if (fixPara.indexOf("slope") >= 0) {
				double[] paras = parseFixPara(fixPara, "slope");
				sl = paras[0];
				sr = paras[1];
			}

			if (fixPara.indexOf("log_ac50") >= 0) {
				double[] paras = parseFixPara(fixPara, "log_ac50");
				xl = paras[0];
				xr = paras[1];
			}
		}

		// from qhtsclient
		if (FIXED_AC50 != null && FIXED_AC50.length() > 0) {
			double[] ds = parseConstraintPara(FIXED_AC50);
			xl = ds[0];
			xr = ds[1];
		}

		if (FIXED_Y0 != null && FIXED_Y0.length() > 0) {
			double[] ds = parseConstraintPara(FIXED_Y0);
			y0l = ds[0];
			y0r = ds[1];
		}

		if (FIXED_YINF != null && FIXED_YINF.length() > 0) {
			double[] ds = parseConstraintPara(FIXED_YINF);
			yl = ds[0];
			yr = ds[1];
		}

		if (FIXED_SLOPE != null && FIXED_SLOPE.length() > 0) {
			double[] ds = parseConstraintPara(FIXED_SLOPE);
			sl = ds[0];
			sr = ds[1];
		}

		if (MIN_YINF != null) {
			double td = Double.parseDouble(MIN_YINF);
			if (yl < td)
				yl = td;
		}

		if (MAX_YINF != null) {
			double td = Double.parseDouble(MAX_YINF);
			if (yr > td)
				yr = td;
		}

		if (YINF_MORE_80) {
			if (yr > 0.0) {
				yl = 80.0;

				if (yr < 100.0)
					yr = 100.0;
			} else {
				if (yl > -100.0)
					yl = -100.0;

				yr = -80;
			}
		}

		for (double y0 = y0l; y0 <= y0r; y0 += y0d) {
			for (double yinf = yl; yinf <= yr; yinf += yd) {
				if (Y0_LESS_THAN_YINF && Math.abs(y0) >= Math.abs(yinf))
					continue;

				for (double x05 = xl; x05 <= xr; x05 += xd) {
					for (double slope = sl; slope <= sr; slope += delta) {
						if (Math.abs(slope) < MIN_SLOPE)
							continue;

						// for curve class > 2.1 or < -2.1
						if (slope > 5.0)
							delta = 5.0;

						double dev = calcHillDeviation(x05, y0, yinf, slope, flags, ws, xs, ys);

						if (dev < dev_min) {
							dev_min = dev;
							slope_min = slope;
							x05_min = x05;
							yinf_min = yinf;

							y0_min = y0;
						}

						delta = Math.abs(slope) * 0.1;
						if (delta < 0.1)
							delta = 0.1;

						no++;
					}
				}
			}
		}

		double dev_const = calcConstantDeviation(flags, ws, xs, ys);

		double r2 = 0.0;

		// for curve class 5, compare flat fit with no fit
		if (CURVE_CLASS5) {
			double dev0 = 9 * xs.length * CLASSIFICATION_SD * CLASSIFICATION_SD;

			r2 = 1.0 - dev_min / dev0;
		} else {
			r2 = 1.0 - dev_min / (1.0 + dev_const);
		}

		if (r2 < R2)
			return null;

		double[] values = { x05_min, yinf_min, slope_min, r2, dev_const, dev_min, y0_min };

		return values;
	}

	 
	// 4p fit
	public static double[] hillFitFast_0(double y0l, double y0r, double y0d, double yl, double yr, double yd, double xl,
			double xr, double xd, double sl, double sr, double sd, boolean[] flags, double[] ws, double[] xs,double[] ys, String fixPara) {
		double y0_min = 0.0;
		double yinf_min = 0.0;
		double slope_min = 0.0;
		double x05_min = 0.0;
		double dev_min = 100000000.0;
		double delta = sd;

		int no = 0;

		// from console client
		if (fixPara != null) {
			if (fixPara.indexOf("y0") >= 0) {
				double[] paras = parseFixPara(fixPara, "y0");
				y0l = paras[0];
				y0r = paras[1];
			}

			if (fixPara.indexOf("yinf") >= 0) {
				double[] paras = parseFixPara(fixPara, "yinf");
				yl = paras[0];
				yr = paras[1];
			}

			if (fixPara.indexOf("slope") >= 0) {
				double[] paras = parseFixPara(fixPara, "slope");
				sl = paras[0];
				sr = paras[1];
			}

			if (fixPara.indexOf("log_ac50") >= 0) {
				double[] paras = parseFixPara(fixPara, "log_ac50");
				xl = paras[0];
				xr = paras[1];
			}
		}

		// from qhtsclient
		if (FIXED_AC50 != null && FIXED_AC50.length() > 0) {
			double[] ds = parseConstraintPara(FIXED_AC50);
			xl = ds[0];
			xr = ds[1];
		}

		// any parameter fixed?
		if (FIXED_Y0 != null && FIXED_Y0.length() > 0) {
			double[] ds = parseConstraintPara(FIXED_Y0);
			y0l = ds[0];
			y0r = ds[1];
		}

		if (FIXED_YINF != null && FIXED_YINF.length() > 0) {
			double[] ds = parseConstraintPara(FIXED_YINF);
			yl = ds[0];
			yr = ds[1];
		}

		if (FIXED_SLOPE != null && FIXED_SLOPE.length() > 0) {
			double[] ds = parseConstraintPara(FIXED_SLOPE);
			sl = ds[0];
			sr = ds[1];
		}

		if (YINF_MORE_80) {
			if (yr > 0.0) {
				yl = 80.0;

				if (yr < 100.0)
					yr = 100.0;
			} else {
				if (yl > -100.0)
					yl = -100.0;

				yr = -80;
			}
		}

		for (double x05 = xl; x05 <= xr; x05 += xd) {
			for (double y0 = y0l; y0 <= y0r; y0 += y0d) {
				for (double yinf = yl; yinf <= yr; yinf += yd) {
					if (Y0_LESS_THAN_YINF && Math.abs(y0) >= Math.abs(yinf))
						continue;

					double s0 = calcSlope(x05, y0, yinf, flags, xs, ys);

					if (s0 < MIN_SLOPE)
						s0 = MIN_SLOPE;

					if (s0 > MAX_SLOPE)
						s0 = MAX_SLOPE;

					// at s0
					double dev0 = calcHillDeviation(x05, y0, yinf, s0, flags, ws, xs, ys);
					if (dev0 < dev_min) {
						dev_min = dev0;
						slope_min = s0;
						x05_min = x05;
						yinf_min = yinf;
						y0_min = y0;
					}

					if (slopeFlag == false)
						continue;

					// move right
					sr = 1.5 * s0;
					delta = 0.05 * s0;

					for (double slope = s0 + delta; slope <= sr; slope += delta) {
						double dev = calcHillDeviation(x05, y0, yinf, slope, flags, ws, xs, ys);
						if (dev > dev0)
							break;

						dev0 = dev;

						if (dev < dev_min) {
							dev_min = dev;
							slope_min = slope;
							x05_min = x05;
							yinf_min = yinf;
							y0_min = y0;
						}

						no++;
					}

					// move left
					sl = 0.5 * s0;
					for (double slope = s0 - delta; slope >= sl; slope -= delta) {
						double dev = calcHillDeviation(x05, y0, yinf, slope, flags, ws, xs, ys);
						if (dev > dev0)
							break;

						dev0 = dev;

						if (dev < dev_min) {
							dev_min = dev;
							slope_min = slope;
							x05_min = x05;
							yinf_min = yinf;
							y0_min = y0;
						}

						no++;
					}
				}
			}
		}

		double dev_const = calcConstantDeviation(flags, ws, xs, ys);

		double r2 = 0.0;

		// for curve class 5, compare flat fit with no fit
		if (CURVE_CLASS5) {
			double dev0 = xs.length * CLASSIFICATION_SD * CLASSIFICATION_SD;

			r2 = 1.0 - dev_const / dev0;
		} else {
			r2 = 1.0 - dev_min / (1.0 + dev_const);
		}

		// double r2 = 1.0-dev_min/dev_const;
		// if (r2 < R2)
		// return null;

		double[] values = { x05_min, yinf_min, slope_min, r2, dev_const, dev_min, y0_min };

		return values;
	}  


	// 3p fit
	public static double[] hillFit3p(double y0l, double y0r, double y0d, double yl, double yr, double yd, double xl,
			double xr, double xd, double sl, double sr, double sd, boolean[] flags, double[] ws, double[] xs,double[] ys, String fixPara) {
		double y0_min = 0.0;
		double yinf_min = 0.0;
		double slope_min = 0.0;
		double x05_min = 0.0;
		double dev_min = 100000000.0;
		double delta = sd;

		int no = 0;

		if (fixPara != null) {
			if (fixPara.indexOf("y0") >= 0) {
				double[] paras = parseFixPara(fixPara, "y0");
				y0l = paras[0];
				y0r = paras[1];
			}

			if (fixPara.indexOf("yinf") >= 0) {
				double[] paras = parseFixPara(fixPara, "yinf");
				yl = paras[0];
				yr = paras[1];
			}

			if (fixPara.indexOf("slope") >= 0) {
				double[] paras = parseFixPara(fixPara, "slope");
				sl = paras[0];
				sr = paras[1];
			}

			if (fixPara.indexOf("log_ac50") >= 0) {
				double[] paras = parseFixPara(fixPara, "log_ac50");
				xl = paras[0];
				xr = paras[1];
			}
		}

		for (double y0 = y0l; y0 <= y0r; y0 += y0d) {
			for (double yinf = yl; yinf <= yr; yinf += yd) {
				for (double x05 = xl; x05 <= xr; x05 += xd) {
					double slope = 1.0;

					double dev = calcHillDeviation(x05, y0, yinf, slope, flags, ws, xs, ys);

					if (dev < dev_min) {
						dev_min = dev;
						slope_min = slope;
						x05_min = x05;
						yinf_min = yinf;

						y0_min = y0;
					}

					no++;
				}
			}
		}

		double dev_const = calcConstantDeviation(flags, ws, xs, ys);

		// double r2 = 1.0-dev_min/dev_const;

		double r2 = 0.0;

		// for curve class 5, compare flat fit with no fit
		if (CURVE_CLASS5) {
			double dev0 = xs.length * CLASSIFICATION_SD * CLASSIFICATION_SD;

			r2 = 1.0 - dev_const / dev0;
		} else {
			r2 = 1.0 - dev_min / (1.0 + dev_const);
		}

		if (r2 < R2)
			return null;

		double[] values = { x05_min, yinf_min, slope_min, r2, dev_const, dev_min, y0_min };

		return values;
	}

	
	public static double[] hillFitFast_1(double ymin, double ymax, boolean[] flags, double[] ws, double[] xs, double[] ys, String fixPara) {
		
		if (checkRange(ymin, ymax, flags, ys) == false)
			return null;

		double[] fitValues = lineFit(flags, xs, ys);

		if (fitValues == null) {
			double[] fitValues1 = hillFitFast_2(ymin, ymax, flags, ws, xs, ys, fixPara, true);
			double[] fitValues2 = hillFitFast_2(ymin, ymax, flags, ws, xs, ys, fixPara, false);

			// double[] values = {x05_min, yinf_min, slope_min, r2, dev_const, dev_min,
			// y0_min};

			if (fitValues1 == null && fitValues2 == null)
				return null;
			else if (fitValues1 == null && fitValues2 != null)
				return fitValues2;
			else if (fitValues1 != null && fitValues2 == null)
				return fitValues1;
			else if (fitValues1[3] > fitValues2[3])
				return fitValues1;
			else
				return fitValues2;
		}

		if (fitValues[0] > 0)
			return hillFitFast_2(ymin, ymax, flags, ws, xs, ys, fixPara, true);
		else
			return hillFitFast_2(ymin, ymax, flags, ws, xs, ys, fixPara, false);
	}

	
	public static double[] hillFitFast_2(double ymin, double ymax, boolean[] flags, double[] ws, double[] xs, double[] ys, String fixPara, boolean upFlag) {
		
		if (FIT_TYPE.equals("P5")) {
			double[] fitValues = lineFit(flags, xs, ys);

			if (fitValues[0] > 0)
				return hillFitFast5P(ymin, ymax, flags, ws, xs, ys, fixPara, true);
			else
				return hillFitFast5P(ymin, ymax, flags, ws, xs, ys, fixPara, false);
		}

		if (checkRange(ymin, ymax, flags, ys) == false)
			return null;

		double xl = -10.0;
		if (xs[0] * 1.2 < xl)
			xl = xs[0] * 1.2;

		double xr = -2.0;
		double xd = 0.5;

		double y0l = -150.0;
		double y0r = 150.0;
		double y0d = 5.0;

		double yl = -150.0;
		double yr = 150.0;
		double yd = 5.0;

		double sl = MIN_SLOPE;
		double sr = MAX_SLOPE;
		double sd = 0.2;

		double max_y = -Double.MAX_VALUE;
		double min_y = Double.MAX_VALUE;

		for (int i = 0; i < ys.length; i++)
			if (flags[i]) {
				if (ys[i] > max_y)
					max_y = ys[i];

				if (ys[i] < min_y)
					min_y = ys[i];
			}

		if (upFlag) {
			if (min_y < 0) {
				y0l = Y0_INF_COEF * min_y;
				y0r = 0.0;
			} else {
				y0l = 0.0;
				y0r = Y0_INF_COEF * min_y;
			}

			if (max_y > 0.0) {
				yl = 0.0;
				yr = Y0_INF_COEF * max_y;

				// if (fastFlag && yr < 100.0)
				// yr = 100.0;
			} else {
				yr = 0.0;
				yl = Y0_INF_COEF * max_y;
			}

			sl = MIN_SLOPE;
			sr = MAX_SLOPE;
		} else {
			if (max_y < 0) {
				y0l = Y0_INF_COEF * max_y;
				y0r = 0.0;
			} else {
				y0l = 0.0;
				y0r = Y0_INF_COEF * max_y;
			}

			if (min_y < 0.0) {
				yl = Y0_INF_COEF * min_y;
				yr = 0;

				// if (fastFlag && yl > -100.0)
				// yl = -100.0;
			} else {
				yr = Y0_INF_COEF * min_y;
				yl = 0;
			}

			sl = MIN_SLOPE;
			sr = MAX_SLOPE;
		}

		// y0d = delta_y;
		// yd = delta_y;

		double delta_y = 0.05 * (max_y - min_y);

		if (delta_y < 1e-6)
			delta_y = 0.05 * Math.abs(min_y);

		y0d = 0.05 * (y0r - y0l);
		y0d = (y0d < 2.0) ? 2.0 : y0d;
		// y0d = (y0d < delta_y)?delta_y:y0d;

		yd = 0.05 * (yr - yl);
		yd = (yd < 2.0) ? 2.0 : yd;
		yd = (yd < delta_y) ? delta_y : yd;

		if (delta_y < 1.0 && y0d > 1.0)
			y0d = 1.0;

		if (delta_y < 1.0 && yd > 1.0)
			yd = 1.0;

		//System.out.println(y0l + ":" + y0l + ":" + yl + ":" + yr + ":" + xl + ":" + xr);

		// System.out.println(p4Fit+":"+min_y+":"+max_y+":"+delta_y);
		double[] fitValues = null;

		if (p4Fit && FIT_TYPE.equals("P4"))
			fitValues = hillFit(y0l, y0r, y0d, yl, yr, yd, xl, xr, xd, sl, sr, sd, flags, ws, xs, ys, fixPara);
		else
			fitValues = hillFit3p(y0l, y0r, y0d, yl, yr, yd, xl, xr, xd, sl, sr, sd, flags, ws, xs, ys, fixPara);

		if (fitValues == null)
			return null;

		xl = fitValues[0] - 1.0;
		xr = fitValues[0] + 1.0;
		xd = EC50_GRID_SIZE;

		y0l = fitValues[6] - 2 * y0d;
		y0r = fitValues[6] + 2 * y0d;

		yl = fitValues[1] - 2 * yd;
		yr = fitValues[1] + 2 * yd;

		y0d = 0.1 * y0d;
		y0d = (y0d < 0.5) ? 0.5 : y0d;

		yd *= 0.1;
		yd = (yd < 0.5) ? 0.5 : yd;

		if (fitValues[2] - 0.5 < 0.1)
			sl = 0.1;
		else
			sl = fitValues[2] - 0.5;

		sr = fitValues[2] + 0.5;
		sd = 0.1;

		if (p4Fit && FIT_TYPE.equals("P4"))
			fitValues = hillFit(y0l, y0r, y0d, yl, yr, yd, xl, xr, xd, sl, sr, sd, flags, ws, xs, ys, fixPara);
		else
			fitValues = hillFit3p(y0l, y0r, y0d, yl, yr, yd, xl, xr, xd, sl, sr, sd, flags, ws, xs, ys, fixPara);

		return fitValues;
	}

 
}
