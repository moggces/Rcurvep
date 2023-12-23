package gov.nih.ncats.ifx.qhts.curvefitting.algorithm;
 
import gov.nih.ncats.ifx.qhts.utils.*;
 
public class CurveClass2 {

	private boolean useMask = true;

	private boolean allowBell = true;
	private double sdPct = 10.0;
	private double sdFactor = 3.0;
	
	private double robustPct = 6*sdPct;
	
	private double r2Cutoff = HillConstants.R2_CLASS1_CUTOFF;
	private double asymThresh = 0.75;
	private boolean relativeActivity = true;

	private double pvalueCutoff = 0.05;
   

	public void setSDPct(double sd) {
		sdPct = sd;
		robustPct = 6*sdPct;
	}
      
	
	public double calcCurveClass(double[] fitParas, double[] concs, double[] acts, boolean[] maskFlags, double pvalue) throws Exception { 
		return calcCurveClass(fitParas, 1, concs, acts, maskFlags, pvalue);
	}


	public void setSdFactor(double sdFactor) {
		this.sdFactor = sdFactor;
	}


	public void setRelativeActivity(boolean relativeActivity) {
		this.relativeActivity = relativeActivity;
	}


	public double calcCurveClass(double[] fitParas1, int dn, double[] concs, double[] acts1, boolean[] maskFlags,double pvalue) throws Exception { 

		double[] fitParas = null;
		double[] acts = null;

		// shift or not
		if (NCGCConstants.ACTIVITY_DELTA > 0.0) {
			int no = acts1.length;

			acts = new double[no];
			for (int i = 0; i < no; i++) {
				acts[i] = acts1[i] - NCGCConstants.ACTIVITY_DELTA;
			}

			if (fitParas1 != null) {
				no = fitParas1.length;

				fitParas = new double[no];
				for (int i = 0; i < no; i++) {
					fitParas[i] = fitParas1[i];
				}

				fitParas[1] -= NCGCConstants.ACTIVITY_DELTA;
				fitParas[6] -= NCGCConstants.ACTIVITY_DELTA;
			}
		} else if (HillConstants.SHIFT_ACTIVITY && fitParas1 != null) {
			double zeroAct = fitParas1[6];

			int no = acts1.length;

			acts = new double[no];
			for (int i = 0; i < no; i++) {
				acts[i] = acts1[i] - zeroAct;
			}

			no = fitParas1.length;

			fitParas = new double[no];
			for (int i = 0; i < no; i++) {
				fitParas[i] = fitParas1[i];
			}

			fitParas[1] -= zeroAct;
			fitParas[6] -= zeroAct;
		} else if (HillConstants.SHIFT_ACTIVITY) {
			int no = acts1.length;

			double median = PlateUtil.calcMedian(acts1);

			acts = new double[no];
			for (int i = 0; i < no; i++) {
				acts[i] = acts1[i] - median;
			}

			fitParas = fitParas1;
		} else {
			int no = acts1.length;

			acts = new double[no];
			for (int i = 0; i < no; i++) {
				acts[i] = acts1[i];
			}

			fitParas = fitParas1;
		}

		int no = acts.length;
		// double pctCutoff = sdPct*sdFactor;
		double pctCutoff = HillConstants.CLASSIFICATION_SD * HillConstants.CLASSIFICATION_SD_FACTOR;

		// before relative, check special case
		double maxAct = -Double.MAX_VALUE;
		double minAct = Double.MAX_VALUE;

		for (int i = 0; i < no; i++) {
			double act = acts[i];

			// if (!useMask || maskFlags == null || maskFlags[i])
			if (maskFlags == null || maskFlags[i]) {
				maxAct = Math.max(maxAct, act);
				minAct = Math.min(minAct, act);
			}
		}

		if (HillConstants.AUTO_SUPER_ACTIVE && pvalue > 0.05) {
			// super actives
			int no1 = 0;
			int no2 = 0;
			for (int i = 0; i < acts.length; i++) {
				if (acts[i] < -robustPct)
					no1++;

				if (acts[i] < -6 * sdPct)
					no2++;
			}

			if (1.0 * no1 / acts.length >= HillConstants.SUPER_ACTIVE_RATIO) {
				return -1.01;
			} else if (1.0 * no2 / acts.length >= HillConstants.SUPER_ACTIVE_RATIO && isSuperActive(acts)) {
				return -1.02;
			}
		}

		if (no > 1 && minAct > pctCutoff && maxAct > pctCutoff && maxAct - minAct < pctCutoff)
			return 5.0;

		if (no > 1 && minAct < -pctCutoff && maxAct < -pctCutoff && maxAct - minAct < pctCutoff)
			return 5.0;

		if (pvalue > 0.05 && Math.abs(minAct) < pctCutoff && Math.abs(maxAct) < pctCutoff)
			return 4.0;

		if (fitParas == null)
			return calcCurveClassWithoutFit(fitParas, dn, concs, acts, maskFlags, pvalue);
		else
			return calcCurveClassWithFit(fitParas, dn, concs, acts, maskFlags, pvalue);
	} 

	// with fit
	public double calcCurveClassWithFit(double[] fitParas, int dn, double[] concs, double[] acts, boolean[] maskFlags,double pvalue) throws Exception {
		 
		return calcCurveClassWithFit2(fitParas, dn, concs, acts, maskFlags, pvalue);
	}

	
	// with fit
	public double calcCurveClassWithoutFit(double[] fitParas, double[] concs, double[] acts, boolean[] maskFlags,double pvalue) throws Exception {
			
		return calcCurveClassWithoutFit(fitParas, 1, concs, acts, maskFlags, pvalue);
	}
	
		
	// with fit; 1.* and 2.* based upon max act
	public double calcCurveClassWithFit1(double[] fitParas, int dn, double[] concs, double[] acts, boolean[] maskFlags,double pvalue) throws Exception {
		
		double pctCutoff = sdPct * sdFactor;

		double ac50 = fitParas[0];
		double zeroAct = fitParas[6];
		double infAct = fitParas[1];
		double r2 = fitParas[3];

		// do not work for shifted curves
		// if (Math.abs(zeroAct) > sdPct && Math.abs(infAct) < 2*sdPct)
		// return 4.0;

		if (Math.abs(zeroAct - infAct) < HillConstants.RANGE_CUTOFF * sdPct)
			return 4.0;

		// curve go up or down?
		int sign = (zeroAct < infAct) ? 1 : -1;

		// special case: carryover; manual curation needed
		// if (Math.abs(zeroAct) > Math.abs(infAct) && Math.abs(zeroAct) > pctCutoff &&
		// pvalue < 0.05)
		if (Math.abs(zeroAct) > 6 * sdPct) // && Math.abs(infAct-zeroAct) < pctCutoff)
			return 5.0;

		// inhibitor
		if (zeroAct > infAct && 0.5 * (zeroAct + infAct) > pctCutoff)
			return 5.0;

		// activator
		if (zeroAct < infAct && 0.5 * (zeroAct + infAct) < -pctCutoff)
			return 5.0;

		// if (Math.abs(zeroAct)-Math.abs(infAct) > pctCutoff && pvalue < 0.05)
		// return 5.0;

		int no = acts.length;

		// use relative?
		if (relativeActivity && pvalue < 0.05 && Math.abs(zeroAct) > sdPct) {
			for (int i = 0; i < no; i++) {
				acts[i] = acts[i] - zeroAct;
			}

			infAct -= zeroAct;
			zeroAct = 0.0;
		}

		// some stats
		double maxAct = -Double.MAX_VALUE;
		double minAct = Double.MAX_VALUE;

		int uSDF = 0;
		int uSD = 0;
		int dSD = 0;
		int dSDF = 0;
		int lastPoint = -1;

		double lastAct = 0.0;
		double lastAct2 = 0.0;

		for (int i = 0; i < no; i++) {
			double act = acts[i];

			if (!useMask || maskFlags == null || maskFlags[i]) {
				lastAct2 = lastAct;

				lastAct = act;
				lastPoint = i;

				maxAct = Math.max(maxAct, act);
				minAct = Math.min(minAct, act);

				if (act > sdPct) {
					uSD++;

					if (act > pctCutoff)
						uSDF++;
				} else if (act < -sdPct) {
					dSD++;

					if (act < -pctCutoff)
						dSDF++;
				}

			}
		}

		// inactive; should not be used
		// if (pvalue < 0.05 && Math.max(Math.abs(maxAct),Math.abs(minAct)) < pctCutoff)
		// return 4.0;

		// has asymptote? More than one point within 80% of max efficacy
		int pts = 0;
		for (int i = 0; i < no; i++) {
			if (!useMask || maskFlags == null || maskFlags[i]) {
				if (concs[i] > ac50 && sign * acts[i] > asymThresh * sign * infAct)
					pts++;
			}
		}

		// 1.*
		if (pts >= 2) {
			if (pvalue < pvalueCutoff) {
				if (sign > 0) {
					if (maxAct >= robustPct)
						return 1.1;
					else if (maxAct >= pctCutoff)
						return 1.2;
					// else if (maxAct >= actCutoff24Factor*sdPct)
					else if (maxAct >= HillConstants.C124_FACTOR * sdPct)
						return 1.4;
				} else {
					if (minAct <= -robustPct)
						return -1.1;
					else if (minAct <= -pctCutoff)
						return -1.2;
					// else if (minAct <= -actCutoff24Factor*sdPct)
					else if (minAct <= -HillConstants.C124_FACTOR * sdPct)
						return -1.4;
				}
			} else {
				if (sign > 0) {
					if (maxAct >= robustPct)
						return 1.3;
					else if (maxAct >= pctCutoff)
						return 1.4;
				} else {
					if (minAct <= -robustPct)
						return -1.3;
					else if (minAct <= -pctCutoff)
						return -1.4;
				}
			}
		}

		// 3: second point < half of last data point
		if (sign > 0) {
			if (uSDF == 1 && lastPoint > 0 && maxAct == lastAct && lastAct2 < sdPct) // 0.33*lastAct)
				return 3.0;
		}

		if (sign < 0) {
			if (dSDF == 1 && lastPoint > 0 && minAct == lastAct && lastAct2 > -sdPct) // 0.33*lastAct)
				return -3.0;
		}

		// 2.*
		if (pvalue < pvalueCutoff) {
			if (sign > 0) {
				if (Math.min(maxAct, infAct) >= robustPct)
					return 2.1;
				else if (Math.min(maxAct, infAct) >= pctCutoff)
					return 2.2;
				else if (Math.min(maxAct, infAct) >= HillConstants.C124_FACTOR * sdPct)
					return 2.4;
			} else {
				if (Math.max(minAct, infAct) <= -robustPct)
					return -2.1;
				else if (Math.max(minAct, infAct) <= -pctCutoff)
					return -2.2;
				else if (Math.max(minAct, infAct) <= -HillConstants.C124_FACTOR * sdPct)
					return -2.4;
			}
		} else {
			if (sign > 0) {
				if (Math.min(maxAct, infAct) >= robustPct)
					return 2.3;
				else if (Math.min(maxAct, infAct) >= pctCutoff)
					return 2.4;
			} else {
				if (Math.max(minAct, infAct) <= -robustPct)
					return -2.3;
				else if (Math.max(minAct, infAct) <= -pctCutoff)
					return -2.4;
			}
		}

		return 4.0;
	}

	// with fit 1.* and 2.* based upon yinf
	public double calcCurveClassWithFit2(double[] fitParas, int dn, double[] concs, double[] acts, boolean[] maskFlags,double pvalue) throws Exception {
		
		//System.out.println("\n// with fit 1.* and 2.* based upon yinf");
		//System.out.println("357");
		
		
		double pctCutoff = HillConstants.CLASSIFICATION_SD * HillConstants.CLASSIFICATION_SD_FACTOR;

		double ac50 = fitParas[0];
		double zeroAct = fitParas[6];
		double infAct = fitParas[1];
		double r2 = fitParas[3];

		// do not work for shifted curves
		// if (Math.abs(zeroAct) > sdPct && Math.abs(infAct) < 2*sdPct)
		// return 4.0;

		if (Math.abs(zeroAct - infAct) < HillConstants.RANGE_CUTOFF * sdPct)
			return 4.0;

		// curve go up or down?
		int sign = (zeroAct < infAct) ? 1 : -1;

		// special case: carryover; manual curation needed
		// if (Math.abs(zeroAct) > Math.abs(infAct) && Math.abs(zeroAct) > pctCutoff &&
		// pvalue < 0.05)
		if (Math.abs(zeroAct) > 6 * sdPct) // && Math.abs(infAct-zeroAct) < pctCutoff)
			return 5.0;

		// inhibitor
		if (zeroAct > infAct && 0.5 * (zeroAct + infAct) > pctCutoff)
			return 5.0;

		// activator
		if (zeroAct < infAct && 0.5 * (zeroAct + infAct) < -pctCutoff)
			return 5.0;

		// if (Math.abs(zeroAct)-Math.abs(infAct) > pctCutoff && pvalue < 0.05)
		// return 5.0;

		int no = acts.length;

		// use relative?
		if (relativeActivity && pvalue < 0.05 && Math.abs(zeroAct) > sdPct) {
			for (int i = 0; i < no; i++) {
				acts[i] = acts[i] - zeroAct;
			}

			infAct -= zeroAct;
			zeroAct = 0.0;
		}

		// some stats
		double maxAct = -Double.MAX_VALUE;
		double minAct = Double.MAX_VALUE;

		int uSDF = 0;
		int uSD = 0;
		int dSD = 0;
		int dSDF = 0;
		int lastPoint = -1;

		double lastAct = 0.0;
		double lastAct2 = 0.0;

		for (int i = 0; i < no; i++) {
			double act = acts[i];

			if (!useMask || maskFlags == null || maskFlags[i]) {
				lastAct2 = lastAct;

				lastAct = act;
				lastPoint = i;

				maxAct = Math.max(maxAct, act);
				minAct = Math.min(minAct, act);

				if (act > sdPct) {
					uSD++;

					if (act > pctCutoff)
						uSDF++;
				} else if (act < -sdPct) {
					dSD++;

					if (act < -pctCutoff)
						dSDF++;
				}

			}
		}

		// inactive; should not be used
		// if (pvalue < 0.05 && Math.max(Math.abs(maxAct),Math.abs(minAct)) < pctCutoff)
		// return 4.0;

		// has asymptote? More than one point within 80% of max efficacy
		int pts = 0;
		for (int i = 0; i < no; i++) {
			if (!useMask || maskFlags == null || maskFlags[i]) {
				if (concs[i] > ac50 && sign * acts[i] > asymThresh * sign * infAct)
					pts++;
			}
		}

		// 1.*
		if (pts >= 2) {
			if (pvalue < pvalueCutoff) {
				if (sign > 0) {
					if (Math.min(maxAct, infAct) >= robustPct)
						return 1.1;
					else if (Math.min(maxAct, infAct) >= pctCutoff)
						return 1.2;
					else if (Math.min(maxAct, infAct) >= HillConstants.C124_FACTOR * sdPct)
						return 1.4;
				} else {
					if (Math.max(minAct, infAct) <= -robustPct)
						return -1.1;
					else if (Math.max(minAct, infAct) <= -pctCutoff)
						return -1.2;
					else if (Math.max(minAct, infAct) <= -HillConstants.C124_FACTOR * sdPct)
						return -1.4;
				}
			} else {
				if (sign > 0) {
					if (Math.min(maxAct, infAct) >= robustPct)
						return 1.3;
					else if (Math.min(maxAct, infAct) >= pctCutoff)
						return 1.4;
				} else {
					if (Math.max(minAct, infAct) <= -robustPct)
						return -1.3;
					else if (Math.max(minAct, infAct) <= -pctCutoff)
						return -1.4;
				}
			}
		}

		// 3: second point < half of last data point
		if (sign > 0) {
			if (uSDF == 1 && lastPoint > 0 && maxAct == lastAct && lastAct2 < sdPct) // 0.33*lastAct)
				return 3.0;
		}

		if (sign < 0) {
			if (dSDF == 1 && lastPoint > 0 && minAct == lastAct && lastAct2 > -sdPct) // 0.33*lastAct)
				return -3.0;
		}

		// 2.*
		if (pvalue < pvalueCutoff) {
			if (sign > 0) {
				if (Math.min(maxAct, infAct) >= robustPct)
					return 2.1;
				else if (Math.min(maxAct, infAct) >= pctCutoff)
					return 2.2;
				else if (Math.min(maxAct, infAct) >= HillConstants.C124_FACTOR * sdPct)
					return 2.4;
			} else {
				if (Math.max(minAct, infAct) <= -robustPct)
					return -2.1;
				else if (Math.max(minAct, infAct) <= -pctCutoff)
					return -2.2;
				else if (Math.max(minAct, infAct) <= -HillConstants.C124_FACTOR * sdPct)
					return -2.4;
			}
		} else {
			if (sign > 0) {
				if (Math.min(maxAct, infAct) >= robustPct)
					return 2.3;
				else if (Math.min(maxAct, infAct) >= pctCutoff)
					return 2.4;
			} else {
				if (Math.max(minAct, infAct) <= -robustPct)
					return -2.3;
				else if (Math.max(minAct, infAct) <= -pctCutoff)
					return -2.4;
			}
		}

		return 4.0;
	}

	// no fit
	public double calcCurveClassWithoutFit(double[] fitParas, int dn, double[] concs, double[] acts, boolean[] maskFlags, double pvalue) throws Exception {
		 
		int no = acts.length;

		double median = PlateUtil.calcMedian(acts);
		if (Math.abs(median) > 6 * sdPct)
			return 5.0;

		double pctCutoff = sdPct * sdFactor;

		double maxAct = -Double.MAX_VALUE;
		double minAct = Double.MAX_VALUE;

		int uSDF = 0;
		int uSD = 0;
		int dSD = 0;
		int dSDF = 0;
		int lastPoint = -1;

		double lastAct = 0.0;
		double lastAct2 = 0.0;

		for (int i = 0; i < no; i++) {
			double act = acts[i];

			if (!useMask || maskFlags == null || maskFlags[i]) {
				lastAct2 = lastAct;

				lastAct = act;
				lastPoint = i;

				maxAct = Math.max(maxAct, act);
				minAct = Math.min(minAct, act);

				if (act > sdPct) {
					uSD++;

					if (act > pctCutoff)
						uSDF++;
				} else if (act < -sdPct) {
					dSD++;

					if (act < -pctCutoff)
						dSDF++;
				}

			}
		}

		// up or down
		int sign = 0;

		if (maxAct < pctCutoff && minAct < -pctCutoff)
			sign = -1;
		else if (maxAct > pctCutoff && minAct > -pctCutoff)
			sign = 1;
		else if (uSDF > dSDF)
			sign = 1;
		else if (dSDF > uSDF)
			sign = -1;
		else if (uSD > dSD)
			sign = 1;
		else
			sign = -1; // greater probability of freak up than down

		// replicate
		if (dn > 1) {
			int delta = concs.length / dn;

			double maxc = -10.0;
			double minc = 10.0;

			for (int i = 0; i < dn; i++) {
				double[] xx = new double[delta];
				double[] yy = new double[delta];

				for (int j = 0; j < delta; j++) {
					xx[j] = concs[j * dn + i];
					yy[j] = acts[j * dn + i];
				}

				double cc = calcCurveClassWithoutFit(fitParas, 1, xx, yy, maskFlags, pvalue);

				maxc = Math.max(maxc, cc);
				minc = Math.min(minc, cc);
			}

			if (sign > 0)
				return minc;
			else
				return maxc;
		}

		// If a curve has no fit and only has the last unmasked point >3SD (all other
		// points are <3SD), then we call it class 3, other wise class 4.

		// one point
		if (sign > 0) {
			if (uSDF == 1 && maxAct == lastAct)
				return 3.0;
		}

		if (sign < 0) {
			if (dSDF == 1 && minAct == lastAct)
				return -3.0;
		}

		return 4.0;
	}

	// for inhibition, if first point is the highest and second point second
	// highest, it is not super active
	public boolean isSuperActive(double[] ys) {
		boolean flag1 = true;
		boolean flag2 = true;

		for (int i = 1; i < ys.length; i++) {
			if (ys[i] > ys[0]) {
				flag1 = false;
				break;
			}
		}

		for (int i = 2; i < ys.length; i++) {
			if (ys[i] > ys[1]) {
				flag2 = false;
				break;
			}
		}

		if (flag1 && flag2)
			return false;

		return true;
	}
	
}
