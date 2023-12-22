package gov.nih.ncats.ifx.qhts.curvefitting.algorithm;
 
import java.util.*;
import gov.nih.ncats.ifx.qhts.utils.*; 
import static gov.nih.ncats.ifx.qhts.utils.HillConstants.*; 
import org.apache.commons.math3.stat.regression.SimpleRegression;

public class AlgotithmUtils {

	public static double LN10 = Math.log(10.0);
	
	
	
	public static String composeMaskFlag(boolean[] flags){
		
           if (flags == null)
               return null;

           StringBuffer buffer = new StringBuffer();
           for (int i=0; i<flags.length; i++){
               if (flags[i])
                   buffer.append("0 ");
               else
                   buffer.append("1 ");
           }

           return buffer.toString().trim();
    }


	public static double[] parseFixPara(String para, String key) {
		if (para == null)
			return null;

		int k = para.indexOf(key);
		int l = para.indexOf(";", k + 1);
		if (l < 0)
			l = para.length();

		String tok = para.substring(k, l);
		String[] cols = tok.split(",");

		double[] values = new double[2];
		values[0] = Double.parseDouble(cols[1]);

		if (cols.length == 2)
			values[1] = Double.parseDouble(cols[1]);
		else
			values[1] = Double.parseDouble(cols[2]);

		return values;
	}
	

	public static double[] parseConstraintPara(String para) {

		double[] ds = new double[2];

		StringTokenizer toks = new StringTokenizer(para, ", ;");
		if (toks.countTokens() == 1) {
			ds[0] = ds[1] = Double.parseDouble(toks.nextToken());
		} else {
			ds[0] = Double.parseDouble(toks.nextToken());
			ds[1] = Double.parseDouble(toks.nextToken());
		}

		return ds;
	}
	

	public static double calcSlope(double x05, double y0, double yinf, boolean[] flags, double[] xs, double[] ys) {
		
		double[] x = new double[xs.length];
		double[] y = new double[xs.length];

		int no = 0;
		ArrayList<String> list = new ArrayList<>();

		for (int i = 0; i < xs.length; i++) {
			if (flags != null && flags[i] == false)
				continue;

			if (Math.abs(ys[i] - y0) < 0.01)
				continue;

			double td1 = (yinf - ys[i]) / (ys[i] - y0);
			if (td1 < 0.0)
				continue;

			td1 = Math.log10(td1);

			double td2 = x05 - xs[i];

			if (no > 0 && Math.abs(td1 - y[no - 1]) < 0.01) {
				if (no > 1) {
					double slope = calcSlope(no, x, y);
					list.add(String.valueOf(slope));
				}
				no = 0;
			}

			// System.out.printf("%5d %8.3f %8.3f\n",i,td1,td2);

			x[no] = td2;
			y[no] = td1;

			no++;
		}

		if (no > 1) {
			double slope = calcSlope(no, x, y);
			list.add(String.valueOf(slope));
		}

		if (list.size() == 0)
			return -1;

		return calcMedian(list);
		//return PlateUtil.calcMedian(list);
	}

	
	public static double calcMedian(ArrayList<String> list){
		
        double[] values = new double[list.size()];

        for (int k=0; k<list.size(); k++){
            values[k] = Double.parseDouble((String)list.get(k));
        }

        Arrays.sort(values);

        return values[values.length/2];
    }

	
	public static double calcSlope(int no, double[] x, double[] y) {
		
		double xm = 0;
		double ym = 0;

		int k = 0;
		for (int i = 0; i < no; i++) {
			xm += x[i];
			ym += y[i];

			k++;
		}

		xm /= k;
		ym /= k;

		double lxx = 0;
		////double lyy = 0;
		double lxy = 0;
		for (int i = 0; i < no; i++) {
			lxx += (x[i] - xm) * (x[i] - xm);
			// lyy += (y[i]-ym)*(y[i]-ym);
			lxy += (x[i] - xm) * (y[i] - ym);
		}

		double b = lxy / lxx;
		// double a = ym-b*xm;
		// double r = lxy/Math.sqrt(lxx*lyy);

		return b;
	}

	
	public static double calcHillDeviation(double x05, double y0, double yinf, double slope, boolean[] flags,
			double[] ws, double[] xs, double[] ys) {
		double dev = 0.0;

		for (int i = 0; i < xs.length; i++)
			if (flags == null || flags[i]) {

				double y = y0 + (yinf - y0) / (1.0 + Math.exp(LN10 * slope * (x05 - xs[i])));
				double delta = (y - ys[i]);

				if (ws != null)
					dev += ws[i] * delta * delta;

				else
					dev += delta * delta;
			}

		return dev;
	}

	
	public static double calcHillDeviation(double x05, double y0, double yinf, double slope, double sym,
			boolean[] flags, double[] ws, double[] xs, double[] ys) {
		double dev = 0.0;

		for (int i = 0; i < xs.length; i++)
			if (flags == null || flags[i]) {

				double y = y0 + (yinf - y0) / Math.pow(1.0 + Math.exp(LN10 * slope * (x05 - xs[i])), sym);
				double delta = (y - ys[i]);

				if (ws != null)
					dev += ws[i] * delta * delta;

				else
					dev += delta * delta;
			}

		return dev;
	}

	
	public static ArrayList<String> calcHillDeviationAsList(double x05, double y0, double yinf, double slope,
			boolean[] flags, double[] ws, double[] xs, double[] ys) {
		ArrayList<String> list = new ArrayList<>();

		for (int i = 0; i < xs.length; i++) {

			double y = y0 + (yinf - y0) / (1.0 + Math.exp(LN10 * slope * (x05 - xs[i])));

			double delta = (y - ys[i]);

			list.add(String.valueOf(delta));
		}

		return list;
	}

	
	public static double calcConstantDeviation(boolean[] flags, double[] ws, double[] xs, double[] ys) {

		double mean = 0.0;
		int no = 0;

		for (int i = 0; i < xs.length; i++)

			if (flags == null || flags[i]) {
				mean += ys[i];
				no++;
			}

		mean /= no;
		double dev = 0.0;

		for (int i = 0; i < xs.length; i++)

			if (flags == null || flags[i]) {

				double delta = (ys[i] - mean);

				if (ws != null)
					dev += ws[i] * delta * delta;
				else
					dev += delta * delta;
			}

		return dev;
	}

	
	public static double[] calcHillFitConstant(double[] xs, double[] ys) {

		double[] yy = new double[xs.length];

		double mean = 0.0;
		for (int i = 0; i < xs.length; i++)
			mean += ys[i];
		mean /= xs.length;

		for (int i = 0; i < xs.length; i++) {
			yy[i] = mean;
		}

		return yy;
	}
	

	public static boolean[] parseMaskFlag(String maskFlag) {

		if (maskFlag == null)
			return null;

		String[] cols = maskFlag.split(" ");
		boolean[] flags = new boolean[cols.length];

		for (int i = 0; i < cols.length; i++) {
			if (cols[i].equals("1"))
				flags[i] = false;
			else
				flags[i] = true; // true means not masked
		}

		return flags;
	}
	

	public static boolean checkRange(boolean[] flags, double[] ys) {
		double y_min = 1000000.0;

		double y_max = -1000000.0;

		for (int i = 0; i < ys.length; i++)

			if (flags == null || flags[i]) {

				if (ys[i] < y_min)
					y_min = ys[i];

				if (ys[i] > y_max)
					y_max = ys[i];
			}

		if (y_max - y_min < MIN_Y_RANGE)
			return false;

		return true;
	}

	
	public static boolean checkRange(double ymin, double ymax, boolean[] flags, double[] ys) {

		double y_min = Double.MAX_VALUE;

		double y_max = -Double.MAX_VALUE;

		for (int i = 0; i < ys.length; i++)

			if (flags != null && flags[i]) {
				if (ys[i] < y_min)
					y_min = ys[i];

				if (ys[i] > y_max)
					y_max = ys[i];
			}

		double range = y_max - y_min;

		if (y_min > SUPER_Y || y_max < -SUPER_Y) {
			return true;
		}

		// MIN_Y_RANGE = sd/2
		// if ( 100.0*range/Math.max(Math.abs(y_min), Math.abs(y_max)) < MIN_Y_RANGE)
		if (range < MIN_Y_RANGE) {
			// System.out.println("range check: "+(100.0*range/Math.max(Math.abs(y_min),
			// Math.abs(y_max))));
			return false;
		}

		return true;
	}

	
	public static boolean isOnePoint(double[] ys) {
		int no = 0;

		for (int i = 0; i < ys.length; i++) {
			if (Math.abs(ys[i]) < Y0)
				no++;
		}

		if (no < ys.length - 1)
			return false;

		if (ys.length < 3)
			return false;

		if (Math.abs(ys[ys.length - 1]) < Y0 && Math.abs(ys[ys.length - 2]) < Y0)
			return true;

		return false;
	}

	
	public static ArrayList<String> getMaskList(boolean[] flags, double[] xs, double[] ys) {

		ArrayList<String> list = new ArrayList<>();

		for (int i = 0; i < flags.length; i++) {
			if (flags[i] == false)
				list.add("0-" + String.valueOf(i));
		}

		return list;
	}

	
	public static double[] lineFit(boolean flags[], double[] xs, double ys[]) {
		int no = xs.length;

		double xm = 0;
		double ym = 0;
		int nn = 0;

		double minY = 0.0;
		double maxY = 0.0;

		for (int i = 0; i < no; i++)
			if (flags[i]) {
				xm += xs[i];
				ym += ys[i];

				minY = Math.min(minY, ys[i]);
				maxY = Math.max(maxY, ys[i]);

				nn++;
			}

		xm /= nn;
		ym /= nn;

		double lxx = 0;
		double lyy = 0;
		double lxy = 0;

		for (int i = 0; i < no; i++)
			if (flags[i]) {
				lxx += (xs[i] - xm) * (xs[i] - xm);
				lyy += (ys[i] - ym) * (ys[i] - ym);
				lxy += (xs[i] - xm) * (ys[i] - ym);
			}

		double b = lxy / lxx;
		double a = ym - b * xm;
		double r = lxy / Math.sqrt(lxx * lyy);

		double[] values = { b, a, r };

		if (b > 0.0 && minY < -CLASSIFICATION_SD * CLASSIFICATION_SD_FACTOR)
			return null;

		if (b < 0.0 && maxY > CLASSIFICATION_SD * CLASSIFICATION_SD_FACTOR)
			return null;

		return values;
	}

	
	public static double[] linearFit2(boolean flags[], double[] xs, double ys[]) {
		return lineFit(flags, xs, ys);
	}

	
	public static double dist(double xk, double yk, double xl, double yl) {

		double dy = 0.1 * (yk - yl);
		double dx = xk - xl;
		return Math.sqrt(dx * dx + dy * dy);
	}

	
	public static double dist(double[] xs, double[] ys, int k, int l) {

		double dy = 0.1 * (ys[k] - ys[l]);
		double dx = xs[k] - xs[l];
		return Math.sqrt(dx * dx + dy * dy);
	}
	

	public static double calcTheta(double[] xs, double[] ys, int k, int l, int m) {

		double dkl = dist(xs, ys, k, l);
		double dlm = dist(xs, ys, l, m);
		double dkm = dist(xs, ys, k, m);
		double cos = (dkl * dkl + dlm * dlm - dkm * dkm) / (2 * dkl * dlm);

		// System.out.println(cos);
		return 57.0 * Math.acos(cos);
	}

	
	public static double calcTheta(double xk, double yk, double xl, double yl, double xm, double ym) {

		double dkl = dist(xk, yk, xl, yl);
		double dlm = dist(xl, yl, xm, ym);
		double dkm = dist(xk, yk, xm, ym);
		double cos = (dkl * dkl + dlm * dlm - dkm * dkm) / (2 * dkl * dlm);
		// System.out.println(cos);

		return 57.0 * Math.acos(cos);
	}

	
	public static HashMap<String, String> maskInfoMap(double xs, double ys, double t, String reason) {

		HashMap<String, String> map = new HashMap<>();

		map.put("sample_conc", String.valueOf(xs));
		map.put("sample_resp", String.valueOf(ys));
		map.put("mask_value", String.valueOf(t));
		map.put("mask_reason", reason);

		return map;
	}

	
	// ec conversion; log unit
	public static double eccalc(double ec50, double slope, double pct) {
		return ec50 - (1 / slope) * Math.log10((100 - pct) / pct);
	}

	
	// ic or absolute
	public static double iccalc(double y0, double yinf, double ec50, double slope, double pct) {
		return ec50 - (1 / slope) * Math.log10((yinf - y0) / (pct - y0) - 1);
	}
	

	public static double[] calcHillFitCurve(double x05, double y0, double yinf, double slope, double[] xs) {

		double[] ys = new double[xs.length];

		for (int i = 0; i < xs.length; i++)
			ys[i] = y0 + (yinf - y0) / (1.0 + Math.exp(LN10 * slope * (x05 - xs[i])));

		return ys;
	}
	

	public static double calcHillResponse(double x05, double y0, double yinf, double slope, double xs) {

		return y0 + (yinf - y0) / (1.0 + Math.exp(LN10 * slope * (x05 - xs)));
	}
	

	public static double[] calcHillFitCurve(double x05, double y0, double yinf, double slope, double sym, double[] xs) {

		double[] ys = new double[xs.length];

		for (int i = 0; i < xs.length; i++)
			ys[i] = y0 + (yinf - y0) / Math.pow(1.0 + Math.exp(LN10 * slope * (x05 - xs[i])), sym);

		return ys;
	}
	

	public static boolean[] maskAngle(double[] xs, double[] ys, int[][] ds, double[] ms) {

		int n = xs.length - 1;

		boolean[] flags = new boolean[xs.length];

		for (int i = 0; i < xs.length; i++)
			flags[i] = true;

		if (P0 > 0.0 && Math.abs(ys[0]) > P0)
			flags[0] = false;

		for (int i = 1; i < n; i++) {
			double theta = 0.0;

			if (i == 1)
				theta = calcTheta(xs[i - 1], 0.0, xs[i], ys[i], xs[i + 1], ys[i + 1]);

			else
				theta = calcTheta(xs, ys, i - 1, i, i + 1);
			// System.out.println(i+":"+theta);

			ds[i][(int) theta]++;
		}

		return flags;
	}

	
	public static boolean[] maskLine(double[] xs, double[] ys) {

		boolean[] flag1 = new boolean[xs.length];

		for (int i = 0; i < xs.length; i++)
			flag1[i] = true;

		double[] fitValues = HtsUtil.lineFit(xs, ys);

		// System.out.println(fitValues[0] + ":" + fitValues[1] + ":" + fitValues[2]);

		double dev = 0.0;

		for (int i = 0; i < xs.length; i++) {
			double delta = ys[i] - fitValues[1] - fitValues[0] * xs[i];
			dev += delta * delta;
		}

		dev = Math.sqrt(dev / xs.length);

		boolean[] flag2 = new boolean[xs.length];

		for (int i = 0; i < xs.length; i++) {
			double y0 = fitValues[1] - fitValues[0] * xs[i];

			if (Math.abs(ys[i] - y0) > 3 * dev)
				flag2[i] = false;
			else
				flag2[i] = true;

			////double ratio = Math.abs(ys[i] - y0) / dev;
			//System.out.println(i + ":" + y0 + ":" + dev + ":" + (ys[i] - y0) + ":" + ratio + ":" + flag2[i]);
		}

		return flag2;
	}
	
	
	public static boolean checkBellShape(double[] xs, double[] ys1, double sdf, double sd) {
		double[] ys = new double[ys1.length];
		for (int i = 0; i < ys1.length; i++)
			ys[i] = ys1[i]; // Math.abs(ys1[i]);

		// check bell curve peak > sdf*sd
		double maxY = 0.0;

		// check position
		int pos = 0;

		for (int i = ys.length / BELL_CHECK_START; i < ys.length; i++) {
			double yy = ys[i];

			if (BELL_ABS)
				yy = Math.abs(ys[i]);

			if (yy > maxY) {
				maxY = yy;
				pos = i;
			}
		}

		if (pos <= 1 || maxY < sdf * sd)
			return false;

		// no of significant points
		int no1 = 0;

		for (int i = 0; i < ys.length; i++) {
			double yy = ys[i];

			if (BELL_ABS)
				yy = Math.abs(ys[i]);

			if (yy >= sd)
				no1++;
		}

		if (no1 < 2)
			return false;

		if (BELL_CHECK_REGRESSION) {
			// second check the data points after
			SimpleRegression sr = new SimpleRegression();

			int no = 0;
			for (int j = pos + 1; j < ys.length; j++) {
				sr.addData(xs[j], ys[j]);
				no++;
			}

			////double r = sr.getR();
			double s = sr.getSlope();
			////double c = sr.getIntercept();
			double p = 1.0;

			try {
				if (no > 2)
					p = sr.getSignificance();
			} catch (Exception ex) {
				ex.printStackTrace();
			}

			if (p < 0.05 && s > 0.0)
				return false;

			if (p < 0.05 && s < 0.0)
				return true;

			// first check a supporting point; plus the peak > 3*sd
			// left side
			// if (pos <= 1 || ys[pos] < 3*sd || ys[pos-1] < ys[pos]/3)
			if (ys[pos - 1] < ys[pos] / 3)
				return false;
		}

		return true;
	}

}
