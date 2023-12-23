package gov.nih.ncats.ifx.qhts.curvefitting.algorithm;

import java.util.*;

import static gov.nih.ncats.ifx.qhts.curvefitting.algorithm.AlgotithmUtils.*;
import static gov.nih.ncats.ifx.qhts.utils.HillConstants.*;


public class Masking {
	
	
	public static boolean[] mask(double[] xs, double[] ys) {
		
		double sd = CLASSIFICATION_SD;
		double sdf = CLASSIFICATION_SD_FACTOR;

		return maskDif(xs, ys, null, null, null, 1, sd, sdf); 
	}


	public static boolean[] maskDif(double[] xs, double[] ys, int[][] ds, double[] ms, HashMap[] maps, int dn, double sd, double sdf) {
		// replicate
		if (dn > 1) {
			boolean[] flags = new boolean[xs.length];

			int delta = xs.length / dn;

			for (int i = 0; i < dn; i++) {
				double[] xx = new double[delta];
				double[] yy = new double[delta];

				for (int j = 0; j < delta; j++) {
					xx[j] = xs[j * dn + i];
					yy[j] = ys[j * dn + i];
				}

				boolean[] fs = maskDif(xx, yy, null, null, null, 1, sd, sdf);

				for (int j = 0; j < delta; j++)
					flags[j * dn + i] = fs[j];
			}

			return flags;
		}
		

		boolean[] flags = new boolean[xs.length];

		for (int i = 0; i < xs.length; i++)
			flags[i] = true;

		if (xs.length < 3 || MASK_FLAG == false)
			return flags;

		// one point masking

		int pos = -1;
		int no = 0;

		for (int i = 0; i < ys.length; i++) {
			if (Math.abs(ys[i]) >= sd * sdf) {
				no++;
				pos = i;
			}
		}

		if (BELL_MASK == false && no == 1 && pos < ys.length / 2) {
			flags[pos] = false;
			return flags;
		}

		// calc mean and std
		no = 0;

		double mean = 0.0;
		double minY = Double.MAX_VALUE;
		double maxY = -Double.MAX_VALUE;

		for (int i = 0; i < xs.length; i++) {
			maxY = Math.max(maxY, ys[i]);
			minY = Math.min(minY, ys[i]);

			mean += ys[i];

			no++;
		}
		mean /= no;

		double range = maxY - minY;
		if (range < 100.0)
			range = 100.0;

		double dev = 0.0;

		for (int i = 0; i < xs.length; i++) {
			double delta = ys[i] - mean;
			dev += delta * delta;
		}

		dev = Math.sqrt(dev / no);

		no = xs.length - 1;

		// if (Math.abs(ys[0]) > P0 && Math.abs(ys[0])-Math.abs(ys[1]) > MIN_Y_RANGE)
		// make it more stringent
		if (Math.abs(ys[0]) > P0 && Math.abs(ys[0]) - Math.abs(ys[1]) > 2 * sd) {
			//System.out.println("abs(y0) > abs(y1)");

			flags[0] = false;

			if (maps != null) 
				maps[0] = maskInfoMap(xs[0], ys[0], ys[0], "abs(y0) > abs(y1):" + String.valueOf(ys[0]) + " " + String.valueOf(ys[1]));
		}

		ArrayList<String> list = new ArrayList<>();

		boolean prevVFlag = false;

		for (int i = 1; i < xs.length; i++) {
			
			// check against deviation from mean; make it stringent
			double ratio = Math.abs(ys[i] - mean) / dev;

			if (TR > 0.0 && ratio > TR) {
				//System.out.println("(resp-mean)/std > Pref.TR");

				list.add(String.valueOf(i) + ":" + String.valueOf(ys[i]));

				if (maps != null)
					maps[i] = maskInfoMap(xs[i], ys[i], ratio, "(resp-mean)/std > Pref.TR:" + String.valueOf(ratio) + " " + String.valueOf(TR));

				continue;
			}

			// shape like v
			boolean vFlag = false;
			if (i < no && (ys[i] - ys[i - 1]) * (ys[i] - ys[i + 1]) > 0.0 && Math.abs(ys[i] - ys[i - 1]) > V_DEPTH
					&& Math.abs(ys[i] - ys[i + 1]) > V_DEPTH)
				vFlag = true;

			// check against angle and delta-the difference from the extrapolated line
			double t1 = 0.0;
			if (i > 1)
				t1 = calcTheta(xs, ys, i - 2, i - 1, i);

			double t2 = 0.0;
			if (i < no)
				t2 = calcTheta(xs, ys, i - 1, i, i + 1);

			double t = 0.0;
			if (vFlag == false)
				t = Math.max(t1, t2);
			else
				t = Math.min(t1, t2);

			// this and previous
			double delta0 = Math.abs(ys[i] - ys[i - 1]);
			double delta1 = Double.MAX_VALUE;
			double delta2 = Double.MAX_VALUE;

			if (i == no) {
				delta1 = Math.abs(ys[no] - (xs[no] - xs[no - 2]) * (ys[no - 1] - ys[no - 2]) / (xs[no - 1] - xs[no - 2])
						- ys[no - 2]);
			} else {
				delta1 = Math.abs(
						ys[i] - (xs[i] - xs[i - 1]) * (ys[i + 1] - ys[i - 1]) / (xs[i + 1] - xs[i - 1]) - ys[i - 1]);
				delta2 = Math.abs(ys[i] - ys[i + 1]);
			}

			double delta = Double.MAX_VALUE;
			delta = Math.min(delta0, delta1);
			delta = Math.min(delta2, delta);

			//System.out.println(i + ":" + ":" + delta0 + ":" + delta1 + ":" + delta2 + ":" + t1 + ":" + t2);

			if (prevVFlag && i == no) {
				if (delta * 100.0 / range > TPN && TPN > 0.0) {
					//System.out.println(delta + ":" + TPN + ":" + "Last point's delta > Pref.TPN");

					list.add(String.valueOf(i) + ":" + String.valueOf(ys[i]));

					if (maps != null)
						maps[i] = maskInfoMap(xs[i], ys[i], delta, "Last point's delta distance > Pref.TPN:" + String.valueOf(delta) + " "+ String.valueOf(TPN));
				}

				continue;
			}

			if (i != no && delta * 100.0 / range > TP && TP > 0.0) {
				//System.out.println(delta + ":" + TP + ":" + "Current point's delta > Pref.TP");

				list.add(String.valueOf(i) + ":" + String.valueOf(ys[i]));

				if (maps != null)
					maps[i] = maskInfoMap(xs[i], ys[i], delta, "Current point's delta distance > Pref.TP:"+ String.valueOf(delta * 100.0 / range) + " " + String.valueOf(TP));

				continue;
			} else if (vFlag && t < THETA && delta * 100.0 / range > TPV && THETA > 0.0 && TPV > 0.0) {
				//System.out.println(i + ":" + t + ":" + THETA + ":" + delta + ":" + TPN + ":"+ "Current point's delta > Pref.TPV");

				list.add(String.valueOf(i) + ":" + String.valueOf(ys[i]));

				if (maps != null)
					maps[i] = maskInfoMap(xs[i], ys[i], t, "Current point's delta angle > Pref.TPV:" + String.valueOf(t) + " "+ String.valueOf(THETA));
			}

			if (ds != null) {
				if (delta > 499)
					delta = 499;

				ds[i][(int) delta]++;
			}

			prevVFlag = vFlag;
		}

		// bell masking
		if (BELL_MASK) {
			boolean bellFlag = checkBellShape(xs, ys, sdf, sd);

			if (bellFlag) {
				// find maxY
				maxY = 0.0;
				int maxI = -1;

				for (int i = ys.length / BELL_CHECK_START; i < ys.length; i++) {
					double yy = ys[i];

					if (BELL_ABS)
						yy = Math.abs(ys[i]);

					if (yy > maxY) {
						maxY = yy;
						maxI = i;
					}
				}

				// assign masked point
				for (int j = 0; j < list.size(); j++) {
					String line = (String) list.get(j);
					int k = line.indexOf(":");

					int l = Integer.parseInt(line.substring(0, k));
					flags[l] = false;
				}

				// mask bell shape
				for (int i = ys.length - 1; i >= ys.length / BELL_CHECK_START; i--) {
					double yy = ys[i];

					if (BELL_ABS)
						yy = Math.abs(ys[i]);

					// unmask if masked
					if (Math.abs(yy - maxY) < 1.0e-6) {
						int k = i / dn;

						for (int j = k * dn; j < k * dn + dn; j++) {
							if (j >= flags.length)
								continue;

							flags[j] = true;

							if (maps != null)
								maps[j] = null;
						}

						break;
					}

					if (flags[i] == false)
						continue;

					double td = ys[maxI] * ys[i];

					if ((td > 0.0 && yy < maxY * 0.80) || td < 0.0) {
						flags[i] = false;

						if (maps != null)
							maps[i] = maskInfoMap(xs[i], ys[i], ys[i], "bell y<maxY:" + String.valueOf(ys[i]) + " " + String.valueOf(maxY));
					}
				}

				return flags;
			}
		}

		no = xs.length / MASK;

		for (int i = 0; i < no; i++) {
			if (list.size() == 0)
				break;

			int maxj = -1;
			int maxl = -1;

			double maxd = -Double.MAX_VALUE;

			for (int j = 0; j < list.size(); j++) {
				String line = (String) list.get(j);
				int k = line.indexOf(":");

				int l = Integer.parseInt(line.substring(0, k));
				double d = Double.parseDouble(line.substring(k + 1));

				if (d > maxd) {
					maxd = d;
					maxl = l;
					maxj = j;
				}
			}

			flags[maxl] = false;
			list.remove(maxj);

		}

		for (int i = 0; i < xs.length; i++) {
			if (flags[i] && maps != null)
				maps[i] = null;
		}

		if (USER_MASK != null && USER_MASK.length() > 0) {
			int l = flags.length;

			StringTokenizer toks = new StringTokenizer(USER_MASK, " ,");

			while (toks.hasMoreTokens()) {
				String tok = toks.nextToken();

				if (tok.startsWith("-")) {
					int k = Integer.parseInt(tok.substring(1));
					flags[l - 1 - k] = false;
				} else {
					int k = Integer.parseInt(tok);
					flags[k] = false;
				}

			}
		}
		return flags;
	}

	
	public static boolean[] maskDif1(double[] xs, double[] ys, int[][] ds, double[] ms){

		boolean[] flags = new boolean[xs.length];

		for (int i = 0; i < xs.length; i++)
			flags[i] = true;

		if (xs.length < 3)
			return flags;

		int no = 0;

		double mean = 0.0;

		for (int i = 0; i < xs.length; i++){

			mean += ys[i];
			no++;
		}

		mean /= no;

		double dev = 0.0;

		for (int i = 0; i < xs.length; i++){

			double delta = ys[i] - mean;
			dev += delta * delta;
		}

		dev = Math.sqrt(dev / no);

		if (P0 > 0.0 && Math.abs(ys[0]) > P0){
			//System.out.println("ABS of y0 >Pref.P0");

			flags[0] = false;
		}

		no = xs.length - 1;

		ArrayList list = new ArrayList();

		for (int i = 1; i < xs.length; i++){

			if (Math.abs(ys[i]) > PMAX){
				list.add(String.valueOf(i) + ":" + String.valueOf(ys[i]));

				continue;
			}

			// based upon deviation from mean

			double ratio = Math.abs(ys[i] - mean) / dev;

			//System.out.println(i + ":" + ratio);

			if (ratio > TR){
				//System.out.println("Deviation from mean > Pref.TR");

				list.add(String.valueOf(i) + ":" + String.valueOf(ys[i]));

				continue;
			}

			// based upon delta

			double delta0 = Math.abs(ys[i] - ys[i - 1]);

			double delta2 = 10000000.0;

			double delta1 = 10000000.0;

			// shape like v

			boolean vFlag = false;

			double t1 = 0.0;

			if (i > 1)
				t1 = calcTheta(xs, ys, i - 2, i - 1, i);

			double t2 = 0.0;

			if (i < no){

				if ((ys[i] - ys[i - 1]) * (ys[i] - ys[i + 1]) > 0.0)
					vFlag = true;

				t2 = calcTheta(xs, ys, i - 1, i, i + 1);
			}

			double t = 0.0;

			if (vFlag == false)
				t = Math.max(t1, t2);
			else
				t = Math.min(t1, t2);

			if (i == no){
				delta1 = Math.abs(ys[no] - (xs[no] - xs[no - 2]) * (ys[no - 1] - ys[no - 2]) / (xs[no - 1] - xs[no - 2])- ys[no - 2]);

			}else{
				delta1 = Math.abs(ys[i] - (xs[i] - xs[i - 1]) * (ys[i + 1] - ys[i - 1]) / (xs[i + 1] - xs[i - 1]) - ys[i - 1]);

				delta2 = Math.abs(ys[i] - ys[i + 1]);
			}

			double delta = 1000000.0;

			delta = Math.min(delta0, delta1);

			delta = Math.min(delta2, delta);

			//System.out.println(i + ":" + ":" + delta0 + ":" + delta1 + ":" + delta2 + ":" + t1 + ":" + t2);

			if (i == no || i == no - 1){

				if (ys[i] - ys[i - 1] > 0.0 && delta > TPV){
					//System.out.println(delta + " First and last point's delta > Pref.TPV");

					list.add(String.valueOf(i) + ":" + String.valueOf(ys[i]));

				}else if (delta > TPN){
					//System.out.println("Current point's delta > Pref.TPN");

					list.add(String.valueOf(i) + ":" + String.valueOf(ys[i]));
				}

			}else if (delta > TP){
				//System.out.println("Current point's delta > Pref.TP");

				list.add(String.valueOf(i) + ":" + String.valueOf(ys[i]));
				
			}else if (t < THETA && delta > TPV){
				//System.out.println("Current point's delta > Pref.TPV2");

				list.add(String.valueOf(i) + ":" + String.valueOf(ys[i]));
			}

			if (ds != null){

				if (delta > 499)
					delta = 499;

				ds[i][(int) delta]++;
			}
		}

		no = xs.length / MASK;

		for (int i = 0; i < no; i++){

			if (list.size() == 0)
				break;

			int maxj = -1;

			int maxl = -1;

			double maxd = -100000.0;

			for (int j = 0; j < list.size(); j++){

				String line = (String) list.get(j);

				int k = line.indexOf(":");

				int l = Integer.parseInt(line.substring(0, k));

				double d = Double.parseDouble(line.substring(k + 1));

				if (d > maxd){
					maxd = d;
					maxl = l;
					maxj = j;
				}
			}

			flags[maxl] = false;
			list.remove(maxj);
		}
		// System.out.println(i+":"+mean+":"+dev+":"+ys[i]+":"+(ys[i]-mean)+":"+ratio+":"+ratio2+":"+flag2[i]);

		return flags;
	}


	public static boolean[] maskDif2(double[] xs, double[] ys, int[][] ds, double[] ms) {

		boolean[] flags = new boolean[xs.length];

		for (int i = 0; i < xs.length; i++)
			flags[i] = true;

		if (xs.length < 3)
			return flags;

		// calc mean and std

		int no = 0;
		double mean = 0.0;
		for (int i = 0; i < xs.length; i++){

			mean += ys[i];
			no++;
		}

		mean /= no;
		double dev = 0.0;

		for (int i = 0; i < xs.length; i++){
			double delta = ys[i] - mean;
			dev += delta * delta;
		}

		dev = Math.sqrt(dev / no);
		// check first point
		// if (false && P0 >0.0 && Math.abs(ys[0]) > P0)

		if (PMAX > 0.0 && Math.abs(ys[0]) > PMAX) {
			//System.out.println("Abs(y0)>Pref.PMAX");
			flags[0] = false;
		}

		no = xs.length - 1;
		ArrayList list = new ArrayList();
		boolean prevVFlag = false;

		for (int i = 1; i < xs.length; i++){

			// check against absolute value
			if (PMAX > 0.0 && Math.abs(ys[i]) > PMAX){
				
				System.out.println("Abs(y)>Pref.PMAX");
				list.add(String.valueOf(i) + ":" + String.valueOf(ys[i]));
				
				continue;
			}

			// check against deviation from mean
			double ratio = Math.abs(ys[i] - mean) / dev;

			if (TR > 0.0 && ratio > TR) {
				System.out.println("Deviation from mean > Pref.TR");
				list.add(String.valueOf(i) + ":" + String.valueOf(ys[i]));

				continue;
			}

			// shape like v

			boolean vFlag = false;
			if (i < no && (ys[i] - ys[i - 1]) * (ys[i] - ys[i + 1]) > 0.0)
				vFlag = true;

			// check against angle and delta-the difference from the extrapolated line
			double t1 = 0.0;
			
			if (i > 1)
				t1 = calcTheta(xs, ys, i - 2, i - 1, i);

			double t2 = 0.0;

			if (i < no)
				t2 = calcTheta(xs, ys, i - 1, i, i + 1);

			double t = 0.0;

			if (vFlag == false)
				t = Math.max(t1, t2);
			else
				t = Math.min(t1, t2);

			double delta0 = Math.abs(ys[i] - ys[i - 1]);

			double delta1 = 10000000.0;
			double delta2 = 10000000.0;

			if (i == no){
				delta1 = Math.abs(ys[no] - (xs[no] - xs[no - 2]) * (ys[no - 1] - ys[no - 2]) / (xs[no - 1] - xs[no - 2])
						- ys[no - 2]);

			}else{

				delta1 = Math.abs(ys[i] - (xs[i] - xs[i - 1]) * (ys[i + 1] - ys[i - 1]) / (xs[i + 1] - xs[i - 1]) - ys[i - 1]);
				delta2 = Math.abs(ys[i] - ys[i + 1]);
			}

			double delta = 1000000.0;

			delta = Math.min(delta0, delta1);
			delta = Math.min(delta2, delta);

			System.out.println(i + ":" + ":" + delta0 + ":" + delta1 + ":" + delta2 + ":" + t1 + ":" + t2);

			if (prevVFlag && i == no){

				if (delta > TPN && TPN > 0.0){
					System.out.println(delta + ":" + TPN + ":" + "Last point's delta > Pref.TPN");
					list.add(String.valueOf(i) + ":" + String.valueOf(ys[i]));
				}

				continue;
			}

			if (i != no && delta > TP && TP > 0.0){
				System.out.println(delta + ":" + TP + ":" + "Current point's delta > Pref.TP");
				list.add(String.valueOf(i) + ":" + String.valueOf(ys[i]));

				continue;
			}else if (vFlag && t < THETA && delta > TPV && THETA > 0.0 && TPV > 0.0){

				System.out.println(i + ":" + t + ":" + THETA + ":" + delta + ":" + TPN + ":"
							+ "Current point's delta > Pref.TPV");

				list.add(String.valueOf(i) + ":" + String.valueOf(ys[i]));
			}

			if (ds != null){

				if (delta > 499)
					delta = 499;

				ds[i][(int) delta]++;
			}

			prevVFlag = vFlag;
		}

		no = xs.length / MASK;
		for (int i = 0; i < no; i++){

			if (list.size() == 0)
				break;

			int maxj = -1;
			int maxl = -1;
			double maxd = -100000.0;

			for (int j = 0; j < list.size(); j++){

				String line = (String) list.get(j);
				int k = line.indexOf(":");
				int l = Integer.parseInt(line.substring(0, k));
				double d = Double.parseDouble(line.substring(k + 1));

				if (d > maxd){
					maxd = d;
					maxl = l;
					maxj = j;
				}
			}

			flags[maxl] = false;
			list.remove(maxj);
		}
		// System.out.println(i+":"+mean+":"+dev+":"+ys[i]+":"+(ys[i]-mean)+":"+ratio+":"+ratio2+":"+flag2[i]);

		return flags;
	}


	public static boolean[] maskDif3(double[] xs, double[] ys, int[][] ds, double[] ms, HashMap<String, String>[] maps) {

		boolean[] flags = new boolean[xs.length];

		for (int i = 0; i < xs.length; i++)
			flags[i] = true;

		if (xs.length < 3)
			return flags;

		// calc mean and std

		int no = 0;
		double mean = 0.0;

		for (int i = 0; i < xs.length; i++){
			mean += ys[i];
			no++;
		}

		mean /= no;
		double dev = 0.0;

		for (int i = 0; i < xs.length; i++){
			double delta = ys[i] - mean;
			dev += delta * delta;
		}

		dev = Math.sqrt(dev / no);

		// check first point

		// if (false && P0 >0.0 && Math.abs(ys[0]) > P0)
		if (PMAX > 0.0 && Math.abs(ys[0]) > PMAX) {
			//System.out.println("Abs(y0)>Pref.PMAX");

			flags[0] = false;

			if (maps != null) 
				maps[0] = maskInfoMap(xs[0], ys[0], ys[0], "abs(y0)>Pref.PMAX:" + String.valueOf(ys[0]) + " " + String.valueOf(PMAX));
		}

		no = xs.length - 1;

		ArrayList list = new ArrayList();

		boolean prevVFlag = false;

		for (int i = 1; i < xs.length; i++) {

			// check against absolute value

			if (PMAX > 0.0 && Math.abs(ys[i]) > PMAX) {
				//System.out.println("Abs(y)>Pref.PMAX");

				list.add(String.valueOf(i) + ":" + String.valueOf(ys[i]));

				if (maps != null) 
					maps[i] = maskInfoMap(xs[i], ys[i], ys[i], "abs(y)>Pref.PMAX:" + String.valueOf(ys[i]) + " " + String.valueOf(PMAX)); 

				continue;

			}
			// check against deviation from mean

			double ratio = Math.abs(ys[i] - mean) / dev;

			if (TR > 0.0 && ratio > TR) {
				//System.out.println("(resp-mean)/std > Pref.TR");
				list.add(String.valueOf(i) + ":" + String.valueOf(ys[i]));

				if (maps != null) 
					maps[i] = maskInfoMap(xs[i], ys[i], ratio, "(resp-mean)/std > Pref.TR:" + String.valueOf(ratio) + " " + String.valueOf(TR)); 
				
				continue;
			}

			// shape like v
			boolean vFlag = false;
			if (i < no && (ys[i] - ys[i - 1]) * (ys[i] - ys[i + 1]) > 0.0)
				vFlag = true;

			// check against angle and delta-the difference from the extrapolated line

			double t1 = 0.0;

			if (i > 1)
				t1 = calcTheta(xs, ys, i - 2, i - 1, i);

			double t2 = 0.0;

			if (i < no)
				t2 = calcTheta(xs, ys, i - 1, i, i + 1);

			double t = 0.0;

			if (vFlag == false)
				t = Math.max(t1, t2);
			else
				t = Math.min(t1, t2);

			double delta0 = Math.abs(ys[i] - ys[i - 1]);
			double delta1 = 10000000.0;
			double delta2 = 10000000.0;

			if (i == no)
				delta1 = Math.abs(ys[no] - (xs[no] - xs[no - 2]) * (ys[no - 1] - ys[no - 2]) / (xs[no - 1] - xs[no - 2])- ys[no - 2]);

			else{
				delta1 = Math.abs(ys[i] - (xs[i] - xs[i - 1]) * (ys[i + 1] - ys[i - 1]) / (xs[i + 1] - xs[i - 1]) - ys[i - 1]);
				delta2 = Math.abs(ys[i] - ys[i + 1]);
			}

			double delta = 1000000.0;
			delta = Math.min(delta0, delta1);
			delta = Math.min(delta2, delta);

			//System.out.println(i + ":" + ":" + delta0 + ":" + delta1 + ":" + delta2 + ":" + t1 + ":" + t2);

			if (prevVFlag && i == no) {
				if (delta > TPN && TPN > 0.0) {
					//System.out.println(delta + ":" + TPN + ":" + "Last point's delta > Pref.TPN");

					list.add(String.valueOf(i) + ":" + String.valueOf(ys[i]));

					if (maps != null)
						maps[i] = maskInfoMap(xs[i], ys[i], delta, "Last point's delta distance > Pref.TPN:" + String.valueOf(delta) + " "+ String.valueOf(TPN));
				}

				continue;
			}

			if (i != no && delta > TP && TP > 0.0){
				//System.out.println(delta + ":" + TP + ":" + "Current point's delta > Pref.TP");

				list.add(String.valueOf(i) + ":" + String.valueOf(ys[i]));

				if (maps != null)
					maps[i] = maskInfoMap(xs[i], ys[i], delta, "Current point's delta distance > Pref.TP:" + String.valueOf(delta) + " "+ String.valueOf(TP));

				continue;

			}else if (vFlag && t < THETA && delta > TPV && THETA > 0.0 && TPV > 0.0){
				//System.out.println(i + ":" + t + ":" + THETA + ":" + delta + ":" + TPN + ":"+ "Current point's delta > Pref.TPV");

				list.add(String.valueOf(i) + ":" + String.valueOf(ys[i]));

				if (maps != null) 
					maps[i] = maskInfoMap(xs[i], ys[i], t, "Current point's delta angle > Pref.TPV:" + String.valueOf(t) + " "+ String.valueOf(THETA));
			}

			if (ds != null){
				if (delta > 499)
					delta = 499;

				ds[i][(int) delta]++;
			}

			prevVFlag = vFlag;
		}

		no = xs.length / MASK;

		for (int i = 0; i < no; i++){

			if (list.size() == 0)
				break;

			int maxj = -1;
			int maxl = -1;
			double maxd = -100000.0;

			for (int j = 0; j < list.size(); j++){

				String line = (String) list.get(j);
				int k = line.indexOf(":");
				int l = Integer.parseInt(line.substring(0, k));
				double d = Double.parseDouble(line.substring(k + 1));

				if (d > maxd){
					maxd = d;
					maxl = l;
					maxj = j;
				}
			}

			flags[maxl] = false;
			list.remove(maxj);
		}

		for (int i = 0; i < xs.length; i++) {
			if (flags[i] && maps != null)
				maps[i] = null;
		}
		// System.out.println(i+":"+mean+":"+dev+":"+ys[i]+":"+(ys[i]-mean)+":"+ratio+":"+ratio2+":"+flag2[i]);
		return flags;
	}
	
}
