/*
 * HillStat.java
 *
 * Created on February 21, 2007, 11:12 AM
 *
 * To change this template, choose Tools | Template Manager
 * and open the template in the editor.
 */

package gov.nih.ncats.ifx.qhts.utils;
 
import java.util.ArrayList;
import org.apache.commons.math3.distribution.FDistribution;

/**
 *
 * @author yw84380
 */
public class HillStat
{
    public static double LN10 = Math.log(10.0);
    public static double DELTA = 0.1;

    public static boolean debug = false;

    public static void main(String[] args)
    {


        double f = Stat.calcF(5.0, 5.0, 10.967);


        System.out.println(f);
    }


    public HillStat()
    {

    }


    //We dont use masked points at high concentrations, i.e.
    //we would only include masked points at concentrations lower than the first unmasked point.
    //This way the p-values will make more sense for bell-shaped curves.
    public static double calcPValue(double y0min, double ymin, double xmin, double smin, double[] xs, double[] ys, boolean[] flags)
    {
        boolean[] flags2 = null;
        int no = xs.length;

        if (HillConstants.USE_MASKED_POINT_OPTION.equals("No for All"))
        {
            flags2 = flags;

            no = 0;
            for (int i=0; i<flags.length; i++)
            {
                if (flags[i])
                    no++;
            }
        }
        else if (HillConstants.USE_MASKED_POINT_OPTION.equals("No for Bell"))
        {
            //find first unmasked point
            int k = -1;
            for (int i=flags.length-1; i>=0; i--)
            {
                if (flags[i])
                {
                    k = i;
                    break;
                }
            }

            if (k > 0)
            {
                no = k+1;

                flags2 = new boolean[flags.length];

                for (int i=0; i<=k; i++)
                    flags2[i] = true;

                for (int i=k+1; i<flags.length; i++)
                    flags2[i] = false;
            }
        }

        if (no <= 3)
            return 1.0;

        double bestSS = calcHillDeviation(xmin, y0min, ymin, smin, flags2, null, xs, ys);

        double flatSS = calcConstantDeviation(flags2, null, xs, ys );

        //F test
        double ratio = (flatSS/(no-1))/(bestSS/(no-3)); //should be -4; to be user friendly haha

        double pvalue = Stat.calcF((double)(no-1), (double)(no-3), ratio);

        return pvalue;
     }




    public static double calcHillDeviation(double x05,double y0, double yinf,double slope,double sym,
                    boolean[] flags,double[] ws,double[] xs, double[] ys)
    {
        double dev = 0.0;

        for (int i=0; i<xs.length; i++)
        if (flags == null || flags[i])
        {

                double y = y0+(yinf-y0)/Math.pow(1.0+Math.exp(LN10*slope*(x05-xs[i])), sym);

                double delta = (y-ys[i]);



                if (ws != null)

                        dev += ws[i]*delta*delta;

                else

                        dev += delta*delta;

        }



        return dev;

    }


    public static double calcPValue(double y0min, double ymin, double xmin, double smin, double sym, double[] xs, double[] ys, boolean[] flags)
    {
        int no = xs.length;

        if (no <= 3)
            return 1.0;

        double bestSS = calcHillDeviation(xmin, y0min, ymin, smin, sym, flags, null, xs, ys);

        double flatSS = calcConstantDeviation(flags, null, xs, ys );

        //F test
        double ratio = (flatSS/(no-1))/(bestSS/(no-3)); //should be -4; to be user friendly haha

        double pvalue = Stat.calcF((double)(no-1), (double)(no-3), ratio);

        return pvalue;
    }


    public static double calcSS(ArrayList list)
    {
        double ss = 0.0;

        for (int i=0; i<list.size(); i++)
        {
            double td = Double.parseDouble((String)list.get(i));

            ss += td*td;
        }

        return ss;
    }






    public static double calcPValueKW(double y0min, double ymin, double xmin, double smin, double[] xs, double[] ys, boolean[] flags)
    {
        double[] hs = calcHillDeviations(xmin, y0min, ymin, smin, flags, null, xs, ys);
        double[] fs = calcConstantDeviations(flags, null, xs, ys);


        double[][] xys = {xs, ys};

        double p = Stat.kwRankTest(xys);

        return p;
     }



    public static boolean calcPValueMWU(double y0min, double ymin, double xmin, double smin, double[] xs, double[] ys, boolean[] flags)
    {
        double[] hs = calcHillDeviations(xmin, y0min, ymin, smin, flags, null, xs, ys);
        double[] fs = calcConstantDeviations(flags, null, xs, ys);

        return Stat.mwuTest(hs, fs);
    }



    public static double[] calcCI(boolean p4Flag, double y0min, double ymin, double xmin, double smin, double p, double[] xs, double[] ys, boolean[] flags)
    {
        p4Flag = false;

        double dfp = 3.0;
        if (p4Flag)
            dfp = 4.0;

        double dfnp = xs.length-dfp;

        double bestSS = calcHillDeviation(xmin, y0min, ymin, smin, flags, null, xs, ys);

        if (dfnp <= 0)
        {
            double[] values3 ={y0min, y0min, ymin, ymin, xmin, xmin, smin,smin, bestSS, -bestSS};
            return values3;
        }

        FDistribution fd = new FDistribution(dfp, dfnp);

        //inverse F
        double f = fd.inverseCumulativeProbability(1.0-p);


        double allowedSS = bestSS*(f*dfp/dfnp+1);

        //first estimate
        double y0l = y0min-Math.max(0.5*y0min, 20.0);
        double y0r = y0min+Math.max(0.5*y0min, 20.0);
        double y0d = 2.0;

        if (p4Flag == false)
        {
            y0l = y0min;
            y0r = y0min;
        }

        double yl = ymin-Math.max(0.5*ymin, 20.0);
        double yr = ymin+Math.max(0.5*ymin, 20.0);
        double yd = 2.0;

        double xl = xmin-4.0;
        double xr = xmin+4.0;
        double xd = 0.1;

        double sl = smin-2.0;
        if (sl < 0.2)
            sl = 0.2;

        double sr = smin+2.0;
        double sd = 0.2;

        if (debug)
            System.out.println(p4Flag+"\t"+bestSS+"\t"+allowedSS);

        double[] values = calcCI(bestSS, allowedSS, p4Flag, y0min, ymin, xmin, smin, y0l, y0r, y0d, yl, yr, yd, xl, xr, xd, sl, sr, sd, flags, xs, ys);

        if (values == null)
        {
            double[] values3 ={y0min,y0min,ymin,ymin,xmin,xmin,smin,smin,bestSS,allowedSS};
            return values3;
        }

        if (debug)
        {
            String info = Util.formatDouble(y0min, 2, 2)+"\t";
            info += Util.formatDouble(values[0], 2, 2)+"\t";
            info += Util.formatDouble(values[1], 2, 2)+"\t";
            info += Util.formatDouble(ymin, 2, 2)+"\t";
            info += Util.formatDouble(values[2], 2, 2)+"\t";
            info += Util.formatDouble(values[3], 2, 2)+"\t";
            info += Util.formatDouble(xmin, 2, 2)+"\t";
            info += Util.formatDouble(values[4], 2, 2)+"\t";
            info += Util.formatDouble(values[5], 2, 2)+"\t";
            info += Util.formatDouble(smin, 2, 2)+"\t";
            info += Util.formatDouble(values[6], 2, 2)+"\t";
            info += Util.formatDouble(values[7], 2, 2);
            System.out.println(info);
        }

        return values;
    }


    public static double[] calcCI(double bestSS, double allowedSS, boolean p4Flag, double y0min, double ymin, double xmin, double smin,
            double y0l, double y0r, double y0d,
            double yl,double yr, double yd,
            double xl,double xr,double xd,
            double sl,double sr,double sd,
            boolean[] flags, double[] xs, double[] ys)
    {
            double y0_left = 0.0;
            double yinf_left = 0.0;
            double slope_left = 0.0;
            double x05_left = 0.0;

            double y0_right = 0.0;
            double yinf_right = 0.0;
            double slope_right = 0.0;
            double x05_right = 0.0;

            boolean flag = false;

            double maxSS = 0.0;

            for (double yinf=yl; yinf<=yr; yinf+=yd)
            {
                for (double x05=xl; x05<=xr; x05+=xd)
                {
                    for (double slope=sl;slope<=sr;slope+=sd)
                    {
                        double y0 = y0l;

                        while (y0 <= y0r)
                        {
                            double dev = calcHillDeviation(x05,y0,yinf,slope,flags,null,xs,ys);

                            if (dev > maxSS)
                                maxSS = dev;

                            if (Math.abs(dev-allowedSS)<DELTA*allowedSS)
                            {
                                flag = true;

                                if (y0 > y0min && y0-y0min > y0_right)
                                    y0_right = y0-y0min;

                                if (y0 < y0min && y0min-y0 > y0_left)
                                    y0_left = y0min-y0;

                                if (yinf > ymin && yinf-ymin > yinf_right)
                                    yinf_right = yinf-ymin;

                                if (yinf < ymin && ymin-yinf > yinf_left)
                                    yinf_left = ymin-yinf;

                                if (x05 > xmin && x05-xmin > x05_right)
                                    x05_right = x05-xmin;

                                if (x05 < xmin && xmin-x05 > x05_left)
                                    x05_left = xmin-x05;

                                if (slope > smin && slope-smin > slope_right)
                                    slope_right = slope-smin;

                                if (slope < smin && smin-slope > slope_left)
                                    slope_left = smin-slope;
                            }

                            y0 = y0+y0d;

                            if (p4Flag == false)
                                break;
                        }
                    }
                }
            }

            if (debug)
                System.out.println(maxSS);

            if (flag == false)
                return null;

            double[] values ={y0min-y0_left, y0min+y0_right, ymin-yinf_left, ymin+yinf_right, xmin-x05_left,xmin+x05_right, smin-slope_left,smin+slope_right, bestSS, allowedSS};

            return values;
    }



    public static double calcHillDeviation(double x05,double y0, double yinf,double slope,
                    boolean[] flags,double[] ws,double[] xs, double[] ys)
    {
            double dev = 0.0;

            for (int i=0; i<xs.length; i++)
            if (flags == null || flags[i])
            {
                    double y = y0+(yinf-y0)/(1.0+Math.exp(LN10*slope*(x05-xs[i])));
                    double delta = (y-ys[i]);

                    if (ws != null)
                            dev += ws[i]*delta*delta;
                    else
                            dev += delta*delta;
            }

            return dev;
    }



    public static double[] calcHillDeviations(double x05,double y0, double yinf,double slope,
                    boolean[] flags,double[] ws,double[] xs, double[] ys)
    {
            double[] ds = new double[xs.length];

            for (int i=0; i<xs.length; i++)
            if (flags == null || flags[i])
            {
                    double y = y0+(yinf-y0)/(1.0+Math.exp(LN10*slope*(x05-xs[i])));
                    double delta = (y-ys[i]);

                    ds[i] = Math.abs(delta);
            }

            return ds;
    }



    public static double[] calcCI1(boolean p4Flag, double y0min, double ymin, double xmin, double smin, double p, double[] xs, double[] ys)
    {
        p4Flag = false;

        double dfp = 3.0;
        if (p4Flag)
            dfp = 4.0;

        double dfnp = xs.length-dfp;

        FDistribution fd = new FDistribution(dfp, dfnp);

        //inverse F
        double f = fd.inverseCumulativeProbability(1.0-p);

        double bestSS = calcHillDeviation(xmin, y0min, ymin, smin, null, null, xs, ys);

        double allowedSS = bestSS*(f*dfp/dfnp+1);

        //first estimate
        double y0l = y0min-20.0;
        if (y0l < 0.0)
            y0l = 0.0;

        double y0r = y0min+20.0;
        double y0d = 5.0;

        if (p4Flag == false)
        {
            y0l = y0min;
            y0r = y0min;
        }

        double yl = ymin-20.0;
        if (yl < 0.0)
            yl = 0.0;

        double yr = ymin+20.0;
        double yd = 5.0;

        double xl = xmin-2.0;
        double xr = xmin+2.0;
        double xd = 0.5;

        double sl = smin-2.0;
        if (sl < 0.1)
            sl = 0.1;

        double sr = smin+2.0;
        double sd = 0.5;

        double coef = 1.0;

        double[] values = null;

        System.out.println(p4Flag+"\t"+bestSS+"\t"+allowedSS);

        boolean noFlag = false;

        while (true)
        {
            double coef2 = 1.0+coef*0.1;

            values = calcCI(bestSS, allowedSS, p4Flag, y0min, ymin, xmin, smin, coef2*y0l, coef2*y0r, y0d/coef, coef2*yl, coef2*yr, yd/coef, coef2*xl, coef2*xr, xd/coef, sl, sr, sd, null, xs, ys);

            if (values == null)
            {
                coef = coef +1.0;

                System.out.println(coef);

                if (coef > 10.0)
                {
                    noFlag = true;
                    break;
                }
                continue;
            }

            String info = Util.formatDouble(y0min, 2, 2)+"\t";
            info += Util.formatDouble(values[0], 2, 2)+"\t";
            info += Util.formatDouble(values[1], 2, 2)+"\t";
            info += Util.formatDouble(ymin, 2, 2)+"\t";
            info += Util.formatDouble(values[2], 2, 2)+"\t";
            info += Util.formatDouble(values[3], 2, 2)+"\t";
            info += Util.formatDouble(xmin, 2, 2)+"\t";
            info += Util.formatDouble(values[4], 2, 2)+"\t";
            info += Util.formatDouble(values[5], 2, 2)+"\t";
            info += Util.formatDouble(smin, 2, 2)+"\t";
            info += Util.formatDouble(values[6], 2, 2)+"\t";
            info += Util.formatDouble(values[7], 2, 2);
            System.out.println(info);


            break;
        }

        if (noFlag)
        {
            double[] values3 ={y0min, y0min, ymin, ymin, xmin, xmin, smin,smin, 0.0, 0.0};
            return values3;
        }

        //refine
        y0l = values[0];
        y0r = values[1];
        y0d = 0.5;

        yl = values[2];
        yr = values[3];
        yd = 0.5;

        xl = values[4];
        xr = values[5];
        xd = 0.05;

        sl = values[6];
        sr = values[7];
        sd = 0.1;

        double[] values2 = calcCI(bestSS, allowedSS, p4Flag, y0min, ymin, xmin, smin, y0l, y0r, y0d, yl, yr, yd, xl, xr, xd, sl, sr, sd, null, xs, ys);

        if (values2 != null)
            values = values2;

        String info = Util.formatDouble(y0min, 2, 2)+"\t";
        info += Util.formatDouble(values[0], 2, 2)+"\t";
        info += Util.formatDouble(values[1], 2, 2)+"\t";
        info += Util.formatDouble(ymin, 2, 2)+"\t";
        info += Util.formatDouble(values[2], 2, 2)+"\t";
        info += Util.formatDouble(values[3], 2, 2)+"\t";
        info += Util.formatDouble(xmin, 2, 2)+"\t";
        info += Util.formatDouble(values[4], 2, 2)+"\t";
        info += Util.formatDouble(values[5], 2, 2)+"\t";
        info += Util.formatDouble(smin, 2, 2)+"\t";
        info += Util.formatDouble(values[6], 2, 2)+"\t";
        info += Util.formatDouble(values[7], 2, 2);
        System.out.println(info);



        return values;
    }


    	public static double calcConstantDeviation(boolean[] flags,double[] ws,double[] xs, double[] ys)
	{
		double mean = 0.0;

		int no = 0;

		for (int i=0; i<xs.length; i++)
		if (flags == null || flags[i])
		{
			mean += ys[i];

			no++;
		}

		mean /= no;

		double dev = 0.0;


		for (int i=0; i<xs.length; i++)
		if (flags == null || flags[i])
		{
			double delta = (ys[i]-mean);

			if (ws != null)
				dev += ws[i]*delta*delta;
			else
				dev += delta*delta;
		}


		return dev;

	}



        public static double[] calcConstantDeviations(boolean[] flags,double[] ws,double[] xs, double[] ys)
	{
		double mean = 0.0;

		int no = 0;

		for (int i=0; i<xs.length; i++)
		if (flags == null || flags[i])
		{
			mean += ys[i];

			no++;
		}

		mean /= no;


                double[] ds = new double[xs.length];

		for (int i=0; i<xs.length; i++)
		if (flags == null || flags[i])
		{
			double delta = (ys[i]-mean);
                        ds[i] = Math.abs(delta);
		}


		return ds;

	}

}
