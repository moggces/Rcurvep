/*

	parse sw.xml and put them into oracle database
 */

package  gov.nih.ncats.ifx.qhts.utils;

import org.apache.commons.math3.distribution.ChiSquaredDistribution;

public class Stat
	extends Object
{

        private static  double[] cof  = null;
        private static double stp= 2.50662827465;
        private static double half = 0.5;
        private static double one = 1.0;
        private static double fpf = 5.5;


        private static int[][] mwuTable;

      static
      {
          cof = new double[6];
          cof[0] = 76.18009173;
          cof[1] = -86.50532033;
          cof[2] = 24.014909822;
          cof[3] = -1.231739516;
          cof[4] = 0.00120858003;
          cof[5] = -0.536382e-5;



          String mwuText = "2-8 0 0 0 0 1 1 1 1 1 2 2 2 2;";
          mwuText += "3-5 0 1 1 2 2 3 3 4 4 5 5 6 6 7 7 8;";
          mwuText += "4-7 3 4 4 5 6 7 8 9 10 11 11 12 13 13;";
          mwuText += "5-2 0 1 2 2 3 5 6 7 8 9 10 12 13 14 15 17 18 19 20;";
          mwuText += "6-6 5 6 8 10 11 13 14 16 17 19 21 22 24 25 27;";
          mwuText += "7-7 8 10 12 14 16 18 20 22 24 26 28 30 32 34;";
          mwuText += "8-8 13 15 17 19 22 24 26 29 31 34 36 38 41;";
          mwuText += "9-9 17 20 23 26 28 31 34 37 39 42 45 48;";
          mwuText += "10-10 23 26 29 33 36 39 42 45 48 52 55;";
          mwuText += "11-11 30 33 37 40 44 47 51 55 58 62;";
          mwuText += "12-12 37 41 45 49 53 57 61 65 69;";
          mwuText += "13-13 45 50 54 59 63 67 72 76;";
          mwuText += "14-14 55 59 64 67 74 78 83;";
          mwuText += "15-15 64 70 75 80 85 90;";
          mwuText += "16-16 75 81 86 92 98;";
          mwuText += "17-17 87 93 99 105;";
          mwuText += "18-18 99 106 112;";
          mwuText += "19-19 113 119;";
          mwuText += "20-20 127";

          mwuTable = new int[20][20];
          String[] lines = mwuText.split(";");
          for (int i=0; i<lines.length; i++)
          {
              String[] cols = lines[i].split(" ");
              String[] toks = cols[0].split("-");

              int r = Integer.parseInt(toks[0])-1;
              int c = Integer.parseInt(toks[1])-1;

              for (int j=0; j<cols.length-1; j++)
              {
                  mwuTable[r][c+j] = Integer.parseInt(cols[j+1]);
              }
          }

          for (int i=0; i<20; i++)
          for (int j=0; j<i; j++)
          {
              mwuTable[i][j] = mwuTable[j][i];
          }
      }



        public static void main(String[] args)
        {


            double[] xs = {28,27,33,23,24,17,25,23,31,23,23,22,15};
            double[] ys = {14,20,16,13,18,21,23,20,14,20,20,14,11,16,13};

            double[][] xys = {xs, ys};

            double p = Stat.kwRankTest(xys);

            System.out.println(p);
        }



        public static boolean mwuTest(double[] xs, double[] ys)
        {
            int n1 = xs.length;
            int n2 = ys.length;

            //convert xs into ranks
            int[] rxs = new int[n1];
            int[] rys = new int[n2];

            int r = 1;

            while (true)
            {
                double min = Double.MAX_VALUE;

                int s = -1;
                int p = -1;

                for (int i=0; i<n1; i++)
                {
                    if (rxs[i] == 0)
                    {
                        if (xs[i] <= min)
                        {
                            min = xs[i];
                            s = 1;
                            p = i;
                        }
                    }
                }

                for (int i=0; i<n2; i++)
                {
                    if (rys[i] == 0)
                    {
                       if (ys[i] <= min)
                        {
                            min = ys[i];
                            s = 2;
                            p = i;
                        }
                    }
                }

                if (s < 0)
                    break;

                if (s == 1)
                    rxs[p] = r;
                else
                    rys[p] = r;

                r++;
            }

            //average r1
            int sumr = 0;
            for (int i=0; i<rxs.length; i++)
            {
                sumr += rxs[i];
                //System.out.println(rxs[i]+":"+sumr);
            }

            //System.out.println(sumr);


            int u1 = n1*n2+n1*(n1+1)/2-sumr;
            int u2 = n1*n2-u1;

            int u = Math.min(u1, u2);

            System.out.println(u+":"+mwuTable[n1-1][n2-1]);

            if (u < mwuTable[n1-1][n2-1])
                return true;
            else
                return false;

        }


        public static double kwRankTest(double[][] xs)
        {
            //total sample no
            int n = 0;
            for (int i=0; i<xs.length; i++)
                n+= xs[i].length;

            //total average rank
            double meanR = 0.5*(n+1);


            //convert xs into ranks
            int[][] rs = new int[xs.length][];
            for (int i=0; i<xs.length; i++)
                 rs[i] = new int[xs[i].length];

            int r = 1;

            while (true)
            {
                int s = -1;
                int p = -1;

                double min = Double.MAX_VALUE;

                for (int i=0; i<xs.length; i++)
                for (int j=0; j<xs[i].length; j++)
                {
                    if (rs[i][j] == 0)
                    {
                        if (xs[i][j] <= min)
                        {
                            min = xs[i][j];
                            s = i;
                            p = j;
                        }
                    }
                }

                if (s < 0)
                    break;

                rs[s][p] = r;

                r++;
            }

            //average r of each group
            double[] meanRs = new double[rs.length];

            for (int i=0; i<rs.length; i++)
            {
                meanRs[i] = 0.0;
                for (int j=0; j<rs[i].length; j++)
                    meanRs[i] += rs[i][j];

                meanRs[i] /= rs[i].length;
            }

            //k test
            double dk = 0.0;

            for (int i=0; i<rs.length; i++)
            {
                double delta = meanRs[i]-meanR;
                dk += rs[i].length*delta*delta;
            }
            dk *= 12.0/(n*(n+1));

            ChiSquaredDistribution cd = new ChiSquaredDistribution(xs.length-1);

            double p = cd.cumulativeProbability(dk);
            //double p = cd.probability(dk);

            return 1-p;
        }


	public static double[] getP(double td)
	{
		double[] ps = new double[2];

		if (td < 1.642)
		{
			ps[0] = 0.2;
			ps[1] = 0.2;
		}
		else if (td >= 1.642 && td < 2.706)
		{
			ps[0] = 0.1;
			ps[1] = 0.2;
		}
		else if (td >= 2.706 && td < 3.841)
		{
			ps[0] = 0.05;
			ps[1] = 0.1;
		}
		else if (td >= 3.841 && td < 5.412)
		{
			ps[0] = 0.02;
			ps[1] = 0.05;
		}
		else if (td >= 5.412 && td < 6.635)
		{
			ps[0] = 0.01;
			ps[1] = 0.02;
		}
		else if (td >= 6.635 && td < 10.827)
		{
			ps[0] = 0.001;
			ps[1] = 0.01;
		}
		else
		{
			ps[0] = 0.001;
			ps[1] = 0.001;
		}

		return ps;
	}


	public static double[] linecorr(double[] x, double[] y)
		throws Exception
	{
		return linecorr(x, y, 0);
	}



	public static double[] linecorr(double[] x, double[] y, int delta)
		throws Exception
	{
		int no = x.length-delta;

		double xm = 0;
		double ym = 0;

		for (int i=delta; i<no; i++)
		{
			xm += x[i];
			ym += y[i-delta];
		}

		xm /= no;
		ym /= no;

		double lxx = 0;
		double lyy = 0;
		double lxy = 0;
		for (int i=delta; i<no; i++)
		{
			lxx += (x[i]-xm)*(x[i]-xm);
			lyy += (y[i-delta]-ym)*(y[i-delta]-ym);
			lxy += (x[i]-xm)*(y[i-delta]-ym);
		}

		double[] values = new double[3];
		values[0] = lxy/Math.sqrt(lxx*lyy);
		values[1] = lxy/lxx;
		values[2] = ym-values[1]*xm;

		return values;
	}


	public static boolean[] data_mask(double[] x, double[] y)
		throws Exception
	{
		int no = x.length;

		boolean[] flags = new boolean[no];

		if (x.length == 2)
		{
			for (int i=0; i<no; i++)
				flags[i] = true;

			return flags;
		}

		double[] grad = new double[no];

		for (int i=0; i<no-1; i++)
		{
			grad[i] = (y[i+1]-y[i])/(x[i+1]-x[i]);
		}

		grad[no-1]= grad[no-2];

		double[] sta = calcMeanStd(grad);

		for (int i=0; i<no; i++)
		{
			flags[i] = true;

			if (Math.abs(grad[i]-sta[0])>3*sta[1])
				flags[i] = false;

			System.out.println(sta[0]+":"+sta[1]+":"+grad[i]+":"+Math.abs(grad[i]-sta[0])/sta[1]);

			if (flags[i] == false)
			{
				System.out.println(sta[0]+":"+sta[1]+":"+grad[i]+":"+Math.abs(grad[i]-sta[0])/sta[1]);
				System.exit(1);
			}
		}

		return flags;
	}


	public static double[] linecorr(double[] x, double[] y, boolean flag)
		throws Exception
	{
		int no = x.length;

		boolean[] flags = null ; //data_mask(x, y);

		double xm = 0;
		double ym = 0;

		no = 0;
		for (int i=0; i<x.length; i++)
		{
			if (flags != null && flags[i] == false)
				continue;

			xm += x[i];
			ym += y[i];

			no++;
		}

		xm /= no;
		ym /= no;

		double lxx = 0;
		double lyy = 0;
		double lxy = 0;

		for (int i=0; i<x.length; i++)
		{
			if (flags != null && flags[i] == false)
				continue;

			lxx += (x[i]-xm)*(x[i]-xm);
			lyy += (y[i]-ym)*(y[i]-ym);
			lxy += (x[i]-xm)*(y[i]-ym);
		}

		double[] values = new double[3];
		values[0] = lxy/Math.sqrt(lxx*lyy);
		values[1] = lxy/lxx;
		values[2] = ym-values[1]*xm;

		return values;
	}



	public static double[] linefit(double[] fits, double[] x)
		throws Exception
	{
		double[] y = new double[x.length];

		for (int i=0; i<x.length; i++)
		{
			y[i] = x[i]*fits[1]+fits[2];
		}

		return y;
	}


	public static double linefit(double[] fits, double x)
		throws Exception
	{
		return x*fits[1]+fits[2];
	}


	public static double[] calcMeanStd(double[] x)
		throws Exception
	{
		return calcMeanStd(x, x.length);
	}


	public static double[] calcMeanStd(double[] x, int no)
		throws Exception
	{
		if (no<0 || no>x.length)
			no = x.length;

		double xm = 0;
		for (int i=0; i<no; i++)
			xm += x[i];

		xm /= no;

		double std = 0.0;
		for (int i=0; i<no; i++)
			std += (x[i]-xm)*(x[i]-xm);

		std = Math.sqrt(std/no);

		double[] values = new double[2];
		values[0] = xm;
		values[1] = std;

		return values;
	}



        public static double sigOfCorr(double r, int n)
        {
            double tiny = 1.0e-6;

            double div_1 = 0.0;
            double div_2 = 0.0;
            double df = 0.0;
            double tt = 0.0;
            double pr = 0.0;

            if( n > 2 )
            {
                // Fishers z transformation of pcc

                div_1 = ( (1.0+r)+tiny );
                div_2 = ( (1.0-r)+tiny );

                df = n-2;
                tt = r*Math.sqrt( df/( div_1*div_2 ) );
                pr = betai( 0.5*df, 0.5, df/(df+tt*tt) );
          }
          else
          {
            // n<=2
                pr = 1.00;
          }

        return pr;
   }


    // incomplete beta function from numerical recepeies
    public static double betai(double a, double b, double x)
    {
      double bt = 0.0;
      double pr = 0.0;

      if ( x < 0.0 || x > 1.0)
      {
        return 1.0;
      }

      if ( x==0.0 || x==1.0)
      {
         bt = 0.0;
      }
      else
      {
         bt = Math.exp(gammln(a+b)-gammln(a)-gammln(b) + a*Math.log(x)+b*Math.log(1.0-x));
      }

      if( x < (a+1.0)/(a+b+2.0) )
        pr = bt*betacf(a,b,x)/a;
      else
        pr = 1.0-bt*betacf(b,a,1.0-x)/b;

      return pr;
   }


    // continude fraction fo incomplete beta function
   public static double betacf(double a, double b, double x)
   {
      int itmax= 100;
      double eps = 3.0e-7;
      double am, bm, az, bz, qab, qap, qam;
      int em, tem;
      double ap, bp, d, app, bpp, aold;

      am = 1.0;
      bm = 1.0;
      az = 1.0;
      qab = a+b;
      qap = a+1.0;
      qam = a-1.0;
      bz = 1.0-qab*x/qap;

      for(int m=1; m<=itmax; m++)
      {
         em = m;
         tem = em+em;
         d = em*(b-m)*x/((qam+tem)*(a+tem));
         ap = az+d*am;
         bp = bz+d*bm;
         d = -(a+em)*(qab+em)*x/((a+tem)*(qap+tem));
         app = ap+d*az;
         bpp = bp+d*bz;
         aold = az;
         am = ap/bpp;
         bm = bp/bpp;
         az = app/bpp;
         bz = 1.0;

         if( Math.abs(az-aold) < eps*Math.abs(az))
             return az;

       }

       return az;
   }


    // from numerical recepies
   public static double gammln(double xx)
   {
      double result = 0.0;

        double x = 0.0;
        double tmp = 0.0;
        double ser = 0.0;

      x = xx - one;
      tmp = x+fpf;
      tmp = (x+half)*Math.log(tmp)-tmp;
      ser = one;

      for(int j=0; j<6; j++)
      {
        x = x + one;
        ser = ser + cof[j]/x;
      }

      result = tmp + Math.log(stp*ser);

      return result;
   }


        public static double calcF(double df1, double df2, double f)
        {
            double x = df2/(df1*f+df2);

            double pi = Math.PI;
            double pj2 = Math.PI/2;

            if (((int)(df1))%2 == 0)
            {
                return calcF(1-x,df2,df1+df2-4,df2-2)*Math.pow(x,df2/2);
            }

            if (((int)(df2))%2 == 0)
            {
                return 1.0-calcF(x,df1,df1+df2-4,df1-2)*Math.pow(1-x,df1/2);
            }


            double tan = Math.atan(Math.sqrt(df1*f/df2));
            double a = tan/pj2;
            double sat = Math.sin(tan);
            double cot=Math.cos(tan);

            if(df2>1)
            {
                a=a+sat*cot*calcF(cot*cot,2,df2-3,-1)/pj2;
            }

            if(df1==1)
            {
                return 1-a;
            }

            double c=4*calcF(sat*sat,df2+1,df1+df2-4,df2-2)*sat*Math.pow(cot,df2)/pi;

            if(df2==1)
            {
                return 1-a+c/2;
            }

            int k=2;

            while(k<=(df2-1)/2)
            {
                c=c*k/(k-.5);
                k=k+1;
            }

            return 1-a+c;
    }


    public static double calcF(double q, double i, double j, double b)
    {
        double zz=1;
        double z=zz;
        double k=i;

        while (k<=j)
        {
            zz=zz*q*k/(k-b);
            z=z+zz;
            k=k+2;
        }

        return z;
    }



    public static double calcGammaLn(double xx)
    {
            int j;
            double x, tmp, ser;
            double cof[] = { 76.18009173, -86.50532033, 24.01409822, -1.231739156, 0.120858003e-2, -0.536382e-5 };

            x = xx - 1.;
            tmp = x + 5.5;
            tmp -= (x + 0.5) * Math.log(tmp);
            ser = 1.;

            for (j=0; j<6; ++j)
            {
                    x += 1.;
                    ser += cof[j] / x;
            }

            return(-tmp + Math.log(2.50662827465 * ser));
    }



    public static double calcGamma(double x)
    {
        return Math.exp(calcGammaLn(x));
    }


    public static int factorial(int n)
    {
        int k = 1;
        for (int i=1; i<=n; i++)
            k *= i;

        return k;
    }



    public static double calcBeta(int p, int q)
    {
        return (double)(factorial(p-1)*factorial(q-1))/factorial(p+q-1);
    }




       public static double  doAnova(double[][] xs)
            throws Exception
        {
            //grand mean

            double[] means = new double[xs.length];
            int no = 0;
            double mean = 0.0;

            for (int i=0; i<xs.length; i++)
            {
                means[i] = 0.0;

                for (int j=0; j<xs[i].length; j++)
                {
                    mean += xs[i][j];
                    means[i] += xs[i][j];
                    no++;
                }

                means[i] /= xs[i].length;
            }
            mean /= no;

            //sst
            double sst = 0.0;
            for (int i=0; i<xs.length; i++)
            for (int j=0; j<xs[i].length; j++)
                sst += (xs[i][j]-mean)*(xs[i][j]-mean);

            //ssb
            double ssb = 0.0;
            for (int i=0; i<xs.length; i++)
            {
                ssb += xs[i].length*(means[i]-mean)*(means[i]-mean);
            }


            //ss within groups
            double sgg = sst-ssb;


            ssb = ssb/(xs.length-1);
            sgg = sgg/(no-xs.length);

            double f = ssb/sgg;

            double p = Stat.calcF(xs.length-1, no-xs.length, f);

            //System.out.println(sgg+":"+ssb+":"+sst+":"+f+":"+p);

            return p;
        }


}

