package gov.nih.ncats.ifx.qhts.utils;

import java.util.*;
import java.sql.*;


public class  PlateUtil
{
    public double TINY = 1.0e-10;


    public static int colNo = 48;
    public static int rowNo = 32;

    private static HashMap cacheMap = new HashMap();

    public static double OUTLIER_PCT_CTL = 0.15;

    public static int CMPD_START_COL = 4;

    public static double OUTLIER_PCT = 0.10;
	public static char	OUTLIER_DIRECTION = ' ';

    public static double CONTROL_OUTLIER_DELTA = 3.0;


    private static PreparedStatement plateConcIS;
    private static PreparedStatement plateSampleIS;

    private static PreparedStatement plateConc384IS;
    private static PreparedStatement plateSample384IS;

    private static PreparedStatement plateConc96IS;
    private static PreparedStatement plateSample96IS;

    public static void main(String[] args)
    {

        try
        {
	    OUTLIER_PCT = 0.2;
	    OUTLIER_DIRECTION = '-';

            ArrayList list = new ArrayList();
            for (int i=0; i<10; i++)
                list.add(String.valueOf(i));

            list = removeOutlier(list);

            for (int i=0; i<list.size(); i++)
                System.out.println(list.get(i));

        }
        catch (Exception ex)
        {
            ex.printStackTrace();
        }

    }



    

    public static void setRowCol(int rn, int cn)
    {
        rowNo = rn;
        colNo = cn;
	
        //if (rowNo == 16)
        //    setCmpdStartCol(2);
    }




    public static void setCmpdStartCol(int col)
    {
        CMPD_START_COL = col;
    }



    public static void setOutlierPct(double td)
    {
        OUTLIER_PCT = td;
    }



    public static double getOutlierPct()
    {
        return OUTLIER_PCT;
    }


     public static ArrayList removeOutlier1(ArrayList list)
     {
        double[] values = new double[list.size()];

        for (int i=0; i<values.length; i++)
            values[i] = Double.parseDouble((String)list.get(i));

        Arrays.sort(values);

        ArrayList list2 = new ArrayList();
        for (int i=0; i<values.length; i++)
            list2.add(String.valueOf(values[i]));

        int no = list2.size();

		//top
        int nt = (int) (no * 0.375);
		int nb = (int) (no*0.375);

		for (int i=0; i<nt; i++)
				list2.remove(list2.size()-1);

		for (int i=0; i<nb; i++)
			list2.remove(0);


        return list2;
    }



     public static ArrayList removeOutlier(ArrayList list)
     {
        double[] values = new double[list.size()];

        for (int i=0; i<values.length; i++)
            values[i] = Double.parseDouble((String)list.get(i));

        Arrays.sort(values);

        ArrayList list2 = new ArrayList();
        for (int i=0; i<values.length; i++)
            list2.add(String.valueOf(values[i]));

        int no = list2.size();

        int n = (int) (no * OUTLIER_PCT);

		if (OUTLIER_DIRECTION == ' ')
		{
			for (int i=0; i<n; i++)
				list2.remove(0);

			for (int i=0; i<n; i++)
				list2.remove(list2.size()-1);
		}
		else if (OUTLIER_DIRECTION == '+')
		{
			for (int i=0; i<n; i++)
				list2.remove(list2.size()-1);
		}
		else if (OUTLIER_DIRECTION == '-')
		{
			for (int i=0; i<n; i++)
				list2.remove(0);
		}

        return list2;
    }


    public static double calcMean(double[][] data, int startCol, boolean[][] flags)
    {
        ArrayList list = new ArrayList();

        for (int i = 0; i < rowNo; i++)
        for (int j = startCol; j < colNo; j++) {
            if (flags != null && flags[i][j] == false)
                continue;

            if (Double.isNaN(data[i][j]))
                continue;

            list.add(String.valueOf(data[i][j]));
        }

        list = removeOutlier(list);

        int no = list.size();

        double[] values = new double[no];

        for (int i = 0; i < no; i++)
            values[i] = Double.parseDouble((String) list.get(i));


        double mean = 0.0;
        for (int i=0; i <no; i++)
            mean += values[i];

        return mean/no;
    }



    public static double calcMean(double[][] data, int col)
    {

        double mean = 0.0;
        for (int i=0; i <data.length; i++)
            mean += data[i][col];

        return mean/data.length;
    }



    public static double calcMean(double[][] data, String[][] labels, String label, boolean[][] flags)
    {

        ArrayList list = new ArrayList();

        for (int i=0;i<rowNo;i++)

        //for (int j=0;j<CMPD_START_COL;j++)
        for (int j=0;j<labels[i].length;j++)

        {
            if (flags != null && flags[i][j] == false)
                continue;

            if (Double.isNaN(data[i][j]))
                continue;

            if (labels[i][j].equals(label))
            {
                list.add(String.valueOf(data[i][j]));
            }

        }

        list = removeOutlier(list);


        double[] values = new double[list.size()];

        for (int i=0; i<values.length; i++)
            values[i] = Double.parseDouble((String)list.get(i));


        double mean = 0.0;
        for (int i=0; i<values.length; i++)
            mean += values[i];


        return mean/values.length;

    }





    public static double calcMedian(double[][] data, int startCol)
    {
        return calcMedian(data, startCol, null);
    }



    public static double calcMedian(double[][] data, int startCol, boolean[][] flags)
    {

        ArrayList list = new ArrayList();


        for (int i=0;i<rowNo;i++)
        for (int j=startCol;j<colNo;j++)
        {
            if (flags != null && flags[i][j] == false)
                continue;

            if (Double.isNaN(data[i][j]))
                continue;

            list.add(String.valueOf(data[i][j]));
        }

	if (list.size() == 0)
	{
	    return Double.NaN;
	}

        double[] values = new double[list.size()];

        for (int i=0; i<values.length; i++)
            values[i] = Double.parseDouble((String)list.get(i));


        Arrays.sort(values);


        return values[values.length/2];
    }





    public static double calcMedianCol(double[][] data, int c)
    {

        ArrayList list = new ArrayList();


        for (int i=0;i<rowNo;i++)
        {
            if (Double.isNaN(data[i][c]))
                continue;

            list.add(String.valueOf(data[i][c]));
        }


        double[] values = new double[list.size()];

        for (int i=0; i<values.length; i++)
            values[i] = Double.parseDouble((String)list.get(i));


        Arrays.sort(values);


        return values[values.length/2];
    }


    
    public static double calcMedianCol(double[][] data, int c, int r1, int r2)
    {
        ArrayList list = new ArrayList();

        for (int i=r1;i<r2;i++)
        {
            if (Double.isNaN(data[i][c]))
                continue;

            list.add(String.valueOf(data[i][c]));
        }


        double[] values = new double[list.size()];

        for (int i=0; i<values.length; i++)
            values[i] = Double.parseDouble((String)list.get(i));


        Arrays.sort(values);


        return values[values.length/2];
    }


    
    
    

    public static double calcMedianCol(double[][] data, int rowNo, int c)
    {

        ArrayList list = new ArrayList();


        for (int i=0;i<rowNo;i++)
        {
            if (Double.isNaN(data[i][c]))
                continue;

            list.add(String.valueOf(data[i][c]));
        }


        double[] values = new double[list.size()];

        for (int i=0; i<values.length; i++)
            values[i] = Double.parseDouble((String)list.get(i));


        Arrays.sort(values);


        return values[values.length/2];
    }


    
    

    public static double calcStdCol(double[][] data, int rowNo, int c)
    {
        double median = calcMedianCol(data, rowNo, c);
	
        double std = 0.0;

        for (int i=0; i<rowNo; i++)
            std += (data[i][c]-median)*(data[i][c]-median);

        return Math.sqrt(std/(data.length-1));
    }


    

    public static double calcStdCol(double[][] data, int c, double median)
    {
        double std = 0.0;

        for (int i=0; i<data.length; i++)
            std += (data[i][c]-median)*(data[i][c]-median);

        return Math.sqrt(std/(data.length-1));
    }




    public static double calcMedian(double[] values, boolean[] flags)
    {
        ArrayList list = new ArrayList();


        for (int i=0;i<values.length;i++)
        {
            if (flags == null || flags[i])
                list.add(String.valueOf(values[i]));

        }

        return calcMedian(list);
    }





    public static double calcMedian(double[][] data, String[][] labels, String label, boolean[][] flags)
    {

        ArrayList list = new ArrayList();


        for (int i=0;i<rowNo;i++)
        for (int j=0;j<labels[i].length;j++)
        {
            if (flags != null && flags[i][j] == false)
                continue;

            if (Double.isNaN(data[i][j]))
                continue;

            if (labels[i][j].equals(label))
            {
                list.add(String.valueOf(data[i][j]));
            }

        }

        if (list.size() == 0)
            return Double.NaN;

		list = removeOutlier(list);

        double[] values = new double[list.size()];

        for (int i=0; i<values.length; i++)
            values[i] = Double.parseDouble((String)list.get(i));

		Arrays.sort(values);

		if (values.length % 2 == 0)
        {
            int mid = values.length/2;

            return (values[mid]+values[mid-1])/2;
        }

        return values[values.length/2];
    }



    public static double calcMedian(ArrayList list)
    {
        double[] values = new double[list.size()];

        for (int k=0; k<list.size(); k++)
        {
            values[k] = Double.parseDouble((String)list.get(k));
        }

        Arrays.sort(values);

        return values[values.length/2];
    }





    public static double calcMedian(double[] values)
    {
        double[] values2 = new double[values.length];

        for (int i=0;i<values.length;i++)
        {
            values2[i] = values[i];
        }

        Arrays.sort(values2);

        int k = values2.length/2;

        return values2[k];
    }






    public static double calcMedianToken(double[][] data, String[][] labels, String tok, boolean[][] flags)
            throws Exception
    {
        double dm = 0.0;


        if (HtsUtil.isNumber(tok))
        {
            dm = calcMedian(data, Integer.parseInt(tok), flags);
        }
        else if (HtsUtil.isDecimal(tok))
        {
            dm = Double.parseDouble(tok);
        }
        else if (tok.length() == 1)
        {
            dm = calcMedian(data, labels, tok, flags);
        }
        else
        {
            throw new Exception("Wrong format: "+tok);
        }

        return dm;
    }





    public static double calcMedian(double[][] data, String[][] labels, String label, double d0, boolean[][] flags)
    {
        ArrayList list = new ArrayList();

        for (int i=0;i<rowNo;i++)
        for (int j=0;j<labels[i].length;j++)
        {
            if (flags != null && flags[i][j] == false)
                continue;

            if (Double.isNaN(data[i][j]))
                continue;

            if (labels[i][j].equals(label))
            {
                list.add(String.valueOf(Math.abs(data[i][j]-d0)));
            }
        }


        double[] values = new double[list.size()];

        for (int i=0; i<values.length; i++)
            values[i] = Double.parseDouble((String)list.get(i));

        Arrays.sort(values);


        return values[values.length/2];

    }





    public static double calcMedian(double[][] data, int startCol, double d0, boolean[][] flags)
    {

        ArrayList list = new ArrayList();


        for (int i=0;i<rowNo;i++)
        for (int j=startCol;j<colNo;j++)

        {
            if (flags != null && flags[i][j] == false)
                continue;

            if (Double.isNaN(data[i][j]))
                continue;

            list.add(String.valueOf(Math.abs(data[i][j]-d0)));

        }



        double[] values = new double[list.size()];

        for (int i=0; i<values.length; i++)

            values[i] = Double.parseDouble((String)list.get(i));



        Arrays.sort(values);



        return values[values.length/2];

    }







    //for correction and kinetics; no removeOutlier
    public static double calcStd(double[] data, double median)
    {
        double std = 0.0;

        for (int i=0; i<data.length; i++)
            std += (data[i]-median)*(data[i]-median);

        return Math.sqrt(std/(data.length-1));
    }



    //for plate rendering; do not remove outlier
    public static double calcStd(int startCol, double[][] data, double median)
            throws Exception
    {
        double std = 0.0;

        for (int i=0;i<rowNo;i++)
        for (int j=startCol;j<colNo;j++)
        {
            std += (data[i][j]-median)*(data[i][j]-median);
        }

        int no = rowNo*(colNo-startCol)-1;

        return Math.sqrt(std/no);
    }






    //for correction; do not apply removeOutlier
    public static double calcStd(ArrayList list, double median, int delta)
        throws Exception
    {
        double[] values = new double[list.size()];

        for (int i=0; i<values.length; i++)
            values[i] = Double.parseDouble((String)list.get(i));

        Arrays.sort(values);

        double std = 0.0;

        for (int i=delta; i<values.length-delta; i++)
        {
            std += (values[i]-median)*(values[i]-median);
        }


        return Math.sqrt(std/(values.length-2*delta-1));
    }




    public static double calcStd(double[][] data, String[][] labels, String label, double d0, double pct)
    {
        return calcStd(data, labels, label, d0, null);

    }



    public static double calcStd(double[][] data, String[][] labels, String label, double d0, boolean[][] flags)
    {

        ArrayList list = new ArrayList();



        for (int i=0;i<rowNo;i++)
        //for (int j=0;j<CMPD_START_COL;j++)
        for (int j=0;j<labels[i].length;j++)
        {
            if (flags != null && flags[i][j] == false)
                continue;

            if (Double.isNaN(data[i][j]))
                continue;

            if (labels[i][j].equals(label))
            {
                double td = data[i][j]-d0;

                list.add(String.valueOf(td*td));
            }

        }

        list = removeOutlier(list);


        double[] values = new double[list.size()];

        for (int i=0; i<values.length; i++)
            values[i] = Double.parseDouble((String)list.get(i));


        double std = 0.0;
        for (int i=0; i<values.length; i++)
            std += values[i];


        return Math.sqrt(std/(values.length-1));
    }




    public static double calcRobustStd(double[][] data, String[][] labels, String label, double d0, boolean[][] flags)
    {
        ArrayList list = new ArrayList();

        for (int i=0;i<rowNo;i++)
        for (int j=0;j<labels[i].length;j++)
        {
            if (flags != null && flags[i][j] == false)
                continue;

            if (Double.isNaN(data[i][j]))
                continue;

            if (labels[i][j].equals(label))
            {
                double td = Math.abs(data[i][j]-d0);

                list.add(String.valueOf(td));
            }
        }

        //list = removeOutlier(list);

        double[] values = new double[list.size()];

        for (int i=0; i<values.length; i++)
            values[i] = Double.parseDouble((String)list.get(i));

        Arrays.sort(values);

        return values[values.length/2];
    }




    public static double calcStd(double[][] data, int startCol, double d0)
    {
        return calcStd(data, startCol, d0, null);
    }




    public static double calcStd(double[][] data, int startCol, double d0, boolean[][] flags)
    {

        ArrayList list = new ArrayList();


        for (int i=0;i<rowNo;i++)
        for (int j=startCol;j<colNo;j++)
        {
            if (flags != null && flags[i][j] == false)
                continue;

            if (Double.isNaN(data[i][j]))
                continue;

            double td = data[i][j]-d0;

            list.add(String.valueOf(td*td));

        }


        list = removeOutlier(list);


        double[] values = new double[list.size()];

        for (int i=0; i<values.length; i++)
            values[i] = Double.parseDouble((String)list.get(i));


        double std = 0.0;
        for (int i=0; i<values.length; i++)
            std += values[i];


        return Math.sqrt(std/(values.length-1));
    }






    public static double[] calcMinMax(double[][] data, String[][] labels, String label, boolean[][] flags, double pct)
    {
        ArrayList list = new ArrayList();

        for (int i=0;i<rowNo;i++)
        for (int j=0;j<labels[i].length;j++)
        {
            if (flags != null && flags[i][j] == false)
                continue;

            if (Double.isNaN(data[i][j]))
                continue;

            if (labels[i][j].equals(label))
            {
                list.add(String.valueOf(data[i][j]));
            }

        }

        if (list.size() == 0)
            return null;

        double[] values = new double[list.size()];

        for (int i=0; i<values.length; i++)
            values[i] = Double.parseDouble((String)list.get(i));

        Arrays.sort(values);

	int no = (int)(values.length*pct);
	//min
	list = new ArrayList();
	for (int i=0; i<no; i++)
	{
		list.add(String.valueOf(values[i]));
	}
	double min = calcMedian(list);

	//max
	list = new ArrayList();
	int len = values.length;
	for (int i=0; i<no; i++)
	{
		list.add(String.valueOf(values[len-i-1]));
	}
	double max = calcMedian(list);

        double[] retValues = {min, max};

        return retValues;
    }




    public static double[] calcMinMax(double[][] data, String[][] labels, String label, boolean[][] flags)
    {
        ArrayList list = new ArrayList();

        for (int i=0;i<rowNo;i++)
        for (int j=0;j<labels[i].length;j++)
        {
            if (flags != null && flags[i][j] == false)
                continue;

            if (Double.isNaN(data[i][j]))
                continue;

            if (labels[i][j].equals(label))
            {
                list.add(String.valueOf(data[i][j]));
            }

        }

        if (list.size() == 0)
            return null;

        double[] values = new double[list.size()];

        for (int i=0; i<values.length; i++)
            values[i] = Double.parseDouble((String)list.get(i));

        Arrays.sort(values);


        double[] retValues = {values[0], values[values.length-1]};

        return retValues;
    }


}


