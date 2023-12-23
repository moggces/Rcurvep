package gov.nih.ncats.ifx.qhts.utils;

import java.io.*;
import java.util.*;
import java.sql.*;

import java.sql.Timestamp;
import java.text.*;

import org.apache.commons.math3.stat.inference.TestUtils;

public class  HtsUtil
{
	public double TINY = 1.0e-10;

    public static double LN10 = Math.log(10.0);

    public static int colNo = 48;
    public static int rowNo = 32;

    //public static SQLClient sqlClient;

    private static HashMap plateDataISMap = new HashMap();

    //private static PreparedStatement plateDataIS;
    private static PreparedStatement plateLayerIS;

    public static String[] rowLabels = {"A","B","C","D","E","F","G","H","I","J","K","L","M","N","O","P","Q","R","S","T","U","V","W","X","Y","Z","AA","AB","AC","AD","AE","AF"};

    private static HashMap layerNameMap = new HashMap();

    public static HashMap rowIndexMap;

	static
	{

		try
		{
			//String[] args2 = {"ncgcprobe.nhgri.nih.gov","1521","probedb","hts","ncgc"};
			//sqlClient = new SQLClient(args2);
			//sqlClient.connect();


                        rowIndexMap = new HashMap();
                        for (int i=0; i<rowLabels.length; i++)
                            rowIndexMap.put(rowLabels[i], String.valueOf(i+1));

		}
		catch (Exception ex)
		{
			ex.printStackTrace();
		}
	}



        public static String rc2well(int row, int col)
        {
            return rowLabels[row-1]+String.valueOf(col);
        }



        public static String rc2well(String row, String col)
        {
            return rc2well(Integer.parseInt(row), Integer.parseInt(col));

        }


        public static void setRowCol(int rn, int cn)
        {
            rowNo = rn;
            colNo = cn;
        }


        public static double[] list2array(ArrayList list)
        {
            double[] values = new double[list.size()];

            for (int i=0; i<list.size(); i++)
            {
                values[i] = Double.parseDouble((String)list.get(i));
            }

            return values;
        }


        public static double tTest(ArrayList list1, ArrayList list2)
                throws Exception
        {

            double[] values1 = list2array(list1);
            double[] values2 = list2array(list2);

            double p = TestUtils.pairedTTest(values1, values2);

            return p;
        }




        public static double[][] flipPlate(double[][] data)
        {
            int rn = data.length;
            int cn = data[0].length;

            double[][] data2 = new double[rn][cn];

            for (int r=0; r<rn; r++)
            for (int c=0; c<cn; c++)
            {
                data2[r][c] = data[rn-r-1][cn-c-1];
            }

            return data2;
        }



        public static String[][] flipPlate(String[][] data)
        {
            int rn = data.length;
            int cn = data[0].length;

            String[][] data2 = new String[rn][cn];

            for (int r=0; r<rn; r++)
            for (int c=0; c<cn; c++)
            {
                data2[r][c] = data[rn-r-1][cn-c-1];
            }

            return data2;
        }





	public static long toSecond(String ts)
	{
		return Timestamp.valueOf(ts).getTime()/1000;
	}


	public static double toHour(String ts)
	{
		return 0.001*Timestamp.valueOf(ts).getTime()/3600;
	}


	public static double calcMedian(int startCol, double[][] data)
		throws Exception
	{
		int no = rowNo*(colNo-startCol);

		double[] values = new double[no];

		int k = 0;
		for (int i=0;i<rowNo;i++)
		for (int j=startCol;j<colNo; j++)
		{
			values[k] = data[i][j];
			k++;
		}

		Arrays.sort(values);

		return values[values.length/2];
	}


	public static String fitname2id(String fit)
	{
		fit = fit.toLowerCase();

		if (fit.indexOf("smart") >= 0)
			return "0";
		else if (fit.indexOf("consant") >= 0)
			return "1";
		else if (fit.indexOf("hill") >= 0)
			return "2";
		else if (fit.indexOf("robust") >= 0)
			return "3";
		else
			return "-1";
	}


	public static String[] map2keys(HashMap map)
	{
		Object[] objs = map.keySet().toArray();
		String[] keys = new String[objs.length];
		for (int i=0; i<keys.length;i++)
			keys[i] = (String)objs[i];
		return keys;

	}


	public static String flag2id(String flag)
	{
		if (flag != null || flag.toLowerCase().equals("false"))
			return "0";
		else if (flag != null || flag.toLowerCase().equals("true"))
			return "1";
		else return "0";
	}


	public static String[] parseTokens(String line, String del)
		throws Exception
	{
		StringTokenizer tokens = new StringTokenizer(line,del);

		String[] values = new String[tokens.countTokens()];
		int i=0;
		while (tokens.hasMoreTokens())
		{
			values[i] = tokens.nextToken();
			i++;
		}

		return values;
	}



	public static void sortConcLabels(String[] labels)
	{

		Comparator comp = new Comparator()
		{
			public int compare(Object o1, Object o2)
			{
				double d1 = parseConcentration((String)o1);
				double d2 = parseConcentration((String)o2);

				if (d1 > d2)
					return 1;
				else
					return -1;
			}
		};


		Arrays.sort(labels, comp);
	}


	public static double[][] parsePlateData(String text)
		throws Exception
	{
		return parsePlateData(text, 0);
	}



	public static double[][] parsePlateData(String text, int skipLines)
		throws Exception
	{
		StringReader sr = new StringReader(text);
		BufferedReader  br = new BufferedReader(sr);

		for (int i=0; i<skipLines; i++)
			br.readLine();

		double[][] values = new double[rowNo][colNo];
		int r = 0;

		while (true)
		{
			String line = br.readLine();
			if (line == null)
				break;

			if (line.trim().length() == 0)
				continue;

			StringTokenizer tokens = new StringTokenizer(line, " \t");
			int c = 0;
			while (tokens.hasMoreTokens())
			{
				values[r][c] = Double.parseDouble(tokens.nextToken());
				c++;
			}

			r++;
		}

		br.close();
		sr.close();

		return values;
	}



	public static char[][] readPlateDef(String defFile)
		throws Exception
	{
		FileReader fr = new FileReader(new File(defFile));
		BufferedReader  br = new BufferedReader(fr);

		char[][] values = new char[rowNo][colNo];

		int r = 0;

		boolean flag = false;

		while (true)
		{
			String line = br.readLine();
			if (line == null)
				break;

			if (line.trim().length() == 0)
				continue;

			if (line.indexOf("WellTypeDisposition")> 0)
			{
				r = 0;
				flag = true;
				continue;
			}

			if (flag == false)
				continue;

			values[r] = line.trim().toCharArray();


			r++;

			if (line.indexOf("\">") > 0)
				break;
		}

		br.close();
		fr.close();

		return values;
	}




	public static ArrayList readPlateFile(String dataFile)
		throws Exception
	{
		FileReader fr = new FileReader(new File(dataFile));
		BufferedReader  br = new BufferedReader(fr);

		ArrayList list = new ArrayList();

		StringBuffer buffer = new StringBuffer();

		boolean flag = false;

		while (true)
		{
			String line = br.readLine();
			if (line == null)
				break;

			if (line.trim().length() == 0)
				continue;

			if (flag && line.indexOf(";Plate")>= 0)
			{
				if (buffer != null && buffer.length() > 0)
				{
					list.add(buffer);
					buffer = new StringBuffer();
				}

				continue;
			}

			if (line.indexOf(";Plate")>= 0)
			{
				flag = true;
				continue;
			}

			buffer.append(line+"\n");
		}

		br.close();
		fr.close();

		if (buffer != null && buffer.length() > 0)
		{
			list.add(buffer);
		}

		return list;
	}

	public static HashMap readSectionData(StringBuffer dataBuffer, String[] labels)
		throws Exception
	{
		StringReader sr = new StringReader(dataBuffer.toString());
		BufferedReader  br = new BufferedReader(sr);

		HashMap map = new HashMap();

		StringBuffer buffer = new StringBuffer();

		String prevLabel = null;
		String currLabel = null;

		while (true)
		{
			String line = br.readLine();
			if (line == null)
				break;

			if (line.trim().length() == 0)
				continue;

			if (line.indexOf("Identification barcode")>= 0)
			{
				int k = line.lastIndexOf(" ");
				String barcode = line.substring(k).trim();
				map.put("barcode", barcode);
			}

			boolean flag = false;

			for (int i=0; i<labels.length; i++)
			{
				if (line.indexOf(labels[i]) >= 0)
				{
					flag = true;
					currLabel = labels[i];

					if (buffer != null && buffer.length() > 0)
					{
						map.put(prevLabel, buffer);
						prevLabel = currLabel;
						buffer = new StringBuffer();
						buffer.append(line+"\n");
					}
					else
					{
						prevLabel = currLabel;
						buffer.append(line+"\n");
					}

					break;
				}
			}

			if (flag)
				continue;

			buffer.append(line+"\n");
		}

		br.close();
		sr.close();

		map.put(prevLabel, buffer);

		return map;
	}


	public static HashMap[] readCondoseoMap(String mapFile, int skipLines, int[] indexes)
		throws Exception
	{
		FileReader fr = new FileReader(new File(mapFile));
		BufferedReader  br = new BufferedReader(fr);

		HashMap[] maps = new HashMap[indexes.length/2];
		for (int i=0; i<maps.length; i++)
			maps[i] = new HashMap();


		for (int j=0; j<skipLines; j++)
		{
			String line = br.readLine();

			StringTokenizer tokens = new StringTokenizer(line,"\t");

			String[] values = new String[tokens.countTokens()];
			int i=0;
			while (tokens.hasMoreTokens())
			{
				values[i] = tokens.nextToken();
				i++;
			}

			for (i=0; i<indexes.length/2; i++)
			{
				int k = indexes[2*i];
				int l = indexes[2*i+1];

				//System.out.println(line+":"+(2*i+1)+":"+l+":"+values[l]);

				maps[i].put("header",values[l]);
			}
		}


		while (true)
		{
			String line = br.readLine();
			if (line == null)
				break;

			if (line.trim().length() == 0)
				continue;

			if (line.startsWith("#"))
				continue;

			StringTokenizer tokens = new StringTokenizer(line,"\t");

			String[] values = new String[tokens.countTokens()];
			int i=0;
			while (tokens.hasMoreTokens())
			{
				values[i] = tokens.nextToken();
				i++;
			}

			for (i=0; i<indexes.length/2; i++)
			{
				int k = indexes[2*i];
				int l = indexes[2*i+1];

				maps[i].put(values[k],values[l]);
			}
		}

		br.close();
		fr.close();

		return maps;
	}


	public static String[] readHeaders(String dataFile)
		throws Exception
	{
		FileReader fr = new FileReader(new File(dataFile));
		BufferedReader  br = new BufferedReader(fr);

		String line = br.readLine();

		StringTokenizer tokens = new StringTokenizer(line,"\t");

		String[] values = new String[tokens.countTokens()];
		int i=0;
		while (tokens.hasMoreTokens())
		{
			values[i] = tokens.nextToken();
			i++;
		}

		return values;
	}

	public static int conc2index(String conc1, String[] headers)
	{
		double d1 = Double.parseDouble(conc1);

		double minR = 10000.0;
		int minI = -1;

		for (int j=5; j<headers.length-1; j++)
		{
			double d2 = label2conc(headers[j]);
			double r = d1/d2;

			if (Math.abs(r-1.0) < minR)
			{
				minR = Math.abs(r-1.0);
				minI = j-5;
			}
		}

		return minI;
	}


	public static double label2conc(String conc2)
	{
		double d2 = 0.0;

		conc2 = conc2.trim();

		int l = conc2.length();

		String unit = conc2.substring(l-2,l).toLowerCase();
		String amount = conc2.substring(0, l-2);

		//System.out.println(unit+":"+amount);

		if (unit.indexOf("nm")>=0 )
			d2 = Double.parseDouble(amount);
		else if (unit.indexOf("um") >= 0)
			d2 = Double.parseDouble(amount)*1000.0;
		else if (unit.indexOf("mm") >= 0)
			d2 = Double.parseDouble(amount)*1000000.0;
		else
		{
			System.out.println("Unknown concentration unit!");
			System.exit(1);
		}

		return d2;
	}



	public static double parseConcentration(String text)
	{
		text = text.trim().toLowerCase();

		String value = text.substring(0,text.length()-2);

		double d = Double.parseDouble(value);

		if (text.endsWith("nm"))
			return d*0.001;
		else if (text.endsWith("mm"))
			return 1000.0*d;
		else
			return d;


	}








        public static double parseConcentrationMolar(String text)

	{

		text = text.trim().toLowerCase();



		String value = text.substring(0,text.length()-2);



		double d = Double.parseDouble(value);



		if (text.endsWith("nm"))

			return d*1e-9;

		else if (text.endsWith("mm"))

			return d*1e-3;

		else if (text.endsWith("um"))

			return d*1e-6;

                else if (text.endsWith("m"))

                    return Double.parseDouble(text.substring(0,text.length()-1));
                else
                    return Double.parseDouble(value);

	}



        public static String parseNumber(String text)
        {
            StringBuffer buffer = new StringBuffer();
            for (int i=0; i<text.length(); i++)
            {
                char c = text.charAt(i);
                if (Character.isDigit(c) || c == '.')
                    buffer.append(c);
            }

            return buffer.toString();
        }


        public int getDigitStart(String text)
        {
            for (int i=0; i<text.length(); i++)
            {
                char c = text.charAt(i);
                if (Character.isDigit(c))
                    return i;
            }

            return -1;

        }







	public static int[] parseIndexes(String text)
	{
		StringTokenizer tokens = new StringTokenizer(text,",;");
		int[] indexes = new int[tokens.countTokens()];

		int k = 0;
		while (tokens.hasMoreTokens())
		{
			indexes[2*k] = Integer.parseInt(tokens.nextToken());
			indexes[2*k+1] = Integer.parseInt(tokens.nextToken());
			k++;
		}

		return indexes;
	}

        public static String formatConcentration(String conc, String pat)
        {
            double tc = Double.parseDouble(conc);

            return formatDouble(tc, pat);
        }


        public static String formatConcentration(String conc)
        {
            double tc = Double.parseDouble(conc);

            return formatDouble(tc, "#0.0000");
        }


        public static String formatConcentration(double tc)
        {
            return formatDouble(tc, "#0.0000");
        }


        public  static  String  formatDouble(String conc, String digits)
	{
		double td = Double.parseDouble(conc);

                int no = Integer.parseInt(digits);

                //"0.###E0"

                String pattern = "0.";
                for (int i=0; i<no; i++)
                    pattern += "#";
                pattern += "E0";

		DecimalFormat formatter = new DecimalFormat(pattern);


                return formatter.format(td);
        }


	public  static  String  formatDouble(double d, int right)
	{

		DecimalFormat format = new DecimalFormat("#,##0.00000000");
		FieldPosition f = new FieldPosition(0);
		StringBuffer s = new StringBuffer();

		String dString = format.format(d, s, f).toString();

		int pos = dString.indexOf('.');

		if (pos < 0)
			return dString;

		String	leftString = dString.substring(0, pos);

		int delta = dString.length() - pos;

		String rightString = "";

		if (delta > right)
			rightString = dString.substring(pos+1, pos+1+right);
		else
			rightString = dString.substring(pos+1);

		rightString = "."+rightString;

		return leftString + rightString;
	}



        public static double reformatDouble(double d, String pattern)
        {
            return Double.parseDouble(formatDouble(d, pattern));
        }




	public  static  String  formatDouble(double d, String pattern)
	{
            DecimalFormat myFormatter = new DecimalFormat(pattern);
            String output = myFormatter.format(d);

            return output;
	}





	public static String combineNames(String name1,String name2)
	{
		String name = "";
		if (name1 != null)
			name = name+" "+name1;

		if (name2 != null)
			name = name+" "+name2;

		return name;
	}

	public  static  String  expLogAC50(String value)
	{
		if (value == null || value.trim().length() == 0)
			return "";

		double d = Double.parseDouble(value);
		d = Math.exp(d*Math.log(10.0));

		DecimalFormat myFormatter = new DecimalFormat("0.#####E0");
		String output = myFormatter.format(d);

		return output;
	}


	public static double[] lineFit(double[] xs, double ys[])
	{
	    	int no = xs.length;

	    	boolean[] flags = new boolean[no];

	    	for (int i=0; i<no; i++)
	    		flags[i] = true;

	    	return lineFit(flags,xs,ys);
	}



    public static double[] lineFit(boolean flags[], double[] xs, double ys[])
    {
    	int no = xs.length;

	    double xm = 0;
	    double ym = 0;

	    for (int i=0; i<no; i++)
	    if (flags[i])
	    {
	            xm += xs[i];
	            ym += ys[i];
	    }

	    xm /= no;
	    ym /= no;

	    double lxx = 0;
	    double lyy = 0;
	    double lxy = 0;
	    for (int i=0; i<no; i++)
	    if (flags[i])
	    {
	            lxx += (xs[i]-xm)*(xs[i]-xm);
	            lyy += (ys[i]-ym)*(ys[i]-ym);
	            lxy += (xs[i]-xm)*(ys[i]-ym);
	    }

	    double b = lxy/lxx;
	    double a = ym-b*xm;
	    double r = lxy/Math.sqrt(lxx*lyy);

	    double[] values = {b,a,r};

	    return values;
	}


	public static String removeControls(String text)
	{
		StringBuffer buffer = new StringBuffer(text.trim());
		for (int i=0; i<buffer.length(); i++)
		{
			char ch = buffer.charAt(i);

			if (ch == '-' || ch == '_' || Character.isLetterOrDigit(ch))
				continue;

			buffer.deleteCharAt(i);
				i--;
		}

		return buffer.toString();
	}


	public static ArrayList getFileList(String dir, String filter, String filter2)
		throws Exception
	{
		ArrayList fileList = new ArrayList();

		final String pattern = filter;
		final String pattern2 = filter2;

		FileFilter fileFilter = new FileFilter()
		{
			public boolean accept(File file)
			{
				String name = file.getName();

				boolean flag1 = pattern == null || name.matches(pattern);
				boolean flag2 = pattern2 == null || name.matches(pattern2) == false;

				if (flag1 && flag2)
					return true;
				else
					return false;
			}
		};


		getFileList(new File(dir), fileFilter, fileList);

		return fileList;
	}


	public static void getFileList(File file, FileFilter filter, ArrayList fileList)
		throws Exception
	{

		if (file.isFile() && file.length() > 0)
			fileList.add(file);
		else if (file.isFile() == false)
		{
			File[] files = file.listFiles(filter);

			for (int i=0; i<files.length; i++)
				getFileList(files[i], filter, fileList);
		}
	}


        //Q4
	public static String convertBarcode(String barcode)
		throws Exception
	{
            if (barcode == null)
                return null;

            barcode = barcode.trim();

            int len = barcode.length();

            if (len < 2)
		    return barcode;
            
            if (Character.isDigit(barcode.charAt(len-1)) && barcode.charAt(len-2) == 'Q')
            {
                return convertBarcodeNoQ(barcode.substring(0,len-2))+barcode.substring(len-2);
            }
            
            return convertBarcodeNoQ(barcode);
        }
        
            
	public static String convertBarcodeNoQ(String barcode)
		throws Exception
	{
            if (barcode == null)
                return null;

            barcode = barcode.trim();

            int len = barcode.length();

            if (len < 2)
		    return barcode;
                //throw new Exception("Missing barcode: "+barcode);

            //if (len != 8)
            //	throw new Exception("Barcode length :"+len);

            int k = 1;
            while (k<len)
            {
                char c  = barcode.charAt(k);

                if (c != '0' && Character.isDigit(c))
                            break;

                k++;
            }

            if (k == len)
                    return barcode;

            String barcode2 = barcode;

            String s1 = barcode.substring(k);

            long l = Long.parseLong(s1);

            if (l%2 == 0)
            {
                    String s2 = String.valueOf(l-1);

                    int d = s1.length()-s2.length();

                    barcode2 = barcode.substring(0,k);

                    for (int i=0; i<d; i++)
                            barcode2 = barcode2+"0";

                    barcode2 = barcode2+s2;
            }

            //System.out.println(barcode+":"+barcode2);

            return barcode2;
	}


	public static String[] splitStringNumber(String value)
	{
            if (value == null)
                return null;

            value = value.trim();

            int len = value.length();

            int k = 1;
            while (k<len)
            {
                char c  = value.charAt(k);

                if (Character.isDigit(c))
		    break;

                k++;
            }

	    if (k == len)
	    {
		String[] vals = {value, null};
		
		return vals;
	    }

           
            String[] vals = {value.substring(0, k), value.substring(k)};
	    
	    return vals;
	}





	public static String array2in(String[][] values, int start, int end, String tag, int col, int inFlag)
		throws Exception
	{
            StringBuffer buffer = null;

            if (inFlag == 0)
                buffer = new StringBuffer(" in (");
            else
                buffer = new StringBuffer(" not in (");

            buffer.append(tag);
            buffer.append(values[start][col]);
            buffer.append(tag);

            for (int i=start+1; i<values.length && i<end; i++)
            {
                buffer.append(",");
                buffer.append(tag);
                buffer.append(values[i][col]);
                buffer.append(tag);
            }

            buffer.append(")");

            return buffer.toString();
	}



	
	public static String array2in(ArrayList list, int start, int end, String tag, int inFlag)
		throws Exception
	{
           StringBuffer buffer = null;

            if (inFlag == 0)
                buffer = new StringBuffer(" in (");
            else
                buffer = new StringBuffer(" not in (");

            buffer.append(tag);
            buffer.append(list.get(start));
            buffer.append(tag);

            for (int i=start+1; i<list.size() && i<end; i++)
            {
                buffer.append(",");
                buffer.append(tag);
                buffer.append(list.get(i));
                buffer.append(tag);
            }

            buffer.append(")");

            return buffer.toString();
	}


	public static String array2in(String[] values, int start, int end, String tag, int inFlag)
		throws Exception
	{
           StringBuffer buffer = null;

            if (inFlag == 0)
                buffer = new StringBuffer(" in (");
            else
                buffer = new StringBuffer(" not in (");

            buffer.append(tag);
            buffer.append(values[start]);
            buffer.append(tag);

            for (int i=start+1; i<values.length && i<end; i++)
            {
                buffer.append(",");
                buffer.append(tag);
                buffer.append(values[i]);
                buffer.append(tag);
            }

            buffer.append(")");

            return buffer.toString();
	}


        public static String array2in(String[][] values, String tag, int col, String colName)
		throws Exception
	{
            return array2in(values, tag, col, colName, 0);
        }


        public static String array2in(String[][] values, String tag, int col, String colName, int inFlag)
		throws Exception
	{
            if (values.length < 1000)
                return colName+" "+array2in(values, tag, col, inFlag);

            StringBuffer buffer = new StringBuffer(" ( ");

            int n = values.length/1000;
            if (values.length % 1000 != 0)
                n++;

            for (int i=0; i<n; i++)
            {
                int start = 1000*i;
                int end = 1000*i+1000;
                if (end > values.length)
                    end = values.length;

                String in = array2in(values, start, end, tag, col, inFlag);

                if (i == 0)
                    buffer.append(colName +" "+in+" ");
                else
                    buffer.append(" or "+colName+" "+in);
            }

            buffer.append(" ) ");

            return buffer.toString();
	}




    public static String query2in(String text, String tag, String colName, int inFlag)
        throws Exception
    {
        if (text == null || text.trim().length() == 0)
		{
			if (inFlag == 0)
				return colName+" in ";
			else
				return colName+" not in ";
		}

        StringTokenizer tokens = new StringTokenizer(text, " ,;\t\r\n");

        String[] values = new String[tokens.countTokens()];

        int k = 0;

        while (tokens.hasMoreTokens())
        {
            values[k] = tokens.nextToken();
            k++;
        }

        return array2in(values, tag, colName, inFlag);

    }




        public static String array2in(String[] values, String tag, String colName, int inFlag)
		throws Exception
	{
            if (values.length < 1000)
                return colName+" "+array2in(values, tag, inFlag);

            StringBuffer buffer = new StringBuffer(" ( ");

            int n = values.length/1000;
            if (values.length % 1000 != 0)
                n++;

            for (int i=0; i<n; i++)
            {
                int start = 1000*i;
                int end = 1000*i+1000;
                if (end > values.length)
                    end = values.length;

                String in = array2in(values, start, end, tag, inFlag);

                if (i == 0)
                    buffer.append(colName +" "+in+" ");
                else
                    buffer.append(" or "+colName+" "+in);
            }

            buffer.append(" ) ");

            return buffer.toString();
	}



	public static String array2in(String[][] values, String tag, int col)
		throws Exception
	{
            return array2in(values, tag, col, 0);
        }


	public static String array2in(String[][] values, String tag, int col, int inFlag)
		throws Exception
	{
                StringBuffer buffer = null;

                if (inFlag == 0)
                    buffer = new StringBuffer(" in (");
                else
                    buffer = new StringBuffer(" not in (");

		buffer.append(tag);
		buffer.append(values[0][col]);
		buffer.append(tag);

		for (int i=1; i<values.length; i++)
		{
			buffer.append(",");
			buffer.append(tag);
			buffer.append(values[i][col]);
			buffer.append(tag);
		}

		buffer.append(")");

		return buffer.toString();
	}




        public static String string2in(String text, String del, String tag, int inFlag)
		throws Exception
	{
                StringTokenizer tokens = new StringTokenizer(text, del);

                StringBuffer buffer = null;

                if (inFlag == 0)
                    buffer = new StringBuffer(" in (");
                else
                    buffer = new StringBuffer(" not in (");

		buffer.append(tag);
		buffer.append(tokens.nextToken());
		buffer.append(tag);

		while (tokens.hasMoreTokens())
                {
			buffer.append(",");
			buffer.append(tag);
			buffer.append(tokens.nextToken());
			buffer.append(tag);
		}

		buffer.append(")");

		return buffer.toString();
	}






        public static String array2in(Object[] values, String tag, int inFlag)
		throws Exception
	{
                StringBuffer buffer = null;

                if (inFlag == 0)
                    buffer = new StringBuffer(" in (");
                else
                    buffer = new StringBuffer(" not in (");

		buffer.append(tag);
		buffer.append((String)values[0]);
		buffer.append(tag);

		for (int i=1; i<values.length && i<999; i++)
		{
			buffer.append(",");
			buffer.append(tag);
			buffer.append((String)values[i]);
			buffer.append(tag);
		}

		buffer.append(")");

		return buffer.toString();
	}




        public static String array2in(String[] values, String tag, int inFlag)
		throws Exception
	{
                StringBuffer buffer = null;

                if (inFlag == 0)
                    buffer = new StringBuffer(" in (");
                else
                    buffer = new StringBuffer(" not in (");

		buffer.append(tag);
		buffer.append(values[0]);
		buffer.append(tag);

		for (int i=1; i<values.length && i<999; i++)
		{
			buffer.append(",");
			buffer.append(tag);
			buffer.append(values[i]);
			buffer.append(tag);
		}

		buffer.append(")");

		return buffer.toString();
	}






        public static String list2in(ArrayList list, String tag, int end, int inFlag)
		throws Exception
	{
            return list2in(list, tag, 0, end, inFlag);

        }

        public static String list2in(ArrayList list, String tag, int start, int end, int inFlag)
		throws Exception
	{
		StringBuffer buffer = null;

                if (inFlag == 0)
                    buffer = new StringBuffer(" in (");
                else
                    buffer = new StringBuffer(" not in (");

		buffer.append(tag);
		buffer.append((String)list.get(start));
		buffer.append(tag);

                if (end < 0)
                    end = list.size();

		for (int i=start+1; i<end; i++)
		{
                    buffer.append(",");
                    buffer.append(tag);
                    buffer.append((String)list.get(i));
                    buffer.append(tag);
		}

		buffer.append(")");

		return buffer.toString();
	}


        public static String list2in(ArrayList list, String tag, String colName, int inFlag)
		throws Exception
	{
            if (list.size() < 1000)
                return colName+" "+list2in(list, tag, 0, list.size(), inFlag);

            StringBuffer buffer = new StringBuffer(" ( ");

            int n = list.size()/1000;
            if (list.size() % 1000 != 0)
                n++;

            for (int i=0; i<n; i++)
            {
                int start = 1000*i;
                int end = 1000*i+1000;
                if (end > list.size())
                    end = list.size();

                String in = list2in(list, tag, start, end, inFlag);

                if (i == 0)
                    buffer.append(colName +" "+in+" ");
                else
                    buffer.append(" or "+colName+" "+in);
            }

            buffer.append(" ) ");

            return buffer.toString();
        }



        public static String list2in(ArrayList list, String tag, String colName, int inFlag, boolean mapFlag)
		throws Exception
	{
            if (mapFlag == false)
                return list2in(list, tag, colName, inFlag);

            HashMap map = new HashMap();

            for (int i=0; i<list.size(); i++)
                map.put(list.get(i), null);

            ArrayList list2 = new ArrayList();
            Object[] keys = map.keySet().toArray();

            for (int i=0; i<keys.length; i++)
                list2.add(keys[i]);

            return list2in(list2, tag, colName, inFlag);
        }


        public static boolean isNumber(String text)
		throws Exception
	{
		boolean flag = true;

		for (int i=0; i<text.length(); i++)
                {
                    char c = text.charAt(i);

                    flag = flag & (Character.isDigit(c) | c =='-' | c == '.' ) ;
                }



		return flag;
	}


        public static boolean isDecimal(String text)
		throws Exception
	{
                if (text.length()<1) return false;

                char c = text.charAt(0);
                if (text.length() == 1 && !Character.isDigit(c))
                    return false;

                else if (!(Character.isDigit(c) | c =='-' | c == '.' ))
                    return false;

		for (int i=1; i<text.length(); i++)
                {
                    c = text.charAt(i);
                    if (!(Character.isDigit(c) | c =='.'))
                        return false;
                }

		return true;
	}





    public static String[] parseRC(String well)
    {
        int k = 0;

        for (int i=0; i<well.length(); i++)
        {
            if (Character.isDigit(well.charAt(i)))
            {
                k = i;
                break;
            }
        }

        String[] rc = {well.substring(0, k), well.substring(k)};

        return rc;
    }


	public static int[] parseWell(String well)
	{
		String[] toks = parseRC(well);

		int[] rc = new int[2];

		String r = (String)rowIndexMap.get(toks[0]);

		rc[0] = Integer.parseInt(r)-1;
		rc[1] = Integer.parseInt(toks[1])-1;

		return rc;
	}

    public static boolean[][] createIndexList(int row, int col, String toks)
    {
		boolean[][] flags = new boolean[row][col];

		StringTokenizer tokens = new StringTokenizer(toks, " ,;");
		while (tokens.hasMoreTokens())
		{
			String tok = tokens.nextToken();

			int[][] rcs = createIndexList(tok);

			for (int r=0; r<rcs.length; r++)
			for (int c=0; c<rcs[r].length; c++)
			{
				flags[r][c] = true;
			}
		}

		return flags;
	}


    public static int[][] createIndexList(String tok)
    {

        String[] rcs = tok.split("-");
        String[] rc1 = parseRC(rcs[0]);
        String[] rc2 = parseRC(rcs[1]);

        String s1 = (String)rowIndexMap.get(rc1[0]);
        String s2 = (String)rowIndexMap.get(rc2[0]);

        int r1 = Integer.parseInt(s1);
        int r2 = Integer.parseInt(s2);



        int c1 = Integer.parseInt(rc1[1])-1;
        int c2 = Integer.parseInt(rc2[1])-1;

        int d = ((r2-r1)+1)*((c2-c1)+1);

        int[][] indexes = new int[d][2];

        int k = 0;

        for (int r=r1; r<=r2; r++)
        for (int c=c1; c<=c2; c++)
        {
            indexes[k][0] = r;
            indexes[k][1] = c;
            k++;
        }

        return indexes;
    }



    //ec conversion; log unit
    public static double eccalc(double ec50, double slope, double pct)
    {
        return ec50-(1/slope)*Math.log10((100-pct)/pct);
    }


    //ic  or absolute
    public static double iccalc(double y0, double yinf, double ec50, double slope, double pct)
    {

        //return ec50-(1/slope)*Math.log10((yinf-y0)/pct-1);
        return ec50-(1/slope)*Math.log10((yinf-y0)/(pct-y0)-1);
    }

}


