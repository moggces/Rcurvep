package  gov.nih.ncats.ifx.qhts.utils;

/*

	parse sw.xml and put them into oracle database
 */


import java.io.*;
import java.sql.*;
import java.util.*;
import java.net.*;
import java.text.*;

import java.util.Date;

public class Util
	extends Object
{
	public static Integer I = new Integer(1);

	public static HashMap monthMap;
        public static HashMap monthMap2;
        public static HashMap monthMap3;

	private static String[] patches;

	static
	{
			patches = new String[10];

			for (int i=1; i<10; i++)
			{
				patches[i] = new String();
				for (int j=0; j<i; j++)
					patches[i] = patches[i]+" ";
			}


            monthMap = new HashMap();
            monthMap.put("jan", new Integer(1));
            monthMap.put("feb", new Integer(2));
            monthMap.put("mar", new Integer(3));
            monthMap.put("apr", new Integer(4));
            monthMap.put("may", new Integer(5));
            monthMap.put("jun", new Integer(6));
            monthMap.put("jul", new Integer(7));
            monthMap.put("aug", new Integer(8));
            monthMap.put("sep", new Integer(9));
            monthMap.put("oct", new Integer(10));
            monthMap.put("nov", new Integer(11));
            monthMap.put("dec", new Integer(12));

            monthMap2 = new HashMap();
            monthMap2.put("jan", "1");
            monthMap2.put("feb", "2");
            monthMap2.put("mar", "3");
            monthMap2.put("apr", "4");
            monthMap2.put("may", "5");
            monthMap2.put("jun", "6");
            monthMap2.put("jul", "7");
            monthMap2.put("aug", "8");
            monthMap2.put("sep", "9");
            monthMap2.put("oct", "10");
            monthMap2.put("nov", "11");
            monthMap2.put("dec", "12");

            monthMap3 = new HashMap();
            monthMap3.put("jan", "01");
            monthMap3.put("feb", "02");
            monthMap3.put("mar", "03");
            monthMap3.put("apr", "04");
            monthMap3.put("may", "05");
            monthMap3.put("jun", "06");
            monthMap3.put("jul", "07");
            monthMap3.put("aug", "08");
            monthMap3.put("sep", "09");
            monthMap3.put("oct", "10");
            monthMap3.put("nov", "11");
            monthMap3.put("dec", "12");
	}

	public static void main(String[] args)
	{
		try
		{
			sscanf("");
			if (true) System.exit(1);
			
			String str = "2005-07-06 16:29:56.0";

			System.out.println(Util.oracleDateToSystemTime(str));

			str = "2005-07-06 16:30:56.0";

			//System.out.println(Util.oracleDateToSystemTime(str));

                        String[] toks = parseLineCSV("", ",");
                        
                        for (int i=0; i<toks.length; i++)
                        {
                            System.out.println(toks[i]);
                        }
			//String pattern = "#0.00000000";
			//DecimalFormat myFormatter = new DecimalFormat(pattern);
			//String output = myFormatter.format(Double.parseDouble(args[0]));
			//System.out.println(args[0] + " " + pattern + " " + output);

			//System.out.println(Util.getLocalIP());
		}
		catch (Exception ex)
		{
			ex.printStackTrace();
		}
	}

	private Util()
	{
	}



	public static double[] sscanf(String text)
	{
		text = "-79.842 -20.813-100.380 4.5 -5.6 7.8";
		
		StringTokenizer toks = new StringTokenizer(text, " -", true);
		
		ArrayList<String> list1 = new ArrayList<String>();
		
		while (toks.hasMoreTokens())
		{
			String tok = toks.nextToken();
			
			if (tok.equals(" "))
				continue;
			
			if (tok.equals("-") && toks.hasMoreTokens())
			{
				list1.add(tok+toks.nextToken());
			}
			
			list1.add(tok);
		}
		
		ArrayList<String> list2 = new ArrayList<String>();
		
		int i = 0;
		while (i < list1.size())
		{
			if (list1.get(i).equals("-"))
			{
				i++;
				continue;
			}
			
			list2.add(list1.get(i));
			i++;
		}
		
		double[] ds = new double[list2.size()];
		
		for (i=0; i<ds.length; i++)
		{
			ds[i] = Double.parseDouble(list2.get(i));
			
			System.out.println(ds[i]);
		}
		
		return ds;
	}

	
	public static int countChar(String line, char ch)

	{

		int no = 0;



		for (int i=0; i<line.length(); i++)

		{

			if (line.charAt(i) == ch)

			no++;

		}



		return no;

	}



        public static String[] splitLine(String line, String del)

		{

			if (line == null || line.trim().length() == 0)

				return null;



			String row = line;



			ArrayList list = new ArrayList();



			int start = 0;

			int end = row.indexOf(del, start);



			while (end >= 0 && end <= row.length()-1)

			{

				String str = row.substring(start, end).trim();

				if (str.startsWith("\""))

				str = str.substring(1);



				if (str.endsWith("\""))

				str = str.substring(0, str.length()-1);



				list.add(str);



				start = end+del.length();;

				end = row.indexOf(del, start);



				int count = 0;

				while (end >= 0 && end <= row.length()-1)

				{

					str = row.substring(start, end).trim();



					int no = countChar(str,'\"');



					if (no%2 == 0)

					break;



					end = row.indexOf(del, end+1);

				}


			}



		end = row.indexOf(del, start);

		if (end < 0)

		list.add(row.substring(start).trim());

		else

		list.add(row.substring(start, end).trim());



		String[] fields = new String[list.size()];



		for (int i=0; i<list.size(); i++)

		{

			String value = (String)list.get(i);

			if (value.toLowerCase().equals("null"))

			value = null;

			else if (value.toLowerCase().equals("nan"))

			value = null;

			else

			fields[i] = value;

		}



	return fields;

	}


        
	public static String[] parseLineCSV(String line, String del)
            throws Exception
	{                
            if (line == null || line.trim().length() == 0)
                return null;

            String row = line;

            ArrayList list = new ArrayList();

            int dqn = 0;
            
            int start = 0;
            int end = row.indexOf(del, start);

            while (end >= 0 && end <= row.length()-1)
            {
                String tok = row.substring(start, end).trim();
                
                if (tok.startsWith("\""))
                {
                    end = row.indexOf("\"",end+1);
                    
                    if (end > 0)
                    {
                        tok = row.substring(start+1, end).trim();
                    }
                }
                
                list.add(tok);

                start = end+del.length();;
                end = row.indexOf(del, start);
            }

            end = row.indexOf(del, start);
            if (end < 0)
                    list.add(row.substring(start).trim());
            else
                    list.add(row.substring(start, end).trim());

            String[] fields = new String[list.size()];

            for (int i=0; i<list.size(); i++)
            {
                    String value = (String)list.get(i);
                    if (value.toLowerCase().equals("null"))
                            value = null;
                    else if (value.toLowerCase().equals("nan"))
                            value = null;
                    else
                            fields[i] = value;
            }

            return fields;        
        }

        
	public static String[] parseLine(String line, String del)
		throws Exception
	{
            if (line == null || line.trim().length() == 0)
                return null;
            
            String row = line;

            ArrayList list = new ArrayList();

            int start = 0;
            int end = row.indexOf(del, start);

            while (end >= 0 && end <= row.length()-1)
            {
                    list.add(row.substring(start, end).trim());

                    start = end+del.length();;
                    end = row.indexOf(del, start);
            }

            end = row.indexOf(del, start);
            if (end < 0)
                    list.add(row.substring(start).trim());
            else
                    list.add(row.substring(start, end).trim());

            String[] fields = new String[list.size()];

            for (int i=0; i<list.size(); i++)
            {
                    String value = (String)list.get(i);
                    if (value.toLowerCase().equals("null"))
                            value = null;
                    else if (value.toLowerCase().equals("nan"))
                            value = null;
                    else
                            fields[i] = value;
            }

            return fields;
	}



	public static HashMap calcWordMap(String[] rows, String del, String del2)
	{
		HashMap map = new HashMap();

		StringBuffer buffer = new StringBuffer();
		for (int i=0; i<rows.length; i++)
		{
			if (rows[i] == null)
				continue;

			buffer.append(rows[i]);
			buffer.append(del);
		}

		StringTokenizer tokens = new StringTokenizer(buffer.toString(), del);
		while (tokens.hasMoreTokens())
		{
			String token = tokens.nextToken();
			int k = token.indexOf(del2);
			if (k > 0)
				token = token.substring(0,k);

			Object obj = map.get(token);

			if (obj == null)
				map.put(token, new Integer(1));
			else
			{
				map.put(token, new Integer(1+((Integer)obj).intValue()));
			}
		}

		return map;
	}


	public static HashMap row2map(String row)
		throws Exception
	{
		StringBuffer buffer = new StringBuffer(row);

		int count = 0;
		for (int j=0; j<buffer.length(); j++)
		{
			if (buffer.charAt(j) == '"')
				count++;


			if (count%2 == 0 && buffer.charAt(j) == ' ')
				buffer.setCharAt(j,';');
		}

		StringTokenizer tokens = new StringTokenizer(buffer.toString(),";");

		HashMap map = new HashMap();

		while (tokens.hasMoreTokens())
		{
			String token = tokens.nextToken();
			map.put(token, I);
		}

		return map;
	}


	public static String cleanText(String name, String pat)
		throws Exception
	{
		StringTokenizer tokens = new StringTokenizer(name,pat);

		StringBuffer buffer = new StringBuffer();
		while (tokens.hasMoreTokens())
		{
			buffer.append(tokens.nextToken());
		}

		return buffer.toString();

	}

        


	public static String cleanControl(String name)
		throws Exception
	{
            
            if (true) 
                return name.replaceAll("[^\\x00-\\x7F]", "");
            
		char[] cs = name.toCharArray();
		
		StringBuffer buffer = new StringBuffer();
		
		for (int i=0; i<cs.length; i++)
		{
		   
		    
		    buffer.append(cs[i]);
		}

		return buffer.toString();
	}

	
	public static boolean isNumber(String text)
		throws Exception
	{
		text = cleanText(text," .,:;-_");

		boolean flag = true;

		for (int i=0; i<text.length(); i++)
			flag = flag & Character.isDigit(text.charAt(i));

		return flag;
	}





	public static ArrayList split(String line, String pat)
		throws Exception
	{
		int len = pat.length();

		String row = line;

		ArrayList list = new ArrayList();

		int start = 0;
		int end = row.indexOf(pat, start);

		while (end >= 0 && end <= row.length()-1)
		{
			if (end > start)
				list.add(row.substring(start+len, end).trim());

			start = end+len;
			end = row.indexOf(pat, start);
		}

		if (end < 0)
			list.add(row.substring(start).trim());
		else
			list.add(row.substring(start, end).trim());

		return list;
	}


	public static String empty2null(String text)
	{
	    if (text != null && text.trim().length() == 0)
		return null;
	    
	    return text;
	}
	
	
	public static Calendar julianToDate(double jd)
	{

		double z, f, a, b, c, d, e, m, aux;
		Date date = new Date();
		jd += 0.5;
		z = Math.floor(jd);
		f = jd - z;

		if (z >= 2299161.0)
		{
		  a = Math.floor((z - 1867216.25) / 36524.25);
		  a = z + 1 + a - Math.floor(a / 4);
		}
		else
		{
		  a = z;
		}

		b = a + 1524;
		c = Math.floor((b - 122.1) / 365.25);
		d = Math.floor(365.25 * c);
		e = Math.floor((b - d) / 30.6001);
		aux = b - d - Math.floor(30.6001 * e) + f;

		Calendar calendar = new GregorianCalendar();
		calendar.setTime(date);
		calendar.set(Calendar.DAY_OF_MONTH, (int) aux);
		aux = ((aux - calendar.get(Calendar.DAY_OF_MONTH)) * 24);
		calendar.set(Calendar.HOUR_OF_DAY, (int) aux);
		calendar.set(Calendar.MINUTE, (int) ((aux - calendar.get(Calendar.HOUR_OF_DAY)) * 60));

		if (e < 13.5)
		{
		  m = e - 1;
		}
		else
		{
		  m = e - 13;
		}
		// Se le resta uno al mes por el manejo de JAVA, donde los meses empiezan en 0.
		calendar.set(Calendar.MONTH, (int) m - 1);

		if (m > 2.5)
		{
		  calendar.set(Calendar.YEAR, (int) (c - 4716));
		}
		else
		{
		  calendar.set(Calendar.YEAR, (int) (c - 4715));
		}


		return calendar;
	}


	public static String  toDateString()
	{
		GregorianCalendar gDate = new GregorianCalendar();
		String year = String.valueOf(gDate.get(Calendar.YEAR));

		String month = "";
		int m = gDate.get(Calendar.MONTH)+1;
		if (m < 10)
			month = "0";
		month += String.valueOf(m);

		String day = "";
		int d = gDate.get(Calendar.DAY_OF_MONTH);
		if (d < 10)
			day = "0";
		day += String.valueOf(d);

		String dateString = month+"-"+day+"-"+year;

		return "to_date('"+dateString+"','MM-DD-YYYY')";
	}


	public static String  toDateString(long date)
	{
		GregorianCalendar gDate = new GregorianCalendar();
		gDate.setTime(new java.util.Date(date));
		String year = String.valueOf(gDate.get(Calendar.YEAR));

		String month = "";
		int m = gDate.get(Calendar.MONTH)+1;
		if (m < 10)
			month = "0";
		month += String.valueOf(m);

		String day = String.valueOf(gDate.get(Calendar.DAY_OF_MONTH));
		String dateString = month+"-"+day+"-"+year;

		return "to_date('"+dateString+"','MM-DD-YYYY')";
	}


	public static int getMonth(String m)
		throws Exception
	{
		Integer M = (Integer)monthMap.get(m.toLowerCase());

		if (M == null)
			return 0;

		return M.intValue();
	}


        public static String getMonth2(String m)
        {
            return (String)monthMap2.get(m.toLowerCase());
        }


	public void getLocalTime()
	{
		GregorianCalendar cal =(GregorianCalendar)GregorianCalendar.getInstance();
		cal.get(Calendar.DAY_OF_MONTH);
		cal.get(Calendar.YEAR);
		cal.get(Calendar.MONTH);
	}


	public static String  getLocalIP()
		throws Exception
	{
		return InetAddress.getLocalHost().getHostAddress();
		//return InetAddress.getByName().getHostAddress();
	}


	public static void getByteUrlContent(String urlStr, String user, String password, String name)
		throws Exception
	{
		byte[]          buffer = new byte[1000000];

		FileOutputStream fos = new FileOutputStream(name);
		BufferedOutputStream bos = new BufferedOutputStream(fos);

		URL url = new URL(urlStr);
		URLConnection con = url.openConnection();

		if (user != null && password != null)
		{
			String login = user + ":" + password;

			// Encode String
			String encoding = Base64Converter.encode (login.getBytes());
			con.setRequestProperty  ("Authorization", "Basic " + encoding);
		}

		con.setDoInput(true);
		con.setAllowUserInteraction(false);
		con.setDoOutput(false);


		//read data

		int  count = 0;

		BufferedInputStream bin = new BufferedInputStream(con.getInputStream());

		for (;;)
		{
			int size = bin.read(buffer);

			if (size == -1)
				break;

			bos.write(buffer, 0, size);
		}

		bin.close();
		fos.close();
		bos.close();
	}

	public static String formatSeq(String str)
	{
		StringBuffer buffer = new StringBuffer();

		int k = 0;
		if (str.length()%60 == 0)
			k = str.length()/60;
		else
			k = str.length()/60+1;

		for (int i=0; i<k; i++)
		{
			int l = (i+1)*60;
			if (l > str.length())
				l = str.length();

			buffer.append(str.substring(i*60, l));
			buffer.append("\n");
		}

		return buffer.toString();
	}


	public static String escapeString(String str)
	{
            return str;

            /*
		StringBuffer buffer = new StringBuffer(str);

		for (int i=buffer.length()-1; i>0; i--)
		if (buffer.charAt(i) == '_')
		{
			buffer.insert(i,'\\');
		}

		return buffer.toString();
             */
	}



	public static Object[] execCommand(String cmd)
		throws Exception
	{
		Object[] objects = new Object[3];

		Process process = null;

		process = Runtime.getRuntime().exec(cmd);

		process.waitFor();

		objects[0] = new Integer(process.exitValue());
		objects[1] = new String(parseStream(process.getInputStream()));
		objects[2] = new String(parseStream(process.getErrorStream()));

		process.destroy();

		return objects;
	}


	public static Object[] execCommand(String[] cmd)
		throws Exception
	{
		Object[] objects = new Object[3];

		Process process = null;

		process = Runtime.getRuntime().exec(cmd);

		process.waitFor();

		objects[0] = new Integer(process.exitValue());
		objects[1] = new String(parseStream(process.getInputStream()));
		objects[2] = new String(parseStream(process.getErrorStream()));

		process.destroy();

		return objects;
	}


	public static byte[] parseStream(InputStream sin)
		throws Exception
	{
		Integer Size = null;
		int size = 0;
		int total = 0;

		Vector sizeVector = new Vector();
		Vector bufferVector = new Vector();

		while (true)
		{
			byte[] 	buffer = new byte[4096];

			size = sin.read(buffer);

			if (size == -1)
				break;

			sizeVector.addElement(new Integer(size));
			bufferVector.addElement(buffer);
		}

		sin.close();

		total = 0;
		for (int i=0; i<sizeVector.size(); i++)
		{
			Size = (Integer)sizeVector.elementAt(i);
			size = Size.intValue();
			total += size;
		}

		byte[] output = new byte[total];
		total = 0;

		for (int i=0; i<sizeVector.size(); i++)
		{
			Size = (Integer)sizeVector.elementAt(i);
			size = Size.intValue();

			byte[] buffer = (byte[])bufferVector.elementAt(i);
                	System.arraycopy(buffer, 0, output, (int)total, size);

			total += size;
		}


		return output;
	}



	public  static  String  formatDouble(double  d, int left, int right)
	{
		String dString = String.valueOf(d);

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

		delta = left - leftString.length();

		String patch = "";
		for (int i=0; i<delta; i++)
			patch += " ";

		return patch + leftString + rightString;
	}


        public  static  String  formatDouble(String ds, String pattern)
	{
                double d = Double.parseDouble(ds);
		DecimalFormat myFormatter = new DecimalFormat(pattern);
		String output = myFormatter.format(d);

		return output;
	}



	public  static  String  formatDouble(double d, String pattern)
	{
		DecimalFormat myFormatter = new DecimalFormat(pattern);
		String output = myFormatter.format(d);

		return output;
	}


 	public  static  String  formatDouble(double d)
	{
		String pattern = "0.#####E0";

		DecimalFormat myFormatter = new DecimalFormat(pattern);

		return myFormatter.format(d);
	}


	public  static  String  formatDouble(double d, int right)
	{

		String pattern = "#0.";
		for (int i=0; i<right; i++)
			pattern = pattern+"0";

		DecimalFormat myFormatter = new DecimalFormat(pattern);
		String output = myFormatter.format(d);

		if (true)
			return output;

		//DecimalFormat format = new DecimalFormat("#,##0.00000000");
		//FieldPosition f = new FieldPosition(0);
		//StringBuffer s = new StringBuffer();
		//String dString = format.format(Double.parseDouble(d), s, f).toString();

		String dString = String.valueOf(d);

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

	public static String formatSeq(String q0, String q1, String s0, String s1,
				String qSeq, String sSeq, String mSeq)
	{
		int qFrom = Integer.parseInt(q0);
		int qTo = Integer.parseInt(q1);

		int sFrom = Integer.parseInt(s0);
		int sTo = Integer.parseInt(s1);

		StringBuffer buffer = new StringBuffer();

		int start = 0;
		int end = 0;

		int len = qSeq.length()/60+1;
		for (int i=0; i<len; i++)
		{
			int k = i*60;
			int l = (i+1)*60;
			if (l > qSeq.length())
				l = qSeq.length();

			//quer
			if (qFrom < qTo)
			{
				start = qFrom + k;
				end = qFrom + l;
			}
			else
			{
				start = qFrom - k;
				end = qFrom - l;
			}


			buffer.append("Query:   ");
			buffer.append(formatInteger(start, 10));
			buffer.append(qSeq.substring(k, l));
			buffer.append(" ");
			buffer.append(formatInteger(end, 10));
			buffer.append("\n");

			//middle line
			buffer.append("         ");
			buffer.append("          ");
			buffer.append(mSeq.substring(k, l));
			buffer.append("\n");

			//subject
			if (sFrom < sTo)
			{
				start = sFrom + k;
				end = sFrom + l;
			}
			else
			{
				start = sFrom - k;
				end = sFrom - l;
			}

			buffer.append("Subject: ");
			buffer.append(formatInteger(start, 10));
			buffer.append(sSeq.substring(k, l));
			buffer.append(" ");
			buffer.append(formatInteger(end, 10));
			buffer.append("\n");
			buffer.append("\n");
		}

		return buffer.toString();
	}


	public static String formatInteger(int k, int n)
	{
		String str = String.valueOf(k);

		return str + patches[n-str.length()];
	}



	public static String capFirst(String txt)
		throws Exception
	{
		return txt.substring(0,1).toUpperCase()+txt.substring(1).toLowerCase();
	}

	public static String capStr(String desc)
	{
		if (desc.equals("id"))
			return "ID";

		char c1 = desc.charAt(0);
		return Character.toUpperCase(c1) + desc.substring(1);
	}


	public static String toDesc(String desc)
	{
		return desc(desc, true);
	}


	public static String toXmlDesc(String desc)
	{
		return desc(desc, false);
	}


	public static String desc(String desc, boolean spaceFlag)
	{
		String space = "";
		if (spaceFlag)
			space = " ";

		int k = desc.indexOf("_");
		if (k < 0)
			return capStr(desc);

		int l = desc.indexOf("_", k+1);
		if (l < 0)
		{
			String str1 = desc.substring(0, k);
			String str2 = desc.substring(k+1);

			return capStr(str1)+space+capStr(str2);
		}

		String str1 = desc.substring(0, k);
		String str2 = desc.substring(k+1, l);
		String str3 = desc.substring(l+1);

		return capStr(str1)+space+capStr(str2)+space+capStr(str3);
	}

	public static String[] parseTokens(String line, String del)
		throws Exception
	{
		StringTokenizer tokens = new StringTokenizer(line, del);

		String[] toks = new String[tokens.countTokens()];

		int i = 0;

		while (tokens.hasMoreTokens())
		{
			toks[i] = tokens.nextToken();
			i++;
		}

		return toks;
	}


	public static String parseBarcode(String name)
	{
		//barcode
		int k1 = name.indexOf(" ");
		int k2 = name.indexOf("-");
		int p = name.indexOf(".");
		
		if (k1 > 0 && k1 < p)
			p = k1;
		
		if (k2 > 0 && k2 < p)
			p = k2;
		
		String barcode = name.substring(0, p);

		return barcode;
	}
	
	
	public static String parseToken(String line, String pre, String suf)
		throws Exception
	{
		int k = line.indexOf(pre);
		if (k < 0)
			return null;

		int l = line.indexOf(suf,k+pre.length()+1);
		if (l<0)
			return null;

		return line.substring(k+pre.length(),l);
	}


	public static long	oracleDateToSystemTime(String date)
	{
		//2005-07-06 16:29:56.0
		//Month value is 0-based. e.g., 0 for January.

		StringTokenizer tokens = new StringTokenizer(date,"- :.");
		int year = Integer.parseInt(tokens.nextToken());
		int month = Integer.parseInt(tokens.nextToken()) -1;
		int day = Integer.parseInt(tokens.nextToken());
		int hour = Integer.parseInt(tokens.nextToken());
		int min = Integer.parseInt(tokens.nextToken());
		int sec = Integer.parseInt(tokens.nextToken());

		System.out.println(year+":"+month+":"+day+":"+hour+":"+min+":"+sec);

		GregorianCalendar gc = new GregorianCalendar(year,month,day,hour,min, sec);

		return gc.getTimeInMillis()/1000;
    }


}

