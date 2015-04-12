package bin;
import java.lang.*;
import java.*;
import java.io.*;
import java.util.regex.*;
import java.util.ArrayList;
import java.util.*;
/*
Code to be transform tomtom q-values

java Tomtom_parse
*/

public class Tomtom_parse
{
	public static int parse(String tomtom_out,String mcl_file) throws IOException
	{
	try
	{
		//File containing fasta format of genome
		BufferedReader b = new BufferedReader(new FileReader(tomtom_out));
		BufferedWriter buf = new BufferedWriter(new FileWriter(mcl_file));
		String lines="";
		b.readLine();
		while ((lines = b.readLine())!=null)
		{	
			String[] result = lines.split("\t");
			double temp = Math.log10(Double.parseDouble(result[5]))*-1;
			buf.write(result[0]+"\t"+result[1]+"\t"+temp+"\n");
		}
		b.close();buf.flush();buf.close();
	}
	catch(Exception e)
	{e.printStackTrace(); return -1;}
	return 0;
	}
}
