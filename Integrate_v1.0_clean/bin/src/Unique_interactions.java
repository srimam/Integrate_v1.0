package bin;
import java.lang.*;
import java.*;
import java.io.*;
import java.util.regex.*;
import java.util.ArrayList;
import java.util.*;
/*
Code to parse DISTILLER output...

usage example:java Unique_interactions multiple_motif_clusters_withArrays_operonExt.txt result.txt
*/

public class Unique_interactions
{
	public static void filter(String in, String out) throws IOException
	{
	try
	{
		//DISTILLER output
		BufferedReader b = new BufferedReader(new FileReader(in));

		//write out results
		BufferedWriter buf = new BufferedWriter(new FileWriter(out));

		ArrayList<String> interactions = new ArrayList<String>();
		ArrayList<String> genes = new ArrayList<String>();
		ArrayList<String> motifs = new ArrayList<String>();
		
		Pattern p3 = Pattern.compile("\t");
		String lines = "";
		int counter=0,counter2=0;
		while ((lines = b.readLine())!=null)
		{
			String[] result = p3.split(lines);
			if(!interactions.contains(result[0] +" "+result[1]))
			{
				interactions.add(result[0] +" "+result[1]);
				buf.write(lines);buf.newLine();
				counter++;
			}
			if(!motifs.contains(result[0]))
			{
				counter2++;
				motifs.add(result[0]);
			}
		}
		System.out.println("No. of interactions: "+counter);
		System.out.println("No. of motifs: "+counter2);
		buf.flush();buf.close();b.close();
	}
	catch(Exception e)
	{e.printStackTrace(); }
	}
}
