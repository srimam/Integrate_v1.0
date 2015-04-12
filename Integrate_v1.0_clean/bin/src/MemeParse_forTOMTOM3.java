package bin;
import java.lang.*;
import java.*;
import java.io.*;
import java.util.regex.*;
import java.util.ArrayList;
import java.util.*;
/*
Code to parse meme.txt files for clusters that make it through distiller analysis for Tomtom comparisons

usage example: MemeParse_forTOMTOM3
*/

public class MemeParse_forTOMTOM3
{
	public static int parse(String distiller_results,String meme_out) throws IOException
	{
		try
		{
		//Read in DISTILLER results
		BufferedReader b1 = new BufferedReader(new FileReader(distiller_results));
		BufferedWriter buf2 = new BufferedWriter(new FileWriter("Results/DISTILLER/DISTILLER_motifs_TOMTOM.txt"));
		Pattern p = Pattern.compile(" ");
		ArrayList<String> clusters = new ArrayList<String>();

		String lines="";
		while((lines=b1.readLine())!=null)
		{
			String[] result = lines.split("\t");
			if(!clusters.contains(result[0].trim()))//obtain unique list of clusters
				clusters.add(result[0]);
		}

		//Parsing meme.txt file...
		String background = "";
		lines="";
		for(int i=0;i<clusters.size();i++)
		{
			//meme.txt file
			BufferedReader b = new BufferedReader(new FileReader(meme_out+clusters.get(i)+"/meme.txt"));
			int motif=1;
			while ((lines = b.readLine())!=null)
			{
				if(lines.length()>10 && lines.substring(0,10).equals("Background"))
				{
					background = b.readLine();
				}
				if(lines.length()> 25 && lines.substring(0,18).equals("letter-probability"))
				{
					buf2.write("MEME version 4\n\nALPHABET= ACGT\n\nstrands: + -\n\nBackground letter frequencies (from meme_bfile):\n"+background+"\n\n\nMOTIF "+clusters.get(i)+"\n"+lines+"\n");
					while ((lines = b.readLine())!=null)
					{
						if(lines.charAt(0)=='-')
							break;
						buf2.write(lines+"\n");
					}
					buf2.newLine();buf2.newLine();
				}
			}
			b.close();
		}
		buf2.flush();buf2.close();
	}
	catch(Exception e)
	{e.printStackTrace(); return -1;}
	return 0;
	}
}
