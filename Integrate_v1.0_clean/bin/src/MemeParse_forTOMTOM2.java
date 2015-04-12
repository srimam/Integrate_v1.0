package bin;
import java.lang.*;
import java.*;
import java.io.*;
import java.util.regex.*;
import java.util.ArrayList;
import java.util.*;
/*
Code to parse meme.txt files from phylogenetic footprinting analysis for Tomtom comparisons

usage example: MemeParse_forTOMTOM2
*/

public class MemeParse_forTOMTOM2
{
	public static int parse(int number_of_groups, String meme_out, String output) throws IOException
	{
		try
		{
		BufferedWriter buf2 = new BufferedWriter(new FileWriter(output));
		Pattern p = Pattern.compile(" ");

		//Parsing meme.txt file...
		String background = "",lines = "";
		for(int i=1;i<number_of_groups+1;i++)
		{
			//meme.txt file
			BufferedReader b = new BufferedReader(new FileReader(meme_out+"/Cluster_"+i+"/meme.txt"));
			while ((lines = b.readLine())!=null)
			{
				if(lines.length()>10 && lines.substring(0,10).equals("Background"))
				{
					background = b.readLine();
				}
				if(lines.length()> 25 && lines.substring(0,18).equals("letter-probability"))
				{
					buf2.write("MEME version 4\n\nALPHABET= ACGT\n\nstrands: + -\n\nBackground letter frequencies (from meme_bfile):\n"+background+"\n\n\nMOTIF Cluster_"+i+"\n"+lines+"\n");
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
