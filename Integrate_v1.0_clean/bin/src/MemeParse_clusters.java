package bin;
import java.lang.*;
import java.*;
import java.io.*;
import java.util.regex.*;
import java.util.ArrayList;
import java.util.*;
/*
Code to obtain genes within clusters after meme analysis

java MemeParse_clusters
*/

public class MemeParse_clusters
{
	public static int parse(int size, String output) throws IOException
	{
	try
	{	
		BufferedWriter buf = new BufferedWriter(new FileWriter(output));
		String lines = "",test="";

		for(int i=1;i<size+1;i++)
		{
			BufferedReader b = new BufferedReader(new FileReader("Results/Motif_finding/meme_out3/Cluster_"+i+"/meme.txt"));
			ArrayList<String> temp = new ArrayList<String>();
			int test1=0;
			String temp1="Cluster_"+i+":";	
			while ((lines = b.readLine())!=null)
			{
				if(lines.contains("Motif 1 block diagrams"))
				{
					b.readLine();b.readLine();b.readLine();
					while ((lines = b.readLine())!=null)
					{
						if(lines.contains("-----------------------")) test1++;
						lines=lines.replace("                                 ","\t").replace("                                ","\t").replace("                               ","\t").replace("                              ","\t").replace("                             ","\t").replace("                            ","\t").replace("                           ","\t").replace("                          ","\t").replace("                         ","\t").replace("                        ","\t").replace("                       ","\t").replace("                      ","\t").replace("                     ","\t").replace("                    ","\t").replace("                   ","\t").replace("                  ","\t").replace("                 ","\t").replace("                ","\t").replace("               ","\t").replace("              ","\t").replace("             ","\t").replace("            ","\t").replace("           ","\t").replace("          ","\t").replace("         ","\t").replace("        ","\t").replace("       ","\t").replace("      ","\t").replace("     ","\t").replace("    ","\t").replace("   ","\t").replace("  ","\t").replace(" ","\t");
						String[] result = lines.split("\t");
						if(!temp.contains(result[0].trim()) && test1==0)
						{
							temp.add(result[0].trim());
							temp1=temp1+"\t"+result[0].trim();
						}
					}
				}
			}
			if(temp.size()>0)
				buf.write(temp1+"\n");
			b.close();
		}
		buf.flush();buf.close();
	}
	catch(Exception e)
	{e.printStackTrace(); return -1;}
	return 0;
	}

}
