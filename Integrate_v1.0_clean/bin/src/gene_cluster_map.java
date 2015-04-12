package bin;
import java.lang.*;
import java.*;
import java.io.*;
import java.util.regex.*;
import java.util.ArrayList;
import java.util.*;
/*
Code to be used generate two column file of geneIDs and clusterIDs

java gene_cluster_map
*/

public class gene_cluster_map
{
	public static void format() throws IOException
	{
	try
	{
		//File containing clusters from footprinting analysis
		BufferedReader b = new BufferedReader(new FileReader("Results/Clusters.txt"));
		BufferedReader b1 = new BufferedReader(new FileReader("Results/Miscellaneous/Final_clusters.txt"));

		BufferedWriter buf = new BufferedWriter(new FileWriter("Results/Miscellaneous/Phylo_Clusters.txt"));
		
		ArrayList<String> clusters = new ArrayList<String>();
		ArrayList<String> goID = new ArrayList<String>();
		String lines = "";
		while ((lines = b1.readLine())!=null)
		{	
			clusters.add(lines.trim());
		}
		while ((lines = b.readLine())!=null)
		{	
			String[] result = lines.split("\t");
			if(clusters.contains(result[0].trim().replace(":","")))
			{
				for(int i=1;i<result.length;i++)
					buf.write(result[0].trim().replace(":","")+"\t"+result[i].replace("old_","").trim()+"\n");
			}
		}
		b.close();buf.flush();buf.close();
	}
	catch(Exception e)
	{e.printStackTrace(); }
	}
}
