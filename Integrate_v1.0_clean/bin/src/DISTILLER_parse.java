package bin;
import java.lang.*;
import java.*;
import java.io.*;
import java.util.regex.*;
import java.util.ArrayList;
import java.util.*;
/*
Code to parse DISTILLER output...

DISTILLER_parse
*/

public class DISTILLER_parse
{
	public static int parse() throws IOException
	{
	try
	{
		//DISTILLER output
		BufferedReader b = new BufferedReader(new FileReader("Results/DISTILLER/DISTILLER_modules.txt"));

		//File index to geneID mapping
		BufferedReader b2 = new BufferedReader(new FileReader("Results/DISTILLER/Index_to_geneID.txt"));

		//File index to clusterID mapping
		BufferedReader b3 = new BufferedReader(new FileReader("Results/DISTILLER/Index_to_clusterID.txt"));

		//write out results
		BufferedWriter buf = new BufferedWriter(new FileWriter("Results/DISTILLER/DISTILLER_modules_processed.txt"));

		ArrayList<String> ID = new ArrayList<String>();
		ArrayList<String> geneID = new ArrayList<String>();
		ArrayList<String> IDclust = new ArrayList<String>();
		ArrayList<String> clusterID = new ArrayList<String>();
		ArrayList<String> genes = new ArrayList<String>();
		ArrayList<String> motifs = new ArrayList<String>();
		ArrayList<String> uniqueMotifs = new ArrayList<String>();
		ArrayList<String> uniqueGenes = new ArrayList<String>();

		Pattern p3 = Pattern.compile("\t");
		Pattern p4 = Pattern.compile("\\[");
		Pattern p5 = Pattern.compile(",");

		String lines = "";
		while ((lines = b2.readLine())!=null)//obtain genes to indices mapping
		{
			String[] result = p3.split(lines);
			ID.add(result[0]);
			geneID.add(result[1]);
		}
		
		lines = "";
		while ((lines = b3.readLine())!=null)//obtain cluster to indices mapping
		{
			String[] result = p3.split(lines);
			IDclust.add(result[0]);
			clusterID.add(result[1]);
		}

		lines="";
		while ((lines = b.readLine())!=null)//parse DISTILLER output
		{
			String[] result3 = {""};
			if(lines.length()>10 && lines.substring(0,5).equals("items"))
			{
				String[] result = p4.split(lines);
				String[] result2 = p5.split(result[1].replace("] + 1;",""));
				for(int i=0;i<result2.length;i++)
				{
					genes.add(geneID.get(ID.indexOf(result2[i].trim())+1));
					if(!uniqueGenes.contains(geneID.get(ID.indexOf(result2[i].trim())+1)))
						uniqueGenes.add(geneID.get(ID.indexOf(result2[i].trim())+1));
				}
				lines = b.readLine();
				if(lines.length()>10 && lines.substring(0,7).equals("tidset1"))
				{
					result = p4.split(lines);
					result2 = p5.split(result[1].replace("] + 1;",""));
					for(int i=0;i<result2.length;i++)
					{
						motifs.add(clusterID.get(IDclust.indexOf(result2[i].trim())+1));
						if(!uniqueMotifs.contains(clusterID.get(IDclust.indexOf(result2[i].trim())+1)))
							uniqueMotifs.add(clusterID.get(IDclust.indexOf(result2[i].trim())+1));
						
					}
				}
				lines = b.readLine();
				if(lines.length()>10 && lines.substring(0,10).equals("boxtidset1"))
				{
					result3 = p4.split(lines);
				}
			}
			for(int j=0;j<motifs.size();j++)
			{
				for(int k=0;k<genes.size();k++)
				{
					buf.write(motifs.get(j)+"\t"+genes.get(k)+"\t"+result3[1].replace("] + 1;","")+"\n");
				}
			}
			genes.clear();
			motifs.clear();
		}
		System.out.println("Number of unique motifs: "+ uniqueMotifs.size());
		System.out.println("Number of unique genes: "+ uniqueGenes.size());
		buf.flush();buf.close();
		b.close();b2.close();b3.close();
		
		Unique_interactions.filter("Results/DISTILLER/DISTILLER_modules_processed.txt","Results/DISTILLER/DISTILLER_modules_processed_UNIQUE.txt");
	}
	catch(Exception e)
	{e.printStackTrace(); return -1;}
	return 0;
	}
}
