package bin;
import java.lang.*;
import java.*;
import java.io.*;
import java.util.regex.*;
import java.util.ArrayList;
import java.util.*;
/*
Code to parse meme.txt files from phylogenetic footprinting analysis

usage example: java MemeParse_forMast

Written by Saheed Imam
Institute for Systems biology
*/

public class MemeParse_forMast
{
	public static int parse(String target_abbr, int number_of_groups,String mast_bash_outfile, String gene_motif_map,String meme_out,String bg,String pvalue) throws IOException
	{
		int counter=1;
		try
		{
		//write out all identified motifs
		BufferedWriter buf2 = new BufferedWriter(new FileWriter(mast_bash_outfile));//mast commands
		BufferedWriter buf3 = new BufferedWriter(new FileWriter(gene_motif_map));
		
		Pattern p = Pattern.compile(" ");
		Pattern p1 = Pattern.compile("\\b"+target_abbr+"[|]\\w+\\p{Graph}+\\b");

		//Parsing meme.txt file...
		String background = "",motifName="",group="",lines = "";
		File fasta = new File("Results/Motif_finding/mast_motifs");
		deleteFolder.delete(fasta);
		fasta.mkdir();
		fasta = new File("Results/Motif_finding/mast_out");
		fasta.delete();
		fasta.mkdir();
		buf2.write("#!/bin/bash\n");
		buf2.write("# Mast bash script\n");
		for(int i=1;i<number_of_groups+1;i++)
		{
			//meme.txt file
			BufferedReader b = new BufferedReader(new FileReader(meme_out+"/group_"+i+"/meme.txt"));
			int motif=1;
			int test=0;
			lines="";
			motifName="";
			while ((lines = b.readLine())!=null)
			{
				if(lines.contains("COMMAND LINE SUMMARY")) test++;
				if(test>0)
				{
					Matcher m = p1.matcher(lines);
					if(m.find())
					{
						motifName = m.group()+"_"+motif;
					}
					if(lines.length()>10 && lines.substring(0,8).equals("DATAFILE"))
					{
						String[] result = p.split(lines);
						group= result[1];
					}
					if(lines.length()>10 && lines.substring(0,10).equals("Background"))
					{
						background = b.readLine();
					}
					if(lines.length()> 25 && lines.substring(0,8).equals("log-odds") && !motifName.equals(""))
					{
						buf2.write("mast Results/Motif_finding/mast_motifs/Motif_"+counter + " Data/IGRs/"+target_abbr+"_IGRs.txt -bfile "+bg +" -mt "+pvalue+" -ev 0.1 -hit_list >Results/Motif_finding/mast_out/Motif_"+counter+"\n");buf2.flush();
						BufferedWriter buf = new BufferedWriter(new FileWriter("Results/Motif_finding/mast_motifs/Motif_"+counter));//file containing log odds matrix
						buf.write("MEME version 4\n\nALPHABET= ACGT\n\nstrands: + -\n\nBackground letter frequencies (from meme_bfile):\n"+background+"\n\n\nMOTIF "+motifName+"\n"+lines+"\n");
						buf3.write("Motif_"+counter+"\t"+motifName+"\tgroup_"+i+"\n");buf3.flush();
						while ((lines = b.readLine())!=null)
						{
							if(lines.charAt(0)=='-')
							{ 
								motif++;
								counter++;
								break;
							}
							buf.write(lines+"\n");
						}
						buf.newLine();buf.flush();buf.close();
						motifName="";
					}
				}
			}
			b.close();
		}
		buf2.flush();buf2.close();buf3.flush();buf3.close();
	}
	catch(Exception e)
	{e.printStackTrace(); return -1;}
	return counter-1;
	}
}
