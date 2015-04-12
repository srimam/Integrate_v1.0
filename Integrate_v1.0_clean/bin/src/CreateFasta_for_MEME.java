package bin;
import java.lang.*;
import java.*;
import java.io.*;
import java.util.regex.*;
import java.util.ArrayList;
import java.util.*;
/*
Code to create fasta for running MEME

Written by Saheed Imam
Institute for Systems biology
*/

public class CreateFasta_for_MEME
{
	public static int file_prep(String IGRs, String gene_mRNA_map,String ortho_groups, String target, String bg) throws IOException
	{
		int counter = 1;
		try
		{
		//Fasta file with all input sequences
		BufferedReader b2 = new BufferedReader(new FileReader(IGRs));
		//mapping of gene IDs and mRNA names
		BufferedReader b1 = new BufferedReader(new FileReader(gene_mRNA_map));

		//Text file with OrthoMCL results
		BufferedReader b = new BufferedReader(new FileReader(ortho_groups));
		//Mapping file
		BufferedWriter buf3 = new BufferedWriter(new FileWriter("Results/Miscellaneous/Group_map.txt"));

		ArrayList<String> ID = new ArrayList<String>();
		ArrayList<String> sequences = new ArrayList<String>();
		ArrayList<String> orthologs = new ArrayList<String>();
		ArrayList<String> orthologs2 = new ArrayList<String>();
		ArrayList<String> geneID = new ArrayList<String>();
		ArrayList<String> mRNAname = new ArrayList<String>();
		ArrayList<String> temp = new ArrayList<String>();

		String lines = "";
		while ((lines = b2.readLine())!=null)//reading in input sequences...
		{
			if(lines.charAt(0)=='>')
			{
				ID.add(lines.trim().replace(">",""));
				sequences.add(b2.readLine());
			}
		}
		System.out.println("Finished reading sequences...");
		lines = "";
		while ((lines = b1.readLine())!=null)
		{
			String[] result = lines.split("\t");
			geneID.add(result[0].trim());
			mRNAname.add(result[1].trim());
		}

		int sizeTest = 0;lines=""; 
		File fasta = new File("Results/Motif_finding/FASTA");
		if(!fasta.exists())
		{
			fasta.mkdir();
			System.out.println("Making FASTA directory");
		}
		else
		{
			System.out.println("Overwriting FASTA directory");
			deleteFolder.delete(fasta);
			fasta.mkdir();
		}

		fasta = new File("Results/Motif_finding/meme_out");
		deleteFolder.delete(fasta);
		fasta.mkdir();

		BufferedWriter buf1 = new BufferedWriter(new FileWriter("bin/meme_bash"));
		buf1.write("#!/bin/bash\n");
		buf1.write("# Meme bash script\n");
		while ((lines = b.readLine())!=null)
		{	
			int test=0;
			orthologs.clear();
			orthologs2.clear();
			sizeTest=0;
			String[] result = lines.split(" ");
			for(int i=1;i<result.length;i++)//skip group ID
			{
				String[] result2 = result[i].split("[|]");
				orthologs.add(result2[1].trim());
				orthologs2.add(result[i].trim());
			}
			//check if there are at least 4 promoters of adequate size for each orthologs set and if this includes one for the target organism
			for(int i = 0; i<orthologs.size();i++)
			{
				if(mRNAname.indexOf(orthologs.get(i))>=0)
				{
					if(ID.contains(geneID.get(mRNAname.indexOf(orthologs.get(i)))))
					{
						sizeTest++;
						if(result[i+1].contains(target)) test++;
					}
				}
			}
			//if(sizeTest>=4 && test>0 && test <3)//if there are at least 4 proteins in a group and no more than 2 paralogs for the target organism
			if(sizeTest>=4 && test>0)//if there are at least 4 proteins in a group including at least one from target
			{
				buf3.write(result[0]+"\tgroup_"+counter+"\n");buf3.flush();
				//write out fasta file for each group of orthologs
				BufferedWriter buf = new BufferedWriter(new FileWriter("Results/Motif_finding/FASTA/group_"+counter+".fa"));
				buf1.write("meme Results/Motif_finding/FASTA/group_"+counter + ".fa -oc Results/Motif_finding/meme_out/group_" + counter + " -dna  -mod zoops -evt 0.01 -nmotifs 3 -maxw 20 -bfile "+bg+" -revcomp -maxsize 1000000\n");buf1.flush();
				counter++;
				for(int i = 0; i<orthologs.size();i++)
				{
					if(mRNAname.indexOf(orthologs.get(i))>=0 && ID.contains(geneID.get(mRNAname.indexOf(orthologs.get(i)))))
					{
						String[] temp2 = orthologs2.get(i).split("[|]");
						buf.write(">"+temp2[0]+"|"+ geneID.get(mRNAname.indexOf(orthologs.get(i)))+"\n");
						buf.write(sequences.get(ID.indexOf(geneID.get(mRNAname.indexOf(orthologs.get(i)))))+"\n");
					}
				}
				buf.flush();
				buf.close();
			}
		}
		System.out.println("Eligible groups : "+(counter-1));
		b.close();b2.close();buf1.flush();buf1.close();buf3.flush();buf3.close();
	}
	catch(Exception e)
	{e.printStackTrace(); return -1;}
	return counter-1;
	}

}
