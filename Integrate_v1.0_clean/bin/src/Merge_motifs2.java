package bin;
import java.lang.*;
import java.*;
import java.io.*;
import java.util.regex.*;
import java.util.ArrayList;
import java.util.*;
/*
Code to parse results from mcl analysis, create fasta files and use meme for motif detection

usage example: java Merge_motifs2
*/

public class Merge_motifs2
{
	public static int merge(String mclOut,String fasta_directory, String target,String bg, String output) throws IOException
	{
		int counter = 1;
		try
		{
			BufferedReader b = new BufferedReader(new FileReader(mclOut));
			BufferedWriter buf1 = new BufferedWriter(new FileWriter(output));
			ArrayList<String> geneID_igr = new ArrayList<String>();
			ArrayList<String> seq = new ArrayList<String>();
			
			String lines="";
			while((lines=b.readLine())!=null)//read each line of the mcl results (i.e. each set of clusters
			{
				String[] result = lines.split("\t");
				int promoter_count=0;
				for(int i=0;i<result.length;i++)
				{
					BufferedReader b1 = new BufferedReader(new FileReader("Results/Motif_finding/FASTA2/"+result[i].trim()));
					String lines1="";
					while((lines1=b1.readLine())!=null)
					{
						if(!geneID_igr.contains(lines1.trim().replace(">","")) && promoter_count<=1500)
						{
							geneID_igr.add(lines1.trim().replace(">",""));
							seq.add(b1.readLine().trim());
							promoter_count++;
						}
						else b1.readLine();
					}
					b1.close();
				}
				BufferedWriter buf = new BufferedWriter(new FileWriter(fasta_directory+"/Cluster_"+counter));
				buf1.write("Cluster_"+counter+":\t");
				for(int i=0;i<geneID_igr.size();i++)
				{
					buf.write(">"+geneID_igr.get(i)+"\n"+seq.get(i)+"\n");
					buf1.write(geneID_igr.get(i)+"\t");
				}
				counter++;
				buf1.newLine();
				buf.flush();buf.close();
				geneID_igr.clear();
				seq.clear();
			}
			buf1.flush();buf1.close();

			//Run meme analysis on sequences
			BufferedReader reader;
			System.out.println("Running MEME..."); 
			for(int i=1;i<counter;i++)
			{
				Runtime runtime = Runtime.getRuntime();
				Process process = runtime.exec("meme "+fasta_directory+"/Cluster_"+ i+" -oc Results/Motif_finding/meme_out3/Cluster_"+i+" -dna  -mod zoops -evt 0.01 -nmotifs 1 -maxw 20 -bfile "+bg+" -revcomp -maxsize 1000000");
				reader = new BufferedReader(new InputStreamReader(process.getErrorStream()));
				lines="";
				while((lines=reader.readLine())!=null)
				{
				}
				process.waitFor();
				runtime.gc();
			}
		}
		catch(Exception e)
		{e.printStackTrace(); return -1;}
		return counter-1;
	}
}
