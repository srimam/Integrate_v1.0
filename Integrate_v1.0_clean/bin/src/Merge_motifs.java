package bin;
import java.lang.*;
import java.*;
import java.io.*;
import java.util.regex.*;
import java.util.ArrayList;
import java.util.*;
/*
Code to parse results from mcl analysis, create fasta files and use meme for motif detection

usage example: java Merge_motifs
*/

public class Merge_motifs
{
	public static int merge(String mclOut,String fasta_directory, String target,String bg, String output) throws IOException
	{
		int counter = 1;
		try
		{
			BufferedReader b = new BufferedReader(new FileReader(mclOut));
			BufferedReader b2 = new BufferedReader(new FileReader("Data/IGRs/"+target+"_IGRs.txt"));
			BufferedWriter buf1 = new BufferedWriter(new FileWriter(output));
			ArrayList<String> geneID = new ArrayList<String>();
			ArrayList<Integer> start = new ArrayList<Integer>();
			ArrayList<Integer> end = new ArrayList<Integer>();
			ArrayList<String> geneID_igr = new ArrayList<String>();
			ArrayList<String> seq = new ArrayList<String>();
			
			String lines="";
			BufferedReader b1;
			while((lines=b2.readLine())!=null)//get sequences for target organism
			{
				geneID_igr.add(lines.trim().replace(">",""));
				seq.add(b2.readLine().trim());
			}
			lines="";
			while((lines=b.readLine())!=null)//read each line of the mcl results (i.e. each set of clusters
			{
				String[] result = lines.split("\t");
				int promoter_count=0;
				for(int i=0;i<result.length;i++)
				{
					b1 = new BufferedReader(new FileReader("Results/Motif_finding/mast_out/"+result[i].trim()));
					String lines1="";
					while((lines1=b1.readLine())!=null)
					{
						if(!lines1.contains("#"))
						{
							String temp = lines1.replace("      ","\t").replace("     ","\t").replace("    ","\t").replace("   ","\t").replace("  ","\t").replace(" ","\t");
							String[] result1 = temp.split("\t");
							if(Double.parseDouble(result1[4])>750.0 && !geneID.contains(result1[0].trim()) && promoter_count<=1500)// only consider motifs with a match score of upto 750 and only consider up to 1500 sequences so MEME doesn't take forever
							{
								geneID.add(result1[0].trim());
								start.add(Integer.parseInt(result1[2])-5);
								end.add(Integer.parseInt(result1[3])+5);
								promoter_count++;
							}
						}
					}
					b1.close();
				}
				if(geneID.size()>1)
				{
					BufferedWriter buf = new BufferedWriter(new FileWriter(fasta_directory+"/Cluster_"+counter));
					buf1.write("Cluster_"+counter+":\t");
					for(int i=0;i<geneID.size();i++)
					{
						int length = seq.get(geneID_igr.indexOf(geneID.get(i))).length();
						int begin=start.get(i)-1;
						int finish=end.get(i);
						if(begin<0) begin=0;
						if(finish>length) finish=length;
						buf.write(">"+geneID.get(i)+"\n"+seq.get(geneID_igr.indexOf(geneID.get(i))).substring(begin,finish)+"\n");
						buf1.write(geneID.get(i)+"\t");
					}
					counter++;
					buf1.newLine();
					buf.flush();buf.close();
				}
				geneID.clear();
				start.clear();
				end.clear();
			}
			buf1.flush();buf1.close();

			//Run meme analysis on sequences
			BufferedReader reader;
			System.out.println("Running MEME..."); 
			for(int i=1;i<counter;i++)
			{
				Runtime runtime = Runtime.getRuntime();
				Process process = runtime.exec("meme "+fasta_directory+"/Cluster_"+ i+" -oc Results/Motif_finding/meme_out2/Cluster_"+i+" -dna  -mod zoops -evt 0.01 -nmotifs 1 -maxw 20 -bfile "+bg+" -revcomp -maxsize 1000000");
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
