package bin;
import java.lang.*;
import java.*;
import java.io.*;
import java.util.regex.*;
import java.util.ArrayList;
import java.util.*;
/*
Code to sort predicted interactions for AUPR analysis

BY

Written by Saheed Imam
Institute for Systems biology
*/

public class AUPR_preprocessing
{
	public static int sortInteractionList() throws IOException
	{
	try
	{

		//Final modules
		BufferedReader b = new BufferedReader(new FileReader("Results/Integrated_modules_final.txt"));
		
		//Predicted regulators
		BufferedReader b1 = new BufferedReader(new FileReader("Results/Miscellaneous/Top10_regulators.txt"));

		BufferedWriter buf = new BufferedWriter(new FileWriter("Results/Final_sorted_interactions.txt"));
		
		ArrayList<String> clusterID_unique = new ArrayList<String>();
		ArrayList<String> cluster_gene_ID = new ArrayList<String>();
		ArrayList<String> cluster_gene_ID_meme = new ArrayList<String>();
		ArrayList<String> geneID = new ArrayList<String>();
		ArrayList<Double> score = new ArrayList<Double>();
		ArrayList<String> clusterID_reg = new ArrayList<String>();
		ArrayList<Double> reg_score = new ArrayList<Double>();
		ArrayList<String> TF = new ArrayList<String>();
		ArrayList<String> TF_target = new ArrayList<String>();
		ArrayList<String> cluster_final = new ArrayList<String>();
		ArrayList<Double> final_score = new ArrayList<Double>();
		String lines = "",test="";
		while ((lines = b.readLine())!=null)//modules
		{	
			String[] result = lines.split("\t");
			cluster_gene_ID.add(result[0]+"\t"+result[1]);
			if(!clusterID_unique.contains(result[0]))
				clusterID_unique.add(result[0]);		
		}
		b.close();
		lines="";
		while ((lines = b1.readLine())!=null)//regulator predictions
		{	
			String[] result = lines.split("\t");
			if(!clusterID_reg.contains(result[0].trim()))//top regulator data
			{
				clusterID_reg.add(result[0].trim());
				reg_score.add(Double.parseDouble(result[2].trim()));
				TF.add(result[1].trim());
				//buf.write(result[0]+"\t"+result[1]+"\t"+result[2]+"\n");
			}		
		}

		for(int i=0;i<clusterID_unique.size();i++)
		{
			b = new BufferedReader(new FileReader("Results/Motif_finding/meme_out3/"+clusterID_unique.get(i)+"/meme.txt"));
			ArrayList<String> temp = new ArrayList<String>();
			int test1=0;
			while ((lines = b.readLine())!=null)
			{
				if(lines.contains("Combined block diagrams: non-overlapping sites"))
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
							score.add(Math.log10(Double.parseDouble(result[1].trim()))*-1);
							cluster_gene_ID_meme.add(clusterID_unique.get(i) +"\t"+result[0].trim());
							//buf.write(clusterID_unique.get(i)+"\t"+result[0].trim()+"\t"+(Math.log10(Double.parseDouble(result[1].trim()))*-1)+"\n");
						}
					}
				}
			}
			b.close();
		}

		double interaction_score=0.0;
		for(int i=0;i<cluster_gene_ID.size();i++)
		{
			
			if(cluster_gene_ID_meme.indexOf(cluster_gene_ID.get(i))>=0)
			{
				interaction_score=Math.sqrt((reg_score.get(clusterID_reg.indexOf(cluster_gene_ID.get(i).split("\t")[0]))*score.get(cluster_gene_ID_meme.indexOf(cluster_gene_ID.get(i)))));
				TF_target.add(TF.get(clusterID_reg.indexOf(cluster_gene_ID.get(i).split("\t")[0]))+"\t"+cluster_gene_ID.get(i).split("\t")[1]);
				cluster_final.add(cluster_gene_ID.get(i).split("\t")[0]);
				final_score.add(interaction_score);
				//buf.write(TF.get(clusterID_reg.indexOf(cluster_gene_ID.get(i).split("\t")[0]))+"\t"+cluster_gene_ID.get(i).split("\t")[1]+"\t"+interaction_score +"\n");
			}
			else
			{
				//buf.write(TF.get(clusterID_reg.indexOf(cluster_gene_ID.get(i).split("\t")[0]))+"\t"+cluster_gene_ID.get(i).split("\t")[1]+"\t"+interaction_score+"\n");
				TF_target.add(TF.get(clusterID_reg.indexOf(cluster_gene_ID.get(i).split("\t")[0]))+"\t"+cluster_gene_ID.get(i).split("\t")[1]);
				cluster_final.add(cluster_gene_ID.get(i).split("\t")[0]);
				final_score.add(interaction_score);
			}
		}
		ArrayList<Double> final_score_sorted = new ArrayList<Double>(final_score);
		Collections.sort(final_score_sorted);
		ArrayList<String> holder = new ArrayList<String>();
		for(int i=final_score_sorted.size()-1;i>=0;i--)
		{
			for(int j=0;j<final_score.size();j++)
			{
				if(final_score_sorted.get(i)==final_score.get(j) && !holder.contains(TF_target.get(j)+"\t"+cluster_final.get(j)))
				{
					holder.add(TF_target.get(j)+"\t"+cluster_final.get(j));
					buf.write(TF_target.get(j)+"\t"+final_score.get(j)+"\t"+cluster_final.get(j)+"\n");
				}
			}
		}
		b1.close();buf.flush();buf.close();
	}
	catch(Exception e)
	{e.printStackTrace(); return -1;}
	return 0;
	}

}
