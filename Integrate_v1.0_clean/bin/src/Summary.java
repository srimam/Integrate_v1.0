package bin;
import java.lang.*;
import java.*;
import java.io.*;
import java.util.regex.*;
import java.util.ArrayList;
import java.util.*;
import java.math.*;
/*
Code to create summary of comparative genomics analysis

BY

Saheed Imam
Institute for Systems Biology
*/

public class Summary
{

	public static void summarize()
	{
		try
		{
			//File containing final clusters
			BufferedReader b = new BufferedReader(new FileReader("Results/Clusters.txt"));
			//File containing tomtom RegulonDB results
			BufferedReader b2 = new BufferedReader(new FileReader("Results/Motif_finding/tomtom_out2/tomtom.txt"));

			BufferedWriter buf = new BufferedWriter(new FileWriter("Summary_phylo.html"));
			BufferedWriter buf2 = new BufferedWriter(new FileWriter("Results/Miscellaneous/Final_clusters.txt"));

			Pattern p = Pattern.compile("E-value = ");		
			ArrayList<String> goID = new ArrayList<String>();
			ArrayList<String> goDescription = new ArrayList<String>();
			ArrayList<String> pvalue = new ArrayList<String>();
			ArrayList<String> ClusterID = new ArrayList<String>();
			ArrayList<String> regDB = new ArrayList<String>();
			ArrayList<String> evalue = new ArrayList<String>();
			ArrayList<String> cluster_ID = new ArrayList<String>();
			ArrayList<String> clusters = new ArrayList<String>();
			ArrayList<BigDecimal> evalues = new ArrayList<BigDecimal>();
			ArrayList<String> holder = new ArrayList<String>();

			String lines = "";
			File gsea_data = new File("Results/GSEA_analysis_phylo.txt");
			if(gsea_data.exists())
			{
				//File containing GSEA results
				BufferedReader b1 = new BufferedReader(new FileReader("Results/GSEA_analysis_phylo.txt"));
				b1.readLine();
				while ((lines = b1.readLine())!=null)//GSEA analysis
				{
					if(!lines.trim().equals(""))
					{
					
						String[] result=lines.split("\t");
						goID.add(result[0].trim());
						goDescription.add(result[1].trim());
						pvalue.add(result[4].trim());
						ClusterID.add(result[5].trim().replace(":","").trim());
					}
				}
				b1.close();
			}
			lines="";
			b2.readLine();
			while ((lines = b2.readLine())!=null)//tomtom analysis
			{
				String[] result=lines.split("\t");
				regDB.add(result[0].trim());
				cluster_ID.add(result[1].trim());
				evalue.add(result[5].trim());
			}
			
			lines="";			
			while ((lines = b.readLine())!=null)
			{	
				String motifScore="";
				String[] result = lines.split(":");
				BufferedReader b3 = new BufferedReader(new FileReader("Results/Motif_finding/meme_out3/"+result[0].trim()+"/meme.txt"));
				String lines2="";
				while((lines2 = b3.readLine())!=null)
				{
					Matcher m = p.matcher(lines2);
					if(m.find())
					{
						String[] result1=lines2.split(" ");
						motifScore=result1[result1.length-1];
					}
				}
				if(!motifScore.equals(""))
				{
					clusters.add(lines);
					evalues.add(new BigDecimal(motifScore.trim()));
				}
			}
			ArrayList<BigDecimal> top_evalues = new ArrayList<BigDecimal>(evalues);
			Collections.sort(top_evalues);
			buf.write("<!DOCTYPE html><html>\n<table border=\"2\" bordercolor=\"grey\" align=\"center\">\n");
			buf.write("<tr>\n<th>Cluster ID</th>\n<th>Genes in cluster</th>\n<th>Motif</th>\n<th>Motif_RC</th>\n<th>Motif_Evalue</th>\n<th>Enriched functional categories</th>\n<th>Best hit RegulonDB</th>\n</tr>");
			int counts=0;
			//while ((lines = b.readLine())!=null)
			for(int x=0;x<top_evalues.size();x++)
			{	
				for(int y=0;y<evalues.size();y++)
				{
					String[] result = clusters.get(y).split(":");
					if(top_evalues.get(x)==evalues.get(y) && !holder.contains(result[0].trim()) && top_evalues.get(x).doubleValue() <0.0001)
					{
						counts++;	
						String motifScore="";
						BufferedReader b3 = new BufferedReader(new FileReader("Results/Motif_finding/meme_out3/"+result[0].trim()+"/meme.txt"));
						holder.add(result[0].trim());
						String lines2="";
						while((lines2 = b3.readLine())!=null)
						{
							Matcher m = p.matcher(lines2);
							if(m.find())
							{
								String[] result1=lines2.split(" ");
								motifScore=result1[result1.length-1];
							}
						}
						b3.close();
						if(!motifScore.equals(""))
						{buf2.write(result[0].trim()+"\n");
						buf.write("<tr>\n<th>"+result[0].trim()+"</th>\n"+"<td>"+result[1].trim().replace("\t",", ")+"</td>\n");
						buf.write("<td><img src=\"Results/Motif_finding/meme_out3/"+result[0].trim()+"/logo1.png\" alt=\"Smiley face\" height=150 width=300></td>\n");
						buf.write("<td><img src=\"Results/Motif_finding/meme_out3/"+result[0].trim()+"/logo_rc1.png\" alt=\"Smiley face\" height=150 width=300></td>\n");	
						buf.write("<td>"+motifScore+"</td>");				
						buf.write("<td>");
						int tracker=0;
						int counter=0;
						ArrayList<String> temp_ID = new ArrayList<String>();
						ArrayList<String> temp_Desc = new ArrayList<String>();
						ArrayList<String> temp_val = new ArrayList<String>();
						for(int i=0;i<ClusterID.size();i++)
						{
							tracker=0;
							if(ClusterID.get(i).equals(result[0]))
							{
								int loop=temp_ID.size();
								for(int j=0;j<loop;j++)
								{
									if(Double.parseDouble(pvalue.get(i))<Double.parseDouble(temp_val.get(j)) || Double.parseDouble(pvalue.get(i))==Double.parseDouble(temp_val.get(j)))
									{ 
										temp_ID.add(j,goID.get(i));
										temp_Desc.add(j,goDescription.get(i));
										temp_val.add(j,pvalue.get(i));
										tracker++;
										j=loop;
									}
								}
							}
							if(ClusterID.get(i).equals(result[0]) && tracker==0)
							{
								temp_ID.add(goID.get(i));
								temp_Desc.add(goDescription.get(i));
								temp_val.add(pvalue.get(i));
							}
						}
						for(int i=0;i<temp_ID.size();i++)
						{
							BigDecimal bd = new BigDecimal(temp_val.get(i));
							bd = bd.round(new MathContext(2));
							buf.write("<b>"+temp_ID.get(i)+"</b><br>"+temp_Desc.get(i)+"<br>(p="+bd+")<br>");
							if(i==1) i=temp_ID.size();
							counter++;
						}
				
						if(counter==0) buf.write("None");
						buf.write("</td>\n");
						buf.write("<td>");

						ArrayList<String> temp_DB = new ArrayList<String>();
						temp_val = new ArrayList<String>();
						for(int i=0;i<cluster_ID.size();i++)
						{
							tracker=0;
							if(cluster_ID.get(i).equals(result[0]))
							{
								int loop=temp_DB.size();
								for(int j=0;j<loop;j++)
								{
									if(Double.parseDouble(evalue.get(i))<Double.parseDouble(temp_val.get(j)))
									{
										temp_DB.add(j,regDB.get(i));
										temp_val.add(j,evalue.get(i));
										j=loop;
										tracker++;
									}
								}
							}
							if(cluster_ID.get(i).equals(result[0]) && tracker==0)
							{
								temp_DB.add(regDB.get(i));
								temp_val.add(evalue.get(i));
							}
						}
						for(int i=0;i<temp_DB.size();i++)
						{
							BigDecimal bd = new BigDecimal(temp_val.get(i));
							bd = bd.round(new MathContext(2));
							buf.write("<b>"+temp_DB.get(i)+"</b> <br>(qval="+bd+")<br>");
							if(i==1) i=temp_DB.size();
						}
						buf.write("</td>\n");
						buf.write("</tr>\n");
						}
				
				
					}
				}
				if(counts==400)//select no more than top 400 motifs
					break;
			}
			buf.write("</table></html>");
			b.close();b2.close();buf.flush();buf.close();buf2.flush();buf2.close();
	}
	catch(Exception e)
	{e.printStackTrace(); }
	}

}
