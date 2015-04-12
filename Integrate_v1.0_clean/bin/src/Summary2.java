package bin;
import java.lang.*;
import java.*;
import java.io.*;
import java.util.regex.*;
import java.util.ArrayList;
import java.util.*;
import java.math.*;
/*
Code to create summary of overall analysis

BY

Saheed Imam
Institute for Systems Biology
*/

public class Summary2
{

	public static void summarize()
	{
		try
		{
			//File containing final clusters
			BufferedReader b = new BufferedReader(new FileReader("Results/Integrated_modules_final.txt"));
			//File containing tomtom RegulonDB results
			BufferedReader b2 = new BufferedReader(new FileReader("Results/Motif_finding/tomtom_out2/tomtom.txt"));
			//File containing regulator predictions
			BufferedReader b4 = new BufferedReader(new FileReader("Results/Miscellaneous/Top10_regulators.txt"));

			BufferedWriter buf = new BufferedWriter(new FileWriter("Summary_final.html"));

			Pattern p = Pattern.compile("E-value = ");		
			ArrayList<String> goID = new ArrayList<String>();
			ArrayList<String> goDescription = new ArrayList<String>();
			ArrayList<String> pvalue = new ArrayList<String>();
			ArrayList<String> ClusterID = new ArrayList<String>();
			ArrayList<String> regDB = new ArrayList<String>();
			ArrayList<String> evalue = new ArrayList<String>();
			ArrayList<String> cluster_ID = new ArrayList<String>();
			ArrayList<String> cluster_holder = new ArrayList<String>();
			ArrayList<String> clusterID_reg = new ArrayList<String>();
			ArrayList<String> TFs = new ArrayList<String>();
			ArrayList<String> score = new ArrayList<String>();
			ArrayList<String> holder = new ArrayList<String>();
			
			String lines = "";
			File gsea_data = new File("Results/GSEA_analysis_final.txt");
			if(gsea_data.exists())
			{
				//File containing GSEA results
				BufferedReader b1 = new BufferedReader(new FileReader("Results/GSEA_analysis_final.txt"));
				b1.readLine();
				while ((lines = b1.readLine())!=null)//parse enrichment results
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
			while ((lines = b2.readLine())!=null)//parse tomtom results
			{
				String[] result=lines.split("\t");
				regDB.add(result[0].trim());
				cluster_ID.add(result[1].trim());
				evalue.add(result[5].trim());
			}

			lines="";
			while ((lines = b4.readLine())!=null)//parse regulator list
			{
				String[] result=lines.split("\t");
				clusterID_reg.add(result[0].trim());
				TFs.add(result[1].trim());
				score.add(result[2].trim());
			}
			lines="";
			String init="",temp="";
			int track=0;
			while ((lines = b.readLine())!=null)//parse cluster results file
			{
				String[] result = lines.split("\t");
				//System.out.println(temp);
				if(init.equals(""))
					init=result[0].trim();
				if(!init.equals(result[0].trim()))
				{
					cluster_holder.add(temp);
					track=0;
					temp="";
					init=result[0].trim();
				}
				if(track==0)
				{
					temp=result[0].trim()+":\t"+result[1].trim();
					track++;
				}
				else
					temp+="\t"+result[1].trim();
			}
			cluster_holder.add(temp);
			buf.write("<!DOCTYPE html><html>\n<table border=\"2\" bordercolor=\"grey\" align=\"center\">\n");
			buf.write("<tr>\n<th>Cluster ID</th>\n<th>Genes in cluster</th>\n<th>Motif</th>\n<th>Motif_Evalue</th>\n<th>Enriched functional categories</th>\n<th>Top predicted regulators</th>\n<th>Best hit RegulonDB</th>\n</tr>");
			for(int x=0;x<cluster_holder.size();x++)//summarize data for each cluster
			{	
				String motifScore="";
				String[] result = cluster_holder.get(x).split(":");
				BufferedReader b3 = new BufferedReader(new FileReader("Results/Motif_finding/meme_out3/"+result[0].replace(":","").trim()+"/meme.txt"));
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
				{
				buf.write("<tr>\n<th>"+result[0].trim()+"</th>\n"+"<td>"+result[1].trim().replace("\t",", ")+"</td>\n");
				buf.write("<td><img src=\"Results/Motif_finding/meme_out3/"+result[0].trim()+"/logo1.png\" alt=\"Smiley face\" height=150 width=300></td>\n");	
				buf.write("<td>"+motifScore+"</td>");				
				buf.write("<td>");
				int tracker=0;
				int counter=0;
				ArrayList<String> temp_ID = new ArrayList<String>();
				ArrayList<String> temp_Desc = new ArrayList<String>();
				ArrayList<BigDecimal> temp_val = new ArrayList<BigDecimal>();
				for(int i=0;i<ClusterID.size();i++)//enrichment data
				{
					if(ClusterID.get(i).equals(result[0]))
					{
						temp_ID.add(goID.get(i));
						temp_Desc.add(goDescription.get(i));
						temp_val.add(new BigDecimal(pvalue.get(i)));
					}
				}
				ArrayList<BigDecimal> temp_val2 = new ArrayList<BigDecimal>(temp_val);
				Collections.sort(temp_val);
				for(int i=0;i<temp_ID.size();i++)//show top 2 enriched gene sets
				{
					BigDecimal bd = temp_val.get(i).round(new MathContext(2));
					if(holder.contains(temp_ID.get(temp_val2.indexOf(temp_val.get(i)))))	
						buf.write("<b>"+temp_ID.get(temp_val2.lastIndexOf(temp_val.get(i)))+"</b><br>"+temp_Desc.get(temp_val2.lastIndexOf(temp_val.get(i)))+"<br>(p="+bd+")<br>");
					else
					{
						holder.add(temp_ID.get(temp_val2.indexOf(temp_val.get(i))));
						buf.write("<b>"+temp_ID.get(temp_val2.indexOf(temp_val.get(i)))+"</b><br>"+temp_Desc.get(temp_val2.indexOf(temp_val.get(i)))+"<br>(p="+bd+")<br>");
					}
					if(i==1) i=temp_ID.size();
					counter++;
				}
				
				if(counter==0) buf.write("None");
				buf.write("</td>\n");

				buf.write("<td>");
				int y = 0;
				for(int i=0;i<clusterID_reg.size();i++)//regulator data
				{
					if(clusterID_reg.get(i).equals(result[0]))
					{
						BigDecimal bd = new BigDecimal(score.get(i));
						bd = bd.round(new MathContext(2));
						buf.write("<b>"+TFs.get(i)+"</b> <br>(R_score="+bd+")<br>");
						y++;
					}
					if(y==3) 
						i=clusterID_reg.size();
				}
				buf.write("</td>\n");

				buf.write("<td>");
				ArrayList<String> temp_DB = new ArrayList<String>();
				temp_val = new ArrayList<BigDecimal>();
				tracker=0;
				for(int i=0;i<cluster_ID.size();i++)//tomtom data
				{
					if(cluster_ID.get(i).equals(result[0]) && tracker==0)
					{
						temp_DB.add(regDB.get(i));
						temp_val.add(new BigDecimal(evalue.get(i)));
						tracker++;
					}
					else if(cluster_ID.get(i).equals(result[0]) && tracker>0)
					{
						int loop=temp_DB.size();
						for(int j=0;j<loop;j++)
						{
							if(Double.parseDouble(evalue.get(i))<temp_val.get(j).doubleValue())
							{
								temp_DB.add(j,regDB.get(i));
								temp_val.add(j,new BigDecimal(evalue.get(i)));
							}
						}
					}
				}
				for(int i=0;i<temp_DB.size();i++)//show two best reg DB hits if available
				{
					BigDecimal bd = temp_val.get(i).round(new MathContext(2));
					buf.write("<b>"+temp_DB.get(i)+"</b> <br>(qval="+bd+")<br>");
					if(i==1) i=temp_DB.size();
				}
				buf.write("</td>\n");
				buf.write("</tr>\n");
				}	
			}
			buf.write("</table></html>");
			b.close();b2.close();b4.close();buf.flush();buf.close();
	}
	catch(Exception e)
	{e.printStackTrace(); }
	}

}
