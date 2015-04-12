package bin;
import java.lang.*;
import java.*;
import java.io.*;
import java.util.regex.*;
import java.util.ArrayList;
import java.util.*;
/*
Code to operon extend predicted interactions


Written by Saheed Imam
Institute for Systems biology
*/

public class DISTILLER_extend_corr_operons
{
	public static int extend(String expression, String operonStructure) throws IOException
	{
	try
	{
		//Read in expression matrix
		BufferedReader b = new BufferedReader(new FileReader("Data/Expression/"+expression));

		//Read in file with interactions
		BufferedReader b2 = new BufferedReader(new FileReader("Results/DISTILLER/DISTILLER_modules_processed_UNIQUE_seedExtended.txt"));

		//Read in file with operon structure
		BufferedReader b4 = new BufferedReader(new FileReader("Data/TFs.txt"));

		//Read in cluster from footprinting analysis
		BufferedReader b5 = new BufferedReader(new FileReader("Results/Miscellaneous/Phylo_Clusters.txt"));

		//write out results
		BufferedWriter buf = new BufferedWriter(new FileWriter("Results/Miscellaneous/temp_modules_final.txt"));

		ArrayList<String> conditions = new ArrayList<String>();
		ArrayList<String> geneID = new ArrayList<String>();
		ArrayList<String> motif_ID = new ArrayList<String>();
		ArrayList<String> column1 = new ArrayList<String>();
		ArrayList<String> column2 = new ArrayList<String>();
		ArrayList<String> operons = new ArrayList<String>();
		ArrayList<String> geneID3 = new ArrayList<String>();
		ArrayList<String> TFs = new ArrayList<String>();
		ArrayList<ArrayList<String>> TFs_per_cluster = new ArrayList<ArrayList<String>>();
		ArrayList<String> motif_ID2 = new ArrayList<String>();
		ArrayList<String> geneID_temp = new ArrayList<String>();

		Pattern p3 = Pattern.compile("\t");
		Pattern p5 = Pattern.compile(",");
		
		String lines = "";
		b.readLine();//skip header
		int rows=0,cols=0;
		while ((lines = b.readLine())!=null)//identify all genes with available expression data
		{
			String[] result = lines.split("\t");
			if(cols==0)
				cols=result.length-1;
			geneID3.add(result[0]);
			rows++;
		}
		double[][] expression_mat = new double[rows][cols];//hold expression data
		b.close();
		b = new BufferedReader(new FileReader("Data/Expression/"+expression));
		lines="";
		int rowCount=0;
		b.readLine();
		String all_cols="0";
		while ((lines = b.readLine())!=null)//place expression data into matrix
		{
			String[] result = lines.split("\t");
			for(int i=1;i<result.length;i++)
			{
				expression_mat[rowCount][i-1]=Double.parseDouble(result[i]);
				if(rowCount==0)
					all_cols=all_cols+","+(i-1);
			}
			rowCount++;
		}

		lines="";
		while ((lines = b4.readLine())!=null)//TF list
		{
			String[] result = lines.split("\t");
			TFs.add(result[0].trim());
		}

		lines="";
		String initialmotif="";
		ArrayList<String> holder = new ArrayList<String>();
		while ((lines = b5.readLine())!=null)//phylo results
		{
			String[] result = lines.split("\t");
			if(initialmotif.equals("")) initialmotif=result[0];
			if(TFs.contains(result[1].trim())) holder.add(result[1].trim());
			if(!initialmotif.equals(result[0]))
			{
				TFs_per_cluster.add(holder);
				motif_ID2.add(initialmotif);
				initialmotif=result[0];
				holder = new ArrayList<String>();
			}
		}
		TFs_per_cluster.add(holder);
		motif_ID2.add(initialmotif);

		lines="";
		initialmotif="";
		while ((lines = b2.readLine())!=null)//DISTILLER results updated with TFs that might have been dropped in DISTILLER analysis due not so tight expression
		{
			String[] result = lines.split("\t");

			if(initialmotif.equals("")) initialmotif=result[0];
			if(initialmotif!=result[0])
			{
				for(int i=0;i<TFs_per_cluster.get(motif_ID2.indexOf(initialmotif)).size();i++)
				{
					if(!geneID_temp.contains(TFs_per_cluster.get(motif_ID2.indexOf(initialmotif)).get(i)))
					{
						geneID.add(TFs_per_cluster.get(motif_ID2.indexOf(initialmotif)).get(i));
						motif_ID.add(initialmotif);
						conditions.add(all_cols);
					}
				}
				initialmotif=result[0];
			}
			motif_ID.add(result[0]);
			geneID.add(result[1]);
			conditions.add(result[2]);
			geneID_temp.add(result[1]);
		}

		lines="";
		String test="";
		int counter=0;
		File operons_struct = new File("Data/Operons/"+operonStructure);
		if(operons_struct.exists())
		{
			//Read in file with operon structure
			BufferedReader b3 = new BufferedReader(new FileReader("Data/Operons/"+operonStructure));
			while ((lines = b3.readLine())!=null)//operon file
			{
				String[] result =  lines.split("\t");
				if(result[6].trim().equals("TRUE"))
				{
					if(!test.equals(result[2].trim()))
					{
						counter++;
						test=result[3].trim();
					}
					column1.add(result[2]);
					column2.add(result[3]);
					operons.add("operon_"+counter);
					test=result[3].trim();
				}
			}
			b3.close();
		}
		String operon="";
		ArrayList<String> temp = new ArrayList<String>();
		for(int i=0;i<geneID.size();i++)
		{
			temp.clear();
			buf.write(motif_ID.get(i) + "\t" + geneID.get(i) + "\t" + conditions.get(i));buf.newLine();buf.flush();
			if(column1.contains(geneID.get(i)) || column2.contains(geneID.get(i)))//if gene in an operon
			{
				if(column1.indexOf(geneID.get(i))>=0)
					operon = operons.get(column1.indexOf(geneID.get(i)));
				else
					operon = operons.get(column2.indexOf(geneID.get(i)));
				
				for(int j = 0;j<column1.size();j++)
				{
					if(operons.get(j).equals(operon))
					{
						if(!temp.contains(column1.get(j)))
							temp.add(column1.get(j));
						if(!temp.contains(column2.get(j)))
							temp.add(column2.get(j));
					}
				}			
				for(int j=0;j<temp.size();j++)
				{
					if(geneID3.contains(temp.get(j)) && geneID3.contains(geneID.get(i)) && !temp.get(j).equals(geneID.get(i)))//if expression data exists for both genes
					{
						String[] result = p5.split(conditions.get(i));
						double sum=0.0;
						double r=0.0;
						double product=0.0;
						double sumX=0.0,sumY=0.0;
						for(int k=0;k<result.length;k++)
						{
							int row1 = geneID3.indexOf(temp.get(j));
							int row2 = geneID3.indexOf(geneID.get(i));
							int col = Integer.parseInt(result[k].trim());
							sumX+=expression_mat[row1][col];
							sumY+=expression_mat[row2][col];
						}
						double meanX =sumX/result.length;
						double meanY =sumY/result.length;
						double X=0.0,Y=0.0;						
						for(int k=0;k<result.length;k++)
						{
							int row1 = geneID3.indexOf(temp.get(j));
							int row2 = geneID3.indexOf(geneID.get(i));
							int col = Integer.parseInt(result[k].trim());
							X+=Math.pow((expression_mat[row1][col]-meanX),2);
							Y+=Math.pow((expression_mat[row2][col]-meanY),2);
						}
						double sdX = Math.sqrt((X/(result.length-1)));
						double sdY = Math.sqrt((Y/(result.length-1)));
						for(int k=0;k<result.length;k++)
						{
							int row1 = geneID3.indexOf(temp.get(j));
							int row2 = geneID3.indexOf(geneID.get(i));
							int col = Integer.parseInt(result[k].trim());
							double x = (expression_mat[row1][col]-meanX)/sdX;
							double y = (expression_mat[row2][col]-meanY)/sdY;
							product=x*y;
							sum+=product;
						}
						r=sum/(result.length-1);
						if(Math.abs(r)>=0.5)
						{
							buf.write(motif_ID.get(i) + "\t" + temp.get(j) + "\t" + conditions.get(i));buf.newLine();buf.flush();
						}
					}
				}
			}
		}		
		buf.flush();buf.close();
		b.close();b2.close();b4.close();b5.close();

		Unique_interactions.filter("Results/Miscellaneous/temp_modules_final.txt","Results/Integrated_modules_final.txt");
	}
	catch(Exception e)
	{e.printStackTrace(); return -1;}
	return 0;
	}
}
