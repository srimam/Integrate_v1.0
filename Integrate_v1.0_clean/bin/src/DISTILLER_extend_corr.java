package bin;
import java.lang.*;
import java.*;
import java.io.*;
import java.util.regex.*;
import java.util.ArrayList;
import java.util.*;
/*
Code to extend DISTILLER output by identifying other genes with motifs which are highly correlated with those already in seed modules from DISTILLER analysis


Written by Saheed Imam
Institute for Systems biology
*/

public class DISTILLER_extend_corr
{
	public static int extend(String expression) throws IOException
	{
	try
	{
		//Read in expression matrix
		BufferedReader b = new BufferedReader(new FileReader("Data/Expression/"+expression));

		//Read in file with unique DISTILLER results
		BufferedReader b2 = new BufferedReader(new FileReader("Results/DISTILLER/DISTILLER_modules_processed_UNIQUE.txt"));

		//Read in cluster from footprinting analysis
		BufferedReader b3 = new BufferedReader(new FileReader("Results/Miscellaneous/Phylo_Clusters.txt"));

		//write out results
		BufferedWriter buf = new BufferedWriter(new FileWriter("Results/DISTILLER/DISTILLER_modules_processed_seedExtended.txt"));

		//write out results
		BufferedWriter buf2 = new BufferedWriter(new FileWriter("Results/DISTILLER/Cluster_genes_without_expression_data.txt"));

		//write out results
		BufferedWriter buf3 = new BufferedWriter(new FileWriter("Results/DISTILLER/Cluster_genes_with_low_Correlation.txt"));

		ArrayList<String> conditions = new ArrayList<String>();
		ArrayList<String> geneID = new ArrayList<String>();
		ArrayList<String> motif_ID = new ArrayList<String>();
		ArrayList<String> geneID2 = new ArrayList<String>();
		ArrayList<String> motif_ID2 = new ArrayList<String>();
		ArrayList<String> geneID3 = new ArrayList<String>();

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
		b = new BufferedReader(new FileReader("Data/Expression/"+expression));
		lines="";
		int rowCount=0;
		b.readLine();
		while ((lines = b.readLine())!=null)//place expression data into matrix
		{
			String[] result = lines.split("\t");
			for(int i=1;i<result.length;i++)
				expression_mat[rowCount][i-1]=Double.parseDouble(result[i]);
			rowCount++;
		}

		lines="";
		while ((lines = b2.readLine())!=null)//DISTILLER results
		{
			String[] result = lines.split("\t");
			motif_ID.add(result[0]);
			geneID.add(result[1]);
			conditions.add(result[2]);
		}

		lines="";
		while ((lines = b3.readLine())!=null)//phylo results
		{
			String[] result = lines.split("\t");
			motif_ID2.add(result[0]);
			geneID2.add(result[1]);
		}

		String initialMotif="";
		ArrayList<String> temp = new ArrayList<String>();
		for(int i=0;i<geneID2.size();i++)//loop through phylo clusters
		{
			int counter=0;
			if(!initialMotif.equals(motif_ID2.get(i)))
			{
				temp.clear();
				initialMotif=motif_ID2.get(i);
				for(int j=0;j<geneID.size();j++)
				{
					if(motif_ID.get(j).equals(motif_ID2.get(i)))
					{
						temp.add(geneID.get(j));
						buf.write(motif_ID.get(j) + "\t" + geneID.get(j) + "\t" + conditions.get(j));buf.newLine();buf.flush();
					}
				}
			}
			if(!temp.contains(geneID2.get(i)))//if genes not in DISTILLER seed module
			{
				if(geneID3.contains(geneID2.get(i)))//if expression data exists for gene within a given cluster
				{
				for(int j=0;j<geneID.size();j++)//loop through all DISTILLER seed module genes
				{
					if(motif_ID.get(j).equals(motif_ID2.get(i)))//if genes belong to same cluster
					{
						String[] result = conditions.get(j).split(",");
						double sum=0.0;
						double r=0.0;
						double product=0.0;
						double sumX=0.0,sumY=0.0;
						for(int k=0;k<result.length;k++)
						{
							int row1 = geneID3.indexOf(geneID.get(j));
							int row2 = geneID3.indexOf(geneID2.get(i));
							int col = Integer.parseInt(result[k].trim());
							sumX+=expression_mat[row1][col];
							sumY+=expression_mat[row2][col];
						}
						double meanX =sumX/result.length;
						double meanY =sumY/result.length;
						double X=0.0,Y=0.0;
						
						for(int k=0;k<result.length;k++)
						{
							int row1 = geneID3.indexOf(geneID.get(j));
							int row2 = geneID3.indexOf(geneID2.get(i));
							int col = Integer.parseInt(result[k].trim());
							X+=Math.pow((expression_mat[row1][col]-meanX),2);
							Y+=Math.pow((expression_mat[row2][col]-meanY),2);
						}
						double sdX = Math.sqrt((X/(result.length-1)));
						double sdY = Math.sqrt((Y/(result.length-1)));
						for(int k=0;k<result.length;k++)
						{
							int row1 = geneID3.indexOf(geneID.get(j));
							int row2 = geneID3.indexOf(geneID2.get(i));
							int col = Integer.parseInt(result[k].trim());
							double x = (expression_mat[row1][col]-meanX)/sdX;
							double y = (expression_mat[row2][col]-meanY)/sdY;
							product=x*y;
							sum+=product;
						}
						r=sum/(result.length-1);
						if(Math.abs(r)>=0.5)
						{
							buf.write(motif_ID.get(j) + "\t" + geneID2.get(i) + "\t" + conditions.get(j)+ "\t" +r);buf.newLine();buf.flush();
							counter++;
						}
					}
				}
					if(counter==0)//if gene not correlated to any in seed module
						buf3.write(motif_ID2.get(i) + "\t" + geneID2.get(i)+"\n");buf3.flush();
				}
				else//phylo identified but no expression data
					buf2.write(motif_ID2.get(i) + "\t"+ geneID2.get(i)+"\n");buf2.flush();
			}
		}
		
		buf.flush();buf2.flush();buf3.flush();
		buf.close();buf2.close();buf3.close();
		b.close();b2.close();b3.close();
		
		Unique_interactions.filter("Results/DISTILLER/DISTILLER_modules_processed_seedExtended.txt","Results/DISTILLER/DISTILLER_modules_processed_UNIQUE_seedExtended.txt");
	}
	catch(Exception e)
	{e.printStackTrace(); return -1;}
	return 0;
	}
}
