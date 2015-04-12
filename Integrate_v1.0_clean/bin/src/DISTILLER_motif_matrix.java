package bin;
import java.lang.*;
import java.*;
import java.io.*;
import java.util.regex.*;
import java.util.ArrayList;
import java.util.*;
/*
Code to generate binary motif matrix for DISTILLER

Written by Saheed Imam
Institute for Systems biology
*/

public class DISTILLER_motif_matrix
{
	public static int build(String expression) throws IOException
	{
	try
	{
		//Read in cluster information
		BufferedReader b1 = new BufferedReader(new FileReader("Results/Miscellaneous/Phylo_Clusters.txt"));
		//Read in expression matrix just for the genes
		BufferedReader b = new BufferedReader(new FileReader("Data/Expression/"+expression));
		
		//write out DISTILLER motif matrix
		BufferedWriter buf = new BufferedWriter(new FileWriter("bin/DISTILLER-V2/Data/DISTILLER_matrix.txt"));

		//write out normalized expression matrix
		BufferedWriter buf1 = new BufferedWriter(new FileWriter("bin/DISTILLER-V2/Data/Norm_expression_matrix.txt"));

		//write out index to geneID map
		BufferedWriter buf2 = new BufferedWriter(new FileWriter("Results/DISTILLER/Index_to_geneID.txt"));

		//write out index to clusterID map
		BufferedWriter buf3 = new BufferedWriter(new FileWriter("Results/DISTILLER/Index_to_clusterID.txt"));
		
		ArrayList<String> geneID = new ArrayList<String>();
		ArrayList<String> geneID2 = new ArrayList<String>();
		ArrayList<String> clusterID = new ArrayList<String>();
		ArrayList<String> clusterID2 = new ArrayList<String>();

		Matcher m;

		String lines = "";
		b.readLine();//skip header
		int rows=0,cols=0;
		while ((lines = b.readLine())!=null)//identify all genes with available expression data
		{
			String[] result = lines.split("\t");
			if(cols==0)
				cols=result.length-1;
			if(!geneID2.contains(result[0]))
				geneID2.add(result[0]);
			else
			{
				System.out.println("Duplicate expression values for "+result[0]);
				System.out.println("Exiting run as only one set of values expected per gene");
				return -1;
			}	
			rows++;
		}

		double[][] expression_matrix = new double[rows][cols];//hold expression data
		b = new BufferedReader(new FileReader("Data/Expression/"+expression));
		lines="";
		int rowCount=0;
		b.readLine();
		while ((lines = b.readLine())!=null)//place expression data into matrix
		{
			String[] result = lines.split("\t");
			for(int i=1;i<result.length;i++)
			{
				expression_matrix[rowCount][i-1]=Double.parseDouble(result[i]);
				if(rowCount==0)//only do this once (writing header for expression matrix)
					buf1.write("\t"+i);
			}
			rowCount++;
		}
		buf1.newLine();
	
		System.out.println("Normalizing expression matrix...");
		for(int i=0;i<rows;i++)
		{
			buf1.write(String.valueOf(i+1));
			double[] mean_sd = Mean_SD.compute(expression_matrix[i]);
			for(int j=0;j<cols;j++)
				buf1.write("\t"+(expression_matrix[i][j]-mean_sd[0])/mean_sd[1]);
			buf1.newLine();
		}

		lines="";
		while ((lines = b1.readLine())!=null)//obtain genes to cluster mapping
		{
			String[] result = lines.split("\t");
			geneID.add(result[1]);
			clusterID.add(result[0]);
			if(!clusterID2.contains(result[0]))
				clusterID2.add(result[0]);
		}
		
		int clusterCount=1;
		for(int i=0;i<clusterID2.size();i++)
		{
			buf.write("\t"+clusterCount);
			buf3.write(clusterCount+"\t"+clusterID2.get(i)+"\n");
			clusterCount++;
		}
		buf.newLine();
		int geneCount=1;
		for(int i=0;i<geneID2.size();i++)
		{
			buf.write(String.valueOf(geneCount));
			buf2.write(geneCount+"\t"+geneID2.get(i)+"\n");
			for(int j=0;j<clusterID2.size();j++)
			{
				int counter=0;
				for(int k=0;k<geneID.size();k++)
				{
					if(geneID2.get(i).equals(geneID.get(k)) && clusterID.get(k).equals(clusterID2.get(j)))
						counter++;
				}
				if(counter>0)
					buf.write("\t"+1);
				else
					buf.write("\t"+0);
			}
			buf.newLine();
			geneCount++;
		}
		buf.flush();b1.close();b.close();buf.close();buf1.flush();buf1.close();buf2.flush();buf2.close();buf3.flush();buf3.close();
	}
	catch(Exception e)
	{e.printStackTrace(); return -1;}
	return 0;
	}
}
