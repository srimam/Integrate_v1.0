package bin;
import java.lang.*;
import java.*;
import java.io.*;
import java.util.regex.*;
import java.util.ArrayList;
import java.util.*;
/*
Code to calculate the TF-cluster average correlations and pvalues

java regulator_predict_corr
*/

public class regulator_predict_corr implements Runnable
{
	public static String gff2,expression2;
	public void run()
	{
		try{
		this.compute(gff2,expression2);
		System.out.println("Started corr thread");
		}
		catch(IOException e){e.printStackTrace();}
	}
	public regulator_predict_corr(String gff, String expression)
	{
		gff2=gff;expression2=expression;
	}
	public static int compute(String gff, String expression) throws IOException
	{
	try
	{
		gff2=gff;
		expression2=expression;
		//Read in interaction network
		BufferedReader b2 = new BufferedReader(new FileReader("Results/DISTILLER/DISTILLER_modules_processed_UNIQUE_seedExtended.txt"));

		//Read in list of TFs
		BufferedReader b3 = new BufferedReader(new FileReader("Data/TFs.txt"));

		//Read in expression data file
		BufferedReader b4 = new BufferedReader(new FileReader("Data/Expression/"+expression));

		//write out
		BufferedWriter buf = new BufferedWriter(new FileWriter("Results/Miscellaneous/Distance_correlation_results.txt"));
		
		ArrayList<String> start = new ArrayList<String>();
		ArrayList<String> geneID2 = new ArrayList<String>();
		ArrayList<String> geneID = new ArrayList<String>();
		ArrayList<String> motif = new ArrayList<String>();
		ArrayList<String> conditions = new ArrayList<String>();
		ArrayList<String> geneID3 = new ArrayList<String>();
		ArrayList<String> TF_geneID = new ArrayList<String>();
		ArrayList<String> chr = new ArrayList<String>();

		Pattern p5 = Pattern.compile(",");
		String lines="";

		b4.readLine();//skip header
		int rows=0,cols=0;
		while ((lines = b4.readLine())!=null)//identify all genes with available expression data
		{
			String[] result = lines.split("\t");
			if(cols==0)
				cols=result.length-1;
			geneID3.add(result[0]);
			rows++;
		}
		String[][] expression_mat = new String[rows][cols];//hold expression data [number of genes][number of exp.]
		b4.close();
		b4 = new BufferedReader(new FileReader("Data/Expression/"+expression));
		lines="";
		int rowCount=0;
		b4.readLine();
		while ((lines = b4.readLine())!=null)//place expression data into matrix
		{
			String[] result = lines.split("\t");
			for(int i=1;i<result.length;i++)
				expression_mat[rowCount][i-1]=result[i];
			rowCount++;
		}

		String[] gff_result = gff_parse.parse("Data/Gffs/"+gff);
		String[] gff_result2 = gff_result[0].split("\n");
		for(int i=1;i<gff_result2.length;i++)//skip header and start from 2nd line
		{
			String[] result = gff_result2[i].split("\t");
			geneID2.add(result[0]);
			chr.add(result[4]);
			if(result[5].equals("+"))
				start.add(result[2]);
			else
				start.add(result[3]);
		}

		lines="";
		while ((lines = b2.readLine())!=null)//DISTILLER file
		{
			String[] result = lines.split("\t");
			motif.add(result[0]);
			geneID.add(result[1]);
			conditions.add(result[2]);
		}

		lines="";
		while ((lines = b3.readLine())!=null)//TF file
		{
			String[] result = lines.split("\t");
			TF_geneID.add(result[0].trim());
		}

		//Compute TF-target correlations and distances pvalues
		String initialMotif="";
		ArrayList<String> motif_all = new ArrayList<String>();
		ArrayList<String> candidate_TF_all = new ArrayList<String>();
		ArrayList<Double> aveCor_all = new ArrayList<Double>();
		ArrayList<Double> minCor_all = new ArrayList<Double>();
		ArrayList<Double> maxCor_all = new ArrayList<Double>();
		System.out.println("Computing correlation scores... This may take a while...");
		for(int i=0;i<geneID.size();i++)
		{
			if(!initialMotif.equals(motif.get(i)))
			{
				initialMotif=motif.get(i);
				
				ArrayList<Double> correlations = new ArrayList<Double>();
				for(int j=0;j<TF_geneID.size();j++)
				{
					correlations.clear();
					for(int l=0;l<geneID.size();l++)
					{
						if(initialMotif.equals(motif.get(l)) && geneID3.contains(geneID.get(l)) && geneID3.contains(TF_geneID.get(j)) && geneID2.contains(TF_geneID.get(j)) && !TF_geneID.get(j).equals(geneID.get(l)))
						{
						String[] result = p5.split(conditions.get(l));
						double sum=0.0;
						double r=0.0;
						double product=0.0;
						double sumX=0.0,sumY=0.0;
						for(int k=0;k<result.length;k++)
						{
							int row1 = geneID3.indexOf(TF_geneID.get(j));
							int row2 = geneID3.indexOf(geneID.get(l));
							int col = Integer.parseInt(result[k].trim());
							sumX=sumX+Double.parseDouble(expression_mat[row1][col]);
							sumY=sumY+Double.parseDouble(expression_mat[row2][col]);
						}
						double meanX =sumX/result.length;
						double meanY =sumY/result.length;
						double X=0.0,Y=0.0;
						
						for(int k=0;k<result.length;k++)
						{
							int row1 = geneID3.indexOf(TF_geneID.get(j));
							int row2 = geneID3.indexOf(geneID.get(l));
							int col = Integer.parseInt(result[k].trim());
							X=X+Math.pow((Double.parseDouble(expression_mat[row1][col])-meanX),2);
							Y=Y+Math.pow((Double.parseDouble(expression_mat[row2][col])-meanY),2);
						}
						double sdX = Math.sqrt((X/(result.length-1)));
						double sdY = Math.sqrt((Y/(result.length-1)));
						for(int k=0;k<result.length;k++)
						{
							int row1 = geneID3.indexOf(TF_geneID.get(j));
							int row2 = geneID3.indexOf(geneID.get(l));
							int col = Integer.parseInt(result[k].trim());
							double x = (Double.parseDouble(expression_mat[row1][col])-meanX)/sdX;
							double y = (Double.parseDouble(expression_mat[row2][col])-meanY)/sdY;
							product=x*y;
							sum=sum+product;
						}
						r=sum/(result.length-1);
						correlations.add(r);
						}
					}
						double min=1.0,max=0.0,corrSum=0.0,counter3=0.0;
						for(int l=0;l<correlations.size();l++)
						{
							if(min>Math.abs(correlations.get(l)))
								min=Math.abs(correlations.get(l));
							if(max<Math.abs(correlations.get(l)))
								max=Math.abs(correlations.get(l));
							corrSum=corrSum+Math.abs(correlations.get(l));
							counter3++;
						}
						Collections.sort(correlations);
						if(geneID3.contains(TF_geneID.get(j)) && geneID2.contains(TF_geneID.get(j)))
						{
							if((corrSum/counter3)==Double.NaN)
								aveCor_all.add(0.0);
							else
								aveCor_all.add(corrSum/counter3);
							minCor_all.add(min);
							maxCor_all.add(max);
							motif_all.add(motif.get(i));
							candidate_TF_all.add(TF_geneID.get(j));	
						}
				}
			}
		}
		b2.close();b3.close();b4.close();
		System.out.println(motif_all.size()+"\t"+candidate_TF_all.size()+"\t"+maxCor_all.size()+"\t"+minCor_all.size()+"\t"+aveCor_all.size());

		//calculate pvalues
		ArrayList<Double> aveCor_pvalues = new ArrayList<Double>();
		
		Random rand = new Random();
		for(int j=0;j<aveCor_all.size();j++)
		{
			double counts=1.0;
			for(int q=0;q<1000;q++)
			{
				if(aveCor_all.get(rand.nextInt(aveCor_all.size()))>aveCor_all.get(j))
							counts++;
			}
			aveCor_pvalues.add(counts/1000.0);
		}

		buf.write("Motif\tCandidtate TF\tMax. corr.\tMin. corr.\tAve. corr\tpvalues\n");
		for(int j=0;j<motif_all.size();j++)
		{
			buf.write(motif_all.get(j) + "\t" + candidate_TF_all.get(j) + "\t" +maxCor_all.get(j) + "\t"+minCor_all.get(j)+"\t"+aveCor_all.get(j)+"\t"+aveCor_pvalues.get(j));
			buf.newLine();buf.flush();
		}
		buf.flush();buf.close();
	}
	catch(Exception e)
	{e.printStackTrace(); return -1;}
	return 0;
	}
}
