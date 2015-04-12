package bin;
import java.lang.*;
import java.*;
import java.io.*;
import java.util.regex.*;
import java.util.ArrayList;
import java.util.*;
/*
Code to retrieve protein sequences of TFs for Pfam analysis

java regulator_predict_merge
*/

public class regulator_predict_merge
{
	public static int merge() throws IOException
	{
	try
	{
		//Read in DBD results
		BufferedReader b = new BufferedReader(new FileReader("Results/Miscellaneous/DBD_results.txt"));

		//Read in correlation/distance results
		BufferedReader b1 = new BufferedReader(new FileReader("Results/Miscellaneous/Distance_correlation_results.txt"));

		//Read in correlation/distance results
		BufferedReader b2 = new BufferedReader(new FileReader("Results/Miscellaneous/Distance_results.txt"));

		//write out results
		BufferedWriter buf = new BufferedWriter(new FileWriter("Results/Miscellaneous/Top10_regulators.txt"));

		ArrayList<String> DBD_gene_clusterID = new ArrayList<String>();
		ArrayList<Double> DBD_score = new ArrayList<Double>();
		ArrayList<String> corr_gene_clusterID = new ArrayList<String>();
		ArrayList<Double> corr_score = new ArrayList<Double>();
		ArrayList<String> Dist_clusterID = new ArrayList<String>();
		ArrayList<Double> dist_score = new ArrayList<Double>();
		
		System.out.println("Calculating R_score...");
		String lines = "";
		while ((lines = b.readLine())!=null)//DBD data
		{
			String[] result = lines.split("\t");
			DBD_gene_clusterID.add(result[0].trim()+"\t"+result[1].trim());
			DBD_score.add(Math.log10(Double.parseDouble(result[3].trim()))*-1);
		}

		lines="";
		b1.readLine();
		while ((lines = b1.readLine())!=null)//correlation data
		{
			String[] result=lines.split("\t");
			corr_gene_clusterID.add(result[0].trim()+"\t"+result[1].trim());
			corr_score.add(Math.log10(Double.parseDouble(result[5].trim()))*-1);
		}

		lines="";
		b2.readLine();
		while ((lines = b2.readLine())!=null)//distance data
		{
			String[] result=lines.split("\t");
			Dist_clusterID.add(result[0].trim()+"\t"+result[1].trim());
			dist_score.add(Math.log10(Double.parseDouble(result[4].trim()))*-1);
		}

		ArrayList<Double> R_score = new ArrayList<Double>();
		ArrayList<String> ID = new ArrayList<String>();
		String init="";
		for(int i=0;i<DBD_gene_clusterID.size();i++)//Sum up R score
		{
			String[] result = DBD_gene_clusterID.get(i).split("\t");
			if(init.equals(""))
				init = result[0];
			if(!init.equals(result[0]))
			{
				ArrayList<Double> track = new ArrayList<Double>(R_score);
				Collections.sort(R_score);
				for(int j=track.size()-1;j>track.size()-6;j--)
					buf.write(ID.get(track.indexOf(R_score.get(j)))+"\t"+R_score.get(j)+"\n");
				R_score = new ArrayList<Double>();
				ID = new ArrayList<String>();
				init = result[0];
			}
			ID.add(DBD_gene_clusterID.get(i));
			if(corr_gene_clusterID.contains(DBD_gene_clusterID.get(i)) && Dist_clusterID.contains(DBD_gene_clusterID.get(i)))
				R_score.add(DBD_score.get(i)+corr_score.get(corr_gene_clusterID.indexOf(DBD_gene_clusterID.get(i)))+dist_score.get(Dist_clusterID.indexOf(DBD_gene_clusterID.get(i))));
			else if(corr_gene_clusterID.contains(DBD_gene_clusterID.get(i)) && !Dist_clusterID.contains(DBD_gene_clusterID.get(i)))
				R_score.add(DBD_score.get(i)+corr_score.get(corr_gene_clusterID.indexOf(DBD_gene_clusterID.get(i))));
			else if(Dist_clusterID.contains(DBD_gene_clusterID.get(i)) && !corr_gene_clusterID.contains(DBD_gene_clusterID.get(i)))
				R_score.add(DBD_score.get(i)+dist_score.get(Dist_clusterID.indexOf(DBD_gene_clusterID.get(i))));
			else
				R_score.add(DBD_score.get(i));
			if(i==DBD_gene_clusterID.size()-1)
			{
				ArrayList<Double> track = new ArrayList<Double>(R_score);
				Collections.sort(R_score);
				for(int j=track.size()-1;j>track.size()-6;j--)
					buf.write(ID.get(track.indexOf(R_score.get(j)))+"\t"+R_score.get(j)+"\n");
			}
		}
		buf.flush();buf.close();b.close();b1.close();
	}
	catch(Exception e)
	{e.printStackTrace(); return -1;}
	return 0;
	}
}
