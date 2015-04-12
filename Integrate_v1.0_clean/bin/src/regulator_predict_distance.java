package bin;
import java.lang.*;
import java.*;
import java.io.*;
import java.util.regex.*;
import java.util.ArrayList;
import java.util.*;
/*
Code to calculate the TF-cluster minimum distances and pvalues

java regulator_predict_distance
*/

public class regulator_predict_distance implements Runnable
{
	public static String gff2;
	public void run()
	{
		try
		{
			this.compute(gff2);
		}
		catch(IOException e)
		{e.printStackTrace();}
	}
	public regulator_predict_distance(String gff)
	{
		gff2=gff;
	}
	public static int compute(String gff) throws IOException
	{
	try
	{
		gff2=gff;
		//Read in interaction network
		BufferedReader b2 = new BufferedReader(new FileReader("Results/DISTILLER/DISTILLER_modules_processed_UNIQUE_seedExtended.txt"));

		//Read in list of TFs
		BufferedReader b3 = new BufferedReader(new FileReader("Data/TFs.txt"));

		//write out
		BufferedWriter buf = new BufferedWriter(new FileWriter("Results/Miscellaneous/Distance_results.txt"));
		
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
		ArrayList<String> cluster_members = new ArrayList<String>();
		ArrayList<String> candidate_TF = new ArrayList<String>();
		ArrayList<Double> distance = new ArrayList<Double>();

		ArrayList<String> motif_all = new ArrayList<String>();
		ArrayList<String> candidate_TF_all = new ArrayList<String>();
		ArrayList<Double> distance_all = new ArrayList<Double>();
		ArrayList<Double> averageDim_all = new ArrayList<Double>();
		System.out.println("Computing distance scores... This may take a while...");
		for(int i=0;i<geneID.size();i++)
		{
			if(!initialMotif.equals(motif.get(i)))
			{
				candidate_TF.clear();
				cluster_members.clear();
				distance.clear();
				initialMotif=motif.get(i);
				for(int j=0;j<geneID.size();j++)
				{
					if(initialMotif.equals(motif.get(j)))
					{
						cluster_members.add(geneID.get(j));
					}
				}
				for(int j=0;j<TF_geneID.size();j++)
				{
					for(int k=0;k<cluster_members.size();k++)
					{	
						if(geneID2.contains(TF_geneID.get(j)) && geneID2.contains(cluster_members.get(k)))
						{
						if(chr.get(geneID2.indexOf(TF_geneID.get(j))).equals(chr.get(geneID2.indexOf(cluster_members.get(k)))))//if TF on the same chromosome as target
						{	
							int begin=0,end=0;
							double counter=0;
							if(Integer.parseInt(start.get(geneID2.indexOf(TF_geneID.get(j))))>Integer.parseInt(start.get(geneID2.indexOf(cluster_members.get(k)))))
							{
								begin = Integer.parseInt(start.get(geneID2.indexOf(cluster_members.get(k))));
								end = Integer.parseInt(start.get(geneID2.indexOf(TF_geneID.get(j))));
							}
							else
							{	
								begin = Integer.parseInt(start.get(geneID2.indexOf(TF_geneID.get(j))));
								end = Integer.parseInt(start.get(geneID2.indexOf(cluster_members.get(k))));
							}
							for(int l=0;l<start.size();l++)
							{
								if(Integer.parseInt(start.get(l))>begin && Integer.parseInt(start.get(l))<=end && chr.get(l).equals(chr.get(geneID2.indexOf(TF_geneID.get(j)))))
									counter++;
							}
							
							if(candidate_TF.contains(TF_geneID.get(j)))
							{
								double oldDmin = distance.get(candidate_TF.indexOf(TF_geneID.get(j)));
								double newDmin= counter;
								if(newDmin<oldDmin)
									distance.set(candidate_TF.indexOf(TF_geneID.get(j)),newDmin);
							}
							else
							{
								candidate_TF.add(TF_geneID.get(j));
								candidate_TF_all.add(TF_geneID.get(j));
								distance.add(counter);
							}
						}
						else
						{
							if(!candidate_TF.contains(TF_geneID.get(j)))
							{
								candidate_TF.add(TF_geneID.get(j));
								candidate_TF_all.add(TF_geneID.get(j));
								distance.add(Double.POSITIVE_INFINITY);
							}
						}
						}
					}
				}
				double counter2=0;
				double distanceSum=0.0;
				for(int j=0;j<distance.size();j++)
				{
					if(distance.get(j)!=Double.POSITIVE_INFINITY)
					{
						distanceSum=distanceSum+distance.get(j);
						counter2++;
					}
					distance_all.add(distance.get(j));
				}
				double averageDim=distanceSum/counter2;
				for(int j=0;j<candidate_TF.size();j++)
				{
					motif_all.add(motif.get(i));
					averageDim_all.add(averageDim);
				}
			}
		}
		b2.close();b3.close();
		System.out.println(motif_all.size()+"\t"+candidate_TF_all.size()+"\t"+averageDim_all.size()+"\t"+distance_all.size());

		//calculate pvalues
		ArrayList<Double> pvalues = new ArrayList<Double>();
		
		Random rand = new Random();
		for(int j=0;j<distance_all.size();j++)
		{
			if(distance_all.get(j)!=Double.POSITIVE_INFINITY)
			{
				double counts=1.0;
				for(int q=0;q<1000;q++)
				{
					if(distance_all.get(rand.nextInt(distance_all.size()))<distance_all.get(j))
						counts++;
				}
				pvalues.add(counts/1000.0);
			}
			else
				pvalues.add(1.0);
		}

		buf.write("Motif\tCandidtate TF\tDmin\tNorm Dim\tpvalues\n");
		for(int j=0;j<motif_all.size();j++)
		{
			buf.write(motif_all.get(j) + "\t" + candidate_TF_all.get(j) + "\t" +distance_all.get(j)+"\t"+(distance_all.get(j)/averageDim_all.get(j))+"\t"+pvalues.get(j));
			buf.newLine();buf.flush();
		}
		buf.flush();buf.close();
	}
	catch(Exception e)
	{e.printStackTrace(); return -1;}
	return 0;
	}
}
