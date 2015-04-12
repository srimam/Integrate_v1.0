package bin;
import java.lang.*;
import java.*;
import java.io.*;
import java.util.regex.*;
import java.util.ArrayList;
import java.util.*;
/*
Code to calculate a score for each TF for a given cluster based on cluster motif(s) similarity to Ecoli motif of TF with similar DBD 

usage example:java TF_Domain_family
*/

public class TF_Domain_family implements Runnable
{

	public static String answer2;
	public void run()
	{
		try
		{
			this.compute(answer2);
		}
		catch(IOException e)
		{e.printStackTrace();}
	}
	public TF_Domain_family(String answer)
	{
		answer2=answer;
	}

	public static int compute(String answer) throws IOException
	{
	try
	{
		//Read in DISTILLER Results
		BufferedReader b = new BufferedReader(new FileReader("Results/Integrated_modules_final.txt"));

		//Read in Tomtom result
		BufferedReader b1 = new BufferedReader(new FileReader("Results/DISTILLER/tomtom_out/tomtom.txt"));

		//Read in Ecoli RegDB TFs pfam results
		BufferedReader b2 = new BufferedReader(new FileReader("bin/Pfam_TF_results.txt"));

		//Read in pfam results
		BufferedReader b3 = new BufferedReader(new FileReader("Data/"+answer));
		
		//Read in TFs
		BufferedReader b4 = new BufferedReader(new FileReader("Data/TFs.txt"));

		//write out results
		BufferedWriter buf = new BufferedWriter(new FileWriter("Results/Miscellaneous/DBD_results.txt"));
		
		ArrayList<String> motifID_distiller = new ArrayList<String>();
		ArrayList<String> motifID_tomtom = new ArrayList<String>();
		ArrayList<String> Ecoli_ID_tomtom = new ArrayList<String>();
		ArrayList<String> Ecoli_ID_Pfam = new ArrayList<String>();
		ArrayList<String> pvalue = new ArrayList<String>();
		ArrayList<String> family_ID_ecoli = new ArrayList<String>();
		ArrayList<String> family_ID_target = new ArrayList<String>();
		ArrayList<String> target_ID_Pfam = new ArrayList<String>();
		ArrayList<String> target_ID_TF = new ArrayList<String>();
		
		Pattern p4 = Pattern.compile("[|]");
		
		String lines="";
		//DISTILLER file
		while ((lines = b.readLine())!=null)
		{
			String[] result = lines.split("\t");
			if(!motifID_distiller.contains(result[0]))
				motifID_distiller.add(result[0].trim());
		}

		lines="";
		//tomtom file
		b1.readLine();//discard header
		while ((lines = b1.readLine())!=null)
		{
			String[] result = lines.split("\t");
			motifID_tomtom.add(result[0].trim());
			Ecoli_ID_tomtom.add(result[1].trim());
			pvalue.add(result[5].trim());
		}

		lines="";
		//Ecoli Pfam results
		b2.readLine();
		while ((lines = b2.readLine())!=null)
		{
			String[] result = lines.split("\t");
			String[] result1 = p4.split(lines);
			family_ID_ecoli.add(result[1].trim());
			Ecoli_ID_Pfam.add(result1[0].trim());
		}
		lines="";

		//All target TFs Pfam results
		b3.readLine();//discard header
		while ((lines = b3.readLine())!=null)
		{
			if(lines.trim().length()>0 && lines.charAt(0)!='#')
			{
				lines = lines.replace("                 ","\t").replace("                ","\t").replace("               ","\t").replace("              ","\t").replace("             ","\t").replace("            ","\t").replace("           ","\t").replace("          ","\t").replace("         ","\t").replace("        ","\t").replace("       ","\t").replace("      ","\t").replace("     ","\t").replace("    ","\t").replace("   ","\t").replace("  ","\t").replace(" ","\t");
				String[] result =lines.split("\t");
				if(family_ID_ecoli.contains(result[6].trim()))
				{
					family_ID_target.add(result[6].trim());
					target_ID_Pfam.add(result[0].trim());
				}
			}
		}
		
		//Target genome TFs
		while ((lines = b4.readLine())!=null)
		{
			String[] result = lines.split("\t");
			target_ID_TF.add(result[0].trim());
		}
		System.out.println("Calculating DBD scores...");
		//obtain all scores
		ArrayList<Double> all_scores = new ArrayList<Double>();
		for(int i=0;i<motifID_distiller.size();i++)//for each cluster
		{
			Pattern p = Pattern.compile(motifID_distiller.get(i));
			for(int j=0;j<target_ID_TF.size();j++)//for each TF
			{
				double max_score=0.0;
				for(int x=0;x<motifID_tomtom.size();x++)//for each motif to DB comparison score
				{
					if(motifID_distiller.get(i).equals(motifID_tomtom.get(x)) && target_ID_Pfam.indexOf(target_ID_TF.get(j))>0 && family_ID_ecoli.get(Ecoli_ID_Pfam.indexOf(Ecoli_ID_tomtom.get(x))).equals(family_ID_target.get(target_ID_Pfam.indexOf(target_ID_TF.get(j)))))//if TF family equals comparison family
					{
						double temp = Math.log10(Double.parseDouble(pvalue.get(x)))*-1;
						if(max_score<temp) max_score=temp;
					}
				}
				all_scores.add(max_score);
			}
		}
		
		//the determine TF score for each clusted based on its DBD
		for(int i=0;i<motifID_distiller.size();i++)//for each cluster
		{
			Pattern p = Pattern.compile(motifID_distiller.get(i));
			ArrayList<Double> score = new ArrayList<Double>();
			for(int j=0;j<target_ID_TF.size();j++)//for each TF
			{
				double max_score=0.0;
				for(int x=0;x<motifID_tomtom.size();x++)//for each motif to DB comparison score
				{
					if(motifID_distiller.get(i).equals(motifID_tomtom.get(x)) && target_ID_Pfam.indexOf(target_ID_TF.get(j))>0 && family_ID_ecoli.get(Ecoli_ID_Pfam.indexOf(Ecoli_ID_tomtom.get(x))).equals(family_ID_target.get(target_ID_Pfam.indexOf(target_ID_TF.get(j)))))//if TF family equals comparison family
					{
						double temp = Math.log10(Double.parseDouble(pvalue.get(x)))*-1;
						if(max_score<temp) max_score=temp;
					}
				}
				score.add(max_score);
			}
			Random rand = new Random();
			ArrayList<Double> pvalues2 = new ArrayList<Double>();
			for(int j=0;j<score.size();j++)
			{
				if(score.get(j)>=0.52)//if qvalue not upto 0.3 don't bother...
				{
					double counts=1.0;
					for(int q=0;q<10000;q++)
					{
						//if(all_scores.get(rand.nextInt(all_scores.size()))>score.get(j))
						if(score.get(rand.nextInt(score.size()))>score.get(j))
							counts++;
					}
					pvalues2.add(counts/10000.0);
				}
				else
					pvalues2.add(1.0);
			}
			for(int j=0;j<target_ID_TF.size();j++)//for each TF
			{
				buf.write(motifID_distiller.get(i)+"\t"+target_ID_TF.get(j)+"\t"+score.get(j)+"\t"+pvalues2.get(j)+"\n");buf.flush();
			}
		}
		buf.flush();buf.close();
		b.close();b1.close();b2.close();b3.close();
		
	}
	catch(Exception e)
	{e.printStackTrace(); return -1;}
	return 0;
	}
}
