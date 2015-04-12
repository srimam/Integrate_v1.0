package bin;
import java.lang.*;
import java.*;
import java.io.*;
import java.util.regex.*;
import java.util.ArrayList;
import java.util.*;
/*
Code to parse and compare mast files to identify overlap in their identified binding sites. This generates input for mcl clustering

usage example: java Motif_bindsite_overlap

Written by Saheed Imam
Institute for Systems biology
*/

public class Motif_bindsite_overlap
{
	public static int overlap(int number_of_motifs,String mclIn, int IGR_length) throws IOException
	{
		try
		{
			BufferedWriter buf = new BufferedWriter(new FileWriter(mclIn));
			ArrayList<Integer> motifID = new ArrayList<Integer>();
			ArrayList<String> geneID = new ArrayList<String>();
			ArrayList<Integer> start = new ArrayList<Integer>();
			ArrayList<Integer> end = new ArrayList<Integer>();
			//get all the motif hits
			int initial=0,counter=0;
			for(int i=1;i<number_of_motifs+1;i++)
			{
				String lines="";
				BufferedReader b = new BufferedReader(new FileReader("Results/Motif_finding/mast_out/Motif_"+i));
				while ((lines = b.readLine())!=null)
				{
					if(!lines.contains("#"))
					{
						String temp = lines.replace("      ","\t").replace("     ","\t").replace("    ","\t").replace("   ","\t").replace("  ","\t").replace(" ","\t");
						String[] result = temp.split("\t");
						if(Double.parseDouble(result[4])>900.0)// only consider motifs with a match score of upto 900 (hits below this are often bogus)
						{
							if(counter==0)
							{
								initial = i;
							}
							motifID.add(i);
							geneID.add(result[0].trim());
							if(result[1].trim().equals("+1"))
							{
								start.add(Integer.parseInt(result[2]));
								end.add(Integer.parseInt(result[3]));
							}
							else
							{
								start.add(IGR_length-Integer.parseInt(result[3]));
								end.add(IGR_length-Integer.parseInt(result[2]));
							}
							counter++;
						}
					}
				}
				b.close();
			}
			//compare motif hits across all motifs and generate mcl input file
			ArrayList<String> first_geneID = new ArrayList<String>();
			ArrayList<Integer> first_start = new ArrayList<Integer>();
			ArrayList<Integer> first_end = new ArrayList<Integer>();
			int tracker = initial;
			for(int i=0;i<motifID.size();i++)
			{
				if(motifID.get(i)==tracker)
				{
					first_geneID.add(geneID.get(i));
					first_start.add(start.get(i));
					first_end.add(end.get(i));
				}
				else
				{
					double overlap=0.0;
					int tracker2=initial;
					for(int j=0;j<motifID.size();j++)
					{
						if(motifID.get(j)==tracker2)
						{
							if(first_geneID.contains(geneID.get(j)))
							{
								if((start.get(j) <first_end.get(first_geneID.indexOf(geneID.get(j)))-5 && start.get(j) >= first_start.get(first_geneID.indexOf(geneID.get(j)))) || (end.get(j) >first_start.get(first_geneID.indexOf(geneID.get(j)))+5 && end.get(j) <=first_end.get(first_geneID.indexOf(geneID.get(j)))))
								{
									overlap++;
								}
							}					
						}
						else
						{
							if(overlap/(double)first_geneID.size()>=0.3)
								buf.write("Motif_"+tracker+"\tMotif_"+tracker2+"\t"+overlap/(double)first_geneID.size()+"\n");
							overlap=0.0;
							tracker2=motifID.get(j);
							if(first_geneID.contains(geneID.get(j)))
							{
								if((start.get(j) <first_end.get(first_geneID.indexOf(geneID.get(j)))-5 && start.get(j) >= first_start.get(first_geneID.indexOf(geneID.get(j)))) || (end.get(j) >first_start.get(first_geneID.indexOf(geneID.get(j)))+5 && end.get(j) <=first_end.get(first_geneID.indexOf(geneID.get(j)))))
								{
									overlap++;
								}
							}
						}
						//if at end of loop 
						if(j==motifID.size()-1)
						{
							if(overlap/(double)first_geneID.size()>=0.3)
								buf.write("Motif_"+tracker+"\tMotif_"+tracker2+"\t"+overlap/(double)first_geneID.size()+"\n");
						}
					}
					tracker=motifID.get(i);
					first_geneID.clear();
					first_start.clear();
					first_end.clear();
					first_geneID.add(geneID.get(i));
					first_start.add(start.get(i));
					first_end.add(end.get(i));
				}
				//if at end of loop 
				if(i==motifID.size()-1)
				{
					double overlap=0.0;
					int tracker2=initial;
					for(int j=0;j<motifID.size();j++)
					{
						if(motifID.get(j)==tracker2)
						{
							if(first_geneID.contains(geneID.get(j)))
							{
								if((start.get(j) <first_end.get(first_geneID.indexOf(geneID.get(j)))-5 && start.get(j) >= first_start.get(first_geneID.indexOf(geneID.get(j)))) || (end.get(j) >first_start.get(first_geneID.indexOf(geneID.get(j)))+5 && end.get(j) <=first_end.get(first_geneID.indexOf(geneID.get(j)))))
								{
									overlap++;
								}
							}					
						}
						else
						{
							if(overlap/(double)first_geneID.size()>=0.3)
								buf.write("Motif_"+tracker+"\tMotif_"+tracker2+"\t"+overlap/(double)first_geneID.size()+"\n");
							tracker2=motifID.get(j);
							overlap=0.0;
							tracker2=motifID.get(j);
							if(first_geneID.contains(geneID.get(j)))
							{
								if((start.get(j) <first_end.get(first_geneID.indexOf(geneID.get(j)))-5 && start.get(j) >= first_start.get(first_geneID.indexOf(geneID.get(j)))) || (end.get(j) >first_start.get(first_geneID.indexOf(geneID.get(j)))+5 && end.get(j) <=first_end.get(first_geneID.indexOf(geneID.get(j)))))
								{
									overlap++;
								}
							}
						}

						if(j==motifID.size()-1)
						{
							if(overlap/(double)first_geneID.size()>=0.3)
								buf.write("Motif_"+tracker+"\tMotif_"+tracker2+"\t"+overlap/(double)first_geneID.size()+"\n");
						}
					}
				}
			}
			buf.flush();buf.close();
		}
		catch(Exception e)
		{e.printStackTrace(); return -1;}
		return 0;
	}
}
