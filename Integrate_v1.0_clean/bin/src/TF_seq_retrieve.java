package bin;
import java.lang.*;
import java.*;
import java.io.*;
import java.util.regex.*;
import java.util.ArrayList;
import java.util.*;
/*
Code to retrieve protein sequences of TFs for Pfam analysis

java TF_seq_retrieve
*/

public class TF_seq_retrieve
{
	public static int fetch(String proteome) throws IOException
	{
	try
	{
		//List of TFs
		BufferedReader b = new BufferedReader(new FileReader("Data/TFs.txt"));

		//Proteome file
		BufferedReader b1 = new BufferedReader(new FileReader("Data/Proteomes/"+proteome));

		//gene to ID mapping
		BufferedReader b2 = new BufferedReader(new FileReader("Results/Miscellaneous/geneID_mRNA-name_all.txt"));

		//write out results
		BufferedWriter buf = new BufferedWriter(new FileWriter("Data/TF_protein_seqs.txt"));

		ArrayList<String> ID = new ArrayList<String>();
		ArrayList<String> seq = new ArrayList<String>();
		ArrayList<String> geneID = new ArrayList<String>();
		ArrayList<String> name = new ArrayList<String>();
		
		String lines = "",sequence="";
		int counter=0;
		while ((lines = b1.readLine())!=null)//read in proteome file
		{
			if(lines.length() > 1 && lines.charAt(0)=='>' && counter==0)
			{
				ID.add(lines);
				counter++;
			}
			else if(lines.length() > 1 && lines.charAt(0)=='>')
			{
				ID.add(lines);
				seq.add(sequence);
				sequence="";
			}
			else
				sequence+=lines.trim();
		}
		seq.add(sequence);

		lines="";
		b2.readLine();
		while ((lines = b2.readLine())!=null)//read in gene to ID map
		{
			String[] result=lines.split("\t");
			geneID.add(result[0]);
			name.add(result[1]);
		}

		lines="";
		while ((lines = b.readLine())!=null)//read in TFs file
		{
			for(int i=0;i<ID.size();i++)
			{
				if(geneID.indexOf(lines.trim())>-1 && ID.get(i).contains(name.get(geneID.indexOf(lines.trim()))))
				{
					buf.write(">"+lines.trim()+"\n"+seq.get(i)+"\n");
				}
			}
		}
		
		buf.flush();buf.close();b.close();b1.close();b2.close();
	}
	catch(Exception e)
	{e.printStackTrace(); return -1;}
	return 0;
	}
}
