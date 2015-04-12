package bin;
import java.lang.*;
import java.*;
import java.io.*;
import java.util.regex.*;
import java.util.ArrayList;
import java.util.*;
/*
Code to parse gff file

java gff_parse
*/

public class gff_parse
{
	public static String[] parse(String gff_file) throws IOException
	{
		String temp="";
		String[] myResult = new String[2];
		temp+="GeneID\tmRNA Name\tStart\tEnd\tchr\torientation\n";
		try
		{
		//GFF file
		BufferedReader b = new BufferedReader(new FileReader(gff_file));

		ArrayList<String> ID = new ArrayList<String>();
		ArrayList<String> ID2 = new ArrayList<String>();
		ArrayList<String> start = new ArrayList<String>();
		ArrayList<String> end = new ArrayList<String>();
		ArrayList<String> chr = new ArrayList<String>();
		ArrayList<String> orientation = new ArrayList<String>();
		String lines="";
		//Parsing gff
		lines="";
		String IDstring = "";
		int count=0;
		while ((lines = b.readLine())!=null)
		{
			if(lines.contains("locus_tag="))
			{ 
				count++;
				break;
			}
		}
		
		b = new BufferedReader(new FileReader(gff_file));		
		if(count==0)//for phytozome files
		{
			while ((lines = b.readLine())!=null)
			{
				String[] result = lines.split("\t");
				if(lines.charAt(0)=='>')//if fasta format sequence included in gff
					break;
				if(lines.charAt(0)!='#' && result[2].trim().equals("gene"))
				{
					String[] result2 = result[8].split(";");
					for(int i=0;i<result2.length;i++)
					{
						String[] result3 = result2[i].split("=");
						if(result3[0].equalsIgnoreCase("Name") && !ID.contains(result2[i].replace("Name=","")))
							IDstring = result2[i].replace("Name=","");
					}
				}
				if(lines.charAt(0)!='#' && result[2].trim().equals("mRNA"))
				{
					String[] result2 = result[8].split(";");
					for(int i=0;i<result2.length;i++)
					{
						String[] result3 = result2[i].split("=");
						if(result3[0].equalsIgnoreCase("Name") && !ID2.contains(result2[i].replace("Name=","")))
						{
							ID2.add(result2[i].replace("Name=",""));
							ID.add(IDstring);
							start.add(result[3].trim());
							end.add(result[4].trim());
							chr.add(result[0].trim());
							orientation.add(result[6].trim());
						}
					}
				}
			}
		}
		else //for NCBI files
		{
			while ((lines = b.readLine())!=null)
			{
				String[] result = lines.split("\t");
				if(lines.charAt(0)=='>')//if fasta format sequence included in gff
					break;
				if(lines.charAt(0)!='#' && result[2].trim().equals("gene"))
				{
					String[] result2 = result[8].split(";");
					for(int i=0;i<result2.length;i++)
					{
						String[] result3 = result2[i].split("=");
						if(result3[0].equalsIgnoreCase("locus_tag") && !ID.contains(result2[i].replace("locus_tag=","")))
						{
							IDstring = result2[i].replace("locus_tag=","");
						}
					}
				}
				if(lines.charAt(0)!='#' && result[2].trim().equals("CDS"))
				{
					String[] result2 = result[8].split(";");
					for(int i=0;i<result2.length;i++)
					{
						String[] result3 = result2[i].split("=");
						if(result3[0].equalsIgnoreCase("Name") && !ID2.contains(result2[i].replace("Name=","")))
						{
							ID2.add(result2[i].replace("Name=",""));
							ID.add(IDstring);
							start.add(result[3].trim());
							end.add(result[4].trim());
							chr.add(result[0].trim());
							orientation.add(result[6].trim());
						}
					}
				}
			}
		}
		int counter=0;
		for(int i=0;i<ID.size();i++)
		{
			temp+=ID.get(i) +"\t"+ID2.get(i)+"\t"+start.get(i)+"\t"+end.get(i)+"\t"+chr.get(i)+"\t"+orientation.get(i)+"\n";
			counter++;
		}
		b.close();
		myResult[0] = temp;
		myResult[1] = String.valueOf(counter);
	}
	catch(Exception e)
	{e.printStackTrace(); myResult[1]="-1"; return myResult;}
	return myResult;
	}

}
