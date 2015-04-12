package bin;
import java.lang.*;
import java.*;
import java.io.*;
import java.util.regex.*;
import java.util.ArrayList;
import java.util.*;
/*
Code to be used to extract intergenomic regions from genome file and gff

java Extract_upstream_regions Creinhardtii_236.fa Creinhardtii_236_gene.gff3 500 result.txt

where 500 is the IGR size 
*/

public class Extract_upstream_regions_revcomp
{
	public static void extract(String genome_file, String gff_file, int IGR_length,String outfile) throws IOException
	{
		try
		{
		//File containing fasta format of genome
		BufferedReader b = new BufferedReader(new FileReader(genome_file));
		//GFF file
		BufferedReader b1 = new BufferedReader(new FileReader(gff_file));

		BufferedWriter buf = new BufferedWriter(new FileWriter(outfile));
		
		ArrayList<String> chr = new ArrayList<String>();
		ArrayList<String> seq = new ArrayList<String>();
		StringBuilder allSeq= new StringBuilder();

		ArrayList<String> chrName = new ArrayList<String>();
		ArrayList<Integer> start = new ArrayList<Integer>();
		ArrayList<Integer> end = new ArrayList<Integer>();
		ArrayList<String> orientation = new ArrayList<String>();
		ArrayList<String> ID = new ArrayList<String>();
		String lines="";
		//Parsing genome file...
		//System.out.println("Extracting IGRs for "+genome_file);
		while ((lines = b.readLine())!=null)
		{	
			if(lines.length() > 1 && lines.charAt(0)=='>')
			{
				String[] result = lines.split(" ");
				String[] result1 = result[0].split("[|]");
				if(result1.length>1)
					chr.add(result1[3].trim());
				else
					chr.add(result[0].trim().replace(">",""));
				//System.out.println(result[0].trim().replace(">","")+" done...");
				if(allSeq.length()!=0)
				{
					seq.add(allSeq.toString());
					allSeq.delete(0,allSeq.length());
				}
			}
			else
			{
				allSeq.append(lines);
			}
		}
		seq.add(allSeq.toString());
		//Parsing gff
		lines="";
		int counter=1;
		while ((lines = b1.readLine())!=null)
		{
			String[] result = lines.split("\t");
			if(lines.charAt(0)=='>')//if fasta format sequence included in gff
				break;
			if(lines.charAt(0)!='#' && result[2].trim().equals("gene"))
			{
				//sort entries by gene start position
				int track=0,j=0;
				for(int i=0;i<chrName.size();i++)
				{	j++;
					if(chrName.get(i).equals(result[0]) && Integer.parseInt(result[3])<start.get(i))
					{
						chrName.add(i,result[0]);
						start.add(i,Integer.parseInt(result[3]));
						end.add(i,Integer.parseInt(result[4]));
						orientation.add(i,result[6]);
						track++;
						i=chrName.size();
					}
				}
				if(track==0)
				{
					chrName.add(result[0]);
					start.add(Integer.parseInt(result[3]));
					end.add(Integer.parseInt(result[4]));
					orientation.add(result[6]);
				}
				String[] result2 = result[8].split(";");
				int count=0;
				for(int i=0;i<result2.length;i++)
				{
					if(result2[i].contains("locus_tag="))//use locus tag
					{
						if(!ID.contains(result2[i].replace("locus_tag=","")))
						{
							if(track==0)
								ID.add(result2[i].replace("locus_tag=",""));
							else
								ID.add(j-1,result2[i].replace("locus_tag=",""));
						}
						else
						{
							if(track==0)
								ID.add(result2[i].replace("locus_tag=","")+"_"+counter);
							else
								ID.add(j-1,result2[i].replace("locus_tag=","")+"_"+counter);
							counter++;
						}
						count++;
					}
				}
				if(count==0)//if locus tag not found use ID
				{
					for(int i=0;i<result2.length;i++)
					{
						if(result2[i].contains("Name="))
						{
							if(!ID.contains(result2[i].replace("Name=","")))
							{
								if(track==0)
									ID.add(result2[i].replace("Name=",""));
								else
									ID.add(j-1,result2[i].replace("Name=",""));
							}
							else
							{
								if(track==0)
									ID.add(result2[i].replace("Name=","")+"_"+counter);
								else
									ID.add(j-1,result2[i].replace("Name=","")+"_"+counter);
								counter++;
							}
						}
					}
				}
			}
		}
		String trial="";
		for(int i=0;i<chrName.size();i++)
		{
			String temp="";
			try{
			if(orientation.get(i).equals("+"))
			{//System.out.println(ID.get(i));
				if((start.get(i)-IGR_length-2)<0 && !trial.equals(chrName.get(i)))//if end of scaffold and is first gene
				{
					temp=seq.get(chr.indexOf(chrName.get(i))).substring(0,start.get(i)-1);
					trial=chrName.get(i);
				}
				else if(i>0 && chrName.get(i).equals(chrName.get(i-1)) && ((start.get(i)-IGR_length)-1)<(end.get(i-1)-1) && end.get(i-1)<start.get(i)-1)//prevent overlap with next gene
					temp=seq.get(chr.indexOf(chrName.get(i))).substring((end.get(i-1)),start.get(i)-1);
				else if(i>0 && chrName.get(i).equals(chrName.get(i-1)) && ((start.get(i)-IGR_length)-1)<(end.get(i-1)-1) && end.get(i-1)>=start.get(i)-1)//if structural genes overlap
				{/*Do nothing*/}
				else
					temp=seq.get(chr.indexOf(chrName.get(i))).substring((start.get(i)-IGR_length-1),start.get(i)-1);
			}
			else
			{
				if((end.get(i)+IGR_length+1)>seq.get(chr.indexOf(chrName.get(i))).length() && i+1<chrName.size() && !trial.equals(chrName.get(i+1)))//if end of scaffold and last gene
					temp=seq.get(chr.indexOf(chrName.get(i))).substring(end.get(i),seq.get(chr.indexOf(chrName.get(i))).length());
				else if((end.get(i)+IGR_length+1)>seq.get(chr.indexOf(chrName.get(i))).length() && i==chrName.size()-1)
					temp=seq.get(chr.indexOf(chrName.get(i))).substring(end.get(i),seq.get(chr.indexOf(chrName.get(i))).length());
				else if(i<chrName.size()-1 && chrName.get(i).equals(chrName.get(i+1)) && (end.get(i)+IGR_length+1)>(start.get(i+1)-1) && end.get(i)<(start.get(i+1)-1)) //prevent overlap with next gene 
					temp=seq.get(chr.indexOf(chrName.get(i))).substring(end.get(i),(start.get(i+1)-1));
				else if(i<chrName.size()-1 && chrName.get(i).equals(chrName.get(i+1)) && (end.get(i)+IGR_length+1)>(start.get(i+1)-1) && end.get(i)>=(start.get(i+1)-1))//if structural genes overlap
				{/*Do nothing*/}
				else
					temp=seq.get(chr.indexOf(chrName.get(i))).substring(end.get(i),(end.get(i)+IGR_length));
				char[] test = temp.toCharArray();
				String temp2="";				
				for(int x=test.length-1;x>=0;x--)
				{
					if(test[x]==('A')) temp2+="T";
					else if(test[x]==('T')) temp2+="A";
					else if(test[x]==('G')) temp2+="C";
					else if(test[x]==('C')) temp2+="G";
					else if(test[x]==('N')) temp2+="N";
					else
					{ 
						//System.out.println("Unknown character encountered in reverse complement... replacing with N");
						temp2+="N";
					}
				}
				temp=temp2;
			}
			if(temp.length()>=40)
				buf.write(">"+ID.get(i)+"\n"+temp+"\n");
			}catch(Exception e){System.out.println(start.get(i)+"\t"+end.get(i)+"\t"+ID.get(i));}
		}
		b.close();b1.close();buf.flush();buf.close();
	}
	catch(Exception e)
	{e.printStackTrace(); }
	}
}
