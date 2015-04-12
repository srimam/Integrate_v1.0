package bin;
import java.lang.*;
import java.*;
import java.io.*;
import java.util.regex.*;
import java.util.ArrayList;
import java.util.Collections;
import java.math.*;

/*
Compute the enrichment of phylo clusters from GO terms

usage: java GO_enrichment
*/

public class GO_enrichment
{
	public static void significance(String clusters, String go_terms, int genomeSize, String results) throws IOException
	{
	try
	{
		Matcher m1;
		//First file should contain your list of clusters formatted appropriately
		BufferedReader b = new BufferedReader(new FileReader(clusters));
		
		//Second file should contain your list of gene sets/go terms formatted appropriately
		BufferedReader b2 = new BufferedReader(new FileReader(go_terms));

		//Name of file to write results
		BufferedWriter buf = new BufferedWriter(new FileWriter(results));	

		String firstLine = "ID"+"\t"+"Gene set Name"+"\t"+"Gene set size"+"\t"+"DE genes"+"\t"+"Pvalue (hyper)";
		buf.write(firstLine);
		buf.newLine();

		String genes,DEg;
		String GID = "";
		int N = genomeSize; //Total number of annotated genes in organism of interest
		int M,K,X;
		int test = 0;
		ArrayList<String> geneSet = new ArrayList<String>();
		ArrayList<String> DEgenes = new ArrayList<String>();
		ArrayList<String> regulators = new ArrayList<String>();
		ArrayList<String> uniqueRegulators = new ArrayList<String>();
		ArrayList<String> genesInCluster = new ArrayList<String>();
		ArrayList<String> goTerm = new ArrayList<String>();
		ArrayList<String> unique_goTerm = new ArrayList<String>();
		ArrayList<String> goTermMember = new ArrayList<String>();
		ArrayList<String> temp = new ArrayList<String>();

		//Save the list of all genes in cluster in DEgenes arraylist
		
		while ((DEg = b.readLine())!=null)
		{	
			String[] result = DEg.split("\t");
			if(DEg.contains(":"))//if phylo clusters
			{
				for(int i=1;i<result.length;i++)
				{
					DEgenes.add(result[i]);//gene IDs
					regulators.add(result[0].replace(":","").trim());//cluster IDs
				}
				if(!uniqueRegulators.contains(result[0])) uniqueRegulators.add(result[0].replace(":","").trim());//unique clusters
			}
			else//if final cluster
			{
				DEgenes.add(result[1]);//gene IDs
				regulators.add(result[0].trim());//cluster IDs
				if(!uniqueRegulators.contains(result[0])) uniqueRegulators.add(result[0].trim());//unique clusters
			}
		}
		String regulator = regulators.get(0);
		b.close();
		while ((genes = b2.readLine())!=null)//parse GOA file with go terms
		{
			if(genes.charAt(0)!='!')
			{
				String[] result = genes.split("\t");
				if(result.length>3)
				{
					if(!unique_goTerm.contains(result[4].trim()))
						unique_goTerm.add(result[4].trim());
					String[] result1 = result[10].split("[|]");
					if(!temp.contains(result[4].trim()+"\t"+result1[result1.length-1].trim()))
					{
						goTerm.add(result[4].trim());
						goTermMember.add(result1[result1.length-1].trim());
						temp.add(result[4].trim()+"\t"+result1[result1.length-1].trim());
					}
				}
				else
				{
					if(!unique_goTerm.contains(result[1].trim()))
						unique_goTerm.add(result[1].trim());
					if(!temp.contains(result[0].trim()+"\t"+result[1].trim()))
					{
						goTerm.add(result[1].trim());
						goTermMember.add(result[0].trim());
						temp.add(result[0].trim()+"\t"+result[1].trim());
					}
				}
			}
		}
		int count=0;
		for(int x = uniqueRegulators.size()-1; x >=0 ;x--)//for each cluster
		{
			//System.out.println(uniqueRegulators.get(x));
			int counter=0;

			for(int y = 0; y < regulators.size(); y++)
			{
				if(regulators.get(y).equals(uniqueRegulators.get(x)))
					genesInCluster.add(DEgenes.get(y));//pool together all genes for that cluster
			}
			K = genesInCluster.size();
			for(int j=0;j<unique_goTerm.size();j++)//for each go term
			{	
				for(int k=0;k<goTerm.size();k++)
				{
					if(goTerm.get(k).equals(unique_goTerm.get(j)))
						geneSet.add(goTermMember.get(k));//pool together all genes per go term
				}
				X = 0;
				M = geneSet.size();
				for(int i=0;i<geneSet.size();i++)
				{
					if(genesInCluster.contains(geneSet.get(i)))//count number of hits per go term
					{
						X++;
					}
				}
				BigDecimal pvalue =GO_enrichment.phyper(N,M,K,X);//determine significance
				if(pvalue.doubleValue() <= 0.001 && X > 1 && M > 2)//only consider situations with at least an overlap of 2 and p val of 0.001
				{
					if(pvalue.doubleValue()<10e-25)
						buf.write(unique_goTerm.get(j)+"\t"+GO_enrichment.getGeneSetName(unique_goTerm.get(j))+"\t"+M+"\t"+X+"\t"+0+"\t"+uniqueRegulators.get(x));//print out results
					else
						buf.write(unique_goTerm.get(j)+"\t"+GO_enrichment.getGeneSetName(unique_goTerm.get(j))+"\t"+M+"\t"+X+"\t"+pvalue+"\t"+uniqueRegulators.get(x));//print out results
					buf.newLine();buf.flush();
					counter++;
				}
				geneSet.clear();
			}
			if(counter>0)
			{
				buf.newLine();
				count++;
			}
			buf.flush();b2.close();
			genesInCluster.clear();
		}
		System.out.println("Number clusters with significant GO terms : "+count);
	}
	catch(Exception e)
	{e.printStackTrace();}
	}
	public static String getGeneSetName(String geneSetID)throws Exception
	{
		Pattern p5 = Pattern.compile("GO:");
		Matcher m2 = p5.matcher(geneSetID);
		String geneSetDescription = "";	
		if(m2.find())
		{
			BufferedReader b = new BufferedReader(new FileReader("bin/GO_descriptions.txt"));
			String line="";
			while((line=b.readLine())!=null)
			{
				String[] results = line.split("\t");
				if(geneSetID.equals(results[0].trim().replace("\"","")))
				{
					geneSetDescription = results[1];
				}
			}
			b.close();
		}
		//else System.out.println("Could not find required gene set description file");
		return geneSetDescription;	
	}
	private static BigDecimal phyper(int N, int m, int k, int x)
	{
		BigDecimal prob= new BigDecimal(0.0);
		int n = N-m;
		BigDecimal nmfact= new BigDecimal("1");
		BigDecimal kfact= new BigDecimal("1");
		BigDecimal mfact= new BigDecimal("1");
		BigDecimal xfact= new BigDecimal("1");
		BigDecimal nfact= new BigDecimal("1");
		BigDecimal k_xfact= new BigDecimal("1");
		BigDecimal choose_mx= new BigDecimal("0");
		BigDecimal choose_nk_x= new BigDecimal("0");
		BigDecimal choose_m_n_k= new BigDecimal("0");
		BigDecimal test1= new BigDecimal("0");
		
		for(int i=n+m;i>(m+n-k);i--)
		{
			nmfact=nmfact.multiply(new BigDecimal(String.valueOf(i)));
		}
		for(int i=k;i>1;i--)
		{
			kfact=kfact.multiply(new BigDecimal(String.valueOf(i)));
		}
		choose_m_n_k=nmfact.divide(kfact);
		for(int j=x;j>=0;j--)
		{
			for(int i=m;i>(m-j);i--)
			{
				mfact=mfact.multiply(new BigDecimal(String.valueOf(i)));
			}
			for(int i=j;i>1;i--)
			{
				xfact=xfact.multiply(new BigDecimal(String.valueOf(i)));
			}
			choose_mx =mfact.divide(xfact);
			for(int i=n;i>(n-(k-j));i--)
			{
				nfact=nfact.multiply(new BigDecimal(String.valueOf(i)));
			}
			for(int i=k-j;i>1;i--)
			{
				k_xfact=k_xfact.multiply(new BigDecimal(String.valueOf(i)));
			}
			choose_nk_x=nfact.divide(k_xfact);
			prob=prob.add(choose_mx.multiply(choose_nk_x.divide(choose_m_n_k,30, RoundingMode.HALF_UP)));
			mfact= new BigDecimal("1");
			xfact= new BigDecimal("1");
			nfact= new BigDecimal("1");
			k_xfact= new BigDecimal("1");
		}
		
		return new BigDecimal(1.0).subtract(prob);
	}
}
