import java.lang.*;
import java.*;
import java.io.*;
import java.util.regex.*;
import java.util.ArrayList;
import java.util.*;
import bin.*;
import java.nio.file.StandardCopyOption.*;
import java.util.concurrent.*;
/*
Controller class that calls all the other classes when needed

usage: java Integrate

BY

Saheed Imam
Institute for Systems Biology
*/

public class Integrate
{

	public static void main(String[] args) throws IOException
	{
	try
	{
		System.out.println("Run started at :"+Calendar.getInstance().getTime());
		//File containing fasta format of genome
		BufferedReader b = new BufferedReader(new FileReader("Parameters.txt"));

		BufferedWriter buf = new BufferedWriter(new FileWriter("bin/commands"));
		
		ArrayList<String> genomes = new ArrayList<String>();
		ArrayList<String> proteomes = new ArrayList<String>();
		ArrayList<String> gff = new ArrayList<String>();
		ArrayList<String> abbr = new ArrayList<String>();
		CreateFasta_for_MEME ext;
		String lines="",username="",password="",target="",goa="",expression="",operon="",pvalue="";
		int IGR_length=500,order=2,format=4;
		buf.write("#!/bin/bash\n");
		buf.write("# Meme bash script\n");

///////////////////Parsing Parameters.txt
		while ((lines = b.readLine())!=null)
		{	
			if(!lines.contains("#"))
			{
				if(lines.contains("Genomic_sequences"))
				{
					String[] result = lines.split("[=]");
					String[] result1 = result[1].split(",");
					for(int i=0;i<result1.length;i++)
					{
						genomes.add(result1[i].trim());
					}
				}
				if(lines.contains("Proteomes"))
				{
					String[] result = lines.split("[=]");
					String[] result1 = result[1].split(",");
					for(int i=0;i<result1.length;i++)
					{
						proteomes.add(result1[i].trim());
					}
				}
				if(lines.contains("GFF"))
				{
					String[] result = lines.split("[=]");
					String[] result1 = result[1].split(",");
					for(int i=0;i<result1.length;i++)
					{
						gff.add(result1[i].trim());
					}
				}
				if(lines.contains("Abbr"))
				{
					String[] result = lines.split("[=]");
					String[] result1 = result[1].split(",");
					for(int i=0;i<result1.length;i++)
					{
						abbr.add(result1[i].trim());
					}
				}
				if(lines.contains("IGR_length"))
				{
					String[] result = lines.split("[=]");
					try
					{
						IGR_length = Integer.parseInt(result[1].trim());
					}
					catch(Exception e)
					{
						System.out.println("No number provided for IGR length... Using default of 500...\n\n");
						IGR_length = 500;
					}
				}
				if(lines.contains("Order"))
				{
					String[] result = lines.split("[=]");
					try
					{
						order = Integer.parseInt(result[1].trim());
					}
					catch(Exception e)
					{
						System.out.println("No number provided for bg distribution order... Using default of 2...\n\n");
						order = 2;
					}
				}
				if(lines.contains("Username"))
				{
					String[] result = lines.trim().split("[=]");
					if(result.length<2)
					{
						System.out.println("Mysql login credentials required");
						return;
					}
					else
						username = result[1].trim(); 	
				}
				if(lines.contains("Password"))
				{
					String[] result = lines.trim().split("[=]");
					if(result.length<2)
					{
						System.out.println("Mysql login credentials required");
						return;
					}
					else
						password = result[1].trim(); 	
				}
				if(lines.contains("Target"))
				{
					String[] result = lines.trim().split("[=]");
					if(result.length<2)
					{
						target = abbr.get(0);
						System.out.println("WARNING: Target organism not specified. Using first organism from Abbr\n\n");
					}
					else
						target = result[1].trim(); 	
				}
				if(lines.contains("GOA"))
				{
					String[] result = lines.trim().split("[=]");
					if(result.length<2)
					{
						System.out.println("No GO annotation file provided. Gene set enrichment will not be computed\n\n");
					}
					else
						goa = result[1].trim(); 	
				}
				if(lines.contains("Expression_data"))
				{
					String[] result = lines.trim().split("[=]");
					if(result.length<2)
					{
						System.out.println("No expression data file provided. Analysis will stopped at the end of comparative genomes analysis\n\n");
					}
					else
						expression = result[1].trim(); 	
				}
				if(lines.contains("Operon_data"))
				{
					String[] result = lines.trim().split("[=]");
					if(result.length<2)
					{
						System.out.println("No operon data file provided. Operon extension will not be performed...\n\n");
					}
					else
						operon = result[1].trim(); 	
				}
				if(lines.contains("pvalue"))
				{
					String[] result = lines.trim().split("[=]");
					if(result.length<2)
					{
						pvalue="0.000001";
					}
					else
						pvalue = result[1].trim(); 	
				}
			}	
		}
///////////////Done parsing Parameters.txt

		if(gff.size()!=proteomes.size()|| genomes.size()!=gff.size()||abbr.size()!=gff.size())
		{
			System.out.println("Missing gff, genome, proteome or abbreviations... Check parameter file carefully!");
			return;
		}
		Runtime runtime;
		Process process;
		BufferedReader reader;
		File fasta;
		Scanner s;
		int test2=0;
		int size=0;

		//Prep for new run
		initialize.init(username,password);

		//1. IGR_extract
		System.out.println("Extracting IGRs...");
		for(int i=0; i<genomes.size();i++)//write commands to extract IGRs
		{
			try
			{
				Extract_upstream_regions_revcomp.extract("Data/Genomes/"+genomes.get(i),"Data/Gffs/"+gff.get(i),IGR_length,"Data/IGRs/"+abbr.get(i)+"_IGRs.txt");
			}
			catch(Exception e){System.out.println("You need to provide one genome, gff and proteome file for each organism"); return;}
		}
		buf.write("cat ");//Concantenate files
		for(int i=0; i<abbr.size();i++)
		{
			buf.write("Data/IGRs/"+abbr.get(i)+"_IGRs.txt ");
		}

		buf.write("> Data/IGRs/Combined_IGRs.txt\n");
		System.out.println("Calculating bg distribution......");
		buf.write("fasta-get-markov -m "+order +" <Data/IGRs/Combined_IGRs.txt >Data/bg/Background_dist_"+order+"_order.txt\n");//Create bg distribution
		buf.write("fasta-get-markov -m 2 <Data/IGRs/"+target+"_IGRs.txt"+" >Data/bg/Background_dist_target_2_order.txt\n");//Create target-specific bg distribution
		

		//2. orthomcl
		System.out.println("Running orthomcl...... (This may take a while depending on the number and sizes of your genomes...)");
		for(int i=0; i<proteomes.size();i++)
		{
			File newFile = new File("Data/Proteomes/"+proteomes.get(i));
			if(!newFile.exists())
			{
				System.out.println("Missining proteome file(s)");
				return;	
			}
			reader = new BufferedReader(new FileReader("Data/Proteomes/"+proteomes.get(i)));
			if(reader.readLine().split("[|]").length>=4) format = 4;
			else format = 1;
			reader.close();
			buf.write("perl bin/orthomclSoftware-v2.0.9/bin/orthomclAdjustFasta "+abbr.get(i)+" Data/Proteomes/"+proteomes.get(i)+" "+format+"\n");
			buf.write("mv "+abbr.get(i)+".fasta Data/Proteomes/compliantFasta/"+abbr.get(i)+".fasta\n");
		}		
		buf.write("perl bin/orthomclSoftware-v2.0.9/bin/orthomclFilterFasta Data/Proteomes/compliantFasta 10 20\n");
		buf.write("mv goodProteins.fasta Data/Proteomes/goodProteins.fasta\n");
		buf.write("mv poorProteins.fasta Data/Proteomes/poorProteins.fasta\n");
		buf.write("formatdb -i Data/Proteomes/goodProteins.fasta\n");
		buf.write("mv formatdb.log Data/Proteomes/formatdb.log\n");
		buf.write("blastall -i Data/Proteomes/goodProteins.fasta -d Data/Proteomes/goodProteins.fasta -p blastp -o Data/Proteomes/all_v_all_proteomes -v 100000 -b 100000 -F  'm S' -m 8 -e 1e-5\n");
		buf.write("perl bin/orthomclSoftware-v2.0.9/bin/orthomclBlastParser Data/Proteomes/all_v_all_proteomes Data/Proteomes/compliantFasta > Data/Proteomes/SimilarSequences.txt\n");
		buf.write("echo \"drop database orthomcl\" | mysql -u"+username+" -p"+password+"\n");
		buf.write("echo \"create database orthomcl\" | mysql -u"+username+" -p"+password+"\n");
		buf.write("mysql -u"+username+" -p"+password+" orthomcl < bin/orthomclSoftware-v2.0.9/bin/orthomclInstallSchema.sql\n");
		buf.write("mysqlimport --verbose -u"+username+" -p"+password+" orthomcl Data/Proteomes/SimilarSequences.txt --local\n");
		buf.write("perl bin/orthomclSoftware-v2.0.9/bin/orthomclPairs bin/orthomcl.config bin/orthomclSoftware-v2.0.9/orthomcl_pairs.log cleanup=yes\n");
		buf.write("perl bin/orthomclSoftware-v2.0.9/bin/orthomclDumpPairsFiles bin/orthomcl.config\n");
		buf.write("mv mclInput Results/Miscellaneous/mclInput\n");
 		buf.write("mcl Results/Miscellaneous/mclInput --abc -I 1.5 -o Results/Miscellaneous/mclOutput\n");
		buf.write("perl bin/orthomclSoftware-v2.0.9/bin/orthomclMclToGroups SI 1000 < Results/Miscellaneous/mclOutput > Results/Miscellaneous/Orthologous_groups.txt\n");
		b.close();buf.flush();buf.close();

		runtime = Runtime.getRuntime();
		process = runtime.exec("sh bin/commands");
		reader = new BufferedReader(new InputStreamReader(process.getErrorStream()));
		lines="";
		while((lines=reader.readLine())!=null)
		{
			//System.out.println(lines);
		}
		process.waitFor();
		runtime.gc();
		reader.close();
		BufferedReader test =  new BufferedReader (new FileReader("Results/Miscellaneous/Orthologous_groups.txt"));
		if(test.readLine()==null)
		{
			System.out.println("Orthomcl analyses failed... Check error output and make necessary modifications");
			return;
		}

		//delete pairs folder
		fasta = new File("pairs");
		deleteFolder.delete(fasta);
				
		//3. gff_parse
		System.out.println("Parsing gff......");
		String temp="";
		for(int i=0; i<gff.size();i++)
		{
			try
			{
				temp+=gff_parse.parse("Data/Gffs/"+gff.get(i))[0];
			}
			catch(Exception e){System.out.println("Error in parsing GFF files"); return;}
		}
		buf = new BufferedWriter(new FileWriter("Results/Miscellaneous/geneID_mRNA-name_all.txt"));
		buf.write(temp);buf.flush();buf.close();
	
		//4. Create fasta files for meme
		System.out.println("Running MEME......");
		int counter = CreateFasta_for_MEME.file_prep("Data/IGRs/Combined_IGRs.txt", "Results/Miscellaneous/geneID_mRNA-name_all.txt","Results/Miscellaneous/Orthologous_groups.txt", target, "Data/bg/Background_dist_"+order+"_order.txt");
		if(counter==-1)
		{
			System.out.println("Error while creating fasta files for meme analysis... Exiting");
			return;
		}


		//5. Run MEME analysis
		runtime = Runtime.getRuntime();
		process = runtime.exec("sh bin/meme_bash");
		reader = new BufferedReader(new InputStreamReader(process.getErrorStream()));
		lines="";
		while((lines=reader.readLine())!=null)
		{
			//System.out.println(lines);
		}
		process.waitFor();
		runtime.gc();
		reader.close();
		//6. Create mast bash file
		int counter2 = MemeParse_forMast.parse(target, counter,"bin/mast_bash", "Results/Miscellaneous/Gene_motif_map.txt","Results/Motif_finding/meme_out", "Data/bg/Background_dist_target_2_order.txt",pvalue);

		if(counter2==-1)
		{
			System.out.println("Error while creating fasta files for mast analysis... Exiting");
			return;
		}

		//7. Run MAST analysis
		System.out.println("Running mast analysis with phylo motifs...");
		runtime = Runtime.getRuntime();
		process = runtime.exec("sh bin/mast_bash");
		reader = new BufferedReader(new InputStreamReader(process.getErrorStream()));
		lines="";
		while((lines=reader.readLine())!=null)
		{
			//System.out.println(lines);
		}
		process.waitFor();
		runtime.gc();
		reader.close();
		//8. Parse MAST results and create new input for Mcl
		test2=Motif_bindsite_overlap.overlap(counter2,"Results/Miscellaneous/mclIn.txt",IGR_length);
		if(test2==-1)
		{
			System.out.println("Error while preping files for motif clustering... Exiting");
			return;
		}

		//9. Cluster motifs with mcl
		runtime = Runtime.getRuntime();
		process = runtime.exec("mcl Results/Miscellaneous/mclIn.txt --abc -I 1.5 -o Results/Miscellaneous/mclOut.txt");
		reader = new BufferedReader(new InputStreamReader(process.getErrorStream()));
		process.waitFor();
		lines="";
		while((lines=reader.readLine())!=null)
		{
			//System.out.println(lines);
		}
		runtime.gc();
		reader.close();
		
		//10. Pool together clustered motifs and build final PSSMs
		fasta = new File("Results/Motif_finding/FASTA2");
		deleteFolder.delete(fasta);
		fasta.mkdir();
		fasta = new File("Results/Motif_finding/meme_out2");
		deleteFolder.delete(fasta);
		fasta.mkdir();
		size = Merge_motifs.merge("Results/Miscellaneous/mclOut.txt","Results/Motif_finding/FASTA2",target,"Data/bg/Background_dist_target_2_order.txt","Results/Miscellaneous/Groups.txt");
		if(size==-1)
		{
			System.out.println("Error while merging motifs... Exiting");
			return;
		}
		
		//further clustering
		test2=MemeParse_forTOMTOM2.parse(size, "Results/Motif_finding/meme_out2","Results/Miscellaneous/All_initial_motifs_TOMTOM.txt");
		if(test2==-1)
		{
			System.out.println("Error while parse meme files for tomtom analysis... Exiting");
			return;
		}
		System.out.println("Running tomtom analysis...");
		runtime = Runtime.getRuntime();
		process = runtime.exec("tomtom -oc Results/Motif_finding/tomtom_out -thresh 0.001 Results/Miscellaneous/All_initial_motifs_TOMTOM.txt Results/Miscellaneous/All_initial_motifs_TOMTOM.txt");
		reader = new BufferedReader(new InputStreamReader(process.getErrorStream()));
		lines="";
		while((lines=reader.readLine())!=null)
		{
		}
		process.waitFor();
		runtime.gc();
		reader.close();

		Tomtom_parse.parse("Results/Motif_finding/tomtom_out/tomtom.txt","Results/Miscellaneous/mclIn2.txt");
		runtime = Runtime.getRuntime();
		process = runtime.exec("mcl Results/Miscellaneous/mclIn2.txt --abc -I 1.5 -o Results/Miscellaneous/mclOut2.txt");
		reader = new BufferedReader(new InputStreamReader(process.getErrorStream()));
		process.waitFor();
		lines="";
		while((lines=reader.readLine())!=null)
		{
		}
		runtime.gc();
		reader.close();

		fasta = new File("Results/Motif_finding/FASTA3");
		deleteFolder.delete(fasta);
		fasta.mkdir();
		fasta = new File("Results/Motif_finding/meme_out3");
		deleteFolder.delete(fasta);
		fasta.mkdir();
		size = Merge_motifs2.merge("Results/Miscellaneous/mclOut2.txt","Results/Motif_finding/FASTA3",target,"Data/bg/Background_dist_target_2_order.txt","Results/Miscellaneous/preClusters.txt");
		if(size==-1)
		{
			System.out.println("Error while merging motifs... Exiting");
			return;
		}
		
		test2=MemeParse_clusters.parse(size,"Results/Clusters.txt");
		if(test2==-1)
		{
			System.out.println("Error while parse meme files to obtain cluster membership... Exiting");
			return;
		}

		//11. Parse MEME for Tomtom and run mcl
		test2=MemeParse_forTOMTOM2.parse(size, "Results/Motif_finding/meme_out3","Results/All_motifs_TOMTOM.txt");
		if(test2==-1)
		{
			System.out.println("Error while parse meme files for tomtom analysis... Exiting");
			return;
		}
		System.out.println("Running tomtom analysis...");
		runtime = Runtime.getRuntime();
		process = runtime.exec("tomtom -oc Results/Motif_finding/tomtom_out2 -thresh 0.1 bin/All_motifs_Ecoli Results/All_motifs_TOMTOM.txt");
		reader = new BufferedReader(new InputStreamReader(process.getErrorStream()));
		lines="";
		while((lines=reader.readLine())!=null)
		{
		}
		process.waitFor();
		runtime.gc();
		reader.close();	

		//12. Run gene set enrichment analysis
		System.out.println("Running gene set enrichment analysis (GO)...");
		if(!goa.equals(""))
		{
			int genomeSize = Integer.parseInt(gff_parse.parse("Data/Gffs/"+gff.get(abbr.indexOf(target)))[1]);
			if(genomeSize==-1)
			{
				System.out.println("Error while parsing gff files for GSEA... Exiting");
				return;
			}
			System.out.println("Number of ORFs (genes) : "+genomeSize);
			GO_enrichment.significance("Results/Clusters.txt", "Data/GO_terms/"+goa, genomeSize, "Results/GSEA_analysis_phylo.txt");
		}

		//13. Print out summary webpage of phylogenetic footprinting analysis
		Summary.summarize();
		gene_cluster_map.format();

		//14. Begin integration with expression data using DISTILLER
		if(!expression.equals(""))//only run DISTILLER analysis if expression matrix provided
		{
			fasta = new File("bin/DISTILLER-V2/dataFolder");
			if(fasta.exists()) deleteFolder.delete(fasta);
			if(!fasta.exists()) fasta.mkdir();

			fasta = new File("bin/DISTILLER-V2/outputInitial.m");
			if(fasta.exists()) fasta.delete();
				
			System.out.println("Building motif matrix for DISTILLER");
			test2 = DISTILLER_motif_matrix.build(expression);
			if(test2==-1)
			{
				System.out.println("Error: Exiting run...");
				return;
			}

			boolean conf=true;
			while(conf)
			{
				int count1=0,count2=0;
				System.out.println("Running DISTILLER...");//need to write code to ensure distiller works well enough
				runtime = Runtime.getRuntime();
				process = runtime.exec("java -jar -Xmx8000m bin/DISTILLER-V2/Miner.jar bin/DISTILLER-V2/input.txt");
				reader = new BufferedReader(new InputStreamReader(process.getErrorStream()));
				lines="";
				while((lines=reader.readLine())!=null)
				{
					
				}
				process.waitFor();
				runtime.gc();
				reader.close();
				reader = new BufferedReader(new FileReader("bin/DISTILLER-V2/outputInitial.m"));
				lines="";
				while((lines=reader.readLine())!=null)
				{
					if(lines.contains("Significances")) count1++;
				}
				

				System.out.println("Running DISTILLER... phase2");
				runtime = Runtime.getRuntime();
				process = runtime.exec("java -jar -Xmx8000m bin/DISTILLER-V2/Filter.jar bin/DISTILLER-V2/inputModuleSelection.txt");
				reader = new BufferedReader(new InputStreamReader(process.getErrorStream()));
				lines="";
				while((lines=reader.readLine())!=null)
				{
					//System.out.println(lines);
				}
				process.waitFor();
				runtime.gc();
				reader.close();

				reader = new BufferedReader(new FileReader("Results/DISTILLER/DISTILLER_modules.txt"));
				lines="";
				while((lines=reader.readLine())!=null)
				{
					if(lines.contains("NaN")) count2++;
				}
				if(count1>=300 && count2==0) conf=false; //if DISTILLER works ok move on else retry 
				reader.close();
			}
			
			
			//parse DISTILLER output
			System.out.println("Parsing DISTILLER results");
			test2=DISTILLER_parse.parse();
			if(test2==-1)
			{
				System.out.println("Error: error occurred parsing DISTILLER results...Exiting run.");
				return;
			}

			//Extend seed modules
			System.out.println("Seed extending...");
			test2=DISTILLER_extend_corr.extend(expression);
			if(test2==-1)
			{
				System.out.println("Error: error occurred during seed extension...Exiting run.");
				return;
			}

			//Operon extend
			System.out.println("Operon extension (only if valid operon file provided...)");
			test2=DISTILLER_extend_corr_operons.extend(expression, operon);
			if(test2==-1)
			{
				System.out.println("Error: error occurred during operon extension...Exiting run.");
				return;
			}

			//Get sequences for all TFs to be used for Pfam analysis by user
			TF_seq_retrieve.fetch(proteomes.get(abbr.indexOf(target)));

			s = new Scanner(System.in);

			System.out.print("Do you have results of your pfam analysis of TF sequences ready, appropriately named and placed in the Data directory? [y/n]: ");
			String answer = s.nextLine();
			
			if(answer.equalsIgnoreCase("y") || answer.equalsIgnoreCase("yes"))
			{
				System.out.println("Proceeding with analysis...");
				test2=MemeParse_forTOMTOM3.parse("Results/Integrated_modules_final.txt", "Results/Motif_finding/meme_out3/");
				if(test2==-1)
				{
					System.out.println("Error: Parsing meme file for tomtom analysis (distiller results)...");
					return;
				}
				System.out.println("Running tomtom analysis...");
				runtime = Runtime.getRuntime();
				process = runtime.exec("tomtom -oc Results/DISTILLER/tomtom_out -thresh 1 Results/DISTILLER/DISTILLER_motifs_TOMTOM.txt bin/All_motifs_Ecoli");
				reader = new BufferedReader(new InputStreamReader(process.getErrorStream()));
				lines="";
				while((lines=reader.readLine())!=null)
				{
				}
				process.waitFor();
				runtime.gc();
				reader.close();

				//compute corr score, dist score and dbd scores
				try
				{
					Thread[] threadHolder=new Thread[3];
					Thread trial = new Thread(new regulator_predict_distance(gff.get(abbr.indexOf(target))));
					trial.start();
					threadHolder[2]=trial;
					trial = new Thread(new regulator_predict_corr(gff.get(abbr.indexOf(target)),expression));
					trial.start();
					threadHolder[1]=trial;
					trial = new Thread(new TF_Domain_family("TF_domains.txt"));
					trial.start();
					threadHolder[0]=trial;
					for(int i=0;i<threadHolder.length;i++)
						threadHolder[i].join();//ensure all thread complete before moving to next step
				}
				catch(Exception e)
				{e.printStackTrace(); System.out.println("Error computing DBD/Correlation/Distances scores...");return;}	

				//compute R score
				test2=regulator_predict_merge.merge();
				if(test2==-1)
				{
					System.out.println("Error: Exiting run...");
					return;
				}

				//Run gene set enrichment analysis
				if(!goa.equals(""))
				{
					int genomeSize = Integer.parseInt(gff_parse.parse("Data/Gffs/"+gff.get(abbr.indexOf(target)))[1]);
					if(genomeSize==-1)
					{
						System.out.println("Error while parsing gff files for GSEA... Exiting");
						return;
					}
					System.out.println("Number of ORFs (genes) : "+genomeSize);
					GO_enrichment.significance("Results/Integrated_modules_final.txt", "Data/GO_terms/"+goa, genomeSize, "Results/GSEA_analysis_final.txt");
				}

				//Build summary report
				System.out.println("Summarizing...");
				Summary2.summarize();
			}
			else
			{
				System.out.println("Exiting run...");
				return;
			}
			AUPR_preprocessing.sortInteractionList();
		}
		
		System.out.println("Run finished at :"+Calendar.getInstance().getTime());

	}
	catch(Exception e)
	{e.printStackTrace(); }
	}
}
