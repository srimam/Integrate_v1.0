package bin;
import java.lang.*;
import java.*;
import java.io.*;
import java.util.regex.*;
import java.util.ArrayList;
import java.util.*;
import bin.*;
/*
Code to initialize folders and orthomcl config file for new run

java initialize

Written by Saheed Imam
Institute for Systems biology
*/

public class initialize
{
	public static int init(String username,String password) throws IOException
	{
		try
		{
		File fasta;
		Scanner s;
		//create Results folder and subfolders
		fasta = new File("Results");
		deleteFolder.delete(fasta);
		fasta.mkdir();

		fasta = new File("Results/Miscellaneous");
		fasta.mkdir();

		fasta = new File("Results/DISTILLER");
		fasta.mkdir();

		fasta = new File("Results/Motif_finding");
		fasta.mkdir();

		fasta = new File("Data");
		if(!fasta.exists()) fasta.mkdir();

		fasta = new File("Data/bg");
		if(!fasta.exists()) fasta.mkdir();

		fasta = new File("Data/Expression");
		if(!fasta.exists()) fasta.mkdir();

		fasta = new File("Data/Genomes");
		if(!fasta.exists()) fasta.mkdir();

		fasta = new File("Data/Gffs");
		if(!fasta.exists()) fasta.mkdir();

		fasta = new File("Data/GO_terms");
		if(!fasta.exists()) fasta.mkdir();

		fasta = new File("Data/IGRs");
		if(!fasta.exists()) fasta.mkdir();

		fasta = new File("Data/Operons");
		if(!fasta.exists()) fasta.mkdir();

		fasta = new File("Data/Proteomes");
		if(!fasta.exists()) fasta.mkdir();

		fasta = new File("Data/Proteomes/compliantFasta");
		if(fasta.exists()) deleteFolder.delete(fasta);

		fasta = new File("Data/Proteomes/all_v_all_proteomes");
		if(fasta.exists()) fasta.delete();

		fasta = new File("Data/Proteomes/SimilarSequences.txt");
		if(fasta.exists()) fasta.delete();

		fasta = new File("Data/Proteomes/goodProteins.fasta");
		if(fasta.exists()) fasta.delete();

		fasta = new File("bin/DISTILLER-V2/outputInitial.m");
		if(fasta.exists()) fasta.delete();

		fasta = new File("Data/Proteomes/compliantFasta");
		if(!fasta.exists()) fasta.mkdir();

		fasta = new File("Data/TFs.txt");
		if(!fasta.exists())
		{
			s = new Scanner(System.in);
			System.out.println("Please provide a list of TFs with appropriate IDs for your organism in a file called \"TFs.txt\" and place this in the \"Data\" directory. Press ENTER when ready. If only interested in phylogenetic footprinting analysis just press ENTER");
			s.nextLine();
		}
		
		//Create orthomcl config file
		BufferedWriter buf =new BufferedWriter(new FileWriter("bin/orthomcl.config"));
		buf.write("# this config assumes a mysql database named 'orthomcl'.  adjust according");buf.newLine();
		buf.write("# to your situation.");buf.newLine();
		buf.write("dbVendor=mysql ");buf.newLine();
		buf.write("dbConnectString=dbi:mysql:orthomcl");buf.newLine();
		buf.write("dbLogin="+username);buf.newLine();
		buf.write("dbPassword="+password);buf.newLine();
		buf.write("similarSequencesTable=SimilarSequences");buf.newLine();
		buf.write("orthologTable=Ortholog");buf.newLine();
		buf.write("inParalogTable=InParalog");buf.newLine();
		buf.write("coOrthologTable=CoOrtholog");buf.newLine();
		buf.write("interTaxonMatchView=InterTaxonMatch");buf.newLine();
		buf.write("percentMatchCutoff=50");buf.newLine();
		buf.write("evalueExponentCutoff=-5");buf.newLine();
		buf.write("oracleIndexTblSpc=NONE");
		buf.flush();
		}
		catch(Exception e)
		{e.printStackTrace(); return -1;}
		return 0;
	}
	
}
