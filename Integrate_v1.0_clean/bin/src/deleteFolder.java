package bin;
import java.lang.*;
import java.*;
import java.util.*;
import java.io.*;
/*
Code to delete entire folders


Written by Saheed Imam
Institute for Systems biology
*/

public class deleteFolder
{
	public static void delete(File f) throws IOException 
	{
		if (f.isDirectory())
		{
			for (File c : f.listFiles())
				delete(c);
		}
		f.delete();
	}	
}
