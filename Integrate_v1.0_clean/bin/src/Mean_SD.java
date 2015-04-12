package bin;
import java.lang.*;
import java.*;
import java.util.*;
/*
Code to compute the mean and standard deviation of an array

java Mean_SD
*/

public class Mean_SD
{
	public static double[] compute(double[] data)
	{
	double mean=0.0,sd=0.0;
	try
	{
		double sum=0.0;
		for(int i=0;i<data.length;i++)
			sum+=data[i];
		mean = sum/(double)data.length;
		sum=0.0;
		for(int i=0;i<data.length;i++)
			sum+=Math.pow((data[i]-mean),2);
		sd = Math.pow((sum/((double)data.length-1)),0.5);
		
	}
	catch(Exception e)
	{e.printStackTrace(); }
	return new double[]{mean,sd};
	}
	
}
