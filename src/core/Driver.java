package core;
import java.util.Scanner;
import tinker.ARC_Parser;
import tinker.XYZ_Parser;

/**
 * The driver is used to execute any implemented parser or analysis operation
 * @author Trevor Fisher
 *
 */
public class Driver {
	
	/**
	 * The method called when the program starts. It implements a 
	 * @param args	Currently unused, but eventually the args can be used to specify which operation
	 * 				to run when calling the executable
	 */
	public static void main(String[] args) {
		Scanner sc = new Scanner(System.in);
		
		//Read in file name of file
		System.out.println("Enter the name of the coordinate file:");
		String fileName = sc.nextLine();
		
		arc_Analysis(fileName);
		//xyz_Analysis(fileName);
		//TODO Implement more calculations here:
		
		sc.close();
	}
	
	/**
	 * Begins an operation to parse a Tinker .arc file and perform calculations specified in the <code>parseFile</code> method.
	 * @param fileName		Provided by the scanner in the driver, if no file path is specified, it defaults
	 * 						to the current directory.
	 */
	public static void arc_Analysis(String fileName)
	{
		BulkSystem simulation = new BulkSystem();
		ARC_Parser parser = new ARC_Parser(fileName, simulation, 10000, 12.5, 0.03);
		parser.parseFile();
	}
	
	/**
	 * OUTDATED: Because the Tinker .arc and .xyz files are the same format, the ARC_RDF function
	 * 			can also parse the .xyz file if the number of systems in the Arc_Parser class is specified
	 * 			to 1. Parses a Tinker .xyz file and performs a radial distribution function calculation.
	 * @param fileName	Provided by the scanner in the driver, if no file path is specified, it defaults
	 * 					to the current directory.	
	 */
	public static void xyz_Analysis(String fileName)
	{

		BulkSystem simulation = new BulkSystem();
		
		XYZ_Parser parser = new XYZ_Parser(fileName, simulation);
		parser.parseFile();
		parser.rDF(simulation, 12.5, 0.25);
	}
}
