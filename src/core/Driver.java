package core;

import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Scanner;

import calculations.AtomDensity;
import calculations.HBondAngleDist;
import calculations.HBondLife;
import calculations.HBondLifeAC;
import calculations.HBondLifeACT;
import calculations.HBondLifeSS;
import calculations.HBondLifeTC;
import calculations.HBondLifeTrajWriter;
import calculations.RadialDF2;
import lammps.DatWriter;
import pdb.PDB_Parser;
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
		//Set up a text input instead of reading from console
		/*try{
		System.setIn(new FileInputStream("input.txt"));
		} catch(FileNotFoundException e) {
			System.err.println("Input file not found.");
			System.exit(0);
		}
		Scanner sc = new Scanner(System.in);*/
		
		//Read in file name of file
		//System.out.println("Enter the name of the coordinate file:");
		//String fileName = sc.nextLine();
		
		//arc_Analysis(fileName);
		//xyz_Analysis(fileName);
		
		//TODO Implement more calculations here:
		//pdbToLmpDat(sc);
		
		//HBondLife hL = new HBondLife();
		//hL.runAnalysis();
		
		//HBondLifeTC hLTC = new HBondLifeTC();
		//hLTC.runAnalysis();
		
		HBondLifeAC hLAC = new HBondLifeAC("2", "2", 3.8, 35, 10);
		hLAC.runAnalysis();
		
		//HBondLifeACT hLACT = new HBondLifeACT();
		//hLACT.runAnalysis();
		
		//HBondLifeSS hLSS = new HBondLifeSS();
		//hLSS.runAnalysis();
		
		//HBondLifeTrajWriter hLTW = new HBondLifeTrajWriter();
		//hLTW.runAnalysis();
		
		
		//RadialDF2 rdf2 = new RadialDF2(6.0, "2", "1");
		//rdf2.runAnalysis();
		
		//HBondAngleDist hbad = new HBondAngleDist(3.85, "2", "1");
		//hbad.runAnalysis();
		
		//AtomDensity ad = new AtomDensity();
		//ad.runAnalysis();
		//sc.close();
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
	
	/**
	 * Converts a .pdb file into a Lammps .dat input file while conserving bond and angle definitions.
	 * Requires both the system .pdb file and individual .pdb files for each molecule containing bond and
	 * angle associations. 
	 * @param sc	The console prompt
	 */
	public static void pdbToLmpDat(Scanner sc)
	{
		PDB_Parser par = new PDB_Parser();
		
		//Parse the blueprint files containing atoms, bonds, angles
		while(true)
		{
			System.out.println("Enter the name of the file containing a single molecule, or type \"DONE\"");
			String bPName = sc.nextLine();
			
			if(bPName.equals("DONE"))
			{
				break;
			}
			
			par.parseBlueprintFile(bPName);
		}
		
		//Parse the file containing the entire system
		String sysName = null;
		boolean badFile = true;
		while(badFile)
		{
			System.out.println("Enter the name of the system file.");
			sysName = sc.nextLine();
			Path path = Paths.get(sysName);
			if(Files.exists(path))
			{
				badFile = false;
			}
			else
			{
				System.err.println("File does not exist. Please enter a valid file name.");
			}
		}
		
		par.parseSystemFile(sysName);
		
		//Write the molecules to a Lammps .dat file
		ArrayList<Molecule> mols = par.getMolecules();		
		DatWriter.writeFile(mols, sc);
	}
}
