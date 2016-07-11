package tinker;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;

import core.Atom;
import core.BulkSystem;

/**
 * A class for defining the structure of Tinker formatted files. Tinker is a MD software created by
 * Jay Ponder Lab, Department of Chemistry, at Washington University, Saint Louis, Missouri.
 * Information about Tinker can be found <a href="http://dasher.wustl.edu/tinker/">here</a>.
 * @see core.BulkSystem
 * @author Trevor Fisher
 *
 */
public abstract class Tinker_Parser 
{
	protected File file;
	protected FileReader fr;
	protected BufferedReader br;
	
	/**The system that the parser will store the atoms into*/
	protected BulkSystem system;
	
	/**
	 * Constructs the readers for the file with the name that matches <code>file</code>. It also takes in a
	 * new or existing <code>BulkSystem</code> to read the new atoms into.
	 * @see core.BulkSystem
	 * @param file		The filename of the Tinker file.
	 * @param system	The <code>BulkSystem</code> to read the atoms from the file into.
	 */
	public Tinker_Parser(String file, BulkSystem system)
	{
		this.file = new File(file);
		
		try {
			this.fr = new FileReader(this.file);
		} 
		catch (FileNotFoundException e) {
			System.err.println("File not found");
		}
		
		this.br = new BufferedReader(this.fr);
		
		this.system = system;
	}
	
	/**
	 * Handles the specific file extension's format for parsing out coordinates and atoms.
	 */
	public abstract void parseFile();
	
	/**
	 * Stores the simulation box dimensions in <code>bs</code>.
	 * @param nextLine		The line containing the box dimensions.
	 * @param bs			The <code>BulkSystem</code> to store the coordinates.
	 */
	public static void parseCoords(String nextLine, BulkSystem bs)
	{
		String[] cutLine = nextLine.split(" ");
		String[] parts = new String[20];
		int j = 0;
		for(int i = 0; i < cutLine.length; i++)
		{
			if(cutLine[i].length() != 0)
			{
				parts[j] = cutLine[i];
				j++;
			}
		}
		
		bs.setXBox(Double.valueOf(parts[0]));
		bs.setYBox(Double.valueOf(parts[1]));
		bs.setZBox(Double.valueOf(parts[2]));
	}
	
	/**
	 * Parses all the containing atom information from a line and stores it in <code>bs</code>
	 * @param nextLine			The line containing the atom information.
	 * @param bs				The <code>BulkSystem</code> to store the coordinates.
	 * @throws BadAtomFormat	If the line does not match the Tinker format, then this is thrown.
	 */
	public static void parseAtom(String nextLine, BulkSystem bs) throws BadAtomFormat
	{
		String[] cutLine = nextLine.split(" ");
		String[] parts = new String[50];
		int j = 0;
		for(int i = 0; i < cutLine.length & i < 50; i++)
		{
			if(!cutLine[i].isEmpty())
			{
				parts[j] = cutLine[i];
				j++;
			}
		}
		try{
		String type = parts[1];
		double xCord = Double.valueOf(parts[2]);
		double yCord = Double.valueOf(parts[3]);
		double zCord = Double.valueOf(parts[4]);
		Atom atom = new Atom(type, xCord, yCord, zCord);
		bs.add(atom);
		}catch (NullPointerException e)
		{
			System.err.println("Null Pointer err - Bad line:\n" + nextLine);
			throw new BadAtomFormat();
		}catch (NumberFormatException e)
		{
			System.err.println("Number Format err - Bad line:\n" + nextLine);
			throw new BadAtomFormat();
		}
	}
}
