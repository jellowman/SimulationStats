package tinker;

import java.io.IOException;
import java.util.Arrays;

import calculations.RadialDF;
import core.BulkSystem;

/**
 * OUTDATED - The <code>ARC_Parser</code> class can handle .xyz files as well since a Tinker
 * .arc file is made up of Tinker .xyz files. This class will no longer be updated, and it is suggested
 * to use the Arc_Parser class with <code>numSystems</code> set to 1.
 * @see ARC_Parser
 * @see Tinker_Parser
 * @author Trevor Fisher
 *
 */
public class XYZ_Parser extends Tinker_Parser
{
	
	public XYZ_Parser(String file, BulkSystem system)
	{
		super(file, system);
	}
	
	/**
	 * Parses a Tinker-formatted xyz file for the atom coordinates and atom types
	 */
	public void parseFile()
	{
		try {
			String nextLine = br.readLine();
			nextLine = nextLine.trim();
			system.setNumAtoms(Integer.valueOf(nextLine));
			
			nextLine = br.readLine();
			parseCoords(nextLine, system);
			
			nextLine = br.readLine();
			while(nextLine != null && !nextLine.isEmpty())
			{
				try {
					parseAtom(nextLine, system);
				} catch (BadAtomFormat e) {
					System.err.println(file.getName() + " has improper syntax");
					System.exit(1);
				}
				nextLine = br.readLine();
			}
		}
		catch (IOException ex)
		{
			System.err.println("Error in parsing file.");
		}
	}
	
	public void rDF(BulkSystem simulation, double cutoff, double interval)
	{
		//Set up RDF calculator
		RadialDF calc1 = new RadialDF(null, null, interval);
		
		int size = (int) Math.floor(cutoff/interval);
		int[] distances = new int[size];
		Arrays.fill(distances, 0);
		double[] normalizedDistances;
		
		calc1.systemRDFCalculation("OW", "OW", simulation, distances, 12.5);

		normalizedDistances = calc1.normalizeRDF(distances);
		
		for(int i = 0; i < normalizedDistances.length; i++)
		{
			System.out.println("Layer:\t" + i*interval + "\t\t" + normalizedDistances[i]);
		}
	}
}
