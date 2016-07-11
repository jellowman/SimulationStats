package tinker;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.Arrays;

import calculations.RadialDF;
import core.BulkSystem;

/**
 * Parses a Tinker .arc file and performs calculations on each configuration contained in the file. v1.1 is updated
 * to match the changes to v1.1 of the <code>RadialDF</code> class.
 * @see Tinker_Parser
 * @see calculations.RadialDF
 * @see core.BulkSystem
 * @author Trevor Fisher
 * @version 1.1
 */
public class ARC_Parser extends Tinker_Parser
{
	private int numSystems=500;
	private double interval=0.25;
	private double cutoff=12.5;
	
	public ARC_Parser(String file, BulkSystem system, int numSystems, double cutoff, double interval) 
	{
		super(file, system);
		this.numSystems = numSystems;
		this.interval = interval;
		this.cutoff = cutoff;
	}
	
	/**
	 * Parses the .arc file and performs analysis on each configuration and post-analysis over the entire set of configurations.
	 * Note the TO-DO sections where additional analysis can be implemented, as well as changes to the parser if a different
	 * ensemble other than NVT was used.
	 */
	public void parseFile() 
	{
		try {
			String nextLine = br.readLine();
			nextLine = nextLine.trim();
			int numAtoms = Integer.parseInt(nextLine);
			system.setNumAtoms(numAtoms);
			//Simulation box size only needs to be parsed once as long as an NVT Ensemble was used
			nextLine = br.readLine();
			parseCoords(nextLine, system);
			
			//TODO Set up initial variables that will need to be carried over each configuration
			
			//Set up RDF calculator
			RadialDF loopCalc = new RadialDF(null, null, interval);
			
			//Set up RDF shells
			int size = (int) Math.floor(cutoff/interval);
			int[] distances = new int[size];
			Arrays.fill(distances, 0);
			
			//How many configurations to iterate through
			
			//Used to determine how many simulations were actually calculated in order
			//to obtain the correct RDF average among all simulations
			int numGoodIterations = numSystems;
			for(int config = 0; config < numSystems && nextLine!=null; config++)
			{
				System.out.println("Calculating system " + (config+1));
				//Make a new system to hold the current configuration
				BulkSystem currentSystem = new BulkSystem();
				
				/* To speed up parsing, the dimensions of the box and number of atoms
				are grabbed from the initial system.
				
				TODO: Must be changed if the simulation did not use the NVT ensemble */
				/*nextLine = br.readLine();
				nextLine = nextLine.trim();
				currentSystem.setNumAtoms(Integer.parseInt(nextLine));
				nextLine = br.readLine();
				parseCoords(nextLine, currentSystem); */
				currentSystem.setNumAtoms(numAtoms);
				currentSystem.setXBox(system.getXBox());
				currentSystem.setYBox(system.getYBox());
				currentSystem.setZBox(system.getZBox());
				
				boolean goodSet=true;
				//NOTE: arc file must be clean (did not stop during a writeout)
				//Otherwise, error handling will discard the current simulation and move onto the next one
				for(int i = 0; i < numAtoms; i++)
				{
					nextLine = br.readLine();
					try{
					parseAtom(nextLine, currentSystem);
					}
					catch(BadAtomFormat bf)
					{
						goodSet=false;
						br.readLine();
						System.err.println("Bad atom in config " + (config+1));
						numGoodIterations--;
						break;
					}
				}
				
				if(goodSet)
				{
					system = currentSystem;
					//TODO add calculations to do on each configuration here:
					
					loopCalc.systemRDFCalculation("OW", "OW", currentSystem, distances, cutoff);
					
					//Read through the next header of the next system (numAtoms then coordinates)
					//TODO Must remove these lines if you are parsing the number of atoms or volume for each configuration!
					br.readLine();
					nextLine = br.readLine();
				}
			}
			
			//TODO add processing after all configurations have been analyzed here:
			
			//Divide the increments by the density and the number of atoms in a system);
			double[] normalizedDistances = loopCalc.normalizeRDF(distances, numGoodIterations);
			
			//Print out RDF results to console
			for(int i = 0; i < normalizedDistances.length; i++)
			{
				DecimalFormat df = new DecimalFormat("#.##");
				System.out.println("Layer:\t" + df.format(i*interval) + "\t\t" + normalizedDistances[i]);
			}
			
		} catch (IOException e) {
			System.out.println("Error in parsing file");
		}
	}
	
}
