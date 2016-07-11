package calculations;

import core.Atom;
import core.BulkSystem;
/**
 * Calculates a radial distribution function for a system, and has support
 * for averaging over many different configurations. Version 1.1 fixes a major bug in the calculation of the distances between atoms
 * where it would only use the 0 indexed central atom after the current atom index reached the current atom size.
 * @author Trevor Fisher
 * @version 1.1
 *
 */
public class RadialDF 
{
	/**Used to determine number of RDF functions calculated*/
	private BulkSystem centerAtoms;
	/**Used for determining atomic density*/
	private BulkSystem otherAtoms;
	/**The radius of each shell used in the RDF*/
	private double increment;
	
	public RadialDF(BulkSystem centr, BulkSystem other, double increment)
	{
		this.centerAtoms = centr;
		this.otherAtoms = other;
		this.increment = increment;
	}
	
	/**
	 * Adds to an array of radial shells extending from <code>atomIndex</code>
	 * @param atomIndex		The atom treated as r1
	 * @param distances		Iterates each section by 1 if a distance falls within that shell radius
	 */
	public void rDFCalculation(int atomIndex, int[] distances)
	{
		//Ensure that the distance between itself is not calculated
		for(int i = 0; i < atomIndex; i++)
		{
			double r = periodicDistance(centerAtoms.get(atomIndex), otherAtoms.get(i));
			int radiusLayer = (int)(r/increment);
			if(radiusLayer < distances.length)
				distances[radiusLayer]++;
		}
		for(int i = atomIndex+1; i < otherAtoms.size(); i++)
		{
			double r = periodicDistance(centerAtoms.get(atomIndex), otherAtoms.get(i));
			int radiusLayer = (int)(r/increment);
			if(radiusLayer < distances.length)
				distances[radiusLayer]++;
		}
	}
	
	/**
	 * Calculates the RDF for each atom in the system.
	 * @param centerAtom		The central atom type for a RDF calculation
	 * @param otherAtom			The type of atom that will be measured for the distance from the centerAtom
	 * @param currentSystem		The current configuration holding all the atoms within it.
	 * @param distances			Keeps track of the number of distance calculations that falls within each RDF shell.
	 * @param cutoff			The maximum distance that the RDF function will calculate.
	 * 							Note that this should be half the length of the smallest system box side at most.
	 */
	public void systemRDFCalculation(String centerAtom, String otherAtom, BulkSystem currentSystem, int[] distances, double cutoff)
	{
		//Clone the system and remove all atoms except for the two to compare
		//BulkSystem systemCopy = (BulkSystem) currentSystem.clone();
		//systemCopy.filterAtomType(atomType);
		BulkSystem centerAtoms = (BulkSystem) currentSystem.clone();
		centerAtoms.clear();
		BulkSystem otherAtoms = new BulkSystem();
		
		for(Atom atom : currentSystem)
		{
			if(atom.getType().equals(centerAtom))
					centerAtoms.add(atom);
			if(atom.getType().equals(otherAtom))
					otherAtoms.add(atom);
		}
		
		this.centerAtoms = centerAtoms;
		this.otherAtoms = otherAtoms;
		
		for(int i = 0; i < centerAtoms.size(); i++)
		{
			//Increments number in each distance shell
			this.rDFCalculation(i, distances);
		}
	}
	
	/**
	 * Normalizes the relative density for each RDF shell to the system density.
	 * Should be used after all configurations have been tallied in <code>distances</code>.
	 * @param distances				Keeps track of the number of distance calculations that falls within each RDF shell.
	 * @param numGoodIterations		The number of configurations that were tallied using <code>systemRDFCalculation</code>.
	 * @return						An array containing the relative density for each RDF shell.
	 */
	public double[] normalizeRDF(int[] distances, int numGoodIterations)
	{
		double[] normalDistances = new double[distances.length];
		double aDensity = getAtomDensity();
		for(int i = 0; i < normalDistances.length; i++)
		{
			double volume = getShellVolume(i+1, increment);
			normalDistances[i] = distances[i]/volume/aDensity/centerAtoms.size()/numGoodIterations;
		}
		
		return normalDistances;
	}
	
	/**
	 * Normalizes the relative density for each RDF shell to the system density. Assumes only one configuration was tallied.
	 * @param distances		Keeps track of the number of distance calculations that falls within each RDF shell.
	 * @return				An array containing the relative density for each RDF shell.
	 */
	public double[] normalizeRDF(int[] distances)
	{
		return normalizeRDF(distances, 1);
	}
	
	/**
	 * Calculates the distance between two atoms, accounting for periodic
	 * conditions.
	 * @param a		Atom a
	 * @param b		Atom b
	 * @return		The radial distance between the two atoms
	 */
	private double periodicDistance(Atom a, Atom b)
	{
		double xDist = b.getXCord() - a.getXCord();
		double yDist = b.getYCord() - a.getYCord();
		double zDist = b.getZCord() - a.getZCord();
		
		xDist = Math.round(xDist/centerAtoms.getXBox()) * centerAtoms.getXBox() - xDist;
		yDist = Math.round(yDist/centerAtoms.getYBox()) * centerAtoms.getYBox() - yDist;
		zDist = Math.round(zDist/centerAtoms.getZBox()) * centerAtoms.getZBox() - zDist;
		double radius = Math.sqrt(xDist*xDist + yDist*yDist + zDist*zDist);
		return radius;
	}
	
	/**
	 * @return	The volume of the simulation box.
	 */
	private double getSimulationVolume()
	{
		return (centerAtoms.getXBox() * centerAtoms.getYBox() * centerAtoms.getZBox());
	}
	
	/**
	 * Calculates the volume of a specific RDF shell.
	 * @param count			The shell ordered in increased distance from the central atom.
	 * @param increment		The radius of each shell.
	 * @return				The volume of an RDF shell.
	 */
	private static double getShellVolume(int count, double increment)
	{
		return (4.0/3 * 3.14159 * (Math.pow(count*increment, 3) - Math.pow(((count-1)*increment), 3)));
	}
	
	/**
	 * @return	The atomic density of a simulation box.
	 */
	private double getAtomDensity()
	{
		return (otherAtoms.size()/getSimulationVolume());
	}

	public void setCenterAtomSystem(BulkSystem system)
	{
		centerAtoms = system;
	}
}
