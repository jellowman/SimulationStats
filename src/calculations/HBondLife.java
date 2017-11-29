package calculations;

import java.util.ArrayList;
import java.util.LinkedHashMap;

import core.Atom;
import core.BulkSystem;
import core.Molecule;
import lammps.OutParser;

public class HBondLife 
{
	/**Holds a list of expired hydrogen bonds with their bond duration*/
	private ArrayList<Integer> lifetimes;
	
	/**Holds info about the simulation box dimensions*/
	private BulkSystem currentSystem;
	
	/**Holds info of the atoms throughout all timesteps (assumes no new atoms added)*/
	private ArrayList<Atom> atoms;
	
	/**Keeps track live hydrogen bonds from previous steps*/
	private LinkedHashMap<String, Bond> hBonds;
	
	/**The maximum length between two oxygen atoms for a hydrogen bond to be defined*/
	private static final double MAX_H_BOND_LENGTH = 3.5;
	
	/**The maximum angle between some atoms for a hydrogen bond to be defined*/
	private static final double H_BOND_ANGLE = Math.toRadians(30.0);
	
	/**The maximum allowable deviation from the H_BOND_ANGLE*/
	private static final double ANGLE_TOLERANCE = 5.0;
	
	public HBondLife() {
		lifetimes = new ArrayList<Integer>();
		currentSystem = null;
		atoms = new ArrayList<Atom>();
		
	}
	
	private class Bond {
		private Atom atom1;
		private Atom atom2;
		private int lifeTime;
		
		public Bond(Atom a1, Atom a2) {
			lifeTime = 0;
			atom1 = a1;
			atom2 = a2;
		}
		
		public int getLifetime() {
			return lifeTime;
		}
		
		public void incrementLifeTime(int inc) {
			lifeTime = lifeTime + inc;
			}
		
		public void setLifeTime(int lt) {
			lifeTime = lt;
		}
		
		
	}
	
	public void runAnalysis(String filename) {
		//Create parser to read in new atomic coordinates for each timestep
		OutParser timesteps = new OutParser(filename);
		
		//Create lists of atoms and molecule groups
		ArrayList<Molecule> molecules = new ArrayList<Molecule>();
		timesteps.buildSystem(atoms, molecules);
		
		//Create separate list for just Oxygen atoms
		ArrayList<Atom> oxygen = new ArrayList<Atom>();
		for(int i = 0; i < atoms.size(); i++) {
			if(atoms.get(i).getType().equals("O")) {
				oxygen.add(atoms.get(i));
			}
		}
		
		while(true) {
			//Update atom coordinates
			if(timesteps.nextTimestep(atoms)) {
				break;
			}
			
			//Build the h bonds that exist from the new timestep
			LinkedHashMap<String, Bond> newHBonds = new LinkedHashMap<String, Bond>();
			for(int i = 0; i < oxygen.size(); i++) {
				for(int j=0; j < oxygen.size(); j++) {
					if(this.isHBonded(oxygen.get(i), oxygen.get(j))) {
						if(i != j) {
							//Build string ID for the unique type of hydrogen bond
							String key = String.valueOf(oxygen.get(i).getID()) + "-" + String.valueOf(oxygen.get(j).getID());
							newHBonds.put(key, new Bond(oxygen.get(i), oxygen.get(j)));
						}
					}
				}
			}
			
			//Check to see if any h bonds are missing from the previous timestep
			String[] keys = (String[]) hBonds.keySet().toArray();
			for(String key : keys) {
				//If it isn't present in the new timestep, output lifetime
				if(!newHBonds.containsKey(key)) {
					lifetimes.add(hBonds.get(key).getLifetime());
				} else { //Otherwise, increment lifetime in new map
					newHBonds.get(key).setLifeTime(hBonds.get(key).getLifetime()+1);
				}
			}
			
			
			hBonds = newHBonds;
			
		} //END While loop
		
		
	}
	
	private boolean isHBonded(Atom a1, Atom a2) {
		//If the bonding specifications are met
		//First check distance criteria
		if(Geometry.PeriodicDistance(a1, a2, currentSystem.getXBox(), currentSystem.getYBox(),
				currentSystem.getZBox()) > MAX_H_BOND_LENGTH) {
			//Then check angle criteria
			if(Math.abs(Geometry.PeriodicAngle(a1, a2, atoms.get(a1.getID()+1), currentSystem.getXBox(),
					currentSystem.getYBox(), currentSystem.getZBox()) - H_BOND_ANGLE) > ANGLE_TOLERANCE) {
				return true;
			}
		}
		//Otherwise, they are not hydrogen bonded
		return false;
	}
}

/*
 * Hydrogen bond life specifics:
 * First: Set up data structures to hold onto bonds between specific IDs
 * Second: Obtain a simulation file containing box dimensions, molecules, atoms, etc..
 * for a given timestep
 * Third: Find out existing H-Bonds in the simulation
 * Fourth: Determine the new H-Bonds, start tracking in the data structure (LinkedHashMap)
 * Fifth: Determine old H-Bonds, update lifetime
 * Sixth: Determine which H-Bonds have broken, take lifetime and add to statistic list (ArrayList)
 * 
 */
