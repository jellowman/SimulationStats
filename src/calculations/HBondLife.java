package calculations;

import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.Scanner;
import java.util.Set;

import core.Atom;
import core.Bin3D;
import core.BulkSystem;
import core.Molecule;
import lammps.OutParser;

public class HBondLife 
{
	/**Holds a list of expired hydrogen bonds with their bond duration*/
	private ArrayList<Integer> lifetimes, lifetimesG, lifetimesO, lifetimesGO;
	
	/**Holds info about the simulation box dimensions*/
	private BulkSystem currentSystem;
	
	/**Holds info of the atoms throughout all timesteps (assumes no new atoms added)*/
	private ArrayList<Atom> atoms;
	
	/**Holds molecules with atoms associated with the above atoms list*/
	private ArrayList<Molecule> molecules;
	
	/**Keeps track live hydrogen bonds from previous steps*/
	private LinkedHashMap<String, Bond> hBonds, hBondsG, hBondsO, hBondsGO;
	
	/**The maximum length between two oxygen atoms for a hydrogen bond to be defined*/
	private static final double MAX_H_BOND_LENGTH = 3.5;
	
	/**The maximum angle between some atoms for a hydrogen bond to be defined*/
	private static final double H_BOND_ANGLE = Math.toRadians(30.0);
	
	/**The maximum allowable deviation from the H_BOND_ANGLE*/
	private static final double ANGLE_TOLERANCE = 5.0;
	
	public HBondLife() {
		lifetimes = new ArrayList<Integer>();
		lifetimesG = new ArrayList<Integer>();
		lifetimesO = new ArrayList<Integer>();
		lifetimesGO = new ArrayList<Integer>();
		currentSystem = new BulkSystem();
		atoms = new ArrayList<Atom>();
		molecules = new ArrayList<Molecule>();
		hBonds = new LinkedHashMap<String, Bond>();
		hBondsG = new LinkedHashMap<String, Bond>();
		hBondsO = new LinkedHashMap<String, Bond>();
		hBondsGO = new LinkedHashMap<String, Bond>();
		
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
	
	public void runAnalysis() {
		Scanner sc = new Scanner(System.in);
		System.out.println("Provide lammps output file for hydrogen bond analysis");
		String fileName = sc.nextLine();
		runAnalysis(fileName);
	}
	public void runAnalysis(String filename) {
		//Create parser to read in new atomic coordinates for each timestep
		OutParser timesteps = new OutParser(filename);
		
		//Create lists of atoms and molecule groups
		//ArrayList<Molecule> molecules = new ArrayList<Molecule>();
		timesteps.buildSystem(atoms, molecules);
		
		//Create separate list for just Oxygen atoms
		ArrayList<Atom> oxygen = new ArrayList<Atom>();
		for(int i = 0; i < atoms.size(); i++) {
			if(atoms.get(i).getType().equals("1") || atoms.get(i).getType().equals("2")) {
				oxygen.add(atoms.get(i));
			}
		}
		//Make list of just Glycerol oxygen
		//Make list of just Water oxygen
		ArrayList<Atom> glyO = new ArrayList<Atom>();
		ArrayList<Atom> watO = new ArrayList<Atom>();
		for(int i = 0; i < oxygen.size(); i++) {
			if(oxygen.get(i).getType().equals("2")) {
				watO.add(oxygen.get(i));
			} else {
				glyO.add(oxygen.get(i));
			}
		}
		
		boolean lastStep = false;
		while(!lastStep) {
			//Update atom coordinates
			lastStep = timesteps.nextTimestep(atoms, currentSystem);
			
			//Build the h bonds that exist from the new timestep
			LinkedHashMap<String, Bond> newHBonds = new LinkedHashMap<String, Bond>();
			LinkedHashMap<String, Bond> newHBondsG = new LinkedHashMap<String, Bond>();
			LinkedHashMap<String, Bond> newHBondsO = new LinkedHashMap<String, Bond>();
			LinkedHashMap<String, Bond> newHBondsGO = new LinkedHashMap<String, Bond>();
			
			
			
			//Update neighbor list---
			Bin3D<Atom> atomBin = new Bin3D<Atom>(MAX_H_BOND_LENGTH, currentSystem.getXMin(),
					currentSystem.getXMax(), currentSystem.getYMin(), currentSystem.getYMax(),
					currentSystem.getZMin(), currentSystem.getZMax());
			for(Atom atom : atoms) {
				atomBin.addFrom3Space(atom, atom.getXCord(), atom.getYCord(), atom.getZCord());
			}
			
			//Perform H-Bonds for every hydrogen bond
			for(int i = 0; i < oxygen.size(); i++) {
				//Get neighbors to oxygen
				ArrayList<Atom> neigh = atomBin.getNeighbors(oxygen.get(i).getXCord(),
						oxygen.get(i).getYCord(), oxygen.get(i).getZCord());
				for(int j = 0; j < neigh.size(); j++) {
					if(oxygen.get(i).getMolID() != neigh.get(j).getMolID() && 
							(neigh.get(j).getType().equals("1") || neigh.get(j).getType().equals("2"))) {
						if(this.isHBonded(oxygen.get(i), neigh.get(j))) {
							//Build string ID for the unique type of hydrogen bond
							String key = String.valueOf(oxygen.get(i).getID()) + "-" + String.valueOf(neigh.get(j).getID());
							newHBonds.put(key, new Bond(oxygen.get(i), neigh.get(j)));
						}
					}
				}
			}
			
			//Perform for glycerol h-bonds
			for(int i = 0; i < glyO.size(); i++) {
				//Get neighbors to oxygen
				ArrayList<Atom> neigh = atomBin.getNeighbors(glyO.get(i).getXCord(),
						glyO.get(i).getYCord(), glyO.get(i).getZCord());
				for(int j = 0; j < neigh.size(); j++) {
					if(glyO.get(i).getMolID() != neigh.get(j).getMolID()) {
						if(neigh.get(j).getType().equals("1")) {
							if(this.isHBonded(glyO.get(i), neigh.get(j))) {
								//Build string ID for the gly-gly hydrogen bond
								String key = String.valueOf(glyO.get(i).getID()) + "-" + String.valueOf(neigh.get(j).getID());
								newHBondsG.put(key, new Bond(glyO.get(i), neigh.get(j)));
							}
						} else if (neigh.get(j).getType().equals("2")) {
							if(this.isHBonded(glyO.get(i), neigh.get(j))) {
								//Build string ID for the gly-water hydrogen bond
								String key = String.valueOf(glyO.get(i).getID()) + "-" + String.valueOf(neigh.get(j).getID());
								newHBondsGO.put(key, new Bond(glyO.get(i), neigh.get(j)));
							}
						}
					}
				}
			}
			
			//Perform for water h-bonds
			for(int i = 0; i < watO.size(); i++) {
				//Get neighbors to oxygen
				ArrayList<Atom> neigh = atomBin.getNeighbors(watO.get(i).getXCord(),
						watO.get(i).getYCord(), watO.get(i).getZCord());
				for(int j = 0; j < neigh.size(); j++) {
					if(watO.get(i).getMolID() != neigh.get(j).getMolID()) {
						if(neigh.get(j).getType().equals("2")) {
							if(this.isHBonded(watO.get(i), neigh.get(j))) {
								//Build string ID for the water-water hydrogen bond
								String key = String.valueOf(watO.get(i).getID()) + "-" + String.valueOf(neigh.get(j).getID());
								newHBondsO.put(key, new Bond(watO.get(i), neigh.get(j)));
							}
						} else if (neigh.get(j).getType().equals("1")) {
							if(this.isHBonded(watO.get(i), neigh.get(j))) {
								//Build string ID for the water-gly hydrogen bond
								String key = String.valueOf(watO.get(i).getID()) + "-" + String.valueOf(neigh.get(j).getID());
								newHBondsGO.put(key, new Bond(watO.get(i), neigh.get(j)));
							}
						}
					}
				}
			}
			
			
			//Check to see if any h bonds are missing from the previous timestep
			//String[] keys = (String[]) hBonds.keySet().toArray();
			/*Set<String> keys = hBonds.keySet();
			for(String key : keys) {
				//If it isn't present in the new timestep, output lifetime
				if(!newHBonds.containsKey(key)) {
					lifetimes.add(hBonds.get(key).getLifetime());
				} else { //Otherwise, increment lifetime in new map
					newHBonds.get(key).setLifeTime(hBonds.get(key).getLifetime()+1);
				}
			}*/
			
			updateBondsList(hBonds, newHBonds, lifetimes);
			updateBondsList(hBondsG, newHBondsG, lifetimesG);
			updateBondsList(hBondsO, newHBondsO, lifetimesO);
			updateBondsList(hBondsGO, newHBondsGO, lifetimesGO);
			hBonds = newHBonds;
			hBondsG = newHBondsG;
			hBondsO = newHBondsO;
			hBondsGO = newHBondsGO;
			
			
		} //END While loop
		
		System.out.println("ALL H-BONDS");
		printLifetimes(lifetimes, hBonds);
		System.out.println("GLY-GLY H-BONDS");
		printLifetimes(lifetimesG, hBondsG);
		System.out.println("WATER-WATER H-BONDS");
		printLifetimes(lifetimesO, hBondsO);
		System.out.println("GLY-WATER H-BONDS");
		printLifetimes(lifetimesGO, hBondsGO);
		
	}
	
	private boolean isHBonded(Atom a1, Atom a2) {
		//If the bonding specifications are met
		//First check distance criteria
		if(Geometry.PeriodicDistance(a1, a2, currentSystem.getXBox(), currentSystem.getYBox(),
				currentSystem.getZBox()) < MAX_H_BOND_LENGTH) {
			//Then check angle criteria
			if(Math.abs(Geometry.PeriodicAngle(a1, a2, atoms.get(a1.getID()+1), currentSystem.getXBox(),
					currentSystem.getYBox(), currentSystem.getZBox()) - H_BOND_ANGLE) > ANGLE_TOLERANCE) {
				return true;
			}
		}
		//Otherwise, they are not hydrogen bonded
		return false;
	}
	
	/*private void calcHBonds(ArrayList<Atom> oGroup, LinkedHashMap<String, Bond> newHBonds, Bin3D<Atom> bin) {
		for(int i = 0; i < oGroup.size(); i++) {
			//Get neighbors to oxygen
			ArrayList<Atom> neigh = bin.getNeighbors(oGroup.get(i).getXCord(),
					oGroup.get(i).getYCord(), oGroup.get(i).getZCord());
			for(int j = 0; j < neigh.size(); j++) {
				if(oGroup.get(i).getMolID() != neigh.get(j).getMolID() && 
						(neigh.get(j).getType().equals("1") || neigh.get(j).getType().equals("2"))) {
					if(this.isHBonded(oGroup.get(i), neigh.get(j))) {
						//Build string ID for the unique type of hydrogen bond
						String key = String.valueOf(oGroup.get(i).getID()) + "-" + String.valueOf(neigh.get(j).getID());
						newHBonds.put(key, new Bond(oGroup.get(i), neigh.get(j)));
					}
				}
			}
		}
	}*/
	
	private static void updateBondsList(LinkedHashMap<String, Bond> oldBonds, LinkedHashMap<String, Bond> newBonds, ArrayList<Integer> lifetimes) {
		Set<String> keys = oldBonds.keySet();
		for(String key : keys) {
			//If it isn't present in the new timestep, output lifetime
			if(!newBonds.containsKey(key)) {
				lifetimes.add(oldBonds.get(key).getLifetime());
			} else { //Otherwise, increment lifetime in new map
				newBonds.get(key).setLifeTime(oldBonds.get(key).getLifetime()+1);
			}
		}
		
	}
	
	private static void printLifetimes(ArrayList<Integer> lifetimes, LinkedHashMap<String, Bond> hBonds) {
		//Find largest size
		int largest = 0;
		for(int i = 0; i < lifetimes.size(); i++) {
			if (lifetimes.get(i) > largest) {
				largest = lifetimes.get(i);
			}
		}
		for(Bond bond : hBonds.values()) {
			if(bond.getLifetime() > largest) {
				largest = bond.getLifetime();
			}
		}
		int[] times = new int[largest+1];
		for(int i = 0; i < times.length; i++) {
			times[i] = 0;
		}
		System.out.println("List of lifetimes:");
		for(int i = 0; i < lifetimes.size(); i++) {
			times[lifetimes.get(i)] += 1;
		}
		for(int i = 0; i < times.length; i++) {
			System.out.println("[" + i + "] - " + times[i]);
			times[i] = 0;
		}
		System.out.println("Still bonded:");
		for(Bond bond : hBonds.values()) {
			times[bond.getLifetime()] += 1;
		}
		for(int i = 0; i < times.length; i++) {
			System.out.println("[" + i + "] - " + times[i]);
		}
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
