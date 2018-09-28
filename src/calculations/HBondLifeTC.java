package calculations;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.Scanner;
import java.util.Set;

import core.Atom;
import core.Bin3D;
import core.BulkSystem;
import core.Molecule;
import lammps.OutParser;

public class HBondLifeTC 
{	
	/**Holds info about the simulation box dimensions*/
	private BulkSystem currentSystem;
	
	/**Holds info of the atoms throughout all timesteps (assumes no new atoms added)*/
	private ArrayList<Atom> atoms;
	
	/**Holds molecules with atoms associated with the above atoms list*/
	private ArrayList<Molecule> molecules;
	
	/**Keeps track live hydrogen bonds from previous steps*/
	private LinkedHashMap<String, HBond> hBonds;
	
	/**Holds all H-Bonds that have ended*/
	private ArrayList<HBond> oldHBonds;
	
	/**The maximum length between two oxygen atoms for a hydrogen bond to be defined*/
	private static final double MAX_H_BOND_LENGTH = 4.2;
	
	/**The maximum angle between some atoms for a hydrogen bond to be defined*/
	private static final double H_BOND_ANGLE = 30.0;
	
	/**The maximum allowable deviation from the H_BOND_ANGLE*/
	private static final double ANGLE_TOLERANCE = 5.0;
	
	/**The maximum number of timesteps where an H-Bond can violate criteria*/
	private static final int H_BOND_TIME_ALLOWANCE = 10;
	
	public HBondLifeTC() {
		currentSystem = new BulkSystem();
		atoms = new ArrayList<Atom>();
		molecules = new ArrayList<Molecule>();
		hBonds = new LinkedHashMap<String, HBond>();
		oldHBonds = new ArrayList<HBond>();
	}
	
	private class HBond {
		private Atom atom1;
		private Atom atom2;
		private int lifeStart;
		private int lifeEnd;
		private int timesViolated;
		private boolean isVisited;
		
		public HBond(Atom a1, Atom a2, int start) {
			lifeStart = start;
			lifeEnd = start;
			atom1 = a1;
			atom2 = a2;
			timesViolated = 0;
			isVisited = true;
		}
		
		public Atom getAcceptorAtom() {
			return atom1;
		}
		
		public Atom getDonatorAtom() {
			return atom2;
		}
		public int getLifetime() {
			return (lifeEnd - lifeStart);
		}
		
		public int getLifeStart() {
			return lifeStart;
		}
		
		public int getLifeEnd()	 {
			return lifeEnd;
		}
		
		public void incrementLifeTime(int inc) {
			lifeEnd = lifeEnd + inc;
			}
		
		public void setEndLifetime(int end)	 {
			lifeEnd = end;
		}
		
		public void resetViolations() {
			timesViolated = 0;
		}
		
		public void incrementViolations() {
			timesViolated += 1;
		}
		
		public int getTimesViolated() {
			return timesViolated;
		}
		
		public void setVisited(boolean visit) {
			isVisited = visit;
		}
		
		public boolean isVisited() {
			return isVisited;
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
		int timestep = 0;
		while(!lastStep) {
			//Update atom coordinates
			lastStep = timesteps.nextTimestep(atoms, currentSystem);
			
			//Update neighbor list---
			Bin3D<Atom> atomBin = new Bin3D<Atom>(MAX_H_BOND_LENGTH, currentSystem.getXMin(),
					currentSystem.getXMax(), currentSystem.getYMin(), currentSystem.getYMax(),
					currentSystem.getZMin(), currentSystem.getZMax());
			for(Atom atom : atoms) {
				atomBin.addFrom3Space(atom, atom.getXCord(), atom.getYCord(), atom.getZCord());
			}
			
			//Set visited flags to false for every existing H-Bond
			resetVisitedHBonds(hBonds);
			
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
							//newHBonds.put(key, new Bond(oxygen.get(i), neigh.get(j)));
							HBond check = hBonds.putIfAbsent(key, new HBond(oxygen.get(i), neigh.get(j), timestep));
							//If not null, H-Bond Existed last timestep, so update endLifetime for existing HBond and put back in
							if(check != null) {
								//hBonds.get(key).setEndLifetime(timestep);
								check.setEndLifetime(timestep);
								check.resetViolations();
								check.setVisited(true);
								//hBonds.put(key, check);
							}
						}
					}
				}
			}
			
			//Expire H-Bonds that have violated H-Bond criteria for too long
			checkVisitedHBonds(hBonds, oldHBonds);
			
			timestep += 1;
		} //END While loop
		
		System.out.println("ALL H-BONDS");
		printLifetimes(oldHBonds, timestep);
		
		
	}
	
	private boolean isHBonded(Atom a1, Atom a2) {
		//If the bonding specifications are met
		//First check distance criteria
		if(Geometry.PeriodicDistance(a1, a2, currentSystem.getXBox(), currentSystem.getYBox(),
				currentSystem.getZBox()) < MAX_H_BOND_LENGTH) {
			//Then check angle criteria
			double angle;
			boolean result = false;
			/*if(a2.getType().equals("2")) {
				angle = Geometry.PeriodicAngle(a2, a1, atoms.get(a2.getID()+1), currentSystem.getXBox(),
						currentSystem.getYBox(), currentSystem.getZBox());
				if(Math.abs(angle) < H_BOND_ANGLE) {
					result = true;
				} else {
					angle = Geometry.PeriodicAngle(a2, a1, atoms.get(a2.getID()+2), currentSystem.getXBox(),
							currentSystem.getYBox(), currentSystem.getZBox());
					if(Math.abs(angle) < H_BOND_ANGLE) {
						result = true;
					}
				}
			} else {*/
				angle = Geometry.PeriodicAngle(a2, a1, atoms.get(a2.getID()+1), currentSystem.getXBox(),
						currentSystem.getYBox(), currentSystem.getZBox());
				if(Math.abs(angle) < H_BOND_ANGLE) {
					result = true;
				}
			//}
			return result;
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
	
	private static void resetVisitedHBonds(LinkedHashMap<String, HBond> hBonds) {
		Set<String> keys = hBonds.keySet();
		for(String key : keys) {
			hBonds.get(key).setVisited(false);
		}
	}
	private static void checkVisitedHBonds(LinkedHashMap<String, HBond> hBonds, ArrayList<HBond> oldHBonds) {
		Iterator<String> keyIterator = hBonds.keySet().iterator();
		/*for(String key : keys) {
			//If not visited, increment violation allowance
			if(!hBonds.get(key).isVisited()) {
				hBonds.get(key).incrementViolations();
				if(hBonds.get(key).getTimesViolated() > H_BOND_TIME_ALLOWANCE) {
					//Retire H-Bond
					oldHBonds.add(hBonds.get(key));
					keys.remove(key);
				}
			}
		}*/
		while(keyIterator.hasNext()) {
			String nextKey = keyIterator.next();
			if(!hBonds.get(nextKey).isVisited()) {
				hBonds.get(nextKey).incrementViolations();
				if(hBonds.get(nextKey).getTimesViolated() > H_BOND_TIME_ALLOWANCE) {
					//Retire H-Bond
					oldHBonds.add(hBonds.get(nextKey));
					keyIterator.remove();
				}
			}
		}
	}
	
	private static void printLifetimes(ArrayList<HBond> retiredBonds, int timesteps) {
		int middleTimestep = timesteps/2;
		System.out.println("Total Timesteps: " + timesteps);
		//Determine H-Bonds that existed during middle timestep and determine average
		int histSize = timesteps/5;
		int[] hist = new int[histSize];
		int[] hist2 = new int[histSize];
		for(int i = 0; i < hist.length; i++) {
			hist[i] = 0;
			hist2[i] = 0;
		}
		double dz = (double)timesteps/histSize;
		int sumTc = 0;
		int numHBondsTc = 0;
		int sumTt = 0;
		int numHBondsTt = 0;
		for(HBond bond : retiredBonds) {
			if((bond.getAcceptorAtom().getType().equals("1")) && (bond.getDonatorAtom().getType().equals("1"))) {
				if((bond.getLifeStart() < middleTimestep) && (bond.getLifeEnd() > middleTimestep)) {
					if(bond.getLifeStart() != 0) {
						sumTc += bond.getLifetime();
						numHBondsTc += 1;
						hist[((int)(bond.getLifetime()/dz))] += 1;
					}
				}
				if(bond.getLifeStart() != 0) {
					sumTt += bond.getLifetime();
					numHBondsTt += 1;
					hist2[((int)(bond.getLifetime()/dz))] += 1;
				}
			}
		}
		
		double avgHBondLifetimeTc = (double)sumTc/numHBondsTc;
		double avgHBondLifetimeTt = (double)sumTt/numHBondsTt;
		System.out.println("Tc Average Lifetime: " + avgHBondLifetimeTc);
		System.out.println("Tt Average Lifetime: " + avgHBondLifetimeTt);
		
		System.out.println("\nConfiguration Averaged Bins\n");
		for(int i=0; i < hist.length; i++) {
			System.out.println("Bin " + i + " - " + hist[i]);
		}
		System.out.println("\nTrajectory Averaged Bins\n");
		for(int i=0; i < hist2.length; i++) {
			System.out.println("Bin " + i + " - " + hist2[i]);
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
