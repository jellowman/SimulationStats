package calculations;

import java.lang.reflect.Array;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.Scanner;
import java.util.Set;

import core.Atom;
import core.Bin3D;
import core.BulkSystem;
import core.Molecule;
import lammps.OutParser;

/**Calculates the dynamic H-Bond distributions Pc(t) and Ptt(t)*/
public class HBondLifeTC 
{	
	/**Holds info about the simulation box dimensions*/
	protected BulkSystem currentSystem;
	
	/**Holds info of the atoms throughout all timesteps (assumes no new atoms added)*/
	protected ArrayList<Atom> atoms;
	
	/**Holds molecules with atoms associated with the above atoms list*/
	protected ArrayList<Molecule> molecules;
	
	/**Keeps track live hydrogen bonds from previous steps*/
	protected LinkedHashMap<String, HBond> hBonds;
	
	/**Holds all H-Bonds that have ended*/
	protected ArrayList<HBond> oldHBonds;
	
	/**The maximum length between two oxygen atoms for a hydrogen bond to be defined*/
	protected static final double MAX_H_BOND_LENGTH = 3.8;
	
	/**Extra distance outside of max H-Bond range for proximity population*/
	 protected static final double EXTRA_DISTANCE = 0.0;
	
	/**The minimum angle between some atoms for a hydrogen bond to be defined*/
	protected static final double H_BOND_ANGLE = 140.0;
	
	//**The maximum allowable deviation from the H_BOND_ANGLE*/
	//private static final double ANGLE_TOLERANCE = 5.0;
	
	/**The maximum number of timesteps where an H-Bond can violate criteria*/
	protected static final int H_BOND_TIME_ALLOWANCE = 3;
	
	/**The timestep between each trajectory*/
	protected static final double TIMESTEP = 20;
	
	protected static final String DONOR = "1";
	protected static final String ACCEPTOR = "1";
	//2 is water, 1 is glycerol
	
	public HBondLifeTC() {
		currentSystem = new BulkSystem();
		atoms = new ArrayList<Atom>();
		molecules = new ArrayList<Molecule>();
		hBonds = new LinkedHashMap<String, HBond>();
		oldHBonds = new ArrayList<HBond>();
	}
	
	protected class HBond {
		private Atom acceptor;
		private Atom donor;
		private Atom proton;
		private int lifeStart;
		private int lifeEnd;
		private int timesViolated;
		private boolean isVisited;
		private int numViolations;
		
		/**
		 * 
		 * @param a1 The atom accepting an H in the H-Bond
		 * @param a2 The atom donating an H in the H-Bond
		 * @param start Start of H-Bond Lifetime
		 */
		public HBond(Atom a1, Atom a2, int start) {
			lifeStart = start;
			lifeEnd = start;
			acceptor = a1;
			donor = a2;
			proton = null;
			timesViolated = 0;
			isVisited = true;
			numViolations = 0;
		}
		
		/**
		 * 
		 * @param acceptor The atom accepting an H in the H-Bond
		 * @param donor The atom donating an H in the H-bond
		 * @param proton The H-atom in the H-Bond
		 * @param start Start of H-Bond Lifetime
		 */
		public HBond(Atom acceptor, Atom donor, Atom proton, int start) {
			this(acceptor, donor, start);
			this.proton = proton;
		}
		
		public Atom getAcceptorAtom() {
			return acceptor;
		}
		
		public Atom getDonatorAtom() {
			return donor;
		}
		
		public Atom getProtonAtom() {
			return proton;
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
		
		public int getNumViolations() {
			return numViolations;
		}
		
		public void incrementTotalViolations() {
			numViolations += 1;
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
		ArrayList<Integer> violationTimes = new ArrayList<Integer>();
		double avgViolationTime = 0;
		int violationCounter = 0;
		ArrayList<Double> bR = new ArrayList<Double>();
		ArrayList<Double> fR = new ArrayList<Double>();
		
		boolean lastStep = false;
		int timestep = 0;
		while(!lastStep) {
			int numFormed = 0;
			int proximityPop = 0;
			int numBroken = 0;
			int bondedPop = 0;
			//Update atom coordinates
			lastStep = timesteps.nextTimestep(atoms, currentSystem);
			
			//Update neighbor list---
			Bin3D<Atom> atomBin = new Bin3D<Atom>(MAX_H_BOND_LENGTH+EXTRA_DISTANCE, currentSystem.getXMin(),
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
							if(oxygen.get(i).getType().equals(ACCEPTOR) && neigh.get(j).getType().equals(DONOR)) {
								//Increment number of H-Bonds
								bondedPop += 1;
							}
							//Build string ID for the unique type of hydrogen bond
							String key = String.valueOf(oxygen.get(i).getID()) + "-" + String.valueOf(neigh.get(j).getID());
							//newHBonds.put(key, new Bond(oxygen.get(i), neigh.get(j)));
							HBond check = hBonds.putIfAbsent(key, new HBond(oxygen.get(i), neigh.get(j), timestep));
							//If not null, H-Bond Existed last timestep, so update endLifetime for existing HBond and put back in
							if(check != null) {
								//hBonds.get(key).setEndLifetime(timestep);
								check.setEndLifetime(timestep);
								if(check.getTimesViolated() != 0) {
									if(check.getAcceptorAtom().getType().equals(ACCEPTOR) && check.getDonatorAtom().getType().equals(DONOR)) {
										violationTimes.add(new Integer(check.getTimesViolated()));
										avgViolationTime += check.getTimesViolated();
										violationCounter += 1;
										check.incrementTotalViolations();
									}
								}
								check.resetViolations();
								check.setVisited(true);
								//hBonds.put(key, check);
							} else {
								//H-Bond is not an existing one, forms new bond
								numFormed += 1;
							}
						} else if(Geometry.PeriodicDistance(oxygen.get(i), neigh.get(j), currentSystem.getXBox(),
								currentSystem.getYBox(), currentSystem.getZBox()) < (MAX_H_BOND_LENGTH + EXTRA_DISTANCE)) {
							//Oxygen are not oriented for a H-Bond, but they are within proximity of each other
							//Also not in temporary violation (Can't be in existing H-Bond population)
							//String key = String.valueOf(oxygen.get(i).getID()) + "-" + String.valueOf(neigh.get(j).getID());
							//HBond check = hBonds.putIfAbsent(key, new HBond(oxygen.get(i), neigh.get(j), timestep));
							//if(check == null) { //If H-Bond didn't exist last timestep, then check is null
								if(oxygen.get(i).getType().equals(ACCEPTOR) && neigh.get(j).getType().equals(DONOR)) {
									proximityPop += 1;
								}
							//}
						}
					}
				}
			}
			
			//Expire H-Bonds that have violated H-Bond criteria for too long, these must be assigned to
			//previous timestep corresponding to the moment they broke (H_BOND_ALLOWANCE earlier)
			numBroken = checkVisitedHBonds(hBonds, oldHBonds);
			//checkVisitedHBonds(hBonds,oldHBonds);
			//Print out H-Bond formation and break rate constants
			System.out.println(numBroken);
			System.out.println(numFormed);
			//bondedPop = hBonds.size();
			double breakRate = (double)numBroken/bondedPop;
			double formationRate = (double)numFormed/proximityPop;
			if(breakRate != 0) {
			bR.add(breakRate);
			fR.add(formationRate);
			}
			//System.out.println("Break rate: " + breakRate + " per " + TIMESTEP + " picoseconds.");
			//System.out.println("Formation rate: " + formationRate + " per " + TIMESTEP + " picoseconds.");
			timestep += 1;
		} //END While loop
		
		double avgBR = 0; double avgFR = 0;
		for(int i = 0; i < bR.size(); i++) {
			avgBR += bR.get(i);
			avgFR += fR.get(i);
		}
		avgBR = avgBR / bR.size(); avgFR = avgFR / fR.size();
		
		System.out.println("Average Break Rate: " + avgBR + " per " + TIMESTEP + " picoseconds.");
		System.out.println("Average Formation Rate: " + avgFR + " per " + TIMESTEP + " picoseconds.");
		printLifetimes(hBonds, oldHBonds, timestep-1);
		
		if(H_BOND_TIME_ALLOWANCE > 0){
		avgViolationTime = avgViolationTime/violationCounter;
		Collections.sort(violationTimes);
		System.out.println("Average time H-Bonds Violated: " + avgViolationTime*TIMESTEP);
		System.out.println("Median time H-Bonds Violated: " + violationTimes.get(violationTimes.size()/2)*TIMESTEP);
		
		int[] vioBin = new int[H_BOND_TIME_ALLOWANCE+1];
		Arrays.fill(vioBin, 0);
		for(Integer vio : violationTimes) {
			vioBin[vio] += 1;
		}
		for(int i = 0; i < vioBin.length; i++) {
			System.out.println("Bin " + (i) + " - " + vioBin[i]);
		}
		}
		
		int numActive = 0;
		double activeMean = 0;
		int numWholeTime = 0;
		Set<String> activeBondKeys = hBonds.keySet();
		for(String key : activeBondKeys) {
			numActive += 1;
			activeMean = activeMean + (timestep - 1 - hBonds.get(key).getLifeStart());
			if(hBonds.get(key).getLifeStart() == 0) {
				numWholeTime += 1;
			}
		}
		activeMean = activeMean / numActive;
		
		System.out.println("Number of active bonds: " + numActive);
		System.out.println("Active bond mean lifetime: " + activeMean*TIMESTEP);
		System.out.println("Number of bonds lasting whole time: " + numWholeTime);
		
	}
	
	/**
	 * 
	 * @param a1 The atom accepting a H
	 * @param a2 The atom donating a H
	 * @return
	 */
	protected boolean isHBonded(Atom a1, Atom a2) {
		//If the bonding specifications are met
		//First check distance criteria
		if(Geometry.PeriodicDistance(a1, a2, currentSystem.getXBox(), currentSystem.getYBox(),
				currentSystem.getZBox()) < MAX_H_BOND_LENGTH) {
			//Then check angle criteria
			double angle;
			boolean result = false;
			if(a2.getType().equals("2")) { //Donator atom is water, must check both H on O
				angle = Geometry.PeriodicAngle(atoms.get(a2.getID()+1), a2, a1, currentSystem.getXBox(),
						currentSystem.getYBox(), currentSystem.getZBox());
				if(Math.abs(angle) > H_BOND_ANGLE) {
					result = true;
				} else {
					angle = Geometry.PeriodicAngle(atoms.get(a2.getID()+2), a2, a1, currentSystem.getXBox(),
							currentSystem.getYBox(), currentSystem.getZBox());
					if(Math.abs(angle) > H_BOND_ANGLE) {
						result = true;
					}
				}
			} else {
			//TODO add if statement to check both H on the water O-type
				angle = Geometry.PeriodicAngle(atoms.get(a2.getID()+1), a2, a1, currentSystem.getXBox(),
						currentSystem.getYBox(), currentSystem.getZBox());
				if(Math.abs(angle) > H_BOND_ANGLE) {
					result = true;
				}
			}
			return result;
		}
		//Otherwise, they are not hydrogen bonded
		return false;
	}
	
	protected static void resetVisitedHBonds(LinkedHashMap<String, HBond> hBonds) {
		Set<String> keys = hBonds.keySet();
		for(String key : keys) {
			hBonds.get(key).setVisited(false);
		}
	}
	protected static int checkVisitedHBonds(LinkedHashMap<String, HBond> hBonds, ArrayList<HBond> oldHBonds) {
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
		int numBroken = 0;
		while(keyIterator.hasNext()) {
			String nextKey = keyIterator.next();
			if(!hBonds.get(nextKey).isVisited()) {
				hBonds.get(nextKey).incrementViolations();
				if(hBonds.get(nextKey).getTimesViolated() > H_BOND_TIME_ALLOWANCE) {
					//Retire H-Bond
					oldHBonds.add(hBonds.get(nextKey));
					keyIterator.remove();
					numBroken += 1;
				}
			}
		}
		return numBroken;
	}
	
	protected static void printLifetimes(LinkedHashMap<String, HBond> hBonds, ArrayList<HBond> retiredBonds, int timesteps) {
		int middleTimestep = timesteps/2;
		System.out.println("Total Timesteps: " + timesteps);
		//Determine H-Bonds that existed during middle timestep and determine average
		//double dz = 0.5/TIMESTEP;
		//int histSize = (int)(timesteps/dz)+1;
		int histSize = timesteps/10;
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
		double avgNumViolations = 0;
		for(HBond bond : retiredBonds) {
			//if((bond.getAcceptorAtom().getID() % 14) != 10 && (bond.getDonatorAtom().getID() % 14) != 10) {
			if((bond.getAcceptorAtom().getType().equals(ACCEPTOR)) && (bond.getDonatorAtom().getType().equals(DONOR))) {
				if((bond.getLifeStart() <= middleTimestep) && (bond.getLifeEnd() >= middleTimestep)) {
					if(bond.getLifeStart() != 0) {
						sumTc += bond.getLifetime();
						numHBondsTc += 1;
						hist[((int)(bond.getLifetime()/dz))] += 1;
						avgNumViolations += bond.getNumViolations();
						//System.out.println("Lifetime: " + bond.getLifetime() + " - " + bond.getNumViolations() + " out of criteria.");
					} else {
						//Increment error counter for HBonds excluded
					}
				}
				if(bond.getLifeStart() != 0) {
					sumTt += bond.getLifetime();
					numHBondsTt += 1;
					hist2[((int)(bond.getLifetime()/dz))] += 1;
					//System.out.println("Lifetime: " + bond.getLifetime() + " - " + bond.getNumViolations() + " out of criteria.");
					//avgNumViolations += bond.getNumViolations();
				}
			}
			//}
		}
		//avgNumViolations = avgNumViolations / numHBondsTt;
		avgNumViolations = avgNumViolations / numHBondsTc;
		double avgHBondLifetimeTc = (double)sumTc/numHBondsTc;
		double avgHBondLifetimeTt = (double)sumTt/numHBondsTt;
		System.out.println("Tc Lifetime: " + (avgHBondLifetimeTc*TIMESTEP));
		System.out.println("Avg. Hbond from Tc: " + (avgHBondLifetimeTc*TIMESTEP/2));
		System.out.println("Tt Lifetime: " + (avgHBondLifetimeTt*TIMESTEP));
		System.out.println("Bonds had " + avgNumViolations + " violations on average.");
		
		System.out.println("\nConfiguration Averaged Bins\n");
		for(int i=0; i < hist.length; i++) {
			//System.out.println("Bin " + i + " - " + hist[i]);
			System.out.println(hist[i]);
		}
		System.out.println("\nTrajectory Averaged Bins\n");
		for(int i=0; i < hist2.length; i++) {
			//System.out.println("Bin " + i + " - " + hist2[i]);
			System.out.println(hist2[i]);
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
