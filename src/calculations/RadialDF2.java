package calculations;

import java.lang.reflect.Array;
import java.text.DecimalFormat;
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

public class RadialDF2 
{	
	/**Holds info about the simulation box dimensions*/
	protected BulkSystem currentSystem;
	
	/**Holds info of the atoms throughout all timesteps (assumes no new atoms added)*/
	protected ArrayList<Atom> atoms;
	
	/**Holds molecules with atoms associated with the above atoms list*/
	protected ArrayList<Molecule> molecules;
	
	/**Interval between each RDF bin**/
	protected static final double INTERVAL = 0.05;
	
	/**Cutoff for RDF**/
	protected double cutoff;
	protected String centerAtom;
	protected String otherAtom;
	
	public RadialDF2(double cutoff, String centerAtom, String otherAtom) {
		currentSystem = new BulkSystem();
		atoms = new ArrayList<Atom>();
		molecules = new ArrayList<Molecule>();
		this.centerAtom = centerAtom;
		this.otherAtom = otherAtom;
		this.cutoff = cutoff;
	}
	
	public void runAnalysis() {
		Scanner sc = new Scanner(System.in);
		System.out.println("Provide lammps output file for RDF analysis");
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
		ArrayList<Atom> centerAtoms = new ArrayList<Atom>();
		ArrayList<Atom> otherAtoms = new ArrayList<Atom>();
		
		for(int i = 0; i < atoms.size(); i++) {
			if(atoms.get(i).getType().equals(centerAtom)) {
				centerAtoms.add(atoms.get(i));
			}
			if(atoms.get(i).getType().equals(otherAtom)) {
				otherAtoms.add(atoms.get(i));
			}
		}
		
		int size = (int) Math.floor(cutoff/INTERVAL);
		int[] distances = new int[size];
		Arrays.fill(distances, 0);
		
		boolean lastStep = false;
		int timestep = 0;
		while(!lastStep) {
			//Update atom coordinates and box dimensions
			lastStep = timesteps.nextTimestep(atoms, currentSystem);
			
			//Update neighbor list---
			Bin3D<Atom> atomBin = new Bin3D<Atom>(cutoff, currentSystem.getXMin(),
					currentSystem.getXMax(), currentSystem.getYMin(), currentSystem.getYMax(),
					currentSystem.getZMin(), currentSystem.getZMax());
			
			for(Atom atom : centerAtoms) {
				atomBin.addFrom3Space(atom, atom.getXCord(), atom.getYCord(), atom.getZCord());
			}
			if(!centerAtom.equals(otherAtom)) {
				for(Atom atom : otherAtoms) {
					atomBin.addFrom3Space(atom, atom.getXCord(), atom.getYCord(), atom.getZCord());
				}
			}
			
			//Perform half-volume RDF to account for glycerol backbone interference
			if(otherAtom.equals("1")) {
				for(int i = 0; i < centerAtoms.size(); i++) {
					//Get neighbors to oxygen
					ArrayList<Atom> neigh = atomBin.getNeighbors(centerAtoms.get(i).getXCord(),
							centerAtoms.get(i).getYCord(), centerAtoms.get(i).getZCord());
					//Determine Carbon atom that glycerol oxygen is bonded to
					int repeatGly = centerAtoms.get(i).getID() % 14;
					Atom carbonAtom = null;
					if(repeatGly == 8) carbonAtom = atoms.get(centerAtoms.get(i).getID()-8);
					else if(repeatGly == 10) carbonAtom = atoms.get(centerAtoms.get(i).getID()-7);
					else if(repeatGly == 12) carbonAtom = atoms.get(centerAtoms.get(i).getID()-7);
					for(int j = 0; j < neigh.size(); j++) {
						if(neigh.get(j).getType().equals(otherAtom)) {
							//if(centerAtoms.get(i).getID() != neigh.get(j).getID()) {
							//if(Geometry.PeriodicAngle(centerAtoms.get(i), neigh.get(j), carbonAtom, currentSystem.getXBox(), currentSystem.getYBox(), currentSystem.getZBox()) > 90) {
								if(centerAtoms.get(i).getMolID() != neigh.get(j).getMolID()) {
									double r = Geometry.PeriodicDistance(centerAtoms.get(i), neigh.get(j), currentSystem.getXBox(),
											currentSystem.getYBox(), currentSystem.getZBox());
									int radiusLayer = (int)(r/INTERVAL);
									if(radiusLayer < distances.length) {
										//if(Geometry.PeriodicAngle(centerAtoms.get(i), neigh.get(j), carbonAtom, currentSystem.getXBox(), currentSystem.getYBox(), currentSystem.getZBox()) > 90) {
											if(Geometry.PeriodicAngle(neigh.get(j), atoms.get(neigh.get(j).getID()+1), centerAtoms.get(i), currentSystem.getXBox(), currentSystem.getYBox(), currentSystem.getZBox()) < 50) {
												distances[radiusLayer] += 1;
											}
										//}
									}
								}
							//}
						}
					}
				}
			} else { //Perform normal RDF calculation
				for(int i = 0; i < centerAtoms.size(); i++) {
					//Get neighbors to oxygen
					ArrayList<Atom> neigh = atomBin.getNeighbors(centerAtoms.get(i).getXCord(),
							centerAtoms.get(i).getYCord(), centerAtoms.get(i).getZCord());
					for(int j = 0; j < neigh.size(); j++) {
						if(neigh.get(j).getType().equals(otherAtom)) {
							//if(centerAtoms.get(i).getID() != neigh.get(j).getID()) {
							if(centerAtoms.get(i).getMolID() != neigh.get(j).getMolID()) {
								double r = Geometry.PeriodicDistance(centerAtoms.get(i), neigh.get(j), currentSystem.getXBox(),
										currentSystem.getYBox(), currentSystem.getZBox());
								int radiusLayer = (int)(r/INTERVAL);
								if(radiusLayer < distances.length) {
									if(Geometry.PeriodicAngle(neigh.get(j), centerAtoms.get(i), atoms.get(neigh.get(j).getID()+1), currentSystem.getXBox(), currentSystem.getYBox(), currentSystem.getZBox()) < 35
											|| Geometry.PeriodicAngle(neigh.get(j), centerAtoms.get(i), atoms.get(neigh.get(j).getID()+2), currentSystem.getXBox(), currentSystem.getYBox(), currentSystem.getZBox()) < 35) {
										distances[radiusLayer] += 1;
									}
								}
							}
						}
					}
				}
			}
			timestep += 1;
		} //END While loop
		
		//Normalize distances
		double[] normalDistances = new double[distances.length];
		double cumulative = 0.0;
		double aDensity = otherAtoms.size()/currentSystem.getXBox()/currentSystem.getYBox()/currentSystem.getZBox();
		System.out.println("Layer: Relative Density - Cumulative number of neighbors");
		for(int i = 0; i < normalDistances.length; i++) {
			double volume = Geometry.getShellVolume(i+1, INTERVAL);
			if(centerAtom.equals("1")) {
				volume = volume / 2.0;
			}
			normalDistances[i] = distances[i]/volume/aDensity/centerAtoms.size()/timestep;
			cumulative += (double)(distances[i])/centerAtoms.size()/timestep;
			DecimalFormat df = new DecimalFormat("#.##");
			System.out.println("Layer: " + df.format(i*INTERVAL) + " - " + normalDistances[i] + " - " + cumulative);
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
