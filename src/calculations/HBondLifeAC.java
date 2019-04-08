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

/**Calculates the standard H-Bond autocorrelation function c(t)*/
public class HBondLifeAC extends HBondLifeTC
{		
	//2 is water, 1 is glycerol
	private String donor;
	private String acceptor;
	private double maxHBondLength;
	private double angleCutoff;
	private double timeInterval;
	
	public HBondLifeAC(String donor, String acceptor, double maxHBondLength, double angleCutoff, double timeInterval) {
		super();
		this.donor = donor;
		this.acceptor = acceptor;
		this.maxHBondLength = maxHBondLength;
		this.angleCutoff = angleCutoff;
		this.timeInterval = timeInterval;
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
		
		ArrayList<Integer> violationTimes = new ArrayList<Integer>();
		double avgViolationTime = 0;
		int violationCounter = 0;
		
		boolean lastStep = false;
		int timestep = 0;
		//Determine H-Bond population at t=0
		lastStep = timesteps.nextTimestep(atoms, currentSystem);
		
		//Update neighbor list---
		Bin3D<Atom> atomBin = new Bin3D<Atom>(maxHBondLength, currentSystem.getXMin(),
				currentSystem.getXMax(), currentSystem.getYMin(), currentSystem.getYMax(),
				currentSystem.getZMin(), currentSystem.getZMax());
		
		for (Atom atom : oxygen) {
			atomBin.addFrom3Space(atom, atom.getXCord(), atom.getYCord(), atom.getZCord());
		}
		
		ArrayList<HBond> initHBonds = new ArrayList<HBond>();
		
		//Determine initial population of H-Bonds
		for(int i = 0; i < oxygen.size(); i++) {
			//Get neighbors to oxygen
			ArrayList<Atom> neigh = atomBin.getNeighbors(oxygen.get(i).getXCord(),
						oxygen.get(i).getYCord(), oxygen.get(i).getZCord());
			for(int j = 0; j < neigh.size(); j++) {
				if(oxygen.get(i).getMolID() != neigh.get(j).getMolID()) {
					if(oxygen.get(i).getType().equals(acceptor) && neigh.get(j).getType().equals(donor)) {
						if(donor.equals("2")) {
							if(this.isHBonded(oxygen.get(i), neigh.get(j), atoms.get(neigh.get(j).getID()+1))) {
								initHBonds.add(new HBond(oxygen.get(i), neigh.get(j), atoms.get(neigh.get(j).getID()+1), timestep));
							}
							if(this.isHBonded(oxygen.get(i), neigh.get(j), atoms.get(neigh.get(j).getID()+2))) {
								initHBonds.add(new HBond(oxygen.get(i), neigh.get(j), atoms.get(neigh.get(j).getID()+2), timestep));
							}
						} else {
							if(this.isHBonded(oxygen.get(i), neigh.get(j), atoms.get(neigh.get(j).getID()+1))) {
								initHBonds.add(new HBond(oxygen.get(i), neigh.get(j), atoms.get(neigh.get(j).getID()+1), timestep));
							}
						}
					}
				}
			}
		}
		
		ArrayList<Integer> cT = new ArrayList<Integer>();
		ArrayList<Integer> nT = new ArrayList<Integer>();
		cT.add(initHBonds.size());
		nT.add(0);
		//Loop through to determine c(t) and n(t) for all timesteps
		timestep += 1;
		while(!lastStep) {
			int numFormed = 0;
			int proximityPop = 0;
			int numBroken = 0;
			int bondedPop = 0;
			//Update atom coordinates
			lastStep = timesteps.nextTimestep(atoms, currentSystem);
			
			//Update neighbor list---
			/*atomBin = new Bin3D<Atom>(maxHBondLength, currentSystem.getXMin(),
					currentSystem.getXMax(), currentSystem.getYMin(), currentSystem.getYMax(),
					currentSystem.getZMin(), currentSystem.getZMax());
			for(Atom atom : atoms) {
				atomBin.addFrom3Space(atom, atom.getXCord(), atom.getYCord(), atom.getZCord());
			}*/
			
			//Examine status of H-Bond for initial population of H-Bonds
			for(int i = 0; i < initHBonds.size(); i++) {
				//If H-Bonded, count towards c(t) population
				if(this.isHBonded(initHBonds.get(i).getAcceptorAtom(), initHBonds.get(i).getDonatorAtom(), initHBonds.get(i).getProtonAtom())) {
					bondedPop += 1;
				} //If not H-Bonded but still in proximity, count towards n(t) population
				else if(Geometry.PeriodicDistance(initHBonds.get(i).getAcceptorAtom(), initHBonds.get(i).getDonatorAtom(), currentSystem.getXBox(),
						currentSystem.getYBox(), currentSystem.getZBox()) < (maxHBondLength)) {
					//Check to see if either atoms in the pair are in a different H-Bond
					//Boolean hBondFlag = false;
					//Do not count bonds where roles reversed
					/*if(this.isHBonded(initHBonds.get(i).getDonatorAtom(), initHBonds.get(i).getAcceptorAtom())) {
						hBondFlag = true;
					}*/
					//Do not count pair if the acceptor atom is accepting from other pairs
					/*ArrayList<Atom> atomBinAcc = atomBin.getNeighbors(initHBonds.get(i).getAcceptorAtom().getXCord(), initHBonds.get(i).getAcceptorAtom().getYCord(),
							initHBonds.get(i).getAcceptorAtom().getZCord());
					for(int j = 0; j < atomBinAcc.size(); j++) {
						if(atomBinAcc.get(j).getType().equals("1") || atomBinAcc.get(j).getType().equals("2")) {
							//if(this.isHBonded(initHBonds.get(i).getAcceptorAtom(), atomBinAcc.get(j)) || this.isHBonded(atomBinAcc.get(j), initHBonds.get(i).getAcceptorAtom())) {
							if(this.isHBonded(initHBonds.get(i).getAcceptorAtom(), atomBinAcc.get(j))) {	
								hBondFlag = true;
							}
						}
					}
					//Do not count pair if the donator atom is donating to other pairs
					ArrayList<Atom> atomBinDon = atomBin.getNeighbors(initHBonds.get(i).getDonatorAtom().getXCord(), initHBonds.get(i).getDonatorAtom().getYCord(),
							initHBonds.get(i).getDonatorAtom().getZCord());
					for(int j = 0; j < atomBinDon.size(); j++) {
						if(atomBinDon.get(j).getType().equals("1") || atomBinDon.get(j).getType().equals("2")) {
							//if(this.isHBonded(initHBonds.get(i).getDonatorAtom(), atomBinDon.get(j)) || this.isHBonded(atomBinDon.get(j), initHBonds.get(i).getDonatorAtom())) {
							if(this.isHBonded(atomBinDon.get(j), initHBonds.get(i).getDonatorAtom())) {	
								hBondFlag = true;
							}
						}
					}*/
					//If hBondFlag is still false, then neither atoms in the pair are in a different H-Bond
					/*if(!hBondFlag) {
						proximityPop += 1;
					}*/
					proximityPop += 1;
				}
			}
			
			cT.add(bondedPop);
			nT.add(proximityPop);
			
			timestep += 1;
		} //END While loop
		
		//Print out c(t) and n(t)
		System.out.println("Timestep(ps) c(t) n(t)");
		for(int i = 0; i < timestep; i++) {
			System.out.println(new DecimalFormat("#.##").format(i*timeInterval) + " " + cT.get(i) + " " + nT.get(i));
		}
		double endCT = (double)cT.get(timestep-1)/cT.get(0);
		System.out.println("End value for c(t) (including transient effects): " + endCT);
	}
	
	/**
	 * 
	 * @param acceptor The atom accepting a proton
	 * @param donor The atom donating a proton
	 * @param proton A proton covalently bonded to the donor atom
	 * @return True if the donor atom is donating the specified proton to the acceptor atom, false otherwise
	 */
	private boolean isHBonded(Atom acceptor, Atom donor, Atom proton) {
		//If the bonding specifications are met
		//First check distance criteria
		if(Geometry.PeriodicDistance(acceptor, donor, currentSystem.getXBox(), currentSystem.getYBox(),
				currentSystem.getZBox()) < maxHBondLength) {
			//Then check angle criteria
			double angle;
			
			angle = Geometry.PeriodicAngle(donor, proton, acceptor, currentSystem.getXBox(),
						currentSystem.getYBox(), currentSystem.getZBox());
			if(Math.abs(angle) < angleCutoff) {
					return true;
			}
		}
		//Otherwise, they are not hydrogen bonded
		return false;
	}
}
