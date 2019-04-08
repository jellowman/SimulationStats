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

/**Calculates an H-Bond autocorrelation function c(t) while estimating the
 * initial proximity population n(t)
 */
public class HBondLifeACT extends HBondLifeTC
{	
	//2 is water, 1 is glycerol
	
	public HBondLifeACT() {
		super();
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
		Bin3D<Atom> atomBin = new Bin3D<Atom>(MAX_H_BOND_LENGTH+EXTRA_DISTANCE, currentSystem.getXMin(),
				currentSystem.getXMax(), currentSystem.getYMin(), currentSystem.getYMax(),
				currentSystem.getZMin(), currentSystem.getZMax());
		
		for (Atom atom : atoms) {
			atomBin.addFrom3Space(atom, atom.getXCord(), atom.getYCord(), atom.getZCord());
		}
		
		ArrayList<HBond> initHBonds = new ArrayList<HBond>();
		LinkedHashMap<String, HBond> hBondMap = new LinkedHashMap<String, HBond>();
		//Determine initial population of H-Bonds
		
		int bondedPop = 0;
		int proximityPop = 0;
		int bothBonded = 0;
		for(int i = 0; i < oxygen.size(); i++) {
			//Get neighbors to oxygen
			ArrayList<Atom> neigh = atomBin.getNeighbors(oxygen.get(i).getXCord(),
						oxygen.get(i).getYCord(), oxygen.get(i).getZCord());
			for(int j = 0; j < neigh.size(); j++) {
				if(oxygen.get(i).getMolID() != neigh.get(j).getMolID() && 
						(neigh.get(j).getType().equals("1") || neigh.get(j).getType().equals("2"))) {
					//Will visit both pair types if both acceptor and donor are the same atom type
					if(oxygen.get(i).getType().equals(ACCEPTOR) && neigh.get(j).getType().equals(DONOR)) {
						if(this.isHBonded(oxygen.get(i), neigh.get(j))) {
							//TODO Add conditional statement to add only 1 combination of pairs for same-type bonds
							//This will allow for dynamics without allowing donor-acceptor role to switch
							HBond newBond = new HBond(oxygen.get(i), neigh.get(j), 0);
							initHBonds.add(newBond);
							String key = String.valueOf(oxygen.get(i).getID()) + "-" + String.valueOf(neigh.get(j).getID());
							hBondMap.put(key, newBond);
							bondedPop += 1;
							if(this.isHBonded(neigh.get(j), oxygen.get(i))) {
								bothBonded += 1;
							}
						} else if(Geometry.PeriodicDistance(oxygen.get(i), neigh.get(j), currentSystem.getXBox(),
								currentSystem.getYBox(), currentSystem.getZBox()) < (MAX_H_BOND_LENGTH + EXTRA_DISTANCE)) {
							//If not bonded, close so add to prox. population
							//Note how if both are same atom types, proximity population double counts pairs.
							HBond newBond = new HBond(oxygen.get(i), neigh.get(j), 0);
							initHBonds.add(newBond);
							String key = String.valueOf(oxygen.get(i).getID()) + "-" + String.valueOf(neigh.get(j).getID());
							hBondMap.put(key, newBond);
							//Check to see if the reverse pair is H-Bonded
							//Add to proximity population if it is not present
							if(!this.isHBonded(neigh.get(j), oxygen.get(i))) {
								proximityPop += 1;
							}
						}
					}
				}
			}
		}
		System.out.println("There are " + bothBonded/2.0 + " instances of atoms donating and accepting simultaneously.");
		
		ArrayList<Integer> cT = new ArrayList<Integer>();
		ArrayList<Integer> nT = new ArrayList<Integer>();
		cT.add(bondedPop);
		nT.add(proximityPop);
		//Loop through to determine c(t) and n(t) for all timesteps
		timestep += 1;
		while(!lastStep) {
			//int numFormed = 0;
			proximityPop = 0;
			//int numBroken = 0;
			bondedPop = 0;
			//Update atom coordinates
			lastStep = timesteps.nextTimestep(atoms, currentSystem);
			
			//Update neighbor list---
			atomBin = new Bin3D<Atom>(MAX_H_BOND_LENGTH+EXTRA_DISTANCE, currentSystem.getXMin(),
					currentSystem.getXMax(), currentSystem.getYMin(), currentSystem.getYMax(),
					currentSystem.getZMin(), currentSystem.getZMax());
			for(Atom atom : atoms) {
				atomBin.addFrom3Space(atom, atom.getXCord(), atom.getYCord(), atom.getZCord());
			}
			
			//Examine status of H-Bond for initial population of H-Bonds
			for(int i = 0; i < initHBonds.size(); i++) {
				//If H-Bonded, count towards c(t) population
				if(this.isHBonded(initHBonds.get(i).getAcceptorAtom(), initHBonds.get(i).getDonatorAtom())) {
					bondedPop += 1;
				} //If not H-Bonded but still in proximity, count towards n(t) population
				else if(Geometry.PeriodicDistance(initHBonds.get(i).getAcceptorAtom(), initHBonds.get(i).getDonatorAtom(), currentSystem.getXBox(),
						currentSystem.getYBox(), currentSystem.getZBox()) < (MAX_H_BOND_LENGTH + EXTRA_DISTANCE)) {
					/*//Check to see if either atoms in the pair are in a different H-Bond
					Boolean hBondFlag = false;
					
					ArrayList<Atom> atomBinAcc = atomBin.getNeighbors(initHBonds.get(i).getAcceptorAtom().getXCord(), initHBonds.get(i).getAcceptorAtom().getYCord(),
							initHBonds.get(i).getAcceptorAtom().getZCord());
					for(int j = 0; j < atomBinAcc.size(); j++) {
						if(atomBinAcc.get(j).getType().equals("1") || atomBinAcc.get(j).getType().equals("2")) {
							if(this.isHBonded(initHBonds.get(i).getAcceptorAtom(), atomBinAcc.get(j)) || this.isHBonded(atomBinAcc.get(j), initHBonds.get(i).getAcceptorAtom())) {
								hBondFlag = true;
							}
						}
					}
					ArrayList<Atom> atomBinDon = atomBin.getNeighbors(initHBonds.get(i).getDonatorAtom().getXCord(), initHBonds.get(i).getDonatorAtom().getYCord(),
							initHBonds.get(i).getDonatorAtom().getZCord());
					for(int j = 0; j < atomBinDon.size(); j++) {
						if(atomBinDon.get(j).getType().equals("1") || atomBinDon.get(j).getType().equals("2")) {
							if(this.isHBonded(initHBonds.get(i).getDonatorAtom(), atomBinDon.get(j)) || this.isHBonded(atomBinDon.get(j), initHBonds.get(i).getDonatorAtom())) {
								hBondFlag = true;
							}
						}
					}
					//If hBondFlag is still false, then neither atoms in the pair are in a different H-Bond
					if(!hBondFlag) {
						proximityPop += 1;
					}*/
					
					//Check to see if opposite pair is H-Bonded
					//Add to proximity population neither pair combination is H-Bonded
					if(!this.isHBonded(initHBonds.get(i).getDonatorAtom(), initHBonds.get(i).getAcceptorAtom())) {
						proximityPop += 1;
					}
				}
			}
			
			cT.add(bondedPop);
			nT.add(proximityPop);
			
			timestep += 1;
		} //END While loop
		
		//Print out c(t) and n(t)
		System.out.println("Timestep(ps) c(t) n(t)");
		for(int i = 0; i < timestep; i++) {
			System.out.println(new DecimalFormat("#.#").format(i*TIMESTEP) + " " + cT.get(i) + " " + nT.get(i));
		}
	}
}
