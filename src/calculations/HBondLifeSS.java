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

import calculations.HBondLifeTC.HBond;
import core.Atom;
import core.Bin3D;
import core.BulkSystem;
import core.Molecule;
import lammps.OutParser;

/**Calculates the fractions of H-Bonds, at a given time, that are within the criteria
 * and violating the criteria. Can be used for initial conditions for c(t) relaxations.
 */
public class HBondLifeSS extends HBondLifeTC
{	
	//2 is water, 1 is glycerol
	
	public HBondLifeSS() {
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
		
		ArrayList<HBond> totalHBonds = new ArrayList<HBond>();
		LinkedHashMap<String, HBond> hBondMap = new LinkedHashMap<String, HBond>();
		//Determine initial population of H-Bonds
		
		ArrayList<Integer> cT = new ArrayList<Integer>();
		ArrayList<Integer> nT = new ArrayList<Integer>();
		ArrayList<Integer> unique = new ArrayList<Integer>();
		
		boolean lastStep = false;
		int timestep = 0;
		//Determine H-Bond population for timestep and add to total if a new one forms.
		//Determine proximity population by pairs that are close but are not bonded, nor the opposing pair
		
		//Loop through to determine c(t) and n(t) for all timesteps
		while(!lastStep) {
			
			int bondedPop = 0;
			int proximityPop = 0;
			int bothBonded = 0;
			
			lastStep = timesteps.nextTimestep(atoms, currentSystem);
			
			//Update neighbor list---
			Bin3D<Atom> atomBin = new Bin3D<Atom>(MAX_H_BOND_LENGTH+EXTRA_DISTANCE, currentSystem.getXMin(),
					currentSystem.getXMax(), currentSystem.getYMin(), currentSystem.getYMax(),
					currentSystem.getZMin(), currentSystem.getZMax());
			
			for (Atom atom : atoms) {
				atomBin.addFrom3Space(atom, atom.getXCord(), atom.getYCord(), atom.getZCord());
			}
			
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
								bondedPop += 1;
								HBond newBond = new HBond(oxygen.get(i), neigh.get(j), timestep);
								String key = String.valueOf(oxygen.get(i).getID()) + "-" + String.valueOf(neigh.get(j).getID());
								HBond check = hBondMap.putIfAbsent(key, newBond);
								//if check is null, then it's a new H-Bond pair not seen before.
								if(check == null) {
									totalHBonds.add(newBond);
								}
								/* Check for instances where both atoms in pair are accepting and donating to each other
								if(this.isHBonded(neigh.get(j), oxygen.get(i))) {
									bothBonded += 1;
								}*/
							} else if(Geometry.PeriodicDistance(oxygen.get(i), neigh.get(j), currentSystem.getXBox(),
									currentSystem.getYBox(), currentSystem.getZBox()) < (MAX_H_BOND_LENGTH + EXTRA_DISTANCE)) {
								//If not bonded but close, see if the pair has been H-Bonded before
								String key = String.valueOf(oxygen.get(i).getID()) + "-" + String.valueOf(neigh.get(j).getID());
								HBond check = hBondMap.get(key);
								//If check is not null, then H-Bond has existed before and can be in prox. pop
								if(check != null) {
									//Check to see if the reverse pair is H-Bonded
									//Add to proximity population if opposite pairing is not currently H-Bonded
									if(!this.isHBonded(neigh.get(j), oxygen.get(i))) {
										proximityPop += 1;
									}
								}
							}
						}
					}
				}
			}
			//System.out.println("There are " + bothBonded/2.0 + " instances of atoms donating and accepting simultaneously.");
			
			
			cT.add(bondedPop);
			nT.add(proximityPop);
			unique.add(totalHBonds.size());
			timestep += 1;
		}
		
		//Print out c(t) and n(t)
		System.out.println("Timestep(ps) c(t) n(t) [Total unique H-Bonds seen] [Num diffused away]");
		for(int i = 0; i < timestep; i++) {
			System.out.println(new DecimalFormat("#.##").format(i*TIMESTEP) + " " + cT.get(i) + " " + nT.get(i) + " " + unique.get(i) + " " + (unique.get(i)-cT.get(i)-nT.get(i)));
		}
	}
}
