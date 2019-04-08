package calculations;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
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

/**Takes a simulation trajectory file and appends a column of data indicating which atoms
 * are in the specified H-Bond
 */
public class HBondLifeTrajWriter extends HBondLifeTC
{	
	//2 is water, 1 is glycerol
	
	public HBondLifeTrajWriter() {
		super();
	}
	
	
	public void runAnalysis(String filename) {
		//Create parser to read in new atomic coordinates for each timestep
		OutParser timesteps = new OutParser(filename);
		
		//Create file writer for modified trajectory file
		try {
		File hBondFile = new File(filename.substring(0, filename.lastIndexOf('.', filename.length()))+"HB"+".lammpstrj");
		FileWriter hBondFR;
		hBondFR = new FileWriter(hBondFile);
		BufferedWriter hBondBR = new BufferedWriter(hBondFR);
		
		//Write first line in trajectory file for loop offset
		hBondBR.write("ITEM: TIMESTEP"); hBondBR.newLine();
		
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
		
		boolean lastStep = false;
		int timestep = 0;
		while(!lastStep) {
			//Hold onto raw line data from Trajectory File
			ArrayList<String> lineList = new ArrayList<String>();
			//Update atom coordinates
			lastStep = timesteps.nextTimestep(atoms, currentSystem, lineList);
			
			//Set all atoms to unbonded
			for(Atom atom : oxygen) {
				atom.setisHBonded(false);
			}
			
			//Update neighbor list---
			Bin3D<Atom> atomBin = new Bin3D<Atom>(MAX_H_BOND_LENGTH+EXTRA_DISTANCE, currentSystem.getXMin(),
					currentSystem.getXMax(), currentSystem.getYMin(), currentSystem.getYMax(),
					currentSystem.getZMin(), currentSystem.getZMax());
			for(Atom atom : atoms) {
				atomBin.addFrom3Space(atom, atom.getXCord(), atom.getYCord(), atom.getZCord());
			}
			
			//Look over whole system for H-Bonds
			for(int i = 0; i < oxygen.size(); i++) {
				//Get neighbors to oxygen
				ArrayList<Atom> neigh = atomBin.getNeighbors(oxygen.get(i).getXCord(),
						oxygen.get(i).getYCord(), oxygen.get(i).getZCord());
				for(int j = 0; j < neigh.size(); j++) {
					if(oxygen.get(i).getMolID() != neigh.get(j).getMolID() && 
							(neigh.get(j).getType().equals("1") || neigh.get(j).getType().equals("2"))) {
						if(this.isHBonded(oxygen.get(i), neigh.get(j))) {
							if(oxygen.get(i).getType().equals(ACCEPTOR) && neigh.get(j).getType().equals(DONOR)) {
								//Indicate that atoms are H-Bonded
								//Only indicate the donating atom
								//atoms.get(oxygen.get(i).getID()).setisHBonded(true);
								atoms.get(neigh.get(j).getID()).setisHBonded(true);
							}
						}
					}
				}
			}
			
			//Write out timestep information to new file
			
			//Write out header information for new timestep
			for(int i = 0; i < 7; i++) {
				hBondBR.write(lineList.get(i)); hBondBR.newLine();
			}
			//Modify column info line to indicate new data addition
			hBondBR.write(lineList.get(7)+"hbond "); hBondBR.newLine();
			//Rest of the lines are atom data except last line unless it is the last timestep
			for(int i = 8; i < lineList.size()-1; i++) {
				//Note index 8 is start of atoms, but first atom should be index 0
				int hBondStatus = 0;
				if(atoms.get(i-8).isHBonded()) {
					hBondStatus = 1;
				}
				hBondBR.write(lineList.get(i)+ " " + hBondStatus + " "); hBondBR.newLine();
			}
			//Check if last line is timestep line or last atom at end of file
			if(lineList.get(lineList.size()-1).charAt(0) == ' ') {
				int hBondStatus = 0;
				if(atoms.get(lineList.size()-9).isHBonded()) {
					hBondStatus = 1;
				}
				//If first character is blank space, then it's an atom
				hBondBR.write(lineList.get(lineList.size()-1) + " " + hBondStatus + " "); hBondBR.newLine();
			} else { //Otherwise, write the timestep line
				hBondBR.write(lineList.get(lineList.size()-1)); hBondBR.newLine();
			}
			
			timestep += 1;
			System.out.println("TIMESTEP: " + timestep);
		//END While Loop
		}
		System.out.println("Done writing H-Bond info to " + hBondFile.getName());
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
}
