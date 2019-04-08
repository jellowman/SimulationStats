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
public class HBondAngleDist extends RadialDF2
{	
	//2 is water, 1 is glycerol
	
	public HBondAngleDist(double cutoff, String centerAtom, String otherAtom) {
		super(cutoff, centerAtom, otherAtom);
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
		
		int bins = 180/5;
		double dz = 5;
		Integer[] angleHist = new Integer[bins];
		
		for(int i=0; i < angleHist.length; i++) {
			angleHist[i] = 0;
		}
		
		
		boolean lastStep = false;
		int timestep = 0;
		//Determine H-Bond population for timestep and add to total if a new one forms.
		//Determine proximity population by pairs that are close but are not bonded, nor the opposing pair
		
		
		while(!lastStep) {
			lastStep = timesteps.nextTimestep(atoms, currentSystem);
			
			//Update neighbor list---
			Bin3D<Atom> atomBin = new Bin3D<Atom>(cutoff, currentSystem.getXMin(),
					currentSystem.getXMax(), currentSystem.getYMin(), currentSystem.getYMax(),
					currentSystem.getZMin(), currentSystem.getZMax());
			
			for (Atom atom : centerAtoms) {
				atomBin.addFrom3Space(atom, atom.getXCord(), atom.getYCord(), atom.getZCord());
			}
			
			if (!centerAtom.equals(otherAtom)) {
				for (Atom atom : otherAtoms) {
					atomBin.addFrom3Space(atom, atom.getXCord(), atom.getYCord(), atom.getZCord());
				}
			}
			
			for(int i = 0; i < centerAtoms.size(); i++) {
				//Get neighbors to centerAtom
				ArrayList<Atom> neigh = atomBin.getNeighbors(centerAtoms.get(i).getXCord(),
							centerAtoms.get(i).getYCord(), centerAtoms.get(i).getZCord());
				for(int j = 0; j < neigh.size(); j++) {
					if(centerAtoms.get(i).getMolID() != neigh.get(j).getMolID()) {
						//Will visit both pair types if both acceptor and donor are the same atom type
						if(neigh.get(j).getType().equals(otherAtom)) {
							if(Geometry.PeriodicDistance(centerAtoms.get(i), neigh.get(j), currentSystem.getXBox(),
									currentSystem.getYBox(), currentSystem.getZBox()) < (cutoff)) {
								//Put in bin based on angle
								//If water molecule is donor, get angle for both H-atoms
								double angle;
								if(centerAtoms.get(i).getType().equals("2")) {
									//if(Geometry.PeriodicDistance(atoms.get(neigh.get(j).getID()+1), oxygen.get(i), currentSystem.getXBox(),
									//currentSystem.getYBox(), currentSystem.getZBox()) < Geometry.PeriodicDistance(neigh.get(j), oxygen.get(i), currentSystem.getXBox(),
									//currentSystem.getYBox(), currentSystem.getZBox())) {
										angle = Geometry.PeriodicAngle(centerAtoms.get(i), atoms.get(centerAtoms.get(i).getID()+1), neigh.get(j), currentSystem.getXBox(),
												currentSystem.getYBox(), currentSystem.getZBox());
										angleHist[(int)(angle/dz)] += 1;
									//}
									
									//if(Geometry.PeriodicDistance(atoms.get(neigh.get(j).getID()+2), oxygen.get(i), currentSystem.getXBox(),
									//		currentSystem.getYBox(), currentSystem.getZBox()) < Geometry.PeriodicDistance(neigh.get(j), oxygen.get(i), currentSystem.getXBox(),
									//		currentSystem.getYBox(), currentSystem.getZBox())) {
										angle = Geometry.PeriodicAngle(centerAtoms.get(i), atoms.get(centerAtoms.get(i).getID()+2), neigh.get(j), currentSystem.getXBox(),
												currentSystem.getYBox(), currentSystem.getZBox());
										angleHist[(int)(angle/dz)] += 1;
									//}
								} else {
									//if(Geometry.PeriodicDistance(atoms.get(neigh.get(j).getID()+1), oxygen.get(i), currentSystem.getXBox(),
									//		currentSystem.getYBox(), currentSystem.getZBox()) < Geometry.PeriodicDistance(neigh.get(j), oxygen.get(i), currentSystem.getXBox(),
									//				currentSystem.getYBox(), currentSystem.getZBox())){
										angle = Geometry.PeriodicAngle(centerAtoms.get(i), atoms.get(centerAtoms.get(i).getID()+1), neigh.get(j), currentSystem.getXBox(),
											currentSystem.getYBox(), currentSystem.getZBox());
										angleHist[(int)(angle/dz)] += 1;
									//}
								}
							}
						}
					}
				}
			}
			
			//System.out.println(timestep);
			timestep += 1;
		}
		
		//Print out angle hist
		//double[] angleNorm = new double[angleHist.length];
		int numDonors = centerAtoms.size();
		
		double integration = 0;
		int degeneracy = 1;
		if(centerAtom == "2") {
			degeneracy = 2;
		}
		for(int i=0; i < angleHist.length; i++) {
			//normalize based on volume
			/*double sin1 = Math.sin((i+1)*dz);
			double sin2 = Math.sin(i*dz);
			double angle12 = Math.asin(sin1/MAX_H_BOND_LENGTH*0.98);
			double angle22 = Math.asin(sin2/MAX_H_BOND_LENGTH*0.98);
			double angle13 = 180 - (i+1)*dz - angle12;
			double angle23 = 180 - i*dz - angle22; 
			
			Integration result * 2pi / 2
			double area = Math.PI * Math.abs(0.98 * MAX_H_BOND_LENGTH * (-Math.cos(Math.toRadians(angle13))+Math.cos(Math.toRadians(angle23))));
			double area = Math.abs((Math.sin(Math.toRadians((i+1)*dz))-Math.sin(Math.toRadians(i*dz)))*0.5*0.98*b)*2*Math.PI;
			angleNorm[i] = angleHist[i] / area;
			System.out.println((i*dz+2.5) + " " + angleNorm[i]);
			*/
			//Normalize to avg. num neighbors
			double part = (double)(angleHist[i])/numDonors/timestep;
			integration += part;
			System.out.println((i*dz+2.5) + " " + part + " " + integration/degeneracy);
		}
		
		System.out.println("There are " + (integration/degeneracy) + " average acceptor neighbors per donor");
		
	}
}
