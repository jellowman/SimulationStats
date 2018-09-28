package calculations;

import java.util.ArrayList;
import java.util.Scanner;

import core.Atom;
import core.BulkSystem;
import core.Molecule;
import lammps.OutParser;

public class AtomDensity {
	
	/**Holds info of the atoms throughout all timesteps (assumes no new atoms added)*/
	private ArrayList<Atom> atoms;
	
	/**Holds molecules with atoms associated with the above atoms list*/
	private ArrayList<Molecule> molecules;
	
	/**Holds info about the simulation box dimensions*/
	private BulkSystem currentSystem;
	
	/**Modifier to determine number of bins. Bin size will be approximately this size.*/
	private static final double BIN_SIZE = 2.0;
	
	public AtomDensity() {
		atoms = new ArrayList<Atom>();
		molecules = new ArrayList<Molecule>();
		currentSystem = new BulkSystem();
	}
	
	public void runAnalysis() {
		Scanner sc = new Scanner(System.in);
		System.out.println("Provide lammps output file for molecule density analysis");
		String fileName = sc.nextLine();
		runAnalysis(fileName);
	}
	
	public void runAnalysis(String filename) {
		//Create parser to read in new atomic coordinates for each timestep
		OutParser timesteps = new OutParser(filename);
				
		//Create lists of atoms and molecule groups
		//ArrayList<Molecule> molecules = new ArrayList<Molecule>();
		timesteps.buildSystem(atoms, molecules);
		
		//Create separate list for just Carbon atoms
		//NOT USED
		ArrayList<Atom> carbon = new ArrayList<Atom>();
		for(int i = 0; i < atoms.size(); i++) {
			if(atoms.get(i).getType().equals("2") || atoms.get(i).getType().equals("2")) {
				carbon.add(atoms.get(i));
			}
		}
		
		int numSteps = 0;
		double[] avgBins = new double[100];
		for(int i=0;i<avgBins.length; i++) {
			avgBins[i] = 0;
		}
		
		boolean lastStep = false;
		while(!lastStep) {
			//Update atom coordinates
			lastStep = timesteps.nextTimestep(atoms, currentSystem);
			
			//Make bins along z-axis to count water and glycerol
			int binsZ = (int) ((currentSystem.getZMax()-currentSystem.getZMin())/BIN_SIZE);
			double dz = (currentSystem.getZMax()-currentSystem.getZMin())/binsZ;
			int[] gBin = new int[binsZ];
			int[] wBin = new int[binsZ];
			
			for(int i = 0; i < gBin.length; i++) {
				gBin[i] = 0;
				wBin[i] = 0;
			}
			
			int[][] molCount = new int[binsZ][molecules.size()+1];
			for(int i = 0; i < binsZ; i++) {
				for(int j = 0; j < molecules.size()+1; j++) {
						molCount[i][j] = 0;
				}
			}
			
			for(int i = 0; i < atoms.size(); i++) {
				double zOffset = atoms.get(i).getZCord() - currentSystem.getZMin();
				int zAdd = (int) (zOffset/dz);
				int zAddd = (zAdd+binsZ)%binsZ;
				//if(atoms.get(i).getMolID() < 146) { //If it's in glycerol
				if(atoms.get(i).getType().equals("1")) { //Treat outer oxygens as COM
					//if(molCount[zAddd][atoms.get(i).getMolID()] == 0) { //If molecule hasn't been counted yet
					if(true) {
						gBin[zAddd] += 1;
						molCount[zAddd][atoms.get(i).getMolID()] = 1;
					}
				//} else {
				} else if(atoms.get(i).getType().equals("2")) {
					if(molCount[zAddd][atoms.get(i).getMolID()] == 0) {
						wBin[zAddd] += 1;
						molCount[zAddd][atoms.get(i).getMolID()] = 1;
					}
				}
			}
			
			for(int i = 0; i < gBin.length; i++) {
				double gWeight;
				if(gBin[i]+wBin[i] != 0) {
					gWeight = (double)(gBin[i])/(gBin[i]+wBin[i]);
				} else {
					gWeight = 0;
				}
				System.out.println(i + " - " + gWeight);
				avgBins[i] += gWeight;
			}
			numSteps += 1;
		}
		
		for(int i = 0; i < avgBins.length; i++) {
			avgBins[i] = avgBins[i] / numSteps;
			System.out.println(i + " - " + avgBins[i]);
		}
	}
	
	
}
