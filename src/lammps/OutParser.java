package lammps;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Scanner;

import core.Atom;
import core.BulkSystem;
import core.Molecule;

public class OutParser 
{
	private FileReader fr;
	private BufferedReader br;
	
	/**FileReader for LAMMPS output file*/
	
	public OutParser(String fileName) {
		try {
		File file = new File(fileName);
		fr = new FileReader(file);
		br = new BufferedReader(fr);
		//Read first line in to have same offset
		br.readLine();
		} catch (IOException ex) {
			
		}
	}
	
	//--Move to DatReader---
	//Reads atoms and molecules from .dat file. Used b/c trajectory file
	//Stores atom and molecule data in passed in Atom and Molecule arrays
	//Does not output bond information
	public void buildSystem(ArrayList<Atom> atoms, ArrayList<Molecule> molecules) {
		Scanner sc = new Scanner(System.in);
		System.out.println("Enter lammps input .dat filename");
		String datName = sc.nextLine();
		File datFile = new File(datName);
		try{
			FileReader datFr = new FileReader(datFile);
			BufferedReader datBr = new BufferedReader(datFr);
			
			//Skip beginning of .dat file to atoms section
			String nextLine = datBr.readLine();
			while(!nextLine.contains("Atoms")) {
				System.out.println(nextLine);
				nextLine = datBr.readLine();
			}
			datBr.readLine();
			
			nextLine = datBr.readLine();
			//System.out.println(nextLine);
			while(!nextLine.isEmpty()) {
				//System.out.println(nextLine);
				String parts[] = nextLine.split("\\s+");
				//Make it 0-indexed
				int id = Integer.valueOf(parts[1]) - 1;
				//double x = Double.valueOf(parts[4]);
				//double y = Double.valueOf(parts[5]);
				//double z = Double.valueOf(parts[6]);
				Atom newAtom = new Atom(id, parts[3], 0.0, 0.0, 0.0);
				newAtom.setMolID(Integer.valueOf(parts[2]));
				atoms.add(newAtom);
				
				//Add new molecule to molecule Array if does not exist
				int molID = Integer.valueOf(parts[2]);
				if(molID == molecules.size()+1) {
					molecules.add(new Molecule());
				} else if(molID > molecules.size()+1) {
					System.err.println("Skipped molID");
					System.exit(1);
				}
				
				molecules.get(molecules.size()-1).addAtom(newAtom);
				nextLine = datBr.readLine();
			}
		} catch(IOException io) {
			System.err.println("IO error in output parser");
		}
		
	}
	
	/**
	 * 
	 * @param atoms
	 * @return True if the end of file is reached, false otherwise
	 */
	public boolean nextTimestep(ArrayList<Atom> atoms, BulkSystem system) {
		try {
			//Read header info
			String nextLine;
			for(int i = 0; i < 4; i++) {
				nextLine = br.readLine();
			}
			
			//Update box dimensions
			nextLine = br.readLine();
			String dims[] = nextLine.split("\\s+");
			//System.out.println(dims[1]);
			system.setXBox(Double.valueOf(dims[1]) - Double.valueOf(dims[0]));
			system.setXMax(Double.valueOf(dims[1]));
			system.setXMin(Double.valueOf(dims[0]));
			nextLine = br.readLine();
			dims = nextLine.split("\\s+");
			system.setYBox(Double.valueOf(dims[1]) - Double.valueOf(dims[0]));
			system.setYMax(Double.valueOf(dims[1]));
			system.setYMin(Double.valueOf(dims[0]));
			nextLine = br.readLine();
			dims = nextLine.split("\\s+");
			system.setZBox(Double.valueOf(dims[1]) - Double.valueOf(dims[0]));
			system.setZMax(Double.valueOf(dims[1]));
			system.setZMin(Double.valueOf(dims[0]));
			
			br.readLine();
			//Update atom coordinates
			nextLine = br.readLine();
			int atomID = -1;
			//System.out.println(nextLine.charAt(0));
			while(nextLine.charAt(0) == ' ') {
				String parts[] = nextLine.split("\\s+");
				//System.out.println(parts.length);
				//System.out.println(nextLine);
				atomID++;
				
				atoms.get(atomID).setXCord(Double.valueOf(parts[2]));
				atoms.get(atomID).setYCord(Double.valueOf(parts[3]));
				atoms.get(atomID).setZCord(Double.valueOf(parts[4]));
						
				nextLine = br.readLine();
				//Check for end of file
				if(nextLine == null) {
					return true;
				}
			}
			return false;
		} catch (IOException io) {
			
		}
		
		return false;
	}
}
