package lammps;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;

import core.Atom;
import core.Molecule;

public class OutParser 
{
	/**FileReader for LAMMPS output file*/
	
	public OutParser(String fileName) {
		try {
		File file = new File(fileName);
		FileReader fr = new FileReader(file);
		BufferedReader br = new BufferedReader(fr);
		} catch (IOException ex) {
			
		}
	}
	
	public void buildSystem(ArrayList<Atom> atoms, ArrayList<Molecule> molecules) {
		
	}
	public boolean nextTimestep(ArrayList<Atom> atoms) {
		return true;
	}
}
