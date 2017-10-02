package lammps;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Scanner;

import core.Atom;
//import core.BulkSystem;
import core.Molecule;
import core.Molecule.Angle;
import core.Molecule.Bond;
import core.Molecule.Dihedral;

public class DatWriter 
{
	/**An offset to ensure that molecules do not intersect with each other over periodic boundaries*/
	private static final int SIMULATION_MARGIN = 2;
	
	/*
	/**
	 * <code>writeFile</code> called on a <code>BulkSystem</code> will generate just atomic coordinates
	 * for a .dat file for lammps input.
	 * @param bs	The system containing all the atoms in a system
	 *
	public writeFile(BulkSystem bs)
	{
		System.out.println("Please provide a name for the Lammps input file.");
	}*/
	
	/**
	 * <code>writeFile</code> called on an <code>ArrayList</code> of type <code>Molecule</code> will generate atomic
	 * coordinates and any additional properties for a lammps .dat input file
	 * such as bonds and angles between atoms.
	 * @param mols	The array of molecules in the system
	 */
	public static void writeFile(ArrayList<Molecule> mols)
	{
		//Initialize the file.
		@SuppressWarnings("resource")
		Scanner sc = new Scanner(System.in);
		System.out.println("Please provide a name for the Lammps input file.");
		String name = sc.nextLine();
		try{
		File file = new File(name);
		FileWriter fr = new FileWriter(file);
		BufferedWriter br = new BufferedWriter(fr);
		
		//Read through the molecules to see how many atoms, bonds, and angles there are
		HashMap<String, Integer> atomTypeToNumType = new HashMap<String, Integer>();
		HashMap<String, Double> atomTypeToMass = new HashMap<String, Double>();
		HashMap<String, Double> atomTypeToCharge = new HashMap<String, Double>();
		writeCBA(br, mols, atomTypeToNumType, atomTypeToMass, atomTypeToCharge);
		
		//Write out Atoms
		writeAtoms(br, mols, atomTypeToNumType, atomTypeToCharge);
		//Write out Bonds
		writeBonds(br, mols);
		//Write out Angles
		writeAngles(br, mols);
		//Write out Dihedrals
		writeDihedrals(br, mols);
		
		System.out.println("Finished!");
		
		br.close();
		fr.close();
		} catch(IOException io)
		{
			
		}
	}
	
	/**
	 * Writes the coordinates, bonds, and angles (CBA) to the file linked to the <code>BufferedWriter</code>.
	 */
	private static void writeCBA(BufferedWriter br, ArrayList<Molecule> mols, HashMap<String, Integer> atomTypeToNumType,
			HashMap<String, Double> atomTypeToMass, HashMap<String, Double> atomTypeToCharge)
	{
		int numAtoms = 0, numBonds = 0, numAngles = 0, numDihedrals = 0;
		double maxX = 0, maxY = 0, maxZ = 0;
		HashSet<String> atomTypes = new HashSet<String>();
		HashSet<Integer> bondTypes = new HashSet<Integer>();
		HashSet<Integer> angleTypes = new HashSet<Integer>();
		HashSet<Integer> dihedralTypes = new HashSet<Integer>();
		
		for(Molecule mol : mols)
		{
			for(Atom atom : mol.getAtoms())
			{
				numAtoms++;
				atomTypes.add(atom.getType());
				if(atom.getXCord() > maxX)
				{
					maxX = atom.getXCord();
				}
				if(atom.getYCord() > maxY)
				{
					maxY = atom.getYCord();
				}
				if(atom.getZCord() > maxZ)
				{
					maxZ = atom.getZCord();
				}
			}
			for(Bond bond : mol.getBonds())
			{
				numBonds++;
				bondTypes.add(bond.getID());
			}
			for(Angle angle : mol.getAngles())
			{
				numAngles++;
				angleTypes.add(angle.getID());
			}
			for(Dihedral dihedral : mol.getDihedrals())
			{
				numDihedrals++;
				dihedralTypes.add(dihedral.getID());
			}
		}
		
		//Write info to header of file
		try {
		br.newLine(); br.newLine();
		br.write("\t\t" + numAtoms + " atoms"); br.newLine();
		br.write("\t\t" + numBonds + " bonds"); br.newLine();
		br.write("\t\t" + numAngles + " angles"); br.newLine();
		br.write("\t\t" + numDihedrals + " dihedrals"); br.newLine();
		br.newLine();
		br.write("\t\t\t" + atomTypes.size() + " atom types"); br.newLine();
		br.write("\t\t\t" + bondTypes.size() + " bond types"); br.newLine();
		br.write("\t\t\t" + angleTypes.size() + " angle types"); br.newLine();
		br.write("\t\t\t" + dihedralTypes.size() + " dihedral types"); br.newLine();
		br.newLine();
		
		//Write max and min simulation coordinates
		br.write(" 0.0000000E+00   " + String.format("%-15.12f", maxX + SIMULATION_MARGIN) + "      xlo xhi"); br.newLine();
		br.write(" 0.0000000E+00   " + String.format("%-15.12f", maxY + SIMULATION_MARGIN) + "      ylo yhi"); br.newLine();
		br.write(" 0.0000000E+00   " + String.format("%-15.12f", maxZ + SIMULATION_MARGIN) + "      zlo zhi"); br.newLine(); br.newLine(); br.newLine();
		
		//Write masses and charge of each atom type
		br.write(" Masses"); br.newLine(); br.newLine();
		
		int atomNum = 1;
		@SuppressWarnings("resource")
		Scanner sc = new Scanner(System.in);
		sc.reset();
		for(String atomType : atomTypes)
		{
			//Link Atom type to a number
			atomTypeToNumType.put(atomType, atomNum);
			atomNum++;
			
			//Link Atom type to a mass
			System.out.println("Provide the mass of the atom type " + atomType);
			Double mass = Double.valueOf(sc.nextLine());
			atomTypeToMass.put(atomType, mass);
			
			//Link Atom type to a charge
			System.out.println("Provide the charge of the atom type " + atomType);
			Double charge = Double.valueOf(sc.nextLine());
			atomTypeToCharge.put(atomType, charge);
			
			br.write(" " + atomTypeToNumType.get(atomType) + "   " + atomTypeToMass.get(atomType)); br.newLine();
		}
		
		br.newLine(); br.newLine();
		
		} catch (IOException e) {
			System.err.println("Error in printing .dat header.");
		}
		
	}
	
	private static void writeAtoms(BufferedWriter br, ArrayList<Molecule> mols, HashMap<String, Integer> atomTypeToNumType,
			HashMap<String, Double> atomTypeToCharge)
	{
		try {
			br.write(" Atoms"); br.newLine(); br.newLine();
			
			int molNum = 1;
			int atomNum = 1;
			for(Molecule mol : mols)
			{
				for(Atom atom : mol.getAtoms())
				{
					String formatted = String.format("%10d %6d %6d %12f %12f %12f %12f", atomNum, molNum, atomTypeToNumType.get(atom.getType()),
							atomTypeToCharge.get(atom.getType()), atom.getXCord(), atom.getYCord(), atom.getZCord());
					br.write(formatted); br.newLine();
					atomNum++;
				}
			molNum++;
			}
		} catch (IOException io) {
			System.err.println("Error writing atoms.");
			io.printStackTrace();
		}
	}
	
	private static void writeBonds(BufferedWriter br, ArrayList<Molecule> mols)
	{
		try{
			br.newLine(); br.write(" Bonds"); br.newLine(); br.newLine();
			
			int bondNum = 1;
			for(Molecule mol : mols)
			{
				for(Bond bond : mol.getBonds())
				{
					String formatted = String.format("%12d %s", bondNum, bond.toString());
					br.write(formatted); br.newLine();
					bondNum++;
				}
			}
		}
		catch (IOException io)
		{
			System.err.println("Error writing bonds.");
			io.printStackTrace();
		}
	}
	
	private static void writeAngles(BufferedWriter br, ArrayList<Molecule> mols)
	{
		try{
			br.newLine(); br.write(" Angles"); br.newLine(); br.newLine();
			
			int angleNum = 1;
			for(Molecule mol : mols)
			{
				for(Angle angle : mol.getAngles())
				{
					String formatted = String.format("%12d %s", angleNum, angle.toString());
					br.write(formatted); br.newLine();
					angleNum++;
				}
			}
		}
		catch (IOException io)
		{
			System.err.println("Error writing angles.");
			io.printStackTrace();
		}
	}
	
	private static void writeDihedrals(BufferedWriter br, ArrayList<Molecule> mols)
	{
		try{
			br.newLine(); br.write(" Dihedrals"); br.newLine(); br.newLine();
			
			int dihedralNum = 1;
			for(Molecule mol : mols)
			{
				for(Dihedral dihedral : mol.getDihedrals())
				{
					String formatted = String.format("%12d %s", dihedralNum, dihedral.toString());
					br.write(formatted); br.newLine();
					dihedralNum++;
				}
			}
		}
		catch (IOException io)
		{
			System.err.println("Error writing dihedrals.");
			io.printStackTrace();
		}
	}
}
