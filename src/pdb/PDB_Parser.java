package pdb;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.LinkedHashMap;

import core.Atom;
import core.MolBlueprint;
import core.Molecule;

public class PDB_Parser 
{
	/*protected FileReader fr;
	protected BufferedReader br;*/
	
	protected LinkedHashMap<String, MolBlueprint> blueprints;
	protected ArrayList<Molecule> molecules;
	
	public PDB_Parser()
	{
		blueprints = new LinkedHashMap<String, MolBlueprint>();
		molecules = new ArrayList<Molecule>();
	}
	
	public void parseBlueprintFile(String filename)
	{
		File file = new File(filename);
		
		try {
			FileReader fr = new FileReader(file);
			BufferedReader br = new BufferedReader(fr);
			
			MolBlueprint blueprint = new MolBlueprint();
			
			blueprint.setKey("");
			
			String nextLine = br.readLine();
			while(nextLine != null)
			{
				String[] parts = nextLine.split("\\s+");
				
				if(parts[0].equals("HETATM"))
				{
					blueprint.addToKey(parts[2]);
				}
				else if (parts[0].equals("CONECT"))
				{
					//Build Bond
					if(parts.length == 3)
					{
						//Assume format from file is second, first
						blueprint.addBondIndex(parts[2], parts[1]);
					}
					//Build Angle
					else if(parts.length == 4)
					{
						//Assume format from file is center, left, right
						blueprint.addAngleIndex(parts[2], parts[1], parts[3]);
					}
				}
				nextLine = br.readLine();
			}
			
			blueprints.put(blueprint.getKey(), blueprint);
			
			System.out.println("Loaded molecule: " + blueprint.getKey());
			System.out.println("Bonds: " + blueprint.getBondIndexes().size());
			System.out.println("Angles: " + blueprint.getAngleIndexes().size());
			
			br.close();
			fr.close();
		} catch (FileNotFoundException e) {
			System.err.println("File not found.");
		} catch (IOException e) {
			System.err.println("IO error.");
		}
	}
	
	public void parseSystemFile(String fileName)
	{
		try {
			File file = new File(fileName);
			FileReader fr = new FileReader(file);
			BufferedReader br = new BufferedReader(fr);
			
			String nextLine = br.readLine();
			String molNum = "1";
			ArrayList<Atom> tempAtomCollection = new ArrayList<Atom>();
			while(nextLine != null)
			{
				String[] parts = nextLine.split("\\s+");
				if(parts[0].equals("HETATM"))
				{
					//Check to see if the new atom is part of a new molecule
					if(parts[3].length() != 1)
					{
						String[] fixedParts = parts.clone();
						parts[4] = fixedParts[3].substring(1, fixedParts[3].length());
						parts[5] = fixedParts[4];
						parts[6] = fixedParts[5];
						parts[7] = fixedParts[6];
					}
					if(!molNum.equals(parts[4]))
					{
						String blueprintKey = "";
						for(Atom atom : tempAtomCollection)
						{
							blueprintKey += atom.getType();
						}
						Molecule mol = new Molecule(blueprints.get(blueprintKey));
						for(Atom atom : tempAtomCollection)
						{
							mol.addAtom(atom);
						}
						
						mol.buildBondsAndAngles();
						molecules.add(mol);
						tempAtomCollection = new ArrayList<Atom>();
						molNum = parts[4];
					}
					
					Atom atom = parseAtom(parts);
					tempAtomCollection.add(atom);
				}
				if(nextLine.equals("END"))
				{
					String blueprintKey = "";
					for(Atom atom : tempAtomCollection)
					{
						blueprintKey += atom.getType();
					}
					Molecule mol = new Molecule(blueprints.get(blueprintKey));
					for(Atom atom : tempAtomCollection)
					{
						mol.addAtom(atom);
					}
					
					mol.buildBondsAndAngles();
					molecules.add(mol);
					tempAtomCollection = new ArrayList<Atom>();
				}
				nextLine = br.readLine();
			}
			
			br.close();
			fr.close();
			
		} catch (FileNotFoundException e) {
			System.err.println("File not found.");
		} catch (IOException e) {
			System.err.println("IO Error");
		}
		
	}
	
	/**
	 * Assumes the line is formatted from a Packmol built pdb file.
	 * @param parts
	 * @return	An atom parsed form the line
	 */
	private Atom parseAtom(String[] parts)
	{
		double x,y,z;
		Atom atom = null;
		String atomType;
		int atomID;
		try{
		atomType = parts[2];
		atomID = Integer.parseInt(parts[1]);
		
		x = Double.parseDouble(parts[5]);
		y = Double.parseDouble(parts[6]);
		z = Double.parseDouble(parts[7]);
		
		atom = new Atom(atomID, atomType, x, y, z);
		}catch (NumberFormatException e)
		{
			//The "A" and the molecule number may be mashed together
			if(parts[3].length() != 1)
			{
				atomType = parts[2];
				atomID = Integer.parseInt(parts[1]);
				
				x = Double.parseDouble(parts[4]);
				y = Double.parseDouble(parts[5]);
				z = Double.parseDouble(parts[6]);
				
				atom = new Atom(atomID, atomType, x, y, z);
			}
		}
		return atom;
	}
	
	public ArrayList<Molecule> getMolecules()
	{
		return molecules;
	}
}
