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
				else if (parts[0].contains("CONECT"))
				{
					//Assign same-type connection ID
					int connecID;
					try {
					connecID = Integer.parseInt(parts[0].substring(6));
					} catch (NumberFormatException e) {
						connecID = 0;
					}
					
					//Build Bond
					if(parts.length == 3)
					{
						//Assume format from file is first, second
						blueprint.addBondIndex(connecID, parts[1], parts[2]);
					}
					//Build Angle
					else if(parts.length == 4)
					{
						//Assume format from file is left, center, right; need to insert as left, center, right
						blueprint.addAngleIndex(connecID, parts[1], parts[2], parts[3]);
					}
					//Build Dihedral
					else if(parts.length == 5) {
						blueprint.addDihedralIndex(connecID, parts[1], parts[2], parts[3], parts[4]);
					}
				}
				nextLine = br.readLine();
			}
			
			blueprints.put(blueprint.getKey(), blueprint);
			
			System.out.println("Loaded molecule: " + blueprint.getKey());
			System.out.println("Bonds: " + blueprint.getBondIndexes().size());
			System.out.println("Angles: " + blueprint.getAngleIndexes().size());
			System.out.println("Dihedrals: " + blueprint.getDihedralIndexes().size());
			
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
				
				if(parts[0].contains("HETATM"))
				{
					parts = fixParts(parts);
					
					//Check to see if the new atom is part of a new molecule
					if(!molNum.equals(parts[5]))
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
						
						mol.buildConnections();
						molecules.add(mol);
						tempAtomCollection = new ArrayList<Atom>();
						molNum = parts[5];
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
					
					mol.buildConnections();
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
		//try{
		atomType = parts[2];
		atomID = Integer.parseInt(parts[1]);
		
		x = Double.parseDouble(parts[parts.length-3]);
		y = Double.parseDouble(parts[parts.length-2]);
		z = Double.parseDouble(parts[parts.length-1]);
		
		atom = new Atom(atomID, atomType, x, y, z);
		/*}catch (NumberFormatException e)
		{
			//The "A" and the molecule number may be mashed together
			if(parts[3].length() != 1)
			{
				atomType = parts[2];
				atomID = Integer.parseInt(parts[1]);
				
				x = Double.parseDouble(parts[parts.length-3]);
				y = Double.parseDouble(parts[parts.length-2]);
				z = Double.parseDouble(parts[parts.length-1]);
				
				atom = new Atom(atomID, atomType, x, y, z);
			}
		}
		catch (ArrayIndexOutOfBoundsException e)
		{
			atomType = parts[2];
			atomID = Integer.parseInt(parts[1]);
			
			x = Double.parseDouble(parts[parts.length-3]);
			y = Double.parseDouble(parts[parts.length-2]);
			z = Double.parseDouble(parts[parts.length-1]);
		}*/
		return atom;
	}
	
	private String[] fixParts(String[] parts)
	{
		//Check to see if the ID was concatenated to the beginning
		if(!parts[0].equals("HETATM"))
		{
			String[] tempParts = new String[parts.length+1];
			tempParts[0] = "HETATM";
			tempParts[1] = parts[0].replace("HETATM", "");
			for(int i = 2; i < tempParts.length; i++)
			{
				tempParts[i] = parts[i-1];
			}
			parts = tempParts;
		}
		
		//Check to see if the molecule ID and the Atom type are concatenated
		if(parts[4].length() != 1)
		{
			String[] tempParts = new String[parts.length+1];
			for(int i = 0; i < 4; i++)
			{
				tempParts[i] = parts[i];
			}
			tempParts[4] = parts[4].substring(0, 1);
			tempParts[5] = parts[4].substring(1);
			for(int i = 6; i < tempParts.length; i++)
			{
				tempParts[i] = parts[i-1];
			}
			parts = tempParts;
		}
		return parts;
	}
	
	public ArrayList<Molecule> getMolecules()
	{
		return molecules;
	}
}
