package core;

import java.util.ArrayList;

/**
 * Holds information on how to create bonds and angles between atoms
 * in a molecule.
 * @author Trevor Fisher
 */
public class MolBlueprint 
{
	protected ArrayList<BondIndex> bondIndexes;
	protected ArrayList<AngleIndex> angleIndexes;
	protected ArrayList<DihedralIndex> dihedralIndexes;
	String key;
	/*
	 * Will be implemented at a later time
	 * protected Dihedral[] dihedrals;
	 */
	
	public MolBlueprint()
	{
		bondIndexes = new ArrayList<BondIndex>();
		angleIndexes = new ArrayList<AngleIndex>();
		dihedralIndexes = new ArrayList<DihedralIndex>();
	}
	
	public void addBondIndex(int connecID, String atom1, String atom2)
	{
		int a1 = Integer.parseInt(atom1);
		int a2 = Integer.parseInt(atom2);
		bondIndexes.add(new BondIndex(connecID, a1, a2));
	}
	
	public void addAngleIndex(int connecID, String atom1, String atom2, String atom3)
	{
		int a1 = Integer.parseInt(atom1);
		int a2 = Integer.parseInt(atom2);
		int a3 = Integer.parseInt(atom3);
		angleIndexes.add(new AngleIndex(connecID, a1, a2, a3));
	}
	
	public void addDihedralIndex(int connecID, String atom1, String atom2, String atom3, String atom4) {
		int a1 = Integer.parseInt(atom1);
		int a2 = Integer.parseInt(atom2);
		int a3 = Integer.parseInt(atom3);
		int a4 = Integer.parseInt(atom4);
		dihedralIndexes.add(new DihedralIndex(connecID, a1, a2, a3, a4));
	}
	
	public ArrayList<BondIndex> getBondIndexes()
	{
		return bondIndexes;
	}
	
	public ArrayList<AngleIndex> getAngleIndexes()
	{
		return angleIndexes;
	}
	
	public ArrayList<DihedralIndex> getDihedralIndexes() {
		return dihedralIndexes;
	}
	
	public void addToKey(String add)
	{
		key += add;
	}
	
	public void setKey(String key)
	{
		this.key = key;
	}
	
	public String getKey()
	{
		return key;
	}
}
