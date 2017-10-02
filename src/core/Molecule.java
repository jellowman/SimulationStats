package core;

import java.util.ArrayList;

public class Molecule 
{
	/**Contains information about relative bond and angle connections*/
	protected MolBlueprint blueprint;
	
	/**All of the atoms in the molecule*/
	protected ArrayList<Atom> atoms;
	
	protected ArrayList<Bond> bonds;
	protected ArrayList<Angle> angles;
	protected ArrayList<Dihedral> dihedrals;
	
	public Molecule(MolBlueprint blueprint)
	{
		this.blueprint = blueprint;
		atoms = new ArrayList<Atom>();
		bonds = new ArrayList<Bond>();
		angles = new ArrayList<Angle>();
		dihedrals = new ArrayList<Dihedral>();
	}
	
	public void addAtom(Atom atom)
	{
		atoms.add(atom);
	}
	
	public class Bond
	{
		private Atom atom1;
		private Atom atom2;
		private int bondID;
		
		public Bond(int bondID, Atom atom1, Atom atom2)
		{
			this.bondID = bondID;
			this.atom1 = atom1;
			this.atom2 = atom2;
		}
		
		public String toString()
		{
			String formatted = String.format("%12d %12d %12d", bondID, atom1.getID(), atom2.getID());
			return formatted;
		}
		
		public int getID()
		{
			return bondID;
		}
	}
	
	public class Angle
	{
		private Atom atom1, atom2, atom3;
		private int angleID;
		
		public Angle(int angleID, Atom atom1, Atom atom2, Atom atom3)
		{
			this.angleID = angleID;
			this.atom1 = atom1;
			this.atom2 = atom2;
			this.atom3 = atom3;
		}
		
		public String toString()
		{
			String formatted = String.format("%12d %12d %12d %12d", angleID, atom1.getID(), atom2.getID(), atom3.getID());
			return formatted;
		}
		
		public int getID()
		{
			return angleID;
		}
	}
	
	public class Dihedral
	{
		private Atom atom1, atom2, atom3, atom4;
		private int dihedralID;
		
		public Dihedral(int dihedralID, Atom atom1, Atom atom2, Atom atom3, Atom atom4)
		{
			this.dihedralID = dihedralID;
			this.atom1 = atom1;
			this.atom2 = atom2;
			this.atom3 = atom3;
		}
		
		public String toString()
		{
			String formatted = String.format("%12d %12d %12d %12d", dihedralID, atom1.getID(), atom2.getID(), atom3.getID(), atom4.getID());
			return formatted;
		}
		
		public int getID()
		{
			return dihedralID;
		}
	}

	public void buildConnections() 
	{
		try{
		for(BondIndex bondIndex : blueprint.getBondIndexes())
		{
			Bond bond = new Bond(bondIndex.getBondID(), atoms.get(bondIndex.getAtom(0)), atoms.get(bondIndex.getAtom(1)));
			bonds.add(bond);
		}
		
		for(AngleIndex angleIndex : blueprint.getAngleIndexes())
		{
			Angle angle = new Angle(angleIndex.getAngleID(), atoms.get(angleIndex.getAtom(0)),
					atoms.get(angleIndex.getAtom(1)), atoms.get(angleIndex.getAtom(2)));
			angles.add(angle);
		}
		
		for(DihedralIndex dihedralIndex : blueprint.getDihedralIndexes())
		{
			Dihedral dihedral = new Dihedral(dihedralIndex.getDihedralID(), atoms.get(dihedralIndex.getAtom(0)),
					atoms.get(dihedralIndex.getAtom(1)), atoms.get(dihedralIndex.getAtom(2)), atoms.get(dihedralIndex.getAtom(3)));
			dihedrals.add(dihedral);
		}
		}
		catch(NullPointerException ne)
		{
			String molID = "";
			for(Atom atom : atoms)
			{
				molID += atom.getType();
			}
			System.err.println("Could not find a blueprint of this molecule: " + molID + "\n"
					+ "Ignoring bonds and angles.");
		}
	}
	
	public ArrayList<Atom> getAtoms()
	{
		return atoms;
	}
	public ArrayList<Bond> getBonds()
	{
		return bonds;
	}
	public ArrayList<Angle> getAngles()
	{
		return angles;
	}
	public ArrayList<Dihedral> getDihedrals()
	{
		return dihedrals;
	}
	
	public String toString()
	{
		return "Atoms: " + atoms.size() + " Bonds: " + bonds.size() + " Angles: " + angles.size() + " Dihedrals: " + dihedrals.size();
		
	}
}
