package core;

/**
 * Creates an angle declaration that specifies which atom indexes in a molecule
 * should be bonded together.
 * @see core.MolBlueprint
 * @author Trevor Fisher
 *
 */
public class DihedralIndex 
{
	/**Assigns a unique angle ID to each angle type*/
	protected static int newDihedralID = 1;
	
	/**This angle's assigned ID*/
	protected int dihedralID;
	
	protected int[] dihedralIndex;
	
	/**
	 * Creates a new angle association between the four specified atom indexes. Blueprint molecule is
	 * one-indexed, but the internal arrays are 0-indexed. Thus, a -1 shift occurs in assignment.
	 * @param atom1		One of the atom indexes at the end of the dihedral
	 * @param atom2		A central atom in the dihedral, attached to atom1
	 * @param atom3		A central atom in the dihedral, attached to atom4
	 * @param atom4		One of the atom indexes at the end of the dihedral
	 */
	public DihedralIndex(int atom1, int atom2, int atom3, int atom4)
	{
		dihedralIndex = new int[4];
		dihedralIndex[0] = atom1-1;
		dihedralIndex[1] = atom2-1;
		dihedralIndex[2] = atom3-1;
		dihedralIndex[3] = atom4-1;
		dihedralID = newDihedralID;
		newDihedralID++;
	}
	
	public int getAtom(int i)
	{
		return dihedralIndex[i];
	}
	
	public int getDihedralID()
	{
		return dihedralID;
	}
}
