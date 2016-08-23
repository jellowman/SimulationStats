package core;

/**
 * Creates a bond declaration that specifies which atom indexes in a molecule
 * should be bonded together.
 * @author Trevor Fisher
 *
 */
public class BondIndex 
{
	/**Assigns a unique bond ID to each angle type*/
	protected static int newBondID = 1;
	
	/**This bond's assigned ID*/
	protected int bondID;
	
	protected int[] bondIndex;
	
	/**
	 * Creates a new bond association between three specified atom indexes
	 * @param atom1		The first atom index
	 * @param atom2		The second atom index
	 */
	public BondIndex(int atom1, int atom2)
	{
		bondIndex = new int[2];
		bondIndex[0] = atom1-1;
		bondIndex[1] = atom2-1;
		bondID = newBondID;
		newBondID++;
	}
	
	public int getAtom(int i)
	{
		return bondIndex[i];
	}
	
	public int getBondID()
	{
		return bondID;
	}
}
