package core;

/**
 * Creates an angle declaration that specifies which atom indexes in a molecule
 * should be bonded together.
 * @see core.MolBlueprint
 * @author Trevor Fisher
 *
 */
public class AngleIndex 
{
	/**Assigns a unique angle ID to each angle type*/
	protected static int newAngleID = 1;
	
	/**This angle's assigned ID*/
	protected int angleID;
	
	protected int[] angleIndex;
	
	/**
	 * Creates a new angle association between the three specified atom indexes
	 * @param atom1		One of the atom indexes at the end of the angle
	 * @param atom2		The central atom index that forms the origin of the angle
	 * @param atom3		One of the atom indexes at the end of the angle
	 */
	public AngleIndex(int atom1, int atom2, int atom3)
	{
		angleIndex = new int[3];
		angleIndex[0] = atom1-1;
		angleIndex[1] = atom2-1;
		angleIndex[2] = atom3-1;
		angleID = newAngleID;
		newAngleID++;
	}
	
	public int getAtom(int i)
	{
		return angleIndex[i];
	}
	
	public int getAngleID()
	{
		return angleID;
	}
}
