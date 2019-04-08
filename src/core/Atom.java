package core;

/**
 * Holds all relevant information contained in an atom. As of version 1.0, quantum properties
 * are not implemented.
 * @author Trevor Fisher
 * @version 1.0
 */
public class Atom 
{
	/**The 3D Cartesian coordinates of the atom in the system.*/
	private double xCord, yCord, zCord;
	
	/**The atom type, used for identification in some calculations. */
	private String type;
	
	/**The numeric ID assigned by the simulation system*/
	private int iD;
	
	/**The molecular ID that this atom is associated with*/
	private int molID;
	
	/**Indicates if the atom is in an H-Bond */
	private boolean isHBonded;
	
	/**
	 * Initializes the atom with coordinates, but will give it a default name of "None".
	 * @param x		The x-coordinate of the atom
	 * @param y		The y-coordinate of the atom
	 * @param z		The z-coordinate of the atom
	 */
	public Atom(double x, double y, double z)
	{
		this(0, "None", x, y, z);
	}
	
	/**
	 * Initializes the atom with the specified identifier and coordinates.
	 * @param type	The identifier for the atom
	 * @param x		The x-coordinate of the atom
	 * @param y		The y-coordinate of the atom
	 * @param z		The z-coordinate of the atom
	 */
	public Atom(String type, double x, double y, double z)
	{
		this(0, type, x, y, z);
	}
	
	/**
	 * Initializes the atom with the specified identifier and coordinates.
	 * @param iD	The numeric number assigned by the simulation system
	 * @param type	The identifier for the atom
	 * @param x		The x-coordinate of the atom
	 * @param y		The y-coordinate of the atom
	 * @param z		The z-coordinate of the atom
	 */
	public Atom(int iD, String type, double x, double y, double z)
	{
		xCord = x;
		yCord = y;
		zCord = z;
		this.type = type;
		this.iD = iD;
	}
	
	/**
	 * Gets the atom identifier
	 * @return		The atom identifier
	 */
	public String getType()
	{
		return this.type;
	}
	
	public double getXCord()
	{
		return xCord;
	}
	public double getYCord()
	{
		return yCord;
	}
	public double getZCord()
	{
		return zCord;
	}
	/**Returns a 0-indexed id of an atom**/
	public int getID()
	{
		return iD;
	}
	public int getMolID() {
		return molID;
	}
	public boolean isHBonded() {
		return isHBonded;
	}
	public void setMolID(int mid) {
		molID = mid;
	}
	public void setXCord(double x) {
		xCord = x;
	}
	public void setYCord(double y) {
		yCord = y;
	}
	public void setZCord(double z) {
		zCord = z;
	}
	public void setCords(double x, double y, double z) {
		xCord = x; yCord = y; zCord = z;
	}
	public void setisHBonded(boolean hBondStatus) {
		isHBonded = hBondStatus;
	}
	
	public String toString()
	{
		String formatted = String.format("%10d %10d %10d", xCord, yCord, zCord);
		return formatted;
	}
}
