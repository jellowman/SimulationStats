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
	public int getID()
	{
		return iD;
	}
	
	public String toString()
	{
		String formatted = String.format("%10d %10d %10d", xCord, yCord, zCord);
		return formatted;
	}
}
