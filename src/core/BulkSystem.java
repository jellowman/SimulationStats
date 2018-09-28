package core;
import java.util.ArrayList;
import java.util.function.Predicate;

/**
 * Holds a collection of atoms in a system with a specified dimension size.
 * @see Atom
 * @author Trevor Fisher
 * @version 1.0
 *
 */
public class BulkSystem extends ArrayList<Atom>
{
	private static final long serialVersionUID = 5667655335676813356L;
	
	/**The dimensions of the system in Cartesian coordinates.*/
	private double xBox, yBox, zBox;
	private double xMax, yMax, zMax;
	private double xMin, yMin, zMin;
	/**The number of atoms in a system.*/
	private int numAtoms;
	
	public BulkSystem() {
		xBox = 0; yBox = 0; zBox = 0;
		numAtoms = 0;
	}
	public void setNumAtoms(int num) {
		numAtoms = num;
	}
	public int getNumAtoms(String type) {
		return numAtoms;
	}
	
	public void setXBox(double x) {
		xBox = x;
	}
	
	public void setYBox(double y) {
		yBox = y;
	}
	
	public void setZBox(double z) {
		zBox = z;
	}

	public double getXBox() {
		return xBox;
	}
	public double getYBox() {
		return yBox;
	}
	public double getZBox() {
		return zBox;
	}
	public void setXMax(double x) {
		xMax = x;
	}
	public void setXMin(double x) {
		xMin = x;
	}
	public void setYMax(double y) {
		yMax = y;
	}
	public void setYMin(double y) {
		yMin = y;
	}
	public void setZMax(double z) {
		zMax = z;
	}
	public void setZMin(double z) {
		zMin = z;
	}
	public double getXMax() {
		return xMax;
	}
	public double getXMin() {
		return xMin;
	}
	public double getYMax() {
		return yMax;
	}
	public double getYMin() {
		return yMin;
	}
	public double getZMax() {
		return zMax;
	}
	public double getZMin() {
		return zMin;
	}
	
	/**
	 * Predicate for filtering atoms with a whitelist specification.
	 * @param <T>	An <code>Atom</code> derived class
	 */
	private class AtomFilterPredicate<T extends Atom> implements Predicate<T>
	{
		/**The desired atom type to keep*/
		String atomType;
		
		public AtomFilterPredicate(String type){
			atomType = type;
		}
		@Override
		public boolean test(Atom t) {
			if(t.getType().equals(atomType)){
				return false;
			} else{
				return true;
			}
		}
	}
	/**
	 * Removes any atoms in the system that do not match the atom identifier given.
	 * @param type	The desired atom to keep.
	 */
	public void filterAtomType(String type)
	{
		AtomFilterPredicate<Atom> typeFilter = new AtomFilterPredicate<Atom>(type);
		this.removeIf(typeFilter);
	}
	
	/**
	 * Predicate for filtering atoms with a blacklist specification.
	 * @param <T>	An <code>Atom</code> derived class
	 */
	private class AtomRemovePredicate<T extends Atom> implements Predicate<T>
	{
		/**The atom type to remove from the system.*/
		String atomType;
		
		public AtomRemovePredicate(String type){
			atomType = type;
		}
		@Override
		public boolean test(Atom t) {
			if(t.getType().equals(atomType)){
				return true;
			} else{
				return false;
			}
		}
	}
	/**
	 * Removes any atoms in the system that match the atom identifier given.
	 * @param type	The atom type to remove from the system.
	 */
	public void removeAtomType(String type)
	{
		AtomRemovePredicate<Atom> typeFilter = new AtomRemovePredicate<Atom>(type);
		this.removeIf(typeFilter);
	}
}
