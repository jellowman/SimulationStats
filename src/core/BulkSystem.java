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
	/**The number of atoms in a system.*/
	private int numAtoms;
	
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
