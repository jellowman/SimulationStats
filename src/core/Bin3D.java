package core;

import java.lang.reflect.Array;
import java.util.ArrayList;

/**
 * Stores objects in bins from periodic 3D Space
 * @author trevorfisher
 *
 * @param <T>
 */
public class Bin3D<T> {
	/**Used to determine offset for bins.*/
	private double minx, miny, minz;
	
	/**Length of each bin in the x, y , and z direction, respectively.*/
	private double dx, dy, dz;
	
	/**Number of bins in the x, y, and z direction, respectively.*/
	private int binsX, binsY, binsZ;
	
	/**3D array containing a list of atoms in the bin*/
	private ArrayList<T>[][][] bin;
	
	/**
	 * Initializes a bin structure for a 3D system
	 * @param minLen	Minimum length for a bin dimension
	 */
	public Bin3D(double minLen, double xmin, double xmax, double ymin,
			double ymax, double zmin, double zmax) {
		minx = xmin;
		miny = ymin;
		minz = zmin;
		
		binsX = (int) ((xmax - xmin)/minLen);
		binsY = (int) ((ymax - ymin)/minLen);
		binsZ = (int) ((zmax - zmin)/minLen);
		
		dx = (xmax-xmin)/binsX;
		dy = (ymax-ymin)/binsY;
		dz = (zmax-zmin)/binsZ;
		
		ArrayList<T> obj = new ArrayList<T>();
		@SuppressWarnings("unchecked")
		ArrayList<T>[][][] a = (ArrayList<T>[][][]) Array.newInstance(obj.getClass(), binsX, binsY, binsZ);
		bin = a;
		for(int i = 0; i < binsX; i++) {
			for(int j = 0; j < binsY; j++) {
				for(int k = 0; k < binsZ; k++){
					bin[i][j][k] = new ArrayList<T>();
				}
			}
		}
	}
	
	/**
	 * Adds object to specified bin
	 * @param obj
	 * @param x
	 * @param y
	 * @param z
	 */
	public void addToBin(T obj, int x, int y, int z) {
		bin[x%binsX][y%binsY][z%binsZ].add(obj);
	}
	
	public void addFrom3Space(T obj, double x, double y, double z) {
		double xOffset = x - minx;
		double yOffset = y - miny;
		double zOffset = z - minz;
		
		int xadd = (int) (xOffset/dx);
		int yadd = (int) (yOffset/dy);
		int zadd = (int) (zOffset/dz);
		
		bin[(xadd+binsX)%binsX][(yadd+binsY)%binsY][(zadd+binsZ)%binsZ].add(obj);
	}
	
	/**
	 * Returns an ArrayList containing all of the neighbors guaranteed to be less than minLen
	 * Also returns other neighbors that could be as far as minLen 2*minLen*sqrt(3)
	 * @param x
	 * @param y
	 * @param z
	 * @return
	 */
	public ArrayList<T> getNeighbors(int x, int y, int z) {
		ArrayList<T> rtn = new ArrayList<T>();
		//Must add objects in the current and 26 surrounding bins
		//System.out.println("Current bin: " + x + "-" + y + "-" + z);
		for(int i = -1; i < 2; i++) {
			for(int j = -1; j < 2; j++) {
				for(int k = -1; k < 2; k++) {
					//System.out.println("\tLooking in bin: " + ((x+i+binsX)%binsX) + "- "+ ((y+j+binsY)%binsY) + "-" + ((z+k+binsZ)%binsZ));
					rtn.addAll(bin[(x+i+binsX)%binsX][(y+j+binsY)%binsY][(z+k+binsZ)%binsZ]);
				}
			}
		}
		return rtn;
	}
	
	public ArrayList<T> getNeighbors(double x, double y, double z) {
		double xOffset = x - minx;
		double yOffset = y - miny;
		double zOffset = z - minz;
		
		int xbin = (int) (xOffset/dx);
		int ybin = (int) (yOffset/dy);
		int zbin = (int) (zOffset/dz);
		return getNeighbors(xbin, ybin, zbin);
	}
	
	@SuppressWarnings("unchecked")
	public void clearBin() {
		bin = (ArrayList<T>[][][]) Array.newInstance(bin.getClass(), binsX, binsY, binsZ);
	}
	
}
