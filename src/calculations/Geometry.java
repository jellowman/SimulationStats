package calculations;

import core.Atom;

/**
 * Contains basic geometry calculations used for many post-analysis code
 * @author trevorfisher
 *
 */
public class Geometry 
{
	/**
	 * Calculates the distance between two atoms, accounting for periodic
	 * conditions.
	 * @param a		Atom a
	 * @param b		Atom b
	 * @return		The radial distance between the two atoms
	 */
	public static double PeriodicDistance(Atom a, Atom b, double xBox, double yBox, double zBox) {	
		double xDist = b.getXCord() - a.getXCord();
		double yDist = b.getYCord() - a.getYCord();
		double zDist = b.getZCord() - a.getZCord();
		
		xDist = Math.round(xDist/xBox) * xBox - xDist;
		yDist = Math.round(yDist/yBox) * yBox - yDist;
		zDist = Math.round(zDist/zBox) * zBox - zDist;
		double radius = Math.sqrt(xDist*xDist + yDist*yDist + zDist*zDist);
		return radius;
	}
	
	public static double PeriodicDistance(double a, double b, double dimBox) {
		double dist = a-b;
		if(dist > dimBox/2) {
			dist = dimBox-dist;
		}
		return dist;
	}
	
	/**
	 * Calculates the angle for a bond.
	 * @param a - An atom specified as the root of the angle
	 * @param b - An atom as a leaf of the angle
	 * @param c - An atom as a leaf of the angle
	 * @return	The angle (in degrees) of the H-A-B angle
	 */
	public static double PeriodicAngle(Atom a, Atom b, Atom c, double xBox, double yBox, double zBox) {
		//Calculate dot product of a-b bond and a-h bond
		double xDistAB = PeriodicDistance(a.getXCord(), b.getXCord(), xBox);
		double xDistAC = PeriodicDistance(a.getXCord(), c.getXCord(), xBox);
		
		double yDistAB = PeriodicDistance(a.getYCord(), b.getYCord(), yBox);
		double yDistAC = PeriodicDistance(a.getYCord(), c.getYCord(), yBox);
		
		double zDistAB = PeriodicDistance(a.getZCord(), b.getZCord(), zBox);
		double zDistAC = PeriodicDistance(a.getZCord(), c.getZCord(), zBox);
		
		double dotP = xDistAB*xDistAC + yDistAB*yDistAC + zDistAB*zDistAC;
		double distMagn  =   PeriodicDistance(a, b, xBox, yBox, zBox);
		distMagn = distMagn * PeriodicDistance(a, c, xBox, yBox, zBox);
		double cosTh = dotP/distMagn;
		return Math.toDegrees(Math.acos(cosTh));
	}
}
