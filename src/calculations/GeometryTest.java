package calculations;

import static org.junit.Assert.*;

import org.junit.Test;

import core.Atom;

public class GeometryTest {

	@Test
	public void testAngle() {
		Atom a1 = new Atom("O", 0.0, 0.0, 0.0);
		Atom a2 = new Atom("H", 1.732, 0.0, 0.0);
		Atom a3 = new Atom("O", 1.732, 1.0, 0.0);
		
		double angle = Geometry.PeriodicAngle(a1, a2, a3, 10, 10, 10);
		System.out.println("Angle: " + angle);
		assertEquals(angle, 30, 2);
	}
	
	@Test
	public void testDistance() {
		double x1 = 1.0;
		double x2 = 8.0;
		double sbox = -1;
		double lbox = 11;
		
		double dist = Geometry.PeriodicDistance(x1, x2, sbox, lbox);
		assertEquals(dist, -5, 0.1);
	}
	
	

}
