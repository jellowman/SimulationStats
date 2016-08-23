package core;

import org.junit.Test;

public class BondIndexTest {

	@Test
	public void test() {
		
		String formatted = String.format("%-10s %10d %10f", "osne", 12098531, -135364.657);
		String formatted2 = String.format("%-10s %10d %10d", "fouwer", -351, 137365);
		System.out.println(formatted);
		System.out.println(formatted2);
	}

}
