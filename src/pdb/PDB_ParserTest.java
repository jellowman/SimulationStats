package pdb;

import static org.junit.Assert.*;

import org.junit.Test;

public class PDB_ParserTest {

	@Test
	public void test() {
		String one = "HETATM    1  O           0      -1.176  -0.169   0.000                       O";
		String[] parts = one.split("\\s+");
		assertEquals(parts.length == 8, true);
		for(String part : parts)
		{
			System.out.println(part);
		}
		
		String two = "CONECT    1    2    3";
		String three = "CONECT    2    1";
		
		parts = two.split("\\s+");
		System.out.println(parts.length);
		
		parts = three.split("\\s+");
		System.out.println(parts.length);
	}
}
