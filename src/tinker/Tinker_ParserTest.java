package tinker;
import static org.junit.Assert.*;

import org.junit.Test;

import core.Atom;
import core.BulkSystem;

public class Tinker_ParserTest {

	@Test
	public void testParser() {
		String numberLine = "     1  OW    -6.191289   -2.188269    4.647389   164     2     3";
		String[] cutLine = numberLine.split(" ");
		String[] parts = new String[50];
		int j = 0;
		for(int i = 0; i < cutLine.length; i++)
		{
			if(cutLine[i].length() != 0)
			{
				parts[j] = cutLine[i];
				j++;
			}
		}
		
		for(String part : parts)
		{
			System.out.println(part);
		}
		
		String type = parts[1];
		double xCord = Double.valueOf(parts[2]);
		double yCord = Double.valueOf(parts[3]);
		double zCord = Double.valueOf(parts[4]);
		Atom atom = new Atom(type, xCord, yCord, zCord);
		
		assertEquals(atom.getType().equals("OW"), true);
		assertEquals(atom.getXCord() == -6.191289, true);
		assertEquals(atom.getYCord() == -2.188269, true);
		assertEquals(atom.getZCord() == 4.647389, true);
	}
	
	@Test
	public void testMethod() throws BadAtomFormat
	{
		String numberLine = "     1  OW    -6.191289   -2.188269    4.647389   164     2     3";
		BulkSystem bs = new BulkSystem();
		Tinker_Parser.parseAtom(numberLine, bs);
	}

}
