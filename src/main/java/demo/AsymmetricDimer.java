package demo;

import org.biojava.bio.structure.Atom;
import org.biojava.bio.structure.ResidueNumber;
import org.biojava.bio.structure.Structure;
import org.biojava.bio.structure.StructureException;
import org.biojava.bio.structure.StructureTools;
import org.biojava.bio.structure.align.StructureAlignment;
import org.biojava.bio.structure.align.StructureAlignmentFactory;
import org.biojava.bio.structure.align.ce.CeMain;
import org.biojava.bio.structure.align.gui.StructureAlignmentDisplay;
import org.biojava.bio.structure.align.model.AFPChain;
import org.biojava.bio.structure.align.util.AFPAlignmentDisplay;
import org.biojava.bio.structure.align.util.AlignmentTools;

public class AsymmetricDimer {

	public static void main(String[] args) {
		// /Users/blivens/dev/capitani/tlpa_superpositions/tlpa_scoI_on_cox_by_tlpa.pdb
		if( args.length != 1 || args[0].equalsIgnoreCase("-h")) {
			System.out.println("usage: [-h] structure");
			System.exit(1);
			return;
		}
		String name = args[0];
		Structure s;
		try {
			s = StructureTools.getStructure(name);
		} catch (Exception e) {
			System.err.println("Error loading structure "+name);
			System.exit(1);
			return;
		}
		
		if(s.size() != 2) {
			System.err.println("Error. Only runs on dimers.");
		}
		
		Atom[] ca1 = StructureTools.getAtomCAArray(s.getChain(0));
		Atom[] ca2 = StructureTools.getAtomCAArray(s.getChain(1));
		
		AFPChain chainAlign;
		try {
			StructureAlignment ce = StructureAlignmentFactory.getAlgorithm(CeMain.algorithmName);
			chainAlign = ce.align(ca1, ca2);
		} catch (StructureException e) {
			e.printStackTrace();
			System.exit(1);
			return;
		}
		
		Atom[] full1 = StructureTools.getAtomCAArray(s);
		
		ResidueNumber[] aligned1 = new ResidueNumber[2*chainAlign.getOptLength()];
		ResidueNumber[] aligned2 = new ResidueNumber[aligned1.length];
		int alignedPos = 0;
		
		int[] blocklens = chainAlign.getOptLen();
		int[][][] optAln = chainAlign.getOptAln();
		for(int block=0;block<chainAlign.getBlockNum();block++) {
			for(int pos=0;pos<blocklens[block];pos++) {
				int atom1 = optAln[block][0][pos];
				int atom2 = optAln[block][1][pos];
				ResidueNumber res1 = ca1[ atom1 ].getGroup().getResidueNumber();
				ResidueNumber res2 = ca2[ atom2 ].getGroup().getResidueNumber();
				aligned1[alignedPos] = res1;
				aligned2[alignedPos] = res2;
				aligned1[chainAlign.getOptLength() + alignedPos ] = res2;
				aligned2[chainAlign.getOptLength() + alignedPos ] = res1;
				
				alignedPos++;
			}
		}
		try {
			AFPChain fullAlignment = AlignmentTools.createAFPChain(ca1, ca2, aligned1, aligned2);
			
			Atom[] full2 = StructureTools.cloneCAArray(full1);
			StructureAlignmentDisplay.display(fullAlignment, full1, full2);
		} catch (StructureException e) {
			e.printStackTrace();
			System.exit(1);
			return;
		}
	}

}
