package demo;

import java.util.LinkedHashMap;
import java.util.Map;

import org.biojava.bio.structure.Atom;
import org.biojava.bio.structure.Structure;
import org.biojava.bio.structure.StructureException;
import org.biojava.bio.structure.StructureTools;
import org.biojava.bio.structure.align.StructureAlignment;
import org.biojava.bio.structure.align.StructureAlignmentFactory;
import org.biojava.bio.structure.align.ce.CeMain;
import org.biojava.bio.structure.align.gui.StructureAlignmentDisplay;
import org.biojava.bio.structure.align.model.AFPChain;
import org.biojava.bio.structure.align.util.AlignmentTools;

/**
 * We expect homodimeric structures to be symmetric. Extending this to
 * heterodimers with related structures, we expect them to be pseudo-symmetric.
 *
 * <p>This script takes a two-chain structure as input. It aligns the two
 * chains individually, then creates a rigid-body superposition of the full
 * complex using this alignment.
 * @author Spencer Bliven
 *
 */
public class AsymmetricDimer {

	public static void main(String[] args) {
		// /Users/blivens/dev/capitani/tlpa_superpositions/tlpa_scoI_on_cox_by_tlpa.pdb
		if( args.length != 1 || args[0].equalsIgnoreCase("-h")) {
			System.out.println("usage: [-h] structureID");
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
		assert(full1.length == ca1.length+ca2.length);
		
		int chain1len = ca1.length;
		
		Map<Integer,Integer> alignmentMap = new LinkedHashMap<Integer,Integer>();
		
		int[] blocklens = chainAlign.getOptLen();
		int[][][] optAln = chainAlign.getOptAln();
		for(int block=0;block<chainAlign.getBlockNum();block++) {
			for(int pos=0;pos<blocklens[block];pos++) {
				int atom1 = optAln[block][0][pos];
				int atom2 = optAln[block][1][pos] + chain1len;
				alignmentMap.put(atom1, atom2);
				alignmentMap.put(atom2, atom1);
			}
		}
		try {
			Atom[] full2 = StructureTools.cloneCAArray(full1);
			AFPChain fullAlign = AlignmentTools.replaceOptAln(chainAlign, full1, full2, alignmentMap);
			AlignmentTools.updateSuperposition(fullAlign, full1, full2);
			StructureAlignmentDisplay.display(fullAlign, full1, full2);
		} catch (StructureException e) {
			e.printStackTrace();
			System.exit(1);
			return;
		}
	}

}
