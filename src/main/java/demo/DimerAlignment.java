package demo;

import java.util.LinkedHashMap;
import java.util.Map;

import org.biojava.bio.structure.Atom;
import org.biojava.bio.structure.Calc;
import org.biojava.bio.structure.Structure;
import org.biojava.bio.structure.StructureException;
import org.biojava.bio.structure.StructureTools;
import org.biojava.bio.structure.align.StructureAlignment;
import org.biojava.bio.structure.align.StructureAlignmentFactory;
import org.biojava.bio.structure.align.ce.CeMain;
import org.biojava.bio.structure.align.gui.StructureAlignmentDisplay;
import org.biojava.bio.structure.align.gui.jmol.StructureAlignmentJmol;
import org.biojava.bio.structure.align.model.AFPChain;
import org.biojava.bio.structure.align.util.AlignmentTools;
import org.biojava.bio.structure.align.util.RotationAxis;
import org.biojava.bio.structure.jama.Matrix;

/**
 * We expect homodimeric structures to be symmetric. Extending this to
 * heterodimers with related structures, we expect them to be pseudo-symmetric.
 *
 * <p>This script takes a two-chain structure as input. It aligns the two
 * chains individually, then creates a rigid-body superposition of the full
 * complex using this alignment.
 * 
 * <p>The script also displays a rotation axis on one of the chains, showing
 * how that chain is rotated away from its optimal alignment by the constraint
 * of the other domain.
 * @author Spencer Bliven
 *
 */
public class DimerAlignment {
	private AFPChain alignment;
	private RotationAxis axis;
	private Atom[] full1;
	private Atom[] full2;



	/**
	 * 
	 * @param aligner
	 * @param ca Concatenated array of CA atoms for both chains. The first `chain1len`
	 *  atoms should come from one chain, with the remainder coming from the second.
	 * @param chain1len Number of atoms in `ca` from the first chain
	 * @return
	 * @throws StructureException
	 */
	public DimerAlignment(StructureAlignment aligner, Structure s) throws StructureException{
		if(s.size() != 2) {
			throw new IllegalArgumentException("Expect two-chain input structure");
		}

		Atom[] ca1 = StructureTools.getAtomCAArray(s.getChain(0));
		Atom[] ca2 = StructureTools.getAtomCAArray(s.getChain(1));

		// Align chains individually
		AFPChain chainAlign = aligner.align(ca1, ca2);

		Matrix chainRotation = chainAlign.getBlockRotationMatrix()[0];
		Atom chainShift = chainAlign.getBlockShiftVector()[0];
		Matrix inverseChainRotation = chainRotation.inverse();

		// Create new alignment with the same aligned residues, but include the full structure

		// Build map of aligned positions
		Map<Integer,Integer> alignmentMap = new LinkedHashMap<Integer,Integer>();
		int[] blocklens = chainAlign.getOptLen();
		int[][][] optAln = chainAlign.getOptAln();
		for(int block=0;block<chainAlign.getBlockNum();block++) {
			for(int pos=0;pos<blocklens[block];pos++) {
				int atom1 = optAln[block][0][pos];
				int atom2 = optAln[block][1][pos] + ca1.length;
				alignmentMap.put(atom1, atom2);
				alignmentMap.put(atom2, atom1);
			}
		}
		// Create an AFPChain for the full alignment
		this.full1 = StructureTools.getAtomCAArray(s);
		this.full2 = StructureTools.cloneCAArray(full1);
		this.alignment = AlignmentTools.replaceOptAln(chainAlign, full1, full2, alignmentMap);
		
		Matrix fullRotation = alignment.getBlockRotationMatrix()[0];
		Atom fullShift = alignment.getBlockShiftVector()[0];

		// matrices are premultiplied
		// rotation from chainAlign to fullAlign is (x-chainShift)*chainRot^-1*fullRot+fullShift
		// or x*combinedMatrix + (fullShift-chainShift*combinedMatrix)
		Matrix combinedMatrix = inverseChainRotation.times(fullRotation);
		Atom combinedShift = chainShift;
		Calc.rotate(combinedShift, combinedMatrix);
		combinedShift = Calc.subtract(fullShift, combinedShift);
		
		this.axis = new RotationAxis(combinedMatrix, combinedShift);
	}


	public AFPChain getAlignment() {
		return alignment;
	}

	public RotationAxis getAxis() {
		return axis;
	}

	public Atom[] getCA1() {
		return full1;
	}

	public Atom[] getCA2() {
		return full2;
	}

	public static void main(String[] args) {
		// Read structure from input. May be a PDB ID, filename, SCOP domain, etc.
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


		// Only accept dimer input
		if(s.size() != 2) {
			System.err.println("Error. Only runs on dimers.");
		}


		try {
			StructureAlignment ce = StructureAlignmentFactory.getAlgorithm(CeMain.algorithmName);
			DimerAlignment dimer = new DimerAlignment(ce,s);
			Atom[] full1 = dimer.getCA1();
			Atom[] full2 = dimer.getCA2();
			AFPChain fullAlign = dimer.getAlignment();
			// Display
			StructureAlignmentJmol jmol = StructureAlignmentDisplay.display(fullAlign, full1, full2);
			
			RotationAxis axis = dimer.getAxis();
			jmol.evalString(axis.getJmolScript(full1));
		} catch (StructureException e) {
			e.printStackTrace();
			System.exit(1);
			return;
		}
	}
}
