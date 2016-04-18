package srp.evolution.operators;

import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.List;

import beast.core.Description;
import beast.core.Input;
import beast.core.Operator;
import beast.evolution.tree.Node;
import beast.evolution.tree.Tree;
import beast.util.Randomizer;
import srp.evolution.OperationType;
import srp.evolution.haplotypes.Haplotype;
import srp.evolution.haplotypes.HaplotypeModel;

@Description("Swap singel base")
public class HaplotypeBaseSingleOperator extends Operator {

	final public Input<HaplotypeModel> haplotypeModelInput = new Input<>("hapModel", 
			"haplotype model");
	final public Input<Boolean> gaussianInput = new Input<>("gaussian",
			"Gaussian (=true=default) or uniform delta", true);
	final public Input<Boolean> optimiseInput = new Input<>("optimise",
			"flag to indicate that the scale factor is automatically changed in order to achieve a good acceptance rate (default true)",
			true);
	final public Input<Double> limitInput = new Input<>("limit",
			"limit on step size, default disable, "
					+ "i.e. -1. (when positive, gets multiplied by tree-height/log2(n-taxa).",
			-1.0);
	// shadows size
	public static final char[] DNA_CHARS = {'A','C','G','T'};
	
	public static final OperationType OP = OperationType.SINGLE;
	double size;
	private double limit;
	private HaplotypeModel haplotypeModel;

	@Override
	public void initAndValidate() {
		haplotypeModel = haplotypeModelInput.get();
		haplotypeModel.setID("aoeu");
//		limit = limitInput.get();
	}

	/**
	 * Do a probabilistic subtree slide move.
	 *
	 * @return log of Hastings Ratio, or Double.NEGATIVE_INFINITY if proposal
	 *         should not be accepted *
	 */
	@Override
	public double proposal() {
//		System.out.println("Start proposal");
		HaplotypeModel hap = haplotypeModelInput.get(this); // call store() here, how to pass record??
//		System.out.println(hap.getID());
		hap.setOperationType(OP);
		
		
		
		double logq = 0;

		int hapIndex = getNextHapIndex();
		int siteIndex = getNextSiteIndex();
//		int hapIndex = 1;
//		int siteIndex = 1;

		siteIndex %= 4;
		Haplotype haplotype = hap.getHaplotype(hapIndex);
		int oldState = haplotype.getState(siteIndex);
		int newState = getNextState();
//		char newChar = 'A';
//		newChar = '-';
		
		hap.updateState(hapIndex, siteIndex, newState);

		haplotype.setStateAt(siteIndex, newState);
		hap.setOperationRecord(OP, hapIndex, siteIndex);

//		hap.startEditing(this);
		hap.storeState();// call explicitely, should be able to change this
		return logq;
	}


	public int getNextHapIndex() {
		return Randomizer.nextInt(haplotypeModel.getHaplotypeCount());
	}

	public int getNextSiteIndex() {
		return getNextSiteIndex(haplotypeModel.getHaplotypeLength());
	}
	public int getNextSiteIndex(int length){
		return Randomizer.nextInt(length);
	}
	
	public OperationType getOperationType(){
		return OperationType.SINGLE;
	};

	public static char getNextDiffBase(int oldState) {
		int i = oldState;
		do {
			i = Randomizer.nextInt(4);
		} while (i == oldState);
	
		return DNA_CHARS[i];
	}

	public static int getNextState() {
		int i = Randomizer.nextInt(4);
		return i;
	}
	// /**
	// * automatic parameter tuning *
	// */
	// @Override
	// public void optimize(final double logAlpha) {
	// if (optimiseInput.get()) {
	// double delta = calcDelta(logAlpha);
	// delta += Math.log(size);
	// final double f = Math.exp(delta);
	//// double f = Math.exp(delta);
	// if( limit > 0 ) {
	// final Tree tree = treeInput.get();
	// final double h = tree.getRoot().getHeight();
	// final double k = Math.log(tree.getLeafNodeCount()) / Math.log(2);
	// final double lim = (h / k) * limit;
	// if( f <= lim ) {
	// size = f;
	// }
	// } else {
	// size = f;
	// }
	// }
	// }

	// @Override
	// public double getCoercableParameterValue() {
	// return size;
	// }
	//
	// @Override
	// public void setCoercableParameterValue(final double value) {
	// size = value;
	// }
	//
	// @Override
	// public String getPerformanceSuggestion() {
	// final double prob = m_nNrAccepted / (m_nNrAccepted + m_nNrRejected +
	// 0.0);
	// final double targetProb = getTargetAcceptanceProbability();
	//
	// double ratio = prob / targetProb;
	//
	// if (ratio > 2.0) ratio = 2.0;
	// if (ratio < 0.5) ratio = 0.5;
	//
	// final double newDelta = size * ratio;
	//
	// final DecimalFormat formatter = new DecimalFormat("#.###");
	// if (prob < 0.10) {
	// return "Try decreasing size to about " + formatter.format(newDelta);
	// } else if (prob > 0.40) {
	// return "Try increasing size to about " + formatter.format(newDelta);
	// } else return "";
	// }

}
