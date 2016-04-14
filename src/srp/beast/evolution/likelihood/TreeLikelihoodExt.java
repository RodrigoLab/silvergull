package srp.beast.evolution.likelihood;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import beast.core.Description;
import beast.core.Input;
import beast.core.Input.Validate;
import beast.core.util.Log;
import beast.evolution.alignment.Alignment;
import beast.evolution.branchratemodel.StrictClockModel;
import beast.evolution.likelihood.BeagleTreeLikelihood;
import beast.evolution.likelihood.BeerLikelihoodCore;
import beast.evolution.likelihood.BeerLikelihoodCore4;
import beast.evolution.likelihood.TreeLikelihood;
import beast.evolution.likelihood.TreeLikelihood.Scaling;
import beast.evolution.sitemodel.SiteModel;
import beast.evolution.substitutionmodel.SubstitutionModel;
import beast.evolution.tree.Node;
import beast.evolution.tree.Tree;
import beast.evolution.tree.TreeInterface;
import srp.dr.ext.SitePatternsExt;
import srp.evolution.OperationRecord;
import srp.evolution.OperationType;
import srp.evolution.haplotypes.HaplotypeModel;


@Description("Extension of beast2 TreeLikelihood class")
public class TreeLikelihoodExt extends TreeLikelihood {
	
    final public Input<HaplotypeModel> haplotypeInput = new Input<>("hapModel", "sequence data for the beast.tree", Validate.REQUIRED);
    
    BeagleTreeLikelihood beagle;

    /**
     * BEASTObject associated with inputs. Since none of the inputs are StateNodes, it
     * is safe to link to them only once, during initAndValidate.
     */
    SubstitutionModel substitutionModel;

    
    /**
     * memory allocation for probability tables obtained from the SiteModel *
     */
    double[] probabilities;

    int matrixSize;

    /**
     * flag to indicate ascertainment correction should be applied *
     */
    boolean useAscertainedSitePatterns = false;

    /**
     * dealing with proportion of site being invariant *
     */
    double proportionInvariant = 0;
    List<Integer> constantPattern = null;

    ///////////////////////////
    
    final private static boolean DEBUG = false;
    private HaplotypeModel haplotypeModel;
    
    
//    final public Input<Alignment> dataInput = nulnew Input<>("data", "sequence data for the beast.tree", Validate.REQUIRED);
    
    @Override
    public void initAndValidate() {
    	System.out.println("initAndValidate TreeLikelihoodExt with haplotypeModel");
        // sanity check: alignment should have same #taxa as tree
        if (haplotypeInput.get().getTaxonCount() != treeInput.get().getLeafNodeCount()) {
        	System.out.println("aoeu2" + haplotypeInput.get().getTaxonCount() +"\t"+ treeInput.get().getLeafNodeCount());
            throw new IllegalArgumentException("The number of nodes in the tree does not match the number of sequences");
        }
        beagle = null;
        
//        beagle = new BeagleTreeLikelihood();
        System.out.println("Only java core now, NO beagle yet");
//        try {
//        	System.out.println("aoeu4");
////	        beagle.initByName(
////                    "data", dataInput.get(), "tree", treeInput.get(), "siteModel", siteModelInput.get(),
////                    "branchRateModel", branchRateModelInput.get(), "useAmbiguities", m_useAmbiguities.get(), 
////                    "useTipLikelihoods", m_useTipLikelihoods.get(),"scaling", scaling.get().toString());
////	        if (beagle.beagle != null) {
////	            //a Beagle instance was found, so we use it
////	            return;
////	        }
//        } catch (Exception e) {
//			// ignore
//		}
        // No Beagle instance was found, so we use the good old java likelihood core
        beagle = null;

        int nodeCount = treeInput.get().getNodeCount();
        if (!(siteModelInput.get() instanceof SiteModel.Base)) {
        	throw new IllegalArgumentException("siteModel input should be of type SiteModel.Base");
        }
        m_siteModel = (SiteModel.Base) siteModelInput.get();
        m_siteModel.setDataType(haplotypeInput.get().getDataType());
        substitutionModel = m_siteModel.substModelInput.get();

        if (branchRateModelInput.get() != null) {
            branchRateModel = branchRateModelInput.get();
        } else {
            branchRateModel = new StrictClockModel();
        }
        m_branchLengths = new double[nodeCount];
        storedBranchLengths = new double[nodeCount];

        int stateCount = haplotypeInput.get().getMaxStateCount();//TODO: maybe fix at 4?
        int patterns = haplotypeInput.get().getPatternCount();//TODO: maybe fix at hap length?
        if (stateCount == 4) {
            likelihoodCore = new BeerLikelihoodCore4();
        } else {
            likelihoodCore = new BeerLikelihoodCore(stateCount);
        }

        String className = getClass().getSimpleName();

//        Alignment alignment = haplotypeInput.get();//FIXME
//
        Log.info.println(className + "(" + getID() + ") uses " + likelihoodCore.getClass().getSimpleName());
//        Log.info.println("  " + alignment.toString(true));
        // print startup messages via Log.print*

        proportionInvariant = m_siteModel.getProportionInvariant();
        m_siteModel.setPropInvariantIsCategory(false);
        if (proportionInvariant > 0) {
            calcConstantPatternIndices(patterns, stateCount);
        }

        initCore();

        patternLogLikelihoods = new double[patterns];
        m_fRootPartials = new double[patterns * stateCount];
        matrixSize = (stateCount + 1) * (stateCount + 1);
        probabilities = new double[(stateCount + 1) * (stateCount + 1)];
        Arrays.fill(probabilities, 1.0);

		if (haplotypeInput.get().isAscertained) {
			useAscertainedSitePatterns = true;
		}
		

    }
    
    /**
     * Determine indices of m_fRootProbabilities that need to be updates
     * // due to sites being invariant. If none of the sites are invariant,
     * // the 'site invariant' category does not contribute anything to the
     * // root probability. If the site IS invariant for a certain character,
     * // taking ambiguities in account, there is a contribution of 1 from
     * // the 'site invariant' category.
     */
    void calcConstantPatternIndices(final int patterns, final int stateCount) {
        constantPattern = new ArrayList<>();
        for (int i = 0; i < patterns; i++) {
            final int[] pattern = haplotypeInput.get().getPattern(i);
            final boolean[] isInvariant = new boolean[stateCount];
            Arrays.fill(isInvariant, true);
            for (final int state : pattern) {
                final boolean[] isStateSet = haplotypeInput.get().getStateSet(state);
                if (m_useAmbiguities.get() || !haplotypeInput.get().getDataType().isAmbiguousState(state)) {
                    for (int k = 0; k < stateCount; k++) {
                        isInvariant[k] &= isStateSet[k];
                    }
                }
            }
            for (int k = 0; k < stateCount; k++) {
                if (isInvariant[k]) {
                    constantPattern.add(i * stateCount + k);
                }
            }
        }
    }
    @Override
    protected void initCore() {
        final int nodeCount = treeInput.get().getNodeCount();
        likelihoodCore.initialize(
                nodeCount,
                haplotypeInput.get().getPatternCount(),
                m_siteModel.getCategoryCount(),
                true, m_useAmbiguities.get()
        );

        final int extNodeCount = nodeCount / 2 + 1;
        final int intNodeCount = nodeCount / 2;

        if (m_useAmbiguities.get() || m_useTipLikelihoods.get()) {
            setPartials(treeInput.get().getRoot(), haplotypeInput.get().getPatternCount());
        } else {
            setStates(treeInput.get().getRoot(), haplotypeInput.get().getPatternCount());
        }
        hasDirt = Tree.IS_FILTHY;
        for (int i = 0; i < intNodeCount; i++) {
            likelihoodCore.createNodePartials(extNodeCount + i);
        }
    }
    

    /**
     * set leaf states in likelihood core *
     */
    @Override
    protected void setStates(Node node, int patternCount) {
        if (node.isLeaf()) {
            HaplotypeModel data = haplotypeInput.get();
            int i;
            int[] states = new int[patternCount];
            int taxonIndex = getTaxonIndex(node.getID(), data);
            for (i = 0; i < patternCount; i++) {
                int code = data.getPattern(taxonIndex, i);
                int[] statesForCode = data.getDataType().getStatesForCode(code);
                if (statesForCode.length==1)
                    states[i] = statesForCode[0];
                else
                    states[i] = code; // Causes ambiguous states to be ignored.
            }
            likelihoodCore.setNodeStates(node.getNr(), states);

        } else {
            setStates(node.getLeft(), patternCount);
            setStates(node.getRight(), patternCount);
        }
    }
    
    /**
    *
    * @param taxon the taxon name as a string
    * @param data the alignment
    * @return the taxon index of the given taxon name for accessing its sequence data in the given alignment,
    *         or -1 if the taxon is not in the alignment.
    */
    private int getTaxonIndex(String taxon, HaplotypeModel data) {
		int taxonIndex = data.getTaxonIndex(taxon);
		if (taxonIndex == -1) {
			if (taxon.startsWith("'") || taxon.startsWith("\"")) {
				taxonIndex = data.getTaxonIndex(taxon.substring(1, taxon.length() - 1));
			}
			if (taxonIndex == -1) {
				throw new RuntimeException(
						"Could not find sequence " + taxon + " in the alignment");
			}
		}
		return taxonIndex;
	}



	/**
     * set leaf partials in likelihood core *
     */
    @Override
    protected void setPartials(Node node, int patternCount) {
        if (node.isLeaf()) {
        	HaplotypeModel data = haplotypeInput.get();
            int states = data.getDataType().getStateCount();
            double[] partials = new double[patternCount * states];
            int k = 0;
            int taxonIndex = getTaxonIndex(node.getID(), data);
            for (int patternIndex_ = 0; patternIndex_ < patternCount; patternIndex_++) {                
                double[] tipLikelihoods = data.getTipLikelihoods(taxonIndex,patternIndex_);
                if (tipLikelihoods != null) {
                	for (int state = 0; state < states; state++) {
                		partials[k++] = tipLikelihoods[state];
                	}
                }
                else {
                	int stateCount = data.getPattern(taxonIndex, patternIndex_);
	                boolean[] stateSet = data.getStateSet(stateCount);
	                for (int state = 0; state < states; state++) {
	                	 partials[k++] = (stateSet[state] ? 1.0 : 0.0);                
	                }
                }
            }
            likelihoodCore.setNodePartials(node.getNr(), partials);

        } else {
            setPartials(node.getLeft(), patternCount);
            setPartials(node.getRight(), patternCount);
        }
    }
    
    /**
     * Calculate the log likelihood of the current state.
     *
     * @return the log likelihood.
     */
    double m_fScale = 1.01;
    int m_nScale = 0;
    int X = 100;


	
    @Override
    public double calculateLogP() {
        if (beagle != null) {
            logP = beagle.calculateLogP();
            return logP;
        }
        final TreeInterface tree = treeInput.get();

        try {
        	if (traverse(tree.getRoot()) != Tree.IS_CLEAN)
        		calcLogP();
        }
        catch (ArithmeticException e) {
        	return Double.NEGATIVE_INFINITY;
        }
        m_nScale++;
        if (logP > 0 || (likelihoodCore.getUseScaling() && m_nScale > X)) {
//            System.err.println("Switch off scaling");
//            m_likelihoodCore.setUseScaling(1.0);
//            m_likelihoodCore.unstore();
//            m_nHasDirt = Tree.IS_FILTHY;
//            X *= 2;
//            traverse(tree.getRoot());
//            calcLogP();
//            return logP;
        } else if (logP == Double.NEGATIVE_INFINITY && m_fScale < 10 && !scaling.get().equals(Scaling.none)) { // && !m_likelihoodCore.getUseScaling()) {
            m_nScale = 0;
            m_fScale *= 1.01;
            Log.warning.println("Turning on scaling to prevent numeric instability " + m_fScale);
            likelihoodCore.setUseScaling(m_fScale);
            likelihoodCore.unstore();
            hasDirt = Tree.IS_FILTHY;
            traverse(tree.getRoot());
            calcLogP();
            return logP;
        }
        return logP;
    }
    
    
	protected void calcLogP() {
        logP = 0.0;
        if (useAscertainedSitePatterns) {
            final double ascertainmentCorrection = haplotypeInput.get().getAscertainmentCorrection(patternLogLikelihoods);
            for (int i = 0; i < haplotypeInput.get().getPatternCount(); i++) {
                logP += (patternLogLikelihoods[i] - ascertainmentCorrection) * haplotypeInput.get().getPatternWeight(i);
            }
        } else {
            for (int i = 0; i < haplotypeInput.get().getPatternCount(); i++) {
                logP += patternLogLikelihoods[i] * haplotypeInput.get().getPatternWeight(i);
            }
        }
    }
    

    /* Assumes there IS a branch rate model as opposed to traverse() */
    
	protected int traverse(final Node node) {

        int update = (node.isDirty() | hasDirt);

        final int nodeIndex = node.getNr();

        final double branchRate = branchRateModel.getRateForBranch(node);
        final double branchTime = node.getLength() * branchRate;

        // First update the transition probability matrix(ices) for this branch
        //if (!node.isRoot() && (update != Tree.IS_CLEAN || branchTime != m_StoredBranchLengths[nodeIndex])) {
        if (!node.isRoot() && (update != Tree.IS_CLEAN || branchTime != m_branchLengths[nodeIndex])) {
            m_branchLengths[nodeIndex] = branchTime;
            final Node parent = node.getParent();
            likelihoodCore.setNodeMatrixForUpdate(nodeIndex);
            for (int i = 0; i < m_siteModel.getCategoryCount(); i++) {
                final double jointBranchRate = m_siteModel.getRateForCategory(i, node) * branchRate;
                System.out.println(node);
                System.out.println(parent.getHeight());
                		System.out.println(node.getHeight()); 
                				System.out.println(jointBranchRate);
                				System.out.println("prob length2");
                				System.out.println(probabilities.length);
                				System.out.println(substitutionModel);
                substitutionModel.getTransitionProbabilities(node, parent.getHeight(), node.getHeight(), jointBranchRate, probabilities);
                //System.out.println(node.getNr() + " " + Arrays.toString(m_fProbabilities));
                likelihoodCore.setNodeMatrix(nodeIndex, i, probabilities);
            }
            update |= Tree.IS_DIRTY;
        }

        // If the node is internal, update the partial likelihoods.
        if (!node.isLeaf()) {

            // Traverse down the two child nodes
            final Node child1 = node.getLeft(); //Two children
            final int update1 = traverse(child1);

            final Node child2 = node.getRight();
            final int update2 = traverse(child2);

            // If either child node was updated then update this node too
            if (update1 != Tree.IS_CLEAN || update2 != Tree.IS_CLEAN) {

                final int childNum1 = child1.getNr();
                final int childNum2 = child2.getNr();

                likelihoodCore.setNodePartialsForUpdate(nodeIndex);
                update |= (update1 | update2);
                if (update >= Tree.IS_FILTHY) {
                    likelihoodCore.setNodeStatesForUpdate(nodeIndex);
                }

                if (m_siteModel.integrateAcrossCategories()) {
                    likelihoodCore.calculatePartials(childNum1, childNum2, nodeIndex);
                } else {
                    throw new RuntimeException("Error TreeLikelihood 201: Site categories not supported");
                    //m_pLikelihoodCore->calculatePartials(childNum1, childNum2, nodeNum, siteCategories);
                }

                if (node.isRoot()) {
                    // No parent this is the root of the beast.tree -
                    // calculate the pattern likelihoods
                    final double[] frequencies = //m_pFreqs.get().
                            substitutionModel.getFrequencies();

                    final double[] proportions = m_siteModel.getCategoryProportions(node);
                    likelihoodCore.integratePartials(node.getNr(), proportions, m_fRootPartials);

                    if (constantPattern != null) { // && !SiteModel.g_bUseOriginal) {
                        proportionInvariant = m_siteModel.getProportionInvariant();
                        // some portion of sites is invariant, so adjust root partials for this
                        for (final int i : constantPattern) {
                            m_fRootPartials[i] += proportionInvariant;
                        }
                    }

                    likelihoodCore.calculateLogLikelihoods(m_fRootPartials, frequencies, patternLogLikelihoods);
                }

            }
        }
        return update;
    } // traverseWithBRM

	
}