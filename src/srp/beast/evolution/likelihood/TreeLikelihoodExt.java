package srp.beast.evolution.likelihood;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;

import javax.swing.tree.TreeModel;

import beast.core.Description;
import beast.core.Input;
import beast.core.Input.Validate;
import beast.core.util.Log;
import beast.evolution.alignment.TaxonSet;
import beast.evolution.branchratemodel.BranchRateModel;
import beast.evolution.branchratemodel.StrictClockModel;
import beast.evolution.likelihood.BeagleTreeLikelihood;
import beast.evolution.likelihood.BeerLikelihoodCore;
import beast.evolution.likelihood.BeerLikelihoodCore4;
import beast.evolution.likelihood.GenericTreeLikelihood;
import beast.evolution.likelihood.LikelihoodCore;
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
public class TreeLikelihoodExt extends GenericTreeLikelihood {
	
    final public Input<HaplotypeModel> haplotypeInput = new Input<>("hapModel", "sequence data for the beast.tree", Validate.REQUIRED);

    final public Input<Boolean> m_useAmbiguities = new Input<>("useAmbiguities", "flag to indicate that sites containing ambiguous states should be handled instead of ignored (the default)", false);
    final public Input<Boolean> m_useTipLikelihoods = new Input<>("useTipLikelihoods", "flag to indicate that partial likelihoods are provided at the tips", false);
    
    
    public static enum Scaling {none, always, _default};
    final public Input<Scaling> scaling = new Input<>("scaling", "type of scaling to use, one of " + Arrays.toString(Scaling.values()) + ". If not specified, the -beagle_scaling flag is used.", Scaling._default, Scaling.values());
    


    /**
     * calculation engine *
     */
    protected LikelihoodCore likelihoodCore;
    BeagleTreeLikelihood beagle;

    /**
     * BEASTObject associated with inputs. Since none of the inputs are StateNodes, it
     * is safe to link to them only once, during initAndValidate.
     */
    SubstitutionModel substitutionModel;
    protected SiteModel.Base m_siteModel;
    protected BranchRateModel.Base branchRateModel;

    /**
     * flag to indicate the
     * // when CLEAN=0, nothing needs to be recalculated for the node
     * // when DIRTY=1 indicates a node partial needs to be recalculated
     * // when FILTHY=2 indicates the indices for the node need to be recalculated
     * // (often not necessary while node partial recalculation is required)
     */
    protected int hasDirt;

    /**
     * Lengths of the branches in the tree associated with each of the nodes
     * in the tree through their node  numbers. By comparing whether the
     * current branch length differs from stored branch lengths, it is tested
     * whether a node is dirty and needs to be recomputed (there may be other
     * reasons as well...).
     * These lengths take branch rate models in account.
     */
    protected double[] m_branchLengths;
    protected double[] storedBranchLengths;

    /**
     * memory allocation for likelihoods for each of the patterns *
     */
    protected double[] patternLogLikelihoods;
    /**
     * memory allocation for the root partials *
     */
    protected double[] m_fRootPartials;
    /**
     * memory allocation for probability tables obtained from the SiteModel *
     */
    double[] probabilities;

    int matrixSize;

    /**
     * flag to indicate ascertainment correction should be applied *
     */
//    boolean useAscertainedSitePatterns = false;

    /**
     * dealing with proportion of site being invariant *
     */
    double proportionInvariant = 0;
    List<Integer> constantPattern = null;
    protected HaplotypeModel haplotypeModel;
    
    
    final private static boolean DEBUG = false;
    
    
    @Override
    public void initAndValidate() {
    	if (DEBUG) {
		
    		System.out.println("TreeLikelihoodExt initAndValidate() with haplotypeModel");
    	}
    	haplotypeModel = haplotypeInput.get();
        // sanity check: alignment should have same #taxa as tree
        if (haplotypeInput.get().getTaxonCount() != treeInput.get().getLeafNodeCount()) {
        	System.out.println("aoeu2: " + haplotypeInput.get().getTaxonCount() +"\t"+ treeInput.get().getLeafNodeCount());
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
//        beagle = null;

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

        int stateCount = haplotypeInput.get().getMaxStateCount();
        int patterns = haplotypeInput.get().getPatternCount();
        if (stateCount == 4) {
            likelihoodCore = new BeerLikelihoodCore4();
        } else {
        	Log.err.println("StateCount!=4 "+ stateCount + ". This should not happen in HaplotypeModel!");
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

//		if (haplotypeInput.get().isAscertained) {
//			Log.warning.println("WARNING: ascertainment correction is NOT supported!");
////			useAscertainedSitePatterns = true;
//		}
		
		//////////// new method
		int patternCount = patterns;
		states = new int[nodeCount][patternCount];
	    storedStates = new int[nodeCount][patternCount];
//        for (int i = 0; i < nodeCount; i++) {
//            states[i] = new int[patternCount];
//            storedStates[i] = new int[patternCount];
//        }
//		LogManager.getLogManager().reset();
//		final Logger logger = Logger.getLogger("dr.evomodel");
//        logger.setLevel(Level.OFF);
        treeTaxonIndex = new int[haplotypeModel.getTaxonCount()];
        List<String> taxaNames = haplotypeModel.getTaxaNames();
        
        Node[] nodesAsArray = treeInput.get().getNodesAsArray();
    	for (int i = 0; i < treeTaxonIndex.length; i++) {
    		String taxonId = taxaNames.get(i);
//			treeTaxonIndex[i] = treeInput.get().getTaxonIndex(taxonId);
    		treeTaxonIndex[i] = -1;
    		for (int j = 0; j < nodesAsArray.length; j++) {
    			//TODO: Double check the index
                if (taxonId.equals( nodesAsArray[j].getID() )){
                	treeTaxonIndex[i] = j;
                	break;
                }
            }

    	}
    	for (int j = 0; j < nodesAsArray.length; j++) {
    		Node node = nodesAsArray[j];
	    	if (node.isLeaf()) {
	            HaplotypeModel data = haplotypeInput.get();
	            int taxonIndex = getTaxonIndex(node.getID(), data);
//	            System.out.println(taxonIndex);
	    	}
    	}
    	
//    	System.out.println();
//    	System.out.println(Arrays.toString(treeTaxonIndex));
    	
		
    	
		for (int i = 0; i < haplotypeModel.getHaplotypeCount(); i++) {
//    		String taxonId = haplotypeModel.getTaxonId(i);
			int updateExternalNodeIndex = treeTaxonIndex[i];//treeModel.getTaxonIndex(taxonId);
//	        int updateExternalNodeIndex2 = treeTaxonIndex[i];
//			treeTaxonIndex
//			System.out.println(i +"\t"+ updateExternalNodeIndex +"\t"+
//					states.length +"\t"+ states[updateExternalNodeIndex].length);
//    		likelihoodCore.getNodeStates(updateExternalNodeIndex, states[updateExternalNodeIndex]);
//    		likelihoodCore.setNodeStates(updateExternalNodeIndex, states[i]);
    		
    		if (m_useAmbiguities.get() || m_useTipLikelihoods.get()) {
//                setPartials(treeInput.get().getRoot(), haplotypeInput.get().getPatternCount());
//                likelihoodCore.getNodePartials(updateExternalNodeIndex, states[updateExternalNodeIndex]);
            } else {
            	likelihoodCore.getNodeStates(updateExternalNodeIndex, states[updateExternalNodeIndex]);
            }
		}
		haplotypeModel.setOperationType(OperationType.FULL);
//        store();
		

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
        	Log.warning.println("m_useAmbiguities.get() || m_useTipLikelihoods.get()) are not tested/implemented for haplotype model.");
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
    
    protected void setPartials(Node node, int patternCount) {
    	Log.warning.println("setPartials are NOT tested");
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
        	if (traverse(tree.getRoot()) != Tree.IS_CLEAN){
//        		traverse(tree.getRoot());
        		//new
        		
//        		hasDirt = Tree.IS_FILTHY;
//                traverse(tree.getRoot());
        		// end new
        		
        		
        		calcLogP();
        	}
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
//        if (useAscertainedSitePatterns) {
//        	Log.warning.println("WARNING: ascertainment correction is NOT supported!");
//            final double ascertainmentCorrection = haplotypeInput.get().getAscertainmentCorrection(patternLogLikelihoods);
//            for (int i = 0; i < haplotypeInput.get().getPatternCount(); i++) {
//                logP += (patternLogLikelihoods[i] - ascertainmentCorrection) * haplotypeInput.get().getPatternWeight(i);
//            }
//        } else {
            for (int i = 0; i < haplotypeInput.get().getPatternCount(); i++) {
                logP += patternLogLikelihoods[i] * haplotypeInput.get().getPatternWeight(i);
                
            }
//            for (int i = 0; i < haplotypeModel.getHaplotypeCount(); i++) {
//                likelihoodCore.getNodeStates(i, states[i]);
//                System.out.println(Arrays.toString(states[i]));
//			}
//            System.out.println(Arrays.toString(patternLogLikelihoods));
            
//        }
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

	
	

    /** CalculationNode methods **/

    /**
     * check state for changed variables and update temp results if necessary *
     */
    @Override
    protected boolean requiresRecalculation() {
    	if (DEBUG) {
		
    		System.out.println("TreeLikelihoodExt.java requiresRecalculation()");
    	}
        if (beagle != null) {
        	Log.err.println("Beagle calculation is NOT implemented!!");
//            return beagle.requiresRecalculation();
        }
        hasDirt = Tree.IS_CLEAN;
        
//        System.out.println(haplotypeInput.get().somethingIsDirty());
//        System.out.println(haplotypeModel.somethingIsDirty());
//        System.out.println("M:"+m_siteModel.isDirtyCalculation());
        if (haplotypeInput.get().somethingIsDirty()) {
        	updatePatternList();
            hasDirt = Tree.IS_FILTHY;//TODO: check, maybe Tree.IS_DIRTY;
            return true;
        }
        
        if (m_siteModel.isDirtyCalculation()) {
            hasDirt = Tree.IS_DIRTY;
            return true;
        }
        if (branchRateModel != null && branchRateModel.isDirtyCalculation()) {
            //m_nHasDirt = Tree.IS_DIRTY;
            return true;
        }
        return treeInput.get().somethingIsDirty();
    }

    @Override
    public void store() {
        if (beagle != null) {
            beagle.store();
            super.store();
            return;
        }
        if (likelihoodCore != null) {
            likelihoodCore.store();
        }
        
        super.store();
        System.arraycopy(m_branchLengths, 0, storedBranchLengths, 0, m_branchLengths.length);
        
        if (DEBUG) {
			System.out.println("==store() "+this.getClass().getName());
		}
        storeState();
    }

    @Override
    public void restore() {
        if (beagle != null) {
            beagle.restore();
            super.restore();
            return;
        }
        if (likelihoodCore != null) {
            likelihoodCore.restore();
        }
        
        super.restore();
        double[] tmp = m_branchLengths;
        m_branchLengths = storedBranchLengths;
        storedBranchLengths = tmp;
        
        if (DEBUG) {
			System.out.println("==restore() "+this.getClass().getName());
		}
        restoreState();
    }

    /**
     * @return a list of unique ids for the state nodes that form the argument
     */
    @Override
	public List<String> getArguments() {
        return Collections.singletonList(dataInput.get().getID());
    }

    /**
     * @return a list of unique ids for the state nodes that make up the conditions
     */
    @Override
	public List<String> getConditions() {
        return m_siteModel.getConditions();
    }


	// *****************************************
	// Methods to update likeliood
	// ************************************


//	private HaplotypeModel haplotypeModel;
//	private SitePatternsExt sitePatternExt;
//	private AbstractLikelihoodCore likelihoodCoreA;
	private int[] treeTaxonIndex;
//	private NativeNucleotideLikelihoodCore likelihoodCoreExt;
//	private NativeNucleotideLikelihoodCoreExt likelihoodCoreExt2;
	
	protected int[][] states;
    protected int[][] storedStates;

//    
//	public TreeLikelihoodExt(HaplotypeModel haplotypeModel, TreeModel treeModel,
//			SiteModel siteModel, BranchRateModel branchRateModel,
//			TipStatesModel tipStatesModel, boolean useAmbiguities,
//			boolean allowMissingTaxa, boolean storePartials,
//			boolean forceJavaCore, boolean forceRescaling) {
//		
//		
//		super(new SitePatternsExt (haplotypeModel),
//				treeModel, siteModel, branchRateModel,
//				tipStatesModel, useAmbiguities, allowMissingTaxa,
//				storePartials, forceJavaCore, forceRescaling);
//		
//		this.sitePatternExt = (SitePatternsExt) getPatternList(); 
//		this.haplotypeModel = haplotypeModel;
//		addModel(this.haplotypeModel);
////		tempstates = new int[patternCount];
//		likelihoodCoreA = (AbstractLikelihoodCore)likelihoodCore;
////		likelihoodCoreExt = (NativeNucleotideLikelihoodCore) likelihoodCore;
////		likelihoodCoreExt2 = (NativeNucleotideLikelihoodCoreExt) likelihoodCoreA;
//		states = new int[nodeCount][];
//	    storedStates = new int[nodeCount][];
//        for (int i = 0; i < nodeCount; i++) {
//            
//            states[i] = new int[patternCount];
//            storedStates[i] = new int[patternCount];
//        }
////		LogManager.getLogManager().reset();
////		final Logger logger = Logger.getLogger("dr.evomodel");
////        logger.setLevel(Level.OFF);
//        treeTaxonIndex = new int[haplotypeModel.getHaplotypeCount()];
//    	for (int i = 0; i < treeTaxonIndex.length; i++) {
//    		String taxonId = haplotypeModel.getTaxonId(i);
//			treeTaxonIndex[i] = treeModel.getTaxonIndex(taxonId);
//
//    	}
//    	
//		for (int i = 0; i < haplotypeModel.getHaplotypeCount(); i++) {
//    		String taxonId = haplotypeModel.getTaxonId(i);
//			int updateExternalNodeIndex = treeModel.getTaxonIndex(taxonId);
////	        int updateExternalNodeIndex2 = treeTaxonIndex[i];
//
//    		likelihoodCoreA.getNodeStates(updateExternalNodeIndex, states[updateExternalNodeIndex]);
////    		likelihoodCore.setNodeStates(updateExternalNodeIndex, states[i]);
//    		
//		}
//	}
//
////	@Override
////	public void makeDirty() {
////		haplotypeModel.resetOperation();
////		super.makeDirty();
////	
////	}
//
//    

	    public void storeState(){
	
	//    	sitePatternExt.storeState();
	    	
	    	
	        OperationRecord record = haplotypeModel.getOperationRecord();
			int haplotypeIndex = record.getSpectrumIndex();
	        int updateExternalNodeIndex = treeTaxonIndex[haplotypeIndex];
	
	//		updateNode[updateExternalNodeIndex] = true;
			int site;
			
			OperationType operation = record.getOperation();
	        if (DEBUG) {
				System.out.println("====storeState() in TreeLikelihoodExt:\t"+ operation);
			}
			
			switch (operation) {
			case NONE:
				break;
			case SINGLE:
				site = record.getSingleIndex();
				storedStates[updateExternalNodeIndex][site] = states[updateExternalNodeIndex][site];
	//			likelihoodCoreA.getNodeStates(updateExternalNodeIndex, states[updateExternalNodeIndex]);
				break;
			case MULTI:
				int[] sites = record.getAllSiteIndexs();
				for (int s : sites) {
					storedStates[updateExternalNodeIndex][s] = states[updateExternalNodeIndex][s];
				}
	
				break;
			case COLUMN:
				site = record.getSingleIndex();
				for (int h = 0; h < haplotypeModel.getHaplotypeCount(); h++) {
					updateExternalNodeIndex = treeTaxonIndex[h];
					storedStates[updateExternalNodeIndex][site] = states[updateExternalNodeIndex][site];
				}
				
				
				break;
			case RECOMBINATION:
				int[] twoPositions = record.getRecombinationPositionIndex();
				int[] twoHapIndexs = record.getRecombinationSpectrumIndex();
				for (int h : twoHapIndexs) {
					updateExternalNodeIndex = treeTaxonIndex[h];
					for (int s = twoPositions[0]; s < twoPositions[1]; s++) {
						storedStates[updateExternalNodeIndex][s] = states[updateExternalNodeIndex][s];
					}
					
	
				}
				break;
			case FULL:
	//			System.out.println("storeState updatePatternListExt() Full");
	//			updatePatternListExt(sitePatternExt);
	//			for (int i = 0; i < haplotypeModel.getHaplotypeCount(); i++) {
	//	    		String taxonId = haplotypeModel.getTaxonId(i);
	//				updateExternalNodeIndex = treeModel.getTaxonIndex(taxonId);
	////		        int updateExternalNodeIndex2 = treeTaxonIndex[i];
	//
	//	    		likelihoodCoreA.getNodeStates(updateExternalNodeIndex, states[updateExternalNodeIndex]);
	////	    		likelihoodCore.setNodeStates(updateExternalNodeIndex, states[i]);
	//	    		
	//			}
				for (int i = 0; i < states.length; i++) {
					System.arraycopy(states[i], 0, storedStates[i], 0, haplotypeModel.getHaplotypeCount());
					
				}
				break;
			default:
	//			break;
				throw new IllegalArgumentException("Invalid operation type:"+record.getOperation());
			}
	//    	
//	    	super.store();
	    }


	public void restoreState() {
    	
//    	updatePatternListExt();
//		sitePatternExt.updateAlignment(haplotypeModel);
		
        OperationRecord record = haplotypeModel.getOperationRecord();
		int haplotypeIndex = record.getSpectrumIndex();
        int updateExternalNodeIndex = treeTaxonIndex[haplotypeIndex];
        int site;
//		updateNode[updateExternalNodeIndex] = true;
//		int site;
//        for (int i = 0; i < nodeCount; i++) {
//            updateNode[i] = true;
//        }
//		likelihoodCoreA.getNodeStates(updateExternalNodeIndex, tempstates);
        OperationType operation = record.getOperation();
        if (DEBUG ) {
			System.out.println("====restoreState() in TreeLikelihoodExt:\t"+ operation);
		}
		
		switch (operation) {
		case SINGLE:
			site = record.getSingleIndex();
//			System.out.println(patternList.getPatternState(haplotypeIndex, site));
//			System.out
//					.println(haplotypeModel.getState(haplotypeIndex, site));;
//					System.out.println();
//			System.out.println(tempstates.length);
//			tempstates[site] = patternList.getPatternState(haplotypeIndex, site);
			states[updateExternalNodeIndex][site] = storedStates[updateExternalNodeIndex][site];
			likelihoodCore.setNodeStates(updateExternalNodeIndex, states[updateExternalNodeIndex]);
//			setStates(likelihoodCore, patternList, haplotypeIndex, updateExternalNodeIndex);
			break;
		case MULTI:
			int[] sites = record.getAllSiteIndexs();
			for (int s  : sites) {
				states[updateExternalNodeIndex][s] = storedStates[updateExternalNodeIndex][s];
//				System.arraycopy(storedStates[i], 0, states[i], 0, patternCount);
			}
//			likelihoodCore.setNodeStates(updateExternalNodeIndex, tempstates);
			likelihoodCore.setNodeStates(updateExternalNodeIndex, storedStates[updateExternalNodeIndex]);
//			System.arraycopy(storedStates[i], 0, states[i], 0, patternCount);
//			setStates(likelihoodCore, patternList, haplotypeIndex, updateExternalNodeIndex);
			break;
		case COLUMN:
			site = record.getSingleIndex();
			for (int h = 0; h < treeTaxonIndex.length; h++) {
				updateExternalNodeIndex = treeTaxonIndex[h];
				likelihoodCore.setNodeStates(updateExternalNodeIndex, storedStates[updateExternalNodeIndex]);
				states[updateExternalNodeIndex][site] = storedStates[updateExternalNodeIndex][site];
			}
			
			
			break;
		case RECOMBINATION:
			int[] twoPositions = record.getRecombinationPositionIndex();
			int[] twoHapIndexs = record.getRecombinationSpectrumIndex();
			for (int h : twoHapIndexs) {
				updateExternalNodeIndex = treeTaxonIndex[h];
				likelihoodCore.setNodeStates(updateExternalNodeIndex, storedStates[updateExternalNodeIndex]);
//				states[updateExternalNodeIndex][site] = storedStates[updateExternalNodeIndex][site];
//		        updateExternalNodeIndex = treeTaxonIndex[h];
//				updateNode[updateExternalNodeIndex] = true;
//				likelihoodCoreA.getNodeStates(updateExternalNodeIndex, tempstates);
//			
				for (int s = twoPositions[0]; s < twoPositions[1]; s++) {
					states[updateExternalNodeIndex][s] = storedStates[updateExternalNodeIndex][s];
				}
//				likelihoodCore.setNodeStates(updateExternalNodeIndex, tempstates);

			}
			break;
		case FULL:
//			System.out.println("restoreState updatePatternListExt() Full");
//			updatePatternListExt(sitePatternExt);
//			for (int i = 0; i < states.length; i++) {
//				
//				
//			}
			for (int h = 0; h < treeTaxonIndex.length; h++) {
				updateExternalNodeIndex = treeTaxonIndex[h];
				
				
				
					int[] tmp = states[updateExternalNodeIndex];
					states[updateExternalNodeIndex] = storedStates[updateExternalNodeIndex];
					storedStates[updateExternalNodeIndex] = tmp;

				
				
//				System.arraycopy(storedStates[updateExternalNodeIndex], 0, states[updateExternalNodeIndex], 0, haplotypeModel.getHaplotypeLength());
				likelihoodCore.setNodeStates(updateExternalNodeIndex, states[updateExternalNodeIndex]);

			}

			break;
		default:
//			break;
			throw new IllegalArgumentException("Invalid operation type:"+operation);
		}
//		updatePatternListExt();
//    	super.restore();

    }
//    
    public void updatePatternList() {

///////////////////
//    	if (node.isLeaf()) {
//            HaplotypeModel data = haplotypeInput.get();
//            int i;
//            int[] states = new int[patternCount];
//            int taxonIndex = getTaxonIndex(node.getID(), data);
//            for (i = 0; i < patternCount; i++) {
//                int code = data.getPattern(taxonIndex, i);
//                int[] statesForCode = data.getDataType().getStatesForCode(code);
//                if (statesForCode.length==1)
//                    states[i] = statesForCode[0];
//                else
//                    states[i] = code; // Causes ambiguous states to be ignored.
//            }
//            likelihoodCore.setNodeStates(node.getNr(), states);
//
//        } else {
//            setStates(node.getLeft(), patternCount);
//            setStates(node.getRight(), patternCount);
//        }
    	
    	
    	/////////////////
//		sitePatternExt.updateAlignment(haplotypeModel);
    	
        OperationRecord record = haplotypeModel.getOperationRecord();
		int haplotypeIndex = record.getSpectrumIndex();
		int updateExternalNodeIndex = treeTaxonIndex[haplotypeIndex];

//		updateNode[updateExternalNodeIndex] = true;
		int site;

		likelihoodCore.getNodeStates(updateExternalNodeIndex, states[updateExternalNodeIndex]);

		if (DEBUG) {
    		System.out.println("TreeLikelihoodExt.java UpdatePatternList()" +"\t"+ record.getOperation() +"\t"+hasDirt);
    	}
		
		switch (record.getOperation()) {
		case SINGLE:
			site = record.getSingleIndex();
//			System.out.println(patternList.getPatternState(haplotypeIndex, site));
//			System.out
//					.println(haplotypeModel.getState(haplotypeIndex, site));;
//					System.out.println();
//			System.out.println(tempstates.length);
			
			
                int code = haplotypeModel.getPattern(haplotypeIndex, site);
//                int[] statesForCode = data.getDataType().getStatesForCode(code);
//                if (statesForCode.length==1)
//                    states[site] = statesForCode[0];
//                else
//			System.out.println(states[updateExternalNodeIndex][site] + "\t"
//					+ storedStates[updateExternalNodeIndex][site] + "\t" + code);
			states[updateExternalNodeIndex][site] = code; // Causes ambiguous states to be ignored.


			
//			states[updateExternalNodeIndex][site] = patternList.getPatternState(haplotypeIndex, site);
			likelihoodCore.setNodeStates(updateExternalNodeIndex, states[updateExternalNodeIndex]);
			
//			likelihoodCoreExt2.setNodeStatesSite(updateExternalNodeIndex, states[updateExternalNodeIndex], site);
//			System.out.println("SET: "+"\t"+updateExternalNodeIndex +"\t"+site +"\t"+states[updateExternalNodeIndex][site]);
			break;
		case MULTI:
			int[] sites = record.getAllSiteIndexs();
			for (int s : sites) {
//				states[updateExternalNodeIndex][s] = patternList.getPatternState(haplotypeIndex, s);
//				likelihoodCoreExt2.setNodeStatesSite(updateExternalNodeIndex, states[updateExternalNodeIndex], s);
			}
			likelihoodCore.setNodeStates(updateExternalNodeIndex, states[updateExternalNodeIndex]);
			
//			setStates(likelihoodCore, patternList, haplotypeIndex, updateExternalNodeIndex);
			break;
		case COLUMN:
			site = record.getSingleIndex();
			for (int h = 0; h < haplotypeModel.getHaplotypeCount(); h++) {
//				updateExternalNodeIndex = treeTaxonIndex[h];
////				likelihoodCoreA.
//				updateNode[updateExternalNodeIndex] = true;
//				likelihoodCoreA.getNodeStates(updateExternalNodeIndex, states[updateExternalNodeIndex]);
//				states[updateExternalNodeIndex][site] = patternList.getPatternState(h, site);
//				likelihoodCore.setNodeStates(updateExternalNodeIndex, states[updateExternalNodeIndex]);

			}
			
			
			break;
		case RECOMBINATION:
			int[] twoPositions = record.getRecombinationPositionIndex();
			int[] twoHapIndexs = record.getRecombinationSpectrumIndex();
			for (int h : twoHapIndexs) {

		        updateExternalNodeIndex = treeTaxonIndex[h];

//				updateNode[updateExternalNodeIndex] = true;
//				likelihoodCoreA.getNodeStates(updateExternalNodeIndex, states[updateExternalNodeIndex]);
//			
//				for (int s = twoPositions[0]; s < twoPositions[1]; s++) {
//					states[updateExternalNodeIndex][s] = patternList.getPatternState(h, s);
//				}
//				likelihoodCore.setNodeStates(updateExternalNodeIndex, states[updateExternalNodeIndex]);

			}
			break;
		case FULL:
//			System.out.println("updatePatternListExt() Full"); //TODO: redo this
//			updatePatternListExt(sitePatternExt);
			break;
		default:
//			break;
			throw new IllegalArgumentException("Invalid operation type:"+record.getOperation());
		}
//		int[] sites = record.getAllSiteIndexs();
//		((AbstractLikelihoodCore)likelihoodCore).getNodeStates(updateExternalNodeIndex, tempstates);
//		for (int s : sites) {
//			tempstates[s] = patternList.getPatternState(haplotypeIndex, s);
//		}
		

	}
//
	
	
}
