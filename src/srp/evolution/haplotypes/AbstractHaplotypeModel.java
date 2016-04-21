package srp.evolution.haplotypes;


import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;
import java.util.Set;

import beast.core.Input;
import beast.core.Input.Validate;
import beast.core.StateNode;
import beast.core.util.Log;
import beast.evolution.alignment.Alignment;
//import beast.evolution.alignment.PatternList;
import beast.evolution.alignment.Sequence;
import beast.evolution.alignment.TaxonSet;
import beast.evolution.datatype.DataType;
import beast.evolution.datatype.Nucleotide;
import beast.evolution.datatype.StandardData;


public abstract class AbstractHaplotypeModel extends StateNode {
	//Store Haplotype[] and reimplement alignment class 

	public static final DataType DATA_TYPE = new Nucleotide(); 
	
	protected Haplotype[] haplotypes;
	protected int haplotypeCount;
	protected int haplotypeLength;

	
	public AbstractHaplotypeModel(String name, int hapCount, int hapLength) {
//		super(name);

		this.haplotypeCount = hapCount;
		this.haplotypeLength = hapLength;
		haplotypes = new Haplotype[haplotypeCount];
	}


	public Haplotype getHaplotype(int i){
		return haplotypes[i];
	}
	
	public char getHaplotypeCharAt(int hapIndex, int charIndex) {
		return haplotypes[hapIndex].getChar(charIndex);
	}
	
	protected void setHaplotype(int hapIndex, Haplotype haplotype){
		haplotypes[hapIndex] = haplotype;
	}
	
	public int getHaplotypeLength() {
		return haplotypeLength;
	}

	public int getHaplotypeCount() {
		return haplotypeCount;
	}


    /**
     * @return number of states for this siteList
     */
//    @Override
	public int getStateCount() {
        return DATA_TYPE.getStateCount();
    }

    /**
     * Gets the length of the pattern strings which will usually be the
     * same as the number of taxa
     *
     * @return the length of patterns
     */
//    @Override
	public int getPatternLength() {
        return getHaplotypeCount();
    }


    /**
     * @return state at (taxonIndex, patternIndex)
     */
//    @Override
	public int getPatternState(int taxonIndex, int patternIndex) {
    	return 0;//getState(taxonIndex, patternIndex);
    	//TODO: FIX THIS
    }



//    /**
//     * @return the array of pattern weights
//     */
////    @Override
//	public double[] getPatternWeights() {
//        double[] weights = new double[getHaplotypeCount()];
//        for (int i = 0; i < weights.length; i++)
//            weights[i] = 1.0;
//        return weights;
//    }

	
    /**
     * @return the frequency of each state
     */
//    @Override
//	public double[] getStateFrequencies() {
//        return PatternList.Utils.empiricalStateFrequencies(this);
//    }


	
    // **************************************************************
    // Other old IMPLEMENTATION
    // **************************************************************
	
//	/**
//	 * call getHaplotypeString(sequenceIndex);
//	 */
////	@Override
//	public String getAlignedSequenceString(int sequenceIndex) {
//		
//		return haplotypes[sequenceIndex].getSequenceString();
//	}
//
//	/**
//	 * call getHaplotypeString(sequenceIndex);
//	 */
////	@Override
//	public String getUnalignedSequenceString(int sequenceIndex) {
//		
//		return haplotypes[sequenceIndex].getSequenceString();
//	}

	// **************************************************************
	// Loggable beast.core.loggable
	// **************************************************************

    @Override
	public void init(PrintStream out) {
//        Node node = getRoot();
        out.println("#Hapoltype\n");
//        out.println("Begin taxa;");
//        out.println("\tDimensions ntax=" + getLeafNodeCount() + ";");
//        out.println("\t\tTaxlabels");
//        printTaxa(node, out, getNodeCount() / 2);
//        out.println("\t\t\t;");
//        out.println("End;");
//
//        out.println("Begin trees;");
//        out.println("\tTranslate");
//        printTranslate(node, out, getNodeCount() / 2);
//        out.print(";");
    }

    @Override
	public void log(int sample, PrintStream out) {
//        Tree tree = (Tree) getCurrent();
        out.print("Haplotype STATE_" + sample + " = ");

        for (Haplotype h : haplotypes) {
        	out.print(h.getSequenceString());
		}
        
        out.print(";");
    }

    /**
     * @see beast.core.Loggable *
     */
    @Override
	public void close(PrintStream out) {
        out.print("End Haplotype;");
    }

    
    // ***************************************************************
    // Copied from Alignment, don't need most of them
    // ***************************************************************
    

//    static List<String> types = new ArrayList<>();
    /**
     * default data type *
     */
    protected final static String NUCLEOTIDE = "nucleotide";


//    final public Input<List<Sequence>> sequenceInput =
//            new Input<>("sequence", "sequence and meta data for particular taxon", new ArrayList<>(), Validate.OPTIONAL);

//    final public Input<TaxonSet> taxonSetInput =
//            new Input<>("taxa", "An optional taxon-set used only to sort the sequences into the same order as they appear in the taxon-set.", new TaxonSet(), Validate.OPTIONAL);

//    final public Input<Integer> stateCountInput = new Input<>("statecount", "maximum number of states in all sequences");
//    final public Input<String> dataTypeInput = new Input<>("dataType", "data type, one of " + types, NUCLEOTIDE, types.toArray(new String[0]));
//    final public Input<DataType.Base> userDataTypeInput = new Input<>("userDataType", "non-standard, user specified data type, if specified 'dataType' is ignored");
//    final public Input<Boolean> stripInvariantSitesInput = new Input<>("strip", "sets weight to zero for sites that are invariant (e.g. all 1, all A or all unkown)", false);
//    final public Input<String> siteWeightsInput = new Input<>("weights", "comma separated list of weights, one for each site in the sequences. If not specified, each site has weight 1");
//
//    final public Input<Boolean> isAscertainedInput = new Input<>("ascertained", "is true if the alignment allows ascertainment correction, i.e., conditioning the " +
//            "Felsenstein likelihood on excluding constant sites from the alignment", false);
//    /**
//     * Inputs from AscertainedAlignment
//     */
//    final public Input<Integer> excludefromInput = new Input<>("excludefrom", "first site to condition on, default 0", 0);
//    final public Input<Integer> excludetoInput = new Input<>("excludeto", "last site to condition on (but excluding this site), default 0", 0);
//    final public Input<Integer> excludeeveryInput = new Input<>("excludeevery", "interval between sites to condition on (default 1)", 1);

    /**
     * list of sequences in the alignment *
     */
    protected List<Sequence> sequences = new ArrayList<>();

    /**
     * list of taxa names defined through the sequences in the alignment *
     */
    protected List<String> taxaNames = new ArrayList<>();

    /**
     * list of state counts for each of the sequences, typically these are
     * constant throughout the whole alignment.
     */
//    protected List<Integer> stateCounts = new ArrayList<>();

    /**
     * maximum of m_nStateCounts *
     */
    final protected int maxStateCount = DATA_TYPE.getStateCount();;
    
    /**
     * state codes for the sequences *
     */
    protected List<List<Integer>> counts = new ArrayList<>();

    /**
     * data type, useful for converting String sequence to Code sequence, and back *
     */
    final protected DataType m_dataType = DATA_TYPE;

    /**
     * weight over the columns of a matrix *
     */
    protected int[] patternWeight;

    /**
     * weights of sites -- assumed 1 for each site if not specified
     */
//    protected int[] siteWeights = null;

    /**
     * Probabilities associated with each tip of the tree, for use when the
     * characters are uncertain.
     */
    public List<double[][]> tipLikelihoods = new ArrayList<>(); // #taxa x #sites x #states
    private boolean usingTipLikelihoods = false;
    
    /**
     * pattern state encodings *
     */
    protected int [][] sitePatterns; // #patterns x #taxa
    protected int[][] storedSitePatterns;
    /**
     * maps site nr to pattern nr *
     */
    @Deprecated
    protected int[] patternIndex;

    /**
     * From AscertainedAlignment
     */
//    Set<Integer> excludedPatterns;

    /**
     * A flag to indicate if the alignment is ascertained
     */
//    public boolean isAscertained;


    /**
     * Constructor for testing purposes.
     *
     * @param sequences
     * @param dataType
     */
    public static void Alignment(List<Sequence> sequences, String dataType) {
//            for (Sequence sequence : sequences) {
//                sequenceInput.setValue(sequence, this);
//            }
//            dataTypeInput.setValue(dataType, this);
//            initAndValidate();
    }

    @Override
    abstract public void initAndValidate();
    /**
     * Initializes the alignment given the provided list of sequences and no other information.
     * It site weights and/or data type have been previously set up with initAndValidate then they
     * remain in place. This method is used mainly to re-order the sequences to a new taxon order
     * when an analysis of multiple alignments on the same taxa are undertaken.
     *
     * @param sequences
     */
    protected void initializeWithSequenceList(List<Sequence> sequences, boolean log) {
        this.sequences = sequences;
        taxaNames.clear();
//        stateCounts.clear();
        counts.clear();
        System.out.println("initializeWithSequenceList");
        try {
            for (Sequence seq : sequences) {

                counts.add(seq.getSequence(m_dataType));
                if (taxaNames.contains(seq.getTaxon())) {
                    throw new RuntimeException("Duplicate taxon found in alignment: " + seq.getTaxon());
                }
                taxaNames.add(seq.getTaxon());
                tipLikelihoods.add(seq.getLikelihoods());
                // if seq.isUncertain() == false then the above line adds 'null'
	            // to the list, indicating that this particular sequence has no tip likelihood information
                usingTipLikelihoods |= (seq.getLikelihoods() != null);	            

//                if (seq.totalCountInput.get() != null) {
//                    stateCounts.add(seq.totalCountInput.get());
//                } else {
//                    stateCounts.add(m_dataType.getStateCount());
//                }
            }
            if (counts.size() == 0) {
                // no sequence data
                throw new RuntimeException("Sequence data expected, but none found");
            }
        } catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }
        sanityCheckCalcPatternsSetUpAscertainment(log);
    }

    /**
     * Checks that sequences are all the same length, calculates patterns and sets up ascertainment.
     */
    private void sanityCheckCalcPatternsSetUpAscertainment(boolean log) {
        // Sanity check: make sure sequences are of same length
        int length = counts.get(0).size();
        
        if (!(m_dataType instanceof StandardData)) {
            for (List<Integer> seq : counts) {
                if (seq.size() != length) {
                    throw new RuntimeException("Two sequences with different length found: " + length + " != " + seq.size());
                }
                if (seq.size() != haplotypeLength) {
                    throw new RuntimeException("Two haplotypes with different length found: " + haplotypeLength + " != " + seq.size());
                }
            }
        }
        if( counts.size() != haplotypeCount ){
            throw new RuntimeException("Different number of hapoltype count: " + haplotypeCount + " != " + counts.size());
        }
//        if (siteWeights != null && siteWeights.length != length) {
//            throw new RuntimeException("Number of weights (" + siteWeights.length + ") does not match sequence length (" + length + ")");
//        }

        calcPatterns(log);
//        setupAscertainment();
    }

    /**
     * Sorts an alignment by a provided TaxonSet, so that the sequence/taxon pairs in the alignment match the order
     * that the taxa appear in the TaxonSet (i.e. not necessarily alphabetically).
     *
     * @param toSortBy the taxon set that species the order on the taxa.
     */
//    public void sortByTaxonSet(TaxonSet toSortBy) {
//
//        List<Sequence> sortedSeqs = new ArrayList<>();
//        sortedSeqs.addAll(sequences);
//        Collections.sort(sortedSeqs, (Sequence o1, Sequence o2) -> {
//                return Integer.compare(toSortBy.getTaxonIndex(o1.getTaxon()), toSortBy.getTaxonIndex(o2.getTaxon()));
//            }
//        );
//        initializeWithSequenceList(sortedSeqs, false);
//    }

//    void setupAscertainment() {
//        isAscertained = isAscertainedInput.get();
//
//        if (isAscertained) {
//            //From AscertainedAlignment
//        	Log.warning.println("WARNING: ascertainment correction is NOT supported!");
////            int from = excludefromInput.get();
////            int to = excludetoInput.get();
////            int every = excludeeveryInput.get();
////            excludedPatterns = new HashSet<>();
////            for (int i = from; i < to; i += every) {
////                int patternIndex_ = patternIndex[i];
////                // reduce weight, so it does not confuse the tree likelihood
////                patternWeight[patternIndex_] = 0;
////                excludedPatterns.add(patternIndex_);
////            }           
//
//
//        } else {
//        	// sanity check
//            int from = excludefromInput.get();
//            int to = excludetoInput.get();
//            if (from != excludefromInput.defaultValue || to != excludetoInput.defaultValue) {
//            	Log.warning.println("WARNING: excludefrom or excludeto is specified, but 'ascertained' flag is not set to true");
//            	Log.warning.println("WARNING: to suppress this warning, remove the excludefrom or excludeto attributes (if no astertainment correction is required)");
//            	Log.warning.println("WARNING: or set the 'ascertained' flag to true on element with id=" + getID());
//            }
//        }
//
//    } // initAndValidate

    static String getSequence(Alignment data, int taxonIndex) {

        int[] states = new int[data.getPatternCount()];
        for (int i = 0; i < data.getPatternCount(); i++) {
            int[] sitePattern = data.getPattern(i);
            states[i] = sitePattern[taxonIndex];
        }
        try {
            return data.getDataType().state2string(states);
        } catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }
        return null;
    }


    /*
     * assorted getters and setters *
     */
    public List<String> getTaxaNames() {
        if (taxaNames.size() == 0) {
            try {
                initAndValidate();
            } catch (Exception e) {
                e.printStackTrace();
                throw new RuntimeException(e);
            }
        }
        return taxaNames;
    }

//    public List<Integer> getStateCounts() {
//        return stateCounts;
//    }

    /**
     * Returns a List of Integer Lists where each Integer List represents
     * the sequence corresponding to a taxon.  The taxon is identified by
     * the position of the Integer List in the outer List, which corresponds
     * to the nodeNr of the corresponding leaf node and the position of the
     * taxon name in the taxaNames list.
     *
     * @return integer representation of sequence alignment
     */
    public List<List<Integer>> getCounts() {
        return counts;
    }

    public DataType getDataType() {
        return m_dataType;
    }

    /**
     * @return number of taxa in Alignment.
     */
    public int getTaxonCount() {
        //if (taxonsetInput.get() != null) {
        //	return taxonsetInput.get().getTaxonCount();
        //}
//        return taxaNames.size();
    	//FIXME:!! change back to   return taxaNames.size();
    	return haplotypeCount;
    }

//    /**
//     * @return number of taxa in Alignment.
//     * @deprecated Use getTaxonCount() instead.
//     */
//        @Deprecated
//        public int getNrTaxa() {
//            return getTaxonCount();
//        }

    public int getTaxonIndex(String id) {
        return taxaNames.indexOf(id);
    }

    /**
     * @return Number of unique character patterns in alignment.
     */
    public int getPatternCount() {
        return sitePatterns.length;
    }

    public int[] getPattern(int patternIndex_) {
        return sitePatterns[patternIndex_];
    }

    public int getPattern(int taxonIndex, int patternIndex_) {
        return sitePatterns[patternIndex_][taxonIndex];
    }

    /**
     * Retrieve the "weight" of a particular pattern: the number of sites
     * having that pattern.
     *
     * @param patternIndex_ Index into pattern array.
     * @return pattern weight
     */
    public int getPatternWeight(int patternIndex_) {
        return patternWeight[patternIndex_];
    }

    public int getMaxStateCount() {
        return maxStateCount;
    }

    /**
     * Retrieve index of pattern corresponding to a particular site.
     *
     * @param site Index of site.
     * @return Index of pattern.
     */
    @Deprecated
    public int getPatternIndex(int site) {
        return patternIndex[site];
    }

    /**
     * @return Total number of sites in alignment.
     */
    public int getSiteCount() {
        return haplotypeLength;
    }

    /**
     * Retrieve an array containing the number of times each character pattern
     * occurs in the alignment.
     *
     * @return Pattern weight array.
     */
    public int[] getWeights() {
        return patternWeight;
    }


    /**
     * SiteComparator is used for ordering the sites,
     * which makes it easy to identify patterns.
     */
    class SiteComparator implements Comparator<int[]> {
        @Override
		public int compare(int[] o1, int[] o2) {
            for (int i = 0; i < o1.length; i++) {
                if (o1[i] > o2[i]) {
                    return 1;
                }
                if (o1[i] < o2[i]) {
                    return -1;
                }
            }
            return 0;
        }
    } // class SiteComparator


    protected void calcPatterns() {
        calcPatterns(true);
    }

    /**
     * calculate patterns from sequence data
     * *
     */
    private void calcPatterns(boolean log) {
        int taxonCount = counts.size();
        int siteCount = counts.get(0).size();

        // convert data to transposed int array
        int[][] data = new int[siteCount][taxonCount];
        for (int i = 0; i < taxonCount; i++) {
            List<Integer> sites = counts.get(i);
            for (int j = 0; j < siteCount; j++) {
                data[j][i] = sites.get(j);
            }
        }

//        // sort data
//        for (int i = 0; i < data.length; i++) {
//        	System.out.println(Arrays.toString(data[i]));
//		}
//        
        
//        Arrays.sort(data, comparator);
//System.out.println();
//        for (int i = 0; i < data.length; i++) {
//        	System.out.println(Arrays.toString(data[i]));
//		}
//        System.exit(2);
        
        // count patterns in sorted data
        // if (siteWeights != null) the weights are recalculated below
        int patterns = siteCount;
//        int[] weights = new int[siteCount];
//        weights[0] = 1;
//        for (int i = 1; i < siteCount; i++) {
//            if (true || usingTipLikelihoods || comparator.compare(data[i - 1], data[i]) != 0) {
//            	// In the case where we're using tip probabilities, we need to treat each 
//            	// site as a unique pattern, because it could have a unique probability vector.
//                patterns++;
//                data[patterns - 1] = data[i];
//            }
////            weights[patterns - 1]++;
//        }
        
        // reserve memory for patterns
        patternWeight = new int[patterns];
        Arrays.fill(patternWeight, 1);
        
        sitePatterns = new int[patterns][taxonCount];
        storedSitePatterns = new int [patterns][taxonCount];
        for (int i = 0; i < patterns; i++) {
//            patternWeight[i] = weights[i];
            sitePatterns[i] = data[i];
        }

        // find patterns for the sites
        SiteComparator comparator = new SiteComparator();
        patternIndex = new int[siteCount];
        for (int i = 0; i < siteCount; i++) {
            int[] sites = new int[taxonCount];
            for (int j = 0; j < taxonCount; j++) {
                sites[j] = counts.get(j).get(i);
            }
            patternIndex[i] = Arrays.binarySearch(sitePatterns, sites, comparator);
            
        }
//System.out.println(Arrays.toString(patternIndex));
        
//        if (siteWeights != null) {
//            Arrays.fill(patternWeight, 0);
//            for (int i = 0; i < siteCount; i++) {
//                patternWeight[patternIndex[i]] += siteWeights[i];
//            }
//        }
        
        // determine maximum state count
        // Usually, the state count is equal for all sites,
        // though for SnAP analysis, this is typically not the case.
//        maxStateCount = 4;
//        for (int m_nStateCount1 : stateCounts) {
//            maxStateCount = Math.max(maxStateCount, m_nStateCount1);
//        }
        // report some statistics
        if (log && taxaNames.size() < 30) {
            for (int i = 0; i < taxaNames.size(); i++) {
                Log.info.println(taxaNames.get(i) + ": " + counts.get(i).size() + " " );
            }
        }

//        if (stripInvariantSitesInput.get()) {
//            // don't add patterns that are invariant, e.g. all gaps
//            if (log) Log.info.println("Stripping invariant sites");
//
//            int removedSites = 0;
//            for (int i = 0; i < patterns; i++) {
//                int[] pattern = sitePatterns[i];
//                int value = pattern[0];
//                boolean isInvariant = true;
//                for (int k = 1; k < pattern.length; k++) {
//                    if (pattern[k] != value) {
//                        isInvariant = false;
//                        break;
//                    }
//                }
//                if (isInvariant) {
//                    removedSites += patternWeight[i];
//                    patternWeight[i] = 0;
//
//                    if (log) Log.info.print(" <" + value + "> ");
//                }
//            }
//            if (log) Log.info.println(" removed " + removedSites + " sites ");
//        }

        
    } // calcPatterns

    /**
     * @return the total weight of all the patterns (this is the effective number of sites)
     */
    private long getTotalWeight() {
        long totalWeight = 0;
        for (int weight : patternWeight) {
            totalWeight += weight;
        }
        return totalWeight;
    }

    /**
     * Pretty printing of vital statistics of an alignment including id, #taxa, #sites, #patterns and totalweight
     *
     * @param singleLine true if the string should fit on one line
     * @return string representing this alignment
     */
    public String toString(boolean singleLine) {
        long totalWeight = getTotalWeight();
        StringBuilder builder = new StringBuilder();
        builder.append(getClass().getSimpleName() + "(" + getID() + ")");

        if (singleLine) {
            builder.append(": [taxa, patterns, sites] = [" + getTaxonCount() + ", " + getPatternCount());
            builder.append(", " + getTotalWeight() + "]");
        } else {

            long siteCount = getSiteCount();

            builder.append('\n');
            builder.append("  " + getTaxonCount() + " taxa");
            builder.append('\n');
            builder.append("  " + siteCount + (siteCount == 1 ? " site" : " sites") + (totalWeight == getSiteCount() ? "" : " with weight " + totalWeight + ""));
            builder.append('\n');
            if (siteCount > 1) {
                builder.append("  " + getPatternCount() + " patterns");
                builder.append('\n');
            }
        }
        return builder.toString();
    }

    public double[] getTipLikelihoods(int taxonIndex, int patternIndex_) {
    	if (taxonIndex >= tipLikelihoods.size() || tipLikelihoods.get(taxonIndex) == null) { 
    		return null; 
    	} else { 
    		return tipLikelihoods.get(taxonIndex)[patternIndex_];
    	}
    	
    }
    /**
     * returns an array containing the non-ambiguous states that this state represents.
     */
    public boolean[] getStateSet(int state) {
        return m_dataType.getStateSet(state);
//            if (!isAmbiguousState(state)) {
//                boolean[] stateSet = new boolean[m_nMaxStateCount];
//                stateSet[state] = true;
//                return stateSet;
//            } else {
//            }
    }

    boolean isAmbiguousState(int state) {
        return (state >= 0 && state < maxStateCount);
    }

//    //Methods from AscertainedAlignment
//    public Set<Integer> getExcludedPatternIndices() {
//        return excludedPatterns;
//    }
//
//    public int getExcludedPatternCount() {
//        return excludedPatterns.size();
//    }

//    public double getAscertainmentCorrection(double[] patternLogProbs) {
//        double excludeProb = 0, includeProb = 0, returnProb = 1.0;
//
//        for (int i : excludedPatterns) {
//            excludeProb += Math.exp(patternLogProbs[i]);
//        }
//
//        if (includeProb == 0.0) {
//            returnProb -= excludeProb;
//        } else if (excludeProb == 0.0) {
//            returnProb = includeProb;
//        } else {
//            returnProb = includeProb - excludeProb;
//        }
//        return Math.log(returnProb);
//    } // getAscertainmentCorrection

//        /**
//         * Should not be used. No special order of taxa are assumed. Taxa order should be left to user input.
//         */
//        @Deprecated
//        static public void sortByTaxonName(List<Sequence> seqs) {
//            Collections.sort(seqs, (Sequence o1, Sequence o2) -> {
//                    return o1.taxonInput.get().compareTo(o2.taxonInput.get());
//                }
//            );
//        }

    /** 
     * Get String representation of a sequence according to the current datatype
     * @param taxon the name of the taxon to get the sequence from in the alignment
     * @return sequence in String representation
     */
	public String getSequenceAsString(String taxon) {
		int i = getTaxonIndex(taxon);		

		// build up string from underlying data using the current datatype
		int [] states = new int[getSiteCount()];
		for (int k = 0; k < getSiteCount(); k++) {
//			int d = sitePatterns[patternIndex[k]][i];
			int d = sitePatterns[k][i];
			states[k] = d;
		}
		String seq = null;
		try {
			seq = m_dataType.state2string(states);
		} catch (Exception e) {
			e.printStackTrace();
		}
		return seq;
	}


// 
	
	///**
	//* directory to pick up data types from *
	//*/
	//final static String[] IMPLEMENTATION_DIR = {"beast.evolution.datatype"};
	
	/**
	* list of data type descriptions, obtained from DataType classes *
	*/
	//static List<String> types = new ArrayList<>();
	
	//static {
	//  findDataTypes();
	//}
	
	//static public void findDataTypes() {
	//  // build up list of data types
	//  List<String> m_sDataTypes = AddOnManager.find(beast.evolution.datatype.DataType.class, IMPLEMENTATION_DIR);
	//  for (String dataTypeName : m_sDataTypes) {
	//      try {
	//          DataType dataType = (DataType) Class.forName(dataTypeName).newInstance();
	//          if (dataType.isStandard()) {
	//              String description = dataType.getTypeDescription();
	//              if (!types.contains(description)) {
	//                  types.add(description);
	//              }
	//          }
	//      } catch (Exception e) {
	//          // TODO: handle exception
	//      }
	//  }
	//}
    
}