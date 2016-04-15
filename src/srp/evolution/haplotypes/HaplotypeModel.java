package srp.evolution.haplotypes;


import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import org.w3c.dom.Node;

import beast.core.Input;
import beast.core.StateNode;
import beast.core.Input.Validate;
import beast.core.util.Log;
import beast.evolution.alignment.Alignment;
import beast.evolution.alignment.Sequence;
import beast.evolution.datatype.DataType;
import beast.evolution.datatype.Nucleotide;
import beast.evolution.tree.Tree;
import beast.util.AddOnManager;
import beast.util.TreeParser;
import srp.evolution.OperationRecord;
import srp.evolution.OperationType;
import srp.evolution.shortreads.ShortReadMapping;
//import beast.evolution.tree.TreeModel;
//import beast.inference.model.Parameter;
//import beast.util.NumberFormatter;

//public abstract class Base<T> extends StateNode implements Parameter<T> {
//
// What do we need to do in BESAT 2 is
// HapoltypeModel extends StateNode -> mcmc
// haplotypeModel ?extends Alignment -> treeLikelihood, this cant' happen, maybe modify treelikelihood later
public class HaplotypeModel extends AbstractHaplotypeModel  {
	
	public static final DataType DATA_TYPE = new Nucleotide();

//	private static final int NUCLEOTIDE_STATES[] = Nucleotide.NUCLEOTIDE_STATES;
	private static final int STATE_COUNT = DATA_TYPE.getStateCount();
	
	public static final String TAXON_PREFIX = "hap_";
	
	private static final boolean DEBUG = false;

	@Deprecated
	private static final String MODEL_NAME = "HaplotypeModel";

	@Deprecated
	public Input<List<Sequence>> sequenceInput =
            new Input<>("sequence", "sequence and meta data for particular taxon", new ArrayList<>(), Validate.OPTIONAL);

	final public Input<List<Haplotype>> haplotypeInput =
            new Input<>("sequence", "sequence and meta data for particular taxon", new ArrayList<>(), Validate.OPTIONAL);

	
	
	private ShortReadMapping srpMap;
//	private boolean DEBUG = true;

	
	protected OperationRecord operationRecord;

	
	private void initHaplotypes() {
		for (int i = 0; i < haplotypeCount; i++) {
			String taxon = TAXON_PREFIX + i;
			Haplotype haplotype = new Haplotype(taxon, haplotypeLength);
			setHaplotype(i, haplotype);
		}
	}
	public HaplotypeModel(int hapCount, int hapLength) {
		super(MODEL_NAME, hapCount, hapLength);
		initHaplotypes();
//		storeEverything();
	}
	
	
	public HaplotypeModel(Alignment trueAlignment) {
		this(trueAlignment.getTaxonCount(), trueAlignment.getSiteCount());

//		List<String>
		taxaNames = trueAlignment.getTaxaNames();
		for (int i = 0; i < taxaNames.size(); i++) {
			Haplotype haplotype = new Haplotype(trueAlignment.getSequenceAsString(taxaNames.get(i)));
			setHaplotype(i, haplotype);
		}
		
		//HACK:
		sequenceInput = trueAlignment.sequenceInput;
		initAndValidate();
//		storeEverything();
	}
	
	@Override
	public void initAndValidate() {
		System.out.println("initAndValidate HapoltyeModel");
		if (haplotypes.length != haplotypeCount){
			throw new IllegalArgumentException("The number of haplotypes does not match!!");
		}
//

//        if (siteWeightsInput.get() != null) {
//            String str = siteWeightsInput.get().trim();
//            String[] strs = str.split(",");
//            siteWeights = new int[strs.length];
//            for (int i = 0; i < strs.length; i++) {
//                siteWeights[i] = Integer.parseInt(strs[i].trim());
//            }
//        }

        // determine data type, either user defined or one of the standard ones
//        if (userDataTypeInput.get() != null) {
//            m_dataType = userDataTypeInput.get();
//        } else {
//            if (types.indexOf(dataTypeInput.get()) < 0) {
////                throw new IllegalArgumentException("data type + '" + dataTypeInput.get() + "' cannot be found. " +
////                        "Choose one of " + Arrays.toString(types.toArray(new String[0])));
//            }
//            // seems to spend forever in there??
////            List<String> dataTypes = AddOnManager.find(beast.evolution.datatype.DataType.class, IMPLEMENTATION_DIR);
////            for (String dataTypeName : dataTypes) {
////                DataType dataType;
////				try {
////					dataType = (DataType) Class.forName(dataTypeName).newInstance();
////	                if (dataTypeInput.get().equals(dataType.getTypeDescription())) {
////	                    m_dataType = dataType;
////	                    break;
////	                }
////				} catch (InstantiationException | IllegalAccessException | ClassNotFoundException e) {
////					throw new IllegalArgumentException(e.getMessage());
////				}
////            }
//        }
        m_dataType = DATA_TYPE;
        // initialize the sequence list
        if (sequenceInput.get().size() > 0) {
            sequences = sequenceInput.get();
        } else {
            // alignment defined by a map of id -> sequence
        	java.util.Map<String, String> map = null;
            List<String> taxa = new ArrayList<>();
            taxa.addAll(map.keySet());
            sequences.clear();
            for (String key : taxa) {
                String sequence = map.get(key);
                sequences.add(new Sequence(key, sequence));
            }
        }

        // initialize the alignment from the given list of sequences
        initializeWithSequenceList(sequences, true);
//        if (taxonSetInput.get() != null && taxonSetInput.get().getTaxonCount() > 0) {
//            sortByTaxonSet(taxonSetInput.get());
//        }
        Log.info.println(toString(false));
    
		
	}
	
//	public HaplotypeModel(int noOfRecoveredHaplotype, ShortReadMapping srpMap) {
//		this(noOfRecoveredHaplotype, srpMap.getLength());
//		for (int i = 0; i < sequenceCount; i++) {
//			Haplotype haplotype = getHaplotype(i);
//			for (int s = 0; s < sequenceLength; s++) {
//				char newChar = srpMap.getBaseAt(s);
//				haplotype.setCharAt(s, newChar);
//			}
//			
//		}
//	}
	//TODO: What is a good starting point?
	public HaplotypeModel(int noOfRecoveredHaplotype, ShortReadMapping srpMap) {
			this(noOfRecoveredHaplotype, srpMap.getLength());
			this.srpMap = srpMap;
			for (int i = 0; i < haplotypeCount; i++) {
				char[] randHap = this.srpMap.getSemiRandHaplotype2();
//				char[] randHap = this.srpMap.getRandHaplotype();
//				System.out.println(String.valueOf(randHap));
				Haplotype haplotype = getHaplotype(i);
//				for (int s = 0; s < sequenceLength; s++) {
//					haplotype.setCharAt(s, randHap[s]);
//				}
				haplotype.setSequenceString(randHap);
				setHaplotype(i, haplotype);
				
			}
			
	//		storeEverything();
		}

	public  ArrayList<Integer> getGetMapToSrpAt(int pos){
		return srpMap.getMapToSrp(pos);
	}
	
	public ArrayList<int[]> getListOfAvailableChar2(){
		return srpMap.getListOfAvailableChar2();
	}
	public char getNextBaseFromSrpMap(int s){
		char baseAt = srpMap.getBaseAt(s);
		return baseAt;
	}
	public int calculateSPS(){
		return 0;//TODO test this
//		return SPSDist.calculeteSPS(this, this);
		
	}

	public String getHaplotypeString(int i) {
		return getHaplotype(i).getSequenceString();
	}

	
	public int[] getStoredSitePattern(int siteIndex){

		int n = getHaplotypeCount();
    	int[] pattern = new int[n];
    	int[] newPattern = new int[n];
        for (int i = 0; i < n; i++) {
            Haplotype hap = getHaplotype(i);
//            pattern[i] = hap.getStoredState(siteIndex);
//            newPattern[i] = hap.getState(siteIndex);
//            System.out.println(hap.getChar(siteIndex)+"\t"+ hap.getStoredChar(siteIndex));
        }
//        System.out.println(Arrays.toString(pattern) +"\t"+ Arrays.toString(newPattern));
        
        return pattern;

	}
	

	
//	@Override
//	public String toString(){
//
//        NumberFormatter formatter = new NumberFormatter(6);
//
//        StringBuilder buffer = new StringBuilder();
//
////	        boolean countStatistics = !(dataType instanceof Codons) && !(dataType instanceof GeneralDataType);
//
////        if (countStatistics) {
////            buffer.append("Site count = ").append(getSiteCount()).append("\n");
////            buffer.append("Invariant sites = ").append(getInvariantCount()).append("\n");
////            buffer.append("Singleton sites = ").append(getSingletonCount()).append("\n");
////            buffer.append("Parsimony informative sites = ").append(getInformativeCount()).append("\n");
////            buffer.append("Unique site patterns = ").append(getUniquePatternCount()).append("\n\n");
////        }
//        for (int i = 0; i < getHaplotypeCount(); i++) {
//            String name = formatter.formatToFieldWidth(getTaxon(i).getId(), 10);
//            buffer.append(">" + name + "\n");
//            buffer.append(getAlignedSequenceString(i) + "\n");
//        }
//
//        return buffer.toString();
//		
//	}
	
	protected void storeEverything(){
		for (int i = 0; i < getHaplotypeCount(); i++) {
			Haplotype haplotype = getHaplotype(i);
			haplotype.storeState();
		}
		
	}
	
	protected void storeState() {

		OperationType operation = operationRecord.getOperation();
		int haplotypeIndex;
		int siteIndex;
		int[] siteIndexs;
		Haplotype haplotype;
		if (DEBUG) {
			System.out.println("StoreState in HaplotypeModel:\t"+ operation);
		}
		switch (operation) {
		case NONE:
//			System.out.println("Init HaplotypeModel StoreState()");//FIXME: maybe not the best default option
			break;
		case FULL:
			storeEverything();
			break;

		case SINGLE:
			haplotypeIndex = operationRecord.getSpectrumIndex();
			siteIndex = operationRecord.getSingleIndex();
			haplotype = getHaplotype(haplotypeIndex);
			haplotype.storeState(siteIndex);
			break;

		case MULTI:
			haplotypeIndex = operationRecord.getSpectrumIndex();
			siteIndexs = operationRecord.getAllSiteIndexs();
			haplotype = getHaplotype(haplotypeIndex);
			for (int s: siteIndexs){
				haplotype.storeState(s);
			}
			break;

		case COLUMN:
			siteIndex = operationRecord.getSingleIndex();
			for (int i = 0; i < getHaplotypeCount(); i++) {
				haplotype = getHaplotype(i);
				haplotype.storeState(siteIndex);
			}
			break;

		case RECOMBINATION:
			int[] twoSpectrums = operationRecord.getRecombinationSpectrumIndex();
			int[] twoPositions = operationRecord.getRecombinationPositionIndex();
			for (int i : twoSpectrums) {
				haplotype = getHaplotype(i);
				for (int s = twoPositions[0]; s < twoPositions[1]; s++) {
					haplotype.storeState(s);
				}
			}
			break;
		default:
			throw new IllegalArgumentException("Unknown operation type: "
					+ operation + "\tin"
					+ HaplotypeModel.class.getSimpleName());

		}
	}

	
	protected void restoreState() {
	
		OperationType operation = operationRecord.getOperation();
		int hapIndex;
		int siteIndex;
		Haplotype haplotype;
		int[] siteIndexs;
		if (DEBUG) {
			System.out.println("RestoreState in HaplotypeModel:\t" + operation);
		}
		switch (operation) {
		
		case NONE:
			break;
		case FULL:
			for (int i = 0; i < getHaplotypeCount(); i++) {
				haplotype = getHaplotype(i);
				haplotype.restoreState();
			}
			break;

		case SINGLE:
			hapIndex = operationRecord.getSpectrumIndex();
			siteIndex = operationRecord.getSingleIndex();
//			System.out.println(hapIndex +"\t"+ siteIndex);
			haplotype = getHaplotype(hapIndex);
			haplotype.restoreState(siteIndex);
			break;
		
		case MULTI:
			hapIndex = operationRecord.getSpectrumIndex();
			siteIndexs = operationRecord.getAllSiteIndexs();
			haplotype = getHaplotype(hapIndex);
			for (int s: siteIndexs){
				haplotype.restoreState(s);
			}
			break;

		case COLUMN:
	//			spectrumIndex = spectrumOperationRecord.getSpectrumIndex();
			siteIndex = operationRecord.getSingleIndex();
			for (int i = 0; i < getHaplotypeCount(); i++) {
				haplotype = getHaplotype(i);
				haplotype.restoreState(siteIndex);
			}
			break;
			
		case RECOMBINATION:
//			System.err.println("Restore alignment recombination");
			int[] twoSpectrums = operationRecord.getRecombinationSpectrumIndex();
			int[] twoPositions = operationRecord.getRecombinationPositionIndex();
			for (int i : twoSpectrums) {
				haplotype = getHaplotype(i);
				for (int s = twoPositions[0]; s < twoPositions[1]; s++) {
					haplotype.restoreState(s);
				}
			}
			break;
		default:
			throw new IllegalArgumentException("Unknown operation type: " + operation);

		}
	}

	public static HaplotypeModel factory(Alignment shortReads, Alignment trueAlignment){
		
		HaplotypeModel haplotypeModel = new HaplotypeModel(trueAlignment);
		return haplotypeModel;
	}


	public static HaplotypeModel duplicateHaplotypeModel(
			HaplotypeModel oldModel) {
		HaplotypeModel haplotypeModel = new HaplotypeModel(oldModel.getHaplotypeCount(),
				oldModel.getHaplotypeLength());

		for (int i = 0; i < oldModel.getHaplotypeCount(); i++) {
			Haplotype haplotype = Haplotype.duplicateHaplotype(oldModel.getHaplotype(i));
			haplotypeModel.setHaplotype(i, haplotype);
		}
		
		return haplotypeModel;

	}
	
	public static Alignment CreateAlignmentFromHaplotypeModel(HaplotypeModel hm){
		
		ArrayList<Sequence> seqList = new ArrayList<>();
		for (int i = 0; i < hm.getHaplotypeCount(); i++) {
			Sequence s = new Sequence("hap_"+i, hm.getHaplotypeString(i));
			seqList.add(s);
		}
		Alignment ali = new Alignment(seqList, DATA_TYPE.getTypeDescription());
//		ali.
		return ali;
	}
	
//
//	public void simulateSequence(TreeLikelihoodExt treeLikelihood) {
//
//        double substitutionRate = 1000.1;//(getHaplotypeCount()*getHaplotypeLength()) ;
//        double damageRate = 0;
//        SiteModel siteModel = treeLikelihood.getSiteModel();
//        SubstitutionModel substitutionModel = siteModel.getSubstitutionModel();
//        
//        int[] initialSequence = aMap.getConsensusSequenceState();
//        StringBuffer buffer = new StringBuffer();
//        for (int i = 0; i < initialSequence.length; i++) {
//            buffer.append(Nucleotides.INSTANCE.getChar(    initialSequence[i] ));
//        }
//        treeLikelihood.makeDirty();
//        treeLikelihood.getLogLikelihood();
//        System.err.println(buffer.toString());
//
//
//
////    	Arrays.fill(initialSequence, 1);
//        SeqGenExt seqGen = new SeqGenExt(initialSequence, 
//                substitutionRate, substitutionModel, siteModel,
//                damageRate);
//        
//        Tree tree = treeLikelihood.getTreeModel();
//        jebl.evolution.alignments.Alignment jeblAlignment = seqGen.simulate(tree);
//        List<Sequence> sequenceList = jeblAlignment.getSequenceList();
//        for (int j = 0; j < sequenceList.size(); j++) {
//			
//
////        	System.err.println(getHaplotypeString(j));
//			System.out.println(sequenceList.get(j).getString());
//			Haplotype haplotype = getHaplotype(j);
//			haplotype.setSequenceString(sequenceList.get(j).getString());
//			System.out.println(haplotype.getSequenceString());
//			System.out.println();
//		}
////        SimpleAlignment
//		
//		
//	}
//
//	public void simulateSequence(double errorRate, SiteModel siteModel, SubstitutionModel substitutionModel,
//			TreeModel treeModel, ShortReadMapping srpMap) {
//
//        double substitutionRate = errorRate/(getHaplotypeCount()*getHaplotypeLength()) ;
//        System.err.println(substitutionRate);
//        double damageRate = 0;
////        SiteModel siteModel = treeLikelihood.getSiteModel();
////        SubstitutionModel substitutionModel = siteModel.getSubstitutionModel();
//        
//        int[] initialSequence = srpMap.getConsensusSequenceState();
//        StringBuffer buffer = new StringBuffer();
//        for (int i = 0; i < initialSequence.length; i++) {
//            buffer.append(Nucleotides.INSTANCE.getChar(    initialSequence[i] ));
//        }
//
//        SeqGenExt seqGen = new SeqGenExt(initialSequence, 
//                substitutionRate, substitutionModel, siteModel, damageRate);
//        
//        jebl.evolution.alignments.Alignment jeblAlignment = seqGen.simulate(treeModel);
//        List<Sequence> sequenceList = jeblAlignment.getSequenceList();
//        for (int j = 0; j < sequenceList.size(); j++) {
//			Haplotype haplotype = getHaplotype(j);
//			haplotype.setSequenceString(sequenceList.get(j).getString());
//		}
//		
//	}
	
	////////////////////////////////////////////////////
	//// implement functions from StateNode
	///////////////////////////////////////////////////
	
	@Override
	public int getDimension() {
		// TODO Auto-generated method stub
		return 0;
	}

	@Override
	public double getArrayValue() {
		// TODO Auto-generated method stub
		return 0;
	}

	@Override
	public double getArrayValue(int dim) {
		// TODO Auto-generated method stub
		return 0;
	}

	

	@Override
	public void setEverythingDirty(boolean isDirty) {
		// TODO should be a lot more here, maybe some operationRecord things
		setSomethingIsDirty(isDirty);
		
	}

	@Override
	public HaplotypeModel copy() {
		// TODO Auto-generated method stub
		return null;
	}


    /**
     * copy of all values into existing tree *
     */
    @Override
    public void assignTo(final StateNode other) {
		// TODO Auto-generated method stub
//        final Tree tree = (Tree) other;
//        final Node[] nodes = new Node[nodeCount];
//        listNodes(tree.root, nodes);
//        tree.setID(getID());
//        //tree.index = index;
//        root.assignTo(nodes);
//        tree.root = nodes[root.getNr()];
//        tree.nodeCount = nodeCount;
//        tree.internalNodeCount = internalNodeCount;
//        tree.leafNodeCount = leafNodeCount;
    }

    /**
     * copy of all values from existing tree *
     */
    @Override
    public void assignFrom(final StateNode other) {
		// TODO Auto-generated method stub
//        final Tree tree = (Tree) other;
//        final Node[] nodes = new Node[tree.getNodeCount()];//tree.getNodesAsArray();
//        for (int i = 0; i < tree.getNodeCount(); i++) {
//            nodes[i] = newNode();
//        }
//        setID(tree.getID());
//        //index = tree.index;
//        root = nodes[tree.root.getNr()];
//        root.assignFrom(nodes, tree.root);
//        root.parent = null;
//        nodeCount = tree.nodeCount;
//        internalNodeCount = tree.internalNodeCount;
//        leafNodeCount = tree.leafNodeCount;
//        initArrays();
    }

    /**
     * as assignFrom, but only copy tree structure *
     */
    @Override
    public void assignFromFragile(final StateNode other) {
		// TODO Auto-generated method stub
//        final Tree tree = (Tree) other;
//        if (m_nodes == null) {
//            initArrays();
//        }
//        root = m_nodes[tree.root.getNr()];
//        final Node[] otherNodes = tree.m_nodes;
//        final int rootNr = root.getNr();
//        assignFrom(0, rootNr, otherNodes);
//        root.height = otherNodes[rootNr].height;
//        root.parent = null;
//        if (otherNodes[rootNr].getLeft() != null) {
//            root.setLeft(m_nodes[otherNodes[rootNr].getLeft().getNr()]);
//        } else {
//            root.setLeft(null);
//        }
//        if (otherNodes[rootNr].getRight() != null) {
//            root.setRight(m_nodes[otherNodes[rootNr].getRight().getNr()]);
//        } else {
//            root.setRight(null);
//        }
//        assignFrom(rootNr + 1, nodeCount, otherNodes);
    }


	@Override
	public void fromXML(Node node) {
		// TODO Auto-generated method stub

    /**
     * reconstruct tree from XML fragment in the form of a DOM node *
     */

//        final String newick = node.getTextContent();
//        final TreeParser parser = new TreeParser();
//        try {
//            parser.thresholdInput.setValue(1e-10, parser);
//        } catch (Exception e1) {
//            e1.printStackTrace();
//        }
//        try {
//            parser.offsetInput.setValue(0, parser);
//            setRoot(parser.parseNewick(newick));
//        } catch (Exception e) {
//            // TODO Auto-generated catch block
//            e.printStackTrace();
//        }
//        initArrays();
//    	
	}

	@Override
	public int scale(double scale) {
        // nothing to do
        Log.warning.println("Attempt to scale HapoltypeModel " + getID() + "  has no effect");
        return 0;

	}

	@Override
	protected void store() {
		// TODO Auto-generated method stub
		
	}

	@Override
	public void restore() {
		// TODO Auto-generated method stub
		
	}
	
////////////////////////////////////////////////////////////////////////////////////////////
////XXX not used after this
//	public double getLogqFrequency(int oldChar, int newChar){
//		return storedLogqMatrix[NUCLEOTIDE_STATES[oldChar]][NUCLEOTIDE_STATES[newChar]];
//	}
//	
//	public double getLogqFrequencyStates(int oldState, int newState){
//		return storedLogqMatrix[oldState][newState];
//	}
//	
//	private void checkFrequencyParameter(Parameter frequency) {
//
//		for (int i = 0; i < storedFrequency.length; i++) {
//			if(storedFrequency[i]!= frequency.getParameterValue(i)){
//
//				for (int j = i; j < storedFrequency.length; j++) {
//					storedFrequency[j] = frequency.getParameterValue(j);
//					logFreq[j] = Math.log(storedFrequency[j]);
//				}			
//				storedCumSumFrequency[0] = storedFrequency[0];
//				storedCumSumFrequency[1] = storedCumSumFrequency[0]+storedFrequency[1];
//				storedCumSumFrequency[2] = storedCumSumFrequency[1]+storedFrequency[2];
////				storedCumSumFrequency[2] = storedCumSumFrequency[1]+storedFrequency[2];
//
//				for (int j = 0; j < logFreq.length; j++) {
//					for (int k = j+1; k < logFreq.length; k++) {
//						storedLogqMatrix[j][k] = logFreq[j]-logFreq[k];
//						storedLogqMatrix[k][j] = -storedLogqMatrix[j][k];
//					}
////					System.out.println(Arrays.toString(storedLogqMatrix[j]));
//				}
//				
//				break;
//			}
//		}
//	}
//
//	private void checkFrequencyParameterTune(Parameter frequency) {
//
//		for (int i = 0; i < storedFrequency.length; i++) {
//			if(storedFrequency[i]!= frequency.getParameterValue(i)){
//				System.out.println(i);
//long time1 = System.currentTimeMillis();
//for (int t = 0; t < 1e7; t++) {
//	
//				for (int j = i; j < storedCumSumFrequency.length; j++) {
//					storedFrequency[j] = frequency.getParameterValue(j);
//					logFreq[j] = Math.log(storedFrequency[j]);
//				}			
////				storedFrequency = frequency.getParameterValues();
//				storedCumSumFrequency[0] = storedFrequency[0];
//				storedCumSumFrequency[1] = storedCumSumFrequency[0]+storedFrequency[1];
//				storedCumSumFrequency[2] = storedCumSumFrequency[1]+storedFrequency[2];
//				// Too short for a loop?
//				// for (int j = 1; j < INDEX_OF_LAST_VALID_CHARS; j++) {
//				// storedCumSumFrequency[j] = storedCumSumFrequency[j - 1]
//				// + storedFrequency[j];
//				// }
//				
////				double[] logFreq = new double[4];
////				for (int j = 0; j < storedCumSumFrequency.length; j++) {
////					logFreq[j] = Math.log(storedFrequency[j]);
////				}
////				System.err.println(Arrays.toString(storedFrequency));
//				for (int j = 0; j < storedFrequency.length; j++) {
//					for (int k = j+1; k < storedFrequency.length; k++) {
////						storedLogqMatrix[j][k] = Math.log(storedFrequency[j]/storedFrequency[k]);
//						storedLogqMatrix[j][k] = logFreq[j]-logFreq[k];
//						storedLogqMatrix[k][j] = -storedLogqMatrix[j][k];
//					}
////					System.out.println(Arrays.toString(storedLogqMatrix[j]));
//				}
////				System.out.println();
////				for (int j = 0; j < storedFrequency.length; j++) {
////					for (int k = 0; k < storedFrequency.length; k++) {
//////						storedLogqMatrix[j][k] = Math.log(storedFrequency[j]/storedFrequency[k]);
////						storedLogqMatrix[j][k] = logFreq[j]-logFreq[k];
////					}
//////					System.out.println(Arrays.toString(storedLogqMatrix[j]));
////				}
//
////				System.out.println(i +"\t"+ Arrays.toString(storedFrequency));
////				System.err.println(Arrays.toString(storedCumSumFrequency));
//}
//long time2 = System.currentTimeMillis();
//
//System.out.println((time2 - time1) + "\t");
//				
//				break;
//			}
//		}
//	}
//	
//	private double[] logFreq = new double[4];
//	private double[] storedFrequency = new double[4];
//	private double[] storedCumSumFrequency = new double[STATE_COUTN];
//	private double[][] storedLogqMatrix = new double[4][4];
//
//	@Deprecated AlignmentMapping aMap;
//	
//	@Deprecated
//	public void addShortReadMap(ShortReadMapping srpMap) {
//		this.srpMap = srpMap;
//		
//	}
//
//	public ShortReadMapping getShortReadMapping() {
//		
//		return srpMap;
//	}
//
//
	
}