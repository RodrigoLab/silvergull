package srp.evolution.haplotypes;


import java.io.PrintStream;
import java.util.ArrayList;
import java.util.List;

import org.w3c.dom.Node;

import beast.core.StateNode;
import beast.evolution.alignment.Alignment;
import beast.evolution.datatype.DataType;
import beast.evolution.datatype.Nucleotide;
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
	private static final int STATE_COUTN = DATA_TYPE.getStateCount();

	private static final String MODEL_NAME = "HaplotypeModel";
	
	private static final boolean DEBUG = false;

	private ShortReadMapping srpMap;
//	private boolean DEBUG = true;

	public static final String TAXON_PREFIX = "hap_";
	
	protected OperationRecord operationRecord;
	
	private void initHaplotypes() {
		for (int i = 0; i < sequenceCount; i++) {
			String taxon = TAXON_PREFIX + i;
			Haplotype haplotype = new Haplotype(taxon, sequenceLength);
			setHaplotype(i, haplotype);
		}

	}
	
//	public HaplotypeModel(String[] sequences){
//		this(AlignmentUtils.createAlignment(sequences));
//	}

	
	public HaplotypeModel(int hapCount, int hapLength) {
		super(MODEL_NAME, hapCount, hapLength);
		initHaplotypes();
//		storeEverything();
	}
	
	public HaplotypeModel(Alignment trueAlignment) {
		this(trueAlignment.getTaxonCount(), trueAlignment.getSiteCount());

		List<String> taxaNames = trueAlignment.getTaxaNames();
		for (int i = 0; i < taxaNames.size(); i++) {
			Haplotype haplotype = new Haplotype(trueAlignment.getSequenceAsString(taxaNames.get(i)));
			setHaplotype(i, haplotype);
		}
//		storeEverything();
	}
	
	@Override
	public void initAndValidate() {
		// TODO Auto-generated method stub
		
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
			for (int i = 0; i < sequenceCount; i++) {
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

	@Override
	public void assignTo(StateNode other) {
		// TODO Auto-generated method stub
		
	}

	@Override
	public void assignFrom(StateNode other) {
		// TODO Auto-generated method stub
		
	}

	@Override
	public void assignFromFragile(StateNode other) {
		// TODO Auto-generated method stub
		
	}

	@Override
	public void fromXML(Node node) {
		// TODO Auto-generated method stub
		
	}

	@Override
	public int scale(double scale) {
		// TODO Auto-generated method stub
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
	
//	
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