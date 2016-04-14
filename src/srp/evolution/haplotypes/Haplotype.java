package srp.evolution.haplotypes;

import beast.evolution.datatype.DataType;
import beast.evolution.datatype.Nucleotide;
import beast.util.Randomizer;

import java.util.ArrayList;
import java.util.List;

import beast.core.BEASTObject;
import beast.core.Description;
import beast.core.Input;
import beast.evolution.alignment.Alignment;
import beast.evolution.alignment.Sequence;
//import beast.evolution.util.Taxon;
//import beast.math.MathUtils;
import beast.evolution.alignment.Taxon;

@Description("Single haplotype in the HaplotypeModel. ONLY support Nucleotide for now")
public class Haplotype extends BEASTObject {

	final public Input<String> taxonInput = new Input<>("taxon", "name of the taxa",
			Input.Validate.REQUIRED);
	final public Input<String> dataInput = new Input<>("value",
			"sequence data, either encoded as a string, whitespace is ignored.",
			Input.Validate.REQUIRED);
	final public Input<Boolean> uncertainInput = new Input<>("uncertain",
			"if true, sequence is provided as comma separated probabilities for each character, with sites separated by a semi-colons. In this formulation, gaps are coded as 1/K,...,1/K, where K is the number of states in the model.");

	
	
	char[] haplotype;
	char[] storedHaplotype;
	int haplotypeLength;
//	private String sequenceString;

	public static final DataType DATA_TYPE = new Nucleotide();
	public static final int STATE_COUNT = DATA_TYPE.getStateCount();
	
	//WORK with state in BEAST 2
	
	public Haplotype(int haplotypeLength) {

		this.haplotypeLength = haplotypeLength;
		haplotype = new char[this.haplotypeLength];
		storedHaplotype = new char[this.haplotypeLength];
		
		initHaplotype();
	}

	private void initHaplotype() {

		for (int i = 0; i < haplotypeLength; i++) {
//			int state = MathUtils.nextInt(STATE_COUNT);
			int state = Randomizer.nextInt(STATE_COUNT);
//			char newChar = DATA_TYPE.getChar(state);
			haplotype[i] = (char) state;
		}
		storeState();
		
	}

	public Haplotype(String sequence) {
		this(sequence.length());
		setSequenceString(sequence.toUpperCase());
		storeState();
	}


	public Haplotype(String taxon, int haplotypeLength) {
		this(haplotypeLength);
		taxonInput.setValue(taxon, this);

	}

	public Haplotype(Sequence sequence) {
		this(sequence.getTaxon(), sequence.getData());
	}

	public Haplotype(String taxon, String sequence) {
		taxonInput.setValue(taxon, this);
		dataInput.setValue(sequence, this);
		initAndValidate();

	}
	
//	public Haplotype(Alignment alignment) {
//		
//		taxonInput.setValue(alignment.getTaxaNames(), this);
//		dataInput.setValue(alignment.getsInputs(), this);
//		initAndValidate();
//
//	}
	
	@Override
	public void initAndValidate() {
		Sequence x;

		String data = dataInput.get();
        // remove spaces
        data = data.replaceAll("\\s", "");
        List<Integer> sequence = DATA_TYPE.string2state(data);
		
	}
	
	/**
     * @return the taxon of this sequence as a string.
     */
    public final String getTaxon() {
        return taxonInput.get();
    }

    /**
     * @return the data of this sequence as a string.
     */
    public final String getData() {
        return dataInput.get();
    }



	public void setCharAt(int index, int newChar) {
		setCharAt(index, (char) newChar); 
	}
	
	public void setCharAt(int index, char newChar) {
		haplotype[index]=newChar;
	}

	@Deprecated
	public char replaceCharAt(int index, int newChar){
		char oldChar = getChar(index);
		setCharAt(index, (char) newChar);
		return oldChar;
	}
	
//	@Override
	public char getStoredChar(int index) {
        return storedHaplotype[index];
    }
//	public int getStoredState(int index) {
//		return dataType.getState(storedHaplotype[index]);
////        return storedHaplotype[index];
//    }
	// **************************************
	// OVERRIDE ALL (almost all) methods
	// Do NOT call setState()!!
	// ************************************
	
//	@Override
//    public void setDataType(DataType dataType) {
//        this.dataType = DATA_TYPE;
//    }

    /**
     * @return the length of the sequences.
     */
    
	public int getLength() {
        return haplotypeLength;
    }

    /**
     * @return a String containing the sequences.
     */
    
	public String getSequenceString() {
    	return String.valueOf(haplotype);
    	
    }

    /**
     * @return a char containing the state at index.
     */
    
	public char getChar(int index) {
        return haplotype[index];
    }

    /**
     * @return the state at site index.
     */
//    @Override
//	public int getState(int index) {
//        return DATA_TYPE.getState(haplotype[index]);
//    }

    /**
     */
//	public void setState(int index, int state) {
//
//        sequenceString.setCharAt(index, dataType.getChar(state));
//    }

    /**
     * Characters are copied from the sequences into the destination character array dst.
     */
//    @Override
//	public void getChars(int srcBegin, int srcEnd, char[] dst, int dstBegin) {
//        System.arraycopy(haplotype, srcBegin, dst, dstBegin, srcEnd - srcBegin);
//    }

    /**
     * Set the DataType of the sequences.
     */
//    @Override
//	public DataType guessDataType() {
//        DataType guessDataType = DataType.guessDataType(String.valueOf(haplotype));
//        if(guessDataType.getName().equals(DATA_TYPE.getName())){
//        	return DATA_TYPE;
//        }
//        else{
//        	throw new IllegalArgumentException("Only support "+DATA_TYPE.getName()+". Please check your haplotypes");
//        }
//        
//    }

    /**
     * Set the sequences using a string.
     */
    
	public void setSequenceString(String data) {
		List<Integer> sequence;
    	
	    	
	        // remove spaces
        data = data.replaceAll("\\s", "");
        sequence = DATA_TYPE.string2state(data);
		//setSequenceString(sequence.toCharArray());
    }
    
    public void setSequenceString(char[] sequenceArray) {
    	if(sequenceArray.length != haplotypeLength){
			throw new IllegalArgumentException("Invalid sequence length: "
					+ sequenceArray.length
					+ ". Haplotype length must be equal to " + haplotypeLength);
		}
    	else{
			System.arraycopy(sequenceArray, 0, haplotype, 0, haplotypeLength);
    	}
    	
    }


	public void storeState(int index) {
		storedHaplotype[index] = haplotype[index];
//		storeState();
	}
	public void restoreState(int index) {
		haplotype[index] = storedHaplotype[index];
//		restoreState();
	}
	
	public void storeState() {
		System.arraycopy(haplotype, 0, storedHaplotype, 0, haplotypeLength);
	}

	public void restoreState() {
		char[] temp = storedHaplotype;
		storedHaplotype = haplotype;
		haplotype = temp;
	}

//	public static Haplotype duplicateHaplotype(Haplotype oldHaplotype) {
//		
//		String newTaxon = oldHaplotype.getTaxon();
//		Haplotype newHaplotype = new Haplotype(newTaxon, oldHaplotype.getSequenceString());
//
//		return newHaplotype;
//	}

	public String getStoredSequenceString() {
    	return String.valueOf(storedHaplotype);
	}

	public static Haplotype duplicateHaplotype(Haplotype haplotype2) {
		// TODO Auto-generated method stub
		return null;
	}

	



}