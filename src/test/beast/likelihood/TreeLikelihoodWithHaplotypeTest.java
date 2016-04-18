package test.beast.likelihood;


import java.util.ArrayList;

import org.junit.Test;

import beast.core.Operator;
import beast.core.OperatorSchedule;
import beast.core.State;
import beast.core.StateNode;
import beast.core.parameter.RealParameter;
import beast.evolution.alignment.Alignment;
import beast.evolution.alignment.Sequence;
import beast.evolution.datatype.UserDataType;
import beast.evolution.operators.DeltaExchangeOperator;
import beast.evolution.sitemodel.SiteModel;
import beast.evolution.substitutionmodel.BinaryCovarion;
import beast.evolution.substitutionmodel.Blosum62;
import beast.evolution.substitutionmodel.CPREV;
import beast.evolution.substitutionmodel.Dayhoff;
import beast.evolution.substitutionmodel.Frequencies;
import beast.evolution.substitutionmodel.GTR;
import beast.evolution.substitutionmodel.GeneralSubstitutionModel;
import beast.evolution.substitutionmodel.HKY;
import beast.evolution.substitutionmodel.JTT;
import beast.evolution.substitutionmodel.JukesCantor;
import beast.evolution.substitutionmodel.MTREV;
import beast.evolution.substitutionmodel.MutationDeathModel;
import beast.evolution.substitutionmodel.SubstitutionModel;
import beast.evolution.substitutionmodel.WAG;
import beast.evolution.tree.Tree;
import beast.util.Randomizer;
import junit.framework.TestCase;
import srp.beast.evolution.likelihood.TreeLikelihoodExt;
import srp.evolution.haplotypes.Haplotype;
import srp.evolution.haplotypes.HaplotypeModel;
import srp.evolution.operators.HaplotypeBaseSingleOperator;
import test.beast.BEASTTestCase;
import test.beast.evolution.alignment.UncertainAlignmentTest;


public class TreeLikelihoodWithHaplotypeTest extends TestCase {
//
//    public TreeLikelihoodTest() {
//        super();
//    }
//
//    protected TreeLikelihood newTreeLikelihood() {
//    	System.setProperty("java.only","true");
//        return new TreeLikelihood();
//    }

////    @Test
//    public void testJC69Likelihood() throws Exception {
//        // Set up JC69 model: uniform freqs, kappa = 1, 0 gamma categories
//        Alignment data = BEASTTestCase.getAlignment();
//        Tree tree = BEASTTestCase.getTree(data);
//
//        JukesCantor JC = new JukesCantor();
//        JC.initAndValidate();
//
//        SiteModel siteModel = new SiteModel();
//        siteModel.initByName("mutationRate", "1.0", "gammaCategoryCount", 1, "substModel", JC);
//
//        TreeLikelihood likelihood = new TreeLikelihood();
//        likelihood.initByName("data", data, "tree", tree, "siteModel", siteModel);
//        double logP = 0;
//        logP = likelihood.calculateLogP();
//        assertEquals(logP, -1992.2056440317247, BEASTTestCase.PRECISION);
//
//        likelihood.initByName("useAmbiguities", true, "data", data, "tree", tree, "siteModel", siteModel);
//        logP = likelihood.calculateLogP();
//        assertEquals(logP, -1992.2056440317247, BEASTTestCase.PRECISION);
//    }
    
    @Test
    public void testJC69Likelihood2() throws Exception {
        // Set up JC69 model: uniform freqs, kappa = 1, 0 gamma categories
        Alignment dataAli = BEASTTestCase.getAlignment();
        HaplotypeModel data = new HaplotypeModel(dataAli);
        Tree tree = BEASTTestCase.getTree(dataAli);

        JukesCantor JC = new JukesCantor();
        JC.initAndValidate();

        SiteModel siteModel = new SiteModel();
        siteModel.initByName("mutationRate", "1.0", "gammaCategoryCount", 1, "substModel", JC);

        TreeLikelihoodExt likelihood = new TreeLikelihoodExt();
        likelihood.initByName("hapModel", data, "tree", tree, "siteModel", siteModel);
        double logP = 0;
        logP = likelihood.calculateLogP();
        assertEquals(logP, -1992.2056440317247, BEASTTestCase.PRECISION);

        
        
        Operator op = new HaplotypeBaseSingleOperator();
        op.setID("baseSingle");
        op.initByName("hapModel", data, "weight" ,1.0);
        

        OperatorSchedule operatorSchedule = new OperatorSchedule();
//        for (final Operator op : operatorsInput.get()) {
            operatorSchedule.addOperator(op);
//        }


        
    
        
		RealParameter parameter = new RealParameter(new Double[] { 1., 1., 1., 1. });
		ArrayList<StateNode> stateNodeList = new ArrayList<StateNode>();
		
		stateNodeList.add(parameter);
		stateNodeList.add(data);

		
		State state = new State();
		state.initByName("stateNode", stateNodeList );
		state.initialise();
		state.setPosterior(likelihood);
		
		DeltaExchangeOperator d = new DeltaExchangeOperator();
		// An invalid operator should either fail in initByName or make valid
		// proposals
		try {
			d.initByName("parameter", parameter);
		} catch (RuntimeException e) {
			return;
		}
		for (final StateNode stateNode : op.listStateNodes()) {
            System.out.println(stateNode +"\t"+ stateNode.getID());
        }	
		logP = likelihood.calculateLogP();
		System.out.println("LogP: "+ logP);
		for (int i = 0; i < 5; i++) {
			System.out.println("\n\n_Start "+i);
//			logP = likelihood.calculateLogP();
			System.out.println("LogP: "+ logP);
			op.proposal();
			state.storeCalculationNodes();
			state.checkCalculationNodesDirtiness();
			
		
//		d.proposal();    
//		state.storeCalculationNodes();
//		state.checkCalculationNodesDirtiness();
		logP = likelihood.calculateLogP();
		System.out.println("LogP: "+ logP);
		//TODO: FIXME: store/restore haplotypeModel. Previously, it deal with haploytpe.
		//But we calculate everything with sitePatterns now. Either store/restore in Haplotype,
		//Or create and use. storeSitePattern
		
		  if (Randomizer.nextDouble() < 0) {
              // accept
              
              state.acceptCalculationNodes();

              
              op.accept();
              
              //System.err.print(" accept");
          } else {
              // reject
              
                  op.reject(0);
              
              state.restore();
              state.restoreCalculationNodes();
              //System.err.print(" reject");
          }
          state.setEverythingDirty(false);
//          logP = likelihood.calculateLogP();
//  		System.out.println("LogP: -1992.205  "+ logP);
          
		}
		for (int i = 0; i < data.getHaplotypeCount(); i++) {
			System.out.println(data.getHaplotypeString(i));
			
		}
//		likelihood.initAndValidate();
		logP = likelihood.calculateLogP();
		System.out.println("LogP: "+ logP);
        assertEquals(logP, -1992.2056440317247, BEASTTestCase.PRECISION);

            
            
            
        
//        likelihood.initByName("useAmbiguities", true, "hapModel", data, "tree", tree, "siteModel", siteModel);
//        logP = likelihood.calculateLogP();
//        assertEquals(logP, -1992.2056440317247, BEASTTestCase.PRECISION);
//        
        
        
    }
//    @Test
//    public void testK80Likelihood() throws Exception {
//        // Set up K80 model: uniform freqs, kappa = 27.402591, 0 gamma categories
//    	Alignment dataAli = BEASTTestCase.getAlignment();
//        HaplotypeModel data = new HaplotypeModel(dataAli);
//        Tree tree = BEASTTestCase.getTree(dataAli);
//
//        Frequencies freqs = new Frequencies();
//        freqs.initByName("data", dataAli,
//                "estimate", false);
//
//        HKY hky = new HKY();
//        hky.initByName("kappa", "27.40259", "frequencies", freqs);
//
//        SiteModel siteModel = new SiteModel();
//        siteModel.initByName("mutationRate", "1.0", "gammaCategoryCount", 1, "substModel", hky);
//
//        TreeLikelihoodExt likelihood = new TreeLikelihoodExt();
//        likelihood.initByName("hapModel", data, "tree", tree, "siteModel", siteModel);
//
//        double logP = 0;
//        logP = likelihood.calculateLogP();
//        assertEquals(logP, -1856.303048876734, BEASTTestCase.PRECISION);
//
//        likelihood.initByName("useAmbiguities", true, "hapModel", data, "tree", tree, "siteModel", siteModel);
//        logP = likelihood.calculateLogP();
//        assertEquals(logP, -1856.303048876734, BEASTTestCase.PRECISION);
//    }
//
//    @Test
//    public void testHKY85Likelihood() throws Exception {
//        // Set up HKY85 model: estimated freqs, kappa = 29.739445, 0 gamma categories
//    	Alignment dataAli = BEASTTestCase.getAlignment();
//        HaplotypeModel data = new HaplotypeModel(dataAli);
//        Tree tree = BEASTTestCase.getTree(dataAli);
//
//
//        Frequencies freqs = new Frequencies();
//        freqs.initByName("data", dataAli);
//
//        HKY hky = new HKY();
//        hky.initByName("kappa", "29.739445", "frequencies", freqs);
//
//        SiteModel siteModel = new SiteModel();
//        siteModel.initByName("mutationRate", "1.0", "gammaCategoryCount", 1, "substModel", hky);
//
//        TreeLikelihoodExt likelihood = new TreeLikelihoodExt();
//        likelihood.initByName("hapModel", data, "tree", tree, "siteModel", siteModel);
//
//        double logP = 0;
//        logP = likelihood.calculateLogP();
//        assertEquals(logP, -1825.2131708068507, BEASTTestCase.PRECISION);
//
//        likelihood.initByName("useAmbiguities", true, "hapModel", data, "tree", tree, "siteModel", siteModel);
//        logP = likelihood.calculateLogP();
//        assertEquals(logP, -1825.2131708068507, BEASTTestCase.PRECISION);
//    }
//
//
//    @Test
//    public void testHKY85GLikelihood() throws Exception {
//        // Set up HKY85+G model: estimated freqs, kappa = 38.82974, 4 gamma categories, shape = 0.137064
//    	Alignment dataAli = BEASTTestCase.getAlignment();
//        HaplotypeModel data = new HaplotypeModel(dataAli);
//        Tree tree = BEASTTestCase.getTree(dataAli);
//
//
//        Frequencies freqs = new Frequencies();
//        freqs.initByName("data", dataAli);
//
//        HKY hky = new HKY();
//        hky.initByName("kappa", "38.82974", "frequencies", freqs);
//
//        SiteModel siteModel = new SiteModel();
//        siteModel.initByName("mutationRate", "1.0", "gammaCategoryCount", 4,
//                "shape", "0.137064",
//                "proportionInvariant", "0.0",
//                "substModel", hky);
//
//        TreeLikelihoodExt likelihood = new TreeLikelihoodExt();
//        likelihood.initByName("hapModel", data, "tree", tree, "siteModel", siteModel);
//
//        double logP = 0;
//        logP = likelihood.calculateLogP();
//        System.err.println(logP - -1789.7593576610134);
//        assertEquals(logP, -1789.7593576610134, BEASTTestCase.PRECISION);
//
//        likelihood.initByName("useAmbiguities", true, "hapModel", data, "tree", tree, "siteModel", siteModel);
//        logP = likelihood.calculateLogP();
//        assertEquals(logP, -1789.7593576610134, BEASTTestCase.PRECISION);
//    }
//
//    @Test
//    public void testHKY85ILikelihood() throws Exception {
//        // Set up HKY85+I model: estimated freqs, kappa = 38.564672, 0 gamma categories, prop invariant = 0.701211
//    	Alignment dataAli = BEASTTestCase.getAlignment();
//        HaplotypeModel data = new HaplotypeModel(dataAli);
//        Tree tree = BEASTTestCase.getTree(dataAli);
//
//
//        Frequencies freqs = new Frequencies();
//        freqs.initByName("data", dataAli);
//
//        HKY hky = new HKY();
//        hky.initByName("kappa", "38.564672", "frequencies", freqs);
//
//        SiteModel siteModel = new SiteModel();
//        siteModel.initByName("mutationRate", "1.0", "gammaCategoryCount", 1,
//                "shape", "0.137064",
//                "proportionInvariant", "0.701211",
//                "substModel", hky);
//
//        TreeLikelihoodExt likelihood = new TreeLikelihoodExt();
//        likelihood.initByName("hapModel", data, "tree", tree, "siteModel", siteModel);
//
//        double logP = 0;
//        logP = likelihood.calculateLogP();
//        assertEquals(logP, -1789.912401996943, BEASTTestCase.PRECISION);
//
//        likelihood.initByName("useAmbiguities", true, "hapModel", data, "tree", tree, "siteModel", siteModel);
//        logP = likelihood.calculateLogP();
//        assertEquals(logP, -1789.912401996943, BEASTTestCase.PRECISION);
//    }
//
//    @Test
//    public void testHKY85GILikelihood() throws Exception {
//        // Set up HKY85+G+I model: estimated freqs, kappa = 39.464538, 4 gamma categories, shape = 0.587649, prop invariant = 0.486548
//    	Alignment dataAli = BEASTTestCase.getAlignment();
//        HaplotypeModel data = new HaplotypeModel(dataAli);
//        Tree tree = BEASTTestCase.getTree(dataAli);
//
//
//        Frequencies freqs = new Frequencies();
//        freqs.initByName("data", dataAli);
//
//        HKY hky = new HKY();
//        hky.initByName("kappa", "39.464538", "frequencies", freqs);
//
//        SiteModel siteModel = new SiteModel();
//        siteModel.initByName("mutationRate", "1.0", "gammaCategoryCount", 4,
//                "shape", "0.587649",
//                "proportionInvariant", "0.486548",
//                "substModel", hky);
//
//        TreeLikelihoodExt likelihood = new TreeLikelihoodExt();
//        likelihood.initByName("hapModel", data, "tree", tree, "siteModel", siteModel);
//
//        double logP = 0;
//        logP = likelihood.calculateLogP();
//        assertEquals(logP, -1789.639227747059, BEASTTestCase.PRECISION);
//
//        likelihood.initByName("useAmbiguities", true, "hapModel", data, "tree", tree, "siteModel", siteModel);
//        logP = likelihood.calculateLogP();
//        assertEquals(logP, -1789.639227747059, BEASTTestCase.PRECISION);
//    }
//
//
//    @Test
//    public void testGTRLikelihood() throws Exception {
//        // Set up GTR model: no gamma categories, no proportion invariant
//    	Alignment dataAli = BEASTTestCase.getAlignment();
//        HaplotypeModel data = new HaplotypeModel(dataAli);
//        Tree tree = BEASTTestCase.getTree(dataAli);
//
//
//        Frequencies freqs = new Frequencies();
//        freqs.initByName("data", dataAli);
//
//        GTR gsm = new GTR();
//        gsm.initByName("frequencies", freqs);
//
//        SiteModel siteModel = new SiteModel();
//        siteModel.initByName("mutationRate", "1.0", "gammaCategoryCount", 1,
//                "substModel", gsm);
//
//        TreeLikelihoodExt likelihood = new TreeLikelihoodExt();
//        likelihood.initByName("hapModel", data, "tree", tree, "siteModel", siteModel);
//
//        double logP = 0;
//        logP = likelihood.calculateLogP();
//        assertEquals(logP, -1969.145839307625, BEASTTestCase.PRECISION);
//
//        likelihood.initByName("useAmbiguities", false, "hapModel", data, "tree", tree, "siteModel", siteModel);
//        logP = likelihood.calculateLogP();
//        assertEquals(logP, -1969.145839307625, BEASTTestCase.PRECISION);
//    }
//
//    @Test
//    public void testGTRILikelihood() throws Exception {
//        // Set up GTR model: prop invariant = 0.5
//    	Alignment dataAli = BEASTTestCase.getAlignment();
//        HaplotypeModel data = new HaplotypeModel(dataAli);
//        Tree tree = BEASTTestCase.getTree(dataAli);
//
//
//        Frequencies freqs = new Frequencies();
//        freqs.initByName("data", dataAli);
//
//        GeneralSubstitutionModel gsm = new GeneralSubstitutionModel();
//        gsm.initByName("rates", "1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0", "frequencies", freqs);
//
//        SiteModel siteModel = new SiteModel();
//        siteModel.initByName("mutationRate", "1.0", "gammaCategoryCount", 1,
//                "proportionInvariant", "0.5",
//                "substModel", gsm);
//        //siteModel.init("1.0", 1, null, "0.5", gsm);
//
//        TreeLikelihoodExt likelihood = new TreeLikelihoodExt();
//        likelihood.initByName("hapModel", data, "tree", tree, "siteModel", siteModel);
//
//        double logP = 0;
//        logP = likelihood.calculateLogP();
//        assertEquals(logP, -1948.8417455357564, BEASTTestCase.PRECISION);
//
//        likelihood.initByName("useAmbiguities", false, "hapModel", data, "tree", tree, "siteModel", siteModel);
//        logP = likelihood.calculateLogP();
//        assertEquals(logP, -1948.8417455357564, BEASTTestCase.PRECISION);
//    }
//
//    @Test
//    public void testGTRGLikelihood() throws Exception {
//        // Set up GTR model: 4 gamma categories, gamma shape = 0.5
//    	Alignment dataAli = BEASTTestCase.getAlignment();
//        HaplotypeModel data = new HaplotypeModel(dataAli);
//        Tree tree = BEASTTestCase.getTree(dataAli);
//
//
//        Frequencies freqs = new Frequencies();
//        freqs.initByName("data", dataAli);
//
//        GeneralSubstitutionModel gsm = new GeneralSubstitutionModel();
//        gsm.initByName("rates", "1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0", "frequencies", freqs);
//
//        SiteModel siteModel = new SiteModel();
//        siteModel.initByName("mutationRate", "1.0", "gammaCategoryCount", 4,
//                "shape", "0.5",
//                "substModel", gsm);
//
//        TreeLikelihoodExt likelihood = new TreeLikelihoodExt();
//        likelihood.initByName("hapModel", data, "tree", tree, "siteModel", siteModel);
//
//        double logP = 0;
//        logP = likelihood.calculateLogP();
//        assertEquals(logP, -1949.0360143622, BEASTTestCase.PRECISION);
//
//        likelihood.initByName("useAmbiguities", false, "hapModel", data, "tree", tree, "siteModel", siteModel);
//        logP = likelihood.calculateLogP();
//        assertEquals(logP, -1949.0360143622, BEASTTestCase.PRECISION);
//    }
//
//    @Test
//    public void testGTRGILikelihood() throws Exception {
//        // Set up GTR model: 4 gamma categories, gamma shape = 0.5, prop invariant = 0.5
//    	Alignment dataAli = BEASTTestCase.getAlignment();
//        HaplotypeModel data = new HaplotypeModel(dataAli);
//        Tree tree = BEASTTestCase.getTree(dataAli);
//
//
//        Frequencies freqs = new Frequencies();
//        freqs.initByName("data", dataAli);
//
//        GeneralSubstitutionModel gsm = new GeneralSubstitutionModel();
//        gsm.initByName("rates", "1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0", "frequencies", freqs);
//
//        SiteModel siteModel = new SiteModel();
//        siteModel.initByName("mutationRate", "1.0", "gammaCategoryCount", 4,
//                "shape", "0.5",
//                "proportionInvariant", "0.5",
//                "substModel", gsm);
//
//        TreeLikelihoodExt likelihood = new TreeLikelihoodExt();
//        likelihood.initByName("hapModel", data, "tree", tree, "siteModel", siteModel);
//
//        double logP = 0;
//        logP = likelihood.calculateLogP();
//        assertEquals(logP, -1947.5829396144961, BEASTTestCase.PRECISION);
//
//        likelihood.initByName("useAmbiguities", false, "hapModel", data, "tree", tree, "siteModel", siteModel);
//        logP = likelihood.calculateLogP();
//        assertEquals(logP, -1947.5829396144961, BEASTTestCase.PRECISION);
//    }
    
} // class TreeLikelihoodTest
