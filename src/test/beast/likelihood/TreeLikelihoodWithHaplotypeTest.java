package test.beast.likelihood;


import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import org.junit.Ignore;
import org.junit.Test;

import beast.core.BEASTObject;
import beast.core.Distribution;
import beast.core.Input;
import beast.core.Logger;
import beast.core.MCMC;
import beast.core.Operator;
import beast.core.OperatorSchedule;
import beast.core.State;
import beast.core.StateNode;
import beast.core.StateNodeInitialiser;
import beast.core.Logger.LOGMODE;
import beast.core.Logger.SORTMODE;
import beast.core.parameter.RealParameter;
import beast.evolution.alignment.Alignment;
import beast.evolution.alignment.Sequence;
import beast.evolution.datatype.UserDataType;
import beast.evolution.operators.DeltaExchangeOperator;
import beast.evolution.operators.ScaleOperator;
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

//
//	@Ignore
//    public void testJC69Likelihood2() throws Exception {
//        // Set up JC69 model: uniform freqs, kappa = 1, 0 gamma categories
//        Alignment dataAli = BEASTTestCase.getAlignment();
//        HaplotypeModel data = new HaplotypeModel(dataAli);
//        Tree tree = BEASTTestCase.getTree(dataAli);
//
//        JukesCantor JC = new JukesCantor();
//        JC.initAndValidate();
//
//        SiteModel siteModel = new SiteModel();
//        siteModel.initByName("mutationRate", "1.0", "gammaCategoryCount", 1, "substModel", JC);
//
//        TreeLikelihoodExt likelihood = new TreeLikelihoodExt();
//        likelihood.initByName("hapModel", data, "tree", tree, "siteModel", siteModel);
//        double logP = 0;
//        logP = likelihood.calculateLogP();
//        assertEquals(logP, -1992.2056440317247, BEASTTestCase.PRECISION);
//
//        
//        
//        Operator op = new HaplotypeBaseSingleOperator();
//        op.setID("baseSingle");
//        op.initByName("hapModel", data, "weight" ,1.0);
//        
//
//        OperatorSchedule operatorSchedule = new OperatorSchedule();
////        for (final Operator op : operatorsInput.get()) {
//            operatorSchedule.addOperator(op);
////        }
//
//
//        
//    
//        
//		RealParameter parameter = new RealParameter(new Double[] { 1., 1., 1., 1. });
//		ArrayList<StateNode> stateNodeList = new ArrayList<StateNode>();
//		
//		stateNodeList.add(parameter);
//		stateNodeList.add(data);
//
//		
//		State state = new State();
//		state.initByName("stateNode", stateNodeList );
//		state.initialise();
//		state.setPosterior(likelihood);
//		
//		DeltaExchangeOperator d = new DeltaExchangeOperator();
//		// An invalid operator should either fail in initByName or make valid
//		// proposals
//		try {
//			d.initByName("parameter", parameter);
//		} catch (RuntimeException e) {
//			return;
//		}
//		for (final StateNode stateNode : op.listStateNodes()) {
//            System.out.println(stateNode +"\t"+ stateNode.getID());
//        }	
//		logP = likelihood.calculateLogP();
//		System.out.println("LogP: "+ logP);
//		for (int i = 0; i < 5; i++) {
//			System.out.println("\n\n_Start "+i);
////			logP = likelihood.calculateLogP();
//			System.out.println("LogP: "+ logP);
//			op.proposal();
//			state.storeCalculationNodes();
//			state.checkCalculationNodesDirtiness();
//			
//		
////		d.proposal();    
////		state.storeCalculationNodes();
////		state.checkCalculationNodesDirtiness();
//		logP = likelihood.calculateLogP();
//		System.out.println("LogP: "+ logP);
//		//TODO: FIXME: store/restore haplotypeModel. Previously, it deal with haploytpe.
//		//But we calculate everything with sitePatterns now. Either store/restore in Haplotype,
//		//Or create and use. storeSitePattern
//		
//		  if (Randomizer.nextDouble() < 0) {
//              // accept
//              
//              state.acceptCalculationNodes();
//
//              
//              op.accept();
//              
//              //System.err.print(" accept");
//          } else {
//              // reject
//              
//                  op.reject(0);
//              
//              state.restore();
//              state.restoreCalculationNodes();
//              //System.err.print(" reject");
//          }
//          state.setEverythingDirty(false);
////          logP = likelihood.calculateLogP();
////  		System.out.println("LogP: -1992.205  "+ logP);
//          
//		}
//		for (int i = 0; i < data.getHaplotypeCount(); i++) {
//			System.out.println(data.getHaplotypeString(i));
//			
//		}
////		likelihood.initAndValidate();
//		logP = likelihood.calculateLogP();
//		System.out.println("LogP: "+ logP);
//        assertEquals(logP, -1992.2056440317247, BEASTTestCase.PRECISION);
//
//            
//            
//            
//        
////        likelihood.initByName("useAmbiguities", true, "hapModel", data, "tree", tree, "siteModel", siteModel);
////        logP = likelihood.calculateLogP();
////        assertEquals(logP, -1992.2056440317247, BEASTTestCase.PRECISION);
////        
//        
//        
//    }
//    
//    

    @Test
    public void testLikelihoodMCMC() throws Exception {
        // Set up JC69 model: uniform freqs, kappa = 1, 0 gamma categories
        Alignment dataAli = BEASTTestCase.getAlignment();
        HaplotypeModel data = new HaplotypeModel(dataAli);
        Tree tree = BEASTTestCase.getTree(dataAli);

//        JukesCantor JC = new JukesCantor();
//        JC.initAndValidate();
        Frequencies freqs = new Frequencies();
        freqs.initByName("data", dataAli,
              "estimate", false);

//        RealParameter parameter = new RealParameter(new Double[] { 1., 1., 1., 1. });
//		DeltaExchangeOperator d = new DeltaExchangeOperator();
//			d.initByName("parameter", parameter);
        
//        RealParameter kappa = new RealParameter(new Double[] {27.40259});
        RealParameter kappa = new RealParameter(27.40259+"");
        HKY hky = new HKY();
        hky.initByName("kappa", kappa, "frequencies", freqs);

        
        SiteModel siteModel = new SiteModel();
        siteModel.initByName("mutationRate", "1.0", "gammaCategoryCount", 1, "substModel", hky);

        TreeLikelihoodExt likelihood = new TreeLikelihoodExt();
        likelihood.initByName("hapModel", data, "tree", tree, "siteModel", siteModel);
        double logP = 0;
        logP = likelihood.calculateLogP();
        assertEquals(logP, -1856.303048876734, BEASTTestCase.PRECISION);

        
        ArrayList<Operator> opsList = new ArrayList<>();
        Operator op;
        op = new HaplotypeBaseSingleOperator();
        op.setID("baseSingle");
        op.initByName("hapModel", data, "weight" ,1.0);
        opsList.add(op);
//        
        op = new ScaleOperator();
        op.setID("treeScale");
        op.initByName("tree", tree, "weight", 1.0);
        opsList.add(op);
        
        op = new ScaleOperator();
        op.setID("kappaScalp");
        op.initByName("parameter", kappa, "weight", 1.0);
        opsList.add(op);
        
        OperatorSchedule operatorSchedule = new OperatorSchedule();
//        for (final Operator op : operatorsInput.get()) {
        operatorSchedule.addOperator(op);
//        }
        
        
        likelihood.setID("treeLikelihoodExt");
        kappa.setID("kappa");
        Logger logger = new Logger();
        logger.loggersInput.setValue(likelihood, logger);
//        logger.loggersInput.setValue(kappa, logger);
        logger.setInputValue("log", kappa);
//        logger.setInputValue("model", kappa);
//        logger.setInputValue("model", data);
//        System.out.println(logger.modeInput.getTipText());
//        logger.setInputValue("sort", SORTMODE.smart);
        logger.setInputValue("mode", LOGMODE.compound);
        logger.setInputValue("sanitiseHeaders", true);
        logger.initByName("logEvery", 1000);
        


//        final public Input<String> fileNameInput = new Input<>("fileName", "Name of the file, or stdout if left blank");
//
//        final public Input<Integer> everyInput = new Input<>("logEvery", "Number of the samples logged", 1);
//        final public Input<BEASTObject> modelInput = new Input<>("model", "Model to log at the top of the log. " +
//                "If specified, XML will be produced for the model, commented out by # at the start of a line. " +
//                "Alignments are suppressed. This way, the log file documents itself. ");
//        final public Input<LOGMODE> modeInput = new Input<>("mode", "logging mode, one of " + Arrays.toString(LOGMODE.values()), LOGMODE.autodetect, LOGMODE.values());
//        final public Input<SORTMODE> sortModeInput = new Input<>("sort", "sort items to be logged, one of " + Arrays.toString(SORTMODE.values()), SORTMODE.none, SORTMODE.values());
//        final public Input<Boolean> sanitiseHeadersInput = new Input<>("sanitiseHeaders", "whether to remove any clutter introduced by Beauti" , false);

        
        
        MCMC mcmc = new MCMC();
        mcmc.setInputValue("operator", opsList);
        mcmc.initByName("chainLength", 10000,
        				"storeEvery", 10,
        				"preBurnin", 0,
        				"distribution", likelihood,
        				"logger", logger
//        				"operatorschedule", operatorSchedule
//        				"operator"
        				);
//        final public Input<Integer> chainLengthInput =
//                new Input<>("chainLength", "Length of the MCMC chain i.e. number of samples taken in main loop",
//                        Input.Validate.REQUIRED);
//
//        final public Input<State> startStateInput =
//                new Input<>("state", "elements of the state space");
//
//        final public Input<List<StateNodeInitialiser>> initialisersInput =
//                new Input<>("init", "one or more state node initilisers used for determining " +
//                        "the start state of the chain",
//                        new ArrayList<>());
//
//        final public Input<Integer> storeEveryInput =
//                new Input<>("storeEvery", "store the state to disk every X number of samples so that we can " +
//                        "resume computation later on if the process failed half-way.", -1);
//
//        final public Input<Integer> burnInInput =
//                new Input<>("preBurnin", "Number of burn in samples taken before entering the main loop", 0);
//
//
//        final public Input<Integer> numInitializationAttempts =
//                new Input<>("numInitializationAttempts", "Number of initialization attempts before failing (default=10)", 10);
//
//        final public Input<Distribution> posteriorInput =
//                new Input<>("distribution", "probability distribution to sample over (e.g. a posterior)",
//                        Input.Validate.REQUIRED);
//
//        final public Input<List<Operator>> operatorsInput =
//                new Input<>("operator", "operator for generating proposals in MCMC state space",
//                        new ArrayList<>());//, Input.Validate.REQUIRED);
//
//        final public Input<List<Logger>> loggersInput =
//                new Input<>("logger", "loggers for reporting progress of MCMC chain",
//                        new ArrayList<>(), Input.Validate.REQUIRED);
//
//        final public Input<Boolean> sampleFromPriorInput = new Input<>("sampleFromPrior", "whether to ignore the likelihood when sampling (default false). " +
//                "The distribution with id 'likelihood' in the posterior input will be ignored when this flag is set.", false);
//
//        final public Input<OperatorSchedule> operatorScheduleInput = new Input<>("operatorschedule", "specify operator selection and optimisation schedule", new OperatorSchedule());

        mcmc.run();
        
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
