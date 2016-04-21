package test.core;

import java.util.ArrayList;

import org.apache.commons.math3.distribution.LogNormalDistribution;
import org.junit.Test;

import beast.core.Distribution;
import beast.core.Function;
import beast.core.Input;
import beast.core.Logger;
import beast.core.MCMC;
import beast.core.Operator;
import beast.core.OperatorSchedule;
import beast.core.Input.Validate;
import beast.core.Logger.LOGMODE;
import beast.core.parameter.RealParameter;
import beast.core.util.CompoundDistribution;
import beast.evolution.alignment.Alignment;
import beast.evolution.operators.DeltaExchangeOperator;
import beast.evolution.operators.Exchange;
import beast.evolution.operators.ScaleOperator;
import beast.evolution.operators.SubtreeSlide;
import beast.evolution.operators.Uniform;
import beast.evolution.operators.WilsonBalding;
import beast.evolution.sitemodel.SiteModel;
import beast.evolution.substitutionmodel.Frequencies;
import beast.evolution.substitutionmodel.HKY;
import beast.evolution.tree.RandomTree;
import beast.evolution.tree.coalescent.Coalescent;
import beast.evolution.tree.coalescent.ConstantPopulation;
import beast.evolution.tree.coalescent.TreeIntervals;
import beast.math.distributions.LogNormalDistributionModel;
import beast.math.distributions.OneOnX;
import beast.math.distributions.ParametricDistribution;
import beast.math.distributions.Prior;
import junit.framework.TestCase;
import srp.beast.evolution.likelihood.TreeLikelihoodExt;
import srp.core.DataImporter;
import srp.evolution.haplotypes.HaplotypeModel;
import srp.evolution.likelihood.haplotypes.ShortReadsHaplotypeLikelihood;
import srp.evolution.operators.HaplotypeBaseSingleOperator;
import srp.evolution.shortreads.ShortReadMapping;
import test.beast.BEASTTestCase;
public class MCMCHaplotypeTest  extends TestCase {

	@Test
	public void testMMCMC() throws Exception {
/*
 * 
 * isntead nxn dist matirx, build 2nx2n, try dist tree and MDS
 */
//		String dataDir = "/home/sw167/workspaceSrp/snowgoose/srp/unittest/testData/";
//		int runIndex = 51;
//		int totalSamples = 1000;
//		int logInterval = 1000;
//		int noOfTrueHaplotype = 7;
//		int noOfRecoveredHaplotype=7;
//		dataDir += "H7_"+runIndex + File.separator;
//		
//		String dataDir = args[0];
//		int runIndex = Integer.parseInt(args[1]);
//		int totalSamples = Integer.parseInt(args[2]);
//		int logInterval = Integer.parseInt(args[3]);
//		int noOfTrueHaplotype = Integer.parseInt(args[4]);
//		int noOfRecoveredHaplotype= Integer.parseInt(args[5]);

		
		String dataDir;
		int runIndex;
		int totalSamples;
		int logInterval;
		int noOfTrueHaplotype;
		int noOfRecoveredHaplotype;
		boolean randomTree = true;
		boolean randomHaplotype = true;
		String inputReadSuffix = "";
		
		boolean isLocal = false;
//		commandLine = false;
		
//		if(args.length == 7){
//			dataDir = args[0];
//			runIndex = Integer.parseInt(args[1]);
//			totalSamples = Integer.parseInt(args[2]);
//			logInterval = Integer.parseInt(args[3]);
//			noOfTrueHaplotype = Integer.parseInt(args[4]);
//			noOfRecoveredHaplotype= Integer.parseInt(args[5]);
//			inputReadSuffix= (args[6]);
//			if(!inputReadSuffix.equals("ART") && !inputReadSuffix.equals("ART_errFree") && !inputReadSuffix.equals("Shrimp") ){
//				System.out.println("Invalid input: InputSuffix should be one of [ART|ART_errFree|Shrimp]");
//				System.exit(-3);
//			}
//		}
//		
//		else{	
			System.out.println("local parameters");
			isLocal = true;
			dataDir = "/home/steven/workspaceDbHap/silvergull/unittest/testData/";
			runIndex = 0;
//			dataDir += "H10_"+runIndex+"/";
			dataDir += "H5_"+runIndex+"/";
//			dataDir += "H5_001/";
			//TODO: local control
			totalSamples = 100	;
			logInterval  = 1000 ;
			
			randomTree = true;
			randomHaplotype = true;
			
//			randomTree = false;
//			randomHaplotype = false;
//			inputReadSuffix = "ART_errFree";
			inputReadSuffix = "ART";
			noOfTrueHaplotype = 5;
			noOfRecoveredHaplotype=5;
//		}
		
		String hapRunIndex = "H"+noOfTrueHaplotype+"_"+runIndex;
//		String shortReadFile = hapRunIndex +"_Srp.fasta";
		String shortReadFile = hapRunIndex +"_ShortRead_"+inputReadSuffix+".fasta";
		String trueHaplotypeFile = hapRunIndex +"_FullHaplotype.fasta";
//		shortReadFile = trueHaplotypeFile;//TODO Remove later. Full test on this later
		
		String prefix = dataDir+"Result_"+hapRunIndex;
		String logTracerName = prefix+".log";
		String logTreeName = prefix+".trees";
		String logHaplotypeName = prefix+".haplotype";
		String operatorAnalysisFile = prefix+"_operatorAnalysisFile.txt";
		
		System.out.println("Input reads file: "+shortReadFile);
		DataImporter dataImporter = new DataImporter(dataDir);

		Alignment shortReads = dataImporter.importAlignment(shortReadFile);
		ShortReadMapping srpMap = new ShortReadMapping(shortReads);

//		srpMap.summary();
//		for (int i = 0; i < 10; i++) {
//			System.out.println((char)srpMap.nextBaseAt(100));
//		}
//		System.exit(10);
		
		HaplotypeModel haplotypeModel = null;
		
//		if(randomHaplotype){
//			haplotypeModel = new HaplotypeModel(noOfRecoveredHaplotype, shortReads.getSiteCount());
			haplotypeModel = new HaplotypeModel(noOfRecoveredHaplotype, srpMap);
//			haplotypeModel = new HaplotypeModel(noOfRecoveredHaplotype, shortReads.getSiteCount());
//			haplotypeModel.addShortReadMap(srpMap);
//		}
//		else{
//////			String partialHaplotypeName = prefix+".haplotypepartial";
////			String partialHaplotypeName = hapRunIndex+"_FullHaplotype.fasta";
//////			Alignment trueAlignment = dataImporter.importAlignment(trueHaplotypeFile);
////			Alignment trueAlignment = dataImporter.importAlignment(partialHaplotypeName);
////			haplotypeModel = new HaplotypeModel(trueAlignment);
////			haplotypeModel.addShortReadMap(srpMap);
//////			haplotypeModel = dataImporter.importPartialSpectrumFile(partialHaplotypeName);
//		}

		// ShortReadLikelihood
		ShortReadsHaplotypeLikelihood srpLikelihood = new ShortReadsHaplotypeLikelihood(haplotypeModel, srpMap);
		System.out.println("Error rate: "+srpLikelihood.ERROR_RATE);
		System.out.println("SRP Likelihood: " + srpLikelihood.calculateSrpLikelihoodFullMaster());
		// coalescent
//		Parameter popSize = new Parameter.Default(ConstantPopulationModelParser.POPULATION_SIZE, 3000.0, 100, 100000.0);
//
//		// Random treeModel
//		ConstantPopulationModel popModel = new ConstantPopulationModel(popSize, Units.Type.YEARS);
//		TreeModel treeModel = MCMCSetupHelperHaplotype.setupRandomTreeModel(popModel, haplotypeModel, Units.Type.YEARS);
		String realHap = hapRunIndex+"_FullHaplotype.fasta";
		Alignment initAlignment = DataImporter.importAlignment(dataDir, realHap);
		
		
		//Tree likelihood
		RealParameter popSize = new RealParameter(new Double[]{0.3});
		popSize.setID("popSize.t");
		
		ConstantPopulation popFunction = new ConstantPopulation();
        popFunction.setInputValue("popSize", popSize);
        popFunction.setID("ConstantPopulation");

        RandomTree tree = new RandomTree();
		tree.initByName("taxa", initAlignment, "populationModel", popFunction);
		tree.setID("Tree.t");

		@Deprecated HaplotypeModel data = haplotypeModel;
		RealParameter freqParameter = new RealParameter(new Double[]{0.25,0.25,0.25,0.25});
		freqParameter.setID("freqParameter.s");
		
		Frequencies estimatedFreqs = new Frequencies();
		estimatedFreqs.initByName("frequencies", freqParameter);
		estimatedFreqs.setID("estimatedFreqs.s");
//        RealParameter parameter = new RealParameter(new Double[] { 1., 1., 1., 1. });
//		DeltaExchangeOperator d = new DeltaExchangeOperator();
//			d.initByName("parameter", parameter);
        
        RealParameter kappa = new RealParameter(new Double[] {27.40259});
        kappa.setID("kappa.s");
        HKY hky = new HKY();
        hky.initByName("kappa", kappa, "frequencies", estimatedFreqs);
        
        SiteModel siteModel = new SiteModel();
        siteModel.initByName("mutationRate", "1.0", "gammaCategoryCount", 1, "substModel", hky);

//        TreeLikelihood likelihood = new TreeLikelihood();
//        likelihood.initByName("data", initAlignment, "tree", tree, "siteModel", siteModel);
        TreeLikelihoodExt likelihood = new TreeLikelihoodExt();
        likelihood.initByName("hapModel", haplotypeModel, "tree", tree, "siteModel", siteModel);
        likelihood.setID("treeLikelihoodExt");
        
        double logP = 0;
        logP = likelihood.calculateLogP();

        //Prior
        Distribution kappaPrior = new Prior();
        LogNormalDistributionModel lnd = new LogNormalDistributionModel();
        lnd.initByName("M", new RealParameter("1.0"), "S", new RealParameter("1.25"));
        kappaPrior.initByName("x", kappa, "distr", lnd);
        kappaPrior.setID("KappaPrior.s");
        
        
        Distribution CoalescentConstant = new Coalescent();
        CoalescentConstant.initByName("populationModel", popFunction, "treeIntervals", new TreeIntervals(tree));
        CoalescentConstant.setID("CoalescentConstant.t");
        
        Distribution popSizePrior = new Prior();
        OneOnX oox = new OneOnX();
        popSizePrior.initByName("x", popSize, "distr", oox);
        popSizePrior.setID("PopSizePrior.t");
        
        //Posterior
        ArrayList<Distribution> posteriorList = new ArrayList<>();
        posteriorList.add(likelihood);
        posteriorList.add(kappaPrior);
        posteriorList.add(CoalescentConstant);
        posteriorList.add(popSizePrior);
        
        CompoundDistribution posterior = new CompoundDistribution();
        posterior.initByName("distribution", posteriorList);
        posterior.setID("posterior");
        
        
        //Operators
        ArrayList<Operator> opsList = new ArrayList<>();
        Operator op;
        
        op = new ScaleOperator();
        op.initByName("parameter", kappa, "scaleFactor", 0.5, "weight", 0.1);
        op.setID("kappaScale.s");
        opsList.add(op);
        
        op = new ScaleOperator();
        op.initByName("tree", tree, "scaleFactor", 0.5, "weight", 3.0);
        op.setID("CoalescentConstantTreeScaler.t");
        opsList.add(op);
        
        op = new ScaleOperator();
        op.initByName("tree", tree, "rootOnly", true, "scaleFactor", 0.5, "weight", 3.0);
        op.setID("CoalescentConstantTreeRootScaler.t");
        opsList.add(op);
        
        op = new Uniform();
        op.initByName("tree", tree, "weight", 30.0);
        op.setID("CoalescentConstantUniformOperator.t");
        opsList.add(op);
        
        op = new SubtreeSlide();
        op.initByName("tree", tree, "weight", 15.0);
        op.setID("CoalescentConstantSubtreeSlide.t");
        opsList.add(op);
        
        op = new Exchange();
        op.initByName("tree", tree, "weight", 15.0);
        op.setID("CoalescentConstantNarrow.t");
        opsList.add(op);
        
        op = new Exchange();
        op.initByName("tree", tree, "isNarrow", false, "weight", 3.0);
        op.setID("CoalescentConstantWide.t");
        opsList.add(op);
        
        op = new WilsonBalding();
        op.initByName("tree", tree, "weight", 3.0);
        op.setID("CoalescentConstantWilsonBalding.t");
        opsList.add(op);
        
        
        op = new ScaleOperator();
        op.initByName("parameter", popSize, "scaleFactor", 0.75, "weight", 3.0);
        op.setID("PopSizeScaler.t");
        opsList.add(op);
        
        op = new DeltaExchangeOperator();
        op.initByName("parameter", freqParameter, "delta", 0.01, "weight", 0.1);
        op.setID("FrequenciesExchanger.s");
        opsList.add(op);        
        
        
        
        
        op = new HaplotypeBaseSingleOperator();
        op.setID("baseSingle");
        op.initByName("hapModel", haplotypeModel, "weight" ,1.0);
//        opsList.add(op);

        OperatorSchedule operatorSchedule = new OperatorSchedule();
//        for (final Operator op : operatorsInput.get()) {
        operatorSchedule.addOperator(op);
//        }
        

        
        Logger logger = new Logger();
        
//        logger.loggersInput.setValue(likelihood, logger);
        logger.setInputValue("log", posterior);
        logger.loggersInput.setValue(likelihood, logger);
//        logger.setInputValue("log", prior);
//        logger.setInputValue("log", srpLikelihood);
        logger.setInputValue("log", kappa);
        logger.setInputValue("log", popSize);
        logger.setInputValue("log", CoalescentConstant);
        logger.setInputValue("log", freqParameter);
        logger.setInputValue("log", kappaPrior);
        
//        logger.setInputValue("model", kappa);
//        logger.setInputValue("model", data);
//        System.out.println(logger.modeInput.getTipText());
//        logger.setInputValue("sort", SORTMODE.smart);
        logger.setInputValue("mode", LOGMODE.compound);
        logger.setInputValue("sanitiseHeaders", true);
        logger.initByName("logEvery", 1000);
        
        
        MCMC mcmc = new MCMC();
        mcmc.setInputValue("operator", opsList);
        mcmc.initByName("chainLength", 10000,
        				"storeEvery", 10,
        				"preBurnin", 0,
//        				"distribution", likelihood,
        				"distribution", posterior,
        				"logger", logger
//        				"operatorschedule", operatorSchedule
//        				"operator"
        				);
		mcmc.run();
		
		logP = srpLikelihood.calculateLogP();
		System.out.println(logP);
		
	}
}
