package srp.core;

import java.io.File;
import java.io.FileOutputStream;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;

import cern.colt.Arrays;
import srp.dr.ext.TreeLikelihoodExt;
import srp.evolution.haplotypes.HaplotypeLoggerWithTrueHaplotype;
import srp.evolution.haplotypes.HaplotypeModel;
import srp.evolution.likelihood.AbstractShortReadsLikelihood;
import srp.evolution.likelihood.haplotypes.ShortReadsHaplotypeLikelihood;
import srp.evolution.shortreads.ShortReadMapping;
import dr.app.pathogen.TemporalRooting;
import dr.evolution.alignment.Alignment;
import dr.evolution.alignment.Patterns;
import dr.evolution.datatype.Nucleotides;
import dr.evolution.distance.DistanceMatrix;
import dr.evolution.distance.F84DistanceMatrix;
import dr.evolution.tree.NeighborJoiningTree;
import dr.evolution.tree.SimpleTree;
import dr.evolution.tree.Tree;
import dr.evolution.tree.UPGMATree;
import dr.evolution.util.TaxonList;
import dr.evolution.util.Units;
import dr.evomodel.branchratemodel.StrictClockBranchRates;
import dr.evomodel.coalescent.CoalescentLikelihood;
import dr.evomodel.coalescent.ConstantPopulationModel;
import dr.evomodel.sitemodel.GammaSiteModel;
import dr.evomodel.substmodel.FrequencyModel;
import dr.evomodel.substmodel.HKY;
import dr.evomodel.tree.TreeLogger;
import dr.evomodel.tree.TreeModel;
import dr.evomodel.treelikelihood.TreeLikelihood;
import dr.evomodelxml.coalescent.ConstantPopulationModelParser;
import dr.evomodelxml.substmodel.HKYParser;
import dr.inference.loggers.MCLogger;
import dr.inference.loggers.TabDelimitedFormatter;
import dr.inference.mcmc.MCMC;
import dr.inference.mcmc.MCMCOptions;
import dr.inference.model.Likelihood;
import dr.inference.model.Parameter;
import dr.inference.operators.MCMCOperator;
import dr.inference.operators.OperatorSchedule;
import dr.inference.operators.SimpleOperatorSchedule;
import dr.inferencexml.model.CompoundLikelihoodParser;

import org.apache.*;

import beast.evolution.likelihood.BeagleTreeLikelihood;
public class MainMCMCHaplotype {

	public static void main(String[] args) throws Exception {
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

		
		BeagleTreeLikelihood likelihood = new BeagleTreeLikelihood();
        likelihood.initByName("data", data, "tree", tree, "siteModel", siteModel);
        TreeLikelihood likelihood = new TreeLikelihood();
        likelihood.initByName("data", data, "tree", tree, "siteModel", siteModel);
        logP += likelihood.calculateLogP();

        // Now calculate again, with shuffled taxon order
        Collections.shuffle(data.sequenceInput.get());
        likelihood = new TreeLikelihood();
        likelihood.initByName("data", data, "tree", tree, "siteModel", siteModel);
        shuffledLogP += likelihood.calculateLogP();
		
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
		
		if(args.length == 7){
			dataDir = args[0];
			runIndex = Integer.parseInt(args[1]);
			totalSamples = Integer.parseInt(args[2]);
			logInterval = Integer.parseInt(args[3]);
			noOfTrueHaplotype = Integer.parseInt(args[4]);
			noOfRecoveredHaplotype= Integer.parseInt(args[5]);
			inputReadSuffix= (args[6]);
			if(!inputReadSuffix.equals("ART") && !inputReadSuffix.equals("ART_errFree") && !inputReadSuffix.equals("Shrimp") ){
				System.out.println("Invalid input: InputSuffix should be one of [ART|ART_errFree|Shrimp]");
				System.exit(-3);
			}
		}
		
		else{	
			System.out.println("local parameters");
			isLocal = true;
			dataDir = "/home/steven/workspaceSrp/snowgoose/srp/unittest/testData/";
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
		}
		
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

		Alignment shortReads = dataImporter.importShortReads(shortReadFile);
		ShortReadMapping srpMap = new ShortReadMapping(shortReads);

//		srpMap.summary();
//		for (int i = 0; i < 10; i++) {
//			System.out.println((char)srpMap.nextBaseAt(100));
//		}
//		System.exit(10);
		
		HaplotypeModel haplotypeModel = null;
		
		if(randomHaplotype){
//			haplotypeModel = new HaplotypeModel(noOfRecoveredHaplotype, shortReads.getSiteCount());
			haplotypeModel = new HaplotypeModel(noOfRecoveredHaplotype, srpMap);
//			haplotypeModel = new HaplotypeModel(noOfRecoveredHaplotype, shortReads.getSiteCount());
//			haplotypeModel.addShortReadMap(srpMap);
		}
		else{
//			String partialHaplotypeName = prefix+".haplotypepartial";
			String partialHaplotypeName = hapRunIndex+"_FullHaplotype.fasta";
//			Alignment trueAlignment = dataImporter.importAlignment(trueHaplotypeFile);
			Alignment trueAlignment = dataImporter.importAlignment(partialHaplotypeName);
			haplotypeModel = new HaplotypeModel(trueAlignment);
			haplotypeModel.addShortReadMap(srpMap);
//			haplotypeModel = dataImporter.importPartialSpectrumFile(partialHaplotypeName);
		}

		// ShortReadLikelihood
		ShortReadsHaplotypeLikelihood srpLikelihood = new ShortReadsHaplotypeLikelihood(haplotypeModel, srpMap);
		System.out.println("Error rate: "+srpLikelihood.ERROR_RATE);
		
		// coalescent
		Parameter popSize = new Parameter.Default(ConstantPopulationModelParser.POPULATION_SIZE, 3000.0, 100, 100000.0);

		// Random treeModel
		ConstantPopulationModel popModel = new ConstantPopulationModel(popSize, Units.Type.YEARS);
//		TreeModel treeModel = MCMCSetupHelperHaplotype.setupRandomTreeModel(popModel, haplotypeModel, Units.Type.YEARS);
		TreeModel treeModel;
		if(randomTree){
//			treeModel = MCMCSetupHelperHaplotype.setupRandomTreeModel(popModel, haplotypeModel, Units.Type.YEARS);
			
            DistanceMatrix distances = new F84DistanceMatrix(haplotypeModel);
            Tree tree = new UPGMATree(distances);
            treeModel = new TreeModel(tree);			
		}
		else{
//			String partialTreeName = prefix+".treespartial";
//			String partialTreeName = prefix+".trees";
			String partialTreeName = hapRunIndex+"_FullTree.tree";
			System.out.println("load tree: "+ partialTreeName);
			Tree partialPhylogeny = dataImporter.importTree(partialTreeName);
			treeModel = new TreeModel(TreeModel.TREE_MODEL, partialPhylogeny, false);
		}
		
		// Coalescent likelihood
		CoalescentLikelihood coalescent = new CoalescentLikelihood(treeModel,null, new ArrayList<TaxonList>(), popModel);
		coalescent.setId("coalescent");

		// Simulate haplotypes, treeLikelihood
		HashMap<String, Object> parameterList = MCMCSetupHelperHaplotype.setupTreeLikelihoodHaplotypeModel(treeModel, haplotypeModel, srpMap);

		Parameter kappa = (Parameter) parameterList.get("kappa");
		Parameter freqs = (Parameter) parameterList.get("freqs");
		StrictClockBranchRates branchRateModel = (StrictClockBranchRates) parameterList.get("branchRateModel");
		TreeLikelihoodExt treeLikelihood = (TreeLikelihoodExt) parameterList.get("treeLikelihood");
//		TreeLikelihoodExt treeLikelihood = null; 

		// CompoundLikelihood
		HashMap<String, Likelihood> compoundlikelihoods = MCMCSetupHelperHaplotype
				.setupCompoundLikelihood(popSize, kappa, coalescent,
						treeLikelihood, srpLikelihood);
		Likelihood prior = compoundlikelihoods.get(CompoundLikelihoodParser.PRIOR);
		Likelihood likelihood = compoundlikelihoods.get(CompoundLikelihoodParser.LIKELIHOOD);
		Likelihood shortReadLikelihood = compoundlikelihoods.get(AbstractShortReadsLikelihood.SHORT_READ_LIKELIHOOD);
		Likelihood posterior = compoundlikelihoods.get(CompoundLikelihoodParser.POSTERIOR);

		
		Parameter rootHeight = treeModel.getRootHeightParameter();
		rootHeight.setId("rootHeight");
		// Operators
		OperatorSchedule schedule = new SimpleOperatorSchedule();

		MCMCSetupHelperHaplotype.defalutTreeOperators(schedule, treeModel);
		double totalTree = 0;
		for (int i = 0; i < schedule.getOperatorCount(); i++) {
			MCMCOperator op= schedule.getOperator(i);
			totalTree += op.getWeight() ;
		}
		
		
		MCMCSetupHelperHaplotype.defalutOperators(schedule, haplotypeModel, freqs, popSize, kappa);
		double total = 0;
		System.out.println("All Operators:");
		for (int i = 0; i < schedule.getOperatorCount(); i++) {
			MCMCOperator op= schedule.getOperator(i);
			System.out.println(op.getOperatorName() +"\t"+ op.getWeight());
			total += op.getWeight() ;
		}
		System.out.println("\n===========\ntotal Tree Weight: "+totalTree);
		System.out.println("total non-Tree Weight: "+ (total-totalTree));
		

		// MCLogger
		MCLogger[] loggers = new MCLogger[4];
		// log tracer
		loggers[0] = new MCLogger(logTracerName, logInterval, false, 0);
		MCMCSetupHelperHaplotype.addToLogger(loggers[0], 
				posterior, 
				prior, 
				likelihood, 
				shortReadLikelihood,
				rootHeight, 
////				//rateParameter,
				popSize, kappa, coalescent,
				freqs
				);
		// System.out
		loggers[1] = new MCLogger(new TabDelimitedFormatter(System.out), logInterval, true, logInterval*2);
		MCMCSetupHelperHaplotype.addToLogger(loggers[1],
//				
				posterior, prior,
				likelihood, 
				shortReadLikelihood,
				popSize, kappa, coalescent, rootHeight
//				,freqs
				);
		// log Tree
		TabDelimitedFormatter treeFormatter = new TabDelimitedFormatter(
				new PrintWriter(new FileOutputStream(new File(logTreeName))));

		loggers[2] = new TreeLogger(treeModel, branchRateModel, null, null,
				treeFormatter, logInterval, true, true, true, null, null);

		// log Haplotype
		
		Alignment trueAlignment = dataImporter.importAlignment(trueHaplotypeFile);
		HaplotypeModel trueHaplotypeModel = new HaplotypeModel(trueAlignment);
		ShortReadsHaplotypeLikelihood trueSrp = new ShortReadsHaplotypeLikelihood(trueHaplotypeModel, srpMap);
		System.out.println("True haplotype likelihood: "+trueSrp.getLogLikelihood());

//		AlignmentMapping alignmentMapping = new AlignmentMapping(shortReads);
//		ShortReadsHaplotypeLikelihood trueSrp = new ShortReadsHaplotypeLikelihood(HaplotypeModel.factory(shortReads, trueAlignment), srpMap);
//		System.err.println("\'trueShortReadLikelihood\': "+trueSrp.getLogLikelihood());
		loggers[3] = new HaplotypeLoggerWithTrueHaplotype(haplotypeModel,
				trueAlignment, logHaplotypeName, (int) (logInterval
						* totalSamples * 0.1));
		// loggers[3] = new TreeLogger(treeModel, branchRateModel, null, null,
//				treeFormatter, logInterval, true, true, true, null, null);
		
		// MCMC
		MCMCOptions options;
		if(isLocal ){
			options = MCMCSetupHelper.setMCMCOptions(logInterval, totalSamples, 0, 0);
		}
		else{
			options = MCMCSetupHelper.setMCMCOptions(logInterval, totalSamples);
		}
		
		MCMC mcmc = new MCMC("mcmc1");
		mcmc.setShowOperatorAnalysis(true);
		mcmc.setOperatorAnalysisFile(new File(operatorAnalysisFile));
		
		mcmc.init(options, posterior, schedule, loggers);
		mcmc.run();

		{		
			System.out.println(mcmc.getTimer().toString());
			System.out.println("True haplotype Likelihood: "+trueSrp.getLogLikelihood());
			srpLikelihood.makeDirty();
			System.out.println("MCMC haplotype likelihood: "+srpLikelihood.getLogLikelihood());
			
			//		treeModel
			TreeLikelihood reCalTreeLikelihood = new TreeLikelihood(
					haplotypeModel, treeModel, treeLikelihood.getSiteModel(), treeLikelihood.getBranchRateModel(), null,
					false, false, true, false, false);
			System.out.println("Tree treeLikelihood: "+reCalTreeLikelihood.getLogLikelihood());

		}
	
	}


}
