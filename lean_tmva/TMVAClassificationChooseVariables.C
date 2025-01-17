// @(#)root/tmva $Id$
/**********************************************************************************
 * Project   : TMVA - a ROOT-integrated toolkit for multivariate data analysis    *
 * Package   : TMVA                                                               *
 * Root Macro: TMVAClassification                                                 *
 *                                                                                *
 * This macro provides examples for the training and testing of the               *
 * TMVA classifiers.                                                              *
 *                                                                                *
 * As input data is used a toy-MC sample consisting of four Gaussian-distributed  *
 * and linearly correlated input variables.                                       *
 *                                                                                *
 * The methods to be used can be switched on and off by means of booleans, or     *
 * via the prompt command, for example:                                           *
 *                                                                                *
 *    root -l ./TMVAClassification.C\(\"Fisher,Likelihood\"\)                     *
 *                                                                                *
 * (note that the backslashes are mandatory)                                      *
 * If no method given, a default set of classifiers is used.                      *
 *                                                                                *
 * The output file "TMVA.root" can be analysed with the use of dedicated          *
 * macros (simply say: root -l <macro.C>), which can be conveniently              *
 * invoked through a GUI that will appear at the end of the run of this macro.    *
 * Launch the GUI via the command:                                                *
 *                                                                                *
 *    root -l ./TMVAGui.C                                                         *
 *                                                                                *
 **********************************************************************************/

#include <cstdlib>
#include <iostream>
#include <map>
#include <string>
#include <unistd.h>
#include <fstream>
#include <cstdlib>
#include <set>
#include <algorithm>
#include <ctime>

#include "TChain.h"
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TObjString.h"
#include "TSystem.h"
#include "TROOT.h"

#include "TMVAGui.C"
#include "TMVA/Factory.h"
#include "TMVA/Tools.h"

using namespace TMVA;

typedef std::pair<double, std::pair<long long, std::string> > method_stats;
method_stats make_method_stats(long long v, double r, std::string method) {
	     return std::make_pair(r, make_pair(v, method));
}

void TMVAClassificationChooseVariables(int MVA_type, string variable_file)
{
	//---------------------------------------------------------------
	// This loads the library
	Tools::Instance();
	gROOT->LoadMacro("mydyncastcopy.C+");

	// Default MVA methods to be trained + tested
	std::map<std::string,int> Use;

	//
	// --- Neural Networks (all are feed-forward Multilayer Perceptrons)
	Use["MLP"]             = bool(1&MVA_type); // Recommended ANN
	// 
	// --- Boosted Decision Trees
	Use["BDTG"]            = bool(2&MVA_type);
	// ---------------------------------------------------------------

	std::cerr << std::endl;
	std::cerr << "==> Start TMVAClassification" << std::endl;

	// --------------------------------------------------------------------------------------------------

	// Read training and test data  
	TString data_path = "/mnt/xrootdb/alister/MVA_studies/samples/"; 
	TFile *input_sig_el = TFile::Open(data_path + "nominal_el/tprime_650_1M.root");
	TFile *input_sig_mu = TFile::Open(data_path + "nominal_mu/tprime_650_1M.root");
	TFile *input_ttbar_el = TFile::Open(data_path + "nominal_el/ttbar.root");
	TFile *input_ttbar_mu = TFile::Open(data_path + "nominal_mu/ttbar.root");
	TFile *input_wjets_el = TFile::Open(data_path + "nominal_el/wjets.root");
	TFile *input_wjets_mu = TFile::Open(data_path + "nominal_mu/wjets.root");
	TFile *input_zjets_el = TFile::Open(data_path + "nominal_el/zjets.root");
	TFile *input_zjets_mu = TFile::Open(data_path + "nominal_mu/zjets.root");
	TFile *input_singletop_el = TFile::Open(data_path + "nominal_el/singletop.root");
	TFile *input_singletop_mu = TFile::Open(data_path + "nominal_mu/singletop.root");
	TFile *input_diboson_el = TFile::Open(data_path + "nominal_el/diboson.root");
	TFile *input_diboson_mu = TFile::Open(data_path + "nominal_mu/diboson.root");

	// --- Register the training and test trees
	TTree *signal_el  = (TTree*)input_sig_el->Get("mini");
	TTree *signal_mu  = (TTree*)input_sig_mu->Get("mini");
	TTree *ttbar_el = (TTree*)input_ttbar_el->Get("mini");
	TTree *ttbar_mu = (TTree*)input_ttbar_mu->Get("mini");
	TTree *wjets_el = (TTree*)input_wjets_el->Get("mini");
	TTree *wjets_mu = (TTree*)input_wjets_mu->Get("mini");
	TTree *zjets_el = (TTree*)input_zjets_el->Get("mini");
	TTree *zjets_mu = (TTree*)input_zjets_mu->Get("mini");
	TTree *singletop_el = (TTree*)input_singletop_el->Get("mini");
	TTree *singletop_mu = (TTree*)input_singletop_mu->Get("mini");
	TTree *diboson_el = (TTree*)input_diboson_el->Get("mini");
	TTree *diboson_mu = (TTree*)input_diboson_mu->Get("mini");

	// --- Here the preparation phase begins
	std::vector<std::vector<TString> > variables; //each variable set specified in a 4-tuple.
	std::fstream fin(variable_file.c_str(), std::fstream::in);
	std::vector<TString> inp(4); //var, title, unit, type.
	while (fin >> inp[0] >> inp[1] >> inp[2] >> inp[3]) {
		variables.push_back(inp);
	}

	std::set<method_stats> rankings;

	// Create a ROOT output file where TMVA will store ntuples, histograms, etc.
	TString outfileName( "TMVAChosenVariables.root" );
	TFile* outputFile = TFile::Open( outfileName, "RECREATE" );

	// Create the factory object. Later you can choose the methods
	// whose performance you'd like to investigate. The factory is 
	// the only TMVA object you have to interact with
	//
	// The first argument is the base of the name of all the
	// weightfiles in the directory weight/
	//
	// The second argument is the output file for the training results
	// All TMVA output can be suppressed by removing the "!" (not) in
	// front of the "Silent" argument in the option string
	Factory *factory = new Factory( "TMVAClassification", outputFile,
			"!V:Silent:Color:DrawProgressBar:AnalysisType=Classification" );

	std::cerr << std::endl;
	std::cerr << "================================================" << std::endl;

	for (int i = 0; i < variables.size(); i++) {
		const std::vector<TString>& tup = variables[i];
		factory->AddVariable(tup[0], tup[1], tup[2], tup[3][0]);
		std::cerr << "Adding variable: " << tup[1] << std::endl;
	}
	std::cerr << "================================================" << std::endl;

	// global event weights per tree (see below for setting event-wise weights)
	Double_t signalWeight     = 1.0;
	Double_t backgroundWeight = 1.0;

	// You can add an arbitrary number of signal or background trees
	factory->AddSignalTree    ( signal_el,     signalWeight     );
	factory->AddSignalTree    ( signal_mu,     signalWeight     );
	factory->AddBackgroundTree( ttbar_el, backgroundWeight );
	factory->AddBackgroundTree( ttbar_mu, backgroundWeight );
	factory->AddBackgroundTree( wjets_el, backgroundWeight );
	factory->AddBackgroundTree( wjets_mu, backgroundWeight );
	factory->AddBackgroundTree( zjets_el, backgroundWeight );
	factory->AddBackgroundTree( zjets_mu, backgroundWeight );
	factory->AddBackgroundTree( singletop_el, backgroundWeight );
	factory->AddBackgroundTree( singletop_mu, backgroundWeight );
	factory->AddBackgroundTree( diboson_el, backgroundWeight );
	factory->AddBackgroundTree( diboson_mu, backgroundWeight );

	// Set individual event weights (the variables must exist in the original TTree)
	factory->SetSignalWeightExpression    ("weight_1btin_70*weight_70/btweight_70");
	factory->SetBackgroundWeightExpression("weight_1btin_70*weight_70/btweight_70");

	// Apply additional cuts on the signal and background samples (can be different)
	TCut mycuts = "weight_70>0 && analysis_channel==1"; // for example: TCut mycuts = "abs(var1)<0.5 && abs(var2-0.5)<1";
	TCut mycutb = "weight_70>0 && analysis_channel==1"; // for example: TCut mycutb = "abs(var1)<0.5";

	// Tell the factory how to use the training and testing events
	//
	// If no numbers of events are given, half of the events in the tree are used 
	// for training, and the other half for testing:
	//    factory->PrepareTrainingAndTestTree( mycut, "SplitMode=random:!V" );
	// To also specify the number of testing events, use:
	//    factory->PrepareTrainingAndTestTree( mycut,
	//                                         "NSigTrain=3000:NBkgTrain=3000:NSigTest=3000:NBkgTest=3000:SplitMode=Random:!V" );
	factory->PrepareTrainingAndTestTree( mycuts, mycutb,
			"nTrain_Signal=0:nTrain_Background=0:SplitMode=Random:NormMode=NumEvents:!V" );

	// ---- Book MVA methods
	//
	// Please lookup the various method configuration options in the corresponding cxx files, eg:
	// src/MethoCuts.cxx, etc, or here: http://tmva.sourceforge.net/optionRef.html


	// TMVA ANN: MLP (recommended ANN) -- all ANNs in TMVA are Multilayer Perceptrons
	if (Use["MLP"])
		factory->BookMethod( Types::kMLP, "MLP", "!H:!V:NeuronType=tanh:VarTransform=N:NCycles=600:HiddenLayers=N+5:TestRate=5:!UseRegulator" );

	// Boosted decision trees
	// In later versions of TMVA, use UsedBaggedBoost and BaggedSampleFraction
	// Instead of UsedBaggedGrad and GradBaggingFraction
	if (Use["BDTG"])
		factory->BookMethod( Types::kBDT, "BDTG",
				"!H:!V:NTrees=2000::BoostType=Grad:Shrinkage=0.1:UseBaggedGrad:GradBaggingFraction=0.5:nCuts=20:MaxDepth=3:MaxDepth=4" );

	// ---- Now you can tell the factory to train, test, and evaluate the MVAs

	// Train MVAs using the set of training events
	factory->TrainAllMethods();

	// ---- Evaluate all MVAs using the set of test events
	factory->TestAllMethods();

	// ----- Evaluate and compare performance of all configured MVAs
	factory->EvaluateAllMethods();

	// --------------------------------------------------------------

	// Save the output
	outputFile->Close();

	std::cerr << "==> Wrote root file: " << outputFile->GetName() << std::endl;
	std::cerr << "==> TMVAClassification is done!" << std::endl;

	//get efficiency of methods, compare with current bests
	IMethod* i_met;
	MethodBase* method;
	if (Use["MLP"]) {
		i_met = factory->GetMethod("MLP");
		method = dyncast(i_met);
		std::cerr << "MLP Significance: " << method->GetSignificance() <<  std::endl;
		rankings.insert(make_method_stats(variable_choice, method->GetSignificance(), "MLP"));
	}

	if (Use["BDTG"]) {
		i_met = factory->GetMethod("BDTG");
		method = dyncast(i_met);
		std::cerr << "BDTG Significance: " << method->GetSignificance() <<  std::endl;
		rankings.insert(make_method_stats(variable_choice, method->GetSignificance(), "BDTG"));
	}

	delete factory;

	std::set<method_stats>::iterator it;
	std::cout << "Best variables:" << std::endl;
	for (it = rankings.begin(); it != rankings.end(); it++) {
		std::cout << "================================================" << std::endl;
		std::cout << "Significance: " << it->first << std::endl;
		std::cout << "Method Name: " << (it->second).second << std::endl;
		std::cout << "variables used: ";
		for (int i = 0; i < variables.size(); i++) {
			if ((it->second).first & (1LL << i)) std::cout << "[ " << variables[i][1] << " ] ";
		}
		cout << endl;
		std::cout << "================================================" << std::endl;
	}

	// Launch the GUI for the root macros
	//if (!gROOT->IsBatch()) TMVAGui( outfileName );
	gApplication->Terminate(0);
}

