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
#include <cmath>
#include <sstream>

#include "TChain.h"
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TObjString.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TStopwatch.h"

#include "TMVAGui.C"
#include "TMVA/Factory.h"
#include "TMVA/Tools.h"
#include "TMVA/Reader.h"
#include "TMVA/Config.h"
#include "TMVA/Tools.h"
#include "TMVA/MethodCuts.h"

using namespace TMVA;

typedef std::pair<double, std::pair<long long, std::string> > method_stats;
method_stats make_method_stats(long long v, double r, std::string method) {
	return std::make_pair(r, make_pair(v, method));
}


//void TMVAClassification(int MVA_type) {
int main(int argc, char* argv[]) {
	srand(time(NULL));
	//---------------------------------------------------------------
	// This loads the library
	Tools::Instance();

	if (argc == 1) {
		std::cerr << "Filename of input variables needed!" << std::endl;
		return 0;
	}

	//printf("%s\n", argv[1]);
	std::fstream fin(argv[1], std::fstream::in);

	int MVA_type;
	if (argc > 2) MVA_type = atoi(argv[2]);
	else MVA_type = 3;
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


	TFile *input_sig_test_el = TFile::Open(data_path + "nominal_el/tprime_650.root", "READ");
	TFile *input_sig_test_mu = TFile::Open(data_path + "nominal_mu/tprime_650.root", "READ");

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

	TTree *signal_test_el  = (TTree*)input_sig_test_el->Get("mini");
	TTree *signal_test_mu  = (TTree*)input_sig_test_mu->Get("mini");

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
	std::vector<TString> inp(4); //var, title, unit, type.
	while (fin >> inp[0] >> inp[1] >> inp[2] >> inp[3]) {
		variables.push_back(inp);
	}

	std::set<method_stats> rankings;
	int num_used = variables.size();
	for (int num_nodes = 20; num_nodes <= 20; num_nodes++) {
		std::stringstream ss;
		ss << num_nodes;
		std::string str_nodes = ss.str();
		if (num_nodes >= 0) str_nodes = "+"+str_nodes;
		// Create a ROOT output file where TMVA will store ntuples, histograms, etc.
		TString outfileName( "TMVA.root" );
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
				"!V:Silent:Color:DrawProgressBar:Transformations=I;D;P;G,D:AnalysisType=Classification" );

		std::cerr << std::endl;
		std::cerr << "================================================" << std::endl;
		//add a random num_used sized subset of variables to train on.

		for (int i = 0; i < num_used; i++) {
			if (variables[i][1] == "analysis_channel") continue;
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
		if (Use["MLP"]) {
			/*
			std::string layers = "N";
			for (int r = 0; r < num_nodes; r++) {
				layers += ",N";
			}*/
			factory->BookMethod( Types::kMLP, "MLP", "!H:!V:NeuronType=tanh:VarTransform=N:NCycles=700:HiddenLayers=N+"+str_nodes+":TestRate=5:!UseRegulator:SamplingTraining=False:EstimatorType=CE" );
		}
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

		//std::cerr << "==> Wrote root file: " << outputFile->GetName() << std::endl;
		//std::cerr << "==> TMVAClassification is done!" << std::endl;

		delete factory;

		//std::cerr << "Pausing before application stage..." << std::endl;
		//sleep(10);

		//=================================== Begin Application =====================================//

		// Create the Reader Object

		Float_t weight_70, btweight_70, weight_1btin_70;
		Int_t analysis_channel;

		Reader *reader = new TMVA::Reader( "!Color:!Silent" );	
		std::vector<Float_t> var_val(num_used);

		int cnt = 0;
		for (int i = 0; i < num_used; i++) {
			if (variables[i][1] == "analysis_channel") continue;
			const std::vector<TString>& tup = variables[i];
			//std::cerr << "Adding variable: " << tup[1] << std::endl;
			reader->AddVariable( tup[0], &var_val[cnt++] );
		}

		// --- Book the MVA methods
		if (Use["MLP"]) {
			TString methodName = TString("MLP method");
			reader->BookMVA( methodName, "weights/TMVAClassification_MLP.weights.xml" );
		}

		if (Use["BDTG"]) {
			TString methodName = TString("BDTG method");
			reader->BookMVA( methodName, "weights/TMVAClassification_BDTG.weights.xml" );
		}

		TTree* background[10] = {ttbar_el, ttbar_mu, wjets_el, wjets_mu, zjets_el, zjets_mu, 
			singletop_el, singletop_mu, diboson_el, diboson_mu};

		const int cuts = 10;
		double l = 0.99, r = 0.999;
		double del = (r-l)/cuts;

		double max_mlp_sig = 0, max_bdt_sig = 0;
		for (double mva_cut = l; mva_cut < r; mva_cut += del) {
			std::cerr << "trying cut value " << mva_cut << std::endl;
			double tot_bg_mlp = 0, tot_bg_bdtg = 0;
			for (int bg = 0; bg < 10; bg++) {
				// Prepare input tree (this must be replaced by your data source)
				// in this example, there is a toy tree with signal and one with background events
				// we'll later on use only the "signal" events for the test in this example.
				//

				//std::cerr << "--- Testing background sample ---" << std::endl;

				TTree* theTree = background[bg];	

				cnt = 0;
				for (int i = 0; i < num_used; i++) {
					if (variables[i][1] == "analysis_channel") continue;
					const std::vector<TString>& tup = variables[i];
					//std::cerr << "Adding variable: " << tup[1] << std::endl;
					theTree->SetBranchAddress( tup[0], &var_val[cnt++] );
				}

				theTree->SetBranchAddress( "weight_70", &weight_70);
				theTree->SetBranchAddress( "btweight_70", &btweight_70);
				theTree->SetBranchAddress( "weight_1btin_70", &weight_1btin_70);
				theTree->SetBranchAddress( "analysis_channel", &analysis_channel);
				//std::cerr << "--- Processing: " << theTree->GetEntries() << " events" << std::endl;
				//TStopwatch sw;
				//sw.Start();

				Int_t nEvent = theTree->GetEntries();

				for (Long64_t ievt=0; ievt<nEvent; ievt++) {
					//if (ievt%10000 == 0) {
					//	std::cerr << "--- ... Processing event: " << ievt << std::endl;
					//}

					theTree->GetEntry(ievt);
					if (analysis_channel == 0) continue;

					//theTree->Show(ievt);

					Double_t mc_weight = weight_1btin_70*weight_70/btweight_70;
					// --- Return the MVA outputs and fill into histograms
					if (Use["MLP"]) {
						Double_t MVA_MLP = reader->EvaluateMVA( "MLP method" );
						if (MVA_MLP >= mva_cut) {
							//std::cerr << "got one! " << mc_weight << std::endl; 
							tot_bg_mlp += mc_weight;
						}
					}

					if (Use["BDTG"]) {
						Double_t MVA_BDTG = reader->EvaluateMVA("BDTG method");	
						//std::cerr << MVA_BDTG << std::endl;
						if ((MVA_BDTG+1)/2.0 >= mva_cut) {
							//std::cerr << "got one!" << mc_weight << std::endl;
							tot_bg_bdtg += mc_weight;
						}
					}
				}

				theTree->ResetBranchAddresses();
				// Get elapsed time
				//sw.Stop();
				//std::cerr << "--- End of event loop: "; //sw.Print();
				//std::cerr << "==> TMVAClassificationApplication is done!" << std::endl << std::endl;
			}

			TTree* signal[2] = {signal_test_el, signal_test_mu};

			double signal_mlp = 0, signal_bdtg = 0;
			for (int s = 0; s < 2; s++) {
				// Prepare input tree (this must be replaced by your data source)
				// in this example, there is a toy tree with signal and one with background events
				// we'll later on use only the "signal" events for the test in this example.
				//

				//std::cerr << "--- Testing signal sample ---" << std::endl;

				TTree* theTree = signal[s];	

				cnt = 0;
				for (int i = 0; i < num_used; i++) {
					if (variables[i][1] == "analysis_channel") continue;
					const std::vector<TString>& tup = variables[i];
					//std::cerr << "Adding variable: " << tup[1] << std::endl;
					theTree->SetBranchAddress( tup[0], &var_val[cnt++] );
				}

				theTree->SetBranchAddress( "weight_70", &weight_70);
				theTree->SetBranchAddress( "btweight_70", &btweight_70);
				theTree->SetBranchAddress( "weight_1btin_70", &weight_1btin_70);
				theTree->SetBranchAddress( "analysis_channel", &analysis_channel);
				//std::cerr << "--- Processing: " << theTree->GetEntries() << " events" << std::endl;
				//TStopwatch sw;
				//sw.Start();



				Int_t nEvent = theTree->GetEntries();

				for (Long64_t ievt=0; ievt<nEvent; ievt++) {
					//if (ievt%10000 == 0) {
					//	std::cerr << "--- ... Processing event: " << ievt << std::endl;
					//}

					theTree->GetEntry(ievt);
					if (analysis_channel == 0) continue;

					//theTree->Show(ievt);

					Double_t mc_weight = weight_1btin_70*weight_70/btweight_70;
					// --- Return the MVA outputs and fill into histograms
					if (Use["MLP"]) {
						Double_t MVA_MLP = reader->EvaluateMVA( "MLP method" );
						if (MVA_MLP >= mva_cut) {
							//std::cerr << "got one! " << mc_weight << std::endl; 
							signal_mlp += mc_weight;
						}
					}

					if (Use["BDTG"]) {
						Double_t MVA_BDTG = reader->EvaluateMVA("BDTG method");	
						//std::cerr << MVA_BDTG << std::endl;
						if ((MVA_BDTG+1)/2.0 >= mva_cut) {
							//std::cerr << "got one!" << mc_weight << std::endl;
							signal_bdtg += mc_weight;
						}
					}
				}

				// Get elapsed time
				//sw.Stop();
				theTree->ResetBranchAddresses();
				//std::cerr << "--- End of event loop: "; //sw.Print();
				//std::cerr << "==> TMVAClassificationApplication is done!" << std::endl << std::endl;


				//if (ievt % 10000 == 0)
				//	std::cerr << signal_mlp << " " << std::signal_bdtg < endl;
			}

			if (Use["MLP"]) {
				std::cout << "cut at " << mva_cut << " gives " << signal_mlp/sqrt(signal_mlp + tot_bg_mlp) << std::endl;
				max_mlp_sig = std::max(max_mlp_sig, signal_mlp/sqrt(signal_mlp+tot_bg_mlp));
			}

			if (Use["BDTG"]) {
				max_bdt_sig = std::max(max_bdt_sig, signal_bdtg/sqrt(signal_bdtg+tot_bg_bdtg));
			}
		}

		std::cerr << "Num nodes: N" << num_nodes << std::endl;
		std::cerr << "MLP, BDT: " << max_mlp_sig << ", " << max_bdt_sig << std::endl;
		//get efficiency of methods, compare with current bests
		if (Use["MLP"]) {
			std::cerr << "MLP Significance: " << max_mlp_sig <<  std::endl;
			rankings.insert(make_method_stats(num_nodes, max_mlp_sig, "MLP"));
		}

		if (Use["BDTG"]) {
			std::cerr << "BDTG Significance: " << max_bdt_sig <<  std::endl;
			rankings.insert(make_method_stats(num_nodes, max_bdt_sig, "BDTG"));
		}
	}

	std::set<method_stats>::iterator it;
	std::cout << "Best variables:" << std::endl;
	for (it = rankings.begin(); it != rankings.end(); it++) {
		std::cout << "================================================" << std::endl;
		std::cout << "Significance: " << it->first << std::endl;
		std::cout << "Method Name: " << (it->second).second << std::endl;
		std::cout << "Num trees: " << (it->second).first << std::endl;
		std::cout << "Variables used: ";
		for (int i = 0; i < variables.size(); i++) {
			std::cout << "[ " << variables[i][1] << " ] ";
		}
		cout << endl;
		std::cout << "================================================" << std::endl;
	}

	// Launch the GUI for the root macros
	//if (!gROOT->IsBatch()) TMVAGui( outfileName );
	//gApplication->Terminate(0);

}
