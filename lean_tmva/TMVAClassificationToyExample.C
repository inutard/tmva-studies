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

long long random_ksubset(int n, int k) {
    std::vector<int> pool;
    for (int i = 0; i < n; i++) pool.push_back(i);
    long long res = 0;
    while (k--) {
        int cand = rand() % pool.size();
        res |= (1LL << pool[cand]);
        std::swap(pool[cand], pool[pool.size()-1]);
        pool.pop_back();
    }
    return res;
}

typedef std::pair<double, std::pair<long long, std::string> > method_stats;
method_stats make_method_stats(long long v, double r, std::string method) {
    return std::make_pair(r, make_pair(v, method));
}

void TMVAClassificationToyExample()
{
	// The explicit loading of the shared libTMVA is done in TMVAlogon.C, defined in .rootrc
	// if you use your private .rootrc, or run from a different directory, please copy the
	// corresponding lines from .rootrc

	// methods to be processed can be given as an argument; use format:
	//
	// mylinux~> root -l TMVAClassification.C\(\"myMethod1,myMethod2,myMethod3\"\)
	//
	// if you like to use a method via the plugin mechanism, we recommend using
	//
	// mylinux~> root -l TMVAClassification.C\(\"P_myMethod\"\)
	// (an example is given for using the BDT as plugin (see below),
	// but of course the real application is when you write your own
	// method based)

	//---------------------------------------------------------------
	// This loads the library
	Tools::Instance();
    gROOT->LoadMacro("mydyncast.C+");
    
	// Default MVA methods to be trained + tested
	std::map<std::string,int> Use;

	//
	// --- Neural Networks (all are feed-forward Multilayer Perceptrons)
	Use["MLP"]             = 1; // Recommended ANN
	// 
	// --- Boosted Decision Trees
	Use["BDTG"]            = 0;
	// ---------------------------------------------------------------

	std::cout << std::endl;
	std::cout << "==> Start TMVAClassification" << std::endl;

	// --------------------------------------------------------------------------------------------------

	// --- Here the preparation phase begins

    // Read training and test data
	// (it is also possible to use ASCII format as input -> see TMVA Users Guide)
	TString fname = "./tmva_class_example.root";
	TFile *input = TFile::Open( fname );

    // --- Register the training and test trees

	TTree *signal     = (TTree*)input->Get("TreeS");
	TTree *background = (TTree*)input->Get("TreeB");
	
	std::cout << "--- TMVAClassification       : Using input file: " << input->GetName() << std::endl;

    TString tvars[] = {"myvar1 := var1+var2", "myvar2 := var1-var2", "var3", "var4", "var1", "var2"};
    std::vector<TString> variables;
    for (int i = 0; i < 6; i++) variables.push_back(tvars[i]);
	
	std::set<method_stats> rankings;
    const int max_trials = 1;
    for (int num_used = 2; num_used <= variables.size(); num_used++) {
        for (int trial = 0; trial < max_trials; trial++) {
        
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
			        "!V:Silent:Color:DrawProgressBar:AnalysisType=Classification" );

	        std::cout << std::endl;
			std::cout << "================================================" << std::endl;
			//add a random num_used sized subset of variables to train on.
            long long variable_choice = random_ksubset(variables.size(), num_used);
            for (int i = 0; i < variables.size(); i++) {
                if (variable_choice & (1LL<<i)) {
                    factory->AddVariable(variables[i], 'F');
					std::cout << "Adding variable: " << variables[i] << std::endl;
                }
            }
			std::cout << "================================================" << std::endl;            

	        // global event weights per tree (see below for setting event-wise weights)
	        Double_t signalWeight     = 1.0;
	        Double_t backgroundWeight = 1.0;

	        // You can add an arbitrary number of signal or background trees
	        factory->AddSignalTree    ( signal,     signalWeight     );
	        factory->AddBackgroundTree( background, backgroundWeight );

	        factory->SetBackgroundWeightExpression( "weight" );

	        // Apply additional cuts on the signal and background samples (can be different)
	        TCut mycuts = ""; // for example: TCut mycuts = "abs(var1)<0.5 && abs(var2-0.5)<1";
	        TCut mycutb = ""; // for example: TCut mycutb = "abs(var1)<0.5";

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
	        if (Use["BDTG"])
		        factory->BookMethod( Types::kBDT, "BDTG",
				        "!H:!V:NTrees=2000::BoostType=Grad:Shrinkage=0.1:UseBaggedBoost:BaggedSampleFraction=0.5:nCuts=20:MaxDepth=3:MaxDepth=4" );

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

	        std::cout << "==> Wrote root file: " << outputFile->GetName() << std::endl;
	        std::cout << "==> TMVAClassification is done!" << std::endl;

            //get efficiency of methods, compare with current bests
	        IMethod* i_met;
	        MethodBase* method;
	        if (Use["MLP"]) {
		        i_met = factory->GetMethod("MLP");
		        method = dyncast(i_met);
		        //std::cout << "MLP ROC Area: " << method->GetROCIntegral() <<  std::endl;
		        rankings.insert(make_method_stats(variable_choice, method->GetROCIntegral(), "MLP"));
	        }

	        if (Use["BDTG"]) {
		        i_met = factory->GetMethod("BDTG");
		        method = dyncast(i_met);
		        //std::cout << "BDTG ROC Area: " << method->GetROCIntegral() <<  std::endl;
		        rankings.insert(make_method_stats(variable_choice, method->GetROCIntegral(), "BDTG"));
	        }

	        delete factory;
	    }
	}
    
    std::set<method_stats>::iterator it;
    std::cout << "Best variables:" << std::endl;
    for (it = rankings.begin(); it != rankings.end(); it++) {
        std::cout << "================================================" << std::endl;
        std::cout << "ROC Integral: " << it->first << std::endl;
        std::cout << "Method Name: " << (it->second).second << std::endl;
        std::cout << "variables used: ";
        for (int i = 0; i < variables.size(); i++) {
            if ((it->second).first & (1LL << i)) std::cout << "[ " << variables[i] << " ] ";
        }
        cout << endl;
        std::cout << "================================================" << std::endl;
    }
	// Launch the GUI for the root macros
	//if (!gROOT->IsBatch()) TMVAGui( outfileName );
}
