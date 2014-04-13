/**********************************************************************************
 * Project   : TMVA - a Root-integrated toolkit for multivariate data analysis    *
 * Package   : TMVA                                                               *
 * Exectuable: TMVAClassificationApplication                                      *
 *                                                                                *
 * This macro provides a simple example on how to use the trained classifiers     *
 * within an analysis module                                                      *
 **********************************************************************************/

#include <cstdlib>
#include <vector>
#include <iostream>
#include <map>
#include <string>
#include <math.h>

#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TStopwatch.h"

#include "TMVA/Reader.h"
#include "TMVA/Config.h"
#include "TMVA/Tools.h"
#include "TMVA/MethodCuts.h"

using namespace std;

int main( int argc, char** argv )
{

	std::map<std::string,int> nIdenticalResults;

	std::cout << std::endl;
	std::cout << "==> Start TMVAClassificationApplication" << std::endl;

	std::string input_file;
	std::string output_file;
	if(argc>1){
		input_file = std::string(argv[1]);
		output_file = std::string(argv[2]);
	}
	else{
		input_file = "/mnt/xrootdb/alister/MVA_studies/samples/nominal_mu/tprime_650.root";
		output_file = "/home/inutard/training_code_v2/TMVA_output/nominal_mu/tprime_650.root";
	}

	TString weight_file;
	if(argc>3)
		weight_file = std::string(argv[3]);
	else
		weight_file = "weights/TMVAClassification_MLP.weights.xml";

	int nEvts = -1;
	if(argc>4)
		nEvts = atoi(argv[4]);

	cout<<"Input file: "<<input_file<<endl;
	cout<<"Output file: "<<output_file<<endl;
	cout<<"MVA Weights file: "<<weight_file<<endl;
	cout<<"Running over N events (-1 means all): "<<nEvts<<endl;

	// --------------------------------------------------------------------------------------------------

	// --- Create the Reader object

	TMVA::Reader *reader = new TMVA::Reader( "!Color:!Silent" );

	// Create a set of variables and declare them to the reader
	// - the variable names must corresponds in name and type to
	// those given in the weight file(s) that you use


	//Float_t jet_pt[50]; // AL: this was wrong: cannot initialise an array with a non-initialised number

	Float_t	Ht;
	Float_t	b_lep_pt;
	Float_t	b_had_pt;
	Float_t	dR_lnu;
	Float_t	dR_Whad_blep;
	Float_t	dR_Whad_bhad;
	Float_t	dR_lb_lep;
	Float_t	dR_lb_had;

	reader->AddVariable( "Ht", &Ht );
	reader->AddVariable( "b_had_pt", &b_had_pt );
	reader->AddVariable( "dR_lb_lep", &dR_lb_lep );
	reader->AddVariable( "b_lep_pt", &b_lep_pt );
	reader->AddVariable( "dR_Whad_bhad", &dR_Whad_bhad );
	reader->AddVariable( "dR_lnu", &dR_lnu );
	reader->AddVariable( "dR_Whad_blep", &dR_Whad_blep );
	reader->AddVariable( "dR_lb_had", &dR_lb_had );

	// Spectator variables declared in the training have to be added to the reader, too
	Float_t mc_weight;
	Float_t m_reco;
	Int_t runNumber, eventNumber, analysis_channel;
	reader->AddSpectator( "mc_weight := weight_1btin_70*weight_70/btweight_70", &mc_weight ); // SN hard code
	reader->AddSpectator( "m_reco", &m_reco);
	reader->AddSpectator( "runNumber", &runNumber);
	reader->AddSpectator( "eventNumber", &eventNumber);

	// --- Book the MVA methods
	TString methodName = TString("MLP") + TString(" method");
	reader->BookMVA( methodName, weight_file ); 


	// Prepare input tree (this must be replaced by your data source)
	// in this example, there is a toy tree with signal and one with background events
	// we'll later on use only the "signal" events for the test in this example.
	//
	TFile* input = TFile::Open( (TString)input_file, "READ" ); // check if file in local directory exists

	if (!input) {
		std::cout << "ERROR: could not open data file" << std::endl;
		exit(1);
	}
	std::cout << "--- TMVAClassificationApp    : Using input file: " << input->GetName() << std::endl;

	// --- Event loop

	// Prepare the event tree
	// - here the variable names have to corresponds to your tree
	// - you can use the same variables as above which is slightly faster,
	//   but of course you can use different ones and copy the values inside the event loop
	//
	std::cout << "--- Select signal sample" << std::endl;
	TTree* theTree = (TTree*)input->Get("mini");

	Float_t weight_70, btweight_70, weight_1btin_70;
	theTree->SetBranchAddress( "Ht", &Ht );
	theTree->SetBranchAddress( "b_lep_pt", &b_lep_pt );
	theTree->SetBranchAddress( "b_had_pt", &b_had_pt );
	theTree->SetBranchAddress( "dR_lnu", &dR_lnu );
	theTree->SetBranchAddress( "dR_Whad_blep", &dR_Whad_blep );
	theTree->SetBranchAddress( "dR_Whad_bhad", &dR_Whad_bhad );
	theTree->SetBranchAddress( "dR_lb_lep", &dR_lb_lep );
	theTree->SetBranchAddress( "dR_lb_had", &dR_lb_had );   
	theTree->SetBranchAddress( "analysis_channel", &analysis_channel);
	theTree->SetBranchAddress( "weight_70", &weight_70);
	theTree->SetBranchAddress( "btweight_70", &btweight_70);
	theTree->SetBranchAddress( "weight_1btin_70", &weight_1btin_70);
	// Efficiency calculator for cut method
	//Int_t    nSelCutsGA = 0;
	//Double_t effS       = 0.7;

	std::vector<Float_t> vecVar(5); // vector for EvaluateMVA tests

	std::cout << "--- Processing: " << theTree->GetEntries() << " events" << std::endl;
	TStopwatch sw;
	sw.Start();


	//TString target_name = output_file;
	//TFile *target = 0;
	//target = TFile::Open( target_name,"RECREATE" );
	TTree *tree = theTree->CloneTree(0);
	//tree->SetDirectory(target);
	//target->cd();
	Float_t MVA;
	tree->Branch("MVA", &MVA, "MVA/F");


	Int_t nEvent = theTree->GetEntries();
	if(nEvts != -1) nEvent = nEvts;

	int cnt = 0;
	double tot = 0;

	//UInt_t nbin = 100;
	//TH1F *histMLP(0);
	//histMLP = new TH1F( "MVA_MLP",           "MVA_MLP",           nbin, -1.25, 1.5 );
	for (Long64_t ievt=0; ievt<nEvent; ievt++) {


		if (ievt%10000 == 0) {
			std::cout << "--- ... Processing event: " << ievt << std::endl;
			std::cout << tot << endl; 
		}

		theTree->GetEntry(ievt);
		if (analysis_channel == 0) continue;

		//theTree->Show(ievt);

		// --- Return the MVA outputs and fill into histograms
		MVA = reader->EvaluateMVA( "MLP method" );

		if (MVA <= 0.9965) continue;

		cnt++;
		//histMLP->Fill(MVA);

		tot += weight_1btin_70*weight_70/btweight_70;//MVA; //mc_weight;
		/*
		   if (tot > 10.3) {
		   cout << cnt << " " << tot << endl;
		   break;
		   }
		 */
		/*
		   if (cnt >= 4737) {
		   cout << tot << endl;    
		   break; 
		   }
		 */
		// tree->Fill();
	}

	cout << "Final count: " << tot << endl;
	// Get elapsed time
	sw.Stop();
	std::cout << "--- End of event loop: "; sw.Print();

	std::cout << "--- Created root file containing the MVA output histograms" << std::endl;

	//histMLP->Write();
	//target->Write();
	//target->Close();

	delete reader;

	std::cout << "==> TMVAClassificationApplication is done!" << std::endl << std::endl;
}
