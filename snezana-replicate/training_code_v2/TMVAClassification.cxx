/* ////////////////////////////////
Usage:
./TMVAClassification MLP (Cuts) -p /atlas/data1/userdata/snezana/TprimeAnalysis/2011_7TeV_full/TRCR-11-00-00-05/VLQTree/
el resolved loose _1btin
/////////////////////////////////*/
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <map>
#include <string>
#include <cmath>
#include <string>

#include "TChain.h"
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TObjString.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TH1.h"
#include "TLegend.h"
#include "TText.h"
#include "TH2.h"
#include "TMath.h"
#include "TGraph.h"
#include "THStack.h"
#include <TRandom3.h>
#include "tmvaglob.C"

#include "TMVA/Factory.h"
#include "TMVA/DataSetManager.h"
#include "TMVA/DataSetFactory.h"
#include "TMVA/DataInputHandler.h"
#include "TMVA/Tools.h"

using namespace std;

// read input data file with ascii format (otherwise ROOT) ?
Bool_t ReadDataFromAsciiIFormat = kFALSE;

enum leptonChannel {el, mu, elmu};
enum analysisChannel {rest, boosted, resolved, combined};
enum selection {loose, tight};
enum b_tag {_0btin, _1btin, _2btin};
enum HistType { MVAType = 0, ProbaType = 1, RarityType = 2, CompareType = 3 };
struct flags {
  TString path;
  leptonChannel lep;
  analysisChannel analysis_channel;
  b_tag btag;
  TString training_type;
  TString output_path;
};

flags getflags(int argc, char *argv[]);
TCut getBaseCut(flags f);
TString getPath(flags f);
TString getWeight(flags f);
TH1* move_overflow(TH1 *h);

//////////////////////////////////////////////////////////////////

int main( int argc, char** argv )
{

   //---------------------------------------------------------------
   // Default MVA methods to be trained + tested
   std::map<std::string,int> Use;
   //
   // --- Neural Networks (all are feed-forward Multilayer Perceptrons)
   Use["MLP"]             = 1; // Recommended ANN
   Use["MLPBFGS"]         = 0; // Recommended ANN with optional training method
   Use["MLPBNN"]          = 0; // Recommended ANN with BFGS training method and bayesian regulator   

   // --- Cut optimisation
   Use["Cuts"]            = 1;
   Use["CutsD"]           = 0;
   Use["CutsPCA"]         = 0;
   Use["CutsGA"]          = 0;
   Use["CutsSA"]          = 0;

   // --- Boosted Decision Trees
   Use["BDT"]             = 1; // uses Adaptive Boost
   Use["BDTG"]            = 0; // uses Gradient Boost
   Use["BDTB"]            = 0; // uses Bagging
   Use["BDTD"]            = 0; // decorrelation + Adaptive Boost
   Use["BDTF"]            = 0; // allow usage of fisher discriminant for node splitting 

   // ---------------------------------------------------------------

   std::cout << std::endl << "==> Start TMVAClassification" << std::endl;

   bool batchMode(false);
   bool useDefaultMethods(true);

   flags fl = getflags(argc,argv);
   // Select methods (don't look at this code - not of interest)
   int pos = 1;
   if (argc>1) {
     string entry = string(argv[pos]);
     while ((pos<argc) && (entry != "-p")) {
       pos++;
       entry = string(argv[pos]);
     }
   }
   cout << "pos = " << pos << endl;

   for (int i=1; i<pos; i++) {
     std::string regMethod(argv[i]);
     if(regMethod=="-b" || regMethod=="--batch") {
       batchMode=true;
       continue;
     }
     if (Use.find(regMethod) == Use.end()) {
       std::cout << "Method \"" << regMethod << "\" not known in TMVA under this name. Choose among the following:" << std::endl;
       for (std::map<std::string,int>::iterator it = Use.begin(); it != Use.end(); it++) std::cout << it->first << " ";
       std::cout << std::endl;
       return 1;
     }
     useDefaultMethods = false;
   }
   
   if (!useDefaultMethods) {
     for (std::map<std::string,int>::iterator it = Use.begin(); it != Use.end(); it++) it->second = 0;
     for (int i=1; i<argc; i++) {
       std::string regMethod(argv[i]);
       if(regMethod=="-b" || regMethod=="--batch") continue;
       Use[regMethod] = 1;
     }
   }
   
   // --------------------------------------------------------------------------------------------------

   // --- Here the preparation phase begins

   // Create a ROOT output file where TMVA will store ntuples, histograms, etc.
   TString outfileName( "TMVA.root" );
   TFile* outputFile = TFile::Open( outfileName, "RECREATE" );

   TMVA::Factory *factory = new TMVA::Factory( "TMVAClassification", outputFile,
                                               "!V:!Silent:Color:DrawProgressBar:Transformations=I;D;P;G,D:AnalysisType=Classification" );
   
   // Define the signal + background tree
   TString path = getPath(fl);
   TChain* tree = new TChain("mini");
   tree->Add(path + "/nominal_el/tprime_650_1M.root");
   tree->Add(path + "/nominal_mu/tprime_650_1M.root");
   tree->Add(path + "/nominal_el/ttbar.root");
   tree->Add(path + "/nominal_mu/ttbar.root");
   tree->Add(path + "/nominal_el/wjets.root");
   tree->Add(path + "/nominal_mu/wjets.root");
   tree->Add(path + "/nominal_el/zjets.root");
   tree->Add(path + "/nominal_mu/zjets.root");
   tree->Add(path + "/nominal_el/singletop.root");
   tree->Add(path + "/nominal_mu/singletop.root");
   tree->Add(path + "/nominal_el/diboson.root");
   tree->Add(path + "/nominal_mu/diboson.root");

   Float_t	Ht;
   Float_t	b_lep_pt;
   Float_t	b_had_pt;
   Float_t	dR_lnu;
   Float_t	dR_Whad_blep;
   Float_t	dR_Whad_bhad;
   Float_t	dR_lb_lep;
   Float_t	dR_lb_had;
   Int_t 	runNumber;
   Int_t 	eventNumber;

   tree->SetBranchAddress( "Ht", &Ht );
   tree->SetBranchAddress( "b_had_pt", &b_had_pt );
   tree->SetBranchAddress( "dR_lb_lep", &dR_lb_lep );
   tree->SetBranchAddress( "b_lep_pt", &b_lep_pt );
   tree->SetBranchAddress( "dR_Whad_bhad", &dR_Whad_bhad );
   tree->SetBranchAddress( "dR_lnu", &dR_lnu );
   tree->SetBranchAddress( "dR_Whad_blep", &dR_Whad_blep );
   tree->SetBranchAddress( "dR_lb_had", &dR_lb_had );
   tree->SetBranchAddress( "runNumber", &runNumber );
   tree->SetBranchAddress( "eventNumber", &eventNumber );

   // Define the input variables that shall be used for the MVA training
   fstream input_file;
   input_file.open("input_variables.txt");
   string line;
   
   while (getline(input_file, line)) {
     TString var_info[3];
     int itr = 0;
     size_t found;
     found = line.find_first_of(" ");
     while (found!=string::npos) {
       var_info[itr] = (TString)line.substr(0, found);
       line = line.substr(found+1, line.length());
       found = line.find_first_of(" ");
       itr++;
     }
     string var_type = line;
     if (var_info[2] == "nan") var_info[2]="";
     if (var_type == "F")
       factory->AddVariable( var_info[0], var_info[1], var_info[2], 'F');
     else if (var_type == "I")
       factory->AddVariable( var_info[0], var_info[1], var_info[2], 'I');
     else {
       cout << "Unknown variable type" << endl;
       exit(-1);
     }
   }

   // Add so-called "Spectator variables"
   TString mcweight = "mc_weight := " + getWeight(fl);
   factory->AddSpectator( mcweight, "mc_weight", "MeV", 'F' );
   factory->AddSpectator( "m_reco", "m_reco", "MeV", 'F' );
   factory->AddSpectator( "runNumber", "runNumber", "", 'I' );
   factory->AddSpectator( "eventNumber", "eventNumber", "", 'I' );

   // Register the trees to the factory
   // At this step the cuts for splitting trees are defined
   // if channelNumber == 119688 -> signal, else background

   // The splittinf criterion
   TCut train_a = "(int(abs(lep_phi)*125434413.)%2==0)";
   TCut test_a = "(int(abs(lep_phi)*125434413.)%2!=0)";

   TCut cutSigTrain;
   TCut cutSigTest;
   TCut cutBgdTrain;
   TCut cutBgdTest;

   TString nr = "119688";
   TString cut_ss = "(channelNumber == "+nr+")";
   TCut cut_s = (TCut)cut_s;

   if (fl.training_type == "a") {
     // Training A:
     cutSigTrain = getBaseCut(fl) + "(channelNumber == 119688)" + train_a;
     cutSigTest = getBaseCut(fl) + "(channelNumber == 119688)" + test_a;
     cutBgdTrain = getBaseCut(fl) + "(channelNumber != 119688)" + train_a;
     cutBgdTest = getBaseCut(fl) + "(channelNumber != 119688)" + test_a;
   }
   else {
     // Training B:
     cutSigTrain = getBaseCut(fl) + "(channelNumber == 119688)" + test_a;
     cutSigTest = getBaseCut(fl) + "(channelNumber == 119688)" + train_a;
     cutBgdTrain = getBaseCut(fl) + "(channelNumber != 119688)" + test_a;
     cutBgdTest = getBaseCut(fl) + "(channelNumber != 119688)" + train_a;
   }

   factory->AddTree (tree, "Signal", 1.0, cutSigTrain, "train");
   factory->AddTree (tree, "Signal", 1.0, cutSigTest, "test");
   factory->AddTree (tree, "Background", 1.0, cutBgdTrain, "train");
   factory->AddTree (tree, "Background", 1.0, cutBgdTest, "test");
   
   // --- end of tree registration 

   // Set individual event weights (the variables must exist in the original TTree)
   factory->SetSignalWeightExpression    (getWeight(fl));
   factory->SetBackgroundWeightExpression(getWeight(fl));

   ///////////////////////////////////////////////////////////////////////////////////////
   // ---- Book MVA methods
   // Cut optimisation
   if (Use["Cuts"])
      factory->BookMethod( TMVA::Types::kCuts, "Cuts",
                           "!H:!V:FitMethod=MC:EffSel:SampleSize=200000:VarProp=FSmart" );

   if (Use["CutsD"])
      factory->BookMethod( TMVA::Types::kCuts, "CutsD",
                           "!H:!V:FitMethod=MC:EffSel:SampleSize=200000:VarProp=FSmart:VarTransform=Decorrelate" );

   if (Use["CutsPCA"])
      factory->BookMethod( TMVA::Types::kCuts, "CutsPCA",
                            "!H:!V:FitMethod=MC:EffSel:SampleSize=200000:VarProp=FSmart:VarTransform=PCA" );

   if (Use["CutsGA"])
      factory->BookMethod( TMVA::Types::kCuts, "CutsGA",
                           "H:!V:FitMethod=GA:CutRangeMin[0]=-10:CutRangeMax[0]=10:VarProp[1]=FMax:EffSel:Steps=30:Cycles=3:PopSize=400:SC_steps=10:SC_rate=5:SC_factor=0.95" );

   if (Use["CutsSA"])
      factory->BookMethod( TMVA::Types::kCuts, "CutsSA",
                           "!H:!V:FitMethod=SA:EffSel:MaxCalls=150000:KernelTemp=IncAdaptive:InitialTemp=1e+6:MinTemp=1e-6:Eps=1e-10:UseDefaultScale" );

   // TMVA ANN: MLP (recommended ANN) -- all ANNs in TMVA are Multilayer Perceptrons
   if (Use["MLP"])
     factory->BookMethod( TMVA::Types::kMLP, "MLP", "H:!V:NeuronType=tanh:VarTransform=N:NCycles=700:HiddenLayers=N+20:TestRate=5:!UseRegulator:SamplingTraining=False:EstimatorType=CE" );

   if (Use["MLPBFGS"])
      factory->BookMethod( TMVA::Types::kMLP, "MLPBFGS", "H:!V:NeuronType=tanh:VarTransform=N:NCycles=650:HiddenLayers=N+1:TestRate=5:TrainingMethod=BFGS:!UseRegulator:SamplingTraining=False:EstimatorType=CE" );

   if (Use["MLPBNN"])
      factory->BookMethod( TMVA::Types::kMLP, "MLPBNN", "H:!V:NeuronType=tanh:VarTransform=N:NCycles=650:HiddenLayers=N+1:TestRate=5:TrainingMethod=BFGS:UseRegulator:SamplingTraining=False:EstimatorType=CE" ); // BFGS training with bayesian regulators

   // Boosted Decision Trees
   if (Use["BDTG"]) // Gradient Boost
      factory->BookMethod( TMVA::Types::kBDT, "BDTG",
                           "!H:!V:NTrees=1000:BoostType=Grad:Shrinkage=0.10:UseBaggedGrad:GradBaggingFraction=0.5:nCuts=20:NNodesMax=5" );

   if (Use["BDT"])  // Adaptive Boost
      factory->BookMethod( TMVA::Types::kBDT, "BDT",
                           "!H:!V:NTrees=850:nEventsMin=150:MaxDepth=3:BoostType=AdaBoost:AdaBoostBeta=0.5:SeparationType=GiniIndex:nCuts=20:PruneMethod=NoPruning" );

   if (Use["BDTB"]) // Bagging
      factory->BookMethod( TMVA::Types::kBDT, "BDTB",
                           "!H:!V:NTrees=400:BoostType=Bagging:SeparationType=GiniIndex:nCuts=20:PruneMethod=NoPruning" );

   if (Use["BDTD"]) // Decorrelation + Adaptive Boost
      factory->BookMethod( TMVA::Types::kBDT, "BDTD",
                           "!H:!V:NTrees=400:nEventsMin=400:MaxDepth=3:BoostType=AdaBoost:SeparationType=GiniIndex:nCuts=20:PruneMethod=NoPruning:VarTransform=Decorrelate" );

   if (Use["BDTF"])  // Allow Using Fisher discriminant in node splitting for (strong) linearly correlated variables
      factory->BookMethod( TMVA::Types::kBDT, "BDTMitFisher",
                           "!H:!V:NTrees=50:nEventsMin=150:UseFisherCuts:MaxDepth=3:BoostType=AdaBoost:AdaBoostBeta=0.5:SeparationType=GiniIndex:nCuts=20:PruneMethod=NoPruning" );
   
   // --------------------------------------------------------------------------------------------------

   // ---- Now you can optimize the setting (configuration) of the MVAs using the set of training events

   // factory->OptimizeAllMethods("SigEffAt001","Scan");
   // factory->OptimizeAllMethods("ROCIntegral","GA");

   // --------------------------------------------------------------------------------------------------

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

   std::cout << "==> Wrote root file: " << outputFile->GetName() << std::endl
             << "==> TMVAClassification is done!" << std::endl
             << std::endl
             << "==> To view the results, launch the GUI: \"root -l ./TMVAGui.C\"" << std::endl
             << std::endl;

   // Clean up
   delete factory;
   cout << "Output path: " << fl.output_path << endl;
}

//////////////////////////////////////////////////////////////////

flags getflags(int argc, char *argv[]) {
  flags f;

  int pos = 1;
  if (argc>1) {
    string entry = string(argv[pos]);
    while ((pos<argc) && (entry != "-p")) {
      pos++;
      entry = string(argv[pos]);
    }
  }

  if (argc > pos+1) {
    f.path = TString(argv[pos+1]);
  }

  if(argc > pos+2) {
    string lep_channel = string(argv[pos+2]);
    if (lep_channel.find("el") != string::npos)
      f.lep = el;
    else if (lep_channel.find("mu") != string::npos)
      f.lep = mu;
    else if (lep_channel.find("elmu") != string::npos)
      f.lep = elmu;
    else {
      cout << "Invalid lepton channel!" << endl;
      exit(-1);
    }
  }  

  if(argc > pos+3) {
    string an_channel = string(argv[pos+3]);
    if (an_channel.find("rest") != string::npos)
      f.analysis_channel = rest;
    else if (an_channel.find("boosted") != string::npos)
      f.analysis_channel = boosted;
    else if (an_channel.find("resolved") != string::npos)
      f.analysis_channel = resolved;
    else if (an_channel.find("combined") != string::npos)
      f.analysis_channel = combined;
    else {
      cout << "Invalid analysis channel!" << endl;
      exit(-1);
    }
  }

  if(argc > pos+4) {
    string bt =  string(argv[pos+4]);
    if (bt.find("_0btin") != string::npos)
      f.btag = _0btin;
    else if (bt.find("_1btin") != string::npos)
      f.btag = _1btin;
    else if (bt.find("_2btin") != string::npos)
      f.btag = _2btin;
    else {
      cout << "Invalid btag configutation!" << endl;
      exit(-1);
    }
  }

  if(argc > pos+5) {
    string training_type = string(argv[pos+5]);
    f.training_type = training_type;
  }

  if(argc > pos+6) {
    string output_path = string(argv[pos+6]);
    f.output_path = output_path;
  }

  return f;
}

TCut getBaseCut(flags f) {
  TCut cut = "weight_70>0";

  if (f.analysis_channel == rest) cut += "analysis_channel==0";
  else if (f.analysis_channel == boosted) cut += "analysis_channel==1";
  else if (f.analysis_channel == resolved) cut += "analysis_channel==2";
  else if (f.analysis_channel == combined) cut += "analysis_channel!=0";
  else {
    cout << "Invalid analysis channel!" << endl;
    exit(-1);
  }

  cout << "Base cut: " << cut << endl;
  return cut;
 }

TString getPath(flags f) {
  TString path = f.path;
  return path;  
}

TString getWeight(flags f) {
  if (f.btag == _0btin) return "weight_70";
  else if (f.btag == _1btin) return "weight_1btin_70*weight_70/btweight_70";
  else if (f.btag == _2btin) return "weight_2btin_70*weight_70/btweight_70";
}

TH1* move_overflow(TH1 *h)
{
  h->SetBinContent(1, h->GetBinContent(1)+h->GetBinContent(0));
  h->SetBinContent(0, 0);
  h->SetBinContent(h->GetXaxis()->GetNbins(), 
		   h->GetBinContent(h->GetXaxis()->GetNbins())+
		   h->GetBinContent(h->GetXaxis()->GetNbins()+1));
  h->SetBinContent(h->GetXaxis()->GetNbins()+1, 0);
  
  return h;
}
