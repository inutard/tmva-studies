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

#include "tmvaglob.C"

#include "TMVA/Factory.h"
#include "TMVA/DataSetManager.h"
#include "TMVA/DataSetFactory.h"
#include "TMVA/DataInputHandler.h"
//#include "TMVA/EventVectorOfClassesOfTreeType.h"
//#include "TMVA/EvtStatsPerClass.h"
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
  //selection selec;
  b_tag btag;
  TString output_path;
};

// SN hack
Float_t n_tprime = 20.2394; // combined elmu
Float_t n_bgd = 2073.49; // all MC backgrounds combined elmu
Float_t n_tprime_cut = 8.1; // acceptance of the t' 650 sample in my implementation of the tight selection

flags getflags(int argc, char *argv[]);
TCut getBaseCut(flags f);
TString getPath(flags f);
TString getWeight(flags f);
Float_t getNBgd_for_NSig(Float_t nsig, TDirectory* mytitDir );
TH1* move_overflow(TH1 *h);
void mvas( TString fin, HistType htype, Bool_t useTMVAStyle, TString output_path );
void annconvergencetest( TString fin, Bool_t useTMVAStyle, TString output_path  );
void printout_train_test(Float_t nsig);
void PrintOut();
Float_t find_MLP_cut(Float_t nsig, TFile* f );
void find_max_SoverB(TString output_path);


//////////////////////////////////////////////////////////////////

int main( int argc, char** argv )
{

   //---------------------------------------------------------------
   // Default MVA methods to be trained + tested
   std::map<std::string,int> Use;
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
       //cout << line.substr(0, found) << endl;
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
       cout << "Unknown variable type (TMVAClassification.cxx:170)" << endl;
       exit(-1);
     }
   }


   // You can add so-called "Spectator variables", which are not used in the MVA training,
   // but will appear in the final "TestTree" produced by TMVA. This TestTree will contain the
   // input variables, the response values of all trained MVAs, and the spectator variables
   //factory->AddSpectator( "jet_pt[2]", "jet2 pt", "MeV", 'F' );
   TString mcweight = "mc_weight := " + getWeight(fl);
   factory->AddSpectator( mcweight, "mc_weight", "MeV", 'F' );
   factory->AddSpectator( "m_reco", "m_reco", "MeV", 'F' );
   factory->AddSpectator( "runNumber", "runNumber", "", 'I' );

   // Read training and test data   
   TFile *input_sig_el = TFile::Open(getPath(fl) + "/el/nominal/tprime_650_1M.root");
   TFile *input_sig_mu = TFile::Open(getPath(fl) + "/mu/nominal/tprime_650_1M.root");
   TFile *input_ttbar_el = TFile::Open(getPath(fl) + "/el/nominal/ttbar.root");
   TFile *input_ttbar_mu = TFile::Open(getPath(fl) + "/mu/nominal/ttbar.root");
   TFile *input_wjets_el = TFile::Open(getPath(fl) + "/el/nominal/wjets.root");
   TFile *input_wjets_mu = TFile::Open(getPath(fl) + "/mu/nominal/wjets.root");
   TFile *input_zjets_el = TFile::Open(getPath(fl) + "/el/nominal/zjets.root");
   TFile *input_zjets_mu = TFile::Open(getPath(fl) + "/mu/nominal/zjets.root");
   TFile *input_singletop_el = TFile::Open(getPath(fl) + "/el/nominal/singletop.root");
   TFile *input_singletop_mu = TFile::Open(getPath(fl) + "/mu/nominal/singletop.root");
   TFile *input_diboson_el = TFile::Open(getPath(fl) + "/el/nominal/diboson.root");
   TFile *input_diboson_mu = TFile::Open(getPath(fl) + "/mu/nominal/diboson.root");

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
   
   // To give different trees for training and testing, do as follows:
   //    factory->AddSignalTree( signalTrainingTree, signalTrainWeight, "Training" );
   //    factory->AddSignalTree( signalTestTree,     signalTestWeight,  "Test" );
   //
   // --- end of tree registration 

   // Set individual event weights (the variables must exist in the original TTree)
   factory->SetSignalWeightExpression    (getWeight(fl));
   factory->SetBackgroundWeightExpression(getWeight(fl));

   // Apply additional cuts on the signal and background samples (can be different)
   TCut mycuts = getBaseCut(fl);
   TCut mycutb = getBaseCut(fl);

   // Tell the factory how to use the training and testing events
   factory->PrepareTrainingAndTestTree( mycuts, mycutb,
					"nTrain_Signal=0:nTrain_Background=0:SplitMode=Random:NormMode=NumEvents:!V" );

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
     factory->BookMethod( TMVA::Types::kMLP, "MLP", "H:!V:NeuronType=tanh:VarTransform=N:NCycles=650:HiddenLayers=N+1:TestRate=5:!UseRegulator:SamplingTraining=False:EstimatorType=CE" );

   if (Use["MLPBFGS"])
      factory->BookMethod( TMVA::Types::kMLP, "MLPBFGS", "H:!V:NeuronType=tanh:VarTransform=N:NCycles=600:HiddenLayers=N+10:TestRate=5:TrainingMethod=BFGS:!UseRegulator:SamplingTraining=False:EstimatorType=CE" );

   if (Use["MLPBNN"])
      factory->BookMethod( TMVA::Types::kMLP, "MLPBNN", "H:!V:NeuronType=tanh:VarTransform=N:NCycles=600:HiddenLayers=N+5:TestRate=5:TrainingMethod=BFGS:UseRegulator" ); // BFGS training with bayesian regulators

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
   mvas("TMVA.root", CompareType, kTRUE, fl.output_path);
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
      //cout << "getflags: " << pos << " " << entry << endl;
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
    string output_path = string(argv[pos+5]);
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
  //if (f.selec == loose) cut += "loose==1";
  //else if (f.selec == loose) cut += "tight==1";

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

void mvas( TString fin, HistType htype, Bool_t useTMVAStyle, TString output_path )
{
   // set style and remove existing canvas'
   TMVAGlob::Initialize( useTMVAStyle );

   // switches
   const Bool_t Save_Images = kTRUE;

   // checks if file with name "fin" is already open, and if("MethodMLP/MLP") not opens one
   TFile* file = TMVAGlob::OpenFile( fin );  

   // define Canvas layout here!
   Int_t xPad = 1; // no of plots in x
   Int_t yPad = 1; // no of plots in y
   Int_t noPad = xPad * yPad ; 
   const Int_t width = 600;   // size of canvas

   // this defines how many canvases we need
   TCanvas *c = 0;

   // counter variables
   Int_t countCanvas = 0;

   // search for the right histograms in full list of keys
   TIter next(file->GetListOfKeys());
   TKey *key(0);   
   while ((key = (TKey*)next())) {

      if (!TString(key->GetName()).BeginsWith("Method_")) continue;
      if (!gROOT->GetClass(key->GetClassName())->InheritsFrom("TDirectory")) continue;

      TString methodName;
      TMVAGlob::GetMethodName(methodName,key);

      TDirectory* mDir = (TDirectory*)key->ReadObj();

      TIter keyIt(mDir->GetListOfKeys());
      TKey *titkey;
      while ((titkey = (TKey*)keyIt())) {

         if (!gROOT->GetClass(titkey->GetClassName())->InheritsFrom("TDirectory")) continue;

         TDirectory *titDir = (TDirectory *)titkey->ReadObj();
         TString methodTitle;
         TMVAGlob::GetMethodTitle(methodTitle,titDir);

         cout << "--- Found directory for method: " << methodName << "::" << methodTitle << flush;
         TString hname = "MVA_" + methodTitle;
         if      (htype == ProbaType  ) hname += "_Proba";
         else if (htype == RarityType ) hname += "_Rarity";
         TH1* sig = dynamic_cast<TH1*>(titDir->Get( hname + "_S" ));
         TH1* bgd = dynamic_cast<TH1*>(titDir->Get( hname + "_B" ));

         if (sig==0 || bgd==0) {
            if     (htype == MVAType)     
               cout << ":\t mva distribution not available (this is normal for Cut classifier)" << endl;
            else if(htype == ProbaType)   
               cout << ":\t probability distribution not available" << endl;
            else if(htype == RarityType)  
               cout << ":\t rarity distribution not available" << endl;
            else if(htype == CompareType) 
               cout << ":\t overtraining check not available" << endl;
            else cout << endl;
            continue;
         }

         cout << " containing " << hname << "_S/_B" << endl;
         // chop off useless stuff
         sig->SetTitle( Form("TMVA response for classifier: %s", methodTitle.Data()) );
         if      (htype == ProbaType) 
            sig->SetTitle( Form("TMVA probability for classifier: %s", methodTitle.Data()) );
         else if (htype == RarityType) 
            sig->SetTitle( Form("TMVA Rarity for classifier: %s", methodTitle.Data()) );
         else if (htype == CompareType) 
            sig->SetTitle( Form("TMVA overtraining check for classifier: %s", methodTitle.Data()) );
         
         // create new canvas
         TString ctitle = ((htype == MVAType) ? 
                           Form("TMVA response %s",methodTitle.Data()) : 
                           (htype == ProbaType) ? 
                           Form("TMVA probability %s",methodTitle.Data()) :
                           (htype == CompareType) ? 
                           Form("TMVA comparison %s",methodTitle.Data()) :
                           Form("TMVA Rarity %s",methodTitle.Data()));
         
         c = new TCanvas( Form("canvas%d", countCanvas+1), ctitle, 
                          countCanvas*50+200, countCanvas*20, width, (Int_t)width*0.78 ); 
	 c->SetLogy();
    
         // set the histogram style
         TMVAGlob::SetSignalAndBackgroundStyle( sig, bgd );
         
         // normalise both signal and background
         TMVAGlob::NormalizeHists( sig, bgd );
         
         // frame limits (choose judicuous x range)
         Float_t nrms = 10;
         cout << "--- Mean and RMS (S): " << sig->GetMean() << ", " << sig->GetRMS() << endl;
         cout << "--- Mean and RMS (B): " << bgd->GetMean() << ", " << bgd->GetRMS() << endl;
         Float_t xmin = TMath::Max( TMath::Min(sig->GetMean() - nrms*sig->GetRMS(), 
                                               bgd->GetMean() - nrms*bgd->GetRMS() ),
                                    sig->GetXaxis()->GetXmin() );
         Float_t xmax = TMath::Min( TMath::Max(sig->GetMean() + nrms*sig->GetRMS(), 
                                               bgd->GetMean() + nrms*bgd->GetRMS() ),
                                    sig->GetXaxis()->GetXmax() );
	 xmin = -0.2;
	 xmax = 1.2;
         Float_t ymin = 0;
         Float_t maxMult = (htype == CompareType) ? 1.3 : 1.2;
         Float_t ymax = TMath::Max( sig->GetMaximum(), bgd->GetMaximum() )*maxMult;
   
         // build a frame
         Int_t nb = 500;
         TString hFrameName(TString("frame") + methodTitle);
         TObject *o = gROOT->FindObject(hFrameName);
         if(o) delete o;
         TH2F* frame = new TH2F( hFrameName, sig->GetTitle(), 
                                 nb, xmin, xmax, nb, 1e-4, ymax );
         //TH2F* frame = new TH2F( hFrameName, sig->GetTitle(), 
         //                        nb, xmin, xmax, nb, ymin, ymax );
	 //frame->SetMinimum(1e-4);
         frame->GetXaxis()->SetTitle( methodTitle + ((htype == MVAType || htype == CompareType) ? " response" : "") );
         if      (htype == ProbaType  ) frame->GetXaxis()->SetTitle( "Signal probability" );
         else if (htype == RarityType ) frame->GetXaxis()->SetTitle( "Signal rarity" );
         frame->GetYaxis()->SetTitle("(1/N) dN^{ }/^{ }dx");
         TMVAGlob::SetFrameStyle( frame );
   
         // eventually: draw the frame
	 frame->SetMinimum(1e-5);
         frame->Draw();  
    
         c->GetPad(0)->SetLeftMargin( 0.105 );
         frame->GetYaxis()->SetTitleOffset( 1.2 );

         // Draw legend               
         TLegend *legend= new TLegend( c->GetLeftMargin(), 1 - c->GetTopMargin() - 0.12, 
                                       c->GetLeftMargin() + (htype == CompareType ? 0.40 : 0.3), 1 - c->GetTopMargin() );
         legend->SetFillStyle( 1 );
         legend->AddEntry(sig,TString("Signal")     + ((htype == CompareType) ? " (test sample)" : ""), "F");
         legend->AddEntry(bgd,TString("Background") + ((htype == CompareType) ? " (test sample)" : ""), "F");
         legend->SetBorderSize(1);
         legend->SetMargin( (htype == CompareType ? 0.2 : 0.3) );
         legend->Draw("same");

         sig->Draw("samehist");
         bgd->Draw("samehist");
   
         if (htype == CompareType) {
            // if overtraining check, load additional histograms
            TH1* sigOv = 0;
            TH1* bgdOv = 0;

            TString ovname = hname += "_Train";
            sigOv = dynamic_cast<TH1*>(titDir->Get( ovname + "_S" ));
            bgdOv = dynamic_cast<TH1*>(titDir->Get( ovname + "_B" ));
      
            if (sigOv == 0 || bgdOv == 0) {
               cout << "+++ Problem in \"mvas.C\": overtraining check histograms do not exist" << endl;
            }
            else {
               cout << "--- Found comparison histograms for overtraining check" << endl;

               TLegend *legend2= new TLegend( 1 - c->GetRightMargin() - 0.42, 1 - c->GetTopMargin() - 0.12,
                                              1 - c->GetRightMargin(), 1 - c->GetTopMargin() );
               legend2->SetFillStyle( 1 );
               legend2->SetBorderSize(1);
               legend2->AddEntry(sigOv,"Signal (training sample)","P");
               legend2->AddEntry(bgdOv,"Background (training sample)","P");
               legend2->SetMargin( 0.1 );
               legend2->Draw("same");
            }
            // normalise both signal and background
            TMVAGlob::NormalizeHists( sigOv, bgdOv );

            Int_t col = sig->GetLineColor();
            sigOv->SetMarkerColor( col );
            sigOv->SetMarkerSize( 0.7 );
            sigOv->SetMarkerStyle( 20 );
            sigOv->SetLineWidth( 1 );
            sigOv->SetLineColor( col );
	    sigOv = move_overflow(sigOv);
            sigOv->Draw("e1same");
      
            col = bgd->GetLineColor();
            bgdOv->SetMarkerColor( col );
            bgdOv->SetMarkerSize( 0.7 );
            bgdOv->SetMarkerStyle( 20 );
            bgdOv->SetLineWidth( 1 );
            bgdOv->SetLineColor( col );
	    bgdOv = move_overflow(bgdOv);
            bgdOv->Draw("e1same");

            ymax = TMath::Max( (Double_t) ymax, (Double_t) TMath::Max( sigOv->GetMaximum(), bgdOv->GetMaximum() )*maxMult );
            //frame->GetYaxis()->SetLimits( 0, ymax );
      
            // for better visibility, plot thinner lines
            sig->SetLineWidth( 1 );
            bgd->SetLineWidth( 1 );

            // perform K-S test
            cout << "--- Perform Kolmogorov-Smirnov tests" << endl;
            Double_t kolS = sig->KolmogorovTest( sigOv );
            Double_t kolB = bgd->KolmogorovTest( bgdOv );
            cout << "--- Goodness of signal (background) consistency: " << kolS << " (" << kolB << ")" << endl;

            TString probatext = Form( "Kolmogorov-Smirnov test: signal (background) probability = %5.3g (%5.3g)", kolS, kolB );
            TText* tt = new TText( 0.12, 0.74, probatext );
            tt->SetNDC(); tt->SetTextSize( 0.032 ); tt->AppendPad(); 
         }

         // redraw axes
         frame->Draw("sameaxis");

         // text for overflows
         Int_t    nbin = sig->GetNbinsX();
         Double_t dxu  = sig->GetBinWidth(0);
         Double_t dxo  = sig->GetBinWidth(nbin+1);
         TString uoflow = Form( "U/O-flow (S,B): (%.1f, %.1f)%% / (%.1f, %.1f)%%", 
                                sig->GetBinContent(0)*dxu*100, bgd->GetBinContent(0)*dxu*100,
                                sig->GetBinContent(nbin+1)*dxo*100, bgd->GetBinContent(nbin+1)*dxo*100 );
         TText* t = new TText( 0.975, 0.115, uoflow );
         t->SetNDC();
         t->SetTextSize( 0.030 );
         t->SetTextAngle( 90 );
         t->AppendPad();    
   
         // update canvas
	 frame->SetMinimum(1e-5);
         c->Update();

         // save canvas to file

	 cout << "Saving plots in: weights/"+output_path+"/plots" << endl;
         TMVAGlob::plot_logo(1.058);
         if (Save_Images) {
            if      (htype == MVAType)     TMVAGlob::imgconv( c, Form("weights/"+output_path+"/plots/mva_%s",     methodTitle.Data()) );
            else if (htype == ProbaType)   TMVAGlob::imgconv( c, Form("weights/"+output_path+"/plots/proba_%s",   methodTitle.Data()) ); 
            else if (htype == CompareType) TMVAGlob::imgconv( c, Form("weights/"+output_path+"/plots/overtrain_%s", methodTitle.Data()) ); 
            else                           TMVAGlob::imgconv( c, Form("weights/"+output_path+"/plots/rarity_%s",  methodTitle.Data()) ); 
         }
         countCanvas++;
         
      }
      cout << "";
   }
}

void annconvergencetest( TString fin, Bool_t useTMVAStyle, TString output_path )
{
   // set style and remove existing canvas'
   TMVAGlob::Initialize( useTMVAStyle );
  
   // checks if file with name "fin" is already open, and if not opens one
   TFile* file = TMVAGlob::OpenFile( fin );  
   TDirectory *lhdir = file->GetDirectory("Method_MLP/MLP");

   TString jobName = lhdir->GetName();
   Int_t icanvas = -1;
   icanvas++;
   TCanvas* c = new TCanvas( Form("MLPConvergenceTest_%s",jobName.Data()), Form("MLP Convergence Test, %s",jobName.Data()), 
                             100 + (icanvas)*40, 0 + (icanvas+1)*20, 600, 580*0.8  );
  
   TH1* estimatorHistTrain = (TH1*)lhdir->Get( "estimatorHistTrain" );
   TH1* estimatorHistTest  = (TH1*)lhdir->Get( "estimatorHistTest"  );

   Double_t m1  = estimatorHistTrain->GetMaximum();
   Double_t m2  = estimatorHistTest ->GetMaximum();
   Double_t max = TMath::Max( m1, m2 );
   m1  = estimatorHistTrain->GetMinimum();
   m2  = estimatorHistTest ->GetMinimum();
   Double_t min = TMath::Min( m1, m2 );
   estimatorHistTrain->SetMaximum( max + 0.1*(max - min) );
   estimatorHistTrain->SetMinimum( min - 0.1*(max - min) );
   estimatorHistTrain->SetLineColor( 2 );
   estimatorHistTrain->SetLineWidth( 2 );
   estimatorHistTrain->SetTitle( TString("MLP Convergence Test") );
  
   estimatorHistTest->SetLineColor( 4 );
   estimatorHistTest->SetLineWidth( 2 );

   estimatorHistTrain->GetXaxis()->SetTitle( "Epochs" );
   estimatorHistTrain->GetYaxis()->SetTitle( "Estimator" );
   estimatorHistTrain->GetXaxis()->SetTitleOffset( 1.20 );
   estimatorHistTrain->GetYaxis()->SetTitleOffset( 1.65 );

   estimatorHistTrain->Draw();
   estimatorHistTest ->Draw("same");

   // need a legend
   TLegend *legend= new TLegend( 1 - c->GetRightMargin() - 0.45, 1-c->GetTopMargin() - 0.20, 
                                 1 - c->GetRightMargin() - 0.05, 1-c->GetTopMargin() - 0.05 );

   legend->AddEntry(estimatorHistTrain,"Training Sample","l");
   legend->AddEntry(estimatorHistTest,"Test sample","l");
   legend->Draw("same");
   legend->SetMargin( 0.3 );

   c->cd();
   TMVAGlob::plot_logo(); // don't understand why this doesn't work ... :-(
   c->Update();

   TString fname = "weights/"+output_path+"/plots/annconvergencetest";
   TMVAGlob::imgconv( c, fname );
}



