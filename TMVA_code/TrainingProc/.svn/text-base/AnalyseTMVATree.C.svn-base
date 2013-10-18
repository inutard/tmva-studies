// Compile with:
// g++ -o AnalyseTMVATree AnalyseTMVATree.C `root-config --cflags --glibs`
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
#include "TCanvas.h"
#include "TPad.h"
#include "TF1.h"
#include "tmvaglob.C"
#include "TMVA/Tools.h"

TString training = "Training_elmu_07";
TString folder="weights/"+training+"/";

TH1* move_overflow(TH1 *h);
void make_overtraining_plot(TTree* train_tree, TTree* test_tree);
void make_anconvergence_plot(TFile* f);
void printout_train_test(TTree* train_tree, TTree* test_tree);
void cut_scan(TTree* train_tree, TTree* test_tree);
void plot_varibles (TTree* train_tree, TTree* test_tree);


//int main() {
void AnalyseTMVATree() {
  gROOT->LoadMacro("/atlas/users/snezana/testarea/15.6.5.1/atlasstyle-00-02-04/AtlasUtils.C");
  gROOT->LoadMacro("/atlas/users/snezana/testarea/15.6.5.1/atlasstyle-00-02-04/AtlasStyle.C");
  //gStyle->SetAtlasStyle();
  gStyle->SetPalette(1);
  gStyle->SetPaintTextFormat();
  gStyle->SetOptStat(0);

  TFile* f = TFile::Open(folder+"TMVA.root");
  TTree* train_tree = f->Get("TrainTree");
  TTree* test_tree = f->Get("TestTree");

  make_overtraining_plot(train_tree, test_tree);
  make_anconvergence_plot(f);
  printout_train_test(train_tree, test_tree);
  cut_scan(train_tree, test_tree);
  plot_varibles(train_tree, test_tree);
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

void make_overtraining_plot(TTree* train_tree, TTree* test_tree){
  UInt_t nbins = 70;
  Float_t xmin = -0.2;
  Float_t xmax = 1.2;
 
  /////////////////////////////////////////////////////////
  // Make the training/testing sample comparison for signal
  TCanvas* c_overtrain = new TCanvas("overtrain", "overtrain", 800, 800);
  c_overtrain->SetLeftMargin(c_overtrain->GetLeftMargin()-0.05);
  c_overtrain->SetRightMargin(c_overtrain->GetRightMargin()+0.05);
  TPad* pad_1 = new TPad("pad_1", "", 0., 0.3, 1., 1.);
  pad_1->SetBottomMargin(0);
  pad_1->SetLogy();
  pad_1->Draw();
  
  TPad* pad_2 = new TPad("pad_2", "", 0., 0.0, 1., 0.3);
  pad_2->SetTopMargin(0);
  pad_2->SetBottomMargin(pad_2->GetBottomMargin()+0.05);
  pad_2->Draw();

  pad_1->cd();
  // Full training and testing histograms for sig and bgd
  TH1F* train_sig = new TH1F("train_sig", "train_sig", nbins, xmin, xmax);
  train_sig->Sumw2();
  train_tree->Draw("MLP>>train_sig","mc_weight*(classID==0)");
  TH1F* train_bgd = new TH1F("train_bgd", "train_bgd", nbins, xmin, xmax);
  train_bgd->Sumw2();
  train_tree->Draw("MLP>>train_bgd","mc_weight*(classID==1)");
  TH1F* test_sig = new TH1F("test_sig", "test_sig", nbins, xmin, xmax);
  test_sig->Sumw2();
  test_tree->Draw("MLP>>test_sig","mc_weight*(classID==0)");
  TH1F* test_bgd = new TH1F("test_bgd", "test_bgd", nbins, xmin, xmax);
  test_bgd->Sumw2();
  test_tree->Draw("MLP>>test_bgd","mc_weight*(classID==1)");

  // Training and testing histograms for NN>0.9795
  

  // Mooving overflow and underflow
  train_sig = (TH1F*)move_overflow(train_sig);
  train_bgd = (TH1F*)move_overflow(train_bgd);
  test_sig = (TH1F*)move_overflow(test_sig);
  test_bgd = (TH1F*)move_overflow(test_bgd);

  cout << "train_sig->Integral(0,-1): " << train_sig->Integral(0,-1) << endl;
  cout << "test_sig->Integral(0,-1): " << test_sig->Integral(0,-1) << endl;
  cout << "train_bgd->Integral(0,-1): " << train_bgd->Integral(0,-1) << endl;
  cout << "test_bgd->Integral(0,-1): " << test_bgd->Integral(0,-1) << endl;

  // Draw the plots in the pad_1
  test_bgd->SetTitle(training);
  test_bgd->GetXaxis()->SetTitle("MLP");
  test_bgd->GetXaxis()->SetTitleSize(0.03);
  test_bgd->GetXaxis()->SetLabelSize(0.03);
  test_bgd->GetYaxis()->SetTitle("Events");
  test_bgd->GetYaxis()->SetTitleSize(0.05);
  test_bgd->GetYaxis()->SetTitleOffset(1);
  test_bgd->GetYaxis()->SetLabelSize(0.03);
  test_bgd->GetYaxis()->SetRangeUser(1e-3, 2e3);
  test_bgd->SetLineColor(kPink-7);
  test_bgd->SetLineWidth(2);
  train_bgd->SetMarkerStyle(20);
  train_bgd->SetMarkerColor(kPink-6);

  test_bgd->Draw("hist");
  train_bgd->Draw("same,ep");

  train_sig->SetMarkerColor(kAzure-6);
  train_sig->SetMarkerStyle(20);
  test_sig->SetLineWidth(2);
  test_sig->SetLineColor(kAzure-7);
  test_sig->Draw("same,hist");
  train_sig->Draw("same,ep");

  float x = 0.6;
  float y = 0.9;
  TLegend * leg = new TLegend(x,y-0.2,x+0.25,y,NULL,"brNDC");
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetFillStyle(0);
  leg->SetTextSize(0.03);
  leg->AddEntry(train_sig,"training sample S","lep");
  leg->AddEntry(test_sig,"testing sample S","l");
  leg->AddEntry(train_bgd,"training sample B","lep");
  leg->AddEntry(test_bgd,"testing sample B","l");
  leg->Draw();

  //TString fname = "weights/"+output_path+"/plots/SoverB";
  //TString fname = "weights/"+output_path+"/plots/SoverB";
  //TMVAGlob::imgconv( c_SoverB, fname );

  pad_2->cd();
  TH1F* ratio_sig = new TH1F("ratio_sig", "ratio_sig", nbins, xmin, xmax);
  ratio_sig->Sumw2();
  ratio_sig->Add(test_sig);
  ratio_sig->Divide(train_sig);
  ratio_sig->SetMarkerStyle(20);
  ratio_sig->SetMarkerColor(kAzure-6);

  TH1F* ratio_bgd = new TH1F("ratio_bgd", "ratio_bgd", nbins, xmin, xmax);
  ratio_bgd->Sumw2();
  ratio_bgd->Add(test_bgd);
  ratio_bgd->Divide(train_bgd);

  ratio_bgd->GetYaxis()->SetTitle("test/train ");
  ratio_bgd->GetXaxis()->SetTitle("MLP");
  ratio_bgd->GetXaxis()->SetTitleSize(0.1);
  ratio_bgd->GetYaxis()->SetTitleSize(0.1);
  ratio_bgd->GetXaxis()->SetTitleOffset(0.8);
  ratio_bgd->GetYaxis()->SetTitleOffset(0.55);
  ratio_bgd->GetXaxis()->SetLabelSize(0.08);
  ratio_bgd->GetYaxis()->SetLabelSize(0.08);
  ratio_bgd->GetYaxis()->SetRangeUser(0.,2.5);
  ratio_bgd->SetMarkerStyle(20);
  ratio_bgd->SetMarkerColor(kPink-6);
  ratio_bgd->Draw("ep");
	    
  TF1 *unity = new TF1("unity", "1", ratio_sig->GetBinLowEdge(1), ratio_sig->GetBinLowEdge(ratio_sig->GetNbinsX()+1));
  unity->SetLineStyle(2);
  unity->Draw("LSAME");

  x = 0.6;
  y = 0.85;
  leg = new TLegend(x,y-0.15,x+0.25,y,NULL,"brNDC");
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetFillStyle(0);
  leg->SetTextSize(0.08);
  leg->AddEntry(ratio_sig,"Signal","lep");
  leg->AddEntry(ratio_bgd,"Background","lep");
  leg->Draw();

  ratio_sig->Draw("same,ep");
  ratio_bgd->Draw("same,ep");

  Double_t kolS = test_sig->KolmogorovTest( train_sig );
  Double_t kolB = test_bgd->KolmogorovTest( train_bgd );
  TString probatext = Form( "Kolmogorov-Smirnov test: S (B) probability = %5.3g (%5.3g)", kolS, kolB);
  TText* tt = new TText( 0.18, 0.87, probatext );
  tt->SetNDC(); tt->SetTextSize( 0.08 ); tt->AppendPad();
  
  c_overtrain->Print(folder+"plots/Overtraining_Ratio.pdf");
  }

void make_anconvergence_plot(TFile* f){
  /////////////////////////////////////////////////////////
  // Make the convergence plot
  TCanvas* c_convergence = new TCanvas("convergence", "convergence", 800, 800);
  c_convergence->SetLeftMargin(0.05);
  c_convergence->SetRightMargin(0.05);
  c_convergence->cd();
  c_convergence->SetRightMargin(0.05);
  c_convergence->SetLeftMargin(0.15);
  c_convergence->SetLeftMargin(0.05);
  TDirectoryFile* m_mlp = (TDirectoryFile*)f->Get("Method_MLP");
  TDirectoryFile* mlp = (TDirectoryFile*)m_mlp->Get("MLP");

  TH1F* est_train = (TH1F*)mlp->Get("estimatorHistTrain");
  TH1F* est_test = (TH1F*)mlp->Get("estimatorHistTest");

  cout << "est_train->Integral(): " << est_train->Integral() << endl;

  est_train->SetLineColor(kPink-6);
  est_train->SetLineWidth(2);
  est_train->GetXaxis()->SetTitle("Training Cycle");
  est_train->GetYaxis()->SetTitle("Error Function");
  est_train->GetXaxis()->SetTitleSize(0.035);
  est_train->GetYaxis()->SetTitleSize(0.035);
  est_train->GetYaxis()->SetTitleOffset(1.7);
  est_train->GetXaxis()->SetLabelSize(0.035);
  est_train->GetYaxis()->SetLabelSize(0.035);
  est_test->SetLineColor(kAzure-6);
  est_train->Draw("hist");
  est_test->Draw("same,hist");

  float x = 0.6;
  float y = 0.9;
  TLegend * leg = new TLegend(x,y-0.15,x+0.25,y,NULL,"brNDC");
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetFillStyle(0);
  leg->SetTextSize(0.03);
  leg->AddEntry(est_train,"training sample","l");
  leg->AddEntry(est_test,"testing sample","l");
  leg->Draw();
  c_convergence->Print(folder+"plots/convergence_test.pdf");
}

void printout_train_test(TTree* train_tree, TTree* test_tree){

  TH1F* ts_sig = new TH1F("ts_sig", "ts_sig", 25, 0., 2500000.);
  test_tree->Draw("Ht>>ts_sig","mc_weight*((classID==0)&&(MLP>0.9795))");
  TH1F* ts_bgd = new TH1F("ts_bgd", "ts_bgd", 25, 0., 2500000.);
  test_tree->Draw("Ht>>ts_bgd","mc_weight*((classID==1)&&(MLP>0.9795))");

  TH1F* ts_sig_all = new TH1F("ts_sig_all", "ts_sig_all", 25, 0., 2500000.);
  test_tree->Draw("Ht>>ts_sig_all","mc_weight*(classID==0)");
  TH1F* ts_bgd_all = new TH1F("ts_bgd_all", "ts_bgd_all", 25, 0., 2500000.);
  test_tree->Draw("Ht>>ts_bgd_all","mc_weight*(classID==1)");

  TH1F* tr_sig = new TH1F("tr_sig", "tr_sig", 25, 0., 2500000.);
  train_tree->Draw("Ht>>tr_sig","mc_weight*((classID==0)&&(MLP>0.9795))");
  TH1F* tr_bgd = new TH1F("tr_bgd", "tr_bgd", 25, 0., 2500000.);
  train_tree->Draw("Ht>>tr_bgd","mc_weight*((classID==1)&&(MLP>0.9795))");

  TH1F* tr_sig_all = new TH1F("tr_sig_all", "tr_sig_all", 25, 0., 2500000.);
  train_tree->Draw("Ht>>tr_sig_all","mc_weight*(classID==0)");
  TH1F* tr_bgd_all = new TH1F("tr_bgd_all", "tr_bgd_all", 25, 0., 2500000.);
  train_tree->Draw("Ht>>tr_bgd_all","mc_weight*(classID==1)");

  // Printing out numbers for the training tree only
  ts_sig->Add(tr_sig);
  ts_bgd->Add(tr_bgd);

  ts_sig_all->Add(tr_sig_all);
  ts_bgd_all->Add(tr_bgd_all);

  ts_sig = (TH1F*)move_overflow(ts_sig);
  ts_bgd = (TH1F*)move_overflow(ts_bgd);
  ts_sig_all = (TH1F*)move_overflow(ts_sig_all);
  ts_bgd_all = (TH1F*)move_overflow(ts_bgd_all);

  float integral_sig = ts_sig_all->Integral(0, -1);
  float integral_bgd = ts_bgd_all->Integral(0, -1);

  cout << "Number of events with the testing+training tree:" << endl;
  cout << "MVA cut	sig	bgrd:" << endl;
  cout << ts_sig->Integral() << "\t" << ts_bgd->Integral() << endl;

}

void cut_scan(TTree* train_tree, TTree* test_tree) {

  static Int_t nbins = 1200;

  TH1F* test_sig = new TH1F("test_sig", "test_sig", nbins, 0., 1.2);
  test_tree->Draw("MLP>>test_sig","mc_weight*(classID==0)");
  TH1F* test_bgd = new TH1F("test_bgd", "test_bgd", nbins, 0., 1.2);
  test_tree->Draw("MLP>>test_bgd","mc_weight*(classID==1)");

  TH1F* train_sig = new TH1F("train_sig", "train_sig", nbins, 0., 1.2);
  train_tree->Draw("MLP>>train_sig","mc_weight*(classID==0)");
  TH1F* train_bgd = new TH1F("train_bgd", "train_bgd", nbins, 0., 1.2);
  train_tree->Draw("MLP>>train_bgd","mc_weight*(classID==1)");
  
  test_sig->Add(train_sig);
  test_bgd->Add(train_bgd);

  test_sig = (TH1F*)move_overflow(test_sig);
  test_bgd = (TH1F*)move_overflow(test_bgd);

  float integral_sig = test_sig->Integral(0, -1);
  float integral_bgd = test_bgd->Integral(0, -1);

  float MLP_cut[1200];
  float SoverB[1200];
  float significance[1200];
  float max_signif = 0;
  float max_cut_signif = 0;
  float max_SoverB = 0;
  float max_cut_SoverB = 0;
  float vec_sig[1200];
  float vec_bgd[1200];

  cout << "itr		MLP cut		sig		bgd		significance	    S/B	     sig[bin]		bgd[bin]	sinificance[bin]	S/B[bin]:" << endl;

  int itr = 1;
  float n_sg = test_sig->Integral(itr, -1);
  float n_bg = test_bgd->Integral(itr, -1);
  while ((n_sg!=0) && (n_bg!=0) && (itr<nbins)) {
    vec_sig[itr-1] = n_sg;
    vec_bgd[itr-1] = n_bg;
    MLP_cut[itr-1] = (test_sig->GetBinLowEdge(itr)>0) ? test_sig->GetBinLowEdge(itr) : 0;
    SoverB[itr-1] = (n_sg/n_bg)>0 ? n_sg/n_bg : 0;
    significance[itr-1] = (n_sg/sqrt(n_sg+n_bg)>0) ? n_sg/sqrt(n_sg+n_bg) : 0;
    float sg_bin = test_sig->GetBinContent(itr);/
    float bg_bin = test_bgd->GetBinContent(itr);/
    if (MLP_cut[itr-1]>0.8)
      cout << itr  << "\t \t" << MLP_cut[itr-1] << "\t \t" << n_sg << "\t \t" << n_bg << "\t \t" << significance[itr-1] << "\t \t" << SoverB[itr-1] << "\t \t" << sg_bin << "\t \t" << bg_bin << "\t \t" << sg_bin/sqrt(sg_bin+bg_bin)<< "\t \t" << sg_bin/bg_bin << endl;
    if (significance[itr-1]>max_signif) {
      max_signif = significance[itr-1];
      max_cut_signif = MLP_cut[itr-1];
    }
    if (SoverB[itr-1]>max_SoverB) {
      max_SoverB = SoverB[itr-1];
      max_cut_SoverB = MLP_cut[itr-1];
    }

    itr++;
    n_sg = test_sig->Integral(itr, -1);
    n_bg = test_bgd->Integral(itr, -1);
  }

  cout << "Maximal SoverB is " << max_SoverB << " for the MLP cut at " << max_cut_SoverB << endl;
  cout << "Maximal significance is " << max_signif << " for the MLP cut at " << max_cut_signif << endl;

  cout << "bins (187-200): " << "sig: " <<  test_sig->Integral(187, 20) << "bgd: " << test_bgd->Integral(187, 20) << endl;

  TCanvas* c_SoverB = new TCanvas("SoverB", "SoverB", 600, 400);
  c_SoverB->cd();
  TGraph* gr_SoverB = new TGraph(itr, MLP_cut, SoverB);
  gr_SoverB->GetXaxis()->SetTitle("MLP Cut");
  gr_SoverB->GetXaxis()->SetTitleSize(0.03);
  gr_SoverB->GetXaxis()->SetLabelSize(0.03);
  gr_SoverB->GetYaxis()->SetTitle("S/B");
  gr_SoverB->GetYaxis()->SetTitleSize(0.03);
  gr_SoverB->GetYaxis()->SetTitleOffset(1.7);
  gr_SoverB->GetYaxis()->SetLabelSize(0.03);
  gr_SoverB->GetYaxis()->SetRangeUser(0., max_SoverB*1.2);
  gr_SoverB->SetMarkerStyle(20);
  gr_SoverB->SetMarkerColor(kAzure-6);
  gr_SoverB->Draw("AP");
  c_SoverB->Print(folder+"plots/SoverB.pdf");

  TCanvas* c_signif = new TCanvas("signif", "signif", 600, 400);
  c_signif->cd();
  TGraph* gr_signif = new TGraph(itr, MLP_cut, significance);
  gr_signif->GetXaxis()->SetTitle("MLP Cut");
  gr_signif->GetXaxis()->SetTitleSize(0.03);
  gr_signif->GetXaxis()->SetLabelSize(0.03);
  gr_signif->GetYaxis()->SetTitle("S/#sqrt{S+B}");
  gr_signif->GetYaxis()->SetTitleSize(0.03);
  gr_signif->GetYaxis()->SetTitleOffset(1.7);
  gr_signif->GetYaxis()->SetLabelSize(0.03);
  gr_signif->GetYaxis()->SetRangeUser(0., max_signif*1.2);
  gr_signif->SetMarkerStyle(20);
  gr_signif->SetMarkerColor(kPink-6);
  gr_signif->Draw("AP");
  c_signif->Print(folder+"plots/significance.pdf");
  
  TCanvas* c_bg_vs_sg = new TCanvas("bg_vs_sg", "bg_sv_sg", 600, 400);
  c_bg_vs_sg->cd();
  TH1F* dummy = new TH1F("dummy", "dummy", 1, 0, 93);
  dummy->GetYaxis()->SetRangeUser(0., 2100.);
  dummy->Draw();
  
  TGraph* gr_bg_vs_sg = new TGraph(itr-1, vec_sig, vec_bgd);
  gr_bg_vs_sg->GetXaxis()->SetTitle("Signal Events");
  gr_bg_vs_sg->GetXaxis()->SetTitleSize(0.03);
  gr_bg_vs_sg->GetXaxis()->SetLabelSize(0.03);
  gr_bg_vs_sg->GetYaxis()->SetTitle("Background Events");
  gr_bg_vs_sg->GetYaxis()->SetTitleSize(0.03);
  gr_bg_vs_sg->GetYaxis()->SetTitleOffset(1.7);
  gr_bg_vs_sg->GetYaxis()->SetLabelSize(0.03);
  gr_bg_vs_sg->SetMarkerStyle(20);
  gr_bg_vs_sg->SetMarkerColor(kTeal-6);
  gr_bg_vs_sg->Draw("same,AP");
  c_bg_vs_sg->Print(folder+"plots/bgd_vs_sig_eff.pdf");

}

void plot_varibles (TTree* train_tree, TTree* test_tree){

  TString variables[] = {"MLP",
			 "m_reco",
			 "Ht",
			 "b_had_pt",
			 "b_lep_pt",
			 "dR_lb_had",
			 "dR_lb_lep",
			 "dR_Whad_bhad",
			 "dR_Whad_blep",
			 "dR_lnu"};
  int nbins[] = {14, 6, 8, 16, 16, 20, 20, 20, 20, 20};

  float binning[][46] = {{-0.2, -0.1, 0., 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1., 1.1, 1.2},	// MVA, 11 bins
			 {0., 300., 450., 600., 750., 900., 1200.},		// m_reco, 6 bins
			 {0., 100., 200., 300., 400., 500., 600., 700., 800.}, // Ht, 12 bins
			 {0., 25., 50., 75., 100., 125., 150., 175., 200., 225., 250., 275., 300., 325., 350., 375., 400.}, // b_had_pt, 16 bins
			 {0., 25., 50., 75., 100., 125., 150., 175., 200., 225., 250., 275., 300., 325., 350., 375., 400.}, // b_lep_pt, 16 bins
			 {0., 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4, 3.6, 3.8, 4.}, // dR_lb_had, 20 bins
			 {0., 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4, 3.6, 3.8, 4.}, // dR_lb_lep, 20 bins
			 {0., 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4, 3.6, 3.8, 4.}, // dR_Whad_bhad, 20 bins
			 {0., 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4, 3.6, 3.8, 4.}, // dR_Whad_blep, 20 bins
			 {0., 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4, 3.6, 3.8, 4.}}; // dR_lnu, 20 bins
  TString scale[] = {"1.", "1000.", "1000.", "1000.", "1000.", "1.", "1.", "1.", "1.", "1."};
  
  for (int i = 0; i<10; ++i) {
    cout << variables[i] << endl;
    TCanvas* c = new TCanvas("overtrain", "overtrain", 800, 800);
    TString hist_name = variables[i];
    
    TCanvas* c = new TCanvas(hist_name, hist_name, 800, 800);
    c->SetLeftMargin(c->GetLeftMargin()-0.05);
    c->SetRightMargin(c->GetRightMargin()+0.05);
    c->SetBottomMargin(c->GetBottomMargin()+0.05);
    TPad* pad_1 = new TPad("pad_1", "", 0., 0.3, 1., 1.);
    pad_1->SetBottomMargin(0);
    pad_1->SetLogy();
    pad_1->Draw();
    
    TPad* pad_2 = new TPad("pad_2", "", 0., 0.0, 1., 0.3);
    pad_2->SetTopMargin(0);
    pad_2->SetBottomMargin(pad_2->GetBottomMargin()+0.05);
    pad_2->Draw();
    
    pad_1->cd();
    TString draw_option = variables[i] + "/" + scale[i];
    train_tree->Draw(draw_option,"mc_weight*(classID==0)");
    
    TH1F* train_sig = new TH1F("train_sig", "train_sig", nbins[i], binning[i]);
    train_sig->Sumw2();
    train_tree->Draw(draw_option+">>train_sig","mc_weight*(classID==0)");
    TH1F* train_bgd = new TH1F("train_bgd", "train_bgd", nbins[i], binning[i]);
    train_bgd->Sumw2();
    train_tree->Draw(draw_option+">>train_bgd","mc_weight*(classID==1)");
    TH1F* test_sig = new TH1F("test_sig", "test_sig", nbins[i], binning[i]);
    test_sig->Sumw2();
    test_tree->Draw(draw_option+">>test_sig","mc_weight*(classID==0)");
    TH1F* test_bgd = new TH1F("test_bgd", "test_bgd", nbins[i], binning[i]);
    test_bgd->Sumw2();
    test_tree->Draw(draw_option+">>test_bgd","mc_weight*(classID==1)");
    
    // Mooving overflow and underflow
    train_sig = (TH1F*)move_overflow(train_sig);
    train_bgd = (TH1F*)move_overflow(train_bgd);
    test_sig = (TH1F*)move_overflow(test_sig);
    test_bgd = (TH1F*)move_overflow(test_bgd);
   
    float y_max = test_bgd->GetMaximum();
    if (test_sig->GetMaximum()>y_max)
      y_max = test_sig->GetMaximum;
    float y_min = test_bgd->GetMinimum();
    if (test_sig->GetMinimum()<y_min)
      y_min = test_sig->GetMinimum();
    if (y_min<0) y_min=0.001
    
 
    // Draw the plots in the pad_1
    test_bgd->SetTitle(training);
    test_bgd->GetXaxis()->SetTitle("MLP");
    test_bgd->GetXaxis()->SetTitleSize(0.03);
    test_bgd->GetXaxis()->SetLabelSize(0.03);
    //test_bgd->GetYaxis()->SetTitle(test_bgd->GetName());
    //test_bgd->GetYaxis()->SetTitle("S/B");
    test_bgd->GetYaxis()->SetTitle("Events");
    test_bgd->GetYaxis()->SetTitleSize(0.05);
    test_bgd->GetYaxis()->SetTitleOffset(1);
    test_bgd->GetYaxis()->SetLabelSize(0.03);
    test_bgd->GetYaxis()->SetRangeUser(0.01, y_max*1.3);
    test_bgd->SetLineColor(kPink-7);
    test_bgd->SetLineWidth(2);
    train_bgd->SetMarkerStyle(20);
    train_bgd->SetMarkerColor(kPink-6);
    
    test_bgd->Draw("hist");
    train_bgd->Draw("same,ep");
    
    train_sig->SetMarkerColor(kAzure-6);
    train_sig->SetMarkerStyle(20);
    test_sig->SetLineWidth(2);
    test_sig->SetLineColor(kAzure-7);
    test_sig->Draw("same,hist");
    train_sig->Draw("same,ep");
    
    float x = 0.6;
    float y = 0.9;
    TLegend * leg = new TLegend(x,y-0.2,x+0.25,y,NULL,"brNDC");
    leg->SetBorderSize(0);
    leg->SetFillColor(0);
    leg->SetFillStyle(0);
    leg->SetTextSize(0.03);
    leg->AddEntry(train_sig,"training sample S","lep");
    leg->AddEntry(test_sig,"testing sample S","l");
    leg->AddEntry(train_bgd,"training sample B","lep");
    leg->AddEntry(test_bgd,"testing sample B","l");
    leg->Draw();
    
    pad_2->cd();
    TH1F* ratio_sig = new TH1F("ratio_sig", "ratio_sig", nbins[i], binning[i]);
    ratio_sig->Sumw2();
    ratio_sig->Add(test_sig);
    ratio_sig->Divide(train_sig);
    ratio_sig->SetMarkerStyle(20);
    ratio_sig->SetMarkerColor(kAzure-6);
    ratio_sig->SetLineColor(kAzure-6);
    
    TH1F* ratio_bgd = new TH1F("ratio_bgd", "ratio_bgd", nbins[i], binning[i]);
    ratio_bgd->Sumw2();
    ratio_bgd->Add(test_bgd);
    ratio_bgd->Divide(train_bgd);
    
    ratio_bgd->GetYaxis()->SetTitle("test/train ");
    ratio_bgd->GetXaxis()->SetTitle(variables[i]);
    ratio_bgd->GetXaxis()->SetTitleSize(0.1);
    ratio_bgd->GetYaxis()->SetTitleSize(0.1);
    ratio_bgd->GetXaxis()->SetTitleOffset(0.8);
    ratio_bgd->GetYaxis()->SetTitleOffset(0.55);
    ratio_bgd->GetXaxis()->SetLabelSize(0.08);
    ratio_bgd->GetYaxis()->SetLabelSize(0.08);
    ratio_bgd->GetYaxis()->SetRangeUser(0.5,1.5);
    ratio_bgd->SetMarkerStyle(20);
    ratio_bgd->SetMarkerColor(kPink-6);
    ratio_bgd->Draw("ep");
    
    TF1 *unity = new TF1("unity", "1", ratio_sig->GetBinLowEdge(1), ratio_sig->GetBinLowEdge(ratio_sig->GetNbinsX()+1));
    unity->SetLineStyle(2);
    unity->Draw("LSAME");
    
    x = 0.6;
    y = 0.85;
    leg = new TLegend(x,y-0.15,x+0.25,y,NULL,"brNDC");
    leg->SetBorderSize(0);
    leg->SetFillColor(0);
    leg->SetFillStyle(0);
    leg->SetTextSize(0.08);
    leg->AddEntry(ratio_sig,"Signal","lep");
    leg->AddEntry(ratio_bgd,"Background","lep");
    leg->Draw();
    
    ratio_sig->Draw("same,ep");
    ratio_bgd->Draw("same,ep");
    
    Double_t kolS = test_sig->KolmogorovTest( train_sig );
    Double_t kolB = test_bgd->KolmogorovTest( train_bgd );
    TString probatext = Form( "Kolmogorov-Smirnov test: S (B) probability = %5.3g (%5.3g)", kolS, kolB);
    TText* tt = new TText( 0.18, 0.87, probatext );
    tt->SetNDC(); tt->SetTextSize( 0.08 ); tt->AppendPad();
    c->Print(folder+"plots/overtraining_"+variables[i]+".pdf");

    delete c;
    delete train_sig;
    delete train_bgd;
    delete test_sig;
    delete test_bgd;
    delete ratio_sig;
    delete ratio_bgd;
  }
}
