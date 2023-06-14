#include "RooRealVar.h"
#include "RooConstVar.h"
#include "RooGaussian.h"
#include "RooArgusBG.h"
#include "RooAddPdf.h"
#include "RooDataSet.h"
#include "RooPlot.h"
#include "TTree.h"
#include "TChain.h"
#include "TFile.h"
#include "PhysicsTools/TagAndProbe/interface/RooCMSShape.h"
#include "RooCrystalBall.cxx"
#include "RooCrystalBall.h"

using namespace RooFit;

void FaillingProbeHistoFitter()
{
  TCanvas *canvas = new TCanvas("cs_data","cs_data",2200,1200);

  TFile* inputBC = TFile::Open("HistoOutputNewSelectionBC.root");
  TFile* inputA = TFile::Open("HistoOutputNewSelectionA.root");

  TH1F* Fail4p5infBC = (TH1F*)inputBC->Get("FailSV3Mass4p5infAllCuts");
  TH1F* Fail4p5infA = (TH1F*)inputA->Get("FailSV3Mass4p5infAllCuts");

  Fail4p5infA->Add(Fail4p5infBC);
  Fail4p5infA->Rebin(2);

  RooRealVar InvMass("InvMass", "InvMass", 2., 4., "GeV/c^{2}");
  RooDataHist *data_hist = new RooDataHist("data_hist", "binned data", InvMass, Fail4p5infA);

  //Roo CMS shape variables
  RooRealVar exp_alpha("#alpha", "alpha", 40.0, 0, 500);
  RooRealVar exp_beta("#beta", "beta", 0.05, 0.0, 4.0);
  RooRealVar exp_gamma("#gamma", "gamma", 0.02, 0.0, 0.3);
  RooRealVar exp_peak("peak", "peak", 3 ,2 ,4 );
  RooRealVar bkg_yield("bkg_yield", "yield of background", 300, 0 , 900);

  //CB function variables
  RooRealVar* dcbMean  = new RooRealVar("dcbMean","dcbMean",3.0, 3.0,3.15);
  RooRealVar* dcbSigma = new RooRealVar("dcbSigma","dcbSigma",0.1, 0., 3.);
  RooRealVar* dcbAlphaL = new RooRealVar("dcbAlphaL","dcbAlphaL",0.57, 30, 50.);
  RooRealVar* dcbNL = new RooRealVar("dcbNL","dcbNL",30., 0., 100.);
  RooRealVar* dcbAlphaR = new RooRealVar("dcbAlphaR","dcbAlphaR",1.1, 0, 60.);
  RooRealVar* dcbNR = new RooRealVar("dcbNR","dcbNR",30., 0., 100.);
  RooRealVar sig_yield("sig_yield", "yield of signal peak", 20, 0, 100);

  dcbMean->setVal(3.05); // Fix JPsi mass
  dcbMean->setConstant(kTRUE); //Fix Jpsi mass
  //dcbAlphaL->setVal(31);
  //dcbAlphaL->setConstant(kTRUE);
  //dcbNL->setVal(100);
  //dcbNL->setConstant(kTRUE);
  //dcbNR->setVal(100);
  //dcbNR->setConstant(kTRUE);
  //sig_yield.setVal(20);
  //sig_yield.setConstant(kTRUE);
  //dcbSigma->setVal(0.05);
  //dcbSigma->setConstant(kTRUE);

  //exp_peak.setVal(3.05);
  //exp_peak.setConstant(kTRUE);
  //exp_alpha.setVal(2.5);
  //exp_alpha.setConstant(kTRUE);
  exp_beta.setVal(2);
  exp_beta.setConstant(kTRUE);
  //exp_gamma.setVal(0.03);
  //exp_gamma.setConstant(kTRUE);
  //Sig and bkg functions
  RooCMSShape exp_pdf("exp_pdf", "bkg shape", InvMass, exp_alpha, exp_beta, exp_gamma, exp_peak);
  RooCrystalBall* dcb = new RooCrystalBall("dcb","dcb",InvMass,*dcbMean,*dcbSigma, *dcbAlphaL, *dcbNL, *dcbAlphaR, *dcbNR);


  RooArgList shapes, yields;
  shapes.add(exp_pdf);
  shapes.add(*dcb);
  yields.add(bkg_yield);
  yields.add(sig_yield);
  RooAddPdf totalPdf("totalPdf", "sum of signal and background PDF's", shapes, yields);

  totalPdf.fitTo(*data_hist, Extended());

  RooPlot* massframe = InvMass.frame();
  data_hist->plotOn(massframe);
  totalPdf.plotOn(massframe, LineColor(kBlue));
  totalPdf.plotOn(massframe,Components(exp_pdf), LineColor(kRed));
  totalPdf.paramOn(massframe,Parameters(sig_yield),Layout(0.75));
  massframe->Draw();
  canvas->SaveAs("FaillingFit.pdf");

}
