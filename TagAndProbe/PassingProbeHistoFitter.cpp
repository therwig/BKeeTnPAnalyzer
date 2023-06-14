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

void PassingProbeHistoFitter()
{
  TCanvas *canvas = new TCanvas("cs_data","cs_data",2200,1200);

  TFile* inputBC = TFile::Open("HistoOutputNewSelectionBC.root");
  TFile* inputA = TFile::Open("HistoOutputNewSelectionA.root");

  TH1F* Pass4p5infBC = (TH1F*)inputBC->Get("PassSV3Mass4p5infAllCuts");
  TH1F* Pass4p5infA = (TH1F*)inputA->Get("PassSV3Mass4p5infAllCuts");

  Pass4p5infA->Add(Pass4p5infBC);
  Pass4p5infA->Rebin(1);

  RooRealVar InvMass("InvMass", "InvMass", 2., 4., "GeV/c^{2}");
  RooDataHist *data_hist = new RooDataHist("data_hist", "binned data", InvMass, Pass4p5infA);



  //DCB J/Psi function variables
  RooRealVar* dcbMean  = new RooRealVar("dcbMean","dcbMean",3.0, 3.0,3.15);
  RooRealVar* dcbSigma = new RooRealVar("dcbSigma","dcbSigma",0.1, 0., 3.);
  RooRealVar* dcbAlphaL = new RooRealVar("dcbAlphaL","dcbAlphaL",0.57, 0, 50.);
  RooRealVar* dcbNL = new RooRealVar("dcbNL","dcbNL",30., 0., 100.);
  RooRealVar* dcbAlphaR = new RooRealVar("dcbAlphaR","dcbAlphaR",1.1, 0, 60.);
  RooRealVar* dcbNR = new RooRealVar("dcbNR","dcbNR",30., 0., 100.);
  RooRealVar sig_yield("sig_yield", "yield of signal peak", 20, 0, 200);

  //DCB J/Psi Prime function variables
  RooRealVar* dcbMeanPrime  = new RooRealVar("dcbMeanPrime","dcbMeanPrime",3.65, 3.0,4.0);
  RooRealVar* dcbSigmaPrime = new RooRealVar("dcbSigmaPrime","dcbSigmaPrime",0.1, 0., 3.);
  RooRealVar* dcbAlphaLPrime = new RooRealVar("dcbAlphaLPrime","dcbAlphaLPrime",0.57, 0, 50.);
  RooRealVar* dcbNLPrime = new RooRealVar("dcbNLPrime","dcbNLPrime",30., 0., 100.);
  RooRealVar* dcbAlphaRPrime = new RooRealVar("dcbAlphaRPrime","dcbAlphaRPrime",1.1, 0, 60.);
  RooRealVar* dcbNRPrime = new RooRealVar("dcbNRPrime","dcbNRPrime",30., 0., 100.);
  RooRealVar sig_yieldPrime("sig_yieldPrime", "yield of signal peak", 20, 0, 200);

  dcbMean->setVal(3.05); // Fix JPsi mass
  dcbMean->setConstant(kTRUE); //Fix Jpsi mass
  dcbMeanPrime->setVal(3.65); // Fix JPsi Prime mass
  dcbMeanPrime->setConstant(kTRUE); //Fix Jpsi Prime mass

  //Sig and bkg functions
  RooCrystalBall* dcb = new RooCrystalBall("dcb","dcb",InvMass,*dcbMean,*dcbSigma, *dcbAlphaL, *dcbNL, *dcbAlphaR, *dcbNR);
  RooCrystalBall* dcbPrime = new RooCrystalBall("dcbPrime","dcbPrime",InvMass,*dcbMeanPrime,*dcbSigmaPrime, *dcbAlphaLPrime, *dcbNLPrime, *dcbAlphaRPrime, *dcbNRPrime);

  RooArgList shapes, yields;
  shapes.add(*dcb);
  //shapes.add(gauss);
  shapes.add(*dcbPrime);
  yields.add(sig_yield);
  yields.add(sig_yieldPrime);
  RooAddPdf totalPdf("totalPdf", "sum of signal and background PDF's", shapes, yields);

  totalPdf.fitTo(*data_hist, Extended());

  RooPlot* massframe = InvMass.frame();
  data_hist->plotOn(massframe);
  totalPdf.plotOn(massframe, LineColor(kBlue));
  totalPdf.plotOn(massframe,Components(*dcbPrime), LineColor(kRed));
  totalPdf.paramOn(massframe,Parameters(sig_yield),Layout(0.7));
  massframe->Draw();
  canvas->SaveAs("PassingFit.pdf");

}
