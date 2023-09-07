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
#include "RooGaussExp.cxx"
#include "RooGaussExp.h"

using namespace RooFit;

void HistoFitterPass()
{
  TCanvas *canvas = new TCanvas("cs_data","cs_data",2200,1200);

  TFile* input = TFile::Open("HistoOutputCheckBins.root");

  TH1F* Pass4p5inf = (TH1F*)input->Get("PassSV3Mass4p5infAllCuts");
  //Pass4p5inf->Rebin(2);

  RooRealVar InvMass("InvMass", "InvMass", 2., 4.5, "GeV/c^{2}");
  RooDataHist *data_hist = new RooDataHist("data_hist", "binned data", InvMass, Pass4p5inf);


  //DCB J/Psi function variables
  RooRealVar* dcbMean  = new RooRealVar("dcbMean","dcbMean",3.0, 2.9,3.2);
  RooRealVar* dcbSigma = new RooRealVar("dcbSigma","dcbSigma",0.1, 0.03, 0.3);
  RooRealVar* dcbAlphaL = new RooRealVar("dcbAlphaL","dcbAlphaL",0.57, 0, 50.);
  RooRealVar* dcbNL = new RooRealVar("dcbNL","dcbNL",30., 0., 100.);
  RooRealVar* dcbAlphaR = new RooRealVar("dcbAlphaR","dcbAlphaR",1.1, 0, 60.);
  RooRealVar* dcbNR = new RooRealVar("dcbNR","dcbNR",30., 0., 100.);
  RooRealVar sig_yield("sig_yield", "yield of signal peak", 20, 0, 1500);

  //DCB J/Psi Prime function variables
  RooRealVar* dcbMeanPrime  = new RooRealVar("dcbMeanPrime","dcbMeanPrime",3.65, 3.5,3.7);
  RooRealVar* dcbSigmaPrime = new RooRealVar("dcbSigmaPrime","dcbSigmaPrime",0.07, 0.03, 0.15);
  RooRealVar* dcbAlphaLPrime = new RooRealVar("dcbAlphaLPrime","dcbAlphaLPrime",0.57, 0, 50.);
  RooRealVar* dcbNLPrime = new RooRealVar("dcbNLPrime","dcbNLPrime",30., 0., 100.);
  RooRealVar* dcbAlphaRPrime = new RooRealVar("dcbAlphaRPrime","dcbAlphaRPrime",1.1, 0, 60.);
  RooRealVar* dcbNRPrime = new RooRealVar("dcbNRPrime","dcbNRPrime",30., 0., 100.);
  RooRealVar sig_yieldPrime("sig_yieldPrime", "yield of signal peak", 20, 0, 500);

  //GaussExp variables
  RooRealVar* GaussExpMean  = new RooRealVar("GaussExpMean","GaussExpMean",3.0, 2.9,3.2);
  RooRealVar* GaussExpSigma  = new RooRealVar("GaussExpSigma","GaussExpSigma",0.01,0.07,3);
  RooRealVar* GaussExpAlpha  = new RooRealVar("GaussExpAlpha","GaussExpAlpha", 2, 0.01, 100.);

  RooRealVar* GaussExpMeanPrime  = new RooRealVar("GaussExpMeanPrime","GaussExpMeanPrime",3.6, 3.45,3.8);
  RooRealVar* GaussExpSigmaPrime  = new RooRealVar("GaussExpSigmaPrime","GaussExpSigmaPrime",0.01,0.07,10);
  RooRealVar* GaussExpAlphaPrime  = new RooRealVar("GaussExpAlphaPrime","GaussExpAlphaPrime", 2, 0.01, 100.);

  // 10 to 20: dcbNLPrime 2. ,2., 10.
  // Below 10: dcbNLPrime 30.,0.,100.

  //Sig and bkg functions
  RooCrystalBall* dcb = new RooCrystalBall("dcb","dcb",InvMass,*dcbMean,*dcbSigma, *dcbAlphaL, *dcbNL, *dcbAlphaR, *dcbNR);
  RooCrystalBall* dcbPrime = new RooCrystalBall("dcbPrime","dcbPrime",InvMass,*dcbMeanPrime,*dcbSigmaPrime, *dcbAlphaLPrime, *dcbNLPrime, *dcbAlphaRPrime, *dcbNRPrime);

  RooGaussExp* GaussExp = new RooGaussExp("GaussExp","GaussExp",InvMass, *GaussExpMean,*GaussExpSigma,*GaussExpAlpha);
  RooGaussExp* GaussExpPrime = new RooGaussExp("GaussExpPrime","GaussExpPrime",InvMass, *GaussExpMeanPrime,*GaussExpSigmaPrime,*GaussExpAlphaPrime);

  RooArgList shapes, yields;
  shapes.add(*dcb);
  shapes.add(*dcbPrime);
  //shapes.add(*GaussExp);
  //shapes.add(*GaussExpPrime);
  yields.add(sig_yield);
  yields.add(sig_yieldPrime);
  RooAddPdf totalPdf("totalPdf", "sum of signal and background PDF's", shapes, yields);

  totalPdf.fitTo(*data_hist, Extended());

  RooPlot* massframe = InvMass.frame();
  data_hist->plotOn(massframe);
  totalPdf.plotOn(massframe, LineColor(kBlue));
  totalPdf.plotOn(massframe,Components(*dcbPrime), LineColor(kRed));
  totalPdf.paramOn(massframe,/*Parameters(sig_yield),*/Layout(0.67,0.99,0.99));

  /*Int_t npar = totalPdf.getParameters(*data_hist)->selectByAttrib("Constant",kFALSE)-> getSize();
  Double_t chi2ndf = massframe->chiSquare(npar);
  std::string numberStr = std::to_string(chi2ndf);
  std::string prefix = "chi2/ndof=";
  const char* numberChar = prefix.append(numberStr).c_str();
  TText *txt = new TText(2.1, 340, numberChar);
  txt->SetTextSize(0.04);
  massframe->addObject(txt);*/

  massframe->Draw();
  canvas->SaveAs("PassingFit.pdf");

}
