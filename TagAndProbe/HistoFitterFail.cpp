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

void HistoFitterFail()
{
  TCanvas *canvas = new TCanvas("cs_data","cs_data",2200,1200);

  TFile* input = TFile::Open("HistoOutputCheckBins.root");

  TH1F* Fail4p5inf = (TH1F*)input->Get("FailSV3Mass4p5infAllCuts");
  Fail4p5inf->Rebin(2);

  RooRealVar InvMass("InvMass", "InvMass", 2., 4.5, "GeV/c^{2}");
  RooDataHist *data_hist = new RooDataHist("data_hist", "binned data", InvMass, Fail4p5inf);

  //Roo CMS shape variables
  RooRealVar exp_alpha("#alpha", "alpha", 40.0, 0, 500);
  RooRealVar exp_beta("#beta", "beta", 0.05, 0.0, 4.0);
  RooRealVar exp_gamma("#gamma", "gamma", 0.02, 0.0, 0.3);
  RooRealVar exp_peak("peak", "peak", 3 ,2 ,4 );
  RooRealVar bkg_yield("bkg_yield", "yield of background", 300, 0 , 6000);

  //DCB function variables
  RooRealVar* dcbMean  = new RooRealVar("dcbMean","dcbMean",3.05, 3.0,3.15);
  RooRealVar* dcbSigma = new RooRealVar("dcbSigma","dcbSigma",0.1, 0., 0.3);
  RooRealVar* dcbAlphaL = new RooRealVar("dcbAlphaL","dcbAlphaL",0.57, 0., 50.);
  RooRealVar* dcbNL = new RooRealVar("dcbNL","dcbNL",30., 0., 100.);
  RooRealVar* dcbAlphaR = new RooRealVar("dcbAlphaR","dcbAlphaR",1.1, 0, 60.);
  RooRealVar* dcbNR = new RooRealVar("dcbNR","dcbNR",30., 0., 100.);
  RooRealVar sig_yield("sig_yield", "yield of signal peak", 20, 0, 2000);

  //CB function variables
  RooRealVar* cbMean  = new RooRealVar("cbMean","cbMean",3.05, 3.0,3.15);
  RooRealVar* cbRes  = new RooRealVar("cbRes","cbRes", 0.05, 0.02, 0.15);
  RooRealVar* cbTail  = new RooRealVar("cbTail","cbTail",1.1, 0.005, 3.);
  RooRealVar* cbN  = new RooRealVar("cbN","cbN",1.1, 0., 100.);

  //GaussExp variables
  RooRealVar* GaussExpMean  = new RooRealVar("GaussExpMean","GaussExpMean",3.0, 3.05,3.1);
  RooRealVar* GaussExpSigma  = new RooRealVar("GaussExpSigma","GaussExpSigma",0.08,0.07,0.1); //ok result is 0.07 to 3
  RooRealVar* GaussExpAlpha  = new RooRealVar("GaussExpAlpha","GaussExpAlpha", 2, 0.1, 100.);

  //Below 10 Mean: 3.0,3.05,3.1 Sigma:0.08,0.07,0.1 Alpha: 2,0.1,100
  //10to20 Mean:   3.0,3.1,3.1 Sigma:0.08,0.045,0.1 Alpha: 5, 7, 100.

  /*cbRes->setVal(0.07);
  cbRes->setConstant(kTRUE);
  cbTail->setVal(0.2);
  cbTail->setConstant(kTRUE);
  cbN->setVal(99.985);
  cbN->setConstant(kTRUE);*/

  //polynom times exp_pdf
  RooRealVar Alpha("Alpha", "Alpha", 2, 0., 10.);
  RooRealVar Beta("Beta", "Beta", 2, 0.,10. );

  //polynom
  RooRealVar p1("p1","coeff #1", 100, -1000., 1000.);
  RooRealVar p2("p2","coeff #2", 10, -100., 100.);
  RooRealVar p3("p3","coeff #3", 10, -1000., 1000.);
  RooRealVar p4("p4","coeff #4", 10, -1000., 1000.);
  RooRealVar p5("p5","coeff #5", 10, -1000., 1000.);

  //LogTimesExp
  RooRealVar par1("par1","par1", 1, 0.1, 10.);
  RooRealVar par2("par2","par2", 1, 0., 5.);

  //Gaussian
  RooRealVar mean("mean","mean",3.05,2.9,3.2);
  RooRealVar sigma("sigma","sigma",0.1,0.07,0.1);

  //dcbMean->setVal(3.05); // Fix JPsi mass
  //dcbMean->setConstant(kTRUE); //Fix Jpsi mass
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
  RooGenericPdf* PolyTimesExp = new RooGenericPdf("PolyExp","PolyExp","TMath::Power(InvMass,Alpha)*TMath::Exp(-1.*Beta*InvMass)",RooArgList(InvMass,Alpha,Beta));
  RooGenericPdf* LogTimesExp = new RooGenericPdf("LogTimesExp","LogTimesExp","TMath::Log(par1*InvMass)*TMath::Exp(-1*par2*InvMass)",RooArgList(InvMass,par1,par2));
  RooPolynomial Poly("poly","poly",InvMass,RooArgList(p1,p2,p3));
  RooCBShape* cb = new RooCBShape("CrystallBall","CrystallBall",InvMass,*cbMean, *cbTail, *cbRes, *cbN);
  RooGaussian* gauss= new RooGaussian("Gaussian","Gaussian",InvMass, mean,sigma);
  RooGaussExp* GaussExp = new RooGaussExp("GaussExp","GaussExp",InvMass, *GaussExpMean,*GaussExpSigma,*GaussExpAlpha);

  RooArgList shapes, yields;
  //shapes.add(exp_pdf);
  //shapes.add(*PolyTimesExp);
  shapes.add(Poly);
  //shapes.add(*dcb);
  //shapes.add(*LogTimesExp);
  //shapes.add(*cb);
  //shapes.add(*gauss);
  shapes.add(*GaussExp);
  yields.add(bkg_yield);
  yields.add(sig_yield);
  RooAddPdf totalPdf("totalPdf", "sum of signal and background PDF's", shapes, yields);

  totalPdf.fitTo(*data_hist, Extended());

  RooPlot* massframe = InvMass.frame();
  data_hist->plotOn(massframe);
  totalPdf.plotOn(massframe, LineColor(kBlue));
  //totalPdf.plotOn(massframe,Components(*PolyTimesExp), LineColor(kRed));
  //totalPdf.plotOn(massframe,Components(exp_pdf), LineColor(kRed));
  totalPdf.plotOn(massframe,Components(Poly), LineColor(kRed));
  //totalPdf.plotOn(massframe,Components(*LogTimesExp), LineColor(kRed));
  totalPdf.paramOn(massframe,/*Parameters(sig_yield),*/Layout(0.75,0.99,0.95));

  /*Int_t npar = totalPdf.getParameters(*data_hist)->selectByAttrib("Constant",kFALSE)-> getSize();
  Double_t chi2ndf = massframe->chiSquare(npar);
  std::string numberStr = std::to_string(chi2ndf);
  std::string prefix = "chi2/ndof=";
  const char* numberChar = prefix.append(numberStr).c_str();
  TText *txt = new TText(2.1, 340, numberChar);
  txt->SetTextSize(0.04);
  massframe->addObject(txt);*/

  massframe->Draw();
  //std::cout<<massframe->chiSquare(10)<<std::endl;
  canvas->SaveAs("FaillingFit.pdf");

}
