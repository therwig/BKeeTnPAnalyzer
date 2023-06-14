#include <iostream>
#include <vector>
#include "TROOT.h"
#include "TChain.h"
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TChain.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TH1F.h"
#include "TLorentzVector.h"
#include "TRandom.h"

void HistoReaderNewSelection()
{
  TFile* input1 = TFile::Open("HistoOutputNewSelectionA.root");
  TFile* input2 = TFile::Open("HistoOutputNewSelectionBC.root");

  TH1F* Fail4p5infdR1 = (TH1F*)input1->Get("FailSV3Mass4p5infdR");
  TH1F* Pass4p5infdR1 = (TH1F*)input1->Get("PassSV3Mass4p5infdR");
  TH1F* Fail4p5infSVprob1 = (TH1F*)input1->Get("FailSV3Mass4p5infSVprob");
  TH1F* Pass4p5infSVprob1 = (TH1F*)input1->Get("PassSV3Mass4p5infSVprob");
  TH1F* Fail4p5infSVPtSum1 = (TH1F*)input1->Get("FailSV3Mass4p5infSVPtSum");
  TH1F* Pass4p5infSVPtSum1 = (TH1F*)input1->Get("PassSV3Mass4p5infSVPtSum");
  TH1F* Fail4p5infMinDeltaR1 = (TH1F*)input1->Get("FailSV3Mass4p5infMinDeltaR");
  TH1F* Pass4p5infMinDeltaR1 = (TH1F*)input1->Get("PassSV3Mass4p5infMinDeltaR");
  TH1F* Fail4p5infAllCuts1 = (TH1F*)input1->Get("FailSV3Mass4p5infAllCuts");
  TH1F* Pass4p5infAllCuts1 = (TH1F*)input1->Get("PassSV3Mass4p5infAllCuts");

  TH1F* Fail4p5infdR2 = (TH1F*)input2->Get("FailSV3Mass4p5infdR");
  TH1F* Pass4p5infdR2 = (TH1F*)input2->Get("PassSV3Mass4p5infdR");
  TH1F* Fail4p5infSVprob2 = (TH1F*)input2->Get("FailSV3Mass4p5infSVprob");
  TH1F* Pass4p5infSVprob2 = (TH1F*)input2->Get("PassSV3Mass4p5infSVprob");
  TH1F* Fail4p5infSVPtSum2 = (TH1F*)input2->Get("FailSV3Mass4p5infSVPtSum");
  TH1F* Pass4p5infSVPtSum2 = (TH1F*)input2->Get("PassSV3Mass4p5infSVPtSum");
  TH1F* Fail4p5infMinDeltaR2 = (TH1F*)input2->Get("FailSV3Mass4p5infMinDeltaR");
  TH1F* Pass4p5infMinDeltaR2 = (TH1F*)input2->Get("PassSV3Mass4p5infMinDeltaR");
  TH1F* Fail4p5infAllCuts2 = (TH1F*)input2->Get("FailSV3Mass4p5infAllCuts");
  TH1F* Pass4p5infAllCuts2 = (TH1F*)input2->Get("PassSV3Mass4p5infAllCuts");

  Fail4p5infdR1->Add(Fail4p5infdR2);
  Pass4p5infdR1->Add(Pass4p5infdR2);
  Fail4p5infSVprob1->Add(Fail4p5infSVprob2);
  Pass4p5infSVprob1->Add(Pass4p5infSVprob2);
  Fail4p5infSVPtSum1->Add(Fail4p5infSVPtSum2);
  Pass4p5infSVPtSum1->Add(Pass4p5infSVPtSum2);
  Fail4p5infMinDeltaR1->Add(Fail4p5infMinDeltaR2);
  Pass4p5infMinDeltaR1->Add(Pass4p5infMinDeltaR2);
  Fail4p5infAllCuts1->Add(Fail4p5infAllCuts2);
  Pass4p5infAllCuts1->Add(Pass4p5infAllCuts2);

  TH1F *Fail4p5infBindR1 = (TH1F*)Fail4p5infdR1->Rebin(2);
  TH1F *Pass4p5infBindR1 = (TH1F*)Pass4p5infdR1->Rebin(2);

  TH1F *Fail4p5infBinSVprob1 = (TH1F*)Fail4p5infSVprob1->Rebin(2);
  TH1F *Pass4p5infBinSVprob1 = (TH1F*)Pass4p5infSVprob1->Rebin(2);

  TH1F *Fail4p5infBinSVPtSum1 = (TH1F*)Fail4p5infSVPtSum1->Rebin(2);
  TH1F *Pass4p5infBinSVPtSum1 = (TH1F*)Pass4p5infSVPtSum1->Rebin(2);

  TH1F *Fail4p5infBinMinDeltaR1 = (TH1F*)Fail4p5infMinDeltaR1->Rebin(2);
  TH1F *Pass4p5infBinMinDeltaR1 = (TH1F*)Pass4p5infMinDeltaR1->Rebin(2);

  TH1F *Fail4p5infBinAllCuts1 = (TH1F*)Fail4p5infAllCuts1->Rebin(2);
  TH1F *Pass4p5infBinAllCuts1 = (TH1F*)Pass4p5infAllCuts1->Rebin(2);

  TCanvas *canva1 = new TCanvas("cs1","cs1",1600,800);
  TCanvas *canva2 = new TCanvas("cs2","cs2",1600,800);
  TCanvas *canva3 = new TCanvas("cs3","cs3",1600,800);
  TCanvas *canva4 = new TCanvas("cs4","cs4",1600,800);
  TCanvas *canva5 = new TCanvas("cs5","cs5",1600,800);
  TCanvas *canva6 = new TCanvas("cs6","cs6",1600,800);
  TCanvas *canva7 = new TCanvas("cs7","cs7",1600,800);
  TCanvas *canva8 = new TCanvas("cs8","cs8",1600,800);
  TCanvas *canva9 = new TCanvas("cs9","cs9",1600,800);
  TCanvas *canva10 = new TCanvas("cs10","cs10",1600,800);


  canva1->cd();
  Fail4p5infBindR1->Draw();
  canva1->SaveAs("SV3Pt5NewSelection/Fail4p5infBindR.pdf");

  canva2->cd();
  Pass4p5infBindR1->Draw();
  canva2->SaveAs("SV3Pt5NewSelection/Pass4p5infBindR.pdf");

  canva3->cd();
  Fail4p5infBinSVprob1->Draw();
  canva3->SaveAs("SV3Pt5NewSelection/Fail4p5infBinSVprob.pdf");

  canva4->cd();
  Pass4p5infBinSVprob1->Draw();
  canva4->SaveAs("SV3Pt5NewSelection/Pass4p5infBinSVprob.pdf");

  canva5->cd();
  Fail4p5infBinSVPtSum1->Draw();
  canva5->SaveAs("SV3Pt5NewSelection/Fail4p5infBinSVPtSum.pdf");

  canva6->cd();
  Pass4p5infBinSVPtSum1->Draw();
  canva6->SaveAs("SV3Pt5NewSelection/Pass4p5infBinSVPtSum.pdf");

  canva7->cd();
  Fail4p5infBinMinDeltaR1->Draw();
  canva7->SaveAs("SV3Pt5NewSelection/Fail4p5infBinMinDeltaR.pdf");

  canva8->cd();
  Pass4p5infBinMinDeltaR1->Draw();
  canva8->SaveAs("SV3Pt5NewSelection/Pass4p5infBinMinDeltaR.pdf");

  canva9->cd();
  Fail4p5infAllCuts1->Draw();
  canva9->SaveAs("SV3Pt5NewSelection/Fail4p5infAllCuts.pdf");

  canva10->cd();
  Pass4p5infBinAllCuts1->Draw();
  canva10->SaveAs("SV3Pt5NewSelection/Pass4p5infBinAllCuts.pdf");



}
