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
#include "DataFormats/Math/interface/deltaR.h"
#include "TMath.h"

float InvariantMass(float pt1,float eta1,float phi1,float pt2,float eta2,float phi2)
{
   TLorentzVector P,P1,P2;
   P1.SetPtEtaPhiM(pt1,eta1,phi1,0);
   P2.SetPtEtaPhiM(pt2,eta2,phi2,0);
   P=P1+P2;
   return P.M();
}
double distancexy(float x,float y)
{
  return sqrt(x*x+y*y);
}

void NewSelectionPlotter()
{

  TFile *f = new TFile("v5_BC.root");

  TTree *tree = (TTree*)f->Get("reco/Events");

  std::vector<float> *ele_pt=0,*ele_eta=0,*ele_phi=0,*tr_pt=0,*tr_eta=0,*tr_phi=0,
  *ElectronMVAEstimatorRun2Fall17NoIsoV2Values=0,*SVpt=0,*SVeta=0,*SVphi=0,*SVmass=0,*SVdistanceToMaxPtPV=0,
  *SVx=0,*SVy=0,*SVz=0,*SVndof=0,*SVchi2=0;
  std::vector<int> *EleTRKref=0,*SVsize=0,*EleSVtrkref=0,*SVcharge=0,*SVchargeSum=0,*ele_charge=0,*tr_charge=0,
  *SVTrEleRef=0;
  int numele,numtr,numsv;
  float PVmaxPtX,PVmaxPtY,PVmaxPtZ;

  //help variables
  std::vector<float> TmpDeltaR;
  float MinDeltaR;

  tree->SetBranchAddress("ele_pt",&ele_pt);
  tree->SetBranchAddress("ele_charge",&ele_charge);
  tree->SetBranchAddress("tr_charge",&tr_charge);
  tree->SetBranchAddress("ele_eta",&ele_eta);
  tree->SetBranchAddress("ele_phi",&ele_phi);
  tree->SetBranchAddress("tr_pt",&tr_pt);
  tree->SetBranchAddress("tr_eta",&tr_eta);
  tree->SetBranchAddress("tr_phi",&tr_phi);
  tree->SetBranchAddress("numele",&numele);
  tree->SetBranchAddress("numtr",&numtr);
  tree->SetBranchAddress("EleTRKref",&EleTRKref);
  tree->SetBranchAddress("ElectronMVAEstimatorRun2Fall17NoIsoV2Values",&ElectronMVAEstimatorRun2Fall17NoIsoV2Values);
  tree->SetBranchAddress("numsv",&numsv);
  tree->SetBranchAddress("SVsize",&SVsize);
  tree->SetBranchAddress("SVcharge",&SVcharge);
  tree->SetBranchAddress("SVchargeSum",&SVchargeSum);
  tree->SetBranchAddress("SVmass",&SVmass);
  tree->SetBranchAddress("SVpt",&SVpt);
  tree->SetBranchAddress("SVeta",&SVeta);
  tree->SetBranchAddress("SVphi",&SVphi);
  tree->SetBranchAddress("EleSVtrkref",&EleSVtrkref);
  tree->SetBranchAddress("SVTrEleRef",&SVTrEleRef);
  tree->SetBranchAddress("SVdistanceToMaxPtPV",&SVdistanceToMaxPtPV);
  tree->SetBranchAddress("SVx",&SVx);
  tree->SetBranchAddress("SVy",&SVy);
  tree->SetBranchAddress("SVz",&SVz);
  tree->SetBranchAddress("PVmaxPtX",&PVmaxPtX);
  tree->SetBranchAddress("PVmaxPtY",&PVmaxPtY);
  tree->SetBranchAddress("PVmaxPtZ",&PVmaxPtZ);
  tree->SetBranchAddress("SVndof",&SVndof);
  tree->SetBranchAddress("SVchi2",&SVchi2);


  TH1F *h_PassSV3Mass4p5infdR = new TH1F("PassSV3Mass4p5infdR","PassSV3Mass4p5infdR",100,2,7);
  TH1F *h_FailSV3Mass4p5infdR = new TH1F("FailSV3Mass4p5infdR","FailSV3Mass4p5infdR",100,2,7);
  TH1F *h_PassSV3Mass4p5infSVprob = new TH1F("PassSV3Mass4p5infSVprob","PassSV3Mass4p5infSVprob",100,2,7);
  TH1F *h_FailSV3Mass4p5infSVprob = new TH1F("FailSV3Mass4p5infSVprob","FailSV3Mass4p5infSVprob",100,2,7);
  TH1F *h_PassSV3Mass4p5infSVPtSum = new TH1F("PassSV3Mass4p5infSVPtSum","PassSV3Mass4p5infSVPtSum",100,2,7);
  TH1F *h_FailSV3Mass4p5infSVPtSum = new TH1F("FailSV3Mass4p5infSVPtSum","FailSV3Mass4p5infdRSVPtSum",100,2,7);
  TH1F *h_PassSV3Mass4p5infMinDeltaR = new TH1F("PassSV3Mass4p5infMinDeltaR","PassSV3Mass4p5infMinDeltaR",100,2,7);
  TH1F *h_FailSV3Mass4p5infMinDeltaR = new TH1F("FailSV3Mass4p5infMinDeltaR","FailSV3Mass4p5infMinDeltaR",100,2,7);
  TH1F *h_PassSV3Mass4p5infAllCuts = new TH1F("PassSV3Mass4p5infAllCuts","PassSV3Mass4p5infAllCuts",100,2,7);
  TH1F *h_FailSV3Mass4p5infAllCuts = new TH1F("FailSV3Mass4p5infAllCuts","FailSV3Mass4p5infAllCuts",100,2,7);

  TFile* output = TFile::Open("HistoOutputNewSelectionBC.root","RECREATE");

  Int_t nentries=(Int_t)tree->GetEntries();

  for(Int_t i=0;i<nentries;i++)
  {

    tree->GetEntry(i);

    for(int j=0;j<numele;j++)
    {
      if(EleSVtrkref->at(j)!=99)
      {
        int helper_up=0;
        int iter = 0;
        while(helper_up<=EleSVtrkref->at(j))
        {
          helper_up+=SVsize->at(iter);
          iter++;
        }

        int helper_dn=helper_up-SVsize->at(iter-1);

        int ChargeSum=SVchargeSum->at(iter-1);
        int SVTrackSize=SVsize->at(iter-1);
        float mass=SVmass->at(iter-1);
        float distance3D=SVdistanceToMaxPtPV->at(iter-1);
        float distance2D=sqrt((SVx->at(iter-1)-PVmaxPtX)*(SVx->at(iter-1)-PVmaxPtX)+(SVy->at(iter-1)-PVmaxPtY)*(SVy->at(iter-1)-PVmaxPtY));


         //Ele and track position are on same place so I can use SV3eta and SV3phi
         TmpDeltaR.clear();
         for(int m=helper_dn;m<helper_up;m++)
         {
           for(int n=helper_dn;n<helper_up;n++)
           {
             if(n!=m)
             {
               TmpDeltaR.push_back(reco::deltaR2(SVeta->at(m),SVphi->at(m),SVeta->at(n),SVphi->at(n)));
             }
           }
         }
         MinDeltaR=*min_element(TmpDeltaR.begin(), TmpDeltaR.end());

        float SVPtSum=0;
        for(int l=helper_dn;l<helper_up;l++)
        {
          if(SVTrEleRef->at(l)==j)
          SVPtSum+=ele_pt->at(j);
          else
          SVPtSum+=SVpt->at(l);
        }


        for(long unsigned int k=0;k<SVpt->size();k++)
        {
          if(SVTrEleRef->at(k)!=99 && SVTrEleRef->at(k)!=j && ele_charge->at(j)==-1*SVcharge->at(k) &&
              ele_pt->at(j)>7. && ElectronMVAEstimatorRun2Fall17NoIsoV2Values->at(j)>0.95 && SVpt->at(k)>5. &&
              k>=helper_dn && k<helper_up && mass>4.5 && SVTrackSize==3 && std::abs(ChargeSum)==1 && reco::deltaR2(ele_eta->at(j),ele_phi->at(j),SVeta->at(k),SVphi->at(k))>0.03)
              h_PassSV3Mass4p5infdR->Fill(InvariantMass(ele_pt->at(j),ele_eta->at(j),ele_phi->at(j),SVpt->at(k),SVeta->at(k),SVphi->at(k)));

          if(SVTrEleRef->at(k)==99 && ele_charge->at(j)==-1*SVcharge->at(k) &&
              ele_pt->at(j)>7. && ElectronMVAEstimatorRun2Fall17NoIsoV2Values->at(j)>0.95 && SVpt->at(k)>5. &&
              k>=helper_dn && k<=helper_up && mass>4.5 && SVTrackSize==3 && std::abs(ChargeSum)==1 && reco::deltaR2(ele_eta->at(j),ele_phi->at(j),SVeta->at(k),SVphi->at(k))>0.03)
              h_FailSV3Mass4p5infdR->Fill(InvariantMass(ele_pt->at(j),ele_eta->at(j),ele_phi->at(j),SVpt->at(k),SVeta->at(k),SVphi->at(k)));

          if(SVTrEleRef->at(k)!=99 && SVTrEleRef->at(k)!=j && ele_charge->at(j)==-1*SVcharge->at(k) &&
              ele_pt->at(j)>7. && ElectronMVAEstimatorRun2Fall17NoIsoV2Values->at(j)>0.95 && SVpt->at(k)>5. &&
              k>=helper_dn && k<helper_up && mass>4.5 && SVTrackSize==3 && std::abs(ChargeSum)==1 && TMath::Prob(SVchi2->at(iter-1),SVndof->at(iter-1))>0.001)
              h_PassSV3Mass4p5infSVprob->Fill(InvariantMass(ele_pt->at(j),ele_eta->at(j),ele_phi->at(j),SVpt->at(k),SVeta->at(k),SVphi->at(k)));

          if(SVTrEleRef->at(k)==99 && ele_charge->at(j)==-1*SVcharge->at(k) &&
              ele_pt->at(j)>7. && ElectronMVAEstimatorRun2Fall17NoIsoV2Values->at(j)>0.95 && SVpt->at(k)>5. &&
              k>=helper_dn && k<=helper_up && mass>4.5 && SVTrackSize==3 && std::abs(ChargeSum)==1 && TMath::Prob(SVchi2->at(iter-1),SVndof->at(iter-1))>0.001)
              h_FailSV3Mass4p5infSVprob->Fill(InvariantMass(ele_pt->at(j),ele_eta->at(j),ele_phi->at(j),SVpt->at(k),SVeta->at(k),SVphi->at(k)));

          if(SVTrEleRef->at(k)!=99 && SVTrEleRef->at(k)!=j && ele_charge->at(j)==-1*SVcharge->at(k) &&
              ele_pt->at(j)>7. && ElectronMVAEstimatorRun2Fall17NoIsoV2Values->at(j)>0.95 && SVpt->at(k)>5. &&
              k>=helper_dn && k<helper_up && mass>4.5 && SVTrackSize==3 && std::abs(ChargeSum)==1 && SVPtSum>3.)
              h_PassSV3Mass4p5infSVPtSum->Fill(InvariantMass(ele_pt->at(j),ele_eta->at(j),ele_phi->at(j),SVpt->at(k),SVeta->at(k),SVphi->at(k)));

          if(SVTrEleRef->at(k)==99 && ele_charge->at(j)==-1*SVcharge->at(k) &&
              ele_pt->at(j)>7. && ElectronMVAEstimatorRun2Fall17NoIsoV2Values->at(j)>0.95 && SVpt->at(k)>5. &&
              k>=helper_dn && k<=helper_up && mass>4.5 && SVTrackSize==3 && std::abs(ChargeSum)==1 && SVPtSum>3.)
              h_FailSV3Mass4p5infSVPtSum->Fill(InvariantMass(ele_pt->at(j),ele_eta->at(j),ele_phi->at(j),SVpt->at(k),SVeta->at(k),SVphi->at(k)));

          if(SVTrEleRef->at(k)!=99 && SVTrEleRef->at(k)!=j && ele_charge->at(j)==-1*SVcharge->at(k) &&
              ele_pt->at(j)>7. && ElectronMVAEstimatorRun2Fall17NoIsoV2Values->at(j)>0.95 && SVpt->at(k)>5. &&
              k>=helper_dn && k<helper_up && mass>4.5 && SVTrackSize==3 && std::abs(ChargeSum)==1 && MinDeltaR>0.03)
              h_PassSV3Mass4p5infMinDeltaR->Fill(InvariantMass(ele_pt->at(j),ele_eta->at(j),ele_phi->at(j),SVpt->at(k),SVeta->at(k),SVphi->at(k)));

          if(SVTrEleRef->at(k)==99 && ele_charge->at(j)==-1*SVcharge->at(k) &&
              ele_pt->at(j)>7. && ElectronMVAEstimatorRun2Fall17NoIsoV2Values->at(j)>0.95 && SVpt->at(k)>5. &&
              k>=helper_dn && k<=helper_up && mass>4.5 && SVTrackSize==3 && std::abs(ChargeSum)==1 && MinDeltaR>0.03)
              h_FailSV3Mass4p5infMinDeltaR->Fill(InvariantMass(ele_pt->at(j),ele_eta->at(j),ele_phi->at(j),SVpt->at(k),SVeta->at(k),SVphi->at(k)));

          if(SVTrEleRef->at(k)!=99 && SVTrEleRef->at(k)!=j && ele_charge->at(j)==-1*SVcharge->at(k) &&
              ele_pt->at(j)>7. && ElectronMVAEstimatorRun2Fall17NoIsoV2Values->at(j)>0.95 && SVpt->at(k)>5. &&
              k>=helper_dn && k<helper_up && mass>4.5 && SVTrackSize==3 && std::abs(ChargeSum)==1 && MinDeltaR>0.03 &&
              reco::deltaR2(ele_eta->at(j),ele_phi->at(j),SVeta->at(k),SVphi->at(k))>0.03 &&
              TMath::Prob(SVchi2->at(iter-1),SVndof->at(iter-1))>0.001 && SVPtSum>3. && MinDeltaR>0.03)
              h_PassSV3Mass4p5infAllCuts->Fill(InvariantMass(ele_pt->at(j),ele_eta->at(j),ele_phi->at(j),SVpt->at(k),SVeta->at(k),SVphi->at(k)));

          if(SVTrEleRef->at(k)==99 && ele_charge->at(j)==-1*SVcharge->at(k) &&
              ele_pt->at(j)>7. && ElectronMVAEstimatorRun2Fall17NoIsoV2Values->at(j)>0.95 && SVpt->at(k)>5. &&
              k>=helper_dn && k<=helper_up && mass>4.5 && SVTrackSize==3 && std::abs(ChargeSum)==1 && MinDeltaR>0.03 &&
              reco::deltaR2(ele_eta->at(j),ele_phi->at(j),SVeta->at(k),SVphi->at(k))>0.03 &&
              TMath::Prob(SVchi2->at(iter-1),SVndof->at(iter-1))>0.001 && SVPtSum>3. && MinDeltaR>0.03)
              h_FailSV3Mass4p5infAllCuts->Fill(InvariantMass(ele_pt->at(j),ele_eta->at(j),ele_phi->at(j),SVpt->at(k),SVeta->at(k),SVphi->at(k)));

        }
      }
    }

   }

   h_PassSV3Mass4p5infdR->Write();
   h_FailSV3Mass4p5infdR->Write();
   h_PassSV3Mass4p5infSVprob->Write();
   h_FailSV3Mass4p5infSVprob->Write();
   h_PassSV3Mass4p5infSVPtSum->Write();
   h_FailSV3Mass4p5infSVPtSum->Write();
   h_PassSV3Mass4p5infMinDeltaR->Write();
   h_FailSV3Mass4p5infMinDeltaR->Write();
   h_PassSV3Mass4p5infAllCuts->Write();
   h_FailSV3Mass4p5infAllCuts->Write();

   output->Close();

}
