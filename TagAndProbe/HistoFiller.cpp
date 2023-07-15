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

void HistoFiller()
{

  TFile *f = new TFile("skim123456ABCD.root");

  TTree *tree = (TTree*)f->Get("reco/Events");

  std::vector<float> *ele_pt=0,*ele_eta=0,*ele_phi=0,
  *ElectronMVAEstimatorRun2Fall17NoIsoV2Values=0,*SVpt=0,*SVeta=0,*SVphi=0,*SVmass=0,
  *SVx=0,*SVy=0,*SVz=0,*SVndof=0,*SVchi2=0;
  std::vector<int> *EleSVtrkref=0,*SVcharge=0,*ele_charge=0, *SVTrEleRef=0;
  int numele,numsv,numGoodSV;

  //help variables
  std::vector<float> TmpDeltaR;
  float MinDeltaR;

  tree->SetBranchAddress("ele_pt",&ele_pt);
  tree->SetBranchAddress("ele_charge",&ele_charge);
  tree->SetBranchAddress("ele_eta",&ele_eta);
  tree->SetBranchAddress("ele_phi",&ele_phi);
  tree->SetBranchAddress("numele",&numele);
  tree->SetBranchAddress("ElectronMVAEstimatorRun2Fall17NoIsoV2Values",&ElectronMVAEstimatorRun2Fall17NoIsoV2Values);
  tree->SetBranchAddress("numGoodSV",&numGoodSV);
  tree->SetBranchAddress("SVcharge",&SVcharge);
  tree->SetBranchAddress("SVmass",&SVmass);
  tree->SetBranchAddress("SVpt",&SVpt);
  tree->SetBranchAddress("SVeta",&SVeta);
  tree->SetBranchAddress("SVphi",&SVphi);
  tree->SetBranchAddress("EleSVtrkref",&EleSVtrkref);
  tree->SetBranchAddress("SVTrEleRef",&SVTrEleRef);
  tree->SetBranchAddress("SVx",&SVx);
  tree->SetBranchAddress("SVy",&SVy);
  tree->SetBranchAddress("SVz",&SVz);
  tree->SetBranchAddress("SVndof",&SVndof);
  tree->SetBranchAddress("SVchi2",&SVchi2);

  TH1F *h_PassSV3Mass4p5infAllCuts = new TH1F("PassSV3Mass4p5infAllCuts","PassSV3Mass4p5infAllCuts",50,2,7);
  TH1F *h_FailSV3Mass4p5infAllCuts = new TH1F("FailSV3Mass4p5infAllCuts","FailSV3Mass4p5infAllCuts",50,2,7);

  TFile* output = TFile::Open("HistoOutput.root","RECREATE");

  Int_t nentries=(Int_t)tree->GetEntries();

  for(Int_t i=0;i<nentries;i++)
  {

    tree->GetEntry(i);

    for(int j=0;j<numele;j++)
    {
      if(EleSVtrkref->at(j)!=99)
      {

        int helper_dn=EleSVtrkref->at(j)/3*3;
        int helper_up=helper_dn+3;
        int iter=EleSVtrkref->at(j)/3;

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
              k>=helper_dn && k<helper_up && reco::deltaR2(ele_eta->at(j),ele_phi->at(j),SVeta->at(k),SVphi->at(k))>0.03 &&
              TMath::Prob(SVchi2->at(iter),SVndof->at(iter))>0.001 && SVPtSum>3. && MinDeltaR>0.03)
              h_PassSV3Mass4p5infAllCuts->Fill(InvariantMass(ele_pt->at(j),ele_eta->at(j),ele_phi->at(j),SVpt->at(k),SVeta->at(k),SVphi->at(k)));

          if(SVTrEleRef->at(k)==99 && ele_charge->at(j)==-1*SVcharge->at(k) &&
              ele_pt->at(j)>7. && ElectronMVAEstimatorRun2Fall17NoIsoV2Values->at(j)>0.95 && SVpt->at(k)>5. &&
              k>=helper_dn && k<=helper_up  && reco::deltaR2(ele_eta->at(j),ele_phi->at(j),SVeta->at(k),SVphi->at(k))>0.03 &&
              TMath::Prob(SVchi2->at(iter),SVndof->at(iter))>0.001 && SVPtSum>3. && MinDeltaR>0.03)
              h_FailSV3Mass4p5infAllCuts->Fill(InvariantMass(ele_pt->at(j),ele_eta->at(j),ele_phi->at(j),SVpt->at(k),SVeta->at(k),SVphi->at(k)));
        }
      }
    }

   }

   h_PassSV3Mass4p5infAllCuts->Write();
   h_FailSV3Mass4p5infAllCuts->Write();


   output->Close();

}
