#include <stdio.h>
#include <memory>
#include <string>
#include "TTree.h"
#include "TFile.h"
#include "TLorentzVector.h"
#include<vector>
#include <math.h>
#include <cmath>
#include "../interface/Def.h"

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "FWCore/Framework/interface/ESHandle.h"

#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/Math/interface/deltaR.h"

#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include <FWCore/ParameterSet/interface/FileInPath.h>

#include "CommonTools/Egamma/interface/ConversionTools.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"

#include "RecoEgamma/EgammaElectronAlgos/interface/ElectronHcalHelper.h"
#include "RecoLocalCalo/HGCalRecAlgos/interface/ClusterTools.h"
#include "RecoEcal/EgammaCoreTools/interface/EcalTools.h"
#include "RecoEgamma/EgammaElectronAlgos/interface/ElectronUtilities.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"
#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/EcalDetId/interface/EBDetId.h"
#include "DataFormats/EcalDetId/interface/EEDetId.h"
#include "DataFormats/HcalDetId/interface/HcalDetId.h"


#include "RecoEgamma/ElectronIdentification/interface/ElectronMVAEstimator.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrackFwd.h"



class RecoSkimAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources> {
public:
  explicit RecoSkimAnalyzer(const edm::ParameterSet&);
  ~RecoSkimAnalyzer();

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  void beginJob() override;
  void analyze(const edm::Event&, const edm::EventSetup&) override;
  void endJob() override;

  edm::EDGetToken electronsToken_;
  edm::EDGetTokenT<reco::TrackCollection> TrackToken;
  edm::EDGetTokenT<edm::ValueMap<float> > mvaV2IsoValuesMapToken;
  edm::EDGetTokenT<edm::ValueMap<float> > mvaV2NoIsoValuesMapToken;
  edm::EDGetTokenT<reco::VertexCollection> SecondaryVertexToken;

  edm::InputTag Electron_;
  edm::InputTag TrackInputTag_;
  edm::InputTag mvaV2IsoValuesMapInputTag_;
  edm::InputTag mvaV2NoIsoValuesMapInputTag_;
  edm::InputTag SecondaryVertexInputTag_;

  TTree *reco_tree;
  std::vector<float> ele_pt,ele_eta,ele_phi,ElectronMVAEstimatorRun2Fall17IsoV2Values,
  ElectronMVAEstimatorRun2Fall17NoIsoV2Values,ele_z,ele_x,ele_y,SVchi2,SVndof, SVpt,SVeta,SVphi,SVcharge,SVmass,SVx,SVy,SVz;
  std::vector<int> ele_charge,EleSVtrkref,EleTrkIsNoNull,sv_tr_size,SVTrEleRef;
  int numele,numGoodSV;
  int evtnum;

};

RecoSkimAnalyzer::RecoSkimAnalyzer(const edm::ParameterSet& iConfig):
Electron_(iConfig.getUntrackedParameter<edm::InputTag>("Electron")),
TrackInputTag_(iConfig.getUntrackedParameter<edm::InputTag>("Track")),
mvaV2IsoValuesMapInputTag_(iConfig.getUntrackedParameter<edm::InputTag>("mvaV2IsoValuesMap")),
mvaV2NoIsoValuesMapInputTag_(iConfig.getUntrackedParameter<edm::InputTag>("mvaV2NoIsoValuesMap")),
SecondaryVertexInputTag_(iConfig.getUntrackedParameter<edm::InputTag>("secondaryVertices"))
{
  electronsToken_    = consumes<edm::View<reco::GsfElectron> >(Electron_);
  TrackToken = consumes<reco::TrackCollection>(TrackInputTag_);
  mvaV2IsoValuesMapToken = consumes<edm::ValueMap<float>>(mvaV2IsoValuesMapInputTag_);
  mvaV2NoIsoValuesMapToken = consumes<edm::ValueMap<float>>(mvaV2NoIsoValuesMapInputTag_);
  SecondaryVertexToken = consumes<reco::VertexCollection>(SecondaryVertexInputTag_);

  edm::Service<TFileService> fs;
  reco_tree = fs->make<TTree>("Events", "Events");

  reco_tree->Branch("numele",&numele);
  reco_tree->Branch("ele_eta",&ele_eta);
  reco_tree->Branch("ele_phi",&ele_phi);
  reco_tree->Branch("ele_pt",&ele_pt);
  reco_tree->Branch("ele_charge",&ele_charge);
  reco_tree->Branch("EleTrkIsNoNull",&EleTrkIsNoNull);
  reco_tree->Branch("ele_x",&ele_x);
  reco_tree->Branch("ele_y",&ele_y);
  reco_tree->Branch("ele_z",&ele_z);

  //V2 MVA scores
  reco_tree->Branch("ElectronMVAEstimatorRun2Fall17IsoV2Values",&ElectronMVAEstimatorRun2Fall17IsoV2Values);
  reco_tree->Branch("ElectronMVAEstimatorRun2Fall17NoIsoV2Values",&ElectronMVAEstimatorRun2Fall17NoIsoV2Values);


  //Secondary Vertex
  reco_tree->Branch("SVchi2",&SVchi2);
  reco_tree->Branch("SVndof",&SVndof);
  reco_tree->Branch("numGoodSV",&numGoodSV);
  reco_tree->Branch("SVpt",&SVpt);
  reco_tree->Branch("SVphi",&SVphi);
  reco_tree->Branch("SVeta",&SVeta);
  reco_tree->Branch("EleSVtrkref",&EleSVtrkref);
  reco_tree->Branch("SVTrEleRef",&SVTrEleRef);
  reco_tree->Branch("SVcharge",&SVcharge);
  reco_tree->Branch("SVmass",&SVmass);
  reco_tree->Branch("SVx",&SVx);
  reco_tree->Branch("SVy",&SVy);
  reco_tree->Branch("SVz",&SVz);

  //Event number
  reco_tree->Branch("evtnum",&evtnum);

}

RecoSkimAnalyzer::~RecoSkimAnalyzer() {

}


void RecoSkimAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  using namespace edm;

  edm::Handle<edm::View<reco::GsfElectron> > electrons;
  iEvent.getByToken(electronsToken_, electrons);

  edm::Handle<reco::TrackCollection> TrackHandle;
  iEvent.getByToken(TrackToken, TrackHandle);

  edm::Handle<edm::ValueMap<float> > mvaV2NoIsoValues;
  iEvent.getByToken(mvaV2NoIsoValuesMapToken,mvaV2NoIsoValues);

  edm::Handle<edm::ValueMap<float> > mvaV2IsoValues;
  iEvent.getByToken(mvaV2IsoValuesMapToken,mvaV2IsoValues);

  edm::Handle<reco::VertexCollection> secondaryVertexHandle;
  iEvent.getByToken(SecondaryVertexToken,secondaryVertexHandle);

  numele=0;
  evtnum=0;

  ele_eta.clear();
  ele_phi.clear();
  ele_pt.clear();
  ele_charge.clear();
  ele_x.clear();
  ele_y.clear();
  ele_z.clear();


  ElectronMVAEstimatorRun2Fall17IsoV2Values.clear();
  ElectronMVAEstimatorRun2Fall17NoIsoV2Values.clear();
  EleTrkIsNoNull.clear();

  SVchi2.clear();
  SVndof.clear();
  SVpt.clear();
  SVphi.clear();
  SVeta.clear();
  EleSVtrkref.clear();
  SVTrEleRef.clear();
  SVcharge.clear();
  SVmass.clear();
  numGoodSV=0;
  SVx.clear();
  SVy.clear();
  SVz.clear();

  //helper variables
  int ChargeHelper;
  TLorentzVector P_sum,P_tmp;

  const reco::TrackCollection* tkColl = TrackHandle.product();


  for(auto sv = secondaryVertexHandle->begin(); sv!= secondaryVertexHandle->end(); ++sv)
  {
      ChargeHelper=0.;
      P_sum.SetPtEtaPhiM(0.,0.,0.,0.);
      int SVtrackSize = 0;
      for(long unsigned int r=0;r<sv->tracksSize();r++) //calcaulte variables needed to select KJPsi events
      {
        if(sv->trackRefAt(r)->pt()>0.7)
        {
          ChargeHelper+=sv->trackRefAt(r)->charge();
          P_tmp.SetPtEtaPhiM(sv->trackRefAt(r)->pt(),sv->trackRefAt(r)->eta(),sv->trackRefAt(r)->phi(),0.);
          P_sum+=P_tmp;
          SVtrackSize++;
        }
      }
      if(P_sum.M()>4.5 && SVtrackSize==3 && std::abs(ChargeHelper)==1) //fill only if satisfied
      {
        for(long unsigned int r=0;r<sv->tracksSize();r++)
        {
          if(sv->trackRefAt(r)->pt()>0.7)
          {
            SVpt.push_back(sv->trackRefAt(r)->pt());
            SVeta.push_back(sv->trackRefAt(r)->eta());
            SVphi.push_back(sv->trackRefAt(r)->phi());
            SVcharge.push_back(sv->trackRefAt(r)->charge());
          }
        }
        SVchi2.push_back(sv->chi2());
        SVndof.push_back(sv->ndof());
        SVx.push_back(sv->x());
        SVy.push_back(sv->y());
        SVz.push_back(sv->z());
        SVmass.push_back(P_sum.M());
        numGoodSV++;
      }
  }


  for(auto el = electrons->begin(); el != electrons->end(); ++el)
  {

    if(el->closestCtfTrackRef().isNonnull() && el->pt()>5.)  //what if the reference to track does not exist ?
    {
      auto ClosestTrackRef = tkColl->begin() + el->closestCtfTrackRef().index();
      int FoundSVmatch=0;
      for(long unsigned int b=0; b<SVpt.size(); b++)
      {
        if((float)ClosestTrackRef->pt()==SVpt.at(b) && (float)ClosestTrackRef->eta()==SVeta.at(b) && (float)ClosestTrackRef->phi()==SVphi.at(b))
        {
          EleSVtrkref.push_back(b);
          FoundSVmatch=1;
        }
      }
      if(FoundSVmatch==0)
      EleSVtrkref.push_back(99);
    }
    else if(el->pt()>5.)
    EleSVtrkref.push_back(99);

    if(el->pt()>5.)
    {
      const auto iter = electrons->ptrAt(numele);
      ele_pt.push_back(el->pt());
      ele_eta.push_back(el->eta());
      ele_phi.push_back(el->phi());
      ele_charge.push_back(el->charge());
      ElectronMVAEstimatorRun2Fall17NoIsoV2Values.push_back((*mvaV2NoIsoValues)[iter]);
      ElectronMVAEstimatorRun2Fall17IsoV2Values.push_back((*mvaV2IsoValues)[iter]);
      EleTrkIsNoNull.push_back(el->closestCtfTrackRef().isNonnull());
      ele_x.push_back(el->vx());
      ele_y.push_back(el->vy());
      ele_z.push_back(el->vz());
      numele++;
    }
  }

  for(long unsigned int b=0; b<SVpt.size();b++)
  {
    int RefSVtoEl=99;
    for(long unsigned c=0; c<EleSVtrkref.size();c++)
    {
      if(EleSVtrkref.at(c)==(int)b)
      RefSVtoEl=c;
    }
    SVTrEleRef.push_back(RefSVtoEl);
  }

  evtnum=iEvent.eventAuxiliary().event();

  if((numele==1 || numele==2) && numGoodSV>=1) //add charge condition ?
  {
    //std::cout<<iEvent.eventAuxiliary().event()<<std::endl;
    reco_tree->Fill();
  }

}


void RecoSkimAnalyzer::beginJob() {

}

void RecoSkimAnalyzer::endJob() {

}

void RecoSkimAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);

}

DEFINE_FWK_MODULE(RecoSkimAnalyzer);
