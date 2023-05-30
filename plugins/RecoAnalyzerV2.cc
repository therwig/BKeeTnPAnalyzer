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

#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "RecoEgamma/EgammaTools/interface/HoECalculator.h"

#include "DataFormats/EcalRecHit/interface/EcalRecHit.h"
#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"
#include "DataFormats/HcalRecHit/interface/HcalRecHitCollections.h"
#include "DataFormats/CaloRecHit/interface/CaloID.h"

#include "RecoEgamma/ElectronIdentification/interface/ElectronMVAEstimator.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrackFwd.h"



class RecoAnalyzerV2 : public edm::one::EDAnalyzer<edm::one::SharedResources> {
public:
  explicit RecoAnalyzerV2(const edm::ParameterSet&);
  ~RecoAnalyzerV2();

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  void beginJob() override;
  void analyze(const edm::Event&, const edm::EventSetup&) override;
  void endJob() override;

  edm::EDGetToken electronsToken_;
  edm::EDGetTokenT<reco::TrackCollection> TrackToken;
  edm::EDGetTokenT<edm::ValueMap<float> > mvaV2IsoValuesMapToken;
  edm::EDGetTokenT<edm::ValueMap<float> > mvaV2NoIsoValuesMapToken;
  edm::EDGetTokenT<reco::GsfTrackCollection> GsfTrackToken;
  edm::EDGetTokenT<reco::VertexCollection> SecondaryVertexToken;
  edm::EDGetTokenT<reco::VertexCollection> PrimaryVertexToken;

  edm::InputTag Electron_;
  edm::InputTag TrackInputTag_;
  edm::InputTag mvaV2IsoValuesMapInputTag_;
  edm::InputTag mvaV2NoIsoValuesMapInputTag_;
  edm::InputTag GsfTrackInputTag_;
  edm::InputTag SecondaryVertexInputTag_;
  edm::InputTag PrimaryVertexInputTag_;

  TTree *reco_tree;
  std::vector<float> ele_pt,ele_eta,ele_phi,trkIsoEle,EleClosestTrackPt,
  tr_pt,tr_eta,tr_phi,ElectronMVAEstimatorRun2Fall17IsoV2Values,ElectronMVAEstimatorRun2Fall17NoIsoV2Values,
  ele_z,ele_x,ele_y,tr_z,tr_x,tr_y,gsf_z,gsf_x,gsf_y, sv_x, sv_y, sv_z,sv_chi2,sv_ndof,gsf_pt,gsf_eta,gsf_phi,
  sv_tr_pt, sv_tr_eta, sv_tr_phi,SV3pt,SV3eta,SV3phi,SV3charge,SV3mass,SVdistanceToMaxPtPV,SV3x,SV3y,SV3z;
  std::vector<int> ele_charge,EleTRKref,EleTrkIsNoNull,ele_gsf_charge,gsf_charge,tr_charge,
  sv_tr_size,EleSV3trkref,SV3TrEleRef;
  int numele, PFnumele,EleCounter,match,numtr,numgsf,numsv,numsv3;
  TLorentzVector P,P0,P1,p,p0,p1;
  float PVmaxPtX,PVmaxPtY,PVmaxPtZ;

  //unsigned long long cachedCaloGeometryID_;
  //edm::ESHandle<CaloGeometry> caloGeometry_;


};

RecoAnalyzerV2::RecoAnalyzerV2(const edm::ParameterSet& iConfig):
Electron_(iConfig.getUntrackedParameter<edm::InputTag>("Electron")),
TrackInputTag_(iConfig.getUntrackedParameter<edm::InputTag>("Track")),
mvaV2IsoValuesMapInputTag_(iConfig.getUntrackedParameter<edm::InputTag>("mvaV2IsoValuesMap")),
mvaV2NoIsoValuesMapInputTag_(iConfig.getUntrackedParameter<edm::InputTag>("mvaV2NoIsoValuesMap")),
GsfTrackInputTag_(iConfig.getUntrackedParameter<edm::InputTag>("GsfTrack")),
SecondaryVertexInputTag_(iConfig.getUntrackedParameter<edm::InputTag>("secondaryVertices")),
PrimaryVertexInputTag_(iConfig.getUntrackedParameter<edm::InputTag>("primaryVertices"))
{
  electronsToken_    = consumes<edm::View<reco::GsfElectron> >(Electron_);
  TrackToken = consumes<reco::TrackCollection>(TrackInputTag_);
  mvaV2IsoValuesMapToken = consumes<edm::ValueMap<float>>(mvaV2IsoValuesMapInputTag_);
  mvaV2NoIsoValuesMapToken = consumes<edm::ValueMap<float>>(mvaV2NoIsoValuesMapInputTag_);
  GsfTrackToken = consumes<reco::GsfTrackCollection>(GsfTrackInputTag_);
  SecondaryVertexToken = consumes<reco::VertexCollection>(SecondaryVertexInputTag_);
  PrimaryVertexToken = consumes<reco::VertexCollection>(PrimaryVertexInputTag_);

  edm::Service<TFileService> fs;
  reco_tree = fs->make<TTree>("Events", "Events");

  reco_tree->Branch("numele",&numele);
  reco_tree->Branch("PFnumele",&PFnumele);
  reco_tree->Branch("ele_eta",&ele_eta);
  reco_tree->Branch("ele_phi",&ele_phi);
  reco_tree->Branch("trkIsoEle",&trkIsoEle);
  reco_tree->Branch("ele_pt",&ele_pt);
  reco_tree->Branch("ele_charge",&ele_charge);
  reco_tree->Branch("EleClosestTrackPt",&EleClosestTrackPt);
  reco_tree->Branch("EleTrkIsNoNull",&EleTrkIsNoNull);
  reco_tree->Branch("ele_x",&ele_x);
  reco_tree->Branch("ele_y",&ele_y);
  reco_tree->Branch("ele_z",&ele_z);
  reco_tree->Branch("ele_gsf_charge",&ele_gsf_charge);

  //V2 MVA scores
  reco_tree->Branch("ElectronMVAEstimatorRun2Fall17IsoV2Values",&ElectronMVAEstimatorRun2Fall17IsoV2Values);
  reco_tree->Branch("ElectronMVAEstimatorRun2Fall17NoIsoV2Values",&ElectronMVAEstimatorRun2Fall17NoIsoV2Values);

  //Track Info
  reco_tree->Branch("tr_pt",&tr_pt);
  reco_tree->Branch("tr_phi",&tr_phi);
  reco_tree->Branch("tr_eta",&tr_eta);
  reco_tree->Branch("numtr",&numtr);
  reco_tree->Branch("EleTRKref",&EleTRKref);
  reco_tree->Branch("tr_x",&tr_x);
  reco_tree->Branch("tr_y",&tr_y);
  reco_tree->Branch("tr_z",&tr_z);
  reco_tree->Branch("tr_charge",&tr_charge);

  //Gsf Info
  reco_tree->Branch("numgsf",&numgsf);
  reco_tree->Branch("gsf_x",&gsf_x);
  reco_tree->Branch("gsf_y",&gsf_y);
  reco_tree->Branch("gsf_z",&gsf_z);
  reco_tree->Branch("gsf_charge",&gsf_charge);
  reco_tree->Branch("gsf_pt",&gsf_pt);
  reco_tree->Branch("gsf_eta",&gsf_eta);
  reco_tree->Branch("gsf_phi",&gsf_phi);

  //Secondary Vertex
  reco_tree->Branch("sv_x",&sv_x);
  reco_tree->Branch("sv_y",&sv_y);
  reco_tree->Branch("sv_z",&sv_z);
  reco_tree->Branch("sv_chi2",&sv_chi2);
  reco_tree->Branch("sv_ndof",&sv_ndof);
  reco_tree->Branch("sv_tr_size",&sv_tr_size);
  reco_tree->Branch("sv_tr_pt",&sv_tr_pt);
  reco_tree->Branch("sv_tr_eta",&sv_tr_eta);
  reco_tree->Branch("sv_tr_phi",&sv_tr_phi);
  reco_tree->Branch("numsv",&numsv);
  reco_tree->Branch("SV3pt",&SV3pt);
  reco_tree->Branch("SV3phi",&SV3phi);
  reco_tree->Branch("SV3eta",&SV3eta);
  reco_tree->Branch("EleSV3trkref",&EleSV3trkref);
  reco_tree->Branch("SV3TrEleRef",&SV3TrEleRef);
  reco_tree->Branch("SV3charge",&SV3charge);
  reco_tree->Branch("SV3mass",&SV3mass);
  reco_tree->Branch("numsv3",&numsv3);
  reco_tree->Branch("SVdistanceToMaxPtPV",&SVdistanceToMaxPtPV);
  reco_tree->Branch("SV3x",&SV3x);
  reco_tree->Branch("SV3y",&SV3y);
  reco_tree->Branch("SV3z",&SV3z);

  //PrimaryVertex
  reco_tree->Branch("PVmaxPtX",&PVmaxPtX);
  reco_tree->Branch("PVmaxPtY",&PVmaxPtY);
  reco_tree->Branch("PVmaxPtZ",&PVmaxPtZ);

}

RecoAnalyzerV2::~RecoAnalyzerV2() {

}

/*float InvariantMass(float pt1,float eta1,float phi1,float pt2,float eta2,float phi2)
{
   TLorentzVector P,P1,P2;
   P1.SetPtEtaPhiM(pt1,eta1,phi1,0);
   P2.SetPtEtaPhiM(pt2,eta2,phi2,0);
   P=P1+P2;
   return P.M();
}*/

void RecoAnalyzerV2::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  using namespace edm;

  /*if (cachedCaloGeometryID_ != iSetup.get<CaloGeometryRecord>().cacheIdentifier()) {
        cachedCaloGeometryID_ = iSetup.get<CaloGeometryRecord>().cacheIdentifier();
        iSetup.get<CaloGeometryRecord>().get(caloGeometry_);
    }*/

  edm::Handle<edm::View<reco::GsfElectron> > electrons;
  iEvent.getByToken(electronsToken_, electrons);


  edm::Handle<reco::TrackCollection> TrackHandle;
  iEvent.getByToken(TrackToken, TrackHandle);

  edm::Handle<edm::ValueMap<float> > mvaV2NoIsoValues;
  iEvent.getByToken(mvaV2NoIsoValuesMapToken,mvaV2NoIsoValues);

  edm::Handle<edm::ValueMap<float> > mvaV2IsoValues;
  iEvent.getByToken(mvaV2IsoValuesMapToken,mvaV2IsoValues);

  edm::Handle<reco::GsfTrackCollection> GsfTrackHandle;
  iEvent.getByToken(GsfTrackToken,GsfTrackHandle);

  edm::Handle<reco::VertexCollection> secondaryVertexHandle;
  iEvent.getByToken(SecondaryVertexToken,secondaryVertexHandle);

  edm::Handle<reco::VertexCollection> primaryVertexHandle;
  iEvent.getByToken(PrimaryVertexToken,primaryVertexHandle);

  //edm::Handle<edm::View<reco::GsfElectron> > electrons;
  //iEvent.getByToken(electronsToken_, electrons);

  /*for (size_t l = 0; l < electrons->size(); ++l){
    const auto el = electrons->ptrAt(l);
    std::cout << (*mvaV2NoIsoValues)[el] << std::endl;}*/

  numele=0;
  PFnumele=0;

  ele_eta.clear();
  ele_phi.clear();
  ele_pt.clear();
  ele_charge.clear();
  trkIsoEle.clear();
  ele_x.clear();
  ele_y.clear();
  ele_z.clear();
  ele_gsf_charge.clear();


  ElectronMVAEstimatorRun2Fall17IsoV2Values.clear();
  ElectronMVAEstimatorRun2Fall17NoIsoV2Values.clear();
  EleClosestTrackPt.clear();
  EleTrkIsNoNull.clear();

  tr_pt.clear();
  tr_eta.clear();
  tr_phi.clear();
  EleTRKref.clear();
  tr_x.clear();
  tr_y.clear();
  tr_z.clear();
  tr_charge.clear();
  numtr=0;

  gsf_x.clear();
  gsf_y.clear();
  gsf_z.clear();
  gsf_charge.clear();
  gsf_pt.clear();
  gsf_eta.clear();
  gsf_phi.clear();
  numgsf=0;

  sv_x.clear();
  sv_y.clear();
  sv_z.clear();
  sv_chi2.clear();
  sv_ndof.clear();
  sv_tr_size.clear();
  sv_tr_pt.clear();
  sv_tr_eta.clear();
  sv_tr_phi.clear();
  SV3pt.clear();
  SV3phi.clear();
  SV3eta.clear();
  EleSV3trkref.clear();
  SV3TrEleRef.clear();
  SV3charge.clear();
  SV3mass.clear();
  numsv=0;
  numsv3=0;
  SVdistanceToMaxPtPV.clear();
  SV3x.clear();
  SV3y.clear();
  SV3z.clear();

  PVmaxPtX=0;
  PVmaxPtY=0;
  PVmaxPtZ=0;

  //helper variables
  float SumPt,TrackPtHelper=0.;
  int index,num5,TrackIndexHelper;
  TLorentzVector P,P1,P2,P3;
  std::vector<float> PVptSum;

  const reco::TrackCollection* tkColl = TrackHandle.product();
  math::XYZPoint pv(primaryVertexHandle->begin()->position());

  PVptSum.clear();
  for(auto pv = primaryVertexHandle->begin(); pv!= primaryVertexHandle->end(); ++pv)
  {
    float PtSum=0;
    for(long unsigned int p=0;p<pv->tracksSize();p++)
    {
      PtSum+=pv->trackRefAt(p)->pt();
    }
    PVptSum.push_back(PtSum);
  }

  auto MaxPtIndex = primaryVertexHandle->begin() + std::distance(PVptSum.begin(),std::max_element(PVptSum.begin(), PVptSum.end()));
  PVmaxPtX = MaxPtIndex->x();
  PVmaxPtY = MaxPtIndex->y();
  PVmaxPtZ = MaxPtIndex->z();


  for(auto sv = secondaryVertexHandle->begin(); sv!= secondaryVertexHandle->end(); ++sv)
  {
    if(sv->tracksSize()==3 && std::fabs(sv->trackRefAt(0)->charge()+sv->trackRefAt(1)->charge()+sv->trackRefAt(2)->charge())==1)
    {
      for(long unsigned int r=0;r<sv->tracksSize();r++)
      {
        SV3pt.push_back(sv->trackRefAt(r)->pt());
        SV3eta.push_back(sv->trackRefAt(r)->eta());
        SV3phi.push_back(sv->trackRefAt(r)->phi());
        SV3charge.push_back(sv->trackRefAt(r)->charge());
      }
      SV3x.push_back(sv->x());
      SV3y.push_back(sv->y());
      SV3z.push_back(sv->z());
      P1.SetPtEtaPhiM(SV3pt.at(0),SV3eta.at(0),SV3phi.at(0),0.);
      P2.SetPtEtaPhiM(SV3pt.at(1),SV3eta.at(1),SV3phi.at(1),0.);
      P3.SetPtEtaPhiM(SV3pt.at(2),SV3eta.at(2),SV3phi.at(2),0.);
      P=P1+P2+P3;
      SV3mass.push_back(P.M());
      SVdistanceToMaxPtPV.push_back(sqrt((sv->x()-PVmaxPtX)*(sv->x()-PVmaxPtX)+(sv->y()-PVmaxPtY)*(sv->y()-PVmaxPtY)+
                                    (sv->z()-PVmaxPtZ)*(sv->z()-PVmaxPtZ)));
      numsv3++;
    }
  }


  for(auto el = electrons->begin(); el != electrons->end(); ++el)
  {

    if(el->closestCtfTrackRef().isNonnull() && el->pt()>5.)  //what if the reference to track does not exist ?
    {
      auto ClosestTrackRef = tkColl->begin() + el->closestCtfTrackRef().index();
      int FoundSVmatch=0;
      for(long unsigned int b=0; b<SV3pt.size(); b++)
      {
        if((float)ClosestTrackRef->pt()==SV3pt.at(b) && (float)ClosestTrackRef->eta()==SV3eta.at(b) && (float)ClosestTrackRef->phi()==SV3phi.at(b))
        {
          EleSV3trkref.push_back(b);
          FoundSVmatch=1;
        }
      }
      if(FoundSVmatch==0)
      EleSV3trkref.push_back(99);
    }
    else if(el->pt()>5.)
    EleSV3trkref.push_back(99);

    if(el->pt()>5)
    {
      const auto iter = electrons->ptrAt(numele);
      ele_pt.push_back(el->pt());
      ele_eta.push_back(el->eta());
      ele_phi.push_back(el->phi());
      ele_charge.push_back(el->charge());
      const auto ele = electrons->ptrAt(numele);
      ElectronMVAEstimatorRun2Fall17NoIsoV2Values.push_back((*mvaV2NoIsoValues)[iter]);
      ElectronMVAEstimatorRun2Fall17IsoV2Values.push_back((*mvaV2IsoValues)[iter]);
      EleTrkIsNoNull.push_back(el->closestCtfTrackRef().isNonnull());
      ele_x.push_back(el->vx());
      ele_y.push_back(el->vy());
      ele_z.push_back(el->vz());
      ele_gsf_charge.push_back(el->gsfTrack()->charge());
      SumPt=0.;
      num5=0;
      TrackIndexHelper=0;
      if(el->closestCtfTrackRef().isNonnull())
      index = el->closestCtfTrackRef().index();  //adding this if(tr_pt(index)==sv_tr_pt->at(k)) tr_pt je ovdi bez cuta
      else
      index=5000; //put random value so code doesn't break
      for(auto tr = tkColl->begin(); tr != tkColl->end(); ++tr)
      {
        //std::cout<<el->pt()<<" "<<el->closestCtfTrackRef().isNonnull()<<" "<<el->closestCtfTrackRef().index()<<" "<<tr->pt()<<std::endl;
        double dR2 = reco::deltaR2(el->eta(),el->phi(),tr->eta(),tr->phi());
        double dz = fabs(tr->vz()-el->vz());
        double dEta = fabs(tr->eta()-el->eta());
        double TrackPt = tr->pt();
        if(dR2>0 && dR2<0.4 && dz<0.1 && TrackPt>2 && dEta>0.005)
        SumPt+=TrackPt;
        if(tr->pt()<5. && TrackIndexHelper<index)
        num5++; //count the number of tracks < 5GeV
        if(TrackIndexHelper==index)
        TrackPtHelper = tr->pt();
        TrackIndexHelper++;
      }
      //std::cout<<el->pt()<<" "<<TrackPtHelper<<" "<<index-num5<<std::endl;
      trkIsoEle.push_back(SumPt);
      if(el->closestCtfTrackRef().isNonnull() && TrackPtHelper>5.)
      EleTRKref.push_back(index-num5);
      else
      EleTRKref.push_back(99);
      EleClosestTrackPt.push_back(TrackPtHelper);
      numele++;
    }
  }

  for(long unsigned int b=0; b<SV3pt.size();b++)
  {
    int RefSV3toEl=99;
    for(long unsigned c=0; c<EleSV3trkref.size();c++)
    {
      if(EleSV3trkref.at(c)==(int)b)
      RefSV3toEl=c;
    }
    SV3TrEleRef.push_back(RefSV3toEl);
  }


  for(auto sv = secondaryVertexHandle->begin(); sv!= secondaryVertexHandle->end(); ++sv)
  {
    sv_x.push_back(sv->x());
    sv_y.push_back(sv->y());
    sv_z.push_back(sv->z());
    sv_chi2.push_back(sv->chi2());
    sv_ndof.push_back(sv->ndof());
    sv_tr_size.push_back(sv->tracksSize());
    for(long unsigned int r=0;r<sv->tracksSize();r++)
    {
      sv_tr_pt.push_back(sv->trackRefAt(r)->pt());
      sv_tr_eta.push_back(sv->trackRefAt(r)->eta());
      sv_tr_phi.push_back(sv->trackRefAt(r)->phi());
    }
    if(sv->tracksSize()==3)
    {
      for(long unsigned int s=0;s<sv->tracksSize();s++)
      {
        SV3pt.push_back(sv->trackRefAt(s)->pt());
        SV3eta.push_back(sv->trackRefAt(s)->eta());
        SV3phi.push_back(sv->trackRefAt(s)->phi());
      }
    }
    numsv++;
    //const auto &ref = sv->trackRefAt(0);
    //std::cout<<sv->tracksSize()<<" "<<ref->pt()<<std::endl;
  }

  for(auto gsf = GsfTrackHandle->begin(); gsf!= GsfTrackHandle->end(); ++gsf)
  {
    gsf_x.push_back(gsf->vx());
    gsf_y.push_back(gsf->vy());
    gsf_z.push_back(gsf->vz());
    gsf_pt.push_back(gsf->pt());
    gsf_eta.push_back(gsf->eta());
    gsf_phi.push_back(gsf->phi());
    gsf_charge.push_back(gsf->charge());
    numgsf++;
  }

  for(auto tr = tkColl->begin(); tr != tkColl->end(); ++tr)
  {
    if(tr->pt()>5)
    {
    tr_pt.push_back(tr->pt());
    tr_eta.push_back(tr->eta());
    tr_phi.push_back(tr->phi());
    tr_x.push_back(tr->vx());
    tr_y.push_back(tr->vy());
    tr_z.push_back(tr->vz());
    tr_charge.push_back(tr->charge());
    //std::cout<<tr->normalizedChi2()<<std::endl;
    numtr++;
    //std::cout<<tr->charge()<<std::endl;
    }
  }


  /*for (auto el = electrons->begin(); el != electrons->end(); ++el)
  {
    if(el->pt()>5)
    {
      const auto iter = electrons->ptrAt(numele);
      ele_pt.push_back(el->pt());
      ele_eta.push_back(el->eta());
      ele_phi.push_back(el->phi());
      ele_charge.push_back(el->charge());
      const auto ele = electrons->ptrAt(numele);
      ElectronMVAEstimatorRun2Fall17NoIsoV2Values.push_back((*mvaV2NoIsoValues)[iter]);
      ElectronMVAEstimatorRun2Fall17IsoV2Values.push_back((*mvaV2IsoValues)[iter]);
      EleTrkIsNoNull.push_back(el->closestCtfTrackRef().isNonnull());
      ele_x.push_back(el->vx());
      ele_y.push_back(el->vy());
      ele_z.push_back(el->vz());
      ele_gsf_charge.push_back(el->gsfTrack()->charge());
      SumPt=0.;
      num5=0;
      TrackIndexHelper=0;
      if(el->closestCtfTrackRef().isNonnull())
      index = el->closestCtfTrackRef().index();  //adding this if(tr_pt(index)==sv_tr_pt->at(k)) tr_pt je ovdi bez cuta
      else
      index=5000; //put random value so code doesn't break
      for(auto tr = tkColl->begin(); tr != tkColl->end(); ++tr)
      {
        //std::cout<<el->pt()<<" "<<el->closestCtfTrackRef().isNonnull()<<" "<<el->closestCtfTrackRef().index()<<" "<<tr->pt()<<std::endl;
        double dR2 = reco::deltaR2(el->eta(),el->phi(),tr->eta(),tr->phi());
        double dz = fabs(tr->vz()-el->vz());
        double dEta = fabs(tr->eta()-el->eta());
        double TrackPt = tr->pt();
        if(dR2>0 && dR2<0.4 && dz<0.1 && TrackPt>2 && dEta>0.005)
        SumPt+=TrackPt;
        if(tr->pt()<5. && TrackIndexHelper<index)
        num5++; //count the number of tracks < 5GeV
        if(TrackIndexHelper==index)
        TrackPtHelper = tr->pt();
        TrackIndexHelper++;
      }
      //std::cout<<el->pt()<<" "<<TrackPtHelper<<" "<<index-num5<<std::endl;
      trkIsoEle.push_back(SumPt);
      if(el->closestCtfTrackRef().isNonnull() && TrackPtHelper>5.)
      EleTRKref.push_back(index-num5);
      else
      EleTRKref.push_back(99);
      EleClosestTrackPt.push_back(TrackPtHelper);
      numele++;
    }

  }*/


  if(numele==1 || numele==2) //add charge condition ?
  {
    reco_tree->Fill();
  }

}


void RecoAnalyzerV2::beginJob() {

}

void RecoAnalyzerV2::endJob() {

}

void RecoAnalyzerV2::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);

}

DEFINE_FWK_MODULE(RecoAnalyzerV2);
