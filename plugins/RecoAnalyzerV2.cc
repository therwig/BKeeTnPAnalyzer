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
  edm::EDGetTokenT<reco::SuperClusterCollection> SuperClusterEBToken;
  edm::EDGetTokenT<reco::SuperClusterCollection> SuperClusterEEToken;
  edm::EDGetTokenT<reco::ConversionCollection> ConversionToken;
  edm::EDGetTokenT<reco::BeamSpot> BeamSpotToken;
  edm::EDGetTokenT<double> RhoToken;
  edm::EDGetTokenT<HBHERecHitCollection> hbheRHcToken;
  edm::EDGetTokenT<reco::TrackCollection> TrackToken;
  edm::EDGetTokenT<reco::VertexCollection> VertexToken;
  edm::EDGetTokenT<edm::ValueMap<float> > mvaV2IsoValuesMapToken;
  edm::EDGetTokenT<edm::ValueMap<float> > mvaV2NoIsoValuesMapToken;
  edm::EDGetTokenT<reco::GsfTrackCollection> GsfTrackToken;
  edm::EDGetTokenT<reco::VertexCollection> SecondaryVertexToken;

  edm::InputTag Electron_;
  edm::InputTag SuperClusterEB_;
  edm::InputTag SuperClusterEE_;
  edm::InputTag Conversion_;
  edm::InputTag BeamSpot_;
  edm::InputTag Rho_;
  edm::InputTag hcalRecHitsInputHBHE_;
  edm::InputTag TrackInputTag_;
  edm::InputTag VertexInputTag_;
  edm::InputTag mvaV2IsoValuesMapInputTag_;
  edm::InputTag mvaV2NoIsoValuesMapInputTag_;
  edm::InputTag GsfTrackInputTag_;
  edm::InputTag SecondaryVertexInputTag_;

  TTree *reco_tree;
  std::vector<float> ele_pt,ele_eta,ele_phi,full5x5_sigmaIetaIeta,dEtaSeed,dPhiIn,HoverE,relIso,Ep,
  el_sc_eta, el_sc_E,el_sc_phi,sc_eta, sc_pt, sc_phi, sc_E,trkIsoSC,trkIsoEle,EleClosestTrackPt,
  tr_pt,tr_eta,tr_phi,ElectronMVAEstimatorRun2Fall17IsoV2Values,ElectronMVAEstimatorRun2Fall17NoIsoV2Values,
  ele_z,ele_x,ele_y,tr_z,tr_x,tr_y,tr_charge,gsf_z,gsf_x,gsf_y, sv_x, sv_y, sv_z,sv_chi2,sv_ndof;
  std::vector<int> ele_charge,ExpMissInnerHits,EleTRKref,EleTrkIsNoNull,ele_gsf_charge,gsf_charge;
  std::vector<bool> PassConversionVeto,CutBasedLoose,CutBasedMedium,CutBasedTight;
  float rho;
  int numele, PFnumele,numSC,EleCounter,match,numtr,numgsf;
  TLorentzVector P,P0,P1,p,p0,p1;

  //unsigned long long cachedCaloGeometryID_;
  //edm::ESHandle<CaloGeometry> caloGeometry_;


};

RecoAnalyzerV2::RecoAnalyzerV2(const edm::ParameterSet& iConfig):
Electron_(iConfig.getUntrackedParameter<edm::InputTag>("Electron")),
SuperClusterEB_(iConfig.getUntrackedParameter<edm::InputTag>("SuperClusterEB")),
SuperClusterEE_(iConfig.getUntrackedParameter<edm::InputTag>("SuperClusterEE")),
Conversion_(iConfig.getUntrackedParameter<edm::InputTag>("Conversions")),
BeamSpot_(iConfig.getUntrackedParameter<edm::InputTag>("BeamSpot")),
Rho_(iConfig.getUntrackedParameter<edm::InputTag>("Rho")),
hcalRecHitsInputHBHE_(iConfig.getUntrackedParameter<edm::InputTag>("HBHERecHit")),
TrackInputTag_(iConfig.getUntrackedParameter<edm::InputTag>("Track")),
VertexInputTag_(iConfig.getUntrackedParameter<edm::InputTag>("Vertex")),
mvaV2IsoValuesMapInputTag_(iConfig.getUntrackedParameter<edm::InputTag>("mvaV2IsoValuesMap")),
mvaV2NoIsoValuesMapInputTag_(iConfig.getUntrackedParameter<edm::InputTag>("mvaV2NoIsoValuesMap")),
GsfTrackInputTag_(iConfig.getUntrackedParameter<edm::InputTag>("GsfTrack")),
SecondaryVertexInputTag_(iConfig.getUntrackedParameter<edm::InputTag>("secondaryVertices"))
{
  electronsToken_    = consumes<edm::View<reco::GsfElectron> >(Electron_);
  SuperClusterEBToken = consumes<reco::SuperClusterCollection>(SuperClusterEB_);
  SuperClusterEEToken = consumes<reco::SuperClusterCollection>(SuperClusterEE_);
  ConversionToken = consumes<reco::ConversionCollection>(Conversion_);
  BeamSpotToken = consumes<reco::BeamSpot>(BeamSpot_);
  RhoToken = consumes<double>(Rho_);
  hbheRHcToken= consumes<HBHERecHitCollection>(hcalRecHitsInputHBHE_);
  TrackToken = consumes<reco::TrackCollection>(TrackInputTag_);
  VertexToken = consumes<reco::VertexCollection>(VertexInputTag_);
  mvaV2IsoValuesMapToken = consumes<edm::ValueMap<float>>(mvaV2IsoValuesMapInputTag_);
  mvaV2NoIsoValuesMapToken = consumes<edm::ValueMap<float>>(mvaV2NoIsoValuesMapInputTag_);
  GsfTrackToken = consumes<reco::GsfTrackCollection>(GsfTrackInputTag_);
  SecondaryVertexToken = consumes<reco::VertexCollection>(SecondaryVertexInputTag_);

  edm::Service<TFileService> fs;
  reco_tree = fs->make<TTree>("Events", "Events");

  reco_tree->Branch("numele",&numele);
  reco_tree->Branch("PFnumele",&PFnumele);
  reco_tree->Branch("ele_eta",&ele_eta);
  reco_tree->Branch("ele_phi",&ele_phi);
  reco_tree->Branch("trkIsoEle",&trkIsoEle);
  reco_tree->Branch("ele_pt",&ele_pt);
  reco_tree->Branch("ele_charge",&ele_charge);
  reco_tree->Branch("full5x5_sigmaIetaIeta",&full5x5_sigmaIetaIeta);
  reco_tree->Branch("dEtaSeed",&dEtaSeed);
  reco_tree->Branch("dPhiIn",&dPhiIn);
  reco_tree->Branch("HoverE",&HoverE);
  reco_tree->Branch("relIso",&relIso);
  reco_tree->Branch("Ep",&Ep);
  reco_tree->Branch("ExpMissInnerHits",&ExpMissInnerHits);
  reco_tree->Branch("PassConversionVeto",&PassConversionVeto);
  reco_tree->Branch("CutBasedLoose",&CutBasedLoose);
  reco_tree->Branch("CutBasedMedium",&CutBasedMedium);
  reco_tree->Branch("CutBasedTight",&CutBasedTight);
  reco_tree->Branch("rho",&rho);
  reco_tree->Branch("EleClosestTrackPt",&EleClosestTrackPt);
  reco_tree->Branch("EleTrkIsNoNull",&EleTrkIsNoNull);
  reco_tree->Branch("ele_x",&ele_x);
  reco_tree->Branch("ele_y",&ele_y);
  reco_tree->Branch("ele_z",&ele_z);
  reco_tree->Branch("ele_gsf_charge",&ele_gsf_charge);

  //V2 MVA scores
  reco_tree->Branch("ElectronMVAEstimatorRun2Fall17IsoV2Values",&ElectronMVAEstimatorRun2Fall17IsoV2Values);
  reco_tree->Branch("ElectronMVAEstimatorRun2Fall17NoIsoV2Values",&ElectronMVAEstimatorRun2Fall17NoIsoV2Values);

  //electron supercluster variables

  reco_tree->Branch("el_sc_eta",&el_sc_eta);
  reco_tree->Branch("el_sc_E",&el_sc_E);
  reco_tree->Branch("el_sc_phi",&el_sc_phi);

  //SC variables
  reco_tree->Branch("numSC",&numSC);
  reco_tree->Branch("sc_eta",&sc_eta);
  reco_tree->Branch("sc_E",&sc_E);
  reco_tree->Branch("sc_pt",&sc_pt);
  reco_tree->Branch("sc_phi",&sc_phi);

  //TrackIsolation
  reco_tree->Branch("trkIsoSC",&trkIsoSC);

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

  //Secondary Vertex
  reco_tree->Branch("sv_x",&sv_x);
  reco_tree->Branch("sv_y",&sv_y);
  reco_tree->Branch("sv_z",&sv_z);
  reco_tree->Branch("sv_chi2",&sv_chi2);
  reco_tree->Branch("sv_ndof",&sv_ndof);

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

  edm::Handle< std::vector<reco::SuperCluster>> superclustersEB;
  iEvent.getByToken(SuperClusterEBToken, superclustersEB);

  edm::Handle< std::vector<reco::SuperCluster>> superclustersEE;
  iEvent.getByToken(SuperClusterEEToken, superclustersEE);

  edm::Handle<std::vector<reco::Conversion>> hConversions;
  iEvent.getByToken(ConversionToken, hConversions);

  edm::Handle<reco::BeamSpot> BeamSpot;
  iEvent.getByToken(BeamSpotToken, BeamSpot);
  const reco::BeamSpot &beamspot = *BeamSpot.product();

  edm::Handle<double> rhoHandle;
  iEvent.getByToken(RhoToken, rhoHandle);

  edm::Handle<HBHERecHitCollection> hcalRecHitsHandle;
  iEvent.getByToken(hbheRHcToken, hcalRecHitsHandle);

  edm::Handle<reco::TrackCollection> TrackHandle;
  iEvent.getByToken(TrackToken, TrackHandle);

  edm::Handle<reco::VertexCollection> VertexHandle;
  iEvent.getByToken(VertexToken, VertexHandle);

  edm::Handle<edm::ValueMap<float> > mvaV2NoIsoValues;
  iEvent.getByToken(mvaV2NoIsoValuesMapToken,mvaV2NoIsoValues);

  edm::Handle<edm::ValueMap<float> > mvaV2IsoValues;
  iEvent.getByToken(mvaV2IsoValuesMapToken,mvaV2IsoValues);

  edm::Handle<reco::GsfTrackCollection> GsfTrackHandle;
  iEvent.getByToken(GsfTrackToken,GsfTrackHandle);

  edm::Handle<reco::VertexCollection> secondaryVertexHandle;
  iEvent.getByToken(SecondaryVertexToken,secondaryVertexHandle);

  //edm::Handle<edm::View<reco::GsfElectron> > electrons;
  //iEvent.getByToken(electronsToken_, electrons);

  /*for (size_t l = 0; l < electrons->size(); ++l){
    const auto el = electrons->ptrAt(l);
    std::cout << (*mvaV2NoIsoValues)[el] << std::endl;}*/

  numSC=0;
  numele=0;
  PFnumele=0;
  rho=0;

  ele_eta.clear();
  ele_phi.clear();
  ele_pt.clear();
  ele_charge.clear();
  trkIsoEle.clear();
  el_sc_eta.clear();
  el_sc_E.clear();
  el_sc_phi.clear();
  ele_x.clear();
  ele_y.clear();
  ele_z.clear();
  ele_gsf_charge.clear();

  sc_eta.clear();
  sc_E.clear();
  sc_pt.clear();
  sc_phi.clear();

  full5x5_sigmaIetaIeta.clear();
  dEtaSeed.clear();
  dPhiIn.clear();
  HoverE.clear();
  relIso.clear();
  Ep.clear();
  ExpMissInnerHits.clear();
  PassConversionVeto.clear();
  CutBasedLoose.clear();
  CutBasedMedium.clear();
  CutBasedTight.clear();
  trkIsoSC.clear();
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
  numgsf=0;

  sv_x.clear();
  sv_y.clear();
  sv_z.clear();
  sv_chi2.clear();
  sv_ndof.clear();

  //helper variables
  float SumPt,TrackPtHelper=0.;
  int t,index,num5,TrackIndexHelper;

  const reco::TrackCollection* tkColl = TrackHandle.product();
  math::XYZPoint pv(VertexHandle->begin()->position());
  rho = *(rhoHandle.product());

  //Test closest track
  /*for (auto el = electrons->begin(); el != electrons->end(); ++el)
  {
    unsigned int num=0;
    if(el->closestCtfTrackRef().isNonnull())
    for(auto tr = tkColl->begin(); tr != tkColl->end(); ++tr)
    {
      if(el->closestCtfTrackRef().index()==num)
      std::cout<<el->pt()<<" "<<el->eta()<<" "<<el->phi()<<" "<<tr->pt()<<" "<<tr->eta()<<" "<<tr->phi()<<std::endl;
      num++;
    }
  }*/

  for(auto sv = secondaryVertexHandle->begin(); sv!= secondaryVertexHandle->end(); ++sv)
  {
    sv_x.push_back(sv->x());
    sv_y.push_back(sv->y());
    sv_z.push_back(sv->z());
    sv_chi2.push_back(sv->chi2());
    sv_ndof.push_back(sv->ndof());
  }

  for(auto gsf = GsfTrackHandle->begin(); gsf!= GsfTrackHandle->end(); ++gsf)
  {
    gsf_x.push_back(gsf->vx());
    gsf_y.push_back(gsf->vy());
    gsf_z.push_back(gsf->vz());
    gsf_charge.push_back(gsf->charge());
    numgsf++;
  }

  for(auto tr = tkColl->begin(); tr != tkColl->end(); ++tr)
  {
    //track electron matching
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

  for (auto sc = superclustersEB->cbegin(); sc != superclustersEB->cend(); ++sc)
  {
    if(sc->energy()*sin(2*atan(exp(-1.0*sc->eta())))>5.)
    {
      SumPt=0;
      for(t = 0; t < numtr; t++)
      {
        double dR2 = reco::deltaR2(sc->eta(),sc->phi(),tr_eta.at(t),tr_phi.at(t));
        double TrackPt = tr_pt.at(t);
        if(dR2>0 && dR2<0.4 && TrackPt>2)
        SumPt+=TrackPt;
      }

      trkIsoSC.push_back(SumPt);
      sc_E.push_back(sc->energy());
      sc_eta.push_back(sc->eta());
      sc_pt.push_back(sc->energy()*sin(2*atan(exp(-1.0*sc->eta())))); //Et=pT, transverse energy is practically equal to transverse momentum
      sc_phi.push_back(sc->phi());
      numSC++;
    }
  }

  for (auto sc = superclustersEE->cbegin(); sc != superclustersEE->cend(); ++sc)
  {
    if(sc->energy()*sin(2*atan(exp(-1.0*sc->eta())))>5.)
    {
      SumPt=0;
      for(t = 0; t < numtr; t++)
      {
        double dR2 = reco::deltaR2(sc->eta(),sc->phi(),tr_eta.at(t),tr_phi.at(t));
        double TrackPt = tr_pt.at(t);
        if(dR2>0 && dR2<0.4 && TrackPt>2)
        SumPt+=TrackPt;
      }

      trkIsoSC.push_back(SumPt);
      sc_E.push_back(sc->energy());
      sc_eta.push_back(sc->eta());
      sc_pt.push_back(sc->energy()*sin(2*atan(exp(-1.0*sc->eta())))); //Et=pT, transverse energy is practically equal to transverse momentum
      sc_phi.push_back(sc->phi());
      numSC++;
    }
  }


  for (auto el = electrons->begin(); el != electrons->end(); ++el)
  {
    if(el->pt()>5)
    {
      const auto iter = electrons->ptrAt(numele);
      ele_pt.push_back(el->pt());
      ele_eta.push_back(el->eta());
      ele_phi.push_back(el->phi());
      el_sc_eta.push_back(el->superCluster()->eta());
      el_sc_phi.push_back(el->superCluster()->phi());
      el_sc_E.push_back(el->superCluster()->energy());
      ele_charge.push_back(el->charge());
      full5x5_sigmaIetaIeta.push_back(el->full5x5_sigmaIetaIeta());
      dEtaSeed.push_back(el->superCluster().isNonnull() && el->superCluster()->seed().isNonnull() ?
      el->deltaEtaSuperClusterTrackAtVtx() - el->superCluster()->eta() + el->superCluster()->seed()->eta() : std::numeric_limits<float>::max());
      dPhiIn.push_back(el->deltaPhiSuperClusterTrackAtVtx());
      HoverE.push_back(el->hadronicOverEm());
      const float ecal_energy_inverse = 1.0/el->ecalEnergy();
      const float eSCoverP = el->eSuperClusterOverP();
      Ep.push_back(std::abs(1.0 - eSCoverP)*ecal_energy_inverse);
      relIso.push_back((float)(el->pfIsolationVariables().sumChargedHadronPt+std::max(float(0.0),el->pfIsolationVariables().sumNeutralHadronEt+el->pfIsolationVariables().sumPhotonEt))/el->pt());
      constexpr auto missingHitType = reco::HitPattern::MISSING_INNER_HITS;
      ExpMissInnerHits.push_back(el->gsfTrack()->hitPattern().numberOfLostHits(missingHitType));
      PassConversionVeto.push_back(!ConversionTools::hasMatchedConversion(*el, *hConversions, beamspot.position()));
      const auto ele = electrons->ptrAt(numele);
      CutBasedLoose.push_back(CutBasedLooseID(full5x5_sigmaIetaIeta[numele],dEtaSeed[numele],dPhiIn[numele],HoverE[numele],Ep[numele],relIso[numele],ExpMissInnerHits[numele],PassConversionVeto[numele],el->superCluster()->energy(),*(rhoHandle.product()),el->pt(),el->superCluster()->eta()));
      CutBasedMedium.push_back(CutBasedMediumID(full5x5_sigmaIetaIeta[numele],dEtaSeed[numele],dPhiIn[numele],HoverE[numele],Ep[numele],relIso[numele],ExpMissInnerHits[numele],PassConversionVeto[numele],el->superCluster()->energy(),*(rhoHandle.product()),el->pt(),el->superCluster()->eta()));
      CutBasedTight.push_back(CutBasedTightID(full5x5_sigmaIetaIeta[numele],dEtaSeed[numele],dPhiIn[numele],HoverE[numele],Ep[numele],relIso[numele],ExpMissInnerHits[numele],PassConversionVeto[numele],el->superCluster()->energy(),*(rhoHandle.product()),el->pt(),el->superCluster()->eta()));
      ElectronMVAEstimatorRun2Fall17NoIsoV2Values.push_back((*mvaV2NoIsoValues)[iter]);
      ElectronMVAEstimatorRun2Fall17IsoV2Values.push_back((*mvaV2IsoValues)[iter]);
      EleTrkIsNoNull.push_back(el->closestCtfTrackRef().isNonnull());
      ele_x.push_back(el->vx());
      ele_y.push_back(el->vy());
      ele_z.push_back(el->vz());
      ele_charge.push_back(el->charge());
      ele_gsf_charge.push_back(el->gsfTrack()->charge());
      SumPt=0.;
      num5=0;
      TrackIndexHelper=0;
      if(el->closestCtfTrackRef().isNonnull())
      index = el->closestCtfTrackRef().index();
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


  /*for(i=0;i<numele;i++)
  {
    for(k=0;k<numtr;k++)
    {
      //if(EleTRKref[i]==k)
      std::cout<<"bla:"<<ele_pt[i]<<" "<<tr_pt[k]<<" "<<EleTRKref[i]<<std::endl;
    }
  }*/

  if(numele==1 || (numele==2 && ele_charge[0]==-1*ele_charge[1]))
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
