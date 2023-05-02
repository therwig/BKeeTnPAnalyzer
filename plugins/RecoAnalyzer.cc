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

class RecoAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources> {
public:
  explicit RecoAnalyzer(const edm::ParameterSet&);
  ~RecoAnalyzer();

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
  edm::EDGetTokenT<edm::ValueMap<float> > eleIdLooseToken;
  edm::EDGetTokenT<edm::ValueMap<float> > eleIdRobustLooseToken;
  edm::EDGetTokenT<edm::ValueMap<float> > eleIdRobustTightToken;
  edm::EDGetTokenT<edm::ValueMap<float> > eleIdTightToken;
  edm::EDGetTokenT<double> RhoToken;
  edm::EDGetTokenT<HBHERecHitCollection> hbheRHcToken;
  edm::EDGetTokenT<reco::TrackCollection> TrackToken;
  edm::EDGetTokenT<reco::VertexCollection> VertexToken;
  edm::EDGetTokenT<edm::ValueMap<float> > mvaV2IsoValuesMapToken;
  edm::EDGetTokenT<edm::ValueMap<float> > mvaV2NoIsoValuesMapToken;

  edm::InputTag Electron_;
  edm::InputTag SuperClusterEB_;
  edm::InputTag SuperClusterEE_;
  edm::InputTag Conversion_;
  edm::InputTag BeamSpot_;
  edm::InputTag eleIdLoose_;
  edm::InputTag eleIdRobustLoose_;
  edm::InputTag eleIdRobustTight_;
  edm::InputTag eleIdTight_;
  edm::InputTag Rho_;
  edm::InputTag hcalRecHitsInputHBHE_;
  edm::InputTag TrackInputTag_;
  edm::InputTag VertexInputTag_;
  edm::InputTag mvaV2IsoValuesMapInputTag_;
  edm::InputTag mvaV2NoIsoValuesMapInputTag_;

  TTree *reco_tree;
  std::vector<float> ele_pt,ele_eta,ele_phi,full5x5_sigmaIetaIeta,dEtaSeed,dPhiIn,HoverE,relIso,Ep,
  el_sc_eta, el_sc_E,el_sc_phi,sc_eta, sc_pt, sc_phi, sc_E, eleIdLoose,eleIdTight,eleIdRobustLoose,eleIdRobustTight,trkIsoSC,trkIsoEle,
  tr_pt,tr_eta,tr_phi, EleSC_mass,EleTrk_mass,ElectronMVAEstimatorRun2Fall17IsoV2Values,ElectronMVAEstimatorRun2Fall17NoIsoV2Values;
  std::vector<int> ele_charge,SCref,ExpMissInnerHits,EleTRKref;
  std::vector<bool> PassConversionVeto,CutBasedLoose,CutBasedMedium,CutBasedTight;
  float rho;
  int numele, PFnumele,numSC,EleCounter,match,numtr;
  TLorentzVector P,P0,P1,p,p0,p1;

  //unsigned long long cachedCaloGeometryID_;
  //edm::ESHandle<CaloGeometry> caloGeometry_;


};

RecoAnalyzer::RecoAnalyzer(const edm::ParameterSet& iConfig):
Electron_(iConfig.getUntrackedParameter<edm::InputTag>("Electron")),
SuperClusterEB_(iConfig.getUntrackedParameter<edm::InputTag>("SuperClusterEB")),
SuperClusterEE_(iConfig.getUntrackedParameter<edm::InputTag>("SuperClusterEE")),
Conversion_(iConfig.getUntrackedParameter<edm::InputTag>("Conversions")),
BeamSpot_(iConfig.getUntrackedParameter<edm::InputTag>("BeamSpot")),
eleIdLoose_(iConfig.getUntrackedParameter<edm::InputTag>("eleIdLoose")),
eleIdRobustLoose_(iConfig.getUntrackedParameter<edm::InputTag>("eleIdRobustLoose")),
eleIdRobustTight_(iConfig.getUntrackedParameter<edm::InputTag>("eleIdRobustTight")),
eleIdTight_(iConfig.getUntrackedParameter<edm::InputTag>("eleIdTight")),
Rho_(iConfig.getUntrackedParameter<edm::InputTag>("Rho")),
hcalRecHitsInputHBHE_(iConfig.getUntrackedParameter<edm::InputTag>("HBHERecHit")),
TrackInputTag_(iConfig.getUntrackedParameter<edm::InputTag>("Track")),
VertexInputTag_(iConfig.getUntrackedParameter<edm::InputTag>("Vertex")),
mvaV2IsoValuesMapInputTag_(iConfig.getUntrackedParameter<edm::InputTag>("mvaV2IsoValuesMap")),
mvaV2NoIsoValuesMapInputTag_(iConfig.getUntrackedParameter<edm::InputTag>("mvaV2NoIsoValuesMap"))
{
  electronsToken_    = consumes<edm::View<reco::GsfElectron> >(Electron_);
  SuperClusterEBToken = consumes<reco::SuperClusterCollection>(SuperClusterEB_);
  SuperClusterEEToken = consumes<reco::SuperClusterCollection>(SuperClusterEE_);
  ConversionToken = consumes<reco::ConversionCollection>(Conversion_);
  BeamSpotToken = consumes<reco::BeamSpot>(BeamSpot_);
  eleIdLooseToken = consumes<edm::ValueMap<float>>(eleIdLoose_);
  eleIdRobustLooseToken = consumes<edm::ValueMap<float>>(eleIdRobustLoose_);
  eleIdRobustTightToken = consumes<edm::ValueMap<float>>(eleIdRobustTight_);
  eleIdTightToken = consumes<edm::ValueMap<float>>(eleIdTight_);
  RhoToken = consumes<double>(Rho_);
  hbheRHcToken= consumes<HBHERecHitCollection>(hcalRecHitsInputHBHE_);
  TrackToken = consumes<reco::TrackCollection>(TrackInputTag_);
  VertexToken = consumes<reco::VertexCollection>(VertexInputTag_);
  mvaV2IsoValuesMapToken = consumes<edm::ValueMap<float>>(mvaV2IsoValuesMapInputTag_);
  mvaV2NoIsoValuesMapToken = consumes<edm::ValueMap<float>>(mvaV2NoIsoValuesMapInputTag_);

  edm::Service<TFileService> fs;
  reco_tree = fs->make<TTree>("Events", "Events");

  reco_tree->Branch("numele",&numele);
  reco_tree->Branch("PFnumele",&PFnumele);
  reco_tree->Branch("ele_eta",&ele_eta);
  reco_tree->Branch("ele_phi",&ele_phi);
  reco_tree->Branch("EleSC_mass",&EleSC_mass);

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
  reco_tree->Branch("eleIdLoose",&eleIdLoose);
  reco_tree->Branch("eleIdTight",&eleIdTight);
  reco_tree->Branch("eleIdRobustLoose",&eleIdRobustLoose);
  reco_tree->Branch("eleIdRobustTight",&eleIdRobustTight);
  reco_tree->Branch("CutBasedLoose",&CutBasedLoose);
  reco_tree->Branch("CutBasedMedium",&CutBasedMedium);
  reco_tree->Branch("CutBasedTight",&CutBasedTight);
  reco_tree->Branch("rho",&rho);

  //V2 MVA scores
  reco_tree->Branch("ElectronMVAEstimatorRun2Fall17IsoV2Values",&ElectronMVAEstimatorRun2Fall17IsoV2Values);
  reco_tree->Branch("ElectronMVAEstimatorRun2Fall17NoIsoV2Values",&ElectronMVAEstimatorRun2Fall17NoIsoV2Values);

  //electron supercluster variables

  reco_tree->Branch("el_sc_eta",&el_sc_eta);
  reco_tree->Branch("el_sc_E",&el_sc_E);
  reco_tree->Branch("el_sc_phi",&el_sc_phi);

  //SC variables
  reco_tree->Branch("SCref",&SCref);
  reco_tree->Branch("numSC",&numSC);
  reco_tree->Branch("sc_eta",&sc_eta);
  reco_tree->Branch("sc_E",&sc_E);
  reco_tree->Branch("sc_pt",&sc_pt);
  reco_tree->Branch("sc_phi",&sc_phi);

  //TrackIsolation
  reco_tree->Branch("trkIsoSC",&trkIsoSC);
  reco_tree->Branch("trkIsoEle",&trkIsoEle);

  //Track Info
  reco_tree->Branch("tr_pt",&tr_pt);
  reco_tree->Branch("tr_phi",&tr_phi);
  reco_tree->Branch("tr_eta",&tr_eta);
  reco_tree->Branch("numtr",&numtr);
  reco_tree->Branch("EleTRKref",&EleTRKref);
  reco_tree->Branch("EleTrk_mass",&EleTrk_mass);
  //reco_tree->Branch("HcalSum",&HcalSum);

  //cachedCaloGeometryID_ = 0;

}

RecoAnalyzer::~RecoAnalyzer() {

}

float InvariantMass(float pt1,float eta1,float phi1,float pt2,float eta2,float phi2)
{
   TLorentzVector P,P1,P2;
   P1.SetPtEtaPhiM(pt1,eta1,phi1,0);
   P2.SetPtEtaPhiM(pt2,eta2,phi2,0);
   P=P1+P2;
   return P.M();
}

void RecoAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
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

  edm::Handle<edm::ValueMap<float> > ele_id_decisions_loose;
  iEvent.getByToken(eleIdLooseToken ,ele_id_decisions_loose);

  edm::Handle<edm::ValueMap<float> > ele_id_decisions_robust_loose;
  iEvent.getByToken(eleIdRobustLooseToken ,ele_id_decisions_robust_loose);

  edm::Handle<edm::ValueMap<float> > ele_id_decisions_robust_tight;
  iEvent.getByToken(eleIdRobustTightToken ,ele_id_decisions_robust_tight);

  edm::Handle<edm::ValueMap<float> > ele_id_decisions_tight;
  iEvent.getByToken(eleIdTightToken ,ele_id_decisions_tight);

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

  //edm::Handle<edm::View<reco::GsfElectron> > electrons;
  //iEvent.getByToken(electronsToken_, electrons);

  /*for (size_t l = 0; l < electrons->size(); ++l){
    const auto el = electrons->ptrAt(l);
    std::cout << (*mvaV2NoIsoValues)[el] << std::endl;}*/

  numSC=0;
  numele=0;
  PFnumele=0;
  rho=0;

  EleSC_mass.clear();
  ele_eta.clear();
  ele_phi.clear();
  ele_pt.clear();
  ele_charge.clear();

  el_sc_eta.clear();
  el_sc_E.clear();
  el_sc_phi.clear();

  SCref.clear();
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
  eleIdLoose.clear();
  eleIdTight.clear();
  eleIdRobustLoose.clear();
  eleIdRobustTight.clear();
  CutBasedLoose.clear();
  CutBasedMedium.clear();
  CutBasedTight.clear();
  trkIsoSC.clear();
  trkIsoEle.clear();
  ElectronMVAEstimatorRun2Fall17IsoV2Values.clear();
  ElectronMVAEstimatorRun2Fall17NoIsoV2Values.clear();
  //HcalSum.clear();

  tr_pt.clear();
  tr_eta.clear();
  tr_phi.clear();
  EleTRKref.clear();
  EleTrk_mass.clear();
  numtr=0;

  float SumPt;
  int i,j,k;

  const reco::TrackCollection* tkColl = TrackHandle.product();
  math::XYZPoint pv(VertexHandle->begin()->position());
  rho = *(rhoHandle.product());


  for (auto el = electrons->begin(); el != electrons->end(); ++el)
  {
    const auto iter = electrons->ptrAt(numele);
    //std::cout << (*mvaV2IsoValues)[iter] <<" "<<(*mvaV2NoIsoValues)[iter]<< std::endl;
    std::vector<float> EleTr_separation;
    std::fill(EleTr_separation.begin(), EleTr_separation.end(), 0);
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
    eleIdLoose.push_back((*ele_id_decisions_loose)[ele]);
    eleIdTight.push_back((*ele_id_decisions_tight)[ele]);
    eleIdRobustLoose.push_back((*ele_id_decisions_robust_loose)[ele]);
    eleIdRobustTight.push_back((*ele_id_decisions_robust_tight)[ele]);
    CutBasedLoose.push_back(CutBasedLooseID(full5x5_sigmaIetaIeta[numele],dEtaSeed[numele],dPhiIn[numele],HoverE[numele],Ep[numele],relIso[numele],ExpMissInnerHits[numele],PassConversionVeto[numele],el->superCluster()->energy(),*(rhoHandle.product()),el->pt(),el->superCluster()->eta()));
    CutBasedMedium.push_back(CutBasedMediumID(full5x5_sigmaIetaIeta[numele],dEtaSeed[numele],dPhiIn[numele],HoverE[numele],Ep[numele],relIso[numele],ExpMissInnerHits[numele],PassConversionVeto[numele],el->superCluster()->energy(),*(rhoHandle.product()),el->pt(),el->superCluster()->eta()));
    CutBasedTight.push_back(CutBasedTightID(full5x5_sigmaIetaIeta[numele],dEtaSeed[numele],dPhiIn[numele],HoverE[numele],Ep[numele],relIso[numele],ExpMissInnerHits[numele],PassConversionVeto[numele],el->superCluster()->energy(),*(rhoHandle.product()),el->pt(),el->superCluster()->eta()));
    ElectronMVAEstimatorRun2Fall17NoIsoV2Values.push_back((*mvaV2NoIsoValues)[iter]);
    ElectronMVAEstimatorRun2Fall17IsoV2Values.push_back((*mvaV2IsoValues)[iter]);
    SumPt=0;
    for(auto tr = tkColl->begin(); tr != tkColl->end(); ++tr)
    {
      if(tr->pt()>5)
      {
      double dR2 = reco::deltaR2(el->eta(),el->phi(),tr->eta(),tr->phi());
      double dz = fabs(tr->vz()-el->vz());
      double dEta = fabs(tr->eta()-el->eta());
      double TrackPt = tr->pt();
      if(dR2>0 && dR2<0.4 && dz<0.1 && TrackPt>2 && dEta>0.005)
      SumPt+=TrackPt;

      EleTr_separation.push_back(dR2);
      //std::cout<<dR2<<" "<<el->eta()<<" "<<el->phi()<<" "<<tr->eta()<<" "<<tr->phi()<<std::endl;
      }
    }
    trkIsoEle.push_back(SumPt);
    numele++;

    auto tmp = std::min_element(std::begin(EleTr_separation), std::end(EleTr_separation));
    int index = std::distance(std::begin(EleTr_separation),tmp);

    if(EleTr_separation.size()>0)
    {
      if(*tmp<0.05)
      EleTRKref.push_back(index);
      else
      EleTRKref.push_back(99);
      //std::cout<<index<<std::endl;
    }

  }

  //float Sum;

  for (auto sc = superclustersEB->cbegin(); sc != superclustersEB->cend(); ++sc)
  {
    EleCounter=0;
    match=99;
    for (auto el = electrons->begin(); el != electrons->end(); ++el)
    {
      if(sc->eta()==el->superCluster()->eta() && sc->phi()==el->superCluster()->phi()) //matching
      match=EleCounter;
      EleCounter++;
    }

    SumPt=0;
    for(auto tr = tkColl->begin(); tr != tkColl->end(); ++tr)
    {
      double dR2 = reco::deltaR2(sc->eta(),sc->phi(),tr->eta(),tr->phi());
      double TrackPt = tr->pt();
      if(dR2>0 && dR2<0.4 && TrackPt>2)
      SumPt+=TrackPt;
    }

    /*Sum=0.;
    for(auto iterator=hcalRecHitsHandle->begin();iterator != hcalRecHitsHandle->end(); ++iterator)
    {
      double HitEnergy = iterator->energy();
      const auto phit = caloGeometry_.product()->getGeometry(HcalDetId(iterator->detid()))->repPos();
      double dR2 = reco::deltaR2(sc->eta(),sc->phi(),phit.eta(),phit.phi());
      if(dR2<0.15) //0<dR2<0.4
      Sum+=HitEnergy;
    }

    HcalSum.push_back(Sum);*/
    trkIsoSC.push_back(SumPt);
    SCref.push_back(match);
    sc_E.push_back(sc->energy());
    sc_eta.push_back(sc->eta());
    sc_pt.push_back(sc->energy()*sin(2*atan(exp(-1.0*sc->eta())))); //Et=pT, transverse energy is practically equal to transverse momentum
    sc_phi.push_back(sc->phi());
    numSC++;
  }


  for (auto sc = superclustersEE->cbegin(); sc != superclustersEE->cend(); ++sc)
  {
    EleCounter=0;
    match=99;

    for (auto el = electrons->begin(); el != electrons->end(); ++el)
    {
      if(sc->eta()==el->superCluster()->eta() && sc->phi()==el->superCluster()->phi()) //matching
      match=EleCounter;
      EleCounter++;
    }

    SumPt=0;
    for(auto tr = tkColl->begin(); tr != tkColl->end(); ++tr)
    {
      double dR2 = reco::deltaR2(sc->eta(),sc->phi(),tr->eta(),tr->phi());
      double TrackPt = tr->pt();
      if(dR2>0 && dR2<0.4 && TrackPt>2)
      SumPt+=TrackPt;
    }

    /*Sum=0.;
    for(auto iterator=hcalRecHitsHandle->begin();iterator != hcalRecHitsHandle->end(); ++iterator)
    {
      double HitEnergy = iterator->energy();
      const auto phit = caloGeometry_.product()->getGeometry(HcalDetId(iterator->detid()))->repPos();
      double dR2 = reco::deltaR2(sc->eta(),sc->phi(),phit.eta(),phit.phi());
      if(dR2<0.15) //0<dR2<0.4
      Sum+=HitEnergy;
    }
    HcalSum.push_back(Sum);*/
    trkIsoSC.push_back(SumPt);
    SCref.push_back(match);
    sc_E.push_back(sc->energy());
    sc_eta.push_back(sc->eta());
    sc_pt.push_back(sc->energy()*sin(2*atan(exp(-1.0*sc->eta())))); //Et=pT, transverse energy is practically equal to transverse momentum
    sc_phi.push_back(sc->phi());
    numSC++;
  }


  //const HBHERecHitCollection* hithbhe_ = hcalRecHitsHandle.product();

  //typename HBHERecHitCollection::const_iterator i = HBHEhits->begin();
  /*for(auto iterator=hcalRecHitsHandle->begin();iterator != hcalRecHitsHandle->end(); ++iterator)
  {
    std::cout<<HcalDetId(iterator->detid()).ieta()<<std::endl;
    const auto phit = caloGeometry_.product()->getGeometry(HcalDetId(iterator->detid()))->repPos();
    std::cout<<phit.eta()<<std::endl;
  }*/


  for(auto tr = tkColl->begin(); tr != tkColl->end(); ++tr)
  {
    //track electron matching
    if(tr->pt()>5)
    {
    tr_pt.push_back(tr->pt());
    tr_eta.push_back(tr->eta());
    tr_phi.push_back(tr->phi());
    numtr++;
    }
  }

  //HoECalculator hoeCalc(caloGeometry_);

  int PassSCcut=0,PassElecut=0;

  if(numele==1 || (numele==2 && ele_charge[0]==-1*ele_charge[1]))
  {
    for(i=0;i<numele;i++)
    {
      if(ele_pt[i]>5)
      PassElecut++;
      for(j=0;j<numSC;j++)
      {
        if(SCref[j]==99)
        EleSC_mass.push_back(InvariantMass(ele_pt[i],ele_eta[i],ele_phi[i],sc_pt[j],sc_eta[j],sc_phi[j]));
        if(sc_pt[j]>5)
        PassSCcut++;
      }
      /*for(k=0;k<numtr;k++)
      {
        if(EleTRKref[i]!=k) //trk not matched el
        EleTrk_mass.push_back(InvariantMass(ele_pt[i],ele_eta[i],ele_phi[i],tr_pt[k],tr_eta[k],tr_phi[k]));
      }*/
    }
    for(k=0;k<numtr;k++)
    {
      if(numele==1)
      {
        if(EleTRKref[0]!=k) //trk not matched el
        EleTrk_mass.push_back(InvariantMass(ele_pt[0],ele_eta[0],ele_phi[0],tr_pt[k],tr_eta[k],tr_phi[k]));
      }
      if(numele==2)
      {
        if(EleTRKref[0]!=k && EleTRKref[1]!=k) //trk not matched el
        {
          EleTrk_mass.push_back(InvariantMass(ele_pt[0],ele_eta[0],ele_phi[0],tr_pt[k],tr_eta[k],tr_phi[k]));
          EleTrk_mass.push_back(InvariantMass(ele_pt[1],ele_eta[1],ele_phi[1],tr_pt[k],tr_eta[k],tr_phi[k]));
        }
      }
    }

    bool NotMatchedTrack = std::all_of(EleTRKref.cbegin(), EleTRKref.cend(), [](int TRreference){ return TRreference==99; });
    bool NotMatchedSC = std::all_of(SCref.cbegin(), SCref.cend(), [](int SCreference){ return SCreference==99; });

    if(PassSCcut==numSC*numele && PassElecut==numele && !NotMatchedTrack && !NotMatchedSC) //dont save the event if electron has 0 matchs with tracks and SCs
    reco_tree->Fill();

  }

}


void RecoAnalyzer::beginJob() {

}

void RecoAnalyzer::endJob() {

}

void RecoAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);

}

DEFINE_FWK_MODULE(RecoAnalyzer);
