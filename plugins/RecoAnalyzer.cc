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

  TTree *reco_tree;
  std::vector<float> ele_pt,ele_eta,ele_phi,full5x5_sigmaIetaIeta,dEtaSeed,dPhiIn,HoverE,relIso,Ep,
  el_sc_eta, el_sc_E,el_sc_phi,sc_eta, sc_pt, sc_phi, sc_E, eleIdLoose,eleIdTight,eleIdRobustLoose,eleIdRobustTight;
  std::vector<int> ele_charge,SCref,ExpMissInnerHits;
  std::vector<bool> PassConversionVeto,CutBasedLoose,CutBasedMedium,CutBasedTight;
  float EleSC_mass,SCSC_mass,rho;
  int numele, PFnumele,numSC,EleCounter,match;
  TLorentzVector P,P0,P1,p,p0,p1;

  unsigned long long cachedCaloGeometryID_;
  edm::ESHandle<CaloGeometry> caloGeometry_;


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
hcalRecHitsInputHBHE_(iConfig.getUntrackedParameter<edm::InputTag>("HBHERecHit"))
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

  edm::Service<TFileService> fs;
  reco_tree = fs->make<TTree>("Events", "Events");

  reco_tree->Branch("numele",&numele);
  reco_tree->Branch("PFnumele",&PFnumele);
  reco_tree->Branch("ele_eta",&ele_eta);
  reco_tree->Branch("ele_phi",&ele_phi);
  reco_tree->Branch("EleSC_mass",&EleSC_mass);
  reco_tree->Branch("SCSC_mass",&SCSC_mass);

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
  //reco_tree->Branch("ElectronMVAEstimatorRun2Fall17IsoV2Values",&ElectronMVAEstimatorRun2Fall17IsoV2Values);
  //reco_tree->Branch("ElectronMVAEstimatorRun2Fall17NoIsoV2Values",&ElectronMVAEstimatorRun2Fall17NoIsoV2Values);

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

  cachedCaloGeometryID_ = 0;

}

RecoAnalyzer::~RecoAnalyzer() {

}

void RecoAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  using namespace edm;

  if (cachedCaloGeometryID_ != iSetup.get<CaloGeometryRecord>().cacheIdentifier()) {
        cachedCaloGeometryID_ = iSetup.get<CaloGeometryRecord>().cacheIdentifier();
        iSetup.get<CaloGeometryRecord>().get(caloGeometry_);
    }

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

  numSC=0;
  numele=0;
  PFnumele=0;
  EleSC_mass=0;
  SCSC_mass=0;
  rho=0;
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


  for (auto sc = superclustersEB->cbegin(); sc != superclustersEB->cend(); ++sc)
  {
    EleCounter=0;
    match=99;
    for (auto it = electrons->begin(); it != electrons->end(); ++it)
    {
      if(sc->eta()==it->superCluster()->eta() && sc->phi()==it->superCluster()->phi()) //matching
      match=EleCounter;
      EleCounter++;
    }
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

    for (auto it = electrons->begin(); it != electrons->end(); ++it)
    {
      if(sc->eta()==it->superCluster()->eta() && sc->phi()==it->superCluster()->phi()) //matching
      match=EleCounter;
      EleCounter++;
    }
    SCref.push_back(match);
    sc_E.push_back(sc->energy());
    sc_eta.push_back(sc->eta());
    sc_pt.push_back(sc->energy()*sin(2*atan(exp(-1.0*sc->eta())))); //Et=pT, transverse energy is practically equal to transverse momentum
    sc_phi.push_back(sc->phi());
    numSC++;
  }

  rho = *(rhoHandle.product());
  const HBHERecHitCollection* hithbhe_ = hcalRecHitsHandle.product();

  for (auto it = electrons->begin(); it != electrons->end(); ++it)
  {
    ele_pt.push_back(it->pt());
    ele_eta.push_back(it->eta());
    ele_phi.push_back(it->phi());
    el_sc_eta.push_back(it->superCluster()->eta());
    el_sc_phi.push_back(it->superCluster()->phi());
    el_sc_E.push_back(it->superCluster()->energy());
    ele_charge.push_back(it->charge());
    full5x5_sigmaIetaIeta.push_back(it->full5x5_sigmaIetaIeta());
    dEtaSeed.push_back(it->superCluster().isNonnull() && it->superCluster()->seed().isNonnull() ?
    it->deltaEtaSuperClusterTrackAtVtx() - it->superCluster()->eta() + it->superCluster()->seed()->eta() : std::numeric_limits<float>::max());
    dPhiIn.push_back(it->deltaPhiSuperClusterTrackAtVtx());
    HoverE.push_back(it->hadronicOverEm());
    const float ecal_energy_inverse = 1.0/it->ecalEnergy();
    const float eSCoverP = it->eSuperClusterOverP();
    Ep.push_back(std::abs(1.0 - eSCoverP)*ecal_energy_inverse);
    relIso.push_back((float)(it->pfIsolationVariables().sumChargedHadronPt+std::max(float(0.0),it->pfIsolationVariables().sumNeutralHadronEt+it->pfIsolationVariables().sumPhotonEt))/it->pt());
    constexpr auto missingHitType = reco::HitPattern::MISSING_INNER_HITS;
    ExpMissInnerHits.push_back(it->gsfTrack()->hitPattern().numberOfLostHits(missingHitType));
    PassConversionVeto.push_back(!ConversionTools::hasMatchedConversion(*it, *hConversions, beamspot.position()));
    const auto el = electrons->ptrAt(numele);
    eleIdLoose.push_back((*ele_id_decisions_loose)[el]);
    eleIdTight.push_back((*ele_id_decisions_tight)[el]);
    eleIdRobustLoose.push_back((*ele_id_decisions_robust_loose)[el]);
    eleIdRobustTight.push_back((*ele_id_decisions_robust_tight)[el]);
    CutBasedLoose.push_back(CutBasedLooseID(full5x5_sigmaIetaIeta[numele],dEtaSeed[numele],dPhiIn[numele],HoverE[numele],Ep[numele],relIso[numele],ExpMissInnerHits[numele],PassConversionVeto[numele],it->superCluster()->energy(),*(rhoHandle.product()),it->pt(),it->superCluster()->eta()));
    CutBasedMedium.push_back(CutBasedMediumID(full5x5_sigmaIetaIeta[numele],dEtaSeed[numele],dPhiIn[numele],HoverE[numele],Ep[numele],relIso[numele],ExpMissInnerHits[numele],PassConversionVeto[numele],it->superCluster()->energy(),*(rhoHandle.product()),it->pt(),it->superCluster()->eta()));
    CutBasedTight.push_back(CutBasedTightID(full5x5_sigmaIetaIeta[numele],dEtaSeed[numele],dPhiIn[numele],HoverE[numele],Ep[numele],relIso[numele],ExpMissInnerHits[numele],PassConversionVeto[numele],it->superCluster()->energy(),*(rhoHandle.product()),it->pt(),it->superCluster()->eta()));
    numele++;
  }

  //typename HBHERecHitCollection::const_iterator i = HBHEhits->begin();
  for(auto iterator=hcalRecHitsHandle->begin();iterator != hcalRecHitsHandle->end(); ++iterator)
  std::cout<<HcalDetId(iterator->detid()).ieta()<<std::endl;

  HoECalculator hoeCalc(caloGeometry_);

  if((numele==1 && numSC==2 && ele_pt[0]>5 && sc_pt[0]>5 && sc_pt[1]>5) || (numele==2 && numSC==2 && ele_pt[0]>5 && ele_pt[1]>5 && sc_pt[0]>5 && sc_pt[1]>5 && ele_charge[0]==-1*ele_charge[1]))
  {
    P0.SetPtEtaPhiM(ele_pt[0],ele_eta[0],ele_phi[0],0);
    if(SCref[0]==0)
    P1.SetPtEtaPhiM(sc_pt[1],sc_eta[1],sc_phi[1],0);
    else if(SCref[0]==1 || SCref[0]==99)
    P1.SetPtEtaPhiM(sc_pt[0],sc_eta[0],sc_phi[0],0);
    P=P0+P1;
    EleSC_mass=P.M();
    p0.SetPtEtaPhiM(sc_pt[0],sc_eta[0],sc_phi[0],0);
    p1.SetPtEtaPhiM(sc_pt[1],sc_eta[1],sc_phi[1],0);
    p=p0+p1;
    SCSC_mass=p.M();
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
