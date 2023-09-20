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

#include "DataFormats/METReco/interface/PFMET.h"
#include "DataFormats/JetReco/interface/PFJet.h"
// vector<reco::PFJet> ak4PFJets

#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include <FWCore/ParameterSet/interface/FileInPath.h>

#include "RecoEgamma/EgammaTools/interface/ConversionTools.h"
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

using std::cout;
using std::endl;


class SosSkimAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources> {
public:
  explicit SosSkimAnalyzer(const edm::ParameterSet&);
  ~SosSkimAnalyzer();

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  void beginJob() override;
  void analyze(const edm::Event&, const edm::EventSetup&) override;
  void endJob() override;

  edm::EDGetToken electronsToken_;
  edm::EDGetToken lowPtElectronsToken_;
  edm::EDGetTokenT<reco::TrackCollection> TrackToken;
  edm::EDGetTokenT<edm::ValueMap<float> > cutvValuesMapToken;
  edm::EDGetTokenT<edm::ValueMap<float> > cutlValuesMapToken;
  edm::EDGetTokenT<edm::ValueMap<float> > cutmValuesMapToken;
  edm::EDGetTokenT<edm::ValueMap<float> > cuttValuesMapToken;
  edm::EDGetTokenT<edm::ValueMap<float> > lowptValuesMapToken;
  edm::EDGetTokenT<edm::ValueMap<float> > mvaV2IsoValuesMapToken;
  edm::EDGetTokenT<edm::ValueMap<float> > mvaV2NoIsoValuesMapToken;
  edm::EDGetTokenT<reco::VertexCollection> PrimaryVertexToken;
  edm::EDGetTokenT<reco::VertexCollection> SecondaryVertexToken;
  edm::EDGetTokenT<std::vector<reco::PFMET> > MetToken;
  edm::EDGetTokenT<std::vector<reco::PFJet> > JetToken;

  edm::InputTag Electron_;
  edm::InputTag LowPtElectron_;
  edm::InputTag TrackInputTag_;
  edm::InputTag cutvValuesMapInputTag_;
  edm::InputTag cutlValuesMapInputTag_;
  edm::InputTag cutmValuesMapInputTag_;
  edm::InputTag cuttValuesMapInputTag_;
  edm::InputTag lowptValuesMapInputTag_;
  edm::InputTag mvaV2IsoValuesMapInputTag_;
  edm::InputTag mvaV2NoIsoValuesMapInputTag_;
  edm::InputTag PrimaryVertexInputTag_;
  edm::InputTag SecondaryVertexInputTag_;
  edm::InputTag MetInputTag_;
  edm::InputTag JetInputTag_;

  TTree *reco_tree;
  // std::vector<float> ele_pt,ele_eta,ele_phi,ElectronMVAEstimatorRun2Fall17IsoV2Values,
  // ElectronMVAEstimatorRun2Fall17NoIsoV2Values,ele_z,ele_x,ele_y,SVchi2,SVndof, SVpt,SVeta,SVphi,SVcharge,SVmass,SVx,SVy,SVz;
  // std::vector<int> ele_charge,EleSVtrkref,EleTrkIsNoNull,sv_tr_size,SVTrEleRef;
  // int numele,numGoodSV;
  // int evtnum;

  // j/psi(ee) cand pair
  float TnP_dr, TnP_mass;
  // sv B(Kee) cand triplet
  float sv_mass, sv_chi2, sv_prob, sv_ndof, sv_charge, sv_pt, sv_eta, sv_phi, sv_x, sv_y, sv_z;
  int sv_nTrack;
  // tag
  float Tag_pt, Tag_eta, Tag_phi, Tag_cutBased, Tag_mvaFall17V2noIso, Tag_mvaFall17V2Iso, Tag_lowPtID;
  int Tag_isLowPt, Tag_charge;
  // probe
  float Probe_pt, Probe_eta, Probe_phi;
  int Probe_charge, Probe_matches_GED, Probe_matches_LowPt;
  float Probe_cutBased, Probe_mvaFall17V2noIso, Probe_mvaFall17V2Iso, Probe_lowPtID;
  // per event
  float TnP_ht, TnP_met, TnP_ipair;
  float K_pt, K_eta, K_phi;

    //
    float pv_x, pv_y, pv_z;
    float K_x, K_y, K_z;
    float Tag_x, Tag_y, Tag_z;
    float Probe_x, Probe_y, Probe_z;
    float dr_K_Probe;
    float dr_K_Tag;
};

SosSkimAnalyzer::SosSkimAnalyzer(const edm::ParameterSet& iConfig):
Electron_(iConfig.getUntrackedParameter<edm::InputTag>("Electron")),
LowPtElectron_(iConfig.getUntrackedParameter<edm::InputTag>("LowPtElectron")),
TrackInputTag_(iConfig.getUntrackedParameter<edm::InputTag>("Track")),
cutvValuesMapInputTag_(iConfig.getUntrackedParameter<edm::InputTag>("cutV")),
cutlValuesMapInputTag_(iConfig.getUntrackedParameter<edm::InputTag>("cutL")),
cutmValuesMapInputTag_(iConfig.getUntrackedParameter<edm::InputTag>("cutM")),
cuttValuesMapInputTag_(iConfig.getUntrackedParameter<edm::InputTag>("cutT")),
lowptValuesMapInputTag_(iConfig.getUntrackedParameter<edm::InputTag>("lowptValuesMap")),
mvaV2IsoValuesMapInputTag_(iConfig.getUntrackedParameter<edm::InputTag>("mvaV2IsoValuesMap")),
mvaV2NoIsoValuesMapInputTag_(iConfig.getUntrackedParameter<edm::InputTag>("mvaV2NoIsoValuesMap")),
PrimaryVertexInputTag_(iConfig.getUntrackedParameter<edm::InputTag>("primaryVertices")),
SecondaryVertexInputTag_(iConfig.getUntrackedParameter<edm::InputTag>("secondaryVertices")),
MetInputTag_(iConfig.getUntrackedParameter<edm::InputTag>("met")),
JetInputTag_(iConfig.getUntrackedParameter<edm::InputTag>("jets"))
{
  electronsToken_    = consumes<edm::View<reco::GsfElectron> >(Electron_);
  lowPtElectronsToken_    = consumes<edm::View<reco::GsfElectron> >(LowPtElectron_);
  TrackToken = consumes<reco::TrackCollection>(TrackInputTag_);
  cutvValuesMapToken = consumes<edm::ValueMap<float>>(cutvValuesMapInputTag_);
  cutlValuesMapToken = consumes<edm::ValueMap<float>>(cutlValuesMapInputTag_);
  cutmValuesMapToken = consumes<edm::ValueMap<float>>(cutmValuesMapInputTag_);
  cuttValuesMapToken = consumes<edm::ValueMap<float>>(cuttValuesMapInputTag_);
  lowptValuesMapToken = consumes<edm::ValueMap<float>>(lowptValuesMapInputTag_);
  mvaV2IsoValuesMapToken = consumes<edm::ValueMap<float>>(mvaV2IsoValuesMapInputTag_);
  mvaV2NoIsoValuesMapToken = consumes<edm::ValueMap<float>>(mvaV2NoIsoValuesMapInputTag_);
  PrimaryVertexToken = consumes<reco::VertexCollection>(PrimaryVertexInputTag_);
  SecondaryVertexToken = consumes<reco::VertexCollection>(SecondaryVertexInputTag_);
  MetToken = consumes<std::vector<reco::PFMET> >(MetInputTag_);
  JetToken = consumes<std::vector<reco::PFJet> >(JetInputTag_);

  edm::Service<TFileService> fs;
  reco_tree = fs->make<TTree>("Events", "Events");

  reco_tree->Branch("TnP_dr",&TnP_dr);
  reco_tree->Branch("TnP_mass",&TnP_mass);
  reco_tree->Branch("sv_mass",&sv_mass);
  reco_tree->Branch("sv_chi2",&sv_chi2);
  reco_tree->Branch("sv_prob",&sv_prob);
  reco_tree->Branch("sv_ndof",&sv_ndof);
  reco_tree->Branch("sv_charge",&sv_charge);
  reco_tree->Branch("sv_pt",&sv_pt);
  reco_tree->Branch("sv_eta",&sv_eta);
  reco_tree->Branch("sv_phi",&sv_phi);
  reco_tree->Branch("sv_nTrack",&sv_nTrack);
  reco_tree->Branch("sv_x",&sv_x);
  reco_tree->Branch("sv_y",&sv_y);
  reco_tree->Branch("sv_z",&sv_z);
  reco_tree->Branch("Tag_pt",&Tag_pt);
  reco_tree->Branch("Tag_eta",&Tag_eta);
  reco_tree->Branch("Tag_phi",&Tag_phi);
  reco_tree->Branch("Tag_cutBased",&Tag_cutBased);
  reco_tree->Branch("Tag_mvaFall17V2noIso",&Tag_mvaFall17V2noIso);
  reco_tree->Branch("Tag_mvaFall17V2Iso",&Tag_mvaFall17V2Iso);
  reco_tree->Branch("Tag_lowPtID",&Tag_lowPtID);
  reco_tree->Branch("Tag_isLowPt",&Tag_isLowPt);
  reco_tree->Branch("Tag_charge",&Tag_charge);
  reco_tree->Branch("Probe_pt",&Probe_pt);
  reco_tree->Branch("Probe_eta",&Probe_eta);
  reco_tree->Branch("Probe_phi",&Probe_phi);
  reco_tree->Branch("Probe_cutBased",&Probe_cutBased);
  reco_tree->Branch("Probe_mvaFall17V2noIso",&Probe_mvaFall17V2noIso);
  reco_tree->Branch("Probe_mvaFall17V2Iso",&Probe_mvaFall17V2Iso);
  reco_tree->Branch("Probe_lowPtID",&Probe_lowPtID);
  reco_tree->Branch("Probe_charge",&Probe_charge);
  reco_tree->Branch("Probe_matches_GED",&Probe_matches_GED);
  reco_tree->Branch("Probe_matches_LowPt",&Probe_matches_LowPt);
  reco_tree->Branch("K_pt",&K_pt);
  reco_tree->Branch("K_eta",&K_eta);
  reco_tree->Branch("K_phi",&K_phi);
  reco_tree->Branch("TnP_ht",&TnP_ht);
  reco_tree->Branch("TnP_met",&TnP_met);
  reco_tree->Branch("TnP_ipair",&TnP_ipair);

  reco_tree->Branch("Probe_x",&Probe_x);
  reco_tree->Branch("Probe_y",&Probe_y);
  reco_tree->Branch("Probe_z",&Probe_z);
  reco_tree->Branch("Tag_x",&Tag_x);
  reco_tree->Branch("Tag_y",&Tag_y);
  reco_tree->Branch("Tag_z",&Tag_z);
  reco_tree->Branch("pv_x",&pv_x);
  reco_tree->Branch("pv_y",&pv_y);
  reco_tree->Branch("pv_z",&pv_z);
  reco_tree->Branch("K_x",&K_x);
  reco_tree->Branch("K_y",&K_y);
  reco_tree->Branch("K_z",&K_z);
  reco_tree->Branch("dr_K_Probe",&dr_K_Probe);
  reco_tree->Branch("dr_K_Tag",&dr_K_Tag);

  //Event number
  // reco_tree->Branch("evtnum",&evtnum);

}

SosSkimAnalyzer::~SosSkimAnalyzer() {

}

bool tracksMatch(const reco::Track a, const reco::Track b){
    return reco::deltaR2(a.eta(),a.phi(),b.eta(),b.phi()) < 0.0001;
}
float getMass(const reco::GsfElectron a, const reco::Track b){
    TLorentzVector ta; ta.SetPtEtaPhiM(a.pt(),a.eta(),a.phi(),0);
    TLorentzVector tb; tb.SetPtEtaPhiM(b.pt(),b.eta(),b.phi(),0);
    return (ta+tb).M();
}

void SosSkimAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  using namespace edm;

  edm::Handle<edm::View<reco::GsfElectron> > electrons;
  iEvent.getByToken(electronsToken_, electrons);

  edm::Handle<edm::View<reco::GsfElectron> > lowPtElectrons;
  iEvent.getByToken(lowPtElectronsToken_, lowPtElectrons);

  // std::cout << " Ele: " << electrons->size() << std::endl;
  // std::cout << " Low: " << lowPtElectrons->size() << std::endl;

  edm::Handle<reco::TrackCollection> TrackHandle;
  iEvent.getByToken(TrackToken, TrackHandle);

  edm::Handle<edm::ValueMap<float> > cutvValues;
  iEvent.getByToken(cutvValuesMapToken,cutvValues);
  edm::Handle<edm::ValueMap<float> > cutlValues;
  iEvent.getByToken(cutlValuesMapToken,cutlValues);
  edm::Handle<edm::ValueMap<float> > cutmValues;
  iEvent.getByToken(cutmValuesMapToken,cutmValues);
  edm::Handle<edm::ValueMap<float> > cuttValues;
  iEvent.getByToken(cuttValuesMapToken,cuttValues);
  edm::Handle<edm::ValueMap<float> > lowptValues;
  iEvent.getByToken(lowptValuesMapToken,lowptValues);

  // cout << "n electrons: " << lowPtElectrons->size() << endl;
  // cout << "n values: " << lowptValues->size() << endl;
  // for(unsigned int i=0;i<lowPtElectrons->size();i++){
  //     //cout << " ele " << (lowPtElectrons->refAt(i)).id() << endl;
  //     ///////// cout << " ele " << (*lowPtElectrons->ptrAt(i)).id() << endl;
  //     const auto iter = lowPtElectrons->ptrAt(i);
  //     const auto val = (*lowptValues)[iter];
  //     cout << " ele " << iter->pt() << " and " << val << endl;      
  // }

  // auto val = lowptValues->begin();
  // for(unsigned int i=0;i<lowptValues->size();i++){
  //     // cout << " val " << (*lowptValues)->ptrAt(i)).id() << endl;
  //     // cout << " val " << val.id() << endl;
  //     val++;
  // }
                      //     Tag_mvaFall17V2Iso = (*mvaV2IsoValues)[electrons->ptrAt(tag_indices.first)];
                      // } else {
                      //     Tag_lowPtID = (*lowptValues)[lowPtElectrons->ptrAt(tag_indices.first)];

  edm::Handle<edm::ValueMap<float> > mvaV2NoIsoValues;
  iEvent.getByToken(mvaV2NoIsoValuesMapToken,mvaV2NoIsoValues);

  edm::Handle<edm::ValueMap<float> > mvaV2IsoValues;
  iEvent.getByToken(mvaV2IsoValuesMapToken,mvaV2IsoValues);

  edm::Handle<reco::VertexCollection> primaryVertexHandle;
  iEvent.getByToken(PrimaryVertexToken,primaryVertexHandle);
  edm::Handle<reco::VertexCollection> secondaryVertexHandle;
  iEvent.getByToken(SecondaryVertexToken,secondaryVertexHandle);

  edm::Handle<std::vector<reco::PFMET> > metHandle;
  iEvent.getByToken(MetToken,metHandle);

  edm::Handle<std::vector<reco::PFJet> > jetHandle;
  iEvent.getByToken(JetToken,jetHandle);

  // edm::ValueMap<float>                  "lowPtGsfElectronID"        ""                "RECO"
  // for(auto el : *electrons){
  //     std::cout << el.electronID("cutBasedElectronID-Fall17-94X-V2-veto") << std::endl;
  // }

  // must have at least one electron (GED or lowpt) and an SV
  if( primaryVertexHandle->size()==0 || secondaryVertexHandle->size()==0 || (electrons->size()==0 && lowPtElectrons->size()==0) ) return;

  // evtnum=0;
  TnP_dr = 0;
  TnP_mass = 0;
  sv_mass = 0;
  sv_chi2 = 0;
  sv_prob = 0;
  sv_ndof = 0;
  sv_charge = -99;
  sv_pt = 0;
  sv_eta = 0;
  sv_phi = 0;
  sv_nTrack = -1;
  sv_x = 0;
  sv_y = 0;
  sv_z = 0;
  Tag_pt = 0;
  Tag_eta = 0;
  Tag_phi = 0;
  Tag_cutBased = -99;
  Tag_mvaFall17V2noIso = -99;
  Tag_mvaFall17V2Iso = -99;
  Tag_lowPtID = -99;
  Tag_isLowPt = -1;
  Tag_charge = -99;
  Probe_pt = 0;
  Probe_eta = 0;
  Probe_phi = 0;
  K_pt = 0;
  K_eta = 0;
  K_phi = 0;
  Probe_cutBased = -99;
  Probe_mvaFall17V2noIso = -99;
  Probe_mvaFall17V2Iso = -99;
  Probe_lowPtID = -99;
  Probe_charge = -99;
  Probe_matches_GED = 0;
  Probe_matches_LowPt = 0;
  TnP_ht = -1;
  TnP_met = -1;
  TnP_ipair = -1;
  pv_x = 0;
  pv_y = 0;
  pv_z = 0;
  K_x = 0;
  K_y = 0;
  K_z = 0;
  Tag_x = 0;
  Tag_y = 0;
  Tag_z = 0;
  Probe_x = 0;
  Probe_y = 0;
  Probe_z = 0;
  dr_K_Probe = 0;
  dr_K_Tag = 0;

  //helper variables
  TLorentzVector P_sum,P_tmp;

  const reco::TrackCollection* tkColl = TrackHandle.product();
  const reco::VertexCollection* svColl = secondaryVertexHandle.product();
  // bool debug=true;
  bool debug=false;

  if(debug)cout<<"== New Event =="<<endl;

  std::map<unsigned int, std::vector<unsigned int>> sv_tracks{};
  std::map<unsigned int, float> sv_massCalc{};
  std::map<unsigned int, int> sv_chargeCalc{};
  for(unsigned int isv=0; isv < svColl->size(); isv++){
      const auto& sv = svColl->at(isv);
      int svCharge=0.;
      P_sum.SetPtEtaPhiM(0.,0.,0.,0.);
      // int SVtrackSize = 0;
      std::vector<unsigned int> track_indices{};
      // search for suitable SV candidates
      if(debug) std::cout << "Checking vertex " << isv << " with " << sv.tracksSize() << " tracks (chi2="<<sv.chi2()<<"), ";
      for(unsigned int r=0;r<sv.tracksSize();r++){
          const auto & tk = sv.trackRefAt(r);
          if(tk->pt()>0.7){
              svCharge+=tk->charge();
              P_tmp.SetPtEtaPhiM(tk->pt(),tk->eta(),tk->phi(),0.);
              P_sum+=P_tmp;
              track_indices.push_back(r);
          }
      }
      if(debug) std::cout << track_indices.size() << " total, with charge " << svCharge << " and mass " << P_sum.M();
      if(P_sum.M()>4.0 && track_indices.size()==3 && std::abs(svCharge)==1){ //fill only if satisfied
      // if(track_indices.size()==3 && std::abs(svCharge)==1){ //fill only if satisfied
          // good_sv_idx.push_back(isv);
          sv_tracks[isv] = track_indices;
          sv_massCalc[isv] = P_sum.M();
          sv_chargeCalc[isv] = svCharge;
          if(debug){
              std::cout << " ITS GOOD! Tk idxs:";
              for(auto idx : track_indices) cout << " " << idx;
              if(debug) std::cout << std::endl;

              for(unsigned int r=0;r<sv.tracksSize();r++){
                  const auto & tk = sv.trackRefAt(r);
                  if(tk->pt()>0.7){                      
                      if(debug) cout<<" Good Tk " << r << " with pt eta phi " << tk->pt() << " " << tk->eta() << " " << tk->phi() << endl;
                  }
              }
          }
      } else if(debug){
          std::cout << std::endl;
      }
  }
  if(sv_tracks.size()==0) return;
  
  int iPair=0;
  float ht25=-1;
  for(const auto& p : sv_tracks){
      const auto& sv = svColl->at(p.first);
      if(debug) cout<<"Considering vertex "<<p.first<<"(chi2="<<sv.chi2()<<") for ele matches"<<endl;;
      // check if any SV tracks match to a TAG electron
      std::vector<std::pair<unsigned int, unsigned int>> ged_tk_tag_indices{};
      std::vector<std::pair<unsigned int, unsigned int>> lowpt_tk_tag_indices{};
      for(const unsigned int track_index : p.second){
          const auto& tk = sv.trackRefAt(track_index);
          if(debug) cout<<" Tk " << track_index << " with pt eta phi " << tk->pt() << " " << tk->eta() << " " << tk->phi() << endl;
          // check for a GED TAG
          for(unsigned int i=0; i < electrons->size(); i++){
              const auto& el = electrons->at(i);
              if(debug) cout<<" - GED " << i << " with pt eta phi " << el.pt() << " " << el.eta() << " " << el.phi();
              if(!el.closestCtfTrackRef().isNonnull()){
                  if(debug) cout<<" noref" << endl;
                  continue;
              }
              const auto &elTkRef = tkColl->begin() + el.closestCtfTrackRef().index();
              if(debug) cout<<" (tkRef " << elTkRef->pt() << " " << elTkRef->eta() << " " << elTkRef->phi()<<")";
              if(tracksMatch(*elTkRef, *tk)) {
                  ged_tk_tag_indices.push_back(std::make_pair(i,track_index));
                  if(debug) cout<<" is match ";
              }
              if(debug) cout<<endl;
          }
          // check for a LowPt TAG
          for(unsigned int i=0; i < lowPtElectrons->size(); i++){
              const auto& el = lowPtElectrons->at(i);
              if(debug) cout<<" - LOW " << i << " with pt eta phi " << el.pt() << " " << el.eta() << " " << el.phi();
              if(!el.closestCtfTrackRef().isNonnull()){
                  if(debug) cout<<" noref" << endl;
                  continue;
              }
              const auto &elTkRef = tkColl->begin() + el.closestCtfTrackRef().index();
              if(debug) cout<<" (tkRef " << elTkRef->pt() << " " << elTkRef->eta() << " " << elTkRef->phi()<<")";
              if(tracksMatch(*elTkRef, *tk)) {
                  lowpt_tk_tag_indices.push_back(std::make_pair(i,track_index));
                  if(debug) cout<<" is match ";
              }
              if(debug) cout<<endl;
          }
      }
      // Look for j/psi(ee) candidates among the other tracks
      for(bool gedTag : {0,1}){
          const auto& el_tk_tag_indices = (gedTag ? ged_tk_tag_indices : lowpt_tk_tag_indices);
          for(const auto tag_indices : el_tk_tag_indices){
              if(gedTag && debug)  cout << " > Checking GED Tags"<<endl;
              if(!gedTag && debug) cout << " > Checking LOW Tags"<<endl;
              const auto &eleColl = (gedTag ? electrons : lowPtElectrons);
              const auto& tagEl = eleColl->at(tag_indices.first);
              const unsigned int iTagTrack = tag_indices.second;
              if(debug)cout<<"  Tag El(" << tag_indices.first << ") Tk(" << iTagTrack<<")" <<endl;
              for(const unsigned int iProbeTrack : p.second){
                  if (iTagTrack==iProbeTrack) continue;
                  const auto& probeTrack = sv.trackRefAt(iProbeTrack);
                  //  check for a j/psi candidate
                  float m = getMass(tagEl, *probeTrack);
                  if (m > 2.0 && m < 4.5 && (tagEl.charge()+probeTrack->charge()==0)){ // FIXME EVENTUALLY
                  // if (m > 0.0 && m < 400.5 && (tagEl->charge()+probeTrack->charge()==0)){ // FIXME EVENTUALLY
                      if(debug)cout<<"  - makes jpsi (m = "<<m<<") with track " << iProbeTrack << endl;
                      // probe track matches GED
                      int iGed=-1;
                      float min_dr=9e9;                      
                      for(unsigned int i=0; i < electrons->size(); i++){
                          const auto& el = electrons->at(i);
                          float dr=9e9;
                          if ( el.closestCtfTrackRef().isNonnull() ){
                              const auto &elTkRef = tkColl->begin() + el.closestCtfTrackRef().index();
                              dr = reco::deltaR(elTkRef->eta(),elTkRef->phi(),probeTrack->eta(),probeTrack->phi());
                          } else {
                              dr = reco::deltaR(el.eta(),el.phi(),probeTrack->eta(),probeTrack->phi());
                          }
                          if (dr < 0.01 && dr < min_dr){
                              min_dr = dr;
                              iGed=i;
                              if(debug)cout<<"    with GED match (dr = "<<dr<<") Ele Idx "<<iGed<< endl;
                              //if(debug)cout<<"    with GED match (dr = "<<dr<<", from " << elTkRef->eta() << " " << elTkRef->phi() << " " << probeTrack->eta() << " " << probeTrack->phi() << ") Ele Idx "<<iGed<< endl;
                          }
                      }
                      // probe track matches LowPt
                      int iLow=-1;
                      min_dr=9e9;
                      for(unsigned int i=0; i < lowPtElectrons->size(); i++){
                          const auto& el = lowPtElectrons->at(i);
                          // const auto &elTkRef = tkColl->begin() + el.closestCtfTrackRef().index();
                          float dr=9e9;
                          if ( el.closestCtfTrackRef().isNonnull() ){
                              const auto &elTkRef = tkColl->begin() + el.closestCtfTrackRef().index();
                              dr = reco::deltaR(elTkRef->eta(),elTkRef->phi(),probeTrack->eta(),probeTrack->phi());
                          } else {
                              dr = reco::deltaR(el.eta(),el.phi(),probeTrack->eta(),probeTrack->phi());
                          }
                          // float dr = reco::deltaR(elTkRef->eta(),elTkRef->phi(),probeTrack->eta(),probeTrack->phi());
                          if (dr < 0.01 && dr < min_dr){
                              min_dr = dr;
                              iLow=i;
                              if(debug)cout<<"    with LOW match (dr = "<<dr<<") Ele Idx "<<iLow<< endl;
                              //if(debug)cout<<"    with LOW match (dr = "<<dr<<", from " << elTkRef->eta() << " " << elTkRef->phi() << " " << probeTrack->eta() << " " << probeTrack->phi() << ") Ele Idx "<<iLow<< endl;
                          }
                      }
                      
                      if(debug)cout<<" ~~ Recording m2 = "<<TnP_mass<<", m3 = "<<sv_massCalc[p.first]<<" with matches GED "<< iGed << " and LOW " << iLow << endl;

                      // fill tree
                      TnP_dr =  reco::deltaR(tagEl.eta(),tagEl.phi(),probeTrack->eta(),probeTrack->phi()) ;
                      TnP_mass = m;
                      sv_mass = sv_massCalc[p.first];
                      sv_chi2 = sv.chi2();
                      sv_ndof = sv.ndof();
                      sv_prob = TMath::Prob(sv_chi2,sv_ndof);
                      sv_charge = sv_massCalc[p.first];
                      // sv_pt = sv.pt();
                      // sv_eta = sv.eta();
                      // sv_phi = sv.phi();
                      sv_nTrack = sv.tracksSize();
                      sv_x = sv.x();
                      sv_y = sv.y();
                      sv_z = sv.z();
                      Tag_pt = tagEl.pt();
                      Tag_eta = tagEl.eta();
                      Tag_phi = tagEl.phi();
                      if(gedTag){
                          const auto iter = electrons->ptrAt(tag_indices.first);
                          const auto valV = (*cutvValues)[iter];
                          const auto valL = (*cutlValues)[iter];
                          const auto valM = (*cutmValues)[iter];
                          const auto valT = (*cuttValues)[iter];
                          // cout << "tag is ged with val " << val << endl;
                          Tag_cutBased = valV+valL+valM+valT; //(*cutvValues)[electrons->ptrAt(tag_indices.first)];
                          Tag_mvaFall17V2noIso = (*mvaV2NoIsoValues)[electrons->ptrAt(tag_indices.first)];
                          Tag_mvaFall17V2Iso = (*mvaV2IsoValues)[electrons->ptrAt(tag_indices.first)];
                      } else {
                          const auto iter = lowPtElectrons->ptrAt(tag_indices.first);
                          const auto val = (*lowptValues)[iter];
                          Tag_lowPtID = val; //(*lowptValues)[lowPtElectrons->ptrAt(tag_indices.first)];
                      }
                      Tag_isLowPt = gedTag?0:1;
                      Tag_charge = tagEl.charge();
                      Tag_x = tagEl.vx();
                      Tag_y = tagEl.vy();
                      Tag_z = tagEl.vz();

                      Probe_pt = probeTrack->pt();
                      Probe_eta = probeTrack->eta();
                      Probe_phi = probeTrack->phi();
                      if(iGed>=0){
                          const auto iter = electrons->ptrAt(iGed);
                          const auto valV = (*cutvValues)[iter];
                          const auto valL = (*cutlValues)[iter];
                          const auto valM = (*cutmValues)[iter];
                          const auto valT = (*cuttValues)[iter];
                          // cout << "tag is ged with val " << val << endl;
                          Probe_cutBased = valV+valL+valM+valT; //
                          //Probe_cutBased = (*cutvValues)[electrons->ptrAt(iGed)];
                          Probe_mvaFall17V2noIso = (*mvaV2NoIsoValues)[electrons->ptrAt(iGed)];
                          Probe_mvaFall17V2Iso = (*mvaV2IsoValues)[electrons->ptrAt(iGed)];
                      }
                      if(iLow>=0){
                          const auto iter = lowPtElectrons->ptrAt(iLow);
                          const auto val = (*lowptValues)[iter];
                          Probe_lowPtID = val; //(*lowptValues)[lowPtElectrons->ptrAt(iLow)];
                      }
                      Probe_charge = probeTrack->charge();
                      Probe_matches_GED = iGed>=0;
                      Probe_matches_LowPt = iLow>=0;
                      Probe_x = probeTrack->vx();
                      Probe_y = probeTrack->vy();
                      Probe_z = probeTrack->vz();

                      pv_x = secondaryVertexHandle->begin()->x();
                      pv_y = secondaryVertexHandle->begin()->y();
                      pv_z = secondaryVertexHandle->begin()->z();
                      
                      // kTrack...
                      int ik=-1;
                      if(int(p.second[0]==iProbeTrack)+int(p.second[1]==iProbeTrack)+int(p.second[1]==iTagTrack)+int(p.second[0]==iTagTrack)==2) ik=p.second[2];
                      if(int(p.second[0]==iProbeTrack)+int(p.second[2]==iProbeTrack)+int(p.second[2]==iTagTrack)+int(p.second[0]==iTagTrack)==2) ik=p.second[1];
                      if(int(p.second[2]==iProbeTrack)+int(p.second[1]==iProbeTrack)+int(p.second[1]==iTagTrack)+int(p.second[2]==iTagTrack)==2) ik=p.second[0];
                      // cout <<"A"<<int(p.second[0]==iProbeTrack)+int(p.second[1]==iProbeTrack)+int(p.second[1]==iTagTrack)+int(p.second[0]==iTagTrack)<<" "<<p.second[2]<<endl;
                      // cout <<"B"<<int(p.second[0]==iProbeTrack)+int(p.second[2]==iProbeTrack)+int(p.second[2]==iTagTrack)+int(p.second[0]==iTagTrack)<<" "<<p.second[1]<<endl;
                      // cout <<"C"<<int(p.second[2]==iProbeTrack)+int(p.second[1]==iProbeTrack)+int(p.second[1]==iTagTrack)+int(p.second[2]==iTagTrack)<<" "<<p.second[0]<<endl;
                      // cout <<"ik " << ik<<endl;
                      if(ik>=0){
                          const auto& kTrack = sv.trackRefAt(ik);
                          K_pt = kTrack->pt();
                          K_eta = kTrack->eta();
                          K_phi = kTrack->phi();
                          K_x = kTrack->vx();
                          K_y = kTrack->vy();
                          K_z = kTrack->vz();
                          dr_K_Probe =  reco::deltaR(probeTrack->eta(), probeTrack->phi(), kTrack->eta(), kTrack->phi());
                          dr_K_Tag   =  reco::deltaR(tagEl.eta(), tagEl.phi(), kTrack->eta(), kTrack->phi());
                          //cout <<Tag_pt << " " << Probe_pt << " " << K_pt <<endl;
                      }
                      
                      // calc ht if we haven't already done so
                      if(ht25<-0.5){
                          ht25=0;
                          for(const auto& j : *jetHandle){
                              if(j.pt()<25) continue;
                              ht25 += j.pt();
                          }
                      }
                      
                      TnP_ht = ht25;
                      TnP_met = metHandle->at(0).pt();
                      TnP_ipair = iPair;

                      reco_tree->Fill();
                      iPair++;
                  }
              }
          }
      }
  }
  return;


  // for(auto sv = secondaryVertexHandle->begin(); sv!= secondaryVertexHandle->end(); ++sv)
  // {
  //     ChargeHelper=0.;
  //     P_sum.SetPtEtaPhiM(0.,0.,0.,0.);
  //     int SVtrackSize = 0;
  //     for(long unsigned int r=0;r<sv->tracksSize();r++) //calcaulte variables needed to select KJPsi events
  //     {
  //       if(sv->trackRefAt(r)->pt()>0.7)
  //       {
  //         ChargeHelper+=sv->trackRefAt(r)->charge();
  //         P_tmp.SetPtEtaPhiM(sv->trackRefAt(r)->pt(),sv->trackRefAt(r)->eta(),sv->trackRefAt(r)->phi(),0.);
  //         P_sum+=P_tmp;
  //         SVtrackSize++;
  //       }
  //     }
  //     if(P_sum.M()>4.5 && SVtrackSize==3 && std::abs(ChargeHelper)==1) //fill only if satisfied
  //     {
  //       for(long unsigned int r=0;r<sv->tracksSize();r++)
  //       {
  //         if(sv->trackRefAt(r)->pt()>0.7)
  //         {
  //           SVpt.push_back(sv->trackRefAt(r)->pt());
  //           SVeta.push_back(sv->trackRefAt(r)->eta());
  //           SVphi.push_back(sv->trackRefAt(r)->phi());
  //           SVcharge.push_back(sv->trackRefAt(r)->charge());
  //         }
  //       }
  //       SVchi2.push_back(sv->chi2());
  //       SVndof.push_back(sv->ndof());
  //       SVx.push_back(sv->x());
  //       SVy.push_back(sv->y());
  //       SVz.push_back(sv->z());
  //       SVmass.push_back(P_sum.M());
  //       numGoodSV++;
  //     }
  // }


  // for(auto el = electrons->begin(); el != electrons->end(); ++el)
  // {

  //   if(el->closestCtfTrackRef().isNonnull() && el->pt()>5.)  //what if the reference to track does not exist ?
  //   {
  //     auto ClosestTrackRef = tkColl->begin() + el->closestCtfTrackRef().index();
  //     int FoundSVmatch=0;
  //     for(long unsigned int b=0; b<SVpt.size(); b++)
  //     {
  //       // TODO: fix the float comparison here
  //       if((float)ClosestTrackRef->pt()==SVpt.at(b) && (float)ClosestTrackRef->eta()==SVeta.at(b) && (float)ClosestTrackRef->phi()==SVphi.at(b))
  //       {
  //         EleSVtrkref.push_back(b);
  //         FoundSVmatch=1;
  //       }
  //     }
  //     if(FoundSVmatch==0)
  //     EleSVtrkref.push_back(99);
  //   }
  //   else if(el->pt()>5.)
  //   EleSVtrkref.push_back(99);

  //   if(el->pt()>5.)
  //   {
  //     const auto iter = electrons->ptrAt(numele);
  //     ele_pt.push_back(el->pt());
  //     ele_eta.push_back(el->eta());
  //     ele_phi.push_back(el->phi());
  //     ele_charge.push_back(el->charge());
  //     ElectronMVAEstimatorRun2Fall17NoIsoV2Values.push_back((*mvaV2NoIsoValues)[iter]);
  //     ElectronMVAEstimatorRun2Fall17IsoV2Values.push_back((*mvaV2IsoValues)[iter]);
  //     EleTrkIsNoNull.push_back(el->closestCtfTrackRef().isNonnull());
  //     ele_x.push_back(el->vx());
  //     ele_y.push_back(el->vy());
  //     ele_z.push_back(el->vz());
  //     numele++;
  //   }
  // }

  // for(long unsigned int b=0; b<SVpt.size();b++)
  // {
  //   int RefSVtoEl=99;
  //   for(long unsigned c=0; c<EleSVtrkref.size();c++)
  //   {
  //     if(EleSVtrkref.at(c)==(int)b)
  //     RefSVtoEl=c;
  //   }
  //   SVTrEleRef.push_back(RefSVtoEl);
  // }

  // evtnum=iEvent.eventAuxiliary().event();

  // if((numele==1 || numele==2) && numGoodSV>=1) //add charge condition ?
  // {
  //   //std::cout<<iEvent.eventAuxiliary().event()<<std::endl;
  //   reco_tree->Fill();
  // }

}


void SosSkimAnalyzer::beginJob() {

}

void SosSkimAnalyzer::endJob() {

}

void SosSkimAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);

}

DEFINE_FWK_MODULE(SosSkimAnalyzer);
