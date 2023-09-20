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
#include "FWCore/Framework/interface/one/EDFilter.h"

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



class SosSkimOnly : public edm::one::EDFilter<edm::one::SharedResources> {
public:
  explicit SosSkimOnly(const edm::ParameterSet&);
  ~SosSkimOnly();

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  void beginJob() override;
  bool filter(edm::Event &, edm::EventSetup const &) override;
  void endJob() override;

  edm::EDGetToken electronsToken_;
  edm::EDGetToken lowPtElectronsToken_;
  edm::EDGetTokenT<reco::TrackCollection> TrackToken;
  edm::EDGetTokenT<edm::ValueMap<float> > cutvValuesMapToken;
  edm::EDGetTokenT<edm::ValueMap<float> > lowptValuesMapToken;
  edm::EDGetTokenT<edm::ValueMap<float> > mvaV2IsoValuesMapToken;
  edm::EDGetTokenT<edm::ValueMap<float> > mvaV2NoIsoValuesMapToken;
  edm::EDGetTokenT<reco::VertexCollection> SecondaryVertexToken;
  edm::EDGetTokenT<std::vector<reco::PFMET> > MetToken;
  edm::EDGetTokenT<std::vector<reco::PFJet> > JetToken;

  edm::InputTag Electron_;
  edm::InputTag LowPtElectron_;
  edm::InputTag TrackInputTag_;
  edm::InputTag cutvValuesMapInputTag_;
  edm::InputTag lowptValuesMapInputTag_;
  edm::InputTag mvaV2IsoValuesMapInputTag_;
  edm::InputTag mvaV2NoIsoValuesMapInputTag_;
  edm::InputTag SecondaryVertexInputTag_;
  edm::InputTag MetInputTag_;
  edm::InputTag JetInputTag_;

    
};

SosSkimOnly::SosSkimOnly(const edm::ParameterSet& iConfig):
Electron_(iConfig.getUntrackedParameter<edm::InputTag>("Electron")),
LowPtElectron_(iConfig.getUntrackedParameter<edm::InputTag>("LowPtElectron")),
TrackInputTag_(iConfig.getUntrackedParameter<edm::InputTag>("Track")),
cutvValuesMapInputTag_(iConfig.getUntrackedParameter<edm::InputTag>("cutV")),
lowptValuesMapInputTag_(iConfig.getUntrackedParameter<edm::InputTag>("lowptValuesMap")),
mvaV2IsoValuesMapInputTag_(iConfig.getUntrackedParameter<edm::InputTag>("mvaV2IsoValuesMap")),
mvaV2NoIsoValuesMapInputTag_(iConfig.getUntrackedParameter<edm::InputTag>("mvaV2NoIsoValuesMap")),
SecondaryVertexInputTag_(iConfig.getUntrackedParameter<edm::InputTag>("secondaryVertices")),
MetInputTag_(iConfig.getUntrackedParameter<edm::InputTag>("met")),
JetInputTag_(iConfig.getUntrackedParameter<edm::InputTag>("jets"))
{
  electronsToken_    = consumes<edm::View<reco::GsfElectron> >(Electron_);
  lowPtElectronsToken_    = consumes<edm::View<reco::GsfElectron> >(LowPtElectron_);
  TrackToken = consumes<reco::TrackCollection>(TrackInputTag_);
  cutvValuesMapToken = consumes<edm::ValueMap<float>>(cutvValuesMapInputTag_);
  lowptValuesMapToken = consumes<edm::ValueMap<float>>(lowptValuesMapInputTag_);
  mvaV2IsoValuesMapToken = consumes<edm::ValueMap<float>>(mvaV2IsoValuesMapInputTag_);
  mvaV2NoIsoValuesMapToken = consumes<edm::ValueMap<float>>(mvaV2NoIsoValuesMapInputTag_);
  SecondaryVertexToken = consumes<reco::VertexCollection>(SecondaryVertexInputTag_);
  MetToken = consumes<std::vector<reco::PFMET> >(MetInputTag_);
  JetToken = consumes<std::vector<reco::PFJet> >(JetInputTag_);


}

SosSkimOnly::~SosSkimOnly() {

}

bool SosSkimOnly::filter(edm::Event& iEvent, const edm::EventSetup& iSetup) {
  using namespace edm;

  edm::Handle<edm::View<reco::GsfElectron> > electrons;
  iEvent.getByToken(electronsToken_, electrons);

  edm::Handle<edm::View<reco::GsfElectron> > lowPtElectrons;
  iEvent.getByToken(lowPtElectronsToken_, lowPtElectrons);

  edm::Handle<reco::TrackCollection> TrackHandle;
  iEvent.getByToken(TrackToken, TrackHandle);

  // edm::Handle<edm::ValueMap<float> > cutvValues;
  // iEvent.getByToken(cutvValuesMapToken,cutvValues);
  // edm::Handle<edm::ValueMap<float> > lowptValues;
  // iEvent.getByToken(lowptValuesMapToken,lowptValues);

  // edm::Handle<edm::ValueMap<float> > mvaV2NoIsoValues;
  // iEvent.getByToken(mvaV2NoIsoValuesMapToken,mvaV2NoIsoValues);

  // edm::Handle<edm::ValueMap<float> > mvaV2IsoValues;
  // iEvent.getByToken(mvaV2IsoValuesMapToken,mvaV2IsoValues);

  edm::Handle<reco::VertexCollection> secondaryVertexHandle;
  iEvent.getByToken(SecondaryVertexToken,secondaryVertexHandle);

  // edm::Handle<std::vector<reco::PFMET> > metHandle;
  // iEvent.getByToken(MetToken,metHandle);

  // edm::Handle<std::vector<reco::PFJet> > jetHandle;
  // iEvent.getByToken(JetToken,jetHandle);

  // must have at least one electron (GED or lowpt) and an SV
  if( secondaryVertexHandle->size()==0 || (electrons->size()==0 && lowPtElectrons->size()==0) ) return 0;

  //helper variables
  TLorentzVector P_sum,P_tmp;

  const reco::TrackCollection* tkColl = TrackHandle.product();
  const reco::VertexCollection* svColl = secondaryVertexHandle.product();

  // std::vector<unsigned int> good_sv_idx{};
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
      for(long unsigned int r=0;r<sv.tracksSize();r++){
          const auto & tk = sv.trackRefAt(r);
          if(tk->pt()>0.7){
              svCharge+=tk->charge();
              P_tmp.SetPtEtaPhiM(tk->pt(),tk->eta(),tk->phi(),0.);
              P_sum+=P_tmp;
              track_indices.push_back(r);
          }
      }
      // if(P_sum.M()>4.5 && track_indices.size()==3 && std::abs(svCharge)==1){ //fill only if satisfied
      if(track_indices.size()==3 && std::abs(svCharge)==1){ //fill only if satisfied
          // good_sv_idx.push_back(isv);
          sv_tracks[isv] = track_indices;
          sv_massCalc[isv] = P_sum.M();
          sv_chargeCalc[isv] = svCharge;
      }
  }
  if(sv_tracks.size()==0) return 0;

  return 1;

  // int iPair=0;
  // float ht25=-1;
  // for(const auto& p : sv_tracks){
  //     const auto& sv = svColl->at(p.first);
  //     // check if any SV tracks match to a TAG electron
  //     std::vector<std::pair<unsigned int, unsigned int>> ged_tk_tag_indices{};
  //     std::vector<std::pair<unsigned int, unsigned int>> lowpt_tk_tag_indices{};
  //     for(const unsigned int track_index : p.second){
  //         const auto& tk = sv.trackRefAt(track_index);
  //         // check for a GED TAG
  //         for(unsigned int i=0; i < electrons->size(); i++){
  //             const auto& el = electrons->at(i);
  //             if(!el.closestCtfTrackRef().isNonnull()) continue;
  //             const auto &elTkRef = tkColl->begin() + el.closestCtfTrackRef().index();
  //             if(tracksMatch(*elTkRef, *tk)) {
  //                 ged_tk_tag_indices.push_back(std::make_pair(i,track_index));
  //             }
  //         }
  //         // check for a LowPt TAG
  //         for(unsigned int i=0; i < lowPtElectrons->size(); i++){
  //             const auto& el = lowPtElectrons->at(i);
  //             if(!el.closestCtfTrackRef().isNonnull()) continue;
  //             const auto &elTkRef = tkColl->begin() + el.closestCtfTrackRef().index();
  //             if(tracksMatch(*elTkRef, *tk)) {
  //                 lowpt_tk_tag_indices.push_back(std::make_pair(i,track_index));
  //             }
  //         }
  //     }
  //     // Look for j/psi(ee) candidates among the other tracks
  //     for(bool gedTag : {0,1}){
  //         const auto& el_tk_tag_indices = (gedTag ? ged_tk_tag_indices : lowpt_tk_tag_indices);
  //         for(const auto tag_indices : el_tk_tag_indices){
  //             const auto &eleColl = (gedTag ? electrons : lowPtElectrons);
  //             const auto& tagEl = eleColl->at(tag_indices.first);
  //             const unsigned int iTagTrack = tag_indices.second;
  //             for(const unsigned int iProbeTrack : p.second){
  //                 if (iTagTrack==iProbeTrack) continue;
  //                 const auto& probeTrack = sv.trackRefAt(iProbeTrack);
  //                 //  check for a j/psi candidate
  //                 float m = getMass(tagEl, *probeTrack);
  //                 if (m > 2.4 && m < 3.8){
  //                     // probe track matches GED
  //                     int iGed=-1;
  //                     float min_dr=9e9;                      
  //                     for(unsigned int i=0; i < electrons->size(); i++){
  //                         const auto& el = electrons->at(i);
  //                         float dr = reco::deltaR(el.eta(),el.phi(),probeTrack->eta(),probeTrack->phi());
  //                         if (dr < 0.01 && dr < min_dr){
  //                             min_dr = dr;
  //                             iGed=i;
  //                         }
  //                     }
  //                     // probe track matches LowPt
  //                     int iLow=-1;
  //                     min_dr=9e9;
  //                     for(unsigned int i=0; i < electrons->size(); i++){
  //                         const auto& el = lowPtElectrons->at(i);
  //                         float dr = reco::deltaR(el.eta(),el.phi(),probeTrack->eta(),probeTrack->phi());
  //                         if (dr < 0.01 && dr < min_dr){
  //                             min_dr = dr;
  //                             iLow=i;
  //                         }
  //                     }
                      

  //                     // fill tree
  //                     TnP_dr =  reco::deltaR(tagEl.eta(),tagEl.phi(),probeTrack->eta(),probeTrack->phi()) ;
  //                     TnP_mass = m;
  //                     sv_mass = sv_massCalc[p.first];
  //                     sv_chi2 = sv.chi2();
  //                     sv_ndof = sv.ndof();
  //                     sv_charge = sv_massCalc[p.first];
  //                     // sv_pt = sv.pt();
  //                     // sv_eta = sv.eta();
  //                     // sv_phi = sv.phi();
  //                     sv_nTrack = sv.tracksSize();
  //                     sv_x = sv.x();
  //                     sv_y = sv.y();
  //                     sv_z = sv.z();
  //                     Tag_pt = tagEl.pt();
  //                     Tag_eta = tagEl.eta();
  //                     Tag_phi = tagEl.phi();
  //                     if(gedTag){
  //                         Tag_cutBased = (*cutvValues)[electrons->ptrAt(tag_indices.first)];
  //                         Tag_mvaFall17V2noIso = (*mvaV2NoIsoValues)[electrons->ptrAt(tag_indices.first)];
  //                         Tag_mvaFall17V2Iso = (*mvaV2IsoValues)[electrons->ptrAt(tag_indices.first)];
  //                     } else {
  //                         Tag_lowPtID = (*lowptValues)[lowPtElectrons->ptrAt(tag_indices.first)];
  //                     }
  //                     Tag_isLowPt = gedTag?0:1;
  //                     Tag_charge = tagEl.charge();
  //                     Probe_pt = probeTrack->pt();
  //                     Probe_eta = probeTrack->eta();
  //                     Probe_phi = probeTrack->phi();
  //                     if(iGed>=0){
  //                         Probe_cutBased = (*cutvValues)[electrons->ptrAt(iGed)];
  //                         Probe_mvaFall17V2noIso = (*mvaV2NoIsoValues)[electrons->ptrAt(iGed)];
  //                         Probe_mvaFall17V2Iso = (*mvaV2IsoValues)[electrons->ptrAt(iGed)];
  //                     }
  //                     if(iLow>=0){
  //                         Probe_lowPtID = (*lowptValues)[lowPtElectrons->ptrAt(iLow)];
  //                     }
  //                     Probe_charge = probeTrack->charge();
  //                     Probe_matches_GED = iGed>=0;
  //                     Probe_matches_LowPt = iLow>=0;
                      
  //                     // calc ht if we haven't already done so
  //                     if(ht25<-0.5){
  //                         ht25=0;
  //                         for(const auto& j : *jetHandle){
  //                             if(j.pt()<25) continue;
  //                             ht25 += j.pt();
  //                         }
  //                     }
                      
  //                     TnP_ht = ht25;
  //                     TnP_met = metHandle->at(0).pt();
  //                     TnP_ipair = iPair;

  //                     reco_tree->Fill();
  //                     iPair++;
  //                 }
  //             }
  //         }
  //     }
  // }
  // return 0;


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


void SosSkimOnly::beginJob() {

}

void SosSkimOnly::endJob() {

}

void SosSkimOnly::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);

}

DEFINE_FWK_MODULE(SosSkimOnly);
