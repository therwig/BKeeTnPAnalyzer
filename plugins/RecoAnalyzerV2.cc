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
  ele_z,ele_x,ele_y,tr_z,tr_x,tr_y,gsf_z,gsf_x,gsf_y,SVchi2,SVndof,gsf_pt,gsf_eta,gsf_phi,
  SVpt,SVeta,SVphi,SVcharge,SVmass,SVdistanceToMaxPtPV,SVx,SVy,SVz;
  std::vector<int> ele_charge,EleTRKref,EleTrkIsNoNull,ele_gsf_charge,gsf_charge,tr_charge,
  sv_tr_size,EleSVtrkref,SVTrEleRef,SVsize,SVchargeSum;
  int numele, PFnumele,EleCounter,match,numtr,numgsf,numsv;
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
  reco_tree->Branch("SVchi2",&SVchi2);
  reco_tree->Branch("SVndof",&SVndof);
  reco_tree->Branch("numsv",&numsv);
  reco_tree->Branch("SVpt",&SVpt);
  reco_tree->Branch("SVphi",&SVphi);
  reco_tree->Branch("SVeta",&SVeta);
  reco_tree->Branch("EleSVtrkref",&EleSVtrkref);
  reco_tree->Branch("SVTrEleRef",&SVTrEleRef);
  reco_tree->Branch("SVcharge",&SVcharge);
  reco_tree->Branch("SVmass",&SVmass);
  reco_tree->Branch("numsv",&numsv);
  reco_tree->Branch("SVdistanceToMaxPtPV",&SVdistanceToMaxPtPV);
  reco_tree->Branch("SVx",&SVx);
  reco_tree->Branch("SVy",&SVy);
  reco_tree->Branch("SVz",&SVz);
  reco_tree->Branch("SVsize",&SVsize);
  reco_tree->Branch("SVchargeSum",&SVchargeSum);

  //PrimaryVertex
  reco_tree->Branch("PVmaxPtX",&PVmaxPtX);
  reco_tree->Branch("PVmaxPtY",&PVmaxPtY);
  reco_tree->Branch("PVmaxPtZ",&PVmaxPtZ);

}

RecoAnalyzerV2::~RecoAnalyzerV2() {

}


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

  SVchi2.clear();
  SVndof.clear();
  SVpt.clear();
  SVphi.clear();
  SVeta.clear();
  EleSVtrkref.clear();
  SVTrEleRef.clear();
  SVcharge.clear();
  SVmass.clear();
  numsv=0;
  numsv=0;
  SVdistanceToMaxPtPV.clear();
  SVx.clear();
  SVy.clear();
  SVz.clear();
  SVsize.clear();
  SVchargeSum.clear();

  PVmaxPtX=0;
  PVmaxPtY=0;
  PVmaxPtZ=0;

  //helper variables
  float SumPt,TrackPtHelper=0.;
  int index,num5,TrackIndexHelper,ChargeHelper;
  TLorentzVector P_sum,P_tmp;
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
      ChargeHelper=0.;
      P_sum.SetPtEtaPhiM(0.,0.,0.,0.);
      int SVtrackSize = 0;
      for(long unsigned int r=0;r<sv->tracksSize();r++)
      {
        if(sv->trackRefAt(r)->pt()>0.7)
        {
          SVpt.push_back(sv->trackRefAt(r)->pt());
          SVeta.push_back(sv->trackRefAt(r)->eta());
          SVphi.push_back(sv->trackRefAt(r)->phi());
          SVcharge.push_back(sv->trackRefAt(r)->charge());
          ChargeHelper+=sv->trackRefAt(r)->charge();
          P_tmp.SetPtEtaPhiM(sv->trackRefAt(r)->pt(),sv->trackRefAt(r)->eta(),sv->trackRefAt(r)->phi(),0.);
          P_sum+=P_tmp;
          SVtrackSize++;
        }
      }
      SVchi2.push_back(sv->chi2());
      SVndof.push_back(sv->ndof());
      SVx.push_back(sv->x());
      SVy.push_back(sv->y());
      SVz.push_back(sv->z());
      SVmass.push_back(P_sum.M());
      SVchargeSum.push_back(ChargeHelper);
      SVsize.push_back(SVtrackSize);
      SVdistanceToMaxPtPV.push_back(sqrt((sv->x()-PVmaxPtX)*(sv->x()-PVmaxPtX)+(sv->y()-PVmaxPtY)*(sv->y()-PVmaxPtY)+
                                    (sv->z()-PVmaxPtZ)*(sv->z()-PVmaxPtZ)));
      numsv++;
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
    numtr++;
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
