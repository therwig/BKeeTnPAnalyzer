#include <memory>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"

#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "FWCore/Common/interface/TriggerResultsByName.h"

using namespace std;

class TriggerAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources> {
public:
  explicit TriggerAnalyzer(const edm::ParameterSet&);
  ~TriggerAnalyzer();

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  void beginJob() override;
  void analyze(const edm::Event&, const edm::EventSetup&) override;
  void endJob() override;

  edm::EDGetTokenT<edm::TriggerResults> triggerBits_;


};

TriggerAnalyzer::TriggerAnalyzer(const edm::ParameterSet& iConfig):
triggerBits_(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("bits")))
{

}

TriggerAnalyzer::~TriggerAnalyzer() {

}


void TriggerAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  using namespace edm;

  edm::Handle<edm::TriggerResults> triggerBits;
  iEvent.getByToken(triggerBits_, triggerBits);

  const edm::TriggerNames &names = iEvent.triggerNames(*triggerBits);
  std::cout << "\n == TRIGGER PATHS= " << std::endl;
  for (unsigned int i = 0, n = triggerBits->size(); i < n; ++i) {
      /*std::cout << "Trigger " << names.triggerName(i) <<
              ", prescale " << triggerPrescales->getPrescaleForIndex(i) <<
              ": " << (triggerBits->accept(i) ? "PASS" : "fail (or not run)")
              << std::endl;*/
      if(triggerBits->accept(i))
      std::cout << "Trigger " << names.triggerName(i)<<std::endl;
  }


}

// ------------ method called once each job just before starting event loop  ------------
void TriggerAnalyzer::beginJob() {
  // please remove this method if not needed
}

// ------------ method called once each job just after ending the event loop  ------------
void TriggerAnalyzer::endJob() {
  // please remove this method if not needed
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void TriggerAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);

  //Specify that only 'tracks' is allowed
  //To use, remove the default given above and uncomment below
  //ParameterSetDescription desc;
  //desc.addUntracked<edm::InputTag>("tracks","ctfWithMaterialTracks");
  //descriptions.addWithDefaultLabel(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(TriggerAnalyzer);
