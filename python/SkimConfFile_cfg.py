import FWCore.ParameterSet.Config as cms
from Configuration.StandardSequences.Reconstruction_cff import*
from Configuration.AlCa.GlobalTag import GlobalTag

process = cms.Process("Reco2") # skim doesn't like the same process name, can modify w/o any issues

process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("Configuration.StandardSequences.GeometryRecoDB_cff")
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

#process.GlobalTag.globaltag = '106X_dataRun2_v35'

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1))
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(10000))
process.MessageLogger.cerr.FwkReport.reportEvery = 10000000

process.source = cms.Source("PoolSource",fileNames = cms.untracked.vstring('file:../TestFiles/CD7C6E2F-3535-B347-A5F9-060D8A63508A.root'))
process.source = cms.Source("PoolSource",fileNames = cms.untracked.vstring('file:/uscms/home/therwig/nobackup/sos/bkeeTnP/local_files/51EDE899-A4CA-A649-B086-F376B94DEA5F.root')) # /SingleMuon/Run2018D-12Nov2019_UL2018-v8/AOD
process.source = cms.Source("PoolSource",fileNames = cms.untracked.vstring('root://cmseos.fnal.gov//store/user/therwig/tmp/230920_BKJPsi_AOD/job33232548AOD/events.0.root'))
process.source = cms.Source("PoolSource",fileNames = cms.untracked.vstring('file:/uscms/home/therwig/nobackup/sos/bkeeTnP/local_files/230920_BKJPsi_AOD.6p7k.skim.root'))

# one that throws an error
process.source = cms.Source("PoolSource",fileNames = cms.untracked.vstring('root://cms-xrd-global.cern.ch//store/data/Run2018D/SingleMuon/AOD/12Nov2019_UL2018-v8/270000/0ACB5BE8-70BE-4D46-97E4-F7D5E685AEA5.root'))

process.reco = cms.EDAnalyzer('SosSkimAnalyzer',
   Electron = cms.untracked.InputTag("gedGsfElectrons"),
   LowPtElectron = cms.untracked.InputTag("lowPtGsfElectrons"),
   Track = cms.untracked.InputTag("generalTracks"),
   lowptValuesMap = cms.untracked.InputTag("lowPtGsfElectronID"),
   cutV = cms.untracked.InputTag("egmGsfElectronIDs:cutBasedElectronID-Fall17-94X-V2-veto"),
   mvaV2IsoValuesMap = cms.untracked.InputTag("electronMVAValueMapProducer:ElectronMVAEstimatorRun2Fall17IsoV2Values"),
   mvaV2NoIsoValuesMap = cms.untracked.InputTag("electronMVAValueMapProducer:ElectronMVAEstimatorRun2Fall17NoIsoV2Values"),
   # mvaV2NoIsoValuesMap = cms.untracked.InputTag("egmGsfElectronIDs:"),
                              # "egmGsfElectronIDs:cutBasedElectronID-Fall17-94X-V2-veto",
                              # "egmGsfElectronIDs:cutBasedElectronID-Fall17-94X-V2-loose",
                              # "egmGsfElectronIDs:cutBasedElectronID-Fall17-94X-V2-medium",
                              # "egmGsfElectronIDs:cutBasedElectronID-Fall17-94X-V2-tight",
   secondaryVertices = cms.untracked.InputTag("inclusiveSecondaryVertices", "", "RECO"),
   met = cms.untracked.InputTag("pfMet", "", "RECO"),
   jets = cms.untracked.InputTag("ak4PFJets", "", "RECO"),
                              )

from EgammaUser.EgammaPostRecoTools.EgammaPostRecoTools import setupEgammaPostRecoSeq,makeEgammaPATWithUserData

setupEgammaPostRecoSeq(process,
                       era='2018-UL',
                       isMiniAOD=False,
                       runVID=True
)

process.TFileService = cms.Service("TFileService", fileName=cms.string("Reco_tree.root"))

process.p = cms.Path(process.egammaPostRecoSeq*process.reco)


# process.out = cms.OutputModule("PoolOutputModule",
#                                fileName = cms.untracked.string("debug.root"),
#                                SelectEvents = cms.untracked.PSet(SelectEvents = cms.vstring("p"))
#                            )
# process.end = cms.EndPath(process.out)
# process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(50))

## add filter
# process.svFilter = cms.EDFilter("CandViewCountFilter", src = cms.InputTag("inclusiveSecondaryVertices","","RECO"),minNumber = cms.uint32(1))
# process.p.insert(0, cms.Sequence(process.svFilter))

if False:
  process.myfilter = cms.EDFilter('SosSkimOnly',
     Electron = cms.untracked.InputTag("gedGsfElectrons"),
     LowPtElectron = cms.untracked.InputTag("lowPtGsfElectrons"),
     Track = cms.untracked.InputTag("generalTracks"),
     lowptValuesMap = cms.untracked.InputTag("lowPtGsfElectronID"),
     cutV = cms.untracked.InputTag("egmGsfElectronIDs:cutBasedElectronID-Fall17-94X-V2-veto"),
     mvaV2IsoValuesMap = cms.untracked.InputTag("electronMVAValueMapProducer:ElectronMVAEstimatorRun2Fall17IsoV2Values"),
     mvaV2NoIsoValuesMap = cms.untracked.InputTag("electronMVAValueMapProducer:ElectronMVAEstimatorRun2Fall17NoIsoV2Values"),
     secondaryVertices = cms.untracked.InputTag("inclusiveSecondaryVertices", "", "RECO"),
     met = cms.untracked.InputTag("pfMet", "", "RECO"),
     jets = cms.untracked.InputTag("ak4PFJets", "", "RECO"),
                                )
  
  process.p = cms.Path(process.myfilter)
  process.out = cms.OutputModule("PoolOutputModule",
                                 fileName = cms.untracked.string("skim.root"),
                                 SelectEvents = cms.untracked.PSet(SelectEvents = cms.vstring("p"))
                             )
  process.end = cms.EndPath(process.out)
  process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1))
