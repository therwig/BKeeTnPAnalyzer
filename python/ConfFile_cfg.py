import FWCore.ParameterSet.Config as cms
from Configuration.StandardSequences.Reconstruction_cff import*
from Configuration.AlCa.GlobalTag import GlobalTag

process = cms.Process("Reco")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("Configuration.StandardSequences.GeometryRecoDB_cff")
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.GlobalTag.globaltag = '106X_dataRun2_v35'

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(10000))
process.MessageLogger.cerr.FwkReport.reportEvery = 100000

process.source = cms.Source("PoolSource",fileNames = cms.untracked.vstring('file:../TestFiles/013C3F9A-271C-A644-B39D-2FD1E9CA9A60.root'))

process.trigger = cms.EDAnalyzer('TriggerAnalyzer',
                              bits = cms.InputTag("TriggerResults","","HLT")
                              )

process.reco = cms.EDAnalyzer('RecoAnalyzerV2',
   Electron = cms.untracked.InputTag("gedGsfElectrons"),
   Track = cms.untracked.InputTag("generalTracks"),
   mvaV2IsoValuesMap = cms.untracked.InputTag("electronMVAValueMapProducer:ElectronMVAEstimatorRun2Fall17IsoV2Values"),
   mvaV2NoIsoValuesMap = cms.untracked.InputTag("electronMVAValueMapProducer:ElectronMVAEstimatorRun2Fall17NoIsoV2Values"),
   GsfTrack = cms.untracked.InputTag("electronGsfTracks"),
   secondaryVertices = cms.untracked.InputTag("inclusiveSecondaryVertices", "", "RECO"),
   primaryVertices = cms.untracked.InputTag("offlinePrimaryVertices")
                              )

from EgammaUser.EgammaPostRecoTools.EgammaPostRecoTools import setupEgammaPostRecoSeq,makeEgammaPATWithUserData

setupEgammaPostRecoSeq(process,
                       era='2018-UL',
                       isMiniAOD=False,
                       runVID=True
)

process.TFileService = cms.Service("TFileService", fileName=cms.string("Reco_tree.root"))

process.p = cms.Path(process.egammaPostRecoSeq*process.reco)
