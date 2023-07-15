import FWCore.ParameterSet.Config as cms
from Configuration.StandardSequences.Reconstruction_cff import*
from Configuration.AlCa.GlobalTag import GlobalTag

process = cms.Process("Reco")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("Configuration.StandardSequences.GeometryRecoDB_cff")
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

#process.GlobalTag.globaltag = '106X_dataRun2_v35'

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1))
process.MessageLogger.cerr.FwkReport.reportEvery = 10000000

process.source = cms.Source("PoolSource",fileNames = cms.untracked.vstring('file:../TestFiles/CD7C6E2F-3535-B347-A5F9-060D8A63508A.root'))

process.reco = cms.EDAnalyzer('RecoSkimAnalyzer',
   Electron = cms.untracked.InputTag("gedGsfElectrons"),
   Track = cms.untracked.InputTag("generalTracks"),
   mvaV2IsoValuesMap = cms.untracked.InputTag("electronMVAValueMapProducer:ElectronMVAEstimatorRun2Fall17IsoV2Values"),
   mvaV2NoIsoValuesMap = cms.untracked.InputTag("electronMVAValueMapProducer:ElectronMVAEstimatorRun2Fall17NoIsoV2Values"),
   secondaryVertices = cms.untracked.InputTag("inclusiveSecondaryVertices", "", "RECO"),
                              )

from EgammaUser.EgammaPostRecoTools.EgammaPostRecoTools import setupEgammaPostRecoSeq,makeEgammaPATWithUserData

setupEgammaPostRecoSeq(process,
                       era='2018-UL',
                       isMiniAOD=False,
                       runVID=True
)

process.TFileService = cms.Service("TFileService", fileName=cms.string("Reco_tree.root"))

process.p = cms.Path(process.egammaPostRecoSeq*process.reco)
