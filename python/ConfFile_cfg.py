import FWCore.ParameterSet.Config as cms
from Configuration.StandardSequences.Reconstruction_cff import*
from Configuration.AlCa.GlobalTag import GlobalTag

process = cms.Process("Reco")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("Configuration.StandardSequences.GeometryRecoDB_cff")
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.GlobalTag.globaltag = '106X_dataRun2_v35'

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100) )
process.MessageLogger.cerr.FwkReport.reportEvery = 1000000

process.source = cms.Source("PoolSource",fileNames = cms.untracked.vstring('file:863EAD78-D853-F348-B135-00B596FB6478.root'))

process.trigger = cms.EDAnalyzer('TriggerAnalyzer',
                              bits = cms.InputTag("TriggerResults","","HLT")
                              )

process.reco = cms.EDAnalyzer('RecoAnalyzer',
   Electron = cms.untracked.InputTag("gedGsfElectrons"),
   SuperClusterEB = cms.untracked.InputTag("particleFlowSuperClusterECAL","particleFlowSuperClusterECALBarrel"),
   SuperClusterEE = cms.untracked.InputTag("particleFlowSuperClusterECAL","particleFlowSuperClusterECALEndcapWithPreshower"),
   Conversions=cms.untracked.InputTag("allConversions"),
   BeamSpot=cms.untracked.InputTag("offlineBeamSpot"),
   eleIdLoose = cms.untracked.InputTag("eidLoose"),
   eleIdRobustLoose = cms.untracked.InputTag("eidRobustLoose"),
   eleIdRobustTight = cms.untracked.InputTag("eidRobustTight"),
   eleIdTight = cms.untracked.InputTag("eidTight"),
   Rho = cms.untracked.InputTag("fixedGridRhoAll"),
   HBHERecHitTag = cms.untracked.InputTag("reducedHcalRecHits:hbhereco")
                              )

process.TFileService = cms.Service("TFileService", fileName=cms.string("Reco_tree.root"))

process.p = cms.Path(process.reco)
