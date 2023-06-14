import CRABClient
from CRABClient.UserUtilities import config

config = config()

#config.General.requestName = 'Bpark_MiniAOD_MC_reliso_UL' # Naming for CRAB and maps
config.General.requestName = 'ParkingBPH1B_AOD_DATA_UL_v5'
config.General.workArea = 'crab_projects' # area where information about crab jobs will be saved. You should create this directory next to this file
config.General.transferOutputs = True

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'ConfFile_cfg.py' # CMS config code we wish to run on Crab jobs


config.Data.inputDataset = '/ParkingBPH1/Run2018B-20Jun2021_UL2018-v1/AOD'
#config.Data.inputDataset = '/BuToKJpsi_Toee_Mufilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/RunIISummer19UL18MiniAODv2-TrkExtra_106X_upgrade2018_realistic_v16_L1v1-v1/MINIAODSIM' # MC
config.Data.inputDBS = 'global'
config.Data.splitting = 'FileBased' # IF we wish to run by files per job
#config.Data.splitting = 'Automatic' # If we wish to optimize files per job (sometimes does not work)
#config.Data.splitting = 'LumiBased'
config.Data.unitsPerJob = 5 #running 20 files per job
#config.Data.lumiMask = 'https://cms-service-dqmdc.web.cern.ch/CAF/certification/Collisions18/13TeV/PromptReco/Cert_314472-325175_13TeV_PromptReco_Collisions18_JSON.txt'
#config.Data.runRange = '275776-275782'
config.Data.publication = True
#config.Data.outputDatasetTag = 'Bpark_MiniAOD_MC_reliso_UL'
config.Data.outputDatasetTag = 'ParkingBPH1B_AOD_DATA_UL_v5' # Naming for CRAB and maps

config.Site.storageSite = 'T3_CH_CERNBOX' # writes the ouput inside your cernbox directoty under the name of 'inputDataset'
config.Data.partialDataset = True
# When this is done, submit CRAB jobs by using command "crab submit -c CrabConfig_mini.py"
