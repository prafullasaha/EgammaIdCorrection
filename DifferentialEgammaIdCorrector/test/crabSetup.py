from CRABClient.UserUtilities import config
config = config()

config.General.workArea = 'crab_projects'
config.General.transferOutputs = True
config.General.transferLogs = True

config.JobType.pluginName  = 'Analysis'
#config.JobType.outputFiles = ['noise_histos.root']


###for data
config.General.requestName = 'DYJetsToLL_M-50'
#config.JobType.psetName = 'ecalvalidationDATAAOD_cfg.py'
config.JobType.psetName = 'DifferentialEgammaIdCorrector_cfg.py'

#config.Data.useParent = True
#config.JobType.maxMemoryMB = 4000

config.Data.inputDataset = '/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3_ext2-v1/MINIAODSIM'



#config.Data.inputDataset = '/ZeroBias/Run2017B-06Jul2017-v2/AOD'
#config.Data.inputDataset = '/SingleNeutrino/PhaseISpring17DR-FlatPU28to62_90X_upgrade2017_realistic_v20-v1/AODSIM'
#config.Data.secondaryInputDataset = '/SingleNeutrino/PhaseISpring17DR-FlatPU28to62_90X_upgrade2017_realistic_v20-v1/GEN_SIM-RECO'
#config.Data.inputDataset = '/MinBias_TuneCUETP8M1_13TeV-pythia8/PhaseIFall16DR-NoPUNZS_90X_upgrade2017_realistic_v6_C1_ext1-v1/AODSIM'

#
#
#config.Data.lumiMask = '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions16/13TeV/ReReco/Final/Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON.txt'
#config.Data.runRange = '283946'
###config.Data.runRange = '275601-275603'


config.Data.inputDBS = 'global'
#config.JobType.inputFiles = ['Summer16_23Sep2016BCDV4_DATA_L2L3Residual_AK8PFchs.txt', 'Summer16_23Sep2016BCDV4_DATA_L2Relative_AK8PFchs.txt']

config.Data.outputDatasetTag = 'Oct08'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 10
config.Data.outLFNDirBase = '/store/user/prsaha/'
config.Data.publication = False

config.Site.storageSite ='T2_IN_TIFR' 
