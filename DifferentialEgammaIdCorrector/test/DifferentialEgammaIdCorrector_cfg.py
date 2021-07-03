import FWCore.ParameterSet.Config as cms
import os
#import EgammaIdCorrection.DifferentialEgammaIdCorrector.DifferentialEgammaIdCorrector_cfi as EgammaIdInp_Corr

InFile = os.environ.get('InFile')
OutFile = os.environ.get('OutFile')

print ("Hi")
print(InFile)
print(OutFile)
print ("bye")
process = cms.Process("EgammaIdCorrection")

process.load("FWCore.MessageLogger.MessageLogger_cfi")

process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32(1000)
process.MessageLogger.cerr = cms.untracked.PSet(threshold = cms.untracked.string("ERROR"))

process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff")
from Configuration.AlCa.GlobalTag_condDBv2 import GlobalTag
#process.GlobalTag = GlobalTag(process.GlobalTag, '94X_mcRun2_asymptotic_v3')
process.GlobalTag = GlobalTag(process.GlobalTag, '94X_dataRun2_v10')
process.load("Configuration.StandardSequences.MagneticField_cff")

process.source = cms.Source("PoolSource",
#    	skipEvents = cms.untracked.uint32(0),
	fileNames = cms.untracked.vstring(InFile))
#"/store/mc/RunIISpring18DRPremix/SingleNeutrino_Spring18/AODSIM/100X_upgrade2018_realistic_v10-v2/40000/86E12AB8-9990-E811-BFEA-FA163E53E163.root"
#MC_16
#"/store/mc/RunIISummer16MiniAODv3/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/PUMoriond17_94X_mcRun2_asymptotic_v3_ext2-v1/270000/FEFC23FC-37C7-E811-97DA-0CC47AA53D86.root"
#"/store/mc/RunIISummer16MiniAODv3/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/PUMoriond17_94X_mcRun2_asymptotic_v3_ext2-v1/100000/8E7A4028-2DBD-E811-B0B1-008CFA1660F8.root",
#"/store/mc/RunIISummer16MiniAODv3/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/PUMoriond17_94X_mcRun2_asymptotic_v3_ext2-v1/270000/FEC9FD38-0EC7-E811-8C4B-00259029E81A.root",
#"/store/mc/RunIISummer16MiniAODv3/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/PUMoriond17_94X_mcRun2_asymptotic_v3_ext2-v1/270000/FE984109-17C7-E811-8D1F-0025901AC0FA.root",
#"/store/mc/RunIISummer16MiniAODv3/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/PUMoriond17_94X_mcRun2_asymptotic_v3_ext2-v1/270000/FE7E48CF-F1C0-E811-9E19-008CFA1660F8.root",
#"/store/mc/RunIISummer16MiniAODv3/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/PUMoriond17_94X_mcRun2_asymptotic_v3_ext2-v1/270000/FCAD4B77-F6C0-E811-99BB-B496910A0554.root",
#"/store/mc/RunIISummer16MiniAODv3/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/PUMoriond17_94X_mcRun2_asymptotic_v3_ext2-v1/270000/FC6CAA61-88BD-E811-AB57-0CC47A57CB62.root",
#"/store/mc/RunIISummer16MiniAODv3/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/PUMoriond17_94X_mcRun2_asymptotic_v3_ext2-v1/270000/FC5A7451-59C7-E811-81C2-002590FD5A48.root",
#"/store/mc/RunIISummer16MiniAODv3/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/PUMoriond17_94X_mcRun2_asymptotic_v3_ext2-v1/270000/FC03A7B5-7FC0-E811-99BA-B496910A9A2C.root",
#"/store/mc/RunIISummer16MiniAODv3/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/PUMoriond17_94X_mcRun2_asymptotic_v3_ext2-v1/270000/FADA5029-2CC7-E811-BECD-0025907D1D6C.root",
#"/store/mc/RunIISummer16MiniAODv3/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/PUMoriond17_94X_mcRun2_asymptotic_v3_ext2-v1/270000/F855E2C9-30C7-E811-A621-0025901AC0F8.root",
#"/store/mc/RunIISummer16MiniAODv3/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/PUMoriond17_94X_mcRun2_asymptotic_v3_ext2-v1/270000/F6A1640F-65C7-E811-BA45-0025901F8740.root",
#"/store/mc/RunIISummer16MiniAODv3/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/PUMoriond17_94X_mcRun2_asymptotic_v3_ext2-v1/270000/F6989592-7DC7-E811-96B7-0CC47A0AD3BC.root"

#MC_17
#"/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14_ext1-v1/90000/FEFADA2F-8C44-E811-914F-B496910A8618.root"

#Data_16
#"/store/data/Run2016B/DoubleEG/MINIAOD/17Jul2018_ver2-v1/50000/F2DAD002-F38B-E811-82DA-008CFA197AEC.root",
#"/store/data/Run2016B/DoubleEG/MINIAOD/17Jul2018_ver2-v1/50000/F019830A-0B8C-E811-9A6E-008CFA1979B0.root"
#Data_17
#"/store/data/Run2017B/DoubleEG/MINIAOD/31Mar2018-v1/80000/FAC7DC8A-3737-E811-8BA7-6CC2173DC380.root",
#"/store/data/Run2017B/DoubleEG/MINIAOD/31Mar2018-v1/80000/F89A88BC-3837-E811-8017-0017A4771048.root"
#	)
#)

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(1000)
)

process.load("EgammaIdCorrection.DifferentialEgammaIdCorrector.DifferentialEgammaIdCorrector_cfi")
import EgammaIdCorrection.DifferentialEgammaIdCorrector.DifferentialEgammaIdCorrector_cfi as EgammaIdInp_Corr
EgammaIdInp_Corr.setup_DifferentialPhoIdInputsCorrection( process )
process.TFileService = cms.Service("TFileService",
#    fileName = cms.string('output.root')
    fileName = cms.string(OutFile)
)

process.p = cms.Path(
     process.differentialEgammaIdCorrection
    )
