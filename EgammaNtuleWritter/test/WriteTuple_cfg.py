import FWCore.ParameterSet.Config as cms
import os
#import EgammaIdCorrection.DifferentialEgammaIdCorrector.DifferentialEgammaIdCorrector_cfi as EgammaIdInp_Corr
process = cms.Process("EgammaIdCorrection")

process.load("FWCore.MessageLogger.MessageLogger_cfi")

process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32(10000)
process.MessageLogger.cerr = cms.untracked.PSet(threshold = cms.untracked.string("ERROR"))

process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff")
from Configuration.AlCa.GlobalTag_condDBv2 import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '94X_mcRun2_asymptotic_v3')
process.load("Configuration.StandardSequences.MagneticField_cff")

process.source = cms.Source("PoolSource",
    	skipEvents = cms.untracked.uint32(0),

	fileNames = cms.untracked.vstring(
                                "/store/mc/RunIISummer16MiniAODv3/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/PUMoriond17_94X_mcRun2_asymptotic_v3_ext2-v1/270000/FEFC23FC-37C7-E811-97DA-0CC47AA53D86.root")
)

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(100)
)

process.load("EgammaIdCorrection.DifferentialEgammaIdCorrector.DifferentialEgammaIdCorrector_cfi")
import EgammaIdCorrection.DifferentialEgammaIdCorrector.DifferentialEgammaIdCorrector_cfi as EgammaIdInp_Corr
EgammaIdInp_Corr.setup_DifferentialPhoIdInputsCorrection( process )

process.dump = cms.EDAnalyzer("EventContentAnalyzer")
process.mytreewriter = cms.EDAnalyzer('EgammaNtuleWritter',
                                      PhoCorrectionContainer=cms.InputTag("differentialEgammaIdCorrection", "NewPhotonCorrections"),
				      rhoFixedGridCollection = cms.InputTag('fixedGridRhoAll')
                                      )


#process.out = cms.OutputModule("PoolOutputModule",
#                               fileName = cms.untracked.string('myOutputFile.root'),
#                               outputCommands = cms.untracked.vstring('drop *',
#                                                                      "keep *_reducedEgamma_*_*",
#                                                                      "keep *_differentialEgammaIdCorrection_*_*",
#                                                                      "keep *_slimmedOOTPhotons_*_*",
#                                                                      "keep *_slimmedPhotons_*_*")
#                               )

process.TFileService = cms.Service("TFileService",
    fileName = cms.string('EGOutputTree.root')
)

process.p = cms.Path(process.differentialEgammaIdCorrection*process.dump*process.mytreewriter)
#process.p = cms.Path(process.differentialEgammaIdCorrection*process.dump)
#process.e = cms.EndPath(process.out)

# customisation of the process.
#print process.dumpPython()
