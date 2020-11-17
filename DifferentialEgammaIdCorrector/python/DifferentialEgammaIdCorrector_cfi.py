import FWCore.ParameterSet.Config as cms

globalVariables = cms.PSet(
    rho =  cms.InputTag('fixedGridRhoAll'),
    vertexes = cms.InputTag("offlineSlimmedPrimaryVertices"),

    puReWeight=cms.bool(False),
    puBins=cms.vdouble(),
    dataPu=cms.vdouble(),
    mcPu=cms.vdouble(),
    puInfo=cms.InputTag("slimmedAddPileupInfo"),
    extraFloats=cms.PSet(),
    )
differentialEgammaIdCorrection = cms.EDProducer("DifferentialEgammaIdCorrector",

	Electrons = cms.InputTag("slimmedElectrons"),
	Photons = cms.InputTag("slimmedPhotons"),

	rhoFixedGridCollection = cms.InputTag('fixedGridRhoAll'),
	reducedBarrelRecHitCollection = cms.InputTag('reducedEgamma','reducedEBRecHits'),
        reducedEndcapRecHitCollection = cms.InputTag('reducedEgamma','reducedEERecHits'),
        reducedPreshowerRecHitCollection = cms.InputTag('reducedEgamma','reducedESRecHits'),
	effAreasConfigFile = cms.FileInPath("RecoEgamma/PhotonIdentification/data/Spring16/effAreaPhotons_cone03_pfPhotons_90percentBased.txt"),
	globalVariables = globalVariables
	)


import json
import os
import EgammaIdCorrection.DifferentialEgammaIdCorrector.dumperConfigTools as cfgTools


def setup_DifferentialPhoIdInputsCorrection( process ):

    process.RandomNumberGeneratorService = cms.Service("RandomNumberGeneratorService",
                                                       differentialEgammaIdCorrection = cms.PSet(
                                                           initialSeed=cms.untracked.uint32(
                                                               90)
                                                       )
                                                   )

    process.load("EgammaIdCorrection.DifferentialEgammaIdCorrector.DifferentialEgammaIdCorrector_cfi")
    # print(process)
    
#    setattr(differentialEgammaIdCorrection, "photonIdMVAweightfile_EB", cms.FileInPath(str(metaConditions['flashggPhotons']['photonIdMVAweightfile_EB'])))
#    setattr(differentialEgammaIdCorrection, "photonIdMVAweightfile_EE", cms.FileInPath(str(metaConditions['flashggPhotons']['photonIdMVAweightfile_EE'])))
#    setattr(differentialEgammaIdCorrection, "effAreasConfigFile", cms.FileInPath(str(metaConditions['flashggPhotons']['effAreasConfigFile'])))
#    setattr(differentialEgammaIdCorrection, "is2017", cms.bool(metaConditions['flashggPhotons']['is2017']))
    corrections_summary = {}
    with open(os.path.expandvars("$CMSSW_BASE/src/EgammaIdCorrection/DifferentialEgammaIdCorrector/data/PhoIdInputsCorrections/corrections_summary_2017.json")) as json_file:
        corrections_summary = json.load(json_file)

        #---Shower shapes
    for var in corrections_summary['shower_shapes'].keys():
        for subdet in ['EB', 'EE']:
            ss_summary = corrections_summary['shower_shapes'][var][subdet]
            setattr(differentialEgammaIdCorrection, var+'_corrector_config_'+subdet, 
                    cms.PSet(
                        variables = cms.VPSet(),
                        weights = cms.FileInPath(str(ss_summary['weights'])),
                        regr_output_scaling = cms.string('x[0]*(%f)+(%f)' % (ss_summary['scale'], ss_summary['center'])),
                        regression = cms.bool(True),
                        classifier = cms.string('BDTG_{}_{}'.format(var, subdet))
                    )
                )        
            xgb_config = getattr(differentialEgammaIdCorrection, var+'_corrector_config_'+subdet)
            # print(list(metaConditions['PhoIdInputCorrections']['SS_variables']))
	    SS_variables = ["f0 := pt",
                          "f1 := superCluster.eta",
                          "f2 := phi",
                          "f3 := global.rho",
                          "f4 := superCluster.phiWidth",
                          "f5 := sieip",
                          "f6 := s4",
                          "f7 := full5x5_r9",
                          "f8 := full5x5_sigmaIetaIeta",
                          "f9 := superCluster.etaWidth"
                         ]
            cfgTools.addVariables(xgb_config.variables, [str(st) for st in SS_variables])
                                   

