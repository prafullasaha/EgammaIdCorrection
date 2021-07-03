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
#	PFCandidatesTag=cms.untracked.InputTag('packedPFCandidates'),	
	PFCandidatesTag=cms.InputTag('packedPFCandidates'),
	vertexes = cms.InputTag("offlineSlimmedPrimaryVertices"),
	rhoFixedGridCollection = cms.InputTag('fixedGridRhoAll'),
	reducedBarrelRecHitCollection = cms.InputTag('reducedEgamma','reducedEBRecHits'),
        reducedEndcapRecHitCollection = cms.InputTag('reducedEgamma','reducedEERecHits'),
        reducedPreshowerRecHitCollection = cms.InputTag('reducedEgamma','reducedESRecHits'),
	effAreasConfigFile = cms.FileInPath("RecoEgamma/PhotonIdentification/data/Spring16/effAreaPhotons_cone03_pfPhotons_90percentBased.txt"),
	globalVariables = globalVariables,

        correctShowerShapes = cms.bool(True),
        correctIsolations = cms.bool(True),
	useVtx0ForNeutralIso = cms.bool(True),	

        phoIso_corrector_config_EB = cms.PSet(),
        phoIso_corrector_config_EE = cms.PSet(),
        chIso_corrector_config_EB = cms.PSet(),
        chIso_corrector_config_EE = cms.PSet(),
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
    with open(os.path.expandvars("$CMSSW_BASE/src/EgammaIdCorrection/DifferentialEgammaIdCorrector/data/PhoIdInputsCorrections/corrections_summary_2016.json")) as json_file:
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
                                   
        #---Isolations
    for subdet in ['EB', 'EE']:
        #-photon isolation
        iso_summary = corrections_summary['isolations']['phoIso'][subdet]
        setattr(differentialEgammaIdCorrection, 'phoIso_corrector_config_'+subdet, 
                cms.PSet(
                    clf_data = cms.PSet(
                        variables = cms.VPSet(),
                        weights = cms.FileInPath(str(iso_summary['peak_tail_clfs']['weights_data'])),
                        classifier = cms.string('BDTG_{}_{}_p0t_data'.format("phoIso", subdet))
                    ),
                    clf_mc = cms.PSet(
                        variables = cms.VPSet(),
                        weights = cms.FileInPath(str(iso_summary['peak_tail_clfs']['weights_mc'])),
                        classifier = cms.string('BDTG_{}_{}_p0t_mc'.format("phoIso", subdet))
                    ),
                    peak2tail = cms.PSet(
                        variables = cms.VPSet(),
                        weights = cms.FileInPath(str(iso_summary['peak2tail']['weights'])),
                        regression = cms.bool(True),
                        classifier = cms.string('BDTG_{}_{}_p2t'.format("phoIso", subdet))
                    ),
                    morphing = cms.PSet(
                        variables = cms.VPSet(),
                        weights = cms.FileInPath(str(iso_summary['morphing']['weights'])),
                        regression = cms.bool(True),
                        classifier = cms.string('BDTG_{}_{}'.format("phoIso", subdet)),
                        regr_output_scaling = cms.string('x[0]*(%f)+(%f)' % (iso_summary['morphing']['scale'], iso_summary['morphing']['center']))
                    )
                )
        )
        phoIso_corrector_config = getattr(differentialEgammaIdCorrection, 'phoIso_corrector_config_'+subdet)
        cfgTools.addVariables(phoIso_corrector_config.clf_data.variables, [str(x) for x in iso_summary['peak_tail_clfs']['variables']])
        cfgTools.addVariables(phoIso_corrector_config.clf_mc.variables, [str(x) for x in iso_summary['peak_tail_clfs']['variables']])
        cfgTools.addVariables(phoIso_corrector_config.peak2tail.variables, [str(x) for x in iso_summary['peak2tail']['variables']])
        cfgTools.addVariables(phoIso_corrector_config.morphing.variables, [str(x) for x in iso_summary['morphing']['variables']])    

        #-charged isolations
        iso_summary = corrections_summary['isolations']['chIso'][subdet]
        setattr(differentialEgammaIdCorrection, 'chIso_corrector_config_'+subdet, 
                cms.PSet(
                    clf_data = cms.PSet(
                        variables = cms.VPSet(),
                        weights = cms.FileInPath(str(iso_summary['peak_tail_clfs']['weights_data'])),
                        multiclass = cms.bool(True),
                        classifier = cms.string('BDTG_{}_{}_3Cat_data'.format("chIso", subdet))
                    ),
                    clf_mc = cms.PSet(
                        variables = cms.VPSet(),
                        weights = cms.FileInPath(str(iso_summary['peak_tail_clfs']['weights_mc'])),
                        multiclass = cms.bool(True),
                        classifier = cms.string('BDTG_{}_{}_3Cat_mc'.format("chIso", subdet))
                    ),
                    chIso_peak2tail = cms.PSet(
                        variables = cms.VPSet(),
                        weights = cms.FileInPath(str(iso_summary['chIso_peak2tail']['weights'])),
                        regression = cms.bool(True),
                        classifier = cms.string('BDTG_{}_{}_p2t'.format("chIso", subdet))
                    ),
                    chIso_morphing = cms.PSet(
                        variables = cms.VPSet(),
                        weights = cms.FileInPath(str(iso_summary['chIso_morphing']['weights'])),
                        classifier = cms.string('BDTG_{}_{}'.format("chIso", subdet)),
                        regression = cms.bool(True),
                        regr_output_scaling = cms.string('x[0]*(%f)+(%f)' % (iso_summary['chIso_morphing']['scale'], 
                                                                             iso_summary['chIso_morphing']['center']))
                    ),
                    chIsoWorst_peak2tail = cms.PSet(
                        variables = cms.VPSet(),
                        weights = cms.FileInPath(str(iso_summary['chIsoWorst_peak2tail']['weights'])),
                        regression = cms.bool(True),
                        classifier = cms.string('BDTG_{}_{}_p2t'.format("chIsoWorst", subdet))
                    ),
                    chIsoWorst_morphing = cms.PSet(
                        variables = cms.VPSet(),
                        weights = cms.FileInPath(str(iso_summary['chIsoWorst_morphing']['weights'])),
                        classifier = cms.string('BDTG_{}_{}'.format("chIsoWorst", subdet)),
                        regression = cms.bool(True),
                        regr_output_scaling = cms.string('x[0]*(%f)+(%f)' % (iso_summary['chIsoWorst_morphing']['scale'], 
                                                                             iso_summary['chIsoWorst_morphing']['center']))
                    )
                )
            )
        chIso_corrector_config = getattr(differentialEgammaIdCorrection, 'chIso_corrector_config_'+subdet)
        cfgTools.addVariables(chIso_corrector_config.clf_data.variables, [str(x) for x in iso_summary['peak_tail_clfs']['variables']])
        cfgTools.addVariables(chIso_corrector_config.clf_mc.variables, [str(x) for x in iso_summary['peak_tail_clfs']['variables']])
        cfgTools.addVariables(chIso_corrector_config.chIso_peak2tail.variables, [str(x) for x in iso_summary['chIso_peak2tail']['variables']])
        cfgTools.addVariables(chIso_corrector_config.chIso_morphing.variables, [str(x) for x in iso_summary['chIso_morphing']['variables']])    
        cfgTools.addVariables(chIso_corrector_config.chIsoWorst_peak2tail.variables, [str(x) for x in iso_summary['chIsoWorst_peak2tail']['variables']])
        cfgTools.addVariables(chIso_corrector_config.chIsoWorst_morphing.variables, [str(x) for x in iso_summary['chIsoWorst_morphing']['variables']])    
    
