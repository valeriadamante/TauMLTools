import FWCore.ParameterSet.Config as cms
def update(process, graphPath, normalizationDict, rateWP):

    from HLTrigger.Configuration.customizeHLTforPatatrack import customizeHLTforPatatrackTriplets
    process = customizeHLTforPatatrackTriplets(process)
    process.load("TauMLTools.Production.l2TauNNTag_cfi")
    process.l2TauNNTag.processName = cms.string('MLProva')
    process.l2TauNNTag.l1taus=cms.InputTag('hltGtStage2Digis','Tau')
    process.l2TauNNTag.ecalInputs =cms.VInputTag("hltEcalRecHit:EcalRecHitsEB", "hltEcalRecHit:EcalRecHitsEE")
    process.l2TauNNTag.hbheInput = cms.InputTag("hltHbhereco")
    process.l2TauNNTag.hoInput = cms.InputTag("hltHoreco")
    process.l2TauNNTag.pataVertices = cms.InputTag("hltTrimmedPixelVertices")
    process.l2TauNNTag.pataTracks = cms.InputTag("hltPixelTracks")
    process.l2TauNNTag.graphPath = cms.string(graphPath)
    process.l2TauNNTag.normalizationDict = cms.string(normalizationDict)
    thWp = {
        'opt_threshold_3': 0.16547765582618013,
        'opt_threshold_4': 0.11188014969047799,
        'opt_threshold_5': 0.07665744051155343,
    }
    process.l2TauNNTag.discr_threshold = cms.double(thWp[rateWP])
    #process.MLPathTest = cms.Path(process.HLTBeginSequence + process.hltL1sDoubleTauBigOR , process.HLTDoLocalPixelTask, process.HLTRecoPixelTracksTask, process.HLTRecopixelvertexingTask)
    process.MLPathTest = cms.Path(process.HLTBeginSequence + process.hltL1sDoubleTauBigOR + process.HLTDoCaloSequence +process.l2TauNNTag, process.HLTDoLocalPixelTask, process.HLTRecoPixelTracksTask, process.HLTRecopixelvertexingTask)

    process.schedule = cms.Schedule(*[ process.HLTriggerFirstPath, process.MLPathTest, process.HLTriggerFinalPath, process.endjob_step ], tasks=[process.patAlgosToolsTask])

    process.options.wantSummary = cms.untracked.bool(True)
    return process
