import FWCore.ParameterSet.Config as cms
def update(process, graphPath, inputTensorName, outputTensorName):

    from HLTrigger.Configuration.customizeHLTforPatatrack import customizeHLTforPatatrackTriplets
    process = customizeHLTforPatatrackTriplets(process)
    process.L2MLTestFilter = cms.EDFilter("L2TauNNTag",
        processName = cms.string('MLProva'),
        l1taus=cms.InputTag('hltGtStage2Digis','Tau'),
        ecalInputs =cms.VInputTag("hltEcalRecHit:EcalRecHitsEB", "hltEcalRecHit:EcalRecHitsEE"),
        hbheInput = cms.InputTag("hltHbhereco"),
        hfInput = cms.InputTag("hltHfreco"),
        hoInput = cms.InputTag("hltHoreco"),
        pataVertices = cms.InputTag("hltTrimmedPixelVertices"),
        pataTracks = cms.InputTag("hltPixelTracks")
    )
    process.L2MLTestFilter.graphPath = cms.string(graphPath)
    process.L2MLTestFilter.inputTensorName = cms.string(inputTensorName)
    process.L2MLTestFilter.outputTensorName = cms.string(outputTensorName)
    #process.MLPathTest = cms.Path(process.HLTBeginSequence + process.hltL1sDoubleTauBigOR , process.HLTDoLocalPixelTask, process.HLTRecoPixelTracksTask, process.HLTRecopixelvertexingTask)
    process.MLPathTest = cms.Path(process.HLTBeginSequence + process.hltL1sDoubleTauBigOR + process.HLTDoCaloSequence +process.L2MLTestFilter, process.HLTDoLocalPixelTask, process.HLTRecoPixelTracksTask, process.HLTRecopixelvertexingTask)

    process.schedule = cms.Schedule(*[ process.HLTriggerFirstPath, process.MLPathTest, process.HLTriggerFinalPath, process.endjob_step ], tasks=[process.patAlgosToolsTask])

    process.options.wantSummary = cms.untracked.bool(False)
    return process
