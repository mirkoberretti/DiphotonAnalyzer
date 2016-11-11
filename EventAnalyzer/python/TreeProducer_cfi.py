import FWCore.ParameterSet.Config as cms

treeProducer = cms.EDAnalyzer('TreeProducer',
    sqrtS = cms.double(13.e3),
    metLabel = cms.InputTag('slimmedMETs'),
    diphotonLabel = cms.InputTag('flashggDiPhotons'),
    protonLabel = cms.InputTag('flashggProtons'),
    vertexLabel = cms.InputTag('offlineSlimmedPrimaryVertices'),
    electronLabel = cms.InputTag('flashggElectrons'),
    muonLabel = cms.InputTag('flashggMuons'),
    jetLabel = cms.InputTag('flashggJets'),
    minPtSinglePhoton = cms.double(50.),
    minR9SinglePhoton = cms.double(0.94),
    maxEtaSinglePhoton = cms.double(2.5),
    minMassDiPhoton = cms.double(500.),
    outputFilename = cms.string('output.root'),
)
