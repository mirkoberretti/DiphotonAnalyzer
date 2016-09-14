import FWCore.ParameterSet.Config as cms

treeProducer = cms.EDAnalyzer('TreeProducer',
    metLabel = cms.InputTag('slimmedMETs'),
    protonLabel = cms.InputTag('flashggProtons'),
    photonLabel = cms.InputTag('flashggRandomizedPhotons'),
    minPtSinglePhoton = cms.double(50.),
    minMassDiPhoton = cms.double(500.),
    outputFilename = cms.string('output.root'),
)
