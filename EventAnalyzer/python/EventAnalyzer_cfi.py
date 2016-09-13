import FWCore.ParameterSet.Config as cms

eventAnalyzer = cms.EDAnalyzer('EventAnalyzer',
    metLabel = cms.InputTag('slimmedMETs'),
    diphotonwithprotonLabel = cms.InputTag('flashggDiProtonsDiPhotons'),
    diphotonLabel = cms.InputTag('flashggDiPhotons'),
    minPtSinglePhoton = cms.double(50.),
    minMassDiPhoton = cms.double(500.),
)
