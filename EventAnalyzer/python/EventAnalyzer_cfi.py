import FWCore.ParameterSet.Config as cms

eventAnalyzer = cms.EDAnalyzer('EventAnalyzer',
    sqrtS = cms.double(13.e3),
    metLabel = cms.InputTag('slimmedMETs'),
    diphotonwithprotonLabel = cms.InputTag('flashggDiProtonsDiPhotons'),
    diphotonLabel = cms.InputTag('flashggDiPhotons'),
    minPtSinglePhoton = cms.double(50.),
    minR9SinglePhoton = cms.double(0.94),
    maxEtaSinglePhoton = cms.double(2.5),
    minMassDiPhoton = cms.double(500.),
)
