import FWCore.ParameterSet.Config as cms

process = cms.Process("analyzer")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
#        'file:myfile.root'
'/store/group/phys_higgs/cmshgg/lforthom/flashgg/pps_run2016/Moriond16WSFinal-106-g90923ae/DoubleEG/pps_run2016-Moriond16WSFinal-106-g90923ae-v0-Run2016B-PromptReco-v2/160624_003754/0000/myMicroAODOutputFile_1.root',
'/store/group/phys_higgs/cmshgg/lforthom/flashgg/pps_run2016/Moriond16WSFinal-106-g90923ae/DoubleEG/pps_run2016-Moriond16WSFinal-106-g90923ae-v0-Run2016B-PromptReco-v2/160624_003754/0000/myMicroAODOutputFile_2.root',
    )
)

process.load('DiphotonAnalyzer.EventAnalyzer.EventAnalyzer_cfi')

# set some parameters to the run
process.eventAnalyzer.minPtSinglePhoton = cms.double(75.)
process.eventAnalyzer.minMassDiPhoton = cms.double(500.)

# prepare the output file
process.TFileService = cms.Service("TFileService", 
    fileName = cms.string("histo.root"),
    closeFileFast = cms.untracked.bool(True)
)

process.p = cms.Path(process.eventAnalyzer)
