# DiphotonAnalyzer

Simple EDAnalyzer for flashgg output

### Prerequisites
- Flashgg
- ...

### (Simplified set of) instructions
- Clone this directory inside your $CMSSW_BASE/src directory containing the flashgg installation (with cmsenv already run)
- scram b
- Edit EventAnalyzer/test/EventAnalyzer_cfg.py accordingly (cuts, input files, any future feature...)
- cmsRun EventAnalyzer/test/EventAnalyzer_cfg.py
