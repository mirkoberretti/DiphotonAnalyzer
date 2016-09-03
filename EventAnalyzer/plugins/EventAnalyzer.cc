// -*- C++ -*-
//
// Package:    DiphotonAnalyzer/EventAnalyzer
// Class:      EventAnalyzer
// 
/**\class EventAnalyzer EventAnalyzer.cc DiphotonAnalyzer/EventAnalyzer/plugins/EventAnalyzer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Laurent Forthomme
//         Created:  Thu, 01 Sep 2016 16:57:56 GMT
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/PatCandidates/interface/MET.h"
#include "flashgg/DataFormats/interface/DiPhotonCandidate.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "TH2.h"

//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<> and also remove the line from
// constructor "usesResource("TFileService");"
// This will improve performance in multithreaded jobs.

class EventAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
   public:
      explicit EventAnalyzer(const edm::ParameterSet&);
      ~EventAnalyzer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

      // ----------member data ---------------------------

      edm::EDGetTokenT< edm::View<flashgg::DiPhotonCandidate> > diphotonToken_;
      edm::EDGetTokenT< edm::View<pat::MET> > metToken_;
      double singlePhotonMinPt_;
      double photonPairMinMass_;

      TH2D* hDiphotonPtVsMET_;
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
EventAnalyzer::EventAnalyzer(const edm::ParameterSet& iConfig) :
  diphotonToken_(consumes< edm::View<flashgg::DiPhotonCandidate> >(iConfig.getParameter<edm::InputTag>("diphotonLabel"))),
  metToken_     (consumes< edm::View<pat::MET> >                  (iConfig.getParameter<edm::InputTag>("metLabel"))),
  singlePhotonMinPt_(iConfig.getParameter<double>("minPtSinglePhoton")),
  photonPairMinMass_(iConfig.getParameter<double>("minMassDiPhoton"))
{
  //now do what ever initialization is needed
  usesResource("TFileService");
  edm::Service<TFileService> fs;
  hDiphotonPtVsMET_ = fs->make<TH2D>("diphoton_pt_vs_met", "Diphoton p_{T} / MET", 100, 0., 100., 100, 0., 100.);
}


EventAnalyzer::~EventAnalyzer()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
EventAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  #include <iostream> // for debugging purposes

  // first start by fetching the diphoton collection from EDM file
  edm::Handle< edm::View<flashgg::DiPhotonCandidate> > diphotons;
  iEvent.getByToken(diphotonToken_, diphotons);

  unsigned int num_diphoton_cand = 0;
  flashgg::DiPhotonCandidate gamgam_cand;
  for (edm::View<flashgg::DiPhotonCandidate>::const_iterator diph=diphotons->begin(); diph!=diphotons->end(); diph++) {
    if (diph->leadingPhoton()->pt()<singlePhotonMinPt_) continue;
    if (diph->subLeadingPhoton()->pt()<singlePhotonMinPt_) continue;
    if (diph->mass()<photonPairMinMass_) continue;

    std::cout << "---> diphoton " << num_diphoton_cand << "::: mass=" << diph->mass() << ", pt=" << diph->pt() << std::endl
              << "     leading photon: pt=" << diph->leadingPhoton()->pt() << std::endl
              << "  subleading photon: pt=" << diph->subLeadingPhoton()->pt() << std::endl;
    num_diphoton_cand++;
    gamgam_cand = *diph;
  }
  std::cout << num_diphoton_cand << " diphoton candidate(s) in the event" << std::endl;

  // retrieve the missing ET
  edm::Handle< edm::View<pat::MET> > mets;
  iEvent.getByToken(metToken_, mets);

  const edm::View<pat::MET>* metColl = mets.product();
  edm::View<pat::MET>::const_iterator met = metColl->begin();
  std::cout << "met=" << met->sumEt() << std::endl;

  if (num_diphoton_cand==1) { //FIXME study the multiple diphoton cases
    hDiphotonPtVsMET_->Fill(gamgam_cand.pt(), met->sumEt());
  }

}


// ------------ method called once each job just before starting event loop  ------------
void 
EventAnalyzer::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
EventAnalyzer::endJob() 
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
EventAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(EventAnalyzer);
