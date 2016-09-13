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
#include "flashgg/DataFormats/interface/DiProtonDiPhotonCandidate.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "TH1.h"
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

      edm::EDGetTokenT< edm::View<flashgg::DiPhotonCandidate> > diphToken_;
      edm::EDGetTokenT< edm::View<flashgg::DiProtonDiPhotonCandidate> > diphprToken_;
      edm::EDGetTokenT< edm::View<pat::MET> > metToken_;
      double singlePhotonMinPt_;
      double photonPairMinMass_;

      TH1D* hMgg_notag_, *hMgg_tag_;
      TH1D* hPtgg_notag_, *hPtgg_tag_;
      TH1D* hYgg_notag_, *hYgg_tag_;
      TH2D* hMgg_vs_Mpp_, *hYgg_vs_Ypp_;
      TH1D* hMpp_;
      TH2D* hMET_vs_Ptgg_notag, *hMET_vs_Ptgg_tag;
      TH1D* hNum_diph_notag_, *hNum_diph_tag_;
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
  diphToken_  (consumes< edm::View<flashgg::DiPhotonCandidate> >        (iConfig.getParameter<edm::InputTag>("diphotonLabel"))),
  diphprToken_(consumes< edm::View<flashgg::DiProtonDiPhotonCandidate> >(iConfig.getParameter<edm::InputTag>("diphotonwithprotonLabel"))),
  metToken_   (consumes< edm::View<pat::MET> >                          (iConfig.getParameter<edm::InputTag>("metLabel"))),
  singlePhotonMinPt_(iConfig.getParameter<double>("minPtSinglePhoton")),
  photonPairMinMass_(iConfig.getParameter<double>("minMassDiPhoton"))
{
  //now do what ever initialization is needed
  usesResource("TFileService");
  edm::Service<TFileService> fs;
  hMgg_notag_ = fs->make<TH1D>("diphoton_m_notag", "Diphoton m (no proton tagging)", 750, 500., 2000.);
  hPtgg_notag_ = fs->make<TH1D>("diphoton_pt_notag", "Diphoton p_{T} (no proton tagging)", 160, 0., 400.);
  hYgg_notag_ = fs->make<TH1D>("diphoton_y_notag", "Diphoton rapidity (no proton tagging)", 20, -5., 5.);
  hMgg_tag_ = fs->make<TH1D>("diphoton_m_tag", "Diphoton m (proton tagging)", 750, 500., 2000.);
  hPtgg_tag_ = fs->make<TH1D>("diphoton_pt_tag", "Diphoton p_{T} (proton tagging)", 160, 0., 400.);
  hYgg_tag_ = fs->make<TH1D>("diphoton_y_tag", "Diphoton rapidity (proton tagging)", 20, -5., 5.);
  hMgg_vs_Mpp_ = fs->make<TH2D>("mgg_vs_mpp", "Diphoton m / diproton missing M", 300, 500., 2000., 300, 500., 2000.);
  hYgg_vs_Ypp_ = fs->make<TH2D>("ygg_vs_ypp", "Diphoton rapidity / diproton rapidity", 20, -5., 5., 20, -5., 5.);
  hMpp_ = fs->make<TH1D>("diproton_m", "Diproton missing mass", 750, 500., 2000.);
  hMET_vs_Ptgg_notag = fs->make<TH2D>("diphoton_pt_vs_met_notag", "Diphoton p_{T} / MET (no proton tagging)", 160, 0., 400., 160, 0., 400.);
  hMET_vs_Ptgg_tag = fs->make<TH2D>("diphoton_pt_vs_met_tag", "Diphoton p_{T} / MET (proton tagging)", 160, 0., 400., 160, 0., 400.);
  hNum_diph_notag_ = fs->make<TH1D>("num_diphoton_evt_notag", "Number of diphoton candidates per event (no proton tagging)", 10, 0., 10.);
  hNum_diph_tag_ = fs->make<TH1D>("num_diphoton_evt_tag", "Number of diphoton candidates per event (proton tagging)", 10, 0., 10.);
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

  // retrieve the missing ET
  edm::Handle< edm::View<pat::MET> > mets;
  iEvent.getByToken(metToken_, mets);

  const edm::View<pat::MET>* metColl = mets.product();
  edm::View<pat::MET>::const_iterator met = metColl->begin();

 /*if (num_diphoton_cand==1) { //FIXME study the multiple diphoton cases
    hDiphotonPtVsMET_->Fill(gamgam_cand.pt(), met->sumEt());
  }*/

  // fetch the diphoton collection from EDM file
  edm::Handle< edm::View<flashgg::DiPhotonCandidate> > diph;
  iEvent.getByToken(diphToken_, diph);

  unsigned int num_cand_notag = 0;
  for ( unsigned int i=0; i<diph->size(); i++ ) {
    edm::Ptr<flashgg::DiPhotonCandidate> pc = diph->ptrAt( i );

    if ( pc->leadingPhoton()->pt()<singlePhotonMinPt_ ) continue;
    if ( pc->subLeadingPhoton()->pt()<singlePhotonMinPt_ ) continue;
    if ( pc->mass()<photonPairMinMass_ ) continue;

    if ( fabs( pc->leadingPhoton()->superCluster()->eta())>=1.4442 and fabs( pc->leadingPhoton()->superCluster()->eta() )<=1.566 or fabs( pc->leadingPhoton()->superCluster()->eta() )>=2.5
      or fabs( pc->subLeadingPhoton()->superCluster()->eta())>=1.4442 and fabs( pc->subLeadingPhoton()->superCluster()->eta())<=1.566 or fabs( pc->subLeadingPhoton()->superCluster()->eta())>=2.5 ) continue;

    if ( max( pc->leadingPhoton()->r9(), pc->subLeadingPhoton()->r9() )<0.94 ) continue;

    hMgg_notag_->Fill( pc->mass() );
    hPtgg_notag_->Fill( pc->pt() );
    hYgg_notag_->Fill( pc->rapidity() );

    hMET_vs_Ptgg_notag->Fill( pc->pt(), met->sumEt() );
    num_cand_notag++;
  }
  hNum_diph_notag_->Fill( num_cand_notag );

  edm::Handle< edm::View<flashgg::DiProtonDiPhotonCandidate> > diphpr;
  iEvent.getByToken(diphprToken_, diphpr);

  unsigned int num_gg_cand_tag = 0;
  for ( unsigned int i=0; i<diphpr->size(); i++ ) {
    edm::Ptr<flashgg::DiProtonDiPhotonCandidate> pc = diphpr->ptrAt( i );

    if ( pc->diphoton()->leadingPhoton()->pt()<singlePhotonMinPt_ ) continue;
    if ( pc->diphoton()->subLeadingPhoton()->pt()<singlePhotonMinPt_ ) continue;
    if ( pc->diphoton()->mass()<photonPairMinMass_ ) continue;

    if ( fabs( pc->diphoton()->leadingPhoton()->superCluster()->eta())>=1.4442 and fabs( pc->diphoton()->leadingPhoton()->superCluster()->eta() )<=1.566 or fabs( pc->diphoton()->leadingPhoton()->superCluster()->eta() )>=2.5
      or fabs( pc->diphoton()->subLeadingPhoton()->superCluster()->eta())>=1.4442 and fabs( pc->diphoton()->subLeadingPhoton()->superCluster()->eta())<=1.566 or fabs( pc->diphoton()->subLeadingPhoton()->superCluster()->eta())>=2.5 ) continue;

    if ( max( pc->diphoton()->leadingPhoton()->r9(), pc->diphoton()->subLeadingPhoton()->r9() )<0.94 ) continue;
    //std::cout << "protons: " << pc->diphoton()->mass() << " // " << pc->diproton()->mass() << std::endl;

    hMgg_tag_->Fill( pc->diphoton()->mass() );
    hPtgg_tag_->Fill( pc->diphoton()->pt() );
    hYgg_tag_->Fill( pc->diphoton()->rapidity() );

    hMET_vs_Ptgg_tag->Fill( pc->diphoton()->pt(), met->sumEt() );

    hMgg_vs_Mpp_->Fill( pc->diproton()->mass(), pc->diphoton()->mass() );
    hYgg_vs_Ypp_->Fill( pc->diproton()->rapidity(), pc->diphoton()->rapidity() );

    hMpp_->Fill( pc->diproton()->mass() );

    num_gg_cand_tag++;

  }
  hNum_diph_tag_->Fill( num_gg_cand_tag );

  /*flashgg::DiPhotonCandidate gamgam_cand;
  for (edm::View<flashgg::DiPhotonCandidate>::const_iterator diph=diphotons->begin(); diph!=diphotons->end(); diph++) {

    std::cout << "---> diphoton " << num_diphoton_cand << "::: mass=" << diph->mass() << ", pt=" << diph->pt() << std::endl
              << "     leading photon: pt=" << diph->leadingPhoton()->pt() << std::endl
              << "  subleading photon: pt=" << diph->subLeadingPhoton()->pt() << std::endl;
    num_diphoton_cand++;
    gamgam_cand = *diph;
  }*/
  //std::cout << num_diphoton_cand << " diphoton candidate(s) in the event" << std::endl;


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
