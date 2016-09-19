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

#include "DiphotonAnalyzer/EventAnalyzer/interface/SelectionUtils.h"

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

      double sqrtS_;
      double singlePhotonMinPt_, singlePhotonMaxEta_, singlePhotonMinR9_;
      double photonPairMinMass_;

      TH1D* hMgg_notag_, *hMgg_tag_;
      TH1D* hPtgg_notag_, *hPtgg_tag_;
      TH1D* hYgg_notag_, *hYgg_tag_;
      TH1D* hXi1gg_notag_, *hXi1gg_tag_;
      TH1D* hXi2gg_notag_, *hXi2gg_tag_;
      TH2D* hMgg_vs_Mpp_, *hYgg_vs_Ypp_, *hXi1gg_vs_Xi1pp_, *hXi2gg_vs_Xi2pp_;
      TH1D* hMpp_;
      TH1D* hXi1_, *hXi2_;
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
  sqrtS_             (iConfig.getParameter<double>("sqrtS")),
  singlePhotonMinPt_ (iConfig.getParameter<double>("minPtSinglePhoton")),
  singlePhotonMaxEta_(iConfig.getParameter<double>("maxEtaSinglePhoton")),
  singlePhotonMinR9_ (iConfig.getParameter<double>("minR9SinglePhoton")),
  photonPairMinMass_ (iConfig.getParameter<double>("minMassDiPhoton"))
{
  //now do what ever initialization is needed
  usesResource("TFileService");
  edm::Service<TFileService> fs;
  hMgg_notag_ = fs->make<TH1D>("diphoton_m_notag", "Diphoton m (no proton tagging)", 750, 500., 2000.);
  hPtgg_notag_ = fs->make<TH1D>("diphoton_pt_notag", "Diphoton p_{T} (no proton tagging)", 160, 0., 400.);
  hYgg_notag_ = fs->make<TH1D>("diphoton_y_notag", "Diphoton rapidity (no proton tagging)", 40, -2.5, 2.5);
  hXi1gg_notag_ = fs->make<TH1D>("diphoton_xi1_notag", "#xi_{1} reconstructed from diphoton (no proton tagging)", 100, 0., 0.5);
  hXi2gg_notag_ = fs->make<TH1D>("diphoton_xi2_notag", "#xi_{2} reconstructed from diphoton (no proton tagging)", 100, 0., 0.5);
  hMgg_tag_ = fs->make<TH1D>("diphoton_m_tag", "Diphoton m (proton tagging)", 750, 500., 2000.);
  hPtgg_tag_ = fs->make<TH1D>("diphoton_pt_tag", "Diphoton p_{T} (proton tagging)", 160, 0., 400.);
  hYgg_tag_ = fs->make<TH1D>("diphoton_y_tag", "Diphoton rapidity (proton tagging)", 40, -2.5, 2.5);
  hXi1gg_tag_ = fs->make<TH1D>("diphoton_xi1_tag", "#xi_{1} reconstructed from diphoton (proton tagging)", 100, 0., 0.5);
  hXi2gg_tag_ = fs->make<TH1D>("diphoton_xi2_tag", "#xi_{2} reconstructed from diphoton (proton tagging)", 100, 0., 0.5);
  hMgg_vs_Mpp_ = fs->make<TH2D>("mgg_vs_mpp", "Diphoton m / diproton missing M", 300, 500., 2000., 300, 500., 2000.);
  hYgg_vs_Ypp_ = fs->make<TH2D>("ygg_vs_ypp", "Diphoton rapidity / diproton rapidity", 100, -5., 5., 100, -5., 5.);
  hXi1gg_vs_Xi1pp_ = fs->make<TH2D>("xi1gg_vs_xi1pp", "#xi_{1} (diphoton) / #xi_{1} (proton)", 100, 0., 0.5, 100, 0., 0.5);
  hXi2gg_vs_Xi2pp_ = fs->make<TH2D>("xi2gg_vs_xi2pp", "#xi_{2} (diphoton) / #xi_{2} (proton)", 100, 0., 0.5, 100, 0., 0.5);
  hMpp_ = fs->make<TH1D>("diproton_m", "Diproton missing mass", 750, 500., 2000.);
  hXi1_ = fs->make<TH1D>("proton1_xi", "#xi_{1} (RP proton)", 100, 0., 0.5);
  hXi2_ = fs->make<TH1D>("proton2_xi", "#xi_{2} (RP proton)", 100, 0., 0.5);
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

  // -----
  // collect the diphoton collection (no proton tag required here)
  // -----

  edm::Handle< edm::View<flashgg::DiPhotonCandidate> > diph;
  iEvent.getByToken(diphToken_, diph);

  unsigned int num_cand_notag = 0;
  for ( unsigned int i=0; i<diph->size(); i++ ) {
    edm::Ptr<flashgg::DiPhotonCandidate> pc = diph->ptrAt( i );

    if ( pc->leadPhotonId()<-0.9 ) continue;
    if ( pc->subLeadPhotonId()<-0.9 ) continue;

    if ( !passSinglePhotonCuts( pc->leadingPhoton() ) ) continue;
    if ( !passSinglePhotonCuts( pc->subLeadingPhoton() ) ) continue;

    if ( fabs( pc->leadingPhoton()->eta() )>=singlePhotonMaxEta_ or fabs( pc->subLeadingPhoton()->eta() )>=singlePhotonMaxEta_ ) continue;
    if ( pc->leadingPhoton()->pt()<singlePhotonMinPt_ or pc->subLeadingPhoton()->pt()<singlePhotonMinPt_ ) continue;
    if ( pc->leadingPhoton()->r9()<singlePhotonMinR9_ or pc->subLeadingPhoton()->r9()<singlePhotonMinR9_ ) continue;

    if ( pc->mass()<photonPairMinMass_ ) continue;

    hMgg_notag_->Fill( pc->mass() );
    hPtgg_notag_->Fill( pc->pt() );
    hYgg_notag_->Fill( pc->rapidity() );

    float xi1, xi2;
    computeXiReco( sqrtS_, pc.get(), &xi1, &xi2 );
    hXi1gg_notag_->Fill( xi1 );
    hXi2gg_notag_->Fill( xi2 );

    hMET_vs_Ptgg_notag->Fill( pc->pt(), met->sumEt() );
    num_cand_notag++;
  }
  hNum_diph_notag_->Fill( num_cand_notag );

  // -----
  // collect the diphoton + diproton tag
  // -----

  edm::Handle< edm::View<flashgg::DiProtonDiPhotonCandidate> > diphpr;
  iEvent.getByToken(diphprToken_, diphpr);

  unsigned int num_gg_cand_tag = 0;
  for ( unsigned int i=0; i<diphpr->size(); i++ ) {
    edm::Ptr<flashgg::DiProtonDiPhotonCandidate> pc = diphpr->ptrAt( i );

    if ( pc->diphoton()->leadPhotonId()<-0.9 ) continue;
    if ( pc->diphoton()->subLeadPhotonId()<-0.9 ) continue;

    if ( !passSinglePhotonCuts( pc->diphoton()->leadingPhoton() ) ) continue;
    if ( !passSinglePhotonCuts( pc->diphoton()->subLeadingPhoton() ) ) continue;

    if ( fabs( pc->diphoton()->leadingPhoton()->eta() )>=singlePhotonMaxEta_ or fabs( pc->diphoton()->subLeadingPhoton()->eta() )>=singlePhotonMaxEta_ ) continue;
    if ( pc->diphoton()->leadingPhoton()->pt()<singlePhotonMinPt_ or pc->diphoton()->subLeadingPhoton()->pt()<singlePhotonMinPt_ ) continue;
    if ( pc->diphoton()->leadingPhoton()->r9()<singlePhotonMinR9_ or pc->diphoton()->subLeadingPhoton()->r9()<singlePhotonMinR9_ ) continue;

    if ( pc->diphoton()->mass()<photonPairMinMass_ ) continue;

    hMgg_tag_->Fill( pc->diphoton()->mass() );
    hPtgg_tag_->Fill( pc->diphoton()->pt() );
    hYgg_tag_->Fill( pc->diphoton()->rapidity() );

    float xi1, xi2;
    computeXiReco( sqrtS_, pc->diphoton(), &xi1, &xi2 );
    hXi1gg_tag_->Fill( xi1 );
    hXi2gg_tag_->Fill( xi2 );

    hXi1gg_vs_Xi1pp_->Fill( pc->diproton()->proton1()->xi(), xi1 );
    hXi2gg_vs_Xi2pp_->Fill( pc->diproton()->proton2()->xi(), xi2 );

    hMET_vs_Ptgg_tag->Fill( pc->diphoton()->pt(), met->sumEt() );

    hMgg_vs_Mpp_->Fill( pc->diproton()->mass(), pc->diphoton()->mass() );
    hYgg_vs_Ypp_->Fill( pc->diproton()->rapidity(), pc->diphoton()->rapidity() );

    hXi1_->Fill( pc->diproton()->proton1()->xi() );
    hXi2_->Fill( pc->diproton()->proton2()->xi() );

    hMpp_->Fill( pc->diproton()->mass() );

    num_gg_cand_tag++;

  }
  hNum_diph_tag_->Fill( num_gg_cand_tag );

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
