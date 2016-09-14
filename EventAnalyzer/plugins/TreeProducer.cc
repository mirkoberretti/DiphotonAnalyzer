// -*- C++ -*-
//
// Package:    DiphotonAnalyzer/EventAnalyzer
// Class:      TreeProducer
// 
/**\class TreeProducer TreeProducer.cc DiphotonAnalyzer/EventAnalyzer/plugins/TreeProducer.cc

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

#include "TFile.h"
#include "TTree.h"
#include "TLorentzVector.h"

//
// class declaration
//

#define MAX_PHOTON 20
#define MAX_DIPHOTON 5
#define MAX_PROTON 10
#define MAX_DIPROTON 5

class TreeProducer : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
  public:
    explicit TreeProducer(const edm::ParameterSet&);
    ~TreeProducer();

    static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


  private:
    virtual void beginJob() override;
    virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
    virtual void endJob() override;

    // ----------member data ---------------------------

    edm::EDGetTokenT< edm::View<flashgg::Photon> > photonToken_;
    edm::EDGetTokenT< edm::View<flashgg::Proton> > protonToken_;
    edm::EDGetTokenT< edm::View<pat::MET> > metToken_;
    double singlePhotonMinPt_;
    double photonPairMinMass_;

    std::string filename_;
    TFile* file_;
    TTree* tree_;

    // --- tree components ---
    unsigned int fPhotonNum;
    float fPhotonPt[MAX_PHOTON], fPhotonEta[MAX_PHOTON], fPhotonEtaSC[MAX_PHOTON], fPhotonPhi[MAX_PHOTON];
    float fPhotonR9[MAX_PHOTON];

    unsigned int fProtonNum;
    float fProtonXi[MAX_PROTON];
    unsigned int fProtonSide[MAX_PROTON];

    unsigned int fDiphotonNum;
    unsigned int fDiphotonPhoton1[MAX_DIPHOTON], fDiphotonPhoton2[MAX_DIPHOTON];
    float fDiphotonM[MAX_DIPHOTON], fDiphotonY[MAX_DIPHOTON];
    float fDiphotonPt[MAX_DIPHOTON], fDiphotonDphi[MAX_DIPHOTON];

    unsigned int fDiprotonNum;
    unsigned int fDiprotonProton1[MAX_DIPROTON], fDiprotonProton2[MAX_DIPROTON];
    float fDiprotonM[MAX_DIPROTON], fDiprotonY[MAX_DIPROTON];

    float fMET;

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
TreeProducer::TreeProducer(const edm::ParameterSet& iConfig) :
  photonToken_( consumes< edm::View<flashgg::Photon> >( iConfig.getParameter<edm::InputTag>( "photonLabel" ) ) ),
  protonToken_( consumes< edm::View<flashgg::Proton> >( iConfig.getParameter<edm::InputTag>( "protonLabel") ) ),
  metToken_   ( consumes< edm::View<pat::MET> >       ( iConfig.getParameter<edm::InputTag>( "metLabel") ) ),
  singlePhotonMinPt_( iConfig.getParameter<double>( "minPtSinglePhoton" ) ),
  photonPairMinMass_( iConfig.getParameter<double>( "minMassDiPhoton" ) ),
  filename_         ( iConfig.getParameter<std::string>( "outputFilename" ) ),
  file_( 0 ), tree_( 0 )
{
  //now do what ever initialization is needed
  file_ = new TFile( filename_.c_str(), "recreate" );
}


TreeProducer::~TreeProducer()
{
 
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)

  if ( file_ ) {
    file_->Write();
    delete file_;
  }
  if ( tree_ ) delete tree_;

}


//
// member functions
//

// ------------ method called for each event  ------------
void
TreeProducer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  #include <iostream> // for debugging purposes

  // fetch the diphoton collection from EDM file
  edm::Handle< edm::View<flashgg::Photon> > photons;
  iEvent.getByToken(photonToken_, photons);

  fPhotonNum = 0;
  for ( unsigned int i=0; i<photons->size(); i++ ) {

    edm::Ptr<flashgg::Photon> photon = photons->ptrAt( i );
    fPhotonPt[i] = photon->pt();
    fPhotonEta[i] = photon->eta();
    fPhotonEtaSC[i] = photon->superCluster()->eta();
    fPhotonPhi[i] = photon->phi();
    fPhotonR9[i] = photon->r9();

    fPhotonNum++;
  }

  fDiphotonNum = 0;
  TLorentzVector ph1, ph2;
  for ( unsigned int i=0; i<fPhotonNum; i++ ) {

    // single photons' quality cuts before diphoton candidates matching
    if ( ( fabs( fPhotonEtaSC[i] )>=1.4442 and fabs( fPhotonEtaSC[i] )<=1.566 ) or fabs( fPhotonEtaSC[i] )>=2.5 ) continue;
    if ( fPhotonR9[i]<0.94 ) continue;
    if ( fPhotonPt[i]<singlePhotonMinPt_ ) continue;
    ph1.SetPtEtaPhiM( fPhotonPt[i], fPhotonEta[i], fPhotonPhi[i], 0. );

    for ( unsigned int j=i+1; j<fPhotonNum; j++ ) {
      if ( ( fabs( fPhotonEtaSC[j] )>=1.4442 and fabs( fPhotonEtaSC[j] )<=1.566 ) or fabs( fPhotonEtaSC[j] )>=2.5 ) continue;
      if ( fPhotonR9[j]<0.94 ) continue;
      if ( fPhotonPt[j]<singlePhotonMinPt_ ) continue;
      ph2.SetPtEtaPhiM( fPhotonPt[j], fPhotonEta[j], fPhotonPhi[j], 0. );

      if ( (ph1+ph2).M()<photonPairMinMass_ ) continue;

      fDiphotonPhoton1[fDiphotonNum] = i;
      fDiphotonPhoton2[fDiphotonNum] = j;
      fDiphotonNum++;
    
    }
  }

  edm::Handle< edm::View<flashgg::Proton> > protons;
  iEvent.getByToken(protonToken_, protons);

  fProtonNum = 0;
  for ( unsigned int i=0; i<protons->size(); i++ ) {
    edm::Ptr<flashgg::Proton> proton = protons->ptrAt( i );
    fProtonXi[i] = proton->xi();
    fProtonSide[i] = proton->side();
    fProtonNum++;
  }

  for ( unsigned int i=0; i<fProtonNum; i++ ) {
    for ( unsigned int j=i+1; j<fProtonNum; j++ ) {
      if ( fProtonSide[i]==fProtonSide[j] ) continue;
      fDiprotonM[fDiprotonNum] = 13.e3*sqrt( fProtonXi[i]*fProtonXi[j] );
      fDiprotonY[fDiprotonNum] = log( fProtonXi[j]/fProtonXi[i] )/2.;
      fDiprotonNum++;
    }
  }


  /*unsigned int num_gg_cand_tag = 0;
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
  hNum_diph_tag_->Fill( num_gg_cand_tag );*/

  if ( fDiphotonNum<1 ) return;

  // retrieve the missing ET
  edm::Handle< edm::View<pat::MET> > mets;
  iEvent.getByToken(metToken_, mets);

  const edm::View<pat::MET>* metColl = mets.product();
  edm::View<pat::MET>::const_iterator met = metColl->begin();
  fMET = met->sumEt();

  tree_->Fill();
}


// ------------ method called once each job just before starting event loop  ------------
void 
TreeProducer::beginJob()
{
  tree_ = new TTree( "evt", "diphoton analyzer tree" );

  tree_->Branch( "num_photon", &fPhotonNum, "num_photon/i" );
  tree_->Branch( "photon_pt", fPhotonPt, "photon_pt[num_photon]/F" );
  tree_->Branch( "photon_eta", fPhotonEta, "photon_eta[num_photon]/F" );
  tree_->Branch( "photon_eta_sc", fPhotonEtaSC, "photon_eta_sc[num_photon]/F" );
  tree_->Branch( "photon_phi", fPhotonPhi, "photon_phi[num_photon]/F" );
  tree_->Branch( "photon_r9", fPhotonR9, "photon_r9[num_photon]/F" );

  tree_->Branch( "num_proton", &fProtonNum, "num_proton/i" );
  tree_->Branch( "proton_xi", fProtonXi, "proton_xi[num_proton]/F" );
  tree_->Branch( "proton_side", fProtonSide, "proton_side[num_proton]/i" );

  tree_->Branch( "num_diphoton", &fDiphotonNum, "num_diphoton/i" );
  tree_->Branch( "diphoton_photon1", fDiphotonPhoton1, "diphoton_photon1[num_diphoton]/i" );
  tree_->Branch( "diphoton_photon2", fDiphotonPhoton2, "diphoton_photon2[num_diphoton]/i" );
  tree_->Branch( "diphoton_mass", fDiphotonM, "diphoton_mass[num_diphoton]/F" );
  tree_->Branch( "diphoton_rapidity", fDiphotonY, "diphoton_rapidity[num_diphoton]/F" );
  tree_->Branch( "diphoton_pt", fDiphotonPt, "diphoton_pt[num_diphoton]/F" );
  tree_->Branch( "diphoton_dphi", fDiphotonDphi, "diphoton_dphi[num_diphoton]/F" );

  tree_->Branch( "num_diproton", &fDiprotonNum, "num_diproton/i" );
  tree_->Branch( "diproton_proton1", fDiprotonProton1, "diproton_proton1[num_diproton]/i" );
  tree_->Branch( "diproton_proton2", fDiprotonProton2, "diproton_proton2[num_diproton]/i" );
  tree_->Branch( "diproton_mass", fDiprotonM, "diproton_mass[num_diproton]/F" );
  tree_->Branch( "diproton_rapidity", fDiprotonY, "diproton_rapidity[num_diproton]/F" );

  tree_->Branch( "met", &fMET );


}

// ------------ method called once each job just after ending the event loop  ------------
void 
TreeProducer::endJob() 
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
TreeProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(TreeProducer);
