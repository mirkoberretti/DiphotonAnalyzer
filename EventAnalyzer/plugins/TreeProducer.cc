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
//         Created:  Tue, 13 Sep 2016 03:57:43 GMT
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

#include "flashgg/DataFormats/interface/Proton.h"
#include "flashgg/DataFormats/interface/DiPhotonCandidate.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/VertexReco/interface/Vertex.h"

#include "DiphotonAnalyzer/EventAnalyzer/interface/SelectionUtils.h"

#include "TFile.h"
#include "TTree.h"
#include "TLorentzVector.h"

//
// class declaration
//

#define MAX_PROTON 10
#define MAX_DIPROTON 5
#define MAX_DIPHOTON 5

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

    void clearTree();

    edm::EDGetTokenT< edm::View<flashgg::DiPhotonCandidate> > diphotonToken_;
    edm::EDGetTokenT< edm::View<flashgg::Proton> > protonToken_;
    edm::EDGetTokenT< edm::View<pat::MET> > metToken_;
    edm::EDGetTokenT< edm::View<reco::Vertex> > vtxToken_;
    double sqrtS_;
    double singlePhotonMinPt_, singlePhotonMaxEta_, singlePhotonMinR9_;
    double photonPairMinMass_;
    std::string filename_;

    TFile* file_;
    TTree* tree_;

    // --- tree components ---

    unsigned int fBX, fRun, fLumiSection;
    unsigned long long fEventNum;

    unsigned int fProtonNum;
    float fProtonXi[MAX_PROTON];
    unsigned int fProtonSide[MAX_PROTON];

    unsigned int fDiprotonNum;
    float fDiprotonM[MAX_DIPROTON], fDiprotonY[MAX_DIPROTON];

    unsigned int fDiphotonNum;
    float fDiphotonPt1[MAX_DIPHOTON], fDiphotonPt2[MAX_DIPHOTON];
    float fDiphotonEta1[MAX_DIPHOTON], fDiphotonEta2[MAX_DIPHOTON];
    float fDiphotonPhi1[MAX_DIPHOTON], fDiphotonPhi2[MAX_DIPHOTON];
    float fDiphotonR91[MAX_DIPHOTON], fDiphotonR92[MAX_DIPHOTON];
    float fDiphotonM[MAX_DIPHOTON], fDiphotonY[MAX_DIPHOTON];
    float fDiphotonPt[MAX_DIPHOTON], fDiphotonDphi[MAX_DIPHOTON];

    float fMET;

    unsigned int fVertexNum;
    unsigned int fDiphotonVertexTracks;
    float fDiphotonVertexX, fDiphotonVertexY, fDiphotonVertexZ;
    float fDiphotonNearestDist;

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
  diphotonToken_( consumes< edm::View<flashgg::DiPhotonCandidate> >( iConfig.getParameter<edm::InputTag>( "diphotonLabel" ) ) ),
  protonToken_  ( consumes< edm::View<flashgg::Proton> >           ( iConfig.getParameter<edm::InputTag>( "protonLabel") ) ),
  metToken_     ( consumes< edm::View<pat::MET> >                  ( iConfig.getParameter<edm::InputTag>( "metLabel") ) ),
  vtxToken_     ( consumes< edm::View<reco::Vertex> >              ( iConfig.getParameter<edm::InputTag>( "vertexLabel" ) ) ),
  sqrtS_             ( iConfig.getParameter<double>( "sqrtS")),
  singlePhotonMinPt_ ( iConfig.getParameter<double>( "minPtSinglePhoton" ) ),
  singlePhotonMaxEta_( iConfig.getParameter<double>( "maxEtaSinglePhoton" ) ),
  singlePhotonMinR9_ ( iConfig.getParameter<double>( "minR9SinglePhoton" ) ),
  photonPairMinMass_ ( iConfig.getParameter<double>( "minMassDiPhoton" ) ),
  filename_          ( iConfig.getParameter<std::string>( "outputFilename" ) ),
  file_( 0 ), tree_( 0 )
{
  //now do what ever initialization is needed
  file_ = new TFile( filename_.c_str(), "recreate" );
  file_->cd();

  tree_ = new TTree( "ntp", "diphoton ntuple" );
}


TreeProducer::~TreeProducer()
{
 
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)

  if ( file_ ) {
    file_->Write();
    file_->Close();
    delete file_;
  }
  if ( tree_ ) delete tree_;

}


//
// member functions
//

void
TreeProducer::clearTree()
{
  fProtonNum = 0;
  for ( unsigned int i=0; i<MAX_PROTON; i++ ) {
    fProtonXi[i] = 0.;
    fProtonSide[i] = 2; //invalid
  }

  fDiprotonNum = 0;
  for ( unsigned int i=0; i<MAX_DIPROTON; i++ ) {
    fDiprotonM[i] = fDiprotonY[i] = 0.;
  }

  fDiphotonNum = 0;
  for ( unsigned int i=0; i<MAX_DIPHOTON; i++ ) {
    fDiphotonPt1[i] = fDiphotonPt2[i] = 0.;
    fDiphotonEta1[i] = fDiphotonEta2[i] = 0.;
    fDiphotonPhi1[i] = fDiphotonPhi2[i] = 0.;
    fDiphotonR91[i] = fDiphotonR92[i] = 0.;
    fDiphotonM[i] = fDiphotonY[i] = fDiphotonPt[i] = fDiphotonDphi[i] = 0.;
  }

  fMET = 0.;

  fVertexNum = 0;
  fDiphotonVertexTracks = 0;
  fDiphotonVertexX = fDiphotonVertexY = fDiphotonVertexZ = 0.;
  fDiphotonNearestDist = 0.;

}

// ------------ method called for each event  ------------
void
TreeProducer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  #include <iostream> // for debugging purposes

  clearTree();

  // Run and BX information
  fBX = iEvent.bunchCrossing();
  fRun = iEvent.id().run();
  fLumiSection = iEvent.luminosityBlock();
  fEventNum = iEvent.id().event();

  // fetch the proton collection from EDM file
  edm::Handle< edm::View<flashgg::Proton> > protons;
  iEvent.getByToken(protonToken_, protons);

  fProtonNum = fDiprotonNum = 0;
  for ( unsigned int i=0; i<protons->size(); i++ ) {
    edm::Ptr<flashgg::Proton> proton = protons->ptrAt( i );

    fProtonXi[i] = proton->xi();
    fProtonSide[i] = proton->side();

    for ( unsigned int j=i+1; j<protons->size(); j++ ) {
      edm::Ptr<flashgg::Proton> proton2 = protons->ptrAt( j );
      if ( proton2->side()==proton->side() ) continue;
      fDiprotonM[fDiprotonNum] = sqrtS_*sqrt( proton->xi()*proton2->xi() );
      fDiprotonY[fDiprotonNum] = log( proton2->xi()/proton->xi() )/2.;

      fDiprotonNum++;
    }

    fProtonNum++;
  }

  // fetch the diphoton collection from EDM file
  edm::Handle< edm::View<flashgg::DiPhotonCandidate> > diphotons;
  iEvent.getByToken(diphotonToken_, diphotons);

  fDiphotonNum = 0;
  edm::Ptr<reco::Vertex> diphoton_vtx;

  for ( unsigned int i=0; i<diphotons->size(); i++ ) {
    edm::Ptr<flashgg::DiPhotonCandidate> diphoton = diphotons->ptrAt( i );

    if ( diphoton->leadPhotonId()<-0.9 ) continue;
    if ( diphoton->subLeadPhotonId()<-0.9 ) continue;

    if ( !passSinglePhotonCuts( diphoton->leadingPhoton() ) ) continue;
    if ( !passSinglePhotonCuts( diphoton->subLeadingPhoton() ) ) continue;

    if ( fabs( diphoton->leadingPhoton()->eta() )>=singlePhotonMaxEta_ or fabs( diphoton->subLeadingPhoton()->eta() )>=singlePhotonMaxEta_ ) continue;
    if ( diphoton->leadingPhoton()->pt()<singlePhotonMinPt_ or diphoton->subLeadingPhoton()->pt()<singlePhotonMinPt_ ) continue;
    if ( diphoton->leadingPhoton()->r9()<singlePhotonMinR9_ or diphoton->subLeadingPhoton()->r9()<singlePhotonMinR9_ ) continue;

    if ( diphoton->mass()<photonPairMinMass_ ) continue;

    diphoton_vtx = diphoton->vtx();
    fDiphotonVertexTracks = diphoton_vtx->tracksSize();
    fDiphotonVertexX = diphoton_vtx->x();
    fDiphotonVertexY = diphoton_vtx->y();
    fDiphotonVertexZ = diphoton_vtx->z();

    fDiphotonPt1[fDiphotonNum] = diphoton->leadingPhoton()->pt();
    fDiphotonEta1[fDiphotonNum] = diphoton->leadingPhoton()->eta();
    fDiphotonPhi1[fDiphotonNum] = diphoton->leadingPhoton()->phi();
    fDiphotonR91[fDiphotonNum] = diphoton->leadingPhoton()->r9();

    fDiphotonPt2[fDiphotonNum] = diphoton->subLeadingPhoton()->pt();
    fDiphotonEta2[fDiphotonNum] = diphoton->subLeadingPhoton()->eta();
    fDiphotonPhi2[fDiphotonNum] = diphoton->subLeadingPhoton()->phi();
    fDiphotonR92[fDiphotonNum] = diphoton->subLeadingPhoton()->r9();

    fDiphotonM[fDiphotonNum] = diphoton->mass();
    fDiphotonY[fDiphotonNum] = diphoton->rapidity();
    fDiphotonPt[fDiphotonNum] = diphoton->pt();

    float dphi = diphoton->leadingPhoton()->phi()-diphoton->subLeadingPhoton()->phi();
    while ( dphi<-TMath::Pi() ) dphi += 2.*TMath::Pi();
    while ( dphi> TMath::Pi() ) dphi -= 2.*TMath::Pi();
    fDiphotonDphi[fDiphotonNum] = dphi;

    std::cout << fDiphotonPt1[fDiphotonNum] << " --- " << fDiphotonPt2[fDiphotonNum] << " --- " << fDiphotonM[fDiphotonNum] << std::endl;

    fDiphotonNum++;
  }

  if ( fDiphotonNum<1 ) return;

std::cout << "---> " << fDiphotonNum << std::endl;
  // retrieve the missing ET
  edm::Handle< edm::View<pat::MET> > mets;
  iEvent.getByToken( metToken_, mets );

  const edm::View<pat::MET>* metColl = mets.product();
  edm::View<pat::MET>::const_iterator met = metColl->begin();
  fMET = met->sumEt();

  edm::Handle< edm::View<reco::Vertex> > vertices;
  iEvent.getByToken( vtxToken_, vertices );
  fVertexNum = vertices->size();

  float mindist_vtx = 999.;
  for ( unsigned int i=0; i<vertices->size(); i++ ) {
    edm::Ptr<reco::Vertex> vtx = vertices->ptrAt( i );
    if ( diphoton_vtx->isValid() and diphoton_vtx->position()==vtx->position() ) continue; // found the diphoton vertex

    const float vtx_dist = sqrt( pow( diphoton_vtx->x()-vtx->x(), 2 )+pow( diphoton_vtx->y()-vtx->y(), 2 )+pow( diphoton_vtx->z()-vtx->z(), 2 ) );
    if ( vtx_dist<mindist_vtx ) mindist_vtx = vtx_dist;
  }
  fDiphotonNearestDist = mindist_vtx;

  tree_->Fill();
}


// ------------ method called once each job just before starting event loop  ------------
void 
TreeProducer::beginJob()
{
  tree_->Branch( "run_id", &fRun, "run_id/i");
  tree_->Branch( "lumisection", &fLumiSection, "lumisection/i");
  tree_->Branch( "bunch_crossing", &fBX, "bunch_crossing/i");
  tree_->Branch( "event_number", &fEventNum, "event_number/l");

  tree_->Branch( "num_proton", &fProtonNum, "num_proton/i" );
  tree_->Branch( "proton_xi", fProtonXi, "proton_xi[num_proton]/F" );
  tree_->Branch( "proton_side", fProtonSide, "proton_side[num_proton]/i" );

  tree_->Branch( "num_diproton", &fDiprotonNum, "num_diproton/i" );
  tree_->Branch( "diproton_mass", fDiprotonM, "diproton_mass[num_diproton]/i" );
  tree_->Branch( "diproton_rapidity", fDiprotonY, "diproton_rapidity[num_diproton]/i" );

  tree_->Branch( "num_diphoton", &fDiphotonNum, "num_diphoton/i" );
  tree_->Branch( "diphoton_pt1", fDiphotonPt1, "diphoton_pt1[num_diphoton]/F" );
  tree_->Branch( "diphoton_pt2", fDiphotonPt2, "diphoton_pt2[num_diphoton]/F" );
  tree_->Branch( "diphoton_eta1", fDiphotonEta1, "diphoton_eta1[num_diphoton]/F" );
  tree_->Branch( "diphoton_eta2", fDiphotonEta2, "diphoton_eta2[num_diphoton]/F" );
  tree_->Branch( "diphoton_phi1", fDiphotonPhi1, "diphoton_phi1[num_diphoton]/F" );
  tree_->Branch( "diphoton_phi2", fDiphotonPhi2, "diphoton_phi2[num_diphoton]/F" );
  tree_->Branch( "diphoton_r91", fDiphotonR91, "diphoton_r91[num_diphoton]/F" );
  tree_->Branch( "diphoton_r92", fDiphotonR92, "diphoton_r92[num_diphoton]/F" );
  tree_->Branch( "diphoton_mass", fDiphotonM, "diphoton_mass[num_diphoton]/F" );
  tree_->Branch( "diphoton_rapidity", fDiphotonY, "diphoton_rapidity[num_diphoton]/F" );
  tree_->Branch( "diphoton_pt", fDiphotonPt, "diphoton_pt[num_diphoton]/F" );
  tree_->Branch( "diphoton_dphi", fDiphotonDphi, "diphoton_dphi[num_diphoton]/F" );

  tree_->Branch( "num_vertex", &fVertexNum, "num_vertex/i" );
  tree_->Branch( "diphoton_vertex_tracks", &fDiphotonVertexTracks, "diphoton_vertex_tracks/i" );
  tree_->Branch( "diphoton_vertex_x", &fDiphotonVertexX, "diphoton_vertex_x/F" );
  tree_->Branch( "diphoton_vertex_y", &fDiphotonVertexY, "diphoton_vertex_y/F" );
  tree_->Branch( "diphoton_vertex_z", &fDiphotonVertexZ, "diphoton_vertex_z/F" );
  tree_->Branch( "diphoton_vertex_nearestvtxdist", &fDiphotonNearestDist, "diphoton_vertex_nearestvtxdist/F" );

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
