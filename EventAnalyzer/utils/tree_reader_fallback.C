#include "Canvas.h"

#define out_path "/afs/cern.ch/user/l/lforthom/www/private/twophoton/"
//#define default_ntp_file "run2016BC_17sep.root"
#define default_ntp_file "run2016BC_28sep.root"

void plot_3hists( const char* name, const char* top_label, TH1* h, TH1* h_1tag, TH1* h_2tag )
{
  Canvas c( name, top_label, true );
  h_1tag->Sumw2();
  h_2tag->Sumw2();
  h->Draw();
  h->SetLineColor( kBlack );
  //h->SetMarkerStyle( 20 );
  h_1tag->Draw("same");
  h_1tag->SetLineColor( kBlack );
  h_1tag->SetMarkerStyle( 20 );
  h_1tag->SetMarkerColor( kRed+1 );
  h_2tag->Draw("same");
  h_2tag->SetLineColor( kBlack );
  h_2tag->SetMarkerStyle( 24 );
  h_2tag->SetMarkerColor( kGreen+2 );
  h_2tag->SetFillStyle( 3005 );
  h_2tag->SetFillColor( kGreen+2 );
  c.AddLegendEntry( h, "All diphotons", "f" );
  c.AddLegendEntry( h_1tag, "#geq 1 proton tag", "ep" );
  c.AddLegendEntry( h_2tag, "#geq 2 proton tags", "ep" );
  c.Prettify( h );
  c.RatioPlot( h, h_1tag, h_2tag, 0., 0.55 );
  //c.Save( "png", out_path );
  c.Save( "pdf", out_path );
  c.Save( "pdf", out_path );
  //c.Save( "png", out_path );
}

void plot_balances( const char* name, const char* top_label, TH2* h2, TH2* h2_mtag=0, TH2* h2_ytag=0, const float& mass_resol=-1., bool abs_unc=false )
{
  Canvas c( name, top_label );
  h2->Draw( "p" );
  c.DrawDiagonal( h2, 0., mass_resol, abs_unc );
  //h2->Draw( "p same" );
  h2->SetMarkerStyle( 24 );
  c.SetLegendY1( 0.18 );
  c.AddLegendEntry( h2, "All candidates", "p" );
  if ( h2_mtag ) {
    h2_mtag->Draw( "p same" );
    h2_mtag->SetMarkerStyle( 20 );
    h2_mtag->SetMarkerSize( .95 );
    h2_mtag->SetMarkerColor( kRed+1 );
    c.AddLegendEntry( h2_mtag, "Mass matching", "p" );
  }
  if ( h2_ytag ) {
    h2_ytag->Draw( "p same" );
    h2_ytag->SetMarkerStyle( 31 );
    h2_ytag->SetMarkerSize( .7 );
    h2_ytag->SetMarkerColor( kGreen+2 );
    c.AddLegendEntry( h2_ytag, "Rapidity matching", "p" );
  }
  c.Prettify( h2 );
  c.Save( "pdf", out_path );
  c.Save( "png", out_path );
}

void tree_reader_fallback( TString file=default_ntp_file )
{
  TFile f(file);
  if ( !f.IsOpen() ) return;

  const float sqrt_s = 13.e3;

  TTree* tr = dynamic_cast<TTree*>( f.Get( "ntp" ) );
  // general quantities
  unsigned int run_id, lumisection;
  unsigned long long event_number;
  tr->SetBranchAddress( "run_id", &run_id );
  tr->SetBranchAddress( "lumisection", &lumisection );
  tr->SetBranchAddress( "event_number", &event_number );
  // diphoton quantities
  unsigned int num_diphoton;
  float diphoton_pt1[5], diphoton_eta1[5], diphoton_phi1[5], diphoton_pt2[5], diphoton_eta2[5], diphoton_phi2[5];
  float diphoton_pt[5], diphoton_mass[5], diphoton_rapidity[5], diphoton_dphi[5];
  tr->SetBranchAddress( "num_diphoton", &num_diphoton );
  tr->SetBranchAddress( "diphoton_pt1", diphoton_pt1 );
  tr->SetBranchAddress( "diphoton_eta1", diphoton_eta1 );
  tr->SetBranchAddress( "diphoton_phi1", diphoton_phi1 );
  tr->SetBranchAddress( "diphoton_pt2", diphoton_pt2 );
  tr->SetBranchAddress( "diphoton_eta2", diphoton_eta2 );
  tr->SetBranchAddress( "diphoton_phi2", diphoton_phi2 );
  tr->SetBranchAddress( "diphoton_pt", diphoton_pt );
  tr->SetBranchAddress( "diphoton_dphi", diphoton_dphi );
  tr->SetBranchAddress( "diphoton_mass", diphoton_mass );
  tr->SetBranchAddress( "diphoton_rapidity", diphoton_rapidity );
  unsigned int diphoton_vertex_vtx1mmdist[5], diphoton_vertex_vtx2mmdist[5], diphoton_vertex_vtx5mmdist[5], diphoton_vertex_vtx1cmdist[5];
  float diphoton_vertex_nearestvtxdist[5];
  unsigned int diphoton_vertex_tracks[5];
  tr->SetBranchAddress( "diphoton_vertex_vtx1mmdist", diphoton_vertex_vtx1mmdist );
  tr->SetBranchAddress( "diphoton_vertex_vtx2mmdist", diphoton_vertex_vtx2mmdist );
  tr->SetBranchAddress( "diphoton_vertex_vtx5mmdist", diphoton_vertex_vtx5mmdist );
  tr->SetBranchAddress( "diphoton_vertex_vtx1cmdist", diphoton_vertex_vtx1cmdist );
  tr->SetBranchAddress( "diphoton_vertex_nearestvtxdist", diphoton_vertex_nearestvtxdist );
  tr->SetBranchAddress( "diphoton_vertex_tracks", diphoton_vertex_tracks );
  // proton quantities
  unsigned int num_proton;
  unsigned int proton_side[5];
  float proton_xi[5];
  tr->SetBranchAddress( "num_proton", &num_proton );
  tr->SetBranchAddress( "proton_side", proton_side );
  tr->SetBranchAddress( "proton_xi", proton_xi );
  // diproton quantities
  unsigned int num_diproton;
  float diproton_mass[5], diproton_rapidity[5];
  tr->SetBranchAddress( "num_diproton", &num_diproton );
  tr->SetBranchAddress( "diproton_mass", diproton_mass );
  tr->SetBranchAddress( "diproton_rapidity", diproton_rapidity );
  // vertex quantities
  unsigned int num_vertex;
  tr->SetBranchAddress( "num_vertex", &num_vertex );
  // other quantities
  float met, met_phi;
  tr->SetBranchAddress( "met", &met );
  tr->SetBranchAddress( "met_phi", &met_phi );

  //                JW
  // Electron Quantities
  unsigned int num_electron;
  float electron_pt[10], electron_eta[10], electron_phi[10], electron_energy[10];
  tr->SetBranchAddress( "num_electron", &num_electron );
  tr->SetBranchAddress( "electron_pt", electron_pt );
  tr->SetBranchAddress( "electron_eta", electron_eta );
  tr->SetBranchAddress( "electron_phi", electron_phi );
  tr->SetBranchAddress( "electron_energy", electron_energy );
  // Muon Quantities
  unsigned int num_muon;
  float muon_pt[10], muon_eta[10], muon_phi[10], muon_energy[10];
  tr->SetBranchAddress("num_muon",&num_muon );
  tr->SetBranchAddress("muon_pt", muon_pt );
  tr->SetBranchAddress("muon_eta",muon_eta );
  tr->SetBranchAddress("muon_phi",muon_phi );
  tr->SetBranchAddress("muon_energy", muon_energy );
  //Jet Quantities
  unsigned int num_jet;
  float jet_pt[10], jet_eta[10], jet_phi[10], jet_energy[10];
  tr->SetBranchAddress("num_jet",&num_jet );
  tr->SetBranchAddress("jet_pt", jet_pt );
  tr->SetBranchAddress("jet_eta",jet_eta );
  tr->SetBranchAddress("jet_phi",jet_phi );
  tr->SetBranchAddress("jet_energy", jet_energy );
  tr->SetBrnachAddress("jet_mass", jet_mass );
  //
  TH1D* h_num_proton = new TH1D( "num_proton", "Number of protons reconstructed in event\\Events", 6, 0., 6. );
  TH1D* h_mpp_over_mgg = new TH1D( "mpp_over_mgg", "m_{pp}^{missing} / m_{#gamma#gamma} for double-tag events\\Events\\?.2f", 30, -2., 4. ),
       *h_ypp_minus_ygg = new TH1D( "ypp_minus_ygg", "y_{pp}^{missing} - y_{#gamma#gamma} for double-tag events\\Events\\?.2f", 50, -2.5, 2.5 );
  TH1D* h_met = new TH1D( "met", "Missing E_{T}\\Events\\GeV?.0f", 48, 0., 240. ),
       *h_met_1tag = (TH1D*)h_met->Clone( "met_1tag" ),
       *h_met_2tag = (TH1D*)h_met->Clone( "met_2tag" );
  TH1D* h_diphoton_pt = new TH1D( "diphoton_pt", "Diphoton p_{T}\\Events\\GeV?.0f", 40, 0., 400. ),
       *h_diphoton_pt_1tag = (TH1D*)h_diphoton_pt->Clone( "diphoton_pt_1tag" ),
       *h_diphoton_pt_2tag = (TH1D*)h_diphoton_pt->Clone( "diphoton_pt_2tag" );
  TH1D* h_diphoton_leadpt = new TH1D( "leadphoton_pt", "Leading photon p_{T}\\Events\\GeV?.0f", 35, 50., 750. ),
       *h_diphoton_leadpt_1tag = (TH1D*)h_diphoton_leadpt->Clone( "leadphoton_pt_1tag" ),
       *h_diphoton_leadpt_2tag = (TH1D*)h_diphoton_leadpt->Clone( "leadphoton_pt_2tag" );
  TH1D* h_diphoton_subleadpt = new TH1D( "subleadphoton_pt", "Subleading photon p_{T}\\Events\\GeV?.0f", 35, 50., 750. ),
       *h_diphoton_subleadpt_1tag = (TH1D*)h_diphoton_subleadpt->Clone( "subleadphoton_pt_1tag" ),
       *h_diphoton_subleadpt_2tag = (TH1D*)h_diphoton_subleadpt->Clone( "subleadphoton_pt_2tag" );
  TH1D* h_diphoton_leadeta = new TH1D( "leadphoton_eta", "Leading photon #eta\\Events\\?.3f", 40, -2.5, 2.5 ),
       *h_diphoton_leadeta_1tag = (TH1D*)h_diphoton_leadeta->Clone( "leadphoton_eta_1tag" ),
       *h_diphoton_leadeta_2tag = (TH1D*)h_diphoton_leadeta->Clone( "leadphoton_eta_2tag" );
  TH1D* h_diphoton_subleadeta = new TH1D( "subleadphoton_eta", "Subleading photon #eta\\Events\\?.3f", 40, -2.5, 2.5 ),
       *h_diphoton_subleadeta_1tag = (TH1D*)h_diphoton_subleadeta->Clone( "subleadphoton_eta_1tag" ),
       *h_diphoton_subleadeta_2tag = (TH1D*)h_diphoton_subleadeta->Clone( "subleadphoton_eta_2tag" );
  TH1D* h_diphoton_dphi = new TH1D( "diphoton_dphi", "Diphoton 1-|#Delta#phi/#pi|\\Events\\?.2f", 50, 0., 1. ),
       *h_diphoton_dphi_1tag = (TH1D*)h_diphoton_dphi->Clone( "diphoton_dphi_1tag" ),
       *h_diphoton_dphi_2tag = (TH1D*)h_diphoton_dphi->Clone( "diphoton_dphi_2tag" );
  TH1D* h_diphoton_mass = new TH1D( "diphoton_mass", "Diphoton mass\\Events\\GeV?.1f", 50, 500., 2000. ),
       *h_diphoton_mass_1tag = (TH1D*)h_diphoton_mass->Clone( "diphoton_mass_1tag" ),
       *h_diphoton_mass_2tag = (TH1D*)h_diphoton_mass->Clone( "diphoton_mass_2tag" );
  TH1D* h_diphoton_rap = new TH1D( "diphoton_rap", "Diphoton rapidity\\Events\\?.1f", 40, -2.5, 2.5 ),
       *h_diphoton_rap_1tag = (TH1D*)h_diphoton_rap->Clone( "diphoton_rap_1tag" ),
       *h_diphoton_rap_2tag = (TH1D*)h_diphoton_rap->Clone( "diphoton_rap_2tag" );
  TH1D* h_diphoton_ntrk = new TH1D( "diphoton_ntrk", "Number of tracks on diphoton vertex\\Events\\?.1f", 20, 0., 20. ),
       *h_diphoton_ntrk_1tag = (TH1D*)h_diphoton_ntrk->Clone( "diphoton_ntrk_1tag" ),
       *h_diphoton_ntrk_2tag = (TH1D*)h_diphoton_ntrk->Clone( "diphoton_ntrk_2tag" );
  TH1D* h_num_vtx = new TH1D( "num_vtx", "Number of primary vertices in event\\Events", 40, 0., 40. ),
       *h_num_vtx_1tag = (TH1D*)h_num_vtx->Clone( "num_vtx_1tag" ),
       *h_num_vtx_2tag = (TH1D*)h_num_vtx->Clone( "num_vtx_2tag" );
  TH1D* h_num_vtx_1mm = new TH1D( "num_vtx_1mm", "Number of primary vertices near diphoton vertex\\Events", 12, 0., 12. ),
       *h_num_vtx_2mm = (TH1D*)h_num_vtx_1mm->Clone( "num_vtx_2mm" ),
       *h_num_vtx_5mm = (TH1D*)h_num_vtx_1mm->Clone( "num_vtx_5mm" ),
       *h_num_vtx_1cm = (TH1D*)h_num_vtx_1mm->Clone( "num_vtx_1cm" );
  TH1D* h_diphoton_closestvtx = new TH1D( "diphoton_closestvtx", "Distance diphoton/nearest vertex\\Events\\mm?.1f", 25, 0., 2.5 ),
       *h_diphoton_closestvtx_1tag = (TH1D*)h_diphoton_closestvtx->Clone( "diphoton_closestvtx_1tag" ),
       *h_diphoton_closestvtx_2tag = (TH1D*)h_diphoton_closestvtx->Clone( "diphoton_closestvtx_2tag" );
  TH2D* h_met_vs_pt = new TH2D( "met_vs_pt", "Missing E_{T} (GeV)\\Diphoton p_{T} (GeV)", 40, 0., 400., 40, 0., 400. ),
       *h_met_vs_pt_2tag = (TH2D*)h_met_vs_pt->Clone( "met_vs_pt_2tag" ),
       *h_metx_vs_mety = new TH2D( "metx_vs_mety", "#slash{E}_{T,x} (GeV)\\#slash{E}_{T,y} (GeV)", 50, -125., 125., 50, -125., 125. ),
       *h_metx_vs_mety_2tag = (TH2D*)h_metx_vs_mety->Clone( "metx_vs_mety_2tag" );
  TH2D* h_ygg_vs_ypp = new TH2D( "ygg_vs_ypp", "Diphoton rapidity\\Diproton rapidity", 600, -3., 3., 600, -3., 3. ),
       *h_ygg_vs_ypp_candm = (TH2D*)h_ygg_vs_ypp->Clone("ygg_vs_ypp_candm"),
       *h_ygg_vs_ypp_candy = (TH2D*)h_ygg_vs_ypp->Clone("ygg_vs_ypp_candy");
  TH2D* h_yggmet_vs_ypp = new TH2D( "yggmet_vs_ypp", "Diphoton + #slash{E}_{T} rapidity\\Diproton rapidity", 600, -3., 3., 600, -3., 3. ),
       *h_yggmet_vs_ypp_candm = (TH2D*)h_yggmet_vs_ypp->Clone("yggmet_vs_ypp_candm"),
       *h_yggmet_vs_ypp_candy = (TH2D*)h_yggmet_vs_ypp->Clone("yggmet_vs_ypp_candy");
  TH2D* h_mgg_vs_mpp = new TH2D( "mgg_vs_mpp", "Diphoton mass (GeV)\\Diproton missing mass (GeV)", 1750, 250., 2000., 1750, 250., 2000. ),
       *h_mgg_vs_mpp_candm = (TH2D*)h_mgg_vs_mpp->Clone("mgg_vs_mpp_candm"),
       *h_mgg_vs_mpp_candy = (TH2D*)h_mgg_vs_mpp->Clone("mgg_vs_mpp_candy");
  TH2D* h_mggmet_vs_mpp = new TH2D( "mggmet_vs_mpp", "Diphoton + #slash{E}_{T} mass (GeV)\\Diproton missing mass (GeV)", 1750, 250., 2000., 1750, 250., 2000. ),
       *h_mggmet_vs_mpp_candm = (TH2D*)h_mggmet_vs_mpp->Clone("mggmet_vs_mpp_candm"),
       *h_mggmet_vs_mpp_candy = (TH2D*)h_mggmet_vs_mpp->Clone("mggmet_vs_mpp_candy");
  TH2D* h_xi1gg_vs_xi1pp = new TH2D( "xi1gg_vs_xi1pp", "#xi_{1} from diphoton system\\Proton #xi_{1}", 1000, 0., 0.5, 1000, 0., 0.5 ),
       *h_xi1gg_vs_xi1pp_candm = (TH2D*)h_xi1gg_vs_xi1pp->Clone( "xi1gg_vs_xi1pp_candm" ),
       *h_xi1gg_vs_xi1pp_candy = (TH2D*)h_xi1gg_vs_xi1pp->Clone( "xi1gg_vs_xi1pp_candy" ),
       *h_xi1ggmet_vs_xi1pp = (TH2D*)h_xi1gg_vs_xi1pp->Clone( "xi1ggmet_vs_xi1pp" ),
       *h_xi2gg_vs_xi2pp = new TH2D( "xi2gg_vs_xi2pp", "#xi_{2} from diphoton system\\Proton #xi_{2}", 1000, 0., 0.5, 1000, 0., 0.5 ),
       *h_xi2gg_vs_xi2pp_candm = (TH2D*)h_xi2gg_vs_xi2pp->Clone( "xi2gg_vs_xi2pp_candm" ),
       *h_xi2gg_vs_xi2pp_candy = (TH2D*)h_xi2gg_vs_xi2pp->Clone( "xi2gg_vs_xi2pp_candy" ),
       *h_xi2ggmet_vs_xi2pp = (TH2D*)h_xi2gg_vs_xi2pp->Clone( "xi2ggmet_vs_xi2pp" );

  ofstream events_list( "events_list_2016BC.txt" );

  const double rel_err_xi = 0.15; // 15% error on xi determination

  unsigned int num_evts_notag = 0, num_evts_with_tag = 0;
  TLorentzVector pho1, pho2, electron, muon, jet;
  // tree readout stage
  for ( unsigned int i=0; i<tr->GetEntries(); i++ ) {
    tr->GetEntry( i );
    // dump the list of events in a text file

    h_num_proton->Fill( num_proton );
    unsigned int num_1tag = 0, num_2tag = 0;
    float xi_prot1 = -1., xi_prot2 = -1.;
    for ( unsigned int j=0; j<num_proton; j++ ) {
      for ( unsigned int k=j+1; k<num_proton; k++ ) {
         if ( proton_side[j]==proton_side[k] ) continue;
         if ( proton_side[j]==0 ) {
           xi_prot1 = proton_xi[j];
           xi_prot2 = proton_xi[k];
         }
         else {
           xi_prot1 = proton_xi[k];
           xi_prot2 = proton_xi[j];
         }
         num_2tag++;
      }
      num_1tag++;
    }

    float max_diproton_mass = -1., max_diproton_mass_rap = -999.;
    for ( unsigned int j=0; j<num_diproton; j++ ) {
      if ( diproton_mass[j]>max_diproton_mass ) {
        max_diproton_mass = diproton_mass[j];
        max_diproton_mass_rap = diproton_rapidity[j];
      }
    }


    //          JW
    for ( unsigned int j=0; j<num_electron; j++ ) {
      electron.SetPtEtaPhiM( electron_pt[j], electron_eta[j], electron_phi[j], 0.000510998928 );
    }
    for ( unsigned int j=0; j<num_muon; j++ ) {
      muon.SetPtEtaPhiM( muon_pt[j], muon_eta[j], muon_phi[j], .1056583715 );
    }
    for ( unsigned int j=0; j<num_jet; j++ ) {
      jet.SetPtEtaPhiM( jet_pt[j], jet_eta[j], jet_phi[j], jet_mass[j] );
    }
    //

    for ( unsigned int j=0; j<num_diphoton; j++ ) {
      const float xi_reco1 = ( diphoton_pt1[j] * exp( -diphoton_eta1[j] ) + diphoton_pt2[j] * exp( -diphoton_eta2[j] ) )/sqrt_s,
                  xi_reco2 = ( diphoton_pt1[j] * exp(  diphoton_eta1[j] ) + diphoton_pt2[j] * exp(  diphoton_eta2[j] ) )/sqrt_s;
      const float xi_reco1_withmet = xi_reco1 + met/sqrt_s,
                  xi_reco2_withmet = xi_reco2 + met/sqrt_s;

      //cout << xi_reco1 << ", " << xi_reco2 << endl;
      h_diphoton_pt->Fill( diphoton_pt[j] );
      h_diphoton_mass->Fill( diphoton_mass[j] );
      h_diphoton_rap->Fill( diphoton_rapidity[j] );
      h_diphoton_closestvtx->Fill( diphoton_vertex_nearestvtxdist[j] );
      h_diphoton_ntrk->Fill( diphoton_vertex_tracks[j] );
      h_diphoton_dphi->Fill( 1-fabs( diphoton_dphi[j]/TMath::Pi() ) );
      h_diphoton_leadpt->Fill( diphoton_pt1[j] );
      h_diphoton_subleadpt->Fill( diphoton_pt2[j] );
      h_diphoton_leadeta->Fill( diphoton_eta1[j] );
      h_diphoton_subleadeta->Fill( diphoton_eta2[j] );

      const float met_x = met*cos( met_phi ),
                  met_y = met*sin( met_phi );
      pho1.SetPtEtaPhiM( diphoton_pt1[j], diphoton_eta1[j], diphoton_phi1[j], 0. );
      pho2.SetPtEtaPhiM( diphoton_pt2[j], diphoton_eta2[j], diphoton_phi2[j], 0. );
 
      //                  JW
      const TLorentzVector lv_met( met_x, met_y, 0., met ),
                           dipho_met = pho1+pho2+electron+muon+jet+lv_met;
      //cout << dipho_met.M() << " <---> " << diphoton_mass[j] << endl;
      const float diphoton_plus_met_mass = dipho_met.M(),
                  diphoton_plus_met_rap = dipho_met.Rapidity();

      h_met_vs_pt->Fill( met, diphoton_pt[j] );
      h_metx_vs_mety->Fill( met_x, met_y );

      if ( num_1tag>0 ) {
        h_diphoton_pt_1tag->Fill( diphoton_pt[j] );
        h_diphoton_mass_1tag->Fill( diphoton_mass[j] );
        h_diphoton_rap_1tag->Fill( diphoton_rapidity[j] );
        h_diphoton_closestvtx_1tag->Fill( diphoton_vertex_nearestvtxdist[j] );
        h_diphoton_ntrk_1tag->Fill( diphoton_vertex_tracks[j] );
        h_diphoton_dphi_1tag->Fill( 1-fabs( diphoton_dphi[j]/TMath::Pi() ) );
        h_diphoton_leadpt_1tag->Fill( diphoton_pt1[j] );
        h_diphoton_subleadpt_1tag->Fill( diphoton_pt2[j] );
        h_diphoton_leadeta_1tag->Fill( diphoton_eta1[j] );
        h_diphoton_subleadeta_1tag->Fill( diphoton_eta2[j] );
      }
      if ( num_2tag>0 ) {
        h_diphoton_pt_2tag->Fill( diphoton_pt[j] );
        h_diphoton_mass_2tag->Fill( diphoton_mass[j] );
        h_diphoton_rap_2tag->Fill( diphoton_rapidity[j] );
        h_diphoton_closestvtx_2tag->Fill( diphoton_vertex_nearestvtxdist[j] );
        h_diphoton_ntrk_2tag->Fill( diphoton_vertex_tracks[j] );
        h_diphoton_dphi_2tag->Fill( 1-fabs( diphoton_dphi[j]/TMath::Pi() ) );
        h_diphoton_leadpt_2tag->Fill( diphoton_pt1[j] );
        h_diphoton_subleadpt_2tag->Fill( diphoton_pt2[j] );
        h_diphoton_leadeta_2tag->Fill( diphoton_eta1[j] );
        h_diphoton_subleadeta_2tag->Fill( diphoton_eta2[j] );

        h_xi1gg_vs_xi1pp->Fill( xi_reco1, xi_prot1 );
        h_xi2gg_vs_xi2pp->Fill( xi_reco2, xi_prot2 );

        h_xi1ggmet_vs_xi1pp->Fill( xi_reco1_withmet, xi_prot1 );
        h_xi2ggmet_vs_xi2pp->Fill( xi_reco2_withmet, xi_prot2 );

        h_met_vs_pt_2tag->Fill( met, diphoton_pt[j] );
        h_metx_vs_mety_2tag->Fill( met_x, met_y );

        events_list << run_id << ":" << lumisection << ":" << event_number << endl;
      }

      if ( num_diproton>0 ) {
        h_mgg_vs_mpp->Fill( diphoton_mass[j], max_diproton_mass );
        h_mggmet_vs_mpp->Fill( diphoton_plus_met_mass, max_diproton_mass );
        h_ygg_vs_ypp->Fill( diphoton_rapidity[j], max_diproton_mass_rap );
	h_yggmet_vs_ypp->Fill( diphoton_plus_met_rap, max_diproton_mass_rap );
        h_mpp_over_mgg->Fill( max_diproton_mass/diphoton_mass[j] );
        h_ypp_minus_ygg->Fill( max_diproton_mass_rap - diphoton_rapidity[j] );
        //if ( fabs( diphoton_mass[j]-max_diproton_mass )<max_diproton_mass*rel_err_xi ) {
        if ( fabs( diphoton_plus_met_mass-max_diproton_mass )<max_diproton_mass*rel_err_xi ) {
          h_mgg_vs_mpp_candm->Fill( diphoton_mass[j], max_diproton_mass );
          h_mggmet_vs_mpp_candm->Fill( diphoton_plus_met_mass, max_diproton_mass );
          h_ygg_vs_ypp_candm->Fill( diphoton_rapidity[j], max_diproton_mass_rap );
	  h_yggmet_vs_ypp_candm->Fill( diphoton_plus_met_rap, max_diproton_mass_rap );
          h_xi1gg_vs_xi1pp_candm->Fill( xi_reco1, xi_prot1 );
          h_xi2gg_vs_xi2pp_candm->Fill( xi_reco2, xi_prot2 );

          //if ( fabs( diphoton_rapidity[j]-max_diproton_mass_rap )<rel_err_xi/sqrt( 2. ) ) {
          if ( fabs( diphoton_plus_met_rap-max_diproton_mass_rap )<rel_err_xi/sqrt( 2. ) ) {
            /*TLorentzVector p1, p2;
            p1.SetPtEtaPhiM( diphoton_pt1[j], diphoton_eta1[j], diphoton_phi1[j], 0. );
            p2.SetPtEtaPhiM( diphoton_pt2[j], diphoton_eta2[j], diphoton_phi2[j], 0. );
            cout << "---------> " << (p1+p2).M() << endl;*/

            cout << " ---> event: " << run_id << ":" << lumisection << ":" << event_number << endl;
            cout << "      diphoton: pt=" << diphoton_pt[j] << ", mass=" << diphoton_mass[j] << ", rapidity=" << diphoton_rapidity[j] << ", dphi=" << diphoton_dphi[j] << endl
                 << "      diproton: mass=" << max_diproton_mass << ", rapidity=" << max_diproton_mass_rap << endl
                 << "      MET: " << met << endl;
            cout << "      single photon: pt=" << diphoton_pt1[j] << ", eta=" << diphoton_eta1[j] << ", phi=" << diphoton_phi1[j] << endl
                 << "                     pt=" << diphoton_pt2[j] << ", eta=" << diphoton_eta2[j] << ", phi=" << diphoton_phi2[j] << endl;
          }
        }
        //if ( fabs( diphoton_rapidity[j]-max_diproton_mass_rap )<rel_err_xi/sqrt( 2. ) ) {
        if ( fabs( diphoton_plus_met_rap-max_diproton_mass_rap )<rel_err_xi/sqrt( 2. ) ) {
          h_mgg_vs_mpp_candy->Fill( diphoton_mass[j], max_diproton_mass );
          h_mggmet_vs_mpp_candy->Fill( diphoton_plus_met_mass, max_diproton_mass );
          h_ygg_vs_ypp_candy->Fill( diphoton_rapidity[j], max_diproton_mass_rap );
	  h_yggmet_vs_ypp_candy->Fill( diphoton_plus_met_rap, max_diproton_mass_rap );
          h_xi1gg_vs_xi1pp_candy->Fill( xi_reco1, xi_prot1 );
          h_xi2gg_vs_xi2pp_candy->Fill( xi_reco2, xi_prot2 );
        }
        num_evts_with_tag++;
      }

      h_num_vtx_1mm->Fill( diphoton_vertex_vtx1mmdist[j] );
      h_num_vtx_2mm->Fill( diphoton_vertex_vtx2mmdist[j] );
      h_num_vtx_5mm->Fill( diphoton_vertex_vtx5mmdist[j] );
      h_num_vtx_1cm->Fill( diphoton_vertex_vtx1cmdist[j] );

      num_evts_notag++;
    }
    h_num_vtx->Fill( num_vertex );
    h_met->Fill( met );
    if ( num_1tag>0 ) {
      h_num_vtx_1tag->Fill( num_vertex );
      h_met_1tag->Fill( met );
    }
    if ( num_2tag>0 ) {
      h_num_vtx_2tag->Fill( num_vertex );
      h_met_2tag->Fill( met );
    }

  }

  // plotting stage
  cout << "events: " << num_evts_notag << ", with tag: " << num_evts_with_tag << endl;

  const float lumi_b = 5.060924481910, // fb-1
              lumi_c = 1.490748474431, // fb-1
              lumi_g = 3.742171002882; // fb-1

  const float lumi = lumi_b+lumi_c;
  /*float lumi = 0.; string run_name;
  switch ( run ) {
    case 'B': lumi = lumi_b; run_name = "B"; break;
    case 'C': lumi = lumi_c; run_name = "C"; break;
    case '0': default: lumi = lumi_b+lumi_c; run_name = "BC"; break;
  }*/

  gStyle->SetOptStat( 0 );

  const string top_label_str = Form( "CMS+CTPPS Preliminary 2016, #sqrt{s} = 13 TeV, L = %.2f fb^{-1}", lumi );
  const char* top_label = top_label_str.c_str();

  {
    plot_3hists( "diphoton_mass", top_label, h_diphoton_mass, h_diphoton_mass_1tag, h_diphoton_mass_2tag );
    plot_3hists( "diphoton_pt", top_label, h_diphoton_pt, h_diphoton_pt_1tag, h_diphoton_pt_2tag );
    plot_3hists( "diphoton_lead_pt", top_label, h_diphoton_leadpt, h_diphoton_leadpt_1tag, h_diphoton_leadpt_2tag );
    plot_3hists( "diphoton_sublead_pt", top_label, h_diphoton_subleadpt, h_diphoton_subleadpt_1tag, h_diphoton_subleadpt_2tag );

    /*h_diphoton_leadeta->GetYaxis()->SetRangeUser(0., 55.);
    h_diphoton_subleadeta->GetYaxis()->SetRangeUser(0., 55.);*/
    h_diphoton_leadeta->SetMinimum( 0. );
    h_diphoton_subleadeta->SetMinimum( 0. );

    plot_3hists( "diphoton_lead_eta", top_label, h_diphoton_leadeta, h_diphoton_leadeta_1tag, h_diphoton_leadeta_2tag );
    plot_3hists( "diphoton_sublead_eta", top_label, h_diphoton_subleadeta, h_diphoton_subleadeta_1tag, h_diphoton_subleadeta_2tag );
    plot_3hists( "diphoton_dphi", top_label, h_diphoton_dphi, h_diphoton_dphi_1tag, h_diphoton_dphi_2tag );
    plot_3hists( "diphoton_rapidity", top_label, h_diphoton_rap, h_diphoton_rap_1tag, h_diphoton_rap_2tag );
    plot_3hists( "diphoton_closest_vtx", top_label, h_diphoton_closestvtx, h_diphoton_closestvtx_1tag, h_diphoton_closestvtx_2tag );
    plot_3hists( "diphoton_vtx_numtracks", top_label, h_diphoton_ntrk, h_diphoton_ntrk_1tag, h_diphoton_ntrk_2tag );
    plot_3hists( "num_vertex", top_label, h_num_vtx, h_num_vtx_1tag, h_num_vtx_2tag );
    plot_3hists( "event_met", top_label, h_met, h_met_1tag, h_met_2tag );

    cout << "total candidates: " << h_mgg_vs_mpp->Integral() << endl
         << " -> with mass matching: " << h_mgg_vs_mpp_candm->Integral() << endl
         << " -> with rapiditiy matching: " << h_mgg_vs_mpp_candy->Integral() << endl;

    plot_balances( "mass_balance", top_label, h_mgg_vs_mpp, h_mgg_vs_mpp_candm, h_mgg_vs_mpp_candy, rel_err_xi );
    plot_balances( "mass_balance_withmet", top_label, h_mggmet_vs_mpp, h_mggmet_vs_mpp_candm, h_mggmet_vs_mpp_candy, rel_err_xi );
    plot_balances( "rapidity_balance", top_label, h_ygg_vs_ypp, h_ygg_vs_ypp_candm, h_ygg_vs_ypp_candy, rel_err_xi/sqrt( 2. ), true );
    plot_balances( "rapidity_balance_withmet", top_label, h_yggmet_vs_ypp, h_yggmet_vs_ypp_candm, h_yggmet_vs_ypp_candy, rel_err_xi/sqrt( 2. ), true );
    plot_balances( "xi1_balance", top_label, h_xi1gg_vs_xi1pp, h_xi1gg_vs_xi1pp_candm, h_xi1gg_vs_xi1pp_candy, rel_err_xi );
    plot_balances( "xi2_balance", top_label, h_xi2gg_vs_xi2pp, h_xi2gg_vs_xi2pp_candm, h_xi2gg_vs_xi2pp_candy, rel_err_xi );
    plot_balances( "xi1_balance_withmet", top_label, h_xi1ggmet_vs_xi1pp, 0, 0, rel_err_xi );
    plot_balances( "xi2_balance_withmet", top_label, h_xi2ggmet_vs_xi2pp, 0, 0, rel_err_xi );

  }

  {
    Canvas c( "met_x_vs_y", top_label );
    h_metx_vs_mety->Draw( "colz" );
    h_metx_vs_mety_2tag->Draw( "p same" );
    h_metx_vs_mety_2tag->SetMarkerStyle( 24 );
    h_metx_vs_mety_2tag->SetMarkerColor( kRed );
    c.SetLegendX1( 0.15 );
    c.SetLegendY1( 0.15 );
    c.AddLegendEntry( h_metx_vs_mety_2tag, "#geq 2 proton tags", "p" );
    c.Prettify( h_metx_vs_mety );
    c.Save( "pdf", out_path );
    c.Save( "png", out_path );
  }
  {
    Canvas c( "diphoton_pt_vs_met", top_label );
    h_met_vs_pt->Draw( "colz" );
    h_met_vs_pt_2tag->Draw( "p same" );
    h_met_vs_pt_2tag->SetMarkerStyle( 24 );
    h_met_vs_pt_2tag->SetMarkerColor( kRed );
    c.AddLegendEntry( h_met_vs_pt_2tag, "#geq 2 proton tags", "p" );
    c.Prettify( h_met_vs_pt );
    c.Save( "pdf", out_path );
    c.Save( "png", out_path );
  }
  {
    Canvas c( "num_close_vertex", top_label );
    h_num_vtx_1mm->Sumw2();
    h_num_vtx_2mm->Sumw2();
    h_num_vtx_5mm->Sumw2();
    h_num_vtx_1cm->Sumw2();
    h_num_vtx_1mm->Draw();
    h_num_vtx_2mm->Draw( "same" );
    h_num_vtx_2mm->SetMarkerColor( kRed+1 );
    h_num_vtx_5mm->Draw( "same" );
    h_num_vtx_5mm->SetMarkerColor( kGreen+2 );
    h_num_vtx_1cm->Draw( "same" );
    h_num_vtx_1cm->SetMarkerColor( kBlue+1 );
    h_num_vtx_1mm->SetMarkerStyle( 20 );
    h_num_vtx_2mm->SetMarkerStyle( 21 );
    h_num_vtx_5mm->SetMarkerStyle( 22 );
    h_num_vtx_1cm->SetMarkerStyle( 23 );
    c.AddLegendEntry( h_num_vtx_1mm, "at 1 mm distance" );
    c.AddLegendEntry( h_num_vtx_2mm, "at 2 mm distance" );
    c.AddLegendEntry( h_num_vtx_5mm, "at 5 mm distance" );
    c.AddLegendEntry( h_num_vtx_1cm, "at 1 cm distance" );
    c.Prettify( h_num_vtx_1mm );
    c.SetLogy();
    c.Save( "pdf", out_path );
    c.Save( "png", out_path );
  }
  {
    Canvas c( "num_proton", top_label );
    h_num_proton->Sumw2();
    h_num_proton->Draw();
    h_num_proton->SetMarkerStyle( 20 );
    c.Prettify( h_num_proton );
    c.Save( "pdf", out_path );
    c.Save( "png", out_path );
  }
  {
    Canvas c( "mass_ratio", top_label );
    h_mpp_over_mgg->Sumw2();
    h_mpp_over_mgg->Draw();
    h_mpp_over_mgg->SetMarkerStyle( 20 );
    c.Prettify( h_mpp_over_mgg );
    c.Save( "pdf", out_path );
    c.Save( "png", out_path );
  }
  {
    Canvas c( "rapidity_difference", top_label );
    h_ypp_minus_ygg->Sumw2();
    h_ypp_minus_ygg->Draw();
    h_ypp_minus_ygg->SetMarkerStyle( 20 );
    c.Prettify( h_ypp_minus_ygg );
    c.Save( "pdf", out_path );
    c.Save( "png", out_path );
  }

}

//  LocalWords:  SetMarkerStyle
