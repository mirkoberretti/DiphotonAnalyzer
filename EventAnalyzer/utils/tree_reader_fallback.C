#include "Canvas.h"

void plot_3hists( const char* name, const char* top_label, TH1* h, TH1* h_1tag, TH1* h_2tag )
{
  Canvas c( name, top_label, true );
  h_1tag->Sumw2();
  h_2tag->Sumw2();
  h->Draw();
  //h->SetMarkerStyle( 20 );
  h_1tag->Draw("same");
  h_1tag->SetMarkerStyle( 20 );
  h_1tag->SetMarkerColor( kRed+1 );
  h_2tag->Draw("same");
  h_2tag->SetMarkerStyle( 24 );
  h_2tag->SetMarkerColor( kGreen+2 );
  h_2tag->SetFillStyle( 3005 );
  h_2tag->SetFillColor( kGreen+2 );
  c.AddLegendEntry( h, "All diphotons", "f" );
  c.AddLegendEntry( h_1tag, "#geq 1 proton tag", "ep" );
  c.AddLegendEntry( h_2tag, "#geq 2 proton tags", "ep" );
  c.Prettify( h );
  c.RatioPlot( h, h_1tag, h_2tag, 0., 0.55 );
  c.Save( "png", "/afs/cern.ch/user/l/lforthom/www/private/twophoton/" );
  c.Save( "pdf", "/afs/cern.ch/user/l/lforthom/www/private/twophoton/" );
}

void tree_reader_fallback( TString file="run2016BB_17sep.root" )
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
  float met;
  tr->SetBranchAddress( "met", &met );

  TH1D* h_num_proton = new TH1D( "num_proton", "Number of protons reconstructed in event\\Events", 6, 0., 6. );
  TH1D* h_met = new TH1D( "met", "Missing E_{T}\\Events\\GeV?.0f", 48, 0., 240. ),
       *h_met_1tag = (TH1D*)h_met->Clone( "met_1tag" ),
       *h_met_2tag = (TH1D*)h_met->Clone( "met_2tag" );
  TH1D* h_diphoton_pt = new TH1D( "diphoton_pt", "Diphoton p_{T}\\Events\\GeV?.0f", 40, 0., 400. ),
       *h_diphoton_pt_1tag = (TH1D*)h_diphoton_pt->Clone( "diphoton_pt_1tag" ),
       *h_diphoton_pt_2tag = (TH1D*)h_diphoton_pt->Clone( "diphoton_pt_2tag" );
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
  TH2D* h_met_vs_pt = new TH2D( "met_vs_pt", "Missing E_{T} (GeV)\\Diphoton p_{T} (GeV)", 40, 0., 400., 40, 0., 400. );
  TH2D* h_ygg_vs_ypp = new TH2D( "ygg_vs_ypp", "Diphoton rapidity\\Diproton rapidity", 600, -3., 3., 600, -3., 3. ),
       *h_ygg_vs_ypp_cand = (TH2D*)h_ygg_vs_ypp->Clone("ygg_vs_ypp_cand");
  TH2D* h_mgg_vs_mpp = new TH2D( "mgg_vs_mpp", "Diphoton mass (GeV)\\Diproton missing mass (GeV)", 1750, 250., 2000., 1750, 250., 2000. ),
       *h_mgg_vs_mpp_cand = (TH2D*)h_mgg_vs_mpp->Clone("mgg_vs_mpp_cand");
  TH2D* h_xi1gg_vs_xi1pp = new TH2D( "xi1gg_vs_xi1pp", "#xi_{1} from diphoton system\\Proton #xi_{1}", 500, 0., 0.5, 500, 0., 0.5 ),
       *h_xi1gg_vs_xi1pp_cand = (TH2D*)h_xi1gg_vs_xi1pp->Clone( "xi1gg_vs_xi1pp_cand" ),
       *h_xi2gg_vs_xi2pp = new TH2D( "xi2gg_vs_xi2pp", "#xi_{2} from diphoton system\\Proton #xi_{2}", 500, 0., 0.5, 500, 0., 0.5 ),
       *h_xi2gg_vs_xi2pp_cand = (TH2D*)h_xi2gg_vs_xi2pp->Clone( "xi2gg_vs_xi2pp_cand" );

  ofstream events_list( "events_list_2016BC.txt" );

  // tree readout stage
  for ( unsigned int i=0; i<tr->GetEntries(); i++ ) {
    tr->GetEntry( i );
    // dump the list of events in a text file
    events_list << run_id << ":" << lumisection << ":" << event_number << endl;

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

    for ( unsigned int j=0; j<num_diphoton; j++ ) {
      const float xi_reco1 = ( diphoton_pt1[j] * exp( -diphoton_eta1[j] ) + diphoton_pt2[j] * exp( -diphoton_eta2[j] ) )/sqrt_s,
                  xi_reco2 = ( diphoton_pt1[j] * exp(  diphoton_eta1[j] ) + diphoton_pt2[j] * exp(  diphoton_eta2[j] ) )/sqrt_s;
      //cout << xi_reco1 << ", " << xi_reco2 << endl;
      h_diphoton_pt->Fill( diphoton_pt[j] );
      h_diphoton_mass->Fill( diphoton_mass[j] );
      h_diphoton_rap->Fill( diphoton_rapidity[j] );
      h_diphoton_closestvtx->Fill( diphoton_vertex_nearestvtxdist[j] );
      h_diphoton_ntrk->Fill( diphoton_vertex_tracks[j] );
      h_diphoton_dphi->Fill( 1-fabs( diphoton_dphi[j]/TMath::Pi() ) );
      if ( num_1tag>0 ) {
        h_diphoton_pt_1tag->Fill( diphoton_pt[j] );
        h_diphoton_mass_1tag->Fill( diphoton_mass[j] );
        h_diphoton_rap_1tag->Fill( diphoton_rapidity[j] );
        h_diphoton_closestvtx_1tag->Fill( diphoton_vertex_nearestvtxdist[j] );
        h_diphoton_ntrk_1tag->Fill( diphoton_vertex_tracks[j] );
        h_diphoton_dphi_1tag->Fill( 1-fabs( diphoton_dphi[j]/TMath::Pi() ) );
      }
      if ( num_2tag>0 ) {
        h_diphoton_pt_2tag->Fill( diphoton_pt[j] );
        h_diphoton_mass_2tag->Fill( diphoton_mass[j] );
        h_diphoton_rap_2tag->Fill( diphoton_rapidity[j] );
        h_diphoton_closestvtx_2tag->Fill( diphoton_vertex_nearestvtxdist[j] );
        h_diphoton_ntrk_2tag->Fill( diphoton_vertex_tracks[j] );
        h_diphoton_dphi_2tag->Fill( 1-fabs( diphoton_dphi[j]/TMath::Pi() ) );

        h_xi1gg_vs_xi1pp->Fill( xi_reco1, xi_prot1 );
        h_xi2gg_vs_xi2pp->Fill( xi_reco2, xi_prot2 );
      }

      if ( num_diproton>0 ) {
        h_mgg_vs_mpp->Fill( diphoton_mass[j], max_diproton_mass );
        h_ygg_vs_ypp->Fill( diphoton_rapidity[j], max_diproton_mass_rap );
        if ( fabs( diphoton_mass[j]-max_diproton_mass )<diphoton_mass[j]*0.15 ) {
          h_mgg_vs_mpp_cand->Fill( diphoton_mass[j], max_diproton_mass );
          h_ygg_vs_ypp_cand->Fill( diphoton_rapidity[j], max_diproton_mass_rap );
          h_xi1gg_vs_xi1pp_cand->Fill( xi_reco1, xi_prot1 );
          h_xi2gg_vs_xi2pp_cand->Fill( xi_reco2, xi_prot2 );
if ( fabs( diphoton_rapidity[j]-max_diproton_mass_rap )<0.1 ) {

TLorentzVector p1, p2;
p1.SetPtEtaPhiM( diphoton_pt1[j], diphoton_eta1[j], diphoton_phi1[j], 0. );
p2.SetPtEtaPhiM( diphoton_pt2[j], diphoton_eta2[j], diphoton_phi2[j], 0. );
cout << "---------> " << (p1+p2).M() << endl;

cout << " ---> event: " << run_id << ":" << lumisection << ":" << event_number << endl;
cout << "      diphoton: pt=" << diphoton_pt[j] << ", mass=" << diphoton_mass[j] << ", dphi=" << diphoton_dphi[j] << "\t" << met << endl;
cout << "      single photon: pt=" << diphoton_pt1[j] << ", eta=" << diphoton_eta1[j] << ", phi=" << diphoton_phi1[j] << endl
     << "                     pt=" << diphoton_pt2[j] << ", eta=" << diphoton_eta2[j] << ", phi=" << diphoton_phi2[j] << endl;
}
        }
      }

      h_num_vtx_1mm->Fill( diphoton_vertex_vtx1mmdist[j] );
      h_num_vtx_2mm->Fill( diphoton_vertex_vtx2mmdist[j] );
      h_num_vtx_5mm->Fill( diphoton_vertex_vtx5mmdist[j] );
      h_num_vtx_1cm->Fill( diphoton_vertex_vtx1cmdist[j] );

      h_met_vs_pt->Fill( met, diphoton_pt[j] );
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

  const float lumi_b = 5.060924481910, // fb-1
              lumi_c = 1.490748474431; // fb-1

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
    plot_3hists( "diphoton_dphi", top_label, h_diphoton_dphi, h_diphoton_dphi_1tag, h_diphoton_dphi_2tag );
    plot_3hists( "diphoton_rapidity", top_label, h_diphoton_rap, h_diphoton_rap_1tag, h_diphoton_rap_2tag );
    plot_3hists( "diphoton_closest_vtx", top_label, h_diphoton_closestvtx, h_diphoton_closestvtx_1tag, h_diphoton_closestvtx_2tag );
    plot_3hists( "diphoton_vtx_numtracks", top_label, h_diphoton_ntrk, h_diphoton_ntrk_1tag, h_diphoton_ntrk_2tag );
    plot_3hists( "num_vertex", top_label, h_num_vtx, h_num_vtx_1tag, h_num_vtx_2tag );
    plot_3hists( "event_met", top_label, h_met, h_met_1tag, h_met_2tag );
  }

  {
    Canvas c( "mass_balance", top_label );
    h_mgg_vs_mpp->Draw( "p" );
    h_mgg_vs_mpp->SetMarkerStyle( 24 );
    h_mgg_vs_mpp_cand->Draw( "p same" );
    h_mgg_vs_mpp_cand->SetMarkerStyle( 20 );
    h_mgg_vs_mpp_cand->SetMarkerColor( kRed+1 );
    c.Prettify( h_mgg_vs_mpp );
    c.DrawDiagonal( h_mgg_vs_mpp );
    c.Save( "png" );
    c.Save( "pdf" );
  }
  {
    Canvas c( "xi1_balance", top_label );
    h_xi1gg_vs_xi1pp->Draw( "p" );
    h_xi1gg_vs_xi1pp->SetMarkerStyle( 24 );
    c.Prettify( h_xi1gg_vs_xi1pp );
    h_xi1gg_vs_xi1pp_cand->Draw( "p same" );
    h_xi1gg_vs_xi1pp_cand->SetMarkerStyle( 20 );
    h_xi1gg_vs_xi1pp_cand->SetMarkerColor( kRed+1 );
    c.DrawDiagonal( h_xi1gg_vs_xi1pp );
    c.Save( "png" );
    c.Save( "pdf" );
  }
  {
    Canvas c( "xi2_balance", top_label );
    h_xi2gg_vs_xi2pp->Draw( "p" );
    h_xi2gg_vs_xi2pp->SetMarkerStyle( 24 );
    c.Prettify( h_xi2gg_vs_xi2pp );
    h_xi2gg_vs_xi2pp_cand->Draw( "p same" );
    h_xi2gg_vs_xi2pp_cand->SetMarkerStyle( 20 );
    h_xi2gg_vs_xi2pp_cand->SetMarkerColor( kRed+1 );
    c.DrawDiagonal( h_xi2gg_vs_xi2pp );
    c.Save( "png" );
    c.Save( "pdf" );
  }
  {
    Canvas c( "rapidity_balance", top_label );
    h_ygg_vs_ypp->Draw( "p" );
    h_ygg_vs_ypp->SetMarkerStyle( 24 );
    c.Prettify( h_ygg_vs_ypp );
    h_ygg_vs_ypp_cand->Draw( "p same" );
    h_ygg_vs_ypp_cand->SetMarkerStyle( 20 );
    h_ygg_vs_ypp_cand->SetMarkerColor( kRed+1 );
    c.DrawDiagonal( h_ygg_vs_ypp );
    c.Save( "png" );
    c.Save( "pdf" );
  }
  {
    Canvas c( "diphoton_pt_vs_met", top_label );
    h_met_vs_pt->Draw( "colz" );
    c.Prettify( h_met_vs_pt );
    c.Save( "png" );
    c.Save( "pdf" );
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
    c.Save( "png" );
    c.Save( "pdf" );
  }
  {
    Canvas c( "num_proton", top_label );
    h_num_proton->Sumw2();
    h_num_proton->Draw();
    h_num_proton->SetMarkerStyle( 20 );
    c.Prettify( h_num_proton );
    c.Save( "png" );
    c.Save( "pdf" );
  }

}
