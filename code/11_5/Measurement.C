/*
 * Project:        Exercise 11.5
 * File:           Measurement.C
 * Author:         Ivo van Vulpen, Aart Heijboer
 * Version (date): 1.0 (23.06.2013)
 *
 * Copyright (C) 2013, Ivo van Vulpen, Aart Heijboer
 * All rights reserved.
 *
 * Description:
 * A code skeleton for the measurement part. For more details, see README file.
 *
 * This code is distributed with the solution manual to the book
 *
 * Data Analysis in High Energy Physics: A Practical Guide to Statistical Methods,
 * Wiley-VCH (2013),
 * O. Behnke, K. Kroeninger, G. Schott, Th. Schoerner-Sadenius (editors)
 */

#include "TArrow.h"
#include "TBox.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH3D.h"

#include "TLatex.h"
#include "TLegend.h"
#include "TLegendEntry.h"
#include "TLine.h"
#include "TMath.h"
#include "TRandom3.h"
#include "TROOT.h"
#include "TString.h"
#include "TStyle.h"

#include "Pdfs.h"

#include <fstream>
//#include "Math/QuantFuncMathCore.h"
//#include "Math/ProbFuncMathCore.h"  // for ROOT::Math::gaussian_cdf
 
#include <vector> 
#include <iostream>
using namespace std;
void AddText( double txt_x = 0.50, double txt_y = 0.50, const char * txt = "dummy", double txt_size = 0.045, double txt_angle = 0., const char * Alignment = "left", int UseNormalizedSize = 1, int txt_color =1 );


TH1D* hist_mass_muon;


TH2D swapxy( const TH2& src )
{  
  int nbinx = src.GetNbinsX();
  int nbiny = src.GetNbinsY();
  double x1 = src.GetXaxis()->GetBinLowEdge(1);
  double x2 = src.GetXaxis()->GetBinLowEdge(nbinx)+src.GetXaxis()->GetBinWidth(nbinx);
  double y1 = src.GetYaxis()->GetBinLowEdge(1);
  double y2 = src.GetYaxis()->GetBinLowEdge(nbiny)+src.GetYaxis()->GetBinWidth(nbiny);  

  TH2D dst ( (TString(src.GetName())+"_swapped").Data(),
	     (TString(src.GetTitle())+"_swapped").Data(),
	     nbiny,y1,y2,nbinx,x1,x2 ); 
  
  for (int i=1; i<= nbinx; i++ ) 
    for (int j=1; j<= nbiny; j++ ) 
      {
	dst.SetBinContent( j, i, src.GetBinContent(i,j ) );
      }
  return dst;
}



//--------------------
//--  global variables
//--------------------
TRandom3 *R = new TRandom3(0);
double lumi_data                          = 1000.;  // pb^-1
double systematic_uncertainty_luminosity  = 0.00;   // 5% uncertainty on the luminosity
int Iprint_global                         = 0;      // print info to user in every step
int I_teststatistic_definition_global     = 1;      // definition test-statistic (1=Tevatron, 2=LHC)


//-- Global histograms
TH1D* h_mass_data                    = 0; // Data
TH1D* h_template_bgr_SM              = 0; // Standard Model templates
TH1D* h_template_sig_SM              = 0; //    ,,
TH1D* h_template_bgr_scaled          = 0; // Scaled templates
TH1D* h_template_sig_scaled          = 0; //    ,,

TH1D* h_template_bgr_normalized          = 0; // Normalized templates
TH1D* h_template_sig_normalized          = 0; //    ,,

TH1D* h_template_sig_plus_bgr_scaled = 0; //    ,,
TH1D* h_cumulative_bgr               = 0; // Cumulative scaled templates
TH1D* h_cumulative_sig_plus_bgr      = 0; //    ,,
TH1D* h_template_mzprime_hypothesis  = 0;


void Init_Templates( double m_zprime = 250., double scalefactor_bgr_xs = 1., double scalefactor_sig_xs = 1);
//void Plot_MassDistribution( double m_zprime, double scalefactor_bgr = 1., double scalefactor_sig = 1.);

struct Evt
{
  double mass;
  double sigma;
  double angle;
};

struct DataSet : public vector<Evt>
{


  DataSet(int size=0) 
  { 
    resize( size ); 
  }
  
  TH1D mass_histogram()
  {
    TH1D r("hmass","hmass",100,0,600);
    for (unsigned i=0; i<size(); i++) r.Fill( at(i).mass );
    return r;
  }

  void dump()
  {
    for (unsigned i=0; i<size(); i++)
      {
	Evt& evt = at(i);
	printf("%-3d mass = %4.3f +/- %4.3f GeV \n", i, evt.mass, evt.sigma );
      }
  }

};



DataSet GenerateDataSet( int Nsig, int Nbgr , double central_signal_mass = 250 )
{
  // generate the background events

  DataSet r(0);
  Evt evt;
   
  for (int i=0; i< Nbgr; i++ )
    {
      evt.mass  = get_f_mass_bgr()->GetRandom();
      evt.sigma = get_f_sigma_background()->GetRandom();
      r.push_back( evt );
    }
  
  for (int i=0; i< Nsig; i++)
    {
      evt.sigma = get_f_sigma_signal()->GetRandom();
      evt.mass  = get_f_mass_sig( central_signal_mass , evt.sigma) ->GetRandom();
      r.push_back( evt );
    }

  // todo reduce the data-sample by contricting the fit range
  

  return r;
}

double unbinned_likelihood_mass_elec( DataSet& dataset, double Nsig, double Nbgr , double mass_hyp, bool use_uncertainties = true )
{
  double NLL = 0;
  for(unsigned i= 0; i< dataset.size() ; i++ )
    {
      Evt& evt = dataset.at(i);

      double P_mass;
      if ( use_uncertainties )
	{
	  P_mass = 
	    Nbgr * pdf_for_mass_background_elec( evt.mass ) +
	    Nsig * pdf_for_mass_signal_elec( evt.mass, mass_hyp ) ;
	}      
      if ( P_mass < 0 ) // can happen for negative nsig
	{
	  return 0;
	}
      
      NLL -= log( P_mass );
    }

  NLL += Nsig+Nbgr;
  return NLL;
}

double unbinned_likelihood_mass( DataSet& dataset, double Nsig, double Nbgr , double mass_hyp, bool use_uncertainties = true )
{
  double NLL = 0;
  
  //TF1* pdf_mass_bgr = get_f_mass_bgr();
  //TF1* pdf_mass_sig = get_f_mass_sig(mass_hyp);

  for(unsigned i= 0; i< dataset.size() ; i++ )
    {
      Evt& evt = dataset.at(i);

      double P_mass;
      if ( use_uncertainties )
	{
	  P_mass = 
	    Nbgr * pdf_for_mass_background( evt.mass )                  * pdf_for_sigma_background( evt.sigma ) +
	    Nsig * pdf_for_mass_signal( evt.mass, mass_hyp, evt.sigma ) * pdf_for_sigma_signal( evt.sigma );
	}
      else 
	{
	  P_mass = 
	    Nbgr * pdf_for_mass_background( evt.mass ) +
	    Nsig * pdf_for_mass_signal( evt.mass, mass_hyp );
	}
      
      if ( P_mass < 0 ) // can happen for negative nsig
	{
	  return 0;
	}
      
      NLL -= log( P_mass );
    }

  NLL += Nsig+Nbgr;
  return NLL;
}

TH1D get_profile_1d_y( TH2D& h) // profile away the x-axis
{ 
  TH1D r = *h.ProjectionY(Form("%s_pry",h.GetName() ) );
  
  r.Scale(0);
    
  for (int iy =1 ; iy<= h.GetNbinsY(); iy++)
    {
      double nll_min = 1e10;

      for (int ix =1 ; ix<= h.GetNbinsX(); ix++)
	{
	  int bin = h.GetBin(ix,iy);
	  double nll = h.GetBinContent( bin );
	  if ( nll < nll_min ) nll_min = nll;	  
	}
      r.SetBinContent( iy, nll_min );
    }
  return r;
}

TH1D get_profile_1d_x( TH2D& h, bool get_yvals=false ) // profile away the y-axis
{
 
  TH1D r = *h.ProjectionX(  get_yvals?  Form("%s_prx",h.GetName() ):Form("%s_yvals",h.GetName() ) );
  r.Scale(0);
    
  for (int ix =1 ; ix<= h.GetNbinsX(); ix++)
    {
      double nll_min = 1e10;
      double y_min = 0;

      for (int iy =1 ; iy<= h.GetNbinsY(); iy++)
	{
	  int bin = h.GetBin(ix,iy);
	  double nll = h.GetBinContent( bin );
	  if ( nll < nll_min ) 
	    {
	      nll_min = nll;	  
	      y_min   = h.GetYaxis()->GetBinCenter( iy );
	    }
	}

      if (get_yvals ) {
	r.SetBinContent( ix, y_min );
	cout << ix << " " << y_min << endl;
      }
      else            r.SetBinContent( ix, nll_min );
    }
  return r;
}




vector<double> analyze_likelihood( const TH1& h )
{
  const double delta_log_lik = 1; // 1 is valid if the histogram is -2 log lik

  vector<double> r(3);
  r[1] = r[2] = 0;
    
  int b = h.GetMinimumBin();
  r[0] = h.GetXaxis()->GetBinCenter( b );
  
  double y0 = h.GetBinContent( b );
  
  // get left point where we pas delta-log-lik = 1
  for ( unsigned i=b ; i>0 ; i--)
    {
      if ( h.GetBinContent( i ) > y0 + delta_log_lik ) 
	{
	  double excess =  h.GetBinContent( i ) - ( y0 + delta_log_lik ); 
	  double dy     =  h.GetBinContent( i )- h.GetBinContent( i-1 );
	  double f = excess / dy;

	  r[1] = h.GetBinCenter( i ) - f * dy;
	  break;
	}
    }
  for (unsigned i =b ; i <= h.GetNbinsX() ; i++)
    {
      if ( h.GetBinContent( i ) > y0 + delta_log_lik ) 
	{
	  double excess =  h.GetBinContent( i ) - ( y0 + delta_log_lik ); 
	  double dy     =  h.GetBinContent( i )- h.GetBinContent( i+1 );
	  double f = excess / dy;

	  r[2] = h.GetBinCenter( i ) - f * dy;
	  break;
	}
    }

  return r;
}



TString format_measurement( vector<double>& v )
{
  return Form("%4.3f +/- %4.3f / %4.3f", v[0],v[1]-v[0],v[2]-v[0] );
}

// profile away the x-axis
TH2D get_profile_2d_yz( TH3D& h , bool set_minimum_to_zero = true)
{
  TH2D r ("likhist2d","likhist;N_{sig};Mass",
	  h.GetNbinsY(), h.GetYaxis()->GetXmin() ,  h.GetYaxis()->GetXmax() , 
	  h.GetNbinsZ(), h.GetZaxis()->GetXmin() ,  h.GetZaxis()->GetXmax() );
  
  double global_min = 1e10;
  
  for (int iy =1 ; iy<= h.GetNbinsY(); iy++)
    for (int iz =1 ; iz<= h.GetNbinsZ(); iz++)
      {
	double nll_min = 1e10;
	
	for (int ix =1 ; ix<= h.GetNbinsX(); ix++)
	  {
	    int bin = h.GetBin(ix,iy,iz);
	    double nll = h.GetBinContent( bin );
	    if ( nll < nll_min ) nll_min = nll;
	  }
	
	if (nll_min < global_min ) global_min = nll_min;
	int rbin = r.GetBin(iy,iz);
	r.SetBinContent( rbin, nll_min );
      }

  if (set_minimum_to_zero)
    {
      for (int ix = 1; ix<= r.GetNbinsX();  ix++)
	for (int iy = 1; iy<= r.GetNbinsY(); iy++)
	  {
	    int b = r.GetBin(ix,iy);
	    r.SetBinContent( b, r.GetBinContent( b ) - global_min );
	  }
    }
  
  return r;
}


TH2D get_profile_2d_xy( TH3D& h , bool set_minimum_to_zero = true)
{
  TH2D r ("likhist2d","likhist;N_{bg};N_{sig}",
	  h.GetNbinsX(), h.GetXaxis()->GetXmin() ,  h.GetXaxis()->GetXmax() , 
	  h.GetNbinsY(), h.GetYaxis()->GetXmin() ,  h.GetYaxis()->GetXmax() );
  
  double global_min = 1e10;
  
  for (int ix =1 ; ix<= h.GetNbinsX(); ix++)
    for (int iy =1 ; iy<= h.GetNbinsY(); iy++)
      {
	double nll_min = 1e10;
	
	for (int iz =1 ; iz<= h.GetNbinsZ(); iz++)
	  {
	    int bin = h.GetBin(ix,iy,iz);
	    double nll = h.GetBinContent( bin );
	    if ( nll < nll_min ) nll_min = nll;
	  }
	
	if (nll_min < global_min ) global_min = nll_min;
	int rbin = r.GetBin(ix,iy);
	r.SetBinContent( rbin, nll_min );
      }

  if (set_minimum_to_zero)
    {
      for (int ix = 1; ix<= r.GetNbinsX();  ix++)
	for (int iy = 1; iy<= r.GetNbinsY(); iy++)
	  {
	    int b = r.GetBin(ix,iy);
	    r.SetBinContent( b, r.GetBinContent( b ) - global_min );
	  }
    }
  
  return r;
 
}


TH2D get_profile_2d( TH3D& h , bool set_minimum_to_zero = true)
{
  return get_profile_2d_yz( h, set_minimum_to_zero );
  
}




/* return -2 log lik */

TH3D get_likelihood_histogram( DataSet& dataset , 
			       bool shift_minimum_to_zero = false , 
			       bool use_uncertainties = true,
			       double fixmass = -1,
			       bool elec = false )
{
  int N = 1;

  //TH3D r("likhist","likhist;N_{bgr};N_{sig};Mass",
  //	 20/N,150,250 ,  //nbgr
  // 45/N, 10,55,    //nsig
  //	 60/N, 235,265); //mass

  double mmax = 280;
  double mmin = 220;
  int mn = 180;

  int Nb = 20;
  int Ns = 20;

  if (fixmass > 0 ) 
    {
      mmax = fixmass+1;
      mmin = fixmass;
      mn   = 1;

      Nb = 200;
      Ns = 800;
      cout << " fixing mass to " << fixmass << endl;
    }

  double   Nb_min = 200;
  double   Nb_max = 400;
  if (elec) 
    {
      Nb_min = 0;
      Nb_max = 100;
    }

  TH3D r("likhist","likhist;N_{bgr};N_{sig};Mass",
	 Nb,Nb_min,Nb_max ,  //nbgr
	 Ns, 0, 80,    //nsig
	 mn, mmin, mmax); //mass

  int ntot = r.GetNbinsX() * r.GetNbinsY() *  r.GetNbinsZ();

  cout << " computing likelihood grid ( " << ntot << " points ) , did : xxxxxxxx " << flush ;
 

  int j=0;
  for (int ix =1 ; ix<= r.GetNbinsX(); ix++)
    for (int iy =1 ; iy<= r.GetNbinsY(); iy++)
      for (int iz =1 ; iz<= r.GetNbinsZ(); iz++)
	{
	  if (j % 100 ==0) printf("\b\b\b\b\b\b\b\b%8d",j);
	  j++;

	  double nbgr = r.GetXaxis()->GetBinCenter( ix );
	  double nsig = r.GetYaxis()->GetBinCenter( iy );
	  double mass = r.GetZaxis()->GetBinCenter( iz );

	  double nll;
	  if (elec)  nll  = 2 * unbinned_likelihood_mass_elec ( dataset ,nsig, nbgr , mass , use_uncertainties);
	  else       nll  = 2 * unbinned_likelihood_mass ( dataset ,nsig, nbgr , mass , use_uncertainties);
	  
	  int bin = r.GetBin(ix,iy,iz);
	  r.SetBinContent( bin, nll);
	}
  
  double m   = r.GetMinimum();
  double max = r.GetMaximum();

  if ( shift_minimum_to_zero )
    {
      for (int ix =1 ; ix<= r.GetNbinsX(); ix++)
  	for (int iy =1 ; iy<= r.GetNbinsY(); iy++)   
	  for (int iz =1 ; iz<= r.GetNbinsZ(); iz++)
  	  {
  	    int bin = r.GetBin(ix,iy,iz);
  	    if ( r.GetBinContent( bin ) != 0 )
  	      {
  		r.SetBinContent( bin, r.GetBinContent(bin) - m );
  	      }
  	    else 
  	      {
  		r.SetBinContent( bin, max - m );
  	      }
  	  }
    }

  return r;
}


//=========================
void makeplots1(){
//=========================

  gStyle->SetPalette(1);

  Init_Templates( 250, 1, 0.1 );

  
  //  gRandom->SetSeed( 11 );
  gRandom->SetSeed(1111 );
  DataSet data;


  // =======================================
  // Plot [1]: data and likelihood functions
  // =======================================
  TCanvas * canvas0 = new TCanvas( "canvas0","Standard Canvas",600,400);  
  canvas0->SetLeftMargin(0.125);
  canvas0->SetBottomMargin(0.125); 

  TH1D h_mass_pdf_sig("h_mass_pdf_sig","h_pdf_sig",600,0,600);
  TH1D h_mass_pdf_bgr("h_mass_pdf_bgr","h_pdf_bgr",600,0,600);
  
  for (int seed = 213; seed < 300; seed++)
    {
      gRandom->SetSeed( seed );
      data = GenerateDataSet( 35, 300 );
      TH1D HH = data.mass_histogram();
      HH.Rebin(2);
      
      double w = HH.GetBinWidth(1);
      
      fill_hist( h_mass_pdf_sig , pdf_for_mass_signal250 );
      fill_hist( h_mass_pdf_bgr , pdf_for_mass_background );
      h_mass_pdf_sig.Scale(35*w);
      h_mass_pdf_bgr.Scale(300*w);
      
      //-- draw histograms
      h_mass_pdf_sig.SetLineWidth(3); h_mass_pdf_sig.SetLineColor(4); h_mass_pdf_sig.SetLineStyle(1);
      h_mass_pdf_bgr.SetLineWidth(3); h_mass_pdf_bgr.SetLineColor(1);
      h_mass_pdf_bgr.SetAxisRange(0.,38.,"Y");
      h_mass_pdf_bgr.SetAxisRange(50.,500.,"X");
      h_mass_pdf_bgr.SetFillColor(kYellow);//DrawCopy("L");
      h_mass_pdf_bgr.DrawCopy("L");
      h_mass_pdf_sig.DrawCopy("Lsame");      
      HH.DrawCopy("Esame");

      // opsmuk
      AddText( 0.900, 0.035, "di-muon mass (GeV)",0.06, 0.,"right"); // X-axis
      AddText( 0.040, 0.900, "Number of events" ,0.06,90.,"right");  // Y-axis

      AddText( 0.620, 0.800, "background pdf",0.06, 0.,"left",1,1);
      AddText( 0.620, 0.700, "Z' signal pdf",0.06, 0.,"left",1,4);

      TLine *l_bgr  = new TLine(300., 33.,325.,33.);
      TLine *l_sig  = new TLine(300., 28.5,325.,28.5);
      l_bgr->SetLineWidth(3); l_bgr->SetLineStyle(1); l_bgr->SetLineColor(1); l_bgr->Draw("same");
      l_sig->SetLineWidth(3); l_sig->SetLineStyle(1); l_sig->SetLineColor(4); l_sig->Draw("same");

      gPad->Update();

      cout << " pres y to redraw seed, anything else to proceed with this dataset " << endl;
      string x ="";
      //cin >> x;
      if (x != "y") break;
    }
    canvas0->SaveAs("./Plots/meas_massdist.eps");
    canvas0->SaveAs("./Plots/meas_massdist.png");

  // =======================================================
  // Plot [2]: sigma distributions for signal and background
  // =======================================================
    TCanvas * canvas1 = new TCanvas( "canvas1","Standard Canvas",600,400);  
    canvas1->SetLeftMargin(0.125);
    canvas1->SetBottomMargin(0.125); 

    TH1F h_sigma_pdf_sig("h_sigma_pdf_sig","h_sigma_pdf_sig",1000,0.,30.);
    TH1F h_sigma_pdf_bgr("h_sigma_pdf_bgr","h_sigma_pdf_bgr",1000,0.,30.);
  
    fill_hist( h_sigma_pdf_sig , pdf_for_sigma_signal );
    fill_hist( h_sigma_pdf_bgr , pdf_for_sigma_background );

    h_sigma_pdf_sig.SetLineWidth(2);  h_sigma_pdf_sig.SetLineColor(4); h_sigma_pdf_sig.SetLineStyle(1);
    h_sigma_pdf_bgr.SetLineWidth(2);  h_sigma_pdf_bgr.SetLineColor(1);

    h_sigma_pdf_bgr.SetFillColor(kYellow);//DrawCopy("L");
    h_sigma_pdf_bgr.DrawCopy("Lhist");
    h_sigma_pdf_sig.DrawCopy("Lsame");
  
    AddText( 0.900, 0.035, "mass uncertainty #sigma (GeV)",0.06, 0.,"right");       // X-axis
    AddText( 0.040, 0.900, "Probability" ,0.06,90.,"right");  // Y-axis
  
    AddText( 0.250, 0.8125, "background pdf",0.06, 0.,"left",1,1);
    AddText( 0.250, 0.700, "Z' signal pdf",0.06, 0.,"left",1,4);

    TLine *l_bgr  = new TLine(3.,0.078, 4.35, 0.078);
    TLine *l_sig  = new TLine(3.,0.065, 4.35, 0.065);
    l_bgr->SetLineWidth(3); l_bgr->SetLineStyle(1); l_bgr->SetLineColor(1); l_bgr->Draw("same");
    l_sig->SetLineWidth(3); l_sig->SetLineStyle(1); l_sig->SetLineColor(4); l_sig->Draw("same");

    canvas1->Print("./Plots/sigma_distribution.eps");
    canvas1->Print("./Plots/sigma_distribution.png");
  
  


  // ================================
  // Plot [3]: 2d likelihood contours (mass measurement)
  // ================================
    
    TCanvas * canvas2 = new TCanvas( "canvas2","Standard Canvas",500,500);  
    canvas2->SetLeftMargin(0.155);
    canvas2->SetBottomMargin(0.125); 

    TH3D h  = get_likelihood_histogram( data, true, true );
    TH2D h2 = get_profile_2d( h );
    h2.SetName("h2");
    TH3D g  = get_likelihood_histogram( data, true, false );
    TH2D g2 = get_profile_2d( g );

    double conts[] = {1,4,9,16}; 
    h2.SetContour(3,conts);
    g2.SetContour(3,conts);
    g2.SetLineStyle(2);
    g2.SetLineWidth(4);
    h2.SetXTitle("");
    h2.SetYTitle("");

    h2.DrawCopy("cont1");
    g2.DrawCopy("cont1,same");
 
 
    AddText( 0.900, 0.035, "N_{sig}",0.06, 0.,"right");  // X-axis
    AddText( 0.040, 0.900, "M_{Z'} (GeV)" ,0.06,90.,"right");  // Y-axis

    AddText( 0.875, 0.23, "using mass resolution",0.05, 0.,"right");
    AddText( 0.875, 0.18, "#Delta(-2log(L))= 1,4,9",0.05, 0.,"right");
  
    AddText( 0.875, 0.81, "#font[12]{disregarding mass resolution}",0.05, 0.,"right");
    AddText( 0.875, 0.76, "#font[12]{#Delta(-2log(L))= 1,4,9}",0.05, 0.,"right");
  
    TArrow *a1 = new TArrow(40., 230., 36., 240., 0.025,"|>"); a1->Draw();
    TArrow *a2 = new TArrow(40., 230., 30., 246., 0.025,"|>"); a2->Draw();
    TArrow *a3 = new TArrow(40., 230., 20., 248., 0.025,"|>"); a3->Draw();

    TLine *l_yes = new TLine(31., 228., 35., 228.);
    TLine *l_no  = new TLine(17., 273., 21., 273.);
    l_yes->SetLineWidth(3); l_yes->SetLineStyle(1);  l_yes->Draw("same");
    l_no->SetLineWidth(3);  l_no->SetLineStyle(2);   l_no->Draw("same");


    canvas2->Print("./Plots/meas_contour.eps");
    canvas2->Print("./Plots/meas_contour.png");

  // =======================================================
  // plot [4] profile likelihood explanation (left plot, 2d)
  // =======================================================

  TCanvas* canvas3 = new TCanvas("canvas3","canvas3",500,500);
  canvas3->SetLeftMargin(0.155);
  canvas3->SetBottomMargin(0.125); 

  TH3D h250    = get_likelihood_histogram( data, true, true , 250);
  TH2D h250xy  = get_profile_2d_xy( h250 );
  TH2D h250yx  = swapxy( h250xy);
  
  h250yx.SetName("h250yx");
  h250yx.SetContour(3,conts);
  h250yx.DrawCopy("cont1");
  
  TH1D hp = get_profile_1d_x( h250yx, true );
  hp.DrawCopy("Csame");
  
  TH1D* hpx = (TH1D*) hp.Clone("hpx");
  for (int i=1; i<= hpx->GetNbinsX(); i++)
    {
      hpx->SetBinContent(i, 300);
    }
  hpx->SetLineStyle(2);
  hpx->Draw("same");


  AddText( 0.900, 0.035, "N_{sig}",0.06, 0.,"right");  // X-axis
  AddText( 0.040, 0.900, "N_{bgr}" ,0.06,90.,"right");  // Y-axis
  
  AddText( 0.875, 0.240, "mass fixed to true value",0.05, 0.,"right");
  AddText( 0.875, 0.185, "#Delta(-2log(L))= 1,4,9",0.05, 0.,"right");
    
  canvas3->Print("./Plots/meas_lik2d_nsig_nbgr.eps");
  canvas3->Print("./Plots/meas_lik2d_nsig_nbgr.png");
  

  // ========================================================
  // plot [5] profile likelihood explanation (right plot, 1d)
  // ========================================================

  TCanvas* canvas4 = new TCanvas("canvas4","canvas4",500,500);
  canvas4->SetLeftMargin(0.155);
  canvas4->SetBottomMargin(0.125); 
     
  TH1D hlik = get_profile_1d_x( h250yx, false );
  int iy = h250yx.GetYaxis()->FindBin(300);
  TH1D* hlik2 = h250yx.ProjectionX("hlik2",iy,iy);  
  double MM = hlik2->GetMinimum();
  for (int i=1; i<= hlik2->GetNbinsX(); i++)
    {
      hlik2->SetBinContent(i, hlik2->GetBinContent(i) - MM);
    }
  hlik.SetAxisRange(10.,70.,"X");
  hlik.SetAxisRange( 0., 6.,"Y");
  hlik.DrawCopy("L");
  hlik2->SetLineStyle(2);
  hlik2->DrawCopy("Lsame");

  vector<double> mx = analyze_likelihood( hlik );
 
  //-- opsmuk  
  TLine *l01 = new TLine( 10., 1.0 , 70., 1.0 );  
  l01->SetLineStyle(1);   l01->Draw("same");
 
  TLine *l02 = new TLine( mx[1], 1.0 ,mx[1], 0.0 ); 
  TLine *l03 = new TLine( mx[2], 1.0 ,mx[2], 0.0 ); 
  l02->Draw("same");
  l03->Draw("same");
  
  AddText( 0.35, 0.650, Form("N_{sig} = %3.2f",mx[0]),0.05,0.,"left");
  AddText( 0.565, 0.650+0.030, Form(" + %3.1f",mx[0]-mx[1]),0.04,0.,"left");
  AddText( 0.565, 0.650-0.010, Form(" -%3.1f",mx[2]-mx[1]),0.04,0.,"left");


  AddText( 0.900, 0.035, "N_{sig}",0.06, 0.,"right");  // X-axis
  AddText( 0.040, 0.900, "-2log(L)" ,0.06,90.,"right");  // Y-axis
 
  canvas4->Print("./Plots/meas_lik1d_nsig.eps");
  canvas4->Print("./Plots/meas_lik1d_nsig.png");
  
  

  // ==================================
  // plot [6] 1d likelihood versus mass
  // ==================================
  
  TCanvas* canvas5 = new TCanvas("canvas5","canvas5",500,500);
  canvas5->SetLeftMargin(0.155);
  canvas5->SetBottomMargin(0.125); 

  TH1D h_lik_mass1 = get_profile_1d_y( h2 );
  TH1D h_lik_mass2 = get_profile_1d_y( g2 );

  hist_mass_muon = (TH1D*) h_lik_mass1.Clone("hist_mass_muon");


  h_lik_mass1.GetXaxis()->SetRangeUser(234,266);
  
  h_lik_mass1.SetMaximum(10);
  h_lik_mass1.DrawCopy("L");
  h_lik_mass2.SetLineStyle(2);
  h_lik_mass2.DrawCopy("Lsame");


  vector<double> m1 = analyze_likelihood( h_lik_mass1 ); 
  vector<double> m2 = analyze_likelihood( h_lik_mass2 );
  
  TLine *l201 = new TLine(234, 1.0 , 266, 1.0 ); 
  l201->Draw("same");

  //l01->SetLineStyle(1); 
  //l01->Draw("same");
 
  TLine *l102 = new TLine( m1[1], 1.0 ,m1[1], 0.0 ); 
  TLine *l103 = new TLine( m1[2], 1.0 ,m1[2], 0.0 ); 

  TLine *l202 = new TLine( m2[1], 1.0 ,m2[1], 0.0 ); 
  TLine *l203 = new TLine( m2[2], 1.0 ,m2[2], 0.0 ); 
  
  l202->SetLineStyle(2);
  l203->SetLineStyle(2);
  
  l102->Draw("same");
  l103->Draw("same");
  l202->Draw("same");
  l203->Draw("same");

  TString meas1 = Form("with res: M_{Z'} = %3.2f#splitline{+%3.1f}{-%3.1f}", m1[0],m1[0]-m1[1],m1[2]-m1[0] );
  TString meas2 = Form("without : M_{Z'} = %3.2f#splitline{+%3.1f}{-%3.1f}", m2[0],m2[0]-m2[1],m2[2]-m2[0] );

  AddText( 0.320,     0.825,       "using mass resolution:"    , 0.04,0.,"left");
  AddText( 0.34,      0.750,       Form("M_{Z'} = %3.1f",m1[0]), 0.04,0.,"left");
  AddText( 0.34+0.17, 0.750+0.030, Form(" +%3.1f",m1[0]-m1[1]), 0.03,0.,"left");
  AddText( 0.34+0.17, 0.750-0.010, Form(" -%3.1f",m1[2]-m1[0]) , 0.03,0.,"left");

  AddText( 0.89,     0.375,       Form("#font[12]{disregarding}") ,0.04,0.,"right");
  AddText( 0.89,     0.325,       Form("#font[12]{mass resolution:}") ,0.04,0.,"right");
  AddText( 0.65,      0.250,       Form("M_{Z'} = %3.1f",m2[0]),0.04,0.,"left");
  AddText( 0.65+0.17, 0.250+0.030, Form(" +%3.1f",m2[0]-m2[1]),0.03,0.,"left");
  AddText( 0.65+0.17, 0.250-0.010, Form(" -%3.1f",m2[2]-m2[0]),0.03,0.,"left");



  //  AddText( 0.35, 0.700, meas2.Data() ,0.05,0.,"left");

  AddText( 0.900, 0.035, "Z' mass (GeV)",0.06, 0.,"right");       // X-axis
  AddText( 0.040, 0.900, "-2 #Delta log(L)" ,0.06,90.,"right");  // Y-axis

  canvas5->Print("./Plots/meas_mass_lik.eps");
  canvas5->Print("./Plots/meas_mass_lik.png");



  // ==================================
  // plot [6] 1d likelihood versus mass
  // ==================================
  
  TCanvas* canvas6 = new TCanvas("canvas6","canvas6",500,500);
  canvas6->SetLeftMargin(0.155);
  canvas6->SetBottomMargin(0.125); 

  TH1D h_lik_nsig1 = get_profile_1d_x( h2 );
  TH1D h_lik_nsig2 = get_profile_1d_x( g2 );
  
  h_lik_nsig1.SetMaximum(10);
  h_lik_nsig1.DrawCopy("L");
  h_lik_nsig2.SetLineStyle(2);
  h_lik_nsig2.DrawCopy("Lsame");

  canvas6->Print("./Plots/meas_nsig_lik.eps");
  canvas6->Print("./Plots/meas_nsig_lik.png");

  return;

} // end makeplots1()



void gen_experiment(ostream& out = cout)
{
  DataSet data = GenerateDataSet( 20, 300 );
  TH3D g  = get_likelihood_histogram( data, true, true );
  TH2D g2 = get_profile_2d( g );
  
  TH1D h_lik_mass1 = get_profile_1d_y( g2 );
  TH1D h_lik_nsig1 = get_profile_1d_x( g2 );

  //g2.DrawCopy("colz");
  //gPad->Update();

  vector<double> m = analyze_likelihood( h_lik_mass1 );
  vector<double> m2 = analyze_likelihood( h_lik_nsig1 );
  
  out << " mass : " << format_measurement( m ) << "   ";
  out << " nsig : " << format_measurement( m2 ) << endl;
  out << flush;

  cout << " mass : " << format_measurement( m ) << "   ";
  cout << " nsig : " << format_measurement( m2 ) << endl;
}





//=================
void run_pes()
//=================
{

  gStyle->SetPalette(1);

  Init_Templates( 250, 1, 2.5 );
  
  ofstream fout("measurement_2013.txt");
  
  int N = 10000;

  for (int i=0; i< N ; i++ )
    {
      cout << i << " of " << N << endl;
      gen_experiment( fout );
    }
  
  fout.close();
}


//=================
void makeplots2(){
//=================

  // ===============================================================
  // print pull distribution from mass-production of toy experiments
  // todo: - add text for the fit parameters
  // ===============================================================

  ifstream f("measurment_250_10k.txt");
  
  TString dum;
  double mass, nsig;
  double errm_up, errm_dwn, errn_up, errn_dwn;

  
  TH1F* hm = new TH1F("hm","hm",40,230, 270 );
  TH1F* hp = new TH1F("hp","hp",60,-3, 3 );
  
 
  while(f)
    {
      string dum;
      f >> dum >> dum >> mass >> dum >> errm_up >> dum >> errm_dwn;
      f >> dum >> dum >> nsig >> dum >> errn_up >> dum >> errn_dwn;

      errm_up*=-1;
      cout << mass << " " << errm_up << " " << errm_dwn << endl;

      double e = errm_up;
      if ( mass < 250 ) e = errm_dwn;

      hm->Fill( mass );
      hp->Fill( ( mass - 250 ) / (e) );
    }
  
  hp->Rebin(2);

  //------------------------------
  //-- plot [1] reconstructed mass
  //------------------------------
  TCanvas* canvas1 = new TCanvas("canvas1","canvas1",500,500);
  canvas1->SetLeftMargin(0.175);
  canvas1->SetBottomMargin(0.125); 
  canvas1->cd();

  hm->DrawCopy();

  AddText( 0.900, 0.035, "di-muon mass"          ,0.06, 0.,"right");  // X-axis
  AddText( 0.040, 0.900, "Number of experiments" ,0.06,90.,"right");  // Y-axis

  canvas1 -> SaveAs("./Plots/massreco.eps");
  canvas1 -> SaveAs("./Plots/massreco.png");

  //----------------
  //-- plot [2] pull
  //----------------
  TCanvas* canvas2 = new TCanvas("canvas2","canvas2",500,500);
  canvas2->SetLeftMargin(0.155);
  canvas2->SetBottomMargin(0.125); 
  canvas2->cd();
  hp->Fit("gaus");
  hp->DrawCopy();

  double fitpargaus[3] = {0.,0.,0.}; 
  TF1 *fgaus = new TF1("fgaus","gaus",-3.,3.);
  hp->Fit(fgaus,"RQ");
  for(int i=0 ; i< 3 ; i++) { 
      fitpargaus[i]= fgaus->GetParameter(i);  
      cout << fitpargaus[i] << endl;
  }

  double mean = fitpargaus[1];
  double rms  = fitpargaus[2];
        
  TH1D* h_pull = (TH1D*)hp->Clone("h_pull");
  h_pull->Draw();
     
  AddText( 0.900, 0.035, "Pull",0.06, 0.,"right");                     // X-axis
  AddText( 0.040, 0.900, "Number of experiments" ,0.06,90.,"right");  // Y-axis

  //-- opsmuk
  //AddText( 0.75, 0.83, "pull",0.06, 0.,"left");

  AddText( 0.40, 0.30, Form("mean = %5.3f",mean),0.05, 0.,"left");
  AddText( 0.40, 0.22, Form("RMS = %6.3f",rms),0.05, 0.,"left");


  canvas2 -> SaveAs("./Plots/pull.eps");
  canvas2 -> SaveAs("./Plots/pull.png");


} // end makeplots2()







//===========================================================================================
void Init_Templates( double m_zprime, double scalefactor_bgr_xs, double scalefactor_sig_xs ){
//===========================================================================================

  //---------------------------------------------------------------
  // Goal: [1] obtain SM template histograms
  //       [2] construct scaled template histograms  
  //       [3] construct cumulative from scaled template histograms  
  //
  // Note: The histograms filled are all global 
  //---------------------------------------------------------------

  // ----------------
  // [1] SM templates
  // ----------------
  printf("\n Init_Templates():\n");
  printf("    [1] Initialising the SM templates\n");

  if(gROOT->FindObject("h_mass_data")){
     gROOT->FindObject("h_mass_data")->Delete();
  }
  if(gROOT->FindObject("h_template_bgr_SM")){
     gROOT->FindObject("h_template_bgr_SM")->Delete();
  }
  if(gROOT->FindObject("h_template_sig_SM")){
     gROOT->FindObject("h_template_sig_SM")->Delete();
  }
  if(gROOT->FindObject("h_template_mzprime_hypothesis")){
     gROOT->FindObject("h_template_mzprime_hypothesis")->Delete();
  }

  TFile *file1 = 0;

  //-- open file and get bgr/sig/data histograms
  TDirectory* dir = gDirectory; 
  file1 = new TFile("DataSamples.root", "READ");
  dir->cd();
  h_mass_data         = (TH1D*) file1->Get("h_mass_data")->Clone("h_mass_data");     
  h_template_bgr_SM   = (TH1D*) file1->Get("h_mass_bgr")->Clone("h_template_bgr_SM"); 
  TString histname;
  histname = Form("h_mass_sig_%d",(int)m_zprime);
  h_template_sig_SM   = (TH1D*) file1->Get(histname)->Clone("h_template_sig_SM");

  //-- rebin
  h_mass_data->Rebin(5);
  h_template_bgr_SM->Rebin(5);
  h_template_sig_SM->Rebin(5);
            
  //-- clean up
  h_template_bgr_SM->SetXTitle("");
  h_template_bgr_SM->SetYTitle("");
  file1->Close();

  // ------------------------------------------------------------------------
  // [2] create scaled templates from SM templates with help of scale factors
  // ------------------------------------------------------------------------
  printf("    [2] Constructing scaled templates\n");

  if(gROOT->FindObject("h_template_bgr_scaled")){
     gROOT->FindObject("h_template_bgr_scaled")->Delete();
  }
  if(gROOT->FindObject("h_template_sig_scaled")){
     gROOT->FindObject("h_template_sig_scaled")->Delete();
  }
  if(gROOT->FindObject("h_template_sig_plus_bgr_scaled")){
     gROOT->FindObject("h_template_sig_plus_bgr_scaled")->Delete();
  }

  //-- define s+b histograms (just copy of SM templates)
  h_template_bgr_scaled          = (TH1D*) h_template_bgr_SM->Clone("h_template_bgr_scaled");    
  h_template_sig_scaled          = (TH1D*) h_template_sig_SM->Clone("h_template_sig_scaled");    

  h_template_bgr_normalized          = (TH1D*) h_template_bgr_SM->Clone("h_template_bgr_normalized");    
  h_template_sig_normalized          = (TH1D*) h_template_sig_SM->Clone("h_template_sig_normalized");   
  
  h_template_bgr_normalized -> Scale( 1.0 / h_template_bgr_normalized->Integral() );
  h_template_sig_normalized -> Scale( 1.0 / h_template_sig_normalized->Integral() );

  h_template_sig_plus_bgr_scaled = (TH1D*) h_template_sig_SM->Clone("h_template_sig_plus_bgr_scaled");  
  h_template_bgr_scaled->Reset();
  h_template_sig_scaled->Reset();
  h_template_sig_plus_bgr_scaled->Reset();
  
  //-- scale the histograms
  int Nbins = h_template_bgr_SM->GetNbinsX();    
  double Nbgr_bin = 0.;
  double Nsig_bin = 0.;
  for(int i_bin = 1; i_bin<=Nbins; i_bin++){
    Nbgr_bin = scalefactor_bgr_xs * h_template_bgr_SM->GetBinContent(i_bin);
    Nsig_bin = scalefactor_sig_xs * h_template_sig_SM->GetBinContent(i_bin);
    h_template_bgr_scaled->SetBinContent( i_bin, Nbgr_bin );
    h_template_sig_scaled->SetBinContent( i_bin, Nsig_bin );
    h_template_sig_plus_bgr_scaled->SetBinContent( i_bin, Nbgr_bin + Nsig_bin );
  }  

  // ------------------------------------------------------
  // [3] create cumulative distributions for b-only and s+b
  // -------------------------------------------------------
  printf("    [3] Constructing cumulative distributions for b and s+b from scaled templates\n\n");

  //-- clean up
  if(gROOT->FindObject("h_cumulative_bgr")){
     gROOT->FindObject("h_cumulative_bgr")->Delete();
  }
  if(gROOT->FindObject("h_cumulative_sig_plus_bgr")){
     gROOT->FindObject("h_cumulative_sig_plus_bgr")->Delete();
  }

  //-- define cumulative histograms
  h_cumulative_bgr          = (TH1D*)h_template_bgr_scaled->Clone("h_cumulative_bgr");
  h_cumulative_sig_plus_bgr = (TH1D*)h_template_sig_plus_bgr_scaled->Clone("h_cumulative_sig_plus_bgr");
  h_cumulative_bgr->Reset();
  h_cumulative_sig_plus_bgr->Reset();
  double Integral_tot_bgr           = h_template_bgr_scaled->Integral();
  double Integral_tot_sig_plus_bgr  = h_template_sig_plus_bgr_scaled->Integral();
  double Integral_left_bgr          = 0.;
  double Integral_left_sig_plus_bgr = 0.;
  for(int i_bin = 1; i_bin<=Nbins; i_bin++){
     Integral_left_bgr          += h_template_bgr_scaled->GetBinContent(i_bin);
     Integral_left_sig_plus_bgr += h_template_sig_plus_bgr_scaled->GetBinContent(i_bin);
     h_cumulative_bgr->SetBinContent(i_bin, Integral_left_bgr / Integral_tot_bgr );
     h_cumulative_sig_plus_bgr->SetBinContent(i_bin, Integral_left_sig_plus_bgr / Integral_tot_sig_plus_bgr);     
   }    

  TH1D * h_dummy = new TH1D("h_template_mzprime_hypothesis"," ", 100000,-0.5,999.5);
  h_template_mzprime_hypothesis = (TH1D*) h_dummy->Clone("h_template_mzprime_hypothesis");
  h_template_mzprime_hypothesis->Reset();
  h_template_mzprime_hypothesis->Fill(m_zprime);
  h_dummy->Delete();

  return;
} // end Init_Templates()




//============================================================================================
void AddText( double txt_x, double txt_y, const char * txt, double txt_size, double txt_angle, 
              const char * Alignment, int UseNormalizedSize, int txt_color){
//============================================================================================
  int txt_align = 12;
  if ( !strcmp(Alignment, "left"))   { txt_align = 12; } // left 
  if ( !strcmp(Alignment, "right"))  { txt_align = 32; } // right
  if ( !strcmp(Alignment, "center")) { txt_align = 22; } // center

  TLatex* t1 = new TLatex( txt_x, txt_y, txt);
  if(UseNormalizedSize) {t1->SetNDC(kTRUE);} // <- use NDC coordinate
  t1->SetTextSize(txt_size);
  t1->SetTextAlign(txt_align);
  t1->SetTextAngle(txt_angle);
  t1->SetTextColor(txt_color);
  t1->Draw();

} // end AddText()



//================
void elecplots(){
//================

  gStyle->SetPalette(1);
  Init_Templates( 250, 1, 0.1 );
  gRandom->SetSeed(1111 );
  DataSet elec_data;

  // =======================
  // plot [1]: electron mass
  // =======================

  TCanvas* canvas1 = new TCanvas("canvas1","canvas1",500,500);
  canvas1->SetLeftMargin(0.155);
  canvas1->SetBottomMargin(0.125); 

  TH1D h_mass_pdf_sig("h_mass_pdf_sig","h_pdf_sig",600,0,600);
  TH1D h_mass_pdf_bgr("h_mass_pdf_bgr","h_pdf_bgr",600,0,600);
  fill_hist( h_mass_pdf_sig , pdf_for_mass_signal_elec250 );
  fill_hist( h_mass_pdf_bgr , pdf_for_mass_background_elec );
    
  int seed = 4;
  
  gRandom->SetSeed( seed );
  
  Evt ee;
  for (int i=0; i< 60; i++)
    {
      ee.mass =  h_mass_pdf_sig.GetRandom();
      elec_data.push_back( ee );
    }
  for (int i=0; i< 60; i++)
    {
      ee.mass =  h_mass_pdf_bgr.GetRandom();	 
      elec_data.push_back( ee );
    }
  
  
  TH1D HH = elec_data.mass_histogram();
  HH.Rebin(4);
  
  double w = HH.GetBinWidth(1);
  
  h_mass_pdf_sig.Scale(60*w);
  h_mass_pdf_bgr.Scale(60*w);
  
  h_mass_pdf_sig.SetLineWidth(3);  h_mass_pdf_sig.SetLineColor(4);
  h_mass_pdf_bgr.SetLineWidth(3);  h_mass_pdf_bgr.SetLineColor(1); h_mass_pdf_bgr.SetFillColor(kYellow);
  
  h_mass_pdf_bgr.SetAxisRange(0.,600.,"X");
  h_mass_pdf_bgr.SetAxisRange(0.,38.,"Y");
  h_mass_pdf_bgr.DrawCopy("Lhist");
  h_mass_pdf_sig.DrawCopy("Lsame");
  
  HH.DrawCopy("Esame");
  AddText( 0.900, 0.035, "di-electron mass (GeV)",0.06, 0.,"right");  // X-axis
  AddText( 0.040, 0.900, "Number of events"      ,0.06,90.,"right");  // Y-axis

  AddText( 0.60, 0.800, "background pdf",0.05, 0.,"left",1,1);
  AddText( 0.60, 0.700, "Z' signal pdf",0.05, 0.,"left",1,4);

  TLine *l_bgr  = new TLine(320., 33., 345.,33.);
  TLine *l_sig  = new TLine(320., 28.5,345.,28.5);
  l_bgr->SetLineWidth(3); l_bgr->SetLineStyle(1); l_bgr->SetLineColor(1); l_bgr->Draw("same");
  l_sig->SetLineWidth(3); l_sig->SetLineStyle(1); l_sig->SetLineColor(4); l_sig->Draw("same");


  gPad->Update();
  canvas1->SaveAs("./Plots/meas_elecsample.png");
  canvas1->SaveAs("./Plots/meas_elecsample.eps");


  // ================================================
  // plot [2]: combine electron and muon measurements  
  // ================================================

  TCanvas* canvas2 = new TCanvas("canvas2","canvas2",500,500);
  canvas2->SetLeftMargin(0.155);
  canvas2->SetBottomMargin(0.125); 

  TH3D h  = get_likelihood_histogram( elec_data, true, true, -1, true );
  TH2D h2 = get_profile_2d( h );
    
  double conts[] = {1,4,9,16}; 
  h2.SetContour(3,conts);
  h2.GetXaxis()->SetTitle("");
  h2.GetYaxis()->SetTitle("");
  h2.DrawCopy("cont1");


  TH1D h_lik_mass1 = get_profile_1d_y( h2 );
  h_lik_mass1.GetXaxis()->SetTitle("");
  h_lik_mass1.SetAxisRange(0.,10.,"Y");
  h_lik_mass1.SetAxisRange(225.,266.,"X");

  h_lik_mass1.SetLineColor(2);
  h_lik_mass1.SetLineStyle(2);

  h_lik_mass1.DrawCopy("L");

  vector<double> m1 = analyze_likelihood( h_lik_mass1 ); 
  
  TLine *l201 = new TLine(225., 1.0 , 266., 1.0 ); 
  l201->Draw("same");
 
  TLine *l102 = new TLine( m1[1], 1.0 ,m1[1], 0.0 ); 
  TLine *l103 = new TLine( m1[2], 1.0 ,m1[2], 0.0 ); 
  l102->SetLineStyle(2);
  l102->SetLineColor(2);
  l103->SetLineStyle(2);
  l103->SetLineColor(2);
  //l102->Draw("same");
  //l103->Draw("same");
       
  hist_mass_muon->GetXaxis()->SetTitle("");
  hist_mass_muon->SetLineStyle(3);
  hist_mass_muon->Draw("Lsame");

  TH1D* htot = (TH1D*) hist_mass_muon->Clone("htot");
  htot->Add( &h_lik_mass1 );

  double MM = htot->GetMinimum();
  for (int i=1; i<= htot->GetNbinsX(); i++)
    {
      htot->SetBinContent(i, htot->GetBinContent(i) - MM);
    }

  htot->GetXaxis()->SetTitle("");
  htot -> SetLineStyle(1);
  htot -> SetLineColor(1);
  htot -> SetLineWidth(2);
  htot -> Draw("Lsame");
 
  vector<double> m2 = analyze_likelihood( *htot ); 
  vector<double> m3 = analyze_likelihood( *hist_mass_muon ); 

  TString meas1 = Form("e sample: M_{Z'} = %3.2f#splitline{+%3.1f}{-%3.1f}", m1[0],m1[0]-m1[1],m1[2]-m1[0] );
  TString meas2 = Form("combined: M_{Z'} = %3.2f#splitline{+%3.1f}{-%3.1f}", m2[0],m2[0]-m2[1],m2[2]-m2[0] );
  TString meas3 = Form("#mu sample: M_{Z'} = %3.2f#splitline{+%3.1f}{-%3.1f}", m3[0],m3[0]-m3[1],m3[2]-m3[0] );

  //AddText( 0.18, 0.800, meas1.Data() ,0.05,0.,"left");
  //AddText( 0.18, 0.700, meas3.Data() ,0.05,0.,"left");
  //AddText( 0.18, 0.600, meas2.Data() ,0.05,0.,"left");

  AddText( 0.900, 0.035, "Z' mass (GeV)",0.06, 0.,"right");       // X-axis
  AddText( 0.040, 0.900, "-2 #Delta log(L)" ,0.06,90.,"right");  // Y-axis


  AddText( 0.20,  0.800, Form("e sample:"),0.05,0.,"left");
  AddText( 0.20,  0.725, Form("M_{Z'} = %3.1f",m1[0]),0.04,0.,"left");
  AddText( 0.37,  0.725+0.030, Form(" +%3.1f",m1[0]-m1[1]),0.03,0.,"left");
  AddText( 0.37,  0.725-0.010, Form(" -%3.1f",m1[2]-m1[0]),0.03,0.,"left");

  AddText( 0.20,  0.600, Form("#mu sample:"),0.05,0.,"left");
  AddText( 0.20,  0.525, Form("M_{Z'} = %3.1f",m3[0]),0.04,0.,"left");
  AddText( 0.37,  0.525+0.030, Form(" +%3.1f",m3[0]-m3[1]),0.03,0.,"left");
  AddText( 0.37,  0.525-0.010, Form(" -%3.1f",m3[2]-m3[0]),0.03,0.,"left");

  AddText( 0.20,  0.400, Form("combined:"),0.05,0.,"left");
  AddText( 0.20,  0.325, Form("M_{Z'} = %3.1f",m2[0]),0.04,0.,"left");
  AddText( 0.37,  0.325+0.030, Form(" +%3.1f",m2[0]-m2[1]),0.03,0.,"left");
  AddText( 0.37,  0.325-0.010, Form(" -%3.1f",m2[2]-m2[0]),0.03,0.,"left");

  TArrow *a1 = new TArrow(238., 8.70, 261.,  8.70, 0.025,"|>"); 
  TArrow *a2 = new TArrow(238., 6.07, 241.5, 6.07, 0.025,"|>"); 
  TArrow *a3 = new TArrow(238., 3.50, 246.,  3.50, 0.025,"|>"); 
  a1->SetLineWidth(1); a1->Draw();
  a2->SetLineWidth(1); a2->Draw();
  a3->SetLineWidth(1); a3->Draw();

  canvas2->SaveAs("./Plots/meas_combine.png");
  canvas2->SaveAs("./Plots/meas_combine.eps");

} // end elecplots()
