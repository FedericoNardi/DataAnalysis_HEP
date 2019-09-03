/*
 * Project:        Exercise 11.5
 * File:           Pdfs.h
 * Author:         Ivo van Vulpen, Aart Heijboer
 * Version (date): 1.0 (23.06.2013)
 *
 * Copyright (C) 2013, Ivo van Vulpen, Aart Heijboer
 * All rights reserved.
 *
 * Description:
 * PDF functions for the two parameters, mass, and angle
 * for the signal and background components
 *
 * This code is distributed with the solution manual to the book
 *
 * Data Analysis in High Energy Physics: A Practical Guide to Statistical Methods,
 * Wiley-VCH (2013),
 * O. Behnke, K. Kroeninger, G. Schott, Th. Schoerner-Sadenius (editors)
 */



#ifndef PDFSINCHH
#define PDFSINCHH

#include "TF1.h"
#include "TF2.h"
#include "TMath.h"

//=====
// Mass
//=====

const double maximal_mass = 1000.;

double turn_on_curve(double mass )
{
  return 0.5 + 0.5 *  TMath::Erf( (mass - 150) / 30 );
}


double pdf_for_mass_background_elec( double mass )
{
  const double slope = 2e-3;

  static bool init   = false;
  static double norm = 0;

  if (!init)
    {
      // compute normalization
      double xmin(0), xmax(600);
      int N=1000;
      double dx = ( xmax - xmin ) / N;
      double c = 0;

      for (int i=0; i< N; i++ )
	{
	  double x = i * dx;
	  c += exp( -slope* x ) * turn_on_curve( x+100);
	}
      
      norm = 1.0 / ( c*dx );
      init = true;
    }

  return norm * exp( -slope* mass ) * turn_on_curve( mass+100);
}



double pdf_for_mass_background( double mass )
{
  const double slope = 1e-2;

  static bool init = false;
  static double norm=0;

  if (!init)
    {
      // compute normalization
      double xmin(0), xmax(maximal_mass);
      int N=1000;
      double dx = ( xmax - xmin ) / N;
      double c = 0;

      for (int i=0; i< N; i++ )
	{
	  double x = i * dx;
	  c += exp( -slope* x ) * turn_on_curve( x);
	}
      
      norm = 1.0 / ( c*dx );
      init = true;
    }

  // cout << norm << " " << mass << " " <<  exp( -slope* mass ) << " " <<  turn_on_curve( mass) << endl;

  return norm * exp( -slope* mass ) * turn_on_curve( mass);
}





double pdf_for_mass_background_r( double* pmass , double* par)
{
  return  pdf_for_mass_background(*pmass);
}

//-----------------------------------------
TF1* get_f_mass_bgr(double sigma_mass = 0 )  //-- background: pdf of mass distribution 
// sigma_mass is not used in this case
//-----------------------------------------
{   
  TF1 *f1 = new TF1("myfunc", pdf_for_mass_background_r,0,1000,0);
  return f1;

  // sigma_mass =0; // do-nothing statement to get rid of compiler warning
  // static TF1* f_mass_bgr = 0;
  //if (!f_mass_bgr) 
  //  {
  //   f_mass_bgr = new TF1("f_mass_bgr","[0]*exp(-[1]*x)",0.,maximal_mass);
  //   f_mass_bgr->SetParameter(0, 1);
  //    f_mass_bgr->SetParName(0, "normalization");
  //    f_mass_bgr->SetParameter(1, 1e-2);
  //    // normalize
  //    f_mass_bgr->SetParameter(0, f_mass_bgr->GetParameter(0) / f_mass_bgr->Integral( f_mass_bgr->GetXmin(), f_mass_bgr->GetXmax() ) );
  //  }
  //return f_mass_bgr;
}


//-------------------
TF1* get_f_sigma_signal()  // return the sigma distribution for signal events
//-------------------
{ 
  const double sigma_min = 5;
  const double sigma_max = 25;

  static TF1* f_sigma_signal = 0;
  if (!f_sigma_signal) 
    {
      f_sigma_signal = new TF1("f_sigma_signal","[0]",sigma_min,sigma_max);
      f_sigma_signal->SetParameter(0, 1.0 / ( sigma_max - sigma_min ) );     
    }
  return f_sigma_signal;
}

const double sigma_min = 5;
const double sigma_max = 25;

inline double pdf_for_sigma_background(double sigma)
{ 
  if ( sigma < sigma_min || sigma > sigma_max ) return 0;

  const double N =  2.0 / ( sigma_max*sigma_max - sigma_min*sigma_min);
  return N * sigma;
}
inline double pdf_for_sigma_signal(double sigma)
{ 
  if ( sigma < sigma_min || sigma > sigma_max ) return 0;

  const double N =  1.0 / ( sigma_max - sigma_min);
  return N;
}

 
//-------------------
TF1* get_f_sigma_background()  // return the sigma distribution for signal events
//-------------------
{ 

  static TF1* f_sigma_background = 0; 
  if (!f_sigma_background)  
    {  
      //f_sigma_background = new TF1("f_sigma_background","[0]",sigma_min,sigma_max); 
      //f_sigma_background->SetParameter(0, 1.0 / ( sigma_max - sigma_min ) );     

      f_sigma_background = new TF1("f_sigma_background","[0]*x",sigma_min,sigma_max); 
      f_sigma_background->SetParameter(0, 2.0 / ( sigma_max*sigma_max - sigma_min*sigma_min) );      
    } 
  return f_sigma_background; 
}



//-------------------
TF1* get_f_mass_sig( double mass_central)  //-- signal: pdf of mass distribution - for a resolution of 15 GeV
//-------------------
{ 
  static TF1* f_mass_sig = 0;
  if (!f_mass_sig) 
    {
      f_mass_sig = new TF1("f_mass_sig","gaus",0.,maximal_mass);
      f_mass_sig->SetParameter(0, 10.0);     
      f_mass_sig->SetParameter(2, 17.5);
    }
  f_mass_sig->SetParameter(1, mass_central);
  f_mass_sig->SetParameter(0, f_mass_sig->GetParameter(0) / f_mass_sig->Integral( f_mass_sig->GetXmin(), f_mass_sig->GetXmax() ) );
  return f_mass_sig;
}


//-------------------
TF1* get_f_mass_sig( double mass_central, double sigma )  //-- signal: pdf of mass distribution
// nb: this is P( mass | mass_central, sigma )
//-------------------
{ 
  static TF1* f_mass_sig = 0;
  if (!f_mass_sig) 
    {
      f_mass_sig = new TF1("f_mass_sig","gaus",0.,maximal_mass);
      f_mass_sig->SetParameter(0, 10.0);     
      f_mass_sig->SetParameter(2, 15.0);
    }
  f_mass_sig->SetParameter(1, mass_central);
  f_mass_sig->SetParameter(2, sigma);
  f_mass_sig->SetParameter(0, 1.0 / sqrt( 2*TMath::Pi() * sigma * sigma )  );
  
  return f_mass_sig;
}






inline double pdf_for_mass_signal( double mass, double mass_central, double sigma=15 )
{ 
  const double sqrt2pi = sqrt( 2*TMath::Pi());
  const double z = (mass - mass_central) / sigma ;
  double G = exp(-0.5 * z*z );
  G /= sqrt2pi * sigma;
  return G;
}

inline double pdf_for_mass_signal_elec( double mass, double mass_central, double sigma=17 )
{ 
  const double sqrt2pi = sqrt( 2*TMath::Pi());
  const double z = (mass - mass_central) / sigma ;
  double G = exp(-0.5 * z*z );
  G /= sqrt2pi * sigma;
  return G;
}



inline double pdf_for_mass_signal250( double mass )
{
  return pdf_for_mass_signal( mass, 250 );

}
inline double pdf_for_mass_signal_elec250( double mass )
{
  return pdf_for_mass_signal_elec( mass, 250 );

}


void fill_hist( TH1& h , double (*f)(double) )
{
  for (int i=0; i<= h.GetNbinsX(); i++)
    {
      double x = h.GetBinCenter( i );
      h.SetBinContent(i, f(x) );
    }
}


//-------------------
//TF2* get_f_mass_and_sigma_sig( double mass_central )  //-- signal: pdf of mass distribution - for a resolution of 15 GeV
// p ( mass, sigma | mass_central )/
//-------------------
//{ 
  // todo/
//  exit(1);
//}



#endif
