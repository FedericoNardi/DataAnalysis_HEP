/*
 * Project:        Exercises 11.4
 * File:           PoissonError.C
 * Author:         Ivo van Vulpen
 * Version (date): 1.0 (01.05.2011)
 *
 * Copyright (C) 2011-2013, Ivo van Vulpen
 * All rights reserved.
 *
 * Description:

 * Compute Poisson errors on observed number of events. This code
 * implements a few of the examples discussed in detail in an article
 * by Bob Cousins -
 * Why isn't every physicist a Bayesian?  Robert D. Cousins, (UCLA)
 * UCLA-HEP-94-005, Sep 1994. 27pp.
 * Published in Am.J.Phys.63:398,1995.
 *
 * This code is distributed with the solution manual to the book
 *
 * Data Analysis in High Energy Physics: A Practical Guide to Statistical Methods,
 * Wiley-VCH (2013),
 * O. Behnke, K. Kroeninger, G. Schott, Th. Schoerner-Sadenius (editors)
 */

#include "TH1D.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TROOT.h"
#include "TLatex.h"
#include "Math/ProbFuncMathCore.h"

#include <iostream>
using namespace std;

//-- Functions in this file
void PoissonError(int n_obs, double * ErrorRange, int IErrorType = 4, int Iplot = 0);
void AddText( double txt_x = 0.50, double txt_y = 0.50, const char * txt = "dummy", double txt_size = 0.045, double txt_angle = 0., const char * Alignment = "left", int UseNormalizedSize = 1, int txt_color =1 );



//============================
void ComputeAllErrorRegions(){
//============================

  double ErrorRange[3];
  int n_obs = 3;
  printf("\n");
  printf("  Observed number of events = %d\n\n",n_obs);
  printf("  Error region:\n");
  PoissonError(n_obs,ErrorRange,1); //  Method 1: Classical central
  printf("    Method 1 (classical central):  %5.2f - %5.2f\n",ErrorRange[0],ErrorRange[2]);
  PoissonError(n_obs,ErrorRange,2); //  Method 2: Likelihood ratio
  printf("    Method 2 (likelihood ratio):   %5.2f - %5.2f\n",ErrorRange[0],ErrorRange[2]);
  PoissonError(n_obs,ErrorRange,3); //  Method 3: Bayesian  central (flat prior - left/right integral = 16% )
  printf("    Method 3 (Bayes central - 1):  %5.2f - %5.2f\n",ErrorRange[0],ErrorRange[2]);
  PoissonError(n_obs,ErrorRange,4); //  Method 4: Bayesian  central (flat prior - integrate pdf (equaly prob) 68%)
  printf("    Method 4 (Bayes central - 2):  %5.2f - %5.2f\n",ErrorRange[0],ErrorRange[2]);
  printf("\n");
  
  return;

} // end ComputeAllErrorRegions()


//===========================================================================
void PoissonError(int n_obs, double * ErrorRange, int IErrorType, int Iplot){
//===========================================================================

// ----------------------------------------------------------------------
// Input: o n_obs   (observed number of events)
//        o IErrorType (type of error to be returned to user)
//          IErrorType: 1 = Classical central
//                      2 = Likelihood ratio 
//	                    3 = Bayesian  central (flat prior - left/right integral = 16% )
//                      4 = Bayesian  central (flat prior - integrate pdf (equaly prob) 68%)
//          For detailed discussion see:
//          Why isn't every physicist a Bayesian?
//          Robert D. Cousins, (UCLA) . UCLA-HEP-94-005, Sep 1994. 27pp. 
//          Published in Am.J.Phys.63:398,1995.
//
// Output:  ErrorRange
// Extra:   Iplot:   1/2 for plots on each of the methods for n_obs=4 Hindeloopen
//
// Howto:   root> double ErrorRange[3] (left, central, right)
//          root> PoissonError(n_obs,ErrorRange,1)  // method 1
// ----------------------------------------------------------------------

  //-----------------
  //-- Standard stuff
  //-----------------
  gROOT->Clear();
  gROOT->Delete();
  TCanvas *canvas1 = 0;
  if(Iplot != 0){
     canvas1 = new TCanvas( "canvas1","Standard Canvas",600,400);  
     canvas1->SetLeftMargin(0.125);
     canvas1->SetBottomMargin(0.125); 
  }

  //-- Define variables and prepare histograms
  double Lambda_min = 0.;
  double Lambda_max = n_obs + 6*sqrt( n_obs );
  int    Nsteps_lambda = int(1e4);	           
  TH1D *h_likelihood = new TH1D("h_likelihood", "Likelihood for lambda   ", Nsteps_lambda, Lambda_min,Lambda_max);
  TH1D *h_min2loglikelihood = new TH1D("h_min2loglikelihood", "Likelihood for lambda   ", Nsteps_lambda, Lambda_min,Lambda_max);
  TH1D *h_pdf_full   = new TH1D("h_pdf_full",   "Pdf for lambda          ", Nsteps_lambda, Lambda_min,Lambda_max);

  
  //-- Define a few standard parameters
  double Integral_fraction_1sigma = ROOT::Math::gaussian_cdf(-1.,1.,0.); // roughly 15.87 %


  // ===================================================================================
  // [A] Prepare Likelihood: Loop over hypotheses of mean of the 'true' Poisson (lambda)
  // ===================================================================================
  double lambda = 0.;
  for(Int_t i_bin_lambda = 1; i_bin_lambda <= Nsteps_lambda; i_bin_lambda++){
    lambda = h_likelihood->GetBinCenter(i_bin_lambda);
    if(fabs(lambda)< 1e-9) lambda = 1e-9;
    
    //-- Compute Poisson probability == Likelihood
    Double_t PoissonProbability = TMath::Poisson(n_obs,lambda);
    Double_t LogLikelihood      = -2.*TMath::Log(PoissonProbability); 

    //-- Fill histogram for Likelihood and probability density function
    h_likelihood->Fill(lambda, LogLikelihood);
    h_min2loglikelihood->Fill(lambda, LogLikelihood);
    h_pdf_full->Fill(lambda, PoissonProbability*1);                  // Bayes prior is constant '(*1)'    
  } // end loop over lambda hypotheses

  //-- save characteristic values
  int    bin_central        = h_likelihood->GetMinimumBin();              // bin with smallest -2*Log(likelihood)
  double LogLikelihood_min  = h_likelihood->GetBinContent(bin_central);
  double Lambda_central     = h_likelihood->GetBinCenter(bin_central);



  // ========================
  // [B] Compute error region
  // ========================

  //-- Compute error region using various options
  double Lambda_1sigma_low = -1.;
  double Lambda_1sigma_up  = -1.;

  // -----------------------------------------------------
  // Option 1: Frequentist (Classical central)
  // -----------------------------------------------------
  if(IErrorType == 1){
    int Nobs_max = n_obs+100;  

    for(int i_bin = 1; i_bin < h_pdf_full->GetNbinsX(); i_bin++){ // loop over lambda
      lambda = h_pdf_full->GetBinCenter(i_bin);    
      double Poisson_sum_low = 0.;
      double Poisson_sum_up  = 0.;

      // lower value
      for(int i_obs = n_obs; i_obs < Nobs_max; i_obs++){ // loop from n_obs to infinity
         Poisson_sum_low += TMath::Poisson(i_obs,lambda);
      } 
      if(Poisson_sum_low > Integral_fraction_1sigma && Lambda_1sigma_low < 0.) {
       	 Lambda_1sigma_low = lambda;
      }

      // upper value
      for(int i_obs = 0; i_obs <= n_obs; i_obs++){ // loop from 0 to n_obs
         Poisson_sum_up += TMath::Poisson(i_obs,lambda);
      }
      if(Poisson_sum_up < Integral_fraction_1sigma && Lambda_1sigma_up < 0.) {
         Lambda_1sigma_up = lambda;
       }        
    } // end loop over lambda  
  } // end IErrorType == 1


  // -----------------------------------------------------------------------
  // Option 2: Likelihood ratio
  //           Find region with delta -2log(likelihood) w.r.t minimum <1.
  // -----------------------------------------------------------------------
  if(IErrorType == 2){
    Lambda_1sigma_low = Lambda_central;
    Lambda_1sigma_up  = Lambda_central;

    int i_bin_lambda     = -1;
    double LogLikelihood = -1.;
    int KeepMoving = 0;
  
    //-- left boundary  
    i_bin_lambda = bin_central;
    KeepMoving = (i_bin_lambda == 1) ? 0 : 1;
    while ( KeepMoving ){
      i_bin_lambda--;
      Lambda_1sigma_low = h_likelihood->GetBinCenter(i_bin_lambda);
      LogLikelihood     = h_likelihood->GetBinContent(i_bin_lambda);
      // decide when to stop
      if(LogLikelihood > LogLikelihood_min+1. || i_bin_lambda==1){
        KeepMoving = 0;
      }
      //printf("lambda = %5.2f LL = %5.2f -> BOUNDARY = %5.2f \n",lambda,LogLikelihood,Lambda_1sigma_low);
    }

    //-- right boundary  
    i_bin_lambda = bin_central;
    KeepMoving = (i_bin_lambda == h_likelihood->GetNbinsX()) ? 0 : 1;
    while ( KeepMoving ){
      i_bin_lambda++;
      Lambda_1sigma_up = h_likelihood->GetBinCenter(i_bin_lambda);
      LogLikelihood    = h_likelihood->GetBinContent(i_bin_lambda);
      // decide when to stop
      if(LogLikelihood > LogLikelihood_min+1. || i_bin_lambda==h_likelihood->GetNbinsX()){
        KeepMoving = 0;
      }
      //printf("lambda = %5.2f LL = %5.2f -> BOUNDARY = %5.2f \n",lambda,LogLikelihood,Lambda_1sigma_up);
    }
    
  } // end IErrorType == 2


  // --------------------------------------------------------------
  // Option 3: Bayes central: flat prior
  //           Find left-error (16%) and right-error (16%) from PDF
  // --------------------------------------------------------------
  if(IErrorType == 3){
    double integral_full     = h_pdf_full->Integral();
    double integral_fraction = 0.;

    for(int i_bin = 1; i_bin < h_pdf_full->GetNbinsX(); i_bin++){
       lambda = h_pdf_full->GetBinCenter(i_bin);
       integral_fraction = h_pdf_full->Integral(1,i_bin)/integral_full;  // fraction of PDF < i_bin
       if(integral_fraction > Integral_fraction_1sigma && Lambda_1sigma_low < 0.) {
       	  Lambda_1sigma_low = lambda;
       }
       if(integral_fraction > (1.-Integral_fraction_1sigma) && Lambda_1sigma_up < 0.) {
       	  Lambda_1sigma_up = lambda;
       }
    }  // end loop over bins

  } // end IErrorType == 3


  // --------------------------------------------------------------
  // Option 4: Bayes central: flat prior
  //           Equal probability
  // --------------------------------------------------------------
  if(IErrorType == 4){
  
    int Nbins = h_pdf_full->GetNbinsX();
    double Integral_tot  = h_pdf_full->Integral(1,Nbins);

    //-- start values
    int    i_bin_mostlikely = h_pdf_full->GetMaximumBin();  
    double i_bin_left       = i_bin_mostlikely;
    double i_bin_right      = i_bin_mostlikely;
    double prob_left         = h_pdf_full->GetBinContent(i_bin_mostlikely);
    double prob_right        = h_pdf_full->GetBinContent(i_bin_mostlikely);
  
    int    KeepMoving  = 1;
    int    Direction   = ( i_bin_mostlikely != Nbins ) ? 1 : -1;  // start right unless you are in the right most bin
    double integral    = 0.;

    while( KeepMoving ){      

      //-- Move right
      if( Direction == 1 ){      
	      i_bin_right += 1;
 	      prob_right = h_pdf_full->GetBinContent(i_bin_right); 
      }
      
      //-- Move left
      if( Direction == -1 ){      
	      i_bin_left -= 1;
 	      prob_left = h_pdf_full->GetBinContent(i_bin_left); 
      }
      
      //-- Decide which way to go    
      Direction = (prob_right < prob_left ) ? -1 : 1;
      if( Direction ==  1 && i_bin_right == Nbins ){ Direction = -1; } 
      if( Direction == -1 && i_bin_left  ==  1    ){ Direction =  1; } 
      
      integral = h_pdf_full->Integral(i_bin_left, i_bin_right) / Integral_tot;

      //-- Check if you can stop
      if( integral > (1.-2.*Integral_fraction_1sigma) ){
     	KeepMoving = 0;
      }

     // printf(" %d %d | %5.2f %5.2f | %5.2f %d\n",i_bin_left, i_bin_right, lik_left, lik_right, integral, KeepMoving);      

    } // end KeepMoving

    Lambda_1sigma_low = h_pdf_full->GetBinCenter(i_bin_left);
    Lambda_1sigma_up  = h_pdf_full->GetBinCenter(i_bin_right);
  } // end IErrorType == 4




  //-----------------------
  //-- Plots
  //-----------------------


  //-----------------  
  //-- Plot option 1: classic central
  //-----------------  
  if(IErrorType == 1 && Iplot != 0){
  
     Int_t Nobs_max = n_obs + 25;
     TH1D *h_1   = new TH1D("h_1",   "poisson probability (sec)      ", Nobs_max+1, -0.5,Nobs_max+0.5);
     TH1D *h_up  = new TH1D("h_up",  "poisson probability (n >= nobs)", Nobs_max+1, -0.5,Nobs_max+0.5);
     TH1D *h_low = new TH1D("h_low", "poisson probability (n <= nobs)", Nobs_max+1, -0.5,Nobs_max+0.5);  
     double lambda = 0.00;
     if(Iplot == 1) {
     	lambda = Lambda_1sigma_low;
     }
     if(Iplot == 2) {
     	lambda = Lambda_1sigma_up;
     }
     double PoissonProb_integrated_low = 0.000;
     double PoissonProb_integrated_up  = 0.000;
     for(Int_t i_obs = 0; i_obs <= Nobs_max; i_obs++){
       Double_t PoissonProb = TMath::Poisson(i_obs,lambda);
       h_1->Fill(Double_t(i_obs), PoissonProb);
       if(i_obs <= n_obs){
          h_low->Fill(Double_t(i_obs), PoissonProb);
          PoissonProb_integrated_low += PoissonProb;
       }
       if(i_obs >= n_obs){
          h_up->Fill(Double_t(i_obs), PoissonProb);
          PoissonProb_integrated_up += PoissonProb;
       }
     }
     printf(" Prob_low = %5.2f     Prob_low = %5.2f \n",PoissonProb_integrated_low,PoissonProb_integrated_up);

     h_1->SetAxisRange(-0.5,16.5,"X");
     h_1->Draw();
     AddText( 0.900, 0.035, "Number of events",0.060,0.,"right");             // X-axis
     AddText( 0.050, 0.900, "Poisson probability",0.060,90.,"right");         // Y-axis
     AddText( 0.800, 0.750, Form("#lambda = %5.2f",lambda),0.060,0.,"right");      if(Iplot == 1){
       h_up->SetFillColor(5);   
       h_up->Draw("same");
       h_1->Draw("axis same");
       canvas1->Print("./PoissonError_frequentist_low.gif");
     }
     if(Iplot == 2){
       h_low->SetFillColor(5);   
       h_low->Draw("same");
       h_1->Draw("axis same");
       canvas1->Print("./PoissonError_frequentist_up.gif");
     }
  } // end plot with option 1 (classic central)

  //-----------------
  //-- Plot option 2: likelihood ratio
  //-----------------
   if(IErrorType == 2 && Iplot != 0){
      if(Iplot ==  1){ 
  	     //h_pdf_full->SetAxisRange(0.,20.,"X");
         //h_pdf_full->SetAxisRange(0.,40.,"Y"); 
         h_pdf_full->Draw();         AddText( 0.900, 0.035, "#lambda",0.060,0.,"right");          canvas1->Print("./PoissonError_likelihood_pdf.gif");
      }
      if(Iplot ==  2){ 
	     h_min2loglikelihood->SetAxisRange(0.,20.,"X");
         h_min2loglikelihood->SetAxisRange(0.,40.,"Y"); 
         h_min2loglikelihood->Draw();         AddText( 0.900, 0.035, "#lambda",0.060,0.,"right");          canvas1->Print("./PoissonError_likelihood.gif");
      }
      if(Iplot ==  3){ 
  	     h_min2loglikelihood->SetAxisRange(1.50,8.00,"X");
         h_min2loglikelihood->SetAxisRange(3.00,5.00,"Y"); 
         h_min2loglikelihood->Draw();         AddText( 0.900, 0.035, "#lambda",0.060,0.,"right");          canvas1->Print("./PoissonError_likelihood_zoom.gif");
      }
   }   // end plot 2 (likelihood ratio)


  //-----------------------
  //-- Plot option 3 and 4: Bayes
  //-----------------------
  TH1D *h_pdf_1sigma = (TH1D*)h_pdf_full->Clone("h_pdf_1sigma"); h_pdf_1sigma->Reset();
  if( (IErrorType == 3 || IErrorType == 4) && Iplot != 0){

      //-- prepare histogram for 1 sigma region
      double lambda = 0.;
      double prob   = 0.;
      double integral_left  = 0.;
      double integral_tot   = 0.;
      double integral_right = 0.;
      for(int i_bin = 1; i_bin < h_pdf_full->GetNbinsX(); i_bin++){
        lambda = h_pdf_full->GetBinCenter(i_bin);      
        prob   = h_pdf_full->GetBinContent(i_bin);       
        if(lambda <= Lambda_1sigma_low || lambda >= Lambda_1sigma_up){
          h_pdf_1sigma->SetBinContent(i_bin,prob);
  	    }
        integral_tot += prob;
        if(lambda <= Lambda_1sigma_low){ integral_left  += prob; }
        if(lambda >= Lambda_1sigma_up) { integral_right += prob; }
      }
      integral_left  = integral_left  / integral_tot;
      integral_right = integral_right / integral_tot;
      printf("Integral error regions:\n");
      printf("  Left (%5.3f) and Right (%5.3f)\n",integral_left,integral_right);      
      
      //-- draw plot
      h_pdf_full->Draw();
      h_pdf_1sigma->SetFillColor(5);
      h_pdf_1sigma->Draw("same");
      h_pdf_full->Draw("axis same");
      AddText( 0.900, 0.035, "#lambda",0.060,0.,"right");                // X-axis
      AddText( 0.050, 0.900, "Probability",0.060,90.,"right");   // Y-axis
      AddText( 0.75, 0.800, Form("N(obs) = %d",n_obs),0.050,0.,"left");  

      if(IErrorType == 3){
      	 canvas1->Print("./PoissonError_bayes_leftright.gif");
      }

      if(IErrorType == 4){
      	 canvas1->Print("./PoissonError_bayes_equalprob.gif");
      }

   }   // end plot 3 or 4 (Bayes)




  //-----------------------
  //-- Return range to user
  //-----------------------
  ErrorRange[0] = Lambda_1sigma_low;
  ErrorRange[1] = Lambda_central;
  ErrorRange[2] = Lambda_1sigma_up;

  //----------------------------
  //-- print to screen for users
  //----------------------------
  if(Iplot != 0){
    printf("\n");
    printf(" Nobs = %d\n",n_obs);
    if(IErrorType == 1){ printf(" Method: Classical central     \n");}
    if(IErrorType == 2){ printf(" Method: Likelihood ratio      \n");}
    if(IErrorType == 3){ printf(" Method: Bayesian (central)    \n");}
    if(IErrorType == 4){ printf(" Method: Bayesian (equal prob) \n");}
    printf("  --> Confidence level interval = (%4.2f,%4.2f)\n",Lambda_1sigma_low,Lambda_1sigma_up);
  }
	

 //-- delete stuff
 //delete h_likelihood;
 //delete h_pdf_full;
   
  return;
} // end PoissonError()




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



