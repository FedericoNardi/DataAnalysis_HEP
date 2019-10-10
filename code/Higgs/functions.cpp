#include "TBox.h"
#include "TArrow.h"
#include "TCanvas.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TLine.h"
#include "TMath.h"
#include "TRandom3.h"
#include "TROOT.h"
#include "TLatex.h"
#include "TFile.h"
#include "TLegend.h"
#include "TLegendEntry.h"
#include "Math/QuantFuncMathCore.h"
#include "Math/ProbFuncMathCore.h" // for ROOT::Math::gaussian_cdf

#include <iostream>
using namespace std;

TH1D * GetMassDistribution(int Itype, double scalefactor){

    // Get histogram from file
    TH1D *h_mass = 0;
    TDirectory* dir = gDirectory;
    TFile *file = new TFile("Histograms_fake.root", "READ");
    dir->cd();

    // Higgs 125
    if(Itype == 125){
        h_mass = (TH1D*) file->Get("h_m4l_Higgs125_fake")->Clone("h_mass");

    }

    // Higgs 200
    if(Itype == 125){
        h_mass = (TH1D*) file->Get("h_m4l_Higgs200_fake")->Clone("h_mass");

    }

    // ZZ background SM
    if(Itype == 1){
        h_mass = (TH1D*) file->Get("h_m4l_ZZ_fake")->Clone("h_mass");
    }

    // data
    if(Itype == 2){
        h_mass = (TH1D*) file->Get("h_m4l_ZZ_fake")->Clone("h_mass");
    }

    // scale histograms
    int Nbins = h_mass->GetNbinsX();
    for(int i_bin =1; i_bin < Nbins; i_bin++){
        double mu_bin = h_mass->GetBinContent(i_bin);
        h_mass->SetBinContent( i_bin, scalefactor * mu_bin );
    }

    file->Close();

    // return histograms
    return h_mass;
}

void MassPlot(int Irebin){
    // plot SM+Higgs+data; rebinning only for plot
    gROOT->Clear();
    gROOT->Delete();

    // Prepare canvas
    TCanvas* canvas1 = new TCanvas("canvas1","Standard Canvas",600,400);
    canvas1->SetLeftMargin(0.125);
    canvas1->SetBottomMargin(0.125);
    canvas1->cd();

    // Get histograms from files
    TH1D *h_sig, *h_bgr, *h_data;
    h_sig = GetMassDistribution(125);
    h_bgr = GetMassDistribution(1);
    h_data = GetMassDistribution(2);

    // Rebin histograms
    h_sig->Rebin(Irebin);
    h_bgr->Rebin(Irebin);
    h_data->Rebin(Irebin);

    // Prepare S+B cumulative histograms
    TH1D *h_sig_bgr = (TH1D*) h_bgr->Clone("h_sig_bgr");
    h_sig_bgr->Reset();

    for(int i_bin = 1; i_bin < h_bgr->GetNbinsX(); i_bin++){
        h_sig_bgr->SetBinContent(i_bin, h_sig->GetBinContent(i_bin) + h_bgr->GetBinContent(i_bin));
        cout << " REBINNED HISTOGRAM: bin " << i_bin << ", Ndata = " << (int)h_data->GetBinContent(i_bin) << endl;

    }

    // Prepare histograms and plot them
    double Data_max = h_data->GetBinContent( h_data->GetMaximumBin() );
    double Ymax_plot = 1.10*(Data_max + TMath::Sqrt(Data_max));

    h_sig_bgr->SetFillColor(7);
    h_sig_bgr->SetAxisRange(0., Ymax_plot, "Y");
    h_sig_bgr->SetAxisRange(0., 400., "X");

    h_bgr->SetFillColor(2);

    h_sig_bgr->Draw("hist");

    h_bgr->Draw("same");
    h_bgr->Draw("axis same");
    h_data->Draw("e same");

    // Image format
    AddText(0.900,0.035, "4-lepton invariant mass [GeV]", 0.060, 0., "right"); //x-axis
    AddText(0.040,0.900, Form("Number of events / %3.1f GeV", h_bgr->GetBinWidth(1)), 0.060, 90. , "right"); //y-axis

    TLegend *leg1 = new TLegend(.65,.65,.90,.85);
    leg1->SetBorderSize(0);
    leg1->SetFillColor(0);
    TLegendEntry *leg1a = leg1->AddEntry(h_bgr, "SM(ZZ)", "f");
    leg1a->SetTextSize(0.04);
    TLegendEntry *leg1b = leg1->AddEntry(h_sig_bgr, "Higgs", "f");
    leg1b->SetTextSize(0.04);
    leg1->Draw();

    // prepare output
    canvas1->Print(Form("./MassPlot_rebin%d.jpg", Irebin));

}

void AddText( Double_t txt_x, Double_t txt_y, const char * txt, Double_t txt_size,
              Double_t txt_angle, const char * Alignment, Int_t UseNormalizedSize, Int_t txt_color){
  Int_t txt_align = 12;
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

}
