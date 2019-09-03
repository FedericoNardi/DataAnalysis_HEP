{

  gROOT->ForceStyle();

  //#include <iostream>

// =======
// Startup
// =======
std::cout << "\n Standard RootLogon \n" << std::endl;

// =========
// Libraries
// =========
// nothing now

//=======================
// Pads
//=======================
gROOT->SetStyle("Plain");       //no grey background on plots

//============================
// Statistics window + Fitting
//============================
gStyle->SetOptLogy(0);
gStyle->SetStatColor(0);        // make the thing transparant
gStyle->SetStatStyle(0);        // to view data behind box too!
gStyle->SetStatBorderSize(0);   // remove shade from stat box
gStyle->SetStatFont(62);        // more readable on screen
gStyle->SetStatH(0.05);         // x,y,h,w define size and position
gStyle->SetStatW(0.225);
gStyle->SetStatX(0.85);
gStyle->SetStatY(0.9);
gStyle->SetPadGridY(0);         // NO grid in Y
gStyle->SetPadGridX(0);         // NO grid in X
gStyle->SetOptStat(0);          // No stat box as default
gStyle->SetOptFit(0);           // No fit results

//=======================
//  Title
//=======================
gStyle->SetOptTitle(0);         // No title as default
gStyle->SetTitleBorderSize(0);  // remove shade from title box
gStyle->SetTitleX(0.3);
gStyle->SetTitleW(0.4);
gStyle->SetTitleH(.092);

//=======================
//  Fonts / Axes
//=======================
// use large Times-Roman fonts
gStyle->SetTextFont(132);
gStyle->SetTextSize(0.08);
gStyle->SetLabelFont(132,"x");
gStyle->SetLabelFont(132,"y");
gStyle->SetLabelFont(132,"z");
gStyle->SetLabelSize(0.05,"x");
gStyle->SetTitleSize(0.06,"x");
gStyle->SetLabelSize(0.10,"y");
gStyle->SetTitleSize(0.06,"y");
gStyle->SetLabelSize(0.10,"z");
gStyle->SetTitleSize(0.06,"z");
gStyle->SetLabelSize(0.05,"X");
gStyle->SetLabelSize(0.05,"Y");
gStyle->SetLabelSize(0.05,"Z");
// Offset various axes
gStyle->SetTitleOffset(1.2,"X"); 
gStyle->SetTitleOffset(1.0,"Y");

//=======================
// Markers
//=======================
gStyle->SetMarkerStyle(20);
gStyle->SetMarkerSize(1.0);       // larger than half a pixel ...
gStyle->SetMarkerColor(kBlack);

//=======================
// Histogram itself
//=======================
gStyle->SetLineColor(kBlack);
gStyle->SetLineWidth(1.5);  
gStyle->SetHistLineWidth(2.0);    // fatter line width for histograms
gStyle->SetHistLineColor(kBlack);

}
