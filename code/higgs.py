
import scipy as np
import ROOT
from ROOT import TBox, TArrow, TCanvas, TH1D, TH2D, TLine, TMath, TRandom3, TROOT, TLatex, TFile, TLegend, TLegendEntry, gROOT, gDirectory, kTRUE


# ------- FUNCTIONS ---------

def GetMassDistribution( Itype=1, scalefactor=1.00 ):
	dir = gDirectory
	file = TFile("Histograms_fake.root","READ")
	dir.cd()

	if(Itype==125):
		h_mass = file.Get("h_m4l_Higgs125_fake").Clone("h_mass")
	if(Itype==200):
		h_mass = file.Get("h_m4l_Higgs200_fake").Clone("h_mass")
	if(Itype==1):
		h_mass = file.Get("h_m4l_ZZ_fake").Clone("h_mass")
	if(Itype==2):
		h_mass = file.Get("h_m4l_data_fake").Clone("h_mass")


	# scale histograms
	Nbins = h_mass.GetNbinsX()
	for i in range(0,Nbins):
		mu_bin = h_mass.GetBinContent(i)
		h_mass.SetBinContent(i, scalefactor*mu_bin)

	h_mass.SetDirectory(0)

	file.Close()

	return h_mass

# ------------------------------------------------------------------

def AddText(txt_x=0.50, txt_y=0.50, txt="dummy", txt_size=0.045, txt_angle=0, Alignment="left", UseNormalizedSize=1, txt_color=1):
	txt_align=12

	if Alignment == "left":
		txt_align = 12
	if Alignment == "right":
		txt_align = 32
	if Alignment == "right":
		txt_align = 22

	t1 = TLatex(txt_x, txt_y, txt)
	
	if UseNormalizedSize: 
		t1.SetNDC(kTRUE)
	
	t1.SetTextSize(txt_size)
	t1.SetTextAlign(txt_align)
	t1.SetTextAngle(txt_angle)
	t1.SetTextColor(txt_color)
	t1.Draw()

	return
	
# -----------------------------------------------------------------

def MassPlot(Irebin):
	gROOT.Clear()
	gROOT.Delete()

	# Prepare canvas
	c1 = TCanvas("canvas1","Standard Canvas",600,400)
	c1.SetLeftMargin(0.125)
	c1.SetBottomMargin(0.125)
	c1.cd()

	# Prepare histograms
	# get from file and rebin
	h_sig = GetMassDistribution(125)
	h_sig.Rebin(Irebin)
	h_bgr = GetMassDistribution(1)
	h_bgr.Rebin(Irebin)
	h_data = GetMassDistribution(2)
	h_data.Rebin(Irebin)

	# cumulative S+B
	h_SB = h_bgr.Clone("h_SB")
	h_SB.Reset()
	for i in range(0,h_bgr.GetNbinsX()):
		h_SB.SetBinContent(i, h_sig.GetBinContent(i)+h_bgr.GetBinContent(i) )

		print(" REBINNED HISTOGRAM: bin", i, "Ndata = ",h_data.GetBinContent(i),"/n")

	# prepare and plot histograms
	Data_max = h_data.GetBinContent(h_data.GetMaximumBin())
	Ymax_plot = 1.20*(Data_max + 0.1*(Data_max))
	h_SB.SetFillColor(7)
	h_SB.SetAxisRange(0.,Ymax_plot,"Y")
	h_SB.SetAxisRange(0.,400.,"X")
	h_bgr.SetFillColor(2)
	h_SB.Draw("hist")
	h_bgr.Draw("same")
	h_bgr.Draw("axis same")
	h_data.Draw("e same")

	# image format
	AddText(0.900,0.035,"4-lepton invariant mass [GeV]",0.060,0.,"right")
	AddText(0.040,0.900, "Number of events / "+str(h_bgr.GetBinWidth(1))+" GeV", 0.060, 90., "right")
	leg1 = TLegend(0.65,0.65,0.90,0.85)
	leg1.SetBorderSize(0); leg1.SetFillColor(0);
	leg1a = leg1.AddEntry(h_bgr, "SM(ZZ)", "f"); leg1a.SetTextSize(0.04);
	leg1b = leg1.AddEntry(h_SB, "Higgs", "f"); leg1b.SetTextSize(0.04);
	leg1.Draw()

	# save image
	c1.Print("Plots/MassPlot_rebin"+str(Irebin)+".pdf")

	return

# ----------------------------------------------------------------------------

def IntegratePoissonFromRight(mu, Nobs):
	integral = 1.
	for i in range(0,int(Nobs)):
		integral -= TMath.Poisson(i, mu)
	return integral

# --------------------------------------------------------------------------

def SignificanceOpt(Lumi_scalefactor = 1.00):
	print("------- Significance Optimization -------/n")
	
	# Prepare histograms
	print("/n Info: Mass distribution in the 4 lepton channel/n")
	
	h_sig = GetMassDistribution(125, Lumi_scalefactor)
	h_bgr = GetMassDistribution(1, Lumi_scalefactor)
	h_data = GetMassDistribution(2, Lumi_scalefactor)


	# Compute significance for various mass windows
	# Define histograms for mass windows to try
	h_masswindow = TH1D("h_masswindow", "", 250, 0.,25.) #full width \in [0,25]
	h_masswindow_exp = TH1D("h_masswindow_expected", "", 250, 0.,25.)
	h_masswindow_obs = TH1D("h_masswindow_observed", "", 250, 0.,25.)	

	# Loop over various mass windows
	for i in range( 0,h_masswindow.GetNbinsX() ):
		masswindow_fullwidth = h_masswindow.GetBinCenter(i)

		print("	Trying as mass window: ",masswindow_fullwidth," GeV\n")
		
		# determine events in the mass window for each event type
		Ndata = h_data.Integral( h_data.FindBin(125-0.5*masswindow_fullwidth),h_data.FindBin(125+0.5*masswindow_fullwidth) )
		Nbgr =  h_bgr.Integral( h_bgr.FindBin(125-0.5*masswindow_fullwidth),h_bgr.FindBin(125+0.5*masswindow_fullwidth) )
		Nsig = h_sig.Integral( h_sig.FindBin(125-0.5*masswindow_fullwidth),h_sig.FindBin(125+0.5*masswindow_fullwidth) )

		print("Number of events: ", Nsig, "\n")
		
		if( (Nbgr+Nsig)<1 ): continue 

		# Computing EXPECTED significance
		pvalue_exp = IntegratePoissonFromRight(Nbgr, Nsig+Nbgr)
		significance_exp = ROOT.Math.gaussian_quantile_c(pvalue_exp,1)
		h_masswindow_exp.SetBinContent(i, significance_exp)

		# Computing OBSERVED significance
		pvalue_obs = IntegratePoissonFromRight(Nbgr, Ndata)
		significance_obs = ROOT.Math.gaussian_quantile_c(pvalue_obs,1)
		h_masswindow_obs.SetBinContent(i, significance_obs)
		
	# Getting maximum
	MaxBin_exp = h_masswindow_exp.GetMaximumBin()
	OptimumSign_exp = h_masswindow_exp.GetBinContent(MaxBin_exp)
	OptMassWindow_exp = h_masswindow_exp.GetXaxis().GetBinCenter(MaxBin_exp)

	print("Luminosity factor: ", Lumi_scalefactor,"\n")
	print("EXPECTED significance - Optimal mass window\n")
	print("	Expected significance: ",OptimumSign_exp,"\n")
	print("	Mass window: ", OptMassWindow_exp, "\n")
		
	if( abs(Lumi_scalefactor-1.0)<0.01 ):
	
		MaxBin_obs = h_masswindow_obs.GetMaximumBin()
		OptimumSign_obs = h_masswindow_obs.GetBinContent(MaxBin_obs)
		OptMassWindow_obs = h_masswindow_obs.GetXaxis().GetBinCenter(MaxBin_obs)

		print("Luminosity factor: ", Lumi_scalefactor,"\n")
		print("OBSERVED significance - Optimal mass window\n")
		print("	Observed significance: ",OptimumSign_exp,"\n")
		print("	Mass window: ", OptMassWindow_exp, "\n")

	# plot and save histograms
	canvas1 = TCanvas( "canvas1","Standard Canvas",600,400);  
  	canvas1.SetLeftMargin(0.125);
  	canvas1.SetBottomMargin(0.125); 
  	canvas1.cd(); 
 	
	h_masswindow_exp.SetLineColor(1);
	h_masswindow_exp.SetLineWidth(2);
	h_masswindow_obs.SetLineColor(4);
	h_masswindow_obs.SetLineWidth(2);

	h_masswindow_exp.SetAxisRange(-1.,6.,"Y");
	h_masswindow_exp.Draw("l");
	if( abs(Lumi_scalefactor-1.00)<0.01 ):
		h_masswindow_obs.Draw("l same");

	# axes
	AddText( 0.900, 0.035, "Mass window GeV",0.060, 0.,"right"); # X-axis
	AddText( 0.040, 0.900, "Significance" ,0.060,90.,"right");   # Y-axis    
	AddText( 0.225, 0.825, "Luminosity scalefactor = "+str(Lumi_scalefactor),0.050, 0.,"left");            

	AddText( 0.700, 0.200, "Expected significance",0.050, 0.,"right",1,1);                        
	if( abs(Lumi_scalefactor-1.00)<0.01 ):
		AddText( 0.700, 0.300, "Observed significance",0.050, 0.,"right",1,4);                        

	canvas1.Print("Plots/Significance_Optimization_lumiscalefactor"+str(Lumi_scalefactor)+".pdf");

	return OptimumSign_exp;

# ---------------------------------------------------------------------------

# =========================================================================
		
# MAIN - RUNNING THE FUNCTIONS

# MassPlot(20)

i=1.0
sign = []
while i<5.5:
	sign = SignificanceOpt(i)
	i += 0.1



	
	
	
	

	

	
