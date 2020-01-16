
import numpy as np
import ROOT
import ctypes
from ROOT import TBox, TArrow, TCanvas, TH1D, TH2D, TLine, TMath, TRandom3, TROOT, TLatex, TFile, TLegend, TLegendEntry, gROOT, gDirectory, kTRUE, kFALSE, TMarker


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
	if Alignment == "center":
		txt_align = 22

	t1 = TLatex(txt_x, txt_y, txt ) 
	
	if UseNormalizedSize: 
		t1.SetNDC(kTRUE)

	#canvas.cd()
	
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
	c1 = TCanvas("c1","Standard Canvas",600,400)
	c1.SetLeftMargin(0.125)
	c1.SetBottomMargin(0.125)
	#c1.SetLogy()
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
	for i in range(1,h_bgr.GetNbinsX()+1):
		h_SB.SetBinContent(i, h_sig.GetBinContent(i)+h_bgr.GetBinContent(i) )

		print(" REBINNED HISTOGRAM: bin", i, "Ndata = ",h_data.GetBinContent(i),"\n")



	h_fit = h_SB.Clone("h_fit")
	h_fit.Reset()

	for i in range( 0,h_bgr.GetNbinsX() ):
		if( h_SB.GetBinCenter(i)>=150 and h_fit.GetBinCenter(i)>=400 ):	h_fit.SetBinContent(i, h_SB.GetBinContent(i) )



	# prepare and plot histograms
	Data_max = h_data.GetBinContent(h_data.GetMaximumBin())
	Ymax_plot = 1.20*(Data_max + 0.1*(Data_max))
	h_SB.SetFillColor(7)
	h_SB.SetAxisRange(1,Ymax_plot,"Y")
	h_SB.SetAxisRange(0,400.,"X")
	h_bgr.SetFillColor(2)

	h_SB.SetTitle("4-LEPTON INVARIANT MASS DISTRIBUTION")


	h_SB.GetXaxis().SetTitle("m [GeV]")
	h_SB.SetTitleSize(0.05,"X")
	h_SB.SetTitleOffset(0.9,"X")


	h_SB.GetYaxis().SetTitle("counts")
	h_SB.SetTitleSize(0.05,"Y")
	h_SB.SetTitleOffset(0.7,"Y")
	

	h_SB.Draw("hist")
	h_bgr.Draw("same")
	h_bgr.Draw("axis same")
	h_data.Draw("e same")

	# image format
	#AddText(0.900,0.035,"4-lepton invariant mass [GeV]",0.060,0.,"right")
	#AddText(0.040,0.900, "Number of events / "+str(h_bgr.GetBinWidth(1))+" GeV", 0.060, 90., "right")
	leg1 = TLegend(0.70,0.70,0.90,0.85)
	leg1.SetBorderSize(1); leg1.SetFillColor(0); #leg1.SetO
	leg1a = leg1.AddEntry(h_bgr, "SM(ZZ)", "f"); leg1a.SetTextSize(0.04);
	leg1b = leg1.AddEntry(h_SB, "Higgs", "f"); leg1b.SetTextSize(0.04);
	leg1.Draw()

	# save image
	c1.Print("Plots/MassPlot_rebin"+str(Irebin)+".pdf")

	h_fit = h_SB.Clone("h_fit")
	h_fit.Reset()

	for i in range( 0,h_bgr.GetNbinsX() ):
		if( h_SB.GetBinCenter(i)>=150 and h_fit.GetBinCenter(i)<=400 ):	h_fit.SetBinContent(i, h_SB.GetBinContent(i) )

	# printing plot for sidebandfit

	c2 = TCanvas("c2","Standard Canvas",600,400)
	c2.SetLeftMargin(0.125)
	c2.SetBottomMargin(0.125)
	c2.cd()

	h_SB.SetAxisRange(1,20,"Y")
	h_fit.SetAxisRange(0,400.,"X")

	h_fit.SetFillColor(8);

	h_SB.SetFillColor(0); h_SB.SetLineColor(9); h_SB.SetTitle("SIDE-BAND FIT REGION")
	h_bgr.SetFillColor(0); h_bgr.SetLineColor(1)

	h_SB.Draw("hist")
	h_bgr.Draw("hist same")
	h_fit.Draw("hist same")
	h_data.Draw("e same")
	leg1.Draw()

	c2.Print("Plots/Sideband_region.pdf")



	return

# ----------------------------------------------------------------------------

def IntegratePoissonFromRight(mu, Nobs):
	integral = 1.
	for i in range(0,int(Nobs)):
		integral -= TMath.Poisson(i, mu)
	return integral

# --------------------------------------------------------------------------

def SignificanceOpt(Lumi_scalefactor = 1.00, plot = "n",speak="") :
	if speak!='mute': print("------- Significance Optimization -------")
	
	# Prepare histograms
	if speak!='mute': print(" Info: Mass distribution in the 4 lepton channel")
	
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

		if speak!='mute': print("	Trying as mass window: ",masswindow_fullwidth," GeV\n")
		
		# determine events in the mass window for each event type
		Ndata = h_data.Integral( h_data.FindBin(125-0.5*masswindow_fullwidth),h_data.FindBin(125+0.5*masswindow_fullwidth) )
		Nbgr =  h_bgr.Integral( h_bgr.FindBin(125-0.5*masswindow_fullwidth),h_bgr.FindBin(125+0.5*masswindow_fullwidth) )
		Nsig = h_sig.Integral( h_sig.FindBin(125-0.5*masswindow_fullwidth),h_sig.FindBin(125+0.5*masswindow_fullwidth) )

		if speak!='mute': print("Number of events: ", Nsig, "\n")
		
		if( (Nbgr+Nsig)<1 ): continue 

		# Computing EXPECTED significance
		pvalue_exp = IntegratePoissonFromRight(Nbgr, Nsig+Nbgr)
		significance_exp = ROOT.Math.gaussian_quantile_c(pvalue_exp,1)
		h_masswindow_exp.SetBinContent(i, significance_exp)

		# Computing OBSERVED significance
		pvalue_obs = IntegratePoissonFromRight(Nbgr, Ndata)
		significance_obs = ROOT.Math.gaussian_quantile_c(pvalue_obs,1)
		h_masswindow_obs.SetBinContent(i, significance_obs)
		if(pvalue_obs<0): continue

		
	# Getting maximum
	MaxBin_exp = h_masswindow_exp.GetMaximumBin()
	OptimumSign_exp = h_masswindow_exp.GetBinContent(MaxBin_exp)
	OptMassWindow_exp = h_masswindow_exp.GetXaxis().GetBinCenter(MaxBin_exp)

	if speak!='mute':
		print("Luminosity factor: ", Lumi_scalefactor,"\n")
		print("EXPECTED significance - Optimal mass window\n")
		print("	Expected significance: ",OptimumSign_exp,"\n")
		print("	Mass window: ", OptMassWindow_exp, "\n")
		
	if( abs(Lumi_scalefactor-1.0)<0.01 ):
	
		MaxBin_obs = h_masswindow_obs.GetMaximumBin()
		OptimumSign_obs = h_masswindow_obs.GetBinContent(MaxBin_obs)
		OptMassWindow_obs = h_masswindow_obs.GetXaxis().GetBinCenter(MaxBin_obs)

		if speak!='mute':
			print("Luminosity factor: ", Lumi_scalefactor,"\n")
			print("OBSERVED significance - Optimal mass window\n")
			print("	Observed significance: ",OptimumSign_obs,"\n")
			print("	Mass window: ", OptMassWindow_obs, "\n")

	# plot and save histograms
	if(plot == "y"):
		canvas1 = TCanvas( "canvas1","Standard Canvas",600,400);  
		canvas1.SetLeftMargin(0.125);
		canvas1.SetBottomMargin(0.125); 
		canvas1.cd(); 
	
		h_masswindow_exp.SetLineColor(1);
		h_masswindow_exp.SetLineWidth(2);
		h_masswindow_obs.SetLineColor(4);
		h_masswindow_obs.SetLineWidth(2);

		h_masswindow_exp.SetAxisRange(-1.,6.,"Y");
		h_masswindow_exp.SetStats(kFALSE)

		h_masswindow_exp.SetTitle("OPTIMAL MASS WINDOW")

		h_masswindow_exp.GetXaxis().SetTitle("mass window [GeV]")
		h_masswindow_exp.GetYaxis().SetTitle("significance Z")

		leg1 = TLegend(0.70,0.70,0.90,0.85)
		leg1.SetBorderSize(1); leg1.SetFillColor(0); #leg1.SetO
		leg1a = leg1.AddEntry(h_masswindow_exp, "expected", "l"); leg1a.SetTextSize(0.04);




		h_masswindow_exp.Draw("l");
		if( abs(Lumi_scalefactor-1.00)<0.01 ):
			h_masswindow_obs.Draw("l same");
			leg1b = leg1.AddEntry(h_masswindow_obs, "observed", "l"); leg1b.SetTextSize(0.04);

		leg1.Draw()                        

		canvas1.Print("Plots/Significance_Optimization_lumiscalefactor"+str(Lumi_scalefactor)+".pdf");

	return OptimumSign_exp;

# ---------------------------------------------------------------------------

def SideBandFit(irebin=1):

	print("\n Side-band fit \n\n")

	# prepare histograms
	h_bgr = GetMassDistribution(1)
	h_data = GetMassDistribution(2)

	# rebin
	h_bgr.Rebin(irebin)
	h_data.Rebin(irebin)
	print("INFO: Rebinning histograms with factor",irebin,". Binwidth: ",h_data.GetBinWidth(1))

	# Loop over scalefactor (alpha_bgr)
	h_sf_bgr = TH1D("h_sf_bgr","",100,0.8,1.5)

	for i in range(1,h_sf_bgr.GetNbinsX()+1): 
		sf_bgr = h_sf_bgr.GetBinCenter(i)

		# Loop over bins, compute loglikelihood, save in histogram
		loglik = 0.
		for j in range(1,h_data.GetNbinsX()+1):

			# Signal-free region
			m4lepBin = h_data.GetBinCenter(j)
			NObsBin = h_data.GetBinContent(j)

			if ( m4lepBin>=150. and m4lepBin<=400.): 
				MeanBgdBin = sf_bgr * h_bgr.GetBinContent(j)

				if MeanBgdBin>0.: loglik += TMath.Log( TMath.Poisson(NObsBin,MeanBgdBin) )

		h_sf_bgr.SetBinContent(i,-2.*loglik)


	# Interpret the likelihood distribution

	h_sf_bgr_rescaled = h_sf_bgr.Clone("h_sf_bgr_rescaled")


	# find minimum
	
	MinBin = h_sf_bgr.GetMinimumBin()
	Minimum = h_sf_bgr.GetBinContent(MinBin)
	BestSF_bgd = h_sf_bgr.GetBinCenter(MinBin)

	# Rescale and find \pm 1\sigma errors
	for i in range(1,h_sf_bgr.GetNbinsX()+1):
		h_sf_bgr_rescaled.SetBinContent(i,Minimum-h_sf_bgr.GetBinContent(i))

	LeftLim = h_sf_bgr.GetBinCenter( h_sf_bgr_rescaled.FindFirstBinAbove(-1)-1 )
	RightLim = h_sf_bgr.GetBinCenter( h_sf_bgr_rescaled.FindLastBinAbove(-1)+1 )

	LeftError = BestSF_bgd - LeftLim
	RightError = RightLim - BestSF_bgd

	# Print summary
	print("   ----------\n","   Result fit: \n","   ----------","Background scale factor from sideband fit: ",BestSF_bgd," - ",LeftError," + ",RightError)

	# Plot histogram
	canvas1 = TCanvas("canvas1","Standard Canvas",600,400)
	canvas1.SetLeftMargin(0.175)
	canvas1.SetBottomMargin(0.125)
	canvas1.cd()
	h_sf_bgr_rescaled.SetStats(kFALSE)

	h_sf_bgr_rescaled.SetTitle("Background scale factor")

	h_sf_bgr_rescaled.GetXaxis().SetTitle("scale factor")
	h_sf_bgr_rescaled.GetYaxis().SetTitle("2 log L")
	#h_sf_bgr_rescaled.SetAxisRange(1,-1	,"Y")
	h_sf_bgr_rescaled.Draw("C")

	CenterLine = TLine(BestSF_bgd,0,BestSF_bgd,-31)
	CenterLine.SetLineColor(2)
	LeftLine = TLine(LeftLim,Minimum-h_sf_bgr.GetBinContent( h_sf_bgr_rescaled.FindFirstBinAbove(-1)-1 ),LeftLim,-31)
	LeftLine.SetLineColor(40)
	RightLine = TLine(RightLim,Minimum-h_sf_bgr.GetBinContent( h_sf_bgr_rescaled.FindLastBinAbove(-1)+1 ),RightLim,-31)
	RightLine.SetLineColor(40)

	CenterLine.Draw("same axis")
	LeftLine.Draw("same axis")
	RightLine.Draw("same axis")

	canvas1.Print("Plots/SideBandFit.pdf")


	# Find expected background
	bgr = h_bgr.Integral(h_bgr.FindBin(125.-0.5*7.15),h_bgr.FindBin(125+0.5*7.15))
	sig = GetMassDistribution(125).Integral(h_bgr.FindBin(125.-0.5*7.15),h_bgr.FindBin(125+0.5*7.15))
	obs = h_data.Integral(h_bgr.FindBin(125.-0.5*7.15),h_bgr.FindBin(125+0.5*7.15))
	print("BACKGROUND - without rescaling	: ", bgr)
	print("Best scalefactor: ",BestSF_bgd)
	print("BACKGROUND - with rescaling	: ", BestSF_bgd*bgr," - ",LeftError*bgr,' + ',RightError*bgr)
	print("SIGNAL EVENTS:	",sig)
	print("OBSERVED EVENTS:	",obs)

	return

# ---------------------------------------------------------

def ExpectedSignificance_ToyMC(n_MC,lumi=1.,profile='n'):
	masswindow = 7.15
	signal = GetMassDistribution(125,lumi)
	bgr = GetMassDistribution(1,lumi)

	Delta_bgd=0

	mean_bgd = bgr.Integral(bgr.FindBin(125-0.5*masswindow),bgr.FindBin(125+0.5*masswindow))
	if profile=='y': Delta_bgd = 0.5*(0.06+0.07)*mean_bgd; mean_bgd*=1.104; 

	mean_sig = signal.Integral(signal.FindBin(125-0.5*masswindow),signal.FindBin(125+0.5*masswindow))

	print("----------------------------")
	print("	Background events:	",mean_bgd,"	+/-	",Delta_bgd)
	print(" Signal events:	",mean_sig)
	print("----------------------------")


	sign = 0.

	rand = TRandom3(42)

	count = 0.

	# h_sign = TH1D("h_sign","",50,0,5)

	# Loop over MC cycles
	for i in range(0,n_MC):
		Mean_bgd = rand.Gaus(mean_bgd,Delta_bgd)
		N_bgd = rand.Poisson(mean_bgd)
		N_sig = rand.Poisson(mean_sig)

		# Calculate p-values
		pvalue = IntegratePoissonFromRight(N_bgd,N_bgd+N_sig)
		if(pvalue<=0 or pvalue>=1.): continue

		significance = ROOT.Math.gaussian_quantile_c(pvalue,1)
		sign += significance
		# print(sign,"	",pvalue)
		# h_sign.Fill(significance)
		count += 1.

	sign /= count

	# cvs = TCanvas("cvs","Standard Canvas",600,400)
	# cvs.cd()
	# h_sign.Draw()
	# cvs.Print("Plots/Significante_toy.pdf")


	print("Expected significance after rescaling:	",sign)
	
	return sign

# -------------------------------------------------

def GetTestStatistics( h_mass_data, h_temp_bgd, h_temp_sig ):
	# Compute likelihood
	Loglik_bgr = 0.
	Loglik_sb = 0.

	for i in range(1,h_mass_data.GetNbinsX()+1):
		# \mu = 0 (no signal)
		Loglik_bgr += TMath.Log( TMath.Poisson( h_mass_data.GetBinContent(i),h_temp_bgd.GetBinContent(i) ) ) 
		# \mu = 1 (signal + background)
		Loglik_sb += TMath.Log( TMath.Poisson( h_mass_data.GetBinContent(i),h_temp_sig.GetBinContent(i)+h_temp_bgd.GetBinContent(i) ) )

	# Get likelihood ratio
	X = 2*( Loglik_bgr-Loglik_sb )
	return X 


def GenerateToyDataset( h_mass_temp ):
	h_toy = h_mass_temp.Clone("h_toy")
	h_toy.Reset()

	rand = TRandom3(0)

	# Loop over bins
	for i in range(1,h_mass_temp.GetNbinsX()+1):
		Nbin = rand.Poisson( h_mass_temp.GetBinContent(i) )
		h_toy.SetBinContent( i,Nbin )

	return h_toy

# ---------------------------------------------------

def ToyDataSet( type="bgr" ):
	if type == "bgr":
		h_temp = GetMassDistribution(1)
		h_mass_temp = GenerateToyDataset(h_temp)

	if type == "sb":
		h_bgr = GetMassDistribution(1)
		h_sig = GetMassDistribution(125)

		# Prepare cumulative S+B
		h_temp = h_bgr.Clone()
		for i in range(1,h_temp.GetNbinsX()+1):
			h_temp.SetBinContent( i, h_bgr.GetBinContent(i)+h_sig.GetBinContent(i) )
		h_mass_temp = GenerateToyDataset(h_temp)

	if type == "data":
		h_temp = GetMassDistribution(2)
		h_mass_temp = GenerateToyDataset(h_temp)

	return h_mass_temp

# ----------------------------------------------

def Quantiles(hist):
	frac_1sigma = ROOT.Math.gaussian_cdf(-1.,1.,0.)
	frac_2sigma = ROOT.Math.gaussian_cdf(-2.,1.,0.)
	probs = np.array( [frac_2sigma,frac_1sigma,0.5,1-frac_1sigma, 1-frac_2sigma] )

	Xvalues = np.zeros(5)
	Xvalues_out = np.zeros(5)

	# Extract quantiles
	hist.GetQuantiles(5,Xvalues,probs)

	for i in range(0,5):
		Xvalues_out[i]=Xvalues[i]

	return Xvalues_out



def TestStatisticsDistribution( N_MC ):
	h_bgr = GetMassDistribution(1)
	h_sig = GetMassDistribution(125)
	h_data = GetMassDistribution(2)

	h_test_bgr = TH1D("test_bgr","",200,-30.,30.)
	h_test_sb = TH1D("test_sb","",200,-30.,30.)

	rand = TRandom3()

	for i in range(1,N_MC+1):
		h_toy_bgr = ToyDataSet("bgr")
		h_toy_sb = ToyDataSet("sb")

		h_test_bgr.Fill( GetTestStatistics(h_toy_bgr,h_bgr,h_sig) )
		h_test_sb.Fill( GetTestStatistics(h_toy_sb,h_bgr,h_sig) )

		if i%100==0: print("step ",i," of ",N_MC )
			

	TestStat_data = GetTestStatistics(h_data,h_bgr,h_sig)
	print("--------------------")
	print("Test Statistics value for real data: ",TestStat_data)
	print("--------------------")

	# Computing p-values

	# DISCOVERY - CLb (b-only hypothesis) 
	Nbins = h_test_bgr.GetNbinsX()
	# Median b-only
	Median_bgr = Quantiles(h_test_bgr)[2]
	CLb_bgr = h_test_bgr.Integral(h_test_bgr.FindBin(Median_bgr),Nbins)/h_test_bgr.Integral()
	
	# Median s+b
	Median_sb = Quantiles(h_test_sb)[2]
	CLb_sb = h_test_bgr.Integral(h_test_bgr.FindBin(Median_sb),Nbins)/h_test_bgr.Integral()

	# Data
	CLb_data = h_test_bgr.Integral(h_test_bgr.FindBin(TestStat_data),h_test_bgr.GetNbinsX())/h_test_bgr.Integral()

	# Print a log
	print("==============================")
	print("b-only hypothesis")
	print("-----------------------------------------")
	print("	CLb 		pvalue 		significance 	")
	print("bgr:	",CLb_bgr,"	",1-CLb_bgr,"	",ROOT.Math.gaussian_quantile_c(1-CLb_bgr,1))
	print("S+B:	",CLb_sb,"	",1-CLb_sb,"	",ROOT.Math.gaussian_quantile_c(1-CLb_sb,1))
	print("Data:	",CLb_data,"	",1-CLb_data,"	",ROOT.Math.gaussian_quantile_c(1-CLb_data,1))
	print("-----------------------------------------")


	# CCL_{s+b} (S+B hypothesis)

	CLsb_bgr = h_test_sb.Integral(h_test_bgr.FindBin(Median_bgr),h_test_sb.GetNbinsX())/h_test_sb.Integral()
	
	# # Median s+b
	CLsb_sb = h_test_sb.Integral(h_test_sb.FindBin(Median_sb),h_test_sb.GetNbinsX())/h_test_sb.Integral()

	# # Data
	CLsb_data = h_test_sb.Integral(h_test_sb.FindBin(TestStat_data),h_test_sb.GetNbinsX())/h_test_sb.Integral()

	print("\ns+b hypothesis")
	print("-----------------------------------------")
	print("	CL_sb 		pvalue 		significance 	")
	print("bgr:	",CLsb_bgr,"	",1-CLsb_bgr,"	",ROOT.Math.gaussian_quantile_c(1-CLsb_bgr,1))
	print("S+B:	",CLb_sb,"	",1-CLb_sb,"	",ROOT.Math.gaussian_quantile_c(1-CLsb_sb,1))
	print("Data:	",CLsb_data,"	",1-CLsb_data,"	",ROOT.Math.gaussian_quantile_c(1-CLsb_data,1))
	print("-----------------------------------------")

	cvs = TCanvas("cvs","Standard Canvas",600,400)
	Data_max = h_test_bgr.GetBinContent(h_test_bgr.GetMaximumBin())
	Ymax_plot = Data_max + np.sqrt(Data_max)
	h_test_bgr.SetAxisRange(0.,Ymax_plot,"Y")
	line = TLine(TestStat_data,0,TestStat_data,Ymax_plot)

	cvs.cd()

	h_test_bgr.Draw("L")
	h_test_sb.Draw("same L")
	line.Draw("same")
	cvs.Print("Plots/TestStat.pdf")

# -------------------------------------------------------

def MuFit(Nbins,irebin=1.):
	# Get histrograms
	h_bgr = GetMassDistribution(1)
	h_data = GetMassDistribution(2)
	h_sig = GetMassDistribution(125)

	h_bgr.Rebin(irebin)
	h_data.Rebin(irebin)
	h_sig.Rebin(irebin)

	h_sf = TH2D("scalefactor","title",Nbins,0.5,2.5,Nbins,0.,6.)

	for i in range(1,h_sf.GetNbinsX()+1):
		for j in range(1,h_sf.GetNbinsY()+1):

			sf_bgr = h_sf.GetXaxis().GetBinCenter(i)
			sf_sig = h_sf.GetYaxis().GetBinCenter(j)

			# Loop over bins, compute likelihood
			loglik = 0.

			for iDataBin in range(1,h_data.GetNbinsX()+1):

				m4lepBin = h_data.GetBinCenter(iDataBin)
				NObsBin = h_data.GetBinContent(iDataBin)

				MeanBin = sf_bgr*h_bgr.GetBinContent(iDataBin) + sf_sig*h_sig.GetBinContent(iDataBin)

				if (MeanBin>0): loglik += TMath.Log( TMath.Poisson(NObsBin,MeanBin) )
				
				#print("Alpha = ",sf_bgr,"	Mu = ",sf_sig,"	likelihood = ",loglik)

			h_sf.SetBinContent(i,j,-2.*loglik)

	# Get best SF parameters

	x=ctypes.c_int(0)
	y=ctypes.c_int(0)
	z=ctypes.c_int(0)
	h_sf.GetBinXYZ(h_sf.GetMinimumBin(),x,y,z)
	Minimum = h_sf.GetBinContent(x.value,y.value)
	best_alpha = h_sf.GetXaxis().GetBinCenter(x.value)
	best_mu = h_sf.GetYaxis().GetBinCenter(y.value)
	print(Minimum,best_alpha,best_mu)
	

	# Rescale histogram
	for i in range(1,h_sf.GetNbinsX()+1):
		for j in range(1,h_sf.GetNbinsY()+1):
			h_sf.SetBinContent( i,j, h_sf.GetBinContent(i,j)-Minimum )


	canvas = TCanvas("canvas","Standard Canvas",600,400)
	canvas.cd()
	Min = TMarker(best_alpha,best_mu,29)
	Min.SetMarkerSize(2)
	h_sf.Draw("COLZ")
	Min.Draw()
	canvas.Print("Plots/MuFit.pdf")


# ------------------------------------

def PoissonError(nObs, ErrorType, plot=""):
	if plot=="plot": 
		canvas1 = TCanvas("canvas1","Standard Canvas",600,400)
		canvas1.SetLeftMargin(0.125)
		canvas1.SetBottomMargin(0.125)

	# Prepare histograms
	LambdaMin = 0.
	LambdaMax = nObs + 6*np.sqrt(nObs)
	Nsteps = 10000
	h_likelihood = TH1D( "h_likelihood","lambda likelihood",Nsteps,LambdaMin,LambdaMax)
	h_2loglik = TH1D("h_2loglik","lambda likelihood",Nsteps,LambdaMin,LambdaMax)
	h_pdf_full = TH1D("h_pdf_full","PDF for lambda 		",Nsteps,LambdaMin,LambdaMax)

	# Standard parameters
	Integral_fraction_1sigma = ROOT.Math.gaussian_cdf(-1,1,0)

	# Loop over hypotheses of mean of the 'true' Poisson (lambda)
	for iBin in range(1,Nsteps+1):
		Lambda = h_likelihood.GetBinCenter(iBin)
		# if np.abs(Lambda)<1e-9: Lambda = 1e-9

		# Compute Poisson probability
		PoissonProb = TMath.Poisson(nObs,Lambda)
		LogLikelihood = -2.*TMath.Log( PoissonProb )

		# Fill histogram for Likelihood and PDF
		h_likelihood.Fill(Lambda, PoissonProb)
		h_2loglik.Fill(Lambda, LogLikelihood)
		h_pdf_full.Fill(Lambda, PoissonProb*1)

	# Get characteristic values
	bin_central = h_2loglik.GetMinimumBin()
	LoglikMin = h_2loglik.GetBinContent(bin_central)
	Lambda_central = h_2loglik.GetBinCenter(bin_central)

	# Compute error region

	Lambda_1sigma_low = -1.
	Lambda_1sigma_up = -1.

	# FREQUENTIST 
	if ErrorType=="Classical Central":
		NobsMax = nObs+100

		for ibin in range(h_pdf_full.GetNbinsX()+1): # loop over lambda
			Lambda = h_pdf_full.GetBinCenter(ibin)
			PoissonSum_low = 0.
			PoissonSum_up = 0.

			# lower value
			for i in range(nObs,NobsMax+1):
				PoissonSum_low += TMath.Poisson(i, Lambda)

			if (PoissonSum_low>Integral_fraction_1sigma and Lambda_1sigma_low<0.):
				Lambda_1sigma_low = Lambda

			# upper value
			for i in range(0,nObs+1):
				PoissonSum_up += TMath.Poisson(i,Lambda)

			if (PoissonSum_up<Integral_fraction_1sigma and Lambda_1sigma_up < 0.):
				Lambda_1sigma_up = Lambda

	if ErrorType=="Likelihood Ratio":
		for i in range(1,h_2loglik.GetNbinsX()+1):
			if (h_2loglik.GetBinCenter(i)<LoglikMin and h_2loglik.GetBinContent(i)-LoglikMin>=1.):  Lambda_1sigma_low=h_2loglik.GetBinCenter(i)
			if (h_2loglik.GetBinCenter(i)>LoglikMin and h_2loglik.GetBinContent(i)-LoglikMin<=1.):  Lambda_1sigma_up=h_2loglik.GetBinCenter(i)


	if ErrorType=="Bayes Central": 
		# Work on likelihood as PDF
		Integral = h_likelihood.Integral()
		for i in range(1,h_likelihood.GetNbinsX()+1):
			h_likelihood.SetBinContent(i,h_likelihood.GetBinContent(i)/Integral)

		# Sum over bins until reached CL
		Sum_low = 0.
		Sum_up = 0.
		for i in range(1,bin_central+1):
			if Sum_low<=Integral_fraction_1sigma:
				Sum_low += h_likelihood.GetBinContent(i)
				Lambda_1sigma_low = h_likelihood.GetBinCenter(i)
		for i in range(h_likelihood.GetNbinsX(),bin_central,-1):
			if Sum_up<=Integral_fraction_1sigma:
				Sum_up += h_likelihood.GetBinContent(i)
				Lambda_1sigma_up = h_likelihood.GetBinCenter(i)


	if ErrorType=="Bayes HDI":
		# Work on likelihood as PDF
		Integral = h_likelihood.Integral()
		for i in range(1,h_likelihood.GetNbinsX()+1):
			h_likelihood.SetBinContent(i,h_likelihood.GetBinContent(i)/Integral)

		Area = 0.	
		RightIndex=bin_central
		while Area<=(1-2*Integral_fraction_1sigma):
			RightIndex+=1
			# find corresponding bin on left side
			for i in range(1,bin_central):
				if h_likelihood.GetBinContent(i)<=h_likelihood.GetBinContent(RightIndex): 
					LeftIndex=i
			Area = h_likelihood.Integral( LeftIndex,RightIndex )

		Lambda_1sigma_low = h_likelihood.GetBinCenter(LeftIndex)
		Lambda_1sigma_up = h_likelihood.GettingBinCenter(RightIndex)

	return [Lambda_1sigma_low, Lambda_1sigma_up]





#	if plot=="plot":
#		NobsMax = nObs+25

#		h1 = TH1D("h1", "poisson probability (sec)	",NobsMax+1,-0.5,NobsMax+0.5)
#		h_up = TH1D("h_up", "poisson probability (n >=nobs)", NobsMax+1,-0.5,NobsMax+0.5)
#		h_low = TH1D("h_low", "poisson probability (n >=nobs)", NobsMax+1,-0.5,NobsMax)













# =========================================================================
		






	
	
	
	

	

	
