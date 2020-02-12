
import numpy as np
import multiprocessing as mp
import time
import ROOT
import ctypes
from ROOT import gStyle, kRed, TGraph, TMultiGraph, TPad, TBox, TArrow, TCanvas, TH1D, TH2D, TLine, TMath, TRandom3, TROOT, TLatex, TFile, TLegend, TLegendEntry, gROOT, gDirectory, kTRUE, kFALSE, TMarker


# ------- FUNCTIONS ---------

def GetMassDistribution( Itype=1, scalefactor=1.00 ):
	dir = gDirectory
	file = TFile("Data/Histograms_fake.root","READ")
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
	for i in range(1,Nbins+1):
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

def IntegratePoisson(mu, Nobs):
	integral = 0.
	for i in range(0,int(Nobs)):
		integral += TMath.Poisson(i, mu)
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

	print("Expected significance after rescaling:	",sign)



def ExpectedSignificance_ToyMC(mean_bgd, Delta_bgd, mean_sig, n_MC):
	gROOT.Clear()
	gROOT.Delete()

	# Define count histograms
	h_Nbgr = TH1D("h_Nbgr","Background events",500,-0.5,499.5)

	# Initialize seed
	rand = TRandom3()
	
	# Generate toy datasets
	for i in range(1,n_MC+1):
		mean_bgr = rand.Gaus(mean_bgd, Delta_bgd)
		mean_sb = rand.Gaus(mean_bgd, Delta_bgd)+mean_sig
		
		h_Nbgr.Fill(rand.Poisson(mean_bgr))

	# Calculate p-values
	pvalue = h_Nbgr.Integral(h_Nbgr.FindBin(mean_bgd+mean_sig),h_Nbgr.GetNbinsX())/h_Nbgr.Integral()
	significance = ROOT.Math.gaussian_quantile_c(pvalue,1)
	#print(pvalue,"	",significance)

	print('Expected significance after rescaling:	',significance)
	
	return significance

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


def FillTestHistogram(templ,N_MC,signal=1.0):
	h_bgr = GetMassDistribution(1)
	h_sig = GetMassDistribution(125)

	print(signal)

	h_test = TH1D("test_"+templ,"",200,-30.,30.)
	for i in range(1,N_MC+1):
		h_toy = ToyDataSet(templ)

		h_test.Fill( GetTestStatistics(h_toy, h_bgr, h_sig) )

		if i%1000==0: print("step ",i," of ",N_MC )

	return h_test




def TestStatisticsDistribution( N_MC, signal=1.0 ):

	pool = mp.Pool(mp.cpu_count())

	templ = ["bgr","sb"]

	[h_test_bgr, h_test_sb] = [pool.apply( FillTestHistogram, args=(i, N_MC, signal) ) for i in  templ]

	pool.close()

	# saving MC dataset
	dir = gDirectory
	dir.cd()
	myfile = TFile("toy_"+str(signal)+".root","RECREATE")
	h_test_bgr.Write()
	h_test_sb.Write()
	myfile.Close()

	return			

def Pvalues(signal=1.0):

	h_data = GetMassDistribution(2)
	h_bgr = GetMassDistribution(1)
	h_sig = GetMassDistribution(125)


	dir = ROOT.gDirectory
	file = ROOT.TFile("toy_"+str(signal)+".root","READ")
	dir.cd()

	h_test_bgr = file.Get("test_bgr").Clone("h_test_bgr")
	h_test_sb = file.Get("test_sb").Clone("h_test_sb")



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
	print("	CL_sb 		pvalue 		CL_s 	")
	print("bgr:	",CLsb_bgr,"	",1-CLsb_bgr,"	",CLsb_bgr/CLb_bgr )
	print("S+B:	",CLsb_sb,"	",1-CLsb_sb,"	",CLsb_sb/CLb_sb )
	print("Data:	",CLsb_data,"	",1-CLsb_data,"	",CLsb_data/CLb_data )
	print("-----------------------------------------")

	if signal!=1.: 
		cvs = TCanvas("cvs","Standard Canvas",600,400)
		Data_max = h_test_bgr.GetBinContent(h_test_bgr.GetMaximumBin())
		Ymax_plot = Data_max + np.sqrt(Data_max)
		h_test_bgr.SetAxisRange(0.,Ymax_plot,"Y")
		line = TLine(TestStat_data,0,TestStat_data,Ymax_plot)

		cvs.cd()
		h_test_bgr.SetTitle("TEST STATISTIC DISTRIBUTION")
		h_test_bgr.SetLineColor(9)
		h_test_sb.SetLineColor(8)
		line.SetLineColor(2)
		h_test_bgr.GetXaxis().SetTitle("test statistic values")
		h_test_bgr.GetYaxis().SetTitle("# toy experiments")
		h_test_bgr.SetStats(kFALSE)

		h_test_bgr.Draw("L")
		h_test_sb.Draw("same L")
		line.Draw("same")
		leg1 = TLegend(0.70,0.70,0.90,0.85)
		leg1.SetBorderSize(1); leg1.SetFillColor(0);
		leg1a = leg1.AddEntry(h_test_bgr, "b- only", "l"); leg1a.SetTextSize(0.04);
		leg1b = leg1.AddEntry(h_test_sb, "s+b", "l"); leg1b.SetTextSize(0.04);
		leg1.Draw()
		cvs.Print("Plots/TestStat"+str(signal)+".pdf")

	file.Close()

	return CLsb_data, CLsb_data/CLb_data

# -------------------------------------------------------

def MuFit(Nbins,irebin=1.):
	# Get histrograms
	h_bgr = GetMassDistribution(1)
	h_data = GetMassDistribution(2)
	h_sig = GetMassDistribution(125)

	h_bgr.Rebin(irebin)
	h_data.Rebin(irebin)
	h_sig.Rebin(irebin)

	h_sf = TH2D("scalefactor","title",Nbins,0.5,2.,Nbins,0.,5.)

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
	h_sf.SetStats(kFALSE)
	h_sf.SetTitle("HIGGS PRODUCTION CROSS SECTION - PARAMETERS")
	h_sf.GetXaxis().SetTitle(r"#alpha")
	h_sf.GetYaxis().SetTitle(r"#mu")
	h_sf.Draw("COLZ")

	h_sigma = h_sf.Clone("h_sigma")
	h_sigma.Reset()
	for i in range(1,h_sigma.GetNbinsX()+1):
		for j in range(1,h_sigma.GetNbinsY()+1):
			if( h_sf.GetBinContent(i,j)<=1.): h_sigma.SetBinContent(i,j, 1.)
	h_sigma.SetMarkerColorAlpha(kRed,0.40)

	h_sigma.SetMarkerSize(10)
	h_sigma.Draw("same L")
	Min.Draw()
	leg1 = TLegend(0.65,0.85,0.85,0.75)
	leg1.SetBorderSize(1); leg1.SetFillColor(0);
	leg1a = leg1.AddEntry(Min, r"optimal (#alpha,#mu)", "p"); leg1a.SetTextSize(0.04);
	leg1.Draw()
	canvas.Print("Plots/MuFit.pdf")


# ------------------------------------

def PoissonError(nObs, ErrorType, plot=""):

	# Prepare histograms
	LambdaMin = 0.
	LambdaMax = nObs + 6*np.sqrt(nObs)
	Nsteps = 1000

	h_likelihood = TH1D( "h_likelihood","",Nsteps,LambdaMin,LambdaMax )
	h_2loglik = TH1D( "h_2loglik","",Nsteps,LambdaMin,LambdaMax )
	h_pdf_full = TH1D( "h_pdf_full","",Nsteps,LambdaMin,LambdaMax )

	IntFraction =ROOT.Math.gaussian_cdf(-1,1,0)

	# loop over possible Lambda values
	for iBin in  range(1,Nsteps+1):
		Lambda=h_likelihood.GetBinCenter(iBin)

		PoissonProb = TMath.Poisson(nObs,Lambda)
		LogLikelihood = -2.*TMath.Log(PoissonProb)

		h_likelihood.Fill( Lambda,PoissonProb )
		h_2loglik.Fill( Lambda, LogLikelihood )
		h_pdf_full.Fill( Lambda, PoissonProb*1. )

	# get characteristic values
	bin_central = h_2loglik.GetMinimumBin()
	LoglikMin = h_2loglik.GetBinContent(bin_central)
	Lambda_central = h_2loglik.GetBinCenter(bin_central)

	LambdaLow = -1.
	LambdaUp = -1.


	if ErrorType=="ClassicalCentral": # Frequentist
		NobsMax = nObs+100

		for iBin in range(1,h_pdf_full.GetNbinsX()+1):
			Lambda = h_pdf_full.GetBinCenter(iBin)
			PoissonSumLow = 0.
			PoissonSumUp = 0.

			for i in range(nObs,NobsMax+1):
				PoissonSumLow += TMath.Poisson(i,Lambda)
			
			if( PoissonSumLow>IntFraction and LambdaLow<0 ):
				LambdaLow = Lambda

			for i in range(0,nObs+1):
				PoissonSumUp += TMath.Poisson(i,Lambda)

			if( PoissonSumUp<IntFraction and LambdaUp<0 ):
				LambdaUp = Lambda

		cvs = TCanvas("Standard Canvas","",600,600)

	
		cvs.Divide(1,2)

		cvs.cd(1)

		hLow = TH1D("hLow","CLASSICAL CENTRAL",1000,0,17)
		for i in range(1,hLow.GetNbinsX()+1):
			hLow.SetBinContent(i, TMath.Poisson( hLow.GetBinCenter(i),LambdaLow ))
		hLow.SetLineColor(1)
		hLow.GetXaxis().SetTitle(r"#mu")

		h_conf = hLow.Clone("h_conf")
		for i in range(1,h_conf.GetNbinsX()+1):
			if(h_conf.GetBinCenter(i)<=nObs): h_conf.SetBinContent(i, 0)
		h_conf.SetFillColorAlpha(9,0.5)
		h_conf.SetLineWidth(1)
		h_conf.SetLineColor(1)


		LowLine = TLine(LambdaLow,0,LambdaLow,hLow.GetBinContent(hLow.FindBin(LambdaLow)))

		hLow.SetStats(kFALSE)


		hLow.Draw("l")
		h_conf.Draw("l same")
		LowLine.Draw()

		text = TLatex()
		text.SetTextSize(0.5)
		text.SetTextColor(1)
		text.DrawLatex(0.5,0.8, "LOL")

		cvs.cd(2)

		hUp = TH1D("hUp","",1000,0,17)
		for i in range(1,hUp.GetNbinsX()+1):
			hUp.SetBinContent(i, TMath.Poisson( hLow.GetBinCenter(i),LambdaUp ))
		hUp.SetLineColor(1)
		hUp.GetXaxis().SetTitle(r"#mu")

		h_conf1 = hUp.Clone("h_conf")
		for i in range(1,h_conf.GetNbinsX()+1):
			if(h_conf.GetBinCenter(i)>=nObs): h_conf1.SetBinContent(i, 0)
		h_conf1.SetFillColorAlpha(9,0.5)
		h_conf1.SetLineWidth(1)
		h_conf1.SetLineColor(1)

		UpLine = TLine(LambdaUp,0,LambdaUp,hUp.GetBinContent(hUp.FindBin(LambdaUp)))

		hUp.SetStats(kFALSE)


		hUp.Draw("l")
		h_conf1.Draw("l same")
		UpLine.Draw()

		cvs.cd()

		leg1 = TLegend(0.35,0.70,0.90,0.85)
		leg1.SetBorderSize(1); leg1.SetFillColor(0);
		h_conf.SetMarkerColor(0)
		leg1a = leg1.AddEntry(h_conf, r"(#mu_{low}, #mu_{up}) = ("+str(round(LambdaLow,2))+","+str(round(LambdaUp,2))+")","p"); leg1a.SetTextSize(0.04);
		leg1.Draw()

		cvs.Print("Plots/Errors_"+ErrorType+".eps")



	if ErrorType=="LikelihoodRatio":
		for i in range(1,h_2loglik.GetNbinsX()+1):
			if (h_2loglik.GetBinCenter(i)<LoglikMin and h_2loglik.GetBinContent(i)-LoglikMin>=1.):  LambdaLow=h_2loglik.GetBinCenter(i)
			if (h_2loglik.GetBinCenter(i)>LoglikMin and h_2loglik.GetBinContent(i)-LoglikMin<=1.):  LambdaUp=h_2loglik.GetBinCenter(i)

		cvs = TCanvas("Standard Canvas","",600,400)
		cvs.cd()

		LowLine = TLine(LambdaLow,2,LambdaLow,h_2loglik.GetBinContent(h_2loglik.FindBin(LambdaLow))); LowLine.SetLineWidth(1)
		UpLine = TLine(LambdaUp,2,LambdaUp,h_2loglik.GetBinContent(h_2loglik.FindBin(LambdaUp))); UpLine.SetLineWidth(1)
		ObsLine = TLine(nObs,2,nObs,h_2loglik.GetBinContent(h_2loglik.FindBin(nObs))); ObsLine.SetLineWidth(1); ObsLine.SetLineColor(2); ObsLine.SetLineStyle(7)

		h_2loglik.SetFillColor(0)
		h_2loglik.SetStats(kFALSE)
		h_2loglik.SetTitle("LIKELIHOOD RATIO")
		h_2loglik.SetAxisRange(0.8,8,"X")
		h_2loglik.SetAxisRange(2,7,"Y")
		h_2loglik.Draw("hist lp")
		h_2loglik.GetXaxis().SetTitle(r"#mu")

		h_conf = h_2loglik.Clone("h_conf")
		for i in range(1,h_conf.GetNbinsX()+1):
			if(h_conf.GetBinCenter(i)>LambdaLow and h_conf.GetBinCenter(i)<LambdaUp): h_conf.SetBinContent(i,0)

		h_conf.SetFillColor(9)
		h_conf.Draw("hist same")

		leg1 = TLegend(0.55,0.70,0.90,0.85)
		leg1.SetBorderSize(1); leg1.SetFillColor(0);
		h_conf.SetMarkerColor(0)
		leg1a = leg1.AddEntry(h_conf, r"(#mu_{low}, #mu_{up}) = ("+str(round(LambdaLow,2))+","+str(round(LambdaUp,2))+")","p"); leg1a.SetTextSize(0.04);
		leg1.Draw()

		LowLine.Draw()
		UpLine.Draw()
		ObsLine.Draw()
		

		cvs.Print("Plots/Errors_"+ErrorType+".eps")




	if ErrorType=="BayesCentral": 
		# Work on likelihood as PDF
		Integral = h_likelihood.Integral()
		for i in range(1,h_likelihood.GetNbinsX()+1):
			h_likelihood.SetBinContent(i,h_likelihood.GetBinContent(i)/Integral)

		# Sum over bins until reached CL
		Sum_low = 0.
		Sum_up = 0.
		for i in range(1,bin_central+1):
			if Sum_low<=IntFraction:
				Sum_low += h_likelihood.GetBinContent(i)
				LambdaLow = h_likelihood.GetBinCenter(i)
		for i in range(h_likelihood.GetNbinsX(),bin_central,-1):
			if Sum_up<=IntFraction:
				Sum_up += h_likelihood.GetBinContent(i)
				LambdaUp = h_likelihood.GetBinCenter(i)

		cvs = TCanvas("Standard Canvas","",600,400)
		cvs.cd()

		LowLine = TLine(LambdaLow,0,LambdaLow,h_likelihood.GetBinContent(h_likelihood.FindBin(LambdaLow))); LowLine.SetLineWidth(1)
		UpLine = TLine(LambdaUp,0,LambdaUp,h_likelihood.GetBinContent(h_likelihood.FindBin(LambdaUp))); UpLine.SetLineWidth(1)
		ObsLine = TLine(nObs,0,nObs,h_likelihood.GetBinContent(h_2loglik.FindBin(nObs))); ObsLine.SetLineWidth(1); ObsLine.SetLineColor(2); ObsLine.SetLineStyle(7)


		h_likelihood.SetFillColor(0)
		h_likelihood.SetStats(kFALSE)
		h_likelihood.SetTitle("BAYES CENTRAL")
		#h_likelihood.SetAxisRange(0.8,8,"X")
		#h_likelihood.SetAxisRange(2,7,"Y")
		h_likelihood.Draw("hist lp")
		h_likelihood.GetXaxis().SetTitle(r"#mu")

		h_conf = h_likelihood.Clone("h_conf")
		for i in range(1,h_conf.GetNbinsX()+1):
			if(h_conf.GetBinCenter(i)>LambdaLow and h_conf.GetBinCenter(i)<LambdaUp): h_conf.SetBinContent(i,0)

		h_conf.SetFillColor(9)
		h_conf.Draw("hist same")

		leg1 = TLegend(0.55,0.70,0.90,0.85)
		leg1.SetBorderSize(1); leg1.SetFillColor(0);
		h_conf.SetMarkerColor(0)
		leg1a = leg1.AddEntry(h_conf, r"(#mu_{low}, #mu_{up}) = ("+str(round(LambdaLow,2))+","+str(round(LambdaUp,2))+")","p"); leg1a.SetTextSize(0.04);
		leg1.Draw()

		LowLine.Draw()
		UpLine.Draw()
		ObsLine.Draw()
		

		cvs.Print("Plots/Errors_"+ErrorType+".eps")


	if ErrorType=="BayesHDI":
		# Work on likelihood as PDF
		Integral = h_likelihood.Integral()
		for i in range(1,h_likelihood.GetNbinsX()+1):
			h_likelihood.SetBinContent(i,h_likelihood.GetBinContent(i)/Integral)

		Area = 0.	
		RightIndex=bin_central
		while Area<=(1-2*IntFraction):
			RightIndex+=1
			# find corresponding bin on left side
			for i in range(1,bin_central):
				if h_likelihood.GetBinContent(i)<=h_likelihood.GetBinContent(RightIndex): 
					LeftIndex=i
			Area = h_likelihood.Integral( LeftIndex,RightIndex )

		LambdaLow = h_likelihood.GetBinCenter(LeftIndex)
		LambdaUp = h_likelihood.GetBinCenter(RightIndex)

		cvs = TCanvas("Standard Canvas","",600,400)
		cvs.cd()

		LowLine = TLine(LambdaLow,0,LambdaLow,h_likelihood.GetBinContent(h_likelihood.FindBin(LambdaLow))); LowLine.SetLineWidth(1)
		UpLine = TLine(LambdaUp,0,LambdaUp,h_likelihood.GetBinContent(h_likelihood.FindBin(LambdaUp))); UpLine.SetLineWidth(1)
		ObsLine = TLine(nObs,0,nObs,h_likelihood.GetBinContent(h_2loglik.FindBin(nObs))); ObsLine.SetLineWidth(1); ObsLine.SetLineColor(2); ObsLine.SetLineStyle(7)

		h_likelihood.SetFillColor(0)
		h_likelihood.SetStats(kFALSE)
		h_likelihood.SetTitle("BAYES SHORTEST")
		h_likelihood.Draw("hist lp")
		h_likelihood.GetXaxis().SetTitle(r"#mu")

		h_conf = h_likelihood.Clone("h_conf")
		for i in range(1,h_conf.GetNbinsX()+1):
			if(h_conf.GetBinCenter(i)>LambdaLow and h_conf.GetBinCenter(i)<LambdaUp): h_conf.SetBinContent(i,0)

		h_conf.SetFillColor(9)
		h_conf.Draw("hist same")

		leg1 = TLegend(0.55,0.70,0.90,0.85)
		leg1.SetBorderSize(1); leg1.SetFillColor(0);
		h_conf.SetMarkerColor(0)
		leg1a = leg1.AddEntry(h_conf, r"(#mu_{low}, #mu_{up}) = ("+str(round(LambdaLow,2))+","+str(round(LambdaUp,2))+")","p"); leg1a.SetTextSize(0.04);
		leg1.Draw()

		LowLine.Draw()
		UpLine.Draw()
		ObsLine.Draw()
		

		cvs.Print("Plots/Errors_"+ErrorType+".eps")


	return [LambdaLow, LambdaUp]









	return












# =========================================================================
		






	
	
	
	

	

	
