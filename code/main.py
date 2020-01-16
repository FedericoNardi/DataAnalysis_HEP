import higgs
import numpy as np
import matplotlib.pyplot as plt

# MAIN - RUNNING THE FUNCTIONS

# higgs.MassPlot(20)

# Plotting expected significance as a function of mass window 
# higgs.SignificanceOpt(1.00, "y")

def IncreaseLumi(method="",npoints=50):
	lumi = np.linspace(1,6.5,npoints)
	sign = np.zeros(lumi.size)
	for i in range(lumi.size):
		if method=="toy":
			sign[i] = higgs.ExpectedSignificance_ToyMC(1000000,lumi=lumi[i],profile='y')
		else:
			sign[i] = higgs.SignificanceOpt(lumi[i],speak='mute')
		if (sign[i]>=5.0): print("	Significance:	",sign[i],"		Scalefactor:	",lumi[i])
	return sign


#sign= IncreaseLumi(npoints=50)
#plt.plot(lumi, sign,'-',linewidth=1,markersize=0,marker='x')
#plt.plot([1,6],[5.0,5.0],'--',linewidth=1)
#plt.grid(linestyle='--',linewidth=1)
#plt.title(r"SIGNIFICANCE OPTIMIZATION")
#plt.xlabel(r"luminosity factor")
#plt.ylabel(r"significance $Z$ $[ \sigma ]$")
#plt.savefig("Plots/Significance_lumi.pdf")
#plt.show()

# higgs.SideBandFit()

#print("Unscaled:	Nbgr = 4.64, Nsig = 5.41")
#higgs.ExpectedSignificance_ToyMC(1000000)

#print("Unscaled:	Nbgr = 5.12 +/- 0.30 , Nsig = 5.41")
#higgs.ExpectedSignificance_ToyMC(1000000,profile='y')

# sign_toy = IncreaseLumi(method="toy", npoints=7)
# sign = IncreaseLumi(npoints=7)
# print(sign_toy)
# plt.plot(np.linspace(1,6.5,7), sign_toy,'-',linewidth=1,markersize=7,marker='x')
# plt.plot(np.linspace(1,6.5,7), sign,'-',linewidth=1,markersize=0,marker='x')
# plt.plot([1,6],[5.0,5.0],'--',linewidth=1)
# plt.grid(linestyle='--',linewidth=1)
# plt.title(r"SIGNIFICANCE OPTIMIZATION")
# plt.xlabel(r"luminosity factor")
# plt.ylabel(r"significance $Z$ $[ \sigma ]$")
# ax=plt.gca()
# ax.legend(["unscaled","scaled"])
# plt.savefig("Plots/Significance_lumi_scale.pdf")
# plt.show()


h_data = higgs.GetMassDistribution(2)
h_sig = higgs.GetMassDistribution(125)
h_bgd = higgs.GetMassDistribution(1)

print("TestStatistics = ",higgs.GetTestStatistics(h_data, h_bgd, h_sig))

higgs.TestStatisticsDistribution( 10000 )

# higgs.MuFit(500,10)

# print(higgs.PoissonError(3,"Bayes HDI"))