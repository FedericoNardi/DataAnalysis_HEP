import higgs
import numpy as np
import matplotlib.pyplot as plt

# MAIN - RUNNING THE FUNCTIONS

#higgs.MassPlot(20)

# lumi = np.linspace(1,6,20)
# sign = np.zeros(lumi.size)
# for i in range(lumi.size):
# 	sign[i] = higgs.SignificanceOpt(lumi[i])
#  	# if  sign[i]>=5.0 : break 

# plt.plot(lumi, sign,'.--',linewidth=1,markersize=8)
# plt.plot([1,6],[5.0,5.0],linewidth=1)
# plt.grid(linestyle='--',linewidth=1)
# plt.savefig("Plots/Significance_lumi.pdf")
# plt.show()

# higgs.SideBandFit(10)

higgs.ExpectedSignificance_ToyMC(7.10,0.41,5.96,1000000)

#h_data = higgs.GetMassDistribution(2)
#h_sig = higgs.GetMassDistribution(125)
#h_bgd = higgs.GetMassDistribution(1)

#print("TestStatistics = ",higgs.GetTestStatistics(h_data, h_bgd, h_sig))

# higgs.TestStatisticsDistribution( 10000 )

# higgs.MuFit(500,10)

# print(higgs.PoissonError(3,"Bayes HDI"))
