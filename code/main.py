import higgs

# MAIN - RUNNING THE FUNCTIONS

#higgs.MassPlot(20)

# lumi = np.linspace(1,6,20)
# sign = np.zeros(lumi.size)
# for i in range(lumi.size):
# 	sign[i] = SignificanceOpt(lumi[i])
# 	# if  sign[i]>=5.0 : break 

# plt.plot(lumi, sign,'.')
# plt.plot([1,6],[5.0,5.0],linewidth=1)
# plt.grid(linestyle='--',linewidth=1)
# plt.show()

# higgs.SideBandFit(10)

# higgs.ExpectedSignificance_ToyMC(6.42,0.,5.96,1e6,"bootstrap")

#h_data = higgs.GetMassDistribution(2)
#h_sig = higgs.GetMassDistribution(125)
#h_bgd = higgs.GetMassDistribution(1)

#print("TestStatistics = ",higgs.GetTestStatistics(h_data, h_bgd, h_sig))

# higgs.TestStatisticsDistribution( 10000 )

higgs.MuFit(50,10)