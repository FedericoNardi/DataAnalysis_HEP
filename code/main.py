import higgs

# MAIN - RUNNING THE FUNCTIONS

higgs.MassPlot(20)

# lumi = np.linspace(1,6,20)
# sign = np.zeros(lumi.size)
# for i in range(lumi.size):
# 	sign[i] = SignificanceOpt(lumi[i])
# 	# if  sign[i]>=5.0 : break 

# plt.plot(lumi, sign,'.')
# plt.plot([1,6],[5.0,5.0],linewidth=1)
# plt.grid(linestyle='--',linewidth=1)
# plt.show()

# SideBandFit()

# higgs.ExpectedSignificance_ToyMC(6.42,0.,5.96,1e6,"bootstrap")