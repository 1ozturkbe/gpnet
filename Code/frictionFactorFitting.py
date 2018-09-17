from gpfit.fit import fit
import numpy as np
from numpy import logspace, log, meshgrid
# fixed initial guess for fitting
Re = logspace(3, 8, 50)
relRough = logspace(-5, -1, 50)
Re, relRough = meshgrid(Re, relRough)
Re = Re.flatten()
relRough = relRough.flatten()

f = 0.25/log(relRough/3.7 + 5.74/Re**0.9)**2 

x = [Re, relRough]
x = np.log(x)
y = np.log(f)
K = 3
cMA, errorMA = fit(x, y, K, "MA")
cSMA, errorSMA = fit(x, y, K, "SMA")
cISMA, errorISMA = fit(x, y, K, "ISMA")
print "MA RMS Error: %.5g" % errorMA
print "SMA RMS Error: %.5g" % errorSMA
print "ISMA RMS Error: %.5g" % errorISMA
