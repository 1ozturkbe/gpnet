from gpfit.fit import fit
import numpy as np
from numpy import logspace, log, random, vstack, hstack
# fixed initial guess for fitting
Re = list(logspace(3, 8, 1001))
Re = hstack(Re)
Re = Re.astype(np.float)

relRough = list(logspace(-5, -1, 1001))
relRough = hstack(relRough)
relRough = relRough.astype(np.float)

f = list(0.25/log(relRough/3.7 + 5.74/Re**0.9)**2)
f = hstack(f)
f = f.astype(np.float)

x = log(Re)
y = log(relRough)
z = log(f)

xData = [Re, relRough]
xData = np.log(xData)

K = 3
cMA, errorMA = fit(xData, z, K, "MA")
cSMA, errorSMA = fit(xData, z, K, "SMA")
cISMA, errorISMA = fit(xData, z, K, "ISMA")
print "MA RMS Error: %.5g" % errorMA
print "SMA RMS Error: %.5g" % errorSMA
print "ISMA RMS Error: %.5g" % errorISMA
