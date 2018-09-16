from gpfit.fit import fit
from numpy import logspace, log, random, vstack
# fixed initial guess for fitting
random.seed(33404)
Re = logspace(3, 8, 1001)
relRough = logspace(-5, -1, 1001)
f = 0.25/log(relRough/3.7 + 5.74/Re**0.9)**2
x = log(Re)
y = log(relRough)
z = log(f)
xData = vstack((x, y))
K = 3
cMA, errorMA = fit(xData, z, K, "MA")
cSMA, errorSMA = fit(xData, z, K, "SMA")
cISMA, errorISMA = fit(xData, z, K, "ISMA")
print "MA RMS Error: %.5g" % errorMA
print "SMA RMS Error: %.5g" % errorSMA
print "ISMA RMS Error: %.5g" % errorISMA
