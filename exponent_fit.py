import numpy as np
import math
import matplotlib
import matplotlib.pyplot as plt
import scipy.special

#x = np.linspace(0, 6)
x = np.arange(0,6,0.1,dtype=np.float)

scale_factor = 0.6
k = 0.75
out = 1.0 - np.exp(-(x/scale_factor)**k)
plt.plot(x, out, 'r-', linewidth=2)

weibull_mean=scale_factor*scipy.special.gamma(1+1./k)
print 'Weibull Mean: ', weibull_mean

out = 1.0 - np.exp(-(x/weibull_mean))
plt.plot(x, out, 'b--', linewidth=2)

plt.show()
