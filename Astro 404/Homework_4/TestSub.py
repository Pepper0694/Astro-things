import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches


perc3, a3 = np.loadtxt('145582795929122.dat', unpack = True, usecols=[0,1])
#absorb3 = ((3**8)/(a3*(10**9)))
plt.plot(perc3, ((3*10**8)/(a3*(10**9))), color= 'deeppink', linewidth=2)
plt.show()
