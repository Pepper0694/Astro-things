import numpy as np
import matplotlib.pyplot as plt


mp= 1.6726*(10**-27)
rho= 10**(-6)
me= 9.11*(10**-31)
k= 1.38*(10**-23)
h= 6.626*(10**-34)
chi= 2.179*(10**-18)

b = (mp/rho)*(((2*np.pi*me*k)/(h**2))**(3/2))


def ion_fraction(t):
    a= (b*(np.exp(-chi/(k*t)))*(t**(3/2)))
    x = -a + np.sqrt(((a)**2) + (4*a))
    return x/2


t = np.arange(5000., 25000.)

plt.figure()
plt.xlabel("Temp (K)")
plt.ylabel("Ion Fraction (N(II)/N(H)")
plt.title("Hydroden Ionization")
plt.plot(t, ion_fraction(t), 'b')
plt.show()
