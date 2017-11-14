import numpy as np

Msun = 1.9885*(10**30)
Lsun = 3.827*(10**27)
c = 3*(10**8)

Llower = 10**(-4.3)
Lupper = 10**(6.006)

tlow = (.007*.072*Msun*(c**2))/(Llower*Lsun)
tupper = (.1*.007*85*Msun*(c**2))/(Llower*Lsun)


print(tlow)
print(tupper)
