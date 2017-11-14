import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches

lmin, lmax = 10**-12, 10**3
fig=plt.figure()

#Top Subplot Wave
sp1 = fig.add_subplot(211,axisbg='grey')
sp1.set_xlabel("Wavelength (m)")
sp1.set_xticks(np.arange(lmax, lmin, 100))
sp1.set_xlim(lmax, lmin)
sp1.set_xscale('log')
sp1.yaxis.set_visible(False)

#Bottom Subplot Atms Absorb.
sp2 = fig.add_subplot(212, axisbg='gray')
sp2.set_title("Atmospheric Absorbtion")
sp2.set_xlabel("Wavelength (m)")
sp2.set_ylabel("Transparency")
sp2.set_xticks(np.arange(lmax, lmin, 100))
sp2.set_xlim(lmax, lmin)
sp2.set_ylim(0,1)
sp2.set_xscale('log')
sp2.yaxis.set_visible(True)

#Atmos Absorb. graph
    #NIR
a1, perc1 = np.loadtxt('mktrans_zm_10_10.dat', unpack = True)
absorb1 = a1*10**-6
sp2.plot(absorb1,perc1, color = 'deeppink', linewidth=.05)
    #MIR
a2, perc2 = np.loadtxt('mktrans_nq_10_10.dat', unpack = True)
absorb2 = a2*10**-6
sp2.plot(absorb2,perc2, color= 'deeppink', linewidth=.05)
    #Submm
a3, perc3 = np.loadtxt('145582795929122.dat', unpack = True, usecols=[0,1])
#absorb3 = ((3**8)/(a3*(10**9)))
sp2.plot(a3*(3*10**-3), 1-perc3, color= 'deeppink', linewidth=.5)
    #Optical
a4, perc4 = np.loadtxt('SCUBA_450.txt', unpack = True, usecols=[0,1])
sp2.plot((10**-9)*a4, perc4, color= 'blue', linewidth=1)

    #Text Stuff for Opacity
sp2.add_patch(patches.Rectangle((.9*10**1,00 ), 10**3, 1, color="black"))
sp2.text(9*10**1, 0.5, 'Opaque', style='oblique',
        verticalalignment='center', horizontalalignment='center',
        color='deeppink', fontsize=12)

sp2.add_patch(patches.Rectangle((.5*10**-4,00 ), 1.8*10**-3, 1, color="black"))
sp2.text(3*10**-4, 0.5, 'Opaque', style='oblique',
        verticalalignment='center', horizontalalignment='center',
        color='deeppink', fontsize=12)

sp2.add_patch(patches.Rectangle((.5*10**-12,00 ), 5*10**-8, 1, color="black"))
sp2.text(5*10**-10, 0.5, 'Opaque', style='oblique',
        verticalalignment='center', horizontalalignment='center',
        color='deeppink', fontsize=12)





#Wavelength Labels

#Radiowave

sp1.add_patch(patches.Rectangle((10**-1,0.75 ), 10**3, .25, color="black"))
sp1.text(10**1, 0.825, 'Radiowave', style='oblique',
        verticalalignment='center', horizontalalignment='center',
        color='white', fontsize=12)

#Microwave

sp1.add_patch(patches.Rectangle((10**-5,0.5 ), 10**-1, .25, color="orange"))
sp1.text(10**-3, 0.625, 'Microwave', style='oblique',
        verticalalignment='center', horizontalalignment='center',
        color='black', fontsize=12)

#Ifrared

sp1.add_patch(patches.Rectangle((7*10**-7,0.25 ), 10**-3, .25, color="pink"))
sp1.text(5*10**-5, 0.325, 'Infrared', style='oblique',
        verticalalignment='center', horizontalalignment='center',
        color='black', fontsize=12)


#Visible Spectrum

#sp1.add_patch(patches.Rectangle((3.8*10**-7,0.25 ), 7*10**-7, .25, color="white"))
sp1.add_patch(patches.Rectangle((620*10**-9,0.25 ), 750*10**-9, .25, color="red"))
sp1.add_patch(patches.Rectangle((590*10**-9,0.25 ), 620*10**-9, .25, color="orange"))
sp1.add_patch(patches.Rectangle((570*10**-9,0.25 ), 590*10**-9, .25, color="yellow"))
sp1.add_patch(patches.Rectangle((495*10**-9,0.25 ), 570*10**-9, .25, color="green"))
sp1.add_patch(patches.Rectangle((450*10**-9,0.25 ), 495*10**-9, .25, color="blue"))
sp1.add_patch(patches.Rectangle((380*10**-9,0.25 ), 450*10**-9, .25, color="violet"))
sp1.text(5*10**-7, 0.125, 'Visible', style='oblique',
        verticalalignment='center', horizontalalignment='center',
        color='black', fontsize=12)

#Ultraviolet

sp1.add_patch(patches.Rectangle((10*10**-9,0.25 ), 3.8*10**-7, .25, color="purple"))
sp1.text(5*10**-8, 0.325, 'Ultraviolet', style='oblique',
        verticalalignment='center', horizontalalignment='center',
        color='white', fontsize=8)

#Soft X-rays

sp1.add_patch(patches.Rectangle((10**-10,0.5 ), 4*10**-8, .25, color="yellow"))
sp1.text(2*10**-9, 0.625, 'Soft X-Rays', style='oblique',
        verticalalignment='center', horizontalalignment='center',
        color='black', fontsize=12)

#Hard X-rays

sp1.add_patch(patches.Rectangle((10**-12,0.75 ), 5*10**-10, .25, color="blue"))
sp1.text(2*10**-11, 0.825, 'Hard X-Rays', style='oblique',
        verticalalignment='center', horizontalalignment='center',
        color='black', fontsize=12)

#Gamma

sp1.add_patch(patches.Rectangle((10**-12,0 ), 10**-11, .25, color="Green"))
sp1.text(6*10**-12, 0.325, 'Gamma', style='oblique',
        verticalalignment='center', horizontalalignment='center',
        color='black', fontsize=12)




plt.subplots_adjust(hspace=.5)



plt.show()
