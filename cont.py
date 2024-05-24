import numpy as np
import matplotlib.pyplot as plt
import pyneb as pn

def Muv(flux):
    # flux in micro Jansky
    return 23.9-2.5*np.log10(flux)

def lflux(muv):
    # flux in micro Jansky
    return (23.9-muv)/2.5

def Halpha(T):
   #erg/s
   #vol in pc assuming n_e=1
   vol=1.0
   n_e=1.0
   lamb=656.46
   c=3e10
   h_p=6.626176e-27
   hnu=h_p*c/(lamb*1e-7)
   return (3./5.)*hnu*n_e**2*alpha(T)*vol

def alpha(T):
   lam= 2.0*157807.0/T;
   rec= (2.753e-14*lam**1.5/(1.0 + (lam/2.740)**0.407)**2.242)
   return rec

temp0=2e4
z_0=6.14

#Hx/Hbeta
balm_dec=np.array([2.86,1,0.468,0.259,0.159,0.105,0.0731,0.0530])/2.86 
m=np.array([3,4,5,6,7,8,9,10])
lam=3645.0682*(m**2/(m**2-2**2)) 
C = pn.Continuum()
wl = np.arange(950, 10000, 1)
#erg/s/cm3/A
cont = C.get_continuum(tem=temp0, den=1e0, He1_H=0.08, He2_H=0.02, wl=wl, HI_label=None)
#cont=cont*4.*np.pi

lamha=6564.46
v=145.0
wid=v/3e5
print(lam[0],lamha,lam[0]*wid)
print(lam)

#Eros: 2.42e-19 cgs -> LHalpha=
#Flux erg/s/cm^2
mu_tot=23.1*2.0
mu_tot=1.0
# Tc2
F1=2.94e-19/mu_tot
#Halpha1=F1*(4.*np.pi*d_L**2)
#250-300 pc
#vol1=300.0**3*100.0
norm=F1/Halpha(temp0)

balmer1=wl*0.0
for i in range(0,8):
  balmer1=np.where((wl < lam[i]*(1+wid)), balm_dec[i]*Halpha(temp0),balmer1)
  balmer1=np.where((wl > lam[i]*(1-wid)), balmer1, 0)
ew=balmer1/cont

fig, axs = plt.subplots(3)
axs[0].plot(wl, 1.0+ew)
axs[0].set_yscale('linear')
axs[0].set_ylabel(r'EW ($\AA$)')

balmer=wl*0.0
for i in range(0,8):
  balmer=np.where((wl < lam[i]*(1+wid)), balm_dec[i]*Halpha(temp0)/(2.*lam[i]*wid),balmer)
  balmer=np.where((wl > lam[i]*(1-wid)), balmer, 0)
flux=norm*(cont+balmer)
axs[1].plot(wl, flux)
axs[1].set_ylim(1e-22,3e-21)
axs[1].set_yscale('linear')
axs[1].set_ylabel(r'Flux (ergs/s/cm2/$\AA$)')

#erg/s/cm3
#cont=cont*wl
#erg/s/cm3/Hz --> micro Jansky
nu=3e10/(wl*1e-8)
flux=flux*1e29*wl/nu

#axs[2].plot(wl, Muv(flux))
#axs[2].set_yscale('linear')
#axs[2].set_ylabel(r'AB Magnitude')
#axs[2].set_ylim(35,29)
axs[2].plot(wl, (flux))
axs[2].set_yscale('log')
axs[2].set_ylabel(r'Flux micro Jy')
axs[2].set_ylim(5e-5,5e-2)
axs[2].set_xlabel(r'$\lambda (\AA)$')

print('H-alpha EW =', max(ew))
print(' Flux in micro Jansky cont=',10**lflux(31.))
print(' Flux in micro Jansky Ha line=',max(flux))
r=max(flux)/10**lflux(31.)
print('ratio =', r, 'd lambda_obsa', 1000.*1e-4/r, 'micro m', 'd v_rest', 3e5*1000*1e-4/r/2./4.69)

#EW for emission lines: EW=2*delta Lambda * (F_line/F_cont-1)~

plt.show()
