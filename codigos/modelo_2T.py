###########################################################
#	       	Modelo con dos capas con T diferente
###########################################################


import math
import numpy as np
import matplotlib.pyplot as plt
import pyfits as pf
import sys
import scipy.optimize as sopt

tau_1_i = 0.5
tau_2_i = 1
T_1_i = 30 #K
T_2_i = 3 #K
v_1_i = -35 #km/s
v_2_i = -30 #km/s
dv_1_i = 4 #km/s
dv_2_i = 1 #km/s
#nu_0 = 89.18853*10**9# Hz
c = 300000 #km/s
v_lsr = -36.4


def B(T,nu):
   #poner T en K y freq en Hz
   h = 6.62606957*10**(-34) #Jxs
   k = 1.3806488*10**(-23)  #JxK-1
   
   B = (h*nu/k) * (np.exp((h*nu)/(k*T)) - 1)**(-1)
   
   return B

def tau_v(tau_gas,v_gas,v_lsr,dv,v):

	tau_v = tau_gas * np.exp(-(v-v_lsr-v_gas)**2/(2*(dv**2)))
	
	return tau_v
	
#v = np.linspace(-50,-15,200)
#v_lsr = 0

#defino la funcion teorica que voy a ajustar al espectro
#def teorico(v,P):
def teorico(nu, P):
   
   #tau_v_1 = tau_v(tau_1,v_1,v_lsr,dv_1,v)
   #tau_v_2 = tau_v(tau_2,v_2,v_lsr,dv_2,v)
   
   #nu = nu_0*(v/c)
   
   tau_v_1 , tau_v_2 , T_1 , T_2 = P

   return B(2.7,nu) * np.exp(-tau_v_1) * np.exp(-tau_v_2) + B(T_1,nu) * (1-np.exp(-tau_v_1)) * np.exp(-tau_v_2) + B(T_2,nu) * (1-np.exp(-tau_v_2))

	
#defino el chi2 que quiero minimizar
	
def chi2(P,nu,spec):

   teorico1 = teorico(nu,P)

   suma = 0
   for i in range(np.shape(nu):

                  resta = ((spec[i]-teorico1[i])**2)/teorico1[i]
                  suma = suma + resta

   return suma

#minimizo con minimos cuadrados
def fit(nu,spec,P):

	resultado = sopt.leastsq(chi2,P,args=(nu,spec), full_output=True)
	return resultado
	

#extraigo espectro del cubo

im_hco = '/Volumes/Elise_ALMA/cleaned_m_uvrange/G305_HCOp.ms.uv.comb_12m1_ACA1_statwt_nopb.3.image.fits'

fitsheader = pf.getheader(im_hco)
nrows,ncols = fitsheader['NAXIS1'],fitsheader['NAXIS2']
nchan = fitsheader['NAXIS3']
arch = pf.open(im_hco)
data = arch[0].data.copy()

spec = data[0,:,nrows/2,ncols/2]

step = fitsheader['CDELT3']
center_vel1 = fitsheader['CRVAL3'] #m/s
center_pix1 = fitsheader['CRPIX3'] #m/s
vi = (center_vel1/1e3-(center_pix1*step1/1e3)+step1/1e3)
nu_0 = fitsheader['RESTFRQ'] #Hz

v = np.zeros(nchan)

for i in range(nchan):
   v[i] = vi*1000 + i*step #en m/s

nu = nu_0*(v/(1000*c))

#ajuste de minimos cuadrados, con los parametros iniciales.
ajuste = fit(nu, spec, [tau_v(tau_1_i , v_1_i, v_lsr , dv_1_i , v), tau_v(tau_2_i , v_2_i, v_lsr , dv_2_i , v), T_1_i , T_2_i]
)            

