"""Physical constants (cgs unless stated) and unit conversions used by arcwind."""
import numpy as np

C_KMS = 299792.458                 #: speed of light (km/s)
C_CGS = 2.99792458e10              #: speed of light (cm/s)
G_CGS = 6.67430e-8                 #: gravitational constant (cgs)
MSUN_G = 1.98847e33                #: solar mass (g)
AU_CM = 1.495978707e13             #: astronomical unit (cm)
AU_KM = 1.495978707e8              #: astronomical unit (km)
PC_CM = 3.0857e18                  #: parsec (cm)
YR_S = 3.15576e7                   #: Julian year (s)
K_B = 1.380649e-16                 #: Boltzmann constant (erg/K)
H_P = 6.62607015e-27               #: Planck constant (erg s)
M_H = 1.6735575e-24                #: hydrogen-atom mass (g)
ARCSEC_PER_RAD = 206264.806        #: arcsec per radian
GM_SUN_AU_KMS2 = G_CGS * MSUN_G / AU_CM / 1e10   #: G Msun in au (km/s)^2 (~887.1)
FWHM2SIG = 1.0 / (2.0 * np.sqrt(2.0 * np.log(2.0)))   #: Gaussian sigma / FWHM
