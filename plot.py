import matplotlib.pyplot as plt
import numpy as np
from astropy import constants as const
from astropy import units as u
from astropy.analytic_functions import blackbody_lambda
from scipy.integrate import quad
from numpy import trapz

'''Plot of flux vs. wavelength for a black body of a given temperature'''

temperature = 10000.0 * u.K

# Wien's displacement law
wavemax = (const.b_wien / temperature).to(u.AA) 

waveset = np.logspace(0, np.log10(10 * wavemax.value), num=1000) * u.AA
with np.errstate(all='ignore'):
    flux = blackbody_lambda(waveset, temperature)

# Converting waveset and flux to SI units
SI_convert_m = u.meter
waveset_SI = waveset.to(SI_convert_m)
wavemax_SI = wavemax.to(SI_convert_m)

SI_convert = u.J / (u.second * (u.meter ** 2) * u.meter * u.sr)
flux_SI = flux.to(SI_convert)

print "We have considered a black body at a temperature of %s" %temperature
print "Wavelength at the maximum value of flux is %s or %s" %(wavemax_SI, wavemax)
print "Total radiated flux = %s" %(const.sigma_sb * (temperature ** 4))

# Plotting flux vs wavelength
fig, ax = plt.subplots(figsize=(8, 5))
ax.plot(waveset_SI.value, flux_SI.value)
ax.axvline(wavemax_SI.value, ls='--')
ax.get_yaxis().get_major_formatter().set_powerlimits((0, 1))
ax.set_xlabel(r'$\lambda$ ({0})'.format(waveset_SI.unit))
ax.set_ylabel(r'$B_\lambda(T)$ ({0})'.format(flux_SI.unit))
ax.set_title('Blackbody, T = {0}'.format(temperature))
plt.show()


'''
# U, B, V and R magnitudes of a star (assuming it is a black body), 
# that is 1 solar-radius in diameter at a distance of 10 parsecs away from the Earth
'''

R = const.R_sun / 2.0
r = 10.0 * const.pc


def blackbody_filter(wave_centre, wave_bandwidth, temperature):
	# Funtion for calculating flux of different filters recieved by detectors on the earth

	waveset_filter = np.linspace((wave_centre - wave_bandwidth / 2), (wave_centre + wave_bandwidth / 2), num=100) * u.AA
	waveset_filter_SI = waveset_filter.to(SI_convert_m)

	with np.errstate(all='ignore'):
	    flux_filter = blackbody_lambda(waveset_filter, temperature)

	# Amount of flux recieved by a detector on Earth
	flux_earth = flux_filter * (R / r) ** 2

	flux_filter_SI = flux_filter.to(SI_convert)
	flux_earth_SI = flux_filter_SI * (R / r) ** 2
	
	#area under curve gives the total radiation recieved per steradian
	area = trapz(flux_earth_SI.value, waveset_filter_SI.value) * np.pi * (u.W / u.meter ** 2)

	return area


def magnitude_star(flux_object, flux_reference_object):
	# Function for calculating magnitude of star 

	# According to the Johnson System star Vega is considered as reference for all filters
	magnitude_reference = 0
	magnitude = magnitude_reference - 2.5 * np.log10(flux_object.value / flux_reference_object.value)
	return magnitude

# Filter values taken from 
# http://www.astrophysicsspectator.com/topics/observation/MagnitudesAndColors.html

# For Ultraviolet filter

wavelength_u = 3735.0 * u.AA
wave_bandwidth_u = 485.0 * u.AA
flux_earth_u = 4.22 * 0.485 * 10.0 ** (-9) * (u.W / u.meter ** 2)

flux_u = blackbody_filter(wavelength_u.value, wave_bandwidth_u.value, temperature)
magnitude_u = magnitude_star(flux_u, flux_earth_u)

print "The magnitude for ultraviolet filter is %s" %magnitude_u

# For Blue filter

wavelength_b = 4443.0 * u.AA
wave_bandwidth_b = 831.0 * u.AA
flux_earth_b = 6.2 * 0.831 * 10.0 ** (-9) * (u.W/u.meter ** 2)

flux_b = blackbody_filter(wavelength_b.value, wave_bandwidth_b.value, temperature)
magnitude_b = magnitude_star(flux_b, flux_earth_b)

print "The magnitude for Blue filter is %s" %magnitude_b

# For Visual filter

wavelength_v = 5483.0 * u.AA
wave_bandwidth_v = 827.0 * u.AA
flux_earth_v = 3.55 * 0.827 * 10.0 ** (-9) * (u.W / u.meter ** 2)

flux_v = blackbody_filter(wavelength_v.value, wave_bandwidth_v.value, temperature)
magnitude_v = magnitude_star(flux_v, flux_earth_v)

print "The magnitude for Visual filter is %s" %magnitude_v

# For Red filter

wavelength_r = 6855.0 * u.AA
wave_bandwidth_r = 1742.0 * u.AA
flux_earth_r = 1.795 * 1.742 * 10 ** (-9) * (u.W / u.meter ** 2)

flux_r = blackbody_filter(wavelength_r.value, wave_bandwidth_r.value, temperature)
magnitude_r = magnitude_star(flux_r, flux_earth_r)

print "The magnitude for Red filter is %s" %magnitude_r
