import healpy
import os
import numpy
import scipy
import astropy.io.fits
import time
import sys
import multiprocessing
import Dust_util_2comp
import matplotlib.pyplot
import pandas
import pysm3
import pysm3.units
import shutil

root_dir = os.path.abspath('.')

# transmission spectra of WMAP
# wavenumber in GHz
steps = 60
# 23 GHz
freq_023  = numpy.arange(23-0.01*steps, 23+0.01*steps, 0.01)
trans_023 = numpy.zeros(2*steps, dtype=float) + 1
# 33 GHz
freq_033  = numpy.arange(33-0.01*steps, 33+0.01*steps, 0.01)
trans_033 = numpy.zeros(2*steps, dtype=float) + 1
# 41 GHz
freq_041  = numpy.arange(41-0.01*steps, 41+0.01*steps, 0.01)
trans_041 = numpy.zeros(2*steps, dtype=float) + 1
# 61 GHz
freq_061  = numpy.arange(61-0.01*steps, 61+0.01*steps, 0.01)
trans_061 = numpy.zeros(2*steps, dtype=float) + 1

# transmission spectra of Planck LFI
# wavenumber in GHz
# transmission is normalized to have an integral of 1 for LFI
# https://wiki.cosmos.esa.int/planck-legacy-archive/index.php/The_RIMO
LFI_RIMO  = astropy.io.fits.open(root_dir+"/PR3-2018/LFI_RIMO_R3.31.fits")
# 30 GHz
freq_030  = LFI_RIMO["bandpass_030"].data.field("wavenumber")
trans_030 = LFI_RIMO["bandpass_030"].data.field("transmission")
# 44 GHz
freq_044  = LFI_RIMO["bandpass_044"].data.field("wavenumber")
trans_044 = LFI_RIMO["bandpass_044"].data.field("transmission")
# 70 GHz
freq_070  = LFI_RIMO["bandpass_070"].data.field("wavenumber")
trans_070 = LFI_RIMO["bandpass_070"].data.field("transmission")
# transmission spectra of Planck HFI
# wavenumber in cm^-1
# transmission is normalized to 1 at the maximum for HFI
# https://wiki.cosmos.esa.int/planck-legacy-archive/index.php/The_RIMO
HFI_RIMO  = astropy.io.fits.open(root_dir+"/PR3-2018/HFI_RIMO_R3.00.fits")
# converse wavenumber in cm^-1 to frequency in GHz by multiplying by 29979245800/1e9
# 100 GHz
freq_100  = 29979245800/1e9*HFI_RIMO["BANDPASS_F100"].data.field("wavenumber")
trans_100 = HFI_RIMO["BANDPASS_F100"].data.field("transmission")
# 143 GHz
freq_143  = 29979245800/1e9*HFI_RIMO["BANDPASS_F143"].data.field("wavenumber")
trans_143 = HFI_RIMO["BANDPASS_F143"].data.field("transmission")
# 217 GHz
freq_217  = 29979245800/1e9*HFI_RIMO["BANDPASS_F217"].data.field("wavenumber")
trans_217 = HFI_RIMO["BANDPASS_F217"].data.field("transmission")
# 353 GHz
freq_353  = 29979245800/1e9*HFI_RIMO["BANDPASS_F353"].data.field("wavenumber")
trans_353 = HFI_RIMO["BANDPASS_F353"].data.field("transmission")
# 545 GHz
freq_545  = 29979245800/1e9*HFI_RIMO["BANDPASS_F545"].data.field("wavenumber")
trans_545 = HFI_RIMO["BANDPASS_F545"].data.field("transmission")
# 857 GHz
freq_857  = 29979245800/1e9*HFI_RIMO["BANDPASS_F857"].data.field("wavenumber")
trans_857 = HFI_RIMO["BANDPASS_F857"].data.field("transmission")


#--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*

# When convolving signal with transmission spectra, integration should not start from the minimum value, as it would introduce significant computational errors.
# Configuration related to integration over frequency for HFI channels
# The selection of lower limit and upper limit of integeration ensures that the transmission larger than 1e-6.
# There are some negative values in trans_143,
# (num = [229 230 231 232 233 234 235 236 252 253 254 255 256 257 258 275 276 277 278 279])
# we have excluded them.
# function freq_low_limit() gives the index of the lower integration limit for each channel
def freq_low_limit(freq_str):
    if freq_str == "100":
        return 159
    elif freq_str == "143":
        return 305
    elif freq_str == "217":
        return 351
    elif freq_str == "353":
        return 586
    elif freq_str == "545":
        return 707
    elif freq_str == "857":
        return 1085
    else:
        return 1+5

# function freq_high_limit() gives the index of the upper integration limit for each channel
def freq_high_limit(freq_str):
    if freq_str == "100":
        return 389
    elif freq_str == "143":
        return 492
    elif freq_str == "217":
        return 816
    elif freq_str == "353":
        return 1238
    elif freq_str == "545":
        return 2408
    elif freq_str == "857":
        return 5529
    else:
        return eval("freq_"+freq_str).shape[0]-1-5

#--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*
# Functions for unit conversion: x_trans(), b_prime_RJ(), b_prime_CMB()

# Planck constant in SI unit
h_Planck = 6.62607015e-34
# Speed of light in SI unit
c_light = 299792458
# Boltzmann constant in SI unit
k_Boltzman = 1.380649e-23
# Average temperature of CMB in SI unit
T_0 = 2.7255
# Reference: A&A 571, A9 (2014) P12

# x_trans = h \nu/kT
def x_trans(freq):
# freq in GHz
    freq = float(freq)
    return h_Planck*1e9*freq/k_Boltzman/T_0

# From Rayleigh-Jeans unit to SI unit
def b_prime_RJ(freq):
# freq in GHz
    freq = float(freq)
    return 2*k_Boltzman*(1e9*freq)**2/c_light**2

# From CMB unit to SI unit
def b_prime_CMB(freq):
# freq in GHz
    freq = float(freq)
    x = x_trans(freq)
    return 2*k_Boltzman*(1e9*freq)**2/c_light**2*x**2*numpy.exp(x)/(numpy.exp(x)-1)**2

#--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*
# FWHM from Planck 2018, I. A&A 641, A1(2020), Table 4
# beam_fwhm() gives FWHM for each channel, in radians
def beam_fwhm(freq_str):
    beam_dict = {
        "023": 0.88*60,
        "030": 32.29,
        "033": 0.66*60,
        "041": 0.51*60,
        "044": 27.94,
        "061": 0.35*60,
        "070": 13.08,
        "100": 9.66,
        "143": 7.22,
        "217": 4.90,
        "353": 4.92,
        "545": 4.67,
        "857": 4.22
    }
    beam = beam_dict.get(freq_str, 0)
    return beam / 60 * numpy.pi / 180

def rebeam(sky_map, Nside):
    N_side = healpy.get_nside(sky_map)
    map_0  = healpy.ud_grade(sky_map, nside_out = 4096)
    # 12 N_side^2 = (lmax + 1)^2
    alm    = healpy.map2alm(map_0, lmax = int(2*numpy.sqrt(3)*N_side))
    map_1  = healpy.alm2map(alm, nside=4096)
    sky_map= healpy.ud_grade(map_1, nside_out = Nside)
    return sky_map


#--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*
# Before performing statistics, mask the regions contaminated by compact sources.
def compact_sourece_mask(freq_str, Nside, fwhm1, fwhm2, fwhm3):
    # Compact source mask for single frequency channel
    # FWHM of compact sources from COM_PCCS_030/044/070_R2.04.fits (PCCS2), COM_PCCS_100/143/217/353/545/857_R2.01.fits (PCCS2), COM_PCCS_100/143/217/353/545/857-excluded_R2.01.fits (PCCS2E)
    # PCCS2 covers most of the sky and allows the user to produce subsamples at higher reliabilities than the target 80% integral reliability of the catalogue. 
    # PCCS2E contains sources detected in sky regions where the diffuse emission makes it difficult to quantify the reliability of the detections. 
    # A&A 594, A26 (2016)
    # https://wiki.cosmos.esa.int/planckpla2015/index.php/Catalogues

    if freq_str=="030" or freq_str=="044" or freq_str=="070":
        mask = numpy.zeros(12*Nside**2, dtype=float) + 1
        compact_catalogue = astropy.io.fits.open(root_dir+"/PR2-2015/COM_PCCS_"+freq_str+"_R2.04.fits")
        galactic_longitude= compact_catalogue[1].data.field("GLON") # in unit of degree
        galactic_latitude = compact_catalogue[1].data.field("GLAT") # in unit of degree
        flux              = compact_catalogue[1].data.field("DETFLUX") # "DETFLUX": Flux density of source as determined by detection method, in unit of MJy
        FWHM              = compact_catalogue[1].data.field("GAU_FWHM_EFF") # "GAU_FWHM_EFF": Gaussian fit effective FWHM, in unit arcmin

    elif freq_str=="100" or freq_str=="143" or freq_str=="217" or freq_str=="353" or freq_str=="545" or freq_str=="857":
        mask = numpy.zeros(12*Nside**2, dtype=float) + 1
        compact_catalogue_1 = astropy.io.fits.open(root_dir+"/PR2-2015/COM_PCCS_"+freq_str+"_R2.01.fits")
        # PCCS2, the sources have been detected in regions of the sky where it is possible to estimate the reliability of the detections, 
        # either statistically or by using external catalogues
        galactic_longitude_1= compact_catalogue_1[1].data.field("GLON") # in unit of degree
        galactic_latitude_1 = compact_catalogue_1[1].data.field("GLAT") # in unit of degree
        flux_1              = compact_catalogue_1[1].data.field("DETFLUX") # "DETFLUX": Flux density of source as determined by detection method, in unit of MJy
        FWHM_1              = compact_catalogue_1[1].data.field("GAU_FWHM_EFF") # "GAU_FWHM_EFF": Gaussian fit effective FWHM, in unit arcmin
        compact_catalogue_2 = astropy.io.fits.open(root_dir+"/PR2-2015/COM_PCCS_"+freq_str+"-excluded_R2.01.fits")
        # PCCS2E, the detected sources are located in regions of the sky where it is not possible to make an estimate of their reliability. 
        galactic_longitude_2= compact_catalogue_2[1].data.field("GLON") # in unit of degree
        galactic_latitude_2 = compact_catalogue_2[1].data.field("GLAT") # in unit of degree
        flux_2              = compact_catalogue_2[1].data.field("DETFLUX") # Flux density of source as determined by detection method
        FWHM_2              = compact_catalogue_2[1].data.field("GAU_FWHM_EFF") # GAU_FWHM_EFF, Gaussian fit effective FWHM, in unit arcmin
        # For 100/143/217/353/545/857 GHz, combine PCCS2 and PCCS2E
        galactic_longitude  = numpy.concatenate((galactic_longitude_1,galactic_longitude_2), axis=0)
        galactic_latitude   = numpy.concatenate((galactic_latitude_1, galactic_latitude_2 ), axis=0)
        flux                = numpy.concatenate((flux_1, flux_2), axis=0)
        FWHM                = numpy.concatenate((FWHM_1, FWHM_2), axis=0)

    compact_source_data = numpy.array([galactic_longitude, galactic_latitude, flux, FWHM])
    #   Transpose: arrange by a certain row
    # No Transpose: arrange by a certain column
    # numpy.argsort(-compact_source_data[2,:]) returns the indices in reversely sorted order, according to the value of DETFLUX (2nd row of compact_source_data)
    compact_source_data = compact_source_data.T[numpy.argsort(-compact_source_data[2,:])].T
    # "number" is the number of comapact sources in each frequency channel
    number = compact_source_data.shape[1]
    # In each channel,
    # mask radius = 3 FWHM for sources with the top 10% intensity
    for i in range(0, int(0.1*number), 1):
        # Direction of the i-th compact source
        # lonlat = True, input angles are assumed to be longitude and latitude in degree
        vector = healpy.ang2vec(compact_source_data[0][i], compact_source_data[1][i], lonlat=True)
        # Mask radius from arcmin to radians
        ipix_disc = healpy.query_disc(nside=Nside, vec=vector, radius=fwhm1*compact_source_data[3][i]/60*numpy.pi/180)
        mask[ipix_disc] = 0.0
    # mask radius = 2 FWHM for sources with the top 10%-30% intensity
    for i in range(int(0.1*number), int(0.3*number), 1):
        vector = healpy.ang2vec(compact_source_data[0][i], compact_source_data[1][i], lonlat=True)
        ipix_disc = healpy.query_disc(nside=Nside, vec=vector, radius=fwhm2*compact_source_data[3][i]/60*numpy.pi/180)
        mask[ipix_disc] = 0.0
    # mask radius = 1.5 FWHM for others
    for i in range(int(0.3*number), number, 1):
        vector = healpy.ang2vec(compact_source_data[0][i], compact_source_data[1][i], lonlat=True)
        ipix_disc = healpy.query_disc(nside=Nside, vec=vector, radius=fwhm3*compact_source_data[3][i]/60*numpy.pi/180)
        mask[ipix_disc] = 0.0
    return mask

def mask_CO(Nside):
# Reading of HFI detector at each channel, in unit of muK_CMB or MJy_sr, and the output with Nside, FWHM = 9.66'
    Nside = int(Nside)
    # High-resolution CO(2-1) line from COMMANDER, in K_RJ km/s, with FWHM = 7.5'
    CO21 = healpy.read_map(root_dir+"/PR2-2015/COM_CompMap_CO21-commander_2048_R2.00.fits", h=False, field=0)
    CO21 = healpy.ud_grade(CO21, Nside)
    CO_mask = numpy.zeros(12*Nside**2, dtype=float)
    # Zero out high-signal regions (>1 K_RJ km/s) of CO(2-1) map, ApJ 798, 88 (2015) P3
    pixel_CO21_low = numpy.where(CO21 < 0.4)[0]
    CO_mask[pixel_CO21_low] = 1
    CO_mask = healpy.ud_grade(CO_mask, Nside)
    CO_mask[CO_mask < 0.9] = 0.0
    return CO_mask

def Statistic_mask(Nside, fwhm1, fwhm2, fwhm3):
    # Mask out compact sources
    compact_source = compact_sourece_mask("030", Nside, fwhm1, fwhm2, fwhm3)
    for freq_str in ["044", "070", "100", "143", "217", "353", "545", "857"]:
        compact_source = compact_source * compact_sourece_mask(freq_str, Nside, fwhm1, fwhm2, fwhm3)
    # Mask out galactic plane, with 80% sky coverage
    # field = 0 -- 20% sky coverage
    # field = 1 -- 40% sky coverage
    # field = 2 -- 60% sky coverage
    # field = 3 -- 70% sky coverage
    # field = 4 -- 80% sky coverage
    # field = 5 -- 90% sky coverage
    galaxy = healpy.read_map(root_dir+"/PR2-2015/HFI_Mask_GalPlane-apo0_2048_R2.00.fits", h=False, field=4)
    galaxy = galaxy.astype(numpy.float64)
    galaxy = healpy.ud_grade(galaxy, Nside)
    galaxy[galaxy < 0.9] = 0.0
    # Mask out CO
    CO = mask_CO(Nside)
    total_mask = compact_source * galaxy * CO
    return total_mask

#--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*
# Map of CMB anisotropies
def CMB_map(freq_str, unit, Nside):
    Nside = int(Nside)
    # Inpainted I intensity map of CMB (SMICA)
    CMB_map = healpy.read_map(root_dir+"/PR3-2018/COM_CMB_IQU-smica_2048_R3.00_full.fits", h=False, field=5) # in unit of K_cmb

    if unit == "muK_CMB":
        CMB_map = 1e6*CMB_map
    elif unit == "MJy_sr":
        part1 = 0
        part2 = 0
        # Transmission spectrum of detectors
        freq_array = eval("freq_"+freq_str)
        trans_array= eval("trans_"+freq_str)
        for i in range(freq_low_limit(freq_str), freq_high_limit(freq_str), 1):
            freq = freq_array[i]
            delta_freq = 0.5*(freq_array[i+1] - freq_array[i-1])
            trans = trans_array[i]
            # part1 is numerator
            part1 = part1 + b_prime_CMB(freq)*trans*delta_freq
            # part2 is denominator
            part2 = part2 + (freq/float(freq_str))**(-1)*trans*delta_freq
        # For CMB anisotropies, T_CMB is independent with frequency
        # From W m^-2 Hz^-1 sr^-1 to MJy/sr
        CMB_map = 1e-6*1e26*CMB_map*part1/part2
    beam1 = healpy.gauss_beam(fwhm= 5/60*numpy.pi/180, pol=False)
    beam2 = healpy.gauss_beam(fwhm=60/60*numpy.pi/180, pol=False)
    alm_CMB = healpy.map2alm(CMB_map)
    alm_CMB = healpy.smoothalm(alm_CMB, beam_window=beam2/beam1)
    CMB_map = healpy.alm2map(alm_CMB, nside=Nside)
    return CMB_map

#--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*
# Map of free-free

# FWHM = 60', Planck 2015, X. A&A 594, A10 (2016), Table 5
# Emission Measure of free-free in pc cm^-6
free_free_EM = healpy.read_map(root_dir+"/PR2-2015/COM_CompMap_freefree-commander_0256_R2.00.fits", h=False, field=0, dtype=numpy.float64)
# Temperature of electrons, in unit of K.
free_free_Te = healpy.read_map(root_dir+"/PR2-2015/COM_CompMap_freefree-commander_0256_R2.00.fits", h=False, field=3, dtype=numpy.float64)
T4 = free_free_Te/1e4

def gaunt_ff(freq):
# Gaunt factor of free-free, Planck 2015, X. A&A 594, A10 (2016), Table 4
# freq in GHz, is \nu9 in Planck 2015, X. A&A 594, A10 (2016), Table 4
    freq = float(freq)
    gaunt = numpy.log(numpy.exp(5.960 - numpy.sqrt(3)/numpy.pi*numpy.log(freq/T4**1.5)) + numpy.exp(1))
    return gaunt

def tau_ff(freq):
# Optical depth of free-free, Planck 2015, X. A&A 594, A10 (2016), Table 4
# freq in GHz, is \nu9 in Planck 2015, X. A&A 594, A10 (2016), Table 4
    freq = float(freq)
    tau = 0.05468/free_free_Te**1.5/freq**2*free_free_EM*gaunt_ff(freq)
    return tau

def K_RJ_ff(freq):
# Brightness temperature of free-free in unit of K_RJ, Planck 2015, X. A&A 594, A10 (2016), Table 4
# freq in GHz
    freq = float(freq)
    s_ff = free_free_Te*(1-numpy.exp(-tau_ff(freq)))
    return s_ff

def intensity_ff(freq):
# Intensity of free-free at freq in W m^-2 Hz^-1 sr^-1
    freq = float(freq)
    return b_prime_RJ(freq)*K_RJ_ff(freq)

def reading_ff(freq_str, unit, Nside):
# Reading of HFI detector at each channel, in unit of muK_CMB or MJy_sr, and the output with Nside
    Nside = int(Nside)
    # Nside = 256 for COM_CompMap_freefree-commander_0256_R2.00.fits
    part1 = numpy.zeros(12*256**2, dtype=numpy.float64)
    part2 = 0
    # Transmission spectrum of detectors
    freq_array = eval("freq_"+freq_str)
    trans_array= eval("trans_"+freq_str)
    for i in range(freq_low_limit(freq_str), freq_high_limit(freq_str), 1):
        freq = freq_array[i]
        delta_freq = 0.5*(freq_array[i+1] - freq_array[i-1])
        trans = trans_array[i]
        # part1 is numerator
        part1 = part1 + intensity_ff(freq)*trans*delta_freq
        if unit == "muK_CMB":
            # part2 is denominator
            part2 = part2 + b_prime_CMB(freq)*trans*delta_freq
        elif unit == "MJy_sr":
            part2 = part2 + (freq/float(freq_str))**(-1)*trans*delta_freq
    if unit == "muK_CMB":
        # From K_CMB to muK_CMB
        map_ff = 1e6*part1/part2
    elif unit == "MJy_sr":
        # From W m^-2 Hz^-1 sr^-1 to MJy/sr
        map_ff = 1e-6*1e26*part1/part2
    # Up_grade from nside=256 to nside=Nside
#    map_ff = healpy.ud_grade(map_ff, Nside)
    map_ff = rebeam(map_ff, Nside)
    return map_ff

#--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*
# Map of synchrotron

_cached_synchro_I = None

#def synch_template(Nside):
#    global _cached_synchro_I
#    if _cached_synchro_I is None:
#        # Extension = "SYNC-TEMP" is the template for synchrotron radiation (Intensity VS Frequency)
#        synchro_fits = astropy.io.fits.open(root_dir+"/PR2-2015/COM_CompMap_Synchrotron-commander_0256_R2.00.fits")
#        # Reference Frequency of synchrotron Stokes I is 408.0 MHz
#        synchro_freq_ref = 0.408
#        
#        # map of synchrotron radiation in unit of muK_RJ
#        synchro_I = healpy.read_map(root_dir+"/PR2-2015/COM_CompMap_Synchrotron-commander_0256_R2.00.fits", h=False, field=0) # in muK_RJ, at 408 MHz
#        # From muK_RJ to SI unit
#        synchro_I = 1e-6*b_prime_RJ(synchro_freq_ref)*synchro_I
##        alm_synchro = healpy.map2alm(synchro_I)
##        synchro_I   = healpy.alm2map(alm_synchro, nside=Nside)
#        synchro_I   = healpy.ud_grade(synchro_I, Nside)
#                
#        galaxy_mask = healpy.read_map(root_dir+"/PR2-2015/HFI_Mask_GalPlane-apo0_2048_R2.00.fits", h=False, field=3)
#        galaxy_mask = healpy.ud_grade(galaxy_mask, Nside)
#        galaxy_mask[galaxy_mask < 0.9] = 0.0
#        synchro_mask = compact_sourece_mask("030", Nside, 2.5, 2, 1.5) * galaxy_mask
#        synchro_I_masked = synchro_I.copy()
#        synchro_I_masked[synchro_mask == 0] = healpy.UNSEEN
#        monopole, dipole = healpy.fit_dipole(synchro_I_masked)
#        npix = healpy.nside2npix(Nside)
#        directions = numpy.array(healpy.pix2vec(Nside, numpy.arange(npix))).T
#        dipole_map = numpy.dot(directions, dipole)
#        
#        synchro_I = synchro_I - monopole - dipole_map
#        _cached_synchro_I = synchro_I
#    return _cached_synchro_I

def synch_template(Nside):
    global _cached_synchro_I
    if _cached_synchro_I is None:
        # Extension = "SYNC-TEMP" is the template for synchrotron radiation (Intensity VS Frequency)
        synchro_fits = astropy.io.fits.open(root_dir+"/PR2-2015/COM_CompMap_Synchrotron-commander_0256_R2.00.fits")
        # Reference Frequency of synchrotron Stokes I is 408.0 MHz
        synchro_freq_ref = 0.408
    
        # map of synchrotron radiation in unit of muK_RJ
        synchro_I = healpy.read_map(root_dir+"/PR2-2015/COM_CompMap_Synchrotron-commander_0256_R2.00.fits", h=False, field=0) # in muK_RJ, at 408 MHz

        # https://doi.org/10.1051/0004-6361/201525659, Table 1
        monopole = 8.9e6
        dipole   = numpy.array([3.2e6, 0.7e6, -0.8e6])
        npix = healpy.nside2npix(256)
        directions = numpy.array(healpy.pix2vec(256, numpy.arange(npix))).T
        dipole_map = numpy.dot(directions, dipole)        
        synchro_I = synchro_I - monopole - dipole_map

        # From muK_RJ to SI unit
        synchro_I = 1e-6*b_prime_RJ(synchro_freq_ref)*synchro_I
#        synchro_I = healpy.ud_grade(synchro_I, Nside)
        synchro_I = rebeam(synchro_I, Nside)
        _cached_synchro_I = synchro_I    
    return _cached_synchro_I

def clear_cache():
    global _cached_synchro_I
    _cached_synchro_I = None

def reading_synch(model, freq_str, unit, Nside):
    if model == "Planck":
        # Extension = "SYNC-TEMP" is the template for synchrotron radiation (Intensity VS Frequency)
        synchro_fits = astropy.io.fits.open(root_dir+"/PR2-2015/COM_CompMap_Synchrotron-commander_0256_R2.00.fits")
        # Reference Frequency of synchrotron Stokes I is 408.0 MHz
        synchro_freq_ref = 0.408
        synchro_freq      = synchro_fits["SYNC-TEMP"].data.field("NU") # GHz
        synchro_intensity = synchro_fits["SYNC-TEMP"].data.field("I") # W Hz^-1 m^-2 sr^-1
        # synchro_template is the spline interpolation obtained based on synchro_freq and synchro_intensity, 
        # representing the relationship between synchrotron intensity and frequency (in log-log).
        synchro_template = scipy.interpolate.splrep(numpy.log10(synchro_freq), numpy.log10(synchro_intensity), k=2)
        # synchro_intensity_ref is the synchrotron intensity at 408 MHz in SI unit
        synchro_intensity_ref = 10**scipy.interpolate.splev(numpy.log10(synchro_freq_ref), synchro_template)

        Nside = int(Nside)
        map_synchro = synch_template(Nside)
        part1 = 0
        part2 = 0
        # Transmission spectrum of HFI detectors
        freq_array = eval("freq_"+freq_str)
        trans_array= eval("trans_"+freq_str)
        for i in range(freq_low_limit(freq_str), freq_high_limit(freq_str), 1):
            freq = freq_array[i]
            delta_freq = 0.5*(freq_array[i+1] - freq_array[i-1])
            trans = trans_array[i]
            # synchro_intensity is synchrotron intensity at freq in SI unit
            synchro_intensity = 10**scipy.interpolate.splev(numpy.log10(freq), synchro_template)
            # Suppose that the spectrum of synchrotron is identical alone all directions
            # part1 is numerator
            part1 = part1 + synchro_intensity/synchro_intensity_ref*trans*delta_freq
            if unit == "muK_CMB":
                # part2 is denominator
                part2 = part2 + b_prime_CMB(freq)*trans*delta_freq
            elif unit == "MJy_sr":
                part2 = part2 + (freq/float(freq_str))**(-1)*trans*delta_freq
        if unit == "muK_CMB":
            # From K_CMB to muK_CMB
            map_synchro = 1e6*map_synchro*part1/part2
        elif unit == "MJy_sr":
            # From W m^-2 Hz^-1 sr^-1 to MJy/sr
            map_synchro = 1e-6*1e26*map_synchro*part1/part2
        clear_cache()
    else:
        step = 1
        if model == "s1":
            sky = pysm3.Sky(nside=512, preset_strings=["s1"])
        elif model == "s2":
            sky = pysm3.Sky(nside=512, preset_strings=["s2"])
        elif model == "s3":
            sky = pysm3.Sky(nside=512, preset_strings=["s3"])
        elif model == "s5":
            sky = pysm3.Sky(nside=512, preset_strings=["s5"])
        elif model == "s7":
            sky = pysm3.Sky(nside=512, preset_strings=["s7"])
        # part1 is numerator
        part1 = numpy.zeros(12*512**2, dtype=numpy.float64)
        # part2 is denominator
        part2 = 0
        # Transmission spectrum of Planck detectors
        freq_array = eval("freq_"+freq_str)
        trans_array= eval("trans_"+freq_str)
        for i in range(freq_low_limit(freq_str), freq_high_limit(freq_str), step):
            freq = freq_array[i]
            delta_freq = 0.5*(freq_array[i+step] - freq_array[i-step])
            trans = trans_array[i]
            synch_intensity = sky.get_emission(freq * pysm3.units.GHz) # in unit of muK_RJ
            synch_intensity = pysm3.apply_smoothing_and_coord_transform(synch_intensity, fwhm=60/60*pysm3.units.deg) # FWHM = 60 arcmin
            synch_intensity = synch_intensity.value
            synch_intensity = synch_intensity[0]
            synch_intensity = b_prime_RJ(freq) * 1e-6 * synch_intensity # in SI unit
            part1 = part1 + synch_intensity*trans*delta_freq
            if unit == "muK_CMB":
                part2 = part2 + b_prime_CMB(freq)*trans*delta_freq
            elif unit == "MJy_sr":
                part2 = part2 + (freq/float(freq_str))**(-1)*trans*delta_freq
        if unit == "muK_CMB":
            # From K_CMB to muK_CMB
            map_synchro = 1e6*part1/part2
        elif unit == "MJy_sr":
            # From W m^-2 Hz^-1 sr^-1 to MJy/sr
            map_synchro = 1e-6*1e26*part1/part2
#        map_synchro = healpy.ud_grade(map_synchro, Nside)
        map_synchro = rebeam(map_synchro, Nside)
    return map_synchro

##--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*
# Map of 94/100 GHz line

def reading_xline(freq_str, unit, Nside):
# Reading of HFI detector at each channel, in unit of muK_CMB, and the output with Nside
# Non-zero at 100 GHz
    Nside = int(Nside)
    if unit == "MJy_sr":
        map_xline = numpy.zeros(12*Nside**2, dtype=float)
    elif unit == "muK_CMB":
        if freq_str == "100":
            map_xline = healpy.read_map(root_dir+"/PR2-2015/COM_CompMap_xline-commander_0256_R2.00.fits", h=False, field=0) # muK_CMB, FWHM = 60.0'
            # field = 0: Intensity map
            # field = 1: Mean Intensity
            map_xline = rebeam(map_xline, Nside)
        else:
            map_xline = numpy.zeros(12*Nside**2, dtype=float)
    return map_xline

#--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*
# thermal dust
def reading_dust(model, freq_str, unit, Nside):
# in unit of muK_CMB or MJy_sr, with FWHM = 9.66', and the output with Nside
    Nside = int(Nside)
    step = 5
    if model == "Planck13" or model == "Planck15-G" or model == "Irfan19":
        if model == "Planck13": # FWHM = 5', A&A 571, A11 (2014)
            # Optical depth of dust at 353 GHz from Planck 2013 release
            dust_tau = healpy.read_map(root_dir+"/PR1-2013/HFI_CompMap_ThermalDustModel_2048_R1.20.fits", h=False, field=0)
            # Dust emission spectral index from Planck 2013
            dust_index = healpy.read_map(root_dir+"/PR1-2013/HFI_CompMap_ThermalDustModel_2048_R1.20.fits", h=False, field=6)
            # Dust equilibrium temperature of dust from Planck 2013 release, in unit of K
            dust_temperature = healpy.read_map(root_dir+"/PR1-2013/HFI_CompMap_ThermalDustModel_2048_R1.20.fits", h=False, field=4)
        elif model == "Planck15-G":
            # Optical depth at 353 GHz from Planck 2015, FWHM = 5'
            dust_tau = healpy.read_map(root_dir+"/PR2-2015/COM_CompMap_Dust-GNILC-Model-Opacity_2048_R2.01.fits", h=False, field=0)
            # Dust emission spectral index from Planck 2015, FWHM = 5'
            dust_index = healpy.read_map(root_dir+"/PR2-2015/COM_CompMap_Dust-GNILC-Model-Spectral-Index_2048_R2.01.fits", h=False, field=0)
            # Dust equilibrium temperature from Planck 2015, in unit of K, FWHM = 5'
            dust_temperature = healpy.read_map(root_dir+"/PR2-2015/COM_CompMap_Dust-GNILC-Model-Temperature_2048_R2.01.fits", h=False, field=0)
        elif model == "Irfan19":
            # Optical depth at 353 GHz from A&A 623, A21 (2019), FWHM = 5', https://cdsarc.cds.unistra.fr/ftp/J/A+A/623/A21/, nside = 2048
            tau = healpy.read_map(root_dir+"/Irfan19/tau.fits", h=False, field=0)
            # Dust emission spectral index from A&A 623, A21 (2019), FWHM = 5'
            beta = healpy.read_map(root_dir+"/Irfan19/beta.fits", h=False, field=0)
            # Dust equilibrium temperature from A&A 623, A21 (2019), in unit of K, FWHM = 5'
            temp = healpy.read_map(root_dir+"/Irfan19/temp.fits", h=False, field=0)
            # There are NaN values in the model from Irfan19, we need to inpaint them.
            dust_model_Irfan19 = numpy.array([beta, tau, temp])
            N_side = 2048
            # If phi is not given or None, theta is interpreted as pixel number, otherwise, theta, phi are angular coordinates in radians.
            neighbour_array = healpy.get_all_neighbours(nside=N_side, theta=numpy.arange(0, 12*N_side**2, 1), phi=None).T
            N_iter = 2000
            # N_iter is the count of iteration for inpainting.
            pixel_with_nan = numpy.where(numpy.isnan(temp))[0]
            # Get the neighboring pixels of the pixels that need to be inpainted
            neighbours = neighbour_array[pixel_with_nan]
            # valid_neighbours are the valid neighboring pixels
            # Since in the healpix algorithm, a pixel has 8 or 7 neighbours, when a pixel has 7 neighbours, healpy.get_all_neighbours() returns a -1 to mark the invalid neighbour
            valid_neighbours = neighbours >= 0
            # Count the neighboring pixels
            counts = numpy.sum(valid_neighbours, axis=1)
            dust_model_Irfan19[:, pixel_with_nan] = 0
            for ii in range(0, N_iter):
                # Sum the values of the valid neighboring pixels
                sums = numpy.sum(dust_model_Irfan19[:, neighbours] * valid_neighbours, axis=2)
                # Calculate the mean value of the neighboring pixels for each masked pixel
                mean_values = sums / counts
                # Update the values with NaN
                dust_model_Irfan19[:, pixel_with_nan] = mean_values
            beta, tau, temp = dust_model_Irfan19
            dust_tau = tau
            dust_index = beta
            dust_temperature = temp
            healpy.write_map(root_dir+"/results_dust/Irfan19_tau.fits", tau, nest=False,overwrite="True",dtype=numpy.float64)
            healpy.write_map(root_dir+"/results_dust/Irfan19_beta.fits",beta,nest=False,overwrite="True",dtype=numpy.float64)
            healpy.write_map(root_dir+"/results_dust/Irfan19_temp.fits",temp,nest=False,overwrite="True",dtype=numpy.float64)
        # Reference frequency for single component model
        freq_d = 353
        # A&A 596, A109 (2016) Equation (8)
        # part1 is numerator
        part1 = numpy.zeros(12*2048**2, dtype=numpy.float64)
        # part2 is denominator
        part2 = 0
        # Transmission spectrum of HFI detectors
        freq_array = eval("freq_"+freq_str)
        trans_array= eval("trans_"+freq_str)
        for i in range(freq_low_limit(freq_str), freq_high_limit(freq_str), step):
            freq = freq_array[i]
            delta_freq = 0.5*(freq_array[i+step] - freq_array[i-step])
            trans = trans_array[i]
            x_trans_freq = h_Planck*(1e9*freq)/k_Boltzman/dust_temperature
            Black_body_spectrum = 2*h_Planck*(1e9*freq)**3/c_light**2/(numpy.exp(x_trans_freq)-1)
            part1 = part1 + (freq/freq_d)**dust_index * Black_body_spectrum * trans*delta_freq
            if unit == "muK_CMB":
                part2 = part2 + b_prime_CMB(freq)*trans*delta_freq
            elif unit == "MJy_sr":
                part2 = part2 + (freq/float(freq_str))**(-1)*trans*delta_freq
        if unit == "muK_CMB":
            # From K_CMB to muK_CMB
            dust_map = 1e6*dust_tau*part1/part2
        elif unit == "MJy_sr":
            # From W m^-2 Hz^-1 sr^-1 to MJy/sr
            dust_map = 1e-6*1e26*dust_tau*part1/part2

    elif model == "Planck15-C":
        # Reference map at 545 GHz from Planck 2015 Commander, FWHM = 60', in unit of muK_RJ
        map_ref = healpy.read_map(root_dir+"/PR2-2015/COM_CompMap_dust-commander_0256_R2.00.fits", h=False, field=0)
        # in unit of W m^-2 Hz^-1 sr^-1
        map_ref = b_prime_RJ(545) * 1e-6 * map_ref
        # Dust emission spectral index from Planck 2015 Commander, FWHM = 60'
        dust_index = healpy.read_map(root_dir+"/PR2-2015/COM_CompMap_dust-commander_0256_R2.00.fits", h=False, field=6)
        # Dust equilibrium temperature from Planck 2015, in unit of K, FWHM = 60'
        dust_temperature = healpy.read_map(root_dir+"/PR2-2015/COM_CompMap_dust-commander_0256_R2.00.fits", h=False, field=3)
        # Reference frequency for single component model
        freq_d = 545
        # part1 is numerator
        part1 = numpy.zeros(12*256**2, dtype=numpy.float64)
        # part2 is denominator
        part2 = 0
        # Transmission spectrum of HFI detectors
        freq_array = eval("freq_"+freq_str)
        trans_array= eval("trans_"+freq_str)
        for i in range(freq_low_limit(freq_str), freq_high_limit(freq_str), step):
            freq = freq_array[i]
            delta_freq = 0.5*(freq_array[i+step] - freq_array[i-step])
            trans = trans_array[i]
            x_trans_freq = h_Planck*(1e9*freq)/k_Boltzman/dust_temperature
            Black_body_spectrum = 2*h_Planck*(1e9*freq)**3/c_light**2/(numpy.exp(x_trans_freq)-1)
            part1 = part1 + (freq/freq_d)**dust_index * Black_body_spectrum * trans*delta_freq
            if unit == "muK_CMB":
                part2 = part2 + b_prime_CMB(freq)*trans*delta_freq
            elif unit == "MJy_sr":
                part2 = part2 + (freq/float(freq_str))**(-1)*trans*delta_freq
        x_trans_d = h_Planck*(1e9*freq_d)/k_Boltzman/dust_temperature
        Black_body_spectrum_d = 2*h_Planck*(1e9*freq_d)**3/c_light**2/(numpy.exp(x_trans_d)-1)
        part1 = part1 / Black_body_spectrum_d
        if unit == "muK_CMB":
            # From K_CMB to muK_CMB
            dust_map = 1e6*map_ref*part1/part2
        elif unit == "MJy_sr":
            # From W m^-2 Hz^-1 sr^-1 to MJy/sr
            dust_map = 1e-6*1e26*map_ref*part1/part2

    elif model == "Meisner":
    # two-component MBB model from ApJ, 798:88, 2015, http://planck.skymaps.info/
        # part1 is numerator
        part1 = numpy.zeros(12*2048**2, dtype=numpy.float64)
        # part2 is denominator
        part2 = 0
        # Transmission spectrum of HFI detectors
        freq_array = eval("freq_"+freq_str)
        trans_array= eval("trans_"+freq_str)
        for i in range(freq_low_limit(freq_str), freq_high_limit(freq_str), step):
            freq = freq_array[i]
            delta_freq = 0.5*(freq_array[i+step] - freq_array[i-step])
            trans = trans_array[i]
            dust_intensity = Dust_util_2comp.getval_2comp(nu=float(freq))
            part1 = part1 + dust_intensity*trans*delta_freq
            if unit == "muK_CMB":
                part2 = part2 + b_prime_CMB(freq)*trans*delta_freq
            elif unit == "MJy_sr":
                part2 = part2 + (freq/float(freq_str))**(-1)*trans*delta_freq
        part1 = healpy.reorder(part1, n2r=True)
        if unit == "muK_CMB":
            # From K_CMB to muK_CMB
            dust_map = 1e6*part1/part2
        elif unit == "MJy_sr":
            # From W m^-2 Hz^-1 sr^-1 to MJy/sr
            dust_map = 1e-6*1e26*part1/part2

    elif model == "DL07":
        DL07_freq_dict = {
            "353" : 2,
            "545" : 1,
            "857" : 0
        }
        DL07_field = DL07_freq_dict.get(freq_str)
        # MJy/sr at 353, 545, 857 GHz. FWHM = 5 arcmin, A&A 586, A132 (2016), P26, Appendix D
        dust_map = healpy.read_map(root_dir+"/Planck_DL07/COM_CompMap_Dust-DL07-ModelFluxes_2048_R2.00.fits", field = DL07_field)

    elif model == "GNILC":
        # MJy/sr at 353, 545, 857 GHz. FWHM = 5 arcmin
        dust_map = healpy.read_map(root_dir+"/PR2-2015/COM_CompMap_Dust-GNILC-F"+freq_str+"_2048_R2.00.fits")

    elif model == "SRoll":
        # K_CMB at 100, 143, 217, 353, 545 GHz, MJy/sr at 857 GHz. nside = 128 smoothed at 1 degree.
        # https://sroll20.ias.u-psud.fr/sroll22_dustmodel.html
        dust_map = healpy.read_map(root_dir+"/SRoll/"+freq_str+"GHz_DUSTMODEL_INTENSITY.fits")
        # Unit conversion
        factor_SI = 0
        factor_CMB= 0
        # Transmission spectrum of HFI detectors
        freq_array = eval("freq_"+freq_str)
        trans_array= eval("trans_"+freq_str)
        for i in range(freq_low_limit(freq_str), freq_high_limit(freq_str), step):
            freq = freq_array[i]
            delta_freq = 0.5*(freq_array[i+step] - freq_array[i-step])
            trans = trans_array[i]
            factor_SI = factor_SI + (freq/float(freq_str))**(-1)*trans*delta_freq
            factor_CMB= factor_CMB+ b_prime_CMB(freq)*trans*delta_freq
        factor_SI = 1/factor_SI
        factor_CMB= 1/factor_CMB
        if unit == "MJy_sr" and (freq_str == "100" or freq_str == "143" or freq_str == "217" or freq_str == "353" or freq_str == "545"):
            # from K_CMB to MJy/sr
            dust_map =  1e20*factor_SI/factor_CMB*dust_map
        if unit == "muK_CMB" and freq_str == "857":
            # from MJy/sr to K_CMB
            dust_map = 1e-20*factor_CMB/factor_SI*dust_map
        if unit == "muK_CMB":
            dust_map = 1e6*dust_map

    elif model == "HD17-d5" or model == "HD17-d7" or model == "HD17-d8" or model == "3D":
        # in unit of muK_RJ
        # model=hensley_draine_2017
        if model == "HD17-d5":
            sky = pysm3.Sky(nside=512, preset_strings=["d5"])
        elif model == "HD17-d7":
            sky = pysm3.Sky(nside=512, preset_strings=["d7"])
        elif model == "HD17-d8":
            sky = pysm3.Sky(nside=512, preset_strings=["d8"])
        elif model == "3D":
            sky = pysm3.Sky(nside=512, preset_strings=["d12"])
        # part1 is numerator
        part1 = numpy.zeros(12*512**2, dtype=numpy.float64)
        # part2 is denominator
        part2 = 0
        # Transmission spectrum of HFI detectors
        freq_array = eval("freq_"+freq_str)
        trans_array= eval("trans_"+freq_str)
        for i in range(freq_low_limit(freq_str), freq_high_limit(freq_str), step):
            freq = freq_array[i]
            delta_freq = 0.5*(freq_array[i+step] - freq_array[i-step])
            trans = trans_array[i]
            dust_intensity = sky.get_emission(freq * pysm3.units.GHz) # in unit of muK_RJ
            dust_intensity = pysm3.apply_smoothing_and_coord_transform(dust_intensity, fwhm=20/60*pysm3.units.deg)
            dust_intensity = dust_intensity.value
            dust_intensity = dust_intensity[0]
            dust_intensity = b_prime_RJ(freq) * 1e-6 * dust_intensity # in SI unit
            part1 = part1 + dust_intensity*trans*delta_freq
            if unit == "muK_CMB":
                part2 = part2 + b_prime_CMB(freq)*trans*delta_freq
            elif unit == "MJy_sr":
                part2 = part2 + (freq/float(freq_str))**(-1)*trans*delta_freq
        if unit == "muK_CMB":
            # From K_CMB to muK_CMB
            dust_map = 1e6*part1/part2
        elif unit == "MJy_sr":
            # From W m^-2 Hz^-1 sr^-1 to MJy/sr
            dust_map = 1e-6*1e26*part1/part2
        dust_map = rebeam(dust_map, Nside)

    # Re-beam to FWHM = 1 degree
    beam_dict = {
        "Planck13"  : 5,
        "Planck15-C": 60,
        "Planck15-G": 5,
        "Irfan19"   : 5,
        "Meisner"   : 6.1,
        "DL07"      : 5,
        "GNILC"     : 5,
        "SRoll"     : 60,
        "HD17-d5"   : 20, # pysm3.apply_smoothing_and_coord_transform(dust_intensity, fwhm=20/60*pysm3.units.deg)
        "HD17-d7"   : 20,
        "HD17-d8"   : 20,
        "3D"        : 20
    }

    beam1 = healpy.gauss_beam(fwhm=beam_dict.get(model)/60*numpy.pi/180, pol=False)
    beam2 = healpy.gauss_beam(fwhm=                  60/60*numpy.pi/180, pol=False)
    alm_dust = healpy.map2alm(dust_map)
    alm_dust = healpy.smoothalm(alm_dust, beam_window=beam2/beam1)
    dust_map = healpy.alm2map(alm_dust, nside=Nside)
    return dust_map
#--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*
# Anomalous Microwave Emission / Spinning Dust

# Extension = "SPINNING-DUST-TEMP" is the template for AME (Intensity VS Frequency)
AME_fits = astropy.io.fits.open(root_dir+"/PR2-2015/COM_CompMap_AME-commander_0256_R2.00.fits")
# Reference Frequency, in GHz
AME_freq_0_1 = 22.8
AME_freq_0_2 = 41.0
AME_freq_p_0 = 30.0
AME_freq_p_1 = healpy.read_map(root_dir+"/PR2-2015/COM_CompMap_AME-commander_0256_R2.00.fits", hdu=1, field=3)
AME_freq_p_2 = 33.35

AME_I_1 = healpy.read_map(root_dir+"/PR2-2015/COM_CompMap_AME-commander_0256_R2.00.fits", hdu=1, field=0) # FWHM = 60 arcmin, in muK_RJ
AME_I_2 = healpy.read_map(root_dir+"/PR2-2015/COM_CompMap_AME-commander_0256_R2.00.fits", hdu=2, field=0) # FWHM = 60 arcmin, in muK_RJ
# From muK_RJ to SI unit
AME_I_1 = 1e-6 * b_prime_RJ(AME_freq_0_1) * AME_I_1
AME_I_2 = 1e-6 * b_prime_RJ(AME_freq_0_2) * AME_I_2

AME_freq      = AME_fits["SPINNING-DUST-TEMP"].data.field("NU")   # GHz
AME_intensity = AME_fits["SPINNING-DUST-TEMP"].data.field("J_NU") # Jy sr^-1 cm^2 / H
# AME_template is the spline interpolation obtained based on AME_freq and AME_intensity, representing the relationship between AME intensity and frequency (in log-log).
AME_template    = scipy.interpolate.splrep(numpy.log10(AME_freq), numpy.log10(AME_intensity), k=2)
# AME_intensity_1 is the intensity of first  component at AME_freq_0_1, that is, $f_sd(\nu_0^1 \nu_{p0} / \nu_p^1)$ in A&A,594,A10(2016), Table 4
# AME_intensity_2 is the intensity of second component at AME_freq_0_2, that is, $f_sd(\nu_0^2 \nu_{p0} / \nu_p^2)$ in A&A,594,A10(2016), Table 4
AME_intensity_1 = 10**scipy.interpolate.splev(numpy.log10(AME_freq_0_1 * AME_freq_p_0 / AME_freq_p_1), AME_template)
AME_intensity_2 = 10**scipy.interpolate.splev(numpy.log10(AME_freq_0_2 * AME_freq_p_0 / AME_freq_p_2), AME_template)

def reading_AME(freq_str, unit, Nside):
    Nside = int(Nside)
    # nside = 256 for COM_CompMap_AME-commander_0256_R2.00.fits
    comp1_part1 = numpy.zeros(12*256**2, dtype=numpy.float64)
    comp2_part1 = numpy.zeros(12*256**2, dtype=numpy.float64)
    part2 = 0
    # Transmission spectrum of detectors
    freq_array = eval("freq_"+freq_str)
    trans_array= eval("trans_"+freq_str)
    for i in range(freq_low_limit(freq_str), freq_high_limit(freq_str), 1):
        freq = freq_array[i]
        delta_freq = 0.5*(freq_array[i+1] - freq_array[i-1])
        trans = trans_array[i]
        # part1 is numerator
        # $f_sd(\nu \nu_{p0} / \nu_p^1)$ in A&A,594,A10(2016), Table 4
        comp1_part1 = comp1_part1 + 10**scipy.interpolate.splev(numpy.log10(freq * AME_freq_p_0 / AME_freq_p_1), AME_template) * trans * delta_freq
        # $f_sd(\nu \nu_{p0} / \nu_p^2)$ in A&A,594,A10(2016), Table 4
        comp2_part1 = comp2_part1 + 10**scipy.interpolate.splev(numpy.log10(freq * AME_freq_p_0 / AME_freq_p_2), AME_template) * trans * delta_freq
        if unit == "muK_CMB":
            # part2 is denominator
            part2 = part2 + b_prime_CMB(freq)*trans*delta_freq
        elif unit == "MJy_sr":
            part2 = part2 + (freq/float(freq_str))**(-1)*trans*delta_freq
    map_AME = AME_I_1 * comp1_part1 / AME_intensity_1 / part2 + AME_I_2 * comp2_part1 / AME_intensity_2 / part2
    if unit == "muK_CMB":
        # From K_CMB to muK_CMB
        map_AME = 1e6 * map_AME
    elif unit == "MJy_sr":
        # From W m^-2 Hz^-1 sr^-1 to MJy/sr
        map_AME = 1e-6 * 1e26 * map_AME
    # Up_grade from nside=256 to nside=Nside
#    map_AME = healpy.ud_grade(map_AME, Nside)
    map_AME = rebeam(map_AME, Nside)
    return map_AME


#--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*
# Function Mosaic() is for statistics
def Mosaic(component, name, freq_str_1, freq_str_2, unit, nside1, nside2, smooth_degree, disk_degree, correlation_threshold):
# nside1 for output map
# nside2 for input map
# smooth_degree is the FWHM for smoothing, in degree
# disk_degree is the radius of the disk for statistics, in degree
# 0-th component of the output is "ratio": R(model) or R'(data)
# 1-st component of the output is "correlation": C_xy(model) or C_{x'y'}(data)
# 2-nd component of the output is used to mark points with a correlation of not less than 80% (only for mid-foreground data)
# name = "data"/"model_XX"
    smooth_degree= float(smooth_degree)
    disk_degree  = float(disk_degree)
    # mid_map_1 and mid_map_2 are (data/model) mid-foreground maps in neighbourinng frequency channels
    mid_map_1 = healpy.read_map(root_dir+"/results_"+component+"/"+component+"_"+name+"_"+unit+"_"+freq_str_1+"GHz.fits", nest=False)
    mid_map_2 = healpy.read_map(root_dir+"/results_"+component+"/"+component+"_"+name+"_"+unit+"_"+freq_str_2+"GHz.fits", nest=False)
    # Smooth mid_map_1 with beam = smooth_degree
    # Smooth mid_map_2 with beam = smooth_degree
    beam1     = healpy.gauss_beam(fwhm=            1*numpy.pi/180, pol=False)
    beam2     = healpy.gauss_beam(fwhm=smooth_degree*numpy.pi/180, pol=False)
    alm_mid   = healpy.map2alm(mid_map_1)
    alm_mid   = healpy.smoothalm(alm_mid, beam_window=beam2/beam1)
    mid_map_1 = healpy.alm2map(alm_mid, nside=nside2)
    alm_mid   = healpy.map2alm(mid_map_2)
    alm_mid   = healpy.smoothalm(alm_mid, beam_window=beam2/beam1)
    mid_map_2 = healpy.alm2map(alm_mid, nside=nside2)
    mask = Statistic_mask(nside2, 2, 1.5, 1)
    # points in ipix_disc form a sample set
    # ipix_disc with anchor point as the center and disk_degree as the radius
    # Initial values are set to  - 1.63750e+30, bad value for healpix, and will be masked by healpy.mollview.
    anchor_point    = numpy.zeros((3, 12*nside1**2), dtype=numpy.float64)
    anchor_point[:] = healpy.UNSEEN
    func = lambda theta : 2*numpy.pi*numpy.sin(theta)
    # "area" is the area of ipix_disc
    area = scipy.integrate.quad(func, 0, numpy.radians(disk_degree))[0]
    # "threshold" is the number of data points contained in 30% of the ipix_disc
    threshold = 0.3*12*nside2**2*area/(4*numpy.pi)
    # for ipix in anchor_point
    for ipix in range(0, 12*nside1**2, 1):
        angular   = healpy.pix2ang(nside1, ipix, nest=False, lonlat=False)
        vector    = healpy.ang2vec(angular[0], angular[1])
        ipix_disc = healpy.query_disc(nside=nside2, vec=vector, radius=numpy.radians(disk_degree))
        ipix_num  = numpy.where(mask[ipix_disc]==1)[0]
        # ipix_area is the un-masked ipix_disc, area of the valid region
        ipix_area = ipix_disc[ipix_num]
        # if the area of valid region larger than 30% of ipix_disc, then the anchor point is valid
        if ipix_area.shape[0] > threshold:
            cov_matrix = numpy.cov(mid_map_1[ipix_area], mid_map_2[ipix_area], ddof=1)
            # the divisor used in the calculation is N - ddof, where N represents the number of elements.
            # numpy.cov(A, B)[0][0]: variance of data set A
            # numpy.cov(A, B)[1][1]: variance of data set B
            # numpy.cov(A, B)[0][1], numpy.cov(A, B)[1][0]: covariance of data set A and B
            # ratio is R(model) or R'(data)
            ratio = cov_matrix[0][1] / cov_matrix[0][0]
            # "correlation" is the correlation between mid_map_1 and mid_map_2
            correlation = cov_matrix[0][1] / numpy.sqrt(cov_matrix[0][0]*cov_matrix[1][1])
            anchor_point[0][ipix] = ratio
            anchor_point[1][ipix] = correlation
            if name == "data":
                if correlation >= correlation_threshold:
                    # mark points with a correlation of not less than correlation_threshold (only for data maps)
                    anchor_point[2][ipix] = 1
    numpy.save(root_dir+"/results_"+component+"/scatter_points_"+name+"_"+str(freq_str_1)+"GHz_"+str(freq_str_2)+"GHz_"+unit+".npy", anchor_point)
    return 1


def component_map(component, name, freq_str, unit):
    mask_map_030 = healpy.read_map(root_dir+"/PR2-2015/LFI_Mask_PointSrc_2048_R2.00.fits", h=False, hdu=1, field=0)
    mask_map_044 = healpy.read_map(root_dir+"/PR2-2015/LFI_Mask_PointSrc_2048_R2.00.fits", h=False, hdu=2, field=0)
    mask_map_070 = healpy.read_map(root_dir+"/PR2-2015/LFI_Mask_PointSrc_2048_R2.00.fits", h=False, hdu=3, field=0)
    mask_map_100 = healpy.read_map(root_dir+"/PR2-2015/HFI_Mask_PointSrc_2048_R2.00.fits", h=False, hdu=1, field=0)
    mask_map_143 = healpy.read_map(root_dir+"/PR2-2015/HFI_Mask_PointSrc_2048_R2.00.fits", h=False, hdu=1, field=1)
    mask_map_217 = healpy.read_map(root_dir+"/PR2-2015/HFI_Mask_PointSrc_2048_R2.00.fits", h=False, hdu=1, field=2)
    mask_map_353 = healpy.read_map(root_dir+"/PR2-2015/HFI_Mask_PointSrc_2048_R2.00.fits", h=False, hdu=1, field=3)
    mask_map_545 = healpy.read_map(root_dir+"/PR2-2015/HFI_Mask_PointSrc_2048_R2.00.fits", h=False, hdu=1, field=4)
    mask_map_857 = healpy.read_map(root_dir+"/PR2-2015/HFI_Mask_PointSrc_2048_R2.00.fits", h=False, hdu=1, field=5)

    component_map = healpy.read_map(root_dir+"/results_"+component+"/"+component+"_"+name+"_"+unit+"_"+freq_str+"GHz.fits")
    Nside = healpy.get_nside(component_map)
    if component == "residual":
        mask_map = Statistic_mask(Nside, 2, 1.5, 1)
    else:
        mask_map = mask_map_030 * mask_map_044 * mask_map_070 * mask_map_100 * mask_map_143 * mask_map_217 * mask_map_353 * mask_map_545 * mask_map_857
        mask_map = healpy.ud_grade(mask_map, Nside)
    if unit == "muK_CMB":
        if freq_str == "545":
            component_map = 1e-3 * component_map
            unit_str = "$10^3\\times\\mu\\mathrm{K_{CMB}}$"
        elif freq_str == "857":
            component_map = 1e-6 * component_map
            unit_str = "$10^6\\times\\mu\\mathrm{K_{CMB}}$"
        else:
            unit_str = "$\\mu\\mathrm{K_{CMB}}$"
    elif unit == "MJy_sr":
        unit_str = "$\\mathrm{MJy\,sr^{-1}}$"
    component_map[mask_map == 0] = healpy.UNSEEN
    component_name = {
        "synch"     : "Synchrotron",
        "AME"       : "AME",
        "free_free" : "Free-free",
        "residual"  : "Residual",
        "dust"      : "Thermal dust"
    }
    prefix = component_name.get(component, "Unknown")
    if name == "data":
        title = f"{prefix} data map\n{freq_str}GHz, {unit_str}"
    else:
        title = f"{prefix} model {name.lstrip('model_')} map\n{freq_str}GHz, {unit_str}"
    MIN = numpy.percentile(component_map[component_map != healpy.UNSEEN], 10)
    MAX = numpy.percentile(component_map[component_map != healpy.UNSEEN], 90)
    if abs(MIN) < 0.001:
        MIN = 0
    matplotlib.pyplot.figure(dpi=200, figsize=(4, 3))
#    healpy.projview(component_map, cmap="Spectral_r", unit=unit_str, norm="none", min=0, max=MAX, format="%.4g", extend="neither", hold=True)
    healpy.mollview(component_map, cmap="Spectral_r", norm="none", min=MIN, max=MAX, format="%.3g", nlocs=3, bgcolor="white", hold=True)
    matplotlib.pyplot.title(title, fontsize=11)
    matplotlib.pyplot.savefig(root_dir+"/figure_"+component+"/"+component+"_map_"+name+"_"+unit+"_"+freq_str+"GHz.pdf", bbox_inches="tight")
    matplotlib.pyplot.close()

# Sky maps for R(model), R'(data), C(model), C'(data), where C is correlation
def map_plot(component, data_name, model_name, freq_name, unit):
    # freq_name = 01: 023*/ 030 GHz
    # freq_name = 02: 030 / 033*GHz
    # freq_name = 03: 033*/ 041*GHz
    # freq_name = 04: 041*/ 044 GHz
    # freq_name = 05: 044 / 061*GHz
    # freq_name = 06: 030 / 044 GHz
    # freq_name = 07: 061*/ 070 GHz
    # freq_name = 08: 044 / 070 GHz
    # freq_name = 09: 070 / 100 GHz
    # freq_name = 10: 100 / 143 GHz
    # freq_name = 11: 143 / 217 GHz
    # freq_name = 12: 217 / 353 GHz
    # freq_name = 13: 353 / 545 GHz
    # freq_name = 14: 545 / 857 GHz
    freq_info = {
        "01": ("023GHz", "030GHz", "$23^*$ - 30 GHz"),
        "02": ("030GHz", "033GHz", "30 - $33^*$ GHz"),
        "03": ("033GHz", "041GHz", "$33^*$ - $41^*$ GHz"),
        "04": ("041GHz", "044GHz", "$41^*$ - 44 GHz"),
        "05": ("044GHz", "061GHz", "44 - $61^*$ GHz"),
        "06": ("030GHz", "044GHz", "30 - 44 GHz"),
        "07": ("061GHz", "070GHz", "$61^*$ - 70 GHz"),
        "08": ("044GHz", "070GHz", "44 - 70 GHz"),
        "09": ("070GHz", "100GHz", "70 - 100 GHz"),
        "10": ("100GHz", "143GHz", "100 - 143 GHz"),
        "11": ("143GHz", "217GHz", "143 - 217 GHz"),
        "12": ("217GHz", "353GHz", "217 - 353 GHz"),
        "13": ("353GHz", "545GHz", "353 - 545 GHz"),
        "14": ("545GHz", "857GHz", "545 - 857 GHz")
    }

    if freq_name in freq_info:
        f1, f2, str_name = freq_info[freq_name]
        mosaic_model = numpy.load(root_dir+"/results_"+component+"/scatter_points_"+model_name+"_"+f1+"_"+f2+"_"+unit+".npy")
        mosaic_data  = numpy.load(root_dir+"/results_"+component+"/scatter_points_"+data_name +"_"+f1+"_"+f2+"_"+unit+".npy")
    if unit == "muK_CMB":
        unit_str = "$\\mu\\mathrm{K_{CMB}}$"
    elif unit=="MJy_sr":
        unit_str = "$\\mathrm{MJy\,sr^{-1}}$"

    # R (model)
    # invalid_points: marked with - 1.63750e+30
    valid_points = numpy.where(mosaic_model[0] > -10000000)[0]
    MIN = numpy.percentile(mosaic_model[0][valid_points],10)
    MAX = numpy.percentile(mosaic_model[0][valid_points],90)
    matplotlib.pyplot.figure(dpi=100, figsize=(3.6,3.6))
    healpy.mollview(mosaic_model[0], nlocs=3, nest=False, min=MIN, max=MAX, cmap="Spectral_r", format="%.4g", hold=True)
    matplotlib.pyplot.title("$R$ of "+str_name+" ("+unit_str+")\n model = "+model_name.lstrip("model_"), fontsize=13)
    matplotlib.pyplot.savefig(root_dir+"/figure_"+component+"/R_"+model_name+"_"+unit+"_"+freq_name+".pdf", bbox_inches="tight")
    matplotlib.pyplot.close()

    # R' (data)
    valid_points = numpy.where(mosaic_data[0] > -10000000)[0]
    MIN = numpy.percentile(mosaic_data[0][valid_points],10)
    MAX = numpy.percentile(mosaic_data[0][valid_points],90)
    matplotlib.pyplot.figure(dpi=100, figsize=(3.6,3.6))
    healpy.mollview(mosaic_data[0], nlocs=3, nest=False, min=MIN, max=MAX, cmap="Spectral_r", format="%.4g", hold=True)
    matplotlib.pyplot.title("$R'$ of "+str_name+" ("+unit_str+")\n data", fontsize=13)
    matplotlib.pyplot.savefig(root_dir+"/figure_"+component+"/R_"+data_name+"_"+unit+"_"+freq_name+".pdf", bbox_inches="tight")
    matplotlib.pyplot.close()

    # Correlation between two neighbouring frequency channels (data)
    valid_points = numpy.where(mosaic_data[1] > -10000000)[0]
    MIN = numpy.percentile(mosaic_data[1][valid_points],10)
    MAX = numpy.percentile(mosaic_data[1][valid_points],90)
    matplotlib.pyplot.figure(dpi=100, figsize=(3.6,3.6))
    healpy.mollview(mosaic_data[1], nlocs=3, nest=False, min=MIN, max=MAX, cmap="Spectral_r", format="%.3g", hold=True)
    matplotlib.pyplot.title("$C'$ of "+str_name, fontsize=13)
    matplotlib.pyplot.savefig(root_dir+"/figure_"+component+"/C_"+data_name+"_"+unit+"_"+freq_name+".pdf", bbox_inches="tight")
    matplotlib.pyplot.close()
    return 1

# Regions the points in the scatter plots
def region(component, data_name, freq_name, unit):
    # freq_name = 01: 023*/ 030 GHz
    # freq_name = 02: 030 / 033*GHz
    # freq_name = 03: 033*/ 041*GHz
    # freq_name = 04: 041*/ 044 GHz
    # freq_name = 05: 044 / 061*GHz
    # freq_name = 06: 030 / 044 GHz
    # freq_name = 07: 061*/ 070 GHz
    # freq_name = 08: 044 / 070 GHz
    # freq_name = 09: 070 / 100 GHz
    # freq_name = 10: 100 / 143 GHz
    # freq_name = 11: 143 / 217 GHz
    # freq_name = 12: 217 / 353 GHz
    # freq_name = 13: 353 / 545 GHz
    # freq_name = 14: 545 / 857 GHz
    freq_info = {
        "01": ("023GHz", "030GHz", "$23^*$ - 30 GHz"),
        "02": ("030GHz", "033GHz", "30 - $33^*$ GHz"),
        "03": ("033GHz", "041GHz", "$33^*$ - $41^*$ GHz"),
        "04": ("041GHz", "044GHz", "$41^*$ - 44 GHz"),
        "05": ("044GHz", "061GHz", "44 - $61^*$ GHz"),
        "06": ("030GHz", "044GHz", "30 - 44 GHz"),
        "07": ("061GHz", "070GHz", "$61^*$ - 70 GHz"),
        "08": ("044GHz", "070GHz", "44 - 70 GHz"),
        "09": ("070GHz", "100GHz", "70 - 100 GHz"),
        "10": ("100GHz", "143GHz", "100 - 143 GHz"),
        "11": ("143GHz", "217GHz", "143 - 217 GHz"),
        "12": ("217GHz", "353GHz", "217 - 353 GHz"),
        "13": ("353GHz", "545GHz", "353 - 545 GHz"),
        "14": ("545GHz", "857GHz", "545 - 857 GHz")
    }

    if freq_name in freq_info:
        f1, f2, str_name = freq_info[freq_name]
        mosaic_data  = numpy.load(root_dir+"/results_"+component+"/scatter_points_"+data_name +"_"+f1+"_"+f2+"_"+unit+".npy")

    # total_pixel_list is index of the points without masked (valid region)
    total_pixel_list = numpy.where(mosaic_data[1]>0)[0]
    # pixel_list is index of the points plotted in the scatter plot with C'>correlation_threshold (reliable region, in green)
    pixel_list = numpy.where(mosaic_data[2]==1)[0]
    # fraction_reliable = reliable region / valid region
    fraction_reliable = numpy.array(pixel_list).shape[0]/numpy.array(total_pixel_list).shape[0]
    matplotlib.pyplot.figure(dpi=100, figsize=(4,4))
    healpy.mollview(mosaic_data[2], cmap="jet", min=0, max=2, cbar=False, hold=True)
    matplotlib.pyplot.title("Reliable region, "+str_name+"\n $f_\mathrm{rel} = $"+"{:.1f}".format(100*fraction_reliable)+"%", fontsize=13)
    matplotlib.pyplot.savefig(root_dir+"/figure_"+component+"/Region_"+data_name+"_"+unit+"_"+freq_name+".pdf", bbox_inches="tight")
    matplotlib.pyplot.close()

# Scatter plots for dust (data) and the sky map for the points in the scatter plots
# If threshold = 5, data points with the lowest 5% probability density values will be excluded.
# If linea_fit = False, No linear fitting. 
def scatter_plot(component, data_name, model_name, freq_name, unit, threshold, linear_fit):
    # freq_name = 01: 023*/ 030 GHz
    # freq_name = 02: 030 / 033*GHz
    # freq_name = 03: 033*/ 041*GHz
    # freq_name = 04: 041*/ 044 GHz
    # freq_name = 05: 044 / 061*GHz
    # freq_name = 06: 030 / 044 GHz
    # freq_name = 07: 061*/ 070 GHz
    # freq_name = 08: 044 / 070 GHz
    # freq_name = 09: 070 / 100 GHz
    # freq_name = 10: 100 / 143 GHz
    # freq_name = 11: 143 / 217 GHz
    # freq_name = 12: 217 / 353 GHz
    # freq_name = 13: 353 / 545 GHz
    # freq_name = 14: 545 / 857 GHz
    freq_info = {
        "01": ("023GHz", "030GHz", "$23^*$ - 30 GHz"),
        "02": ("030GHz", "033GHz", "30 - $33^*$ GHz"),
        "03": ("033GHz", "041GHz", "$33^*$ - $41^*$ GHz"),
        "04": ("041GHz", "044GHz", "$41^*$ - 44 GHz"),
        "05": ("044GHz", "061GHz", "44 - $61^*$ GHz"),
        "06": ("030GHz", "044GHz", "30 - 44 GHz"),
        "07": ("061GHz", "070GHz", "$61^*$ - 70 GHz"),
        "08": ("044GHz", "070GHz", "44 - 70 GHz"),
        "09": ("070GHz", "100GHz", "70 - 100 GHz"),
        "10": ("100GHz", "143GHz", "100 - 143 GHz"),
        "11": ("143GHz", "217GHz", "143 - 217 GHz"),
        "12": ("217GHz", "353GHz", "217 - 353 GHz"),
        "13": ("353GHz", "545GHz", "353 - 545 GHz"),
        "14": ("545GHz", "857GHz", "545 - 857 GHz")
    }

    if freq_name in freq_info:
        f1, f2, str_name = freq_info[freq_name]
        mosaic_model = numpy.load(root_dir+"/results_"+component+"/scatter_points_"+model_name+"_"+f1+"_"+f2+"_"+unit+".npy")
        mosaic_data  = numpy.load(root_dir+"/results_"+component+"/scatter_points_"+data_name +"_"+f1+"_"+f2+"_"+unit+".npy")
    if unit == "muK_CMB":
        unit_str = "$\\mu\\mathrm{K_{CMB}}$"
    elif unit=="MJy_sr":
        unit_str = "$\\mathrm{MJy\,sr^{-1}}$"

    pixel_list = numpy.where(mosaic_data[2]==1)[0]
    # R (model) from the points in reliable region
    model= mosaic_model[0][pixel_list]
    # R' (data) from the points in reliable region
    data = mosaic_data[0][pixel_list]
    x = numpy.arange(0, 2000, 0.1)
    scatters = numpy.array([model, data])
    density  = scipy.stats.gaussian_kde(scatters)(scatters)
    scatters = scatters[:, density >= numpy.percentile(density, threshold)]
    model    = scatters[0]
    data     = scatters[1]

    slope, intercept = numpy.polyfit(model, data, 1)
    data_fit         = slope * x + intercept
    line_name        = "y = "+"{:+.3f}".format(slope)+" x "+"{:+.3f}".format(intercept)

    matplotlib.pyplot.figure(dpi=150, figsize=(7,7))
    matplotlib.pyplot.scatter(model, data, s=2, marker="+", color="black")
    matplotlib.pyplot.plot(x, x, label="$R' = R$", linewidth=3)
    # Linear fitting of R' vs R
    if linear_fit == False:
        pass
    else:
        matplotlib.pyplot.plot(x, data_fit, label=line_name, linewidth=2)
    MIN = min(min(model), min(data))
    MAX = max(max(model), max(data))
    matplotlib.pyplot.xlim(xmin=MIN-0.1*(MAX-MIN), xmax=MIN+1.0*(MAX-MIN))
    matplotlib.pyplot.ylim(ymin=MIN-0.0*(MAX-MIN), ymax=MIN+1.2*(MAX-MIN))
    if model_name == "model":
        matplotlib.pyplot.xlabel(r"$R$ (model)", fontsize=22)
    else:
        matplotlib.pyplot.xlabel(r"$R$ (model = "+model_name.lstrip("model_")+")", fontsize=22)
    matplotlib.pyplot.ylabel(r"$R'$ (data)", fontsize=21)
    matplotlib.pyplot.legend(fontsize=21, loc="upper left")
    matplotlib.pyplot.title(str_name+" ("+unit_str+")", fontsize=23)
    matplotlib.pyplot.savefig(root_dir+"/figure_"+component+"/Scatter_"+model_name+"_"+unit+"_"+freq_name+".png", bbox_inches="tight")
    matplotlib.pyplot.close()
