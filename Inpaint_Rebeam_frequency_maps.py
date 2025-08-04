# 1. Masking primordial sky maps with point source masks. 
# 2. Iterating 2000 times to in-paint the holes after masking. 
# 3. Rebeam the inpainted maps to FWHM = 1 degree. 

import healpy
import os
import numpy
import scipy
import astropy.io.fits
import time
import shutil
import sys
import multiprocessing
import matplotlib.pyplot

print("Inpain_Rebeam_frequency_maps.py, Starting")
root_dir = os.path.abspath('.')

# File name for each channel
def sky_map_file_name(freq_str):
    file_names = {
        "023": "wmap_band_imap_r9_9yr_K_v5.fits",
        # in mK_RJ
        "030": "LFI_SkyMap_030-BPassCorrected_1024_R3.00_full.fits",
        # in K_cmb
        "033": "wmap_band_imap_r9_9yr_Ka_v5.fits",
        # in mK_RJ
        "041": "wmap_band_imap_r9_9yr_Q_v5.fits",
        # in mK_RJ
        "044": "LFI_SkyMap_044-BPassCorrected_1024_R3.00_full.fits",
        # in K_cmb
        "061": "wmap_band_imap_r9_9yr_V_v5.fits",
        # in mK_RJ
        "070": "LFI_SkyMap_070-BPassCorrected_1024_R3.00_full.fits",
        # in K_cmb
        "100": "HFI_SkyMap_100_2048_R3.01_full.fits",
        # in K_cmb
        "143": "HFI_SkyMap_143_2048_R3.01_full.fits",
        # in K_cmb
        "217": "HFI_SkyMap_217_2048_R3.01_full.fits",
        # in K_cmb
        "353": "HFI_SkyMap_353_2048_R3.01_full.fits",
        # in K_cmb
        "545": "HFI_SkyMap_545_2048_R3.01_full.fits",
        # in MJy/sr
        "857": "HFI_SkyMap_857_2048_R3.01_full.fits"
        # in MJy/sr
    }
    # If freq_str exists in file_names, it returns the corresponding value,
    # If freq_str doesn't exist in file_names, it returns the default value 0.
    return file_names.get(freq_str, 0)

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

#--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*
# inpainting
def inpainted_WMAP_map(freq_str, N_iter):
# Inpaint the points where mask = 0
    # Nside is the nside of WMAP map and mask map
    Nside = 512
    # If phi is not given or None, theta is interpreted as pixel number, otherwise, theta, phi are angular coordinates in radians.
    neighbour_array = healpy.get_all_neighbours(nside=Nside, theta=numpy.arange(0, 12*Nside**2, 1), phi=None).T
    # N_iter is the iteration count
    N_iter = int(N_iter)
    # sky_map is WMAP map
    sky_map = healpy.read_map(root_dir+"/WMAP/"+sky_map_file_name(freq_str), h=False, hdu=1, field=0)
    # mask_map
    if freq_str == "023":
        mask_map = healpy.read_map(root_dir+"/PR2-2015/LFI_Mask_PointSrc_2048_R2.00.fits", h=False, hdu=1, field=0)
    if freq_str == "033":
        mask_map = healpy.read_map(root_dir+"/PR2-2015/LFI_Mask_PointSrc_2048_R2.00.fits", h=False, hdu=1, field=0)
    if freq_str == "041":
        mask_map = healpy.read_map(root_dir+"/PR2-2015/LFI_Mask_PointSrc_2048_R2.00.fits", h=False, hdu=2, field=0)
    if freq_str == "061":
        mask_map = healpy.read_map(root_dir+"/PR2-2015/LFI_Mask_PointSrc_2048_R2.00.fits", h=False, hdu=3, field=0)
    # nside of mask_map: from 2048 to 512
    mask_map = healpy.ud_grade(mask_map, Nside)
    mask_map[mask_map < 0.9] = 0.0
    # values in mask.fits are Boolean
    # Boolean values to floats
    mask_map = mask_map.astype(numpy.float64)

    SkyMap = sky_map*mask_map
    # Filter out the points with mask = 0.
    ipix_mask = numpy.where(mask_map < 0.9)[0]
    # Get the neighboring pixels of the pixels that need to be inpainted
    neighbours = neighbour_array[ipix_mask]
    # valid_neighbours are the valid neighboring pixels
    # Since in the healpix algorithm, a pixel has 8 or 7 neighbours, when a pixel has 7 neighbours, healpy.get_all_neighbours() returns a -1 to mark the invalid neighbour
    valid_neighbours = neighbours >= 0
    # Count the neighboring pixels
    counts = numpy.sum(valid_neighbours, axis=1)

    for ii in range(0, N_iter):
        # Sum the values of the valid neighboring pixels
        sums = numpy.sum(SkyMap[neighbours] * valid_neighbours, axis=1)
        # Calculate the mean value of the neighboring pixels for each masked pixel
        mean_values = sums / counts
        # Update the SkyMap values for the pixels that need to be inpainted
        SkyMap[ipix_mask] = mean_values
    return SkyMap

def inpainted_LFI_map(freq_str, N_iter):
# Inpaint the points where mask = 0
    # Nside is the nside of LFI map and mask map
    Nside = 1024
    # If phi is not given or None, theta is interpreted as pixel number, otherwise, theta, phi are angular coordinates in radians.
    neighbour_array = healpy.get_all_neighbours(nside=Nside, theta=numpy.arange(0, 12*Nside**2, 1), phi=None).T
    # N_iter is the iteration count
    N_iter = int(N_iter)
    # sky_map is LFI map
    sky_map = healpy.read_map(root_dir+"/PR3-2018/"+sky_map_file_name(freq_str), h=False, hdu=1, field=0)
    # mask_map
    mask_map_dict = {
        "030": 1,
        "044": 2,
        "070": 3
    }
    mask_hdu = mask_map_dict.get(freq_str)
    mask_map = healpy.read_map(root_dir+"/PR2-2015/LFI_Mask_PointSrc_2048_R2.00.fits", h=False, hdu=mask_hdu, field=0)
    mask_map = mask_map.astype(numpy.float64)
    # nside of mask_map: from 2048 to 1024
    mask_map = healpy.ud_grade(mask_map, Nside)
    mask_map[mask_map < 0.9] = 0.0
    # values in mask.fits are Boolean
    # Boolean values to floats

    SkyMap = sky_map*mask_map
    # Filter out the points with mask = 0.
    ipix_mask = numpy.where(mask_map < 0.9)[0]
    # Get the neighboring pixels of the pixels that need to be inpainted
    neighbours = neighbour_array[ipix_mask]
    # valid_neighbours are the valid neighboring pixels
    # Since in the healpix algorithm, a pixel has 8 or 7 neighbours, when a pixel has 7 neighbours, healpy.get_all_neighbours() returns a -1 to mark the invalid neighbour
    valid_neighbours = neighbours >= 0
    # Count the neighboring pixels
    counts = numpy.sum(valid_neighbours, axis=1)

    for ii in range(0, N_iter):
        # Sum the values of the valid neighboring pixels
        sums = numpy.sum(SkyMap[neighbours] * valid_neighbours, axis=1)
        # Calculate the mean value of the neighboring pixels for each masked pixel
        mean_values = sums / counts
        # Update the SkyMap values for the pixels that need to be inpainted
        SkyMap[ipix_mask] = mean_values
    # For 030, 044, 070 GHz, output in muK_CMB
    SkyMap = 1e6*SkyMap
    return SkyMap

def inpainted_HFI_map(freq_str, N_iter):
    # Inpaint the points where mask = 0
    # Nside is the nside of HFI map and mask map
    Nside = 2048
    # If phi is not given or None, theta is interpreted as pixel number, otherwise, theta, phi are angular coordinates in radians.
    neighbour_array = healpy.get_all_neighbours(nside=Nside, theta=numpy.arange(0, 12*Nside**2, 1), phi=None).T
    # N_iter is the iteration count
    N_iter = int(N_iter)
    # sky_map is HFI map
    sky_map = healpy.read_map(root_dir+"/PR3-2018/"+sky_map_file_name(freq_str), h=False, hdu=1, field=0)
    # mask_map
    mask_map_dict = {
        "100": 0,
        "143": 1,
        "217": 2,
        "353": 3,
        "545": 4,
        "857": 5
    }
    mask_field = mask_map_dict.get(freq_str)
    mask_map = healpy.read_map(root_dir+"/PR2-2015/HFI_Mask_PointSrc_2048_R2.00.fits", h=False, hdu=1, field=mask_field)
    # values in mask.fits are Boolean
    # Boolean values to floats
    mask_map = mask_map.astype(numpy.float64)

    SkyMap = sky_map*mask_map
    # Filter out the points with mask = 0.
    ipix_mask = numpy.where(mask_map < 0.9)[0]
    # Get the neighboring pixels of the pixels that need to be inpainted
    neighbours = neighbour_array[ipix_mask]
    # valid_neighbours are the valid neighboring pixels
    # Since in the healpix algorithm, a pixel has 8 or 7 neighbours, when a pixel has 7 neighbours, healpy.get_all_neighbours() returns a -1 to mark the invalid neighbour
    valid_neighbours = neighbours >= 0
    # Count the neighboring pixels
    counts = numpy.sum(valid_neighbours, axis=1)

    for ii in range(0, N_iter):
        # Sum the values of the valid neighboring pixels
        sums = numpy.sum(SkyMap[neighbours] * valid_neighbours, axis=1)
        # Calculate the mean value of the neighboring pixels for each masked pixel
        mean_values = sums / counts
        # Update the SkyMap values for the pixels that need to be inpainted
        SkyMap[ipix_mask] = mean_values
    if freq_str in {"100", "143", "217", "353"}:
    # For 100, 143, 217, 353 GHz, output in muK_CMB
        SkyMap = 1e6*SkyMap
    return SkyMap


#--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*
# Smooth and rebeam the inpainted HFI maps to FWHM = 60'
def rebeamed_inpainted_sky_map(freq_str, N_iter, Nside):
    Nside = int(Nside)
    if freq_str in {"023", "033", "041", "061"}:
        SkyMap = inpainted_WMAP_map(freq_str, N_iter)
    elif freq_str in {"030", "044", "070"}:
        SkyMap = inpainted_LFI_map(freq_str, N_iter)
    elif freq_str in {"100", "143", "217", "353", "545", "857"}:
        SkyMap = inpainted_HFI_map(freq_str, N_iter)
    # There are bad values in Planck frequency maps, they are set to - 1.63750e+30. This step is taken to remove bad values.
    bad_pixel = numpy.where(SkyMap<-1e+10)[0]
    if freq_str in {"023", "033", "041", "061"}:
        Nside_temp = 512
    elif freq_str in {"030", "044", "070"}:
        Nside_temp = 1024
    elif freq_str in {"100", "143", "217", "353", "545", "857"}:
        Nside_temp = 2048
    for pixel in bad_pixel:
        neighbour = healpy.get_all_neighbours(nside=Nside_temp, theta=pixel, phi=None)
        SkyMap[pixel] = numpy.mean(SkyMap[neighbour])
    beam1 = healpy.gauss_beam(fwhm = beam_fwhm(freq_str))
    beam2 = healpy.gauss_beam(fwhm = 60/60*numpy.pi/180)
    alm_sky = healpy.map2alm(SkyMap)
    alm_sky = healpy.smoothalm(alm_sky, beam_window=beam2/beam1)
    SkyMap = healpy.alm2map(alm_sky, nside=Nside)

    if freq_str in {"023", "033", "041", "061"}:
    # For 023, 033, 041 GHz, in mK_RJ
        healpy.write_map(
            root_dir+"/results/inpainted_rebeamed_map_"+freq_str+"GHz_Nside_"+str(Nside)+".fits", 
            SkyMap, nest = False, coord = "G", column_units = "mK_RJ", overwrite="True", dtype=numpy.float64)
    elif freq_str in {"030", "044", "070", "100", "143", "217", "353"}:
    # For 030, 044, 070, 100, 143, 217, 353 GHz, in muK_CMB
        healpy.write_map(
            root_dir+"/results/inpainted_rebeamed_map_"+freq_str+"GHz_Nside_"+str(Nside)+".fits", 
            SkyMap, nest = False, coord = "G", column_units = "muK_CMB", overwrite="True", dtype=numpy.float64)
    elif freq_str in {"545", "857"}:
    # For 545, 857 GHz, in MJy/sr
        healpy.write_map(
            root_dir+"/results/inpainted_rebeamed_map_"+freq_str+"GHz_Nside_"+str(Nside)+".fits", 
            SkyMap, nest = False, coord = "G", column_units = "MJy_sr", overwrite="True", dtype=numpy.float64)
    return SkyMap

if os.path.exists(root_dir+"/results/"):
   shutil.rmtree(root_dir+"/results/")
os.mkdir(root_dir+"/results/")

print("Inpainting: ")
Nside = 512
N_iter = int(sys.argv[1])
start_time = time.time()
arguments = [("023", N_iter, Nside), 
             ("030", N_iter, Nside), 
             ("033", N_iter, Nside), 
             ("041", N_iter, Nside), 
             ("044", N_iter, Nside), 
             ("061", N_iter, Nside), 
             ("070", N_iter, Nside), 
             ("100", N_iter, Nside), 
             ("143", N_iter, Nside), 
             ("217", N_iter, Nside), 
             ("353", N_iter, Nside), 
             ("545", N_iter, Nside), 
             ("857", N_iter, Nside)]
multiprocessing.Pool(processes=13).starmap(rebeamed_inpainted_sky_map, arguments)
end_time = time.time()
print("Inpaint_Rebeam_frequency_maps.py, Success!")
print("TIME: ", format( (end_time-start_time)/60, '.2f'),"minutes")
