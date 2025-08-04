# Free-free intensity sky maps based on Planck data and model

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
import public_function

root_dir = os.path.abspath('.')

#--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*
# Free-free intensity sky maps based on Planck observed data, with color correction
def reading_data_free_free(freq_str, unit, Nside):
    Nside = int(Nside)
    rebeamed_inpainted_sky_map = healpy.read_map(root_dir+"/results/inpainted_rebeamed_map_"+freq_str+"GHz_Nside_"+str(Nside)+".fits")
    # Unit conversion
    factor_SI = 0
    factor_CMB= 0
    factor_RJ = 0
    # Transmission spectrum of detectors
    freq_array = eval("public_function.freq_"+freq_str)
    trans_array= eval("public_function.trans_"+freq_str)
    for i in range(public_function.freq_low_limit(freq_str), public_function.freq_high_limit(freq_str), 1):
        freq = freq_array[i]
        delta_freq = 0.5*(freq_array[i+1] - freq_array[i-1])
        trans = trans_array[i]
        factor_SI = factor_SI + (freq/float(freq_str))**(-1)     *trans*delta_freq
        factor_CMB= factor_CMB+ public_function.b_prime_CMB(freq)*trans*delta_freq
        factor_RJ = factor_RJ + public_function.b_prime_RJ(freq) *trans*delta_freq
    factor_SI = 1/factor_SI
    factor_CMB= 1/factor_CMB
    factor_RJ = 1/factor_RJ
    if unit == "MJy_sr" and freq_str in {"030", "044", "070", "100", "143", "217", "353"}:
        # from muK_CMB to MJy/sr
        rebeamed_inpainted_sky_map = 1e-6 *1e20*factor_SI /factor_CMB*rebeamed_inpainted_sky_map
    if unit == "muK_CMB" and freq_str in {"545", "857"}:
        # from MJy/sr to muK_CMB
        rebeamed_inpainted_sky_map = 1e-20*1e6 *factor_CMB/factor_SI *rebeamed_inpainted_sky_map
    if unit == "muK_CMB" and freq_str in {"023", "033", "041", "061"}:
        # from mK_RJ to muK_CMB
        rebeamed_inpainted_sky_map = 1e-3 *1e6 *factor_CMB/factor_RJ *rebeamed_inpainted_sky_map
    free_free_map = rebeamed_inpainted_sky_map - public_function.CMB_map(freq_str, unit, Nside) - public_function.reading_dust("Planck15-G", freq_str, unit, Nside) - public_function.reading_synch("Planck", freq_str, unit, Nside) - public_function.reading_AME(freq_str, unit, Nside) - public_function.reading_xline(freq_str, unit, Nside)
    healpy.write_map(root_dir+"/results_free_free/free_free_data_"+unit+"_"+freq_str+"GHz.fits", free_free_map, 
                     nest = False, coord = "G", column_units = unit, overwrite="True", dtype=numpy.float64)
    return 1

#--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*
# Free-free intensity sky maps based on Planck model, with color correction
def reading_model_free_free(freq_str, unit, Nside):
    free_free_map = public_function.reading_ff(freq_str, unit, Nside)
    healpy.write_map(root_dir+"/results_free_free/free_free_model_"+unit+"_"+freq_str+"GHz.fits", free_free_map, 
                     nest = False, coord = "G", column_units = unit, overwrite="True", dtype=numpy.float64)
    return 1

#--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*
Nside = 512

if os.path.exists(root_dir+"/results_free_free/"):
   shutil.rmtree(root_dir+"/results_free_free/")
os.mkdir(root_dir+"/results_free_free/")

print("free-free map from observed data: ")
start_time = time.time()
arguments = [("030", "muK_CMB", Nside), 
             ("044", "muK_CMB", Nside), 
             ("070", "muK_CMB", Nside), 
             ("100", "muK_CMB", Nside), 
             ("143", "muK_CMB", Nside), 
             ("217", "muK_CMB", Nside)]
multiprocessing.Pool(processes=6).starmap(reading_data_free_free, arguments)

for freq_str in ["030", "044", "070", "100", "143", "217"]:
    if os.path.isfile(root_dir+"/results_free_free/free_free_data_muK_CMB_"+freq_str+"GHz.fits"):
        pass
    else:
        reading_data_free_free(freq_str, "muK_CMB", Nside)
end_time = time.time()
print('TIME: ', format((end_time-start_time)/60, '.2f'), 'min')


print("free-free map from model: ")
start_time = time.time()
arguments = [("030","muK_CMB",Nside), 
             ("044","muK_CMB",Nside), 
             ("070","muK_CMB",Nside), 
             ("100","muK_CMB",Nside), 
             ("143","muK_CMB",Nside), 
             ("217","muK_CMB",Nside)]
multiprocessing.Pool(processes=6).starmap(reading_model_free_free, arguments)

for freq_str in ["030", "044", "070", "100", "143", "217"]:
    if os.path.isfile(root_dir+"/results_free_free/free_free_model_muK_CMB_"+freq_str+"GHz.fits"):
        pass
    else:
        reading_model_free_free(freq_str, "muK_CMB", Nside)
end_time = time.time()
print('TIME: ', format((end_time-start_time)/60, '.2f'), 'min')
