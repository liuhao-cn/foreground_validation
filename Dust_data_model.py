# Thermal dust intensity sky maps based on Planck data and models

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
# Thermal dust intensity sky maps based on Planck observed data, with color correction
def reading_data_dust(freq_str, unit, Nside):
    Nside = int(Nside)
    rebeamed_inpainted_sky_map = healpy.read_map(root_dir+"/results/inpainted_rebeamed_map_"+freq_str+"GHz_Nside_"+str(Nside)+".fits")
    # Unit conversion
    factor_SI = 0
    factor_CMB= 0
    # Transmission spectrum of HFI detectors
    freq_array = eval("public_function.freq_"+freq_str)
    trans_array= eval("public_function.trans_"+freq_str)
    for i in range(public_function.freq_low_limit(freq_str), public_function.freq_high_limit(freq_str), 1):
        freq = freq_array[i]
        delta_freq = 0.5*(freq_array[i+1] - freq_array[i-1])
        trans = trans_array[i]
        factor_SI = factor_SI + (freq/float(freq_str))**(-1)     *trans*delta_freq
        factor_CMB= factor_CMB+ public_function.b_prime_CMB(freq)*trans*delta_freq
    factor_SI = 1/factor_SI
    factor_CMB= 1/factor_CMB
    if unit == "MJy_sr" and freq_str in {"100", "143", "217", "353"}:
        # from muK_CMB to MJy/sr
        rebeamed_inpainted_sky_map = 1e-6 *1e20*factor_SI /factor_CMB*rebeamed_inpainted_sky_map
    if unit == "muK_CMB" and freq_str in {"545", "857"}:
        # from MJy/sr to muK_CMB
        rebeamed_inpainted_sky_map = 1e-20*1e6 *factor_CMB/factor_SI *rebeamed_inpainted_sky_map
    dust_map = rebeamed_inpainted_sky_map - public_function.CMB_map(freq_str, unit, Nside) - public_function.reading_ff(freq_str, unit, Nside) - public_function.reading_synch("Planck", freq_str, unit, Nside) - public_function.reading_xline(freq_str, unit, Nside) - public_function.reading_AME(freq_str, unit, Nside)
    healpy.write_map(root_dir+"/results_dust/dust_data_"+unit+"_"+freq_str+"GHz.fits", dust_map, 
                     nest = False, coord = "G", column_units = unit, overwrite="True", dtype=numpy.float64)
    return dust_map

#--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*
# Thermal dust intensity sky maps based on models, with color correction
def reading_model_dust(model, freq_str, unit, Nside):
    dust_map = public_function.reading_dust(model, freq_str, unit, Nside)
    healpy.write_map(root_dir+"/results_dust/dust_model_"+model+"_"+unit+"_"+freq_str+"GHz.fits", dust_map, 
                     nest = False, coord = "G", column_units = unit, overwrite="True", dtype=numpy.float64)
    return 1
#--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*

Nside = 512

if os.path.exists(root_dir+"/results_dust/"):
    shutil.rmtree(root_dir+"/results_dust/")
os.mkdir(root_dir+"/results_dust/")

print("Dust map from observed data: ")
start_time = time.time()
arguments = [("100", "muK_CMB", Nside), 
             ("143", "muK_CMB", Nside), 
             ("217", "muK_CMB", Nside), 
             ("353", "muK_CMB", Nside), 
             ("545", "muK_CMB", Nside), 
             ("857", "muK_CMB", Nside),
             ("353", "MJy_sr",  Nside),
             ("545", "MJy_sr",  Nside),
             ("857", "MJy_sr",  Nside)]
multiprocessing.Pool(processes=9).starmap(reading_data_dust, arguments)

for freq_str in ["100", "143", "217", "353", "545", "857"]:
    if os.path.isfile(root_dir+"/results_dust/dust_data_muK_CMB_"+freq_str+"GHz.fits"):
        pass
    else:
        reading_data_dust(freq_str, "muK_CMB", Nside)
end_time = time.time()
print('TIME: ', format((end_time-start_time)/60, '.2f'), 'min')



print("Dust map from model (3D): ")
start_time = time.time()
arguments = [("3D",      "100", "muK_CMB", Nside), 
             ("3D",      "143", "muK_CMB", Nside), 
             ("3D",      "217", "muK_CMB", Nside), 
             ("3D",      "353", "muK_CMB", Nside), 
             ("3D",      "545", "muK_CMB", Nside), 
             ("3D",      "857", "muK_CMB", Nside)]
multiprocessing.Pool(processes=6).starmap(reading_model_dust, arguments)

for freq_str in ["100", "143", "217", "353", "545", "857"]:
    if os.path.isfile(root_dir+"/results_dust/dust_model_3D_muK_CMB_"+freq_str+"GHz.fits"):
        pass
    else:
        reading_model_dust("3D", freq_str, "muK_CMB", Nside)
end_time = time.time()
print('TIME: ', format((end_time-start_time)/60, '.2f'), 'min')



print("Dust map from model (Planck 2013 single MBB model): ")
start_time = time.time()
arguments = [("Planck13", "100", "muK_CMB", Nside), 
             ("Planck13", "143", "muK_CMB", Nside), 
             ("Planck13", "217", "muK_CMB", Nside), 
             ("Planck13", "353", "muK_CMB", Nside), 
             ("Planck13", "545", "muK_CMB", Nside), 
             ("Planck13", "857", "muK_CMB", Nside)]
multiprocessing.Pool(processes=6).starmap(reading_model_dust, arguments)

for freq_str in ["100", "143", "217", "353", "545", "857"]:
    if os.path.isfile(root_dir+"/results_dust/dust_model_Planck13_muK_CMB_"+freq_str+"GHz.fits"):
        pass
    else:
        reading_model_dust("Planck13", freq_str, "muK_CMB", Nside)
end_time = time.time()
print('TIME: ', format((end_time-start_time)/60, '.2f'), 'min')



print("Dust map from model (Planck 2015 single MBB model, GNILC): ")
start_time = time.time()
arguments = [("Planck15-G", "100", "muK_CMB", Nside), 
             ("Planck15-G", "143", "muK_CMB", Nside), 
             ("Planck15-G", "217", "muK_CMB", Nside), 
             ("Planck15-G", "353", "muK_CMB", Nside), 
             ("Planck15-G", "545", "muK_CMB", Nside), 
             ("Planck15-G", "857", "muK_CMB", Nside)]
multiprocessing.Pool(processes=6).starmap(reading_model_dust, arguments)

for freq_str in ["100", "143", "217", "353", "545", "857"]:
    if os.path.isfile(root_dir+"/results_dust/dust_model_Planck15-G_muK_CMB_"+freq_str+"GHz.fits"):
        pass
    else:
        reading_model_dust("Planck15-G", freq_str, "muK_CMB", Nside)
end_time = time.time()
print('TIME: ', format((end_time-start_time)/60, '.2f'), 'min')



print("Dust map from model (Planck 2015 single MBB model, commander): ")
start_time = time.time()
arguments = [("Planck15-C", "100", "muK_CMB", Nside), 
             ("Planck15-C", "143", "muK_CMB", Nside), 
             ("Planck15-C", "217", "muK_CMB", Nside), 
             ("Planck15-C", "353", "muK_CMB", Nside), 
             ("Planck15-C", "545", "muK_CMB", Nside), 
             ("Planck15-C", "857", "muK_CMB", Nside)]
multiprocessing.Pool(processes=6).starmap(reading_model_dust, arguments)

for freq_str in ["100", "143", "217", "353", "545", "857"]:
    if os.path.isfile(root_dir+"/results_dust/dust_model_Planck15-C_muK_CMB_"+freq_str+"GHz.fits"):
        pass
    else:
        reading_model_dust("Planck15-C", freq_str, "muK_CMB", Nside)
end_time = time.time()
print('TIME: ', format((end_time-start_time)/60, '.2f'), 'min')



print("Dust map from model (Melis O. Irfan et al., 2019): ")
start_time = time.time()
arguments = [("Irfan19", "100", "muK_CMB", Nside), 
             ("Irfan19", "143", "muK_CMB", Nside), 
             ("Irfan19", "217", "muK_CMB", Nside), 
             ("Irfan19", "353", "muK_CMB", Nside), 
             ("Irfan19", "545", "muK_CMB", Nside), 
             ("Irfan19", "857", "muK_CMB", Nside)]
multiprocessing.Pool(processes=6).starmap(reading_model_dust, arguments)

for freq_str in ["100", "143", "217", "353", "545", "857"]:
    if os.path.isfile(root_dir+"/results_dust/dust_model_Irfan19_muK_CMB_"+freq_str+"GHz.fits"):
        pass
    else:
        reading_model_dust("Irfan19", freq_str, "muK_CMB", Nside)
end_time = time.time()
print('TIME: ', format((end_time-start_time)/60, '.2f'), 'min')



print("Dust map from model (Meisner two MBB model): ")
start_time = time.time()
arguments = [("Meisner", "100", "muK_CMB", Nside), 
             ("Meisner", "143", "muK_CMB", Nside), 
             ("Meisner", "217", "muK_CMB", Nside), 
             ("Meisner", "353", "muK_CMB", Nside), 
             ("Meisner", "545", "muK_CMB", Nside), 
             ("Meisner", "857", "muK_CMB", Nside)]
multiprocessing.Pool(processes=6).starmap(reading_model_dust, arguments)

for freq_str in ["100", "143", "217", "353", "545", "857"]:
    if os.path.isfile(root_dir+"/results_dust/dust_model_Meisner_muK_CMB_"+freq_str+"GHz.fits"):
        pass
    else:
        reading_model_dust("Meisner", freq_str, "muK_CMB", Nside)
end_time = time.time()
print('TIME: ', format((end_time-start_time)/60, '.2f'), 'min')



print("Dust map from model (DL07): ")
start_time = time.time()
arguments = [("DL07", "353", "MJy_sr", Nside), 
             ("DL07", "545", "MJy_sr", Nside), 
             ("DL07", "857", "MJy_sr", Nside)]
multiprocessing.Pool(processes=3).starmap(reading_model_dust, arguments)

for freq_str in ["353", "545", "857"]:
    if os.path.isfile(root_dir+"/results_dust/dust_model_DL07_MJy_sr_"+freq_str+"GHz.fits"):
        pass
    else:
        reading_model_dust("DL07", freq_str, "MJy_sr",  Nside)
end_time = time.time()
print('TIME: ', format((end_time-start_time)/60, '.2f'), 'min')



print("Dust map from model (GNILC): ")
start_time = time.time()
arguments = [("GNILC", "353", "MJy_sr", Nside), 
             ("GNILC", "545", "MJy_sr", Nside), 
             ("GNILC", "857", "MJy_sr", Nside)]
multiprocessing.Pool(processes=3).starmap(reading_model_dust, arguments)

for freq_str in ["353", "545", "857"]:
    if os.path.isfile(root_dir+"/results_dust/dust_model_GNILC_MJy_sr_"+freq_str+"GHz.fits"):
        pass
    else:
        reading_model_dust("GNILC", freq_str, "MJy_sr",  Nside)
end_time = time.time()
print('TIME: ', format((end_time-start_time)/60, '.2f'), 'min')



print("Dust map from model (SRoll): ")
start_time = time.time()
arguments = [("SRoll", "100", "muK_CMB", Nside), 
             ("SRoll", "143", "muK_CMB", Nside), 
             ("SRoll", "217", "muK_CMB", Nside), 
             ("SRoll", "353", "muK_CMB", Nside), 
             ("SRoll", "545", "muK_CMB", Nside), 
             ("SRoll", "857", "muK_CMB", Nside)]
multiprocessing.Pool(processes=6).starmap(reading_model_dust, arguments)

for freq_str in ["100", "143", "217", "353", "545", "857"]:
    if os.path.isfile(root_dir+"/results_dust/dust_model_SRoll_muK_CMB_"+freq_str+"GHz.fits"):
        pass
    else:
        reading_model_dust("SRoll", freq_str, "muK_CMB", Nside)
end_time = time.time()
print('TIME: ', format((end_time-start_time)/60, '.2f'), 'min')



print("Dust map from model (HD17-d5, HD17-d7): ")
start_time = time.time()
arguments = [("HD17-d5", "100", "muK_CMB", Nside), 
             ("HD17-d5", "143", "muK_CMB", Nside), 
             ("HD17-d5", "217", "muK_CMB", Nside), 
             ("HD17-d5", "353", "muK_CMB", Nside), 
             ("HD17-d5", "545", "muK_CMB", Nside), 
             ("HD17-d5", "857", "muK_CMB", Nside), 
             ("HD17-d7", "100", "muK_CMB", Nside), 
             ("HD17-d7", "143", "muK_CMB", Nside), 
             ("HD17-d7", "217", "muK_CMB", Nside), 
             ("HD17-d7", "353", "muK_CMB", Nside), 
             ("HD17-d7", "545", "muK_CMB", Nside), 
             ("HD17-d7", "857", "muK_CMB", Nside)]
multiprocessing.Pool(processes=12).starmap(reading_model_dust, arguments)

for freq_str in ["100", "143", "217", "353", "545", "857"]:
    if os.path.isfile(root_dir+"/results_dust/dust_model_HD17-d5_muK_CMB_"+freq_str+"GHz.fits"):
        pass
    else:
        reading_model_dust("HD17-d5", freq_str, "muK_CMB", Nside)
for freq_str in ["100", "143", "217", "353", "545", "857"]:
    if os.path.isfile(root_dir+"/results_dust/dust_model_HD17-d7_muK_CMB_"+freq_str+"GHz.fits"):
        pass
    else:
        reading_model_dust("HD17-d7", freq_str, "muK_CMB", Nside)
end_time = time.time()
print('TIME: ', format((end_time-start_time)/60, '.2f'), 'min')
