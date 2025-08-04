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
import public_function

print("Statistics.py, Starting. ")
root_dir = os.path.abspath('.')

start_time = time.time()

#--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*
# Mosaic(component, name, freq_str_1, freq_str_2, unit, nside1, nside2, smooth_degree, disk_degree, correlation_threshold)
# with color correction
# smooth_degree = 2 degree, disk_degree = 6 degree

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

arguments = [
("dust", "data",             "100", "143", "muK_CMB", 64, 512, 2, 6, 0.95),
("dust", "data",             "143", "217", "muK_CMB", 64, 512, 2, 6, 0.95),
("dust", "data",             "217", "353", "muK_CMB", 64, 512, 2, 6, 0.95),
("dust", "data",             "353", "545", "muK_CMB", 64, 512, 2, 6, 0.95),
("dust", "data",             "545", "857", "muK_CMB", 64, 512, 2, 6, 0.95),
("dust", "data",             "353", "545", "MJy_sr",  64, 512, 2, 6, 0.95),
("dust", "data",             "545", "857", "MJy_sr",  64, 512, 2, 6, 0.95)]
multiprocessing.Pool(processes=7).starmap(public_function.Mosaic, arguments)
arguments = [
("dust", "model_Planck13",   "100", "143", "muK_CMB", 64, 512, 2, 6, 0.95),
("dust", "model_Planck13",   "143", "217", "muK_CMB", 64, 512, 2, 6, 0.95),
("dust", "model_Planck13",   "217", "353", "muK_CMB", 64, 512, 2, 6, 0.95),
("dust", "model_Planck13",   "353", "545", "muK_CMB", 64, 512, 2, 6, 0.95),
("dust", "model_Planck13",   "545", "857", "muK_CMB", 64, 512, 2, 6, 0.95)]
multiprocessing.Pool(processes=5).starmap(public_function.Mosaic, arguments)
arguments = [
("dust", "model_Planck15-G", "100", "143", "muK_CMB", 64, 512, 2, 6, 0.95),
("dust", "model_Planck15-G", "143", "217", "muK_CMB", 64, 512, 2, 6, 0.95),
("dust", "model_Planck15-G", "217", "353", "muK_CMB", 64, 512, 2, 6, 0.95),
("dust", "model_Planck15-G", "353", "545", "muK_CMB", 64, 512, 2, 6, 0.95),
("dust", "model_Planck15-G", "545", "857", "muK_CMB", 64, 512, 2, 6, 0.95)]
multiprocessing.Pool(processes=5).starmap(public_function.Mosaic, arguments)
arguments = [
("dust", "model_Planck15-C", "100", "143", "muK_CMB", 64, 512, 2, 6, 0.95),
("dust", "model_Planck15-C", "143", "217", "muK_CMB", 64, 512, 2, 6, 0.95),
("dust", "model_Planck15-C", "217", "353", "muK_CMB", 64, 512, 2, 6, 0.95),
("dust", "model_Planck15-C", "353", "545", "muK_CMB", 64, 512, 2, 6, 0.95),
("dust", "model_Planck15-C", "545", "857", "muK_CMB", 64, 512, 2, 6, 0.95)]
multiprocessing.Pool(processes=5).starmap(public_function.Mosaic, arguments)
arguments = [
("dust", "model_Irfan19",    "100", "143", "muK_CMB", 64, 512, 2, 6, 0.95),
("dust", "model_Irfan19",    "143", "217", "muK_CMB", 64, 512, 2, 6, 0.95),
("dust", "model_Irfan19",    "217", "353", "muK_CMB", 64, 512, 2, 6, 0.95),
("dust", "model_Irfan19",    "353", "545", "muK_CMB", 64, 512, 2, 6, 0.95),
("dust", "model_Irfan19",    "545", "857", "muK_CMB", 64, 512, 2, 6, 0.95)]
multiprocessing.Pool(processes=5).starmap(public_function.Mosaic, arguments)
arguments = [
("dust", "model_Meisner",    "100", "143", "muK_CMB", 64, 512, 2, 6, 0.95),
("dust", "model_Meisner",    "143", "217", "muK_CMB", 64, 512, 2, 6, 0.95),
("dust", "model_Meisner",    "217", "353", "muK_CMB", 64, 512, 2, 6, 0.95),
("dust", "model_Meisner",    "353", "545", "muK_CMB", 64, 512, 2, 6, 0.95),
("dust", "model_Meisner",    "545", "857", "muK_CMB", 64, 512, 2, 6, 0.95)]
multiprocessing.Pool(processes=5).starmap(public_function.Mosaic, arguments)
arguments = [
("dust", "model_SRoll",      "100", "143", "muK_CMB", 64, 512, 2, 6, 0.95),
("dust", "model_SRoll",      "143", "217", "muK_CMB", 64, 512, 2, 6, 0.95),
("dust", "model_SRoll",      "217", "353", "muK_CMB", 64, 512, 2, 6, 0.95),
("dust", "model_SRoll",      "353", "545", "muK_CMB", 64, 512, 2, 6, 0.95),
("dust", "model_SRoll",      "545", "857", "muK_CMB", 64, 512, 2, 6, 0.95)]
multiprocessing.Pool(processes=5).starmap(public_function.Mosaic, arguments)
arguments = [
("dust", "model_GNILC",      "545", "857", "MJy_sr",  64, 512, 2, 6, 0.95),
("dust", "model_GNILC",      "353", "545", "MJy_sr",  64, 512, 2, 6, 0.95),
("dust", "model_GNILC",      "353", "857", "MJy_sr",  64, 512, 2, 6, 0.95), 
("dust", "model_DL07",       "545", "857", "MJy_sr",  64, 512, 2, 6, 0.95),
("dust", "model_DL07",       "353", "545", "MJy_sr",  64, 512, 2, 6, 0.95),
("dust", "model_DL07",       "353", "857", "MJy_sr",  64, 512, 2, 6, 0.95),
("dust", "model_3D",         "100", "143", "muK_CMB", 64, 512, 2, 6, 0.95),
("dust", "model_3D",         "143", "217", "muK_CMB", 64, 512, 2, 6, 0.95),
("dust", "model_3D",         "217", "353", "muK_CMB", 64, 512, 2, 6, 0.95),
("dust", "model_3D",         "353", "545", "muK_CMB", 64, 512, 2, 6, 0.95),
("dust", "model_3D",         "545", "857", "muK_CMB", 64, 512, 2, 6, 0.95),
("dust", "model_HD17-d5",    "100", "143", "muK_CMB", 64, 512, 2, 6, 0.95),
("dust", "model_HD17-d5",    "143", "217", "muK_CMB", 64, 512, 2, 6, 0.95),
("dust", "model_HD17-d5",    "217", "353", "muK_CMB", 64, 512, 2, 6, 0.95),
("dust", "model_HD17-d5",    "353", "545", "muK_CMB", 64, 512, 2, 6, 0.95),
("dust", "model_HD17-d5",    "545", "857", "muK_CMB", 64, 512, 2, 6, 0.95),
("dust", "model_HD17-d7",    "100", "143", "muK_CMB", 64, 512, 2, 6, 0.95),
("dust", "model_HD17-d7",    "143", "217", "muK_CMB", 64, 512, 2, 6, 0.95),
("dust", "model_HD17-d7",    "217", "353", "muK_CMB", 64, 512, 2, 6, 0.95),
("dust", "model_HD17-d7",    "353", "545", "muK_CMB", 64, 512, 2, 6, 0.95),
("dust", "model_HD17-d7",    "545", "857", "muK_CMB", 64, 512, 2, 6, 0.95)]
multiprocessing.Pool(processes=21).starmap(public_function.Mosaic, arguments)

arguments = [
("synch", "data",            "023", "030", "muK_CMB", 64, 512, 2, 6, 0.90),
("synch", "data",            "030", "033", "muK_CMB", 64, 512, 2, 6, 0.90),
("synch", "data",            "033", "041", "muK_CMB", 64, 512, 2, 6, 0.90),
("synch", "data",            "041", "044", "muK_CMB", 64, 512, 2, 6, 0.90)]
multiprocessing.Pool(processes=4).starmap(public_function.Mosaic, arguments)
arguments = [
("synch", "model_s1",        "023", "030", "muK_CMB", 64, 512, 2, 6, 0.90),
("synch", "model_s1",        "030", "033", "muK_CMB", 64, 512, 2, 6, 0.90),
("synch", "model_s1",        "033", "041", "muK_CMB", 64, 512, 2, 6, 0.90),
("synch", "model_s1",        "041", "044", "muK_CMB", 64, 512, 2, 6, 0.90),
("synch", "model_s2",        "023", "030", "muK_CMB", 64, 512, 2, 6, 0.90),
("synch", "model_s2",        "030", "033", "muK_CMB", 64, 512, 2, 6, 0.90),
("synch", "model_s2",        "033", "041", "muK_CMB", 64, 512, 2, 6, 0.90),
("synch", "model_s2",        "041", "044", "muK_CMB", 64, 512, 2, 6, 0.90),
("synch", "model_s3",        "023", "030", "muK_CMB", 64, 512, 2, 6, 0.90),
("synch", "model_s3",        "030", "033", "muK_CMB", 64, 512, 2, 6, 0.90),
("synch", "model_s3",        "033", "041", "muK_CMB", 64, 512, 2, 6, 0.90),
("synch", "model_s3",        "041", "044", "muK_CMB", 64, 512, 2, 6, 0.90),
("synch", "model_s5",        "023", "030", "muK_CMB", 64, 512, 2, 6, 0.90),
("synch", "model_s5",        "030", "033", "muK_CMB", 64, 512, 2, 6, 0.90),
("synch", "model_s5",        "033", "041", "muK_CMB", 64, 512, 2, 6, 0.90),
("synch", "model_s5",        "041", "044", "muK_CMB", 64, 512, 2, 6, 0.90),
("synch", "model_s7",        "023", "030", "muK_CMB", 64, 512, 2, 6, 0.90),
("synch", "model_s7",        "030", "033", "muK_CMB", 64, 512, 2, 6, 0.90),
("synch", "model_s7",        "033", "041", "muK_CMB", 64, 512, 2, 6, 0.90),
("synch", "model_s7",        "041", "044", "muK_CMB", 64, 512, 2, 6, 0.90)]
multiprocessing.Pool(processes=28).starmap(public_function.Mosaic, arguments)

arguments = [
("AME", "data",              "030", "044", "muK_CMB", 64, 512, 2, 6, 0.90),
("AME", "data",              "044", "070", "muK_CMB", 64, 512, 2, 6, 0.90),
("AME", "data",              "070", "100", "muK_CMB", 64, 512, 2, 6, 0.90)]
multiprocessing.Pool(processes=3).starmap(public_function.Mosaic, arguments)
arguments = [
("AME", "model",             "030", "044", "muK_CMB", 64, 512, 2, 6, 0.90),
("AME", "model",             "044", "070", "muK_CMB", 64, 512, 2, 6, 0.90),
("AME", "model",             "070", "100", "muK_CMB", 64, 512, 2, 6, 0.90)]
multiprocessing.Pool(processes=3).starmap(public_function.Mosaic, arguments)

arguments = [
("free_free", "data",        "030", "044", "muK_CMB", 64, 512, 2, 6, 0.90),
("free_free", "data",        "044", "070", "muK_CMB", 64, 512, 2, 6, 0.90)]
multiprocessing.Pool(processes=2).starmap(public_function.Mosaic, arguments)
arguments = [
("free_free", "model",       "030", "044", "muK_CMB", 64, 512, 2, 6, 0.90),
("free_free", "model",       "044", "070", "muK_CMB", 64, 512, 2, 6, 0.90)]
multiprocessing.Pool(processes=2).starmap(public_function.Mosaic, arguments)

end_time = time.time()
print("Statistics, Succeed ! \n TIME: ", format( (end_time-start_time)/60, '.2f'),"min")