#--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*
# Codes for plots in the article

import healpy
import os
import numpy
import scipy
import astropy.io.fits
import time
import shutil
import sys
import matplotlib.pyplot
import multiprocessing
import public_function

root_dir = os.path.abspath('.')

start_time = time.time()

if os.path.exists(root_dir+"/figure_dust/"):
    shutil.rmtree(root_dir+"/figure_dust/")
os.mkdir(root_dir+"/figure_dust/")

if os.path.exists(root_dir+"/figure_synch/"):
    shutil.rmtree(root_dir+"/figure_synch/")
os.mkdir(root_dir+"/figure_synch/")

if os.path.exists(root_dir+"/figure_AME/"):
    shutil.rmtree(root_dir+"/figure_AME/")
os.mkdir(root_dir+"/figure_AME/")

if os.path.exists(root_dir+"/figure_free_free/"):
    shutil.rmtree(root_dir+"/figure_free_free/")
os.mkdir(root_dir+"/figure_free_free/")

if os.path.exists(root_dir+"/figure_residual/"):
    shutil.rmtree(root_dir+"/figure_residual/")
os.mkdir(root_dir+"/figure_residual/")

#--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--
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
print("Plots: ")

arguments = [("dust", "data",             "100", "muK_CMB"), 
             ("dust", "data",             "143", "muK_CMB"), 
             ("dust", "data",             "217", "muK_CMB"), 
             ("dust", "data",             "353", "muK_CMB"), 
             ("dust", "data",             "545", "muK_CMB"), 
             ("dust", "data",             "857", "muK_CMB")]
multiprocessing.Pool(processes=6).starmap(public_function.component_map, arguments)

arguments = [("dust", "model_Planck13",   "100", "muK_CMB"), 
             ("dust", "model_Planck13",   "143", "muK_CMB"), 
             ("dust", "model_Planck13",   "217", "muK_CMB"), 
             ("dust", "model_Planck13",   "353", "muK_CMB"), 
             ("dust", "model_Planck13",   "545", "muK_CMB"), 
             ("dust", "model_Planck13",   "857", "muK_CMB")]
multiprocessing.Pool(processes=6).starmap(public_function.component_map, arguments)

arguments = [("dust", "model_Planck15-G", "100", "muK_CMB"), 
             ("dust", "model_Planck15-G", "143", "muK_CMB"), 
             ("dust", "model_Planck15-G", "217", "muK_CMB"), 
             ("dust", "model_Planck15-G", "353", "muK_CMB"), 
             ("dust", "model_Planck15-G", "545", "muK_CMB"), 
             ("dust", "model_Planck15-G", "857", "muK_CMB")]
multiprocessing.Pool(processes=6).starmap(public_function.component_map, arguments)

arguments = [("dust", "model_Planck15-C", "100", "muK_CMB"), 
             ("dust", "model_Planck15-C", "143", "muK_CMB"), 
             ("dust", "model_Planck15-C", "217", "muK_CMB"), 
             ("dust", "model_Planck15-C", "353", "muK_CMB"), 
             ("dust", "model_Planck15-C", "545", "muK_CMB"), 
             ("dust", "model_Planck15-C", "857", "muK_CMB")]
multiprocessing.Pool(processes=6).starmap(public_function.component_map, arguments)

arguments = [("dust", "model_Irfan19",    "100", "muK_CMB"), 
             ("dust", "model_Irfan19",    "143", "muK_CMB"), 
             ("dust", "model_Irfan19",    "217", "muK_CMB"), 
             ("dust", "model_Irfan19",    "353", "muK_CMB"), 
             ("dust", "model_Irfan19",    "545", "muK_CMB"), 
             ("dust", "model_Irfan19",    "857", "muK_CMB")]
multiprocessing.Pool(processes=6).starmap(public_function.component_map, arguments)

arguments = [("dust", "model_Meisner",    "100", "muK_CMB"), 
             ("dust", "model_Meisner",    "143", "muK_CMB"), 
             ("dust", "model_Meisner",    "217", "muK_CMB"), 
             ("dust", "model_Meisner",    "353", "muK_CMB"), 
             ("dust", "model_Meisner",    "545", "muK_CMB"), 
             ("dust", "model_Meisner",    "857", "muK_CMB")]
multiprocessing.Pool(processes=6).starmap(public_function.component_map, arguments)

arguments = [("dust", "model_SRoll",      "100", "muK_CMB"), 
             ("dust", "model_SRoll",      "143", "muK_CMB"), 
             ("dust", "model_SRoll",      "217", "muK_CMB"), 
             ("dust", "model_SRoll",      "353", "muK_CMB"), 
             ("dust", "model_SRoll",      "545", "muK_CMB"), 
             ("dust", "model_SRoll",      "857", "muK_CMB")]
multiprocessing.Pool(processes=6).starmap(public_function.component_map, arguments)

arguments = [("dust", "model_3D",         "100", "muK_CMB"), 
             ("dust", "model_3D",         "143", "muK_CMB"), 
             ("dust", "model_3D",         "217", "muK_CMB"), 
             ("dust", "model_3D",         "353", "muK_CMB"), 
             ("dust", "model_3D",         "545", "muK_CMB"), 
             ("dust", "model_3D",         "857", "muK_CMB")]
multiprocessing.Pool(processes=6).starmap(public_function.component_map, arguments)

arguments = [("dust", "model_HD17-d5",    "100", "muK_CMB"), 
             ("dust", "model_HD17-d5",    "143", "muK_CMB"), 
             ("dust", "model_HD17-d5",    "217", "muK_CMB"), 
             ("dust", "model_HD17-d5",    "353", "muK_CMB"), 
             ("dust", "model_HD17-d5",    "545", "muK_CMB"), 
             ("dust", "model_HD17-d5",    "857", "muK_CMB")]
multiprocessing.Pool(processes=6).starmap(public_function.component_map, arguments)

arguments = [("dust", "model_HD17-d7",    "100", "muK_CMB"), 
             ("dust", "model_HD17-d7",    "143", "muK_CMB"), 
             ("dust", "model_HD17-d7",    "217", "muK_CMB"), 
             ("dust", "model_HD17-d7",    "353", "muK_CMB"), 
             ("dust", "model_HD17-d7",    "545", "muK_CMB"), 
             ("dust", "model_HD17-d7",    "857", "muK_CMB")]
multiprocessing.Pool(processes=6).starmap(public_function.component_map, arguments)

for freq_name in ["10", "11", "12", "13", "14"]:
    public_function.region(      "dust",     "data", freq_name, "muK_CMB")
    public_function.map_plot(    "dust",     "data", "model_Planck13",   freq_name, "muK_CMB")
    public_function.map_plot(    "dust",     "data", "model_Planck15-G", freq_name, "muK_CMB")
    public_function.map_plot(    "dust",     "data", "model_Planck15-C", freq_name, "muK_CMB")
    public_function.map_plot(    "dust",     "data", "model_Irfan19",    freq_name, "muK_CMB")
    public_function.map_plot(    "dust",     "data", "model_Meisner",    freq_name, "muK_CMB")
    public_function.map_plot(    "dust",     "data", "model_SRoll",      freq_name, "muK_CMB")
    public_function.map_plot(    "dust",     "data", "model_3D",         freq_name, "muK_CMB")
    public_function.map_plot(    "dust",     "data", "model_HD17-d5",    freq_name, "muK_CMB")
    public_function.map_plot(    "dust",     "data", "model_HD17-d7",    freq_name, "muK_CMB")
    public_function.scatter_plot("dust",     "data", "model_Planck13",   freq_name, "muK_CMB", threshold=3, linear_fit=True)
    public_function.scatter_plot("dust",     "data", "model_Planck15-G", freq_name, "muK_CMB", threshold=3, linear_fit=True)
    public_function.scatter_plot("dust",     "data", "model_Planck15-C", freq_name, "muK_CMB", threshold=3, linear_fit=True)
    public_function.scatter_plot("dust",     "data", "model_Irfan19",    freq_name, "muK_CMB", threshold=3, linear_fit=True)
    public_function.scatter_plot("dust",     "data", "model_Meisner",    freq_name, "muK_CMB", threshold=3, linear_fit=True)
    public_function.scatter_plot("dust",     "data", "model_SRoll",      freq_name, "muK_CMB", threshold=3, linear_fit=True)
    public_function.scatter_plot("dust",     "data", "model_3D",         freq_name, "muK_CMB", threshold=3, linear_fit=True)
    public_function.scatter_plot("dust",     "data", "model_HD17-d5",    freq_name, "muK_CMB", threshold=3, linear_fit=True)
    public_function.scatter_plot("dust",     "data", "model_HD17-d7",    freq_name, "muK_CMB", threshold=3, linear_fit=True)

for freq_name in ["13", "14"]:
#    public_function.scatter_plot("dust",     "data", "model_GNILC",      freq_name, "MJy_sr",  threshold=3, linear_fit=True)
    public_function.scatter_plot("dust",     "data", "model_DL07",       freq_name, "MJy_sr",  threshold=3, linear_fit=True)

for freq_str in ["023", "030", "033", "041", "044"]:
    public_function.component_map("synch", "data",     freq_str, "muK_CMB")
    public_function.component_map("synch", "model_s1", freq_str, "muK_CMB")
    public_function.component_map("synch", "model_s2", freq_str, "muK_CMB")
    public_function.component_map("synch", "model_s3", freq_str, "muK_CMB")
    public_function.component_map("synch", "model_s5", freq_str, "muK_CMB")
    public_function.component_map("synch", "model_s7", freq_str, "muK_CMB")
for freq_name in ["01", "02", "03", "04"]:
    public_function.region(      "synch",     "data", freq_name, "muK_CMB")
    public_function.map_plot(    "synch",     "data", "model_s1",        freq_name, "muK_CMB")
    public_function.map_plot(    "synch",     "data", "model_s2",        freq_name, "muK_CMB")
    public_function.map_plot(    "synch",     "data", "model_s3",        freq_name, "muK_CMB")
    public_function.map_plot(    "synch",     "data", "model_s5",        freq_name, "muK_CMB")
    public_function.map_plot(    "synch",     "data", "model_s7",        freq_name, "muK_CMB")
    public_function.scatter_plot("synch",     "data", "model_s1",        freq_name, "muK_CMB", threshold=0, linear_fit=True)
    public_function.scatter_plot("synch",     "data", "model_s2",        freq_name, "muK_CMB", threshold=0, linear_fit=True)
    public_function.scatter_plot("synch",     "data", "model_s3",        freq_name, "muK_CMB", threshold=0, linear_fit=True)
    public_function.scatter_plot("synch",     "data", "model_s5",        freq_name, "muK_CMB", threshold=0, linear_fit=True)
    public_function.scatter_plot("synch",     "data", "model_s7",        freq_name, "muK_CMB", threshold=0, linear_fit=True)

for freq_str in ["030", "044", "070"]:
    public_function.component_map("AME",      "data",  freq_str, "muK_CMB")
    public_function.component_map("AME",      "model", freq_str, "muK_CMB")
public_function.region(           "AME",      "data",                    "06",      "muK_CMB")
public_function.region(           "AME",      "data",                    "08",      "muK_CMB")
public_function.scatter_plot(     "AME",      "data", "model",           "06",      "muK_CMB", threshold=0, linear_fit=True)
public_function.scatter_plot(     "AME",      "data", "model",           "08",      "muK_CMB", threshold=0, linear_fit=False)

for freq_str in ["030", "044", "070"]:
    public_function.component_map("free_free", "data",  freq_str, "muK_CMB")
    public_function.component_map("free_free", "model", freq_str, "muK_CMB")
for freq_name in ["06", "08"]:
    public_function.region("free_free",  "data",          freq_name, "muK_CMB")
#    public_function.map_plot("free_free",  "data", "model", freq_name, "muK_CMB")
    public_function.scatter_plot("free_free", "data", "model",           freq_name, "muK_CMB", threshold=0, linear_fit=False)

for freq_str in ["030", "044", "070", "100", "143", "217", "353", "545", "857"]:
    public_function.component_map("residual",  "data",  freq_str, "muK_CMB")

def ratio_residual(freq_str):
    residual = healpy.read_map(root_dir+"/results_residual/residual_data_muK_CMB_"+freq_str+"GHz.fits")
    dust     = public_function.reading_dust("Planck15-G", freq_str, "muK_CMB", 512)
    ratio    = residual / dust
    mask_map = public_function.Statistic_mask(512, 2, 1.5, 1)
    ratio[mask_map == 0] = healpy.UNSEEN
    MIN      = numpy.percentile(ratio[mask_map != 0], 10)
    MAX      = numpy.percentile(ratio[mask_map != 0], 90)
    matplotlib.pyplot.figure(dpi=300, figsize=(4, 3))
    healpy.mollview(ratio, cmap="Spectral_r", format="%.3g", nlocs=3, min=MIN, max=MAX, hold=True)
    matplotlib.pyplot.title("residual / dust ("+freq_str+" GHz)", fontsize=11)
    matplotlib.pyplot.savefig(root_dir+"/figure_residual/ratio_"+freq_str+".png", bbox_inches="tight")
    matplotlib.pyplot.close()
    
for freq_str in ["353", "545", "857"]:
    ratio_residual(freq_str)

end_time = time.time()
print("Plots, Succeed ! \n TIME: ", format( (end_time-start_time)/60, '.2f'),"min")