# 0. Introduction
This project is used to test the goodness of Galactic foreground (thermal dust emission, synchrotron, free-free emission, and AME) models in CMB experiments. 

More technical details can be found in the following paper: 

# 1. Dependencies and installation
## 1.1 Dependencies
The main dependencies and their versions are listed below: 

`pip         version: 25.1.1`
`Python      version: 3.9.18`
`PySM3       version: 3.4.2`
`healpy      version: 1.17.3`
`numpy       version: 1.26.4`
`scipy       version: 1.13.1`
`astropy     version: 6.0.1`
`matplotlib  version: 3.9.4`
`pandas      version: 2.3.1`

The computer configuration should be no less than 12 cores and 64GB of memory. 

## Installation

```
git clone https://github.com/XXX
cd XXX
bash configure.sh
bash file.sh
```

The default behavior of `configure.sh` is a full installation of python packages. 

### Installation of python packages


The default behavior of `file.sh` is a full installation of promordial .fits file. 
 
### Download .fits files 

# 2. Running
```
bash run.sh
```
## Step-1
```
python3 -u Inpaint_Rebeam_frequency_maps.py 2000
```

1. Masking HFI maps with point source masks. 
2. Iterating 2000 times to in-paint the holes after masking. 

## Step-2
```
python3 -u Dust_data_model.py
python3 -u Synchrotron_data_model.py
python3 -u free_free_data_model.py
python3 -u AME_data_model.py
python3 -u Residual.py

```
1. Isolating single foreground component data maps from Planck maps and WMAP maps. 
2. Calculating single foreground component model maps. 

## Step-3
```
python Statistics.py
```
Calculating statistical quantities, including linear ratio R (for model), R' (for data), and correlation coefficient C' (for data). 

## Step-5
```
python Plots.py
```
Plotting. 

# 3. Citation

2015 ApJ 798 88
http://planck.skymaps.info
https://doi.org/10.1093/mnras/stx949
http://doi.org/10.21105/joss.03783
https://pysm3.readthedocs.io/en/latest/
https://doi.org/10.1051/0004-6361/201525659