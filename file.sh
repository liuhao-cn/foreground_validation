#!/bin/bash

# Reference
# 2015 ApJ 798 88
# http://planck.skymaps.info
# https://doi.org/10.1093/mnras/stx949
# http://doi.org/10.21105/joss.03783
# https://pysm3.readthedocs.io/en/latest/
# https://doi.org/10.1051/0004-6361/201525659
# http://sroll20.ias.u-psud.fr

# ------------------------------------
# ------------------------------------
# WMAP

if [[ -d 'WMAP' ]]
then
    echo 'You have directory WMAP'
else
    mkdir WMAP
fi

if [[ -f WMAP/wmap_band_imap_r9_9yr_K_v5.fits ]]; then
    echo 'You have file wmap_band_imap_r9_9yr_K_v5.fits'
else
    while true; do
        wget -cP WMAP/ 'https://lambda.gsfc.nasa.gov/data/map/dr5/skymaps/9yr/raw/wmap_band_imap_r9_9yr_K_v5.fits'
        if [[ $(md5sum WMAP/wmap_band_imap_r9_9yr_K_v5.fits | awk '{print $1}') == '50229a5cc6dd3933be56bf334c892b5a' ]]; then
            break
        fi
    done
fi

if [[ -f WMAP/wmap_band_imap_r9_9yr_Ka_v5.fits ]]
then
    echo 'You have file wmap_band_imap_r9_9yr_Ka_v5.fits'
else
    while true; do
        wget -cP WMAP/ 'https://lambda.gsfc.nasa.gov/data/map/dr5/skymaps/9yr/raw/wmap_band_imap_r9_9yr_Ka_v5.fits'
        if [[ $(md5sum WMAP/wmap_band_imap_r9_9yr_Ka_v5.fits | awk '{print $1}') == '4a0afb08d4ad25e765d8e5128b23cfdf' ]]; then
            break
        fi
    done
fi

if [[ -f WMAP/wmap_band_imap_r9_9yr_Q_v5.fits ]]
then
    echo 'You have file wmap_band_imap_r9_9yr_Q_v5.fits'
else
    while true; do
        wget -cP WMAP/ 'https://lambda.gsfc.nasa.gov/data/map/dr5/skymaps/9yr/raw/wmap_band_imap_r9_9yr_Q_v5.fits'
        if [[ $(md5sum WMAP/wmap_band_imap_r9_9yr_Q_v5.fits | awk '{print $1}') == 'e7fbf02b2cfa126abac9a5e7df229b0e' ]]; then
            break
        fi
    done
fi

if [[ -f WMAP/wmap_band_imap_r9_9yr_V_v5.fits ]]
then
    echo 'You have file wmap_band_imap_r9_9yr_V_v5.fits'
else
    while true; do
        wget -cP WMAP/ 'https://lambda.gsfc.nasa.gov/data/map/dr5/skymaps/9yr/raw/wmap_band_imap_r9_9yr_V_v5.fits'
        if [[ $(md5sum WMAP/wmap_band_imap_r9_9yr_V_v5.fits | awk '{print $1}') == 'adf22cccfd04d0f75852427d5840ed45' ]]; then
            break
        fi
    done
fi

if [[ -f WMAP/wmap_processing_mask_K_yr2_r9_9yr_v5.fits ]]
then
    echo 'You have file wmap_processing_mask_K_yr2_r9_9yr_v5.fits'
else
    while true; do
        wget -cP WMAP/ 'https://lambda.gsfc.nasa.gov/data/map/dr5/ancillary/masks/wmap_processing_mask_K_yr2_r9_9yr_v5.fits'
        if [[ $(md5sum WMAP/wmap_processing_mask_K_yr2_r9_9yr_v5.fits | awk '{print $1}') == '9d1689eb051c72538876d49bbaec3b9b' ]]; then
            break
        fi
    done
fi

if [[ -f WMAP/wmap_processing_mask_Ka_r9_9yr_v5.fits ]]
then
    echo 'You have file wmap_processing_mask_Ka_r9_9yr_v5.fits'
else
    while true; do
        wget -cP WMAP/ 'https://lambda.gsfc.nasa.gov/data/map/dr5/ancillary/masks/wmap_processing_mask_Ka_r9_9yr_v5.fits'
        if [[ $(md5sum WMAP/wmap_processing_mask_Ka_r9_9yr_v5.fits | awk '{print $1}') == '32c31034d623ce1c1cf31a85cab4f035' ]]; then
            break
        fi
    done
fi

if [[ -f WMAP/wmap_processing_mask_Q_r9_9yr_v5.fits ]]
then
    echo 'You have file wmap_processing_mask_Q_r9_9yr_v5.fits'
else
    while true; do
        wget -cP WMAP/ 'https://lambda.gsfc.nasa.gov/data/map/dr5/ancillary/masks/wmap_processing_mask_Q_r9_9yr_v5.fits'
        if [[ $(md5sum WMAP/wmap_processing_mask_Q_r9_9yr_v5.fits | awk '{print $1}') == '54d54806239af1969849fe83c085a9b9' ]]; then
            break
        fi
    done
fi

if [[ -f WMAP/wmap_processing_mask_V_r9_9yr_v5.fits ]]
then
    echo 'You have file wmap_processing_mask_V_r9_9yr_v5.fits'
else
    while true; do
        wget -cP WMAP/ 'https://lambda.gsfc.nasa.gov/data/map/dr5/ancillary/masks/wmap_processing_mask_V_r9_9yr_v5.fits'
        if [[ $(md5sum WMAP/wmap_processing_mask_V_r9_9yr_v5.fits | awk '{print $1}') == 'a0a43b7c98782b4fdae31e0957c00779' ]]; then
            break
        fi
    done
fi

# ------------------------------------
# Planck 2013 release

if [[ -d 'PR1-2013' ]]
then
    echo 'You have directory PR1-2013'
else
    mkdir PR1-2013
fi

# ------------------------------------
# Thermal dust model (Planck 2013 Commander)
if [[ -f PR1-2013/HFI_CompMap_ThermalDustModel_2048_R1.20.fits ]]
then
    echo 'You have file HFI_CompMap_ThermalDustModel_2048_R1.20.fits'
else
    while true; do
        wget -cP PR1-2013/ 'https://irsa.ipac.caltech.edu/data/Planck/release_1/all-sky-maps/maps/HFI_CompMap_ThermalDustModel_2048_R1.20.fits'
        if [[ $(md5sum PR1-2013/HFI_CompMap_ThermalDustModel_2048_R1.20.fits | awk '{print $1}') == '8d804f4e64e709f476a63f0dfed1fd11' ]]; then
            break
        fi
    done
fi

# ------------------------------------
# ------------------------------------
# Planck 2015 release

if [[ -d 'PR2-2015' ]]
then
    echo 'You have directory PR2-2015'
else
    mkdir PR2-2015
fi

# ------------------------------------
# Foreground: synchrotron, free-free, CO21
if [[ -f PR2-2015/COM_CompMap_Synchrotron-commander_0256_R2.00.fits ]]
then
    echo 'You have file COM_CompMap_Synchrotron-commander_0256_R2.00.fits'
else
    while true; do
        wget -cP PR2-2015/ 'https://irsa.ipac.caltech.edu/data/Planck/release_2/all-sky-maps/maps/component-maps/foregrounds/COM_CompMap_Synchrotron-commander_0256_R2.00.fits'
        if [[ $(md5sum PR2-2015/COM_CompMap_Synchrotron-commander_0256_R2.00.fits | awk '{print $1}') == '9f67a32409ce073c8aca7ec63f1973da' ]]; then
            break
        fi
    done
fi

if [[ -f PR2-2015/COM_CompMap_freefree-commander_0256_R2.00.fits ]]
then
    echo 'You have file COM_CompMap_freefree-commander_0256_R2.00.fits'
else
    while true; do
        wget -cP PR2-2015/ 'https://irsa.ipac.caltech.edu/data/Planck/release_2/all-sky-maps/maps/component-maps/foregrounds/COM_CompMap_freefree-commander_0256_R2.00.fits'
        if [[ $(md5sum PR2-2015/COM_CompMap_freefree-commander_0256_R2.00.fits | awk '{print $1}') == '78c08bee765f0951ddcd1fe4046d8cf8' ]]; then
            break
        fi
    done
fi

if [[ -f PR2-2015/COM_CompMap_CO21-commander_2048_R2.00.fits ]]
then
    echo 'You have file COM_CompMap_CO21-commander_2048_R2.00.fits'
else
    while true; do
        wget -cP PR2-2015/ 'https://irsa.ipac.caltech.edu/data/Planck/release_2/all-sky-maps/maps/component-maps/foregrounds/COM_CompMap_CO21-commander_2048_R2.00.fits'
        if [[ $(md5sum PR2-2015/COM_CompMap_CO21-commander_2048_R2.00.fits | awk '{print $1}') == 'd533d6949766ef3087e8046e047025a1' ]]; then
            break
        fi
    done
fi

if [[ -f PR2-2015/COM_CompMap_CO-commander_0256_R2.00.fits ]]
then
    echo 'You have file COM_CompMap_CO-commander_0256_R2.00.fits'
else
    while true; do
        wget -cP PR2-2015/ 'https://irsa.ipac.caltech.edu/data/Planck/release_2/all-sky-maps/maps/component-maps/foregrounds/COM_CompMap_CO-commander_0256_R2.00.fits'
        if [[ $(md5sum PR2-2015/COM_CompMap_CO-commander_0256_R2.00.fits | awk '{print $1}') == '3a481db8de0b216f56790f59b7c5f591' ]]; then
            break
        fi
    done
fi

if [[ -f PR2-2015/COM_CompMap_xline-commander_0256_R2.00.fits ]]
then
    echo 'You have file COM_CompMap_xline-commander_0256_R2.00.fits'
else
    while true; do
        wget -cP PR2-2015/ 'https://irsa.ipac.caltech.edu/data/Planck/release_2/all-sky-maps/maps/component-maps/foregrounds/COM_CompMap_xline-commander_0256_R2.00.fits'
        if [[ $(md5sum PR2-2015/COM_CompMap_xline-commander_0256_R2.00.fits | awk '{print $1}') == '0f25975c198dd0fc0d251a887578bc63' ]]; then
            break
        fi
    done
fi

if [[ -f PR2-2015/COM_CompMap_AME-commander_0256_R2.00.fits ]]
then
    echo 'You have file COM_CompMap_AME-commander_0256_R2.00.fits'
else
    while true; do
        wget -cP PR2-2015/ 'https://irsa.ipac.caltech.edu/data/Planck/release_2/all-sky-maps/maps/component-maps/foregrounds/COM_CompMap_AME-commander_0256_R2.00.fits'
        if [[ $(md5sum PR2-2015/COM_CompMap_AME-commander_0256_R2.00.fits | awk '{print $1}') == 'd02e0d8cf4d0cedd8f3e628a2545b87a' ]]; then
            break
        fi
    done
fi

# ------------------------------------
# Masks
if [[ -f PR2-2015/HFI_Mask_GalPlane-apo0_2048_R2.00.fits ]]
then
    echo 'You have file HFI_Mask_GalPlane-apo0_2048_R2.00.fits'
else
    while true; do
        wget -cP PR2-2015/ 'https://irsa.ipac.caltech.edu/data/Planck/release_2/ancillary-data/masks/HFI_Mask_GalPlane-apo0_2048_R2.00.fits'
        if [[ $(md5sum PR2-2015/HFI_Mask_GalPlane-apo0_2048_R2.00.fits | awk '{print $1}') == 'b4a81e555055fe8cca607f078b6dab99' ]]; then
            break
        fi
    done
fi

if [[ -f PR2-2015/HFI_Mask_PointSrc_2048_R2.00.fits ]]
then
    echo 'You have file HFI_Mask_PointSrc_2048_R2.00.fits'
else
    while true; do
        wget -cP PR2-2015/ 'https://irsa.ipac.caltech.edu/data/Planck/release_2/ancillary-data/masks/HFI_Mask_PointSrc_2048_R2.00.fits'
        if [[ $(md5sum PR2-2015/HFI_Mask_PointSrc_2048_R2.00.fits | awk '{print $1}') == 'd26a34315fb8ca1df9014cddb2c7f40f' ]]; then
            break
        fi
    done
fi

if [[ -f PR2-2015/LFI_Mask_PointSrc_2048_R2.00.fits ]]
then
    echo 'You have file LFI_Mask_PointSrc_2048_R2.00.fits'
else
    while true; do
        wget -cP PR2-2015/ 'https://irsa.ipac.caltech.edu/data/Planck/release_2/ancillary-data/masks/LFI_Mask_PointSrc_2048_R2.00.fits'
        if [[ $(md5sum PR2-2015/LFI_Mask_PointSrc_2048_R2.00.fits | awk '{print $1}') == 'd07eb114dc73ee43a799226177f5e60d' ]]; then
            break
        fi
    done
fi

# ------------------------------------
# Compact source catalog
if [[ -f PR2-2015/COM_PCCS_030_R2.04.fits ]]
then
    echo 'You have file COM_PCCS_030_R2.04.fits'
else
    while true; do
        wget -cP PR2-2015/ 'https://irsa.ipac.caltech.edu/data/Planck/release_2/catalogs/COM_PCCS_030_R2.04.fits'
        if [[ $(md5sum PR2-2015/COM_PCCS_030_R2.04.fits | awk '{print $1}') == '8af8803720dfdaa591591fce35a33046' ]]; then
            break
        fi
    done
fi

if [[ -f PR2-2015/COM_PCCS_044_R2.04.fits ]]
then
    echo 'You have file COM_PCCS_044_R2.04.fits'
else
    while true; do
        wget -cP PR2-2015/ 'https://irsa.ipac.caltech.edu/data/Planck/release_2/catalogs/COM_PCCS_044_R2.04.fits'
        if [[ $(md5sum PR2-2015/COM_PCCS_044_R2.04.fits | awk '{print $1}') == '7967b3080d547ba2dcca5db401f924fa' ]]; then
            break
        fi
    done
fi

if [[ -f PR2-2015/COM_PCCS_070_R2.04.fits ]]
then
    echo 'You have file COM_PCCS_070_R2.04.fits'
else
    while true; do
        wget -cP PR2-2015/ 'https://irsa.ipac.caltech.edu/data/Planck/release_2/catalogs/COM_PCCS_070_R2.04.fits'
        if [[ $(md5sum PR2-2015/COM_PCCS_070_R2.04.fits | awk '{print $1}') == 'bb7327d83abf8012d67ee6080cf48274' ]]; then
            break
        fi
    done
fi

if [[ -f PR2-2015/COM_PCCS_100-excluded_R2.01.fits ]]
then
    echo 'You have file COM_PCCS_100-excluded_R2.01.fits'
else
    while true; do
        wget -cP PR2-2015/ 'https://irsa.ipac.caltech.edu/data/Planck/release_2/catalogs/COM_PCCS_100-excluded_R2.01.fits'
        if [[ $(md5sum PR2-2015/COM_PCCS_100-excluded_R2.01.fits | awk '{print $1}') == '06138e9b04d719f556d2b14057d5b8e2' ]]; then
            break
        fi
    done
fi

if [[ -f PR2-2015/COM_PCCS_100_R2.01.fits ]]
then
    echo 'You have file COM_PCCS_100_R2.01.fits'
else
    while true; do
        wget -cP PR2-2015/ 'https://irsa.ipac.caltech.edu/data/Planck/release_2/catalogs/COM_PCCS_100_R2.01.fits'
        if [[ $(md5sum PR2-2015/COM_PCCS_100_R2.01.fits | awk '{print $1}') == '55463c1055af5cf34e348ac88f9c39da' ]]; then
            break
        fi
    done
fi

if [[ -f PR2-2015/COM_PCCS_143-excluded_R2.01.fits ]]
then
    echo 'You have file COM_PCCS_143-excluded_R2.01.fits'
else
    while true; do
        wget -cP PR2-2015/ 'https://irsa.ipac.caltech.edu/data/Planck/release_2/catalogs/COM_PCCS_143-excluded_R2.01.fits'
        if [[ $(md5sum PR2-2015/COM_PCCS_143-excluded_R2.01.fits | awk '{print $1}') == '0e48ead13f813e00075a3dc84f5458c1' ]]; then
            break
        fi
    done
fi

if [[ -f PR2-2015/COM_PCCS_143_R2.01.fits ]]
then
    echo 'You have file COM_PCCS_143_R2.01.fits'
else
    while true; do
        wget -cP PR2-2015/ 'https://irsa.ipac.caltech.edu/data/Planck/release_2/catalogs/COM_PCCS_143_R2.01.fits'
        if [[ $(md5sum PR2-2015/COM_PCCS_143_R2.01.fits | awk '{print $1}') == '43b48f6a45f5ba1b19b8cc6827c09231' ]]; then
            break
        fi
    done
fi

if [[ -f PR2-2015/COM_PCCS_217-excluded_R2.01.fits ]]
then
    echo 'You have file COM_PCCS_217-excluded_R2.01.fits'
else
    while true; do
        wget -cP PR2-2015/ 'https://irsa.ipac.caltech.edu/data/Planck/release_2/catalogs/COM_PCCS_217-excluded_R2.01.fits'
        if [[ $(md5sum PR2-2015/COM_PCCS_217-excluded_R2.01.fits | awk '{print $1}') == '48407343ee676c3282dfd651fd1ae61b' ]]; then
            break
        fi
    done
fi

if [[ -f PR2-2015/COM_PCCS_217_R2.01.fits ]]
then
    echo 'You have file COM_PCCS_217_R2.01.fits'
else
    while true; do
        wget -cP PR2-2015/ 'https://irsa.ipac.caltech.edu/data/Planck/release_2/catalogs/COM_PCCS_217_R2.01.fits'
        if [[ $(md5sum PR2-2015/COM_PCCS_217_R2.01.fits | awk '{print $1}') == '24db4f98910bc3dbbc19b6e6c76cd3f2' ]]; then
            break
        fi
    done
fi

if [[ -f PR2-2015/COM_PCCS_353-excluded_R2.01.fits ]]
then
    echo 'You have file COM_PCCS_353-excluded_R2.01.fits'
else
    while true; do
        wget -cP PR2-2015/ 'https://irsa.ipac.caltech.edu/data/Planck/release_2/catalogs/COM_PCCS_353-excluded_R2.01.fits'
        if [[ $(md5sum PR2-2015/COM_PCCS_353-excluded_R2.01.fits | awk '{print $1}') == 'd4a894d3800c2bccd429eed6d0706f92' ]]; then
            break
        fi
    done
fi

if [[ -f PR2-2015/COM_PCCS_353_R2.01.fits ]]
then
    echo 'You have file COM_PCCS_353_R2.01.fits'
else
    while true; do
        wget -cP PR2-2015/ 'https://irsa.ipac.caltech.edu/data/Planck/release_2/catalogs/COM_PCCS_353_R2.01.fits'
        if [[ $(md5sum PR2-2015/COM_PCCS_353_R2.01.fits | awk '{print $1}') == '18fd480c02fe3aefb3b30cefabc034ba' ]]; then
            break
        fi
    done
fi

if [[ -f PR2-2015/COM_PCCS_545-excluded_R2.01.fits ]]
then
    echo 'You have file COM_PCCS_545-excluded_R2.01.fits'
else
    while true; do
        wget -cP PR2-2015/ 'https://irsa.ipac.caltech.edu/data/Planck/release_2/catalogs/COM_PCCS_545-excluded_R2.01.fits'
        if [[ $(md5sum PR2-2015/COM_PCCS_545-excluded_R2.01.fits | awk '{print $1}') == 'b706e352f522fdfaeaa71fad0b455d11' ]]; then
            break
        fi
    done
fi

if [[ -f PR2-2015/COM_PCCS_545_R2.01.fits ]]
then
    echo 'You have file COM_PCCS_545_R2.01.fits'
else
    while true; do
        wget -cP PR2-2015/ 'https://irsa.ipac.caltech.edu/data/Planck/release_2/catalogs/COM_PCCS_545_R2.01.fits'
        if [[ $(md5sum PR2-2015/COM_PCCS_545_R2.01.fits | awk '{print $1}') == '9c9fc809f9886f00307a6eebbe4412e7' ]]; then
            break
        fi
    done
fi

if [[ -f PR2-2015/COM_PCCS_857-excluded_R2.01.fits ]]
then
    echo 'You have file COM_PCCS_857-excluded_R2.01.fits'
else
    while true; do
        wget -cP PR2-2015/ 'https://irsa.ipac.caltech.edu/data/Planck/release_2/catalogs/COM_PCCS_857-excluded_R2.01.fits'
        if [[ $(md5sum PR2-2015/COM_PCCS_857-excluded_R2.01.fits | awk '{print $1}') == 'fa89e8fc4fd2868169c36b17d53bda05' ]]; then
            break
        fi
    done
fi

if [[ -f PR2-2015/COM_PCCS_857_R2.01.fits ]]
then
    echo 'You have file COM_PCCS_857_R2.01.fits'
else
    while true; do
        wget -cP PR2-2015/ 'https://irsa.ipac.caltech.edu/data/Planck/release_2/catalogs/COM_PCCS_857_R2.01.fits'
        if [[ $(md5sum PR2-2015/COM_PCCS_857_R2.01.fits | awk '{print $1}') == 'ac1fe918e02d5a984dff2c1f4395b3b6' ]]; then
            break
        fi
    done
fi

# ------------------------------------
# Thermal dust model (Planck 2015, GNILC)
if [[ -f PR2-2015/COM_CompMap_Dust-GNILC-Model-Opacity_2048_R2.01.fits ]]
then
    echo 'You have file COM_CompMap_Dust-GNILC-Model-Opacity_2048_R2.01.fits'
else
    while true; do
        wget -cP PR2-2015/ 'https://irsa.ipac.caltech.edu/data/Planck/release_2/all-sky-maps/maps/component-maps/foregrounds/COM_CompMap_Dust-GNILC-Model-Opacity_2048_R2.01.fits'
        if [[ $(md5sum PR2-2015/COM_CompMap_Dust-GNILC-Model-Opacity_2048_R2.01.fits | awk '{print $1}') == 'fc385c2ee5e82edf039cbca6e82d6872' ]]; then
            break
        fi
    done
fi

if [[ -f PR2-2015/COM_CompMap_Dust-GNILC-Model-Spectral-Index_2048_R2.01.fits ]]
then
    echo 'You have file COM_CompMap_Dust-GNILC-Model-Spectral-Index_2048_R2.01.fits'
else
    while true; do
        wget -cP PR2-2015/ 'https://irsa.ipac.caltech.edu/data/Planck/release_2/all-sky-maps/maps/component-maps/foregrounds/COM_CompMap_Dust-GNILC-Model-Spectral-Index_2048_R2.01.fits'
        if [[ $(md5sum PR2-2015/COM_CompMap_Dust-GNILC-Model-Spectral-Index_2048_R2.01.fits | awk '{print $1}') == '3bbe9869872184bfe9d00e61d448c8aa' ]]; then
            break
        fi
    done
fi

if [[ -f PR2-2015/COM_CompMap_Dust-GNILC-Model-Temperature_2048_R2.01.fits ]]
then
    echo 'You have file COM_CompMap_Dust-GNILC-Model-Temperature_2048_R2.01.fits'
else
    while true; do
        wget -cP PR2-2015/ 'https://irsa.ipac.caltech.edu/data/Planck/release_2/all-sky-maps/maps/component-maps/foregrounds/COM_CompMap_Dust-GNILC-Model-Temperature_2048_R2.01.fits'
        if [[ $(md5sum PR2-2015/COM_CompMap_Dust-GNILC-Model-Temperature_2048_R2.01.fits | awk '{print $1}') == '5ff018bd30e911b8c809116c0f80e842' ]]; then
            break
        fi
    done
fi

# ------------------------------------
# Thermal dust from GNILC
if [[ -f PR2-2015/COM_CompMap_Dust-GNILC-F353_2048_R2.00.fits ]]
then
    echo 'You have file COM_CompMap_Dust-GNILC-F353_2048_R2.00.fits'
else
    while true; do
        wget -cP PR2-2015/ 'https://irsa.ipac.caltech.edu/data/Planck/release_2/all-sky-maps/maps/component-maps/foregrounds/COM_CompMap_Dust-GNILC-F353_2048_R2.00.fits'
        if [[ $(md5sum PR2-2015/COM_CompMap_Dust-GNILC-F353_2048_R2.00.fits | awk '{print $1}') == 'ae18f38139dcd426df1e7fc6c93e7670' ]]; then
            break
        fi
    done
fi

if [[ -f PR2-2015/COM_CompMap_Dust-GNILC-F545_2048_R2.00.fits ]]
then
    echo 'You have file COM_CompMap_Dust-GNILC-F545_2048_R2.00.fits'
else
    while true; do
        wget -cP PR2-2015/ 'https://irsa.ipac.caltech.edu/data/Planck/release_2/all-sky-maps/maps/component-maps/foregrounds/COM_CompMap_Dust-GNILC-F545_2048_R2.00.fits'
        if [[ $(md5sum PR2-2015/COM_CompMap_Dust-GNILC-F545_2048_R2.00.fits | awk '{print $1}') == '96de58c0df7a82feebdcf61d4bc5dd24' ]]; then
            break
        fi
    done
fi

if [[ -f PR2-2015/COM_CompMap_Dust-GNILC-F857_2048_R2.00.fits ]]
then
    echo 'You have file COM_CompMap_Dust-GNILC-F857_2048_R2.00.fits'
else
    while true; do
        wget -cP PR2-2015/ 'https://irsa.ipac.caltech.edu/data/Planck/release_2/all-sky-maps/maps/component-maps/foregrounds/COM_CompMap_Dust-GNILC-F857_2048_R2.00.fits'
        if [[ $(md5sum PR2-2015/COM_CompMap_Dust-GNILC-F857_2048_R2.00.fits | awk '{print $1}') == '6dafe091af608429a2c567737e472347' ]]; then
            break
        fi
    done
fi

# ------------------------------------
# Thermal dust model (Planck 2015, Commander)
if [[ -f PR2-2015/COM_CompMap_dust-commander_0256_R2.00.fits ]]
then
    echo 'You have file COM_CompMap_dust-commander_0256_R2.00.fits'
else
    while true; do
        wget -cP PR2-2015/ 'https://irsa.ipac.caltech.edu/data/Planck/release_2/all-sky-maps/maps/component-maps/foregrounds/COM_CompMap_dust-commander_0256_R2.00.fits'
        if [[ $(md5sum PR2-2015/COM_CompMap_dust-commander_0256_R2.00.fits | awk '{print $1}') == '84913c676f5ce834fb54fa353bdaecaa' ]]; then
            break
        fi
    done
fi

# ------------------------------------
# ------------------------------------
# Planck 2018 release

if [[ -d 'PR3-2018' ]]
then
    echo 'You have directory PR3-2018'
else
    mkdir PR3-2018
fi

# ------------------------------------
# CMB SMICA
if [[ -f PR3-2018/COM_CMB_IQU-smica_2048_R3.00_full.fits ]]
then
    echo 'You have file COM_CMB_IQU-smica_2048_R3.00_full.fits'
else
    while true; do
        wget -cP PR3-2018/ 'https://irsa.ipac.caltech.edu/data/Planck/release_3/all-sky-maps/maps/component-maps/cmb/COM_CMB_IQU-smica_2048_R3.00_full.fits'
        if [[ $(md5sum PR3-2018/COM_CMB_IQU-smica_2048_R3.00_full.fits | awk '{print $1}') == 'ee2fc49a2eb70c2eca0d582e4aae5d05' ]]; then
            break
        fi
    done
fi

# ------------------------------------
# synchrotron polarization SMICA
if [[ -f PR3-2018/COM_CompMap_QU-synchrotron-smica_2048_R3.00_full.fits ]]
then
    echo 'You have file COM_CompMap_QU-synchrotron-smica_2048_R3.00_full.fits'
else
    while true; do
        wget -cP PR3-2018/ 'https://irsa.ipac.caltech.edu/data/Planck/release_3/all-sky-maps/maps/component-maps/foregrounds/COM_CompMap_QU-synchrotron-smica_2048_R3.00_full.fits'
        if [[ $(md5sum PR3-2018/COM_CompMap_QU-synchrotron-smica_2048_R3.00_full.fits | awk '{print $1}') == '7db4e6452836fe38e17ee0165afa1619' ]]; then
            break
        fi
    done
fi

# ------------------------------------
# Dust polarization SMICA
#if [[ -f PR3-2018/COM_CompMap_QU-thermaldust-smica_2048_R3.00_full.fits ]]
#then
#    echo 'You have file COM_CompMap_QU-thermaldust-smica_2048_R3.00_full.fits'
#else
#    echo 'You have not file COM_CompMap_QU-thermaldust-smica_2048_R3.00_full.fits, we will download it. '
#    wget -cP PR3-2018/ 'https://irsa.ipac.caltech.edu/data/Planck/release_3/all-sky-maps/maps/component-maps/foregrounds/COM_CompMap_QU-thermaldust-smica_2048_R3.00_full.fits'
#fi

# ------------------------------------
# Dust polarization COMMANDER
if [[ -f PR3-2018/COM_CompMap_QU-thermaldust-commander_2048_R3.00_full.fits ]]
then
    echo 'You have file COM_CompMap_QU-thermaldust-commander_2048_R3.00_full.fits'
else
    while true; do
        wget -cP PR3-2018/ 'https://irsa.ipac.caltech.edu/data/Planck/release_3/all-sky-maps/maps/component-maps/foregrounds/COM_CompMap_QU-thermaldust-commander_2048_R3.00_full.fits'
        if [[ $(md5sum PR3-2018/COM_CompMap_QU-thermaldust-commander_2048_R3.00_full.fits | awk '{print $1}') == '725008422f35706c89de105d52805ebd' ]]; then
            break
        fi
    done
fi

# ------------------------------------
# Dust polarization GNILC
if [[ -f PR3-2018/COM_CompMap_IQU_thermaldust-gnilc-unires_2048_R3.00.fits ]]
then
    echo 'You have file COM_CompMap_IQU_thermaldust-gnilc-unires_2048_R3.00.fits'
else
    while true; do
        wget -cP PR3-2018/ 'https://irsa.ipac.caltech.edu/data/Planck/release_3/all-sky-maps/maps/component-maps/foregrounds/COM_CompMap_IQU-thermaldust-gnilc-unires_2048_R3.00.fits'
        if [[ $(md5sum PR3-2018/COM_CompMap_IQU-thermaldust-gnilc-unires_2048_R3.00.fits | awk '{print $1}') == 'bf0a5adb8a1828614850c217751818da' ]]; then
            break
        fi
    done
fi

# ------------------------------------
# LFI sky maps
if [[ -f PR3-2018/LFI_SkyMap_030-BPassCorrected_1024_R3.00_full.fits ]]
then
   echo 'You have file LFI_SkyMap_030-BPassCorrected_1024_R3.00_full.fits'
else
    while true; do
        wget -cP PR3-2018/ 'https://irsa.ipac.caltech.edu/data/Planck/release_3/all-sky-maps/maps/LFI_SkyMap_030-BPassCorrected_1024_R3.00_full.fits'
        if [[ $(md5sum PR3-2018/LFI_SkyMap_030-BPassCorrected_1024_R3.00_full.fits | awk '{print $1}') == 'aee80baf3cc7ffff6911a7b984e8087b' ]]; then
            break
        fi
    done
fi

if [[ -f PR3-2018/LFI_SkyMap_044-BPassCorrected_1024_R3.00_full.fits ]]
then
   echo 'You have file LFI_SkyMap_044-BPassCorrected_1024_R3.00_full.fits'
else
    while true; do
        wget -cP PR3-2018/ 'https://irsa.ipac.caltech.edu/data/Planck/release_3/all-sky-maps/maps/LFI_SkyMap_044-BPassCorrected_1024_R3.00_full.fits'
        if [[ $(md5sum PR3-2018/LFI_SkyMap_044-BPassCorrected_1024_R3.00_full.fits | awk '{print $1}') == '0503588ddf774ee6597adeddd42d7e6f' ]]; then
            break
        fi
    done
fi

if [[ -f PR3-2018/LFI_SkyMap_070-BPassCorrected_1024_R3.00_full.fits ]]
then
   echo 'You have file LFI_SkyMap_070-BPassCorrected_1024_R3.00_full.fits'
else
    while true; do
        wget -cP PR3-2018/ 'https://irsa.ipac.caltech.edu/data/Planck/release_3/all-sky-maps/maps/LFI_SkyMap_070-BPassCorrected_1024_R3.00_full.fits'
        if [[ $(md5sum PR3-2018/LFI_SkyMap_070-BPassCorrected_1024_R3.00_full.fits | awk '{print $1}') == 'b81d8a2cded20bcf2e0ad9770929c145' ]]; then
            break
        fi
    done
fi

# ------------------------------------
# HFI sky maps
if [[ -f PR3-2018/HFI_SkyMap_100_2048_R3.01_full.fits ]]
then
   echo 'You have file HFI_SkyMap_100_2048_R3.01_full.fits'
else
    while true; do
        wget -cP PR3-2018/ 'https://irsa.ipac.caltech.edu/data/Planck/release_3/all-sky-maps/maps/HFI_SkyMap_100_2048_R3.01_full.fits'
        if [[ $(md5sum PR3-2018/HFI_SkyMap_100_2048_R3.01_full.fits | awk '{print $1}') == '2bb9d1c09f2e75ad1af014b60486e00e' ]]; then
            break
        fi
    done
fi

if [[ -f PR3-2018/HFI_SkyMap_143_2048_R3.01_full.fits ]]
then
    echo 'You have file HFI_SkyMap_143_2048_R3.01_full.fits'
else
    while true; do
        wget -cP PR3-2018/ 'https://irsa.ipac.caltech.edu/data/Planck/release_3/all-sky-maps/maps/HFI_SkyMap_143_2048_R3.01_full.fits'
        if [[ $(md5sum PR3-2018/HFI_SkyMap_143_2048_R3.01_full.fits | awk '{print $1}') == '4151e22bea78a318ea19a8040625f502' ]]; then
            break
        fi
    done
fi

if [[ -f PR3-2018/HFI_SkyMap_217_2048_R3.01_full.fits ]]
then
    echo 'You have file HFI_SkyMap_217_2048_R3.01_full.fits'
else
    while true; do
        wget -cP PR3-2018/ 'https://irsa.ipac.caltech.edu/data/Planck/release_3/all-sky-maps/maps/HFI_SkyMap_217_2048_R3.01_full.fits'
        if [[ $(md5sum PR3-2018/HFI_SkyMap_217_2048_R3.01_full.fits | awk '{print $1}') == '03101e030eb246e307c8e5e81c3cffd1' ]]; then
            break
        fi
    done
fi

if [[ -f PR3-2018/HFI_SkyMap_353_2048_R3.01_full.fits ]]
then
    echo 'You have file HFI_SkyMap_353_2048_R3.01_full.fits'
else
    while true; do
        wget -cP PR3-2018/ 'https://irsa.ipac.caltech.edu/data/Planck/release_3/all-sky-maps/maps/HFI_SkyMap_353_2048_R3.01_full.fits'
        if [[ $(md5sum PR3-2018/HFI_SkyMap_353_2048_R3.01_full.fits | awk '{print $1}') == '4d9fe8caeb3fe31e1732664bcaefa16a' ]]; then
            break
        fi
    done
fi

if [[ -f PR3-2018/HFI_SkyMap_545_2048_R3.01_full.fits ]]
then
    echo 'You have file HFI_SkyMap_545_2048_R3.01_full.fits'
else
    while true; do
        wget -cP PR3-2018/ 'https://irsa.ipac.caltech.edu/data/Planck/release_3/all-sky-maps/maps/HFI_SkyMap_545_2048_R3.01_full.fits'
        if [[ $(md5sum PR3-2018/HFI_SkyMap_545_2048_R3.01_full.fits | awk '{print $1}') == '4f031bdf33def2ab50aab30fc2308952' ]]; then
            break
        fi
    done
fi

if [[ -f PR3-2018/HFI_SkyMap_857_2048_R3.01_full.fits ]]
then
    echo 'You have file HFI_SkyMap_857_2048_R3.01_full.fits'
else
    while true; do
        wget -cP PR3-2018/ 'https://irsa.ipac.caltech.edu/data/Planck/release_3/all-sky-maps/maps/HFI_SkyMap_857_2048_R3.01_full.fits'
        if [[ $(md5sum PR3-2018/HFI_SkyMap_857_2048_R3.01_full.fits | awk '{print $1}') == '468b924c6a8777b3ac4ab282ed61bd24' ]]; then
            break
        fi
    done
fi

# ------------------------------------
# Transmission spectra
if [[ -f PR3-2018/LFI_RIMO_R3.31.fits ]]
then
    echo 'You have file LFI_RIMO_R3.31.fits'
else
    while true; do
        wget -cP PR3-2018/ 'https://irsa.ipac.caltech.edu/data/Planck/release_3/ancillary-data/LFI_RIMO_R3.31.fits'
        if [[ $(md5sum PR3-2018/LFI_RIMO_R3.31.fits | awk '{print $1}') == 'e21e3afd705fc7676b103c838dd22e70' ]]; then
            break
        fi
    done
fi

if [[ -f PR3-2018/HFI_RIMO_R3.00.fits ]]
then
    echo 'You have file HFI_RIMO_R3.00.fits'
else
    while true; do
        wget -cP PR3-2018/ 'https://irsa.ipac.caltech.edu/data/Planck/release_3/ancillary-data/HFI_RIMO_R3.00.fits'
        if [[ $(md5sum PR3-2018/HFI_RIMO_R3.00.fits | awk '{print $1}') == '32ac7c5e153264e5c83914ac3fe7da07' ]]; then
            break
        fi
    done
fi

# ------------------------------------
# ------------------------------------
# Thermal dust model from Melis O. Irfan et al. A&A 623, A21 (2019)
if [[ -d 'Irfan19' ]]
then
    echo 'You have directory Irfan19'
else
    mkdir Irfan19
fi

if [[ -f Irfan19/beta.fits ]]
then
    echo 'You have file beta.fits'
else
    while true; do
        wget -cP Irfan19/ 'https://cdsarc.cds.unistra.fr/ftp/J/A+A/623/A21/fits/beta.fits'
        if [[ $(md5sum Irfan19/beta.fits | awk '{print $1}') == '11304bc5e573ccd3f793a646dd9dbbbb' ]]; then
            break
        fi
    done
fi

if [[ -f Irfan19/tau.fits ]]
then
    echo 'You have file tau.fits'
else
    while true; do
        wget -cP Irfan19/ 'https://cdsarc.cds.unistra.fr/ftp/J/A+A/623/A21/fits/tau.fits'
        if [[ $(md5sum Irfan19/tau.fits | awk '{print $1}') == '2297934041f9f7314257a842cc6e334b' ]]; then
            break
        fi
    done
fi

if [[ -f Irfan19/temp.fits ]]
then
    echo 'You have file temp.fits'
else
    while true; do
        wget -cP Irfan19/ 'https://cdsarc.cds.unistra.fr/ftp/J/A+A/623/A21/fits/temp.fits'
        if [[ $(md5sum Irfan19/temp.fits | awk '{print $1}') == '029bdda19485295957e40460054cc2c5' ]]; then
            break
        fi
    done
fi

# ------------------------------------
# ------------------------------------
if [[ -d '2MBB_Meisner' ]]
then
    echo '2MBB_Meisner'
else
    mkdir 2MBB_Meisner
fi

# ------------------------------------
# Two MBB model from ApJ 798 88 (2015) http://planck.skymaps.info/
if [[ -f 2MBB_Meisner/planck_2comp.fits ]]
then
    echo 'You have file planck_2comp.fits'
else
    while true; do
        wget -cP 2MBB_Meisner/ 'https://faun.rc.fas.harvard.edu/ameisner/planckdust/planck_2comp.fits'
        if [[ $(md5sum 2MBB_Meisner/planck_2comp.fits | awk '{print $1}') == 'ec382f711c8669c78966caac15acde69' ]]; then
            break
        fi
    done
fi

# ------------------------------------
# ------------------------------------
if [[ -d 'Planck_DL07' ]]
then
    echo 'Planck_DL07'
else
    mkdir Planck_DL07
fi

# ------------------------------------
# DL07 model from Planck intermediate results XXIX
if [[ -f Planck_DL07/COM_CompMap_Dust-DL07-ModelFluxes_2048_R2.00.fits ]]
then
    echo 'You have file COM_CompMap_Dust-DL07-ModelFluxes_2048_R2.00.fits'
else
    while true; do
        wget -cP Planck_DL07/ 'https://irsa.ipac.caltech.edu/data/Planck/release_2/all-sky-maps/maps/component-maps/foregrounds/COM_CompMap_Dust-DL07-ModelFluxes_2048_R2.00.fits'
        if [[ $(md5sum Planck_DL07/COM_CompMap_Dust-DL07-ModelFluxes_2048_R2.00.fits | awk '{print $1}') == 'f576a1f9b4ac4fb02d6d5f89a78966e2' ]]; then
            break
        fi
    done
fi

if [[ -f Planck_DL07/COM_CompMap_Dust-DL07-Parameters_2048_R2.00.fits ]]
then
    echo 'You have file COM_CompMap_Dust-DL07-Parameters_2048_R2.00.fits'
else
    while true; do
        wget -cP Planck_DL07/ 'https://irsa.ipac.caltech.edu/data/Planck/release_2/all-sky-maps/maps/component-maps/foregrounds/COM_CompMap_Dust-DL07-Parameters_2048_R2.00.fits'
        if [[ $(md5sum Planck_DL07/COM_CompMap_Dust-DL07-Parameters_2048_R2.00.fits | awk '{print $1}') == '8aa87e9b9659f9a3d0195559e935a080' ]]; then
            break
        fi
    done
fi

# ------------------------------------
# ------------------------------------
if [[ -d 'SRoll' ]]
then
    echo 'SRoll'
else
    mkdir SRoll
fi


# ------------------------------------
# 
if [[ -f SRoll/100GHz_DUSTMODEL_INTENSITY.fits ]]
then
    echo 'You have file 100GHz_DUSTMODEL_INTENSITY.fits'
else
    while true; do
        wget -cP SRoll/ ftp://sroll20.ias.u-psud.fr/sroll22_data/100GHz_DUSTMODEL_INTENSITY.fits
        if [[ $(md5sum SRoll/100GHz_DUSTMODEL_INTENSITY.fits | awk '{print $1}') == '7aa1712fbf5636669487c6d33881053d' ]]; then
            break
        fi
    done
fi

if [[ -f SRoll/143GHz_DUSTMODEL_INTENSITY.fits ]]
then
    echo 'You have file 143GHz_DUSTMODEL_INTENSITY.fits'
else
    while true; do
        wget -cP SRoll/ ftp://sroll20.ias.u-psud.fr/sroll22_data/143GHz_DUSTMODEL_INTENSITY.fits
        if [[ $(md5sum SRoll/143GHz_DUSTMODEL_INTENSITY.fits | awk '{print $1}') == '177ad68d51387f87475e841ca3c3ae05' ]]; then
            break
        fi
    done
fi

if [[ -f SRoll/217GHz_DUSTMODEL_INTENSITY.fits ]]
then
    echo 'You have file 217GHz_DUSTMODEL_INTENSITY.fits'
else
    while true; do
        wget -cP SRoll/ ftp://sroll20.ias.u-psud.fr/sroll22_data/217GHz_DUSTMODEL_INTENSITY.fits
        if [[ $(md5sum SRoll/217GHz_DUSTMODEL_INTENSITY.fits | awk '{print $1}') == 'c5f42c67c78c8ae928d3cea4b370d78e' ]]; then
            break
        fi
    done
fi

if [[ -f SRoll/353GHz_DUSTMODEL_INTENSITY.fits ]]
then
    echo 'You have file 353GHz_DUSTMODEL_INTENSITY.fits'
else
    while true; do
        wget -cP SRoll/ ftp://sroll20.ias.u-psud.fr/sroll22_data/353GHz_DUSTMODEL_INTENSITY.fits
        if [[ $(md5sum SRoll/353GHz_DUSTMODEL_INTENSITY.fits | awk '{print $1}') == 'ed0f81d4cd2ef3f3929adcb04404fdf9' ]]; then
            break
        fi
    done
fi

if [[ -f SRoll/545GHz_DUSTMODEL_INTENSITY.fits ]]
then
    echo 'You have file 545GHz_DUSTMODEL_INTENSITY.fits'
else
    while true; do
        wget -cP SRoll/ ftp://sroll20.ias.u-psud.fr/sroll22_data/545GHz_DUSTMODEL_INTENSITY.fits
        if [[ $(md5sum SRoll/545GHz_DUSTMODEL_INTENSITY.fits | awk '{print $1}') == 'df18a5f2ae483eb69523cc09899f14e6' ]]; then
            break
        fi
    done
fi

if [[ -f SRoll/857GHz_DUSTMODEL_INTENSITY.fits ]]
then
    echo 'You have file 857GHz_DUSTMODEL_INTENSITY.fits'
else
    while true; do
        wget -cP SRoll/ ftp://sroll20.ias.u-psud.fr/sroll22_data/857GHz_DUSTMODEL_INTENSITY.fits
        if [[ $(md5sum SRoll/857GHz_DUSTMODEL_INTENSITY.fits | awk '{print $1}') == '1ffcd9c19ef36669345a131e4416b9b0' ]]; then
            break
        fi
    done
fi
