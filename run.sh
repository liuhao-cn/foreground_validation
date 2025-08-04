#!/bin/bash
# Run Python scripts in the original directory
echo -e "\nRunning Python scripts..."
source myenv/bin/activate
python3 -u Inpaint_Rebeam_frequency_maps.py 2000
python3 -u Dust_data_model.py
python3 -u Synchrotron_data_model.py
python3 -u free_free_data_model.py
python3 -u AME_data_model.py
python3 -u Residual.py
python3 -u Statistics.py
python3 -u Plots.py