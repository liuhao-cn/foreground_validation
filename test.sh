#!/bin/bash
# Run Python scripts in the original directory
echo -e "\nRunning Python scripts..."
source myenv/bin/activate
# Verify installation
echo -e "\nVerifying installation:"
python -c "import pip;             print('pip             version:', pip.__version__)"
python -c "import sys;             print('Python          version:', sys.version)"
python -c "import pysm3;           print('PySM3           version:', pysm3.__version__)"
python -c "import healpy;          print('healpy          version:', healpy.__version__)"
python -c "import numpy;           print('numpy           version:', numpy.__version__)"
python -c "import scipy;           print('scipy           version:', scipy.__version__)"
python -c "import astropy;         print('astropy         version:', astropy.__version__)"
python -c "import matplotlib;      print('matplotlib      version:', matplotlib.__version__)"
python -c "import pandas;          print('pandas          version:', pandas.__version__)"