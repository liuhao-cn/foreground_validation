#!/bin/bash

# Save the original working directory
ORIGINAL_DIR=$(pwd)

# Reference
# 2015 ApJ 798 88
# http://planck.skymaps.info
# https://doi.org/10.1093/mnras/stx949
# http://doi.org/10.21105/joss.03783
# https://pysm3.readthedocs.io/en/latest/
# https://doi.org/10.1051/0004-6361/201525659

# Completely remove the old virtual environment
rm -rf myenv

# Check for available Python 3.9+ versions
if command -v python3.9 &> /dev/null; then
    PYTHON_CMD="python3.9"
elif command -v python3.10 &> /dev/null; then
    PYTHON_CMD="python3.10"
elif command -v python3 &> /dev/null; then
    # Check if default python3 version is 3.9 or higher
    PYTHON_VERSION=$(python3 -c "import sys; print(sys.version_info.major, sys.version_info.minor)" 2>/dev/null)
    if [[ $PYTHON_VERSION =~ ^3\ ([9-9]|[1-9][0-9])$ ]]; then
        PYTHON_CMD="python3"
    else
        # Python 3.9+ not found - install to local directory
        echo "Python 3.9+ not found. Installing to ~/.local/python3.9..."
        
        # Save current directory before changing
        CURRENT_DIR=$(pwd)
        
        # Create temporary directory
        TEMP_DIR=$(mktemp -d)
        cd "$TEMP_DIR"
        
        # Download and compile Python 3.9
        wget -q https://www.python.org/ftp/python/3.9.18/Python-3.9.18.tgz
        tar xzf Python-3.9.18.tgz
        cd Python-3.9.18
        
        # Configure and install
        ./configure --prefix="$HOME/.local/python3.9" --enable-optimizations --with-ensurepip=install
        make -s -j$(nproc)
        make -s install
        
        # Clean up temporary files and return to original directory
        cd "$CURRENT_DIR"
        rm -rf "$TEMP_DIR"
        
        # Set path to newly installed Python
        PYTHON_CMD="$HOME/.local/python3.9/bin/python3.9"
        
        # Add to PATH environment variable
        if ! grep -q "\.local/python3\.9/bin" ~/.bashrc; then
            echo 'export PATH="$HOME/.local/python3.9/bin:$PATH"' >> ~/.bashrc
            source ~/.bashrc
        fi
    fi
else
    echo "Error: Python 3 not found"
    exit 1
fi

# Create virtual environment using detected Python
echo "Creating virtual environment with $PYTHON_CMD..."
"$PYTHON_CMD" -m venv myenv

# Activate environment
source myenv/bin/activate

# Upgrade pip and install required packages
echo "Installing required packages..."
python -m pip install --upgrade pip setuptools --quiet
pip install numpy matplotlib healpy scipy astropy numba toml pysm3 pandas --quiet

# Verify installation
echo -e "\nVerifying installation:"
python -c "import sys; print('Python version:', sys.version)"
python -c "import pysm3; print('PySM3 version:', pysm3.__version__)"

echo -e "\nVirtual environment created successfully!"
echo "Activation command: source myenv/bin/activate"

# Return to original directory before running Python scripts
cd "$ORIGINAL_DIR"

# Verify files exist in current directory
echo -e "\nFiles in current directory ($(pwd)):"
ls -l *.py