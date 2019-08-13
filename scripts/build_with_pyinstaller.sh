#!/bin/bash

# build executable with pyinstaller
# conda install --yes -c conda-forge pyinstaller

rm -r dist
rm -r build

LD_LIBRARY_PATH=/home/ksamuk/anaconda3/lib pyinstaller --onefile pixy.py