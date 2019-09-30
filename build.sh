#!/bin/bash

# python -m nuitka -o pixy --nofollow-imports --plugin-enable=numpy --plugin-enable=qt-plugins --plugin-enable=torch --plugin-enable=pylint-warnings  $RECIPE_DIR/src/pixy.py

mkdir -p $PREFIX/bin

cp $RECIPE_DIR/src/pixy.py  $PREFIX/bin/pixy.py 

