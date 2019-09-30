#!/bin/bash

# python -m nuitka -o pixy --nofollow-imports --plugin-enable=numpy --plugin-enable=qt-plugins --plugin-enable=torch --plugin-enable=pylint-warnings  $RECIPE_DIR/src/pixy.py

version=$1

mkdir -p archive
mkdir -p release

cp pixy.py release/pixy.py
cp meta.yaml release/meta.yaml 
cp build.sh release/build.sh
cp -r src release/src

tar czfv ${version}.tar.gz release 
mv ${version}.tar.gz archive/${version}.tar.gz 

rm -r release