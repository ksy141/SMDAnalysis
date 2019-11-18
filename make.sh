#!/bin/bash

cd compute
python density_help_setup.py build_ext --inplace
rm -rf build *.c *.html
cd ..

