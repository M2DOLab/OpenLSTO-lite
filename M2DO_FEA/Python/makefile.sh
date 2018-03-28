#!/bin/bash 

python FEAsetup.py build
cp build/lib*/pyBind.so ./

#cp pyBind.so ~/Dropbox/packages/topOpt_MDO/FEM2D_v2/run_SIMP
#cp pyBind.so ~/Dropbox/packages/topOpt_MDO/FEM2D_v2
#cp pyBind.so ~/Desktop/WIP/python_lsm
# python test.py
