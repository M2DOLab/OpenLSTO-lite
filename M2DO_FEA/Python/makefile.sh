#!/bin/bash 

python FEAsetup.py build
cp build/lib*/pyBind.so ./

# python test.py