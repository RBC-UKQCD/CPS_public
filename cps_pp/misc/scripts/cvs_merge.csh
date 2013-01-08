#!/bin/tcsh
set branch=v5_0_12_nuc3pt_production_v2
echo cvs update -j $branch $1
cvs update -j $branch $1
echo cvs diff -r $branch $1
cvs diff -r $branch $1
