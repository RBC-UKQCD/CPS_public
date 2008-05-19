#!/bin/tcsh
set branch=v4_9_17-pica_1
echo cvs update -j $branch $1
cvs update -j $branch $1
echo cvs diff -r $branch $1
cvs diff -r $branch $1
