#!/bin/tcsh
sed --quiet -e '1 p'  $1 >temp
echo BGL_DIR = noarch >>temp
sed --quiet -e '2,$ p'  $1 >>temp
echo $1
cat temp
cp temp $1
