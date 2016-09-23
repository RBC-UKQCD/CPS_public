# $Id: mkdocs.sh,v 1.5 2004/09/02 16:51:48 zs Exp $

# Runs doxygen to produce the reference manual and the user guide and 
# installdox to get the cross-references right.

# This script is really meant to be called from Makefile_cps and so the
# CPS directory structure is assumed e.g. the ref and usr directories 
# and the source are in all the right places.
# However it can still be useful to use
# this as a stand-alone script. The bottom line is that this script expects 
# to be called from its own directory or its parent directory.


if [ `basename $PWD` = 'doc' ]
then
    docdir=`pwd -P`
    topdir=`echo $docdir | sed -e 's|/doc$||'`
else
    topdir=`pwd -P` 
    docdir=$topdir/doc
fi

if [ $# -eq 1 ] 
then
    srcdir=$1
else
    srcdir="."
fi

tempdoxcfg=/tmp/dox$$

# The version number:
# I assume that the tag is of the form <something>4_8_0<something>

version=`echo '$Name: v5_0_16_hantao_io_test_v7 $' | sed -e 's/[^1234567890_]//g' -e 'y/_/./'`

for here in usr ref 
do
	sed -e "s/^PROJECT_NUMBER.*/PROJECT_NUMBER = $version/" $docdir/$here/doxygen.cfg > $tempdoxcfg
	mv $tempdoxcfg $docdir/$here/doxygen.cfg
done

sed -e "s|^STRIP_FROM_PATH.*|STRIP_FROM_PATH = $srcdir/|" -e "s|^INPUT *=.*|INPUT = ../../config.h $srcdir/include $srcdir/src|" $docdir/ref/doxygen.cfg > $tempdoxcfg
mv $tempdoxcfg $docdir/ref/doxygen.cfg



# run doxygen  

for here in usr ref 
do
	if [ $here = usr ]
	then
		there=ref
	else
		there=usr
	fi

	cd $docdir/$here
	doxygen doxygen.cfg		
done

# Hack if this doxygen version cannot strip paths  
# for f in `ls $docdir/ref/html/*.html`
# do
# 	sed -e s|$topdir/||g $f > /tmp/painful
# 	mv /tmp/painful $f
# done

# run doxygen again to pick up all cross references

cd $docdir/usr 
doxygen doxygen.cfg		 

# run installdox to edit the cross references 

cd html
perl installdox -l doxtag@../../ref/html


# Make the examples

cd $docdir/examples
for d in t_asqtad_pion s_spect t_gauge_rot
do
	sed -e "s|^OUTPUT_DIRECTORY.*|OUTPUT_DIRECTORY = $d|" \
	-e "s|^INPUT .*|INPUT = $srcdir/tests/$d|" doxygen.cfg > $tempdoxcfg
	doxygen $tempdoxcfg
done
rm $tempdoxcfg

echo "$0: Finished!"








