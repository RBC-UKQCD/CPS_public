# $Id: mkdocs.sh,v 1.2 2003-07-30 16:31:11 zs Exp $
# Runs doxygen to produce the reference manual and the user guide and 
# installdox to get the cross-references right.

# This script is really meant to be called from Makefile_cps and so the
# CPS directory structure is assumed e.g. the ref and usr directories 
# and the source are in all the right places.
# However it can still be useful to use
# this as a stand-alone script. The bottom line is that this script expects 
# to be called from its own directory or its parent directory.


if [ `basename $PWD` == 'doc' ]
then
    docdir=`pwd -P`
    topdir=`echo $docdir | sed -e 's:/doc$::'`
else
    topdir=`pwd -P` 
    docdir=$topdir/doc
fi




# this is not supported for doxygen pre 1.2.17 or thereabouts

sed -e "s:^STRIP_FROM_PATH.*:STRIP_FROM_PATH $topdir/:" $docdir/ref/doxygen.cfg > /tmp/tempdox
mv /tmp/tempdox $docdir/ref/doxygen.cfg

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

# hack if this doxygen version cannot strip paths  

for f in `ls $docdir/ref/html/*.html`
do
	sed -e s:$topdir/::g $f > /tmp/painful
	mv /tmp/painful $f
done

# run doxygen again to pick up all cross references

cd $docdir/usr 
doxygen doxygen.cfg		 

# run installdox to edit the cross references 

for here in usr ref  
do
	if [ $here = usr ]
	then
		there=ref
	else
		there=usr
	fi
	cd $docdir/$here/html	    
	perl installdox -l doxtag@../../$there/html
done	

echo "$0: Finished!"








