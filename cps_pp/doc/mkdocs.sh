# Runs doxygen to produce the reference manual and the user guide and 
# installdox to get the cross-references right.

# This script is really meant to be called from Makefile_cps and so the
# CPS directory structure is assumed e.g. the ref and usr directories 
# and the source are in all the right places. However, we can allow it to 
# be called as  a stand-alone script from the doc directory.

topdir=`echo $PWD | sed -e 's:/doc$::'`
docdir=$topdir/doc

# this does not always work :(

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

# doxygen doesn't always manage to strip paths so this is a failsafe hack

for f in $docdir/ref/html/files.html `ls $docdir/ref/html/*8[Ch]*.html`
do
	sed -e s:$topdir/:: $f > /tmp/painful
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








