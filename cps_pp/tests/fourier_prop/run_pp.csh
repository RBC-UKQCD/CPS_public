
set curdir = `pwd`
set hostdir = /host1/1MB/host/lishu

set config = 5
while ( -e $hostdir/ckpoint/ckpoint_lat.dsp.full.$config )
	echo "Running with Config # $config"
	qrun QCDOC.x ckpoint/ckpoint_lat.dsp.full.$config
	
	echo "Saving data files"
	cd $hostdir
	foreach x ( *.dat )
		mv $x data/$x.$config
	end
	cd $curdir

	@ config += 5
end



