
set curdir = `pwd`
set hostdir = /host1/1MB/host/lishu

set config = 2325
while ( -e $hostdir/ckpoint/Ckpoint_lat.dsp.full.$config )
	echo "Running with Config # $config"
	qrun QCDOC.x ckpoint/Ckpoint_lat.dsp.full.$config
	
	echo "Saving data files"
	cd $hostdir
	foreach x ( *.dat )
		mv $x data_npr/$x.$config
	end
	cd $curdir

	@ config += 5
end



