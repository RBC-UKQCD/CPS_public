# Note that this shebang line may not work everywhere, as some 
# Unixies put perl elsewhere, e.g. /usr/bin/perl.
#
# More portable to call this script via:
#   perl regression.pl
#
#  The location of perl could be put in by the autoconf system
#--------------------------------------------------------------------

#
# Configurations: These are platform/config-dependant:
#

#
# General defs of versions and trusted versions:
#

# directory with good test data in it
$test_dir = "v5_0_0" ;
@test_names = ( "f_hmd" , "f_hmd_dwfso", "g_hb" , "g_hmd" , "s_spect" ,  "w_spect" , "f_wilson_eig" ,  "xi_spect_gsum" , "f_stag_pbp" ) ; 

##@test_names = ( "f_stag_pbp" , "f_wilson_pbp" , "f_clover_pbp" , "f_dwf_pbp" , "f_dwfso_pbp" ); 
#------

$error_tol = 0.001; 


if( 'yes' eq 'yes' ) {
# This makes it run autoconf version:
  $machine = 'i686-pc-linux-gnu';
  $parallel = 'yes';
  $compiler = '@CC@';
  $executable = "regression_cps.x";
} else {
# This is right for the parallel QCDSP version:
  $machine = "qcdsp";
  $parallel = "yes";
  $compiler = "tcpp";
  $executable = "qcdsp.out";
}


#
# Configurations: These should be relevant always:
#

#The name of the shell script this will create:
$shellscript = "regression.sh";
$output_dir = "regressions";


# Variations:
#---------------------------------
if( $machine eq "qcdsp" && $compiler eq "tcpp" ) {
    $exec_prefix = "qrun ./";
} else {
    $exec_prefix = "./";
}

#--------------------------------------------------------------------

# Open the file which will contain the shell script which will run the tests:
open SHOUT,">$shellscript" or die "Could not open $shellscript!\n";

# Loop over each test in the series:

# clean up any old output
    print SHOUT 'rm '.$output_dir.'/*.dat '."\n";

foreach $tni (@test_names){
    print SHOUT '#------------------------------------'."\n";

# Construct the names of the files to store the results:
	$res_file = "../".$output_dir."/stdio.".$tni.".dat";
	$dat_file = "../".$output_dir."/data.".$tni.".dat";

# Construct the name of the file which holds the original results:

# If we are on qcdsp, reset the machine before beginning the next run:
	if( $machine eq "qcdsp" && $compiler eq "tcpp" ) {
	    print SHOUT "echo Rebooting the board...\n";
	    print SHOUT "qreset_boot > qreset.log\n";
	}

# State which test is about to run, run it, and store the results:
	print SHOUT "cd ".$tni."\n";
# Remove any files associated with previous runs (i.e. *.dat)
	print SHOUT 'rm -f '.$exec_cmd.' *.dat'."\n";
#       compile the code (do not echo the commands)
	print SHOUT 'make -s -f Makefile_regression'."\n";
#       check that the executable has been created
	print SHOUT 'if test ! -f '.$exec_cmd.' ; then'."\n";
	print SHOUT 'echo '.$tni." COMPILATION FAILURE\n";
	print SHOUT "cd ..\n";

	print SHOUT "else\n";


# Run the program:
	$exec_cmd = $exec_prefix.$executable;
##	print SHOUT "echo Running $exec_cmd...\n";
	print SHOUT "$exec_cmd > $res_file\n";
# Grab all of the *.dat output files and sling them into a single file:
	print SHOUT "perl ../combine_files.pl *.dat > $dat_file\n";

# Inform the user about the directory being looked at:
##	print "[$tni]: '$exec_cmd'.\n";

# Diff the output and report if the output differs from that expected:
# Drop back to the test directory:
	print SHOUT "cd ..\n";
	print SHOUT "cd regressions ; perl check_data.pl $tni  $error_tol  $test_dir\n" ;
	print SHOUT "cd ..\n" ;

	print SHOUT "fi\n";

# Loop over to the next test:
}

print SHOUT "echo ------------------------------\n";
print SHOUT "echo DISCLAIMER\n";
print SHOUT "echo Please also test this code on \n";
print SHOUT "echo a physical system before using in \n";
print SHOUT "echo production runs \n";
print SHOUT "echo ------------------------------\n";



# close the script file
close SHOUT;

exit;
