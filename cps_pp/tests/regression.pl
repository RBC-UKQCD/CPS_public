#include<config.h>
CPS_START_NAMESPACE
#!/usr/local/bin/perl
# Note that this shebang line may not work everywhere, as some 
# Unixies put perl elsewhere, e.g. /usr/bin/perl.
#
#--------------------------------------------------------------------

#
# Configurations: These are platform/config-dependant:
#

# General defs of versions and trusted versions:
$version = "4_1_0";
$trusted_version = "4_0_0";


if( 'yes' eq 'yes' ) {
# This makes it run autoconf version:
  $machine = 'sparc-sun-solaris2.8';
  $parallel = 'yes';
  $compiler = 'mpCC';
  $executable = "main.gnu.out";
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
#$output_dir = "Qregressions";
$output_dir = "regressions";

# There are three sets of tests...
#---------------------------------
$test_str[0] = "PsiBarPsi tests:";
$test_dirs[0] = "f_stag_pbp f_wilson_pbp f_clover_pbp f_dwf_pbp f_dwfso_pbp";

#---------------------------------
$test_str[1] = "General tests without stdio output:";
$test_dirs[1] = "f_hmd f_hmd_dwfso g_hb g_hmd s_spect w_spect f_wilson_eig xi_spect_gsum";
#[DIES] threept

#---------------------------------
#$test_str[2] = "General tests with stdio output:";
#$test_dirs[2] = "";
#[DIES] fix_gauge

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

# Loop over all series of tests:
for( $tsi = 0; $tsi < scalar(@test_dirs); $tsi++ ) {

# Print out the title of this series of tests:
    print "$test_str[$tsi]\n";

# Loop over each test in the series:
    @test_names = split " ", $test_dirs[$tsi];
    for( $tni = 0; $tni < scalar(@test_names); $tni++ ) {
	print SHOUT '#------------------------------------'."\n";

# Construct the names of the files to store the results:
	$res_file = "../".$output_dir."/stdio.".$test_names[$tni].".dat";
	$dat_file = "../".$output_dir."/data.".$test_names[$tni].".dat";

# Construct the name of the file which holds the original results:

# If we are on qcdsp, reset the machine before beginning the next run:
	if( $machine eq "qcdsp" && $compiler eq "tcpp" ) {
	    print SHOUT "echo Rebooting the board...\n";
	    print SHOUT "qreset_boot > qreset.log\n";
	}

# State which test is about to run, run it, and store the results:
	print SHOUT "echo In ".$test_names[$tni]."\n";
	print SHOUT "cd ".$test_names[$tni]."\n";
# Remove any files associated with previous runs (i.e. *.dat)
	print SHOUT 'rm *.dat'."\n";
# Run the program:
	$exec_cmd = $exec_prefix.$executable;
	print SHOUT "echo Running $exec_cmd...\n";
	print SHOUT "$exec_cmd > $res_file\n";
# Grab all of the *.dat output files and sling them into a single file:
	print SHOUT "more *.dat > $dat_file\n";

# Inform the user about the directory being looked at:
	print "[$tni]: '$exec_cmd'.\n";

# Diff the output and report if the output differs from that expected:
	print SHOUT "echo Checking the output...\n";

# Drop back to the test directory:
	print SHOUT "cd ..\n";

# Loop over to the next test:
	}

# Loop over to the next test series:
}

# Tell the user how to use the script:
print "To run the test programs, use 'source $shellscript'\n";

# close the script file
close SHOUT;

exit;

CPS_END_NAMESPACE
