#
# Routines to run on different platforms
#
#


$qcdoc_final = sub { 
# closing XTERM hangs because it waits for the application to end
system("kill -9 ". $pid_xterm ); 

};


$qcdoc_run = sub { 
    my $exec = $executable ;

    my $log = $exec."_log" ;
    my $out = $exec."_stdout" ;
    my $out_err = $exec."_stderr" ;


    my $pid = open(QCSH, "| ".$qcsh )  or die("Could not open $qcsh") ;
    print QCSH "qinit ".$MACH."  >& $log \n"   ;
    print QCSH "qpartition_connect -p 0  >& $log\n"   ;
    print QCSH "qreset_boot  >& $log\n"   ;
    print QCSH "qdiscover  >& $log\n"   ;
    print QCSH "qpartition_remap $qpartition_remap  >& $log\n"   ;
    print QCSH "qrun ".$exec." >& $log\n"   ;
    print QCSH "qnodes_print -b stdout >& ".$out."\n"   ;
    print QCSH "qnodes_print -b stderr >& ".$out_err."\n"   ;
    print QCSH "qdetach\n"   ;
    close (QCSH) ;

};


$qcdoc_setup = sub { 
    $MACH = "NONE" ; 
    $qpartition_remap = "-X0 -Y1 -Z2 -T3 -S4 -W5" ; 
    my $qtmp = $qpartition_remap ; 
    my $exec = $executable ;


    if( -f "params.pl" )
    {
	require 'params.pl' ; 
    }
    else
    {
	die "ERROR: Need to get the machine name from params.pl \n"  ; 
    }

    if( $MACH eq "NONE" )
    {
	die "ERROR: Need to set \$MACH in params.pl" . 
	    "e.g. \$MACH =   \"ssbp4/08node/mach03\" ;\n" ; 
    }

    if( $qtmp eq $qpartition_remap)
    {
	warn "WARNING: \$qpartition_remap not set in params.pl\n" .  
	     "Using single node qpartition_remap\n " .  
	     "\$qpartition_remap = $qpartition_remap\n" ; 
    }


    my $log = $exec."_log" ;
    my $out = $exec."_stdout" ;
    my $out_err = $exec."_stderr" ;


    print "The logfile is written to ".$log."\n" ; 
    print "Standard output is written to ".$out."\n" ; 
    print "Standard error is written to ".$out_err."\n" ; 


    $pid_xterm = 1 + open(XTERM, "| xterm -T \"Window for daemon output\" -e qdaemon --machine ".$MACH)  or die("ERROR: Could not open daemon xterm") ;

    sleep (6) ;

    $qcsh = "qcsh" ;

};



#
# collect subroutines with machines
#


%machine_names = (
	 QCDOC => {
                       setup         => $qcdoc_setup ,
                       run      => $qcdoc_run , 
                       exec      => "QCDOC.x" , 
                       final         => $qcdoc_final ,
                       output    => "/space/$ENV{'USER'}/" , 
		   },
	  NOARCH    => {
                       setup         => $do_nothing ,
                       run      => $scalar_run , 
                       final         => $do_nothing ,
                       exec      => "NOARCH.x" , 
                       output    => "" , 
		   },
	 );





# required for require
1 ; 
