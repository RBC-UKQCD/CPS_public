#
# List of tests and subroutines associated with these tests
#
#


%test_names = (
	 w_spect => {
                       get_xml         => $combine_file ,
                       check_data      => "YES" , 
		   },
	  t_asqtad_pion    => {
                       get_xml         => $do_nothing ,
                       check_data      => "YES" , 
		   },
	 t_gauge_rot    => {
                     get_xml         => $do_nothing ,
                     check_data      => "YES" , 
		   },
	 asqtad_hmc    => {
                       get_xml         => $combine_file ,
                      check_data      => "YES" , 
	   },
##	  asqtad_hmd_r   => {
##                       get_xml         => $do_nothing ,
##                       check_data      => "NO" , 
##		   },
##
	 );





# required for require
1 ; 
