#!/usr/bin/perl
#
# Perl script to compare an output file
# against some simple output.
#
#


if( $#ARGV  != 2 )
{
    print "Usage: $0 test_name   tolerance directory_with_test_data\n";
    exit ;
}



#  load in the correct data
$test_name = $ARGV[0] ;
$tol =   $ARGV[1] ;
$test_dir = $ARGV[2] ;

# create the correct filenames

$correct_file = $test_dir."/data.".$test_name.".test" ; 
$data_file    = "data.".$test_name.".dat" ; 
$verbose_log  = "data.".$test_name.".checklog" ; 

if( ! -f $correct_file )
{
    print $test_name . "  NO_CORRECT_OUTPUT \n" ;
    exit 0; 
}



open(OUT,">$verbose_log")  || die "I can not open the file $verbose_log" ;

print OUT "Test = " . $test_name . "\n" ; 
print OUT "Tolerance = " . $tol . "\n" ; 

print OUT "\n" ; 

print OUT "Previous data in the file: ".$correct_file."\n" ; 
print OUT "Data to be tested in the file: ".$data_file."\n" ; 



&load_check_data($correct_file); 
&check_load_check_data();


&check_data($data_file,$test_name) ;  

exit 0 ; 
# ---------- end of top level program ----------

#
#  loop through the data file
#

sub check_data
{
    my ($logfile,$tag) = @_ ; 
    
#Open the data file
    open(LOG,$logfile)  || die "I can not open the file $logfile" ;
    print OUT "\nStarting to check the data\n" ;

    my @all_lines = <LOG> ; 

    my $pass_fail = "FAIL" ; 

    my $pt = 0 ; 
    foreach $line (@all_lines) 
    {
	chop($line) ; 

	foreach my $name (keys  %tag_store)
	{
	    if( $name eq $line )
	    {
		print OUT "\n".$name. "   MATCH\n" ; 
		foreach my $my_tag (@{$tag_store{$name}} )
		{
		    ($row_from_name, $column, $value) = split(':',$my_tag) ; 


		    my $target_pt   = $pt + $row_from_name ;
		    my $target_line = $all_lines[$target_pt] ; 
		    chop ($target_line) ; 
		    print OUT $pt . "  " . $row_from_name ."  " .$target_pt . "\n" ;
		    print OUT $target_line . "\n" ;
		    @line_contents = split(' ',$target_line);
		    my $test_value = $line_contents[$column ] ;
		    print OUT "Value from test= ". $test_value ."\n";
		    print OUT "Value from previous/good data = ". $value ."\n";
		    my $diff = abs($test_value - $value)  ;
		    print OUT "Difference = ".$diff . "\n" ; 

		    my $rel = $diff / $test_value ; 

		    if( $rel > $tol ) 
		    {
			print OUT  "FAIL\n" ;
			$pass_fail = "FAIL" ;
		    }
		    else
		    {
			print  OUT "PASS\n" ;
			$pass_fail = "PASS" ;
		    }

		} # end of loop over test data for this particle

	    }

	}  # end loop over names

    $pt += 1 ;
    } # end loop over lines in file


#
#  write to standard output
#
print $tag . "    "  . $pass_fail  ."\n" ; 

} #  sub check_data



#  Load the data to check the file against
#
#
#

sub load_check_data
{
    my ($correct_data) = @_ ; 

    open(GOOD,$correct_data)  || die "I can not open the file $correct_data" ;
    my $line ;

    while(  $line = <GOOD> )
    {
	chop ($line) ; 

	my $name, $column, $row_from_name, $value ;
	($name, $row_from_name, $column, $value) = split(':',$line) ; 
	my $tag = $row_from_name .":". $column .":". $value ;

	# store the test data
	push @{ $tag_store{$name} }, $tag ;


    }

    close(GOOD) ;

}  # sub load_check_data


##
##  Check the loaded correct data
##
##

sub check_load_check_data
{

    my $name , $my_tag ;

    print OUT "The correct data loaded in: \n" ;
    print OUT "Name\t column:row:value \t\n";
    foreach $name (keys  %tag_store)
    {
	print OUT $name . "  "  ; 
	foreach $my_tag (@{$tag_store{$name}} )
	{
	    print OUT "  \t" . $my_tag ;
	}
	print OUT "\n" ; 
    }

}
