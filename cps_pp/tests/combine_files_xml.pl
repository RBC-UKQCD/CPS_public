#!/usr/bin/perl
#

if( $#ARGV  < 0 )
{
 die( "Need files to combine \n");
}


#
#  xml header information
#
print "<?xml version=\"1.0\"?>\n" ; 

print "<CpsDataStore>\n" ;

#
#  loop over each file
#


foreach $dat (@ARGV)
{
print "<".$dat.">\n";

$is_ok = 1 ; 

$count_row = 0  ; 
open(FFF,"$dat") || die("Could not open $dat") ; 
while($line = <FFF>) 
{
#    print $line ; 
    @line_contents = split(' ',$line) ; 
    $len = $#line_contents ;

    if( $count_row == 0 ) 
    {
	for($ii = 0 ; $ii <= $len ; ++$ii)
	{
	    $store[$ii]  = $line_contents[$ii] ;
	}
	$len_store = $len ; 
    }
    else
    {
	for($ii = 0 ; $ii <= $len ; ++$ii)
	{
            $tmp = $store[$ii] ;
	    $store[$ii]  = $tmp. "   " . $line_contents[$ii] ;
	}

    }
	if( $len_store != $len  )
	{
	    $is_ok = 0 ; 
	}


    ++$count_row ;

}


print "<fileContent>\n" ;
if( $is_ok == 1  )
{
    for($jj = 0 ; $jj <= $len ; ++$jj)
    {
	$elem = "Column" . $jj ; 
	if( $jj == $len ) 
	{
	    $elem = "last" ; 
	}
	
	print "<".$elem.">\n" ; 
	print $store[$jj] . "\n" ; 
	print "</".$elem.">\n" ; 
    }
}
else
{
    print "<combine_files_xml_ERROR>" ; 
    print "All row must have same number of columns\n" ;
    print "</combine_files_xml_ERROR>" ; 
}




print "</fileContent>\n" ;
print "</".$dat.">\n";

close(FFF) ; 

}


print "</CpsDataStore>\n" ;

exit 0 ; 

