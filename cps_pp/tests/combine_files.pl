#!/usr/bin/perl
#

if( $#ARGV  < 0 )
{
 die( "Need files to combine \n");
}


foreach $dat (@ARGV)
{
print "--------------------\n" ;
	print $dat."\n";
	print "--------------------\n" ;

	open(FFF,"$dat") || die("Could not open $dat") ; 
        while($line = <FFF>) 
{
print $line ; 
}

	close(FFF) ; 

}

exit 0 ; 

