#
#  Script to write a list of directories that
#  contain source code. 
#
#  At the moment this is required to do
#  source level debugging with the riscwatch
#  debugger. 
#
#  To load the output file into the riscwatch
#  debugger, use:
#     exec  debug_path 
#
#

require "find.pl" ;

$dir = "src" ;
&update_dir($dir) ; 

$dir = "include" ;
&update_dir($dir) ; 


#
# write the directories to a file
# that can be included in the riscwatxh deugger.
#

$script_name = "debug_path" ;
open(DEBUG,">$script_name") || die("Could not open $script_name")  ; 
##print DEBUG "setenv SEARCH_PATH " ;
foreach $dir ( keys %dir_store ) 
{
print DEBUG "srchpath add " . $dir . "\n" ; 
}
print DEBUG "\n" ; 
close(DEBUG); 

exit(0) ; 

#
#  Recurse down the $dir (first argument) directory
#  and store the directory name if the file ends
#  in .C .h or .S
#

sub update_dir
{
    my ($dir) = @_ ; 
    local $dir_target ; 
    print "Seaching for source code in \"".$dir."\" directory\n" ; 

    my $source = $dir ;
    if( ! -d $source )
    {
        print "ERROR can not find " .$source. "\n"  ; 
	exit(1)  ;
    }
    &find($source) ; 

}




sub wanted
{
    if( !($dir =~ /CVS/ || $_ =~ /CVS/)  )
    {
	$here = `pwd` ; 
	chop ($here) ; 
	
	if( -f $_ && 
	    ( $_ =~ /.C$/  || $_ =~ /.h$/  || $_ =~ /.S$/  )   
	    )
	{
	    $dir_store{$here} = "1" ; 
	}

    }


}


# required for require
1 ; 
