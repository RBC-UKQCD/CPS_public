#include<config.h>
CPS_START_NAMESPACE
#!/usr/bin/perl/
# My very own perl script to bludgeon the qcdsp comments into a form
# that doxygen can use.
# Major issues:
#  - convert single line comments // into //!
#  - convert multi-line comments (/* */ OR // \n //) into /*! */
#  - move post-comments into pre-comments (urgh).
#  - Clip useless bits from comments, such as ------------lines.
#
# Also, modified NOT to filter codes that already contain doxygen comments.
#
# usage:
#  anjCommentFilter.pl <infile> > <outfile-via-stdout>
#
# Anj 2001.
#

# define code options:
use constant OUTPUT_MUNGY => 0; # if 1, output lots of stuff for debugging.

# parse CLargs:
if( scalar(@ARGV) != 1 ) {
  die "usage: anjCommentFilter.pl <infile> [> <outfile-via-stdout>]\n";
}
$infile = $ARGV[0];

# define the types into which the file will be split:
use constant EMPTY      => -2;
use constant CODE       => 0;
use constant COMMENT    => 1;
use constant NEWLINE    => 2;
use constant SL_COMM    => 3;
use constant ML_COMM    => 4;
use constant CODE_SPLIT => 5;
use constant DOX_COMM   => 10;
use constant C_START    => 100;
use constant C_END      => 101;

# define possible comment tags:
$cSt[0] = quotemeta('/*'); $cEn[0] = quotemeta('*/');
$cSt[1] = quotemeta('//'); $cEn[1] = quotemeta("\n");
$ictot = 2;
# build any-match string:
$cStart = '(';
$cEnd = '(';
$cAllButNl = '('; $alli = 0;
for($ci = 0; $ci < $ictot; $ci++ ) {
    if ($ci > 0 ) {
	$seperator = "|";
    } else {
	$seperator ="";
    }
    $cStart = $cStart.$seperator.$cSt[$ci];
    $cEnd = $cEnd.$seperator.$cEn[$ci];
    # all but newline:
    if ($alli > 0 ) {
	$seperator = "|";
    } else {
	$seperator ="";
    }
    if( $cSt[$ci] ne quotemeta("\n") ) {
	$cAllButNl = $cAllButNl.$seperator.$cSt[$ci];
	$alli++;
    }
    if ($alli > 0 ) {
	$seperator = "|";
    } else {
	$seperator ="";
    }
    if( $cEn[$ci] ne quotemeta("\n") ) {
	$cAllButNl = $cAllButNl.$seperator.$cEn[$ci];
	$alli++;
    }
}
$cStart = $cStart.')';
$cEnd = $cEnd.')';
$cAllButNl = $cAllButNl.')';
#print STDERR "cStart = $cStart\n";
#print STDERR "cEnd = $cEnd\n";
#print STDERR "cAllButNl = $cAllButNl\n";

# define single and multi-line output tags:
$sOutStart  = "//! "; $sOutEnd = "";
$mOutStart  = "/*! "; $mOutEnd = "*/";

# open dubug-info file is necc.
open FOUT,">munged.mungy.out" or die "Could not open dumpfile!\n" if(OUTPUT_MUNGY==1);

# check the file for valid doxygen comments:
$doxycomm = &check_for_doxycomm($infile);
# If any doxycomments were found, spew out the original file:
if( $doxycomm == 1 ) {
  open DCIN,"$infile" or die("Could not open $infile!\n");
  @original_file = <DCIN>;
  print @original_file;
  close DCIN;
  print STDERR "Doxygen comment identifier(s) found in $infile - no filtering performed.\n";
  exit;
}

# parse file, loading into an array of strings of the type detailed above:
&parse_source($infile);

# parser-debug output
if(OUTPUT_MUNGY==1) {
    for( $iel = 0; $iel < $ftot; $iel++ ) {
	print FOUT "M:l=$iel id=$ftype[$iel]:$fdat[$iel]";
    }
}

# munge the comments:
&munge_file();

# output the munged file:
for( $iel = 0; $iel < $ftot; $iel++ ) {
    print FOUT "l=$iel id=$ftype[$iel]:$fdat[$iel]" if(OUTPUT_MUNGY==1);
    print "$fdat[$iel]";
}

# close debug file:
close FOUT if(OUTPUT_MUNGY==1);

exit;

#----------------------------------------------------------------------
sub munge_file {

    #convert COMMENTs into single line and multi-line COMMENTs.
    $incomment = 0; $comstart = -2;
    for( $iel = 0; $iel < $ftot; $iel++ ) {
	#detect comment start:
	if( $ftype[$iel] == COMMENT ) {
	    $comstart = $iel if( $incomment == 0 );
	    $incomment++;
	}

	# munge comments:
	if ( $incomment > 0 && $ftype[$iel] != COMMENT ) {
	    # single line: add the right bits:
	    if( $iel == $comstart+1 ) {
		if( $fdat[$comstart]=~/[\n]$/ ) {
		    $fdat[$comstart]=~s/[\n]$//;
		    $endofstr = "\n";
		} else {
		    $endofstr = "";
		}
		$fdat[$comstart] = $sOutStart.$fdat[$comstart].$sOutEnd.$endofstr;
	    # multi-line: add the right bits:
	    } else {
		# cram all the comment lines into one string:
		for( $il = $comstart+1; $il < $iel; $il++ ) {
		    $fdat[$comstart] = $fdat[$comstart].$fdat[$il];
		    $fdat[$il] = ""; $ftype[$il] = EMPTY;
		}
		if( $fdat[$comstart]=~/[\n]$/ ) {
		    $fdat[$comstart]=~s/[\n]$//;
		    $endofstr = "\n";
		} else {
		    $endofstr = "";
		}
		$fdat[$comstart] = $mOutStart.$fdat[$comstart].$mOutEnd.$endofstr;
	    }
	    $comstart = -2;
	    $incomment = 0;
	}

    }

    # Look for post-comments to move:
    for( $iel = 0; $iel < $ftot; $iel++ ) {
	# grab the next three-line chunk, avoiding empty lines;
	$iblk = 0; $iline = $iel;
	while( $iblk < 3 && $iline < $ftot ) {
	    if( $ftype[$iline] != EMPTY ) {
		$blk[$iblk] = $iline;
		$iblk++;
	    }
	    $iline++;
	}
	# if NOT EOF and is a backasswards comment...
	if( $iblk == 3 && $ftype[$blk[0]] == CODE &&
	    $ftype[$blk[1]] == COMMENT &&
	    $ftype[$blk[2]] == NEWLINE ) {
	    #...then point comment backwards:
	    $fdat[$blk[1]]=~s/![^!]/!</;
	} 
	#point comment backwards for one liners as well:
	if( $iblk >= 2 && $ftype[$blk[0]] == CODE &&
	    $ftype[$blk[1]] == COMMENT && !($ftype[$blk[0]]=~/$[\n]/) ) {
	    $fdat[$blk[1]]=~s/![^!]/!</;
	}

    }
}

#-------------------------------------
sub remove_comments_from {
    my ($instr) = @_;

    $outstr=$instr;
    $outstr =~s/($cAllButNl)//g;

    return $outstr;
}

#----------------------------------------------------------------------
sub parse_source {
    my ( $fname ) = @_;

    open FIN,"$fname" or die "Could not open $fname!\n";

    $ini = 0;
    $incomment = 0;
    $comtype = -1;
    $foundtype = C_START;
    while (<FIN>) {
        # to deal with combo lines:
	if( $_=~/(\S{1,}\s*)$cStart/ ) {
	    @lbits = split "$cStart|$cEnd",$_;
	} else {
	    @lbits = ();
	    $lbits[0] = $_;
	    #@lbits = split "\n",$_;
	}
	
	# loop over sections:
	for($ilb = 0; $ilb < scalar(@lbits); $ilb++) {
	    next if( $lbits[$ilb] eq "" );
	    # look for any of the comment specifiers in the current line:
	    for($ic = 0; $ic < $ictot; $ic++ ) {
		if( $lbits[$ilb]=~/$cSt[$ic]/ ) {
		    $foundtype = C_START;
		    $incomment = 1;
		    $comtype = $ic;
		}
	    }
	    if( $comtype != -1 && $lbits[$ilb]=~/$cEn[$comtype]/ ) {
		$foundtype = C_END;
		$comtype = -1;
	    }

	    # if comment spec. the either (DOX)COMMENT, continuing COMMENT
	    #  OR a combination of things on one line:
	    if( $incomment == 1 ) {
		$fdat[$ini] = $lbits[$ilb];
		$ftype[$ini] = COMMENT;
		# Strip out unwanted things:
		$fdat[$ini] = &remove_comments_from($fdat[$ini]);
		 # spaces:
		 $fdat[$ini]=~s/ {3,}//g;
		 # lines like -------...
		 $fdat[$ini]=~s/-{3,}//g;
		 # lines like *******...
		 $fdat[$ini]=~s/\*{3,}//g;
		 # lines like =======...
		 $fdat[$ini]=~s/={3,}//g;
		$ini++;
	  # if no comment spec. then either NEWLINE, CODE or DOX COMMENT:
	    } else {
		if( $lbits[$ilb]=~/^\s*[\n]/ ) {
		    $ftype[$ini] = NEWLINE;
		    $fdat[$ini] = $lbits[$ilb];
		    $ini++;
		} else {
		    $ftype[$ini] = CODE;
		    $fdat[$ini] = $lbits[$ilb];
		    $ini++;
		}
	    }
	    # End-Of-Comment
	    $incomment = 0 if( $foundtype == C_END );
	}
    }
    
    close FIN;

    # set total elements:
    $ftot = scalar(@fdat);

}

# $doxycomm = &check_for_doxycomm($infile);
# Check a given file for the presence of doxygen comments:
#----------------------------------------------------------------------
sub check_for_doxycomm {
 my ($dcinfile) = @_;
 
 open DCIN,"$dcinfile" or die("Could not open $dcinfile!\n");
 
 $foundoxycomm = 0;
 while ( <DCIN> ) { 
   if( $_=~/\/(\/|\*)<?\!/ ) {
     $foundoxycomm = 1;
     #print STDERR "Found: $dcinfile\n";
     last;
   }
 }

 close DCIN;

 return $foundoxycomm;
}

CPS_END_NAMESPACE
