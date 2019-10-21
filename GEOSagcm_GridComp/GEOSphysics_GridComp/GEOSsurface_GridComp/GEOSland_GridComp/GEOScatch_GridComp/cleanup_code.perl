#!/usr/bin/perl

$"='';
$linelength=80;
$tab=4;
$editfile=$ARGV[0];


@spacebefore=('call ', '\=\s*GEOS\_', '\=\s*ESMF\_');
@spaceboth=('end ','^\s*where','^\s*contains', '^\s*private', '^\s*public', '^\s*subroutine', '^\s*if.*then', '^\s*do');
@spaceafter=('VERIFY_');

# --------------------------------------------------------------
# 1. Consolidation:  this part of the program compresses all
#    contents of code into an array with no extra spaces
# --------------------------------------------------------------

open OLDFILE, "< $editfile";
@linelist=<OLDFILE>;
close OLDFILE;

$bstr="0";
$backup=$editfile.'_safechange_'.$bstr;
while (-e $backup) {
    $backup=$editfile.'_safechange_'.$bstr;
    $bstr++;
}

open BACKFILE, "> $backup";

$was=$";
$"='';
$bfile="@linelist";
print BACKFILE $bfile;
close BACKFILE;
$"=$was;

# collapse multiple line commands

@newlist=();
$lineindex=0;

foreach $ll (@linelist) {
    
    # stip off leading blanks before
    # any non-whitespace character

#    $ll=~s/^\s*(\S+)/$1/;

    $ll=~s/^\s*//;
    $thisline.=$ll;

    # Don't touch comments
    if ($ll!~m/^\s*\!/) {

	# Shorten any spaces of two or more
	$thisline=~s/  +/ /g;

	# take out any spaces around nonword chars
	$thisline=~s/[ ]*(\W)[ ]*/$1/g;

	# Remove any tabs 
	$thisline=~s/\t*//g;

	# if & present, remove it and \n 
	$thisline=~s/\s*&\s*\n//;

	if ($ll=~m/\&\W*\n/) {
	    next;
	}
    }

    # remove spaces before newline
    $thisline=~s/\s*\n/\n/;

    if ($thisline ne '') {
	push @newlist, $thisline;
    }
    $thisline='';
}
#print "@newlist";

# --------------------------------------------------------------
# 2. Indentation:  The lines are padded with the correct 
#    indentation that adds the correct tabs as specified
#    Comments stick out to the left by one tab space
# --------------------------------------------------------------

@formatlist=();
$indent=0;
$currentlength=$linelength;

foreach $ll (@newlist) {

    if ($ll=~m/end\s*subroutine/ || $ll=~m/end\s*do/ || $ll=~m/end\s*if/) {
	$indent-=$tab;
	$currentlength+=$tab;
	$status='-o-';
    } else { 
	$status='';
    }
    
    if ($ll=~m/^\!/) {
	$cindent=$indent-$tab;
	$newstring=' 'x$cindent.$ll;
    } else {
	$newstring=' 'x$indent.$ll;
    }
    
    push @formatlist, $newstring;

    if ($ll=~m/^(subroutine|do|if)/ )   {
	$indent+=$tab;
	$currentlength-=$tab;
    }    
}

# --------------------------------------------------------------
# 3. Spacing: Add carriage returns for greater readability
# --------------------------------------------------------------

@spacelist=();

$"='|';
$both="@spaceboth";
$after="@spaceafter";
$before="@spacebefore";

$"='';
$prev='';
$cur="0";
foreach $ll (@formatlist) {
  
    if (($ll=~m/^\s*\!/) && ($prev!~m/^\s*\!/)) {
	if ($cur eq "0") {
	    push @spacelist, "\n";
	}
	push @spacelist, "$ll";
	$cur="0";
	$prev=$ll;
	next;
     } elsif (($ll!~m/^\s*\!/) && ($prev=~m/^\s*\!/)) {
        if ($cur eq "0") {
            push @spacelist, "\n";
        }
         push @spacelist, "$ll";
        $cur="0";
	$prev=$ll;
	next;
    }

    if ($ll=~m/$both/) {
	if ($cur eq "0") {
	    push @spacelist, "\n";
	}
	push @spacelist, "$ll";
	push @spacelist, "\n";
	$cur="1";
    } elsif ($ll=~m/$after/) {
	push @spacelist, "$ll";
	push @spacelist, "\n";
	$cur="1";
    } elsif ($ll=~m/$before/) {
	if ($cur eq "0") {
	    push @spacelist, "\n";
	}
	push @spacelist, "$ll";
	$cur="0";
    } else {
       push @spacelist, "$ll";
       $cur="0";
    }
    $prev=$ll;
}

# --------------------------------------------------------------
# 4. breaking down long lines
# two cases dealt with: math and call statements
# --------------------------------------------------------------

@squishlist=();
foreach $ll (@spacelist) {
    if ($ll=~m/^\s*call\s+\b\w+\b\s*\(|^\s*subroutine\s+\b\w+\b|^\s*\S+\=(GEOS|ESMF)\_\S+\(/) {
	

#    if ($ll=~m/\((\s*\W+\s*\,)+/) {
	@p1=split /\(/, $ll;
	$head=shift @p1;
       
	# strip off last space in header
	$head=~s/\s*$//;

	$back=join '(', @p1;
	$back=~s/\s*\)\s*\n//;
	
	@p2=split /\,/, $back;
	
	# strip off any extra space in parameter list
	@params=();
	$openp=0;
	$closep=0;
	$cur='';
	$count=0;
	@newparams=();
	foreach $param (@p2) {
	    $param=~s/^\s*(\S+)\s*$/$1/g;
	    if ($count > 1) {
		$cur.=', '.$param;
	    } else {
		$cur.=$param;
	    }

	    $count++;
	    $openp=($cur=~m/\(/);
	    $closep=$cur=~m/\)/;
	    $same=($openp && $closep) || ((!$openp) && (!$closep)); 
	    
	    if ($same) { 
		push @newparams, $cur;	
		$cur='';
		$count=0;
	    }
	}
	@params=@newparams;

	$"=', ';
	$plist="@params";
	$"='';
	$oneline="$head ".'( '."$plist".' )'."\n";

	if (length($oneline) <= $linelength) {
	    
	    push @squishlist, $oneline;

	} else {

	    # find out how much of an indent there was
	    
	    @newinsert=();

	    $chopped=$head;
	    $chopped=~s/^\s*(\S+)/$1/;
	    $bsp=length($head)-length($chopped);
	    $pad=' 'x$bsp;
	    $padtab=' 'x$bsp.' 'x$tab;
		
	    $bigger=0;
	    $amtleft=$linelength-length($head)-5;
	    foreach $v (@params) {
		$try=length($v);
		if ($try >= $amtleft) {
		    $bigger=1;
		    last;
		}
	    }
	  
	    $finalp=pop(@params);

	    if ($bigger) {

		$sz=$linelength-length($head)-4;
		$sp=' 'x$sz;
		
		$keepline="$head".' ( '."$sp".'&'."\n";
		push @squishlist, $keepline;
		
		foreach $v (@params) {
		    $sz=$linelength-length($v)-length($padtab)-2;
		    $sp=' 'x$sz;
		    $keepline="$padtab".$v.$sp.',&'."\n"; 
		    push @squishlist, $keepline;
		}
		$keepline="$padtab".$finalp.' '.')'."\n"; 
		push @squishlist, $keepline;

	    } else {
		
		$v=@params[0];
		
		$sz=$linelength-length($head)-length($v)-5;
		$sp=' 'x$sz;
		$keepline="$head".' ( '."$params[0]"."$sp".',&'."\n";

		push @squishlist, $keepline;

		shift(@params);

		foreach $v (@params) {
		    $sz=$linelength-length($head)-length($v)-5;
		    $sp=' 'x$sz;
		    $szz=length($head)+3;
		    $spp=' 'x$szz;
		    $keepline="$spp".$v.$sp.',&'."\n"; 
		    push @squishlist, $keepline;
		}
		
		$keepline="$spp".$finalp.' '.')'."\n"; 
		push @squishlist, $keepline;
		
	    }
	}
	
    } else {
	push @squishlist, $ll;
    }
}
$"='';

@dashlist=();

foreach $ll (@squishlist) {
    if ($ll=~m/^\s*\!\s*\-+\s*\n/) {

	$chopped=$ll;
	$chopped=~s/^\s*//;
	$amtleft=length($ll)-length($chopped);
	$indent=' 'x$amtleft;
	$numdashes=$linelength-$amtleft-2;
	$dashes='-'x$numdashes;
	$newline=$indent.'! '.$dashes."\n";
	push @dashlist, $newline;
    } else {
	push @dashlist, $ll;
    }
}

open NEWFILE, "> $editfile";
print NEWFILE "@dashlist";
close NEWFILE;

exit;

sub padline {
    $myscalar=$_[0];
    $len=$_[1];
    $siz=length $myscalar;
    $leftover=$len-$siz-2; 
    $addon=' 'x$leftover;
    $newline="$myscalar".$addon.',&'."\n";
    return $newline;
}

sub insert_here {

    $tag=$_[0];
    $oldfile=$_[1];
    $newtext="! BEGIN_$tag\n"."$_[2]\n"."! END_$tag\n";
    $oldfile=~s/\!\W*?BEGIN\_$tag\W*?\n.*\!\W*?END\_$tag\W*?\n/$newtext/s;
    return $oldfile;

}


