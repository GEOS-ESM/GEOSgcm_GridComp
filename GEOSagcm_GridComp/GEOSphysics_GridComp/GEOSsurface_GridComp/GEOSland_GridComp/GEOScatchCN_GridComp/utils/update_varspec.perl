#!/usr/bin/perl

# ------------------------------------------------------------
# variable spec attributes
# ------------------------------------------------------------

$vlocation='GEOS_VLocationNone';
$averaging_interval='ACCUMINT';
$refresh_interval='MY_STEP';
$linelength=65;

$divider='! '.'-'x($linelength-2)."\n";

$insert1='INSERT_VARSPEC_CODE';
$insert2='INSERT_POINTER_DECLARATIONS_CODE';
$insert3='INSERT_RETRIEVE_POINTERS_CODE';

# ------------------------------------------------------------
# Create the varspecs
# ------------------------------------------------------------

open IMPORTFILE, "< import.txt";
@importlist=<IMPORTFILE>;
$str="@importlist";
close IMPORTFILE;
#print $str;

open INTERNALTFILE, "< internal.txt";
@internallist=<INTERNALTFILE>;
$str="@internallist";
close INTERNALTFILE;
#print $str;

open EXPORTFILE, "< export.txt";
@exportlist=<EXPORTFILE>;
$str="@exportlist";
close EXPORTFILE;
#print $str;

@varspec_code=();
@declarations_code=();
@retrieval_code=();

push @varspec_code,  "\n";
push @varspec_code,  "$divider";
push @varspec_code,  "! IMPORT VARSPEC\n";
push @varspec_code,  "$divider";
push @varspec_code,  "\n";

push @declarations_code,  "\n";
push @declarations_code,  "$divider";
push @declarations_code,  "! IMPORT Pointers\n";
push @declarations_code,  "$divider";
push @declarations_code,  "\n";

push @retrieval_code,  "\n";
push @retrieval_code,  "$divider";
push @retrieval_code,  "! IMPORT Pointers\n";
push @retrieval_code,  "$divider";
push @retrieval_code,  "\n";

foreach $line (@importlist) {

    $line=~s/(.*)\n/$1/;
    @parts=split /\b\s+\b/, $line;
    $shortname=$parts[0];
    $lcshortname=lc $shortname;
    $longname=$parts[1];
    $units=$parts[2];
    $dimensiontag=$parts[3];
    $dimensiontag=~s/\s*$//;

    push @varspec_code,  padline('call GEOS_StateAddSpec(STATE', $linelength);
    push @varspec_code,  padline("  LONG_NAME          = '$longname'", $linelength);
    push @varspec_code,  padline("  UNITS              = '$units'",$linelength);
    push @varspec_code,  padline("  SHORT_NAME         = '$shortname'",$linelength);
    push @varspec_code,  padline("  DIMS               = GEOS_Dims$dimensiontag", $linelength);	
    push @varspec_code,  padline("  VLOCATION          = $vlocation", $linelength);
    push @varspec_code,  padline("  AVERAGING_INTERVAL = $averaging_interval", $linelength);
    push @varspec_code,  padline("  REFRESH_INTERVAL   = $refresh_interval", $linelength);
    push @varspec_code,  padline("  IMPORT             = .true.", $linelength);
    push @varspec_code,  '  RC=STATUS  ) '."\n";
    push @varspec_code,  "\n";
    push @varspec_code,  'VERIFY_(STATUS)'."\n";
    push @varspec_code,  "\n";

    push @retrieval_code, padline("call GET_POINTER(IMPORT",$linelength);
    push @retrieval_code, padline("                 $lcshortname",$linelength);
    push @retrieval_code, padline("                 '$shortname'",$linelength);
    push @retrieval_code,         '                 RC=STATUS  ) '."\n";
    push @retrieval_code, "\n";
    push @retrieval_code, 'VERIFY_(STATUS)'."\n";
    push @retrieval_code, "\n";

    if ($dimensiontag eq 'TileOnly') {
	push @declarations_code, '  real, dimension(:),   pointer :: '.$lcshortname."\n";
    } else {
	push @declarations_code, '  real, dimension(:,:), pointer :: '.$lcshortname."\n";
    }

}

push @varspec_code,  "\n";
push @varspec_code,  "$divider";
push @varspec_code,  "! INTERNAL VARSPEC\n";
push @varspec_code,  "$divider";
push @varspec_code,  "\n";

push @declarations_code,  "\n";
push @declarations_code,  "$divider";
push @declarations_code,  "! INTERNAL Pointers\n";
push @declarations_code,  "$divider";
push @declarations_code,  "\n";

push @retrieval_code,  "\n";
push @retrieval_code,  "$divider";
push @retrieval_code,  "! INTERNAL Pointers\n";
push @retrieval_code,  "$divider";
push @retrieval_code,  "\n";

foreach $line (@internallist) {

    $line=~s/(.*)\n/$1/;
    @parts=split /\b\s+\b/, $line;
    $shortname=$parts[0];
    $lcshortname=lc $shortname;
    $longname=$parts[1];
    $units=$parts[2];
    $dimensiontag=$parts[3];
    $dimensiontag=~s/\s*$//;

    push @varspec_code,  padline('call GEOS_StateAddSpec(STATE', $linelength);
    push @varspec_code,  padline("  LONG_NAME          = '$longname'", $linelength);
    push @varspec_code,  padline("  UNITS              = '$units'",$linelength);
    push @varspec_code,  padline("  SHORT_NAME         = '$shortname'",$linelength);
    push @varspec_code,  padline("  DIMS               = GEOS_Dims$dimensiontag", $linelength);	
    push @varspec_code,  padline("  VLOCATION          = $vlocation", $linelength);
    push @varspec_code,  padline("  INTERNAL           = .true.", $linelength);
    push @varspec_code,  '  RC=STATUS  ) '."\n";
    push @varspec_code,  "\n";
    push @varspec_code,  'VERIFY_(STATUS)'."\n";
    push @varspec_code,  "\n";

    push @retrieval_code, padline("call GET_POINTER(INTERNAL",$linelength);
    push @retrieval_code, padline("                 $lcshortname",$linelength);
    push @retrieval_code, padline("                 '$shortname'",$linelength);
    push @retrieval_code,         '                 RC=STATUS  ) '."\n";
    push @retrieval_code, "\n";
    push @retrieval_code, 'VERIFY_(STATUS)'."\n";
    push @retrieval_code, "\n";
    
    if ($dimensiontag eq 'TileOnly') {
	push @declarations_code, '  real, dimension(:),   pointer :: '.$lcshortname."\n";
    } else {
	push @declarations_code, '  real, dimension(:,:), pointer :: '.$lcshortname."\n";
    }

}

push @varspec_code,  "\n";
push @varspec_code,  "$divider";
push @varspec_code,  "! EXPORT VARSPEC\n";
push @varspec_code,  "$divider";
push @varspec_code,  "\n";

push @declarations_code,  "\n";
push @declarations_code,  "$divider";
push @declarations_code,  "! EXPORT Pointers\n";
push @declarations_code,  "$divider";
push @declarations_code,  "\n";

push @retrieval_code,  "\n";
push @retrieval_code,  "$divider";
push @retrieval_code,  "! EXPORT Pointers\n";
push @retrieval_code,  "$divider";
push @retrieval_code,  "\n";

foreach $line (@exportlist) {

    $line=~s/(.*)\n/$1/;
    @parts=split /\b\s+\b/, $line;
    $shortname=$parts[0];
    $lcshortname=lc $shortname;
    $longname=$parts[1];
    $units=$parts[2];
    $dimensiontag=$parts[3];
    $dimensiontag=~s/\s*$//;

    push @varspec_code,  padline('call GEOS_StateAddSpec(STATE', $linelength);
    push @varspec_code,  padline("  LONG_NAME          = '$longname'", $linelength);
    push @varspec_code,  padline("  UNITS              = '$units'",$linelength);
    push @varspec_code,  padline("  SHORT_NAME         = '$shortname'",$linelength);
    push @varspec_code,  padline("  DIMS               = GEOS_Dims$dimensiontag", $linelength);	
    push @varspec_code,  padline("  VLOCATION          = $vlocation", $linelength);
    push @varspec_code,  padline("  EXPORT             = .true.", $linelength);
    push @varspec_code,  '  RC=STATUS  ) '."\n";
    push @varspec_code,  "\n";
    push @varspec_code,  'VERIFY_(STATUS)'."\n";
    push @varspec_code,  "\n";

    push @retrieval_code, padline("call GET_POINTER(EXPORT",$linelength);
    push @retrieval_code, padline("                 $lcshortname",$linelength);
    push @retrieval_code, padline("                 '$shortname'",$linelength);
    push @retrieval_code, padline("                  alloc=.true.",$linelength);
    push @retrieval_code,         '                 RC=STATUS  ) '."\n";
    push @retrieval_code, "\n";
    push @retrieval_code, 'VERIFY_(STATUS)'."\n";
    push @retrieval_code, "\n";

    if ($dimensiontag eq 'TileOnly') {
	push @declarations_code, '  real, dimension(:),   pointer :: '.$lcshortname."\n";
    } else {
	push @declarations_code, '  real, dimension(:,:), pointer :: '.$lcshortname."\n";
    }
}

open GRIDCOMP, "< GEOS_CatchGridComp.F90";
@catchlist=<GRIDCOMP>;
$oldfile="@catchlist";
close GRIDCOMP;

$oldfile=insert_here($insert1, $oldfile, "@varspec_code");
$oldfile=insert_here($insert2, $oldfile, "@declarations_code");
$oldfile=insert_here($insert3, $oldfile, "@retrieval_code");

#print "\n$newstring";

open GRIDCOMP, "> GEOS_CatchGridComp.F90";
print GRIDCOMP "$oldfile";
close GRIDCOMP;

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


