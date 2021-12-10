#!/usr/bin/env perl
#
#  This script assumes that data in the .raw file is in a v33 format.
#
#  Read in command line arguments. This script requires at least 2 arguments and
#  optionally several more. The two required arguments are the name of the PSS/E
#  .RAW file that is used to generate the contingencies and cost parameter files
#  and the second argument is the root name that will be used for the output
#  files.
#
#  The optional directives are
#
#  -GC  (Single generator contingencies)
#  -LC  (Single line contingencies)
#  -GGC (Two generator contingencies)
#  -GLC (A generator and a line contingecies)
#  -LLC (Two line contingencies)
#
#  Each of these directives can be used to control the number of contingencies
#  of a given type that are added to the contingency file.
#
#  The allowable arguments for these directives are
#
#  all
#  none
#  some number
#
#  If "all is specified, then all possible contingencies of that type will be
#  generated, if "none" is specified then no contingencies of that type will be
#  generated and if a number is specified, then that number of contingencies of
#  that type is generated. By default, the -GC and -LC options are set to "all"
#  and the -GGC, -GLC and -LLC options are set to none
#
$argc    = @ARGV;
$rawfile = "";
if ( $argc > 0 ) {
  $rawfile = $ARGV[0];
}
else {
  print("No Mat Power file specified. Exiting...\n");
  exit(0);
}
print "Mat Power FILE NAME: $rawfile\n";
$rootname = "";
if ( $argc > 1 ) {
  $rootname = $ARGV[1];
}
print "EXPORT ROOT: $rootname\n";
#
#  Parse arguments to find out how many contingencies to create. If no arguments
#  create all N-1 generator and line contingencies and no N-2 contingencies.
#
$allGC  = 1;
$allLC  = 1;
$allGGC = 0;
$allGLC = 0;
$allLLC = 0;
$maxGC  = 0;
$maxLC  = 0;
$maxGGC = 0;
$maxGLC = 0;
$maxLLC = 0;

for ( $i = 0 ; $i < $argc ; $i++ ) {
  if ( $ARGV[$i] eq "-GC" && defined( $ARGV[ $i + 1 ] ) ) {
    $arg = lc( $ARGV[ $i + 1 ] );
    if ( $arg eq "all" ) {
      $allGC = 1;
    }
    elsif ( $arg eq "none" ) {
      $allGC = 0;
    }
    else {
      $allGC = 0;
      $maxGC = $arg;
    }
  }
  elsif ( $ARGV[$i] eq "-LC" && defined( $ARGV[ $i + 1 ] ) ) {
    $arg = lc( $ARGV[ $i + 1 ] );
    if ( $arg eq "all" ) {
      $allLC = 1;
    }
    elsif ( $arg eq "none" ) {
      $allLC = 0;
    }
    else {
      $allLC = 0;
      $maxLC = $arg;
    }
  }
  elsif ( $ARGV[$i] eq "-GGC" && defined( $ARGV[ $i + 1 ] ) ) {
    $arg = lc( $ARGV[ $i + 1 ] );
    if ( $arg eq "all" ) {
      $allGGC = 1;
    }
    elsif ( $arg eq "none" ) {
      $allGGC = 0;
    }
    else {
      $allGGC = 0;
      $maxGGC = $arg;
    }
  }
  elsif ( $ARGV[$i] eq "-GLC" && defined( $ARGV[ $i + 1 ] ) ) {
    $arg = lc( $ARGV[ $i + 1 ] );
    if ( $arg eq "all" ) {
      $allGLC = 1;
    }
    elsif ( $arg eq "none" ) {
      $allGLC = 0;
    }
    else {
      $allGLC = 0;
      $maxGLC = $arg;
    }
  }
  elsif ( $ARGV[$i] eq "-LLC" && defined( $ARGV[ $i + 1 ] ) ) {
    $arg = lc( $ARGV[ $i + 1 ] );
    if ( $arg eq "all" ) {
      $allLLC = 1;
    }
    elsif ( $arg eq "none" ) {
      $allLLC = 0;
    }
    else {
      $allLLC = 0;
      $maxLLC = $arg;
    }
  }
}
#
#  Open file to write out cost data
#
$count = 1;
#
#  Open file for contingencies
#
$filename = "$rootname\.con";
#
#  Parse Mat Power format file for information on generators and lines
#
open( MP,          "$rawfile" );
open( CONTINGENCY, ">$filename" );
$is_gen           = 0;
$is_branch        = 0;
$is_bus           = 0;
$is_gen_parsed    = 0;
$is_branch_parsed = 0;
$is_bus_parsed    = 0;
@nminus2          = ();
$count2           = 0;
$gcnt             = 0;
$lcnt             = 0;
$count_gc         = 0;
$count_lc         = 0;
$buscnt           = 0;
%genbus           = {};
@fromBus          = ();
@toBus            = ();
@lineTags         = ();
%nlines           = {};
%linecnt          = {};

while (<MP>) {
  $line = $_;
  if ( $line =~ /\]\;/ && $is_gen == 1 ) {
    $is_gen        = 0;
    $is_gen_parsed = 1;
    print "Closing generators: $line";
  }
  if ( $line =~ /\]\;/ && $is_branch == 1 ) {
    $is_branch        = 0;
    $is_branch_parsed = 1;
    print "Closing lines: $line";
  }
  if ( $line =~ /\]\;/ && $is_bus == 1 ) {
    $is_bus        = 0;
    $is_bus_parsed = 1;
    print "Closing buses: $line";
  }
  if ( $is_gen == 1 ) {
    @tokens = split( ' ', $line );
    $busid  = $tokens[0];
    $busid =~ s/^\s*//;
    $busid =~ s/\s*$//;
    $active = $tokens[7];
    $active =~ s/^\s*//;
    $active =~ s/\s*$//;
    if ( $active > 0 ) {
      if ( defined( $genbus{$busid} ) ) {
        $genbus{$busid}++;
        $tag = $genbus{$busid};
      }
      else {
        $genbus{$busid} = 1;
        $tag = 1;
      }
      if ( $allGC == 1 || ( $count_gc < $maxGC ) ) {
        print CONTINGENCY "CONTINGENCY G\_$busid\_$tag\_$count\n";
        print CONTINGENCY "REMOVE UNIT $tag FROM BUS $busid\n";
        print CONTINGENCY "END\n";
        $contingency = "GENERATOR\_$busid\_$tag";
        $nminus2[$count2] = $contingency;
        $count++;
        $count_gc++;
        $count2++;
      }
      $gcnt++;
    }
  }
 #
 # Just gather line information for this pass. Construct contingencies on second
 # pass
 #
  if ( $is_branch == 1 ) {
    @tokens = split( ' ', $line );
    $bus1   = $tokens[0];
    $bus1 =~ s/^\s*//;
    $bus1 =~ s/\s*$//;
    $bus2 = $tokens[1];
    $bus2 =~ s/^\s*//;
    $bus2 =~ s/\s*$//;
    $active = $tokens[10];
    $active =~ s/^\s*//;
    $active =~ s/\s*$//;

    if ( $active > 0 ) {
      if ( defined( $linecnt{$bus1} ) ) {
        $linecnt{$bus1}++;
      }
      else {
        $linecnt{$bus1} = 1;
      }
      if ( defined( $linecnt{$bus2} ) ) {
        $linecnt{$bus2}++;
      }
      else {
        $linecnt{$bus2} = 1;
      }
      $qline = "$bus1\_$bus2";
      if ( defined( $nlines{$qline} ) ) {
        $nlines{$qline}++;
        $tag = $nlines{$qline};
      }
      else {
        $nlines{$qline} = 1;
        $tag = 1;
      }
      $fromBus[$lcnt]  = $bus1;
      $toBus[$lcnt]    = $bus2;
      $lineTags[$lcnt] = $tag;
      $lcnt++;
    }
  }
  if ( $is_bus == 1 ) {
    $buscnt++;
  }
  if ( $line =~ /mpc\.gen/ && $is_gen_parsed == 0 ) {
    $is_gen = 1;
    print "Starting line: $line";
  }
  if ( $line =~ /mpc\.branch/ && $is_branch_parsed == 0 ) {
    $is_branch = 1;
    print "Starting line: $line";
  }
  if ( $line =~ /mpc\.bus/ && $is_bus_parsed == 0 ) {
    $is_bus = 1;
    print "Starting bus: $line";
  }
}
print "Buses: $buscnt Generators: $gcnt Active Lines: $lcnt\n";
#
#  create linked list of neighbors
#
@top      = ();
@next     = ();
@nghbrs   = ();
@nghbrtag = ();
for ( $ibus = 0 ; $ibus < $buscnt ; $ibus++ ) {
  $top[$ibus] = -1;
}
$ncnt = 0;
for ( $iline = 0 ; $iline < $lcnt ; $iline++ ) {
  $fbus             = $fromBus[$iline];
  $tbus             = $toBus[$iline];
  $tmp              = $top[ $fbus - 1 ];
  $top[ $fbus - 1 ] = $ncnt;
  $next[$ncnt]      = $tmp;
  $nghbrs[$ncnt]    = $tbus - 1;
  $nghbrtag[$ncnt]  = $lineTags[ $tbus - 1 ];
  $ncnt++;
  $tmp              = $top[ $tbus - 1 ];
  $top[ $tbus - 1 ] = $ncnt;
  $next[$ncnt]      = $tmp;
  $nghbrs[$ncnt]    = $fbus - 1;
  $nghbrtag[$ncnt]  = $lineTags[ $fbus - 1 ];
  $ncnt++;
}
print "Total neighbor count: $ncnt\n";
#
#  Set up line contingencies
#
for ( $iline = 0 ; $iline < $lcnt ; $iline++ ) {
  $bus1   = $fromBus[$iline];
  $bus2   = $toBus[$iline];
  $tag    = $lineTags[$iline];
  $status = 1;
  if ( !defined( $linecnt{$bus1} ) || $linecnt{$bus1} < 2 ) {
    $status = 0;
  }
  if ( !defined( $linecnt{$bus2} ) || $linecnt{$bus2} < 2 ) {
    $status = 0;
  }
  if ( $status == 1 ) {
    if ( $allGC == 1 || ( $count_lc < $maxLC ) ) {
 #
 #  Check to see if contingency splits network into two pieces. Loop over
 #  neighbors recursively until you run out of new neighbors. If total number of
 #  buses in cluster is less than the total number of buses in the original
 #  network, go to next line outage
 #
      $bcnt   = 0;
      @newbus = ();
      @found  = ();
      for ( $ibus = 0 ; $ibus < $buscnt ; $ibus++ ) {
        $found[$ibus] = 0;
      }
      #
      #  Find first bus
      #
      $ifirst = 0;
      while ( $ifirst + 1 == $bus1 || $ifirst + 1 == $bus2 ) {
        $ifirst++;
      }
      $found[$ifirst] = 1;
      $bcnt           = 1;
      $nbus           = 0;
      $ibus           = $top[$ifirst];
      while ( $ibus >= 0 ) {
        $j             = $nghbrs[$ibus];
        $newbus[$nbus] = $j;
        $found[$j]     = 1;
        $bcnt++;
        $nbus++;
        $ibus = $next[$ibus];
      }
      while ( $nbus > 0 ) {
        @oldbus = ();
        for ( $i = 0 ; $i < $nbus ; $i++ ) { $oldbus[$i] = $newbus[$i]; }
        $newcnt = 0;
        @newbus = ();
        for ( $i = 0 ; $i < $nbus ; $i++ ) {
          $j = $top[ $oldbus[$i] ];
          while ( $j >= 0 ) {
            if (
              !(
                (
                     $oldbus[$i] + 1 == $bus1
                  && $nghbrs[$j] + 1 == $bus2
                  && $nghbrtag[$j] == $tag
                )
                || ( $oldbus[$i] + 1 == $bus2
                  && $nghbrs[$j] + 1 == $bus1
                  && $nghbrtag[$j] == $tag )
              )
              )
            {
              $ob = $oldbus[$i] + 1;
              $nb = $nghbrs[$j] + 1;
              if ( $found[ $nghbrs[$j] ] == 0 ) {
                $bcnt++;
                $found[ $nghbrs[$j] ] = 1;
                $newbus[$newcnt] = $nghbrs[$j];
                $newcnt++;
              }
            }
            $j = $next[$j];
          }
        }
        $nbus = $newcnt;
      }
      if ( $bcnt == $buscnt ) {
        print CONTINGENCY "CONTINGENCY L\_$bus1\-$bus2\_$tag\_$count\n";
        print CONTINGENCY "OPEN BRANCH FROM BUS $bus1 TO $bus2 CIRCUT $tag\n";
        print CONTINGENCY "END\n";
        $contingency = "LINE\_$bus1\_$bus2\_$tag";
        $nminus2[$count2] = $contingency;
        $count_lc++;
        $count2++;
        $count++;
      }
      else {
        print
          "Rejecting line contingency $bus1-$bus2 cnt: $bcnt nbus: $buscnt\n";
      }
    }
  }
}

#
# print out N-2 contingencies
#
$ncont = @nminus2;
$ngg   = 0;
$ngl   = 0;
$nll   = 0;
for ( $ii = 0 ; $ii < $ncont ; $ii++ ) {
  for ( $jj = $ii + 1 ; $jj < $ncont ; $jj++ ) {
    $type  = "";
    @buses = ();
    @gtags = ();
    @ltags = ();
    $nbus  = 0;
    $ngtag = 0;
    $nltag = 0;
    if ( $nminus2[$ii] =~ /GENERATOR/ ) {
      $type = "G";
      if ( $nminus2[$ii] =~ /GENERATOR\_(.*)\_(.*)/ ) {
        $buses[$nbus] = $1;
        $nbus++;
        $gtags[$ngtag] = $2;
        $ngtag++;
      }
    }
    elsif ( $nminus2[$ii] =~ /LINE/ ) {
      $type = "L";
      if ( $nminus2[$ii] =~ /LINE\_(.*)\_(.*)\_(.*)/ ) {
        $buses[$nbus] = $1;
        $nbus++;
        $buses[$nbus] = $2;
        $nbus++;
        $ltags[$nltag] = $3;
        $nltag++;
      }
    }
    if ( $nminus2[$jj] =~ /GENERATOR/ ) {
      $type .= "G";
      if ( $nminus2[$jj] =~ /GENERATOR\_(.*)\_(.*)/ ) {
        $buses[$nbus] = $1;
        $nbus++;
        $gtags[$ngtag] = $2;
        $ngtag++;
      }
    }
    elsif ( $nminus2[$jj] =~ /LINE/ ) {
      $type .= "L";
      $contingency = "LINE\_$bus1\_$bus2\_$tag";
      if ( $nminus2[$jj] =~ /LINE\_(.*)\_(.*)\_(.*)/ ) {
        $buses[$nbus] = $1;
        $nbus++;
        $buses[$nbus] = $2;
        $nbus++;
        $ltags[$nltag] = $3;
        $nltag++;
      }
    }
    #
    #  Check for line contingences that split network into two pieces
    #
    $check = 1;
    if ( $type ne "GG" ) {
      $bcnt   = 0;
      @newbus = ();
      @found  = ();
      for ( $ibus = 0 ; $ibus < $buscnt ; $ibus++ ) {
        $found[$ibus] = 0;
      }
      #
      #  Find first bus
      #
      $ifirst = 0;
      $check1 = 0;
      while ( $check1 == 0 ) {
        $ck = 1;
        for ( $i = 0 ; $i < $nbus ; $i++ ) {
          if ( $ifirst + 1 == $buses[$i] ) { $ck = 0; }
        }
        if ( $ck == 0 ) {
          $ifirst++;
        }
        else {
          $check1 = 1;
        }
      }
      $found[$ifirst] = 1;
      $bcnt           = 1;
      $nnbus          = 0;
      $ibus           = $top[$ifirst];
      while ( $ibus >= 0 ) {
        $j              = $nghbrs[$ibus];
        $newbus[$nnbus] = $j;
        $found[$j]      = 1;
        $bcnt++;
        $nnbus++;
        $ibus = $next[$ibus];
      }
      while ( $nnbus > 0 ) {
        @oldbus = ();
        for ( $i = 0 ; $i < $nnbus ; $i++ ) { $oldbus[$i] = $newbus[$i]; }
        $newcnt = 0;
        @newbus = ();
        for ( $i = 0 ; $i < $nnbus ; $i++ ) {
          $j = $top[ $oldbus[$i] ];
          while ( $j >= 0 ) {
            $ob     = $oldbus[$i] + 1;
            $nb     = $nghbrs[$j] + 1;
            $check1 = 0;
            if ( $type eq "GL" ) {
              if (
                !(
                  (
                       $oldbus[$i] + 1 == $buses[1]
                    && $nghbrs[$j] + 1 == $buses[2]
                    && $nghbrtag[$j] == $ltags[0]
                  )
                  || ( $oldbus[$i] + 1 == $buses[2]
                    && $nghbrs[$j] + 1 == $buses[1]
                    && $nghbrtag[$j] == $ltags[0] )
                )
                )
              {
                $check1 = 1;
              }
            }
            elsif ( $type eq "LL" ) {
              $ck1 = 0;
              if (
                !(
                  (
                       $oldbus[$i] + 1 == $buses[0]
                    && $nghbrs[$j] + 1 == $buses[1]
                    && $nghbrtag[$j] == $ltags[0]
                  )
                  || ( $oldbus[$i] + 1 == $buses[1]
                    && $nghbrs[$j] + 1 == $buses[0]
                    && $nghbrtag[$j] == $ltags[0] )
                )
                )
              {
                $ck1 = 1;
              }
              $ck2 = 0;
              if (
                !(
                  (
                       $oldbus[$i] + 1 == $buses[2]
                    && $nghbrs[$j] + 1 == $buses[3]
                    && $nghbrtag[$j] == $ltags[1]
                  )
                  || ( $oldbus[$i] + 1 == $buses[3]
                    && $nghbrs[$j] + 1 == $buses[2]
                    && $nghbrtag[$j] == $ltags[1] )
                )
                )
              {
                $ck2 = 1;
              }
              if ( $ck1 && $ck2 ) {
                $check1 = 1;
              }
            }
            elsif ( $type eq "GG" ) {
              $check1 = 1;
            }
            if ($check1) {
              if ( $found[ $nghbrs[$j] ] == 0 ) {
                $bcnt++;
                $found[ $nghbrs[$j] ] = 1;
                $newbus[$newcnt] = $nghbrs[$j];
                $newcnt++;
              }
            }
            $j = $next[$j];
          }
        }
        $nnbus = $newcnt;
      }
      if ( $bcnt != $buscnt ) {
        $check = 0;
      }
    }
    if ($check) {
      if ( $type eq "GG" ) {
        if ( $allGGC == 1 || $ngg < $maxGGC ) {
          print CONTINGENCY
"CONTINGENCY GG_$buses[0]\_$gtags[0]\_$buses[1]\_$gtags[1]_$count\n";
          print CONTINGENCY "REMOVE UNIT $gtags[0] FROM BUS $buses[0]\n";
          print CONTINGENCY "REMOVE UNIT $gtags[1] FROM BUS $buses[1]\n";
          print CONTINGENCY "END\n";
          $count++;
          $ngg++;
        }
      }
      elsif ( $type eq "GL" ) {
        if ( $allGLC == 1 || $ngl < $maxGLC ) {
          print CONTINGENCY
"CONTINGENCY GL_$buses[0]\_$gtags[0]\_$buses[1]-$buses[2]\_$ltags[0]_$count\n";
          print CONTINGENCY "REMOVE UNIT $gtags[0] FROM BUS $buses[0]\n";
          print CONTINGENCY
            "OPEN BRANCH FROM BUS $buses[1] TO $buses[2] CIRCUT $ltags[0]\n";
          print CONTINGENCY "END\n";
          $count++;
          $ngl++;
        }
      }
      elsif ( $type eq "LL" ) {
        if ( $allLLC == 1 || $nll < $maxLLC ) {
          print CONTINGENCY
"CONTINGENCY LL_$buses[0]-$buses[1]\_$ltags[0]\_$buses[2]-$buses[3]\_$ltags[1]_$count\n";
          print CONTINGENCY
            "OPEN BRANCH FROM BUS $buses[0] TO $buses[1] CIRCUT $ltags[0]\n";
          print CONTINGENCY
            "OPEN BRANCH FROM BUS $buses[2] TO $buses[3] CIRCUT $ltags[1]\n";
          print CONTINGENCY "END\n";
          $count++;
          $nll++;
        }
      }
    }
    else {
      if ( $type eq "GL" ) {
        print
          "Rejecting GL contingency $buses[1]-$buses[2] circuit: $ltags[0]\n";
      }
      if ( $type eq "LL" ) {
        print
"Rejecting LL contingency $buses[0]-$buses[1] circuit: $ltags[0] $buses[2]-$buses[3] circuit: $ltags[1]\n";
      }
    }
  }
}
close(CONTINGENCY);
