#!/usr/bin/env perl

# Use `perltidy -i=2 -b cmake-format.pl` to format this file after editing. Must
# install the perl module Perl::Tidy to format.
#
# Asher Mancinelli <asher.mancinelli@pnnl.gov>
#
# Apply cmake format to all source files in ExaGO, or check for formatting.

use strict;
use warnings;
use v5.16;

use FindBin qw($Bin);
use lib "$Bin/lib";
use ExaGO;
use Getopt::Std;
$Getopt::Std::STANDARD_HELP_VERSION = 1;
our ( $opt_h, $opt_i, $opt_v );

getopts('ihv');

  if ($opt_h) {
    say "cmake-format utility script for ExaGO";
    say "Usage: $0";
    say "\t-h: print this help message";
    say "\t-i: perform formatting in-place.";
    say "\t\tOtherwise, the script verifies that files are formatted.";
    say "\t-v: verbose output";
    print
"\nSet the environment variable CMAKEFORMAT to the cmake format executable you ";
    say "would like to use, if you have multiple.";
    say "\nYou most likely just need to run:";
    say "\n\t\$ $0 -i\n";
    say "from the root ExaGO source directory before committing your code.";
    exit 1;
  }

my $ret = ExaGO::CMakeFormat::tool( $opt_v, $opt_i );

exit $ret;
