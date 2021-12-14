#!/usr/bin/env perl

# Use `perltidy -i=2 -b file-naming-conventions.pl` to format this file after editing. Must
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

getopts('hv');

if ($opt_h) {
  say "Utility script for checking file naming conventions in ExaGO";
  say "Usage: $0";
  say "\t-h: print this help message";
  say "\t-v: verbose output";
  say "\nYou most likely just need to run:";
  say "\n\t\$ $0\n";
  say "from the root ExaGO source directory before committing your code.";
  exit 1;
}

my $ret = ExaGO::FileNamingConventions::tool( $opt_v );

exit $ret;
