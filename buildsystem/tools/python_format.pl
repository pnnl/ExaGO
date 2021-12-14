#!/usr/bin/env perl

# Use `perltidy -i=2 -b` on this file to format this file after editing. Must
# install the perl module Perl::Tidy to format.
#
# Asher Mancinelli <asher.mancinelli@pnnl.gov>
#
# Apply clang format to all source files in ExaGO, or check for formatting.

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
    say "clang-format utility script for ExaGO
Usage: $0
\t-h: print this help message
\t-i: perform formatting in-place.
\t\tOtherwise, the script verifies that files are formatted.
\t-v: verbose output
\nSet the environment variable AUTOPEP8 to the clang format executable you
would like to use, if you have multiple.
\nYou most likely just need to run:
\n\t\$ $0 -i\n
from the root ExaGO source directory before committing your code.";
    return 1;
  }

my $ret = ExaGO::PythonFormat::tool( $opt_v, $opt_i );

exit $ret;
