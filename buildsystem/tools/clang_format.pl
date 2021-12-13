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

my $ret = ExaGO::ClangFormat::tool( $opt_h, $opt_v, $opt_i );

exit $ret;
