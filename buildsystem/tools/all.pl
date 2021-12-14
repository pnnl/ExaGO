#!/usr/bin/env perl

# See buildsystem/tools/lib/ExaGO.pm for documentation

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
  say "Run all ExaGO code quality tools";
  say "Usage: ./all.pl [options]";
  say "\t-h: show this help message";
  say "\t-v: show verbose output";
  say "\t-i: allow tools to perform edits to source code in-place";
  exit 1;
}

my $ret = 0;

say "-- Running clang-format tool";
$ret += ExaGO::ClangFormat::tool( $opt_v, $opt_i );

say "-- Running cmake-format tool";
$ret += ExaGO::CMakeFormat::tool( $opt_v, $opt_i );

say "-- Running file naming convention tool";
$ret += ExaGO::FileNamingConventions::tool( $opt_v );

say "-- Running python formatting tool";
$ret += ExaGO::PythonFormat::tool( $opt_v, $opt_i );

exit $ret;
