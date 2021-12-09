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
use File::Find;
use Sys::Hostname;
use Getopt::Std;
$Getopt::Std::STANDARD_HELP_VERSION = 1;
our ( $opt_h, $opt_i, $opt_v );

getopts('ihv');

if ($opt_h) {
  say "cmake-format utility script for ExaGO";
  say "Usage: ./cmake-format.pl";
  say "\t-h: print this help message";
  say "\t-i: perform formatting in-place.";
  say "\t\tOtherwise, the script verifies that files are formatted.";
  say "\t-v: verbose output";
  print "\nSet the environment variable CF to the cmake format executable you ";
  say "would like to use, if you have multiple.";
  say "\nYou most likely just need to run:";
  say "\n\t\$ ./scripts/cmake-format.pl -i\n";
  say "from the root ExaGO source directory before committing your code.";
  exit 1;
}

my $root = "$Bin/../..";
my $host = hostname;

# Try to find cmake-format executable
my $cf = '';

if ( exists( $ENV{'CF'} ) ) {    # Use env var CF if found
  $cf = $ENV{'CF'};
}
else {    # else just look in PATH
  $cf = `which cmake-format`;
  chomp($cf);
}
`which $cf`;
die "No cmake-format executable could be found" unless 0 == $?;

my $ver = `$cf --version`;
chomp($ver);

say "Found cmake-format version $ver at '$cf'";

my @dirs = (
  "$root/src",                 "$root/include",
  "$root/tests/functionality", "$root/tests/interfaces",
  "$root/tests/unit",          "$root/buildsystem",
  "$root/CMakeLists.txt",
);

find( \&wanted, @dirs );
my @fails;

# If verbose, show all files. Otherwise, keep on one line.
sub printoneline {
  my $s = join " ", @_;
  print "$s" . ( $opt_v ? "\n" : "\c[[K\r" );
}

sub wanted {
  if (/\.cmake$|^CMakeLists\.txt$/s) {
    my $f   = $File::Find::name;
    my $cmd = "$cf " . ( $opt_i ? "--in-place" : "--check" ) . " $f";
    `$cmd`;
    my $ec = $?;
    if ($opt_i) {
      printoneline "(cmake-format) Formatting " . ( $opt_v ? $f : $_ );
    }
    else {
      printoneline "(cmake-format) Checking " . ( $opt_v ? $f : $_ );
      if ($ec) { push @fails, $f; }
    }
  }
}

my $ret = scalar @fails;
if ( $ret and ( not $opt_i ) ) {
  say "The following files are not formatted:\c[[K";
  say join "\n", @fails;
}
say "Returning $ret\c[[K";
exit $ret;
