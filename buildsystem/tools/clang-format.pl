#!/usr/bin/env perl

# Use `perltidy -i=2 -b clang-format.pl` to format this file after editing. Must
# install the perl module Perl::Tidy to format.
#
# Asher Mancinelli <asher.mancinelli@pnnl.gov>
#
# Apply clang format to all source files in ExaGO, or check for formatting.

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
  say "clang-format utility script for ExaGO";
  say "Usage: ./clang-format.pl";
  say "\t-h: print this help message";
  say "\t-i: perform formatting in-place.";
  say "\t\tOtherwise, the script verifies that files are formatted.";
  say "\t-v: verbose output";
  print "\nSet the environment variable CF to the clang format executable you ";
  say "would like to use, if you have multiple.";
  say "\nYou most likely just need to run:";
  say "\n\t\$ ./scripts/clang-format.pl -i\n";
  say "from the root ExaGO source directory before committing your code.";
  exit 1;
}

my $root = "$Bin/../..";
my $host = hostname;

# Try to find clang-format executable
my $cf = '';

if ( exists( $ENV{'CF'} ) ) {    # Use env var CF if found
  $cf = $ENV{'CF'};
}
elsif ( $host =~ /newell/s ) {    # Use full path on newell if no CF env var
  $cf = '/share/apps/llvm/12.0.0/newell/bin/clang-format';
}
else {    # else just look in PATH
  $cf = `which clang-format`;
  chomp($cf);
}
`which $cf`;
die "No clang-format executable could be found" unless 0 == $?;

my $ver = ( split " ", `$cf --version` )[-1];

say "Found clang-format version $ver at '$cf'";

my @dirs = (
  "$root/src",                 "$root/include",
  "$root/tests/functionality", "$root/tests/interfaces",
  "$root/tests/unit"
);

find( \&wanted, @dirs );
my @fails;

# If verbose, show all files. Otherwise, keep on one line.
sub printoneline {
  my $s = join " ", @_;
  print "$s" . ( $opt_v ? "\n" : "\c[[K\r" );
}

sub wanted {
  if (/\.cpp$|\.hpp$|\.h$|\.c$|\.cu$/s) {
    my $f = $File::Find::name;
    my $cmd =
      "$cf -style=file " . ( $opt_i ? "-i" : "--dry-run --Werror" ) . " $f";
    `$cmd`;
    my $ec = $?;
    if ($opt_i) {
      printoneline "(clang-format) Formatting " . ( $opt_v ? $f : $_ );
    }
    else {
      printoneline "(clang-format) Checking " . ( $opt_v ? $f : $_ );
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
