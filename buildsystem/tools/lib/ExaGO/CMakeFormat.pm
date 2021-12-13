package ExaGO::CMakeFormat;

use strict;
use warnings;
use v5.16;

use File::Find;
use ExaGO::Env;
use Exporter;
our @ISA    = qw( Exporter );
our @EXPORT = qw(tool);

sub tool {
  my $help    = shift @_;
  my $verbose = shift @_;
  my $inplace = shift @_;

  if ($help) {
    say "cmake-format utility script for ExaGO";
    say "Usage: ./cmake-format.pl";
    say "\t-h: print this help message";
    say "\t-i: perform formatting in-place.";
    say "\t\tOtherwise, the script verifies that files are formatted.";
    say "\t-v: verbose output";
    print
"\nSet the environment variable CMAKEFORMAT to the cmake format executable you ";
    say "would like to use, if you have multiple.";
    say "\nYou most likely just need to run:";
    say "\n\t\$ ./scripts/cmake-format.pl -i\n";
    say "from the root ExaGO source directory before committing your code.";
    exit 1;
  }

  my $tool = "cmake-format";

  # Try to find cmake-format executable
  my $cf = '';

  if ( exists( $ENV{'CMAKEFORMAT'} ) ) {    # Use env var CMAKEFORMAT if found
    $cf = $ENV{'CMAKEFORMAT'};
  }
  elsif ( $host =~ /newell/s ) {
    &module('load', 'python/miniconda3.8');
    $ENV{'PYTHONPATH'} = '/qfs/projects/exasgd/src/cmake_format_newell';
    $cf = '/qfs/projects/exasgd/src/cmake_format_newell/bin/cmake-format';
  }
  else {                                    # else just look in PATH
    $cf = `which cmake-format`;
    chomp($cf);
  }
  `which $cf`;
  if ($?) {
    say
"ERROR: no cmake-format executable could be found. CMake code will not be linted.";
    return 1;
  }

  my $ver = `$cf --version`;
  chomp($ver);

  say "Found cmake-format version $ver at '$cf'";

  my @dirs = (
    "$root/src",                 "$root/include",
    "$root/tests/functionality", "$root/tests/interfaces",
    "$root/tests/unit",          "$root/buildsystem",
    "$root/CMakeLists.txt",
  );

  my @fails;

  find(
    sub {
      if (/\.cmake$|^CMakeLists\.txt$/s) {
        my $f   = $File::Find::name;
        my $cmd = "$cf " . ( $inplace ? "--in-place" : "--check" ) . " $f";
        `$cmd`;
        my $ec = $?;
        if ($inplace) {
          printoneline $verbose, "($tool) Formatting " . ( $verbose ? $f : $_ );
        }
        else {
          printoneline $verbose, "($tool) Checking " . ( $verbose ? $f : $_ );
          if ($ec) { push @fails, $f; }
        }
      }
    },
    @dirs
  );

  my $ret = scalar @fails;
  if ( $ret and ( not $inplace ) ) {
    say "The following files are not formatted:$el";
    say join "\n", @fails;
  }
  say "Returning $ret $el";
  return $ret;
}
