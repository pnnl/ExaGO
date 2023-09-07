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
  my $verbose = shift @_;
  my $inplace = shift @_;

  my $tool = "cmake-format";

  # Try to find cmake-format executable
  my $cf = '';

  if ( exists( $ENV{'CMAKEFORMAT'} ) ) {    # Use env var CMAKEFORMAT if found
    $cf = $ENV{'CMAKEFORMAT'};
  }
  elsif ( $host =~ /newell/s ) {
    &module( 'load', 'python/miniconda3.8' );
    $ENV{'PYTHONPATH'} = '/qfs/projects/exasgd/src/cmake_format_newell';
    $cf = '/qfs/projects/exasgd/src/cmake_format_newell/bin/cmake-format';
  }
  elsif ( $host =~ /marianas/s || $host =~ /deception/s ) {
    &module( 'load', 'python/miniconda3.8' );
    $ENV{'PYTHONPATH'} = '/qfs/projects/exasgd/marianas/format-tools';
    $cf = '/qfs/projects/exasgd/marianas/format-tools/cmake-format';
  }
  else {
    $cf = `which cmake-format`;
    chomp($cf);
  }
  `which $cf`;
  if ($?) {
    say "ERROR: no cmake-format executable could be found. CMake code will not "
      . "be linted.";
    return 1;
  }

  my $ver = `$cf --version`;
  chomp($ver);

  my $config = `$cf --dump-config`;
  if ($verbose) {
    say "$config";
  }

  say "Found cmake-format version $ver at '$cf'";

  my @dirs = (
    "$root/src",                 "$root/include",
    "$root/tests/functionality", "$root/tests/interfaces",
    "$root/tests/unit",          "$root/buildsystem",
    "$root/CMakeLists.txt",      "$root/interfaces",
    "$root/applications",        "$root/interfaces/python"
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
