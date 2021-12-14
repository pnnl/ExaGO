package ExaGO::ClangFormat;

use strict;
use warnings;
use v5.16;

use File::Find;
use ExaGO::Env;
use Exporter;
our @ISA    = qw( Exporter );
our @EXPORT = qw(tool);

# Files that match this pattern will be skipped
my @skip = ('IpoptAdapter.hpp');

sub tool {
  my $verbose = shift @_;
  my $inplace = shift @_;

  my $tool = "clang-format";

  # Try to find clang-format executable
  my $cf = '';

  if ( exists( $ENV{'CLANGFORMAT'} ) ) {    # Use env var CLANGFORMAT if found
    $cf = $ENV{'CLANGFORMAT'};
  }
  elsif ( $host =~ /newell/s ) {
    $cf =
'/qfs/projects/exasgd/newell/clang+llvm-10.0.1-powerpc64le-linux-rhel-7.4/bin/clang-format';
  }
  else {                                    # else just look in PATH
    $cf = `which clang-format`;
    chomp($cf);
  }
  `which $cf`;
  if ($?) {
    say "ERROR: no clang-format executable could be found. C++ source will not "
      . "be linted.";
    return 1;
  }

  my $ver = ( split " ", `$cf --version` )[-1];
  if ( $ver !~ /^10\./ ) {
    say
"ERROR: this script requires clang-format with major version 10, but got $ver
Please navigate to https://releases.llvm.org/download.html, download clang 10, 
and set the environment variable 'CLANGFORMAT=/path/to/clang-format' to the 
appropriate clang-format executable";
    exit 1;
  }

  say "Found clang-format version $ver at '$cf'";

  my @dirs = (
    "$root/src",              "$root/include",
    "$root/interfaces",       "$root/tests/functionality",
    "$root/tests/interfaces", "$root/tests/unit"
  );

  my @fails;

  find(
    sub {
      foreach my $skip_ (@skip) {
        if (/$skip_/s) {
          return;
        }
      }
      if (/\.cpp$|\.hpp$|\.h$|\.c$|\.cu$/s) {
        my $f = $File::Find::name;
        my $cmd =
          "$cf -style=file "
          . ( $inplace ? "-i" : "--dry-run --Werror" ) . " $f";
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

1;
