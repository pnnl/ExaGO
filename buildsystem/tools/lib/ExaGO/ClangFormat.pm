package ExaGO::ClangFormat;

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
    say "clang-format utility script for ExaGO";
    say "Usage: ./clang-format.pl";
    say "\t-h: print this help message";
    say "\t-i: perform formatting in-place.";
    say "\t\tOtherwise, the script verifies that files are formatted.";
    say "\t-v: verbose output";
    print
"\nSet the environment variable CLANGFORMAT to the clang format executable you ";
    say "would like to use, if you have multiple.";
    say "\nYou most likely just need to run:";
    say "\n\t\$ ./scripts/clang-format.pl -i\n";
    say "from the root ExaGO source directory before committing your code.";
    return 1;
  }

  my $tool = "clang-format";

  # Try to find clang-format executable
  my $cf = '';

  if ( exists( $ENV{'CLANGFORMAT'} ) ) {    # Use env var CLANGFORMAT if found
    $cf = $ENV{'CLANGFORMAT'};
  }
  elsif ( $host =~ /newell/s )
  {    # Use full path on newell if no CLANGFORMAT env var
    $cf = '/share/apps/llvm/12.0.0/newell/bin/clang-format';
  }
  else {    # else just look in PATH
    $cf = `which clang-format`;
    chomp($cf);
  }
  `which $cf`;
  if ($?) {
    say
"ERROR: no clang-format executable could be found. C++ source will not be linted.";
    return 1;
  }

  my $ver = ( split " ", `$cf --version` )[-1];

  say "Found clang-format version $ver at '$cf'";

  my @dirs = (
    "$root/src",                 "$root/include",
    "$root/tests/functionality", "$root/tests/interfaces",
    "$root/tests/unit"
  );

  my @fails;

  find(
    sub {
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
