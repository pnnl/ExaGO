package ExaGO::PythonFormat;

use strict;
use warnings;
use v5.16;

use File::Find;
use ExaGO::Env;
use Exporter;
our @ISA    = qw( Exporter );
our @EXPORT = qw(tool);

# Files that match this pattern will be skipped
my @skip = ('pycache');

sub tool {
  my $verbose = shift @_;
  my $inplace = shift @_;

  my $pycodestyle = '';
  my $autopep8    = '';

  # First check environment variables for appropriate paths
  if ( exists( $ENV{'PYCODESTYLE'} ) ) {
    $pycodestyle = $ENV{'PYCODESTYLE'};
  }
  if ( exists( $ENV{'AUTOPEP8'} ) ) {
    $autopep8 = $ENV{'AUTOPEP8'};
  }

  # If we're on a supported system, the required programs are already installed
  if ( $host =~ /newell/s ) {
    &module( 'load', 'python/miniconda3.8' );
    $ENV{'PYTHONPATH'} =
'/qfs/projects/exasgd/src/pycodestyle_newell:/qfs/projects/exasgd/src/autopep8_newell';
    $autopep8 = '/qfs/projects/exasgd/src/autopep8_newell/bin/autopep8';
    $pycodestyle =
      '/qfs/projects/exasgd/src/pycodestyle_newell/bin/pycodestyle';
  }

  # If either program is not found, we'll just search on the PATH
  if ( not $pycodestyle ) {
    $pycodestyle = `which pycodestyle`;
    chomp($pycodestyle);
  }
  if ( not $autopep8 ) {
    $autopep8 = `which autopep8`;
    chomp($autopep8);
  }

  # Make sure our programs have really been found
  `which $pycodestyle`;
  if ($?) {
    say
"ERROR: no pycodestyle executable could be found. Python source will not be linted.";
    return 1;
  }
  `which $autopep8`;
  if ($?) {
    say
"ERROR: no autopep8 executable could be found. Python source will not be linted.";
  }

  my $pcsver = `$pycodestyle --version`;
  my $ap8ver = `$autopep8 --version`;
  chomp($pcsver);
  chomp($ap8ver);

  say "Found pycodestyle version $pcsver at $pycodestyle";
  say "Found autopep8 version $ap8ver at $autopep8";

  my $tool = "python-format";

  my @dirs = ( "$root/tests/interfaces", "$root/interfaces" );

  my @fails;

  find(
    sub {
      if (/\.py$/s) {

        # Strip prefix to root source directory from filename before performing
        # replacements
        my $f        = $File::Find::name;
        my $fullpath = $f;
        $f =~ s{$root/}{}g;
        foreach my $skip_ (@skip) {
          if ( $f =~ /$skip_/ ) {
            return;
          }
        }
        my $cmd =
          ( $inplace ? "$autopep8 --in-place" : "$pycodestyle" ) . " $fullpath";
        my $res = `$cmd`;
        my $ec = $?;
        if ($inplace) {
          printoneline $verbose, "($tool) Formatting " . ( $verbose ? $f : $_ );
        }
        else {
          printoneline $verbose, "($tool) Checking " . ( $verbose ? $f : $_ );
          if ($ec) {
            say "\n$res";
            push @fails, $f; 
          }
        }
      }
    },
    @dirs
  );

  my $ret = scalar @fails;
  if ($ret) {
    say "The following files are not formatted:$el";
    say join "\n", @fails;
  }
  say "Returning $ret $el";

  return $ret;
}
