package ExaGO::Env;

use strict;
use warnings;
use v5.16;

use Exporter;
our @ISA    = qw( Exporter );
our @EXPORT = qw($root $host $el printoneline module);

use Cwd 'abs_path';
use Sys::Hostname;
use File::Basename;
our ( $root, $host, $el );

my ( $name, $path, $suffix ) = fileparse(__FILE__);

if ( -e '/usr/share/Modules/init/perl.pm' ) {
  require '/usr/share/Modules/init/perl.pm';
}

# Root of ExaGO source tree
$root = abs_path("$path/../../../..");

# Host name used to load system-specific varaibles
$host = hostname;

# End-of-line character that clears all other output
$el = "\c[[K";

# If verbose, show all files. Otherwise, keep on one line.
sub printoneline {
  my $verbose = shift @_;
  my $s       = join " ", @_;
  print "$s" . ( $verbose ? "\n" : "$el\r" );
}
