package ExaGO::FileNamingConventions;

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

  if ($help) {
    say "Utility script for checking file naming conventions in ExaGO";
    say "Usage: ./file-nameing-conventions.pl";
    say "\t-h: print this help message";
    say "\t-i: perform formatting in-place.";
    say "\t\tOtherwise, the script verifies that files are formatted.";
    say "\t-v: verbose output";
    print
      "\nSet the environment variable CF to the cmake format executable you ";
    say "would like to use, if you have multiple.";
    say "\nYou most likely just need to run:";
    say "\n\t\$ ./scripts/file-nameing-conventions.pl -i\n";
    say "from the root ExaGO source directory before committing your code.";
    exit 1;
  }

  my $tool = "file-naming-conventions";

  my @dirs = (
    "$root/src",                 "$root/include",
    "$root/tests/functionality", "$root/tests/interfaces",
    "$root/tests/unit", "$root/applications"
  );

  my @fails;

  find(
    sub {
      # Strip prefix to root source directory from filename before performing
      # replacements
      my $f = $File::Find::name;
      $f =~ s{$root/}{}g;
      printoneline $verbose, "($tool) Checking " . ( $verbose ? $f : $_ );

      # Check for dashes where we expect underscores
      if (/-/s) {
        say "P003: found dash '-' where expected underscore '_'. " .
          "Consider the following replacement:$el";
        my $cmd = "mv $f";
        $f =~ s{-}{_}g;
        $cmd = "$cmd"." $f";
        say "\t\$ $cmd\n";
        push @fails, $cmd;
      }

      # Check for naming of source files
      if ( /\.cpp$|\.hpp$|\.h$|\.c$|\.cu$/s and /[A-Z]/s ) {

        # Check for pascal-cased filenames
        if ( $f =~ /[A-Z]/s ) {
          $f =~ s/.+\/..\/..\///g;

          say "P003: possible pascal cased filename. " .
            "Consider the following replacement:$el";
          my $cmd = "mv $f";

          # Replace leading character with lowercase character
          $f =~ s/\/([A-Z])/\/\L$1/g;

          # Replace all other chars with underscore followed by lowercase char
          $f =~ s/([A-Z])/_\L$1/g;

          $cmd = "$cmd "."$f";

          say "\t\$ $cmd\n";
          push @fails, $cmd;
        }
      }
    },
    @dirs
  );

  my $ret = scalar @fails;
  if ($ret) {
    say "These are the suggested changes:$el";
    say join "\n", @fails;
    say "Please visit docs/DeveloperGuidelines.md if you need to review " .
      "conventions.";
  }
  say "Returning $ret $el";

  return $ret;
}
