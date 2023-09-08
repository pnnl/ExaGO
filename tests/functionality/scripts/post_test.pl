#!/usr/bin/env perl

use Term::ANSIColor;

$argc = $ARGV[0];

sub get_sorted_files {
  my $path = shift;
  opendir my($dirh), $path or die "can't opendir $path: $!";
  my %hash = map { $_ => (stat($_))[9] || undef} # avoid empty list
              map { "$path/$_"} # Need full paths for sorting
              readdir $dirh;
  closedir $dirh;
  return %hash;
}

my %files = get_sorted_files($argc);
@filekeys = keys %files;


$warnings = 0;
foreach my $key (sort {$files{$a} <=> $files{$b}} @filekeys) {
  if ($key =~ /\.warning\s*$/) {
    $warnings++;
  }
}
if ($warnings > 0) {
  print "The following tests generated warnings:\n";
}
foreach my $key (sort {$files{$a} <=> $files{$b}} @filekeys) {
#  print "FOUND FILE $key\n";
  if ($key =~ /\.warning\s*$/) {
    open (FILE,$key);
    $cnt = 0;
    while (<FILE>) {
#      print "$_";
      if ($cnt == 1) {
        $line = $_;
#        $line =~ s/:\s*$//;
#        if ($line =~ /FUNCTIONALITY/) {
          printf "         %s\n", colored($line, 'yellow');
#        }
      }
      $cnt++;
    }
    close(FILE);
#    system("rm $key");
  }
}
if ($warnings > 0) {
  print "Total tests generating warnings: $warnings\n";
  print "More information can be found on warnings by\n";
  print "looking at .warning files located in directory\n";
  printf "   %s\n", colored($argc, 'magenta');
}
