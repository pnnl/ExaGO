#!/usr/bin/env perl
$argc = @ARGV;
if (argc > 0) {
  $directory = $ARGV[0];
  $cmd = "rm $directory/*.warning";
  system($cmd);
}
