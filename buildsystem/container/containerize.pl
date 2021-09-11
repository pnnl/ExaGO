#!/usr/bin/perl
use v5.16;
use strict;
use warnings;
use Getopt::Std;

sub main::HELP_MESSAGE {
  print "-f <file> spack environment to containerize\n";
  print "-c        should this environment generate the build cache as well?\n";

  exit;
}

our($opt_f, $opt_c);
getopts("cf:");

$opt_f = "./spack.yaml" if not defined $opt_f;
die "Could not find spack environment" unless -f $opt_f;

open(FH, '<', $opt_f) or die $!;

# Line continuation
sub bslash {
  my $str = shift;
  print "$str \\\n";
}

# Continue on the next line, contingent on this line's success
sub andand {
  my $str = shift;
  print "$str && \\\n";
}

# End with a posix "true" command
sub endwrap { print "  :\n\n"; }

print "FROM spack/ubuntu-bionic:latest\n\n";

# Create the environment in the container
andand("RUN mkdir /opt/spack-environment");
while (<FH>) {
  chomp();
  andand("  echo \"$_\"");
}
endwrap();

# Base commands used in every build
andand("RUN cd /opt/spack-environment");

# Load secrets to the minio instance
andand("  . /kaniko/s3env.sh");
andand("  set -x");
andand("  spack mirror add minio s3://spack");
andand("  spack mirror add local file:///cache");
andand("  spack mirror add e4s https://cache.e4s.io");
andand("  spack buildcache keys -it");
andand("  spack env activate .");
andand("  spack install --fail-fast");

if ($opt_c) {
  # Generate local buildcache
  andand("  spack gpg init");
  bslash("  spack gpg create 'Asher Mancinelli' 'ashermancinelli\@gmail.com'");
  andand("   --comment 'Key use to rebuild buildcache in CI'");
  andand("  mkdir /cache");
  bslash("  for ii in \$(spack find --format \"yyy {version} /{hash}\" |");
  bslash("        grep -v -E \"^(develop^master)\" |");
  bslash("        grep \"yyy\" |");
  bslash("        cut -f3 -d\" \");");
  bslash("  do");
  bslash("    spack buildcache create -af -d /cache --only=package \$ii ;");
  andand("  done");
  andand("  spack buildcache sync --src-directory /cache --dest-mirror-url s3://spack");

  # Get gpg key id of key used to sign the buildcache packages
  bslash("  keyid=\$(spack gpg list |");
  bslash("    perl -e '");
  bslash("      my \$pub=0;");
  bslash("      while (<>) {");
  bslash("        if (\$pub == 1) { s/\\s+//g; print; exit; };");
  bslash("        if (/^pub/) { \$pub=1; }");
  andand("      }')");

  # Export that key so we can upload it later
  andand("  spack gpg export \"\${keyid}.pub\"");

  # Install minio client
  andand("  curl https://dl.min.io/client/mc/release/linux-amd64/mc -s -o ./mc");
  andand("  chmod +x ./mc");
  andand("  ./mc alias set minio \$S3_ENDPOINT_URL \$AWS_ACCESS_KEY_ID \$AWS_SECRET_ACCESS_KEY");
  andand("  ./mc cp \"\${keyid}.pub\" minio/spack/build_cache/_pgp/");
  andand("  (./mc cp minio/spack/build_cache/_pgp/index.json . || echo '{\"keys\": {}}' >> index.json)");

  # Install JQ to update spack buildcache keys index.json
  andand("  curl https://github.com/stedolan/jq/releases/download/jq-1.6/jq-linux64 -L -s -o jq");
  andand("  chmod +x ./jq");
  andand("  (cat index.json | ./jq --arg keyid \"\${keyid}\" '. * { keys: { (\$keyid): {} } }' > updated_index.json)");
  andand("  ./mc cp updated_index.json minio/spack/build_cache/_pgp/index.json");
}

endwrap();