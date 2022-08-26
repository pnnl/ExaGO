
# DL partitions
my @partitions = grep /newell/, `sinfo`;

# Idle nodes
#my @idle_parts = grep /idle/, @partitions;
#
#@idle_parts = reverse @idle_parts;
#
#if ( scalar @idle_parts eq 0 ) {
#
#  # No idle parititons... Just choose the shared one.
#  print "newell_shared";
#}
#else {
#  my @parts = split / /, $idle_parts[0];
#  print $parts[0];
#}

#force newell_shared to use the upgraded node
print "newell_shared"
