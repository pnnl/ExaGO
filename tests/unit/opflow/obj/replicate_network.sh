#!/bin/bash

# Authors:
#   - Cameron Rutherford <cameron.rutherford@pnnl.gov>

# Returning early if no network size if given
if [ $# -eq 0 ]; then
  echo "No network size given. Call this script with a desired network size."
  exit 1
fi

# Initial Network
initial_net=OF-unittest1.m
network_size=$1
output_net=OF-unittestx$network_size\_temp.m

# Patch of replicating the network once
diff OF-unittest1.m OF-unittestx2.m -ed > patch.txt

# Replicate the patch to be of desired size


# Apply the patch after modification
patch --ed --output=$output_net OF-unittest1.m patch.txt

# rm patch.txt

