#!/bin/bash

# Authors:
#   - Cameron Rutherford <cameron.rutherford@pnnl.gov>

# Returning early if no network size if given
if [ $# -eq 0 ]; then
  echo "No network size given. Call this script with a desired network size."
  exit 1
fi

network_size=$1

#   Replicate the network in this dir, and output file

