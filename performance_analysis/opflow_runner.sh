#!/bin/bash

cd $HOME/exago_repo/exago
opflow -netfile datafiles/case_ACTIVSg200.m \
    -opflow_solver HIOP \
    -opflow_model PBPOLRAJAHIOP \
    -hiop_compute_mode GPU \
    -print_output 0 \
    -hiop_verbosity_level 3 \
    -log_view